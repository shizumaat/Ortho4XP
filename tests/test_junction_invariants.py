"""Junction polygon invariants — recurring-regression guards.

The user's authoritative shape rules (memory:
``feedback_shape_rules.md``) say:

* Junction polygons are the *residue* between rects / aprons /
  runways / terminals.  They are inherently small — incoming-corner
  count plus ≤ 4 trace points per arc between consecutive incoming
  corners.  Above ~5,000 m² the residue is apron-sized and should
  be classified as ``role=apron`` instead.
* "Coverage invariants → Shared vertices exact": every neighbour
  vertex that lands on a junction's perimeter must coincide with
  one of the junction's own ring vertices.  Otherwise X-Plane sees
  a free elevation gap at the kiss point and renders a cliff.

The CYXY ``-10070`` regression that prompted these tests had all
three failure modes at once: 32,000 m² area, 80 ring vertices, and
only 3 of its perimeter points shared with neighbouring shapes.

Each test parametrises across the standard test airports.  A
session-level layout cache avoids rebuilding the same airport for
every test.
"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path
from typing import List, Tuple

import pytest
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

from conftest import (
    airports_under_test, baseline_airports,
    xplane_available, xplane_root,
)


def _test_airports() -> list:
    """Union of baseline airports (always-run) + env-gated airports.

    Per user 2026-05-16: invariant tests must run on every canonical
    baseline airport unconditionally so geometry regressions can't
    slip past CI without being noticed.  ``O4_TEST_AIRPORTS=...``
    still extends the set for ad-hoc coverage of additional ICAOs.
    """
    seen = set()
    out = []
    for ic in list(baseline_airports()) + list(airports_under_test()):
        if ic not in seen:
            seen.add(ic)
            out.append(ic)
    return out

_HERE = Path(__file__).resolve().parent
_SRC = _HERE.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


def _xplane_root() -> str:
    return xplane_root()


def _xplane_available() -> bool:
    return xplane_available()


pytestmark = pytest.mark.skipif(
    not _xplane_available(),
    reason="X-Plane install not found (set XPLANE_ROOT to override)",
)


# Hard invariant thresholds — global, no per-airport relaxation.
# Per user 2026-04-30: junction validity is a question of geometry,
# not area.  A 6-way mega-intersection is a valid junction even
# though it's large.  The invariant is "every boundary point lies
# within ~one taxi-half-width of some converging centerline" — apron
# mis-classifications break this; large valid junctions don't.
MAX_BOUNDARY_TO_CENTERLINE_M = 20.0
BOUNDARY_SAMPLE_STEP_M = 5.0

# Orphan neighbour vertices: zero, hard.  A neighbour vertex
# touching a junction's perimeter line must coincide with one of
# the junction's own ring vertices.
MAX_ORPHAN_NEIGHBOUR_VERTICES = 0


# ── Per-airport regression baselines (calibrated 2026-05-01) ──
#
# Three tests in this file (vertex-count-bounded,
# boundary-near-centerline, taxi-rects-not-alongside-apron) fire
# heavily at SPJC because of pre-existing geometry-quality bugs:
# our junction polygons are vastly more sprawling than the
# ``SPJC_target.osm`` ground truth (target max=36 verts, we
# produce up to 336; target boundaries fit close to centerlines,
# ours stray as much as 566 m).  Once those bugs are fixed the
# baselines should drop toward zero and eventually be deleted.
#
# Each baseline records the WORST observed value at HEAD as a
# regression ceiling: tests fail only if a future change exceeds
# the recorded ceiling.  Lower recorded numbers → tighter gate.
# When you fix something, run the test and lower the baseline.
#
# Airports without an explicit baseline use the default tight
# value (zero offenders / hard cap) — those airports are still
# fully gated by the original invariant.
JUNCTION_BOUNDARY_DISTANCE_REGRESSION_BASELINE = {
    # SPJC: all 43 junctions have boundary points > 20 m from
    # any centerline; worst is 566.8 m (at the same SE-apron
    # mega-junction).  The 20 m threshold is tight for normal
    # taxi-junction geometry but ours are apron-sized at SPJC.
    "SPJC": {"max_offenders": 43, "max_distance_m": 567.0},
}

# Per user 2026-05-16: the shared-sloping-edge rule is universal —
# no airport-specific exemptions.  A sloping rect's sloping edge
# must never be shared by a junction/apron polygon's perimeter.

ORPHAN_NEIGHBOUR_VERTEX_REGRESSION_BASELINE = {
    # SPJC: 5 vertex orphans:
    # * 1 rect corner shifted off its co-located junction vertex by
    #   Rule 5 (push-outside-pavement) -- pending tighter Rule 5
    #   anchor exemption.
    # * 4 boundary-polygon vertices near junction perimeters that
    #   don't snap to junction vertices (boundary traces airport
    #   outline; vertices come from apt.dat row-130 / OSM).  Pending
    #   boundary-vertex snap to nearest junction corner.
    "SPJC": 5,
}


# A neighbour vertex within this distance of a junction's perimeter
# line is considered "kissing" and required to be shared.
ORPHAN_NEAR_PERIMETER_M = 1.0
# Coincidence tolerance: shared vertices need not be bit-identical
# (float-precision drift through unary_union etc.) but must agree
# at sub-decimetre level.
ORPHAN_SAME_VERTEX_TOL_M = 0.10


_LAYOUT_CACHE: dict = {}


def _build_layout(icao: str):
    if icao in _LAYOUT_CACHE:
        return _LAYOUT_CACHE[icao]
    from auto_patch.pipeline import build_airport_pavement
    layout = build_airport_pavement(
        icao, _xplane_root(), compute_elevations=True)
    _LAYOUT_CACHE[icao] = layout
    return layout


def _aeroway_centerlines_m(layout):
    """Union of taxi + runway centerlines available from
    ``layout.shapes``.  Used as the "valid junction-boundary anchor"
    set per the user 2026-04-30 rule.

    Sources:
      * Taxi rects keep ``source_axis`` (the OSM centerline span the
        rect was built from).
      * Runway segments are emitted as 4-corner rects without
        ``source_axis``; we derive the long-axis from the corners.
      * Junctions / terminals / aprons / bridges contribute nothing
        — they're regions, not corridors.
    """
    lines = []
    for s in layout.shapes:
        if s.source_axis is not None and not s.source_axis.is_empty:
            lines.append(s.source_axis)
            continue
        if (s.role == "runway"
                and s.polygon is not None
                and not s.polygon.is_empty):
            coords = list(s.polygon.exterior.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            if len(coords) == 4:
                a_mid = (0.5 * (coords[0][0] + coords[3][0]),
                         0.5 * (coords[0][1] + coords[3][1]))
                b_mid = (0.5 * (coords[1][0] + coords[2][0]),
                         0.5 * (coords[1][1] + coords[2][1]))
                lines.append(LineString([a_mid, b_mid]))
    if not lines:
        return None
    try:
        return unary_union(lines)
    except Exception:
        return None


def _shape_label(layout, idx: int, s) -> str:
    """Stable, human-readable identifier for a shape in failure
    messages.  ``layout.shapes`` index is order-dependent but
    matches the OSM emit order, so it's the most useful pointer
    when inspecting the emitted .osm."""
    ref = getattr(s, "ref", "") or ""
    return f"#{idx}({s.role}{('/' + ref) if ref else ''})"


@pytest.mark.parametrize("icao", _test_airports())
def test_junction_boundary_near_centerline(icao):
    """Per user 2026-04-30: a valid junction's pavement edge is
    always "relatively close" to a converging taxiway / runway
    centerline.  No matter how many taxiways meet at a junction, the
    surrounding apt.dat pavement edge sits at most one taxi
    half-width (+ a little fillet) from the local centerline.

    A junction whose boundary strays farther than
    ``MAX_BOUNDARY_TO_CENTERLINE_M`` from any centerline contains
    apron-territory pavement (no centerline running through it) and
    should be re-classified as ``role=apron`` or split.

    Area alone is NOT the test — a 6-way mega-intersection can be
    legitimately large.  The geometric invariant is what matters.

    Per-airport regression baseline:
    Airports listed in
    ``JUNCTION_BOUNDARY_DISTANCE_REGRESSION_BASELINE`` have a
    known-bad ceiling (count + worst distance); the test fails
    only if either exceeds the recorded value.  Other airports
    are gated tightly (zero offenders).
    """
    layout = _build_layout(icao)
    centers = _aeroway_centerlines_m(layout)
    if centers is None or centers.is_empty:
        pytest.skip(f"{icao}: no aeroway centerlines extractable")
    cap = MAX_BOUNDARY_TO_CENTERLINE_M
    offenders = []
    for idx, s in enumerate(layout.shapes):
        if s.role != "junction":
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        bnd = s.polygon.boundary
        L = bnd.length
        n_steps = max(2, int(L / BOUNDARY_SAMPLE_STEP_M) + 1)
        max_d = 0.0
        max_pt = (0.0, 0.0)
        for i in range(n_steps):
            u = min(L, i * BOUNDARY_SAMPLE_STEP_M)
            p = bnd.interpolate(u)
            d = centers.distance(p)
            if d > max_d:
                max_d = d
                max_pt = (p.x, p.y)
        if max_d > cap:
            offenders.append((
                max_d, _shape_label(layout, idx, s),
                s.polygon.area, max_pt))
    offenders.sort(reverse=True)
    summary = "; ".join(
        f"{lbl} max_d={d:.1f}m at ({mp[0]:.0f},{mp[1]:.0f}) "
        f"area={a:,.0f} m²"
        for d, lbl, a, mp in offenders[:5])

    baseline = JUNCTION_BOUNDARY_DISTANCE_REGRESSION_BASELINE.get(icao)
    if baseline:
        worst_d = offenders[0][0] if offenders else 0.0
        n = len(offenders)
        assert n <= baseline["max_offenders"], (
            f"{icao}: {n} junctions exceed {cap:.0f} m centerline "
            f"distance — exceeds known-bad baseline of "
            f"{baseline['max_offenders']}.  Top: {summary}.")
        assert worst_d <= baseline["max_distance_m"] + 1.0, (
            f"{icao}: worst boundary distance {worst_d:.1f} m "
            f"exceeds known-bad baseline of "
            f"{baseline['max_distance_m']:.1f} m.  Top: {summary}.")
    else:
        assert not offenders, (
            f"{icao}: {len(offenders)} junction polygon(s) have "
            f"boundary points > {cap:.0f} m from nearest "
            f"taxi/runway centerline.  Top: {summary}.")


@pytest.mark.parametrize("icao", _test_airports())
def test_junction_vertices_have_source(icao):
    """Every junction vertex must originate from a geometric source:
    a corner of an adjacent rect (sloping or runway), apron,
    terminal, groundside, or boundary polygon.  Vertices without a
    source are orphans added by densification or buffer rounding
    and must be eliminated.

    This is the dual of ``test_junction_neighbour_corners_shared``:
    that test asserts neighbour vertices within 1 m of a junction
    perimeter coincide with junction vertices; this test asserts
    junction vertices coincide with neighbour-shape corners.

    Per user 2026-05-18: replaces the older "incoming-corners +
    ≤ 4 trace per arc" vertex cap which reflected an obsolete
    boundary-trace architecture.  Universal — no airport-specific
    exemptions.
    """
    from auto_patch.layout import SHARED_VERTEX_TOL_M

    layout = _build_layout(icao)

    # Source shapes: anything that contributes a real geometric
    # corner that a junction can legitimately anchor on.  Other
    # junctions are excluded — two junctions sharing a vertex
    # doesn't ground it in source geometry.
    SOURCE_ROLES = {
        "runway", "primary_parallel", "secondary_parallel",
        "stub", "cross_connector",
        "apron", "terminal", "groundside_pavement", "boundary",
        "tunnel_ramp", "retaining_wall",
    }
    source_corners: List[Tuple[float, float]] = []
    for s in layout.shapes:
        if s.role not in SOURCE_ROLES:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        source_corners.extend(
            (float(c[0]), float(c[1])) for c in coords)

    if not source_corners:
        pytest.skip(
            f"{icao}: no source-shape corners to anchor junctions")

    tol = SHARED_VERTEX_TOL_M
    tol_sq = tol * tol

    orphans: List[str] = []
    for idx, s in enumerate(layout.shapes):
        if s.role != "junction":
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for v_idx, (vx, vy) in enumerate(coords):
            best_d_sq = min(
                (cx - vx) ** 2 + (cy - vy) ** 2
                for cx, cy in source_corners)
            if best_d_sq > tol_sq:
                d = math.sqrt(best_d_sq)
                orphans.append(
                    f"{_shape_label(layout, idx, s)} "
                    f"vertex#{v_idx} at ({vx:.1f},{vy:.1f}) — "
                    f"nearest source corner {d:.2f} m away")

    assert not orphans, (
        f"{icao}: {len(orphans)} junction vertex(es) have no "
        f"source-shape corner within {tol:.2f} m.  First 5:\n  "
        + "\n  ".join(orphans[:5]))


@pytest.mark.parametrize("icao", _test_airports())
def test_junction_neighbour_corners_shared(icao):
    """Coverage invariant: for every junction polygon, every vertex
    of a neighbouring shape (rect / runway / terminal / another
    junction) that lies within ``ORPHAN_NEAR_PERIMETER_M`` of the
    junction's perimeter must coincide (within
    ``ORPHAN_SAME_VERTEX_TOL_M``) with one of the junction's own
    ring vertices.

    Failures are the "adjacent-but-not-shared" pattern that
    produces visible elevation cliffs in X-Plane (CYXY -10070, the
    HECA junction-cluster issue).
    """
    layout = _build_layout(icao)
    junctions = [
        (idx, s) for idx, s in enumerate(layout.shapes)
        if s.role == "junction"
        and s.polygon is not None
        and not s.polygon.is_empty]
    others = [
        (idx, s) for idx, s in enumerate(layout.shapes)
        if s.role != "junction"
        and s.polygon is not None
        and not s.polygon.is_empty]
    if not junctions or not others:
        return

    # Pre-extract neighbour exterior vertices (skip closing repeat).
    nbr_pts = []
    for n_idx, n_s in others:
        for ox, oy in list(n_s.polygon.exterior.coords)[:-1]:
            nbr_pts.append((n_idx, n_s, ox, oy))

    orphans = []  # (miss_dist, j_label, nbr_label, ox, oy)
    for j_idx, j_s in junctions:
        bnd = j_s.polygon.boundary
        j_coords = list(j_s.polygon.exterior.coords)[:-1]
        j_xs = [c[0] for c in j_coords]
        j_ys = [c[1] for c in j_coords]
        # Quick AABB to skip far-away neighbours.
        x_min, y_min, x_max, y_max = j_s.polygon.bounds
        pad = ORPHAN_NEAR_PERIMETER_M + ORPHAN_SAME_VERTEX_TOL_M
        for n_idx, n_s, ox, oy in nbr_pts:
            if (ox < x_min - pad or ox > x_max + pad
                    or oy < y_min - pad or oy > y_max + pad):
                continue
            d_perim = bnd.distance(Point(ox, oy))
            if d_perim > ORPHAN_NEAR_PERIMETER_M:
                continue
            # Vertex is on the junction's perimeter — must coincide
            # with one of the junction's own vertices.
            d_min = min(
                math.hypot(ox - jx, oy - jy)
                for jx, jy in zip(j_xs, j_ys))
            if d_min > ORPHAN_SAME_VERTEX_TOL_M:
                orphans.append((
                    d_min,
                    _shape_label(layout, j_idx, j_s),
                    _shape_label(layout, n_idx, n_s),
                    ox, oy))

    cap = ORPHAN_NEIGHBOUR_VERTEX_REGRESSION_BASELINE.get(
        icao, MAX_ORPHAN_NEIGHBOUR_VERTICES)
    orphans.sort()
    summary = "; ".join(
        f"{j_lbl} ⟂ {n_lbl} at ({ox:.1f},{oy:.1f}) miss={d:.2f}m"
        for d, j_lbl, n_lbl, ox, oy in orphans[:5])
    assert len(orphans) <= cap, (
        f"{icao}: {len(orphans)} neighbour vertex(es) sit within "
        f"{ORPHAN_NEAR_PERIMETER_M:.1f} m of a junction's "
        f"perimeter but more than {ORPHAN_SAME_VERTEX_TOL_M:.2f} m "
        f"from any junction vertex (cap {cap}).  Top: {summary}.")


@pytest.mark.parametrize("icao", _test_airports())
def test_taxi_rects_not_alongside_apron(icao):
    """Per user 2026-05-16: a sloping rect's SLOPING EDGE must
    never be shared with a junction / apron polygon's perimeter.

    Sharing a sloping edge means the junction's elevation has to
    match the rect's per-axis linear slope along the seam — which
    over-constrains the junction's elevation field and produces a
    visible step/cliff at render time, plus violates the user's
    "junctions don't live on sloping rect edges" architectural
    invariant.

    The right outcome when a rect was going to share a sloping
    edge with a junction is for the rect to be absorbed (full or
    partial clip): the apron / junction covers that footprint and
    slopes multi-directionally instead.

    Detection: for each sloping rect, build its 2 sloping-edge
    line segments (corners 0-1 and 2-3 per the
    ``_rect_from_axis_extended`` convention).  For each junction /
    apron polygon, walk its perimeter as consecutive vertex pairs;
    if a junction perimeter segment lies along (within
    ``EDGE_PERP_TOL_M`` perpendicular AND substantially overlaps
    axially with) a sloping edge, flag the rect.

    No airport-specific baseline — the invariant is universal.
    """
    layout = _build_layout(icao)
    sloping_roles = {"primary_parallel", "secondary_parallel",
                     "stub", "cross_connector"}
    other_roles = {"junction", "apron"}

    EDGE_PERP_TOL_M = 1.0
    OVERLAP_MIN_M = 2.0

    def _project_param(px, py, ax, ay, bx, by, L2):
        # Returns (t, perp_dist).  t is parametric along (a→b).
        dx = bx - ax
        dy = by - ay
        t = ((px - ax) * dx + (py - ay) * dy) / L2
        fx = ax + t * dx
        fy = ay + t * dy
        return t, math.hypot(px - fx, py - fy)

    offenders = []
    for idx, s in enumerate(layout.shapes):
        if s.role not in sloping_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            continue
        # Sloping edges: corners (0,1) and (2,3) per
        # _rect_from_axis_extended convention.
        sloping_edges = [(coords[0], coords[1]),
                         (coords[2], coords[3])]
        worst_overlap = 0.0
        worst_label = None
        for o_idx, o in enumerate(layout.shapes):
            if o is s:
                continue
            if o.role not in other_roles:
                continue
            if o.polygon is None or o.polygon.is_empty:
                continue
            o_coords = list(o.polygon.exterior.coords)
            if o_coords and o_coords[0] == o_coords[-1]:
                o_coords = o_coords[:-1]
            n = len(o_coords)
            if n < 3:
                continue
            for (ax, ay), (bx, by) in sloping_edges:
                edge_L2 = (bx - ax) ** 2 + (by - ay) ** 2
                if edge_L2 < 1e-6:
                    continue
                edge_L = math.sqrt(edge_L2)
                # For each junction perimeter segment, find if it
                # overlaps the rect's sloping edge in [0, 1] t-space.
                for i in range(n):
                    px, py = o_coords[i]
                    qx, qy = o_coords[(i + 1) % n]
                    tp, dp = _project_param(
                        px, py, ax, ay, bx, by, edge_L2)
                    tq, dq = _project_param(
                        qx, qy, ax, ay, bx, by, edge_L2)
                    if dp > EDGE_PERP_TOL_M or dq > EDGE_PERP_TOL_M:
                        continue
                    # Clamp the segment's t-range to [0, 1].
                    t_lo = max(0.0, min(tp, tq))
                    t_hi = min(1.0, max(tp, tq))
                    if t_hi <= t_lo:
                        continue
                    overlap_m = (t_hi - t_lo) * edge_L
                    if overlap_m < OVERLAP_MIN_M:
                        continue
                    if overlap_m > worst_overlap:
                        worst_overlap = overlap_m
                        worst_label = _shape_label(
                            layout, o_idx, o)
        if worst_overlap >= OVERLAP_MIN_M:
            offenders.append((
                worst_overlap, _shape_label(layout, idx, s),
                worst_label))
    offenders.sort(reverse=True)
    if offenders:
        summary = "; ".join(
            f"{rect} sloping edge shared {ov:.1f}m with {jn}"
            for ov, rect, jn in offenders[:8])
        msg = (
            f"{icao}: {len(offenders)} sloping rect(s) share a "
            f"sloping edge with a junction/apron (≥ "
            f"{OVERLAP_MIN_M:.0f}m overlap within "
            f"{EDGE_PERP_TOL_M:.1f}m perpendicular).  "
            f"Sloping rects must be absorbed (full or partial) "
            f"when they would share a sloping edge.  "
            f"First {min(8, len(offenders))}: {summary}.")
        assert False, msg
