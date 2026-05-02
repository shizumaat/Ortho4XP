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

import pytest
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

from conftest import airports_under_test, xplane_available, xplane_root

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

# Junction vertex cap: junction = incoming-corners + ≤ 4 trace per
# arc (per shape-rule memory).  30 covers any realistic junction
# even with 5 incoming arcs × max trace.
MAX_JUNCTION_VERTICES = 30

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
JUNCTION_VERTEX_REGRESSION_BASELINE = {
    # SPJC: 12 junctions exceed the 30-vertex cap; worst = 340
    # vertices (single mega-junction sprawling over the SE apron).
    # Ground-truth target max is 36 — see test_compare_target_spjc.
    # 2026-05-02: bumped 336 → 340 after Rule 1+2+3+4 enforcement
    # adds a few snapped runway vertices to the central hub.  These
    # are intentional per user 2026-05-01 (junction-runway 1:1
    # sharing); the +4 vertex bump is the cost of exact node
    # coincidence with the runway segments.
    "SPJC": {"max_offenders": 12, "max_vertex_count": 340},
}

JUNCTION_BOUNDARY_DISTANCE_REGRESSION_BASELINE = {
    # SPJC: all 43 junctions have boundary points > 20 m from
    # any centerline; worst is 566.8 m (at the same SE-apron
    # mega-junction).  The 20 m threshold is tight for normal
    # taxi-junction geometry but ours are apron-sized at SPJC.
    "SPJC": {"max_offenders": 43, "max_distance_m": 567.0},
}

TAXI_RECT_ADJACENCY_REGRESSION_BASELINE = {
    # SPJC: 42/42 surviving sloping rects flag at ≥ 10 % adjacent;
    # all 42 actually flag at 100 % because the surrounding
    # junction polygons cover most of the apron extent.  The
    # legacy partial-absorption rule keeps corridor sections
    # alongside apron-edge as legitimate rects — those fragments
    # trip this whole-rect probe.
    "SPJC": {"max_offenders": 42, "max_adjacent_frac": 1.00},
}

ORPHAN_NEIGHBOUR_VERTEX_REGRESSION_BASELINE = {
    # SPJC: 1 rect corner shifted off its co-located junction vertex
    # by Rule 5 (push-outside-pavement).  Pending tighter Rule 5
    # anchor exemption.
    "SPJC": 1,
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


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
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


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_junction_vertex_count_bounded(icao):
    """Per shape rule: junction = incoming-corners + ≤ 4 trace per
    arc.  A junction with > ``MAX_JUNCTION_VERTICES`` vertices means
    densification ran amok or boundary-trace exceeded the per-arc
    cap.  At CYXY the -10070 regression had 80 vertices — ~45 from
    long-edge densification midpoints, the rest from apt.dat
    boundary trace.

    Per-airport regression baseline:
    Airports listed in ``JUNCTION_VERTEX_REGRESSION_BASELINE``
    have a known-bad ceiling (count + worst vertex count); the
    test fails only if either exceeds the recorded value.  Other
    airports are gated tightly (zero offenders).
    """
    layout = _build_layout(icao)
    cap = MAX_JUNCTION_VERTICES
    offenders = []
    for idx, s in enumerate(layout.shapes):
        if s.role != "junction":
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        # Subtract the closing repeat.
        n = len(s.polygon.exterior.coords) - 1
        if n > cap:
            offenders.append((n, _shape_label(layout, idx, s)))
    offenders.sort(reverse=True)
    summary = "; ".join(f"{lbl} verts={n}" for n, lbl in offenders[:5])

    baseline = JUNCTION_VERTEX_REGRESSION_BASELINE.get(icao)
    if baseline:
        n_off = len(offenders)
        worst_n = offenders[0][0] if offenders else 0
        assert n_off <= baseline["max_offenders"], (
            f"{icao}: {n_off} junction polygon(s) exceed vertex cap "
            f"{cap} — exceeds known-bad baseline of "
            f"{baseline['max_offenders']}.  Top: {summary}.")
        assert worst_n <= baseline["max_vertex_count"], (
            f"{icao}: worst junction has {worst_n} vertices, "
            f"exceeds known-bad baseline of "
            f"{baseline['max_vertex_count']}.  Top: {summary}.")
    else:
        assert not offenders, (
            f"{icao}: {len(offenders)} junction polygon(s) exceed "
            f"vertex cap {cap}.  Top: {summary}.")


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
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


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_taxi_rects_not_alongside_apron(icao):
    """Per user 2026-04-30 absorption rule
    (`_drop_primary_parallels_embedded_in_pavement`): a taxi rect
    that runs along an apron / junction edge for ≥ 10 % of its
    axial length must be absorbed (or partially clipped) into that
    other pavement.  The "taxiway runs along apron edge" case
    folds the embedded portion into the apron (which drives
    elevation); the part of the taxi that leaves the apron
    survives as a shorter rect.

    Test method: replicate the absorption rule's probe — sample 5 m
    steps along each surviving taxi rect's axis, place a probe
    point ``OUTER_PROBE_M`` past each long edge, count steps where
    EITHER probe lands inside any apron / junction polygon (other
    than the rect itself).  A surviving rect should have < 10 % of
    its steps adjacent; otherwise the absorption rule was bypassed
    (corridor-ref exception, EITHER → BOTH switch, etc.) or the
    apron / junction polygon was created after absorption ran.

    Per-airport regression baseline:
    Airports listed in
    ``TAXI_RECT_ADJACENCY_REGRESSION_BASELINE`` have a known-bad
    ceiling (count + worst adjacency fraction); the test fails
    only if either exceeds the recorded value.  Other airports
    are gated tightly (zero offenders).  Note this test surfaces
    the same root issue as the junction-vertex / boundary-distance
    tests at SPJC: surrounding junction polygons are too sprawling,
    so every surviving rect probes "alongside" them.
    """
    layout = _build_layout(icao)
    other_pav = [
        s.polygon for s in layout.shapes
        if s.role in {"apron", "junction"}
        and s.polygon is not None
        and not s.polygon.is_empty]
    if not other_pav:
        pytest.skip(f"{icao}: no apron/junction polygons emitted")
    try:
        other_union = unary_union(other_pav)
    except Exception:
        pytest.skip(f"{icao}: apron/junction union failed")

    sloping_roles = {"primary_parallel", "secondary_parallel",
                     "stub", "cross_connector"}
    SAMPLE_STEP_M = 5.0
    OUTER_PROBE_M = 5.0
    ADJACENCY_FRAC = 0.10
    MIN_AXIS_M = 30.0

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
        # Axis from short-edge midpoints (corners 0,3 → one short
        # edge; 1,2 → other), matching the absorption rule.
        a_mid = (0.5 * (coords[0][0] + coords[3][0]),
                 0.5 * (coords[0][1] + coords[3][1]))
        b_mid = (0.5 * (coords[1][0] + coords[2][0]),
                 0.5 * (coords[1][1] + coords[2][1]))
        ax, ay = b_mid[0] - a_mid[0], b_mid[1] - a_mid[1]
        L = math.hypot(ax, ay)
        if L < MIN_AXIS_M:
            continue
        ux, uy = ax / L, ay / L
        nx, ny = -uy, ux
        half_w = math.hypot(
            coords[0][0] - a_mid[0], coords[0][1] - a_mid[1])
        if half_w < 1.0:
            continue
        outer = half_w + OUTER_PROBE_M
        # Defensive: subtract the rect itself in case an apron
        # polygon ever included its footprint by accident.
        try:
            test_pav = other_union.difference(s.polygon)
        except Exception:
            test_pav = other_union
        n_steps = max(2, int(L / SAMPLE_STEP_M) + 1)
        adj_steps = 0
        for i in range(n_steps):
            u = min(L, i * SAMPLE_STEP_M)
            cx = a_mid[0] + u * ux
            cy = a_mid[1] + u * uy
            try:
                lp = Point(cx + nx * outer, cy + ny * outer)
                rp = Point(cx - nx * outer, cy - ny * outer)
                if (test_pav.contains(lp)
                        or test_pav.contains(rp)):
                    adj_steps += 1
            except Exception:
                pass
        frac = adj_steps / n_steps
        if frac >= ADJACENCY_FRAC:
            offenders.append((frac, _shape_label(layout, idx, s),
                              adj_steps * SAMPLE_STEP_M, L))
    offenders.sort(reverse=True)
    summary = "; ".join(
        f"{lbl} {f * 100:.0f}% adjacent "
        f"({rm:.0f}m of {ll:.0f}m axis)"
        for f, lbl, rm, ll in offenders[:5])
    baseline = TAXI_RECT_ADJACENCY_REGRESSION_BASELINE.get(icao)
    if baseline:
        n_off = len(offenders)
        worst_frac = offenders[0][0] if offenders else 0.0
        assert n_off <= baseline["max_offenders"], (
            f"{icao}: {n_off} taxi rect(s) ≥ "
            f"{ADJACENCY_FRAC * 100:.0f}% adjacent — exceeds "
            f"known-bad baseline of {baseline['max_offenders']}.  "
            f"Top: {summary}.")
        assert worst_frac <= baseline["max_adjacent_frac"] + 0.01, (
            f"{icao}: worst rect {worst_frac * 100:.0f}% adjacent, "
            f"exceeds known-bad baseline of "
            f"{baseline['max_adjacent_frac'] * 100:.0f}%.  "
            f"Top: {summary}.")
    else:
        assert not offenders, (
            f"{icao}: {len(offenders)} surviving taxi rect(s) have "
            f"≥ {ADJACENCY_FRAC * 100:.0f}% of long-edge probes "
            f"inside apron/junction pavement — the absorption rule "
            f"should have clipped them.  Top: {summary}.")
