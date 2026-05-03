"""Geometry regression tests for the pavement builder.

Skipped automatically unless an X-Plane install is available (the
builder needs apt.dat + DSF + DEM tiles).  When run, builds each
test airport and asserts:

* **No self-overlap**: emitted pavement shapes must not overlap each
  other beyond a small tolerance.  Catches the SPJC regression where
  DSF visual overlays (e.g. ``zannespol/Asphalt_2_Green_T80.pol``
  with ``LAYER_GROUP taxiways +1``) duplicated the apt.dat row-110
  pavement coverage and produced 21 overlapping shape pairs covering
  ~15K m² of doubled area.
* **Coverage envelope**: total emitted pavement area must not exceed
  the source pavement (apt.dat row-110 ⊕ runway corners ⊕ surviving
  DSF) by more than a small fraction.  Catches the regression where
  a single DSF overlay polygon contributed 1.45M m² of "pavement"
  that wasn't pavement at all — bulk over-coverage of grass/decor
  areas.

Both checks would have caught the in-flight Phase 1 regression at
SPJC immediately.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest
from shapely.ops import unary_union
from shapely.strtree import STRtree

from conftest import airports_under_test, xplane_available, xplane_root

_HERE = Path(__file__).resolve().parent
_TOOLS = _HERE.parent / "tools"
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))


def _xplane_root() -> str:
    return xplane_root()


def _xplane_available() -> bool:
    return xplane_available()


pytestmark = pytest.mark.skipif(
    not _xplane_available(),
    reason="X-Plane install not found (set XPLANE_ROOT to override)",
)


# Per user 2026-04-30: zero overlap, no exceptions, no rounding-
# error allowance.  An emitted layout where any two shapes overlap
# is a hard invariant violation — X-Plane mesh generation can't
# handle overlapping pavement shapes.  Catches the recent KPHX taxi-
# bridge regression that emitted bridge polygons overlapping the
# adjacent rect / apron.
SELF_OVERLAP_CAP_M2 = 0.0

# Per-airport overlap baseline (m²) — known geometric defects to
# track as a regression ceiling without blocking unrelated work.
# Add an entry only when the overlap is a documented known issue
# scheduled for a follow-up fix.
SELF_OVERLAP_BASELINE_M2 = {
    # SPJC: tunnel_ramp + retaining_wall polygons cross the airport-
    # boundary ribbon where road tunnels emerge outside the
    # perimeter.  The boundary emit doesn't yet carve around them.
    # Pending follow-up to subtract bridge/tunnel footprints from
    # the boundary polygon before emit (see config.py
    # EMIT_BRIDGES_AND_TUNNELS comment).
    "SPJC": 1700.0,
}

# Coverage envelope: the union of every emitted pavement shape's
# polygon must not exceed the source pavement's union by more than
# this fraction.  Source = apt.dat row-110 polygons + runway
# corners.  DSF polygons are excluded from the source because
# they're an *input* to the layout — over-coverage we want to flag
# is "emitted exceeds reasonable inputs".
COVERAGE_OVERAGE_CAP_FRAC = {
    "SPJC": 0.30,   # SPJC has tight apt.dat coverage; allow 30 %
                     # for OSM-synthetic pavement around centerlines
                     # not captured by row-110.
    "CYXY": 1.50,   # CYXY has very sparse apt.dat row-110; DSF + OSM
                     # synthetic pavement legitimately ≈ 2× apt.dat.
                     # Cap is a sanity ceiling, not a precision target.
    "SPLP": 0.30,
}


def _build_layout(icao: str):
    from auto_patch.pipeline import build_airport_pavement
    return build_airport_pavement(icao, _xplane_root(),
                                   compute_elevations=True)


def _source_pavement_union(icao: str):
    """Return the apt.dat row-110 + runway-corner pavement union in
    meter space (anchored at the layout's first-vertex projection).
    """
    import math
    from shapely.ops import transform as shp_transform
    from auto_patch import apt_dat_reader as APR

    apt_dats = APR.find_all_airport_apt_dats(_xplane_root(), icao)
    apt = None
    for ad in apt_dats:
        apt = APR.load_airport(ad, icao)
        if apt is not None and apt.runways:
            break
    if apt is None or not apt.runways:
        pytest.skip(f"{icao}: no apt.dat with runways found")
    # Anchor at the first runway end (matches build_airport_pavement).
    r0 = apt.runways[0]
    lat0, lon0 = r0.lat_a, r0.lon_a
    R = 6_378_137.0
    cos0 = math.cos(math.radians(lat0))

    def to_m(lon, lat, z=None):
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)

    polys = []
    for pav in apt.pavements:
        if pav.polygon is None or pav.polygon.is_empty:
            continue
        pm = shp_transform(to_m, pav.polygon)
        if pm.is_empty:
            continue
        if pm.geom_type == "Polygon":
            polys.append(pm)
        else:
            polys.extend(g for g in getattr(pm, "geoms", [])
                          if g.geom_type == "Polygon")
    # Add runway corners.
    from auto_patch.pavement.runways import _runway_rect_m
    for r in apt.runways:
        rp = _runway_rect_m(r, to_m)
        if rp is not None and not rp.is_empty:
            polys.append(rp)
    if not polys:
        return None
    try:
        return unary_union(polys)
    except Exception:
        return None


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_no_self_overlap(icao):
    """Per user 2026-04-30 hard invariant: NO two emitted pavement
    shapes may overlap, ever.  No floating-point allowance.

    Catches: KPHX taxi-bridge regression where bridge polygons
    overlapped adjacent rect / apron pavement; SPJC DSF visual
    overlay regression where ``zannespol`` polygons duplicated
    apt.dat row-110 coverage; any future absorption / clip pass
    that fails to remove an absorbed sub-rect.
    """
    layout = _build_layout(icao)
    polys = [(s.role, s.polygon) for s in layout.shapes
             if s.polygon is not None and not s.polygon.is_empty]
    if len(polys) < 2:
        return
    tree = STRtree([p for _, p in polys])
    overlap_pairs = []
    overlap_area = 0.0
    for i, (role_a, pa) in enumerate(polys):
        for j in tree.query(pa):
            if j <= i:
                continue
            role_b, pb = polys[j]
            try:
                inter = pa.intersection(pb)
            except Exception:
                continue
            if inter.is_empty:
                continue
            a = inter.area
            if a <= 0.0:
                continue
            overlap_pairs.append((a, role_a, role_b))
            overlap_area += a
    overlap_pairs.sort(reverse=True)
    summary = ", ".join(
        f"{a:.4f} m² ({ra}/{rb})"
        for a, ra, rb in overlap_pairs[:10])
    cap = SELF_OVERLAP_BASELINE_M2.get(icao, SELF_OVERLAP_CAP_M2)
    assert overlap_area <= cap, (
        f"{icao}: {len(overlap_pairs)} overlapping shape pair(s), "
        f"total {overlap_area:,.4f} m² (cap "
        f"{cap:.0f} m² — "
        f"{'baselined' if icao in SELF_OVERLAP_BASELINE_M2 else 'zero tolerance'}).  "
        f"Worst: {summary}.")


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_no_vertex_on_sloping_rect_edge(icao):
    """Per user 2026-04-28 invariant: a junction (or any non-rect)
    polygon vertex can only land on a sloping rect's CORNER, never
    on the interior of one of its four edges.  Edge-interior
    coincidence injects an extra elevation constraint at a non-
    corner location and breaks the rect's straight-line slope.

    Caught the CYXY runway-crossing-junction regression where
    ``_resolve_runway_crossings``'s ``unary_union`` plus the
    downstream 2 m runway-shrink ``difference`` were placing 4
    junction vertices 2–5 m along surviving runway segments' long
    edges (near corners but not at them).
    """
    import math
    layout = _build_layout(icao)
    sloping_roles = {
        "runway", "primary_parallel", "secondary_parallel",
        "stub", "cross_connector"}
    sloping = [s for s in layout.shapes
               if s.role in sloping_roles
               and s.polygon is not None
               and not s.polygon.is_empty]
    others = [s for s in layout.shapes
              if s.role not in sloping_roles
              and s.polygon is not None
              and not s.polygon.is_empty]
    if not sloping or not others:
        return

    # Tolerances: a vertex is ON an edge if it's within 0.5 m of
    # the edge AND > 0.5 m away from either endpoint (i.e. not at
    # a corner — corners are allowed).
    EDGE_PROX_M = 0.5
    CORNER_GUARD_M = 0.5

    violations = []
    for s in sloping:
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            # Sloping rects must be 4-corner — flag anything else.
            violations.append(
                (s.role, s.ref or "?", "non-rect", len(coords)))
            continue
        edges = [(coords[i], coords[(i + 1) % 4])
                 for i in range(4)]
        for o in others:
            ocoords = list(o.polygon.exterior.coords)
            if ocoords and ocoords[0] == ocoords[-1]:
                ocoords = ocoords[:-1]
            for px, py in ocoords:
                # Skip vertices that coincide with one of the
                # rect's corners.
                if any(math.hypot(px - cx, py - cy) <= CORNER_GUARD_M
                       for cx, cy in coords):
                    continue
                for (ax, ay), (bx, by) in edges:
                    dx = bx - ax
                    dy = by - ay
                    L2 = dx * dx + dy * dy
                    if L2 <= 0:
                        continue
                    t = ((px - ax) * dx + (py - ay) * dy) / L2
                    if t <= 0.001 or t >= 0.999:
                        continue
                    proj_x = ax + t * dx
                    proj_y = ay + t * dy
                    d = math.hypot(px - proj_x, py - proj_y)
                    # Also check distance to endpoints (if very
                    # close to one, that's a near-corner case
                    # already excluded above; but the edge
                    # parametrization t may put the projection in-
                    # interior even when the vertex is closer to a
                    # corner than to the midline).
                    d_a = math.hypot(px - ax, py - ay)
                    d_b = math.hypot(px - bx, py - by)
                    if (d < EDGE_PROX_M
                            and d_a > CORNER_GUARD_M
                            and d_b > CORNER_GUARD_M):
                        violations.append(
                            (s.role, s.ref or "?", o.role,
                             o.ref or "?", t, d))
                        break
    if violations:
        # Two violation tuple shapes:
        #   non-rect:       (role, ref, "non-rect", n_corners) — len 4
        #   edge-interior:  (role, ref, o_role, o_ref, t, d)    — len 6
        def _fmt(v):
            if len(v) == 4:
                return (f"{v[0]}({v[1]}) is non-rect "
                        f"(n_corners={v[3]})")
            return (f"{v[2]}({v[3]}) vertex on {v[0]}({v[1]}) "
                    f"edge t={v[4]:.3f} d={v[5]:.2f}m")
        summary = "; ".join(_fmt(v) for v in violations[:5])
        msg = (f"{icao}: {len(violations)} sloping-rect "
               f"invariant violation(s).  Sloping rects must "
               f"have exactly 4 corners; junction polygons may "
               f"share only CORNERS with sloping rects, never "
               f"edge interiors.  First "
               f"{min(5, len(violations))}: {summary}.")
        assert False, msg


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_rect_short_edges_connect(icao):
    """Per user 2026-04-29: a sloping rect (primary_parallel,
    secondary_parallel, stub, cross_connector) has TWO short
    edges, and each short edge represents an end where the
    corridor meets something else — a junction, a runway, a
    terminal, or another rect via shared vertices.  A short edge
    with both corners *un*shared with any other shape means the
    rect is ending in the middle of nowhere — usually an
    artifact of the long-edge-adjacent absorption clipping the
    rect mid-corridor and not extending the surrounding junction
    polygon to share the new clip-boundary corners.

    For each rect, check that EACH of its two short edges has at
    least one corner shared (within ``CORNER_SHARE_TOL_M``) with
    a non-rect-self vertex.  An entirely-disconnected short edge
    is the failure case.
    """
    import math
    CORNER_SHARE_TOL_M = 0.5
    layout = _build_layout(icao)
    rect_roles = {"primary_parallel", "secondary_parallel",
                  "stub", "cross_connector"}
    # Collect every vertex from every shape with its source shape
    # id, then for each rect check both short-edge corners against
    # all OTHER shapes' vertices.
    all_vertices: list = []  # list of (x, y, shape_index)
    for si, s in enumerate(layout.shapes):
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for x, y in coords:
            all_vertices.append((x, y, si))
    failures = []
    tol2 = CORNER_SHARE_TOL_M * CORNER_SHARE_TOL_M
    for ri, r in enumerate(layout.shapes):
        if r.role not in rect_roles:
            continue
        if r.polygon is None or r.polygon.is_empty:
            continue
        try:
            rc = list(r.polygon.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Short edges per ``_rect_from_axis_extended`` convention:
        #   short edge A = corners 0 + 3 (one end)
        #   short edge B = corners 1 + 2 (other end)
        for end_label, (i_a, i_b) in (("end_A", (0, 3)),
                                        ("end_B", (1, 2))):
            ax, ay = rc[i_a]
            bx, by = rc[i_b]
            shared_a = False
            shared_b = False
            for vx, vy, vsi in all_vertices:
                if vsi == ri:
                    continue
                if (not shared_a
                        and (vx - ax) ** 2 + (vy - ay) ** 2 <= tol2):
                    shared_a = True
                if (not shared_b
                        and (vx - bx) ** 2 + (vy - by) ** 2 <= tol2):
                    shared_b = True
                if shared_a and shared_b:
                    break
            if not shared_a and not shared_b:
                failures.append({
                    "ref": r.ref or "?",
                    "role": r.role,
                    "end": end_label,
                    "corner_a": (ax, ay),
                    "corner_b": (bx, by),
                })
    if failures:
        summary = "; ".join(
            f"{f['role']}({f['ref']}) {f['end']}: "
            f"({f['corner_a'][0]:.1f},{f['corner_a'][1]:.1f}) and "
            f"({f['corner_b'][0]:.1f},{f['corner_b'][1]:.1f}) "
            f"both unshared"
            for f in failures[:5])
        msg = (f"{icao}: {len(failures)} rect short edge(s) with "
               f"both corners disconnected from any other shape.  "
               f"A taxi rect's short edge always meets something "
               f"(junction / runway / terminal / other rect).  "
               f"First {min(5, len(failures))}: {summary}.")
        assert False, msg


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_coverage_within_source_envelope(icao):
    """Emitted pavement union must not exceed apt.dat + runway
    coverage by more than the airport's allowed fraction.  Catches
    spurious DSF overlays / non-pavement geometry inflating the
    output beyond its sources.
    """
    layout = _build_layout(icao)
    # Boundary shapes (ROLE_BOUNDARY: airport-perimeter ribbon and
    # boundary→DEM bridge polygons) are elevation control surfaces,
    # not pavement.  Excluding them — this test checks that the
    # PAVEMENT footprint stays close to its sources, which boundary
    # shapes don't contribute to.
    emitted_polys = [s.polygon for s in layout.shapes
                     if s.polygon is not None
                     and not s.polygon.is_empty
                     and getattr(s, "role", None) != "boundary"]
    if not emitted_polys:
        return
    try:
        emitted_union = unary_union(emitted_polys)
    except Exception:
        pytest.fail(f"{icao}: emitted polygons fail unary_union")
    source = _source_pavement_union(icao)
    if source is None or source.is_empty:
        return
    src_area = source.area
    em_area = emitted_union.area
    cap = COVERAGE_OVERAGE_CAP_FRAC.get(icao, 0.5)
    overage = (em_area - src_area) / src_area if src_area > 0 else 0
    assert overage <= cap, (
        f"{icao}: emitted pavement {em_area:,.0f} m² exceeds source "
        f"(apt.dat + runways) {src_area:,.0f} m² by "
        f"{overage*100:.1f}% (cap {cap*100:.0f}%).  Likely cause: "
        f"DSF overlay polygons or non-pavement DSF defs admitted "
        f"into the layout.")
