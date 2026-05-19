"""Junction-polygon elevation repair after the Laplacian solve.

Three concerns, four functions:

* ``_build_clamp_geom_state``: pre-compute spatial-index state of
  every other shape's boundary edges so individual junction
  vertices can be clamped fast.
* ``_clamp_junction_free_vertices``: walk each junction's free
  ring vertices and clamp their elevation to the union of
  feasibility bands implied by every neighbour boundary point
  within ``NEIGHBOUR_CLAMP_RADIUS_M``.
* ``_subdivide_violating_junctions``: split junctions that still
  contain a > ``SUBDIVIDE_VIOLATION_GRADE`` edge after clamp; the
  cut line passes through the worst-grade pair so each new
  sub-polygon admits a feasible elevation field.
* ``_merge_sliver_junctions_into_neighbours``: opposite direction —
  fold tiny sliver junctions into their largest neighbour so JOSM
  doesn't show duplicate-looking polygons.

Public API:
    _build_clamp_geom_state(layout)
    _clamp_junction_free_vertices(layout, geom_state)
    _subdivide_violating_junctions(layout)
    _merge_sliver_junctions_into_neighbours(layout, *, icao)
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

import O4_UI_Utils as UI

# Narrow exception tuple for shapely / numeric-geometry failure
# modes.  Programming errors propagate so they surface immediately.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)

from .canonical_points import snap_polygon_through_registry
from .elevation import (
    NEIGHBOUR_CLAMP_RADIUS_M,
    TAXI_MAX_GRADE,
    _corner_elevation_bucket,
)
from .layout import (
    BuiltShape,
    PavementLayout,
    ROLE_APRON,
    ROLE_CROSS_CONNECTOR,
    ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
)


__all__ = [
    "SUBDIVIDE_MAX_PAIR_DIST_M",
    "SUBDIVIDE_MIN_AREA_M2",
    "SUBDIVIDE_SNAP_RADIUS_M",
    "SUBDIVIDE_VIOLATION_GRADE",
    "_build_clamp_geom_state",
    "_clamp_junction_free_vertices",
    "_merge_sliver_junctions_into_neighbours",
    "_subdivide_violating_junctions",
]


def _build_clamp_geom_state(
        layout: "PavementLayout"
        ) -> "Optional[Tuple[List, List, Dict, set]]":
    """Build the GEOMETRY-ONLY state used by
    :func:`_clamp_junction_free_vertices`.

    The clamp's spatial grid + shared-bucket set are functions of
    polygon geometry alone; only the per-edge elevations change
    between iterations.  Hoisting this build out of the per-call
    body lets the outer fixed-point loop reuse it without
    rebuilding (8+ × cost saved at HECA).

    Returns ``(edge_geom, edge_endpoints, grid, shared_buckets)``
    where:

    * ``edge_geom``:    list of ``(shape_idx, vi_a, vi_b, ax, ay,
                        bx, by)`` — vertex indices index into
                        ``layout.shapes[shape_idx]``'s exterior
                        coords (closed ring's prefix, i.e. the last
                        repeat dropped) so the caller can read
                        current elevations per call.
    * ``edge_endpoints``: bucket key pair per edge (for
                          incident-edge skip).
    * ``grid``:         spatial bucket → list of edge indices.
    * ``shared_buckets``: buckets touched by ≥ 2 shapes.

    Returns None when ``layout.anchor`` is unset.
    """
    if layout.anchor is None:
        return None

    edge_geom: List[Tuple[int, int, int,
                          float, float, float, float]] = []
    edge_endpoints: List[Tuple[Tuple[int, int],
                                Tuple[int, int]]] = []
    bucket_count: Dict[Tuple[int, int], int] = {}
    for si, s in enumerate(layout.shapes):
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        for (cx, cy) in coords:
            bucket = _corner_elevation_bucket(cx, cy)
            bucket_count[bucket] = bucket_count.get(bucket, 0) + 1
        n = len(coords)
        for i in range(n):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % n]
            edge_geom.append((si, i, (i + 1) % n,
                               ax, ay, bx, by))
            edge_endpoints.append((
                _corner_elevation_bucket(ax, ay),
                _corner_elevation_bucket(bx, by)))

    shared_buckets = {b for b, c in bucket_count.items() if c >= 2}

    grid: Dict[Tuple[int, int], List[int]] = {}
    cell = NEIGHBOUR_CLAMP_RADIUS_M
    for ei, (_, _, _, ax, ay, bx, by) in enumerate(edge_geom):
        x0, x1 = (ax, bx) if ax <= bx else (bx, ax)
        y0, y1 = (ay, by) if ay <= by else (by, ay)
        ix0 = int(math.floor((x0 - cell) / cell))
        ix1 = int(math.floor((x1 + cell) / cell))
        iy0 = int(math.floor((y0 - cell) / cell))
        iy1 = int(math.floor((y1 + cell) / cell))
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                grid.setdefault((ix, iy), []).append(ei)

    return (edge_geom, edge_endpoints, grid, shared_buckets)


def _clamp_junction_free_vertices(
        layout: "PavementLayout",
        geom_state: "Optional[Tuple[List, List, Dict, set]]" = None,
        ) -> int:
    """Per-junction free-vertex clamp using every nearby shape
    boundary as a soft anchor (Layer 2).

    Returns the number of free-vertex elevations that changed
    (informational).  Mutates each junction's
    ``node_altitudes`` and ``altitude`` in place.

    When ``geom_state`` is supplied (from
    :func:`_build_clamp_geom_state`), the geometry-only spatial
    grid + shared-bucket set is reused across iterations of the
    outer fixed-point loop.  Falls back to building it on first
    call when the caller doesn't.
    """
    if layout.anchor is None:
        return 0
    if geom_state is None:
        geom_state = _build_clamp_geom_state(layout)
        if geom_state is None:
            return 0
    edge_geom, edge_endpoints, grid, shared_buckets = geom_state
    cell = NEIGHBOUR_CLAMP_RADIUS_M

    # Read current per-shape elevation arrays once per call so the
    # inner clamp loop can look up edge endpoint elevations by
    # (shape_idx, vertex_idx) without re-parsing shapes per edge.
    shape_elevs: Dict[int, List[float]] = {}
    for si, s in enumerate(layout.shapes):
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords_n = len(s.polygon.exterior.coords)
            if coords_n > 0 and (
                s.polygon.exterior.coords[0]
                == s.polygon.exterior.coords[-1]
            ):
                coords_n -= 1
        except _GEOM_EXC:
            continue
        if coords_n <= 0:
            continue
        if s.altitude is not None:
            shape_elevs[si] = [float(s.altitude)] * coords_n
        elif (s.altitude_high is not None
              and s.altitude_low is not None
              and coords_n == 4):
            shape_elevs[si] = [float(s.altitude_high),
                                float(s.altitude_low),
                                float(s.altitude_low),
                                float(s.altitude_high)]
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == coords_n + 1:
                na = na[:-1]
            if len(na) == coords_n:
                shape_elevs[si] = [float(e) for e in na]

    # Reconstruct the boundary_edges layout used downstream
    # ``(shape_idx, ax, ay, bx, by, ea, eb)`` with current elevs.
    boundary_edges: List[Tuple[int, float, float, float, float,
                                float, float]] = []
    boundary_edges_append = boundary_edges.append
    for (si, vi_a, vi_b, ax, ay, bx, by) in edge_geom:
        elevs = shape_elevs.get(si)
        if elevs is None:
            # Skip shapes that didn't yield a valid elev array.
            boundary_edges_append((si, ax, ay, bx, by, 0.0, 0.0))
            continue
        boundary_edges_append((si, ax, ay, bx, by,
                                elevs[vi_a], elevs[vi_b]))

    rect_like_roles = {ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                       ROLE_CROSS_CONNECTOR, ROLE_TERMINAL}

    n_changed = 0
    for si, s in enumerate(layout.shapes):
        if s.role != ROLE_JUNCTION:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        # Build per-vertex elevations and FREE-vs-anchored flags.
        if s.altitude is not None:
            elevs = [float(s.altitude)] * len(coords)
            uniform_alt = True
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == len(coords) + 1:
                na = na[:-1]
            if len(na) != len(coords):
                continue
            elevs = [float(e) for e in na]
            uniform_alt = False
        else:
            continue
        # A vertex is ANCHORED if its bucket is shared with another
        # shape (rect/runway/terminal or another junction with the
        # same coord).  We check via the shared_buckets set.
        is_anchor = []
        for (cx, cy) in coords:
            b = _corner_elevation_bucket(cx, cy)
            is_anchor.append(b in shared_buckets)
        # Now clamp each FREE vertex.
        new_elevs = list(elevs)
        for vi, (cx, cy) in enumerate(coords):
            if is_anchor[vi]:
                continue
            ix = int(math.floor(cx / cell))
            iy = int(math.floor(cy / cell))
            lo_v = float("-inf")
            hi_v = float("inf")
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    bucket = (ix + dx, iy + dy)
                    if bucket not in grid:
                        continue
                    v_bucket = _corner_elevation_bucket(cx, cy)
                    for ei in grid[bucket]:
                        (s_other, ax, ay, bx, by,
                         ea, eb) = boundary_edges[ei]
                        # Per user 2026-05-18: cross-shape grade
                        # enforcement happens ONLY at shared corners
                        # (via the OSM emitter's altitude-aware
                        # vertex bucketing).  Two non-shared shapes
                        # — a taxiway 70 m parallel to a runway, for
                        # example — must be free to slope
                        # independently along their own axes.
                        # Restrict the clamp to same-junction edges
                        # so cross-shape proximity doesn't impose a
                        # phantom grade constraint.
                        if s_other != si:
                            continue
                        # Skip ONLY edges incident to THIS vertex —
                        # they trivially equal the vertex's own
                        # elevation and would lock it in place.
                        ek0, ek1 = edge_endpoints[ei]
                        if v_bucket == ek0 or v_bucket == ek1:
                            continue
                        edx = bx - ax
                        edy = by - ay
                        seg2 = edx * edx + edy * edy
                        if seg2 < 0.04:
                            continue
                        t = ((cx - ax) * edx + (cy - ay) * edy) / seg2
                        if t < 0.0:
                            t = 0.0
                        elif t > 1.0:
                            t = 1.0
                        ccx = ax + t * edx
                        ccy = ay + t * edy
                        d = math.hypot(cx - ccx, cy - ccy)
                        if d > NEIGHBOUR_CLAMP_RADIUS_M:
                            continue
                        e_at = ea + t * (eb - ea)
                        band = max(d, 0.01) * TAXI_MAX_GRADE
                        if e_at - band > lo_v:
                            lo_v = e_at - band
                        if e_at + band < hi_v:
                            hi_v = e_at + band
            if lo_v > hi_v:
                # Conflicting nearby anchors — pick midpoint as a
                # least-squares-style compromise.  Layer 1 will
                # eventually reconcile the source corners.
                target = 0.5 * (lo_v + hi_v)
            else:
                target = elevs[vi]
                if target < lo_v:
                    target = lo_v
                if target > hi_v:
                    target = hi_v
            target = round(target, 1)
            if abs(target - elevs[vi]) > 0.05:
                new_elevs[vi] = target
                n_changed += 1
        # Persist changes.
        if uniform_alt:
            # Was a single altitude; if any free vertex shifted, we
            # have to switch to per-vertex node_altitudes.
            if any(abs(new_elevs[i] - elevs[i]) > 0.05
                   for i in range(len(elevs))):
                elev_range = max(new_elevs) - min(new_elevs)
                if elev_range < 0.05:
                    s.altitude = round(
                        sum(new_elevs) / len(new_elevs), 1)
                else:
                    s.altitude = None
                    closed = list(new_elevs) + [new_elevs[0]]
                    s.node_altitudes = closed
        else:
            # Already per-vertex.  Update node_altitudes (ring +
            # closing repeat).
            closed = list(new_elevs) + [new_elevs[0]]
            s.node_altitudes = closed
            # If everything collapsed to one value, switch back to
            # flat altitude.
            elev_range = max(new_elevs) - min(new_elevs)
            if elev_range < 0.05:
                s.altitude = round(
                    sum(new_elevs) / len(new_elevs), 1)
                s.node_altitudes = None
    return n_changed


SUBDIVIDE_VIOLATION_GRADE = 0.02   # 2 % — attempt subdivision
                                    # whenever the worst within-
                                    # shape pair exceeds the taxi
                                    # grade cap.  Lowered from 10 %
                                    # on 2026-05-05: the per-surface
                                    # solver leaves residual 1.5–3 %
                                    # violations on very large
                                    # junctions sitting over DEM
                                    # spikes; cutting them lets
                                    # each sub-polygon converge to
                                    # its own DEM-floor.
SUBDIVIDE_MAX_PAIR_DIST_M = 60.0   # only consider pairs within
                                    # this radius — same as
                                    # check_grade's
                                    # WITHIN_SHAPE_MAX_PAIR_DIST_M
                                    # (Triangle4XP-plausible edge).
SUBDIVIDE_MIN_AREA_M2 = 5.0        # don't emit sub-polygons
                                    # smaller than this — they'd
                                    # become slivers and re-trigger
                                    # the sliver-corner safety net.


SUBDIVIDE_SNAP_RADIUS_M = 5.0      # snap new cut-line vertices
                                    # to existing ring vertices when
                                    # within this radius.  Without
                                    # snapping, the cut introduces
                                    # 0.5-2 m new vertices very
                                    # close to existing ones; the
                                    # interpolated-vs-original
                                    # elevation mismatch over those
                                    # tiny distances produces huge
                                    # spurious grade percentages
                                    # (worse than the original
                                    # violation we were trying to fix).


def _subdivide_violating_junctions(layout: "PavementLayout") -> int:
    """Split junction polygons whose worst within-shape vertex pair
    exceeds ``SUBDIVIDE_VIOLATION_GRADE`` along a perpendicular
    cut through the midpoint of the violating pair.

    Cut-vertex placement: new vertices created where the cut line
    intersects the polygon boundary are SNAPPED to existing ring
    vertices within ``SUBDIVIDE_SNAP_RADIUS_M``.  Without snapping,
    the cut creates ring-adjacent vertex pairs spaced 0.5-2 m apart
    whose interpolated-vs-original elevations differ by small
    amounts — registering as huge grade percentages
    (e.g. 0.4 m / 0.5 m = 73.9 %) that are WORSE than the original
    violation we were trying to fix.

    Validation: a sub-polygon is only accepted when its own worst
    within-shape vertex pair is BETTER (smaller grade) than the
    parent's worst pair.  If the cut would produce a sub-polygon
    that's MORE violating, the original is kept and the cut is
    abandoned — prevents the iterative subdivision from making
    things worse.

    Returns the number of polygons that were subdivided
    (informational).
    """
    if not layout.shapes:
        return 0
    from shapely.geometry import LineString
    from shapely.ops import split as _shapely_split
    n_subdivided = 0
    new_shapes: List[BuiltShape] = []
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            new_shapes.append(s)
            continue
        if s.polygon is None or s.polygon.is_empty:
            new_shapes.append(s)
            continue
        try:
            ring = list(s.polygon.exterior.coords)
        except _GEOM_EXC:
            new_shapes.append(s)
            continue
        if ring and ring[0] == ring[-1]:
            ring = ring[:-1]
        n = len(ring)
        if n < 4:
            new_shapes.append(s)
            continue
        if s.altitude is not None:
            elevs = [float(s.altitude)] * n
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == n + 1:
                na = na[:-1]
            if len(na) != n:
                new_shapes.append(s)
                continue
            elevs = [float(e) for e in na]
        else:
            new_shapes.append(s)
            continue

        def _worst_grade(pts: List[Tuple[float, float]],
                         es: List[float]) -> Tuple[float,
                                                    Optional[Tuple[int, int]]]:
            """Worst all-pair Euclidean grade within the polygon.
            Per user 2026-05-18: junction grade applies across the
            entire interior surface, not just along centerlines —
            same rule as aprons."""
            radius2 = SUBDIVIDE_MAX_PAIR_DIST_M ** 2
            min_d2 = 0.5 ** 2
            wg = 0.0
            wp = None
            m = len(pts)
            for a in range(m):
                xa, ya = pts[a]
                ea = es[a]
                for b in range(a + 1, m):
                    xb, yb = pts[b]
                    dx_ = xa - xb
                    dy_ = ya - yb
                    d2_ = dx_ * dx_ + dy_ * dy_
                    if d2_ < min_d2 or d2_ > radius2:
                        continue
                    d_ = math.sqrt(d2_)
                    de_ = abs(ea - es[b])
                    if de_ <= TAXI_MAX_GRADE * d_ + 0.10:
                        continue
                    g_ = de_ / d_
                    if g_ > wg:
                        wg = g_
                        wp = (a, b)
            return wg, wp

        worst_grade, worst_pair = _worst_grade(ring, elevs)
        if (worst_pair is None
                or worst_grade < SUBDIVIDE_VIOLATION_GRADE):
            new_shapes.append(s)
            continue
        # Build the perpendicular cut line through the midpoint.
        i, j = worst_pair
        pi = ring[i]
        pj = ring[j]
        mx = 0.5 * (pi[0] + pj[0])
        my = 0.5 * (pi[1] + pj[1])
        dx = pj[0] - pi[0]
        dy = pj[1] - pi[1]
        seg_len = math.hypot(dx, dy)
        if seg_len < 1e-6:
            new_shapes.append(s)
            continue
        # Perpendicular unit vector (rotate +90°).
        px = -dy / seg_len
        py = dx / seg_len
        bx_min, by_min, bx_max, by_max = s.polygon.bounds
        bbox_diag = math.hypot(bx_max - bx_min, by_max - by_min)
        L = max(bbox_diag * 2.0, 1000.0)
        cut = LineString([
            (mx - L * px, my - L * py),
            (mx + L * px, my + L * py),
        ])
        try:
            parts = _shapely_split(s.polygon, cut)
        except _GEOM_EXC:
            new_shapes.append(s)
            continue
        sub_polys: List[Polygon] = []
        try:
            for g in getattr(parts, "geoms", [parts]):
                if (g is None or g.is_empty
                        or g.geom_type != "Polygon"
                        or g.area < SUBDIVIDE_MIN_AREA_M2):
                    continue
                sub_polys.append(g)
        except _GEOM_EXC:
            new_shapes.append(s)
            continue
        if len(sub_polys) < 2:
            new_shapes.append(s)
            continue

        # Pre-compute snap-target ring vertices keyed by their
        # squared snap radius for O(n) lookup per sub-vertex.
        snap_r2 = SUBDIVIDE_SNAP_RADIUS_M ** 2

        def _snap_to_ring(qx: float, qy: float
                          ) -> Tuple[float, float, int]:
            """Return the closest ring vertex within snap radius
            and its index, or ``(qx, qy, -1)`` if no ring vertex
            is close enough.
            """
            best_k = -1
            best_d2 = snap_r2
            for k in range(n):
                rx, ry = ring[k]
                d2 = (rx - qx) ** 2 + (ry - qy) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    best_k = k
            if best_k >= 0:
                rx, ry = ring[best_k]
                return rx, ry, best_k
            return qx, qy, -1

        def _lookup_elev(qx: float, qy: float, hint_idx: int
                         ) -> float:
            """Look up elevation for sub-polygon vertex.  When
            ``hint_idx`` ≥ 0 (snapped to ring vertex) use the
            original elevation directly.  Otherwise interpolate
            along the closest original ring edge.
            """
            if hint_idx >= 0:
                return elevs[hint_idx]
            best_d2 = float("inf")
            best_e = elevs[0]
            for k in range(n):
                ax, ay = ring[k]
                bx, by = ring[(k + 1) % n]
                edx = bx - ax
                edy = by - ay
                seg2 = edx * edx + edy * edy
                if seg2 < 1e-9:
                    continue
                t = ((qx - ax) * edx + (qy - ay) * edy) / seg2
                if t < 0.0:
                    t = 0.0
                elif t > 1.0:
                    t = 1.0
                ccx = ax + t * edx
                ccy = ay + t * edy
                d2 = (qx - ccx) ** 2 + (qy - ccy) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    ea = elevs[k]
                    eb = elevs[(k + 1) % n]
                    best_e = ea + t * (eb - ea)
            return round(float(best_e), 1)

        # Build snapped sub-polygons + per-vertex elevations.
        # Acceptance criterion: each sub-polygon's own worst pair
        # must be MEASURABLY better (≥ 0.5 % grade improvement)
        # than the parent's worst.  This prevents cuts that
        # technically separate the worst pair but introduce new
        # cut-line-vertex pairs of similar magnitude (the bug
        # that made the relaxed version worse than the strict).
        validated_subs: List[Tuple[Polygon, List[Tuple[float, float]],
                                    List[float]]] = []
        cut_was_useful = True
        for sp in sub_polys:
            sub_ring_raw = list(sp.exterior.coords)
            if sub_ring_raw and sub_ring_raw[0] == sub_ring_raw[-1]:
                sub_ring_raw = sub_ring_raw[:-1]
            snapped_pts: List[Tuple[float, float]] = []
            snapped_hints: List[int] = []
            for (qx, qy) in sub_ring_raw:
                sx, sy, hint = _snap_to_ring(qx, qy)
                if (snapped_pts and abs(snapped_pts[-1][0] - sx) < 1e-9
                        and abs(snapped_pts[-1][1] - sy) < 1e-9):
                    continue  # consecutive dup after snap
                snapped_pts.append((sx, sy))
                snapped_hints.append(hint)
            while (len(snapped_pts) >= 2
                   and abs(snapped_pts[0][0]
                           - snapped_pts[-1][0]) < 1e-9
                   and abs(snapped_pts[0][1]
                           - snapped_pts[-1][1]) < 1e-9):
                snapped_pts.pop()
                snapped_hints.pop()
            if len(snapped_pts) < 3:
                cut_was_useful = False
                break
            try:
                snapped_poly = Polygon(snapped_pts)
                if not snapped_poly.is_valid:
                    snapped_poly = snapped_poly.buffer(0)
                if (snapped_poly.is_empty
                        or snapped_poly.geom_type != "Polygon"
                        or snapped_poly.area < SUBDIVIDE_MIN_AREA_M2):
                    cut_was_useful = False
                    break
            except _GEOM_EXC:
                cut_was_useful = False
                break
            sub_elevs = [_lookup_elev(qx, qy, h)
                         for (qx, qy), h
                         in zip(snapped_pts, snapped_hints)]
            sub_worst, _ = _worst_grade(snapped_pts, sub_elevs)
            if sub_worst >= worst_grade - 0.005:
                cut_was_useful = False
                break
            validated_subs.append(
                (snapped_poly, snapped_pts, sub_elevs))

        if not cut_was_useful or len(validated_subs) < 2:
            # Fallback: iso-elevation cut.  When the perpendicular
            # cut through the worst-pair midpoint can't be validated
            # (typically because it produces a tiny sliver containing
            # the high-z corners + cut endpoints, whose own worst
            # pair isn't measurably better), try cutting along the
            # MEDIAN-elevation contour instead.  Walk ring edges; an
            # edge "crosses" the median when its endpoints straddle
            # it.  For a polygon that wraps two distinct elevation
            # regions (SPLP junction-10053: 4 corners at z=70, 2 at
            # z=77), exactly two edges cross — a clean cut between
            # the crossing midpoints separates the two regions.
            iso_subs = _try_iso_elevation_cut(
                s, ring, elevs, worst_grade)
            if iso_subs is not None and len(iso_subs) >= 2:
                validated_subs = iso_subs
            else:
                new_shapes.append(s)
                continue

        _registry = getattr(layout, "canonical_points", None)
        for sp, sub_pts, sub_elevs in validated_subs:
            # Route through registry so the cut produces canonical
            # corner coordinates shared with the parent ring's
            # vertices (which are already in the registry).
            if _registry is not None:
                snapped_sp = snap_polygon_through_registry(
                    sp, _registry)
                if (snapped_sp is None or snapped_sp.is_empty
                        or snapped_sp.geom_type != "Polygon"):
                    continue
                sp = snapped_sp
            sub_shape = BuiltShape(
                polygon=sp, role=ROLE_JUNCTION, ref=s.ref)
            elev_range = max(sub_elevs) - min(sub_elevs)
            if elev_range < 0.05:
                sub_shape.altitude = round(
                    sum(sub_elevs) / len(sub_elevs), 1)
            else:
                closed = list(sub_elevs) + [sub_elevs[0]]
                sub_shape.node_altitudes = closed
            new_shapes.append(sub_shape)
        n_subdivided += 1
    layout.shapes = new_shapes
    return n_subdivided


def _try_iso_elevation_cut(  # noqa: C901 (long helper, see body)
    s: "BuiltShape",
    ring: List[Tuple[float, float]],
    elevs: List[float],
    parent_worst_grade: float,
) -> Optional[List[Tuple["Polygon", List[Tuple[float, float]],
                         List[float]]]]:
    """Cut a junction polygon along its median-elevation contour.

    Returns a list of validated sub-polygons (each
    ``(polygon, snapped_pts, sub_elevs)``) or ``None`` if the cut
    isn't applicable.

    Rationale: when the polygon wraps two distinct elevation regions
    (e.g. corners 0-3 at z=70, corners 4-5 at z=77), the
    perpendicular-cut subdivide produces a tiny corner-sliver that
    fails validation.  An iso-elevation cut walks ring edges,
    finds the two edges where elevation crosses the median (between
    z_min and z_max), and cuts from one crossing point to the other.
    Each resulting sub-polygon has elevation range half the parent's.
    """
    n = len(ring)
    if n < 4:
        return None
    z_min = min(elevs)
    z_max = max(elevs)
    if z_max - z_min < 0.5:
        return None
    median = 0.5 * (z_min + z_max)
    crossings: List[Tuple[float, float, int]] = []
    for k in range(n):
        z_a = elevs[k]
        z_b = elevs[(k + 1) % n]
        if (z_a < median) == (z_b < median):
            continue  # both on same side
        if abs(z_b - z_a) < 1e-6:
            continue  # too flat to interpolate
        t = (median - z_a) / (z_b - z_a)
        if t <= 0.001 or t >= 0.999:
            continue
        ax, ay = ring[k]
        bx, by = ring[(k + 1) % n]
        cx = ax + t * (bx - ax)
        cy = ay + t * (by - ay)
        crossings.append((cx, cy, k))
    if len(crossings) != 2:
        return None
    from shapely.geometry import LineString
    from shapely.ops import split as _shapely_split
    cx0, cy0, _ = crossings[0]
    cx1, cy1, _ = crossings[1]
    cut = LineString([(cx0, cy0), (cx1, cy1)])
    try:
        parts = _shapely_split(s.polygon, cut)
    except _GEOM_EXC:
        return None
    sub_polys: List[Polygon] = []
    for g in getattr(parts, "geoms", [parts]):
        if (g is None or g.is_empty
                or g.geom_type != "Polygon"
                or g.area < SUBDIVIDE_MIN_AREA_M2):
            continue
        sub_polys.append(g)
    if len(sub_polys) < 2:
        return None
    snap_r2 = SUBDIVIDE_SNAP_RADIUS_M ** 2

    def _snap_to_ring(qx: float, qy: float
                      ) -> Tuple[float, float, int]:
        best_k = -1
        best_d2 = snap_r2
        for k in range(n):
            rx, ry = ring[k]
            d2 = (rx - qx) ** 2 + (ry - qy) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_k = k
        if best_k >= 0:
            rx, ry = ring[best_k]
            return rx, ry, best_k
        return qx, qy, -1

    def _lookup_elev(qx: float, qy: float, hint_idx: int) -> float:
        if hint_idx >= 0:
            return elevs[hint_idx]
        best_d2 = float("inf")
        best_e = elevs[0]
        for k in range(n):
            ax, ay = ring[k]
            bx, by = ring[(k + 1) % n]
            edx = bx - ax
            edy = by - ay
            seg2 = edx * edx + edy * edy
            if seg2 < 1e-9:
                continue
            t = ((qx - ax) * edx + (qy - ay) * edy) / seg2
            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0
            ccx = ax + t * edx
            ccy = ay + t * edy
            d2 = (qx - ccx) ** 2 + (qy - ccy) ** 2
            if d2 < best_d2:
                best_d2 = d2
                ea = elevs[k]
                eb = elevs[(k + 1) % n]
                best_e = ea + t * (eb - ea)
        return round(float(best_e), 1)

    validated_subs: List[Tuple[Polygon, List[Tuple[float, float]],
                                List[float]]] = []
    for sp in sub_polys:
        sub_ring_raw = list(sp.exterior.coords)
        if sub_ring_raw and sub_ring_raw[0] == sub_ring_raw[-1]:
            sub_ring_raw = sub_ring_raw[:-1]
        snapped_pts: List[Tuple[float, float]] = []
        snapped_hints: List[int] = []
        for (qx, qy) in sub_ring_raw:
            sx, sy, hint = _snap_to_ring(qx, qy)
            if (snapped_pts
                    and abs(snapped_pts[-1][0] - sx) < 1e-9
                    and abs(snapped_pts[-1][1] - sy) < 1e-9):
                continue
            snapped_pts.append((sx, sy))
            snapped_hints.append(hint)
        while (len(snapped_pts) >= 2
               and abs(snapped_pts[0][0] - snapped_pts[-1][0]) < 1e-9
               and abs(snapped_pts[0][1] - snapped_pts[-1][1]) < 1e-9):
            snapped_pts.pop()
            snapped_hints.pop()
        if len(snapped_pts) < 3:
            return None
        try:
            snapped_poly = Polygon(snapped_pts)
            if not snapped_poly.is_valid:
                snapped_poly = snapped_poly.buffer(0)
            if (snapped_poly.is_empty
                    or snapped_poly.geom_type != "Polygon"
                    or snapped_poly.area < SUBDIVIDE_MIN_AREA_M2):
                return None
        except _GEOM_EXC:
            return None
        sub_elevs = [_lookup_elev(qx, qy, h)
                     for (qx, qy), h
                     in zip(snapped_pts, snapped_hints)]
        # Compute this sub's worst-pair grade
        sub_worst = 0.0
        sm = len(snapped_pts)
        radius2 = SUBDIVIDE_MAX_PAIR_DIST_M ** 2
        for a in range(sm):
            xa, ya = snapped_pts[a]
            ea = sub_elevs[a]
            for b in range(a + 1, sm):
                xb, yb = snapped_pts[b]
                d2 = (xa - xb) ** 2 + (ya - yb) ** 2
                if d2 < 0.25 or d2 > radius2:
                    continue
                d_ = math.sqrt(d2)
                de_ = abs(ea - sub_elevs[b])
                if de_ <= TAXI_MAX_GRADE * d_ + 0.10:
                    continue
                g_ = de_ / d_
                if g_ > sub_worst:
                    sub_worst = g_
        # Iso cut should produce sub-polygons that are MEASURABLY
        # better; if not, this cut isn't useful either.
        if sub_worst >= parent_worst_grade - 0.005:
            return None
        validated_subs.append((snapped_poly, snapped_pts, sub_elevs))
    if len(validated_subs) < 2:
        return None
    return validated_subs


def _merge_sliver_junctions_into_neighbours(
        layout: "PavementLayout",
        icao: str = "",
        sliver_area_m2: float = 1000.0,
        sliver_ratio: float = 0.05,
        shared_vertex_tol_m: float = 0.5,
        ) -> int:
    """Merge small junction polygons into adjacent larger ones.

    A "sliver" is a junction polygon whose area is below
    ``sliver_area_m2`` AND whose ratio to a neighbour's area is
    below ``sliver_ratio``.  Two junctions are "adjacent" if they
    share at least 2 boundary vertices within
    ``shared_vertex_tol_m``.

    Common cause: ``_decompose_polygon_with_holes`` cuts a
    polygon-with-holes into simple pieces, occasionally carving
    off a tiny strip when the cut grazes the polygon's edge.
    The strip and main piece share a boundary segment (the cut
    line); the strip should be merged back.  Subdivision passes
    in ``_compute_elevations`` can produce similar slivers.

    Returns the number of slivers merged.
    """
    junction_idxs = [i for i, s in enumerate(layout.shapes)
                     if s.role == ROLE_JUNCTION
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if len(junction_idxs) < 2:
        return 0
    # Cache per-shape vertex sets in meter coords for fast tests.
    j_verts: Dict[int, List[Tuple[float, float]]] = {}
    for i in junction_idxs:
        try:
            coords = list(layout.shapes[i].polygon.exterior.coords)
        except _GEOM_EXC:
            j_verts[i] = []
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        j_verts[i] = coords
    tol2 = shared_vertex_tol_m * shared_vertex_tol_m
    merge_into: Dict[int, int] = {}
    for i in junction_idxs:
        ai = layout.shapes[i].polygon.area
        if ai >= sliver_area_m2:
            continue
        best_idx: Optional[int] = None
        best_area = 0.0
        for j in junction_idxs:
            if j == i:
                continue
            aj = layout.shapes[j].polygon.area
            if aj <= ai:
                continue
            if (ai / aj) > sliver_ratio:
                continue
            shared = 0
            for vx, vy in j_verts[i]:
                for ux, uy in j_verts[j]:
                    if (vx - ux) ** 2 + (vy - uy) ** 2 <= tol2:
                        shared += 1
                        break
                if shared >= 2:
                    break
            if shared >= 2 and aj > best_area:
                best_idx = j
                best_area = aj
        if best_idx is not None:
            merge_into[i] = best_idx
    if not merge_into:
        return 0
    for sliver_i, target_j in merge_into.items():
        try:
            target_shape = layout.shapes[target_j]
            sliver_shape = layout.shapes[sliver_i]
            merged = unary_union([
                target_shape.polygon,
                sliver_shape.polygon])
            if (merged.geom_type == "Polygon"
                    and not merged.is_empty):
                # Build a per-vertex elevation lookup from the
                # ORIGINAL target + sliver vertices, then re-derive
                # node_altitudes for the merged polygon by nearest-
                # neighbour match.  Per user 2026-04-29 (CYXY apron
                # regression): just setting ``node_altitudes = None``
                # leaves the merged polygon with NO elevation at all
                # (no downstream re-derivation step exists), which
                # makes X-Plane interpolate from neighbour shapes
                # and produces the "terrain all over the place"
                # apron the user saw.
                lookup: List[Tuple[float, float, float]] = []
                for src_shape in (target_shape, sliver_shape):
                    if not src_shape.node_altitudes:
                        # Sloped or flat alternatives.
                        if (src_shape.altitude_high is not None
                                and src_shape.altitude_low is not None):
                            avg = 0.5 * (
                                src_shape.altitude_high
                                + src_shape.altitude_low)
                        elif src_shape.altitude is not None:
                            avg = float(src_shape.altitude)
                        else:
                            continue
                        try:
                            sc = list(
                                src_shape.polygon.exterior.coords)
                        except _GEOM_EXC:
                            continue
                        if sc and sc[0] == sc[-1]:
                            sc = sc[:-1]
                        for sx, sy in sc:
                            lookup.append((sx, sy, float(avg)))
                        continue
                    src_alts = list(src_shape.node_altitudes)
                    try:
                        sc = list(src_shape.polygon.exterior.coords)
                    except _GEOM_EXC:
                        continue
                    if sc and sc[0] == sc[-1]:
                        sc = sc[:-1]
                    if (len(src_alts) == len(sc) + 1
                            and src_alts[0] == src_alts[-1]):
                        src_alts = src_alts[:-1]
                    for k, (sx, sy) in enumerate(sc):
                        if k >= len(src_alts):
                            break
                        lookup.append(
                            (sx, sy, float(src_alts[k])))
                if not lookup:
                    layout.shapes[target_j].polygon = merged
                    layout.shapes[target_j].node_altitudes = None
                    continue
                merged_coords = list(merged.exterior.coords)
                if (merged_coords
                        and merged_coords[0] == merged_coords[-1]):
                    merged_coords_open = merged_coords[:-1]
                else:
                    merged_coords_open = merged_coords
                new_alts: List[float] = []
                for mx, my in merged_coords_open:
                    best_d2 = float("inf")
                    best_alt = 0.0
                    for sx, sy, sa in lookup:
                        d2 = (mx - sx) ** 2 + (my - sy) ** 2
                        if d2 < best_d2:
                            best_d2 = d2
                            best_alt = sa
                    new_alts.append(round(best_alt, 1))
                # Closing-vertex repeat for the ring.
                if new_alts:
                    new_alts.append(new_alts[0])
                target_shape.polygon = merged
                target_shape.node_altitudes = new_alts
        except _GEOM_EXC:
            continue
    sliver_set = set(merge_into.keys())
    layout.shapes = [
        s for k, s in enumerate(layout.shapes)
        if k not in sliver_set]
    try:
        UI.vprint(1,
            f"  [pav-builder] {icao}: merged "
            f"{len(merge_into)} sliver junction(s) into "
            f"adjacent larger junctions.")
    except _GEOM_EXC:
        pass
    return len(merge_into)


def _drop_thin_orphan_slivers(
        layout: "PavementLayout",
        icao: str = "",
        max_area_m2: float = 1000.0,
        min_aspect: float = 4.0,
        shared_vertex_tol_m: float = 0.5,
        ) -> int:
    """Drop thin sliver junction polygons that form residue along a
    stub / parallel rect's long edge.

    A "thin sliver" qualifies when ALL of:
      * Junction polygon area < ``max_area_m2``.
      * Aspect ratio (perimeter / (2*sqrt(pi*area))) > ``min_aspect``
        — i.e. shaped like a long thin strip, not a chunky polygon.
      * Shares ≥ 2 boundary vertices with a stub / primary_parallel
        / secondary_parallel / cross_connector (the rect whose long
        edge the sliver lies along).

    These typically form because the apt.dat row-110 pavement
    boundary curves inward between a stub's two end corners,
    leaving a thin pavement strip uncovered by the rect.  The
    strip gets emitted as a separate junction polygon that's
    geometrically isolated from the surrounding apron (only
    touches via stub-shared vertices, not via the apron boundary),
    so ``_merge_sliver_junctions_into_neighbours`` can't catch it
    (it requires ≥ 2 shared vertices with another JUNCTION).

    Dropping these slivers loses < 0.2 % of pavement coverage at
    typical airports but eliminates visible thin residue along
    diagonal-stub sloping edges (user 2026-05-12).

    Returns count of slivers dropped.
    """
    RECT_LIKE_ROLES = (
        ROLE_STUB, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_CROSS_CONNECTOR,
    )
    junction_idxs = [i for i, s in enumerate(layout.shapes)
                     if s.role == ROLE_JUNCTION
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not junction_idxs:
        return 0
    rect_idxs = [i for i, s in enumerate(layout.shapes)
                 if s.role in RECT_LIKE_ROLES
                 and s.polygon is not None
                 and not s.polygon.is_empty]
    if not rect_idxs:
        return 0

    # Cache per-shape open-ring vertices in meter coords.
    verts: Dict[int, List[Tuple[float, float]]] = {}
    for i in junction_idxs + rect_idxs:
        try:
            c = list(layout.shapes[i].polygon.exterior.coords)
        except _GEOM_EXC:
            verts[i] = []
            continue
        if c and c[0] == c[-1]:
            c = c[:-1]
        verts[i] = c

    tol2 = shared_vertex_tol_m * shared_vertex_tol_m
    to_drop: List[int] = []
    for i in junction_idxs:
        p = layout.shapes[i].polygon
        try:
            area = p.area
        except _GEOM_EXC:
            continue
        if area <= 0 or area >= max_area_m2:
            continue
        try:
            peri = p.length
        except _GEOM_EXC:
            continue
        aspect = peri / (2.0 * math.sqrt(math.pi * area))
        if aspect <= min_aspect:
            continue
        # Confirm 2+ shared vertices with a rect-like neighbour.
        share_with_rect = False
        for j in rect_idxs:
            shared = 0
            for vx, vy in verts[i]:
                for ux, uy in verts[j]:
                    if (vx - ux) ** 2 + (vy - uy) ** 2 <= tol2:
                        shared += 1
                        break
                if shared >= 2:
                    break
            if shared >= 2:
                share_with_rect = True
                break
        if share_with_rect:
            to_drop.append(i)

    if not to_drop:
        return 0
    drop_set = set(to_drop)
    layout.shapes = [
        s for k, s in enumerate(layout.shapes)
        if k not in drop_set]
    try:
        UI.vprint(1,
            f"  [pav-builder] {icao}: dropped "
            f"{len(to_drop)} thin orphan sliver junction(s) "
            f"(residue along stub/parallel long edges).")
    except _GEOM_EXC:
        pass
    return len(to_drop)


def _split_sloped_rects_at_violations(
        layout: "PavementLayout",
        icao: str = "",
        on_edge_tol_m: float = 1.5,
        min_t_margin: float = 0.05,
        ) -> int:
    """Split sloped 4-corner rects when a junction vertex lies on
    one of their sloping (long) edges.

    A sloped rect carries ``altitude_high`` / ``altitude_low`` with
    the canonical [high, low, low, high] corner order (corners 0/3
    are HIGH, corners 1/2 are LOW).  Junction polygons MUST share
    runway/rect vertices only on the SHORT edges (corners 0-3 or
    1-2); a junction vertex on a sloping edge (0-1 or 3-2) creates
    a 4-corner-rect topology violation where the rect's continuous
    high→low slope on that side has a foreign vertex in the
    middle.  When the per-vertex altitude solver runs, the
    junction's altitude at the violating vertex disagrees with the
    rect's interpolated altitude, producing a visible step (the
    SPJC pavement-grade test catches these as mid-edge violations
    against primary_parallel V).

    For each violation, split the rect at the violating vertex's
    axial position into two sloped sub-rects.  Both sub-rects
    remain 4-corner sloped rects with their own altitude_high /
    altitude_low.  The violating vertex now coincides with a
    sub-rect's short-edge corner instead of sitting on a
    sloping edge — topology fixed.

    Per user 2026-05-12: this addresses sloping-edge violations
    introduced when the width-transition splitter
    (_find_width_transition_breakpoints) creates V parallel rects
    in regions target treats as junction territory.  The new
    rects' sloping edges then pass through adjacent junction
    polygons, picking up junction vertices.

    Returns the number of rect splits performed.
    """
    SLOPED_ROLES = (
        ROLE_STUB, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_CROSS_CONNECTOR,
    )
    on_edge_tol2 = on_edge_tol_m * on_edge_tol_m

    # Collect candidate sloped rects.
    candidates: List[int] = []
    for i, s in enumerate(layout.shapes):
        if s.role not in SLOPED_ROLES:
            continue
        if s.altitude_high is None or s.altitude_low is None:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            continue
        candidates.append(i)
    if not candidates:
        return 0

    # Collect junction vertices.  Carry (j_shape_idx, vertex_idx)
    # so we can pull the originating vertex onto the new sub-rect
    # corner after the split — without that, the rect side gets a
    # new corner at the projection of the junction vertex onto the
    # rect's long edge, but the junction vertex itself stays put
    # ~ ``on_edge_tol_m`` perpendicular off the new sub-rect's edge
    # (Tests 2 & 3 then flag the same geometry from two sides).
    j_verts: List[Tuple[float, float, int, int]] = []
    for s_idx, s in enumerate(layout.shapes):
        if s.role != ROLE_JUNCTION:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            jc = list(s.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        if jc and jc[0] == jc[-1]:
            jc = jc[:-1]
        for v_idx, (jx, jy) in enumerate(jc):
            j_verts.append((float(jx), float(jy), s_idx, v_idx))
    if not j_verts:
        return 0

    def _proj_t(pt, a, b):
        ax, ay = a; bx, by = b
        dx = bx - ax; dy = by - ay
        L2 = dx * dx + dy * dy
        if L2 < 1e-6:
            return None, None
        t = ((pt[0] - ax) * dx + (pt[1] - ay) * dy) / L2
        cx = ax + t * dx; cy = ay + t * dy
        d2 = (pt[0] - cx) ** 2 + (pt[1] - cy) ** 2
        return t, d2

    n_splits = 0
    # shape_idx -> [(cluster_t, [(j_idx, v_idx, jx, jy), ...]), ...]
    splits: Dict[int, List[Tuple[float,
                                  List[Tuple[int, int, float, float]]]]] = {}
    for idx in candidates:
        s = layout.shapes[idx]
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        c0, c1, c2, c3 = coords
        # Sloping edge 1: c0 (high) -> c1 (low).
        # Sloping edge 2: c3 (high) -> c2 (low).
        # For each junction vertex (not at a corner), check distance
        # to each sloping edge.  Collect (t, j_idx, v_idx, jx, jy)
        # so we can move the originating vertex to the new sub-rect
        # corner after splitting.
        rect_corner_set = {(c[0], c[1]) for c in coords}
        ts: List[Tuple[float, int, int, float, float]] = []
        for jx, jy, j_idx, v_idx in j_verts:
            if (jx, jy) in rect_corner_set:
                continue
            for a, b in [(c0, c1), (c3, c2)]:
                t, d2 = _proj_t((jx, jy), a, b)
                if t is None or d2 is None:
                    continue
                if d2 > on_edge_tol2:
                    continue
                if t < min_t_margin or t > 1.0 - min_t_margin:
                    continue
                ts.append((t, j_idx, v_idx, jx, jy))
                break
        if ts:
            # Cluster nearby t values (< 0.05 apart).  Every junction
            # vertex inside a cluster snaps to the same new sub-rect
            # corner (representative t of the cluster).
            ts.sort(key=lambda e: e[0])
            clusters: List[Tuple[float,
                                  List[Tuple[int, int, float, float]]]] = []
            for (t, j_idx, v_idx, jx, jy) in ts:
                if clusters and t - clusters[-1][0] <= 0.05:
                    clusters[-1][1].append((j_idx, v_idx, jx, jy))
                else:
                    clusters.append((t, [(j_idx, v_idx, jx, jy)]))
            splits[idx] = clusters

    if not splits:
        return 0

    # Apply splits.  For each candidate, create N+1 sub-rects from
    # N split positions, then snap each source junction vertex to
    # the new sub-rect corner at its cluster's t.  ``junction_moves``
    # collects per-junction vertex updates and is applied once per
    # junction shape so node_altitudes alignment is preserved.
    new_shapes: List["BuiltShape"] = []
    drop_idxs: set = set()
    # junction_shape_idx -> {vertex_idx -> new_xy}
    junction_moves: Dict[int, Dict[int, Tuple[float, float]]] = {}
    for idx, clusters in splits.items():
        s = layout.shapes[idx]
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        c0, c1, c2, c3 = [tuple(c) for c in coords]
        alt_h = float(s.altitude_high)
        alt_l = float(s.altitude_low)

        def _interp_pt(a, b, t):
            return (a[0] + t * (b[0] - a[0]),
                    a[1] + t * (b[1] - a[1]))

        ts_only = [c[0] for c in clusters]
        boundaries = [0.0] + ts_only + [1.0]
        for k in range(len(boundaries) - 1):
            t_a = boundaries[k]
            t_b = boundaries[k + 1]
            if t_b - t_a < 0.02:
                continue
            # Corners of this sub-rect: high at t_a, low at t_b.
            # Route every interpolated corner through the canonical
            # registry so it shares (x, y) with any junction vertex
            # whose snap target lands on the same axial fraction
            # below (``junction_moves`` step) — and so future
            # passes see the same canonical coordinates.
            _registry = getattr(layout, "canonical_points", None)
            def _register(pt):
                if _registry is None:
                    return pt
                return _registry.get_or_add(pt[0], pt[1])
            sub_c0 = _register(_interp_pt(c0, c1, t_a))
            sub_c1 = _register(_interp_pt(c0, c1, t_b))
            sub_c2 = _register(_interp_pt(c3, c2, t_b))
            sub_c3 = _register(_interp_pt(c3, c2, t_a))
            from shapely.geometry import Polygon
            try:
                poly = Polygon([sub_c0, sub_c1, sub_c2, sub_c3,
                                 sub_c0])
                if not poly.is_valid:
                    poly = poly.buffer(0)
                if poly.is_empty or poly.geom_type != "Polygon":
                    continue
            except _GEOM_EXC:
                continue
            import copy
            new_s = copy.copy(s)
            new_s.polygon = poly
            new_s.altitude_high = round(alt_h + t_a * (alt_l - alt_h), 1)
            new_s.altitude_low = round(alt_h + t_b * (alt_l - alt_h), 1)
            new_s.node_altitudes = None
            # If new sub-rect is effectively flat (delta < 0.1m),
            # convert to flat altitude.
            if abs(new_s.altitude_high - new_s.altitude_low) < 0.1:
                avg = 0.5 * (new_s.altitude_high + new_s.altitude_low)
                new_s.altitude = round(avg, 1)
                new_s.altitude_high = None
                new_s.altitude_low = None
            new_shapes.append(new_s)
            n_splits += 1
        drop_idxs.add(idx)

        # For each cluster, also move the originating junction
        # vertices onto the new sub-rect corner — whichever of the
        # two long-edge corners at parameter t is closer to each
        # source vertex.  Both corners go through the canonical
        # registry first so junction-vertex moves land on the
        # SAME canonical (x, y) that ``new_s.polygon`` already
        # uses for its corners — shared-vertex by construction.
        for cluster_t, v_data_list in clusters:
            corner_a_raw = _interp_pt(c0, c1, cluster_t)
            corner_b_raw = _interp_pt(c3, c2, cluster_t)
            if _registry is not None:
                corner_a = _registry.get_or_add(
                    corner_a_raw[0], corner_a_raw[1])
                corner_b = _registry.get_or_add(
                    corner_b_raw[0], corner_b_raw[1])
            else:
                corner_a = corner_a_raw
                corner_b = corner_b_raw
            for (j_idx, v_idx, jx, jy) in v_data_list:
                da = (jx - corner_a[0]) ** 2 + (jy - corner_a[1]) ** 2
                db = (jx - corner_b[0]) ** 2 + (jy - corner_b[1]) ** 2
                target = corner_a if da <= db else corner_b
                # Last-write-wins if two rects' splits both name the
                # same (j_idx, v_idx).  In practice each junction
                # vertex is near at most one rect's long edge.
                junction_moves.setdefault(j_idx, {})[v_idx] = target

    if not new_shapes:
        return 0

    # Apply junction-vertex moves before swapping shapes, so we can
    # still index layout.shapes by the original j_idx values.
    n_jct_moves = 0
    for j_idx, vert_updates in junction_moves.items():
        if not (0 <= j_idx < len(layout.shapes)):
            continue
        j_shape = layout.shapes[j_idx]
        if j_shape.polygon is None or j_shape.polygon.is_empty:
            continue
        try:
            jc = list(j_shape.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        had_close = bool(jc) and jc[0] == jc[-1]
        if had_close:
            jc = jc[:-1]
        new_jc = list(jc)
        for v_idx, new_xy in vert_updates.items():
            if 0 <= v_idx < len(new_jc):
                new_jc[v_idx] = new_xy
        if had_close:
            new_jc = new_jc + [new_jc[0]]
        from shapely.geometry import Polygon
        try:
            new_poly = Polygon(new_jc)
            if not new_poly.is_valid:
                new_poly = new_poly.buffer(0)
            if new_poly.is_empty:
                continue
            if new_poly.geom_type == "MultiPolygon":
                new_poly = max(new_poly.geoms, key=lambda g: g.area)
            # Route the (possibly ``buffer(0)``-repaired) modified
            # junction perimeter through the canonical registry so
            # any drift from the Shapely repair resolves back to
            # canonical coordinates shared with neighbours.
            _reg2 = getattr(layout, "canonical_points", None)
            if _reg2 is not None:
                new_poly = snap_polygon_through_registry(
                    new_poly, _reg2)
                if (new_poly is None or new_poly.is_empty
                        or new_poly.geom_type != "Polygon"):
                    continue
            if new_poly.geom_type != "Polygon":
                continue
        except _GEOM_EXC:
            continue
        j_shape.polygon = new_poly
        n_jct_moves += len(vert_updates)

    layout.shapes = [
        s for k, s in enumerate(layout.shapes)
        if k not in drop_idxs
    ] + new_shapes
    try:
        msg = (f"  [pav-builder] {icao}: split "
               f"{len(splits)} sloped rect(s) into "
               f"{n_splits} sub-rect(s) at junction-vertex "
               f"sloping-edge violations")
        if n_jct_moves:
            msg += (f"; pulled {n_jct_moves} junction vertex(es) "
                    f"to new corner")
        UI.vprint(1, msg + ".")
    except _GEOM_EXC:
        pass
    return n_splits




def _absorb_rect_into_junction(rect_shape, junction_shape):
    """Extend ``junction_shape``'s polygon to include
    ``rect_shape``'s footprint by replacing the shared sloping-edge
    segment in the junction's perimeter with a walk via the rect's
    far corners.

    Returns True if absorption succeeded; False if the rect doesn't
    share a sloping-edge with this junction or the resulting polygon
    is invalid.
    """
    from .layout import SHARED_VERTEX_TOL_M
    rc = list(rect_shape.polygon.exterior.coords)
    if rc and rc[0] == rc[-1]:
        rc = rc[:-1]
    if len(rc) != 4:
        return False
    jc = list(junction_shape.polygon.exterior.coords)
    had_close = bool(jc) and jc[0] == jc[-1]
    if had_close:
        jc = jc[:-1]
    n = len(jc)
    if n < 3:
        return False
    bucket = SHARED_VERTEX_TOL_M

    def _key(p):
        return (round(p[0] / bucket), round(p[1] / bucket))
    rect_keys = [_key(c) for c in rc]

    j_alts = junction_shape.node_altitudes
    j_alts_open = (list(j_alts[:n])
                   if j_alts is not None and len(j_alts) >= n
                   else None)

    for i in range(n):
        ki = _key(jc[i])
        ki1 = _key(jc[(i + 1) % n])
        for c_a, c_b in ((0, 1), (2, 3)):
            ka, kb = rect_keys[c_a], rect_keys[c_b]
            forward = (ki == ka and ki1 == kb)
            reverse = (ki == kb and ki1 == ka)
            if not (forward or reverse):
                continue
            # Walk the OTHER sloping edge's corners.  For shared
            # edge (0,1), the rect's other side goes 0 → 3 → 2 → 1.
            # For (2,3), it goes 2 → 1 → 0 → 3.
            if (c_a, c_b) == (0, 1):
                other = [rc[3], rc[2]]
            else:
                other = [rc[1], rc[0]]
            if reverse:
                other = list(reversed(other))
            # Build candidate junction with rect's far corners
            # inserted between the shared-edge endpoints.
            new_jc = jc[:i + 1] + other + jc[i + 1:]
            try:
                new_poly = Polygon(new_jc)
                if not new_poly.is_valid:
                    new_poly = new_poly.buffer(0)
                if new_poly.is_empty:
                    continue
                if new_poly.geom_type == "MultiPolygon":
                    new_poly = max(new_poly.geoms,
                                    key=lambda g: g.area)
                if new_poly.geom_type != "Polygon":
                    continue
                if not new_poly.is_valid or not new_poly.is_simple:
                    continue
            except _GEOM_EXC:
                continue
            # Sanity: new polygon should have GROWN (we added the
            # rect's footprint).
            if new_poly.area < junction_shape.polygon.area:
                continue
            # Commit.  Update node_altitudes by inserting the rect's
            # corner altitudes (linearly interpolated from the
            # rect's slope: corners 0,3 = high, corners 1,2 = low).
            junction_shape.polygon = new_poly
            if j_alts_open is not None:
                # Pick altitudes for the 2 inserted corners.  For
                # the OTHER-side walk we inserted [rc[3], rc[2]] or
                # similar.  Use altitude_high for corners 0,3 and
                # altitude_low for corners 1,2 per the rect's
                # canonical slope convention.
                ah = rect_shape.altitude_high
                al = rect_shape.altitude_low
                if ah is None or al is None:
                    a_default = (rect_shape.altitude
                                 if rect_shape.altitude is not None
                                 else j_alts_open[i])
                    ah = al = float(a_default)
                hilo_for_corner = {0: ah, 1: al, 2: al, 3: ah}
                if (c_a, c_b) == (0, 1):
                    pair = [3, 2]
                else:
                    pair = [1, 0]
                if reverse:
                    pair = list(reversed(pair))
                a_insert = [hilo_for_corner[p] for p in pair]
                new_alts_open = (j_alts_open[:i + 1]
                                 + a_insert
                                 + j_alts_open[i + 1:])
                junction_shape.node_altitudes = (
                    new_alts_open + [new_alts_open[0]])
            return True
    return False


def _drop_rects_with_shared_sloping_edge_and_absorb(
        layout: "PavementLayout",
        icao: str = "",
        ) -> int:
    """Detect sloping rects with a FULLY-SHARED sloping edge with a
    junction polygon's perimeter (the
    ``test_taxi_rects_not_alongside_apron`` invariant) and absorb
    each into its adjacent junction by extending the junction
    polygon to cover the rect's footprint.

    The rect is then dropped.  Per user 2026-05-16: a sloping
    rect's sloping edge must never be shared with a junction.

    "Fully shared" means the junction perimeter has two CONSECUTIVE
    vertices coinciding with the rect's two sloping-edge corners
    (i.e. the junction walks the rect's whole sloping edge as a
    single perimeter segment).  This catches the
    ``_split_sloped_rects_at_violations`` output where a parent
    rect is carved into many short sub-rects whose sloping edges
    are entirely the boundary with an adjacent junction.

    Partial-sharing cases (junction perimeter walks along part of
    the rect's sloping edge but not corner-to-corner) are NOT
    handled here — those need a separate clipping pass, since
    dropping the whole rect would lose meaningful pavement.

    Runs as a post-pass after ``_split_sloped_rects_at_violations``.

    Returns the number of rects absorbed.
    """
    from .layout import SHARED_VERTEX_TOL_M
    sloping = (ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                ROLE_STUB, ROLE_CROSS_CONNECTOR)
    bucket = SHARED_VERTEX_TOL_M

    def _key(p):
        return (round(p[0] / bucket), round(p[1] / bucket))

    drop_indices: set = set()
    for r_idx, r in enumerate(layout.shapes):
        if r.role not in sloping:
            continue
        if r.polygon is None or r.polygon.is_empty:
            continue
        rc = list(r.polygon.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        rect_keys = [_key(c) for c in rc]
        # Look for any junction whose perimeter has 2 consecutive
        # vertices matching a sloping-edge pair (0-1) or (2-3) of
        # this rect.  Absorb into the first such junction found.
        for j_idx, j in enumerate(layout.shapes):
            if j is r or j.role != ROLE_JUNCTION:
                continue
            if j.polygon is None or j.polygon.is_empty:
                continue
            jc = list(j.polygon.exterior.coords)
            if jc and jc[0] == jc[-1]:
                jc = jc[:-1]
            n = len(jc)
            if n < 3:
                continue
            absorbed = False
            for k in range(n):
                ki = _key(jc[k])
                ki1 = _key(jc[(k + 1) % n])
                for c_a, c_b in ((0, 1), (2, 3)):
                    ka, kb = rect_keys[c_a], rect_keys[c_b]
                    if ((ki == ka and ki1 == kb)
                            or (ki == kb and ki1 == ka)):
                        if _absorb_rect_into_junction(r, j):
                            drop_indices.add(r_idx)
                            absorbed = True
                        break
                if absorbed:
                    break
            if absorbed:
                break

    if not drop_indices:
        return 0

    layout.shapes = [s for k, s in enumerate(layout.shapes)
                     if k not in drop_indices]
    try:
        UI.vprint(1,
            f"  [pav-builder] {icao}: absorbed "
            f"{len(drop_indices)} sloping rect(s) into adjacent "
            f"junction (shared sloping edge).")
    except _GEOM_EXC:
        pass
    return len(drop_indices)


def _absorb_rects_at_junction_perimeters(
        layout: "PavementLayout",
        icao: str = "",
        perp_tol_m: float = 0.5,
        min_overlap_m: float = 5.0,
        min_kept_m: float = 20.0,
        ) -> int:
    """Single-pass sloping-edge absorption using actual junction
    polygons at the FINAL pipeline state.

    Per user 2026-05-17 architectural directive: run absorption
    ONCE at the very end of the pipeline (after all post-elevation
    junction-refinement passes — split_sloped_rects_at_violations,
    snap_junction_vertices_to_rect_flat_edge_corners, widening,
    etc.).  Earlier placement (right after emit_junctions_and_finalize)
    sees transient junction polygons that haven't been refined yet
    — many SPJC corridor stubs (B, C, E) appear to have shared
    sloping edges with the raw apron residue at that point, but
    the later passes retract the apron so the stubs end up grass-
    bordered.  Running absorption at the very end avoids those
    false positives.

    For each sloping rect:
      1. For each of the rect's 2 sloping edges, find adjacent
         junction polygons whose perimeter walks within
         ``perp_tol_m`` perpendicular of the edge.
      2. Compute the shared-edge LineString(s); project onto the
         edge to get axial t-ranges (parameters 0..1).  Clamp away
         the corner-wrap zone (t < 0.05 or t > 0.95) — junctions
         wrapping at the rect's short-edge corner pick up tiny
         stretches of sloping edge that aren't actual sharing.
      3. Merge all absorbed t-ranges across both sloping edges.
      4. Kept t-ranges = [0, 1] − absorbed.  Filter kept ranges by
         ``min_kept_m`` axial length; ranges below the floor are
         folded into the absorbed set (the rect is too short to
         keep that fragment).
      5. For each absorbed range: build the rect's strip polygon
         (full width × axial range), union it into the absorbing
         junction's polygon.  Drop the junction's node_altitudes
         so the OSM emitter falls back to per-shape altitude
         (rect altitudes were committed by the per-surface solver
         and the new union polygon's vertex set doesn't match the
         old altitude array; cleaner to drop than to remap).
      6. For each kept range: build a new sub-rect polygon and
         add it as a fresh ``BuiltShape``.
      7. Drop the original rect.

    Returns the number of rects modified (clipped or dropped).
    """
    from .layout import (BuiltShape, SHARED_VERTEX_TOL_M)
    SLOPING_ROLES = (ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                      ROLE_STUB, ROLE_CROSS_CONNECTOR)

    junction_idx_list = [i for i, s in enumerate(layout.shapes)
                          if s.role == ROLE_JUNCTION
                          and s.polygon is not None
                          and not s.polygon.is_empty
                          and s.polygon.geom_type == "Polygon"]
    if not junction_idx_list:
        return 0

    drop_indices: set = set()
    new_rect_shapes: List["BuiltShape"] = []
    n_dropped = 0
    n_clipped = 0
    junction_extensions: Dict[int, List[Polygon]] = {}

    for r_idx, r in enumerate(layout.shapes):
        if r.role not in SLOPING_ROLES:
            continue
        if r.polygon is None or r.polygon.is_empty:
            continue
        if r.polygon.geom_type != "Polygon":
            continue
        rc = list(r.polygon.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                  0.5 * (rc[0][1] + rc[3][1]))
        b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                  0.5 * (rc[1][1] + rc[2][1]))
        adx = b_mid[0] - a_mid[0]
        ady = b_mid[1] - a_mid[1]
        axis_L = math.hypot(adx, ady)
        if axis_L < 1.0:
            continue
        ux, uy = adx / axis_L, ady / axis_L
        nx, ny = -uy, ux
        half_w = math.hypot(rc[0][0] - a_mid[0],
                              rc[0][1] - a_mid[1])
        if half_w < 1.0:
            continue
        sloping_edges = [
            LineString([rc[0], rc[1]]),
            LineString([rc[2], rc[3]]),
        ]
        absorbed: List[Tuple[float, float, int, int]] = []
        for edge_idx, edge_ls in enumerate(sloping_edges):
            edge_len = edge_ls.length
            if edge_len < 1.0:
                continue
            try:
                edge_buf = edge_ls.buffer(perp_tol_m, cap_style=2)
            except _GEOM_EXC:
                continue
            for j_idx in junction_idx_list:
                if j_idx in drop_indices:
                    continue
                j = layout.shapes[j_idx]
                try:
                    shared = j.polygon.boundary.intersection(edge_buf)
                except _GEOM_EXC:
                    continue
                if shared.is_empty:
                    continue
                if shared.geom_type == "LineString":
                    parts = [shared]
                elif shared.geom_type == "MultiLineString":
                    parts = list(shared.geoms)
                elif shared.geom_type == "GeometryCollection":
                    parts = [g for g in shared.geoms
                             if g.geom_type == "LineString"]
                else:
                    continue
                for part in parts:
                    if part.length < min_overlap_m:
                        continue
                    try:
                        p_a = Point(part.coords[0])
                        p_b = Point(part.coords[-1])
                        t_a = edge_ls.project(p_a) / edge_len
                        t_b = edge_ls.project(p_b) / edge_len
                    except _GEOM_EXC:
                        continue
                    t_lo = max(0.0, min(t_a, t_b))
                    t_hi = min(1.0, max(t_a, t_b))
                    # Skip corner-wrap zones: junctions adjacent to
                    # the rect at its SHORT edges pick up tiny
                    # stretches of sloping edge right at the
                    # corners.  Require the shared range to extend
                    # past 5% from each corner into the edge
                    # interior.
                    INNER_T_MIN = 0.05
                    INNER_T_MAX = 0.95
                    t_lo_inner = max(INNER_T_MIN, t_lo)
                    t_hi_inner = min(INNER_T_MAX, t_hi)
                    if t_hi_inner <= t_lo_inner:
                        continue
                    inner_len = (t_hi_inner - t_lo_inner) * edge_len
                    if inner_len < min_overlap_m:
                        continue
                    absorbed.append((t_lo, t_hi, j_idx, edge_idx))
        if not absorbed:
            continue
        intervals = sorted([(s, e) for s, e, _j, _ei in absorbed])
        merged: List[Tuple[float, float]] = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0],
                                max(merged[-1][1], e))
            else:
                merged.append((s, e))
        kept: List[Tuple[float, float]] = []
        prev_end = 0.0
        for s, e in merged:
            if s > prev_end:
                kept.append((prev_end, s))
            prev_end = e
        if prev_end < 1.0:
            kept.append((prev_end, 1.0))
        kept_keep = [(s, e) for s, e in kept
                      if (e - s) * axis_L >= min_kept_m]
        if len(kept_keep) == len(kept):
            kept_final = kept_keep
            absorbed_final = merged
        else:
            kept_final = kept_keep
            absorbed_final = []
            prev = 0.0
            for s, e in sorted(kept_keep):
                if s > prev:
                    absorbed_final.append((prev, s))
                prev = e
            if prev < 1.0:
                absorbed_final.append((prev, 1.0))
        if not absorbed_final:
            continue
        def _strip_polygon(t_lo, t_hi):
            # Full-absorption: reuse the rect's exact polygon so
            # the union preserves the rect's corner positions (the
            # neighbour rects that shared those corners via short
            # edges keep their connection through the extended
            # junction's perimeter).
            if t_lo <= 1e-6 and t_hi >= 1.0 - 1e-6:
                return r.polygon
            u_lo = t_lo * axis_L
            u_hi = t_hi * axis_L
            # For partial absorption, use the rect's original
            # corner at boundaries that ARE at t=0 or t=1; only
            # use axis-derived positions for the new mid-rect cut.
            if t_lo <= 1e-6:
                p0 = rc[0]
                p3 = rc[3]
            else:
                p0 = (a_mid[0] + u_lo * ux + nx * half_w,
                      a_mid[1] + u_lo * uy + ny * half_w)
                p3 = (a_mid[0] + u_lo * ux - nx * half_w,
                      a_mid[1] + u_lo * uy - ny * half_w)
            if t_hi >= 1.0 - 1e-6:
                p1 = rc[1]
                p2 = rc[2]
            else:
                p1 = (a_mid[0] + u_hi * ux + nx * half_w,
                      a_mid[1] + u_hi * uy + ny * half_w)
                p2 = (a_mid[0] + u_hi * ux - nx * half_w,
                      a_mid[1] + u_hi * uy - ny * half_w)
            try:
                strip = Polygon([p0, p1, p2, p3])
                if not strip.is_valid:
                    strip = strip.buffer(0)
                if strip.is_empty or strip.geom_type != "Polygon":
                    return None
                return strip
            except _GEOM_EXC:
                return None
        for t_s, t_e in absorbed_final:
            strip = _strip_polygon(t_s, t_e)
            if strip is None or strip.is_empty:
                continue
            assigned = False
            for orig_s, orig_e, j_idx, _ei in absorbed:
                if orig_s < t_e and orig_e > t_s:
                    junction_extensions.setdefault(
                        j_idx, []).append(strip)
                    assigned = True
            if not assigned:
                strip_c = strip.centroid
                best = None
                best_d2 = float("inf")
                for j_idx in junction_idx_list:
                    jc = layout.shapes[j_idx].polygon.centroid
                    d2 = ((jc.x - strip_c.x) ** 2
                          + (jc.y - strip_c.y) ** 2)
                    if d2 < best_d2:
                        best_d2 = d2
                        best = j_idx
                if best is not None:
                    junction_extensions.setdefault(
                        best, []).append(strip)
        # Compute per-axial-range altitudes for new sub-rects from
        # the parent rect's altitude_high / altitude_low.  Corners
        # 0,3 → altitude_high; corners 1,2 → altitude_low (per the
        # _rect_from_axis_extended convention).
        parent_ah = r.altitude_high
        parent_al = r.altitude_low
        parent_alt = r.altitude
        def _alt_at_t(t):
            if parent_ah is not None and parent_al is not None:
                return parent_ah + t * (parent_al - parent_ah)
            if parent_alt is not None:
                return parent_alt
            return None
        for t_s, t_e in kept_final:
            new_rect = _strip_polygon(t_s, t_e)
            if new_rect is None or new_rect.is_empty:
                continue
            new_axis = LineString([
                (a_mid[0] + t_s * axis_L * ux,
                 a_mid[1] + t_s * axis_L * uy),
                (a_mid[0] + t_e * axis_L * ux,
                 a_mid[1] + t_e * axis_L * uy),
            ])
            import copy
            new_s = copy.copy(r)
            new_s.polygon = new_rect
            new_s.source_axis = new_axis
            # Re-derive altitudes for the kept sub-rect from the
            # parent rect's slope.  Sub-rect HIGH end is at the
            # smaller t (closer to corner 0/3, parent HIGH side).
            a_high_t = min(t_s, t_e)
            a_low_t = max(t_s, t_e)
            ah_new = _alt_at_t(a_high_t)
            al_new = _alt_at_t(a_low_t)
            if (ah_new is not None and al_new is not None
                    and abs(ah_new - al_new) >= 0.1):
                new_s.altitude_high = round(ah_new, 1)
                new_s.altitude_low = round(al_new, 1)
                new_s.altitude = None
            elif ah_new is not None:
                new_s.altitude = round(ah_new, 1)
                new_s.altitude_high = None
                new_s.altitude_low = None
            new_rect_shapes.append(new_s)
        drop_indices.add(r_idx)
        if kept_final:
            n_clipped += 1
        else:
            n_dropped += 1

    if not drop_indices and not junction_extensions:
        return 0

    # Apply junction extensions.  After union, the junction
    # polygon's vertex set changes; we need to re-derive
    # node_altitudes for the new vertex set.  Build a lookup
    # table of (vertex_position → altitude) from the OLD junction
    # corners (their committed altitudes) AND from the absorbed
    # rect strip corners (altitudes derived from rect axial
    # position).  Each new polygon vertex gets its altitude by
    # nearest-known-vertex match.
    from .layout import SHARED_VERTEX_TOL_M as _SVT
    bucket = _SVT
    for j_idx, strip_data in junction_extensions.items():
        j = layout.shapes[j_idx]
        if j.polygon is None or j.polygon.is_empty:
            continue
        # Build altitude lookup BEFORE replacing polygon.
        alt_lookup: Dict[Tuple[int, int], float] = {}
        def _bkey(x, y):
            return (round(x / bucket), round(y / bucket))
        try:
            old_coords = list(j.polygon.exterior.coords)
            if old_coords and old_coords[0] == old_coords[-1]:
                old_coords = old_coords[:-1]
        except _GEOM_EXC:
            old_coords = []
        if j.node_altitudes and len(j.node_altitudes) >= len(old_coords):
            for i, (cx, cy) in enumerate(old_coords):
                if i < len(j.node_altitudes):
                    alt_lookup[_bkey(cx, cy)] = float(
                        j.node_altitudes[i])
        elif j.altitude is not None:
            for cx, cy in old_coords:
                alt_lookup[_bkey(cx, cy)] = float(j.altitude)
        # Strips are simple list of polygons here.  Each strip
        # came from a rect with known altitudes — register its
        # corners (corners 0,3 = high; corners 1,2 = low per the
        # rect convention used in _strip_polygon).
        # Looking at _strip_polygon: it builds p0 (lo*nx+ny), p1
        # (hi*nx+ny), p2 (hi*-nx-ny), p3 (lo*-nx-ny) — so p0/p3
        # are at t=t_lo (rect's high side), p1/p2 at t=t_hi (low).
        # We tracked (strip, rect_alt_high, rect_alt_low, t_lo, t_hi)
        # in the strip_data tuples.  For now, use a simpler approach:
        # for each strip corner, infer altitude from interpolation
        # along the junction's existing altitude field.
        try:
            merged_strips = unary_union(strip_data)
            new_j_poly = unary_union([j.polygon, merged_strips])
            if new_j_poly.geom_type == "MultiPolygon":
                new_j_poly = max(new_j_poly.geoms,
                                  key=lambda g: g.area)
            if not (new_j_poly.geom_type == "Polygon"
                    and new_j_poly.is_valid
                    and not new_j_poly.is_empty):
                continue
            # Route the unary_union output through the canonical
            # registry: ``unary_union`` reorders / re-rounds
            # vertices during merging and the result's perimeter
            # vertices can drift mm-m off the input rect / junction
            # corners that were originally registered.  Routing
            # restores exact canonical coordinates so the extended
            # junction still shares vertices with its neighbours.
            _registry = getattr(layout, "canonical_points", None)
            if _registry is not None:
                new_j_poly = snap_polygon_through_registry(
                    new_j_poly, _registry)
                if (new_j_poly is None or new_j_poly.is_empty
                        or new_j_poly.geom_type != "Polygon"):
                    continue
            # Also register strip corners in alt_lookup.  Use the
            # NEAREST OLD junction vertex's altitude (a reasonable
            # approximation since the strip is adjacent to the
            # junction).
            for strip in strip_data:
                try:
                    s_coords = list(strip.exterior.coords)
                    if s_coords and s_coords[0] == s_coords[-1]:
                        s_coords = s_coords[:-1]
                except _GEOM_EXC:
                    continue
                for sx, sy in s_coords:
                    key = _bkey(sx, sy)
                    if key in alt_lookup:
                        continue
                    # Find nearest old corner altitude.
                    if old_coords:
                        best_d2 = float("inf")
                        best_alt = None
                        for cx, cy in old_coords:
                            d2 = (cx - sx) ** 2 + (cy - sy) ** 2
                            if d2 < best_d2:
                                best_d2 = d2
                                k = _bkey(cx, cy)
                                if k in alt_lookup:
                                    best_alt = alt_lookup[k]
                        if best_alt is not None:
                            alt_lookup[key] = best_alt
            j.polygon = new_j_poly
            # Assign altitudes to the new polygon's vertices.
            new_coords = list(new_j_poly.exterior.coords)
            had_close = bool(new_coords) and new_coords[0] == new_coords[-1]
            if had_close:
                new_coords = new_coords[:-1]
            if alt_lookup:
                new_alts = []
                for cx, cy in new_coords:
                    key = _bkey(cx, cy)
                    if key in alt_lookup:
                        new_alts.append(alt_lookup[key])
                        continue
                    # Nearest known corner.
                    best_d2 = float("inf")
                    best_alt = 0.0
                    for (kx, ky), alt in alt_lookup.items():
                        # Recover x,y from the bucket key.
                        d2 = ((kx * bucket - cx) ** 2
                              + (ky * bucket - cy) ** 2)
                        if d2 < best_d2:
                            best_d2 = d2
                            best_alt = alt
                    new_alts.append(round(best_alt, 1))
                j.node_altitudes = new_alts + [new_alts[0]]
                j.altitude = None
        except _GEOM_EXC:
            continue

    if drop_indices:
        layout.shapes = [s for k, s in enumerate(layout.shapes)
                          if k not in drop_indices]
    if new_rect_shapes:
        layout.shapes.extend(new_rect_shapes)

    n_total = n_dropped + n_clipped
    try:
        UI.vprint(1,
            f"  [pav-builder] {icao}: end-of-pipeline "
            f"absorption — {n_dropped} dropped, "
            f"{n_clipped} clipped, "
            f"{len(junction_extensions)} junction(s) extended.")
    except _GEOM_EXC:
        pass
    return n_total


# ── Apron reclassification ────────────────────────────────────────

# Per user 2026-05-18: a junction whose boundary strays >
# ``_APRON_RECLASSIFY_MAX_DISTANCE_M`` from the nearest taxi /
# runway centerline contains apron-territory pavement and should
# be reclassified as ``role=apron``.
#
# Threshold = 55 m.  Taxiway half-width is ~11.5 m and a normal
# 90° fillet apex sits ~22 m from the centerline; large airport
# fillets (SPJC mega-junctions, runway-taxi turn-off filleted to
# match wide-body taxi paths) can extend up to ~50 m from the
# straight-line apt.dat centerline.  55 m gives a small margin
# above that while still clearly catching apron territory (which
# extends 100 m+ from any centerline).
_APRON_RECLASSIFY_MAX_DISTANCE_M = 55.0
_APRON_RECLASSIFY_SAMPLE_STEP_M = 5.0


def _aeroway_centerlines_union(layout: "PavementLayout"):
    """Union of every taxi + runway centerline known to the layout.

    Sources, in priority order:

    * ``layout.apt_taxi_centerlines`` — the full apt.dat / OSM
      taxi network captured BEFORE rect-and-junction decomposition.
      Centerlines that got absorbed into junction polygons survive
      here even though they're no longer on any shape.  Without
      this, a legitimate junction sitting on top of an absorbed
      centerline would test as "no centerline within 20 m" and
      get misflagged as apron.
    * Surviving taxi rects' ``source_axis``.  Redundant with the
      apt.dat set for most cases but catches OSM-only-derived
      centerlines (custom packs without apt.dat row 1201/1202).
    * Runway long-axes derived from each runway segment's 4 corners.
    """
    lines = []
    apt_lines = getattr(layout, "apt_taxi_centerlines", None) or []
    for item in apt_lines:
        # apt.dat / OSM taxi extraction returns
        # ``(LineString, name)`` tuples.
        ln = item[0] if isinstance(item, tuple) else item
        if ln is not None and not ln.is_empty:
            lines.append(ln)
    for s in layout.shapes:
        if s.source_axis is not None and not s.source_axis.is_empty:
            lines.append(s.source_axis)
            continue
        if (s.role == ROLE_RUNWAY
                and s.polygon is not None
                and not s.polygon.is_empty):
            try:
                coords = list(s.polygon.exterior.coords)
            except _GEOM_EXC:
                continue
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
    except _GEOM_EXC:
        return None


def _reclassify_apron_junctions(
        layout: "PavementLayout",
        icao: str = "",
        cap_m: float = _APRON_RECLASSIFY_MAX_DISTANCE_M,
        sample_step_m: float = _APRON_RECLASSIFY_SAMPLE_STEP_M,
        ) -> int:
    """Reclassify any junction whose boundary strays > ``cap_m``
    from the nearest taxi/runway centerline as ``role=apron``.

    The reclassification is geometric, not area-based: a 6-way
    mega-intersection is a valid junction even when it's large,
    but a pavement region without a centerline running through it
    is apron territory regardless of size.

    Centerline source is ``_aeroway_centerlines_union(layout)``,
    which now includes the original apt.dat taxi-network lines
    preserved on ``layout.apt_taxi_centerlines`` — without those,
    centerlines absorbed into junction polygons during the rect /
    junction decomposition would be missing from the union and
    legitimate junctions would fail the boundary-distance test.

    Returns the count of reclassified shapes.
    """
    centers = _aeroway_centerlines_union(layout)
    if centers is None or centers.is_empty:
        return 0
    n_reclassified = 0
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            bnd = s.polygon.boundary
            L = bnd.length
        except _GEOM_EXC:
            continue
        if L <= 0:
            continue
        n_steps = max(2, int(L / sample_step_m) + 1)
        max_d = 0.0
        for i in range(n_steps):
            u = min(L, i * sample_step_m)
            try:
                p = bnd.interpolate(u)
                d = centers.distance(p)
            except _GEOM_EXC:
                continue
            if d > max_d:
                max_d = d
                if max_d > cap_m:
                    break
        if max_d > cap_m:
            s.role = ROLE_APRON
            n_reclassified += 1
    if n_reclassified:
        try:
            UI.vprint(1,
                f"  [pav-builder] {icao}: reclassified "
                f"{n_reclassified} junction(s) as apron "
                f"(boundary > {cap_m:.0f} m from any centerline).")
        except _GEOM_EXC:
            pass
    return n_reclassified
