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

from .elevation import (
    NEIGHBOUR_CLAMP_RADIUS_M,
    TAXI_MAX_GRADE,
    _corner_elevation_bucket,
)
from .layout import (
    BuiltShape,
    PavementLayout,
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
                        # Skip ONLY edges incident to THIS vertex —
                        # they trivially equal the vertex's own
                        # elevation and would lock it in place.  All
                        # other same-shape boundary edges DO
                        # constrain it (the polygon's far-side
                        # vertices, narrow-waist near-touches, etc.).
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
            """Return (worst_grade, worst_pair_idx) over all
            vertex pairs of ``pts`` within
            SUBDIVIDE_MAX_PAIR_DIST_M and > 0.5 m apart.  The
            distance floor is critical: < 0.5 m pairs are
            essentially the same vertex with rounding noise
            and produce huge spurious grades."""
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

        for sp, sub_pts, sub_elevs in validated_subs:
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


def _try_iso_elevation_cut(
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


