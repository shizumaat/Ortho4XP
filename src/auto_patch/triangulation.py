"""Per-vertex altitude assignment for junction polygons.

Junction polygons are residue-pavement shapes (apt.dat row-110
covered area minus rect / terminal / runway shapes) that need
per-vertex altitudes so X-Plane can triangulate them with smooth
multi-directional slopes.  This module owns that work.

For every junction:
  1. Collect HARD anchors (rect corners, runway segments, terminals)
     touching the junction's boundary.
  2. Drop sliver corners + spike vertices that would produce
     near-degenerate triangles.
  3. Assign each ring vertex an altitude: planar fit when anchors
     are co-planar enough, otherwise 2D-Euclidean smoothing
     (``_smooth_polygon_grid``).
  4. Stash the per-vertex altitudes on the BuiltShape's
     ``node_altitudes``.

Public API:
    _triangulate_junctions(layout, dem, tile_lat, tile_lon)
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from shapely.geometry import Polygon

import O4_UI_Utils as UI

from .elevation import (
    NEIGHBOUR_CLAMP_RADIUS_M,
    SHARED_AGREE_TOL_M,
    TAXI_MAX_GRADE,
    USE_PER_POLYGON_ELEVATION_FIELD,
    _corner_elev_map,
    _corner_elevation_bucket,
    _match_elev,
    _planar_fit,
    _sample_dem,
)
from .elevation_smoothing import _smooth_polygon_grid
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
    SHARED_VERTEX_TOL_M,
)
from .pavement.junctions import (
    _drop_sliver_corners,
    _splice_holes,
)
from .pavement.vertices import _drop_spike_vertices


__all__ = ["_triangulate_junctions"]


def _triangulate_junctions(
    layout: "PavementLayout",
    dem,
    tile_lat: int,
    tile_lon: int,
    m_to_ll,
) -> int:
    """Replace each junction shape with its ear-clip triangulation,
    setting per-triangle ``node_altitudes`` from the corner-elevation
    map (with neighbour-corner, elevation-graph, and DEM fallbacks).

    Returns the number of triangle shapes produced (informational).
    """
    corner_elev = _corner_elev_map(layout)

    # Cross-junction shared-vertex anchoring.  When two adjacent
    # junction polygons share a boundary point (e.g. either side of
    # a decomposition cut, or both edges of a former hole), they
    # MUST agree on its elevation.  Independent per-junction
    # smoothing can otherwise drift them apart, producing a
    # vertical step at the shared edge.
    #
    # Approach: count how many junction polygons reference each
    # SHARED_VERTEX_TOL_M bucket.  Any bucket touched by ≥ 2
    # junctions is "shared"; we anchor those vertices to a
    # deterministic value (the first junction's graph-sampled or
    # corner-derived elevation) so all junctions read the same
    # value when looking up that bucket.
    junction_bucket_count: Dict[Tuple[int, int], int] = {}
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        try:
            ring = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if ring and ring[0] == ring[-1]:
            ring = ring[:-1]
        seen: set = set()
        for (x, y) in ring:
            key = _corner_elevation_bucket(x, y)
            if key in seen:
                continue
            seen.add(key)
            junction_bucket_count[key] = (
                junction_bucket_count.get(key, 0) + 1)
    shared_junction_buckets = {
        k for k, c in junction_bucket_count.items() if c >= 2}
    # Computed lazily and cached: the first junction to look up a
    # shared-bucket vertex computes its elevation; subsequent
    # junctions read the same value.
    # Two-tier cache so we can distinguish where a shared-bucket
    # value came from: rect-corner / near-corner / edge-interp
    # (HARD — same physical pavement element pinned the value, must
    # be preserved) versus graph or DEM (SOFT — sampled from a 1.5 %
    # network-distance-compliant model, but not 2D-Euclidean
    # compliant against other graph samples that may be far on the
    # network but close in 2D).  Junction-local smoothing must be
    # free to move SOFT values into 2D grade compliance; the cross-
    # junction iteration below averages SOFT shared-bucket values
    # across every junction that touches the bucket so junctions
    # still agree on a single shared elevation.
    shared_junction_elev: Dict[Tuple[int, int], float] = {}      # HARD
    shared_junction_elev_soft: Dict[Tuple[int, int], float] = {}  # SOFT

    # Flat list of (cx, cy, elev) for nearest-corner search beyond
    # the 0.5 m bucket.  Catches boundary vertices that landed close
    # to but not exactly on a neighbour corner (e.g. apt.dat
    # boundary-trace points that abut a terminal pad edge).
    rect_like_roles = {ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL,
                       ROLE_STUB, ROLE_CROSS_CONNECTOR}
    NEAR_CORNER_M = 6.0  # search radius for off-bucket corner matches
    NEAR_EDGE_M = 5.0    # search radius for rect-edge interpolation
    near_corner_list: List[Tuple[float, float, float]] = []
    # Each ``neighbour_edges`` entry is one polygon edge of a non-
    # junction shape, paired with the elevations at its two ends:
    # ``(ax, ay, bx, by, e_a, e_b)``.  Junction vertices that fall
    # close to such an edge are anchored to the linearly-interpolated
    # elevation along the edge — guaranteeing the junction triangle
    # meets a sloped rect's sloping edge at the same height the rect is
    # rendering at, rather than at the centerline-graph value (which
    # is offset by the rect's half-width).
    neighbour_edges: List[Tuple[float, float, float, float,
                                float, float]] = []
    for s in layout.shapes:
        if s.role == ROLE_JUNCTION:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        # Per-corner elevations.
        if (s.role in rect_like_roles
                and s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            elevs = [s.altitude_high, s.altitude_low,
                     s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            elevs = [float(s.altitude)] * len(coords)
        else:
            continue
        for (cx, cy), e in zip(coords, elevs):
            near_corner_list.append((cx, cy, float(e)))
        # Edge list (closed ring).
        m = len(coords)
        for i in range(m):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % m]
            ea = float(elevs[i])
            eb = float(elevs[(i + 1) % m])
            neighbour_edges.append((ax, ay, bx, by, ea, eb))

    # Sloping-rect ALL edges (user 2026-05-04 clarification: junction
    # joins with a sloping rect have to match 1:1 — no intermediate
    # nodes on EITHER the sloping edges OR the cross edges.  The
    # rect's 4 corners are the only legal shared vertices.  Pass
    # every edge of every sloping-role rect to the densify-skip
    # guard and the snap-to-corner pass.
    #
    # Per user 2026-05-04 (earlier today): we no longer skip rects
    # with unassigned altitudes here.  Under the per-surface solver
    # path, altitudes are assigned LATER — at densification time,
    # every rect has ``altitude_high is None``, so the previous skip
    # emptied this list and Rule 2 never triggered.  Treat every
    # sloping-role rect as sloping at this stage; if the rect turns
    # out flat after the elevation pass, the junction's corner-only
    # sharing is still valid (the "flat rects allow free
    # densification" rule allows extras but doesn't require them).
    sloping_rect_edges: List[Tuple[float, float, float, float]] = []
    sloping_roles = {ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                     ROLE_STUB, ROLE_CROSS_CONNECTOR}
    for s in layout.shapes:
        if s.role not in sloping_roles:
            continue
        try:
            rc = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        for i in range(4):
            ax, ay = rc[i]
            bx, by = rc[(i + 1) % 4]
            sloping_rect_edges.append(
                (float(ax), float(ay), float(bx), float(by)))

    # Runway boundary edges (Rule 1: densification midpoints must
    # not land within RUNWAY_BOUNDARY_TOL_M of a runway boundary
    # unless they coincide with a runway vertex).  Treated the same
    # way as sloping rect sloping edges via the per-edge wider-tolerance
    # exclusion.
    runway_edges: List[Tuple[float, float, float, float]] = []
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        try:
            rc = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        m = len(rc)
        for i in range(m):
            ax, ay = rc[i]
            bx, by = rc[(i + 1) % m]
            runway_edges.append((float(ax), float(ay),
                                 float(bx), float(by)))

    # Terminal pad edges (user 2026-05-04: junctions adjacent to a
    # terminal must share its boundary node-for-node, no extra mid-
    # edge vertices).  Terminals lack altitude at triangulation time
    # so they're absent from ``neighbour_edges``; we pass their
    # geometry separately as a skip-only list.
    terminal_edges: List[Tuple[float, float, float, float]] = []
    for s in layout.shapes:
        if s.role != ROLE_TERMINAL:
            continue
        try:
            rc = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        m = len(rc)
        for i in range(m):
            ax, ay = rc[i]
            bx, by = rc[(i + 1) % m]
            terminal_edges.append((float(ax), float(ay),
                                    float(bx), float(by)))

    def _edge_interp_elev(x: float, y: float,
                          max_dist: float = NEAR_EDGE_M
                          ) -> Optional[float]:
        """Return the elevation interpolated along the closest
        non-junction polygon edge within ``max_dist`` of (x, y),
        or None if no edge is within range.

        Projects the query point onto each edge's segment, clamps
        to [0, 1], and linearly interpolates between the edge's
        two endpoint elevations.  This matches the elevation X-Plane
        renders for a sloped rect at any point along its sloping edge
        — so junction triangles abutting a sloped rect at this
        point will meet it without a step.
        """
        best_e: Optional[float] = None
        best_d2 = max_dist * max_dist
        for ax, ay, bx, by, ea, eb in neighbour_edges:
            dx = bx - ax
            dy = by - ay
            seg_len2 = dx * dx + dy * dy
            if seg_len2 < 0.04:  # < 0.2 m segment, skip
                continue
            t = ((x - ax) * dx + (y - ay) * dy) / seg_len2
            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0
            cx = ax + t * dx
            cy = ay + t * dy
            d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)
            if d2 < best_d2:
                best_d2 = d2
                best_e = ea + t * (eb - ea)
        return best_e

    def _vertex_elev_anchored(x: float, y: float
                              ) -> Tuple[Optional[float], bool]:
        """Return ``(elev, is_anchor)``.  ``is_anchor`` is True only
        for HARD anchor sources: rect / runway / terminal corner
        elevations (whether matched by exact bucket or off-bucket
        near-corner / edge-interp search), or a value already cached
        from such a source on a prior shared-bucket lookup.  Returns
        is_anchor=False for graph- or DEM-sampled values so the
        boundary smoother and Mode B grid are free to pull them into
        2D Euclidean grade compliance with the real anchors —
        graph-derived values respect 1.5 % over **network distance**
        (centerline graph), which can yield arbitrarily large 2D
        steps between graph nodes that are far on the network but
        close in 2D.  Cross-junction agreement on shared SOFT
        buckets is restored by the cross-junction averaging pass in
        the outer iteration below."""
        bucket = _corner_elevation_bucket(x, y)
        # 1. Bucket lookup — exact shared-vertex match against
        # rect / runway / terminal corners.  HARD.
        e = corner_elev.get(bucket)
        if e is not None:
            return e, True
        # 2a. Cross-junction shared bucket previously locked from a
        # HARD source (rect corner / near-corner / edge-interp).
        e_shared = shared_junction_elev.get(bucket)
        if e_shared is not None:
            return e_shared, True
        # 2b. Cross-junction shared bucket previously sampled from a
        # SOFT source (graph / DEM).  Reuse the value so this call
        # returns deterministically consistent with the first
        # junction that touched the bucket — but mark it FREE so
        # smoothing can still move it.
        e_shared_soft = shared_junction_elev_soft.get(bucket)
        if e_shared_soft is not None:
            return e_shared_soft, False
        # 3. Wider linear search for off-bucket near-corner matches.
        # Treated as HARD: the value comes from a real rect / runway /
        # terminal corner physically nearby.
        best_e: Optional[float] = None
        best_d2 = NEAR_CORNER_M * NEAR_CORNER_M
        for cx, cy, ce in near_corner_list:
            d2 = (cx - x) * (cx - x) + (cy - y) * (cy - y)
            if d2 < best_d2:
                best_d2 = d2
                best_e = ce
        if best_e is not None:
            if bucket in shared_junction_buckets:
                shared_junction_elev[bucket] = best_e
            return best_e, True
        # 4. Rect-edge interpolation — pulls junction vertices
        # pushed 1m off a sloping edge (or any boundary-trace vertex
        # within NEAR_EDGE_M of a rect/runway/terminal edge) onto
        # the rect's slope at the projected position.  HARD.
        e_edge = _edge_interp_elev(x, y)
        if e_edge is not None:
            if bucket in shared_junction_buckets:
                shared_junction_elev[bucket] = e_edge
            return e_edge, True
        # 5. Free sample from graph / DEM.  SOFT — caller's smoother
        # is free to move this value.  Cache shared-bucket values
        # in the SOFT cache so subsequent junctions see the same
        # initial value.
        e_free: Optional[float] = None
        if dem is not None:
            lat, lon = m_to_ll(x, y)
            e_free = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        if e_free is None:
            return None, False
        if bucket in shared_junction_buckets:
            shared_junction_elev_soft[bucket] = e_free
        return e_free, False

    def _smooth_junction_boundary(
        ring: List[Tuple[float, float]],
        elev: List[float],
        is_anchor: List[bool],
    ) -> List[float]:
        """Bounds-propagation smoothing: each anchor vertex
        constrains every other vertex's elevation to lie within
        ``anchor ± dist × TAXI_MAX_GRADE`` (straight-line distance
        through the junction).  Non-anchor vertices clip to the
        intersection of these bands; anchor vertices stay put.

        Then a Laplacian pass smooths non-anchor elevations toward
        their neighbours' average, re-clipped to the bands so
        anchor compliance is preserved.
        """
        n = len(ring)
        e = list(elev)
        if n < 2:
            return e
        INF = float("inf")
        lo = [-INF] * n
        hi = [INF] * n
        # Anchors lock themselves and contribute bands to others.
        for i in range(n):
            if is_anchor[i]:
                lo[i] = e[i]
                hi[i] = e[i]
        anchor_idx = [i for i in range(n) if is_anchor[i]]
        for i in range(n):
            if is_anchor[i]:
                continue
            xi, yi = ring[i]
            for ai in anchor_idx:
                xa, ya = ring[ai]
                d = math.hypot(xi - xa, yi - ya)
                band = d * TAXI_MAX_GRADE
                lo_i = e[ai] - band
                hi_i = e[ai] + band
                if lo_i > lo[i]:
                    lo[i] = lo_i
                if hi_i < hi[i]:
                    hi[i] = hi_i
        # Initial clip into feasibility intervals.
        for i in range(n):
            if is_anchor[i]:
                continue
            if lo[i] > hi[i]:
                # Conflicting anchors — fall back to midpoint of
                # the conflicting bounds.  Will pull this vertex
                # closer to the closest anchor in practice.
                e[i] = 0.5 * (lo[i] + hi[i])
            else:
                if e[i] < lo[i]:
                    e[i] = lo[i]
                if e[i] > hi[i]:
                    e[i] = hi[i]
        # Laplacian-style smoothing pass between non-anchor
        # neighbours, clipped to bands each iteration so the band
        # constraint stays satisfied.  Convergence in ~20 iters
        # for our junction sizes.
        damping = 0.4
        for _ in range(20):
            new_e = list(e)
            max_change = 0.0
            for i in range(n):
                if is_anchor[i]:
                    continue
                # Mean of all OTHER non-anchor + anchor neighbours,
                # weighted by inverse distance (closer pulls more).
                xi, yi = ring[i]
                num = 0.0
                den = 0.0
                for j in range(n):
                    if j == i:
                        continue
                    xj, yj = ring[j]
                    d = math.hypot(xi - xj, yi - yj)
                    if d < 0.5:
                        continue
                    w = 1.0 / d
                    num += e[j] * w
                    den += w
                if den <= 0:
                    continue
                mean = num / den
                target = e[i] + (mean - e[i]) * damping
                if target < lo[i]:
                    target = lo[i]
                if target > hi[i]:
                    target = hi[i]
                if abs(target - e[i]) > max_change:
                    max_change = abs(target - e[i])
                new_e[i] = target
            e = new_e
            if max_change < 1e-3:
                break
        return e

    # ── Per-junction ring cleanup (independent of elevation) ────
    # Build the cleaned ring + Polygon once per junction; reuse
    # across every cross-junction iteration so we don't redo
    # geometric cleanup at each pass.
    junction_cleaned: List[Tuple["BuiltShape",
                                  List[Tuple[float, float]],
                                  Polygon]] = []
    junction_dropped_shapes: List["BuiltShape"] = []
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        try:
            ring = _splice_holes(shape.polygon)
        except Exception:
            try:
                ring = list(shape.polygon.exterior.coords)
                if ring and ring[0] == ring[-1]:
                    ring = ring[:-1]
            except Exception:
                # Couldn't extract a ring at all — preserve the
                # original shape unchanged downstream.
                junction_dropped_shapes.append(shape)
                continue
        if len(ring) < 3:
            continue
        ring = _drop_spike_vertices(ring)
        if len(ring) < 3:
            continue
        ring = _drop_sliver_corners(ring)
        if len(ring) < 3:
            continue
        try:
            poly = Polygon(ring)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if (poly.is_empty or poly.geom_type != "Polygon"
                    or poly.area < 0.5):
                continue
        except Exception:
            continue
        junction_cleaned.append((shape, ring, poly))

    # ── Single-pass per-junction smoothing + shared-vertex
    # reconciliation.  Per-junction smoothing alone moves SOFT
    # (graph-derived) boundary vertices independently in each
    # junction that touches the same shared bucket, breaking the
    # shared-vertex invariant.  After all junctions have smoothed,
    # we walk every shared SOFT bucket and overwrite each junction's
    # value with the cross-junction average; this restores
    # consistency at shared vertices.
    #
    # Per the elevation field plan (2026-04-26), Step 3's iterative
    # variant of this — average then lock then re-smooth — was
    # tried and made within-shape grade worse, because forcing the
    # average at a shared vertex pulls it away from each junction's
    # locally-optimal smooth solution and constrains the next
    # smoothing round.  A single average pass keeps the per-
    # junction smoothing's local optimum and only sacrifices a
    # small grade residual at shared vertices.
    iter_results: List[Tuple["BuiltShape",
                             List[Tuple[float, float]],
                             List[float]]] = []
    for (shape, ring, poly) in junction_cleaned:
        ev_pairs = [_vertex_elev_anchored(x, y) for (x, y) in ring]
        vert_elev_raw: List[Optional[float]] = [p[0] for p in ev_pairs]
        is_anchor_list: List[bool] = [p[1] for p in ev_pairs]
        known = [e for e in vert_elev_raw if e is not None]
        fallback = sum(known) / len(known) if known else 0.0
        vert_elev = [e if e is not None else fallback
                     for e in vert_elev_raw]
        # Mode B (gated): replace free-vertex elevations with
        # bilinear samples from a 2D smoothed field.  Hard anchors
        # keep their corner-derived value.
        if USE_PER_POLYGON_ELEVATION_FIELD:
            try:
                _anchors_b: List[
                    Tuple[float, float, float]] = [
                        (rx, ry, e)
                        for (rx, ry), e, isa in zip(
                            ring, vert_elev, is_anchor_list)
                        if isa]
                sampler_b = _smooth_polygon_grid(
                    poly, _anchors_b, dem,
                    tile_lat, tile_lon, layout.anchor)
                if sampler_b is not None:
                    for k, (rx, ry) in enumerate(ring):
                        if is_anchor_list[k]:
                            continue
                        e_b = sampler_b(rx, ry)
                        if e_b is not None:
                            vert_elev[k] = float(e_b)
            except Exception:
                pass
        # Junction-local 1D boundary smoothing — pull free vertices
        # into anchor-band compliance with the polygon's hard
        # anchors.
        vert_elev = _smooth_junction_boundary(
            ring, vert_elev, is_anchor_list)
        iter_results.append((shape, ring, vert_elev))

    # Cross-junction shared-vertex reconciliation.  Per-junction
    # smoothing moves SOFT (graph-derived) shared boundary vertices
    # independently — at most shared buckets the moves agree across
    # all junctions to within a small tolerance (junctions have
    # similar local boundaries near a shared vertex), but at some
    # buckets the per-junction smoothing diverges and we'd ship a
    # visible cross-shape cliff.  For the divergent buckets, revert
    # to the cached SOFT graph value (which both junctions would
    # have read consistently anyway — same as baseline behaviour
    # for shared buckets).  Where smoothing converged consistently
    # across junctions, keep the smoothed value (this is the
    # within-shape-grade win that the SOFT classification
    # delivers).
    # Cross-junction reconciliation at shared SOFT buckets.  Each
    # junction's per-junction smoother moved its copy of the
    # shared vertex toward its own neighbour mean, so two
    # junctions touching the same bucket can disagree on the
    # final elevation.  For buckets where the disagreement
    # exceeds SHARED_AGREE_TOL_M we average across all junctions
    # and overwrite — restoring the shared-vertex invariant at
    # the cost of a small within-shape grade residual local to
    # that vertex.  For buckets where every junction's smoother
    # already converged to within tolerance, we leave the
    # per-junction values alone — preserving the SOFT
    # classification's within-shape-grade win.  The tolerance is
    # generous (0.10 m, 1 % grade across a 10 m shared edge) so
    # most shared SOFT buckets pass without averaging.  Constant
    # is module-level so both legacy and Mode B paths can use it.
    shared_smooth_samples: Dict[
        Tuple[int, int], List[float]] = {}
    for (shape, ring, vert_elev) in iter_results:
        seen: set = set()
        for (vx, vy), e in zip(ring, vert_elev):
            bucket = _corner_elevation_bucket(vx, vy)
            if bucket in seen:
                continue
            seen.add(bucket)
            if bucket in corner_elev:
                continue
            if bucket not in shared_junction_buckets:
                continue
            if bucket not in shared_junction_elev_soft:
                continue
            shared_smooth_samples.setdefault(
                bucket, []).append(e)
    shared_lock: Dict[Tuple[int, int], float] = {}
    for bucket, vals in shared_smooth_samples.items():
        if len(vals) >= 2 and (
                max(vals) - min(vals) > SHARED_AGREE_TOL_M):
            shared_lock[bucket] = sum(vals) / len(vals)
    if shared_lock:
        for (shape, ring, vert_elev) in iter_results:
            for k, (vx, vy) in enumerate(ring):
                bucket = _corner_elevation_bucket(vx, vy)
                if bucket in shared_lock:
                    vert_elev[k] = shared_lock[bucket]

    new_shapes: List[BuiltShape] = []
    triangle_count = 0
    grade_violations = 0
    # Carry every non-junction shape through unchanged, plus any
    # junction shape that failed pre-loop cleanup.
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            new_shapes.append(shape)
    for shape in junction_dropped_shapes:
        new_shapes.append(shape)
    # Emit each cleaned junction with its post-iteration vert_elev
    # passed through the FLAT/PLANAR/COMPOUND classifier.
    for (shape, ring, vert_elev) in iter_results:
        # ── Surface-complexity classification ───────────────────
        # User 2026-04-25: only triangulate where the surface has
        # a compound slope.  Flat or planar regions can stay as a
        # single polygon — fewer shapes, cleaner OSM output, less
        # work for X-Plane's mesh builder.
        #
        #   * FLAT       — vertex elevations vary by < 0.2 % of the
        #                  polygon's bbox extent (with a floor of
        #                  ``FLAT_ABS_FLOOR_M``).  Emit one polygon
        #                  with a single ``altitude`` tag.
        #   * PLANAR     — every vertex sits within
        #                  ``PLANAR_RESIDUAL_M`` of the best-fit
        #                  plane through them.  Emit one polygon
        #                  with ``node_altitudes`` (X-Plane
        #                  triangulates internally; since the
        #                  surface is planar the result is identical
        #                  to our pre-triangulated mesh).
        #   * COMPOUND   — triangulate as before.
        FLAT_GRADE_THRESHOLD = 0.001      # 0.1 % (so 0.05 m at 50 m
                                           # bbox) — tightened to keep
                                           # the FLAT mean within
                                           # check_grade's
                                           # SHARED_NID_TOLERANCE_M
                                           # (0.15 m) of every anchor.
                                           # Plateau snapping makes
                                           # most flat regions
                                           # already fit this; what
                                           # doesn't falls to PLANAR
                                           # which emits per-vertex
                                           # elevations matching each
                                           # anchor exactly.
        FLAT_ABS_FLOOR_M = 0.05           # 5 cm — half of one stored
                                           # decimal, so the rounded
                                           # mean can't differ from
                                           # any vertex by more than
                                           # rounding noise.
        PLANAR_RESIDUAL_M = 0.30          # 30 cm — vertices
                                           # within this of best-fit
                                           # plane render virtually
                                           # identical to a
                                           # triangulated mesh
        # Polygon bbox extent.
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        extent = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
        flat_thresh = max(FLAT_ABS_FLOOR_M,
                          FLAT_GRADE_THRESHOLD * extent)
        elev_range = max(vert_elev) - min(vert_elev)
        if elev_range < flat_thresh:
            # FLAT — single polygon, single altitude tag.
            try:
                flat_poly = Polygon(ring)
                if not flat_poly.is_valid:
                    flat_poly = flat_poly.buffer(0)
                if (flat_poly.geom_type == "Polygon"
                        and not flat_poly.is_empty
                        and flat_poly.area >= 0.5):
                    new_shape = BuiltShape(
                        polygon=flat_poly,
                        role=ROLE_JUNCTION,
                        ref=shape.ref)
                    new_shape.altitude = round(
                        sum(vert_elev) / len(vert_elev), 1)
                    new_shapes.append(new_shape)
                    continue
            except Exception:
                pass  # fall through to triangulation
        # Best-fit plane: solve ax + by + c = z via 3×3 normal eqns.
        # Classify as PLANAR only when the fit's residuals are tight
        # AND the plane's slope itself is grade-compliant.  A polygon
        # whose vertices lie on a 5 % plane has tiny residuals but
        # X-Plane would render a 5 % cross-grade across the polygon
        # — well above TAXI_MAX_GRADE.  Demote those to COMPOUND so
        # _smooth_junction_boundary's anchor-band clamping pulls the
        # free vertices into compliance instead.
        plane = _planar_fit(ring, vert_elev)
        residuals = plane[3] if plane is not None else None
        if plane is not None:
            slope_mag = math.hypot(plane[0], plane[1])
        else:
            slope_mag = float("inf")
        if (residuals is not None
                and max(residuals) < PLANAR_RESIDUAL_M
                and slope_mag <= TAXI_MAX_GRADE):
            # PLANAR — single polygon, per-vertex altitudes.
            try:
                planar_poly = Polygon(ring)
                if not planar_poly.is_valid:
                    planar_poly = planar_poly.buffer(0)
                if (planar_poly.geom_type == "Polygon"
                        and not planar_poly.is_empty
                        and planar_poly.area >= 0.5):
                    # node_altitudes spans the closed ring; match
                    # the order Polygon(ring).exterior.coords
                    # produced (which is ring + closing first).
                    closed_ring = list(planar_poly.exterior.coords)
                    closed_elev = [
                        round(float(_match_elev(rx, ry, ring,
                                                vert_elev)), 1)
                        for (rx, ry) in closed_ring]
                    new_shape = BuiltShape(
                        polygon=planar_poly,
                        role=ROLE_JUNCTION,
                        ref=shape.ref)
                    new_shape.node_altitudes = closed_elev
                    new_shapes.append(new_shape)
                    continue
            except Exception:
                pass  # fall through to triangulation

        # ── Compound-slope: emit as a single polygon with per-
        # vertex node_altitudes; defer triangulation to Triangle4XP.
        #
        # Earlier versions ear-clipped or Delaunay-triangulated
        # the junction here, emitting one patch.osm way per output
        # triangle.  Each triangle was geometrically determined
        # (3 vertices), giving Triangle4XP nothing to refine —
        # any sliver / steep-plane triangle from our triangulator
        # made it into the rendered mesh untouched.
        #
        # Submitting the WHOLE polygon as a single constraint
        # boundary lets Triangle4XP's quality-refinement pass
        # ([O4_Mesh_Utils.py:670] flag ``-pq``) insert interior
        # Steiners with min-angle ≥ 20°.  Steiner elevations are
        # bilinear-interpolated from the polygon's boundary
        # ``node_altitudes`` at mesh time — exactly how the legacy
        # Ortho4XP pavement smoothing avoided cliffs without any
        # explicit grade enforcement.
        try:
            poly_for_emit = Polygon(ring)
            if not poly_for_emit.is_valid:
                poly_for_emit = poly_for_emit.buffer(0)
            if (poly_for_emit.is_empty
                    or poly_for_emit.geom_type != "Polygon"
                    or poly_for_emit.area < 0.5):
                continue
        except Exception:
            continue
        # node_altitudes spans the closed ring (one value per
        # vertex including the closing-repeat).  Build it from the
        # smoothed vert_elev in shapely's emitted ring order.
        closed = list(poly_for_emit.exterior.coords)
        elev_for_ring: List[float] = []
        for (rx, ry) in closed:
            best_e = vert_elev[0]
            best_d2 = float("inf")
            for (tx, ty), te in zip(ring, vert_elev):
                d2 = (rx - tx) * (rx - tx) + (ry - ty) * (ry - ty)
                if d2 < best_d2:
                    best_d2 = d2
                    best_e = te
            elev_for_ring.append(round(float(best_e), 1))
        new_shape = BuiltShape(
            polygon=poly_for_emit,
            role=ROLE_JUNCTION,
            ref=shape.ref)
        if (max(elev_for_ring) - min(elev_for_ring)) < 0.05:
            # Pre-classifier missed; emit flat.
            new_shape.altitude = round(
                sum(elev_for_ring[:-1]) / len(elev_for_ring[:-1]),
                1)
        else:
            new_shape.node_altitudes = elev_for_ring
        new_shapes.append(new_shape)
        triangle_count += 1

    layout.shapes = new_shapes
    if grade_violations:
        # Surfaced via stderr so the user sees it during the test
        # tool run; not a hard failure.
        try:
            import sys
            UI.vprint(1,
                f"  [pav-builder] WARN: {grade_violations} junction "
                f"triangle(s) exceed {TAXI_MAX_GRADE * 100:.1f}% grade.")
        except Exception:
            pass
    return triangle_count


# Layer 2: free-vertex clamping with neighbour-boundary lookup ─────
#
# After _triangulate_junctions sets every junction's per-vertex
# elevation, walk each junction's boundary and tighten any FREE
# vertex (one whose elevation is NOT pinned to a rect/runway/
# terminal corner or to another junction's shared vertex) so that
# it satisfies grade compliance with every nearby shape boundary
# point — not just the same-polygon anchored vertices considered
# during the first smoothing pass.
#
# Why a separate pass:
#   The first-pass smoothing (_smooth_junction_boundary) only
#   considers anchors WITHIN the same polygon.  Adjacent junction
#   polygons that are 1.5–5 m apart but don't share node IDs
#   (HECA's "adjacent-but-not-shared" pattern) can have their
#   nearby boundary vertices drift to incompatible elevations,
#   producing the visible sunken-area / elevated-plateau cliffs
#   the user reported at SPJC and the 159% within-shape grade
#   violations at HECA.
#
# What this pass does NOT change:
#   * Rect / runway / terminal altitudes — those are derived from
#     the elevation graph and the 1.5 % grade rule along the taxi
#     network; this pass operates only on junction polygons.
#   * Anchored junction vertices (shared with another shape via
#     the OSM nid).  Moving them would break the shared-vertex
#     invariant.
#
# Where it operates:
#   For each junction polygon's free vertex V at coord (x, y) with
#   current elevation e_v, scan every other shape's boundary
#   edges within NEIGHBOUR_CLAMP_RADIUS_M.  Each nearby boundary
#   point E at distance d contributes a feasibility band
#   ``[E_elev − d × TAXI_MAX_GRADE, E_elev + d × TAXI_MAX_GRADE]``.
#   Intersect these bands and clip e_v.  When the intersection is
#   empty (anchors disagree), fall back to the midpoint — Layer 1
#   will eventually reconcile the conflicting anchors.
# (NEIGHBOUR_CLAMP_RADIUS_M is now defined in elevation.py and
# imported at the top of this module.)


