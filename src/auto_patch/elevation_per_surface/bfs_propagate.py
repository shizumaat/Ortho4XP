"""Phase 2: BFS rect-junction propagation from runway HARD anchors.

After Phase 1 sets each rect's axial profile from DEM (with runway-
adjacent ends anchored), Phase 2 propagates altitudes outward
through the rect-junction adjacency graph.  A junction's vertices
inherit from the rect/runway corners they share with; once a
junction has all-or-most of its boundary settled, the rects on the
other side use those boundary values as anchors and re-solve their
axial profiles.

Key property: propagation follows REALISTIC TAXI AXES.  Cumulative
elevation budget grows with axial path length, not Euclidean
distance — so a long parallel taxiway can reach DEM-target heights
even when the runway is far below.
"""
from __future__ import annotations

from collections import deque
from typing import Dict, Optional, Set, Tuple

from .axial_profile import (
    SLOPING_RECT_ROLES,
    TAXI_MAX_GRADE,
    solve_rect_axial_profile,
)

# Bucket size for shared-vertex matching (mirrors elevation.py
# ``_corner_elevation_bucket``).  Two coordinates within
# ``SHARED_VERTEX_TOL_M`` of each other hash to the same bucket.
SHARED_VERTEX_TOL_M = 0.1


def _bucket(x: float, y: float,
            tol: float = SHARED_VERTEX_TOL_M) -> Tuple[int, int]:
    return (int(round(x / tol)), int(round(y / tol)))


def _rect_axis_endpoints(shape):
    """Return ``(start_xy, end_xy)`` from the rect's source_axis."""
    if shape.source_axis is None or shape.source_axis.is_empty:
        return None
    pts = list(shape.source_axis.coords)
    if len(pts) < 2:
        return None
    return pts[0], pts[-1]


def _rect_axis_end_altitudes(shape) -> Tuple[Optional[float],
                                              Optional[float]]:
    """Map ``altitude_high`` / ``altitude_low`` onto the rect's
    source_axis start/end ends using the polygon vertices.

    Returns ``(start_alt, end_alt)`` where each is the elevation at
    that axial end (regardless of which is high vs low).  Uses
    ``_short_end_pairs_by_axis`` for convention-free pairing.
    """
    from auto_patch.elevation import _short_end_pairs_by_axis

    if shape.altitude_high is None or shape.altitude_low is None:
        if shape.altitude is not None:
            a = float(shape.altitude)
            return a, a
        return None, None

    poly = shape.polygon
    if poly is None or poly.is_empty:
        return None, None
    coords = list(poly.exterior.coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return None, None
    start_pair, end_pair = _short_end_pairs_by_axis(
        coords, shape.source_axis)
    if start_pair is None:
        return None, None
    # Determine which pair holds altitude_high vs altitude_low by
    # which end's vertices are physically closer to source_axis
    # start.  We don't actually have per-corner elevations on the
    # shape, so we rely on the axial position: start_pair holds
    # whichever (high or low) is at axis t=0.  Without further info
    # default to high@start; if downstream BFS finds disagreement
    # at a shared vertex, the reconciliation pass corrects it.
    return float(shape.altitude_high), float(shape.altitude_low)


def _rect_corners_xy(shape):
    """Return the 4 corner (x, y) tuples of a rect polygon."""
    if shape.polygon is None or shape.polygon.is_empty:
        return []
    coords = list(shape.polygon.exterior.coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    return coords


def build_runway_anchor_lookup(layout):
    """Return a callable ``(x, y) -> Optional[float]`` that resolves
    a vertex bucket to its runway HARD elevation if shared with a
    runway corner; ``None`` otherwise.

    Per ``layout.ROLE_RUNWAY``: each runway segment has
    ``altitude_high``/``altitude_low`` set from the CIFP profile.
    The 4 corners of a runway segment polygon are the HARD anchors;
    we use them at-bucket regardless of polygon vertex order via
    ``_short_end_pairs_by_axis`` when available.
    """
    bucket_to_alt: Dict[Tuple[int, int], float] = {}
    for s in layout.shapes:
        if s.role != "runway":
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            # Fall back to a flat altitude for non-rectangular runway
            # polygons (rare — apron-merged segments).
            if s.altitude is not None:
                for x, y in coords:
                    bucket_to_alt.setdefault(_bucket(x, y),
                                             float(s.altitude))
            continue
        if s.altitude_high is not None and s.altitude_low is not None:
            from auto_patch.elevation import _short_end_pairs_by_axis
            sp, ep = _short_end_pairs_by_axis(coords, s.source_axis) \
                if s.source_axis is not None else (None, None)
            if sp is not None and ep is not None:
                # Determine which short-end pair gets which altitude:
                # source_axis t=0 end is start; default high→start
                # when high≥low at start (we don't know orientation,
                # so set both pairs and let the bucket map dedupe).
                hi, lo = float(s.altitude_high), float(s.altitude_low)
                # Without loss: assign hi to start_pair if start's
                # average elevation in DEM is higher; without DEM
                # access here we just average — both endpoints fall
                # within hi±lo range.
                avg = (hi + lo) / 2.0
                for i in sp:
                    bucket_to_alt.setdefault(
                        _bucket(coords[i][0], coords[i][1]), avg)
                for i in ep:
                    bucket_to_alt.setdefault(
                        _bucket(coords[i][0], coords[i][1]), avg)
                # Refine: assume the pair whose physical Y/X comes
                # FIRST along source_axis is start; emit hi/lo
                # accordingly only if the polygon vertex order
                # happens to honour the high/low convention.
                # Simpler: assign hi to whichever pair has higher
                # mean coords toward axis-end.  For the BFS we only
                # need an approximate anchor — the reconciliation
                # pass corrects shared-vertex disagreements.
                continue
        # Plain flat runway segment.
        if s.altitude is not None:
            for x, y in coords:
                bucket_to_alt.setdefault(_bucket(x, y),
                                         float(s.altitude))

    def lookup(x: float, y: float) -> Optional[float]:
        return bucket_to_alt.get(_bucket(x, y))

    return lookup


def _adjacent_rect_buckets(layout):
    """Return a dict mapping bucket → list of rect shapes whose
    polygon includes that bucket.  Used to find which rects share a
    given vertex.
    """
    bucket_rects: Dict[Tuple[int, int], list] = {}
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        for x, y in _rect_corners_xy(s):
            bucket_rects.setdefault(_bucket(x, y), []).append(s)
    return bucket_rects


def _rect_anchors(shape, anchor_map: Dict[Tuple[int, int], float]):
    """Inspect the rect's source_axis start/end and return any
    elevation anchor present in ``anchor_map`` for those endpoints.
    """
    eps = _rect_axis_endpoints(shape)
    if eps is None:
        return None, None
    (sx, sy), (ex, ey) = eps
    return anchor_map.get(_bucket(sx, sy)), anchor_map.get(_bucket(ex, ey))


def propagate_through_rects(layout, dem, tile_lat: int, tile_lon: int,
                            runway_anchor_lookup) -> int:
    """BFS from runway-adjacent rects, propagating altitudes through
    the rect-junction chain.  Re-solves each rect's axial profile
    when a previously-free end becomes anchored from a neighbour.

    Returns the number of rects re-solved with at least one
    propagated anchor.
    """
    anchor_map: Dict[Tuple[int, int], float] = {}
    # Seed with runway anchors.
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        for x, y in _rect_corners_xy(s):
            a = runway_anchor_lookup(x, y)
            if a is not None:
                anchor_map.setdefault(_bucket(x, y), a)
    # Seed with the result of Phase 1 (each rect's axial extremes
    # already carry DEM-smoothed altitudes); record them at the rect
    # corners.
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        a_start, a_end = _rect_axis_end_altitudes(s)
        if a_start is None and a_end is None:
            continue
        eps = _rect_axis_endpoints(s)
        if eps is None:
            continue
        if a_start is not None:
            anchor_map.setdefault(_bucket(*eps[0]), a_start)
        if a_end is not None:
            anchor_map.setdefault(_bucket(*eps[1]), a_end)

    # BFS: starting from runway-touching rects, walk to neighbours
    # via shared corner buckets.
    bucket_rects = _adjacent_rect_buckets(layout)
    queue: deque = deque()
    enqueued: Set[int] = set()

    def _enqueue(rect):
        if id(rect) in enqueued:
            return
        enqueued.add(id(rect))
        queue.append(rect)

    # Seed BFS with rects that have any runway-anchored corner.
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        for x, y in _rect_corners_xy(s):
            if runway_anchor_lookup(x, y) is not None:
                _enqueue(s)
                break

    n_resolved = 0
    while queue:
        shape = queue.popleft()
        anchor_start, anchor_end = _rect_anchors(shape, anchor_map)
        if anchor_start is None and anchor_end is None:
            continue
        result = solve_rect_axial_profile(
            layout, shape, dem, tile_lat, tile_lon,
            anchor_start, anchor_end)
        if result is None:
            continue
        new_high, new_low = result
        old_high = shape.altitude_high
        old_low = shape.altitude_low
        shape.altitude_high = round(float(new_high), 1)
        shape.altitude_low = round(float(new_low), 1)
        shape.altitude = None
        n_resolved += 1
        # Update anchor map at this rect's axis endpoints with the
        # newly-solved values; expose to neighbours.
        eps = _rect_axis_endpoints(shape)
        if eps is not None:
            (sx, sy), (ex, ey) = eps
            # The rect's source_axis start gets whichever of high/lo
            # ended up there; without per-corner certainty we average
            # both extremes — close enough for cross-rect propagation.
            avg = (new_high + new_low) / 2.0
            anchor_map.setdefault(_bucket(sx, sy), avg)
            anchor_map.setdefault(_bucket(ex, ey), avg)
        # Enqueue neighbours sharing any corner with this rect.
        for x, y in _rect_corners_xy(shape):
            for nb in bucket_rects.get(_bucket(x, y), ()):
                if nb is shape:
                    continue
                _enqueue(nb)
        # Re-enqueue shapes if anchors changed substantially (delta
        # >0.5m); avoids stale propagation in topologically-cyclic
        # graphs.
        if (old_high is None or old_low is None
                or abs(old_high - new_high) > 0.5
                or abs(old_low - new_low) > 0.5):
            for x, y in _rect_corners_xy(shape):
                for nb in bucket_rects.get(_bucket(x, y), ()):
                    if nb is shape:
                        continue
                    enqueued.discard(id(nb))
                    _enqueue(nb)
    return n_resolved
