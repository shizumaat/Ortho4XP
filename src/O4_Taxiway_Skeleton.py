"""Voronoi-based medial-axis skeleton extraction for taxiway polygons.

Given an apt.dat taxiway mega-polygon like SPJC's "Taxiway V / U / Q
/ R / L / M" (857 929 m², six concatenated taxiways joined through
a central hub), we want to split it into its individual taxiway
centerlines so each can be emitted as a chain of sloping
rectangles via :func:`O4_Taxiway_Rects.build_rects_along_centerline`.

The morphological decomposition in :mod:`O4_Taxiway_Decompose` only
works when the junction hubs are dramatically wider than the strip
arms, which is not the case for most real airports — the hub is
usually only ~1.5× wider than a taxiway, so morphological opening
can't isolate it.  Medial-axis skeletonisation is the correct tool.

Pipeline:

1. **Densify the polygon boundary** every ``densify_step`` metres
   so the Voronoi diagram sees enough constraint points to
   approximate the medial axis.  A 5 m step gives a 5076-point
   input for the SPJC mega-polygon.

2. **Compute the Voronoi diagram** of the dense boundary points
   via :func:`shapely.ops.voronoi_diagram`.  The diagram's edges
   contain *both* the interior ridge segments (≈ the medial axis
   skeleton) *and* the exterior segments escaping to infinity.

3. **Filter to interior edges** — keep only Voronoi edges whose
   BOTH endpoints lie strictly inside the polygon (a tiny inward
   buffer absorbs boundary-touching precision slop).  These
   segments form the skeleton.

4. **Linemerge** the interior edges into continuous paths.
   :func:`shapely.ops.linemerge` stitches touching segments into
   longer ``LineString`` objects; a Y-junction still shows as
   three separate merged paths meeting at a node, which is
   exactly what we want — each arm of the Y is its own taxiway
   branch.

5. **Drop short paths** below ``min_path_length``.  The Voronoi
   skeleton is noisy at boundary concavities and produces many
   tiny spurs (2-20 m); keeping only long paths prunes these to
   the real taxiway centerlines.

6. **Simplify** each kept path with a Ramer-Douglas-Peucker
   tolerance.  The user has asked for ~1 m polygon detail; on the
   centerline we can be looser (default 3 m) because a 3 m
   centerline deviation on a 30 m-wide taxiway is only 10 % of
   the width and well below any grade-rule threshold.

The output is a list of ``LineString`` centerlines — one per
taxiway branch.  Each centerline is in meter-space and directly
consumable by the rect-chain builder.
"""
from __future__ import annotations

from typing import List, Optional

from shapely.geometry import (
    LineString, MultiLineString, MultiPoint, Point, Polygon,
)
from shapely.ops import linemerge, voronoi_diagram


# Tunables
DEFAULT_DENSIFY_STEP_M = 5.0
DEFAULT_MIN_PATH_LENGTH_M = 200.0   # drop short spurs aggressively
DEFAULT_SIMPLIFY_TOL_M = 20.0  # Skeleton centerline simplification
                               # tolerance in meters.  A larger
                               # value produces longer straighter
                               # segments (one rect each) at the
                               # cost of rotation deviation; the
                               # user's polygon-overflow allowance
                               # makes this safe.
_INTERIOR_SHRINK_M = 0.05


def _iter_voronoi_edges(vd) -> List[LineString]:
    """Walk the GeometryCollection returned by
    ``voronoi_diagram(..., edges=True)`` and yield the
    individual LineString edges.
    """
    edges: List[LineString] = []
    if vd is None or vd.is_empty:
        return edges
    containers = (list(vd.geoms) if hasattr(vd, "geoms") else [vd])
    for g in containers:
        if g is None or g.is_empty:
            continue
        if g.geom_type == "LineString":
            edges.append(g)
        elif g.geom_type == "MultiLineString":
            edges.extend(list(g.geoms))
        elif hasattr(g, "geoms"):
            edges.extend(_iter_voronoi_edges(g))
    return edges


def _densify_boundary(polygon: Polygon, step: float) -> List[tuple]:
    """Sample the polygon's exterior ring and every interior ring
    at uniform spacing.  Returns a flat list of ``(x, y)`` tuples.
    """
    pts: List[tuple] = []

    def _walk(ring):
        length = ring.length
        if length <= 0:
            return
        n = max(4, int(length / step))
        for i in range(n):
            pt = ring.interpolate(i * length / n)
            pts.append((pt.x, pt.y))

    _walk(polygon.exterior)
    for hole in polygon.interiors:
        _walk(hole)
    return pts


def extract_centerlines(
    polygon: Polygon,
    densify_step: float = DEFAULT_DENSIFY_STEP_M,
    min_path_length: float = DEFAULT_MIN_PATH_LENGTH_M,
    simplify_tol: float = DEFAULT_SIMPLIFY_TOL_M,
) -> List[LineString]:
    """Return the taxiway centerlines inside ``polygon``.

    Args:
        polygon: meter-space polygon to skeletonise.  May contain
            interior rings.
        densify_step: boundary sampling step in metres.  Smaller =
            smoother skeleton but slower Voronoi.  5 m is a good
            tradeoff for SPJC's 860 k m² mega-polygon.
        min_path_length: drop merged paths shorter than this.  The
            default 50 m discards Voronoi spurs at boundary
            concavities without sacrificing short real connectors
            (which are usually ≥ 80 m at real airports).
        simplify_tol: Ramer-Douglas-Peucker tolerance applied to
            each kept path.

    Returns:
        A list of :class:`shapely.geometry.LineString` centerlines,
        one per taxiway branch, sorted longest-first.  Empty list
        for a polygon too small or too narrow to skeletonise.
    """
    if polygon is None or polygon.is_empty:
        return []
    if polygon.area < 50.0:
        return []

    dense_pts = _densify_boundary(polygon, densify_step)
    if len(dense_pts) < 4:
        return []

    try:
        vd = voronoi_diagram(MultiPoint(dense_pts), edges=True)
    except Exception:
        return []

    edges = _iter_voronoi_edges(vd)
    if not edges:
        return []

    # Interior-edge filter.  An "interior" Voronoi edge is one
    # whose both endpoints lie strictly inside the polygon — that's
    # the standard medial-axis approximation.  A tiny inward
    # shrink absorbs edges that merely touch the boundary.
    try:
        interior_test = polygon.buffer(-_INTERIOR_SHRINK_M)
    except Exception:
        interior_test = polygon
    if interior_test.is_empty or not hasattr(interior_test, "contains"):
        interior_test = polygon

    interior_edges: List[LineString] = []
    for e in edges:
        cc = list(e.coords)
        if len(cc) < 2:
            continue
        p0 = Point(cc[0])
        p1 = Point(cc[-1])
        try:
            if (interior_test.contains(p0)
                    and interior_test.contains(p1)):
                interior_edges.append(e)
        except Exception:
            continue

    if not interior_edges:
        return []

    try:
        merged = linemerge(MultiLineString(interior_edges))
    except Exception:
        merged = None
    if merged is None or merged.is_empty:
        return []

    if merged.geom_type == "LineString":
        paths = [merged]
    elif hasattr(merged, "geoms"):
        paths = [g for g in merged.geoms
                 if g.geom_type == "LineString" and not g.is_empty]
    else:
        return []

    out: List[LineString] = []
    for p in paths:
        if p.length < min_path_length:
            continue
        try:
            simp = p.simplify(simplify_tol, preserve_topology=False)
        except Exception:
            simp = p
        if simp.is_empty or simp.length < min_path_length:
            continue
        out.append(simp)

    out.sort(key=lambda ln: -ln.length)
    return out


def local_half_width(
    polygon: Polygon,
    center: Point,
    tangent: "tuple[float, float]",
    max_reach: float = 60.0,
) -> float:
    """Return the distance from ``center`` to the nearest polygon
    boundary perpendicular to ``tangent``, capped at ``max_reach``.

    Used by the rect-chain builder to size each slab rect to the
    local polygon width.  Returns 0.0 if the centerpoint is not
    inside the polygon.
    """
    if polygon is None or polygon.is_empty:
        return 0.0
    try:
        if not polygon.contains(center):
            return 0.0
    except Exception:
        return 0.0

    tx, ty = tangent
    # Perpendicular unit vector
    import math
    mag = math.hypot(tx, ty)
    if mag < 1e-9:
        return 0.0
    px, py = -ty / mag, tx / mag

    # Shoot a ray in each direction out to max_reach, find the
    # first intersection with the polygon boundary.
    left_end = (center.x + px * max_reach, center.y + py * max_reach)
    right_end = (center.x - px * max_reach, center.y - py * max_reach)
    try:
        left_line = LineString([(center.x, center.y), left_end])
        right_line = LineString([(center.x, center.y), right_end])
        left_seg = left_line.intersection(polygon)
        right_seg = right_line.intersection(polygon)
    except Exception:
        return 0.0

    def _seg_length(seg):
        if seg is None or seg.is_empty:
            return 0.0
        if hasattr(seg, "length"):
            return float(seg.length)
        return 0.0

    left_d = _seg_length(left_seg)
    right_d = _seg_length(right_seg)
    return min(left_d, right_d)
