"""Junction-polygon construction + residue decomposition.

Two responsibilities live in this module:

1. **Construction** — build junction polygons from rect-endpoint
   corner sets discovered along OSM centerlines:
   ``_find_junction_points``, ``_build_junctions_from_rect_endpoints``,
   ``_rect_end_corners``, ``_build_junction_constructive``,
   ``_build_junction_polys_from_corners``.

2. **Decomposition + densification** — turn a residue Polygon
   (apt.dat pavement minus rect / runway / terminal coverage) into
   a set of junction-polygon pieces with hole splicing, sliver
   removal, long-edge densification, and colinear-vertex pruning:
   ``_decompose_polygon_with_holes``, ``_polygon_min_thickness``,
   ``_merge_thin_decomposed_pieces``, ``_splice_holes``,
   ``_polygon_area``, ``_splice_one_hole``,
   ``_densify_long_boundary_edges``, ``_drop_sliver_corners``,
   ``_drop_colinear_boundary_vertices``.

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _build_junction_constructive
    _build_junction_polys_from_corners
    _build_junctions_from_rect_endpoints
    _decompose_polygon_with_holes
    _densify_long_boundary_edges
    _drop_colinear_boundary_vertices
    _drop_sliver_corners
    _find_junction_points
    _merge_thin_decomposed_pieces
    _polygon_area
    _polygon_min_thickness
    _rect_end_corners
    _splice_holes
    _splice_one_hole
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Set, Tuple

from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.ops import linemerge, nearest_points, unary_union

from ..config import (
    JUNCTION_CLUSTER_DIST_M,
    MAX_BOUNDARY_EDGE_M,
    SLIVER_ANGLE_THRESHOLD_DEG,
)
from ..layout import (
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
    SHARED_VERTEX_TOL_M,
)


# Max distance from a junction-ring midpoint to a neighbour-edge
# below which the midpoint sits on a shared boundary and inherits
# the neighbour's edge-interpolated elevation.
SHARED_NEIGHBOUR_EDGE_TOL_M = 5.0

# Max perpendicular distance to neighbours below which a colinear
# boundary vertex may be dropped without breaking shared-vertex
# alignment.
COLINEAR_DROP_M = 3.0


__all__ = [
    "COLINEAR_DROP_M",
    "SHARED_NEIGHBOUR_EDGE_TOL_M",
    "_build_junction_constructive",
    "_build_junction_polys_from_corners",
    "_build_junctions_from_rect_endpoints",
    "_decompose_polygon_with_holes",
    "_densify_long_boundary_edges",
    "_drop_colinear_boundary_vertices",
    "_drop_sliver_corners",
    "_find_junction_points",
    "_merge_thin_decomposed_pieces",
    "_polygon_area",
    "_polygon_min_thickness",
    "_rect_end_corners",
    "_splice_holes",
    "_splice_one_hole",
]


def _decompose_polygon_with_holes(polygon: Polygon,
                                  min_area_m2: float = 50.0,
                                  max_depth: int = 8,
                                  runway_axis_deg: Optional[float] = None,
                                  corner_snap_pts: Optional[
                                      List[Tuple[float, float]]] = None,
                                  corner_snap_tol_m: float = 5.0,
                                  ) -> List[Polygon]:
    """Return a list of simple (no-hole) polygons that tile the same
    area as ``polygon``.

    Strategy: cut the largest interior hole through its centroid in
    the runway-parallel-or-perpendicular direction that yields the
    most-balanced split.  Recurse on each piece.  Pieces with no
    holes are returned as-is.

    Per user 2026-05-04: this function should ideally be a no-op —
    if every apt.dat hole is bordered by a rect/runway/terminal
    upstream, every residue piece is simply connected and the cut
    machinery is unnecessary.  The cut path here is a defensive
    fallback for residue pieces that still have interior holes.
    """
    from shapely.ops import split as _shp_split
    from shapely.geometry import LineString as _LS

    if (polygon.is_empty or polygon.geom_type != "Polygon"):
        return []
    if not polygon.interiors:
        return [polygon]
    interiors = list(polygon.interiors)
    big_interiors = [h for h in interiors
                     if Polygon(h).area >= min_area_m2]
    if not big_interiors:
        return [Polygon(polygon.exterior.coords)]
    if max_depth <= 0:
        clean = Polygon(polygon.exterior.coords, big_interiors)
        spliced_coords = _splice_holes(clean)
        try:
            return [Polygon(spliced_coords).buffer(0)]
        except Exception:
            return [Polygon(polygon.exterior.coords)]
    big_interiors.sort(key=lambda h: -Polygon(h).area)
    hole = big_interiors[0]
    cx = float(hole.centroid.x)
    cy = float(hole.centroid.y)
    minx, miny, maxx, maxy = polygon.bounds
    span = max(maxx - minx, maxy - miny) + 2.0

    angle_rad = 0.0
    if runway_axis_deg is not None:
        bearing_rad = math.radians(runway_axis_deg)
        runway_x_angle = math.pi / 2.0 - bearing_rad
        cand_a = runway_x_angle % math.pi
        cand_b = (runway_x_angle + math.pi / 2.0) % math.pi

        def _min_piece_area(theta):
            ux, uy = math.cos(theta), math.sin(theta)
            line = _LS([(cx - span * ux, cy - span * uy),
                        (cx + span * ux, cy + span * uy)])
            try:
                result = _shp_split(polygon, line)
            except Exception:
                return -1.0
            geoms = (list(getattr(result, "geoms", []))
                     if result.geom_type != "Polygon" else [result])
            areas = [g.area for g in geoms
                     if g.geom_type == "Polygon" and not g.is_empty]
            return min(areas) if len(areas) >= 2 else 0.0

        angle_rad = (cand_a if _min_piece_area(cand_a) >= _min_piece_area(cand_b)
                     else cand_b)
    else:
        cs = [h.centroid for h in big_interiors]
        if len(cs) > 1:
            x_spread = max(c.x for c in cs) - min(c.x for c in cs)
            y_spread = max(c.y for c in cs) - min(c.y for c in cs)
            angle_rad = math.pi / 2.0 if x_spread > y_spread else 0.0

    dx = math.cos(angle_rad)
    dy = math.sin(angle_rad)
    cut = _LS([(cx - span * dx, cy - span * dy),
               (cx + span * dx, cy + span * dy)])

    try:
        result = _shp_split(polygon, cut)
    except Exception:
        # Fallback: emit the exterior with holes dropped.  Should
        # not occur for valid simple geometries.
        return [Polygon(polygon.exterior.coords)]
    pieces: List[Polygon] = []
    geoms = (list(getattr(result, "geoms", []))
             if result.geom_type != "Polygon" else [result])

    # Per user 2026-05-04: snap cut-induced vertices on each piece's
    # boundary to the nearest existing polygon vertex within 5 m
    # that's MORE aligned with the hole-centroid axis (smaller
    # perpendicular distance to the cut line).  This eliminates
    # mid-rect-edge nodes by replacing the cut crossing with an
    # existing on-axis vertex (typically a rect corner that was
    # already on the polygon's boundary via seam injection).
    SNAP_RADIUS_M = 5.0

    def _perp_to_cut(px, py):
        return abs((px - cx) * dy - (py - cy) * dx)

    # Pre-compute the natural crossing points so we can identify
    # cut-induced verts on each piece.
    try:
        natural_inter = polygon.exterior.intersection(cut)
    except Exception:
        natural_inter = None
    natural_pts: List[Tuple[float, float]] = []
    if natural_inter is not None and not natural_inter.is_empty:
        if natural_inter.geom_type == "Point":
            natural_pts = [(natural_inter.x, natural_inter.y)]
        elif natural_inter.geom_type == "MultiPoint":
            natural_pts = [(p.x, p.y) for p in natural_inter.geoms]
        elif natural_inter.geom_type == "GeometryCollection":
            for g in natural_inter.geoms:
                if g.geom_type == "Point":
                    natural_pts.append((g.x, g.y))
    # Existing polygon boundary verts (exterior + interiors) — these
    # are the candidates for snapping.
    pre_cut_verts: List[Tuple[float, float]] = []
    pe = list(polygon.exterior.coords)
    if pe and pe[0] == pe[-1]:
        pe = pe[:-1]
    pre_cut_verts.extend(pe)
    for h in polygon.interiors:
        hv = list(h.coords)
        if hv and hv[0] == hv[-1]:
            hv = hv[:-1]
        pre_cut_verts.extend(hv)

    def _snap_cut_verts(piece: Polygon) -> Polygon:
        """Replace each cut-induced vertex on this piece's boundary
        with the NEAREST existing pre-cut vertex within
        SNAP_RADIUS_M, provided that vertex is itself within
        ``ON_AXIS_TOL_M`` of the cut line (i.e. it's a candidate
        node that's already aligned with the hole-centroid axis).
        Prevents mid-rect-edge cut-induced verts when an existing
        rect corner sits near the cut line."""
        if not natural_pts:
            return piece
        ON_AXIS_TOL_M = 2.0  # max perp distance from cut line for a
                              # candidate to be considered "aligned"
        coords = list(piece.exterior.coords)
        if not coords:
            return piece
        had_close = (coords[0] == coords[-1])
        if had_close:
            coords = coords[:-1]
        modified = False
        for i, (vx, vy) in enumerate(coords):
            is_cut_vert = any(
                math.hypot(vx - nx, vy - ny) < 0.5
                for nx, ny in natural_pts)
            if not is_cut_vert:
                continue
            # Find the nearest existing pre-cut vertex within range
            # that's also on-axis.  Skip candidates farther from the
            # cut line than ON_AXIS_TOL_M — those would degrade
            # alignment.
            best_v: Optional[Tuple[float, float]] = None
            best_d = SNAP_RADIUS_M
            for px, py in pre_cut_verts:
                if abs(px - vx) > SNAP_RADIUS_M or abs(py - vy) > SNAP_RADIUS_M:
                    continue
                d = math.hypot(px - vx, py - vy)
                if d < 0.01:
                    continue  # same vertex
                if d > best_d:
                    continue
                if _perp_to_cut(px, py) > ON_AXIS_TOL_M:
                    continue
                best_d = d
                best_v = (float(px), float(py))
            if best_v is not None:
                coords[i] = best_v
                modified = True
        if not modified:
            return piece
        # Reconstruct the polygon, dedupe consecutive duplicates.
        deduped: List[Tuple[float, float]] = []
        for c in coords:
            if deduped and (math.hypot(c[0] - deduped[-1][0],
                                        c[1] - deduped[-1][1]) < 0.01):
                continue
            deduped.append(c)
        if len(deduped) < 3:
            return piece
        try:
            new_p = Polygon(deduped, [list(h.coords) for h in piece.interiors])
            if new_p.is_valid and not new_p.is_empty:
                return new_p
            fixed = new_p.buffer(0)
            if (fixed.geom_type == "Polygon" and not fixed.is_empty):
                return fixed
        except Exception:
            pass
        return piece

    for g in geoms:
        if g.geom_type != "Polygon" or g.is_empty:
            continue
        if g.area < min_area_m2:
            continue
        g = _snap_cut_verts(g)
        if g.geom_type != "Polygon" or g.is_empty:
            continue
        pieces.extend(_decompose_polygon_with_holes(
            g, min_area_m2=min_area_m2, max_depth=max_depth - 1,
            runway_axis_deg=runway_axis_deg,
            corner_snap_pts=corner_snap_pts,
            corner_snap_tol_m=corner_snap_tol_m))
    # Sliver clean-up: smart-cut alignment eliminates the most
    # egregious wide-band strips (5 m × 67 m, 10 m × 119 m) that
    # appeared with horizontal-only cuts, but recursive splits can
    # still leave narrow corner slivers when a hole's MRR axis
    # nearly parallels a polygon edge.  Merge anything thinner
    # than 12 m into its largest-shared-boundary neighbour so the
    # apron stays one continuous polygon and we don't get
    # cliff-rendering on thin corners.
    MIN_PIECE_THICKNESS_M = 12.0
    pieces = _merge_thin_decomposed_pieces(
        pieces, min_thickness_m=MIN_PIECE_THICKNESS_M)
    return pieces


def _polygon_min_thickness(poly: "Polygon") -> float:
    """Approximate minimum thickness of a polygon: half-width of
    the rotated minimum bounding rectangle.  Fast computation: try
    the polygon's minimum-rotated-rectangle and return the shorter
    side length."""
    try:
        mrr = poly.minimum_rotated_rectangle
        if mrr.is_empty or mrr.geom_type != "Polygon":
            return 0.0
        coords = list(mrr.exterior.coords)
        if len(coords) < 5:
            return 0.0
        sides = []
        for i in range(4):
            ax, ay = coords[i]
            bx, by = coords[i + 1]
            sides.append(math.hypot(bx - ax, by - ay))
        return min(sides)
    except Exception:
        return 0.0


def _merge_thin_decomposed_pieces(
        pieces: "List[Polygon]",
        min_thickness_m: float = 4.0,
        max_iters: int = 50,
        ) -> "List[Polygon]":
    """Merge any piece in ``pieces`` whose minimum-rotated-rectangle
    thickness is less than ``min_thickness_m`` into the neighbouring
    piece sharing the longest boundary.  Used by
    ``_decompose_polygon_with_holes`` to suppress 2 m-thick horizontal
    strips that arise when multiple holes have close y-centroids.
    Returns a possibly-shorter list with thin strips absorbed.
    """
    if not pieces:
        return pieces
    work = list(pieces)
    for _it in range(max_iters):
        # Find the thinnest piece below threshold.
        thin_idx: Optional[int] = None
        thin_thick = float('inf')
        for i, p in enumerate(work):
            if p is None or p.is_empty:
                continue
            t = _polygon_min_thickness(p)
            if t < min_thickness_m and t < thin_thick:
                thin_idx = i
                thin_thick = t
        if thin_idx is None:
            break
        thin = work[thin_idx]
        # Find the neighbour with the longest shared boundary.
        best_j: Optional[int] = None
        best_share = 0.0
        thin_boundary = thin.boundary
        for j, p in enumerate(work):
            if j == thin_idx or p is None or p.is_empty:
                continue
            try:
                shared = thin_boundary.intersection(
                    p.boundary).length
            except Exception:
                shared = 0.0
            if shared > best_share:
                best_share = shared
                best_j = j
        if best_j is None or best_share <= 0.0:
            # No neighbour to merge with; drop the thin piece so
            # it doesn't render as a cliff.
            work[thin_idx] = None
            continue
        try:
            merged = unary_union([thin, work[best_j]])
            if merged.is_empty:
                work[thin_idx] = None
                continue
            if merged.geom_type == "MultiPolygon":
                # Pick the largest piece — the union didn't fully
                # bridge.  Drop the thin one.
                work[thin_idx] = None
                continue
            if merged.geom_type != "Polygon":
                work[thin_idx] = None
                continue
            work[best_j] = merged
            work[thin_idx] = None
        except Exception:
            work[thin_idx] = None
    return [p for p in work if p is not None and not p.is_empty]


# ── Hole splicing ────────────────────────────────────────────────
#
# Defensive fallback for any polygon with holes that slips through
# decomposition (e.g. shapely.ops.split failed).  Splices each
# hole into the exterior via a zero-width bridge so ear-clipping
# can operate on a single ring; sliver triangles along the bridge
# are filtered downstream by the area threshold.


def _splice_holes(polygon: Polygon) -> List[Tuple[float, float]]:
    """Return the vertex list of the spliced single-ring polygon
    (without closing repeat).  Holes are inserted one at a time by
    finding the closest exterior vertex / hole vertex pair and
    threading the hole into the exterior at that bridge.
    """
    ext = list(polygon.exterior.coords)
    if ext and ext[0] == ext[-1]:
        ext = ext[:-1]
    holes_list: List[List[Tuple[float, float]]] = []
    for h in polygon.interiors:
        h_coords = list(h.coords)
        if h_coords and h_coords[0] == h_coords[-1]:
            h_coords = h_coords[:-1]
        if len(h_coords) >= 3:
            holes_list.append(h_coords)
    if not holes_list:
        return ext
    # Process holes from largest to smallest so big holes get the
    # "best" bridge slots; small holes thread into the still-clean
    # remainder.
    holes_list.sort(key=lambda h: -_polygon_area(h))
    ring = list(ext)
    for hole in holes_list:
        ring = _splice_one_hole(ring, hole)
    return ring


def _polygon_area(coords: Sequence[Tuple[float, float]]) -> float:
    s = 0.0
    n = len(coords)
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) * 0.5


def _splice_one_hole(ring: List[Tuple[float, float]],
                     hole: List[Tuple[float, float]]
                     ) -> List[Tuple[float, float]]:
    """Find the closest (ring_vertex, hole_vertex) pair and splice
    the hole into the ring at that bridge.  Hole is walked in its
    native (CW relative to a CCW exterior) direction so the spliced
    ring stays simple."""
    best_i = best_j = 0
    best_d2 = float("inf")
    for i, (rx, ry) in enumerate(ring):
        for j, (hx, hy) in enumerate(hole):
            d2 = (rx - hx) * (rx - hx) + (ry - hy) * (ry - hy)
            if d2 < best_d2:
                best_d2 = d2
                best_i, best_j = i, j
    # Spliced ring:
    #   ring[0..best_i] + hole[best_j..end] + hole[0..best_j]
    #   + ring[best_i..end]
    # The boundary touches ring[best_i] and hole[best_j] twice —
    # this is the bridge corridor.
    spliced: List[Tuple[float, float]] = []
    spliced.extend(ring[: best_i + 1])
    spliced.extend(hole[best_j:])
    spliced.extend(hole[: best_j + 1])
    spliced.extend(ring[best_i:])
    return spliced


# ── Ear-clip triangulation ───────────────────────────────────────
#
# Pure-Python ear-clipping for simple polygons (no holes).  Returns
# index triples (i,j,k) into the input vertex list.  For an N-vertex
# simple polygon, exactly N-2 triangles are produced — the proven
# minimum count when no Steiner points are added.


                             # boundary segment.  Long segments
                             # (e.g. cuts from hole-decomposition)
                             # leave Triangle4XP's quality refinement
                             # unable to fit small triangles between
                             # the cut endpoints; interior Steiners
                             # span the full distance and produce
                             # steep local gradients.  Densifying
                             # the boundary with linearly-interpolated
                             # midpoints gives Triangle4XP closer
                             # boundary anchors to connect to.


SHARED_NEIGHBOUR_EDGE_TOL_M = 5.0  # max distance from a midpoint
                                    # to a neighbour-shape edge to
                                    # treat the midpoint as on a
                                    # shared boundary; matches
                                    # check_grade's edge-step
                                    # search radius


def _densify_long_boundary_edges(
    ring: List[Tuple[float, float]],
    vert_elev: List[float],
    neighbour_edges: List[Tuple[float, float, float, float,
                                 float, float]],
    sloping_rect_edges: Optional[List[Tuple[float, float, float, float]]] = None,
    runway_edges: Optional[List[Tuple[float, float, float, float]]] = None,
    terminal_edges: Optional[List[Tuple[float, float, float, float]]] = None,
) -> Tuple[List[Tuple[float, float]], List[float]]:
    """Insert interpolated midpoints along long ring segments.

    For each candidate midpoint, check distance to the nearest
    rect/runway/terminal edge.  If within
    ``SHARED_NEIGHBOUR_EDGE_TOL_M``, the midpoint sits on a
    shared boundary; use the NEIGHBOUR'S edge-interpolated
    elevation (matches whatever the neighbour renders at that
    point — including segmented-runway piecewise profiles).
    Otherwise use linear interpolation between the segment's
    endpoint elevations.

    Per Rule 2 (user 2026-05-01) + user 2026-05-04 follow-up: if
    ``sloping_rect_edges`` is supplied, skip any midpoint within
    ``SLOPING_EDGE_SNAP_M`` of ANY edge of a sloping rect (both the
    sloping edges parallel to source_axis AND the cross edges
    perpendicular to it).  Junctions must share a sloping rect's
    boundary 1:1 — only the rect's 4 corners are legal shared
    vertices.  Without this guard, the densification midpoints
    appear as intermediate nodes on the cross edge between two
    corner-coincident vertices, then get pushed ~1 m off the edge
    by Rule 5's outward-push pass.

    Per Rule 1 (user 2026-05-01): if ``runway_edges`` is supplied,
    skip any midpoint within ``RUNWAY_BOUNDARY_TOL_M`` of a runway
    boundary edge — junction vertices on the runway must coincide
    with runway vertices, not float between them.

    Per user 2026-05-04: if ``terminal_edges`` is supplied, skip any
    midpoint within ``SHARED_VERTEX_TOL_M`` of a terminal edge.
    Terminals don't yet have an elevation when triangulation runs,
    so they're absent from ``neighbour_edges`` and the generic
    ``_point_on_neighbour`` guard misses them — junctions adjacent
    to a terminal pad would otherwise grow extra mid-edge vertices
    that the terminal itself doesn't have, breaking the seamless
    meld between the two surfaces.

    This guarantees no shared-boundary step is introduced and
    handles both the simple "junction-on-rect-edge" case (linear
    interp matches) and the "junction-on-segmented-runway-edge"
    case (use the runway segment's interp).
    """
    from ..config import SLOPING_EDGE_SNAP_M, RUNWAY_BOUNDARY_TOL_M
    from ..elevation import _corner_elevation_bucket
    n = len(ring)
    if n < 3 or len(vert_elev) != n:
        return ring, vert_elev

    def _interp_at(mx: float, my: float, fallback: float) -> float:
        """Return neighbour-edge interpolation at (mx, my) if a
        neighbour edge passes within tolerance; else fallback."""
        best_d2 = SHARED_NEIGHBOUR_EDGE_TOL_M * SHARED_NEIGHBOUR_EDGE_TOL_M
        best_e: Optional[float] = None
        for ax, ay, bx, by, ea, eb in neighbour_edges:
            dx = bx - ax
            dy = by - ay
            seg2 = dx * dx + dy * dy
            if seg2 < 0.04:
                continue
            t = ((mx - ax) * dx + (my - ay) * dy) / seg2
            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0
            cx = ax + t * dx
            cy = ay + t * dy
            d2 = (mx - cx) * (mx - cx) + (my - cy) * (my - cy)
            if d2 < best_d2:
                best_d2 = d2
                best_e = ea + t * (eb - ea)
        return best_e if best_e is not None else fallback

    def _point_within_of_edge(mx: float, my: float,
                              edges: List[Tuple[float, float, float, float]],
                              tol_m: float) -> bool:
        """True if (mx, my) lies within ``tol_m`` PERPENDICULAR to
        any edge in the supplied list, AND its projection falls
        strictly within the edge segment.  Per user 2026-05-01:
        vertices reaching toward a rect's cross-edge corner (whose
        projection lies past the sloping edge's endpoint) are
        allowed, so we exempt projections at or beyond either
        endpoint."""
        tol2 = tol_m * tol_m
        for ax, ay, bx, by in edges:
            dx = bx - ax
            dy = by - ay
            seg2 = dx * dx + dy * dy
            if seg2 < 0.04:
                continue
            t = ((mx - ax) * dx + (my - ay) * dy) / seg2
            if t <= 0.0 or t >= 1.0:
                continue
            cx = ax + t * dx
            cy = ay + t * dy
            d2 = (mx - cx) * (mx - cx) + (my - cy) * (my - cy)
            if d2 <= tol2:
                return True
        return False

    def _point_on_neighbour(mx: float, my: float) -> bool:
        """True if point (mx, my) lies within SHARED_VERTEX_TOL_M of
        some neighbour edge's interior.

        A densification midpoint that lands this close to a
        rect/runway/terminal edge will be snapped onto that edge by
        X-Plane's vector_map encroachment check, creating a T-
        junction that breaks the neighbour rect's 4-corner slope
        rendering.  Skip these midpoints — the junction edge stays
        un-densified at that point, which is the right call because
        Triangle4XP doesn't need an interior-anchor near a shared
        boundary anyway (the neighbour rect's own subdivision
        provides the anchors).
        """
        tol2 = SHARED_VERTEX_TOL_M * SHARED_VERTEX_TOL_M
        for nax, nay, nbx, nby, _, _ in neighbour_edges:
            ndx = nbx - nax
            ndy = nby - nay
            seg2 = ndx * ndx + ndy * ndy
            if seg2 < 0.04:
                continue
            t = ((mx - nax) * ndx + (my - nay) * ndy) / seg2
            if t < 0.0 or t > 1.0:
                continue
            cx = nax + t * ndx
            cy = nay + t * ndy
            d2 = (mx - cx) * (mx - cx) + (my - cy) * (my - cy)
            if d2 < tol2:
                return True
        return False

    # Pre-compute the bucket-key of every existing ring vertex so we
    # can reject midpoints that would collide with one via
    # ``to_osm``'s SHARED_VERTEX_TOL_M intern.  A midpoint that
    # collides would emit the same OSM nid as a non-adjacent ring
    # vertex, producing a polygon that visits the same node twice
    # — duplicate-consecutive or self-intersection (figure-8) at
    # OSM-write time.  Either crashes X-Plane's mesh builder.
    existing_buckets = {
        _corner_elevation_bucket(x, y) for (x, y) in ring}
    new_ring: List[Tuple[float, float]] = []
    new_elev: List[float] = []
    for i in range(n):
        a = ring[i]
        b = ring[(i + 1) % n]
        ea = vert_elev[i]
        eb = vert_elev[(i + 1) % n]
        new_ring.append(a)
        new_elev.append(ea)
        d = math.hypot(b[0] - a[0], b[1] - a[1])
        if d <= MAX_BOUNDARY_EDGE_M:
            continue
        n_subs = int(math.ceil(d / MAX_BOUNDARY_EDGE_M))
        # Per-arc cap of 4 — at most 3 inserted midpoints per ring
        # edge.  Without this, the long apt.dat boundary edges that
        # bound a residue junction (CYXY's -10070 had 200+ m edges)
        # add 6-8 midpoints each and total junction-vertex count
        # explodes past Triangle4XP's safe ceiling.  4 anchors per
        # arc is enough for the elevation solver to fit a smooth
        # plane; more just adds free vertices that the solver pushes
        # around until they create cliffs.
        if n_subs > 4:
            n_subs = 4
        for k in range(1, n_subs):
            t = k / n_subs
            mx = a[0] + t * (b[0] - a[0])
            my = a[1] + t * (b[1] - a[1])
            mb = _corner_elevation_bucket(mx, my)
            if mb in existing_buckets:
                continue  # would collide with an existing ring nid
            # Skip midpoints that would land near a neighbour edge
            # — X-Plane's vector_map would snap them onto the
            # neighbour, creating a T-junction that breaks the
            # neighbour rect's 4-corner slope rendering.
            if _point_on_neighbour(mx, my):
                continue
            # Per Rule 2 (user 2026-05-01) + 1:1-corner-sharing rule
            # (user 2026-05-04): no junction midpoint near a sloping
            # rect edge — junctions share corners only.  Tightened
            # from SLOPING_EDGE_SNAP_M (20 m) to 10 m on 2026-05-05
            # so densification can place ring anchors on long edges
            # that pass within the 10–20 m corridor of an adjacent
            # F-taxi rect (otherwise large junctions retain unanchored
            # 200 m+ stretches that let DEM spikes through).
            DENSIFY_SLOPING_RECT_TOL_M = 10.0
            if sloping_rect_edges and _point_within_of_edge(
                    mx, my, sloping_rect_edges,
                    DENSIFY_SLOPING_RECT_TOL_M):
                continue
            # Per Rule 1 (user 2026-05-01): no junction vertex within
            # ``RUNWAY_BOUNDARY_TOL_M`` of a runway boundary edge.
            if runway_edges and _point_within_of_edge(
                    mx, my, runway_edges, RUNWAY_BOUNDARY_TOL_M):
                continue
            # Per user 2026-05-04: no junction vertex on a terminal
            # edge interior — terminals lack altitude at triangulation
            # time and thus don't appear in ``neighbour_edges``.
            if terminal_edges and _point_within_of_edge(
                    mx, my, terminal_edges, SHARED_VERTEX_TOL_M):
                continue
            linear_me = ea + t * (eb - ea)
            me = _interp_at(mx, my, linear_me)
            new_ring.append((mx, my))
            new_elev.append(me)
            existing_buckets.add(mb)
    return new_ring, new_elev


COLINEAR_DROP_M = 0.5  # max perpendicular distance to neighbours
                        # for a non-anchor vertex to be removed.
                        # Originally 0.5 m; bumped to 3.0 m on
                        # 2026-04-25 to thin apt.dat curves; reverted
                        # to 0.5 m on 2026-05-05 once the geometry
                        # baseline produced precise pavement curves
                        # the user wants preserved.  Still strips
                        # near-duplicate vertices (sub-metre noise)
                        # without flattening genuine curve detail.


                             # ring vertex to a non-adjacent edge
                             # of the same polygon for the vertex
                             # to count as a "spike" — the ring
                             # ventured out and returned to (or
                             # very near) itself.  Sub-mm spikes
                             # are valid in shapely's eyes but
                             # become hard self-intersections
                             # after .11f OSM-format truncation,
                             # which would crash X-Plane.



def _drop_sliver_corners(
    ring: List[Tuple[float, float]],
) -> List[Tuple[float, float]]:
    """Drop ring vertices whose interior angle is below
    ``SLIVER_ANGLE_THRESHOLD_DEG``.

    A sliver corner is a needle-tip vertex: the polygon comes in
    along one edge, makes a near-180° fold, and goes back out
    almost on top of the incoming edge, leaving a thin wedge.
    Source: residue construction (apt.dat pav − rects − terminals)
    leaves wedges where two rect/terminal edges meet the boundary at
    nearly-collinear angles.  Shapely calls these polygons valid
    (the two long edges are parallel-but-not-equal, no self-
    intersection), but the polygon's needle-tip corner forces
    Triangle4XP downstream to emit at least one corner triangle
    with that interior angle — for sub-2° tips that triangle is
    near-degenerate (NaN normal) and crashes X-Plane's mesh builder.

    Dropping the tip vertex collapses the wedge into a single edge
    between the two flanking vertices.  Coverage cost: the area
    between the tip and the truncation chord (typically < 50 m²).

    Iterates to a fixed point — dropping one tip can expose
    another.  Capped at 8 passes.
    """
    if len(ring) < 4:
        return ring
    cos_thresh = math.cos(math.radians(SLIVER_ANGLE_THRESHOLD_DEG))
    for _ in range(8):
        n = len(ring)
        if n < 4:
            break
        keep = [True] * n
        for i in range(n):
            ax, ay = ring[(i - 1) % n]
            bx, by = ring[i]
            cx, cy = ring[(i + 1) % n]
            v1x, v1y = ax - bx, ay - by
            v2x, v2y = cx - bx, cy - by
            n1 = math.hypot(v1x, v1y)
            n2 = math.hypot(v2x, v2y)
            if n1 < 1e-9 or n2 < 1e-9:
                continue
            cos = (v1x * v2x + v1y * v2y) / (n1 * n2)
            if cos > cos_thresh:
                # Angle = acos(cos) is below threshold.
                keep[i] = False
        new_ring = [r for r, k in zip(ring, keep) if k]
        if len(new_ring) == n:
            break
        ring = new_ring
    return ring


def _drop_colinear_boundary_vertices(
    ring: List[Tuple[float, float]],
    corner_elev: Dict[Tuple[int, int], float],
    shared_junction_buckets: set,
) -> List[Tuple[float, float]]:
    """Remove vertices whose perpendicular distance to the line
    through their immediate neighbours is below COLINEAR_DROP_M,
    EXCEPT vertices that are shared corners (rect/runway/terminal
    or cross-junction shared buckets) — those carry topological
    meaning and must not be dropped.

    Eliminates the "ear-clip can only emit a sliver here" geometry
    that produces visible step artefacts on long thin apron strips.
    """
    from ..elevation import _corner_elevation_bucket
    if len(ring) < 4:
        return ring
    # Iterate to fixed point: dropping one vertex may make a
    # neighbour droppable too.
    for _ in range(8):
        n = len(ring)
        if n < 4:
            break
        keep = [True] * n
        for i in range(n):
            bucket = _corner_elevation_bucket(*ring[i])
            if bucket in corner_elev or bucket in shared_junction_buckets:
                continue  # anchor — keep no matter what
            ax, ay = ring[(i - 1) % n]
            bx, by = ring[(i + 1) % n]
            cx, cy = ring[i]
            # Perpendicular distance from C to line AB.
            dx = bx - ax
            dy = by - ay
            seg_len = math.hypot(dx, dy)
            if seg_len < 0.1:
                continue
            # Cross-product / line-length = perpendicular distance.
            perp = abs((cx - ax) * dy - (cy - ay) * dx) / seg_len
            if perp < COLINEAR_DROP_M:
                keep[i] = False
        new_ring = [c for c, k in zip(ring, keep) if k]
        if len(new_ring) == n:
            break
        ring = new_ring
    return ring



def _find_junction_points(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    osm_centerlines: Optional[List[Tuple[LineString, str]]] = None,
) -> List[Tuple[float, float]]:
    """Identify junction POINTS — OSM nodes shared by ≥ 2 DIFFERENT refs.

    Per user rule: pure bends within one taxi (two same-ref ways
    meeting at a node) are NOT junctions — they emit as adjacent
    same-role rects sharing a vertex line.  Only nodes where
    multiple distinct refs (or an unrefed way + a refed way) meet
    are junction candidates.

    Candidates within ``JUNCTION_CLUSTER_DIST_M`` of each other
    collapse into one cluster; the cluster centroid is the
    junction point.
    """
    from collections import defaultdict
    refs_at_node: Dict[str, set] = defaultdict(set)
    refed_taxi_nodes: set = set()
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        ref = tags.get("ref", "")
        if not ref:
            continue
        refed_taxi_nodes.update(nds)
        for n in nds:
            refs_at_node[n].add(ref)

    # Per user 2026-05-05: any ``aeroway=taxiway`` (refed OR not)
    # AND any ``aeroway=parking_position`` way contributes a
    # connector node where it meets a refed taxiway.  At SPJC F
    # crosses 2 parking_position centerlines that split it into 3
    # rects; CYXY's V↔U pair has 6 short unrefed taxiway connectors
    # marking the chart-level intersections.  No length filter —
    # any aeroway crossing splits.  Apron-only nodes are still
    # filtered because the connector tag is added ONLY when the
    # node is shared with a refed taxiway.
    for wid, nds, tags in ways:
        aw = tags.get("aeroway")
        if aw == "taxiway" and tags.get("ref", ""):
            continue  # already handled in the refed-taxiway pass above
        if aw not in ("taxiway", "parking_position"):
            continue
        for n in nds:
            if n in refed_taxi_nodes:
                refs_at_node[n].add("_conn")

    candidates: List[Tuple[float, float]] = []
    for nid, refs in refs_at_node.items():
        if len(refs) < 2:
            continue
        if nid not in nodes:
            continue
        lat, lon = nodes[nid]
        candidates.append(to_m(lon, lat))

    # ALSO add geometric crossing points between different-ref
    # centerlines (helps SPLP where few OSM nodes are shared).
    if osm_centerlines:
        for i in range(len(osm_centerlines)):
            ls1, ref1 = osm_centerlines[i]
            for j in range(i+1, len(osm_centerlines)):
                ls2, ref2 = osm_centerlines[j]
                if ref1 and ref2 and ref1 == ref2:
                    continue
                if not ls1.intersects(ls2):
                    continue
                try:
                    inter = ls1.intersection(ls2)
                except Exception:
                    continue
                if inter.is_empty:
                    continue
                if inter.geom_type == "Point":
                    candidates.append((inter.x, inter.y))
                elif inter.geom_type == "MultiPoint":
                    for p in inter.geoms:
                        candidates.append((p.x, p.y))
                elif inter.geom_type == "LineString":
                    candidates.append(inter.centroid.coords[0])

    # Cluster within JUNCTION_CLUSTER_DIST_M (greedy single-link)
    clusters: List[List[Tuple[float, float]]] = []
    for pt in candidates:
        placed = False
        for cl in clusters:
            if any(math.hypot(pt[0]-q[0], pt[1]-q[1]) <= JUNCTION_CLUSTER_DIST_M
                   for q in cl):
                cl.append(pt)
                placed = True
                break
        if not placed:
            clusters.append([pt])

    return [(sum(p[0] for p in cl)/len(cl),
             sum(p[1] for p in cl)/len(cl)) for cl in clusters]


def _build_junctions_from_rect_endpoints(
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    merge_dist: float,
    pav_union: Optional[Polygon],
    terminal_union: Optional[Polygon] = None,
) -> List[Polygon]:
    """Build junctions from rect endpoint clusters (user's approach).

    Algorithm:
      1. Collect each rect's 2 axis endpoints + 2 corner vertices at
         each end (total: 2 endpoints × 2 corners = 4 corners per rect).
      2. Cluster axis endpoints by single-link within ``merge_dist``.
      3. For each cluster of ≥ 2 endpoints:
         * If all endpoints share the same ref → same-taxi bend (no
           junction emitted; same-ref rects connect via shared vertices
           handled elsewhere).
         * Else → emit a junction polygon whose vertices are the 2
           corner vertices of each participating rect at the cluster.

    The polygon vertices are ordered angularly around the cluster
    centroid, giving a star polygon that wraps through each rect's
    corner pair.
    """
    if not taxi_rects:
        return []

    # Endpoint records: (rect_idx, end_index, axis_pt, corner_pair, ref)
    endpoints = []
    for i, (rect, axis, role, ref) in enumerate(taxi_rects):
        pairs = _rect_end_corners(rect, axis)
        if len(pairs) < 2:
            continue
        coords = list(axis.coords)
        endpoints.append((i, 0, coords[0], pairs[0], ref))
        endpoints.append((i, 1, coords[-1], pairs[1], ref))

    # Single-link cluster by axis-endpoint proximity
    n = len(endpoints)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i+1, n):
            ax = endpoints[i][2]
            bx = endpoints[j][2]
            if math.hypot(ax[0]-bx[0], ax[1]-bx[1]) <= merge_dist:
                union(i, j)

    clusters: Dict[int, List[int]] = {}
    for i in range(n):
        r = find(i)
        clusters.setdefault(r, []).append(i)

    out: List[Polygon] = []
    for cl in clusters.values():
        if len(cl) < 2:
            continue
        # Unique refs in cluster
        refs = {endpoints[i][4] for i in cl}
        # Unique rect ids (cluster can contain multiple ends of same rect)
        unique_rects = {endpoints[i][0] for i in cl}
        # Pure same-ref bend (only one ref AND only one rect pairs bends) — skip junction
        if len(refs) == 1 and len(unique_rects) <= 2:
            continue
        # Collect all corner vertices
        all_corners = []
        for i in cl:
            c1, c2 = endpoints[i][3]
            all_corners.append(c1)
            all_corners.append(c2)
        if len(all_corners) < 3:
            continue
        # Deduplicate near-identical corners
        uniq: List[Tuple[float, float]] = []
        for c in all_corners:
            if not any(math.hypot(c[0]-u[0], c[1]-u[1]) < 0.1 for u in uniq):
                uniq.append(c)
        if len(uniq) < 3:
            continue
        cx = sum(p[0] for p in uniq) / len(uniq)
        cy = sum(p[1] for p in uniq) / len(uniq)
        ordered = sorted(uniq, key=lambda p: math.atan2(p[1]-cy, p[0]-cx))
        try:
            poly = Polygon(ordered).buffer(0)
        except Exception:
            continue
        if poly.is_empty:
            continue
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon" or poly.area < 100.0:
            continue
        # Don't let junction bleed into a terminal
        if terminal_union is not None:
            try:
                poly = poly.difference(terminal_union)
            except Exception:
                pass
            if poly.is_empty:
                continue
            if poly.geom_type == "MultiPolygon":
                poly = max(poly.geoms, key=lambda g: g.area)
            if poly.geom_type != "Polygon":
                continue
        out.append(poly)
    return out


def _rect_end_corners(rect: Polygon, axis: LineString
                      ) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Return the 2 pairs of corners at the rect's 2 short ends.

    For a 4-corner rect built as [p1+h*perp, p2+h*perp, p2-h*perp,
    p1-h*perp], the "start" end has corners [idx 0, idx 3] and
    the "end" end has corners [idx 1, idx 2].  This pairing
    matters because junction polygons are built from the pair at
    whichever end of the rect meets the junction centroid.
    """
    coords = list(rect.exterior.coords)
    if len(coords) < 5:
        return []
    # idx 0 = p1+perp, idx 1 = p2+perp, idx 2 = p2-perp, idx 3 = p1-perp
    start_pair = (coords[0], coords[3])   # corners at axis start
    end_pair = (coords[1], coords[2])     # corners at axis end
    return [start_pair, end_pair]


def _build_junction_constructive(
    cluster_centroid: Tuple[float, float],
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    pav_union,
    terminal_union,
    max_corner_dist_m: float = 80.0,
    local_disc_radius_m: float = 120.0,
) -> Optional[Polygon]:
    """Constructive junction-polygon build per the user's
    authoritative shape rule
    (memory: feedback_shape_rules):

      1. One vertex per incoming rect corner (outer corners at the
         rect end that meets this junction).
      2. Between consecutive corners belonging to DIFFERENT rects,
         trace apt.dat pavement vertices that lie on the local pav
         boundary arc between them.
      3. Between consecutive corners of the SAME rect, connect
         directly (that's the rect's short edge — the junction's
         inner face on that side).

    Arcs are bounded by clipping pav_union to a disc of radius
    ``local_disc_radius_m`` around the cluster centroid; this keeps
    the boundary walk local and avoids wrap-around issues at large
    multi-component pavements.

    Returns ``None`` if fewer than 2 rect ends cluster here or
    construction fails.
    """
    if pav_union is None or pav_union.is_empty:
        return None
    cx, cy = cluster_centroid
    cp = Point(cx, cy)
    if terminal_union is not None and not terminal_union.is_empty:
        try:
            if terminal_union.contains(cp):
                return None
        except Exception:
            pass

    # Gather rect ends near the cluster centroid.
    rect_ends: List[Tuple[int, Tuple[float, float], Tuple[float, float]]] = []
    for i, (rect, axis, role, ref) in enumerate(taxi_rects):
        coords_ax = list(axis.coords)
        if len(coords_ax) < 2:
            continue
        ax_start = coords_ax[0]
        ax_end = coords_ax[-1]
        d_start = math.hypot(ax_start[0] - cx, ax_start[1] - cy)
        d_end = math.hypot(ax_end[0] - cx, ax_end[1] - cy)
        if min(d_start, d_end) > max_corner_dist_m:
            continue
        pairs = _rect_end_corners(rect, axis)
        if len(pairs) < 2:
            continue
        idx = 0 if d_start <= d_end else 1
        c1, c2 = pairs[idx]
        rect_ends.append((i, c1, c2))

    if len(rect_ends) < 2:
        return None

    # Local pavement: intersect pav with a disc around the cluster.
    try:
        disc = cp.buffer(local_disc_radius_m)
        local_pav = pav_union.intersection(disc)
    except Exception:
        return None
    if local_pav.is_empty:
        return None
    if local_pav.geom_type == "MultiPolygon":
        # Take the component containing (or closest to) the centroid.
        best = None
        best_d = float('inf')
        for g in local_pav.geoms:
            if g.geom_type != "Polygon":
                continue
            d = g.distance(cp)
            if d < best_d:
                best_d = d
                best = g
        if best is None:
            return None
        local_pav = best
    if local_pav.geom_type != "Polygon":
        return None

    # Projection helper: use the exterior ring as a line for param.
    ext_coords = list(local_pav.exterior.coords)
    if len(ext_coords) < 4:
        return None
    # LinearRing is closed (first == last); use as LineString for project()
    ext_ls = LineString(ext_coords)
    ext_length = ext_ls.length
    if ext_length <= 0:
        return None

    # Project each rect outer corner onto the exterior ring.
    corner_data: List[Tuple[float, Tuple[float, float], int]] = []
    for rect_idx, c1, c2 in rect_ends:
        for c in (c1, c2):
            try:
                param = ext_ls.project(Point(c))
            except Exception:
                continue
            proj_pt = ext_ls.interpolate(param)
            # If the projection is far (corner inside pav interior,
            # not on boundary — e.g. the local disc cut through a
            # rect short edge), use the original corner position.
            if proj_pt.distance(Point(c)) > 20.0:
                continue
            corner_data.append((param, c, rect_idx))

    # Need at least 3 corners for a polygon.
    if len(corner_data) < 3:
        return None

    # Sort corners by boundary param.  This orders them around the
    # local pav exterior ring; same-rect corner pairs will typically
    # be adjacent (a rect's 2 outer corners land close together on
    # the boundary).
    corner_data.sort(key=lambda t: t[0])

    # Collect all exterior-ring vertices with their params for arc walks.
    ext_verts = ext_coords[:-1] if ext_coords[0] == ext_coords[-1] else ext_coords
    ext_vert_params: List[Tuple[float, Tuple[float, float]]] = []
    acc = 0.0
    for i, v in enumerate(ext_verts):
        if i > 0:
            acc += math.hypot(v[0] - ext_verts[i-1][0],
                              v[1] - ext_verts[i-1][1])
        ext_vert_params.append((acc, (v[0], v[1])))

    # Build polygon by walking sorted corners.  Between consecutive
    # corners of different rects, insert all exterior-ring vertices
    # whose params lie between the 2 corner params (forward direction,
    # with wrap-around from last to first).
    n = len(corner_data)
    poly_coords: List[Tuple[float, float]] = []
    for i in range(n):
        cur_param, cur_xy, cur_rect = corner_data[i]
        nxt_param, nxt_xy, nxt_rect = corner_data[(i + 1) % n]
        poly_coords.append(cur_xy)
        if cur_rect == nxt_rect:
            # Same-rect: direct connection (rect short edge), no arc.
            continue
        # Different rect: walk exterior ring from cur_param to nxt_param
        # in increasing-param direction (wrap at end).
        if i == n - 1 or nxt_param < cur_param:
            arc_verts = (
                [v for (vp, v) in ext_vert_params if vp > cur_param] +
                [v for (vp, v) in ext_vert_params if vp < nxt_param])
        else:
            arc_verts = [v for (vp, v) in ext_vert_params
                         if cur_param < vp < nxt_param]
        for v in arc_verts:
            poly_coords.append(v)

    if len(poly_coords) < 3:
        return None
    try:
        poly = Polygon(poly_coords).buffer(0)
    except Exception:
        return None
    if poly.is_empty:
        return None
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    if poly.geom_type != "Polygon" or poly.area < 80.0:
        return None

    # Clip to the local pav (in case rect short-edges extend past it)
    try:
        poly = poly.intersection(local_pav)
    except Exception:
        pass
    if poly.is_empty:
        return None
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    if poly.geom_type != "Polygon":
        return None

    # Exclude terminal overlap
    if terminal_union is not None and not terminal_union.is_empty:
        try:
            poly = poly.difference(terminal_union)
        except Exception:
            pass
        if poly.is_empty or poly.geom_type not in ("Polygon", "MultiPolygon"):
            return None
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon":
            return None

    return poly


def _build_junction_polys_from_corners(
    junction_points: List[Tuple[float, float]],
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    pav: Optional[Polygon],
    terminal_union: Optional[Polygon] = None,
    max_corner_dist_m: float = 100.0,
) -> List[Polygon]:
    """Build each junction polygon from the CORNER VERTICES of the
    adjacent rects, per the user's rule:

        "One vertex for each corner vertex of the rects it's joining."

    For each cluster centroid:
      1. Gather rects whose axis start/end is within
         ``max_corner_dist_m`` of the centroid.
      2. For each gathered rect, take the 2 corner vertices at the
         NEARER end (start or end).
      3. Order all corner vertices angularly around the centroid.
      4. Emit as the junction polygon (simple polygon through these
         vertices).

    If fewer than 2 rect ends cluster here, no junction polygon is
    emitted (would be degenerate).

    Junction polys that overlap a terminal aren't emitted; terminal
    boundary is adjacent-direct (per user rule "aprons can join
    directly to taxiways without a junction").
    """
    if pav is None or not junction_points:
        return []

    polys: List[Polygon] = []
    for (cx, cy) in junction_points:
        jpt = Point(cx, cy)
        # Avoid emitting into a terminal
        if terminal_union is not None and terminal_union.contains(jpt):
            continue

        corner_pts: List[Tuple[float, float]] = []
        for rect, axis, role, ref in taxi_rects:
            coords_ax = list(axis.coords)
            if len(coords_ax) < 2:
                continue
            ax_start = coords_ax[0]
            ax_end = coords_ax[-1]
            d_start = math.hypot(ax_start[0]-cx, ax_start[1]-cy)
            d_end = math.hypot(ax_end[0]-cx, ax_end[1]-cy)
            near_d = min(d_start, d_end)
            if near_d > max_corner_dist_m:
                continue
            pairs = _rect_end_corners(rect, axis)
            if not pairs:
                continue
            # Take the pair at the closer axis end
            idx = 0 if d_start <= d_end else 1
            corner_pts.extend(pairs[idx])

        if len(corner_pts) < 3:
            continue
        # Order corners angularly around centroid
        ordered = sorted(corner_pts,
                         key=lambda p: math.atan2(p[1]-cy, p[0]-cx))
        try:
            poly = Polygon(ordered).buffer(0)
        except Exception:
            continue
        if poly.is_empty:
            continue
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon":
            continue
        if poly.area < 80.0:
            continue
        # Clip to pavement so junction doesn't escape the pavement footprint
        poly = poly.intersection(pav)
        if poly.is_empty:
            continue
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon":
            continue
        polys.append(poly)
    return polys
