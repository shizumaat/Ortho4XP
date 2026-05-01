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

from O4_Pavement_Config import (
    JUNCTION_CLUSTER_DIST_M,
    MAX_BOUNDARY_EDGE_M,
    SLIVER_ANGLE_THRESHOLD_DEG,
)
from O4_Pavement_Layout import (
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
                                  max_depth: int = 8
                                  ) -> List[Polygon]:
    """Return a list of simple (no-hole) polygons that tile the
    same area as ``polygon``.  Cuts through each hole's centroid
    along a direction chosen to follow the hole's natural axes (or
    the centroid-spread of multiple holes), then recurses on each
    side.

    Per user 2026-04-28: previous implementation always cut
    horizontally, which produced thin (~2–12 m) strips between
    parallel taxi-rect holes whose y-centroids happened to be
    close.  Those strips were arbitrary geometric artifacts —
    their cut edges didn't follow any real pavement feature —
    and required a downstream merge band-aid to suppress.  The
    new policy:

      * Multiple holes ⇒ cut PERPENDICULAR to the centroid-spread
        direction.  When holes are spread along x (typical
        parallel-taxiway apron), the cut runs vertically BETWEEN
        the holes rather than horizontally between two close
        cut lines.
      * Single hole ⇒ cut along the hole's MRR long axis.  The
        cut is collinear with the rect's natural orientation and
        produces two pieces straddling the rect rather than
        slicing it across.

    Both rules align cut lines with real geometric features
    (rect axes, hole-cluster alignment) instead of an arbitrary
    horizontal direction.  The thin-piece merge fallback is kept
    at a small threshold purely for floating-point sliver clean-up.
    """
    from shapely.ops import split as _shp_split
    from shapely.geometry import LineString as _LS

    if (polygon.is_empty or polygon.geom_type != "Polygon"):
        return []
    if not polygon.interiors:
        return [polygon]
    if max_depth <= 0:
        # Recursion guard — emit the exterior with holes dropped
        # rather than retry forever.  Should never trigger for
        # realistic airport geometry.
        return [Polygon(polygon.exterior.coords)]
    # Pick the largest remaining hole — it'll be on one side of
    # the cut after slicing.
    interiors = list(polygon.interiors)
    interiors.sort(key=lambda h: -Polygon(h).area)
    hole = interiors[0]
    cx = float(hole.centroid.x)
    cy = float(hole.centroid.y)
    minx, miny, maxx, maxy = polygon.bounds

    # Decide cut direction.  The cut is a line through the hole
    # centroid; ``angle_rad`` is the angle of the cut LINE (not
    # the perpendicular), measured from +x.  Default = horizontal.
    angle_rad = 0.0
    if len(interiors) >= 2:
        # Multi-hole: cut perpendicular to the centroid-spread
        # axis so the cut SEPARATES the holes rather than slicing
        # parallel to their alignment.  E.g. holes spread along
        # +x → cut vertically (angle = π/2) between them.
        cs = [h.centroid for h in interiors]
        x_spread = max(c.x for c in cs) - min(c.x for c in cs)
        y_spread = max(c.y for c in cs) - min(c.y for c in cs)
        # We want the cut LINE to be perpendicular to the spread
        # direction.  spread_along_x ⇒ cut vertical (π/2);
        # spread_along_y ⇒ cut horizontal (0).
        if x_spread > y_spread:
            angle_rad = math.pi / 2.0  # vertical
        else:
            angle_rad = 0.0            # horizontal
    else:
        # Single hole: cut along its MRR long-axis direction.
        try:
            mrr = hole.minimum_rotated_rectangle
            if (mrr is not None and not mrr.is_empty
                    and mrr.geom_type == "Polygon"):
                mc = list(mrr.exterior.coords)
                if len(mc) >= 5:
                    sides = []
                    for i in range(4):
                        ax, ay = mc[i]
                        bx, by = mc[i + 1]
                        sides.append(
                            (math.hypot(bx - ax, by - ay),
                             math.atan2(by - ay, bx - ax)))
                    sides.sort(reverse=True)
                    angle_rad = sides[0][1]
        except Exception:
            angle_rad = 0.0

    # Build a cut line through (cx, cy) at angle_rad, extended well
    # past the polygon bounds on both sides.
    span = max(maxx - minx, maxy - miny) + 2.0
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
    for g in geoms:
        if g.geom_type != "Polygon" or g.is_empty:
            continue
        if g.area < min_area_m2:
            continue
        pieces.extend(_decompose_polygon_with_holes(
            g, min_area_m2=min_area_m2, max_depth=max_depth - 1))
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

    This guarantees no shared-boundary step is introduced and
    handles both the simple "junction-on-rect-edge" case (linear
    interp matches) and the "junction-on-segmented-runway-edge"
    case (use the runway segment's interp).
    """
    from O4_Airport_Pavement_Builder import _corner_elevation_bucket  # lazy: avoids circular import
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
            linear_me = ea + t * (eb - ea)
            me = _interp_at(mx, my, linear_me)
            new_ring.append((mx, my))
            new_elev.append(me)
            existing_buckets.add(mb)
    return new_ring, new_elev


COLINEAR_DROP_M = 3.0  # max perpendicular distance to neighbours
                        # for a non-anchor vertex to be removed.
                        # Bumped from 0.5 m → 3.0 m (user 2026-04-25)
                        # to thin apt.dat boundary curves more
                        # aggressively — eliminates most ear-clip
                        # sliver triangles whose 3 anchored vertices
                        # are nearly colinear.  Preserves rect /
                        # runway / terminal corners (sharp 90° turns
                        # are well above this threshold by definition)
                        # and any vertex shared between junctions.


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
    from O4_Airport_Pavement_Builder import _corner_elevation_bucket  # lazy: avoids circular import
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

    # Second pass: unrefed taxi ways act as CONNECTORS between
    # refed taxis at SPJC (e.g. 6 short unrefed ways bridge V to U,
    # marking the chart-level V↔U intersection points that split V
    # into multiple rects).  Only contribute a connector node when
    # it's also on a refed taxi way — this filters pure apron-area
    # markings (whose nodes touch only other unrefed ways).  Length
    # cap filters the long apron-boundary unrefed ways.
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        if tags.get("ref", ""):
            continue
        # Compute unrefed way length
        path_len = 0.0
        prev = None
        for n in nds:
            if n not in nodes:
                continue
            lat, lon = nodes[n]
            cur = to_m(lon, lat)
            if prev is not None:
                path_len += math.hypot(cur[0] - prev[0], cur[1] - prev[1])
            prev = cur
        if path_len > 200.0 or path_len < 20.0:
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
    max_arc_vertices: int = 4,
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
        # Per user rule (2026-04-18): "use a maximum of 4 points
        # between each" rect corner pair.  Sub-sample evenly.
        if len(arc_verts) > max_arc_vertices:
            step = len(arc_verts) / max_arc_vertices
            arc_verts = [arc_verts[int(k * step)]
                         for k in range(max_arc_vertices)]
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
