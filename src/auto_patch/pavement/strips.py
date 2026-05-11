"""Classify apt.dat pavement polygons into taxi and apron shapes.

This is the decomposition layer for the **apron-deformation**
pavement model (documented in STATUS.md § "Active plan").  The
rules, in full, are:

* Each apt.dat pavement polygon becomes **one `Shape`** (no
  buffering, no junction discs, no pre-carving).
* A polygon whose minimum-rotated-rectangle short side is ≤
  ``narrow_width_m`` (default 30 m) AND whose aspect ratio is ≥
  ``simple_strip_aspect`` (default 2) classifies as **taxi** with
  an axial slope direction (MRR midline, extended to the polygon
  boundary so its endpoints sit on shared edges with neighbours).
* Every other pavement polygon is an **apron**, a 2D smoothed
  surface.  Aprons have no axis.
* Adjacent apron polygons are **unioned** into connected components
  to reduce shape count — their shared boundaries don't carry
  elevation differences and the solver treats them as one surface.
* Taxi polygons are kept individually because each has its own
  axial slope direction.
* Runway polygons are the caller's concern (the pipeline subtracts
  `rwy_union_m` before handing polygons here).

Output is a flat tuple of :class:`Shape`.  There is no `Junction`
type in this model; compound-slope transitions are handled by the
downstream emitter as up-to-30 m triangle zones at the ends of
taxi polygons, and by apron surfaces reshaping to meet neighbour
elevations.  See STATUS.md for the full rationale.

The module is pure: no file I/O, no global state.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, Polygon
from shapely.ops import substring, unary_union
from shapely.validation import make_valid

from . import taxiway_skeleton as TS

# Narrow exception tuple for shapely / numeric-geometry failure
# modes.  Programming errors propagate so they surface immediately.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)


# ──────────────────────────────────────────────────────────────────
# Tunables
# ──────────────────────────────────────────────────────────────────
# Pavement up to this short-side width classifies as taxiway;
# anything wider is apron.  Set by user rule: "anything 45 m or
# less wide is a taxiway, anything larger is apron".
NARROW_WIDTH_M = 45.0

# Minimum MRR aspect for a polygon to be a taxi.  Blocky pavement
# (aspect < 2) has no dominant axis and behaves like an apron even
# if narrow.
SIMPLE_STRIP_ASPECT = 2.0

# Cap the perpendicular width probe.  A strip cannot be meaningfully
# wider than this plus a small margin; rays that hit nothing inside
# this cap are clamped.  Kept conservative so a pathological input
# cannot claim apron-scale width.
WIDTH_PROBE_MAX_REACH_M = NARROW_WIDTH_M + 5.0

# Minimum apron connected-component area; smaller residuals from
# unioning are discarded as slivers.
MIN_APRON_AREA_M2 = 10.0

# Minimum Voronoi skeleton branch length to even consider as a
# taxi candidate in a mega-polygon decomposition.  Shorter branches
# are usually boundary-concavity artifacts, not real taxis.
SKELETON_MIN_BRANCH_M = 50.0

# Minimum taxi polygon area.  After a skeleton branch is buffered
# and clipped to the source polygon, if the result is smaller than
# this it's not emitted — probably a fillet already covered by a
# longer neighbouring branch.
MIN_TAXI_AREA_M2 = 80.0

# Padding added to the local-half-width buffer when extracting a
# taxi's polygon from a mega-polygon.  Absorbs small polygon-
# boundary wobbles so the buffer reliably intersects the source.
TAXI_BUFFER_PAD_M = 2.0

# At a skeleton junction (degree ≥ 3), two branches are considered
# "same trunk continuation" when the angle between their tangent
# vectors is below this threshold.  60° catches curved taxis and
# mild bends while still separating cross-branches (which meet at
# ~90° or more).  Smaller values fragment long parallel taxis at
# every skeleton simplification kink; larger values erroneously
# merge cross-connectors into the trunk.
TRUNK_COLLINEAR_ANGLE_DEG = 60.0

# Skeleton node-identity rounding tolerance in meters.  Coordinates
# within this tolerance are treated as the same graph node.  1 m is
# safely below the shortest meaningful skeleton branch.
SKELETON_NODE_TOL_M = 1.0


# Role labels for classified shapes (see classify_shape_roles).
ROLE_PRIMARY_PARALLEL = "primary_parallel"
ROLE_STUB = "stub"
ROLE_SECONDARY_PARALLEL = "secondary_parallel"
ROLE_CROSS_CONNECTOR = "cross_connector"
ROLE_APRON = "apron"

# Role classification tunables.
PRIMARY_LENGTH_FACTOR = 0.25         # must be ≥ this × runway length
PRIMARY_MAX_DISTANCE_M = 300.0       # max distance from runway axis
BEARING_TOL_DEG = 15.0               # ± for "parallel" / "perpendicular"
STUB_MAX_LENGTH_M = 300.0            # stubs are short connectors

# Minimum 1D shared-boundary length for two shapes to count as
# adjacent.  Point-only contacts create ~1 m of inflation-induced
# false-positive "shared" length at the tolerance corners; 2 m is
# safely above that floor and below the shortest realistic pavement
# adjacency.
MIN_SHARED_LENGTH_M = 2.0

# Boundary-matching tolerance.  apt.dat polygons are drawn
# separately and their edges have small coordinate-precision gaps
# between neighbours; two boundaries within this distance count as
# "shared" for adjacency purposes.  0.5 m is a conservative value:
# large enough to absorb apt.dat simplification (~0.1 – 0.3 m
# drift) but small enough that truly separate pavement stays
# separate.
BOUNDARY_TOLERANCE_M = 0.5


# ──────────────────────────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class Shape:
    """One classified pavement polygon.

    Attributes:
        polygon: meter-space Polygon footprint.  The shape's
            boundary is what adjacent shapes share — no buffering,
            no pre-carving.
        kind: ``"taxi"`` or ``"apron"``.
        axis: for taxi, the 1D axial slope direction (MRR midline
            extended to polygon boundary); ``None`` for apron.
        width_m: taxi representative width (median of local widths
            along the axis, capped at ``WIDTH_PROBE_MAX_REACH_M``);
            ``0.0`` for apron.
    """
    polygon: Polygon
    kind: str
    axis: Optional[LineString]
    width_m: float


@dataclass(frozen=True)
class Adjacency:
    """One shared boundary between two Shape polygons.

    Attributes:
        shape_a: smaller index into the shapes tuple.
        shape_b: larger index (``shape_a < shape_b``).
        shared: longest LineString component of the 1D boundary
            intersection.  When the shared boundary is multi-part
            (e.g. a taxi touching two disjoint edges of a ring-
            shaped apron), ``shared`` is the longest component and
            ``length_m`` is the total length across all components.
        length_m: total shared boundary length in meters, summed
            across all 1D components.
    """
    shape_a: int
    shape_b: int
    shared: LineString
    length_m: float


# ──────────────────────────────────────────────────────────────────
# MRR helpers
# ──────────────────────────────────────────────────────────────────
def _mrr_midline_and_width(polygon: Polygon) -> Tuple[Optional[LineString],
                                                      float]:
    """Return ``(midline, short_side_m)`` for the min-area rotated
    rectangle of ``polygon``.  Midline connects the two short-side
    midpoints.  Returns ``(None, 0.0)`` for degenerate inputs.
    """
    if polygon is None or polygon.is_empty:
        return (None, 0.0)
    try:
        mrr = polygon.minimum_rotated_rectangle
    except _GEOM_EXC:
        return (None, 0.0)
    if mrr is None or mrr.is_empty or not hasattr(mrr, "exterior"):
        return (None, 0.0)
    coords = list(mrr.exterior.coords)
    if len(coords) < 5:
        return (None, 0.0)
    e01 = math.hypot(coords[1][0] - coords[0][0],
                     coords[1][1] - coords[0][1])
    e12 = math.hypot(coords[2][0] - coords[1][0],
                     coords[2][1] - coords[1][1])
    if e01 <= e12:
        m1 = ((coords[0][0] + coords[1][0]) * 0.5,
              (coords[0][1] + coords[1][1]) * 0.5)
        m2 = ((coords[2][0] + coords[3][0]) * 0.5,
              (coords[2][1] + coords[3][1]) * 0.5)
        short_s = e01
    else:
        m1 = ((coords[1][0] + coords[2][0]) * 0.5,
              (coords[1][1] + coords[2][1]) * 0.5)
        m2 = ((coords[3][0] + coords[0][0]) * 0.5,
              (coords[3][1] + coords[0][1]) * 0.5)
        short_s = e12
    return (LineString([m1, m2]), float(short_s))


def _mrr_aspect_short(polygon: Polygon) -> Tuple[float, float]:
    """Return ``(aspect, short_side_m)``.  Aspect is long/short;
    returns ``(0, 0)`` for degenerate inputs.
    """
    if polygon is None or polygon.is_empty:
        return (0.0, 0.0)
    try:
        mrr = polygon.minimum_rotated_rectangle
    except _GEOM_EXC:
        return (0.0, 0.0)
    if mrr is None or mrr.is_empty or not hasattr(mrr, "exterior"):
        return (0.0, 0.0)
    coords = list(mrr.exterior.coords)
    if len(coords) < 5:
        return (0.0, 0.0)
    e1 = math.hypot(coords[1][0] - coords[0][0],
                    coords[1][1] - coords[0][1])
    e2 = math.hypot(coords[2][0] - coords[1][0],
                    coords[2][1] - coords[1][1])
    long_s = max(e1, e2)
    short_s = min(e1, e2)
    aspect = (long_s / short_s) if short_s > 0 else 0.0
    return (aspect, short_s)


# ──────────────────────────────────────────────────────────────────
# Axis helpers
# ──────────────────────────────────────────────────────────────────
def _extend_axis_to_polygon_boundary(
    axis: LineString, polygon: Polygon, max_extension_m: float = 50.0,
) -> LineString:
    """Extend each endpoint of ``axis`` along its tangent to reach
    ``polygon``'s boundary.  MRR midlines stop at the short-side
    midpoints which are INSIDE the polygon; extending to the
    boundary ensures the axis terminals sit on the shared edge with
    a neighbour shape.  Failures fall back to the unextended axis.
    """
    if axis is None or axis.is_empty or polygon is None or polygon.is_empty:
        return axis
    cc = list(axis.coords)
    if len(cc) < 2:
        return axis

    def _extend(idx_end, idx_inner):
        ex, ey = cc[idx_end][0], cc[idx_end][1]
        ix, iy = cc[idx_inner][0], cc[idx_inner][1]
        dx, dy = ex - ix, ey - iy
        mag = math.hypot(dx, dy)
        if mag < 1e-9:
            return None
        dx /= mag
        dy /= mag
        target = (ex + dx * max_extension_m, ey + dy * max_extension_m)
        probe = LineString([(ex, ey), target])
        try:
            inside_seg = probe.intersection(polygon)
        except _GEOM_EXC:
            return None
        if inside_seg.is_empty:
            return None
        if inside_seg.geom_type == "LineString":
            seg = inside_seg
        elif hasattr(inside_seg, "geoms"):
            parts = [g for g in inside_seg.geoms
                     if g.geom_type == "LineString" and not g.is_empty]
            if not parts:
                return None
            parts.sort(key=lambda g: -g.length)
            seg = parts[0]
        else:
            return None
        scc = list(seg.coords)
        d0 = math.hypot(scc[0][0] - ex, scc[0][1] - ey)
        d1 = math.hypot(scc[-1][0] - ex, scc[-1][1] - ey)
        return scc[0] if d0 > d1 else scc[-1]

    new0 = _extend(0, 1)
    newN = _extend(-1, -2)
    out = list(cc)
    if new0 is not None:
        out[0] = new0
    if newN is not None:
        out[-1] = newN
    try:
        return LineString(out)
    except _GEOM_EXC:
        return axis


def _local_widths_along_axis(
    polygon: Polygon, axis: LineString, n_samples: int = 7,
) -> List[float]:
    """Return the perpendicular width of ``polygon`` sampled at
    ``n_samples`` equally-spaced interior points of ``axis``.  Rays
    that hit nothing inside ``WIDTH_PROBE_MAX_REACH_M`` are dropped.
    """
    if polygon is None or polygon.is_empty or axis is None:
        return []
    length = axis.length
    if length <= 0:
        return []
    widths: List[float] = []
    for i in range(n_samples):
        t = (i + 0.5) / n_samples
        c = axis.interpolate(t * length)
        t1 = axis.interpolate(max(0.0, t * length - 1.0))
        t2 = axis.interpolate(min(length, t * length + 1.0))
        tx, ty = t2.x - t1.x, t2.y - t1.y
        mag = math.hypot(tx, ty)
        if mag < 1e-9:
            continue
        tx /= mag
        ty /= mag
        hw = TS.local_half_width(
            polygon, c, (tx, ty),
            max_reach=WIDTH_PROBE_MAX_REACH_M,
        )
        if hw > 0:
            widths.append(2.0 * hw)
    return widths


def _median(xs: Sequence[float]) -> float:
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return 0.0
    if n % 2 == 1:
        return xs[n // 2]
    return 0.5 * (xs[n // 2 - 1] + xs[n // 2])


# ──────────────────────────────────────────────────────────────────
# Classification helpers
# ──────────────────────────────────────────────────────────────────
def _is_simple_strip(polygon: Polygon,
                    narrow_width_m: float,
                    simple_strip_aspect: float) -> bool:
    """Fast-path check: is the whole polygon a single narrow strip
    with a clear long axis?  True when MRR short ≤ threshold AND
    aspect ≥ 2.  Polygons that fail this go through the skeleton-
    based per-branch local-width analysis.
    """
    aspect, short_s = _mrr_aspect_short(polygon)
    return (short_s <= narrow_width_m
            and aspect >= simple_strip_aspect)


def _largest_polygon_part(geom) -> Optional[Polygon]:
    """From a Polygon or MultiPolygon, return the largest Polygon
    part; None for empty / non-polygonal inputs."""
    if geom is None or geom.is_empty:
        return None
    if geom.geom_type == "Polygon":
        return geom
    if hasattr(geom, "geoms"):
        parts = [g for g in geom.geoms
                 if g.geom_type == "Polygon" and not g.is_empty]
        if not parts:
            return None
        parts.sort(key=lambda p: -p.area)
        return parts[0]
    return None


def _validate_polygon(polygon: Polygon) -> Optional[Polygon]:
    """Return a valid Polygon.  Runs ``make_valid`` and falls back to
    a zero-width buffer for self-intersecting boundaries.  Returns
    ``None`` when no polygonal part survives cleanup."""
    if polygon is None or polygon.is_empty:
        return None
    if polygon.is_valid:
        return polygon
    try:
        fixed = make_valid(polygon)
    except _GEOM_EXC:
        try:
            fixed = polygon.buffer(0)
        except _GEOM_EXC:
            return None
    return _largest_polygon_part(fixed)


# ──────────────────────────────────────────────────────────────────
# Trunk extraction
# ──────────────────────────────────────────────────────────────────
def _node_key(pt: Tuple[float, float]) -> Tuple[int, int]:
    """Round a 2D coordinate to the skeleton-node grid.  Two
    coordinates with the same key are treated as the same node."""
    tol = SKELETON_NODE_TOL_M
    return (int(round(pt[0] / tol)), int(round(pt[1] / tol)))


def _unit_tangent(p_from: Tuple[float, float],
                  p_to: Tuple[float, float]) -> Tuple[float, float]:
    """Unit vector pointing from ``p_from`` to ``p_to``.  Returns
    ``(0, 0)`` for coincident points."""
    dx = p_to[0] - p_from[0]
    dy = p_to[1] - p_from[1]
    mag = math.hypot(dx, dy)
    if mag < 1e-9:
        return (0.0, 0.0)
    return (dx / mag, dy / mag)


def _extract_trunks(
    centerlines: Sequence[LineString],
    *,
    min_trunk_length_m: float,
    collinear_angle_deg: float = TRUNK_COLLINEAR_ANGLE_DEG,
    preferred_bearings: Optional[Sequence[float]] = None,
    bearing_tol_deg: float = 15.0,
) -> List[LineString]:
    """Merge skeleton LineStrings into trunks by iteratively
    joining the most-collinear edge pair at each junction.

    At a degree-3+ node (where three or more skeleton branches
    meet), Shapely's ``linemerge`` leaves all the branches as
    separate LineStrings because it can't pick the pass-through
    pair.  This routine walks every junction and, if two outgoing
    branches have tangents aligned within ``collinear_angle_deg``
    AND no other pair is more aligned, stitches them into one
    longer LineString.  Repeating until fixed-point reduces the
    degree of each junction by absorbing straight-through passes.

    Args:
        centerlines: linestrings from Voronoi skeleton.
        min_trunk_length_m: drop trunks shorter than this from the
            final result.
        collinear_angle_deg: maximum bend (in degrees) at which a
            pair of branches is still considered a straight-through
            continuation.
        preferred_bearings: optional list of "favoured" compass
            bearings (in degrees, 0 = +Y = north-ish, 90 = +X =
            east-ish).  At a junction where multiple pair choices
            are within ``collinear_angle_deg``, pairs whose merged
            direction aligns within ``bearing_tol_deg`` of any
            preferred bearing win over pairs that don't.  This is
            how we get V, L, A to extract as single trunks at
            cross-junctions with M, Q, R: the runway bearing wins.
        bearing_tol_deg: how close a pair's direction must be to a
            preferred bearing to count as "aligned".

    LineStrings shorter than ``min_trunk_length_m`` are filtered
    from the final result.
    """
    if not centerlines:
        return []

    # Active set of LineStrings, each with computed endpoints and
    # tangents.  We iteratively replace pairs with their merged
    # successor until no collinear pair remains at any junction.
    cos_threshold = math.cos(math.radians(collinear_angle_deg))

    class _Edge:
        __slots__ = ("coords", "key_a", "key_b", "tan_a", "tan_b",
                     "length")
        def __init__(self, coords):
            self.coords = coords
            self.key_a = _node_key(coords[0])
            self.key_b = _node_key(coords[-1])
            self.tan_a = _unit_tangent(coords[0], coords[1])
            self.tan_b = _unit_tangent(coords[-1], coords[-2])
            self.length = sum(math.hypot(coords[i+1][0] - coords[i][0],
                                         coords[i+1][1] - coords[i][1])
                              for i in range(len(coords) - 1))

    active: List[_Edge] = []
    for ls in centerlines:
        if ls is None or ls.is_empty:
            continue
        cc = list(ls.coords)
        if len(cc) < 2:
            continue
        if _node_key(cc[0]) == _node_key(cc[-1]):
            continue   # loop
        active.append(_Edge(cc))

    def _build_index(edges):
        """node_key -> list of (edge_idx, side) where side=0 means
        the key matches key_a of the edge, 1 means key_b."""
        idx: dict = {}
        for i, e in enumerate(edges):
            idx.setdefault(e.key_a, []).append((i, 0))
            idx.setdefault(e.key_b, []).append((i, 1))
        return idx

    # Normalise preferred bearings (compass degrees → unit vectors
    # along both directions since a taxi has no inherent forward).
    pref_vecs: List[Tuple[float, float]] = []
    if preferred_bearings:
        for b in preferred_bearings:
            rad = math.radians(b)
            vx, vy = math.sin(rad), math.cos(rad)
            pref_vecs.append((vx, vy))
    bearing_cos_tol = math.cos(math.radians(bearing_tol_deg))

    def _pair_aligns_preferred(ti, tj):
        """Does the pair (ti, tj) lie along any preferred bearing?
        A pair merges to a line; its direction is ``ti - tj``
        normalised (ti points AWAY from node, tj also; the merged
        line runs in direction ti).  A taxi has no "forward" so
        alignment is absolute-value dot product."""
        if not pref_vecs:
            return False
        # The merged line's direction, in the frame where the pair
        # is near-opposite, is well-approximated by ti.
        mag = math.hypot(ti[0], ti[1])
        if mag < 1e-9:
            return False
        dx, dy = ti[0] / mag, ti[1] / mag
        for vx, vy in pref_vecs:
            if abs(dx * vx + dy * vy) >= bearing_cos_tol:
                return True
        return False

    # Iteratively find the most-collinear pair at any junction and
    # merge them.  Loops bounded to guard against misbehaviour.
    for _ in range(10000):
        index = _build_index(active)
        # Pick the best pair at any junction, preferring preferred-
        # bearing-aligned pairs over mildly-straighter-but-off-axis
        # pairs when both are within collinear threshold.
        best_i = -1
        best_j = -1
        best_score = (0, cos_threshold)   # (aligned_bonus, straight)
        best_node = None
        best_side_i = 0
        best_side_j = 0
        for node_key, members in index.items():
            if len(members) < 2:
                continue
            out_tans = []
            for (ei, side) in members:
                e = active[ei]
                out_tans.append((ei, side,
                                 e.tan_a if side == 0 else e.tan_b))
            n_members = len(out_tans)
            for pi in range(n_members):
                ei, si, ti = out_tans[pi]
                for pj in range(pi + 1, n_members):
                    ej, sj, tj = out_tans[pj]
                    if ei == ej:
                        continue
                    cos_ang = ti[0] * tj[0] + ti[1] * tj[1]
                    straight = -cos_ang
                    if straight <= cos_threshold:
                        continue
                    aligned = (1 if _pair_aligns_preferred(ti, tj)
                               else 0)
                    # Prefer (aligned, straight) lexicographically:
                    # any aligned pair beats any non-aligned pair,
                    # ties broken by straightness.
                    score = (aligned, straight)
                    if score > best_score:
                        best_score = score
                        best_i = ei
                        best_j = ej
                        best_node = node_key
                        best_side_i = si
                        best_side_j = sj
        if best_i < 0:
            break  # no more merges possible

        # Merge edge best_i and best_j into one linestring through
        # best_node.
        ei = active[best_i]
        ej = active[best_j]
        ci = ei.coords
        cj = ej.coords
        # Orient so the merge happens at best_node: the merged
        # linestring runs from the far end of ei through best_node
        # to the far end of ej.
        if best_side_i == 0:   # best_node is at ci[0]
            ci_oriented = list(reversed(ci))   # ends AT best_node
        else:
            ci_oriented = list(ci)              # ends AT best_node
        if best_side_j == 0:   # best_node is at cj[0]
            cj_oriented = list(cj)              # starts AT best_node
        else:
            cj_oriented = list(reversed(cj))    # starts AT best_node
        merged_coords = ci_oriented + cj_oriented[1:]
        # Replace the two edges with one merged edge.
        new_edge = _Edge(merged_coords)
        # Remove the higher-indexed first to keep indices valid.
        hi, lo = max(best_i, best_j), min(best_i, best_j)
        active.pop(hi)
        active.pop(lo)
        active.append(new_edge)

    # Filter by minimum length and return as LineStrings.
    trunks: List[LineString] = []
    for e in active:
        if e.length < min_trunk_length_m:
            continue
        try:
            trunks.append(LineString(e.coords))
        except _GEOM_EXC:
            continue
    trunks.sort(key=lambda ls: -ls.length)
    return trunks


def _split_trunk_by_local_width(
    trunk: LineString,
    polygon: Polygon,
    narrow_width_m: float,
    step_m: float = 5.0,
    min_run_length_m: float = 50.0,
) -> List[Tuple[LineString, float]]:
    """Walk the trunk, sampling perpendicular polygon width at
    ``step_m`` intervals.  Group consecutive samples by whether
    local width is ≤ ``narrow_width_m``.  Return each contiguous
    narrow run as ``(sub_linestring, median_width)``.  Wide runs
    are dropped (they stay in the apron pool via the residual
    carve-out downstream).

    This is what makes SPJC's M taxi correctly split: the M1/M4
    stubs survive as narrow sub-trunks, while the wide portion
    between them falls into the apron residual.
    """
    length = trunk.length
    if length < min_run_length_m:
        return []
    n_samples = max(3, int(length / step_m) + 1)
    samples: List[Tuple[float, float]] = []
    for i in range(n_samples):
        t = min(length, i * step_m)
        pt = trunk.interpolate(t)
        t_near = trunk.interpolate(max(0.0, t - 1.0))
        t_far = trunk.interpolate(min(length, t + 1.0))
        tan = _unit_tangent((t_near.x, t_near.y), (t_far.x, t_far.y))
        if tan == (0.0, 0.0):
            continue
        hw = TS.local_half_width(
            polygon, pt, tan,
            max_reach=WIDTH_PROBE_MAX_REACH_M,
        )
        w = 2.0 * hw if hw > 0 else 0.0
        samples.append((t, w))

    if not samples:
        return []

    runs: List[Tuple[float, float, List[float]]] = []
    run_start: Optional[float] = None
    run_end = 0.0
    run_ws: List[float] = []
    for (t, w) in samples:
        is_narrow = 0 < w <= narrow_width_m
        if is_narrow:
            if run_start is None:
                run_start = t
            run_end = t
            run_ws.append(w)
        else:
            if run_start is not None:
                runs.append((run_start, run_end, run_ws))
                run_start = None
                run_ws = []
    if run_start is not None:
        runs.append((run_start, run_end, run_ws))

    # Merge adjacent narrow runs separated by a SHORT wide gap.
    # At a cross-junction (e.g. V meets Q) the polygon briefly
    # widens to include both taxis' footprints, creating a ~25 m
    # "wide" stretch along the trunk.  Splitting V here would
    # fragment it into pieces on either side of every intersection;
    # we instead treat wide gaps under ``narrow_width_m`` as still
    # belonging to the same trunk.
    merged_runs: List[Tuple[float, float, List[float]]] = []
    for run in runs:
        if not merged_runs:
            merged_runs.append(run)
            continue
        ps, pe, pw = merged_runs[-1]
        t0, t1, ws = run
        if (t0 - pe) <= narrow_width_m:
            merged_runs[-1] = (ps, t1, pw + ws)
        else:
            merged_runs.append(run)

    out: List[Tuple[LineString, float]] = []
    for (t0, t1, ws) in merged_runs:
        if t1 - t0 < min_run_length_m:
            continue
        try:
            sub = substring(trunk, t0, t1)
        except _GEOM_EXC:
            continue
        if sub is None or sub.is_empty or sub.length < min_run_length_m:
            continue
        out.append((sub, _median(ws)))
    return out


def _decompose_mega_polygon(
    polygon: Polygon,
    *,
    narrow_width_m: float,
    skeleton_min_branch_m: float,
    min_taxi_area_m2: float,
    preferred_bearings: Optional[Sequence[float]] = None,
) -> Tuple[List[Shape], List[Polygon]]:
    """Split a non-simple-strip polygon into narrow-arm taxi shapes
    plus apron residuals.

    Algorithm:
      1. Voronoi-skeletonise the polygon into branch centerlines.
      2. For each branch, measure median local perpendicular width
         via :func:`_local_widths_along_axis`.
      3. Keep only branches whose median local width is ≤
         ``narrow_width_m``.  These are the taxi candidates.
      4. Process candidates longest-first; each taxi's polygon is
         ``branch.buffer(width/2 + pad).intersection(source)`` MINUS
         whatever has already been claimed by a longer branch.
      5. The remaining area (source minus every claimed taxi
         polygon) is apron residual, possibly multi-part.

    Returns ``(taxi_shapes, apron_parts)``.  Empty lists when the
    polygon is too small or has no skeleton.
    """
    taxi_shapes: List[Shape] = []
    apron_parts: List[Polygon] = []

    try:
        # Skeleton with a permissive branch-length filter, so
        # trunk extraction can stitch short sub-branches together
        # across T-junctions.  Trunk extraction later enforces the
        # real minimum length on the stitched result.
        raw_centerlines = TS.extract_centerlines(
            polygon,
            min_path_length=min(skeleton_min_branch_m,
                                TS.DEFAULT_MIN_PATH_LENGTH_M) * 0.25,
            simplify_tol=5.0,  # keep junction geometry intact so
                               # trunk stitching finds the T-nodes
        )
    except _GEOM_EXC:
        raw_centerlines = []

    # Merge sub-branches into trunks along the straightest paths
    # through every junction.  This is what makes V, L, A come out
    # as single continuous linestrings at SPJC instead of being
    # split at every runway-stub T-junction.
    trunks = _extract_trunks(
        raw_centerlines,
        min_trunk_length_m=skeleton_min_branch_m,
        preferred_bearings=preferred_bearings,
    )

    # For each trunk: sample local perpendicular widths.  The
    # Voronoi centerline of a mega-polygon runs through local
    # widenings at every junction (fillets, stub-entry flares) so
    # the reported perpendicular width is systematically larger
    # than the bare taxi width.  The acceptance thresholds below
    # account for that:
    #
    # * MEDIAN ≤ narrow_width_m        → accept whole trunk as
    #                                     taxi.  At 45 m threshold
    #                                     this fits a real taxiway
    #                                     plus up to ~15 m of
    #                                     junction flare without
    #                                     needing a fudge factor.
    # * MEDIAN ≤ narrow_width_m * 1.5  → mixed: part narrow, part
    #                                     wide (e.g. M crossing an
    #                                     apron).  Split into
    #                                     narrow-only runs.
    # * MEDIAN > narrow_width_m * 1.5  → reject; predominantly wide
    #                                     (apron territory).
    accept_median_factor = 1.0
    mix_median_factor = 1.5
    candidates: List[Tuple[LineString, float]] = []
    for cl in trunks:
        if cl is None or cl.is_empty or cl.length < skeleton_min_branch_m:
            continue
        widths = _local_widths_along_axis(polygon, cl)
        if not widths:
            continue
        median_w = _median(widths)
        if median_w <= narrow_width_m * accept_median_factor:
            # Accept whole trunk — use median as representative
            # width (capped at narrow_width_m so downstream sizing
            # matches the actual pavement, not the flared Voronoi
            # reading).
            eff_w = min(median_w, narrow_width_m)
            candidates.append((cl, eff_w))
        elif median_w <= narrow_width_m * mix_median_factor:
            # Mixed; split into narrow-only runs.
            subs = _split_trunk_by_local_width(
                cl, polygon, narrow_width_m,
                min_run_length_m=skeleton_min_branch_m,
            )
            candidates.extend(subs)
        # else: trunk is predominantly wide → reject entirely.

    if not candidates:
        # No narrow arms → whole polygon goes to apron.
        apron_parts.append(polygon)
        return (taxi_shapes, apron_parts)

    # Process candidates longest-first; each claims its buffered
    # corridor from the remaining (yet-unclaimed) source polygon.
    candidates.sort(key=lambda c: -c[0].length)
    claimed = None
    for (cl, w) in candidates:
        half_w = max(w / 2.0, 5.0) + TAXI_BUFFER_PAD_M
        try:
            buf = cl.buffer(half_w, cap_style=3, join_style=2)
        except _GEOM_EXC:
            continue
        if buf.is_empty:
            continue
        # Clip to source polygon.
        try:
            region = buf.intersection(polygon)
        except _GEOM_EXC:
            continue
        if region.is_empty:
            continue
        # Subtract previously-claimed area.
        if claimed is not None and not claimed.is_empty:
            try:
                region = region.difference(claimed)
            except _GEOM_EXC:
                pass
        if region.is_empty:
            continue
        # Pick the largest connected component (buffering through
        # narrow passages can create disconnected bits).
        taxi_poly = _largest_polygon_part(region)
        if taxi_poly is None or taxi_poly.area < min_taxi_area_m2:
            continue
        taxi_poly = _validate_polygon(taxi_poly)
        if taxi_poly is None or taxi_poly.area < min_taxi_area_m2:
            continue
        axis = _extend_axis_to_polygon_boundary(cl, taxi_poly)
        taxi_shapes.append(Shape(
            polygon=taxi_poly, kind="taxi",
            axis=axis, width_m=w,
        ))
        if claimed is None:
            claimed = taxi_poly
        else:
            try:
                claimed = unary_union([claimed, taxi_poly])
            except _GEOM_EXC:
                pass

    # Apron residual = source minus all claimed taxi polygons.
    # Buffer the claimed area by a tiny amount before subtracting
    # so precision drift doesn't leave hairline slivers that then
    # count as taxi/apron overlaps.
    if claimed is not None and not claimed.is_empty:
        try:
            claimed_fat = claimed.buffer(0.1, join_style=2)
        except _GEOM_EXC:
            claimed_fat = claimed
        try:
            residual = polygon.difference(claimed_fat)
        except _GEOM_EXC:
            residual = polygon
        if not residual.is_empty:
            if residual.geom_type == "Polygon":
                residual = _validate_polygon(residual)
                if residual is not None:
                    apron_parts.append(residual)
            elif hasattr(residual, "geoms"):
                for g in residual.geoms:
                    if (g.geom_type == "Polygon"
                            and not g.is_empty):
                        valid_g = _validate_polygon(g)
                        if valid_g is not None:
                            apron_parts.append(valid_g)
    else:
        apron_parts.append(polygon)

    return (taxi_shapes, apron_parts)


# ──────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────
def decompose_pavement(
    taxi_polys: Sequence[Polygon],
    apron_polys: Sequence[Polygon],
    *,
    narrow_width_m: float = NARROW_WIDTH_M,
    simple_strip_aspect: float = SIMPLE_STRIP_ASPECT,
    min_apron_area_m2: float = MIN_APRON_AREA_M2,
    skeleton_min_branch_m: float = SKELETON_MIN_BRANCH_M,
    min_taxi_area_m2: float = MIN_TAXI_AREA_M2,
    preferred_bearings: Optional[Sequence[float]] = None,
) -> Tuple[Shape, ...]:
    """Classify apt.dat pavement polygons into taxi and apron Shapes.

    Classification is geometric, not semantic.  Caller's
    ``taxi_polys`` / ``apron_polys`` split is ignored — every
    polygon is reprocessed against:

    1. **Simple-strip fast path:** polygon's MRR short side ≤
       ``narrow_width_m`` AND aspect ≥ ``simple_strip_aspect`` →
       whole polygon is one taxi shape with MRR midline axis.
    2. **Mega-polygon path:** Voronoi-skeletonise, filter branches
       by median LOCAL width ≤ ``narrow_width_m`` and length ≥
       ``skeleton_min_branch_m``.  Each surviving branch claims a
       buffered corridor from the source polygon (longest-first,
       no overlap).  Remaining area is apron residual.

    Apron residuals (from all inputs) are then unioned into
    connected components — adjacent apt.dat apron polygons and
    mega-polygon residuals merge into one shape where they touch.

    Args:
        taxi_polys, apron_polys: meter-space polygons.  Combined
            and reprocessed; the caller's split is ignored.
        narrow_width_m: width cutoff for taxi classification.
        simple_strip_aspect: minimum aspect for the fast path.
        min_apron_area_m2: discard apron components smaller than
            this after unioning.
        skeleton_min_branch_m: minimum skeleton branch length to
            consider as a taxi candidate.
        min_taxi_area_m2: minimum taxi polygon area after buffer-
            clip-subtract.

    Returns:
        Tuple of :class:`Shape`, in emission order (all taxis
        first, then apron components).
    """
    taxi_valid = [p for p in taxi_polys
                  if p is not None and not p.is_empty]
    apron_valid = [p for p in apron_polys
                   if p is not None and not p.is_empty]

    if not taxi_valid and not apron_valid:
        return tuple()

    # Union each class SEPARATELY.  Unioning taxi-classified
    # polygons with apron-classified polygons dissolves the narrow
    # taxi boundary into the apron interior; perpendicular width
    # probes then travel into the apron and report apron-scale
    # widths, so narrow taxis inside mega-polygons get rejected by
    # the local-width filter.  Keeping the two pools separate
    # preserves the taxi-edge geometry the skeleton analysis needs.
    def _poly_components(union_geom):
        if union_geom is None or union_geom.is_empty:
            return []
        if union_geom.geom_type == "Polygon":
            return [union_geom]
        if hasattr(union_geom, "geoms"):
            return [g for g in union_geom.geoms
                    if g.geom_type == "Polygon" and not g.is_empty]
        return []

    try:
        taxi_union_geom = (unary_union(taxi_valid) if taxi_valid
                           else None)
    except _GEOM_EXC:
        taxi_union_geom = None
    try:
        apron_union_geom = (unary_union(apron_valid) if apron_valid
                            else None)
    except _GEOM_EXC:
        apron_union_geom = None

    taxi_components = _poly_components(taxi_union_geom)
    apron_components = _poly_components(apron_union_geom)

    taxi_shapes: List[Shape] = []
    apron_sources: List[Polygon] = list(apron_components)

    for poly in taxi_components:
        if _is_simple_strip(poly, narrow_width_m, simple_strip_aspect):
            midline, short_s = _mrr_midline_and_width(poly)
            if midline is None:
                apron_sources.append(poly)
                continue
            midline = _extend_axis_to_polygon_boundary(midline, poly)
            widths = _local_widths_along_axis(poly, midline)
            w = _median(widths) if widths else short_s
            taxi_shapes.append(Shape(
                polygon=poly, kind="taxi",
                axis=midline, width_m=w,
            ))
        else:
            ts, ar = _decompose_mega_polygon(
                poly,
                narrow_width_m=narrow_width_m,
                skeleton_min_branch_m=skeleton_min_branch_m,
                min_taxi_area_m2=min_taxi_area_m2,
                preferred_bearings=preferred_bearings,
            )
            taxi_shapes.extend(ts)
            apron_sources.extend(ar)

    # Extracted taxi polygons might overlap apron-classified polys
    # (apt.dat sometimes has cross-class overlap).  Subtract the
    # taxi union from the apron pool so every emitted shape is
    # interior-disjoint.
    apron_shapes: List[Shape] = []
    if apron_sources:
        if taxi_shapes:
            try:
                taxi_u = unary_union([t.polygon for t in taxi_shapes])
                taxi_u_fat = taxi_u.buffer(0.1, join_style=2)
                trimmed_sources = []
                for a in apron_sources:
                    try:
                        diff = a.difference(taxi_u_fat)
                    except _GEOM_EXC:
                        diff = a
                    if diff.is_empty:
                        continue
                    if diff.geom_type == "Polygon":
                        trimmed_sources.append(diff)
                    elif hasattr(diff, "geoms"):
                        for g in diff.geoms:
                            if (g.geom_type == "Polygon"
                                    and not g.is_empty):
                                trimmed_sources.append(g)
                apron_sources = trimmed_sources
            except _GEOM_EXC:
                pass
        try:
            apron_union = unary_union(apron_sources) \
                if apron_sources else None
        except _GEOM_EXC:
            apron_union = None
        if apron_union is not None and not apron_union.is_empty:
            parts = _poly_components(apron_union)
            for part in parts:
                if part.area < min_apron_area_m2:
                    continue
                valid = _validate_polygon(part)
                if valid is None or valid.area < min_apron_area_m2:
                    continue
                apron_shapes.append(Shape(
                    polygon=valid, kind="apron",
                    axis=None, width_m=0.0,
                ))

    return tuple(taxi_shapes + apron_shapes)


# ──────────────────────────────────────────────────────────────────
# Adjacency graph
# ──────────────────────────────────────────────────────────────────
def _boundary_line_parts(geom) -> List[LineString]:
    """Flatten a boundary-intersection geometry to its 1D
    LineString components only.  Point-only contacts are dropped.
    """
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "LineString":
        return [geom]
    if hasattr(geom, "geoms"):
        parts: List[LineString] = []
        for g in geom.geoms:
            parts.extend(_boundary_line_parts(g))
        return parts
    return []


def build_adjacency_graph(
    shapes: Sequence[Shape],
    *,
    min_shared_length_m: float = MIN_SHARED_LENGTH_M,
    tolerance_m: float = BOUNDARY_TOLERANCE_M,
) -> Tuple[Adjacency, ...]:
    """Return all pairwise Shape-to-Shape adjacencies.

    Two shapes are adjacent when the portion of shape A's boundary
    that lies within ``tolerance_m`` of shape B's polygon is at
    least ``min_shared_length_m`` long.  Point-only contacts and
    sub-meter slivers don't count.

    The tolerance absorbs apt.dat coordinate-precision gaps.  Real
    apt.dat polygons for touching pavement are typically 0.05 –
    0.3 m apart at their shared edge due to independent drawing /
    simplification; the default 0.5 m catches every realistic case
    without pulling in genuinely separate pavement.

    The returned tuple is canonicalised: each pair appears once,
    with ``shape_a < shape_b``.
    """
    out: List[Adjacency] = []
    n = len(shapes)
    if n < 2:
        return tuple()

    # Inflate each polygon by the FULL tolerance so shape A's
    # boundary is considered "adjacent to" shape B when it lies
    # within ``tolerance_m`` of polygon B.  The "shared boundary"
    # is the part of A's original boundary that falls inside B's
    # inflation — a LineString lying on A, which is what the solver
    # wants.  Using full tolerance (not half) matches the natural
    # reading of "boundaries within tolerance": gap + tolerance
    # drop is a single boundary-to-boundary test.
    inflated = [s.polygon.buffer(tolerance_m, join_style=2)
                for s in shapes]
    envelopes = [p.envelope for p in inflated]

    for i in range(n):
        boundary_i = shapes[i].polygon.boundary
        env_i = envelopes[i]
        for j in range(i + 1, n):
            if not env_i.intersects(envelopes[j]):
                continue
            try:
                inter = boundary_i.intersection(inflated[j])
            except _GEOM_EXC:
                continue
            parts = _boundary_line_parts(inter)
            if not parts:
                continue
            total_len = sum(p.length for p in parts)
            if total_len < min_shared_length_m:
                continue
            parts.sort(key=lambda p: -p.length)
            out.append(Adjacency(
                shape_a=i, shape_b=j,
                shared=parts[0],
                length_m=total_len,
            ))
    return tuple(out)


# ──────────────────────────────────────────────────────────────────
# Role classification
# ──────────────────────────────────────────────────────────────────
def _linestring_bearing_axis(ls: LineString) -> float:
    """Overall compass bearing of a LineString, mod 180 degrees
    (an axis has no inherent forward direction).
    0 ° = +Y (north), 90 ° = +X (east).
    """
    cc = list(ls.coords)
    if len(cc) < 2:
        return 0.0
    dx = cc[-1][0] - cc[0][0]
    dy = cc[-1][1] - cc[0][1]
    if math.hypot(dx, dy) < 1e-6:
        return 0.0
    return math.degrees(math.atan2(dx, dy)) % 180.0


def _axis_bearing_diff(b1: float, b2: float) -> float:
    """Smallest angular distance between two axis-bearings (mod
    180).  Returns a value in [0, 90]."""
    d = abs(b1 - b2) % 180.0
    return min(d, 180.0 - d)


def classify_shape_roles(
    shapes: Sequence[Shape],
    adjacencies: Sequence[Adjacency],
    runway_centerlines: Sequence[LineString],
    *,
    primary_length_factor: float = PRIMARY_LENGTH_FACTOR,
    primary_max_distance_m: float = PRIMARY_MAX_DISTANCE_M,
    bearing_tol_deg: float = BEARING_TOL_DEG,
    stub_max_length_m: float = STUB_MAX_LENGTH_M,
) -> List[str]:
    """Classify each shape into a role for the emission pipeline.

    Roles, in processing order:

    * ``primary_parallel`` — a taxi whose axis is parallel to a
      runway (within ``bearing_tol_deg``), within
      ``primary_max_distance_m`` of the runway, and at least
      ``primary_length_factor`` times the runway length.  These are
      emitted like runway segments (CIFP-anchored sloped rects) in
      the downstream pipeline.
    * ``stub`` — a short taxi (axis length ≤ ``stub_max_length_m``)
      adjacent to a primary_parallel.  Emitted as a single rect
      with a gap at each end.
    * ``cross_connector`` — a taxi whose axis is perpendicular
      (±``bearing_tol_deg`` of 90°) to a runway.  Single rect
      with gaps at each end.
    * ``secondary_parallel`` — any remaining taxi whose axis is
      parallel to a runway but too short or too far to be primary.
      Emitted as a DEM-driven segmented rect chain.
    * ``apron`` — all non-taxi shapes, plus any taxi that didn't
      match any of the above (unusual diagonal taxi arms).

    Returns a list ``roles[i]`` where ``i`` is the index into
    ``shapes``.
    """
    n = len(shapes)
    roles: List[str] = [ROLE_APRON] * n

    runway_bearings = [_linestring_bearing_axis(rw)
                       for rw in runway_centerlines]
    runway_lengths = [rw.length for rw in runway_centerlines]

    # Pass 1 — primary parallels.
    for i, s in enumerate(shapes):
        if s.kind != "taxi" or s.axis is None:
            continue
        axis_bearing = _linestring_bearing_axis(s.axis)
        axis_len = s.axis.length
        for ri, rw_line in enumerate(runway_centerlines):
            rwy_bearing = runway_bearings[ri]
            rwy_len = runway_lengths[ri]
            if _axis_bearing_diff(axis_bearing, rwy_bearing) > bearing_tol_deg:
                continue
            if axis_len < primary_length_factor * rwy_len:
                continue
            try:
                dist = s.polygon.distance(rw_line)
            except _GEOM_EXC:
                continue
            if dist > primary_max_distance_m:
                continue
            roles[i] = ROLE_PRIMARY_PARALLEL
            break

    # Pass 2 — cross-connectors: perpendicular-to-runway taxis.
    # Runs BEFORE secondary-parallel detection so that a perp
    # taxi never gets mis-tagged as parallel.
    for i, s in enumerate(shapes):
        if roles[i] != ROLE_APRON:
            continue
        if s.kind != "taxi" or s.axis is None:
            continue
        axis_bearing = _linestring_bearing_axis(s.axis)
        for rwy_bearing in runway_bearings:
            diff = _axis_bearing_diff(axis_bearing, rwy_bearing)
            if abs(diff - 90.0) <= bearing_tol_deg:
                roles[i] = ROLE_CROSS_CONNECTOR
                break

    # Pass 3 — secondary parallels: remaining parallel-to-runway
    # taxis that didn't qualify as primary.
    for i, s in enumerate(shapes):
        if roles[i] != ROLE_APRON:
            continue
        if s.kind != "taxi" or s.axis is None:
            continue
        axis_bearing = _linestring_bearing_axis(s.axis)
        for rwy_bearing in runway_bearings:
            if _axis_bearing_diff(axis_bearing, rwy_bearing) <= bearing_tol_deg:
                roles[i] = ROLE_SECONDARY_PARALLEL
                break

    # Pass 4 — stubs: short un-classified ("apron") taxis adjacent
    # to any parallel.  Doesn't overwrite cross_connector or
    # parallel classifications — those have their own semantics.
    # Additionally, a perpendicular taxi adjacent to TWO parallels
    # stays as cross_connector; a perpendicular taxi adjacent to
    # only ONE parallel is really a stub (short perpendicular
    # branch reaching from a parallel out to the airport boundary
    # or runway).
    parallel_set = {ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL}
    for i, s in enumerate(shapes):
        if s.kind != "taxi" or s.axis is None:
            continue
        if s.axis.length > stub_max_length_m:
            continue
        n_parallel_neighbors = 0
        for adj in adjacencies:
            if adj.shape_a == i and roles[adj.shape_b] in parallel_set:
                n_parallel_neighbors += 1
            elif adj.shape_b == i and roles[adj.shape_a] in parallel_set:
                n_parallel_neighbors += 1
        # Demote cross_connector to stub when it has EXACTLY one
        # parallel neighbour (a true between-parallels connector
        # has two; zero neighbours = isolated perpendicular that
        # stays a cross_connector).  Apron → stub when any parallel
        # touches.
        if roles[i] == ROLE_CROSS_CONNECTOR and n_parallel_neighbors == 1:
            roles[i] = ROLE_STUB
        elif roles[i] == ROLE_APRON and n_parallel_neighbors >= 1:
            roles[i] = ROLE_STUB

    return roles


def perimeter_coverage(
    shape_idx: int,
    shapes: Sequence[Shape],
    adjacencies: Sequence[Adjacency],
) -> Tuple[float, float]:
    """Return ``(covered_m, total_perimeter_m)`` for a given shape.

    ``covered_m`` is the sum of shared-boundary lengths the shape
    has with every neighbour; ``total_perimeter_m`` is the shape's
    full polygon boundary length.  The difference is the shape's
    "external" boundary — the portion that touches non-shape
    geometry (airport edge, runway surface, or untracked terrain).

    This is primarily a diagnostic tool: a shape with nearly-full
    coverage is well-embedded in the pavement graph; a shape with
    little coverage is either stand-alone or adjacent to geometry
    the graph doesn't model yet (runways, buildings, terrain).
    """
    if shape_idx < 0 or shape_idx >= len(shapes):
        return (0.0, 0.0)
    total = float(shapes[shape_idx].polygon.boundary.length)
    covered = 0.0
    for adj in adjacencies:
        if adj.shape_a == shape_idx or adj.shape_b == shape_idx:
            covered += adj.length_m
    return (covered, total)
