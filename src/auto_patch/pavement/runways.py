"""Runway rect emission, crossing resolution, and shoulder detection.

Owns the *static* runway geometry produced from apt.dat row-100
records: a 4-vertex rect including blast pads, plus the geometric
post-passes that operate on runway shapes after the initial emit:

* Runway-segment elevation sampling (interpolated from
  ``altitude_high`` / ``altitude_low``).
* Runway-crossing resolution (replace pairs of crossing runway
  segments with a single junction polygon).
* Runway-chain bridging (fill head-to-head gaps in a runway
  designation's segment chain so the surrounding apron-junction
  altitude doesn't dominate the runway's own CIFP profile).
* Runway-shoulder detection (absorb adjacent long thin pavement
  polygons into the runway's perpendicular extent).

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _runway_rect_m
    _sample_runway_segment_elev
    _resolve_runway_crossings
    _insert_runway_chain_bridges
    _detect_runway_shoulders
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from ..layout import (
    BuiltShape,
    PavementLayout,
    ROLE_CROSS_CONNECTOR,
    ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
)
from .vertices import _snap_polygon_vertices_to_rect_corners

# Narrow exception tuple for shapely / numeric-geometry failure
# modes.  Programming errors propagate so they surface immediately.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)


__all__ = [
    "_runway_rect_m",
    "_sample_runway_segment_elev",
    "_resolve_runway_crossings",
    "_insert_runway_chain_bridges",
    "_detect_runway_shoulders",
]


def _runway_rect_m(runway, to_m) -> Polygon:
    """4-vertex runway rect spanning end-to-end including blast pads.

    Targets are drawn with blast pads included (matches how the
    runway paint looks in satellite imagery).
    """
    ax, ay = to_m(runway.lon_a, runway.lat_a)
    bx, by = to_m(runway.lon_b, runway.lat_b)
    dx, dy = bx - ax, by - ay
    mag = math.hypot(dx, dy)
    if mag < 1e-9:
        return Polygon()
    ux, uy = dx / mag, dy / mag
    a_extra = runway.blast_a_m or 0.0
    b_extra = runway.blast_b_m or 0.0
    ax2 = ax - ux * a_extra
    ay2 = ay - uy * a_extra
    bx2 = bx + ux * b_extra
    by2 = by + uy * b_extra
    px, py = -uy, ux
    half = runway.width_m / 2.0
    return Polygon([
        (ax2 + px * half, ay2 + py * half),
        (bx2 + px * half, by2 + py * half),
        (bx2 - px * half, by2 - py * half),
        (ax2 - px * half, ay2 - py * half),
    ])


def _sample_runway_segment_elev(
        shape: "BuiltShape",
        x: float,
        y: float) -> "Optional[float]":
    """Linearly interpolate a runway segment's elevation at point
    ``(x, y)`` from its ``altitude_high`` / ``altitude_low``
    (sloped 4-corner rect) or ``altitude`` (flat).

    Returns None when the segment carries no elevation info.

    Convention used by :func:`_compute_elevations` segment emit
    (line ~2742): corners ``[0, 3]`` are the HIGH-elevation short
    edge, corners ``[1, 2]`` are the LOW-elevation short edge.
    Point's projection onto the high → low axis gives the
    interpolation parameter ``t`` (clipped to [0, 1] outside the
    segment's extent).
    """
    if shape.altitude is not None:
        return float(shape.altitude)
    # Per-vertex node_altitudes (e.g. a sloped runway segment that
    # was cut along a tile boundary — see ``tile_cut.py``).
    # Nearest-neighbour vertex sample.
    if shape.node_altitudes and shape.polygon is not None:
        try:
            coords = list(shape.polygon.exterior.coords)
        except _GEOM_EXC:
            coords = []
        n = min(len(coords), len(shape.node_altitudes))
        if n >= 1:
            best_d2 = float("inf")
            best_alt: Optional[float] = None
            for i in range(n):
                cx, cy = coords[i]
                d2 = (x - cx) ** 2 + (y - cy) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    best_alt = shape.node_altitudes[i]
            if best_alt is not None:
                return float(best_alt)
    if (shape.altitude_high is None
            or shape.altitude_low is None
            or shape.polygon is None):
        return None
    try:
        coords = list(shape.polygon.exterior.coords)
    except _GEOM_EXC:
        return None
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) < 4:
        return 0.5 * (float(shape.altitude_high)
                      + float(shape.altitude_low))
    high_mid_x = 0.5 * (coords[0][0] + coords[3][0])
    high_mid_y = 0.5 * (coords[0][1] + coords[3][1])
    low_mid_x = 0.5 * (coords[1][0] + coords[2][0])
    low_mid_y = 0.5 * (coords[1][1] + coords[2][1])
    ax = low_mid_x - high_mid_x
    ay = low_mid_y - high_mid_y
    L2 = ax * ax + ay * ay
    if L2 < 1e-6:
        return 0.5 * (float(shape.altitude_high)
                      + float(shape.altitude_low))
    t = ((x - high_mid_x) * ax + (y - high_mid_y) * ay) / L2
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    return (float(shape.altitude_high) + t
            * (float(shape.altitude_low)
               - float(shape.altitude_high)))


def _resolve_runway_crossings(
        layout: "PavementLayout",
        min_overlap_m2: float = 20.0,
        proximity_buffer_m: float = 2.0,
        ) -> int:
    """When two runway segments cross significantly, replace BOTH
    with a single junction polygon covering their union, with
    per-vertex altitudes interpolated from the source segments.

    Per user 2026-04-27: the existing overlap-clip pass clips one
    runway segment against another when they cross, producing a
    non-rectangular shape (5+ vertices) that still carries
    ``altitude_high`` / ``altitude_low`` tags — but X-Plane's patch
    format only renders 4-corner rects with H/L tagging, so the
    extra vertex breaks the slope rendering.  Crossings need
    multi-directional sloping which only a junction polygon
    (with per-vertex altitudes and triangulated rendering) can
    express.

    Detection uses an STRtree + union-find so transitively-
    overlapping segment groups are resolved together (e.g. CYXY's
    crosswind 02/20 crosses BOTH 14R/32L and 14L/32R; all three
    segment groups merge into one junction at the triple crossing).

    Returns the number of crossing groups resolved.
    """
    rwy_indices = [i for i, s in enumerate(layout.shapes)
                    if s.role == ROLE_RUNWAY
                    and s.polygon is not None
                    and not s.polygon.is_empty]
    if len(rwy_indices) < 2:
        return 0
    rwy_polys = [layout.shapes[i].polygon for i in rwy_indices]
    from shapely.strtree import STRtree
    try:
        tree = STRtree(rwy_polys)
    except _GEOM_EXC:
        return 0

    # Union-find for transitive grouping.
    n = len(rwy_indices)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union_uf(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for ai in range(n):
        pa = rwy_polys[ai]
        try:
            cands = tree.query(pa)
        except _GEOM_EXC:
            continue
        for ci in cands:
            bi = int(ci)
            if bi <= ai:
                continue
            pb = rwy_polys[bi]
            try:
                inter = pa.intersection(pb)
                if inter.is_empty or inter.area < min_overlap_m2:
                    continue
                union_uf(ai, bi)
            except _GEOM_EXC:
                continue

    # Group by root; ignore singleton groups.
    groups: Dict[int, List[int]] = {}
    for i in range(n):
        r = find(i)
        groups.setdefault(r, []).append(i)

    # Pre-compute sloping-rect corner snap targets for the corner-
    # alignment pass below.  Sloping rects = runway + parallel +
    # stub + cross-connector; junction vertices can only land on
    # their corners, never on their edges.
    sloping_roles_for_snap = (ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                              ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                              ROLE_CROSS_CONNECTOR)
    all_sloping_indices = [i for i, s in enumerate(layout.shapes)
                            if s.role in sloping_roles_for_snap
                            and s.polygon is not None
                            and not s.polygon.is_empty]

    drop_set: set = set()
    new_shapes: List[BuiltShape] = []
    n_resolved = 0
    for members in groups.values():
        if len(members) <= 1:
            continue
        seg_shapes = [layout.shapes[rwy_indices[m]] for m in members]
        seg_polys = [s.polygon for s in seg_shapes]
        member_set = {rwy_indices[m] for m in members}
        try:
            union_poly = unary_union(seg_polys)
            if not union_poly.is_valid:
                union_poly = union_poly.buffer(0)
        except _GEOM_EXC:
            continue
        if union_poly.is_empty:
            continue
        if union_poly.geom_type != "Polygon":
            polys = [g for g in getattr(union_poly, "geoms", [])
                      if g.geom_type == "Polygon"]
            polys.sort(key=lambda p: -p.area)
            if not polys:
                continue
            union_poly = polys[0]
        # Snap any union-polygon vertex that lands within 5 m of a
        # surviving sloping-rect corner to that corner.  Without
        # this, ``unary_union``'s intersection points can sit a few
        # metres along a surviving rect's edge — violating the
        # corner-only-junction-vertex invariant.
        other_sloping_polys = [
            layout.shapes[i].polygon
            for i in all_sloping_indices
            if i not in member_set]
        union_poly = _snap_polygon_vertices_to_rect_corners(
            union_poly, other_sloping_polys, snap_tol_m=5.0)
        try:
            coords = list(union_poly.exterior.coords)
        except _GEOM_EXC:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) < 3:
            continue

        # Per-vertex altitudes — INVERSE-DISTANCE-WEIGHTED average
        # across EVERY crossing segment.  Per user 2026-04-29
        # (CYXY -10075 ridge): the previous "buffer-containment"
        # selection picked a single nearest segment per vertex,
        # so two adjacent ring vertices on opposite sides of a
        # cross-runway gap could pick up different source-runway
        # altitudes (e.g. 14R/32L's 696.7 m vs 02/20's 693.7 m),
        # producing local steps as steep as 10 % between
        # neighbours and a +3.2 m ridge across runway 32L/16R.
        # Inverse-distance weighting blends every vertex's
        # altitude smoothly between all crossing runways — the
        # vertex that's right on a 14R/32L edge gets near-100 %
        # weight from that segment (d→0 ⇒ w→∞), while a vertex
        # equidistant between two runways gets a 50/50 average,
        # and the transition between the two regimes is smooth.
        ring_alts: List[Optional[float]] = []
        for (cx, cy) in coords:
            pt = Point(cx, cy)
            weighted_sum = 0.0
            weight_sum = 0.0
            for s in seg_shapes:
                e = _sample_runway_segment_elev(s, cx, cy)
                if e is None:
                    continue
                try:
                    d = s.polygon.distance(pt)
                except _GEOM_EXC:
                    continue
                # ε floor so a vertex exactly on a segment edge
                # still has a finite weight (just very large).
                d_eff = max(d, 0.1)
                w = 1.0 / (d_eff ** 2)
                weighted_sum += float(e) * w
                weight_sum += w
            if weight_sum > 0:
                ring_alts.append(
                    round(weighted_sum / weight_sum, 1))
            else:
                ring_alts.append(None)
        if any(a is None for a in ring_alts):
            continue
        # node_altitudes spans the closed ring.
        closed_alts: List[float] = list(ring_alts) + [ring_alts[0]]
        ref_combined = "+".join(
            s.ref for s in seg_shapes if s.ref)
        new_shape = BuiltShape(
            polygon=union_poly,
            role=ROLE_JUNCTION,
            ref=ref_combined,
            node_altitudes=closed_alts)
        new_shapes.append(new_shape)
        for m in members:
            drop_set.add(rwy_indices[m])
        n_resolved += 1

    if drop_set:
        layout.shapes = [s for i, s in enumerate(layout.shapes)
                          if i not in drop_set]
        layout.shapes.extend(new_shapes)
    return n_resolved


def _insert_runway_chain_bridges(
        layout: "PavementLayout",
        min_gap_m: float = 5.0,
        max_gap_m: float = 500.0,
        ) -> int:
    """For each runway designation, find pairs of surviving
    runway segments whose facing short edges are head-to-head
    with a gap > ``min_gap_m`` and < ``max_gap_m``, and insert a
    bridging RUNWAY shape between them with linearly
    interpolated altitudes.

    Per user 2026-04-30 (CYXY runway 32L/16R ridge): when one
    runway is absorbed into another via the runway-crossing
    resolver (or into a big apron via the apron-merge drop), the
    remaining runway segments leave a gap in the chain.  The
    surrounding apron-junction polygon's altitudes (DEM-derived,
    700+m at CYXY) then dominate the surface at the gap, sitting
    several metres above the runway's CIFP profile (693-697m at
    CYXY).  X-Plane renders this as a ridge crossing the runway.

    Inserting a bridging runway shape with altitudes anchored at
    the surviving segments' adjacent corners restores a runway-
    altitude plate at the gap.  The shape is emitted as
    ``ROLE_RUNWAY`` so downstream snap passes treat it the same
    way they treat a normal runway segment (corner-bucket
    altitude propagation, sloping rect H/L tagging, etc.).

    Detection algorithm (works at any airport):
      1. Group runway segments by their ``ref`` value (CIFP
         designation, e.g. "RW14R/RW32L").
      2. Within each group, sort segments along their shared
         axis direction.
      3. For each consecutive pair in axis order, measure the
         distance between the FIRST segment's "far" short edge
         and the SECOND segment's "near" short edge.  If
         ``min_gap_m < distance < max_gap_m`` AND the bridging
         rect's footprint is empty of any other runway segment
         (so we don't bridge across an active crossing), emit
         a bridge.

    Returns the number of bridge segments inserted.
    """
    rwy_shapes = [s for s in layout.shapes
                   if s.role == ROLE_RUNWAY
                   and s.polygon is not None
                   and not s.polygon.is_empty
                   and s.ref]
    if len(rwy_shapes) < 2:
        return 0
    # Group by ref.
    by_ref: Dict[str, List[BuiltShape]] = {}
    for s in rwy_shapes:
        by_ref.setdefault(s.ref, []).append(s)
    n_inserted = 0
    other_rwy_polys = [s.polygon for s in rwy_shapes]
    try:
        from shapely.strtree import STRtree as _STRtree
        rwy_tree = _STRtree(other_rwy_polys)
    except _GEOM_EXC:
        rwy_tree = None
    new_shapes: List[BuiltShape] = []
    for ref, segs in by_ref.items():
        if len(segs) < 2:
            continue
        # Determine the runway's axis direction from the FIRST
        # segment's altitude_high → altitude_low axis (corners
        # 0,3 = HIGH, 1,2 = LOW).
        s0 = segs[0]
        try:
            c0 = list(s0.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        if c0 and c0[0] == c0[-1]:
            c0 = c0[:-1]
        if len(c0) != 4:
            continue
        h_mid_0 = (0.5 * (c0[0][0] + c0[3][0]),
                   0.5 * (c0[0][1] + c0[3][1]))
        l_mid_0 = (0.5 * (c0[1][0] + c0[2][0]),
                   0.5 * (c0[1][1] + c0[2][1]))
        ax = l_mid_0[0] - h_mid_0[0]
        ay = l_mid_0[1] - h_mid_0[1]
        alen = math.hypot(ax, ay)
        if alen < 1.0:
            continue
        ux = ax / alen
        uy = ay / alen
        # Project each segment's centroid onto the axis to sort.
        def axis_pos(seg: BuiltShape) -> float:
            cx, cy = seg.polygon.centroid.x, seg.polygon.centroid.y
            return (cx - h_mid_0[0]) * ux + (cy - h_mid_0[1]) * uy
        segs_sorted = sorted(segs, key=axis_pos)
        # Walk consecutive pairs.
        for i in range(len(segs_sorted) - 1):
            a, b = segs_sorted[i], segs_sorted[i + 1]
            try:
                ca = list(a.polygon.exterior.coords)
                cb = list(b.polygon.exterior.coords)
            except _GEOM_EXC:
                continue
            if ca and ca[0] == ca[-1]:
                ca = ca[:-1]
            if cb and cb[0] == cb[-1]:
                cb = cb[:-1]
            if len(ca) != 4 or len(cb) != 4:
                continue
            # Find a's "far end" (the short edge facing b) and
            # b's "near end" (the short edge facing a).  Use
            # axis projection: a's far end has higher axis_pos,
            # b's near end has lower axis_pos.
            a_h_mid = (0.5 * (ca[0][0] + ca[3][0]),
                        0.5 * (ca[0][1] + ca[3][1]))
            a_l_mid = (0.5 * (ca[1][0] + ca[2][0]),
                        0.5 * (ca[1][1] + ca[2][1]))
            b_h_mid = (0.5 * (cb[0][0] + cb[3][0]),
                        0.5 * (cb[0][1] + cb[3][1]))
            b_l_mid = (0.5 * (cb[1][0] + cb[2][0]),
                        0.5 * (cb[1][1] + cb[2][1]))
            ah_pos = ((a_h_mid[0] - h_mid_0[0]) * ux
                       + (a_h_mid[1] - h_mid_0[1]) * uy)
            al_pos = ((a_l_mid[0] - h_mid_0[0]) * ux
                       + (a_l_mid[1] - h_mid_0[1]) * uy)
            bh_pos = ((b_h_mid[0] - h_mid_0[0]) * ux
                       + (b_h_mid[1] - h_mid_0[1]) * uy)
            bl_pos = ((b_l_mid[0] - h_mid_0[0]) * ux
                       + (b_l_mid[1] - h_mid_0[1]) * uy)
            # a's far short edge: whichever of (a_h_mid, a_l_mid)
            # has higher axis_pos.  Same for b's near edge: whichever
            # has lower axis_pos.
            if ah_pos > al_pos:
                a_far_corners = (ca[0], ca[3])
                a_far_alt = a.altitude_high
            else:
                a_far_corners = (ca[1], ca[2])
                a_far_alt = a.altitude_low
            if bh_pos < bl_pos:
                b_near_corners = (cb[0], cb[3])
                b_near_alt = b.altitude_high
            else:
                b_near_corners = (cb[1], cb[2])
                b_near_alt = b.altitude_low
            if (a_far_alt is None and a.altitude is not None):
                a_far_alt = a.altitude
            if (b_near_alt is None and b.altitude is not None):
                b_near_alt = b.altitude
            if a_far_alt is None or b_near_alt is None:
                continue
            # Distance between a's far midpoint and b's near
            # midpoint.
            a_far_mid = (0.5 * (a_far_corners[0][0]
                                  + a_far_corners[1][0]),
                          0.5 * (a_far_corners[0][1]
                                  + a_far_corners[1][1]))
            b_near_mid = (0.5 * (b_near_corners[0][0]
                                   + b_near_corners[1][0]),
                           0.5 * (b_near_corners[0][1]
                                   + b_near_corners[1][1]))
            gap = math.hypot(b_near_mid[0] - a_far_mid[0],
                             b_near_mid[1] - a_far_mid[1])
            if gap < min_gap_m or gap > max_gap_m:
                continue
            # Build bridge polygon.  Corners:
            #   0,3 = HIGH end (whichever of a/b has higher alt)
            #   1,2 = LOW end
            if a_far_alt >= b_near_alt:
                # a's far end is HIGH, b's near end is LOW.
                bridge_corners = [
                    a_far_corners[0],   # 0: A-side corner 1 (HIGH)
                    b_near_corners[0],  # 1: B-side corner 1 (LOW)
                    b_near_corners[1],  # 2: B-side corner 2 (LOW)
                    a_far_corners[1],   # 3: A-side corner 2 (HIGH)
                ]
                eh, el = float(a_far_alt), float(b_near_alt)
            else:
                bridge_corners = [
                    b_near_corners[0],
                    a_far_corners[0],
                    a_far_corners[1],
                    b_near_corners[1],
                ]
                eh, el = float(b_near_alt), float(a_far_alt)
            try:
                poly = Polygon(bridge_corners)
                if not poly.is_valid:
                    poly = poly.buffer(0)
                if poly.is_empty or poly.geom_type != "Polygon":
                    continue
            except _GEOM_EXC:
                continue
            # Reject if the bridge polygon overlaps any OTHER
            # runway segment (we'd be bridging across an active
            # crossing).
            overlap_with_other = False
            if rwy_tree is not None:
                try:
                    for hit in rwy_tree.query(poly):
                        idx = (int(hit) if hasattr(hit, "__int__")
                               else hit)
                        if not isinstance(idx, int):
                            continue
                        other = other_rwy_polys[idx]
                        if other is a.polygon or other is b.polygon:
                            continue
                        try:
                            inter = poly.intersection(other)
                            if (not inter.is_empty
                                    and inter.area > 1.0):
                                overlap_with_other = True
                                break
                        except _GEOM_EXC:
                            continue
                except _GEOM_EXC:
                    pass
            if overlap_with_other:
                continue
            shape = BuiltShape(
                polygon=poly, role=ROLE_RUNWAY, ref=ref)
            if abs(eh - el) >= 0.1:
                shape.altitude_high = round(eh, 1)
                shape.altitude_low = round(el, 1)
            else:
                shape.altitude = round(0.5 * (eh + el), 1)
            new_shapes.append(shape)
            n_inserted += 1
    if new_shapes:
        layout.shapes.extend(new_shapes)
    return n_inserted


def _detect_runway_shoulders(
        runway,
        to_m,
        pav_polys: "List[Polygon]",
        max_lat_gap_m: float = 1.0,
        min_self_inside_frac: float = 0.80,
        min_runway_overlap_frac: float = 0.40,
        min_axial_aspect: float = 2.5,
        min_strip_length_m: float = 50.0,
        ) -> "Tuple[float, float, List[int]]":
    """Scan ``pav_polys`` for polygons that are runway pavement
    (shoulders or the runway's own envelope polygon often labelled
    as a "taxiway" by apt.dat) and return the perpendicular extent
    those polygons + the runway itself jointly occupy.

    Per user 2026-04-27: long thin pavement polygons parallel to a
    runway and touching it (within 1 m of the runway edge) are
    SHOULDERS — fold them into the runway emit so the runway
    pavement reflects actual paved area, not just the apt.dat row-100
    designation.  At HECA, runway 05R/23L's apt.dat row-100 width is
    60 m but a row-110 polygon "New Taxiway 1" sits centered on the
    runway at 75.6 m wide (the runway-with-shoulders envelope).
    Without this absorption, the 7.8 m-per-side shoulder strip
    emits as junctions wrapping the runway.  At CYXY's short
    crosswind 02/20 (22.9 m wide), narrow strips "North of 02" and
    "South of 02" act as shoulders extending the runway pavement
    asymmetrically.

    Returns ``(new_left, new_right, absorbed_indices)``:
      * ``new_left``  — most-negative perpendicular offset from the
        runway centerline that emitted runway pavement should reach.
      * ``new_right`` — most-positive perpendicular offset.
      * ``absorbed_indices`` — indices into ``pav_polys`` of the
        polygons folded into the runway.  Caller should remove them
        from ``pav_polys`` (and ``apt_only_pav_polys`` if applicable)
        so they don't re-emit as separate junction shapes.

    When no shoulders are found, returns ``(-runway_half, +runway_half,
    [])``.

    Detection rules per polygon:
      1. Long axis aligned with runway (axial extent ≥
         ``min_axial_aspect`` × perpendicular extent, axial extent ≥
         ``min_strip_length_m``).
      2. Polygon mostly inside the runway's longitudinal extent:
         ``polygon's u-extent inside [0, L] ≥ min_self_inside_frac
         × polygon's total u-extent``.
      3. Polygon covers a significant fraction of the runway's
         length: ``u-overlap with [0, L] ≥ min_runway_overlap_frac
         × runway length L``.  Filters out tiny end-only blast pads
         (e.g. CYXY's '32L' polygon, 79 m at one end of a 2899 m
         runway) while admitting partial-length shoulders that
         cover roughly half the runway (e.g. CYXY's "North of 02"
         covering 47 % of the 02/20 runway).
      4. Perpendicular interval overlaps or touches
         ``[-runway_half, +runway_half]`` within ``max_lat_gap_m``.
      5. Side-extension limit: the polygon's extension PAST the
         runway edge on each side is less than the runway width.
         Filters out large aprons (e.g. CYXY's "Apron 1 and E",
         329 m wide, that nibbles the runway edge by 4 m) while
         admitting both wide envelopes (HECA's runway-with-
         shoulders polygon, 7.8 m past each edge) and narrow
         single-side shoulders (CYXY's "North of 02", 18 m past
         the north edge).
    """
    ax, ay = to_m(runway.lon_a, runway.lat_a)
    bx, by = to_m(runway.lon_b, runway.lat_b)
    udx = bx - ax
    udy = by - ay
    L = math.hypot(udx, udy)
    runway_half = runway.width_m / 2.0
    if L < 1.0:
        return (-runway_half, runway_half, [])
    ux, uy = udx / L, udy / L
    nx, ny = -uy, ux  # perpendicular (rotated 90° CCW from u)

    runway_width = 2.0 * runway_half
    new_left = -runway_half
    new_right = runway_half
    absorbed: List[int] = []
    for idx, pav in enumerate(pav_polys):
        if pav is None or pav.is_empty:
            continue
        try:
            coords = list(pav.exterior.coords)
        except _GEOM_EXC:
            continue
        if not coords:
            continue
        u_proj = [(cx - ax) * ux + (cy - ay) * uy for cx, cy in coords]
        n_proj = [(cx - ax) * nx + (cy - ay) * ny for cx, cy in coords]
        u_min, u_max = min(u_proj), max(u_proj)
        n_min, n_max = min(n_proj), max(n_proj)
        u_extent = u_max - u_min
        n_extent = n_max - n_min
        if u_extent < min_strip_length_m:
            continue
        if u_extent < min_axial_aspect * n_extent:
            continue
        # Polygon must be mostly inside the runway's longitudinal
        # extent.
        u_overlap = min(u_max, L) - max(u_min, 0.0)
        if u_overlap < min_self_inside_frac * u_extent:
            continue
        # Polygon must cover a significant fraction of the runway
        # length (filters out end-only blast pads / stopways).
        if u_overlap < min_runway_overlap_frac * L:
            continue
        # Perpendicular adjacency / overlap.
        if n_max < -runway_half - max_lat_gap_m:
            continue
        if n_min > runway_half + max_lat_gap_m:
            continue
        # Side-extension limit: the polygon's extension PAST the
        # runway edge on each side must be less than the runway
        # width.  Catches both:
        #   - HECA's "New Taxiway 1" envelope (8.7 m past each edge,
        #     well under 60 m runway width — accepted).
        #   - CYXY's "North of 02" shoulder (18 m past the north
        #     edge, under 22.9 m runway width — accepted).
        # Rejects:
        #   - CYXY's "Apron 1 and E" (extends 325 m past the east
        #     edge of 14R/32L, way over the 45.7 m runway width).
        north_ext = max(0.0, n_max - runway_half)
        south_ext = max(0.0, -runway_half - n_min)
        if north_ext > runway_width or south_ext > runway_width:
            continue
        absorbed.append(idx)
        if n_min < new_left:
            new_left = n_min
        if n_max > new_right:
            new_right = n_max
    return (new_left, new_right, absorbed)
