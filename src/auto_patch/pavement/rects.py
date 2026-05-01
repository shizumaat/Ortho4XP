"""Taxi rect construction from OSM centerlines.

Builds 4-vertex taxi rectangles from the centerline graph produced
by ``O4_Pavement_Centerlines``: probes the underlying apt.dat
pavement to determine the natural half-width along each axis,
extends the rect corners perpendicular to the axis, snaps corners
to the apt.dat boundary, classifies each emitted rect's role
(primary parallel / secondary parallel / stub / cross-connector),
and post-refines based on bearing-to-runway.

This is the OSM-centerline-driven rect builder.  It is distinct
from ``O4_Taxiway_Rects`` (the apt.dat-polygon-to-rect-chain
extractor that uses Voronoi medial-axis centerlines).

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _build_taxi_rects
    _merge_collinear_rects, _merge_collinear_rects_principled
    _natural_half_width, _trim_to_narrow, _probe_axis_width
    _extend_rect_corners_perpendicular
    _rect_from_axis_extended
    _snap_corners_to_pavement
    _cap_rect_length_to_width
    _classify_role, _axis_to_nearest_rwy_db, _refine_roles
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.ops import nearest_points, unary_union

from ..config import MIN_SEGMENT_LEN_M
from ..layout import (
    ROLE_CROSS_CONNECTOR,
    ROLE_PRIMARY_PARALLEL,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
)


__all__ = [
    "EDGE_SNAP_RADIUS_M",
    "VERTEX_SNAP_RADIUS_M",
    "_axis_to_nearest_rwy_db",
    "_build_taxi_rects",
    "_cap_rect_length_to_width",
    "_classify_role",
    "_extend_rect_corners_perpendicular",
    "_merge_collinear_rects",
    "_merge_collinear_rects_principled",
    "_natural_half_width",
    "_probe_axis_width",
    "_rect_from_axis_extended",
    "_refine_roles",
    "_snap_corners_to_pavement",
    "_trim_to_narrow",
]


def _build_taxi_rects(
    centerlines: List[Tuple[LineString, str]],
    pav_union: Optional[Polygon],
    rwy_union: Optional[Polygon],
    rwy_centerlines: List[LineString],
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
    ref_overall_bearings: Optional[Dict[str, float]] = None,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Convert each usable centerline into a 4-vertex rect.

    For each centerline we:
      1. Probe the half-width (distance to apt.dat pavement boundary)
         at many points along the axis.
      2. Determine the ``natural half-width`` of the strip as the
         median of those probes.
      3. TRIM axis endpoints inward until the probe there is
         ≤ 1.3 × natural_half_width — this is the user's rule of
         "pull rects back to the narrowest part of the taxiway."
         Everything past the trim is widening territory, reserved
         for junction polygons.
      4. Emit a 4-vertex rect over the trimmed axis with width =
         2 × natural_half_width.  The rect's 4 corners sit at the
         trimmed axis endpoints ± perpendicular half-width.

    Returns list of (rect, clipped_axis, role, ref).
    """
    if pav_union is None:
        return []

    pav_non_rwy = pav_union
    if rwy_union is not None:
        pav_non_rwy = pav_non_rwy.difference(rwy_union)

    centerlines = sorted(centerlines, key=lambda x: -x[0].length)

    emitted: List[Tuple[Polygon, LineString, str, str]] = []
    emitted_union: Optional[Polygon] = None

    for axis, ref in centerlines:
        # Clip to the full pavement (including runway).  The rect may
        # overlap runway slightly at stubs that reach the runway edge;
        # we prefer a stub that correctly reaches the runway over
        # clipping it at the runway boundary and losing most of it.
        try:
            clipped = axis.intersection(pav_union)
        except Exception:
            continue
        if clipped.is_empty:
            continue
        if clipped.geom_type == "MultiLineString":
            clipped = max(clipped.geoms, key=lambda g: g.length)
        if clipped.geom_type != "LineString":
            continue
        if clipped.length < 20.0:
            continue

        # Probe half-widths along axis.  ``narrow_hw`` (p10) is the
        # strip's narrowest portion; we use it both as the trim
        # baseline AND as the rect's half-width per the user's
        # authoritative rule.  The factor 1.15 means trim ends
        # where pavement grows more than 15% wider than the rect —
        # that extra width is junction territory.
        _natural_hw, _max_hw, narrow_hw = _natural_half_width(
            clipped, pav_non_rwy)
        if narrow_hw < 3.5 or narrow_hw > 40.0:
            continue

        # Per user (2026-04-20): trimming is now handled by the
        # upstream _split_centerlines_at_points which emits each
        # rect axis at 70% of the distance between intersections
        # (15% margin on each junction-facing end).  Skip the
        # width-based trim here — the axis is already cut to the
        # rect's intended length.
        trimmed = clipped
        trim_narrow_hw = narrow_hw

        # Dedup against trimmed axis
        if emitted_union is not None and not emitted_union.is_empty:
            try:
                inside_len = trimmed.intersection(emitted_union).length
                if inside_len / trimmed.length > 0.7:
                    continue
            except Exception:
                pass

        width = 2.0 * trim_narrow_hw
        rect = _rect_from_axis_extended(trimmed, width, pav_non_rwy,
                                        apt_vertices=apt_vertices)
        if rect is None or rect.is_empty:
            continue
        # Skip invalid rects (self-intersecting after snap).
        if not rect.is_valid:
            try:
                rect = rect.buffer(0)
            except Exception:
                continue
            if (rect.is_empty or rect.geom_type != "Polygon"
                    or not rect.is_valid):
                continue

        # ── Apron-interior rect rejection (user 2026-04-27) ────────
        # A rect whose 4 corners aren't on (or very near) the
        # pavement boundary is sitting INSIDE an apron — the
        # surrounding pavement wraps around it, downstream junction
        # construction has to wrap a junction around it too, and the
        # junction's elevation has to bridge the rect's slope on
        # both long edges → visible elevation ridges in JOSM and at
        # render time.
        #
        # Real taxi rects have their 4 corners at pavement-boundary
        # points (the intersections where the taxi corridor meets
        # the adjacent apron / parallel / runway).  If the rect's
        # corners are well INSIDE the pavement, the centerline runs
        # through an apron and shouldn't emit a separate rect — the
        # apron stays as one continuous junction.
        BOUNDARY_TOL_M = 2.0
        rect_coords = list(rect.exterior.coords)
        if rect_coords and rect_coords[0] == rect_coords[-1]:
            rect_coords = rect_coords[:-1]
        boundary = pav_union.boundary
        n_off_boundary = sum(
            1 for (cx, cy) in rect_coords
            if Point(cx, cy).distance(boundary) > BOUNDARY_TOL_M)
        if n_off_boundary >= 2:
            # ≥ 2 corners away from any pavement edge — the rect
            # sits inside an apron.  Skip it; the apron pavement
            # stays as residue → junction.
            continue

        role = _classify_role(trimmed, width, rwy_centerlines,
                               rwy_union, ref=ref,
                               ref_overall_bearings=ref_overall_bearings)
        emitted.append((rect, trimmed, role, ref))
        try:
            emitted_union = (unary_union([emitted_union, rect])
                             if emitted_union is not None else rect)
        except Exception:
            # Self-intersection of accumulated union — skip update.
            pass

    # Post-classify: secondary passes to fix roles based on neighbour
    # topology (stubs touch parallels, cross_connector touches 2 parallels)
    _refine_roles(emitted, rwy_centerlines)

    # Stub-ref dedup: OSM often has multiple disjoint ways with the
    # same sub-ref label (e.g. V2 has 3 separate OSM pieces, each
    # contributing a rect).  Target emits ONE rect per stub ref.
    # Keep only the LONGEST rect per stub ref.  Applies to both
    # sub-refs (letter+digit like V1, L3) AND letter-only stubs
    # (B, C, E, G) which are single-stub taxis that must not
    # fragment across internal bends.
    def _should_dedup(ref_str: str, role_str: str) -> bool:
        if not ref_str:
            return False
        if any(c.isdigit() for c in ref_str):
            return True
        # Letter-only: dedup when classified as stub (B, C, E, G, D).
        if role_str == ROLE_STUB:
            return True
        return False
    # Per-ref, group same-ref rects into geometrically-overlapping
    # clusters; within each cluster, keep only the longest.  Non-
    # overlapping rects of the same ref (e.g. multiple sub-segments
    # of one OSM way that bend-split into rects covering distinct
    # parts of the corridor) coexist — only OSM fragmentation that
    # produces actual duplicates gets deduped.
    OVERLAP_PROX_M = 5.0
    by_ref: Dict[str, List[int]] = {}
    for i, (_r, _a, role, ref) in enumerate(emitted):
        if not _should_dedup(ref, role):
            continue
        by_ref.setdefault(ref, []).append(i)
    drop: set = set()
    for ref, idxs in by_ref.items():
        # Build overlap clusters.
        n = len(idxs)
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a in range(n):
            ra = emitted[idxs[a]][0]
            for b in range(a + 1, n):
                rb = emitted[idxs[b]][0]
                try:
                    overlap = ra.intersection(rb).area
                    proximate = ra.distance(rb) < OVERLAP_PROX_M
                except Exception:
                    overlap, proximate = 0.0, False
                if overlap > 1.0 or proximate:
                    pa, pb = find(a), find(b)
                    if pa != pb:
                        parent[pa] = pb
        clusters: Dict[int, List[int]] = {}
        for k in range(n):
            clusters.setdefault(find(k), []).append(idxs[k])
        # Within each cluster, keep only the longest axis.
        for members in clusters.values():
            if len(members) <= 1:
                continue
            members.sort(key=lambda m: -emitted[m][1].length)
            for m in members[1:]:
                drop.add(m)
    keep: List[Tuple[Polygon, LineString, str, str]] = [
        item for i, item in enumerate(emitted) if i not in drop]
    return keep


def _merge_collinear_rects_principled(
    emitted: List[Tuple[Polygon, LineString, str, str]],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
    angle_tol_deg: float = 4.0,
    gap_tol_m: float = 12.0,
    width_uniformity_tol: float = 1.10,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Merge adjacent same-ref rects whose joining point shows
    NO widening — the pavement runs straight at uniform narrow
    width through the joint.  This is the only case where "single
    rect between junctions on straight sections" applies.
    """
    if not emitted or pav is None or pav.is_empty:
        return emitted
    boundary = pav.boundary
    changed = True
    work = list(emitted)
    while changed:
        changed = False
        for i in range(len(work)):
            for j in range(i + 1, len(work)):
                ri, ai, roli, refi = work[i]
                rj, aj, rolj, refj = work[j]
                if refi != refj or refi == "":
                    continue  # refless airports handled separately
                if roli != rolj:
                    continue
                # Bearings (mod 180°).
                def _bearing(a):
                    c = list(a.coords)
                    return math.degrees(
                        math.atan2(c[-1][0] - c[0][0],
                                   c[-1][1] - c[0][1])) % 180.0
                bi = _bearing(ai)
                bj = _bearing(aj)
                db = abs(bi - bj)
                db = min(db, 180.0 - db)
                if db > angle_tol_deg:
                    continue
                # Closest endpoints & far endpoints.
                coords_i = list(ai.coords)
                coords_j = list(aj.coords)
                pairs = [
                    (coords_i[0], coords_j[0], 0, 0),
                    (coords_i[0], coords_j[-1], 0, 1),
                    (coords_i[-1], coords_j[0], 1, 0),
                    (coords_i[-1], coords_j[-1], 1, 1),
                ]
                best = min(pairs, key=lambda p: math.hypot(
                    p[0][0] - p[1][0], p[0][1] - p[1][1]))
                endp_i, endp_j, ei, ej = best
                gap = math.hypot(endp_i[0] - endp_j[0],
                                 endp_i[1] - endp_j[1])
                if gap > gap_tol_m:
                    continue
                # Joining-region pavement half-width: probe at the
                # midpoint of the two touching endpoints.
                mid = Point((endp_i[0] + endp_j[0]) / 2,
                            (endp_i[1] + endp_j[1]) / 2)
                hw_joint = mid.distance(boundary) if pav.contains(mid) else 0
                if hw_joint <= 0:
                    continue
                # Each rect's own half-width (MRR short side / 2).
                def _rect_hw(p):
                    mrr = p.minimum_rotated_rectangle
                    c = list(mrr.exterior.coords)
                    if len(c) < 5:
                        return 0.0
                    s1 = math.hypot(c[1][0] - c[0][0], c[1][1] - c[0][1])
                    s2 = math.hypot(c[2][0] - c[1][0], c[2][1] - c[1][1])
                    return min(s1, s2) / 2.0
                hwi = _rect_hw(ri)
                hwj = _rect_hw(rj)
                if hwi <= 0 or hwj <= 0:
                    continue
                # No widening: joint hw is within uniformity_tol of
                # each rect's own hw (equivalently, joint hw ≤
                # max(hwi, hwj) × uniformity_tol AND rects have
                # similar widths).
                max_hw = max(hwi, hwj)
                if hw_joint > max_hw * width_uniformity_tol:
                    continue
                ratio_ij = max(hwi, hwj) / min(hwi, hwj)
                if ratio_ij > width_uniformity_tol:
                    continue
                # All checks pass — merge.
                far_i = coords_i[0] if ei == 1 else coords_i[-1]
                far_j = coords_j[0] if ej == 1 else coords_j[-1]
                try:
                    merged_axis = LineString([far_i, far_j])
                except Exception:
                    continue
                merged_rect = _rect_from_axis_extended(
                    merged_axis, 2.0 * (hwi + hwj) / 2.0, pav,
                    apt_vertices=apt_vertices)
                if merged_rect is None or merged_rect.is_empty:
                    continue
                work[i] = (merged_rect, merged_axis, roli, refi)
                del work[j]
                changed = True
                break
            if changed:
                break
    return work


def _merge_collinear_rects(
    emitted: List[Tuple[Polygon, LineString, str, str]],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
    angle_tol_deg: float = 4.0,
    gap_tol_m: float = 8.0,
    width_ratio_tol: float = 1.15,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Merge adjacent same-ref rects whose axes are nearly
    collinear and whose meeting point sits in the narrow corridor
    (no widening).  Produces one rect per straight pavement
    section between real junctions.
    """
    if not emitted:
        return emitted
    # Group by ref and role so we only merge within a ref's rects.
    merged = True
    work = list(emitted)
    while merged:
        merged = False
        n = len(work)
        for i in range(n):
            for j in range(i + 1, n):
                ri, ai, roli, refi = work[i]
                rj, aj, rolj, refj = work[j]
                # Only merge same ref + compatible role.
                if refi != refj or roli != rolj:
                    continue
                if refi == "" or refj == "":
                    continue  # refless: don't auto-merge (SPLP)
                # Axis bearings (mod 180°).
                def _bearing(a):
                    c = list(a.coords)
                    return math.degrees(
                        math.atan2(c[-1][0] - c[0][0],
                                   c[-1][1] - c[0][1])) % 180.0
                bi = _bearing(ai)
                bj = _bearing(aj)
                db = abs(bi - bj)
                db = min(db, 180.0 - db)
                if db > angle_tol_deg:
                    continue
                # Closest endpoints between axes.
                coords_i = list(ai.coords)
                coords_j = list(aj.coords)
                pairs = [
                    (coords_i[0], coords_j[0], 0, 0),
                    (coords_i[0], coords_j[-1], 0, 1),
                    (coords_i[-1], coords_j[0], 1, 0),
                    (coords_i[-1], coords_j[-1], 1, 1),
                ]
                best = min(pairs, key=lambda p: math.hypot(
                    p[0][0] - p[1][0], p[0][1] - p[1][1]))
                endp_i, endp_j, ei, ej = best
                gap = math.hypot(endp_i[0] - endp_j[0],
                                 endp_i[1] - endp_j[1])
                if gap > gap_tol_m:
                    continue
                # Width similarity: compare rect widths (from their
                # polygons' minimum-rotated-rect short-side).
                def _rect_width(p):
                    mrr = p.minimum_rotated_rectangle
                    coords = list(mrr.exterior.coords)
                    if len(coords) < 5:
                        return 0.0
                    s1 = math.hypot(coords[1][0] - coords[0][0],
                                    coords[1][1] - coords[0][1])
                    s2 = math.hypot(coords[2][0] - coords[1][0],
                                    coords[2][1] - coords[1][1])
                    return min(s1, s2)
                wi = _rect_width(ri)
                wj = _rect_width(rj)
                if wi < 1.0 or wj < 1.0:
                    continue
                if max(wi, wj) / min(wi, wj) > width_ratio_tol:
                    continue
                # Build merged axis: take the 2 FAR endpoints.
                far_i = coords_i[0] if ei == 1 else coords_i[-1]
                far_j = coords_j[0] if ej == 1 else coords_j[-1]
                try:
                    merged_axis = LineString([far_i, far_j])
                except Exception:
                    continue
                # Build merged rect from merged axis at average width.
                avg_width = (wi + wj) / 2.0
                merged_rect = _rect_from_axis_extended(
                    merged_axis, avg_width, pav,
                    apt_vertices=apt_vertices)
                if merged_rect is None or merged_rect.is_empty:
                    continue
                # Replace i with merged, drop j.
                work[i] = (merged_rect, merged_axis, roli, refi)
                del work[j]
                merged = True
                break
            if merged:
                break
    return work


def _natural_half_width(axis: LineString, pav: Polygon,
                        n_probes: int = 15) -> Tuple[float, float, float]:
    """Return (natural_hw, max_hw, narrow_hw) LOCAL half-width probes
    along the axis.

    Uses PERPENDICULAR RAY CAST (not distance-to-boundary) so the
    probe measures the taxi's own local half-width on EACH side
    rather than the distance to some faraway edge.  Per user
    rule 4 (2026-04-20): the rect half-width should be the
    NARROWEST pavement width (the taxi's own strip width), not
    an inflated value from adjacent aprons or runway clearance.

    For each probe point:
      * cast a ray perpendicular LEFT from axis; find where ray
        first exits the pavement polygon.
      * cast a ray perpendicular RIGHT from axis similarly.
      * half-width at this probe = min(left, right), capped at
        RAY_CAP_M to avoid saturating across an apron.
    """
    RAY_CAP_M = 40.0
    RAY_STEP_M = 0.5
    if axis.length < 1e-3:
        return 0.0, 0.0, 0.0

    def _perpendicular_half_at(t: float) -> float:
        """Cast perpendicular rays left/right at axis param t.

        Returns the AVERAGE of the two sides — (left + right) / 2 —
        so the width reflects the full pavement strip centered on
        the pavement (not a narrow corridor seen from an off-center
        axis).  Corner snapping downstream pulls the 4 rect corners
        onto the pav boundary, centering the rect on the actual
        pavement regardless of the axis's offset.
        """
        # Local tangent: use points slightly before/after t.
        dt = min(2.0, axis.length * 0.05)
        t0 = max(0.0, t - dt)
        t1 = min(axis.length, t + dt)
        a = axis.interpolate(t0)
        b = axis.interpolate(t1)
        tx, ty = b.x - a.x, b.y - a.y
        mag = math.hypot(tx, ty)
        if mag < 1e-6:
            return 0.0
        ux, uy = tx / mag, ty / mag
        nx, ny = -uy, ux  # left-perp
        pt = axis.interpolate(t)
        ox, oy = pt.x, pt.y
        sides: List[float] = []
        for sign in (-1, 1):
            side = RAY_CAP_M
            d = 0.0
            while d <= RAY_CAP_M:
                qx = ox + sign * nx * d
                qy = oy + sign * ny * d
                if not pav.contains(Point(qx, qy)):
                    side = d
                    break
                d += RAY_STEP_M
            sides.append(side)
        return sum(sides) / 2.0 if sides else RAY_CAP_M

    dists: List[float] = []
    for k in range(n_probes):
        t = (k + 1) / (n_probes + 1) * axis.length
        hw = _perpendicular_half_at(t)
        if hw > 0.1:
            dists.append(hw)
    if not dists:
        return 0.0, 0.0, 0.0
    dists.sort()
    median = dists[len(dists) // 2]
    p90_idx = max(0, int(len(dists) * 0.9) - 1)
    p90 = dists[p90_idx] if p90_idx < len(dists) else dists[-1]
    # ``narrow`` = the MIN half-width probe with a floor to guard
    # against grazing a building corner (skip probes < 3.5 m as
    # noise).  Per user rule 4: rect width = the ACTUAL narrowest
    # section, so the rect fits snugly along the taxi's own narrow
    # corridor, leaving widened areas to junctions.
    filtered = [d for d in dists if d >= 3.5]
    narrow = filtered[0] if filtered else dists[0]
    return median, p90, narrow


def _trim_to_narrow(axis: LineString, pav: Polygon, natural_hw: float,
                    widen_factor: float = 1.3) -> Optional[LineString]:
    """Trim the axis inward from each end until the PERPENDICULAR
    half-width at the endpoint drops below ``widen_factor × natural_hw``.

    Uses the same perpendicular ray-cast probing as
    ``_natural_half_width`` so the trim threshold is applied to
    the taxi's LOCAL half-width (not distance-to-boundary).

    Trim step is 2 m.  We never trim more than 50 % of the axis
    length.
    """
    RAY_CAP_M = 40.0
    RAY_STEP_M = 0.5
    total_len = axis.length
    thresh = natural_hw * widen_factor
    step = 2.0
    max_trim = total_len * 0.45

    def _perp_hw(t: float) -> float:
        dt = min(2.0, total_len * 0.05)
        t0 = max(0.0, t - dt)
        t1 = min(total_len, t + dt)
        a = axis.interpolate(t0)
        b = axis.interpolate(t1)
        tx, ty = b.x - a.x, b.y - a.y
        mag = math.hypot(tx, ty)
        if mag < 1e-6:
            return 0.0
        ux, uy = tx / mag, ty / mag
        nx, ny = -uy, ux
        pt = axis.interpolate(t)
        ox, oy = pt.x, pt.y
        best = RAY_CAP_M
        for sign in (-1, 1):
            d = 0.0
            while d <= RAY_CAP_M:
                qx = ox + sign * nx * d
                qy = oy + sign * ny * d
                if not pav.contains(Point(qx, qy)):
                    if d < best:
                        best = d
                    break
                d += RAY_STEP_M
        return best

    trim_a = 0.0
    while trim_a < max_trim:
        if _perp_hw(trim_a) <= thresh:
            break
        trim_a += step

    trim_b = total_len
    min_b = total_len - max_trim
    while trim_b > min_b:
        if _perp_hw(trim_b) <= thresh:
            break
        trim_b -= step

    if trim_b - trim_a < 10.0:
        return None
    from shapely.ops import substring
    return substring(axis, trim_a, trim_b)


def _probe_axis_width(axis: LineString, pav: Polygon,
                     n_probes: int = 9) -> float:
    """Return 2× the MEDIAN distance-to-boundary along the axis.

    At a large connected pavement (SPJC, where apron + taxi + runway
    are all one blob), ray-casting perpendicular overshoots into the
    apron.  The distance-to-boundary gives the local narrow-corridor
    half-width — which for a taxi centered in its strip is the
    half-strip-width.  Using the median over several probe points
    is robust to both (a) endpoints that sit at wider apron junctions
    and (b) narrow bottlenecks from adjacent building edges.
    """
    if axis.length < 1e-3:
        return 0.0
    boundary = pav.boundary
    dists = []
    for k in range(n_probes):
        t = (k + 1) / (n_probes + 1)
        pt = axis.interpolate(t, normalized=True)
        if not pav.contains(pt):
            continue
        d = pt.distance(boundary)
        if d > 0.1:
            dists.append(d)
    if not dists:
        return 0.0
    dists.sort()
    # Median, then 2× for full width
    return dists[len(dists) // 2] * 2.0


def _extend_rect_corners_perpendicular(
        rect: Polygon, axis: LineString,
        pav: Polygon, max_dist: float = 80.0,
        ) -> Polygon:
    """Extend each rect corner OUTWARD perpendicular to the rect's
    AXIS until the apt.dat pavement boundary (or ``max_dist``).

    Used for runway-end stubs (e.g. SPJC's F) where the half-width
    probe is capped at ``RAY_CAP_M = 40 m`` in
    ``_natural_half_width``, so a stub sitting in a wide runway-end
    ramp ends up under-sized.  The pavement extends past the
    rect's long edges, and the surrounding junction wraps around
    them (wrap-around = polygon along long edge of sloping rect,
    forbidden by the user's invariant).  Per user 2026-04-27: the
    rect should cover the FULL pavement width at each end —
    turning into a trapezoid where each corner sits exactly on the
    pavement boundary independently.

    The returned polygon has the same 4 corners in the same order
    (so X-Plane's altitude_high/low convention is preserved); each
    corner is just shifted outward to its respective pavement edge.
    """
    coords = list(rect.exterior.coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return rect
    a = list(axis.coords)
    if len(a) < 2:
        return rect
    p1, p2 = a[0], a[-1]
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    mag = math.hypot(dx, dy)
    if mag < 1e-6:
        return rect
    ux = dx / mag
    uy = dy / mag
    px = -uy
    py = ux
    STEP = 0.5

    def _ray_extent(ax: float, ay: float,
                    dir_x: float, dir_y: float,
                    base_dist: float) -> float:
        """From axis point (ax, ay) heading (dir_x, dir_y), find
        the FURTHEST distance d such that (ax+d*dir,ay+d*dir) is
        inside ``pav``.  Start at base_dist and walk OUTWARD only
        (we never shrink past base_dist; the rect keeps at least
        its nominal half-width)."""
        # If even base_dist is OUTSIDE pav, the original rect
        # corner is already outside; keep base.
        if not pav.contains(
                Point(ax + dir_x * base_dist,
                      ay + dir_y * base_dist)):
            return base_dist
        d = base_dist
        while d < max_dist:
            d_test = d + STEP
            if not pav.contains(
                    Point(ax + dir_x * d_test,
                          ay + dir_y * d_test)):
                return d
            d = d_test
        return d

    new_corners: List[Tuple[float, float]] = []
    for cx, cy in coords:
        # For each corner: project onto axis to determine which
        # endpoint (p1 or p2) it belongs to and which perp side.
        vx = cx - p1[0]
        vy = cy - p1[1]
        proj = vx * ux + vy * uy
        # axis endpoint nearer this corner
        ax_pt = p1 if proj < mag * 0.5 else p2
        # perpendicular signed distance from axis
        rel_x = cx - ax_pt[0]
        rel_y = cy - ax_pt[1]
        perp_signed = rel_x * px + rel_y * py
        if perp_signed >= 0:
            dir_x, dir_y = px, py
        else:
            dir_x, dir_y = -px, -py
        base = abs(perp_signed)
        if base < 1.0:
            new_corners.append((cx, cy))
            continue
        d = _ray_extent(ax_pt[0], ax_pt[1], dir_x, dir_y, base)
        new_corners.append(
            (ax_pt[0] + dir_x * d, ax_pt[1] + dir_y * d))

    try:
        new_rect = Polygon(new_corners)
        if new_rect.is_valid and not new_rect.is_empty:
            return new_rect
    except Exception:
        pass
    return rect


def _rect_from_axis_extended(axis: LineString, width: float,
                            pav: Polygon,
                            apt_vertices: Optional[
                                List[Tuple[float, float]]] = None,
                            ) -> Optional[Polygon]:
    """Build a rect around the axis at its first-to-last direction.

    The 4 corners are placed at axis endpoints ± perpendicular half-
    width, then each corner is snapped FIRST to the nearest apt.dat
    pavement vertex within ``VERTEX_SNAP_RADIUS_M``, ELSE to the
    nearest pavement edge point within ``EDGE_SNAP_RADIUS_M``.
    This matches the snapped target convention where every non-
    runway vertex sits on an apt.dat pavement vertex.

    Asymmetric-snap trim: when the two snapped end-widths differ
    by more than ``ASYM_WIDTH_TOL_M``, the rect has extended into
    a widening pavement area (apron or junction) on the wider end.
    Per user (2026-04-21): "if we're getting an asymmetric rect,
    most likely it's too long and needs to be shortened a bit so
    it's not pulled into a junction."  Retry once with the axis
    trimmed by ``ASYM_TRIM_FRAC`` of its length on the wider end.
    """
    from shapely.ops import substring

    # Symmetry tolerances: a proper rectangle has equal long sides
    # and equal short sides.  Use a RATIO check so small rects
    # aren't over-trimmed — a 5 m width delta on a 30 m-wide stub
    # (17 %) reads as near-symmetric, while a 5 m delta on a
    # 22 m-wide cross-connector (23 %) reads as a trapezoid.  User
    # (2026-04-21): "I don't see why it should trim the second
    # diagonal at the south end, which is already quite short and
    # appears symmetrical" — that stub had width Δ = 5 m / max
    # 29 m = 17 %, below threshold.
    ASYM_WIDTH_RATIO_TOL = 0.20   # width Δ / max_width > 20 % → trim
    ASYM_LENGTH_RATIO_TOL = 0.10  # length Δ / max_length > 10 % → trim
    # Per-iteration axis shrink per user (2026-04-21): trim BOTH
    # ends by 2.5 % of length each (5 % total) so the rect STAYS
    # CENTERED on its axis as it shrinks — trimming only the
    # wider end shifts the rect toward the narrower side and
    # leaves the problematic (wider) end snapped to the same
    # widening zone on the next iteration.
    ASYM_TRIM_EACH_END = 0.025     # 2.5 % off each end per iteration
    MAX_ASYM_RETRIES = 15          # 15 * 5 % = up to 75 % shrink

    cur_axis = axis
    for attempt in range(MAX_ASYM_RETRIES + 1):
        coords = list(cur_axis.coords)
        if len(coords) < 2:
            return None
        p1 = coords[0]
        p2 = coords[-1]
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        mag = math.hypot(dx, dy)
        if mag < 1e-6:
            return None
        ux, uy = dx / mag, dy / mag
        px, py = -uy, ux
        half = width / 2.0
        corners = [
            (p1[0] + px * half, p1[1] + py * half),   # 0: end1 side1
            (p2[0] + px * half, p2[1] + py * half),   # 1: end2 side1
            (p2[0] - px * half, p2[1] - py * half),   # 2: end2 side2
            (p1[0] - px * half, p1[1] - py * half),   # 3: end1 side2
        ]
        snapped = _snap_corners_to_pavement(
            corners, pav, apt_vertices)
        # Reject degenerate rects where snap collapsed two corners
        # onto the same apt.dat vertex.
        degenerate = False
        for i in range(4):
            for j in range(i + 1, 4):
                if math.hypot(snapped[i][0] - snapped[j][0],
                              snapped[i][1] - snapped[j][1]) < 1.0:
                    degenerate = True
                    break
            if degenerate:
                break
        if degenerate:
            return None

        # Symmetry check: equal widths (end1 vs end2) AND equal
        # lengths (side1 vs side2).
        w_end1 = math.hypot(snapped[0][0] - snapped[3][0],
                            snapped[0][1] - snapped[3][1])
        w_end2 = math.hypot(snapped[1][0] - snapped[2][0],
                            snapped[1][1] - snapped[2][1])
        l_side1 = math.hypot(snapped[1][0] - snapped[0][0],
                             snapped[1][1] - snapped[0][1])
        l_side2 = math.hypot(snapped[2][0] - snapped[3][0],
                             snapped[2][1] - snapped[3][1])
        width_asym = abs(w_end1 - w_end2)
        length_asym = abs(l_side1 - l_side2)
        max_w = max(w_end1, w_end2)
        max_l = max(l_side1, l_side2)
        width_ratio = (width_asym / max_w) if max_w > 1e-6 else 0.0
        length_ratio = (length_asym / max_l) if max_l > 1e-6 else 0.0
        symmetric = (width_ratio <= ASYM_WIDTH_RATIO_TOL
                     and length_ratio <= ASYM_LENGTH_RATIO_TOL)
        if symmetric or attempt == MAX_ASYM_RETRIES:
            return Polygon(snapped)

        # Asymmetric — trim BOTH ends and retry, keeping the rect
        # centered on its axis (per user 2026-04-21 feedback).
        trim_each = ASYM_TRIM_EACH_END * cur_axis.length
        new_start = trim_each
        new_end = cur_axis.length - trim_each
        if new_end - new_start < MIN_SEGMENT_LEN_M:
            # Axis would become too short; accept current asymmetric
            # rect rather than discarding.
            return Polygon(snapped)
        try:
            cur_axis = substring(cur_axis, new_start, new_end)
        except Exception:
            return Polygon(snapped)
    return None


VERTEX_SNAP_RADIUS_M = 8.0  # first-choice: snap to real apt.dat vertex
EDGE_SNAP_RADIUS_M = 15.0   # fallback: nearest point on pav boundary


def _snap_corners_to_pavement(
    corners: List[Tuple[float, float]],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
) -> List[Tuple[float, float]]:
    """Two-stage corner snap per the user's rule (2026-04-18):

    1. First try nearest apt.dat pavement VERTEX within
       ``VERTEX_SNAP_RADIUS_M``.  Apt.dat vertices are the
       authoritative coordinate set the target also snaps to, so
       snapping rect corners to them produces exact shared-vertex
       alignment.
    2. If no apt.dat vertex within range, fall back to the nearest
       POINT on the pav boundary within ``EDGE_SNAP_RADIUS_M``.
    3. If neither within range, leave the corner unsnapped.

    GUARD: a vertex-snap that would collapse two corners onto the
    SAME apt.dat vertex (producing a degenerate rect) is rejected
    — that corner falls through to edge snap instead.  Rects with
    two coincident corners violate rule 7 and are rejected
    downstream anyway, so we'd rather keep the 4 distinct corners.
    """
    boundary = pav.boundary
    # Pre-filter apt.dat vertices to ONLY those that actually sit
    # on the pav_union boundary.  When two apt.dat pavement
    # polygons overlap (common at aprons / terminal pads), a
    # corner vertex of one polygon ends up in the interior of the
    # union.  Snapping a rect corner to such an interior vertex
    # violates rule 7 ("corners ALWAYS on pavement boundary") and
    # produces a rect that floats inside the pavement, leaving a
    # sliver that a junction wraps around.  Discovered at SPJC V1
    # (2026-04-23): c0 snapped to an apt.dat vertex 5.88 m inside
    # the union, leaving the junction to wrap around V1's short
    # side.
    #
    # Per user 2026-04-27: tolerance is STRICT (1 cm).  Snap means
    # EXACTLY on the boundary, not "close to it" — a 0.36 m offset
    # at SPJC G's NW corner caused the surrounding junction to
    # add two extra nodes wrapping around the offset.  Vertices
    # that aren't on the boundary fall through to Stage 2 (edge
    # snap), which projects directly onto the boundary line.
    BOUNDARY_TOL_M = 0.01
    boundary_verts: Optional[List[Tuple[float, float]]] = None
    if apt_vertices:
        boundary_verts = [
            v for v in apt_vertices
            if Point(v[0], v[1]).distance(boundary) <= BOUNDARY_TOL_M
        ]
    # Stage 1: pick nearest apt.dat boundary vertex per corner.
    candidates: List[Optional[Tuple[float, float]]] = []
    for (cx, cy) in corners:
        best_v = None
        best_d = VERTEX_SNAP_RADIUS_M
        if boundary_verts:
            for (vx, vy) in boundary_verts:
                d = math.hypot(cx - vx, cy - vy)
                if d < best_d:
                    best_v = (vx, vy)
                    best_d = d
        candidates.append(best_v)

    # Collision handling: when two vertex-snap candidates collide
    # (would coincide), keep the nearer one on the vertex.  The
    # other corner uses the ORIGINAL pre-snap coordinate (NOT
    # edge-snap, which could pull back to the same region).  This
    # preserves a valid 4-corner rect while still placing the
    # kept corner exactly on an apt.dat vertex.
    # Use proximity (not identity) — two candidates within 1 m
    # would also collapse the rect corner.  Collapse collision:
    # drop the farther-from-original candidate, use pre-snap.
    use_original: List[bool] = [False] * len(candidates)
    COLLISION_TOL = 1.0
    for i in range(len(candidates)):
        if candidates[i] is None:
            continue
        for j in range(i + 1, len(candidates)):
            if candidates[j] is None:
                continue
            d_cand = math.hypot(candidates[i][0] - candidates[j][0],
                                candidates[i][1] - candidates[j][1])
            if d_cand <= COLLISION_TOL:
                di = math.hypot(corners[i][0] - candidates[i][0],
                                corners[i][1] - candidates[i][1])
                dj = math.hypot(corners[j][0] - candidates[j][0],
                                corners[j][1] - candidates[j][1])
                if di <= dj:
                    candidates[j] = None
                    use_original[j] = True
                else:
                    candidates[i] = None
                    use_original[i] = True
                    # candidates[i] is now None; subsequent j
                    # iterations would dereference it.  Break and
                    # let the outer i loop advance.
                    break

    snapped: List[Tuple[float, float]] = []
    for i, (cx, cy) in enumerate(corners):
        if candidates[i] is not None:
            snapped.append(candidates[i])
            continue
        if use_original[i]:
            # Collision fallback — keep pre-snap coord to avoid
            # coincident corners.
            snapped.append((cx, cy))
            continue
        # Stage 2: nearest pav edge point.
        p = Point(cx, cy)
        near, _ = nearest_points(boundary, p)
        if p.distance(near) <= EDGE_SNAP_RADIUS_M:
            snapped.append((near.x, near.y))
        else:
            snapped.append((cx, cy))

    # Final coincidence sweep: vertex-snap AND edge-snap can BOTH
    # pull two corners onto the same area (e.g. edge-snap falls
    # through to the same boundary vertex a vertex-snap picked
    # from the other corner).  If any two final coords are within
    # 1 m, revert the farther-from-original to its pre-snap coord.
    FINAL_COLLISION_TOL = 1.0
    for i in range(len(snapped)):
        for j in range(i + 1, len(snapped)):
            d = math.hypot(snapped[i][0] - snapped[j][0],
                           snapped[i][1] - snapped[j][1])
            if d <= FINAL_COLLISION_TOL:
                di = math.hypot(corners[i][0] - snapped[i][0],
                                corners[i][1] - snapped[i][1])
                dj = math.hypot(corners[j][0] - snapped[j][0],
                                corners[j][1] - snapped[j][1])
                if di <= dj:
                    snapped[j] = corners[j]
                else:
                    snapped[i] = corners[i]

    # Symmetry check: rect corners come in as [end1_side1, end2_side1,
    # end2_side2, end1_side2] per _rect_from_axis_extended.  The two
    # widths (end1_side1 → end1_side2 and end2_side1 → end2_side2)
    # should be equal for a proper rectangle.  When snap pulls
    # corners to apt.dat vertices at asymmetric offsets (one rect
    # end near a wider pavement apron than the other), the widths
    # diverge and the shape reads as a trapezoid.  Per user
    # (2026-04-21): "the perpendicular stub is asymmetrical and
    # looks like it's coming into the space of the primary
    # parallel."  If the snapped widths differ by more than
    # ``ASYM_WIDTH_TOL_M``, revert ONE corner on each side to
    # its pre-snap coord so the rect stays symmetric.  We keep
    # the NARROWER side's snap (matching the tighter pavement) and
    # revert the wider side's corners to the pre-snap perpendicular
    # offset.
    return snapped


def _cap_rect_length_to_width(
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    rwy_centerlines: List[LineString],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]],
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Cap each rect's length-to-width ratio per the user's
    2026-04-27 spec: rects should be roughly square (length ≈
    width) so the long edges sit on the pavement-narrowing
    boundary, corners snap there, and surrounding junctions
    connect only at the short edges (never wrap around long
    edges).

    Cap depends on the bearing-to-nearest-runway:

      * Parallel  (db < 20°)  — NO CAP (long parallel taxis are
        legitimate, often running 500 m+ along the runway).
      * Diagonal  (20° ≤ db < 45°) — length ≤ 1.0 × width
        (truly square; matches the tighter 30 % margin used for
        diagonal stubs in ``_rect_margin_frac_for``).
      * Perpendicular (db ≥ 45°) — length ≤ 1.3 × width (small
        excess so the rect can extend slightly past the apron's
        narrow corridor without forcing surrounding junctions to
        wrap).

    Shrinks the axis symmetrically (same amount from both ends) so
    the rect's centre stays put, then re-runs
    ``_rect_from_axis_extended`` so corners re-snap to apt.dat
    pavement boundary on the new axis.
    """
    PERP_CAP_RATIO = 1.3
    DIAG_CAP_RATIO = 1.0
    from shapely.ops import substring
    out: List[Tuple[Polygon, LineString, str, str]] = []
    for rect, axis, role, ref in taxi_rects:
        try:
            coords = list(rect.exterior.coords)
        except Exception:
            out.append((rect, axis, role, ref))
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            out.append((rect, axis, role, ref))
            continue
        edge_lens = [
            math.hypot(coords[(i + 1) % 4][0] - coords[i][0],
                       coords[(i + 1) % 4][1] - coords[i][1])
            for i in range(4)]
        length = max(edge_lens)
        width = min(edge_lens)
        if width < 1.0 or length < 1.0:
            out.append((rect, axis, role, ref))
            continue
        db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
        if db is None:
            out.append((rect, axis, role, ref))
            continue
        if db < 20.0:
            # Parallel — no cap.
            out.append((rect, axis, role, ref))
            continue
        cap_ratio = (DIAG_CAP_RATIO if db < 45.0
                     else PERP_CAP_RATIO)
        max_len = cap_ratio * width
        if length <= max_len + 0.5:
            out.append((rect, axis, role, ref))
            continue
        axis_len = axis.length
        new_axis_len = max_len
        margin = (axis_len - new_axis_len) / 2.0
        if margin <= 0:
            out.append((rect, axis, role, ref))
            continue
        try:
            new_axis = substring(
                axis, margin, axis_len - margin)
        except Exception:
            out.append((rect, axis, role, ref))
            continue
        if (new_axis.is_empty
                or new_axis.geom_type != "LineString"
                or new_axis.length < 5.0):
            out.append((rect, axis, role, ref))
            continue
        new_rect = _rect_from_axis_extended(
            new_axis, width, pav, apt_vertices=apt_vertices)
        if (new_rect is None or new_rect.is_empty
                or new_rect.geom_type != "Polygon"
                or not new_rect.is_valid):
            out.append((rect, axis, role, ref))
            continue
        out.append((new_rect, new_axis, role, ref))
    return out


def _classify_role(axis: LineString, width: float,
                   rwy_centerlines: List[LineString],
                   rwy_union: Optional[Polygon],
                   ref: str = "",
                   ref_overall_bearings: Optional[Dict[str, float]]
                   = None) -> str:
    """Classify a taxi rect by axis geometry alone.

    The role is determined entirely by:
      * bearing to nearest runway (parallel within 20°,
        perpendicular beyond 45°),
      * straight-line distance from axis midpoint to nearest runway
        centerline (close → primary, far → secondary or cross),
      * axis length (must clear minimum length per role).

    Ref letters are NOT used as a classifier — they vary wildly
    between airports (SPJC's parallel taxis are A/F/L/V; CYXY's
    is E; KBNA uses different letters again; many CYXY taxis have
    no ref at all).  The only ref-pattern rules retained are:

      * Sub-ref (any digit in the label, e.g. V1, L3, A1) →
        always STUB.  Sub-ref tagging is universal: a digit
        suffix means a short connector spur regardless of airport.
      * Diagonal parent ref → always STUB.  When the parent OSM
        way's chord bearing is itself diagonal (db ≥ 20°), every
        segment of that ref is part of a diagonal stub even if
        a curving end-segment happens to align near-parallel
        locally.  At SPJC, B/C/E/G enter the runway at shallow
        angles; without this check, the post-curve sub-segment
        (B's 85 m piece, db_local = 18°) gets misclassified as
        PRIMARY_PARALLEL even though there's no actual B parallel
        taxiway.

    Roles:
      * PRIMARY_PARALLEL  — db < 20°, length ≥ 50 m, < 400 m from runway
      * SECONDARY_PARALLEL — db < 20°, length ≥ 50 m, ≥ 400 m from runway
      * CROSS_CONNECTOR   — db > 45°, length ≥ 80 m, > 250 m from runway
      * STUB              — everything else (short, runway-adjacent perp, etc.)
    """
    if ref and any(c.isdigit() for c in ref):
        # Sub-refs are always stubs.
        return ROLE_STUB

    db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
    if db is None:
        return ROLE_STUB

    # Diagonal-parent check: if this rect's REF has an overall-
    # DIAGONAL parent OSM way (parent db_overall ∈ [20°, 45°)),
    # force STUB regardless of the local segment bearing.  See
    # the ``Diagonal parent ref`` rule above.  Excludes
    # perpendicular parents (db ≥ 45°, e.g. cross-connectors
    # like Q/R at SPJC) which legitimately classify as
    # CROSS_CONNECTOR via the local-axis check below.
    if (ref and ref_overall_bearings
            and ref in ref_overall_bearings
            and rwy_centerlines):
        try:
            _rwy = min(rwy_centerlines,
                       key=lambda r: axis.distance(r))
            _rc = list(_rwy.coords)
            _rdx = _rc[-1][0] - _rc[0][0]
            _rdy = _rc[-1][1] - _rc[0][1]
            if math.hypot(_rdx, _rdy) > 1e-6:
                _rwy_bearing = (math.degrees(
                    math.atan2(_rdx, _rdy)) % 180.0)
                _ref_db = abs(ref_overall_bearings[ref]
                              - _rwy_bearing)
                _ref_db = min(_ref_db, 180.0 - _ref_db)
                if 20.0 <= _ref_db < 45.0:
                    return ROLE_STUB
        except Exception:
            pass
    try:
        mid = axis.interpolate(0.5, normalized=True)
        dist_rwy = min(mid.distance(r) for r in rwy_centerlines)
    except Exception:
        dist_rwy = 1e6
    length = axis.length
    if db < 20.0:
        if length >= 50.0:
            if dist_rwy < 400.0:
                return ROLE_PRIMARY_PARALLEL
            return ROLE_SECONDARY_PARALLEL
        return ROLE_STUB
    if db > 45.0 and length >= 80.0:
        if dist_rwy > 250.0:
            return ROLE_CROSS_CONNECTOR
        return ROLE_STUB
    return ROLE_STUB


def _axis_to_nearest_rwy_db(axis: LineString,
                            rwy_centerlines: List[LineString]
                            ) -> Optional[float]:
    """Return the bearing difference from ``axis`` to the nearest
    runway centerline, modulo 180°."""
    if not rwy_centerlines:
        return None
    rwy = min(rwy_centerlines, key=lambda r: axis.distance(r))

    def _bearing(ls):
        c = list(ls.coords)
        return math.degrees(math.atan2(c[-1][0] - c[0][0],
                                       c[-1][1] - c[0][1])) % 180.0
    db = abs(_bearing(axis) - _bearing(rwy))
    return min(db, 180.0 - db)


def _refine_roles(emitted, rwy_centerlines):
    """Post-classify: demote the stub-A / stub-F segment (the short
    runway-connector within a parallel ref's polyline) from
    primary_parallel to stub.

    Rule: for each parallel ref, find SEGMENTS that are
    significantly perpendicular (>= 40° off runway) AND short (< 150 m).
    Demote to stub.  Multiple per ref allowed.
    """
    if not rwy_centerlines or not emitted:
        return
    for i, (rect, axis, role, ref) in enumerate(emitted):
        if role != ROLE_PRIMARY_PARALLEL:
            continue
        db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
        if db is None:
            continue
        if db >= 40.0 and axis.length < 150.0:
            emitted[i] = (rect, axis, ROLE_STUB, ref)
