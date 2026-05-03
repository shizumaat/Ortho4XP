"""Phase 1: per-rect axial DEM-smoothed elevation profile.

For each taxi rect (primary_parallel / secondary_parallel / stub /
cross_connector), compute the rect's two short-end altitudes
(``altitude_high`` / ``altitude_low``) by:

1. Sampling DEM at ``t=0``, ``t=0.5``, ``t=1`` along source_axis.
2. Replacing samples at runway-shared short-ends with the runway
   corner's HARD elevation.
3. Clamping the start→mid and mid→end gradients to 1.5 % per axial
   metre (TAXI_MAX_GRADE).

Produces only the per-rect axial extremes — the connecting BFS pass
(Phase 2) reconciles cross-rect propagation.
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

from .dem_targets import rect_axial_targets

TAXI_MAX_GRADE = 0.015  # 1.5 %, FAA AC 150/5300-13B


SLOPING_RECT_ROLES = (
    "primary_parallel",
    "secondary_parallel",
    "stub",
    "cross_connector",
)


def clamp_three_point_profile(
        e_start: float, e_mid: float, e_end: float,
        d_sm: float, d_me: float,
        max_grade: float = TAXI_MAX_GRADE,
        anchor: str = "none",
) -> Tuple[float, float, float]:
    """Adjust ``(e_start, e_mid, e_end)`` so that
    ``|delta| / segment_length ≤ max_grade`` between consecutive
    samples.

    ``anchor`` selects which sample is HARD (frozen):
      * ``"start"`` — start is HARD, mid/end may move.
      * ``"end"`` — end is HARD, start/mid may move.
      * ``"none"`` — none HARD; mid is the swing point and end-points
        move toward DEM as much as the cap permits.

    Movement strategy: pull toward the unanchored DEM target, but
    clip to anchored ± ``max_grade × distance``.  Two passes
    (start→mid→end and the reverse) ensure both segment caps are
    satisfied without amplifying DEM noise.
    """
    cap_sm = max_grade * d_sm
    cap_me = max_grade * d_me
    s, m, e = e_start, e_mid, e_end
    if anchor == "start":
        # s is fixed.  m clamped to s ± cap_sm; e clamped to m ± cap_me.
        m = _clip_to_anchor(m, s, cap_sm)
        e = _clip_to_anchor(e, m, cap_me)
    elif anchor == "end":
        m = _clip_to_anchor(m, e, cap_me)
        s = _clip_to_anchor(s, m, cap_sm)
    else:
        # Both ends free: enforce caps symmetrically.  Iterate twice
        # to settle both caps.
        for _ in range(2):
            m_min = max(s - cap_sm, e - cap_me)
            m_max = min(s + cap_sm, e + cap_me)
            if m_min > m_max:
                # Caps disagree; collapse mid to the midpoint of the
                # disagreement range (rest of the pipeline handles
                # the residual through cross-shape continuity).
                m = (m_min + m_max) / 2.0
            else:
                m = max(m_min, min(m_max, m))
            s = _clip_to_anchor(s, m, cap_sm)
            e = _clip_to_anchor(e, m, cap_me)
    return s, m, e


def _clip_to_anchor(target: float, anchor: float, cap: float) -> float:
    if target > anchor + cap:
        return anchor + cap
    if target < anchor - cap:
        return anchor - cap
    return target


def solve_rect_axial_profile(
        layout, shape, dem, tile_lat: int, tile_lon: int,
        anchor_start: Optional[float] = None,
        anchor_end: Optional[float] = None,
) -> Optional[Tuple[float, float]]:
    """Compute ``(altitude_high, altitude_low)`` for one rect.

    ``anchor_start`` / ``anchor_end`` are HARD elevations at the
    rect's source_axis start / end (when shared with a runway corner
    or another already-solved surface).  ``None`` means free.

    Returns ``None`` if DEM samples are unavailable (rect remains
    untouched, falling back to upstream initial state).
    """
    targets = rect_axial_targets(
        layout, dem, tile_lat, tile_lon, shape.source_axis)
    if targets is None:
        return None
    e_start, e_mid, e_end = targets
    pts = list(shape.source_axis.coords)
    ax_len = math.hypot(
        pts[-1][0] - pts[0][0], pts[-1][1] - pts[0][1])
    d_sm = ax_len / 2.0
    d_me = ax_len / 2.0
    if anchor_start is not None:
        e_start = float(anchor_start)
    if anchor_end is not None:
        e_end = float(anchor_end)
    if anchor_start is not None and anchor_end is not None:
        # Both ends HARD: only the midpoint is free.
        s = e_start
        e = e_end
        m = _clamp_mid_between(e_mid, s, e, d_sm, d_me)
    elif anchor_start is not None:
        s, m, e = clamp_three_point_profile(
            e_start, e_mid, e_end, d_sm, d_me, anchor="start")
    elif anchor_end is not None:
        s, m, e = clamp_three_point_profile(
            e_start, e_mid, e_end, d_sm, d_me, anchor="end")
    else:
        s, m, e = clamp_three_point_profile(
            e_start, e_mid, e_end, d_sm, d_me, anchor="none")
    return max(s, e), min(s, e)


def _clamp_mid_between(target: float, s: float, e: float,
                       d_sm: float, d_me: float,
                       max_grade: float = TAXI_MAX_GRADE) -> float:
    """Clamp ``target`` (DEM mid sample) into the range allowed by
    cap × distance from both fixed endpoints ``s`` and ``e``.
    Falls back to the midpoint of the disagreement range when the
    constraints are mutually infeasible.
    """
    cap_sm = max_grade * d_sm
    cap_me = max_grade * d_me
    lo = max(s - cap_sm, e - cap_me)
    hi = min(s + cap_sm, e + cap_me)
    if lo > hi:
        return (lo + hi) / 2.0
    return max(lo, min(hi, target))


def _axial_end_anchor(shape, end_pair_indices, runway_anchor_lookup
                       ) -> Optional[float]:
    """Average runway-corner anchors for the two polygon corners at
    one short end of a rect.  Returns ``None`` if neither corner
    has a runway anchor (the axial end is free).
    """
    if shape.polygon is None or shape.polygon.is_empty:
        return None
    coords = list(shape.polygon.exterior.coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    found = []
    for idx in end_pair_indices:
        if 0 <= idx < len(coords):
            x, y = coords[idx]
            a = runway_anchor_lookup(x, y)
            if a is not None:
                found.append(a)
    if not found:
        return None
    return sum(found) / len(found)


def apply_axial_profiles(layout, dem, tile_lat: int, tile_lon: int,
                         runway_anchor_lookup) -> int:
    """Phase 1 main entry: solve every rect's axial profile from DEM,
    using ``runway_anchor_lookup`` to detect runway-shared short ends.

    Anchors are looked up at the rect's POLYGON CORNERS (short-end
    pairs, identified via ``_short_end_pairs_by_axis``), not the
    centerline midpoints — only the polygon corners share buckets
    with runway corners.

    Returns the number of rects updated.
    """
    from auto_patch.elevation import _short_end_pairs_by_axis

    n_updated = 0
    for shape in layout.shapes:
        if shape.role not in SLOPING_RECT_ROLES:
            continue
        if shape.polygon is None or shape.polygon.is_empty:
            continue
        if shape.source_axis is None or shape.source_axis.is_empty:
            continue
        coords = list(shape.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        anchor_start: Optional[float] = None
        anchor_end: Optional[float] = None
        if len(coords) == 4:
            sp, ep = _short_end_pairs_by_axis(coords, shape.source_axis)
            if sp is not None and ep is not None:
                anchor_start = _axial_end_anchor(
                    shape, sp, runway_anchor_lookup)
                anchor_end = _axial_end_anchor(
                    shape, ep, runway_anchor_lookup)
        result = solve_rect_axial_profile(
            layout, shape, dem, tile_lat, tile_lon,
            anchor_start, anchor_end)
        if result is None:
            continue
        high, low = result
        shape.altitude_high = round(float(high), 1)
        shape.altitude_low = round(float(low), 1)
        shape.altitude = None
        n_updated += 1
    return n_updated
