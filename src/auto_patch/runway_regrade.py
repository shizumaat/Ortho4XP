"""Stage A: runway profile regrading against seam HARD anchors.

When a runway crosses one or more integer lat/lon tile-boundary lines,
its CIFP-derived linear profile may not pass through the DEM-anchored
seam altitudes.  This module re-optimises the runway threshold
altitudes to:

  1. Honor the seam altitudes exactly (HARD constraint — seam wins by
     user 2026-05-13 design rule).
  2. Stay grade-compliant per FAA AC 150/5300-13:
       - longitudinal grade ≤ ``grade_cap`` (default 1.5%);
       - vertical curve K-factor ``L = K × |Δg|`` (default K = 305 m,
         ARC Cat C/D).
  3. Minimise deviation from CIFP threshold altitudes.

When constraints can't all be satisfied (very steep DEM, short runway),
relaxation order is: K-factor first (warning), then grade cap
(warning).  Seam altitudes are never relaxed.

The only unknowns are the two threshold altitudes — every interior
seam vertex is HARD.  Interior segment grades / curvature are
determined entirely by the seam altitudes; we can only warn about
violations there.

Public API: ``regrade_runway``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple


__all__ = ["regrade_runway", "RegradeResult", "regrade_runways_in_layout"]


@dataclass
class RegradeResult:
    threshold_A: float
    threshold_B: float
    warnings: List[str]
    # Per-seam altitudes for downstream solver to apply (unchanged from
    # input; included for convenience so callers don't shuttle the
    # seam list separately).
    seam_altitudes: List[Tuple[float, float]]


# FAA Approach Category C/D K-factor for vertical curves.  Smaller K
# applies to lighter ARCs (A/B: K ≈ 76 m); larger to heavy (E: 610 m).
DEFAULT_ARC_K_M = 305.0   # ARC C/D (most jet airports)
DEFAULT_GRADE_CAP = 0.015  # 1.5% longitudinal


def regrade_runway(
    threshold_A_cifp: float,
    threshold_B_cifp: float,
    runway_axis_length: float,
    seam_anchors: List[Tuple[float, float]],
    *,
    grade_cap: float = DEFAULT_GRADE_CAP,
    arc_K_m: float = DEFAULT_ARC_K_M,
) -> RegradeResult:
    """Re-grade a runway against fixed seam altitudes.

    Parameters:
        threshold_A_cifp: CIFP-derived altitude at runway threshold A.
        threshold_B_cifp: same for threshold B.
        runway_axis_length: total runway centerline length (m).
        seam_anchors: list of ``(dist_from_A, dem_alt)`` tuples;
            distance must be in (0, runway_axis_length).  Sorted
            internally.  An empty list short-circuits to CIFP.

    Returns:
        RegradeResult with adjusted threshold altitudes and any
        constraint-relaxation warnings.
    """
    warnings: List[str] = []
    if not seam_anchors:
        return RegradeResult(
            threshold_A=threshold_A_cifp,
            threshold_B=threshold_B_cifp,
            warnings=warnings,
            seam_altitudes=[])

    seams = sorted(seam_anchors, key=lambda s: s[0])
    # Filter out malformed entries (out-of-range distances).
    seams = [(d, a) for (d, a) in seams
             if 0.0 < d < runway_axis_length]
    if not seams:
        return RegradeResult(
            threshold_A=threshold_A_cifp,
            threshold_B=threshold_B_cifp,
            warnings=warnings,
            seam_altitudes=[])

    n = len(seams)
    # Segment distances along the runway: d[0] = A→seam[0],
    # d[i] = seam[i-1]→seam[i] for i in 1..n-1, d[n] = seam[n-1]→B.
    dists: List[float] = [seams[0][0]]
    for i in range(1, n):
        dists.append(seams[i][0] - seams[i - 1][0])
    dists.append(runway_axis_length - seams[-1][0])
    # Altitudes at each interior node (just the seams).
    alts = [s[1] for s in seams]

    # Warn about interior segment grade / K-factor violations — these
    # are determined entirely by seam altitudes and we can't fix them.
    for i in range(n - 1):
        d = dists[i + 1]
        if d > 1e-6:
            g = (alts[i + 1] - alts[i]) / d
            if abs(g) > grade_cap:
                warnings.append(
                    f"interior runway segment between seam {i} and "
                    f"seam {i+1}: grade {g*100:.2f}% > cap "
                    f"{grade_cap*100:.2f}% (seam altitudes immutable)")
    for i in range(1, n - 1):
        d_prev = dists[i]
        d_next = dists[i + 1]
        if d_prev > 1e-6 and d_next > 1e-6:
            g_prev = (alts[i] - alts[i - 1]) / d_prev
            g_next = (alts[i + 1] - alts[i]) / d_next
            available = 2.0 * min(d_prev, d_next)
            dg_max = available / arc_K_m
            if abs(g_next - g_prev) > dg_max:
                warnings.append(
                    f"interior runway PVI at seam {i}: |Δg|"
                    f"={abs(g_next-g_prev)*100:.2f}% > K-factor "
                    f"limit {dg_max*100:.2f}% (curve cannot fit)")

    # ── Threshold A regrade ─────────────────────────────────────────
    d_A = dists[0]
    seam_alt_0 = alts[0]
    # Grade-cap bound: |g_0| ≤ grade_cap where g_0 = (seam_alt_0 - alt_A) / d_A
    a_lo_grade = seam_alt_0 - grade_cap * d_A
    a_hi_grade = seam_alt_0 + grade_cap * d_A
    # K-factor bound at PVI 0:
    if n >= 2:
        # g_1 is fixed (seam-to-seam), so g_0 bounds are linear in alt_A.
        g_1_fixed = (alts[1] - alts[0]) / dists[1] if dists[1] > 1e-6 else 0.0
        available = 2.0 * min(d_A, dists[1])
        dg_max = available / arc_K_m
        # |g_1 - g_0| ≤ dg_max → g_0 in [g_1 - dg_max, g_1 + dg_max]
        # → alt_A = seam_alt_0 - g_0*d_A
        a_lo_K = seam_alt_0 - (g_1_fixed + dg_max) * d_A
        a_hi_K = seam_alt_0 - (g_1_fixed - dg_max) * d_A
    else:
        # Single seam (n=1): K-factor couples alt_A and alt_B.
        # Defer K-factor to the joint pass below.
        a_lo_K, a_hi_K = -math.inf, math.inf
    a_lo = max(a_lo_grade, a_lo_K)
    a_hi = min(a_hi_grade, a_hi_K)
    if a_lo > a_hi:
        warnings.append(
            f"threshold A: grade cap + K-factor infeasible; relaxing K")
        a_lo, a_hi = a_lo_grade, a_hi_grade
    alt_A = max(a_lo, min(a_hi, threshold_A_cifp))

    # ── Threshold B regrade ─────────────────────────────────────────
    d_B = dists[-1]
    seam_alt_last = alts[-1]
    # g_N = (alt_B - seam_alt_last) / d_B; |g_N| ≤ grade_cap
    b_lo_grade = seam_alt_last - grade_cap * d_B
    b_hi_grade = seam_alt_last + grade_cap * d_B
    if n >= 2:
        g_Nm1_fixed = ((alts[-1] - alts[-2]) / dists[-2]
                       if dists[-2] > 1e-6 else 0.0)
        available = 2.0 * min(d_B, dists[-2])
        dg_max = available / arc_K_m
        b_lo_K = seam_alt_last + (g_Nm1_fixed - dg_max) * d_B
        b_hi_K = seam_alt_last + (g_Nm1_fixed + dg_max) * d_B
    else:
        b_lo_K, b_hi_K = -math.inf, math.inf
    b_lo = max(b_lo_grade, b_lo_K)
    b_hi = min(b_hi_grade, b_hi_K)
    if b_lo > b_hi:
        warnings.append(
            f"threshold B: grade cap + K-factor infeasible; relaxing K")
        b_lo, b_hi = b_lo_grade, b_hi_grade
    alt_B = max(b_lo, min(b_hi, threshold_B_cifp))

    # ── Joint K-factor for n=1 ─────────────────────────────────────
    if n == 1:
        seam_alt = alts[0]
        g_0 = (seam_alt - alt_A) / d_A if d_A > 1e-6 else 0.0
        g_1 = (alt_B - seam_alt) / d_B if d_B > 1e-6 else 0.0
        available = 2.0 * min(d_A, d_B)
        dg_max = available / arc_K_m
        if abs(g_1 - g_0) > dg_max:
            # Project (alt_A, alt_B) onto the K-factor constraint
            # boundary, then re-clip to grade-cap bounds.
            # K-factor constraint: g_1 - g_0 = ±dg_max
            #   → (alt_B - seam)/d_B - (seam - alt_A)/d_A = ±dg_max
            #   → alt_A*d_B + alt_B*d_A = seam*(d_A+d_B) ± dg_max*d_A*d_B
            sign = 1.0 if (g_1 - g_0) > 0 else -1.0
            target = (seam_alt * (d_A + d_B)
                      + sign * dg_max * d_A * d_B)
            current = alt_A * d_B + alt_B * d_A
            # Project: move by t * (d_B, d_A) where t solves
            #   (alt_A - t*d_B)*d_B + (alt_B - t*d_A)*d_A = target
            #   alt_A*d_B + alt_B*d_A - t*(d_B² + d_A²) = target
            #   t = (current - target) / (d_B² + d_A²)
            denom = d_B * d_B + d_A * d_A
            if denom > 1e-9:
                t = (current - target) / denom
                alt_A = alt_A - t * d_B
                alt_B = alt_B - t * d_A
            # Re-clip to grade bounds (K-factor projection may push
            # threshold outside grade cap; grade cap wins).
            alt_A = max(a_lo_grade, min(a_hi_grade, alt_A))
            alt_B = max(b_lo_grade, min(b_hi_grade, alt_B))
            # Verify final state; warn if K-factor still violated.
            g_0 = (seam_alt - alt_A) / d_A
            g_1 = (alt_B - seam_alt) / d_B
            if abs(g_1 - g_0) > dg_max * 1.01:
                warnings.append(
                    f"single-seam runway: K-factor not satisfied "
                    f"alongside grade cap (|Δg|={abs(g_1-g_0)*100:.2f}% "
                    f"> K limit {dg_max*100:.2f}%); grade cap kept")

    # Annotate threshold shifts in the warnings.
    da = alt_A - threshold_A_cifp
    db = alt_B - threshold_B_cifp
    if abs(da) > 0.05:
        warnings.append(
            f"threshold A shifted {da:+.2f}m from CIFP "
            f"({threshold_A_cifp:.2f}→{alt_A:.2f})")
    if abs(db) > 0.05:
        warnings.append(
            f"threshold B shifted {db:+.2f}m from CIFP "
            f"({threshold_B_cifp:.2f}→{alt_B:.2f})")

    return RegradeResult(
        threshold_A=alt_A,
        threshold_B=alt_B,
        warnings=warnings,
        seam_altitudes=list(seams))


def regrade_runways_in_layout(
    layout,
    dem,
    tile_lat: int,
    tile_lon: int,
    *,
    grade_cap: float = DEFAULT_GRADE_CAP,
    arc_K_m: float = DEFAULT_ARC_K_M,
) -> int:
    """Stage A integration: for each runway shape in ``layout`` that
    has seam-inserted vertices, regrade its threshold altitudes and
    write the adjusted values back into ``node_altitudes``.

    Assumes ``layout`` has been processed by
    ``seam_anchors.split_pavement_at_seams`` (runway sloped rects
    converted to ``node_altitudes`` with seam vertices inserted in
    the ring) and ``layout._seam_anchor_keys`` is populated.

    Returns the number of runway shapes that were re-graded.
    """
    from .layout import ROLE_RUNWAY, SHARED_VERTEX_TOL_M
    import math

    anchor_keys = getattr(layout, "_seam_anchor_keys", None)
    if not anchor_keys:
        return 0

    bs = 1.0 / SHARED_VERTEX_TOL_M  # bucket scale

    def _is_seam(x: float, y: float) -> bool:
        return (int(round(x * bs)), int(round(y * bs))) in anchor_keys

    n_regraded = 0
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        if not s.node_altitudes:
            continue
        if s.source_axis is None:
            continue
        ring = list(s.polygon.exterior.coords)
        if ring and ring[0] == ring[-1]:
            ring = ring[:-1]
        if len(ring) < 4:
            continue
        n = len(ring)
        alts = list(s.node_altitudes[:n])
        # Find seam vertices and threshold vertices.
        seam_idxs = [i for i, (x, y) in enumerate(ring) if _is_seam(x, y)]
        if not seam_idxs:
            continue  # no seam crossings on this runway
        thresh_idxs = [i for i in range(n) if i not in seam_idxs]
        # Threshold corners are at the H and L ends of the original
        # 4-corner [H, L, L, H] convention.  With seam vertices
        # inserted, the threshold corners are the non-seam vertices.
        # Project each onto source_axis to determine which end each
        # belongs to.
        axis_coords = list(s.source_axis.coords)
        if len(axis_coords) < 2:
            continue
        ax_x0, ax_y0 = axis_coords[0]
        ax_x1, ax_y1 = axis_coords[-1]
        ax_dx = ax_x1 - ax_x0
        ax_dy = ax_y1 - ax_y0
        ax_len2 = ax_dx * ax_dx + ax_dy * ax_dy
        if ax_len2 < 1e-9:
            continue
        ax_len = math.sqrt(ax_len2)

        def _project(x: float, y: float) -> float:
            """Distance along source_axis from start (in metres)."""
            return ((x - ax_x0) * ax_dx + (y - ax_y0) * ax_dy) / ax_len

        # Group threshold vertices by axis end (low t = "A", high t = "B").
        thresh_with_t = [(i, _project(*ring[i])) for i in thresh_idxs]
        if not thresh_with_t:
            continue
        thresh_with_t.sort(key=lambda x: x[1])
        # Bucket into two clusters by t value — there should be 2 of each.
        # Use the midpoint between min and max as the split.
        t_min = thresh_with_t[0][1]
        t_max = thresh_with_t[-1][1]
        t_mid = 0.5 * (t_min + t_max)
        A_idxs = [i for (i, t) in thresh_with_t if t <= t_mid]
        B_idxs = [i for (i, t) in thresh_with_t if t > t_mid]
        if not A_idxs or not B_idxs:
            continue
        # CIFP threshold altitudes: average the alts at each end.
        cifp_A = sum(alts[i] for i in A_idxs) / len(A_idxs)
        cifp_B = sum(alts[i] for i in B_idxs) / len(B_idxs)
        # Seam anchors: project each seam vertex onto axis, sample DEM.
        seam_inputs = []
        seam_dem_alts: List[Tuple[int, float]] = []  # (ring_idx, dem_alt)
        for i in seam_idxs:
            x, y = ring[i]
            t = _project(x, y)
            lat, lon = layout.m_to_ll(x, y)
            try:
                dem_alt = float(dem.alt_strict(
                    (lon - tile_lon, lat - tile_lat)))
            except (IndexError, ValueError, TypeError):
                dem_alt = float("nan")
            if dem_alt != dem_alt or dem_alt == dem.nodata:   # NaN or NODATA
                continue
            seam_inputs.append((t, dem_alt))
            seam_dem_alts.append((i, dem_alt))
        if not seam_inputs:
            continue
        # Distance along axis = t along the axis vector (already in m).
        # Make sure thresholds are at t≈0 and t≈ax_len.
        result = regrade_runway(
            cifp_A, cifp_B, ax_len, seam_inputs,
            grade_cap=grade_cap, arc_K_m=arc_K_m)
        # Write back: thresholds at axial ends, seam vertices at DEM.
        for i in A_idxs:
            alts[i] = round(result.threshold_A, 1)
        for i in B_idxs:
            alts[i] = round(result.threshold_B, 1)
        for i, dem_alt in seam_dem_alts:
            alts[i] = round(dem_alt, 1)
        s.node_altitudes = alts + [alts[0]]  # close ring
        n_regraded += 1
        if result.warnings:
            from O4_UI_Utils import vprint
            for w in result.warnings:
                vprint(1, f"  [pav-builder] runway {s.ref}: {w}")
    return n_regraded
