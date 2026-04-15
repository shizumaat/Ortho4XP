"""Emit one apt.dat taxiway polygon as a chain of sloping rectangles.

A taxiway polygon is a long, thin strip of pavement.  The previous
code triangulated every taxiway along with the apron, producing
dozens of small triangles where a simple chain of 2-3 sloped rects
would describe the same surface more precisely and match the user's
mental model — "almost all taxiways, except for intersections,
should be either flat or a sloping rectangle".

Pipeline:

1. **Fit a min-rotated bounding rectangle** to the polygon.  If the
   polygon fills less than ``min_fit_ratio`` of its MRR (default
   0.80), the shape is too irregular to approximate as a chain of
   rectangles and ``build_taxiway_rects`` returns ``None`` so the
   caller falls back to triangulation.

2. **Densely sample DEM along the long axis.**  Sample spacing is
   set by ``sample_spacing`` (default 20 m).  A 400 m taxiway yields
   21 samples.

3. **Grade-clamp the sample profile.**  Two passes: (a) longitudinal
   cap so ``|z[i+1]-z[i]| ≤ max_grade * dx`` everywhere; (b)
   rate-of-change cap at the FAA taxiway rule of one percent of
   grade change per 30 m, i.e. ``max_dg_per_m = 1/3000``, to stop
   an adjacent-segment kink from violating vertical-curve smoothness.
   The two caps alternate until convergence.

4. **Ramer-Douglas-Peucker simplification** on the clamped (t, z)
   profile with ``fidelity_tol`` tolerance (default 0.3 m).  The
   result is the minimum set of break-points whose piecewise-linear
   reconstruction stays within tolerance of every sample.  For a
   uniformly-sloped taxiway this collapses to ONE segment spanning
   the whole polygon; for a taxiway with a local dip it might
   produce two or three segments.

5. **Emit one sloping rectangle per simplified segment.**  The
   rectangle is aligned with the MRR long axis and has full MRR
   short-side width, so adjacent rectangles share an exact edge.

Each returned :class:`TaxiwayRect` has its 4 meter-space corners
laid out so ``[0, 1]`` is the high-elevation end and ``[2, 3]`` is
the low-elevation end, matching the ``altitude_high``/``altitude_low``
convention used by the legacy runway patch emitter.

This module is pure: no I/O, no shared state.  Tests in
``tests/test_taxiway_rects.py``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple

from shapely.geometry import Polygon


# ──────────────────────────────────────────────────────────────────────
# Tunables
# ──────────────────────────────────────────────────────────────────────

# How densely to sample DEM along the long axis, in meters.  Smaller
# = more accurate capture of local bumps; larger = faster.  20 m is
# roughly one sample per taxiway width — enough to see a terrain
# bump without being noisy.
DEFAULT_SAMPLE_SPACING_M = 20.0

# RDP tolerance for collapsing adjacent samples into a single sloping
# segment.  0.3 m ~ a gentle rise/fall below visual pixel threshold.
DEFAULT_FIDELITY_TOL_M = 0.3

# Fit-ratio threshold below which the polygon is too irregular to
# approximate with an MRR-aligned rect chain.  0.80 handles simple
# taxiway strips with Bezier-curved ends; anything chunkier falls
# through to triangulation.
DEFAULT_MIN_FIT_RATIO = 0.80

# FAA vertical-curve rate-of-change rule for taxiways: one percent of
# grade per 30 m.  This is the looser taxiway rule; runways use a
# stricter 1 % per 305 m.
DEFAULT_MAX_DG_PER_M = 1.0 / 3000.0

# Maximum iterations for the combined longitudinal-cap / rate-of-change
# convergence loop.
_CLAMP_MAX_ITERS = 50


# ──────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class TaxiwayRect:
    """One sloping-rect fragment of a taxiway chain.

    Corners are arranged so ``[0]``, ``[1]`` are at the elev_high
    end and ``[2]``, ``[3]`` are at the elev_low end, matching
    :func:`O4_Auto_Patch._emit_sloped_rect`'s expectation that
    ``altitude_high`` applies to corners 0-1.

    For a flat rect (``elev_high == elev_low``) the corner
    ordering is still consistent (corners 0-1 at one end, 2-3 at
    the other).
    """
    corners_m: Tuple[Tuple[float, float], ...]   # 4 corners
    elev_low: float
    elev_high: float
    # Centerline endpoints at the high and low ends, for
    # convenient re-projection to lat/lon by the caller.
    center_high: Tuple[float, float]
    center_low: Tuple[float, float]
    width_m: float

    @property
    def is_flat(self) -> bool:
        return abs(self.elev_high - self.elev_low) < 0.1

    @property
    def polygon(self) -> Polygon:
        return Polygon(self.corners_m)


# ──────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────
def _long_axis(polygon: Polygon
               ) -> Optional[Tuple[Tuple[float, float],
                                   Tuple[float, float],
                                   float, float,
                                   float, float]]:
    """Return the long-axis geometry of the polygon's min-rotated rect.

    Result is ``(m_a, m_b, ux, uy, long_len, short_len)`` where:

      * ``m_a`` and ``m_b`` are the midpoints of the two short edges
        of the MRR (centerline endpoints).
      * ``(ux, uy)`` is a unit vector from ``m_a`` to ``m_b``.
      * ``long_len`` is the MRR long-side length.
      * ``short_len`` is the MRR short-side length (= taxiway width).

    Returns ``None`` for an empty/degenerate polygon.
    """
    if polygon is None or polygon.is_empty:
        return None
    try:
        mrr = polygon.minimum_rotated_rectangle
    except Exception:
        return None
    if mrr is None or mrr.is_empty or not hasattr(mrr, "exterior"):
        return None
    coords = list(mrr.exterior.coords)
    if len(coords) < 5:
        return None
    c = coords[:4]
    e1_len = math.hypot(c[1][0] - c[0][0], c[1][1] - c[0][1])
    e2_len = math.hypot(c[2][0] - c[1][0], c[2][1] - c[1][1])

    if e1_len >= e2_len:
        # c[0]-c[1] is long.  Short edges are c[3]-c[0] and c[1]-c[2].
        m_a = ((c[3][0] + c[0][0]) / 2.0, (c[3][1] + c[0][1]) / 2.0)
        m_b = ((c[1][0] + c[2][0]) / 2.0, (c[1][1] + c[2][1]) / 2.0)
        long_len, short_len = e1_len, e2_len
    else:
        # c[1]-c[2] is long.  Short edges are c[0]-c[1] and c[2]-c[3].
        m_a = ((c[0][0] + c[1][0]) / 2.0, (c[0][1] + c[1][1]) / 2.0)
        m_b = ((c[2][0] + c[3][0]) / 2.0, (c[2][1] + c[3][1]) / 2.0)
        long_len, short_len = e2_len, e1_len

    dx = m_b[0] - m_a[0]
    dy = m_b[1] - m_a[1]
    length = math.hypot(dx, dy)
    if length < 1e-6:
        return None
    ux, uy = dx / length, dy / length
    return (m_a, m_b, ux, uy, long_len, short_len)


def _clamp_profile(zs: List[float],
                   seg_len: float,
                   max_grade: float,
                   max_dg_per_m: float) -> None:
    """In-place iterative grade + rate-of-change clamp on a 1D DEM
    profile ``zs`` taken at uniform intervals of ``seg_len`` meters.

    Two passes, alternated until convergence or iter cap:

      1. Longitudinal cap — ``|zs[i+1] - zs[i]| ≤ max_grade*seg_len``
         (split the excess symmetrically between the two endpoints).
      2. Rate-of-change — the grade change between adjacent segments
         is at most ``max_dg_per_m * seg_len``, enforced by nudging
         the middle sample of a (i-1, i, i+1) triple toward the
         expected value ``(zs[i-1]+zs[i+1])/2`` when the deviation
         exceeds the budget.

    Modifies ``zs`` in place.
    """
    if len(zs) < 2 or seg_len <= 0:
        return
    max_dz = max_grade * seg_len
    # The grade-change budget over one segment: change in grade is
    # max_dg_per_m * seg_len, translated to elevation it's
    # max_dg_per_m * seg_len**2.
    grade_change_budget = max_dg_per_m * seg_len * seg_len

    for _ in range(_CLAMP_MAX_ITERS):
        changed = False

        # Longitudinal cap: classic hard-cap "range clamp" around
        # the mean of the profile.  Pulling each sample into the
        # tightest [lower, upper] band imposed by both neighbours
        # converges in a single pass and does not drift upward.
        # The band at index i is:
        #   upper = min_j(zs[j] + max_dz*|i-j|)
        #   lower = max_j(zs[j] - max_dz*|i-j|)
        # which we compute incrementally with two sweeps.
        #
        # Forward pass: upper[i] = min(upper[i-1] + max_dz, upper[i])
        # Backward pass: upper[i] = min(upper[i+1] + max_dz, upper[i])
        # (And lower is the mirror.)
        upper = list(zs)
        lower = list(zs)
        for i in range(1, len(zs)):
            upper[i] = min(upper[i], upper[i - 1] + max_dz)
            lower[i] = max(lower[i], lower[i - 1] - max_dz)
        for i in range(len(zs) - 2, -1, -1):
            upper[i] = min(upper[i], upper[i + 1] + max_dz)
            lower[i] = max(lower[i], lower[i + 1] - max_dz)
        for i in range(len(zs)):
            new_z = zs[i]
            if new_z > upper[i]:
                new_z = upper[i]
            if new_z < lower[i]:
                new_z = lower[i]
            if abs(new_z - zs[i]) > 1e-9:
                zs[i] = new_z
                changed = True

        # Rate-of-change cap: the middle of any triple should be
        # close to the linear interpolation of its neighbours.
        for i in range(1, len(zs) - 1):
            expected = (zs[i - 1] + zs[i + 1]) / 2.0
            dev = zs[i] - expected
            if dev > grade_change_budget + 1e-9:
                zs[i] = expected + grade_change_budget
                changed = True
            elif dev < -grade_change_budget - 1e-9:
                zs[i] = expected - grade_change_budget
                changed = True

        if not changed:
            break


def _rdp_simplify_indices(zs: List[float],
                          tolerance: float) -> List[int]:
    """Ramer-Douglas-Peucker on a 1D z profile sampled at uniform
    intervals.  Returns the sorted list of kept indices.

    Starts with ``[0, len-1]``; recursively splits at the sample
    with the largest deviation from the current piecewise-linear
    reconstruction until every deviation is ``≤ tolerance``.
    """
    n = len(zs)
    if n <= 2:
        return list(range(n))

    keep = {0, n - 1}

    def recurse(lo: int, hi: int) -> None:
        if hi - lo <= 1:
            return
        z_lo = zs[lo]
        z_hi = zs[hi]
        span = hi - lo
        worst_i = -1
        worst_err = 0.0
        for i in range(lo + 1, hi):
            frac = (i - lo) / span
            z_fit = z_lo * (1.0 - frac) + z_hi * frac
            err = abs(zs[i] - z_fit)
            if err > worst_err:
                worst_err = err
                worst_i = i
        if worst_err > tolerance and worst_i > 0:
            keep.add(worst_i)
            recurse(lo, worst_i)
            recurse(worst_i, hi)

    recurse(0, n - 1)
    return sorted(keep)


# ──────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────
def build_taxiway_rects(
    polygon: Polygon,
    sample_dem: Callable[[float, float], Optional[float]],
    max_grade: float = 0.015,
    sample_spacing: float = DEFAULT_SAMPLE_SPACING_M,
    fidelity_tol: float = DEFAULT_FIDELITY_TOL_M,
    min_fit_ratio: float = DEFAULT_MIN_FIT_RATIO,
    max_dg_per_m: float = DEFAULT_MAX_DG_PER_M,
) -> Optional[List[TaxiwayRect]]:
    """Build the rect chain for one taxiway polygon.

    Args:
        polygon: classified-taxiway apt.dat pavement, in METER space.
        sample_dem: callable ``(x, y) → elevation_m | None``.
        max_grade: longitudinal grade cap (e.g. 0.015 for 1.5 %).
        sample_spacing: target distance between DEM samples along
            the long axis.
        fidelity_tol: RDP tolerance for collapsing samples.
        min_fit_ratio: reject the polygon and return ``None`` if
            ``polygon.area / MRR.area`` is below this threshold —
            the caller should fall back to triangulation.
        max_dg_per_m: FAA taxiway vertical-curve rate-of-change
            constant (``1/3000`` m⁻¹).

    Returns:
        A list of :class:`TaxiwayRect` covering the polygon's MRR,
        or ``None`` if the polygon is not strip-like enough.
    """
    if polygon is None or polygon.is_empty:
        return None

    try:
        mrr = polygon.minimum_rotated_rectangle
    except Exception:
        return None
    if mrr is None or mrr.is_empty:
        return None
    try:
        fit_ratio = (polygon.area / mrr.area) if mrr.area > 0 else 0.0
    except Exception:
        fit_ratio = 0.0

    ax = _long_axis(polygon)
    if ax is None:
        return None
    m_a, m_b, ux, uy, long_len, short_len = ax

    # Strip-shape gate.  A polygon is "strip-like enough" to model
    # as an MRR-aligned rect chain when:
    #
    #   * its MRR short side lies within the taxiway width envelope
    #     (9 – 45 m), AND
    #   * its aspect ratio is ≥ 3.5, AND
    #   * it fills at least a scaled fraction of its MRR.
    #
    # Scaled fit threshold: a true rectangle has fit = 1.0; real
    # taxiways have gentle curves, splays at runway touch points,
    # and chamfered corners that depress fit to 0.55 – 0.85 even
    # when the underlying shape is clearly a single strip.
    # High-aspect curved strips get a lenient threshold; short
    # near-square "taxiway pads" get the strict one.
    #
    # L-shapes and C-shapes — the dangerous false-positives where
    # an MRR-aligned rect would massively overshoot the polygon —
    # have bbox aspect ≤ ~1.5 because the bounding rectangle wraps
    # both arms of the bend.  The aspect gate rejects them.
    aspect = (long_len / short_len) if short_len > 0 else 0.0
    if short_len < 9.0 or short_len > 45.0:
        return None
    if aspect < 3.5:
        return None
    if aspect >= 15.0:
        adaptive_min = 0.40
    elif aspect >= 8.0:
        adaptive_min = 0.50
    elif aspect >= 5.0:
        adaptive_min = 0.60
    else:
        adaptive_min = 0.70
    if fit_ratio < adaptive_min:
        return None
    # Perpendicular (90° CCW rotation of the long axis).
    px, py = -uy, ux

    # Dense sample count — at least two, so we always have a
    # beginning and an end.
    n_samples = max(2, int(round(long_len / sample_spacing)) + 1)
    seg_len = long_len / (n_samples - 1)

    zs: List[float] = []
    xs_center: List[Tuple[float, float]] = []
    for i in range(n_samples):
        t = i * seg_len
        cx = m_a[0] + ux * t
        cy = m_a[1] + uy * t
        z = sample_dem(cx, cy)
        if z is None:
            # A missing DEM sample falls back to the nearest
            # neighbour we've already seen — gives a continuous
            # profile even on incomplete DEM.
            z = zs[-1] if zs else 0.0
        zs.append(z)
        xs_center.append((cx, cy))

    _clamp_profile(zs, seg_len, max_grade, max_dg_per_m)

    keep_indices = _rdp_simplify_indices(zs, fidelity_tol)
    if len(keep_indices) < 2:
        keep_indices = [0, n_samples - 1]

    half_w = short_len / 2.0

    rects: List[TaxiwayRect] = []
    for k in range(len(keep_indices) - 1):
        i_a = keep_indices[k]
        i_b = keep_indices[k + 1]
        ax_a = xs_center[i_a]
        ax_b = xs_center[i_b]
        z_a = zs[i_a]
        z_b = zs[i_b]

        # Decide which end is "high" for altitude_high corner 0-1.
        if z_a >= z_b:
            c_high, c_low = ax_a, ax_b
            eh, el = z_a, z_b
        else:
            c_high, c_low = ax_b, ax_a
            eh, el = z_b, z_a

        c0 = (c_high[0] + px * half_w, c_high[1] + py * half_w)
        c1 = (c_high[0] - px * half_w, c_high[1] - py * half_w)
        c2 = (c_low[0] - px * half_w, c_low[1] - py * half_w)
        c3 = (c_low[0] + px * half_w, c_low[1] + py * half_w)

        rects.append(TaxiwayRect(
            corners_m=(c0, c1, c2, c3),
            elev_low=el,
            elev_high=eh,
            center_high=c_high,
            center_low=c_low,
            width_m=short_len,
        ))

    return rects
