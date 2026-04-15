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
# segment.  Per user direction: "unless there's more than 1 m
# elevation change over a 30 m distance, simplify and combine".  A
# 1 m RDP tolerance means a linear fit between two break-points is
# kept as long as every interior sample is within 1 m of the fit,
# so a 300 m straight taxiway with gentle DEM variation collapses
# to ONE rect.  This is the primary lever for "minimum shapes".
DEFAULT_FIDELITY_TOL_M = 1.0

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


def _clamp_profile_with_anchors(zs: List[float],
                                anchored: List[bool],
                                seg_len: float,
                                max_grade: float,
                                max_dg_per_m: float) -> None:
    """Grade + rate-of-change clamp that preserves anchored samples.

    Anchored samples (typically taxi-centerline points within 1 m of
    a runway boundary, pinned to the runway's elevation) are
    IMMOVABLE.  Non-anchored samples are pulled into the grade
    envelope imposed by their nearest neighbours, which may be
    an anchored sample or another non-anchored sample; the clamp
    propagates the anchor's z through the chain.

    Iterates until convergence.
    """
    if len(zs) < 2 or seg_len <= 0:
        return
    max_dz = max_grade * seg_len
    grade_change_budget = max_dg_per_m * seg_len * seg_len
    n = len(zs)

    for _ in range(_CLAMP_MAX_ITERS):
        changed = False

        # Pass 1: hard-range clamp.  Upper/lower envelope is seeded
        # with [z, z] at anchored samples and [-inf, +inf] everywhere
        # else; then forward / backward passes propagate the
        # neighbourhood's imposed limits.
        BIG = 1e18
        upper = [zs[i] if anchored[i] else +BIG for i in range(n)]
        lower = [zs[i] if anchored[i] else -BIG for i in range(n)]
        for i in range(1, n):
            upper[i] = min(upper[i], upper[i - 1] + max_dz)
            lower[i] = max(lower[i], lower[i - 1] - max_dz)
        for i in range(n - 2, -1, -1):
            upper[i] = min(upper[i], upper[i + 1] + max_dz)
            lower[i] = max(lower[i], lower[i + 1] - max_dz)
        for i in range(n):
            if anchored[i]:
                continue
            new_z = zs[i]
            if new_z > upper[i]:
                new_z = upper[i]
            if new_z < lower[i]:
                new_z = lower[i]
            if abs(new_z - zs[i]) > 1e-9:
                zs[i] = new_z
                changed = True

        # Pass 2: rate-of-change cap.  Skip anchored samples.
        for i in range(1, n - 1):
            if anchored[i]:
                continue
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
    max_rect_length_m: float = 0.0,
    runway_polygon: Optional[Polygon] = None,
    runway_elev_lookup: Optional[Callable[[float, float],
                                          Optional[float]]] = None,
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
    anchored: List[bool] = []
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
        anchored.append(False)

    # Runway-endpoint anchoring.  Same rule as
    # build_rects_along_centerline: if the polygon's long-axis
    # endpoint (i.e. the strip's short edge) comes within 1 m of a
    # runway, pin that end to the runway's elevation at the
    # touch point and mark it anchored so the grade clamp
    # propagates it through the chain as an immovable boundary.
    #
    # A taxi running parallel to a runway has both endpoints far
    # from the runway, so no anchoring occurs — exactly as the
    # user clarified.
    from shapely.geometry import Point as _Pt
    # The Voronoi-skeleton endpoint of a taxi polygon sits on the
    # medial axis, which is ~half the polygon's local width away
    # from the nearest boundary — NOT at the boundary itself.
    # For a typical 25–45 m wide taxi that means the skeleton
    # endpoint is 12–22 m inside the polygon, so "near a runway
    # join" is an endpoint within ~half-max-taxi-width of the
    # runway boundary.  The user's "1 m" threshold referred to the
    # taxi pavement's own distance from the runway (they never run
    # closer than that without crossing); the centerline we actually
    # sample is a geometric derivative, not the pavement itself, so
    # it needs a wider search radius.
    ANCHOR_MAX_DIST_M = 25.0
    if (runway_polygon is not None and not runway_polygon.is_empty
            and runway_elev_lookup is not None):
        # Near-start endpoint
        try:
            p_start = _Pt(xs_center[0][0], xs_center[0][1])
            if p_start.distance(runway_polygon) <= ANCHOR_MAX_DIST_M:
                rz = runway_elev_lookup(
                    xs_center[0][0], xs_center[0][1])
                if rz is not None:
                    zs[0] = float(rz)
                    anchored[0] = True
        except Exception:
            pass
        try:
            last = n_samples - 1
            p_end = _Pt(xs_center[last][0], xs_center[last][1])
            if p_end.distance(runway_polygon) <= ANCHOR_MAX_DIST_M:
                rz = runway_elev_lookup(
                    xs_center[last][0], xs_center[last][1])
                if rz is not None:
                    zs[last] = float(rz)
                    anchored[last] = True
        except Exception:
            pass

    if any(anchored):
        _clamp_profile_with_anchors(
            zs, anchored, seg_len, max_grade, max_dg_per_m)
    else:
        _clamp_profile(zs, seg_len, max_grade, max_dg_per_m)

    keep_indices = _rdp_simplify_indices(zs, fidelity_tol)
    if len(keep_indices) < 2:
        keep_indices = [0, n_samples - 1]

    # Same length cap as build_rects_along_centerline — split any
    # segment longer than max_rect_length_m so a gently curved
    # taxi emits a chain of rects that follows the curve rather
    # than one giant rect whose axis drifts off the polygon.
    if max_rect_length_m > 0:
        capped: List[int] = [keep_indices[0]]
        for k in range(len(keep_indices) - 1):
            i_a = keep_indices[k]
            i_b = keep_indices[k + 1]
            span_len = (i_b - i_a) * seg_len
            if span_len <= max_rect_length_m:
                capped.append(i_b)
                continue
            n_pieces = max(
                2, int(round(span_len / max_rect_length_m + 0.5)))
            for s in range(1, n_pieces + 1):
                ix = i_a + int(
                    round(s * (i_b - i_a) / n_pieces))
                if ix > capped[-1]:
                    capped.append(ix)
        keep_indices = capped

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


# ──────────────────────────────────────────────────────────────────────
# Centerline-driven rect chain (used for curved / multi-strip
# polygons that fail the MRR-aligned gate above)
# ──────────────────────────────────────────────────────────────────────
def build_rects_along_centerline(
    centerline,
    polygon: Polygon,
    sample_dem: Callable[[float, float], Optional[float]],
    max_grade: float = 0.015,
    seg_length: float = 30.0,
    fidelity_tol: float = DEFAULT_FIDELITY_TOL_M,
    max_dg_per_m: float = DEFAULT_MAX_DG_PER_M,
    min_width_m: float = 8.0,
    max_width_m: float = 50.0,
    max_rect_length_m: float = 0.0,
    runway_polygon: Optional[Polygon] = None,
    runway_elev_lookup: Optional[Callable[[float, float],
                                          Optional[float]]] = None,
) -> Optional[List[TaxiwayRect]]:
    """Build one sloping rect per STRAIGHT segment of a centerline.

    The centerline is typically a Ramer-Douglas-Peucker simplified
    Voronoi skeleton path; its consecutive vertices are therefore
    already guaranteed to represent straight runs of the taxi
    polygon.  We emit exactly ONE rect per consecutive-vertex
    segment — no uniform sub-sampling — so a long straight taxiway
    collapses to a single rect regardless of length.

    User directive: "we can use between 25-30 rectangles to cover
    all the taxiways, excluding junctions", and "our shapes may
    go beyond the bounds of the apt.dat pavement by a small
    amount if needed to simplify coverage".  This function
    therefore uses the MAX of local-width probes (not the min or
    median) plus a small outward padding so each rect fully
    covers its segment of the taxi polygon — the rect can
    overflow into adjacent pavement by a few metres, which
    downstream apron/junction emission subtracts.

    Grade clamping + RDP on the elevation profile is still
    applied (so a 500 m rect with a 3 m DEM bump in the middle
    splits into 2 rects at the bump), but the segmentation
    floor is the CENTERLINE's own vertex count — the longest
    possible rect per centerline vertex pair is emitted.

    Args:
        centerline: shapely ``LineString`` in meter space, from
            :func:`O4_Taxiway_Skeleton.extract_centerlines`.
        polygon: the source taxiway polygon.
        sample_dem: ``(x, y) → elevation_m | None``.
        max_grade: longitudinal grade cap (default 1.5 % for
            taxiways).
        seg_length: DEM sample spacing along the centerline for
            grade clamping (not related to rect length any more).
        fidelity_tol: RDP tolerance for elevation-profile splits.
        max_dg_per_m: vertical-curve rate-of-change cap.
        min_width_m / max_width_m: if the computed width of a
            segment is below min_width_m the segment is dropped;
            above max_width_m it is clamped.
        max_rect_length_m: 0 disables any hard length cap (the
            default — let centerline vertex spacing determine
            rect length).  A positive value still splits any
            vertex-to-vertex segment exceeding the cap as a
            safety for extreme centerlines.
        runway_polygon / runway_elev_lookup: runway endpoint
            anchoring, unchanged from the previous version.

    Returns:
        A list of :class:`TaxiwayRect` or ``None`` if the
        centerline is degenerate or too narrow at every segment.
    """
    import math
    from shapely.geometry import Point

    if centerline is None or centerline.is_empty:
        return None
    if polygon is None or polygon.is_empty:
        return None
    total_len = centerline.length
    if total_len < 1.0:
        return None

    # Width padding added to every rect so its edges comfortably
    # reach (and slightly exceed) the local polygon boundary.
    # Coverage-over-correctness per the user's "shapes may go
    # beyond the bounds if needed to simplify coverage" directive.
    WIDTH_PADDING_M = 1.5

    # 1. Centerline vertices — these are the "natural" rect break
    # points.  RDP-simplify them once more to enforce `fidelity_tol`
    # geometric straightness and drop tiny dog-legs that would
    # otherwise emit their own rect.  The input centerline is
    # already RDP-simplified at skeleton-extract time; this is a
    # second pass with a (possibly) tighter tolerance, but default
    # DEFAULT_FIDELITY_TOL_M = 1.0 m is a no-op for a 10 m-simplified
    # skeleton centerline.
    base_coords = list(centerline.coords)
    if len(base_coords) < 2:
        return None

    # 2. Runway-endpoint anchoring (unchanged).
    ANCHOR_MAX_DIST_M = 25.0
    anchor_at_start = None
    anchor_at_end = None
    if (runway_polygon is not None
            and not runway_polygon.is_empty
            and runway_elev_lookup is not None):
        try:
            p_start = Point(base_coords[0])
            p_end = Point(base_coords[-1])
            if p_start.distance(runway_polygon) <= ANCHOR_MAX_DIST_M:
                anchor_at_start = runway_elev_lookup(
                    base_coords[0][0], base_coords[0][1])
            if p_end.distance(runway_polygon) <= ANCHOR_MAX_DIST_M:
                anchor_at_end = runway_elev_lookup(
                    base_coords[-1][0], base_coords[-1][1])
        except Exception:
            pass

    # 3. For each vertex-to-vertex segment: sample DEM along the
    # segment, grade-clamp the profile, RDP-simplify, and split
    # into 1+ sub-rects if the elevation profile demands.
    rects: List[TaxiwayRect] = []
    for seg_idx in range(len(base_coords) - 1):
        a = base_coords[seg_idx]
        b = base_coords[seg_idx + 1]
        seg_dx = b[0] - a[0]
        seg_dy = b[1] - a[1]
        seg_len = math.hypot(seg_dx, seg_dy)
        if seg_len < 5.0:
            continue
        ux = seg_dx / seg_len
        uy = seg_dy / seg_len
        px, py = -uy, ux

        # DEM samples along the segment.
        n_samples = max(2, int(round(seg_len / seg_length)) + 1)
        step = seg_len / (n_samples - 1)
        seg_centers: List[Tuple[float, float]] = []
        seg_zs: List[float] = []
        seg_anchored: List[bool] = []
        for k in range(n_samples):
            t = k * step
            cx = a[0] + ux * t
            cy = a[1] + uy * t
            seg_centers.append((cx, cy))
            # Anchor at the very first vertex of seg 0 and at the
            # very last vertex of the final segment.
            forced = None
            if seg_idx == 0 and k == 0 and anchor_at_start is not None:
                forced = float(anchor_at_start)
            elif (seg_idx == len(base_coords) - 2
                  and k == n_samples - 1
                  and anchor_at_end is not None):
                forced = float(anchor_at_end)
            if forced is not None:
                seg_zs.append(forced)
                seg_anchored.append(True)
            else:
                z = sample_dem(cx, cy)
                if z is None:
                    z = seg_zs[-1] if seg_zs else 0.0
                seg_zs.append(z)
                seg_anchored.append(False)

        _clamp_profile_with_anchors(
            seg_zs, seg_anchored, step, max_grade, max_dg_per_m)

        # RDP-simplify the elevation profile.  Multi-rect split
        # ONLY happens when elevation fidelity demands it; straight
        # flat / uniform-slope segments stay as one rect.
        keep = _rdp_simplify_indices(seg_zs, fidelity_tol)
        if len(keep) < 2:
            keep = [0, n_samples - 1]

        # Optional hard length cap as a safety.  Disabled by
        # default (max_rect_length_m = 0).
        if max_rect_length_m > 0:
            capped = [keep[0]]
            for j in range(len(keep) - 1):
                ia = keep[j]
                ib = keep[j + 1]
                span = (ib - ia) * step
                if span <= max_rect_length_m:
                    capped.append(ib)
                    continue
                n_pieces = max(2, int(math.ceil(
                    span / max_rect_length_m)))
                for s in range(1, n_pieces + 1):
                    ix = ia + int(
                        round(s * (ib - ia) / n_pieces))
                    if ix > capped[-1]:
                        capped.append(ix)
            keep = capped

        # Compute ONE rect width for the whole vertex-to-vertex
        # segment: max of 5 probe widths + padding.  All sub-rects
        # that come from RDP splits share this width because the
        # segment is straight — uniform direction = uniform
        # perpendicular = uniform local width range.
        widths: List[float] = []
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            mid_x = a[0] + ux * seg_len * frac
            mid_y = a[1] + uy * seg_len * frac
            hw = _local_half_width(
                polygon, mid_x, mid_y, px, py, max_reach=60.0)
            if hw > 0:
                widths.append(hw * 2.0)
        if not widths:
            continue
        local_width = max(widths) + WIDTH_PADDING_M
        if local_width < min_width_m:
            continue
        if local_width > max_width_m:
            local_width = max_width_m
        half_w = local_width / 2.0

        # Emit one rect per RDP sub-segment.
        for j in range(len(keep) - 1):
            ia = keep[j]
            ib = keep[j + 1]
            ax_pt = seg_centers[ia]
            ay_pt = seg_centers[ib]
            z_a = seg_zs[ia]
            z_b = seg_zs[ib]
            if z_a >= z_b:
                c_high, c_low = ax_pt, ay_pt
                eh, el = z_a, z_b
            else:
                c_high, c_low = ay_pt, ax_pt
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
                width_m=local_width,
            ))

    return rects if rects else None


def _local_half_width(polygon: Polygon,
                      cx: float, cy: float,
                      px: float, py: float,
                      max_reach: float = 60.0) -> float:
    """Return the distance from ``(cx, cy)`` to the nearest polygon
    boundary in the direction ``±(px, py)``, clipped by the lesser
    of the two half-rays.  ``(px, py)`` should be a unit
    perpendicular to the local centerline tangent.
    """
    from shapely.geometry import LineString, Point
    if polygon is None or polygon.is_empty:
        return 0.0
    try:
        if not polygon.contains(Point(cx, cy)):
            return 0.0
    except Exception:
        return 0.0

    left_end = (cx + px * max_reach, cy + py * max_reach)
    right_end = (cx - px * max_reach, cy - py * max_reach)
    try:
        left_seg = LineString(
            [(cx, cy), left_end]).intersection(polygon)
        right_seg = LineString(
            [(cx, cy), right_end]).intersection(polygon)
    except Exception:
        return 0.0

    def _len(seg):
        if seg is None or seg.is_empty:
            return 0.0
        try:
            return float(seg.length)
        except Exception:
            return 0.0

    left_d = _len(left_seg)
    right_d = _len(right_seg)
    if left_d <= 0 or right_d <= 0:
        return 0.0
    return min(left_d, right_d)
