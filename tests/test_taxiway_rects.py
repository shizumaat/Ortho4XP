"""Unit tests for O4_Taxiway_Rects."""
import math

from shapely.affinity import rotate
from shapely.geometry import Polygon

import O4_Taxiway_Rects as TR


# ──────────────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────────────
def _strip(width: float, length: float, cx: float = 0.0, cy: float = 0.0
           ) -> Polygon:
    """Axis-aligned rectangle of given width (y-dim) and length (x-dim),
    centred at (cx, cy)."""
    dx = length / 2.0
    dy = width / 2.0
    return Polygon([
        (cx - dx, cy - dy),
        (cx + dx, cy - dy),
        (cx + dx, cy + dy),
        (cx - dx, cy + dy),
    ])


def _flat_dem(z: float):
    def _sampler(x, y):
        return z
    return _sampler


def _linear_x_dem(z0: float, slope: float):
    """DEM z = z0 + slope * x."""
    def _sampler(x, y):
        return z0 + slope * x
    return _sampler


# ──────────────────────────────────────────────────────────────────────
# _long_axis
# ──────────────────────────────────────────────────────────────────────
def test_long_axis_axis_aligned():
    poly = _strip(20, 200)
    result = TR._long_axis(poly)
    assert result is not None
    m_a, m_b, ux, uy, long_len, short_len = result
    assert math.isclose(long_len, 200.0, abs_tol=1e-6)
    assert math.isclose(short_len, 20.0, abs_tol=1e-6)
    # The long axis runs horizontally; unit vector is ±(1, 0).
    assert math.isclose(abs(ux), 1.0, abs_tol=1e-6)
    assert math.isclose(uy, 0.0, abs_tol=1e-6)


def test_long_axis_rotated():
    poly = rotate(_strip(20, 200), 30.0, origin=(0, 0))
    result = TR._long_axis(poly)
    assert result is not None
    m_a, m_b, ux, uy, long_len, short_len = result
    assert math.isclose(long_len, 200.0, abs_tol=1e-3)
    assert math.isclose(short_len, 20.0, abs_tol=1e-3)
    # Unit vector magnitude check
    assert math.isclose(ux * ux + uy * uy, 1.0, abs_tol=1e-6)


def test_long_axis_empty():
    assert TR._long_axis(Polygon()) is None
    assert TR._long_axis(None) is None


# ──────────────────────────────────────────────────────────────────────
# _clamp_profile
# ──────────────────────────────────────────────────────────────────────
def test_clamp_profile_longitudinal_cap():
    # 100 m segments, 1.5% grade → max dz per segment = 1.5 m.
    # Start with a 4 m jump, should be split symmetrically.
    zs = [0.0, 4.0, 4.0]
    TR._clamp_profile(zs, seg_len=100.0, max_grade=0.015,
                      max_dg_per_m=1.0 / 3000.0)
    assert zs[1] - zs[0] <= 1.5 + 1e-6
    assert zs[2] - zs[1] <= 1.5 + 1e-6


def test_clamp_profile_already_compliant():
    zs = [0.0, 1.0, 2.0, 3.0]
    before = list(zs)
    TR._clamp_profile(zs, seg_len=100.0, max_grade=0.015,
                      max_dg_per_m=1.0 / 3000.0)
    # Flat 1% slope — should be unchanged.
    for a, b in zip(before, zs):
        assert math.isclose(a, b, abs_tol=1e-6)


def test_clamp_profile_rate_of_change():
    # A sharp kink: 0,0,0,5,10,10 over 100 m steps.  max_dg_per_m =
    # 1/3000 so grade-change budget over 100 m = 100*100/3000 ≈ 3.33 m.
    # The kink samples (index 2 is 0, index 3 is 5, index 4 is 10)
    # start with a second derivative that violates the rule and should
    # smooth.
    zs = [0.0, 0.0, 0.0, 5.0, 10.0, 10.0]
    TR._clamp_profile(zs, seg_len=100.0, max_grade=0.10,  # relaxed
                      max_dg_per_m=1.0 / 3000.0)
    # After clamping, no triple should exceed the 3.33 m mid-deviation.
    for i in range(1, len(zs) - 1):
        expected = (zs[i - 1] + zs[i + 1]) / 2.0
        assert abs(zs[i] - expected) <= 3.34


# ──────────────────────────────────────────────────────────────────────
# _rdp_simplify_indices
# ──────────────────────────────────────────────────────────────────────
def test_rdp_linear_collapses_to_endpoints():
    zs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    keep = TR._rdp_simplify_indices(zs, tolerance=0.3)
    assert keep == [0, 5]


def test_rdp_kink_keeps_midpoint():
    zs = [0.0, 0.5, 1.0, 5.0, 9.0, 9.5, 10.0]
    keep = TR._rdp_simplify_indices(zs, tolerance=0.3)
    assert 0 in keep
    assert 6 in keep
    assert len(keep) >= 3   # at least one split


def test_rdp_noise_below_tolerance():
    # 0.2 m oscillation on a 10 m rise — all below 0.3 tol.
    zs = [0.0, 1.9, 4.1, 5.9, 8.1, 10.0]
    keep = TR._rdp_simplify_indices(zs, tolerance=0.3)
    assert keep == [0, 5]


def test_rdp_empty_and_tiny():
    assert TR._rdp_simplify_indices([], 0.3) == []
    assert TR._rdp_simplify_indices([1.0], 0.3) == [0]
    assert TR._rdp_simplify_indices([1.0, 2.0], 0.3) == [0, 1]


# ──────────────────────────────────────────────────────────────────────
# build_taxiway_rects — integration
# ──────────────────────────────────────────────────────────────────────
def test_build_rects_flat_taxiway_single_rect():
    # 20 m × 300 m, flat DEM at 10 m.  Expect exactly one rect.
    poly = _strip(20, 300)
    rects = TR.build_taxiway_rects(poly, _flat_dem(10.0))
    assert rects is not None
    assert len(rects) == 1
    r = rects[0]
    assert r.is_flat
    assert math.isclose(r.elev_high, 10.0, abs_tol=1e-6)
    assert math.isclose(r.elev_low, 10.0, abs_tol=1e-6)
    assert math.isclose(r.width_m, 20.0, abs_tol=1e-6)
    # Polygon area approximately matches the strip.
    assert math.isclose(r.polygon.area, 20.0 * 300.0, rel_tol=1e-6)


def test_build_rects_uniform_slope_single_rect():
    # 20 m × 300 m, DEM slope 1% (below max_grade 1.5%).
    poly = _strip(20, 300)
    # Shift polygon to start at x=0 so DEM z = slope*x works cleanly.
    poly = Polygon([(0, -10), (300, -10), (300, 10), (0, 10)])
    rects = TR.build_taxiway_rects(poly, _linear_x_dem(5.0, 0.01))
    assert rects is not None
    # A linear DEM should collapse to ONE rect under RDP.
    assert len(rects) == 1
    r = rects[0]
    assert not r.is_flat
    # Elevation range ≈ slope * length = 0.01 * 300 = 3.0 m.
    assert math.isclose(r.elev_high - r.elev_low, 3.0, abs_tol=0.05)


def test_build_rects_grade_clamped():
    # DEM slope 5% (way above 1.5% cap).  Output should enforce
    # the cap.
    poly = Polygon([(0, -10), (300, -10), (300, 10), (0, 10)])
    rects = TR.build_taxiway_rects(
        poly, _linear_x_dem(5.0, 0.05), max_grade=0.015)
    assert rects is not None
    # Whatever segmentation, the total rise must be ≤ 1.5% × 300 = 4.5m.
    total_rise = sum(r.elev_high - r.elev_low for r in rects
                     if not r.is_flat)
    assert total_rise <= 4.5 + 0.1


def test_build_rects_irregular_polygon_fallback():
    # A star-ish polygon that fits poorly inside its MRR.
    # MRR ≈ 200 × 200 = 40000; polygon is a plus-sign with area ~ 40*200 = 8000.
    poly = Polygon([
        (-20, -100), (20, -100), (20, -20),
        (100, -20), (100, 20), (20, 20),
        (20, 100), (-20, 100), (-20, 20),
        (-100, 20), (-100, -20), (-20, -20),
    ])
    result = TR.build_taxiway_rects(poly, _flat_dem(10.0))
    assert result is None   # fit ratio too low, caller falls back


def test_build_rects_empty_polygon():
    assert TR.build_taxiway_rects(Polygon(), _flat_dem(10.0)) is None


def test_build_rects_corners_have_high_end_first():
    # High end to the right (elevation 10), low end to the left (elev 0).
    poly = Polygon([(0, -10), (300, -10), (300, 10), (0, 10)])
    rects = TR.build_taxiway_rects(
        poly, _linear_x_dem(0.0, 0.01), max_grade=0.015)
    assert rects is not None
    r = rects[0]
    # elev_high end is at x=300 (the right side).
    assert math.isclose(r.center_high[0], 300.0, abs_tol=1.0)
    assert math.isclose(r.center_low[0], 0.0, abs_tol=1.0)
    # Corners 0-1 are at the high end, so their x is ≈ 300.
    assert math.isclose(r.corners_m[0][0], 300.0, abs_tol=1.0)
    assert math.isclose(r.corners_m[1][0], 300.0, abs_tol=1.0)
    assert math.isclose(r.corners_m[2][0], 0.0, abs_tol=1.0)
    assert math.isclose(r.corners_m[3][0], 0.0, abs_tol=1.0)


def test_build_rects_rotated_polygon():
    # Rotated 45°, same taxiway.
    poly = rotate(Polygon([(0, -10), (300, -10), (300, 10), (0, 10)]),
                  45.0, origin=(150, 0))
    rects = TR.build_taxiway_rects(poly, _flat_dem(10.0))
    assert rects is not None
    assert len(rects) == 1
    r = rects[0]
    # Width preserved.
    assert math.isclose(r.width_m, 20.0, abs_tol=1e-3)
    # Polygon area preserved.
    assert math.isclose(r.polygon.area, 6000.0, rel_tol=1e-3)


def test_build_rects_emits_rect_chain_when_profile_bent():
    # DEM with a dip at the middle: 0 rising to 0, with a -5 at center.
    def _dem(x, y):
        if 140 <= x <= 160:
            return -2.0
        return 0.0
    poly = Polygon([(0, -10), (300, -10), (300, 10), (0, 10)])
    rects = TR.build_taxiway_rects(poly, _dem, max_grade=0.015,
                                    fidelity_tol=0.3)
    assert rects is not None
    # After clamping, the dip is smoothed; the chain has 1-3 segments.
    assert 1 <= len(rects) <= 4
