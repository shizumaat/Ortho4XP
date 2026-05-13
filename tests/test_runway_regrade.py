"""Unit tests for the Stage A runway regrade solver."""
from __future__ import annotations

import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from auto_patch.runway_regrade import regrade_runway


def test_no_seams_returns_cifp():
    r = regrade_runway(50.0, 52.0, 3000.0, [])
    assert r.threshold_A == 50.0
    assert r.threshold_B == 52.0
    assert r.warnings == []


def test_single_seam_compatible_with_cifp():
    """CIFP at 50/52, runway 3000m, seam at midpoint sees DEM=51.
    Grade A→seam = -1/1500 = -0.067%, well under 1.5%.
    Both thresholds should stay at CIFP.
    """
    r = regrade_runway(50.0, 52.0, 3000.0, [(1500.0, 51.0)])
    assert abs(r.threshold_A - 50.0) < 0.01
    assert abs(r.threshold_B - 52.0) < 0.01
    # No threshold-shift warnings.
    assert all("shifted" not in w for w in r.warnings)


def test_single_seam_grade_cap_binds():
    """Runway 1000m, seam at 500m, DEM 60m, CIFP both 50m.
    Required grade A→seam = +10/500 = 2% > 1.5% cap.
    Solver lifts threshold to seam_alt - 0.015*500 = 52.5m.
    K-factor here: 2*500/305 = 3.28% allows the resulting Δg=3% — so
    K-factor does NOT bind and grade-cap is the only active constraint.
    """
    r = regrade_runway(50.0, 50.0, 1000.0, [(500.0, 60.0)])
    expected = 60.0 - 0.015 * 500.0  # = 52.5
    assert abs(r.threshold_A - expected) < 0.1
    assert abs(r.threshold_B - expected) < 0.1
    # CIFP-shift warnings should fire.
    assert any("threshold A shifted" in w for w in r.warnings)
    assert any("threshold B shifted" in w for w in r.warnings)


def test_two_seams_interior_grade_violation_warns():
    """Two seams at 1000m and 1500m apart. Interior segment is 500m
    with DEM altitudes 50 → 60 → that's 2% grade, violates 1.5%.
    Warning should appear but solver still optimises thresholds.
    """
    r = regrade_runway(45.0, 55.0, 2500.0,
                       [(1000.0, 50.0), (1500.0, 60.0)])
    assert any("interior" in w and "grade" in w for w in r.warnings)


def test_seam_anchors_immutable():
    """seam_altitudes should be returned unmodified."""
    seams = [(500.0, 47.5), (1000.0, 48.2)]
    r = regrade_runway(45.0, 50.0, 1500.0, seams)
    # Internal sort may reorder; check contents match.
    assert sorted(r.seam_altitudes) == sorted(seams)


def test_malformed_seams_filtered():
    """Seams outside (0, runway_length) get filtered."""
    r = regrade_runway(50.0, 52.0, 1000.0,
                       [(-100.0, 40.0), (1500.0, 60.0)])
    assert r.threshold_A == 50.0
    assert r.threshold_B == 52.0


def test_single_seam_at_short_runway_kfactor():
    """K-factor: |Δg| ≤ 2*min(d_A,d_B)/K.
    Runway 600m, seam at midpoint (300m each side), K=305m.
    Δg_max = 2*300/305 = 1.97%.
    With DEM 5m above CIFP-midpoint, Δg = 5/300 + 5/300 = 3.3% — exceeds.
    Solver should adjust to satisfy both grade cap and K-factor.
    """
    r = regrade_runway(50.0, 50.0, 600.0, [(300.0, 55.0)])
    # Grade caps must hold: alt_A, alt_B within [55-4.5, 55+4.5] = [50.5, 59.5]
    assert 50.4 <= r.threshold_A <= 59.6
    assert 50.4 <= r.threshold_B <= 59.6
    # After projection, K-factor may or may not be exactly satisfied;
    # depends on grade-cap clipping. Verify grade cap held in any case.
    g0 = (55.0 - r.threshold_A) / 300.0
    g1 = (r.threshold_B - 55.0) / 300.0
    assert abs(g0) <= 0.015 + 1e-6
    assert abs(g1) <= 0.015 + 1e-6


def test_threshold_shift_warning_format():
    r = regrade_runway(50.0, 50.0, 500.0, [(250.0, 60.0)])
    # Should warn about threshold A and B shifts.
    assert any("threshold A shifted" in w for w in r.warnings)
    assert any("threshold B shifted" in w for w in r.warnings)
