"""Unit tests for auto_patch.pavement.runway_geometry.

Pure-function geometric helpers — no DEM / shapely / I/O dependencies
beyond reading apt.dat for ``parse_aptdat_runway_widths``.  Tests
hit the math against synthetic CIFP-style inputs.
"""
import math
import os
import tempfile

import pytest

from auto_patch.pavement.runway_geometry import (
    DEFAULT_RUNWAY_WIDTH,
    DEG_TO_M,
    extend_point,
    get_reciprocal,
    pair_runways,
    parse_aptdat_runway_widths,
    runway_corners,
)


# ──────────────────────────────────────────────────────────────────────
# get_reciprocal
# ──────────────────────────────────────────────────────────────────────
def test_get_reciprocal_basic():
    """Heading +18 (mod 36) gives the reciprocal — RW09 → RW27."""
    assert get_reciprocal("RW09") == "RW27"
    assert get_reciprocal("RW27") == "RW09"
    assert get_reciprocal("RW01") == "RW19"
    assert get_reciprocal("RW36") == "RW18"


def test_get_reciprocal_with_lr_suffix():
    """L/R suffixes flip across the centerline; C stays C."""
    assert get_reciprocal("RW16L") == "RW34R"
    assert get_reciprocal("RW16R") == "RW34L"
    assert get_reciprocal("RW16C") == "RW34C"
    # And the reverse direction.
    assert get_reciprocal("RW34R") == "RW16L"


def test_get_reciprocal_wraparound():
    """Headings beyond 36 wrap correctly: RW20 → RW02 (not RW38)."""
    assert get_reciprocal("RW20") == "RW02"
    assert get_reciprocal("RW19") == "RW01"


def test_get_reciprocal_invalid_designator():
    """Malformed designators return None instead of raising."""
    assert get_reciprocal("not-a-runway") is None
    assert get_reciprocal("") is None
    assert get_reciprocal("RW1") is None  # need 2 digits


# ──────────────────────────────────────────────────────────────────────
# pair_runways
# ──────────────────────────────────────────────────────────────────────
def test_pair_runways_simple_pair():
    """A complete pair gets matched exactly once, with the alpha-
    earlier designator first."""
    runways = {
        "RW09": {"lat": 0.0, "lon": 0.0, "elevation_m": 100.0},
        "RW27": {"lat": 0.0, "lon": 0.027, "elevation_m": 105.0},
    }
    pairs = pair_runways(runways)
    assert len(pairs) == 1
    desig_a, data_a, desig_b, data_b = pairs[0]
    assert {desig_a, desig_b} == {"RW09", "RW27"}
    assert data_a["elevation_m"] in {100.0, 105.0}
    assert data_b["elevation_m"] in {100.0, 105.0}


def test_pair_runways_unpaired_threshold():
    """When only one end of a runway has CIFP data (common at minor
    airports), the unpaired threshold is returned with desig_b/data_b
    set to None so callers can branch."""
    runways = {
        "RW16L": {"lat": 1.0, "lon": 1.0, "elevation_m": 200.0},
    }
    pairs = pair_runways(runways)
    assert len(pairs) == 1
    desig_a, data_a, desig_b, data_b = pairs[0]
    assert desig_a == "RW16L"
    assert data_a["elevation_m"] == 200.0
    assert desig_b is None
    assert data_b is None


def test_pair_runways_multiple_runways():
    """Multiple parallel runways at one airport — each L/R/C variant
    pairs only with its own reciprocal."""
    runways = {
        "RW16L": {"lat": 0.0, "lon": 0.0, "elevation_m": 50.0},
        "RW16R": {"lat": 0.0, "lon": 0.001, "elevation_m": 50.0},
        "RW34L": {"lat": 0.027, "lon": 0.001, "elevation_m": 55.0},
        "RW34R": {"lat": 0.027, "lon": 0.0, "elevation_m": 55.0},
    }
    pairs = pair_runways(runways)
    assert len(pairs) == 2
    pair_sets = [{a, b} for a, _, b, _ in pairs]
    assert {"RW16L", "RW34R"} in pair_sets
    assert {"RW16R", "RW34L"} in pair_sets


def test_pair_runways_no_double_counting():
    """A paired threshold appears in exactly one pair, never two."""
    runways = {
        "RW09": {"lat": 0.0, "lon": 0.0},
        "RW27": {"lat": 0.0, "lon": 0.027},
    }
    pairs = pair_runways(runways)
    all_designators = []
    for desig_a, _, desig_b, _ in pairs:
        all_designators.append(desig_a)
        if desig_b is not None:
            all_designators.append(desig_b)
    assert len(all_designators) == len(set(all_designators))


# ──────────────────────────────────────────────────────────────────────
# runway_corners
# ──────────────────────────────────────────────────────────────────────
def test_runway_corners_returns_4_points():
    """Always returns exactly 4 (lat, lon) tuples for a valid runway."""
    corners = runway_corners(0.0, 0.0, 0.027, 0.0, width_m=45.0)
    assert corners is not None
    assert len(corners) == 4
    for c in corners:
        assert len(c) == 2  # (lat, lon)


def test_runway_corners_degenerate_returns_none():
    """A zero-length runway (both endpoints identical) returns None."""
    assert runway_corners(0.0, 0.0, 0.0, 0.0, width_m=45.0) is None


def test_runway_corners_corner_order_high_then_low():
    """Per the patch-altitude convention:
        c0 = high-side at (lat1, lon1)
        c1 = low-side at (lat2, lon2)
        c2 = low-side at (lat2, lon2) (other side)
        c3 = high-side at (lat1, lon1) (other side)

    So c0 and c3 should be roughly at lat1; c1 and c2 at lat2.
    """
    lat1, lon1 = 0.0, 0.0
    lat2, lon2 = 0.027, 0.0  # ~3 km north
    corners = runway_corners(lat1, lon1, lat2, lon2, width_m=45.0)
    c0, c1, c2, c3 = corners
    # c0/c3 close to lat1; c1/c2 close to lat2.
    assert abs(c0[0] - lat1) < 0.001
    assert abs(c3[0] - lat1) < 0.001
    assert abs(c1[0] - lat2) < 0.001
    assert abs(c2[0] - lat2) < 0.001


def test_runway_corners_width_perpendicular_to_axis():
    """The two corners on the same runway end (c0/c3 or c1/c2) must
    be exactly width_m apart."""
    width_m = 60.0
    corners = runway_corners(0.0, 0.0, 0.027, 0.0, width_m=width_m)
    c0, c1, c2, c3 = corners
    # c0 and c3 are at the same end (both at lat1=0); their separation
    # should equal width_m in meters.
    cos_lat = math.cos(0.0)
    dlat_m = (c0[0] - c3[0]) * DEG_TO_M
    dlon_m = (c0[1] - c3[1]) * cos_lat * DEG_TO_M
    sep_m = math.hypot(dlat_m, dlon_m)
    assert abs(sep_m - width_m) < 0.5


# ──────────────────────────────────────────────────────────────────────
# extend_point
# ──────────────────────────────────────────────────────────────────────
def test_extend_point_basic():
    """Extending 100 m past (lat_to, lon_to) along the from→to axis
    should land 100 m further in that direction."""
    # 1° lat ≈ 111120 m, so going from (0, 0) → (0.001, 0) is ~111 m.
    extended = extend_point(0.0, 0.0, 0.001, 0.0, distance_m=100.0)
    # Result should be further north than (0.001, 0).
    assert extended[0] > 0.001
    # Distance from (0.001, 0) to extended should be ~100 m.
    dlat_m = (extended[0] - 0.001) * DEG_TO_M
    dlon_m = (extended[1] - 0.0) * math.cos(0.001 * math.pi / 180) * DEG_TO_M
    sep_m = math.hypot(dlat_m, dlon_m)
    assert abs(sep_m - 100.0) < 1.0


def test_extend_point_zero_length_input_returns_to_point():
    """When the from→to vector is degenerate (< 1 m), no direction
    can be inferred — return (lat_to, lon_to) unchanged."""
    result = extend_point(0.0, 0.0, 0.0, 0.0, distance_m=100.0)
    assert result == (0.0, 0.0)


def test_extend_point_direction_preserved():
    """Extension always moves in the from→to direction, never back."""
    # West-bound runway: extending past should move further west.
    extended = extend_point(0.0, 0.0, 0.0, -0.001, distance_m=50.0)
    assert extended[1] < -0.001  # further west


# ──────────────────────────────────────────────────────────────────────
# parse_aptdat_runway_widths
# ──────────────────────────────────────────────────────────────────────
APT_DAT_FIXTURE = """\
I
1100 Generated by WorldEditor

1     12 0 0 KFAKE Fake Airport
100 45.72 1 0 0.25 1 2 1 16L 33.94195556 -118.40160833 0 0 2 0 0 0 34R 33.96091944 -118.41877778 0 0 2 0 0 0
100 60.00 1 0 0.25 1 2 1 09  33.95000000 -118.42000000 0 0 2 0 0 0 27  33.95000000 -118.40000000 0 0 2 0 0 0
99

1     50 0 0 KOTHER Other Airport
100 30.00 1 0 0.25 1 2 1 04  40.00 -73.00 0 0 2 0 0 0 22  40.05 -73.02 0 0 2 0 0 0
99
"""


def test_parse_aptdat_runway_widths_basic():
    """Reads row-100 widths from the named ICAO's airport block."""
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False) as f:
        f.write(APT_DAT_FIXTURE)
        path = f.name
    try:
        widths = parse_aptdat_runway_widths(path, "KFAKE")
        # Designators normalised to RWxx form.
        assert widths["RW16L"] == 45.72
        assert widths["RW34R"] == 45.72
        assert widths["RW09"] == 60.00
        assert widths["RW27"] == 60.00
        # No pollution from KOTHER's runways.
        assert "RW04" not in widths
        assert "RW22" not in widths
    finally:
        os.unlink(path)


def test_parse_aptdat_runway_widths_other_icao():
    """Selects the right airport block when multiple are present."""
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False) as f:
        f.write(APT_DAT_FIXTURE)
        path = f.name
    try:
        widths = parse_aptdat_runway_widths(path, "KOTHER")
        assert widths["RW04"] == 30.00
        assert widths["RW22"] == 30.00
        assert "RW16L" not in widths
    finally:
        os.unlink(path)


def test_parse_aptdat_runway_widths_missing_file_returns_empty():
    """Non-existent path returns empty dict instead of raising."""
    widths = parse_aptdat_runway_widths(
        "/does/not/exist/apt.dat", "KFAKE")
    assert widths == {}


def test_parse_aptdat_runway_widths_unknown_icao_returns_empty():
    """ICAO not found in file returns empty dict."""
    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False) as f:
        f.write(APT_DAT_FIXTURE)
        path = f.name
    try:
        widths = parse_aptdat_runway_widths(path, "ZZZZ")
        assert widths == {}
    finally:
        os.unlink(path)


def test_parse_aptdat_runway_widths_none_path_returns_empty():
    """None path is also tolerated (don't blow up at config-loading
    time when no apt.dat is configured)."""
    widths = parse_aptdat_runway_widths(None, "KFAKE")
    assert widths == {}


# ──────────────────────────────────────────────────────────────────────
# Module constants
# ──────────────────────────────────────────────────────────────────────
def test_default_runway_width():
    """Default width matches the FAA major-runway typical (45 m)."""
    assert DEFAULT_RUNWAY_WIDTH == 45.0


def test_deg_to_m_constant():
    """One degree of latitude is ~111.12 km."""
    assert 111000 < DEG_TO_M < 112000
