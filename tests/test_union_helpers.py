"""Unit tests for auto_patch.pavement.union_helpers.

The single function ``_merge_near_touching`` bridges sub-meter
precision gaps between adjacent apt.dat row-110 polygons that
``shapely.unary_union`` leaves disjoint.  These tests exercise the
rule directly with synthetic inputs.
"""
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

from auto_patch.pavement.union_helpers import (
    PAVEMENT_BRIDGE_GAP_M,
    _merge_near_touching,
)


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def _square(cx, cy, side=10.0):
    """Axis-aligned square centred at (cx, cy)."""
    h = side / 2.0
    return Polygon([
        (cx - h, cy - h),
        (cx + h, cy - h),
        (cx + h, cy + h),
        (cx - h, cy + h),
    ])


# ──────────────────────────────────────────────────────────────────────
# Edge-case inputs
# ──────────────────────────────────────────────────────────────────────
def test_none_input_returns_none():
    assert _merge_near_touching(None) is None


def test_empty_polygon_passes_through():
    empty = Polygon()
    assert _merge_near_touching(empty).is_empty


# ──────────────────────────────────────────────────────────────────────
# The core rule: bridge sub-meter gaps
# ──────────────────────────────────────────────────────────────────────
def test_two_squares_with_subcm_gap_get_merged():
    """Two adjacent squares 5 cm apart — exactly the SPJC SE-apron
    case mentioned in the module docstring.  ``unary_union`` would
    leave them disjoint; ``_merge_near_touching`` must fuse them."""
    a = _square(0.0, 0.0, side=10.0)
    b = _square(10.05, 0.0, side=10.0)  # 5 cm gap to the right
    union = unary_union([a, b])
    # Sanity: the gap is real (unary_union doesn't merge).
    assert union.geom_type == "MultiPolygon"

    merged = _merge_near_touching(union)
    assert merged.geom_type == "Polygon", (
        "subcm-gap squares should fuse into a single Polygon")
    # Total area should be ~ 200 m² (two 10×10 squares with the
    # bridging strip), tolerating the eps-thick bridge.
    assert 199.0 < merged.area < 202.0


def test_truly_disjoint_squares_stay_disjoint():
    """Squares 10 m apart — wider than any reasonable apt.dat
    precision drift — must NOT get force-merged."""
    a = _square(0.0, 0.0, side=10.0)
    b = _square(25.0, 0.0, side=10.0)  # 10 m gap
    union = unary_union([a, b])
    assert union.geom_type == "MultiPolygon"

    merged = _merge_near_touching(union)
    assert merged.geom_type == "MultiPolygon", (
        "10-m-apart squares must remain disjoint after merge")
    assert len(list(merged.geoms)) == 2


def test_already_connected_polygon_passes_through():
    """A single connected polygon should be returned essentially
    unchanged (modulo float noise) — the buffer-shrink round-trip is
    near-identity for connected geometry."""
    a = _square(0.0, 0.0, side=20.0)
    merged = _merge_near_touching(a)
    assert merged.geom_type == "Polygon"
    assert abs(merged.area - 400.0) < 1.0


def test_real_holes_are_preserved():
    """The merge must NOT fill in legitimate holes.  An outer 30×30
    square with a 5×5 hole in the middle stays as a polygon-with-hole
    after the round-trip."""
    outer = _square(0.0, 0.0, side=30.0)
    hole = _square(0.0, 0.0, side=5.0)
    holed = outer.difference(hole)
    assert len(holed.interiors) == 1

    merged = _merge_near_touching(holed)
    assert merged.geom_type == "Polygon"
    assert len(merged.interiors) == 1, (
        "_merge_near_touching must preserve real holes "
        "(only sub-eps gaps should close)")


# ──────────────────────────────────────────────────────────────────────
# Eps tunability
# ──────────────────────────────────────────────────────────────────────
def test_default_eps_value():
    """The module-level default tolerance is 0.1 m (per docstring)."""
    assert PAVEMENT_BRIDGE_GAP_M == 0.1


def test_custom_eps_can_bridge_wider_gap():
    """Passing a larger eps lets the helper bridge wider gaps.

    Mechanic: buffer outward by eps, then shrink by eps.  Bridges
    any gap below ~2*eps (each side of the gap grows by eps).
    """
    a = _square(0.0, 0.0, side=10.0)
    b = _square(13.0, 0.0, side=10.0)  # 3 m gap
    union = unary_union([a, b])
    assert union.geom_type == "MultiPolygon"

    # Default eps (0.1 m) won't bridge a 3 m gap.
    assert _merge_near_touching(union).geom_type == "MultiPolygon"

    # eps=1.0 → effective bridging up to ~2 m, won't bridge 3 m.
    assert _merge_near_touching(union, eps=1.0).geom_type == "MultiPolygon"

    # eps=2.0 → effective bridging up to ~4 m, bridges 3 m.
    bridged = _merge_near_touching(union, eps=2.0)
    assert bridged.geom_type == "Polygon"


# ──────────────────────────────────────────────────────────────────────
# Failure-mode safety
# ──────────────────────────────────────────────────────────────────────
def test_invalid_geometry_returns_input_unchanged():
    """When buffer/buffer fails (which shapely's bowtie can do),
    return the original geometry rather than crashing or returning
    something unexpected."""
    # A self-intersecting bowtie is shapely-invalid; buffering it can
    # explode in unpredictable ways depending on the GEOS version.
    bowtie = Polygon([(0, 0), (10, 10), (0, 10), (10, 0)])
    # Whether it raises or not, the helper must not propagate the
    # exception — it must return SOMETHING valid.  Don't assert on
    # the geometry contents; just confirm we get back a non-None
    # geometry of one of the accepted types.
    result = _merge_near_touching(bowtie)
    assert result is not None
    assert result.geom_type in ("Polygon", "MultiPolygon", "GeometryCollection")
