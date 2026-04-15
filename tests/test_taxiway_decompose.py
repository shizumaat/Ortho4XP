"""Unit tests for O4_Taxiway_Decompose."""
import math

from shapely.affinity import rotate, translate
from shapely.geometry import Polygon
from shapely.ops import unary_union

import O4_Taxiway_Decompose as TD


def _strip(width, length, cx=0.0, cy=0.0):
    dx, dy = length / 2.0, width / 2.0
    return Polygon([
        (cx - dx, cy - dy), (cx + dx, cy - dy),
        (cx + dx, cy + dy), (cx - dx, cy + dy),
    ])


# ──────────────────────────────────────────────────────────────────────
# Simple strip — no decomposition
# ──────────────────────────────────────────────────────────────────────
def test_simple_strip_returns_unchanged():
    poly = _strip(25, 400)  # 25 × 400 — classic taxiway strip
    result = TD.decompose_multi_taxiway(poly)
    assert result.used_decomposition is False
    assert len(result.strip_polygons) == 1
    assert not result.junction_polygons
    assert math.isclose(
        result.strip_polygons[0].area, 25 * 400, rel_tol=1e-6)


def test_very_narrow_strip_returns_unchanged():
    # 20 m width — well under 2 * 15 m half-width; morphological
    # opening erases it completely, so caller should rect-chain it.
    poly = _strip(20, 300)
    result = TD.decompose_multi_taxiway(poly)
    assert result.used_decomposition is False
    assert len(result.strip_polygons) == 1


def test_empty_polygon():
    result = TD.decompose_multi_taxiway(Polygon())
    assert result.used_decomposition is False
    assert not result.strip_polygons
    assert not result.junction_polygons


def test_none_polygon():
    result = TD.decompose_multi_taxiway(None)
    assert not result.strip_polygons


# ──────────────────────────────────────────────────────────────────────
# T-junction — two strips meeting at a wide hub
# ──────────────────────────────────────────────────────────────────────
def test_t_junction_splits_into_two_branches_and_one_hub():
    # Horizontal strip 400 × 30 at y=0 joined by a vertical stub
    # 30 × 200 at x=0 extending upward.  The overlap region around
    # the origin is ~60 × 30 = "wide hub" after opening.
    horiz = _strip(30, 400)
    vert = Polygon([(-15, 0), (15, 0), (15, 200), (-15, 200)])
    t_poly = unary_union([horiz, vert])

    result = TD.decompose_multi_taxiway(t_poly, half_width_m=15.0)
    assert result.used_decomposition is True
    # Branches + hub together cover the whole polygon within
    # floating-point slop.
    all_covered = unary_union(
        result.strip_polygons + result.junction_polygons)
    assert math.isclose(all_covered.area, t_poly.area, rel_tol=0.02)


def test_cross_junction_splits_into_four_branches():
    # Plus-sign: horizontal 30 × 400 and vertical 30 × 400.
    horiz = _strip(30, 400)
    vert = Polygon([(-15, -200), (15, -200), (15, 200), (-15, 200)])
    plus = unary_union([horiz, vert])

    result = TD.decompose_multi_taxiway(plus, half_width_m=15.0)
    assert result.used_decomposition is True
    # The central square is the hub; the 4 arms are strips.
    assert len(result.strip_polygons) == 4
    assert len(result.junction_polygons) >= 1
    all_covered = unary_union(
        result.strip_polygons + result.junction_polygons)
    assert math.isclose(all_covered.area, plus.area, rel_tol=0.02)


def test_decomposition_is_disjoint_minus_tiny_overlap():
    # Branches are 1 m-extended into the hub for grade-join
    # continuity; the extension should be small compared to branch
    # area.
    horiz = _strip(30, 400)
    vert = Polygon([(-15, -200), (15, -200), (15, 200), (-15, 200)])
    plus = unary_union([horiz, vert])

    result = TD.decompose_multi_taxiway(plus, half_width_m=15.0)
    branch_sum = sum(p.area for p in result.strip_polygons)
    hub_sum = sum(p.area for p in result.junction_polygons)
    # Total is ≈ plus.area, the small overlap of 1m extensions is
    # absorbed into hub_polys via polygon.difference(branches).
    assert math.isclose(
        branch_sum + hub_sum, plus.area, rel_tol=0.01)


# ──────────────────────────────────────────────────────────────────────
# Blob with NO strips — apron-like
# ──────────────────────────────────────────────────────────────────────
def test_large_blob_with_no_branches():
    # 300 × 300 square — wider than 2*15 in every direction, so
    # there are no "strip branches".  Whole polygon is hub.
    poly = _strip(300, 300)
    result = TD.decompose_multi_taxiway(poly, half_width_m=15.0)
    assert result.used_decomposition is True
    # No branches, or only tiny branches below the area threshold.
    assert sum(p.area for p in result.strip_polygons) < 2000.0
    # Hub covers essentially everything.
    assert sum(
        p.area for p in result.junction_polygons) > poly.area * 0.95
