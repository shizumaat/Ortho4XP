"""Unit tests for ``O4_Surface_Mesh``.

These cover the small helpers (``fit_plane``, ``plane_grade``) and
the full ``adaptive_triangulate`` pipeline against synthetic DEM
samplers and known-good polygons.

Each test exercises a specific invariant from the module docstring:

* a perfectly flat region triangulates to the polygon's own corner
  count (2 triangles for a square, etc.) and preserves elevations,
* a uniformly sloped region under the grade limit triangulates to
  the same minimum count,
* an undulating region adds interior points only where DEM departs
  from the current plane,
* anchored vertices are NEVER moved by the grade clamp,
* an over-steep DEM is clamped to the configured ``max_grade``,
* a polygon whose anchors demand an impossible elevation range
  results in an unresolvable violation that the grade clamp leaves
  in place (the caller's reconciliation step is expected to fix it
  by adjusting the contributing building elevations).
"""
from math import sqrt
import pytest
from shapely.geometry import Polygon

import O4_Surface_Mesh as SM


# ──────────────────────────────────────────────────────────────────────
# Test fixtures
# ──────────────────────────────────────────────────────────────────────
SQUARE_100 = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
L_SHAPE = Polygon([
    (0, 0), (100, 0), (100, 60),
    (40, 60), (40, 100), (0, 100),
])


def constant_dem(z):
    """Return a DEM sampler that always reports `z`."""
    return lambda x, y: z


def linear_dem(slope_x, slope_y, intercept):
    """Return a DEM sampler with a linear surface."""
    return lambda x, y: intercept + slope_x * x + slope_y * y


def bumpy_dem(base, bump_height, bump_cx, bump_cy, bump_radius):
    """Return a DEM sampler with a single radial bump."""
    def _f(x, y):
        d = sqrt((x - bump_cx) ** 2 + (y - bump_cy) ** 2)
        return base + max(0.0, bump_height * (1.0 - d / bump_radius))
    return _f


# ──────────────────────────────────────────────────────────────────────
# fit_plane / plane_grade
# ──────────────────────────────────────────────────────────────────────
class TestFitPlane:
    def test_flat_plane(self):
        plane = SM.fit_plane((0, 0, 5), (10, 0, 5), (0, 10, 5))
        assert plane is not None
        a, b, c = plane
        assert abs(a) < 1e-9
        assert abs(b) < 1e-9
        assert abs(c - 5) < 1e-9
        assert SM.plane_grade(plane) == pytest.approx(0.0, abs=1e-9)

    def test_sloped_plane(self):
        # z = 0.01 * x + 0 * y + 10
        plane = SM.fit_plane((0, 0, 10), (100, 0, 11), (0, 100, 10))
        assert plane is not None
        a, b, c = plane
        assert a == pytest.approx(0.01, abs=1e-9)
        assert b == pytest.approx(0.0, abs=1e-9)
        assert c == pytest.approx(10.0, abs=1e-9)
        assert SM.plane_grade(plane) == pytest.approx(0.01, abs=1e-9)

    def test_compound_slope(self):
        # z = 0.01x + 0.02y + 5 → grade = sqrt(0.0001 + 0.0004) = 0.0224
        plane = SM.fit_plane((0, 0, 5), (100, 0, 6), (0, 100, 7))
        assert plane is not None
        assert plane[0] == pytest.approx(0.01)
        assert plane[1] == pytest.approx(0.02)
        assert SM.plane_grade(plane) == pytest.approx(
            sqrt(0.0001 + 0.0004), abs=1e-9)

    def test_colinear_returns_none(self):
        # All three points have the same y → cross product z is zero.
        plane = SM.fit_plane((0, 0, 0), (10, 0, 1), (20, 0, 2))
        assert plane is None


# ──────────────────────────────────────────────────────────────────────
# adaptive_triangulate — minimum shape count
# ──────────────────────────────────────────────────────────────────────
class TestMinimumShapes:
    def test_flat_square(self):
        """A flat square apron over a uniform DEM should triangulate
        to exactly 2 triangles, both at the DEM elevation."""
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], constant_dem(25.0), max_grade=0.010)
        assert len(tris) == 2
        for t in tris:
            for v in t:
                assert v[2] == pytest.approx(25.0, abs=1e-6)

    def test_compliant_uniform_slope(self):
        """A 0.5 % slope (well under 1 % limit) should triangulate
        to 2 triangles with no interior refinement."""
        dem = linear_dem(0.005, 0.0, 25.0)
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], dem, max_grade=0.010)
        assert len(tris) == 2
        zs = [v[2] for t in tris for v in t]
        assert min(zs) == pytest.approx(25.0, abs=0.01)
        assert max(zs) == pytest.approx(25.5, abs=0.01)

    def test_l_shape_minimum(self):
        """An L-shape has 6 exterior vertices.  A flat L-shape should
        triangulate to 4 triangles (n - 2 for a simple polygon)."""
        tris = SM.adaptive_triangulate(
            L_SHAPE, [], constant_dem(20.0), max_grade=0.010)
        # Exact count depends on Delaunay choice but should be small.
        assert 4 <= len(tris) <= 6
        for t in tris:
            for v in t:
                assert v[2] == pytest.approx(20.0, abs=1e-6)


# ──────────────────────────────────────────────────────────────────────
# adaptive_triangulate — fidelity refinement
# ──────────────────────────────────────────────────────────────────────
class TestFidelity:
    def test_bump_is_captured(self):
        """A 1.5 m bump in the middle of a flat square should cause
        at least one interior refinement insert."""
        dem = bumpy_dem(25.0, 1.5, 50, 50, 30)
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], dem, max_grade=0.010,
            fidelity_tol=1.0)
        # Must add at least one interior point to capture the bump.
        assert len(tris) > 2
        # Some vertex should reflect the bump (within tolerance).
        zs = [v[2] for t in tris for v in t]
        assert max(zs) > 25.5

    def test_fidelity_tolerance_respected(self):
        """A 0.3 m bump under fidelity_tol = 1.0 m should NOT trigger
        refinement — the flat-corners triangulation is already close
        enough."""
        dem = bumpy_dem(25.0, 0.3, 50, 50, 30)
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], dem, max_grade=0.010,
            fidelity_tol=1.0)
        assert len(tris) == 2

    def test_max_extra_points_caps_runaway(self):
        """A noisy DEM with huge variance should still respect
        max_extra_points so we don't get a runaway."""
        import random
        rng = random.Random(0)
        def noisy(x, y):
            return 25.0 + rng.uniform(-2.0, 2.0)
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], noisy, max_grade=0.010,
            fidelity_tol=0.5, max_extra_points=10)
        # We started with 4 polygon corners and were allowed 10
        # interior inserts, so at most 14 vertices = ~24 triangles
        # by Delaunay's vertex-count bound.
        assert len(tris) <= 30


# ──────────────────────────────────────────────────────────────────────
# adaptive_triangulate — anchors
# ──────────────────────────────────────────────────────────────────────
class TestAnchors:
    def test_anchored_vertex_is_pinned(self):
        """An anchor at one corner must be honoured exactly in the
        emitted triangulation, even if DEM disagrees and grade
        clamping happens."""
        anchors = [(0, 0, 30.0)]   # one corner anchored 5 m above DEM
        tris = SM.adaptive_triangulate(
            SQUARE_100, anchors, constant_dem(25.0),
            max_grade=0.010)
        # Find the (0, 0) vertex among the triangles.
        found = False
        for t in tris:
            for v in t:
                if abs(v[0]) < 0.5 and abs(v[1]) < 0.5:
                    assert v[2] == pytest.approx(30.0, abs=0.01)
                    found = True
        assert found

    def test_two_anchored_corners_with_compliant_slope(self):
        """Two corners anchored at heights that fit within 1 % grade
        of each other should produce a clean tilt."""
        # Both left corners at 30, both right at 30.5 → 0.5% across 100 m
        anchors = [(0, 0, 30.0), (0, 100, 30.0),
                   (100, 0, 30.5), (100, 100, 30.5)]
        tris = SM.adaptive_triangulate(
            SQUARE_100, anchors, constant_dem(25.0),
            max_grade=0.010)
        # All anchor corners should still be at their anchor values.
        for t in tris:
            for v in t:
                if abs(v[0]) < 0.5 and abs(v[1]) < 0.5:
                    assert v[2] == pytest.approx(30.0, abs=0.01)
                if abs(v[0] - 100) < 0.5 and abs(v[1]) < 0.5:
                    assert v[2] == pytest.approx(30.5, abs=0.01)

    def test_anchor_preserved_across_grade_clamp(self):
        """When DEM disagrees with anchors, the grade clamp must
        smooth the surface WITHOUT moving anchored corners."""
        # Left edge anchored high, right corners at DEM (25).
        anchors = [(0, 0, 30.0), (0, 100, 30.0)]
        tris = SM.adaptive_triangulate(
            SQUARE_100, anchors, constant_dem(25.0),
            max_grade=0.010)
        anchored_zs = []
        for t in tris:
            for v in t:
                if abs(v[0]) < 0.5 and (abs(v[1]) < 0.5 or
                                         abs(v[1] - 100) < 0.5):
                    anchored_zs.append(v[2])
        assert anchored_zs, "no anchored vertices found"
        for z in anchored_zs:
            assert z == pytest.approx(30.0, abs=0.01)


# ──────────────────────────────────────────────────────────────────────
# adaptive_triangulate — grade enforcement
# ──────────────────────────────────────────────────────────────────────
class TestGradeEnforcement:
    def test_steep_dem_clamped(self):
        """A 5 % DEM slope under a 1 % grade limit must be clamped to
        exactly 1 % on every emitted triangle."""
        dem = linear_dem(0.05, 0.0, 10.0)   # 5 %
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], dem, max_grade=0.010)
        for t in tris:
            plane = SM.fit_plane(*t)
            if plane is None:
                continue
            assert SM.plane_grade(plane) == pytest.approx(0.010, abs=1e-3)

    def test_steep_dem_with_relaxed_limit(self):
        """The same steep DEM under a 1.5 % limit should emit a steeper
        triangulation, also clamped to the limit."""
        dem = linear_dem(0.05, 0.0, 10.0)
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], dem, max_grade=0.015)
        for t in tris:
            plane = SM.fit_plane(*t)
            if plane is None:
                continue
            assert SM.plane_grade(plane) == pytest.approx(0.015, abs=1e-3)

    def test_compliant_dem_unchanged(self):
        """A 0.7 % DEM under a 1 % limit should be left exactly as-is
        — no clamping pull at all."""
        dem = linear_dem(0.007, 0.0, 10.0)
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], dem, max_grade=0.010)
        for t in tris:
            plane = SM.fit_plane(*t)
            if plane is None:
                continue
            assert SM.plane_grade(plane) == pytest.approx(0.007, abs=1e-3)

    def test_impossible_anchors_violation_accepted(self):
        """If all 4 corners are anchored at heights that demand a 5 %
        slope, the grade clamp can't fix it (anchors don't move) and
        the violation is accepted.  Caller's reconciliation should
        adjust building elevations later."""
        anchors = [(0, 0, 10.0), (0, 100, 10.0),
                   (100, 0, 15.0), (100, 100, 15.0)]
        tris = SM.adaptive_triangulate(
            SQUARE_100, anchors, constant_dem(12.0),
            max_grade=0.010)
        # We expect at least one triangle with grade > 0.010 because
        # the anchors are physically incompatible with 1 %.
        max_grade_seen = 0.0
        for t in tris:
            plane = SM.fit_plane(*t)
            if plane is None:
                continue
            g = SM.plane_grade(plane)
            if g > max_grade_seen:
                max_grade_seen = g
        assert max_grade_seen > 0.010


# ──────────────────────────────────────────────────────────────────────
# Integration / smoke
# ──────────────────────────────────────────────────────────────────────
class TestSmoke:
    def test_handles_empty_polygon(self):
        empty = Polygon()
        assert SM.adaptive_triangulate(
            empty, [], constant_dem(10.0), 0.010) == []

    def test_handles_no_dem(self):
        """When sample_dem returns None for everything, the algorithm
        should still produce SOME triangulation using zero elevations
        (rather than crashing)."""
        tris = SM.adaptive_triangulate(
            SQUARE_100, [], lambda x, y: None, max_grade=0.010)
        assert len(tris) >= 1
