"""Regression tests for terrain-following elevation (user 2026-05-02
and 2026-05-03).

The current unified Laplacian solver propagates elevation from
runway HARD anchors outward through the pavement graph at
1.5 % / 1.0 % per-edge caps, which OVER-FLATTENS surfaces whose
natural terrain elevation is significantly above the runway.

These tests assert the elevation pipeline produces taxi/apron
elevations within striking distance of the local DEM, not collapsed
toward the runway.

See ``docs/elevation_per_surface_redesign.md`` for the redesign plan.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pytest

from conftest import airports_under_test, xplane_available, xplane_root


_HERE = Path(__file__).resolve().parent
_SRC = _HERE.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


pytestmark = pytest.mark.skipif(
    not xplane_available(),
    reason="X-Plane install not found (set XPLANE_ROOT to override)",
)


_LAYOUT_CACHE: dict = {}


def _build_layout(icao: str):
    if icao in _LAYOUT_CACHE:
        return _LAYOUT_CACHE[icao]
    from auto_patch.pipeline import build_airport_pavement
    layout = build_airport_pavement(
        icao, xplane_root(), compute_elevations=True)
    _LAYOUT_CACHE[icao] = layout
    return layout


def _shape_max_corner_alt(shape) -> Optional[float]:
    """Return the highest elevation that will actually be emitted to
    OSM for this shape.

    Mirrors ``layout.PavementLayout.to_osm`` field priority:
    ``altitude_high`` (sloped rect / runway) > ``node_altitudes``
    (junction) > ``altitude`` (flat).  Avoids reading stale fields
    that don't drive the rendered surface.
    """
    if shape.altitude_high is not None and shape.altitude_low is not None:
        return float(shape.altitude_high)
    if shape.node_altitudes:
        return max(float(a) for a in shape.node_altitudes)
    if shape.altitude is not None:
        return float(shape.altitude)
    return None


def _shapes_in_bbox(layout, x_lo: float, x_hi: float,
                    y_lo: float, y_hi: float):
    """Yield shapes whose polygon centroid falls inside the bbox."""
    for s in layout.shapes:
        if s.polygon is None or s.polygon.is_empty:
            continue
        c = s.polygon.centroid
        if x_lo <= c.x <= x_hi and y_lo <= c.y <= y_hi:
            yield s


# ── CYXY taxi E / SW apron — terrain follows DEM (user 2026-05-02) ─


def test_cyxy_taxi_e_south_apron_follows_terrain():
    """Taxi E at the south edge of CYXY's SW apron sits on natural
    terrain ~717 m (DEM truth).  The connected runway is at ~704 m
    (CIFP HARD).  The pavement chain back to the runway is long
    enough to absorb the ~14 m delta at 1.5 % grade per surface
    along its own axis, so the per-surface elevation pipeline must
    reach ≥ 714 m at the south end of the SW apron region.

    Failure indicates the solver is over-flattening (propagating
    runway elevation across surfaces via Euclidean / graph-distance
    constraints rather than per-axis grade compliance).

    Reference: ``docs/elevation_per_surface_redesign.md``.
    """
    if "CYXY" not in {"CYXY"}:
        pytest.skip("CYXY-specific regression")
    layout = _build_layout("CYXY")

    # SW apron region (south edge): in CYXY meter coordinates the
    # SW apron sits south-east of the runway at roughly
    # (50..400, -1100..-700).  Taxi E's SW stub centroid is at
    # (124, -820); the DEM in this region is 715-718 m.
    bbox = (-100.0, 500.0, -1200.0, -700.0)
    candidates = list(_shapes_in_bbox(layout, *bbox))
    assert candidates, (
        "CYXY: no shapes found in SW apron bbox "
        f"x∈[{bbox[0]},{bbox[1]}] y∈[{bbox[2]},{bbox[3]}]")

    # Look for the highest pavement elevation reached anywhere in
    # the bbox.  Per user 2026-05-02: must reach ≥ 714 m.
    REQUIRED_M = 714.0

    best_alt = float("-inf")
    best_role = None
    for s in candidates:
        a = _shape_max_corner_alt(s)
        if a is None:
            continue
        if a > best_alt:
            best_alt = a
            best_role = s.role

    assert best_alt >= REQUIRED_M, (
        f"CYXY SW apron: highest pavement elevation in region is "
        f"{best_alt:.1f} m on a {best_role}, expected ≥ {REQUIRED_M:.1f} m. "
        f"Likely over-flattening — the per-surface elevation solver is "
        f"propagating runway altitudes across surfaces instead of letting "
        f"taxi/apron follow DEM along their own axes.")


# ── Within-junction grade audit (user 2026-05-03) ─────────────────


JUNCTION_MAX_GRADE = 0.015  # 1.5 %, matches FAA taxiway cap
APRON_MAX_GRADE = 0.015     # 1.5 % (user 2026-05-18): aligned with
                              # the taxi cap so reclassified apron
                              # shapes stay feasible.

# Absolute rounding allowance for the within-junction grade audit.
# Elevations are stored to 0.1 m precision, so a vertex pair can
# have up to 0.10 m of rounding error.  Plus the solver's own
# convergence tolerance (tol_m=0.005) accumulated across iterations
# can leave residual cap-violation under 0.10 m.  We use 0.20 m as
# the absolute slack — captures genuine grade-rule failures
# (typically 0.5 m+ excess on small junctions) while letting
# rounding-noise pairs pass.
ROUNDING_ALLOWANCE_M = 0.20


def _within_polygon_violations(coords, alts, max_grade: float):
    """All-pair Euclidean grade audit on a polygon.  Returns a list
    of (i, j, dist_m, de_m, grade) for pairs whose
    ``|de| > max_grade * d + ROUNDING_ALLOWANCE_M``.
    """
    out = []
    n = len(coords)
    for i in range(n):
        x1, y1 = coords[i]
        for j in range(i + 1, n):
            x2, y2 = coords[j]
            d = math.hypot(x2 - x1, y2 - y1)
            if d < 1.0:
                continue
            de = abs(alts[i] - alts[j])
            allowed = max_grade * d + ROUNDING_ALLOWANCE_M
            if de > allowed:
                out.append((i, j, d, de, de / d))
    return out


def _polygon_alts_and_coords(shape):
    """Return ``(coords_open, alts_open)`` aligned to the polygon's
    open exterior ring.  ``alts_open`` length matches ``coords_open``;
    excludes the closing duplicate vertex.
    """
    if shape.polygon is None or shape.polygon.is_empty:
        return None, None
    coords = list(shape.polygon.exterior.coords)
    ring_closed = coords and coords[0] == coords[-1]
    coords_open = coords[:-1] if ring_closed else coords
    if shape.node_altitudes:
        alts = list(shape.node_altitudes)
        if ring_closed and len(alts) == len(coords):
            alts = alts[:-1]
        if len(alts) != len(coords_open):
            return None, None
        return coords_open, alts
    if shape.altitude is not None:
        return coords_open, [float(shape.altitude)] * len(coords_open)
    return None, None


# Per-airport regression baselines.  Strict: zero within-junction or
# within-apron grade violations.  When this test fails, the elevation
# pipeline is producing a junction or apron whose internal Euclidean
# gradient exceeds the FAA cap — a real geometric defect that JOSM
# inspection of the patch will confirm visually.
WITHIN_SHAPE_GRADE_BASELINE: Dict[str, int] = {
    # Empty by design.  Add airport-specific entries here only as
    # interim ceilings while a fix is staged.
}


@pytest.mark.parametrize("icao", airports_under_test() or [
    pytest.param("(no airports)", marks=pytest.mark.skip(
        reason="set O4_TEST_TILE=lat,lon or O4_TEST_AIRPORTS=ICAO,..."))])
def test_within_junction_grade_compliance(icao):
    """Junctions and aprons must satisfy ≤ 1.5 % Euclidean grade
    between any pair of their vertices.

    Per user 2026-05-03: junctions are multi-directional but still
    capped at 1.5 % from edge to edge in any direction.  An aircraft
    can taxi across a junction surface in any direction, so the cap
    applies to every vertex pair, not just ring-adjacent ones.

    Failure indicates the elevation pipeline is producing a junction
    whose internal slope exceeds the cap — typically because the
    rect anchors on opposite sides of the junction are too far apart
    in elevation for the junction's spatial extent to bridge.  The
    fix lives in Phase 1 (rect-corner bound derivation from HARD
    anchors via junction grade-reach).
    """
    layout = _build_layout(icao)
    violations: List[str] = []
    for s_idx, s in enumerate(layout.shapes):
        if s.role == "junction":
            cap = JUNCTION_MAX_GRADE
        elif s.role == "apron":
            cap = APRON_MAX_GRADE
        else:
            continue
        coords, alts = _polygon_alts_and_coords(s)
        if coords is None or alts is None:
            continue
        bad = _within_polygon_violations(coords, alts, cap)
        for i, j, d, de, grade in bad:
            violations.append(
                f"{s.role}#{s_idx} verts {i}↔{j}: "
                f"de={de:.2f} m / d={d:.1f} m = {grade*100:.1f}% "
                f"(cap {cap*100:.1f}%, alts {alts[i]:.1f}→{alts[j]:.1f})")

    baseline = WITHIN_SHAPE_GRADE_BASELINE.get(icao, 0)
    if len(violations) > baseline:
        msg = (f"{icao}: {len(violations)} within-junction/apron "
               f"grade violations (baseline {baseline}).\n"
               + "First 10:\n  " + "\n  ".join(violations[:10]))
        if len(violations) > 10:
            msg += f"\n  ... and {len(violations) - 10} more"
        pytest.fail(msg)
