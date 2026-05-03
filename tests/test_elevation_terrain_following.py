"""Regression tests for terrain-following elevation (user 2026-05-02).

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

import sys
from pathlib import Path
from typing import Optional, Tuple

import pytest

from conftest import xplane_available, xplane_root


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
