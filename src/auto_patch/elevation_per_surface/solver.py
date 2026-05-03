"""Per-surface elevation solver — top-level entry point.

Delegates to ``unified_jacobi.solve``.  The DEM + tile coords are
accepted for API parity with the legacy unified solver and for
future use (per-vertex DEM seeding) but are currently unused: the
solver warm-starts soft nodes from existing layout altitudes
(rect altitude_high/low, junction node_altitudes, terminal altitude)
which were already populated upstream from DEM.

See ``unified_jacobi`` for the per-axis grade rule and the
three-change derivation from the legacy unified solver.
"""
from __future__ import annotations

from .unified_jacobi import solve as _jacobi_solve


def solve(layout, icao: str,
          dem=None, tile_lat: int = 0, tile_lon: int = 0) -> None:
    """Per-surface phased elevation solve.  Mutates ``layout`` in
    place: writes ``altitude_high``/``altitude_low`` on rects,
    ``node_altitudes`` on junctions, ``altitude`` on terminals and
    aprons.  Runway segments (HARD anchors) are left untouched.

    When ``dem`` is supplied, SOFT nodes are seeded from per-vertex
    DEM samples — necessary so taxi rects/junctions reach the
    DEM-driven elevations the user expects (e.g. CYXY taxi E sits
    on terrain ~717 m, not pulled down to the 705 m runway).
    """
    _jacobi_solve(layout, icao, dem=dem,
                   tile_lat=tile_lat, tile_lon=tile_lon)
