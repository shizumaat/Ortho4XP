"""Backward-compatibility shim for the moved pavement builder.

The original 17,108-line monolith has been carved into focused
``O4_Pavement_*`` modules.  Orchestration now lives in
``O4_Pavement_Pipeline``; geometry / elevation / feature emit
live in their respective modules.

This shim re-exports the public entry point and the handful of
private helpers that other modules + integration tests still
import by their original ``O4_Airport_Pavement_Builder.<name>``
path.  New code should import directly from the focused modules.

Usage (unchanged from before the refactor)::

    from O4_Airport_Pavement_Builder import build_airport_pavement
    layout = build_airport_pavement("SPJC", "/path/to/X-Plane 12")
    layout.to_osm("/tmp/SPJC_auto.osm")
"""
from __future__ import annotations

# Public entry points + types.
from O4_Pavement_Pipeline import build_airport_pavement
from O4_Pavement_Layout import (
    AEROWAY_FOR_ROLE,
    BuiltShape,
    PavementLayout,
    R_EARTH,
    ROLE_APRON,
    ROLE_BOUNDARY,
    ROLE_CROSS_CONNECTOR,
    ROLE_GROUNDSIDE_PAVEMENT,
    ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RETAINING_WALL,
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
    ROLE_TUNNEL_RAMP,
    SHARED_VERTEX_TOL_M,
    _airport_anchor,
    _projection,
)

# Helpers re-exported for the lazy-import pattern used by extracted
# modules (Centerlines, Junctions, Boundary, Groundside, Bridges)
# and by integration tests that still import by the legacy path.
from O4_Pavement_Pipeline import (
    _add_stub_to_runway_bridges,
    _clip_residue_at_stub_long_edges,
    _load_osm_airports,
    _load_osm_big_roads,
    _load_osm_tile,
    _merge_near_touching,
    _osm_tile_path,
    _pick_best_apt_dat_against_osm,
    _score_apt_dat_against_osm,
)
from O4_Pavement_Rects import _natural_half_width
from O4_Pavement_Runways import _runway_rect_m
from O4_Pavement_Elevation import _corner_elevation_bucket
