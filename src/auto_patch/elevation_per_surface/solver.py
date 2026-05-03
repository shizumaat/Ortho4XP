"""Per-surface elevation solver — top-level orchestrator.

Replaces ``auto_patch.elevation._solve_pavement_elevations_unified``
behind a feature flag.  See
``docs/elevation_per_surface_redesign.md`` for the design.
"""
from __future__ import annotations

import math
import time as _time

from .axial_profile import apply_axial_profiles
from .bfs_propagate import (
    build_runway_anchor_lookup,
    propagate_through_rects,
)
from .junction_field import apply_apron_field, apply_junction_field
from .continuity import collect_anchor_map, reconcile_shared_vertices


def solve(layout, icao: str, dem, tile_lat: int, tile_lon: int) -> None:
    """Per-surface phased elevation solve.

    The runway shapes must already carry their CIFP profile in
    ``altitude_high`` / ``altitude_low``.  Terminal flat altitudes
    are set by the existing terminal-altitude rule before this
    runs.  This solver fills in:

      * taxi rect altitude_high / altitude_low (Phase 1 + 2)
      * junction node_altitudes (Phase 3)
      * apron node_altitudes (Phase 4)

    and finally reconciles shared-vertex disagreements.
    """
    t_start = _time.time()
    runway_anchor_lookup = build_runway_anchor_lookup(layout)

    n_phase1 = apply_axial_profiles(
        layout, dem, tile_lat, tile_lon, runway_anchor_lookup)
    n_phase2 = propagate_through_rects(
        layout, dem, tile_lat, tile_lon, runway_anchor_lookup)

    # Phases 3 + 4: junctions and aprons key off rect / runway /
    # terminal corners.  Snapshot the anchor map after rect solves.
    anchor_map = collect_anchor_map(layout)
    n_phase3 = apply_junction_field(
        layout, dem, tile_lat, tile_lon, anchor_map)
    n_phase4 = apply_apron_field(
        layout, dem, tile_lat, tile_lon, anchor_map)

    n_reconciled = reconcile_shared_vertices(layout)

    elapsed = _time.time() - t_start
    try:
        import O4_UI_Utils as UI
        UI.vprint(1,
            f"  [pav-builder] {icao}: per-surface solver "
            f"({elapsed:.2f} s); rects={n_phase1}+{n_phase2} "
            f"junctions={n_phase3} aprons={n_phase4} "
            f"reconciled_buckets={n_reconciled}.")
    except Exception:
        pass
