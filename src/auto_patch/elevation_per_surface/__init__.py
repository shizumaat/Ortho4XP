"""Per-surface elevation pipeline (user 2026-05-02).

Replacement for the unified Laplacian solver in
``auto_patch.elevation._solve_pavement_elevations_unified``.

The unified solver propagates runway elevations across all pavement
through a single graph with per-edge grade caps, which over-flattens
terrain-following surfaces (CYXY taxi E, HECA elevated taxiways).
This package implements a per-surface phased solve that respects
the FAA per-axis grade rule:

* Phase 1 — taxi rects get an axial DEM-smoothed profile capped at
  1.5 %/m.
* Phase 2 — BFS from runway HARD anchors propagates through the
  rect-junction chain along realistic taxi axes only.
* Phase 3 — junctions take boundary vertices from adjacent
  rects/runways; interior vertices follow DEM smoothed to 1.5 %.
* Phase 4 — aprons same as junctions but capped at 1.0 %.
* Phase 5 — terminals flat at DEM-median (existing rule).

See ``docs/elevation_per_surface_redesign.md`` for the full plan.
"""
from .solver import solve

__all__ = ["solve"]
