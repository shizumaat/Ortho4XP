"""Per-surface elevation solver (user 2026-05-02 / 2026-05-03).

Replacement for the legacy unified Laplacian solver in
``auto_patch.elevation._solve_pavement_elevations_unified``.

The legacy solver treats all pavement as one connected graph with
per-edge grade caps and propagates runway HARD anchors through
every connected vertex.  This violates the per-axis FAA grade rule:

* Taxi rect grade applies along source_axis only — never
  perpendicular to an unrelated surface like a parallel runway.
* Junction / apron grade applies in any direction (multi-directional)
  within the polygon's surface, but not across to a different shape.
* Terminals are flat; aprons follow terrain at 1.5 % grade.
* Only CIFP runway corners are immutable HARD anchors — terminals,
  aprons, and taxi rects can all adjust.

The implementation lives in ``unified_jacobi`` as a single damped
Jacobi solver with role-aware edge generation (rects ring-only;
junctions / aprons / terminals ring + all-pair Euclidean).  See
``docs/elevation_per_surface_redesign.md`` for the full design.
"""
from .solver import solve

__all__ = ["solve"]
