"""Post-phase-1 finalisation: phase-2 elevation + feature emit.

Run after ``pipeline.build_airport_pavement`` has emitted the
phase-1 geometry layout (rects, junctions, terminals, runways).
Performs:

* Phase-2 elevation solve via ``elevation._compute_elevations``.
* Geometry repair passes (shared-vertex enforce + overlap clip)
  to handle subdivisions / decompositions the elevation step may
  have introduced.
* Junction-altitude reconciliation (snap-to-rect-corner +
  shared-vertex average + within-junction adjacent-pair grade
  smoother).
* Final WARN summary of within-shape grade violations.
* Feature emit: airport boundary, groundside pavement, boundary→DEM
  bridges; bridges + tunnel portals when ``EMIT_BRIDGES_AND_TUNNELS``.

Public API:
    run_phase2(layout, icao, xplane_root, apt, *, nodes, ways, to_m,
               apron_candidates)
"""
from __future__ import annotations

import math

import O4_UI_Utils as UI
from shapely.errors import GEOSException, TopologicalError

# Narrow exception tuple for feature-emit failures (shapely +
# OSM/DEM file I/O).  Programming errors propagate so they surface
# immediately rather than being silently masked at runtime.
_GEOM_EXC = (OSError, ValueError, TypeError, KeyError,
             IndexError, RuntimeError,
             GEOSException, TopologicalError)

from .boundary import (
    _emit_airport_boundary_shape,
    _emit_boundary_dem_bridge,
)
from .bridges import (
    _emit_taxi_bridges,
    _emit_through_airport_depressed_roads,
    _emit_tunnel_portals,
    _emit_underpass_road_approaches,
    _scenery_has_bridge_objects,
)
from .config import EMIT_BRIDGES_AND_TUNNELS
from .elevation import (
    SHARED_VERTEX_CLUSTER_TOL_M,
    _compute_elevations,
    _drop_overlap_against_fixed_shapes,
    _enforce_shared_vertex_altitudes,
    _load_airport_dem,
    _merge_sliver_junctions_into_neighbours,
    _report_within_shape_violations,
    _smooth_within_junction_adjacent_pair_grade,
    _snap_junction_altitudes_to_rect_corners,
)
from .groundside import (
    _drop_groundside_orphan_junctions,
    _emit_groundside_pavement_dem,
)
from .pavement.vertices import (
    _enforce_shared_vertices,
    _push_junction_vertices_off_taxi_rect_edges,
)


__all__ = ["run_phase2"]


def run_phase2(layout, icao, xplane_root, apt, *,
               nodes, ways, to_m, apron_candidates,
               tile_dem=None,
               current_tile_lat=None, current_tile_lon=None):
    """Phase-2 elevation solve + feature emit.  Mutates layout.

    ``tile_dem`` (when supplied by the tile-pipeline driver) is the
    pre-loaded ``O4_DEM_Utils.DEM`` from Ortho4XP's
    ``smooth_raster_over_airports`` step.  Threaded through so the
    Phase-2 elevation solver and the boundary-shape emit consume
    the SAME smoothed DEM that drives Ortho4XP's flattening,
    without each per-airport pass loading DEM tiles independently.

    ``current_tile_lat`` / ``current_tile_lon`` identify the tile
    being processed by the driver — used to provide correct
    ``(tile_lat, tile_lon)`` context for DEM sampling on cross-tile
    airports (the airport's anchor tile may differ from the tile
    Ortho4XP is currently generating).  Default ``None`` falls back
    to ``floor(layout.anchor)`` for direct test runs.
    """
    _compute_elevations(
        layout, icao, xplane_root, apt,
        osm_nodes=nodes, osm_ways=ways, to_m=to_m,
        apron_candidates_m=apron_candidates,
        tile_dem=tile_dem)
    # Elevation phase can subdivide junctions, decompose holed
    # polygons, and otherwise modify polygon geometry — re-run
    # the shared-vertex collapse + overlap-clip so the
    # invariants survive into the final layout.
    _enforce_shared_vertices(
        layout, tol=SHARED_VERTEX_CLUSTER_TOL_M)
    _drop_overlap_against_fixed_shapes(layout, icao=icao)
    _enforce_shared_vertices(
        layout, tol=SHARED_VERTEX_CLUSTER_TOL_M)
    # The overlap-clip pass introduces new vertices at
    # intersection points that may land on a taxi rect's
    # edge interior (would split the rect's altitude_high/
    # altitude_low convention at render time).  Push any
    # such vertex off — pure geometry, doesn't touch
    # elevations.
    _push_junction_vertices_off_taxi_rect_edges(layout)
    # Per user 2026-05-03: when the per-surface solver runs the
    # legacy chain of post-elevation reconciliation passes
    # (snap-to-rect-corner, shared-vertex agree, adjacent-pair
    # smoother) re-introduces the within-junction grade violations
    # the solver just fixed.  In particular,
    # ``_snap_junction_altitudes_to_rect_corners`` propagates a
    # rect's ``altitude_high``/``altitude_low`` scalars to junction
    # vertices via the legacy [high, low, low, high] convention,
    # which restores the cross-section delta the solver had
    # equalised to flat.  Skip the chain when the flag is on; the
    # solver's writeback already gives every shared bucket one
    # value, so the snaps are redundant and harmful.
    from .elevation import USE_PER_SURFACE_SOLVER
    if not USE_PER_SURFACE_SOLVER:
        # Per user 2026-04-28: junction polygon vertices that
        # coincide with a runway / sloping-rect / terminal corner
        # MUST emit that shape's altitude tag value.  Run this
        # AFTER all the shared-vertex / overlap-clip passes since
        # they can rearrange polygon coords and put node_altitudes
        # out of sync with the rect's emitted altitude tags.
        _snap_junction_altitudes_to_rect_corners(layout)
        # Per user 2026-04-28: junctions sharing a boundary vertex
        # MUST agree on its altitude.  Subdivide / clamp passes can
        # leave sub-metre disagreement at shared buckets — average
        # them so X-Plane doesn't render a tear at the seam.
        _enforce_shared_vertex_altitudes(layout)
        # Re-run the rect-corner snap after the shared-vertex
        # average, since averaging can pull a shared-with-rect
        # bucket away from the rect's tag value.
        _snap_junction_altitudes_to_rect_corners(layout)
        # Per user 2026-04-28: smooth adjacent-vertex pair grade
        # WITHIN each junction.  Above passes enforce shared-
        # vertex agreement (across polygons) and rect-corner
        # alignment (junction ↔ sloping rect), but interior
        # junction vertices can still violate 1.5 % grade against
        # their immediate ring neighbours.  Iterate adjacent-pair
        # smoothing with hard-vertex anchoring to converge those.
        _smooth_within_junction_adjacent_pair_grade(layout)
        # Re-run shared-vertex + rect-corner snaps so any seam
        # vertex the smoother nudged off-target is restored.
        _enforce_shared_vertex_altitudes(layout)
        _snap_junction_altitudes_to_rect_corners(layout)
    # Per user 2026-04-29: merge small junction slivers into
    # adjacent larger junctions.  Polygon-with-holes
    # decomposition + post-elevation subdivisions can carve
    # off small (< 1000 m²) pieces sharing a boundary segment
    # with a much larger neighbour (HECA -10244 = 485 m²
    # adjacent to -10243 = 30 k m²).  Merge them back so
    # JOSM doesn't show two near-duplicate polygons.
    _merge_sliver_junctions_into_neighbours(layout, icao=icao)
    # Per user 2026-05-03: when the per-surface solver is on,
    # geometry passes above (overlap clip, push-off, sliver merge)
    # may have moved or added vertices since the solver finished.
    # Re-run the solver one more time so its constraints (per-axis
    # rect grade, all-pair junction grade, cross-section flatness)
    # apply to the final geometry the OSM emit will write.
    # Per user 2026-05-03: skip the mid-pipeline per-surface solver
    # call and run the final solver AT THE END (below) after every
    # geometry pass has settled.
    # WARN summary moved to pipeline.py after the absolute final
    # per-surface solver pass — the WARN that ran here would
    # report a stale state because the post-finalize junction-rule
    # passes (widen_junctions_to_runway_corners etc.) and the
    # final solver run after this point.
    # Per user 2026-04-28: emit a 5 m-wide ribbon polygon
    # tracing the airport boundary (apt.dat row-130) with
    # per-vertex altitudes clamped to ≤ 3 % grade from the
    # nearest runway within 400 m, falling back to DEM
    # beyond.  Provides the elevation transition between
    # airport pavement and surrounding terrain.
    try:
        _lat0, _lon0 = layout.anchor
        # Use the current tile being processed (from the driver)
        # for DEM sampling context.  When ``tile_dem`` is the DEM
        # for the current tile (-13/-78 for SPLP when generating
        # that tile), ``_sample_dem`` needs ``tile_lat=-13,
        # tile_lon=-78`` to convert lat/lon → DEM local coords
        # correctly.  Anchor-derived tile coords would be wrong
        # for cross-tile airports.
        _tile_lat = (current_tile_lat if current_tile_lat is not None
                     else int(math.floor(_lat0)))
        _tile_lon = (current_tile_lon if current_tile_lon is not None
                     else int(math.floor(_lon0)))
        _dem = _load_airport_dem(_lat0, _lon0, override_dem=tile_dem)
        n_b = _emit_airport_boundary_shape(
            layout, _dem, _tile_lat, _tile_lon)
        if n_b:
            UI.vprint(1,
                f"  [pav-builder] emitted "
                f"{n_b} airport-boundary shape piece(s).")
        # Per user 2026-04-29: re-emit groundside pavement
        # captured before the airside-apron subtraction, with
        # per-vertex DEM altitudes and a 0.1 m gap from the
        # terminal building.  Lets curbside / drop-off /
        # parking pavement render at local terrain elevation
        # (CYXY: terminal cut into hill, ~4 m higher than the
        # airside apron).
        try:
            n_gs = _emit_groundside_pavement_dem(
                layout, _dem, _tile_lat, _tile_lon)
            if n_gs:
                UI.vprint(1,
                    f"  [pav-builder] emitted "
                    f"{n_gs} groundside pavement "
                    f"polygon(s) with DEM altitudes.")
        except _GEOM_EXC:
            pass
        # Per user 2026-04-29 (CYXY -10111 / -10115): drop
        # junction polygons that are connected ONLY to non-
        # airside pavement (groundside polygons or each
        # other), with no shared vertices on any airside
        # rect/terminal/runway.  These slivers got assigned
        # the airside-flat altitude during the Laplacian
        # solver but visually they're not airside — they
        # share an edge with the DEM-following groundside
        # at a different elevation, creating "valleys" that
        # X-Plane renders as cliffs.  Eliminating them lets
        # the boundary ribbon and the groundside polygon
        # define the surface in that area.
        try:
            n_orph = _drop_groundside_orphan_junctions(layout)
            if n_orph:
                UI.vprint(1,
                    f"  [pav-builder] dropped {n_orph} "
                    f"junction(s) sharing vertices with "
                    f"groundside pavement.")
        except _GEOM_EXC:
            pass
        # Then emit DEM-bridge polygons inside the boundary
        # wherever the clamped boundary altitude differs from
        # raw DEM by > 5 m (per user 2026-04-28).
        try:
            n_br = _emit_boundary_dem_bridge(
                layout, _dem, _tile_lat, _tile_lon)
            if n_br:
                UI.vprint(1,
                    f"  [pav-builder] emitted "
                    f"{n_br} boundary→DEM bridge "
                    f"polygon(s).")
        except _GEOM_EXC:
            pass
        # TODO(bridges): re-enable bridge / tunnel emission
        # after core pavement geometry is stable + refactored.
        # See ``EMIT_BRIDGES_AND_TUNNELS`` at module top for
        # the full to-do list.  Each gated call must carve its
        # footprint out of overlapping airside / groundside
        # pavement before emitting, or
        # ``test_no_self_overlap`` will fail.
        if EMIT_BRIDGES_AND_TUNNELS:
            # Per user 2026-04-29 (KPHX Sky Harbor Blvd): when
            # a public road passes under any airport bridge,
            # its ENTIRE inside-airport stretch must be at
            # apt_elev − 8 m, not just at the bridge crossings.
            # Run BEFORE _emit_tunnel_portals so the latter can
            # skip OSM ways already depressed here.
            _depressed_way_ids: set = set()
            try:
                n_dep, _depressed_way_ids = (
                    _emit_through_airport_depressed_roads(
                        layout, _dem, _tile_lat, _tile_lon,
                        xplane_root=xplane_root, icao=icao))
                if n_dep:
                    UI.vprint(1,
                        f"  [pav-builder] emitted "
                        f"{n_dep} through-airport depressed "
                        f"road segment(s).")
            except _GEOM_EXC:
                _depressed_way_ids = set()
            # Per user 2026-04-29: re-enable tunnel-portal
            # emission.  For each big-roads tunnel crossing the
            # airport boundary, emit a sloped ramp + 3
            # retaining walls that transition the road from
            # outside-DEM elevation down to airport-elevation
            # − 6 m at the portal.  Subtracts tunnel zones from
            # the boundary ribbon and DEM-bridge polygons so
            # they don't overlap.  OSM way ids handled by the
            # through-airport depressed-road emit are excluded
            # so we don't double-emit on Sky-Harbor-style
            # multi-bridge crossings.
            try:
                n_tun = _emit_tunnel_portals(
                    layout, _dem, _tile_lat, _tile_lon,
                    excluded_way_ids=_depressed_way_ids)
                if n_tun:
                    UI.vprint(1,
                        f"  [pav-builder] emitted "
                        f"{n_tun} tunnel-portal cluster(s) "
                        f"(ramp + walls along approach).")
            except _GEOM_EXC:
                pass
            # Per user 2026-04-29: emit retaining walls along
            # taxi bridges (KBNA Taxiway A, KPHX taxis over
            # Sky Harbor Blvd) and road-following approach
            # shapes that descend from outside-DEM down to
            # apt_elev − 8 m under the bridge.  When the
            # scenery pack already includes 3D bridge OBJs
            # (KBNA), the user wants the road to cut straight
            # through and let the OBJ be the bridge — skip our
            # walls and emit the under-bridge flat polygon.
            # When it doesn't (KPHX), keep the terrain flat
            # for the taxiway and only ramp the road up to the
            # bridge edge — emit walls + skip under-bridge.
            try:
                _scn_bridge = _scenery_has_bridge_objects(layout)
            except _GEOM_EXC:
                _scn_bridge = False
            try:
                n_brg = _emit_taxi_bridges(
                    layout, _dem, _tile_lat, _tile_lon,
                    scenery_has_bridge_objects=_scn_bridge)
                if n_brg:
                    UI.vprint(1,
                        f"  [pav-builder] emitted "
                        f"{n_brg} taxi-bridge wall pair(s).")
                elif _scn_bridge:
                    UI.vprint(1,
                        f"  [pav-builder] {icao}: scenery has "
                        f"3D bridge OBJ(s); skipping wall "
                        f"emission.")
            except _GEOM_EXC:
                pass
            try:
                n_app = _emit_underpass_road_approaches(
                    layout, _dem, _tile_lat, _tile_lon,
                    scenery_has_bridge_objects=_scn_bridge)
                if n_app:
                    UI.vprint(1,
                        f"  [pav-builder] emitted underpass-"
                        f"road approaches for {n_app} "
                        f"surface(s)"
                        f"{' (cut through under bridge OBJ)' if _scn_bridge else ' (ramp up to bridge edge)'}.")
            except _GEOM_EXC:
                pass
    except _GEOM_EXC:
        pass

    # The per-surface solver final pass moved to pipeline.py — it
    # must run AFTER pipeline.py's own post-finalize passes
    # (widen_junctions_to_runway_corners, etc.) which can insert new
    # vertices into junctions.  See pipeline.py at end of
    # build_airport_pavement.
