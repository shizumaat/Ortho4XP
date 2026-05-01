"""Pavement elevation: anchors → graph → solver → grade clamp.

Phase-2 elevation work: layer altitudes onto the role-classified
shape topology produced by Phase-1.  Pipeline:

* DEM tile loading + sampling for terrain-derived elevations.
* CIFP threshold anchors for runway profile.
* Apron-side multi-source-Dijkstra grade cone (FAA apron cap).
* Unified Laplacian Jacobi solver over the per-shape elevation
  graph.
* Per-shape grade clamp (taxiway / apron / runway).
* Junction triangulation + free-vertex clamp for the residue
  polygons.
* Sliver-junction merge + violating-junction subdivision.
* Geometric finalization (overlap clip against fixed shapes,
  shared-vertex altitude reconciliation, terminal apron
  re-derivation).

This module is large because the network ↔ solver ↔ grade triple
is tightly coupled — splitting it produces leaky abstractions.
The plan envelope (~900 lines) was optimistic; current size sits
near 3,700 lines.  Future iteration may split into ``_Solver``
and ``_Grade`` siblings if either half can be cleanly separated.

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    Constants:
      APRON_MAX_GRADE, DEM_SUFFIX, ELEVATION_GRID_STEP_M,
      ELEVATION_SMOOTH_CONVERGE_M, ELEVATION_SMOOTH_MAX_ITERS,
      SHARED_AGREE_TOL_M, SUBDIVIDE_MAX_PAIR_DIST_M,
      SUBDIVIDE_SNAP_RADIUS_M, TAXI_ANCHOR_DIST_M, TAXI_MAX_GRADE,
      USE_PER_POLYGON_ELEVATION_FIELD

    Functions (DEM + CIFP):
      _load_airport_dem, _sample_dem, _find_cifp_path

    Functions (main pipeline):
      _compute_elevations, _resample_node_altitudes_nn,
      _apply_geometric_finalization,
      _solve_pavement_elevations_unified

    Functions (per-shape elevation field):
      _smooth_within_junction_adjacent_pair_grade,
      _rederive_terminal_altitude_from_apron_neighbours,
      _enforce_shared_vertex_altitudes,
      _snap_junction_altitudes_to_rect_corners,
      _re_emit_apron_merged_runway_segments,
      _latlon_to_m_local, _orient_rect_for_altitude,
      _planar_fit, _planar_fit_residuals, _match_elev,
      _smooth_polygon_grid

    Functions (corner buckets + clamp + finalization):
      _corner_elevation_bucket, _corner_elev_map,
      _triangulate_junctions, _build_clamp_geom_state,
      _clamp_junction_free_vertices,
      _subdivide_violating_junctions,
      _merge_sliver_junctions_into_neighbours,
      _report_within_shape_violations,
      _drop_overlap_against_fixed_shapes
"""
from __future__ import annotations

import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import O4_UI_Utils as UI

from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon
from shapely.ops import linemerge, nearest_points, unary_union

from . import apt_dat_reader as APR

from .config import (
    EMIT_BRIDGES_AND_TUNNELS,
    RUNWAY_APRON_AREA_RATIO,
    RUNWAY_INSIDE_APRON_FRAC,
)
from .layout import (
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
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
    ROLE_RETAINING_WALL,
    SHARED_VERTEX_TOL_M,
)
from .pavement.vertices import (
    _drop_spike_vertices,
    _enforce_shared_vertices,
    _push_junction_vertices_off_taxi_rect_edges,
    _snap_polygon_vertices_to_rect_corners,
    _validate_shared_vertex_invariant,
)
from .pavement.runways import (
    _insert_runway_chain_bridges,
    _resolve_runway_crossings,
    _sample_runway_segment_elev,
)


__all__ = [
    "APRON_MAX_GRADE",
    "DEM_SUFFIX",
    "ELEVATION_GRID_STEP_M",
    "ELEVATION_SMOOTH_CONVERGE_M",
    "ELEVATION_SMOOTH_MAX_ITERS",
    "SHARED_AGREE_TOL_M",
    "SHARED_VERTEX_CLUSTER_TOL_M",
    "SUBDIVIDE_MAX_PAIR_DIST_M",
    "SUBDIVIDE_SNAP_RADIUS_M",
    "TAXI_ANCHOR_DIST_M",
    "TAXI_MAX_GRADE",
    "USE_PER_POLYGON_ELEVATION_FIELD",
    "_apply_geometric_finalization",
    "_build_clamp_geom_state",
    "_clamp_junction_free_vertices",
    "_compute_elevations",
    "_corner_elev_map",
    "_corner_elevation_bucket",
    "_drop_overlap_against_fixed_shapes",
    "_enforce_shared_vertex_altitudes",
    "_find_cifp_path",
    "_latlon_to_m_local",
    "_load_airport_dem",
    "_match_elev",
    "_merge_sliver_junctions_into_neighbours",
    "_orient_rect_for_altitude",
    "_planar_fit",
    "_planar_fit_residuals",
    "_re_emit_apron_merged_runway_segments",
    "_rederive_terminal_altitude_from_apron_neighbours",
    "_report_within_shape_violations",
    "_resample_node_altitudes_nn",
    "_sample_dem",
    "_smooth_polygon_grid",
    "_smooth_within_junction_adjacent_pair_grade",
    "_snap_junction_altitudes_to_rect_corners",
    "_solve_pavement_elevations_unified",
    "_subdivide_violating_junctions",
    "_triangulate_junctions",
]


# FAA AC 150/5300-13B, Design Group III+ (commercial jetport):
#   Taxiway longitudinal: max 1.5 % grade.
#   Apron: max 1.0 % any direction.
#   Runway longitudinal: max 1.5 % grade (handled by legacy).
TAXI_MAX_GRADE = 0.015
APRON_MAX_GRADE = 0.010   # FAA apron cap, used for the apron-side
                            # multi-source-Dijkstra cone in
                            # ``_pin_apron_targets`` and the post-
                            # pin grade reconciliation.  Tighter
                            # than the taxi cap so apron polygons
                            # don't slope > 1 % between vertices.
TAXI_ANCHOR_DIST_M = 30.0   # snap taxi rect end to runway segment
                             # elevation when within this distance

# ── Per-shape elevation field (Mode A / Mode B) ────────────────────
# Mode A: 1D smoothing along a rect's axis at this spacing.  Mode B:
# 2D grid smoothing within a non-rect pavement polygon at this
# spacing.  Per-step grade cap is ``ELEVATION_GRID_STEP_M ×
# TAXI_MAX_GRADE`` so two adjacent samples (axis or grid) cannot
# differ by more than that amount.
ELEVATION_GRID_STEP_M = 5.0
ELEVATION_SMOOTH_MAX_ITERS = 50
ELEVATION_SMOOTH_CONVERGE_M = 0.01
# Tolerance for cross-junction shared-bucket reconciliation: when
# two junctions touching the same SOFT (graph/DEM-derived) shared
# bucket end up with smoothed values farther apart than this, we
# overwrite both with the average to restore the shared-vertex
# invariant; below this they stay at their per-junction smoothed
# values (preserves the within-shape grade win).  0.10 m is ≈ 1 %
# grade across a 10 m shared edge — well below the visible-cliff
# threshold and below check_grade.py's CROSS-SHAPE 1.5 % bar.
SHARED_AGREE_TOL_M = 0.10
# Wire ``_smooth_polygon_grid`` (Mode B) into the per-junction flow
# instead of the legacy per-vertex anchor lookup + 1D ring smoothing.
# Default OFF (2026-04-26): Mode B's 2D-Euclidean anchor cones impose
# tighter grade constraints than the elevation graph's network-
# distance smoothing, which can cause regressions when rect corner
# anchors that fit network compliance sit outside the 2D cones,
# forcing wild midpoint fallbacks at free cells.  The helper +
# constants are kept so a future iteration can re-engage Mode B
# once rect-corner derivation is reworked to enforce 2D-Euclidean
# grade compliance.
USE_PER_POLYGON_ELEVATION_FIELD = False

# Used by both _build_clamp_geom_state (in this module) and the
# _triangulate_junctions code path (in auto_patch.triangulation).
# Defined here so triangulation can import it without forcing a
# load-order dance.
NEIGHBOUR_CLAMP_RADIUS_M = 5.0

DEM_SUFFIX = ".hgt"
_DEM_CACHE: Dict[Tuple[int, int], object] = {}

def _load_airport_dem(lat0: float, lon0: float):
    """Return an ``O4_DEM_Utils.DEM`` covering the 1° tile that
    contains (lat0, lon0).  Auto-downloads via Ortho4XP's standard
    DEM provider chain when no local .hgt file exists.  Falls back
    to None only when the download itself fails."""
    tile_lat = int(math.floor(lat0))
    tile_lon = int(math.floor(lon0))
    key = (tile_lat, tile_lon)
    if key in _DEM_CACHE:
        return _DEM_CACHE[key]
    hem_ns = "S" if tile_lat < 0 else "N"
    hem_ew = "W" if tile_lon < 0 else "E"
    fname = f"{hem_ns}{abs(tile_lat):02d}{hem_ew}{abs(tile_lon):03d}{DEM_SUFFIX}"
    # Ortho4XP lays out by 10° group.
    group_lat = (tile_lat // 10) * 10
    group_lon = (tile_lon // 10) * 10
    group_dir = (f"{'+' if group_lat >= 0 else '-'}{abs(group_lat):02d}"
                 f"{'+' if group_lon >= 0 else '-'}{abs(group_lon):03d}")
    dem_path = os.path.join("Elevation_data", group_dir, fname)
    try:
        import O4_DEM_Utils as _DEM
        if os.path.isfile(dem_path):
            dem = _DEM.DEM(tile_lat, tile_lon, source=dem_path)
        else:
            # Download via Ortho4XP's default DEM source chain
            # (SRTM/View — DEM.load_data calls ensure_elevation
            # internally).  Same path Ortho4XP's main pipeline
            # uses for the elevation data step.
            try:
                UI.vprint(1,
                    f"  [pav-builder] {fname} missing; downloading "
                    f"DEM tile via Ortho4XP elevation provider...")
            except Exception:
                pass
            dem = _DEM.DEM(tile_lat, tile_lon)
    except Exception as exc:
        try:
            UI.vprint(1,
                f"  [pav-builder] WARN: DEM load/download failed for "
                f"{fname}: {exc}")
        except Exception:
            pass
        _DEM_CACHE[key] = None
        return None
    _DEM_CACHE[key] = dem
    return dem


def _sample_dem(dem, tile_lat: int, tile_lon: int,
                lat: float, lon: float) -> Optional[float]:
    """Sample DEM elevation at (lat, lon).  Returns None if DEM is
    unavailable or out-of-tile."""
    if dem is None:
        return None
    try:
        return float(dem.alt((lon - tile_lon, lat - tile_lat)))
    except Exception:
        return None


def _find_cifp_path(xplane_root: str, icao: str) -> Optional[str]:
    """Locate the CIFP .dat file for an ICAO under the X-Plane
    root.  Returns None if not found."""
    cifp_dir = os.path.join(xplane_root, "Custom Data", "CIFP")
    p = os.path.join(cifp_dir, f"{icao.upper()}.dat")
    if os.path.isfile(p):
        return p
    return None


def _compute_elevations(layout: "PavementLayout", icao: str,
                        xplane_root: str, apt,
                        osm_nodes=None, osm_ways=None,
                        to_m=None,
                        apron_candidates_m: Optional[
                            List[Polygon]] = None) -> None:
    """Phase-2: add altitude tags to runways (segmented), taxi
    rects, and terminal pads.  Junctions / aprons / buildings are
    left un-elevated this iteration.

    Taxi rect elevations come from a grade-compliant elevation
    network built over OSM taxiway centerlines (densified to
    ≤ 30 m edges), anchored at CIFP runway thresholds, and
    post-pass smoothed for rate-of-change compliance
    (FAA 1 %/30 m).  Per user (2026-04-24): real airports
    heavily modify the land, so DEM is a soft preference, not a
    constraint — nodes free from runway anchors follow DEM
    only when no grade-compliance rule forces otherwise.
    """
    lat0, lon0 = layout.anchor
    tile_lat = int(math.floor(lat0))
    tile_lon = int(math.floor(lon0))
    dem = _load_airport_dem(lat0, lon0)

    # Meter-space projection (local — the layout's to_m is not
    # exposed, so reconstruct).
    cos0 = math.cos(math.radians(lat0))

    def m_to_ll(x: float, y: float):
        lon = lon0 + math.degrees(x / (R_EARTH * cos0))
        lat = lat0 + math.degrees(y / R_EARTH)
        return lat, lon

    # ── Segmented runway rectangles (legacy CIFP + DEM) ─────────
    cifp_path = _find_cifp_path(xplane_root, icao)
    runway_segment_chain = []
    if cifp_path is not None and dem is not None:
        try:
            from . import driver as _AP
            from . import cifp_reader as _CIFP
            from .pavement import runway_geometry as _RWY
            cifp_runways = _CIFP.parse_cifp_file(cifp_path)
            if cifp_runways:
                pairs = _RWY.pair_runways(cifp_runways)

                # apt.dat runway geometry (sole source of truth for
                # footprint lat/lon + width) — per legacy contract.
                apt_runway_geom = {}
                for r in apt.runways:
                    key_a = r.desig_a if r.desig_a.startswith("RW") \
                        else "RW" + r.desig_a
                    key_b = r.desig_b if r.desig_b.startswith("RW") \
                        else "RW" + r.desig_b
                    apt_runway_geom[key_a] = (
                        r.lat_a, r.lon_a, r.width_m,
                        r.displaced_a_m, r.blast_a_m)
                    apt_runway_geom[key_b] = (
                        r.lat_b, r.lon_b, r.width_m,
                        r.displaced_b_m, r.blast_b_m)
                runway_widths = {}
                for r in apt.runways:
                    runway_widths[r.desig_a] = r.width_m
                    runway_widths[r.desig_b] = r.width_m

                class _TileStub:
                    pass
                tile = _TileStub()
                tile.lat = tile_lat
                tile.lon = tile_lon
                tile.dem = dem

                _xml, runway_segment_chain = _AP.generate_patch_osm(
                    icao, pairs, runway_widths=runway_widths,
                    tile=tile, apt_runways=apt_runway_geom)
        except Exception:
            runway_segment_chain = []

    new_runway_polys: List[Polygon] = []
    if runway_segment_chain:
        # Drop the single-rect runway shapes; replace with segments.
        old_runways = [s for s in layout.shapes if s.role == ROLE_RUNWAY]
        layout.shapes = [s for s in layout.shapes if s.role != ROLE_RUNWAY]
        ref_fallback = "/".join(sorted(set(
            f"{r.desig_a}/{r.desig_b}" for r in apt.runways)))
        # Legacy generate_patch_osm pads each side by
        # RUNWAY_MARGIN=3 m for imagery coverage; strip that so
        # the segmented runways match the apt.dat width that our
        # Phase-1 junctions/rects were built against.  Keeps the
        # post-elevation layout overlap-free.
        _LEGACY_RUNWAY_MARGIN = 3.0
        for i, seg in enumerate(runway_segment_chain):
            lat_a, lon_a, elev_a, lat_b, lon_b, elev_b, width_m = seg
            width_m = max(1.0, width_m - 2.0 * _LEGACY_RUNWAY_MARGIN)
            ax, ay = _latlon_to_m_local(lat_a, lon_a, lat0, lon0, cos0)
            bx, by = _latlon_to_m_local(lat_b, lon_b, lat0, lon0, cos0)
            length = math.hypot(bx - ax, by - ay)
            if length < 1.0:
                continue
            # Perpendicular half-width offset — always compute
            # relative to the direction from A to B.
            ux = (bx - ax) / length
            uy = (by - ay) / length
            px = -uy * width_m / 2.0
            py = ux * width_m / 2.0
            # X-Plane patch convention: a way's short edge from
            # the last to the first node (way[-2:]) is interpreted
            # as the ``altitude_high`` side; the short edge from
            # way[1] to way[2] is the ``altitude_low`` side.  So
            # corners 0 and 3 must be at the HIGH-elevation end.
            # If B is the higher end, start the ring from B.
            if float(elev_a) >= float(elev_b):
                # A is HIGH: ring starts at A-side corners.
                corners = [
                    (ax + px, ay + py),   # 0: A-left  (HIGH side)
                    (bx + px, by + py),   # 1: B-left  (LOW side)
                    (bx - px, by - py),   # 2: B-right (LOW side)
                    (ax - px, ay - py),   # 3: A-right (HIGH side)
                ]
                eh, el = float(elev_a), float(elev_b)
            else:
                # B is HIGH: reverse — start the ring at B-side.
                # The perpendicular flips sign when direction
                # reverses, so B-left in the reversed walk is the
                # original B-right (and similarly A).
                corners = [
                    (bx - px, by - py),   # 0: B-left  (HIGH side)
                    (ax - px, ay - py),   # 1: A-left  (LOW side)
                    (ax + px, ay + py),   # 2: A-right (LOW side)
                    (bx + px, by + py),   # 3: B-right (HIGH side)
                ]
                eh, el = float(elev_b), float(elev_a)
            poly = Polygon(corners)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if poly.is_empty or poly.geom_type != "Polygon":
                continue
            shape = BuiltShape(
                polygon=poly, role=ROLE_RUNWAY, ref=ref_fallback)
            if abs(eh - el) >= 0.1:
                shape.altitude_high = round(eh, 1)
                shape.altitude_low = round(el, 1)
            else:
                shape.altitude = round((eh + el) / 2.0, 1)
            layout.shapes.append(shape)
            new_runway_polys.append(poly)

        # Drop any newly-emitted runway segment whose footprint is
        # contained inside an apt.dat / DSF pavement polygon that's
        # MUCH LARGER than the segment itself.  Such a polygon is
        # an apron enclosing the runway — the runway physically
        # merges with the surrounding pavement (e.g. CYXY runway 02
        # crosses the south apron at lat 60.7124).  Keeping a
        # separate runway rect there produces a visible rectangular
        # ribbon through the apron.  The junction polygon that
        # fills the apron will smoothly join the LAST surviving
        # runway segment's end at its shared corner.
        #
        # Detection: the segment is ≥ RUNWAY_INSIDE_APRON_FRAC
        # contained inside an apt.dat polygon whose area is
        # ≥ RUNWAY_APRON_AREA_RATIO times the segment area.  A
        # normal runway sits inside a runway-shaped apt.dat
        # polygon that's only marginally larger; an apron-merged
        # runway sits inside a polygon many times its size.
        if apron_candidates_m:
            from shapely.strtree import STRtree
            try:
                index = STRtree(apron_candidates_m)
            except Exception:
                index = None
            kept_shapes: List[BuiltShape] = []
            kept_polys: List[Polygon] = []
            # Track each dropped segment alongside the apron
            # candidate that contained it — used below to clip the
            # hole-fill merge so it can't bleed outside the apron.
            dropped_with_apron: List[Tuple[Polygon, Polygon]] = []
            n_dropped = 0
            for sh in layout.shapes:
                if sh.role != ROLE_RUNWAY:
                    kept_shapes.append(sh)
                    continue
                drop = False
                drop_apron: Optional[Polygon] = None
                if (sh.polygon is not None
                        and not sh.polygon.is_empty):
                    seg_area = sh.polygon.area
                    cand_iter = (index.query(sh.polygon)
                                 if index is not None
                                 else range(len(apron_candidates_m)))
                    for ci in cand_iter:
                        cand = apron_candidates_m[ci]
                        try:
                            if (cand.area
                                    < seg_area
                                    * RUNWAY_APRON_AREA_RATIO):
                                continue
                            inter = sh.polygon.intersection(cand)
                            if (inter.area / seg_area
                                    > RUNWAY_INSIDE_APRON_FRAC):
                                drop = True
                                drop_apron = cand
                                break
                        except Exception:
                            continue
                if drop:
                    n_dropped += 1
                    if (sh.polygon is not None
                            and not sh.polygon.is_empty
                            and drop_apron is not None):
                        dropped_with_apron.append(
                            (sh.polygon, drop_apron))
                    # Per user 2026-04-29 (CYXY runway 32L/16R
                    # ridge): preserve the dropped segment's
                    # altitude info on the layout so a later
                    # pass can imprint the runway's slope onto
                    # any apron-junction polygon that ends up
                    # covering this segment's footprint.  Without
                    # this, the surrounding apron junction's
                    # vertices are 4-5 m higher than the runway
                    # at the same XY, producing a "ridge across
                    # the runway" the user reported.
                    if (sh.polygon is not None
                            and not sh.polygon.is_empty
                            and (sh.altitude_high is not None
                                 or sh.altitude is not None)):
                        if not hasattr(
                                layout,
                                "_apron_merged_runway_drops"):
                            layout._apron_merged_runway_drops = []
                        layout._apron_merged_runway_drops.append(sh)
                    continue
                kept_shapes.append(sh)
                if sh.polygon is not None and not sh.polygon.is_empty:
                    kept_polys.append(sh.polygon)
            if n_dropped:
                try:
                    UI.vprint(1,
                        f"  [pav-builder] {icao}: dropped "
                        f"{n_dropped} runway segment(s) "
                        f"apron-merged.")
                except Exception:
                    pass
                layout.shapes = kept_shapes
                new_runway_polys = kept_polys

                # Per user 2026-04-28: dropping the apron-merged
                # runway segment leaves no hole — the residue
                # computation in ``build_airport_pavement`` already
                # excludes apron-merged regions from the runway-
                # subtraction (see ``_effective_runway_union``), so
                # the surrounding apron junction(s) cover the
                # runway segment's footprint naturally.  Nothing to
                # do here.
                _ = dropped_with_apron  # used only for the log line

        # Resolve runway-runway crossings: when two runway segments
        # overlap significantly (e.g. CYXY's crosswind 02/20 crossing
        # both 14R/32L and 14L/32R), drop both and emit a single
        # junction polygon at the union, with per-vertex altitudes
        # interpolated from the source segments.  Without this, the
        # downstream overlap-clip pass would clip one runway against
        # the other, leaving a 5-vertex shape that still carries
        # ``altitude_high`` / ``altitude_low`` tags — which X-Plane's
        # patch format only renders correctly on 4-corner rects.
        n_crossings = _resolve_runway_crossings(layout)
        if n_crossings:
            try:
                UI.vprint(1,
                    f"  [pav-builder] {icao}: resolved "
                    f"{n_crossings} runway crossing(s) into "
                    f"junction polygon(s).")
            except Exception:
                pass
            new_runway_polys = [
                s.polygon for s in layout.shapes
                if s.role == ROLE_RUNWAY
                and s.polygon is not None
                and not s.polygon.is_empty]
        # Per user 2026-04-30: bridge-segment insertion was
        # tried (option c) but produced worse results — the
        # inserted bridges overlap apron junctions whose
        # altitudes weren't synchronized, creating 10m+ range
        # junctions with 141 % worst-edge grades.  Function
        # ``_insert_runway_chain_bridges`` is left in place but
        # not called pending a different approach.
        if False:
            n_bridges = _insert_runway_chain_bridges(layout)

        # Segmented runway boundaries can drift sub-metre from the
        # original single-rect runway that junctions / rects were
        # built against, leaving tiny overlap slivers.  Subtract
        # the new runway union from junctions (and rects, defensively)
        # to eliminate them.
        if new_runway_polys:
            try:
                new_rwy_union = unary_union(
                    [p.buffer(0) for p in new_runway_polys]
                ).buffer(0)
            except Exception:
                new_rwy_union = None
            if new_rwy_union is not None and not new_rwy_union.is_empty:
                # Two clip regions:
                #
                #   `taxi_clip` — tiny buffer (0.05 m).  For taxi
                #   rects, only kills sub-metre sliver overlap with
                #   the runway while keeping the 4-corner rect
                #   shape intact.
                #
                #   `junction_clip` — 2.0 m outward buffer.  For
                #   junction polygons, shrinks them away from the
                #   runway edge by 2 m (user 2026-04-24: "taxiway
                #   shapes should stop just short of the runway").
                #   Prevents junction boundary vertices from
                #   landing mid-edge on a runway short edge, which
                #   X-Plane's mesh builder would interpret as
                #   splitting the runway's 4-corner slope rect
                #   and break the altitude_high/low rendering.
                try:
                    taxi_clip = new_rwy_union.buffer(0.05)
                except Exception:
                    taxi_clip = new_rwy_union
                try:
                    junction_clip = new_rwy_union.buffer(2.0)
                except Exception:
                    junction_clip = new_rwy_union
                # Rebuild layout.shapes in-place: when the clip
                # produces a MultiPolygon (e.g. a junction that
                # straddled the old runway ends up as two pieces
                # after the new segmented runway with overruns
                # replaces the old single rect), emit EVERY
                # sub-polygon above MIN_JUNCTION_AREA_M2 so no
                # pavement is lost (fixes the missing F/RW34R
                # gap junction, user 2026-04-24).
                from dataclasses import replace as _dc_replace
                new_shapes: List[BuiltShape] = []
                for shape in layout.shapes:
                    if shape.role in (ROLE_RUNWAY, ROLE_TERMINAL):
                        new_shapes.append(shape)
                        continue
                    # Clean sub-polygon before difference — apt.dat
                    # unions can produce polygons with self-kissing
                    # boundaries that trigger "side location" errors
                    # in shapely's overlay.
                    src = shape.polygon
                    if not src.is_valid:
                        try: src = src.buffer(0)
                        except Exception: pass
                    clip_region = (junction_clip
                                   if shape.role == ROLE_JUNCTION
                                   else taxi_clip)
                    try:
                        clipped = src.difference(clip_region)
                    except Exception:
                        try:
                            clipped = src.buffer(0).difference(clip_region)
                        except Exception:
                            new_shapes.append(shape)
                            continue
                    if clipped.is_empty:
                        continue
                    pieces = ([clipped]
                              if clipped.geom_type == "Polygon"
                              else list(getattr(clipped, "geoms", [])))
                    pieces = [p for p in pieces
                              if p.geom_type == "Polygon"
                              and p.area >= 50.0]
                    if not pieces:
                        continue
                    # Keep the shape metadata on the largest piece,
                    # emit any other pieces as new shapes with the
                    # same role/tags.  For junctions this splits
                    # the residue polygon; for rects this almost
                    # never splits (their snap keeps them whole).
                    pieces.sort(key=lambda g: -g.area)
                    shape.polygon = pieces[0]
                    new_shapes.append(shape)
                    for extra in pieces[1:]:
                        new_shapes.append(_dc_replace(
                            shape, polygon=extra, source_axis=None))
                layout.shapes = new_shapes

                # Per user 2026-04-28: junction polygon vertices
                # cannot land on a sloping rect's edge interior —
                # only on corners.  The 2 m runway-shrink difference
                # above can produce boundary intersection points 2 m
                # along a runway rect's edge; snap them to the
                # nearest corner.  Same helper used in
                # ``_resolve_runway_crossings``.
                sloping_rect_polys_for_snap = [
                    s.polygon for s in layout.shapes
                    if s.role in (ROLE_RUNWAY,
                                   ROLE_PRIMARY_PARALLEL,
                                   ROLE_SECONDARY_PARALLEL,
                                   ROLE_STUB,
                                   ROLE_CROSS_CONNECTOR)
                    and s.polygon is not None
                    and not s.polygon.is_empty]
                for shape in layout.shapes:
                    if shape.role != ROLE_JUNCTION:
                        continue
                    if (shape.polygon is None
                            or shape.polygon.is_empty):
                        continue
                    # Exclude this junction's own polygon from snap
                    # candidates (it isn't a sloping rect anyway,
                    # but be defensive).  Snap tolerance 5 m matches
                    # the runway-clip's 2 m buffer plus a small
                    # cushion for Shapely overlay precision.
                    snapped = _snap_polygon_vertices_to_rect_corners(
                        shape.polygon,
                        sloping_rect_polys_for_snap,
                        snap_tol_m=5.0)
                    if snapped is not None and not snapped.is_empty:
                        shape.polygon = snapped

    # ── Terminal pad elevations ─────────────────────────────────
    # Per user 2026-04-28: CIFP runway thresholds are the ONLY
    # truly authoritative elevations.  Everything else, including
    # the terminal altitude, should be derived to satisfy FAA
    # grade rules with the propagated runway / taxi / apron
    # values.  The previous DEM-median rule placed terminals on
    # naturally elevated ground (700.8 m at CYXY) when their
    # actual apron-side neighbours were 6-8 m lower; the apron-
    # pin then had to bridge that gap producing 4-7 % grade
    # violations.
    #
    # Rule (per user clarification):
    #   * Terminal altitude = MAX value that respects
    #     APRON_MAX_GRADE with every nearby HARD anchor (runway /
    #     taxi rect corner) at each terminal corner.  Specifically
    #     at each corner C, max allowed = MIN over nearby
    #     anchors A of (A.elev + APRON_MAX_GRADE × dist(C, A));
    #     terminal altitude = MIN over corners.
    #   * Floor at the highest nearby anchor (don't drop below
    #     the local terrain just because grade allows it).
    #   * Ceiling at the DEM-median (don't raise above natural
    #     ground).
    #   * Fall back to DEM-median if no anchors are available.
    runway_corner_pts: List[Tuple[float, float, float]] = []
    for s in layout.shapes:
        if s.role not in (
                ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                ROLE_CROSS_CONNECTOR):
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            r_coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if r_coords and r_coords[0] == r_coords[-1]:
            r_coords = r_coords[:-1]
        if len(r_coords) != 4:
            continue
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            per = [s.altitude_high, s.altitude_low,
                   s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            per = [float(s.altitude)] * 4
        else:
            continue
        for (x, y), a in zip(r_coords, per):
            runway_corner_pts.append(
                (float(x), float(y), float(a)))

    TERMINAL_NEIGHBOUR_RADIUS_M = 250.0
    INF = float("inf")
    for shape in layout.shapes:
        if shape.role != ROLE_TERMINAL:
            continue
        if shape.polygon is None or shape.polygon.is_empty:
            continue
        # DEM-median ceiling (legacy rule).
        dem_samples: List[float] = []
        try:
            t_corners = list(shape.polygon.exterior.coords)
            if t_corners and t_corners[0] == t_corners[-1]:
                t_corners = t_corners[:-1]
        except Exception:
            t_corners = []
        for x, y in [(shape.polygon.centroid.x,
                      shape.polygon.centroid.y)] + t_corners:
            lat, lon = m_to_ll(x, y)
            e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            if e is not None:
                dem_samples.append(e)
        if dem_samples:
            dem_samples.sort()
            dem_median = dem_samples[len(dem_samples) // 2]
        else:
            dem_median = None
        # Anchor-based max-allowable rule.
        new_alt: Optional[float] = None
        if runway_corner_pts and t_corners:
            per_corner_max: List[float] = []
            for cx, cy in t_corners:
                corner_max = INF
                hits = 0
                for ax, ay, ae in runway_corner_pts:
                    d = math.hypot(cx - ax, cy - ay)
                    if d > TERMINAL_NEIGHBOUR_RADIUS_M:
                        continue
                    allowed = ae + APRON_MAX_GRADE * d
                    if allowed < corner_max:
                        corner_max = allowed
                    hits += 1
                if hits > 0 and corner_max != INF:
                    per_corner_max.append(corner_max)
            if per_corner_max:
                # Terminal altitude must satisfy the MOST
                # restrictive corner.
                max_alt_from_anchors = min(per_corner_max)
                # Floor at the highest nearby anchor: don't
                # drop terminal below local pavement just
                # because grade allows it.
                anchor_floor = -INF
                for ax, ay, ae in runway_corner_pts:
                    # Only consider anchors with at least one
                    # corner within radius.
                    for cx, cy in t_corners:
                        if math.hypot(cx - ax,
                                      cy - ay) <= TERMINAL_NEIGHBOUR_RADIUS_M:
                            if ae > anchor_floor:
                                anchor_floor = ae
                            break
                candidate = max_alt_from_anchors
                if anchor_floor != -INF and candidate < anchor_floor:
                    candidate = anchor_floor
                # Ceiling at DEM-median (don't raise above
                # natural ground).
                if (dem_median is not None
                        and candidate > dem_median):
                    candidate = dem_median
                new_alt = round(float(candidate), 1)
        if new_alt is None and dem_median is not None:
            new_alt = round(float(dem_median), 1)
        if new_alt is None:
            continue
        shape.altitude = new_alt
        try:
            import sys as _sys
            dem_str = (f"{dem_median:.1f}" if dem_median is not None
                       else "n/a")
            UI.vprint(1,
                f"  [pav-builder] terminal({shape.ref or '?'}) "
                f"altitude {new_alt} m (max grade-compliant from "
                f"runway corners within "
                f"{TERMINAL_NEIGHBOUR_RADIUS_M:.0f} m, "
                f"DEM-median ceiling = {dem_str}).")
        except Exception:
            pass

    # ── Phase D: Geometric refinement + unified elevation solve ─
    # Run the geometric polygon-refinement passes (junction edge
    # push, triangulation, clamp, subdivide) interleaved with two
    # unified-Laplacian passes — see
    # ``_apply_geometric_finalization`` for the full sequence.
    # This replaces the previous bottom-up DEM-driven pipeline
    # (centerline graph + propagate_bounds + smooth_rate_of_change
    # + plateau snap + apron-pin + post-pin reconciliation).
    _apply_geometric_finalization(
        layout, icao, dem, tile_lat, tile_lon, m_to_ll)

    # ── Phase E: Diagnostics ────────────────────────────────────
    # NOTE: the within-shape grade WARN is intentionally NOT emitted
    # here.  ``build_airport_pavement`` (the caller) runs a chain of
    # post-elevation passes — ``_enforce_shared_vertices``,
    # ``_snap_junction_altitudes_to_rect_corners``,
    # ``_enforce_shared_vertex_altitudes``,
    # ``_smooth_within_junction_adjacent_pair_grade``, and a final
    # snap/agree round — that meaningfully change per-vertex
    # altitudes after this point.  Reporting here would surface the
    # MID-pipeline state (often 10×–100× worse than the final
    # output) and mislead.  The WARN is emitted from
    # ``build_airport_pavement`` after the smoother converges.


def _resample_node_altitudes_nn(
        new_poly: Polygon,
        old_open: List[Tuple[float, float]],
        old_alts_closed: Optional[List[float]],
        ) -> Optional[List[float]]:
    """Given a new polygon (post geometry edit) and the OLD ring's
    open-form coords + closed-form altitudes, return a fresh
    ``node_altitudes`` list (closed) for ``new_poly`` by nearest-
    neighbour sampling each new ring vertex against the old ring.

    Used wherever a polygon edit (boundary clip, buffer(0) repair,
    push-off, sliver merge, etc.) changes the vertex count and we
    would otherwise have to drop ``node_altitudes`` — a bare drop
    leaves the polygon with no elevation guidance, which X-Plane
    can render as a terrain spike (or, for very large boundary
    polygons, crash on load — see HECA's 1066-vertex airport
    boundary that lost altitudes during tunnel-clip).

    Returns None if the inputs are insufficient to resample.
    """
    if not old_alts_closed or not old_open:
        return None
    if new_poly is None or new_poly.is_empty:
        return None
    src_alts_open = (
        old_alts_closed[:-1]
        if (len(old_alts_closed) == len(old_open) + 1
            and old_alts_closed[0] == old_alts_closed[-1])
        else old_alts_closed[:len(old_open)])
    if not src_alts_open:
        return None
    try:
        new_open = list(new_poly.exterior.coords)
    except Exception:
        return None
    if new_open and new_open[0] == new_open[-1]:
        new_open = new_open[:-1]
    if not new_open:
        return None
    new_alts: List[float] = []
    for nx, ny in new_open:
        best_d2 = float("inf")
        best_a = src_alts_open[0]
        for k, (sx, sy) in enumerate(old_open):
            if k >= len(src_alts_open):
                break
            d2 = (nx - sx) ** 2 + (ny - sy) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_a = src_alts_open[k]
        new_alts.append(round(float(best_a), 1))
    return new_alts + [new_alts[0]]



def _apply_geometric_finalization(
        layout: "PavementLayout",
        icao: str,
        dem,
        tile_lat: int,
        tile_lon: int,
        m_to_ll,
        ) -> None:
    """Run the geometric polygon-refinement passes interleaved
    with two unified-solver passes:

      Phase 1: pre-solve geometry
        * Push junction vertices off taxi-rect edge interiors.
        * Triangulate junctions (each replaced with N-2 ear-clip
          triangles, initial node_altitudes from corner-elev map
          + DEM).

      Phase 2: first unified-solver pass
        Establishes a grade-compliant elevation field on the
        existing geometry.  These are the elevations the
        clamp + subdivide passes use to detect grade violations.

      Phase 3: clamp + subdivide based on real elevations
        * ``_clamp_junction_free_vertices``: tighten free
          boundary vertices to grade-comply with neighbours.
        * ``_subdivide_violating_junctions``: split any junction
          whose worst-pair grade exceeds the subdivision
          threshold along a perpendicular cut.
        * Re-push junction vertices off rect edges (subdivision
          can add new vertices on rect edges).
        * Re-clamp.

      Phase 4: second unified-solver pass
        Re-solves elevations on the refined geometry — produces
        the final FAA-compliant elevation field.
    """
    # Phase 1: pre-solve geometry.
    _push_junction_vertices_off_taxi_rect_edges(layout)
    # Triangulate junctions — initial node_altitudes come from
    # the corner-elev map + DEM fallback.  These are placeholders
    # for the first unified-solver pass below.
    _triangulate_junctions(
        layout, dem, tile_lat, tile_lon, m_to_ll)

    # Phase 2: first unified-solver pass (real elevations on the
    # current geometry — clamp + subdivide need these to detect
    # grade violations, not DEM-noisy fallbacks).
    _solve_pavement_elevations_unified(layout, icao)

    # Phase 3: clamp + subdivide based on the real elevations.
    clamp_geom = _build_clamp_geom_state(layout)
    for _ in range(8):
        n = _clamp_junction_free_vertices(layout, clamp_geom)
        if n == 0:
            break
    for _ in range(4):
        n = _subdivide_violating_junctions(layout)
        if n == 0:
            break
    # Subdivision may introduce vertices on rect edge interiors.
    _push_junction_vertices_off_taxi_rect_edges(layout)
    clamp_geom = _build_clamp_geom_state(layout)
    for _ in range(4):
        n = _clamp_junction_free_vertices(layout, clamp_geom)
        if n == 0:
            break

    # Phase 4: second unified-solver pass (final elevations on
    # refined geometry).
    _solve_pavement_elevations_unified(layout, icao)




def _solve_pavement_elevations_unified(
        layout: "PavementLayout",
        icao: str,
        max_iters: int = 1500,
        tol_m: float = 0.005,
        ) -> None:
    """Unified constrained-Laplacian elevation solver.

    Per user 2026-04-28: replaces the bottom-up DEM-driven
    pipeline (centerline graph + apron-pin + post-pin
    reconciliation + within-junction smoother) with a single
    top-down propagation:

      1. Hard anchors = every CIFP-derived RUNWAY corner (every
         segment in the runway chain has its altitude_high /
         altitude_low set from the FAA-compliant profile, and
         every corner of every segment counts as HARD).
      2. Build the unified pavement graph: every shape's polygon
         ring vertex becomes a node; ring edges + cross-shape
         shared-bucket edges connect them.
      3. Constrained Laplacian solve via Jacobi iteration with
         per-edge grade-cap projection:
           - Each iteration: every non-anchor node moves to the
             length-weighted average of its neighbours, then each
             edge with |Δelev|/length > max_grade pulls its
             endpoints toward each other (or moves the soft one
             toward the anchored one).
           - Per-shape role-specific grade caps:
               runway/runway: 1.5 % (already enforced by HARD)
               taxi rect: 1.5 %
               apron / junction: 1.0 %
               cross-shape: min of the two roles.
           - Terminal corners constrained to be FLAT (all corners
             of one terminal share a single value at every
             iteration).
      4. After convergence, apply the solved elevations back to
         every shape:
           - ROLE_RUNWAY: skip (HARD, already correct).
           - ROLE_PRIMARY_PARALLEL / SECONDARY_PARALLEL / STUB /
             CROSS_CONNECTOR: derive altitude_high/low from the
             rect's 4 corners (avg of corners 0,3 / corners 1,2).
           - ROLE_TERMINAL: avg of all corners → altitude.
           - ROLE_JUNCTION: per-vertex node_altitudes from corner
             elevations.

    The architectural advantage: every step of the previous
    pipeline solves a different subset of the constraints, and
    they sometimes disagree (which is why we kept finding new
    edge cases).  This single solve handles all constraints
    simultaneously.

    Performance: O((V + E) × I) where I = iteration count.
    Typically I ≈ graph_diameter² × log(1/tol) — for CYXY
    (~30-hop diameter) ≈ 900 iterations; for HECA (~150-hop)
    ≈ 22 500.  Each iteration is a single sweep of the edge
    list, ~1 µs/edge; CYXY ≈ 0.5 s, HECA ≈ 5 s.  Compare to the
    old pipeline at ~3.5 s and ~30 s respectively.
    """
    import time as _time
    t_start = _time.time()
    # ── Build node list ─────────────────────────────────────────
    pavement_roles = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_STUB,
        ROLE_CROSS_CONNECTOR, ROLE_TERMINAL, ROLE_JUNCTION,
    }
    bucket_to_idx: Dict[Tuple[int, int], int] = {}
    nodes: List[Tuple[float, float]] = []
    for s in layout.shapes:
        if s.role not in pavement_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for x, y in coords:
            b = _corner_elevation_bucket(x, y)
            if b not in bucket_to_idx:
                bucket_to_idx[b] = len(nodes)
                nodes.append((float(x), float(y)))
    n = len(nodes)
    if n == 0:
        return

    # ── Initial elevations + HARD anchor flags ──────────────────
    elev: List[float] = [0.0] * n
    is_hard: List[bool] = [False] * n
    have_initial: List[bool] = [False] * n

    # CIFP runway corners ⇒ HARD.
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            continue
        if (s.altitude_high is None or s.altitude_low is None):
            continue
        per = [s.altitude_high, s.altitude_low,
               s.altitude_low, s.altitude_high]
        for (x, y), a in zip(coords, per):
            b = _corner_elevation_bucket(x, y)
            if b in bucket_to_idx:
                idx = bucket_to_idx[b]
                # First-writer wins among runway corners (handles
                # adjacent segment shared corners — they should
                # already agree from the FAA profile, but pick
                # one canonically).
                if not is_hard[idx]:
                    elev[idx] = float(a)
                    is_hard[idx] = True
                    have_initial[idx] = True

    if not any(is_hard):
        # No CIFP anchors — can't run the unified solver.
        return

    # Seed soft nodes from their existing layout values when
    # available (gives the solver a warm start).
    for s in layout.shapes:
        if s.role not in pavement_roles or s.role == ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if (s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            per = [s.altitude_high, s.altitude_low,
                   s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            per = [float(s.altitude)] * len(coords)
        elif s.node_altitudes:
            per = [float(a) for a in
                   s.node_altitudes[:len(coords)]]
            if len(per) < len(coords):
                per = list(per) + [per[-1]] * (
                    len(coords) - len(per))
        else:
            continue
        for (x, y), a in zip(coords, per):
            b = _corner_elevation_bucket(x, y)
            if b not in bucket_to_idx:
                continue
            idx = bucket_to_idx[b]
            if is_hard[idx]:
                continue
            if not have_initial[idx]:
                elev[idx] = float(a)
                have_initial[idx] = True

    # Backfill any node still without an initial value via nearest
    # hard anchor's elevation (cheap pass).
    if any(not h for h in have_initial):
        hard_pts: List[Tuple[float, float, float]] = [
            (nodes[i][0], nodes[i][1], elev[i])
            for i in range(n) if is_hard[i]]
        for i in range(n):
            if have_initial[i]:
                continue
            x, y = nodes[i]
            best_d2 = float("inf")
            best_e = 0.0
            for hx, hy, he in hard_pts:
                d2 = (hx - x) * (hx - x) + (hy - y) * (hy - y)
                if d2 < best_d2:
                    best_d2 = d2
                    best_e = he
            elev[i] = best_e
            have_initial[i] = True

    # ── Build edge list with per-edge max grade ────────────────
    # Edge identified by sorted (u, v); max_grade = min over
    # contributing shapes' role caps.
    edge_grade: Dict[Tuple[int, int], float] = {}
    edge_length: Dict[Tuple[int, int], float] = {}

    def _role_grade(role: str) -> float:
        if role == ROLE_RUNWAY:
            return TAXI_MAX_GRADE  # 1.5 %, never tighter than this
        if role in (ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                     ROLE_STUB, ROLE_CROSS_CONNECTOR):
            return TAXI_MAX_GRADE
        # Terminal, junction, apron — apron rule.
        return APRON_MAX_GRADE

    # Per user 2026-04-29: in addition to ring-edge connectivity,
    # add a "spatial-pair" edge between every pair of vertices
    # WITHIN THE SAME SHAPE that are ≤ ``WITHIN_SHAPE_VIOLATION_
    # RADIUS_M`` apart in EUCLIDEAN distance.  The audit /
    # smoother both check spatial distance, not graph distance —
    # so without this, a junction with 100 ring vertices can have
    # vertex pair (i, j) that is 5 m apart spatially but 50 ring-
    # hops away.  The Laplacian's ring-edge cap then permits up
    # to 50 × per-edge cap of cumulative drift across the chain,
    # which the audit reports as a 100 % grade cliff.  Adding
    # the spatial pair as a direct edge constrains the pair the
    # same way the audit will check it.
    spatial_radius_m = WITHIN_SHAPE_VIOLATION_RADIUS_M
    spatial_radius2 = spatial_radius_m * spatial_radius_m
    for s in layout.shapes:
        if s.role not in pavement_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) < 2:
            continue
        gr = _role_grade(s.role)
        m = len(coords)
        # Pre-compute this shape's vertex node indices so we can
        # cheaply emit ring + spatial pairs.
        node_idx: List[Optional[int]] = []
        for x, y in coords:
            b = _corner_elevation_bucket(x, y)
            node_idx.append(bucket_to_idx.get(b))
        # Ring edges (i ↔ i+1).
        for i in range(m):
            ui = node_idx[i]
            uj = node_idx[(i + 1) % m]
            if ui is None or uj is None or ui == uj:
                continue
            x1, y1 = coords[i]
            x2, y2 = coords[(i + 1) % m]
            length = math.hypot(x2 - x1, y2 - y1)
            if length < 0.1:
                continue
            key = (ui, uj) if ui < uj else (uj, ui)
            cur_g = edge_grade.get(key, float("inf"))
            if gr < cur_g:
                edge_grade[key] = gr
            cur_l = edge_length.get(key, length)
            edge_length[key] = min(cur_l, length)
        # Spatial pairs (i, j) with j > i + 1 and Euclidean ≤
        # spatial_radius_m.  Skip pairs already connected as ring
        # edges (handled above) — the dict-min logic would just
        # repeat them.
        for i in range(m):
            xi, yi = coords[i]
            ui = node_idx[i]
            if ui is None:
                continue
            # Start at i+2 to skip the ring-adjacent pair that's
            # already added (and the i,i identity).  Treat the
            # ring-wrap pair (m-1, 0) as already covered too.
            for j in range(i + 2, m):
                # Skip the wrap-around ring edge.
                if i == 0 and j == m - 1:
                    continue
                uj = node_idx[j]
                if uj is None or ui == uj:
                    continue
                xj, yj = coords[j]
                dx = xj - xi
                dy = yj - yi
                d2 = dx * dx + dy * dy
                if d2 > spatial_radius2:
                    continue
                length = math.sqrt(d2)
                if length < 0.1:
                    continue
                key = (ui, uj) if ui < uj else (uj, ui)
                cur_g = edge_grade.get(key, float("inf"))
                if gr < cur_g:
                    edge_grade[key] = gr
                cur_l = edge_length.get(key, length)
                edge_length[key] = min(cur_l, length)

    # Adjacency for Jacobi step.
    adj: List[List[Tuple[int, float, float]]] = [[] for _ in range(n)]
    for (u, v), gr in edge_grade.items():
        L = edge_length[(u, v)]
        adj[u].append((v, L, gr))
        adj[v].append((u, L, gr))

    # ── Per-shape constraint groups ────────────────────────────
    # Terminal corners — flat constraint (all share the same value).
    terminal_groups: List[List[int]] = []
    for s in layout.shapes:
        if s.role != ROLE_TERMINAL:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        idxs: List[int] = []
        for x, y in coords:
            b = _corner_elevation_bucket(x, y)
            if b in bucket_to_idx:
                idxs.append(bucket_to_idx[b])
        if len(idxs) >= 2:
            terminal_groups.append(idxs)

    # ── Constrained Laplacian iteration ────────────────────────
    # Per user 2026-04-28: use a damped Jacobi + multi-sweep cap
    # projection to ensure the per-edge grade-cap wins over the
    # Jacobi neighbour pull.  Without damping, short edges with
    # strong external pull oscillate (Jacobi opens the gap, cap
    # closes it, Jacobi reopens it next iter).  With damping ~0.5
    # and multiple cap sweeps per Jacobi, the equilibrium settles
    # at the cap boundary as required.
    JACOBI_DAMPING = 0.5
    CAP_SWEEPS_PER_ITER = 5
    edge_list = list(edge_grade.keys())
    for it in range(max_iters):
        prev_elev = list(elev)
        # 1) Damped Jacobi neighbour-average step.
        new_elev = list(elev)
        for i in range(n):
            if is_hard[i]:
                continue
            nbrs = adj[i]
            if not nbrs:
                continue
            wsum = 0.0
            wvsum = 0.0
            for v, L, _g in nbrs:
                w = 1.0 / max(L, 0.1)
                wsum += w
                wvsum += w * elev[v]
            if wsum > 0:
                avg = wvsum / wsum
                new_elev[i] = (1.0 - JACOBI_DAMPING) * elev[i] \
                    + JACOBI_DAMPING * avg
        elev = new_elev
        # 2) Multi-sweep edge grade-cap projection.  Repeated
        # until either no edge violates or we hit the per-iter
        # sweep cap; this lets the cap propagate through chains
        # of edges in one outer step.
        for _sweep in range(CAP_SWEEPS_PER_ITER):
            any_proj = False
            for (u, v) in edge_list:
                L = edge_length[(u, v)]
                gr = edge_grade[(u, v)]
                diff = elev[u] - elev[v]
                cap = L * gr
                if abs(diff) <= cap:
                    continue
                excess = abs(diff) - cap
                sign = 1 if diff > 0 else -1
                if is_hard[u] and is_hard[v]:
                    continue
                if is_hard[u]:
                    elev[v] += sign * excess
                    any_proj = True
                elif is_hard[v]:
                    elev[u] -= sign * excess
                    any_proj = True
                else:
                    half = 0.5 * excess * sign
                    elev[u] -= half
                    elev[v] += half
                    any_proj = True
            if not any_proj:
                break
        # 3) Terminal flatness.
        for grp in terminal_groups:
            if not grp:
                continue
            free = [i for i in grp if not is_hard[i]]
            if not free:
                continue
            avg = sum(elev[i] for i in grp) / len(grp)
            for i in free:
                elev[i] = avg
        # 4) Convergence check.
        max_change = 0.0
        for i in range(n):
            if is_hard[i]:
                continue
            d = abs(prev_elev[i] - elev[i])
            if d > max_change:
                max_change = d
        if max_change < tol_m:
            break

    # ── Apply solved elevations back to layout shapes ──────────
    n_terms = 0
    n_rects = 0
    n_junctions = 0
    for s in layout.shapes:
        if s.role not in pavement_roles:
            continue
        if s.role == ROLE_RUNWAY:
            continue  # HARD, already correct
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        ring_closed = (coords and coords[0] == coords[-1])
        coords_open = coords[:-1] if ring_closed else coords
        corner_elevs: List[float] = []
        for x, y in coords_open:
            b = _corner_elevation_bucket(x, y)
            if b in bucket_to_idx:
                corner_elevs.append(elev[bucket_to_idx[b]])
            else:
                corner_elevs.append(float("nan"))
        if any(math.isnan(e) for e in corner_elevs):
            continue
        if s.role == ROLE_TERMINAL:
            avg = sum(corner_elevs) / len(corner_elevs)
            s.altitude = round(float(avg), 1)
            n_terms += 1
        elif s.role in (ROLE_PRIMARY_PARALLEL,
                         ROLE_SECONDARY_PARALLEL,
                         ROLE_STUB, ROLE_CROSS_CONNECTOR):
            if len(corner_elevs) == 4:
                # Sloping rect: corners 0,3 → high; 1,2 → low.
                hi = (corner_elevs[0] + corner_elevs[3]) / 2
                lo = (corner_elevs[1] + corner_elevs[2]) / 2
                if hi < lo:
                    hi, lo = lo, hi
                s.altitude_high = round(float(hi), 1)
                s.altitude_low = round(float(lo), 1)
                s.altitude = None
                n_rects += 1
        elif s.role == ROLE_JUNCTION:
            alts = [round(float(e), 1) for e in corner_elevs]
            if ring_closed:
                alts.append(alts[0])
            s.node_altitudes = alts
            n_junctions += 1

    elapsed = _time.time() - t_start
    try:
        UI.vprint(1,
            f"  [pav-builder] {icao}: unified Laplacian solver "
            f"converged in {it + 1}/{max_iters} iters "
            f"({elapsed:.2f} s); applied to "
            f"{n_terms} terminal(s), {n_rects} rect(s), "
            f"{n_junctions} junction(s).")
    except Exception:
        pass




def _smooth_within_junction_adjacent_pair_grade(
        layout: "PavementLayout",
        max_grade: float = 0.015,
        max_iters: int = 30,
        convergence_m: float = 0.01,
        pair_radius_m: float = 60.0,
        ) -> int:
    """For each junction polygon, iterate over EVERY vertex pair
    within ``pair_radius_m`` (not just immediate ring neighbours).
    When a pair exceeds ``max_grade`` (default 1.5 %, matching the
    FAA taxiway cap), nudge the un-anchored vertex(es) toward the
    grade band.

    Per user 2026-04-28: corner alignment between junctions and
    sloping rects is already enforced by
    ``_snap_junction_altitudes_to_rect_corners`` and
    ``_enforce_shared_vertex_altitudes``, but interior junction
    vertices can still disagree with their immediate neighbours by
    > 1.5 % over short distances (1217 such pairs at CYXY, worst
    17.6 %).  These come from the multi-source apron pin step
    (DEM plane-fit ↔ taxi-graph ↔ apron-pin DEM-clipped-to-ring-
    grade) where adjacent vertices pull from different sources.
    Iterating an adjacent-pair smoother after all the bucket-
    averaging passes have converged the SHARED-vertex constraints
    is the cleanest way to flatten the remaining within-polygon
    grade humps.

    Vertex anchoring rules:
      * A vertex is HARD if its bucket coincides with a sloping
        rect corner (runway / primary_parallel / secondary_parallel
        / stub / cross_connector) — its altitude must equal the
        rect's tag value and CANNOT be modified.
      * Otherwise the vertex is SOFT and may move.

    Pair handling:
      * Both HARD ⇒ skip (constraint unsolvable here).
      * One HARD, one SOFT ⇒ move the soft one to the boundary
        of the hard one's grade band.
      * Both SOFT ⇒ average (preserves volume).

    Returns the number of altitude entries adjusted (cumulative
    across iterations).
    """
    sloping_rect_roles = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_STUB,
        ROLE_CROSS_CONNECTOR,
    }
    # Hard-anchored buckets = sloping rect corner buckets.
    hard_buckets: set = set()
    for s in layout.shapes:
        if s.role not in sloping_rect_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for cx, cy in coords:
            hard_buckets.add(_corner_elevation_bucket(cx, cy))

    n_changed_total = 0
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        if not s.node_altitudes:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        n = len(coords)
        if n < 3 or len(s.node_altitudes) < n:
            continue
        alts = [float(a) for a in s.node_altitudes[:n]]
        # Pre-flag hard vertices.
        is_hard = [
            _corner_elevation_bucket(cx, cy) in hard_buckets
            for cx, cy in coords]
        # Pre-build pairs (i, j, distance) for every pair within
        # pair_radius_m.  For ring-adjacent pairs we always include
        # them (zero-distance edges between the same point are
        # filtered).  For non-adjacent pairs we filter by Euclidean
        # distance — distant pairs in the same polygon don't have
        # a meaningful grade constraint at airport scale.
        pair_radius2 = pair_radius_m * pair_radius_m
        pairs: List[Tuple[int, int, float]] = []
        for i in range(n):
            for j in range(i + 1, n):
                ax, ay = coords[i]
                bx, by = coords[j]
                dx = bx - ax
                dy = by - ay
                d2 = dx * dx + dy * dy
                if d2 > pair_radius2:
                    continue
                d = math.sqrt(d2)
                if d < 0.5:
                    continue
                pairs.append((i, j, d))
        if not pairs:
            continue
        for _it in range(max_iters):
            max_change = 0.0
            for i, j, d in pairs:
                de = abs(alts[i] - alts[j])
                grade = de / d
                if grade <= max_grade:
                    continue
                # Compute the maximum permitted |Δalt| at this
                # separation.
                max_de = max_grade * d
                hi = max(alts[i], alts[j])
                lo = min(alts[i], alts[j])
                hi_idx = i if alts[i] >= alts[j] else j
                lo_idx = j if hi_idx == i else i
                if is_hard[i] and is_hard[j]:
                    # Both anchored — leave alone (the constraint
                    # is unresolvable without moving runway/rect
                    # tags).
                    continue
                if is_hard[hi_idx] and not is_hard[lo_idx]:
                    # Lift the soft (low) vertex up to the band's
                    # lower edge.
                    new_lo = hi - max_de
                    if abs(alts[lo_idx] - new_lo) > convergence_m:
                        max_change = max(
                            max_change,
                            abs(alts[lo_idx] - new_lo))
                        alts[lo_idx] = new_lo
                        n_changed_total += 1
                elif is_hard[lo_idx] and not is_hard[hi_idx]:
                    # Drop the soft (high) vertex down to the
                    # band's upper edge.
                    new_hi = lo + max_de
                    if abs(alts[hi_idx] - new_hi) > convergence_m:
                        max_change = max(
                            max_change,
                            abs(alts[hi_idx] - new_hi))
                        alts[hi_idx] = new_hi
                        n_changed_total += 1
                else:
                    # Both soft — split the violation evenly.
                    avg = (alts[i] + alts[j]) / 2.0
                    excess = (de - max_de) / 2.0
                    new_hi = avg + max_de / 2.0
                    new_lo = avg - max_de / 2.0
                    if abs(alts[hi_idx] - new_hi) > convergence_m:
                        max_change = max(
                            max_change,
                            abs(alts[hi_idx] - new_hi))
                        alts[hi_idx] = new_hi
                        n_changed_total += 1
                    if abs(alts[lo_idx] - new_lo) > convergence_m:
                        max_change = max(
                            max_change,
                            abs(alts[lo_idx] - new_lo))
                        alts[lo_idx] = new_lo
                        n_changed_total += 1
            if max_change < convergence_m:
                break
        # Write back, preserving the closed-ring duplicate at the end.
        s.node_altitudes = [round(a, 1) for a in alts]
        if (len(s.node_altitudes) == n
                and s.polygon.exterior.coords[0]
                == s.polygon.exterior.coords[-1]):
            s.node_altitudes.append(s.node_altitudes[0])
    return n_changed_total


def _rederive_terminal_altitude_from_apron_neighbours(
        layout: "PavementLayout",
        sample_radius_m: float = 250.0,
        ) -> int:
    """Re-derive each terminal's altitude from the HARD-anchored
    sloping rect corners (runway / taxi) within
    ``sample_radius_m`` of the terminal — NOT from DEM under the
    terminal pad and NOT from apron-pinned vertices (which are
    self-referential to the old terminal altitude).

    Per user 2026-04-28: at CYXY the terminal's DEM-median was
    700.8 m, but the apron actually borders runway 02 (694 m) on
    the south and taxiway F (~692 m) on the north.  The terminal
    sits on naturally elevated ground that the apron pavement
    doesn't reach.  Forcing terminal to its DEM altitude forced
    the apron-pin step to bridge a 7 m gap over short distances,
    producing 4-7 % grade violations.

    Sampling only HARD-anchored sloping rect corners (whose
    elevations come from CIFP runway thresholds + chained
    grade-compliant interpolation, NOT from the terminal) breaks
    the self-reference.  The median of those anchors gives a
    terminal altitude consistent with the apron's actual
    elevation range — the apron's grade then flattens naturally
    because all four boundary classes (terminal, runway, taxi,
    junction) end up within a few metres of each other.

    Returns the number of terminal altitudes changed.
    """
    sloping_rect_roles = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_STUB,
        ROLE_CROSS_CONNECTOR,
    }
    hard_pts: List[Tuple[float, float, float]] = []
    for s in layout.shapes:
        if s.role not in sloping_rect_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            continue
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            per = [s.altitude_high, s.altitude_low,
                   s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            per = [float(s.altitude)] * 4
        else:
            continue
        for (x, y), a in zip(coords, per):
            hard_pts.append((float(x), float(y), float(a)))
    if not hard_pts:
        return 0

    n_changed = 0
    radius2 = sample_radius_m * sample_radius_m
    for s in layout.shapes:
        if s.role != ROLE_TERMINAL:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            t_boundary = s.polygon.boundary
        except Exception:
            t_boundary = None
        from shapely.geometry import Point as _P
        nearby: List[Tuple[float, float]] = []  # (distance, elev)
        for px, py, pa in hard_pts:
            try:
                if t_boundary is not None:
                    d = t_boundary.distance(_P(px, py))
                else:
                    d = math.hypot(
                        px - s.polygon.centroid.x,
                        py - s.polygon.centroid.y)
            except Exception:
                continue
            if d * d > radius2:
                continue
            nearby.append((d, pa))
        if not nearby:
            continue
        # Median of the elevations weighted by inverse distance —
        # closer hard anchors carry more weight.  Equivalent to
        # "what elevation do the closest hard anchors agree on?"
        # Take the median for robustness against outliers (e.g. a
        # nearby runway-32L corner at 706 m on the OTHER side of
        # the airport that happens to fall within the radius via
        # straight-line distance).
        nearby.sort()
        # Use just the closest 6 anchors to anchor on local
        # terrain rather than the airport-wide elevation range.
        closest = nearby[:6] if len(nearby) > 6 else nearby
        elevs = sorted(e for _d, e in closest)
        median = elevs[len(elevs) // 2]
        new_alt = round(float(median), 1)
        old_alt = s.altitude
        if old_alt is None or abs(old_alt - new_alt) >= 0.05:
            try:
                UI.vprint(1,
                    f"  [pav-builder] terminal({s.ref or '?'}) "
                    f"altitude {old_alt} → {new_alt} m "
                    f"(median of {len(closest)} closest hard "
                    f"anchors within {sample_radius_m:.0f} m; "
                    f"replaces DEM-median).")
            except Exception:
                pass
            s.altitude = new_alt
            n_changed += 1
    return n_changed


def _enforce_shared_vertex_altitudes(
        layout: "PavementLayout") -> int:
    """For every vertex bucket shared by ≥ 2 shapes, force every
    polygon's per-vertex altitude at that bucket to a single
    canonical value.

    Per user 2026-04-28: junctions sharing a boundary node MUST
    agree on its altitude, otherwise X-Plane renders a tear / step
    at the seam.  Subdivide / clamp / shared-vertex passes can
    leave neighbouring junctions with sub-metre disagreement at
    shared buckets even when the underlying mesh value is
    consistent.

    Policy: take the AVERAGE of the disagreeing altitudes.  Skip
    sloping rect tags (altitude_high / altitude_low / altitude) —
    those are the authoritative source and were already aligned by
    ``_snap_junction_altitudes_to_rect_corners``.

    Returns the number of altitude entries adjusted.
    """
    # Gather per-bucket altitude votes from junction polygons only.
    # (Sloped rect altitudes are tag-level; junctions emit per-vertex.)
    bucket_to_entries: Dict[Tuple[int, int],
                            List[Tuple[int, int, float]]] = {}
    for si, s in enumerate(layout.shapes):
        if s.role != ROLE_JUNCTION:
            continue
        if not s.node_altitudes:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for vi, (cx, cy) in enumerate(coords):
            if vi >= len(s.node_altitudes):
                break
            b = _corner_elevation_bucket(cx, cy)
            bucket_to_entries.setdefault(b, []).append(
                (si, vi, float(s.node_altitudes[vi])))
    n_changed = 0
    for b, entries in bucket_to_entries.items():
        if len(entries) < 2:
            continue
        alts = [e[2] for e in entries]
        spread = max(alts) - min(alts)
        if spread < 0.05:
            continue
        avg = round(sum(alts) / len(alts), 1)
        for si, vi, _e in entries:
            shape = layout.shapes[si]
            if abs(shape.node_altitudes[vi] - avg) < 0.05:
                continue
            shape.node_altitudes[vi] = avg
            # Maintain closed-ring invariant: last == first.
            if (vi == 0 and len(shape.node_altitudes) >= 2):
                shape.node_altitudes[-1] = avg
            n_changed += 1
    return n_changed


def _snap_junction_altitudes_to_rect_corners(
        layout: "PavementLayout",
        interior_proximity_m: float = 1.0,
        ) -> int:
    """For every junction polygon, snap any vertex whose bucket
    coincides with a runway / sloping-rect corner to that rect's
    corresponding altitude tag value (``altitude_high`` for HIGH
    corners 0,3; ``altitude_low`` for LOW corners 1,2; ``altitude``
    for flat shapes).

    Per user 2026-04-29 (CYXY runway-32L ridge): also snap
    junction vertices that lie INSIDE a sloping-rect footprint
    (within ``interior_proximity_m`` of the rect's interior or
    edge) to the rect's INTERPOLATED altitude at that point.
    Without this, a runway-crossing junction polygon whose
    vertices land inside a runway rect can override the rect's
    smooth slope with mesh-interpolated values that are 4-5 m
    off — visible as a ridge crossing the runway surface.

    Without this pass, the smoothing / subdivision / clamping
    passes can leave a junction's ``node_altitudes`` entry at a
    value derived from mesh interpolation rather than the rect's
    EMITTED altitude tag — resulting in a vertical step at the
    shared corner where the rect tag and the junction's per-vertex
    altitude disagree.

    Returns the number of altitude entries adjusted.
    """
    rwy_corner_alt: Dict[Tuple[int, int], float] = {}
    sloping_rect_roles_for_snap = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_STUB,
        ROLE_CROSS_CONNECTOR,
        # Per user 2026-04-28: terminal corners also propagate
        # to apron / junction vertices at the same bucket.  The
        # terminal is FLAT at ``s.altitude`` and the apron must
        # match at every shared corner — otherwise X-Plane
        # renders a step where the apron meets the terminal pad.
        ROLE_TERMINAL,
    }
    # Also collect the FULL sloping-rect shapes for interior-
    # snap (point-in-polygon + interpolated altitude).
    rect_shapes_for_interior: List[BuiltShape] = []
    for s in layout.shapes:
        if s.role not in sloping_rect_roles_for_snap:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        # Sloping rects (4-corner with altitude_high/low) — apply
        # per-corner altitudes.
        if (s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            for i, (cx, cy) in enumerate(coords):
                b = _corner_elevation_bucket(cx, cy)
                e = (s.altitude_high
                     if i in (0, 3)
                     else s.altitude_low)
                # First-writer wins — avoids different runway
                # segments at a shared corner disagreeing about
                # the canonical altitude.
                rwy_corner_alt.setdefault(b, float(e))
            rect_shapes_for_interior.append(s)
        elif s.altitude is not None:
            # Flat shapes (terminal pads, pre-elevation rects).
            # Any number of vertices; all share a single altitude.
            for (cx, cy) in coords:
                b = _corner_elevation_bucket(cx, cy)
                rwy_corner_alt.setdefault(b, float(s.altitude))
            rect_shapes_for_interior.append(s)
    if not rwy_corner_alt and not rect_shapes_for_interior:
        return 0
    # Spatial index for interior-snap probes.
    try:
        from shapely.strtree import STRtree as _STRtree
        if rect_shapes_for_interior:
            interior_tree = _STRtree(
                [s.polygon for s in rect_shapes_for_interior])
        else:
            interior_tree = None
    except Exception:
        interior_tree = None
    n_changed = 0
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        if not s.node_altitudes:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        # node_altitudes spans the closed ring; coords from
        # ``polygon.exterior.coords`` is also closed.  Walk the open
        # ring (drop closing repeat) and update by index.
        if coords and coords[0] == coords[-1]:
            coords_open = coords[:-1]
        else:
            coords_open = coords
        for i, (cx, cy) in enumerate(coords_open):
            if i >= len(s.node_altitudes):
                break
            # Pass 1: corner-bucket snap (prevents shared-corner
            # disagreement between junction and rect tags).
            b = _corner_elevation_bucket(cx, cy)
            target_e = rwy_corner_alt.get(b)
            if target_e is not None:
                if abs(s.node_altitudes[i] - target_e) >= 0.05:
                    s.node_altitudes[i] = round(target_e, 1)
                    n_changed += 1
                continue
            # Pass 2: interior snap — when the junction vertex
            # lies inside a sloping rect's polygon, snap to the
            # rect's interpolated altitude at that location.
            # Only adjust when the disagreement is > 0.5 m so we
            # don't undo the Laplacian solver's small refinements
            # at points that aren't truly inside a rect.
            if interior_tree is None:
                continue
            try:
                _Point = Point  # local alias
                pt = _Point(cx, cy)
                cands = interior_tree.query(pt)
            except Exception:
                cands = []
            best_e: Optional[float] = None
            best_d2 = float("inf")
            for hit in cands:
                ri = int(hit) if hasattr(hit, "__int__") else hit
                if not isinstance(ri, int):
                    continue
                rs = rect_shapes_for_interior[ri]
                try:
                    if rs.polygon.distance(pt) > interior_proximity_m:
                        continue
                except Exception:
                    continue
                e = _sample_runway_segment_elev(rs, cx, cy)
                if e is None:
                    continue
                # Pick the rect whose centroid is closest (deals
                # with overlapping rect candidates — runway +
                # parallel + stub at a complex junction).
                try:
                    rcx = rs.polygon.centroid.x
                    rcy = rs.polygon.centroid.y
                    d2 = (rcx - cx) ** 2 + (rcy - cy) ** 2
                except Exception:
                    d2 = 0.0
                if d2 < best_d2:
                    best_d2 = d2
                    best_e = float(e)
            if best_e is None:
                continue
            if abs(s.node_altitudes[i] - best_e) >= 0.5:
                s.node_altitudes[i] = round(best_e, 1)
                n_changed += 1
        # Maintain closed-ring invariant: last == first.
        if (s.node_altitudes
                and len(s.node_altitudes) >= 2
                and s.node_altitudes[0] != s.node_altitudes[-1]):
            s.node_altitudes[-1] = s.node_altitudes[0]
    return n_changed


def _re_emit_apron_merged_runway_segments(
        layout: "PavementLayout",
        ) -> int:
    """Re-emit each apron-merged-and-dropped runway segment back
    into ``layout.shapes`` AFTER the elevation pipeline has
    finished assigning altitudes to surrounding apron junctions.
    Subtracts the re-emitted footprint from any overlapping
    junction polygon so the two don't conflict, with NN-resample
    of the junction's per-vertex altitudes after the clip.

    Per user 2026-04-29 (CYXY runway 32L/16R ridge): the
    builder drops runway segments that lie inside a much-larger
    apron polygon to avoid emitting a visible rectangular
    ribbon.  The surrounding apron junction then takes over the
    surface in that area — but the junction's altitudes come
    from the Laplacian solver pinning at non-runway corners, so
    the surface at the dropped-segment's footprint can sit 4-5 m
    above the runway's CIFP profile.  Reading along the runway,
    the rendered surface dips into the runway segment, rises
    over the apron-junction-covered void, then dips back into
    the next segment — the user's "ridge across runway".

    Re-emitting the dropped segment with its preserved
    ``altitude``/``altitude_high``/``altitude_low`` puts a real
    runway-altitude plate back into the output.  The
    surrounding apron junction is clipped (subtracted) so it no
    longer covers the runway footprint.  The runway surface is
    now continuous at its profile altitude through the absorbed
    area.

    Returns the number of segments re-emitted.
    """
    drops = list(getattr(
        layout, "_apron_merged_runway_drops", []) or [])
    if not drops:
        return 0
    try:
        from shapely.strtree import STRtree as _STRtree
        # Index the junction polygons so we know which to clip.
        jct_idxs = [i for i, s in enumerate(layout.shapes)
                     if s.role == ROLE_JUNCTION
                     and s.polygon is not None
                     and not s.polygon.is_empty]
        if jct_idxs:
            index = _STRtree([layout.shapes[i].polygon
                               for i in jct_idxs])
        else:
            index = None
    except Exception:
        index = None
        jct_idxs = []
    n_emitted = 0
    for seg in drops:
        if seg.polygon is None or seg.polygon.is_empty:
            continue
        # Subtract the segment's footprint from every junction
        # whose polygon overlaps it.  Use the raw segment polygon
        # (no buffer) so the clip is exactly the runway shape.
        if index is not None:
            try:
                cands = index.query(seg.polygon)
            except Exception:
                cands = []
            for hit in cands:
                ji = int(hit) if hasattr(hit, "__int__") else hit
                if not isinstance(ji, int):
                    continue
                shape_i = jct_idxs[ji]
                target = layout.shapes[shape_i]
                if (target.polygon is None
                        or target.polygon.is_empty):
                    continue
                try:
                    inter_area = target.polygon.intersection(
                        seg.polygon).area
                except Exception:
                    continue
                if inter_area < 1.0:
                    continue
                # Capture old ring + altitudes BEFORE clip.
                try:
                    _old_ring = list(
                        target.polygon.exterior.coords)
                except Exception:
                    _old_ring = []
                if _old_ring and _old_ring[0] == _old_ring[-1]:
                    _old_ring = _old_ring[:-1]
                _old_alts = (list(target.node_altitudes)
                              if target.node_altitudes else None)
                try:
                    new_poly = target.polygon.difference(
                        seg.polygon)
                except Exception:
                    continue
                if (new_poly.is_empty
                        or new_poly.geom_type
                        not in ("Polygon", "MultiPolygon")):
                    continue
                if new_poly.geom_type == "MultiPolygon":
                    new_poly = max(
                        (g for g in new_poly.geoms
                          if g.geom_type == "Polygon"),
                        key=lambda g: g.area, default=None)
                    if new_poly is None or new_poly.is_empty:
                        continue
                target.polygon = new_poly
                resampled = _resample_node_altitudes_nn(
                    new_poly, _old_ring, _old_alts)
                if resampled is not None:
                    target.node_altitudes = resampled
        # Re-add the runway segment.
        layout.shapes.append(seg)
        n_emitted += 1
    return n_emitted


def _latlon_to_m_local(lat: float, lon: float,
                       lat0: float, lon0: float, cos0: float
                       ) -> Tuple[float, float]:
    x = math.radians(lon - lon0) * R_EARTH * cos0
    y = math.radians(lat - lat0) * R_EARTH
    return x, y


def _orient_rect_for_altitude(shape: "BuiltShape",
                              p1: Tuple[float, float],
                              p2: Tuple[float, float],
                              e1: float, e2: float) -> None:
    """Rewrite a 4-corner rect polygon's ring in the X-Plane
    patch convention:

        [n0 high-left, n1 low-left, n2 low-right, n3 high-right]

    where "high" is whichever of ``p1`` / ``p2`` has the larger
    elevation (``e1`` / ``e2``) and "left" / "right" are
    relative to the high→low axis direction.  The way's short
    edges are then:

        way[-2:] = [n3, n0]  = altitude_high short edge
        way[1:3] = [n1, n2]  = altitude_low  short edge

    Earlier pipeline stages (corner snap to pav vertices,
    shared-vertex enforcement) may have permuted the polygon's
    ring order, so this function re-derives the ordering from
    the 4 raw corner positions by classifying each by nearest
    axis endpoint and by left/right of the axis perpendicular.
    Non-4-corner polygons are left alone.
    """
    try:
        coords = list(shape.polygon.exterior.coords)
    except Exception:
        return
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return
    # Classify each corner by nearest axis endpoint.
    p1_corners: List[Tuple[float, float]] = []
    p2_corners: List[Tuple[float, float]] = []
    for c in coords:
        d1 = (c[0] - p1[0]) ** 2 + (c[1] - p1[1]) ** 2
        d2 = (c[0] - p2[0]) ** 2 + (c[1] - p2[1]) ** 2
        (p1_corners if d1 <= d2 else p2_corners).append(c)
    if len(p1_corners) != 2 or len(p2_corners) != 2:
        return
    # Perpendicular for the HIGH→LOW walk direction.
    if e1 >= e2:
        hi_p, lo_p = p1, p2
        hi_corners, lo_corners = p1_corners, p2_corners
    else:
        hi_p, lo_p = p2, p1
        hi_corners, lo_corners = p2_corners, p1_corners
    dx = lo_p[0] - hi_p[0]
    dy = lo_p[1] - hi_p[1]
    ax_len = math.hypot(dx, dy)
    if ax_len < 0.1:
        return
    ux, uy = dx / ax_len, dy / ax_len
    # Left perp when walking high→low = (-uy, ux).
    def _side(c, ref):
        """Positive = left of axis from HIGH end looking LOW."""
        rx, ry = c[0] - ref[0], c[1] - ref[1]
        return rx * (-uy) + ry * ux
    # Sort each endpoint's 2 corners: left first.
    hi_corners = sorted(hi_corners, key=lambda c: -_side(c, hi_p))
    lo_corners = sorted(lo_corners, key=lambda c: -_side(c, lo_p))
    hi_left, hi_right = hi_corners[0], hi_corners[1]
    lo_left, lo_right = lo_corners[0], lo_corners[1]
    # Build ring in legacy convention: high-left → low-left → low-right → high-right.
    new_ring = [hi_left, lo_left, lo_right, hi_right, hi_left]
    try:
        new_poly = Polygon(new_ring, list(shape.polygon.interiors))
        if not new_poly.is_valid:
            new_poly = new_poly.buffer(0)
        if (new_poly.geom_type == "Polygon"
                and not new_poly.is_empty):
            shape.polygon = new_poly
    except Exception:
        pass



# ──────────────────────────────────────────────────────────────────
# Junction-polygon decomposition + densification
# (re-exported from O4_Pavement_Junctions)
# ──────────────────────────────────────────────────────────────────
from .pavement.junctions import (
    _decompose_polygon_with_holes,
    _densify_long_boundary_edges,
    _drop_colinear_boundary_vertices,
    _drop_sliver_corners,
    _merge_thin_decomposed_pieces,
    _polygon_area,
    _polygon_min_thickness,
    _splice_holes,
    _splice_one_hole,
)



def _planar_fit(ring: List[Tuple[float, float]],
                elev: List[float]
                ) -> Optional[Tuple[float, float, float, List[float]]]:
    """Fit a plane ``z = a*x + b*y + c`` to (x, y, z) by least
    squares and return ``(a, b, c, per-vertex residuals)``.  The
    slope magnitude is ``sqrt(a² + b²)`` (rise per metre of horizontal
    travel — directly comparable to ``TAXI_MAX_GRADE``).
    Returns None if the fit is degenerate (colinear xy).
    """
    n = len(ring)
    if n < 3 or len(elev) != n:
        return None
    sxx = sxy = sxc = syy = syc = scc = 0.0
    sxz = syz = szc = 0.0
    for (x, y), z in zip(ring, elev):
        sxx += x * x
        sxy += x * y
        sxc += x
        syy += y * y
        syc += y
        scc += 1.0
        sxz += x * z
        syz += y * z
        szc += z
    det = (sxx * (syy * scc - syc * syc)
           - sxy * (sxy * scc - syc * sxc)
           + sxc * (sxy * syc - syy * sxc))
    if abs(det) < 1e-9:
        return None
    det_a = (sxz * (syy * scc - syc * syc)
             - sxy * (syz * scc - syc * szc)
             + sxc * (syz * syc - syy * szc))
    det_b = (sxx * (syz * scc - syc * szc)
             - sxz * (sxy * scc - syc * sxc)
             + sxc * (sxy * szc - syz * sxc))
    det_c = (sxx * (syy * szc - syz * syc)
             - sxy * (sxy * szc - syz * sxc)
             + sxz * (sxy * syc - syy * sxc))
    a = det_a / det
    b = det_b / det
    c = det_c / det
    residuals = [abs(z - (a * x + b * y + c))
                 for (x, y), z in zip(ring, elev)]
    return (a, b, c, residuals)


def _planar_fit_residuals(ring: List[Tuple[float, float]],
                          elev: List[float]
                          ) -> Optional[List[float]]:
    """Backwards-compat wrapper: residuals only."""
    f = _planar_fit(ring, elev)
    return None if f is None else f[3]


def _match_elev(rx: float, ry: float,
                ring: List[Tuple[float, float]],
                elev: List[float]) -> float:
    """Find the elevation in ``elev`` whose corresponding ring
    vertex is closest to (rx, ry).  Used to map shapely-emitted
    closed-ring coords back to our smoothed elevation array."""
    best_e = elev[0]
    best_d2 = float("inf")
    for (x, y), e in zip(ring, elev):
        d2 = (rx - x) * (rx - x) + (ry - y) * (ry - y)
        if d2 < best_d2:
            best_d2 = d2
            best_e = e
    return best_e




# ── Mode B: per-polygon 2D elevation grid ─────────────────────────
#
# Build a 5 m grid covering a non-rect pavement polygon's bbox,
# pin cells nearest each "hard anchor" (rect/runway/terminal corner
# elevations) to those anchor values, initialize free cells to DEM
# clipped into the local feasibility cone (anchor ± dist × 1.5 %),
# then Laplacian-smooth + grade-cap every adjacent cell pair until
# the field is stable.  Returns a sampler that bilinearly
# interpolates the smoothed grid at any (x, y) the caller asks
# about.  The polygon's boundary vertices then take their
# elevations from this single shared field — which is what
# guarantees within-shape grade compliance for the polygon.
#
# Per the elevation-field plan (2026-04-26) this replaces per-
# vertex independent elevation derivation.  Rect / runway / terminal
# corner anchors stay immutable across smoothing iterations.



# ──────────────────────────────────────────────────────────────────
# 2D polygon-grid smoothing (extracted to auto_patch.elevation_smoothing)
# ──────────────────────────────────────────────────────────────────
from .elevation_smoothing import _smooth_polygon_grid

def _corner_elevation_bucket(x: float, y: float,
                             tol: float = SHARED_VERTEX_TOL_M
                             ) -> Tuple[int, int]:
    """Quantize a meter-space point to a vertex-bucket key.

    Two coordinates within ``tol`` metres of each other (default
    ``SHARED_VERTEX_TOL_M``) hash to the same bucket — used to
    treat vertices on adjacent shapes that should share a node id
    as a single logical point even when their floating-point
    coordinates differ slightly.
    """
    return (int(round(x / tol)), int(round(y / tol)))


def _corner_elev_map(layout: "PavementLayout"
                     ) -> Dict[Tuple[int, int], float]:
    """Return a bucket-keyed elevation lookup for every corner of
    every elevation-bearing non-junction shape.  Uses the same
    bucket size as ``to_osm`` so junction vertices that share a
    node id with a corner will hit the same bucket.

    For sloped rect/runway shapes (altitude_high+altitude_low),
    ring indices 0,3 are the HIGH short edge and 1,2 are the LOW
    short edge — see ``_orient_rect_for_altitude``.  For flat
    polygons (altitude only), every corner gets the single value.
    """
    rect_like_roles = {ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL,
                       ROLE_STUB, ROLE_CROSS_CONNECTOR}
    out: Dict[Tuple[int, int], float] = {}
    for s in layout.shapes:
        if s.role == ROLE_JUNCTION:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        # Sloped 4-corner rect/runway in patch convention.
        if (s.role in rect_like_roles
                and s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            elevs = [s.altitude_high, s.altitude_low,
                     s.altitude_low, s.altitude_high]
            for (cx, cy), e in zip(coords, elevs):
                out.setdefault(
                    _corner_elevation_bucket(cx, cy), float(e))
            continue
        # Flat polygon (terminal, flat rect, or flat runway).
        if s.altitude is not None:
            for (cx, cy) in coords:
                out.setdefault(
                    _corner_elevation_bucket(cx, cy),
                    float(s.altitude))
    return out



# ──────────────────────────────────────────────────────────────────
# Per-vertex junction altitude assignment
# (extracted to auto_patch.triangulation)
# ──────────────────────────────────────────────────────────────────
from .triangulation import _triangulate_junctions

def _build_clamp_geom_state(
        layout: "PavementLayout"
        ) -> "Optional[Tuple[List, List, Dict, set]]":
    """Build the GEOMETRY-ONLY state used by
    :func:`_clamp_junction_free_vertices`.

    The clamp's spatial grid + shared-bucket set are functions of
    polygon geometry alone; only the per-edge elevations change
    between iterations.  Hoisting this build out of the per-call
    body lets the outer fixed-point loop reuse it without
    rebuilding (8+ × cost saved at HECA).

    Returns ``(edge_geom, edge_endpoints, grid, shared_buckets)``
    where:

    * ``edge_geom``:    list of ``(shape_idx, vi_a, vi_b, ax, ay,
                        bx, by)`` — vertex indices index into
                        ``layout.shapes[shape_idx]``'s exterior
                        coords (closed ring's prefix, i.e. the last
                        repeat dropped) so the caller can read
                        current elevations per call.
    * ``edge_endpoints``: bucket key pair per edge (for
                          incident-edge skip).
    * ``grid``:         spatial bucket → list of edge indices.
    * ``shared_buckets``: buckets touched by ≥ 2 shapes.

    Returns None when ``layout.anchor`` is unset.
    """
    if layout.anchor is None:
        return None

    edge_geom: List[Tuple[int, int, int,
                          float, float, float, float]] = []
    edge_endpoints: List[Tuple[Tuple[int, int],
                                Tuple[int, int]]] = []
    bucket_count: Dict[Tuple[int, int], int] = {}
    for si, s in enumerate(layout.shapes):
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        for (cx, cy) in coords:
            bucket = _corner_elevation_bucket(cx, cy)
            bucket_count[bucket] = bucket_count.get(bucket, 0) + 1
        n = len(coords)
        for i in range(n):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % n]
            edge_geom.append((si, i, (i + 1) % n,
                               ax, ay, bx, by))
            edge_endpoints.append((
                _corner_elevation_bucket(ax, ay),
                _corner_elevation_bucket(bx, by)))

    shared_buckets = {b for b, c in bucket_count.items() if c >= 2}

    grid: Dict[Tuple[int, int], List[int]] = {}
    cell = NEIGHBOUR_CLAMP_RADIUS_M
    for ei, (_, _, _, ax, ay, bx, by) in enumerate(edge_geom):
        x0, x1 = (ax, bx) if ax <= bx else (bx, ax)
        y0, y1 = (ay, by) if ay <= by else (by, ay)
        ix0 = int(math.floor((x0 - cell) / cell))
        ix1 = int(math.floor((x1 + cell) / cell))
        iy0 = int(math.floor((y0 - cell) / cell))
        iy1 = int(math.floor((y1 + cell) / cell))
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                grid.setdefault((ix, iy), []).append(ei)

    return (edge_geom, edge_endpoints, grid, shared_buckets)


def _clamp_junction_free_vertices(
        layout: "PavementLayout",
        geom_state: "Optional[Tuple[List, List, Dict, set]]" = None,
        ) -> int:
    """Per-junction free-vertex clamp using every nearby shape
    boundary as a soft anchor (Layer 2).

    Returns the number of free-vertex elevations that changed
    (informational).  Mutates each junction's
    ``node_altitudes`` and ``altitude`` in place.

    When ``geom_state`` is supplied (from
    :func:`_build_clamp_geom_state`), the geometry-only spatial
    grid + shared-bucket set is reused across iterations of the
    outer fixed-point loop.  Falls back to building it on first
    call when the caller doesn't.
    """
    if layout.anchor is None:
        return 0
    if geom_state is None:
        geom_state = _build_clamp_geom_state(layout)
        if geom_state is None:
            return 0
    edge_geom, edge_endpoints, grid, shared_buckets = geom_state
    cell = NEIGHBOUR_CLAMP_RADIUS_M

    # Read current per-shape elevation arrays once per call so the
    # inner clamp loop can look up edge endpoint elevations by
    # (shape_idx, vertex_idx) without re-parsing shapes per edge.
    shape_elevs: Dict[int, List[float]] = {}
    for si, s in enumerate(layout.shapes):
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords_n = len(s.polygon.exterior.coords)
            if coords_n > 0 and (
                s.polygon.exterior.coords[0]
                == s.polygon.exterior.coords[-1]
            ):
                coords_n -= 1
        except Exception:
            continue
        if coords_n <= 0:
            continue
        if s.altitude is not None:
            shape_elevs[si] = [float(s.altitude)] * coords_n
        elif (s.altitude_high is not None
              and s.altitude_low is not None
              and coords_n == 4):
            shape_elevs[si] = [float(s.altitude_high),
                                float(s.altitude_low),
                                float(s.altitude_low),
                                float(s.altitude_high)]
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == coords_n + 1:
                na = na[:-1]
            if len(na) == coords_n:
                shape_elevs[si] = [float(e) for e in na]

    # Reconstruct the boundary_edges layout used downstream
    # ``(shape_idx, ax, ay, bx, by, ea, eb)`` with current elevs.
    boundary_edges: List[Tuple[int, float, float, float, float,
                                float, float]] = []
    boundary_edges_append = boundary_edges.append
    for (si, vi_a, vi_b, ax, ay, bx, by) in edge_geom:
        elevs = shape_elevs.get(si)
        if elevs is None:
            # Skip shapes that didn't yield a valid elev array.
            boundary_edges_append((si, ax, ay, bx, by, 0.0, 0.0))
            continue
        boundary_edges_append((si, ax, ay, bx, by,
                                elevs[vi_a], elevs[vi_b]))

    rect_like_roles = {ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                       ROLE_CROSS_CONNECTOR, ROLE_TERMINAL}

    n_changed = 0
    for si, s in enumerate(layout.shapes):
        if s.role != ROLE_JUNCTION:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        # Build per-vertex elevations and FREE-vs-anchored flags.
        if s.altitude is not None:
            elevs = [float(s.altitude)] * len(coords)
            uniform_alt = True
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == len(coords) + 1:
                na = na[:-1]
            if len(na) != len(coords):
                continue
            elevs = [float(e) for e in na]
            uniform_alt = False
        else:
            continue
        # A vertex is ANCHORED if its bucket is shared with another
        # shape (rect/runway/terminal or another junction with the
        # same coord).  We check via the shared_buckets set.
        is_anchor = []
        for (cx, cy) in coords:
            b = _corner_elevation_bucket(cx, cy)
            is_anchor.append(b in shared_buckets)
        # Now clamp each FREE vertex.
        new_elevs = list(elevs)
        for vi, (cx, cy) in enumerate(coords):
            if is_anchor[vi]:
                continue
            ix = int(math.floor(cx / cell))
            iy = int(math.floor(cy / cell))
            lo_v = float("-inf")
            hi_v = float("inf")
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    bucket = (ix + dx, iy + dy)
                    if bucket not in grid:
                        continue
                    v_bucket = _corner_elevation_bucket(cx, cy)
                    for ei in grid[bucket]:
                        (s_other, ax, ay, bx, by,
                         ea, eb) = boundary_edges[ei]
                        # Skip ONLY edges incident to THIS vertex —
                        # they trivially equal the vertex's own
                        # elevation and would lock it in place.  All
                        # other same-shape boundary edges DO
                        # constrain it (the polygon's far-side
                        # vertices, narrow-waist near-touches, etc.).
                        ek0, ek1 = edge_endpoints[ei]
                        if v_bucket == ek0 or v_bucket == ek1:
                            continue
                        edx = bx - ax
                        edy = by - ay
                        seg2 = edx * edx + edy * edy
                        if seg2 < 0.04:
                            continue
                        t = ((cx - ax) * edx + (cy - ay) * edy) / seg2
                        if t < 0.0:
                            t = 0.0
                        elif t > 1.0:
                            t = 1.0
                        ccx = ax + t * edx
                        ccy = ay + t * edy
                        d = math.hypot(cx - ccx, cy - ccy)
                        if d > NEIGHBOUR_CLAMP_RADIUS_M:
                            continue
                        e_at = ea + t * (eb - ea)
                        band = max(d, 0.01) * TAXI_MAX_GRADE
                        if e_at - band > lo_v:
                            lo_v = e_at - band
                        if e_at + band < hi_v:
                            hi_v = e_at + band
            if lo_v > hi_v:
                # Conflicting nearby anchors — pick midpoint as a
                # least-squares-style compromise.  Layer 1 will
                # eventually reconcile the source corners.
                target = 0.5 * (lo_v + hi_v)
            else:
                target = elevs[vi]
                if target < lo_v:
                    target = lo_v
                if target > hi_v:
                    target = hi_v
            target = round(target, 1)
            if abs(target - elevs[vi]) > 0.05:
                new_elevs[vi] = target
                n_changed += 1
        # Persist changes.
        if uniform_alt:
            # Was a single altitude; if any free vertex shifted, we
            # have to switch to per-vertex node_altitudes.
            if any(abs(new_elevs[i] - elevs[i]) > 0.05
                   for i in range(len(elevs))):
                elev_range = max(new_elevs) - min(new_elevs)
                if elev_range < 0.05:
                    s.altitude = round(
                        sum(new_elevs) / len(new_elevs), 1)
                else:
                    s.altitude = None
                    closed = list(new_elevs) + [new_elevs[0]]
                    s.node_altitudes = closed
        else:
            # Already per-vertex.  Update node_altitudes (ring +
            # closing repeat).
            closed = list(new_elevs) + [new_elevs[0]]
            s.node_altitudes = closed
            # If everything collapsed to one value, switch back to
            # flat altitude.
            elev_range = max(new_elevs) - min(new_elevs)
            if elev_range < 0.05:
                s.altitude = round(
                    sum(new_elevs) / len(new_elevs), 1)
                s.node_altitudes = None
    return n_changed


SUBDIVIDE_VIOLATION_GRADE = 0.10   # 10 % — only attempt
                                    # subdivision when the worst
                                    # vertex pair exceeds this.
                                    # Below 10 % the rendered cliff
                                    # is small (< 1.5 m over 15 m)
                                    # and not worth the polygon-
                                    # split overhead.
SUBDIVIDE_MAX_PAIR_DIST_M = 60.0   # only consider pairs within
                                    # this radius — same as
                                    # check_grade's
                                    # WITHIN_SHAPE_MAX_PAIR_DIST_M
                                    # (Triangle4XP-plausible edge).
SUBDIVIDE_MIN_AREA_M2 = 5.0        # don't emit sub-polygons
                                    # smaller than this — they'd
                                    # become slivers and re-trigger
                                    # the sliver-corner safety net.


SUBDIVIDE_SNAP_RADIUS_M = 5.0      # snap new cut-line vertices
                                    # to existing ring vertices when
                                    # within this radius.  Without
                                    # snapping, the cut introduces
                                    # 0.5-2 m new vertices very
                                    # close to existing ones; the
                                    # interpolated-vs-original
                                    # elevation mismatch over those
                                    # tiny distances produces huge
                                    # spurious grade percentages
                                    # (worse than the original
                                    # violation we were trying to fix).


def _subdivide_violating_junctions(layout: "PavementLayout") -> int:
    """Split junction polygons whose worst within-shape vertex pair
    exceeds ``SUBDIVIDE_VIOLATION_GRADE`` along a perpendicular
    cut through the midpoint of the violating pair.

    Cut-vertex placement: new vertices created where the cut line
    intersects the polygon boundary are SNAPPED to existing ring
    vertices within ``SUBDIVIDE_SNAP_RADIUS_M``.  Without snapping,
    the cut creates ring-adjacent vertex pairs spaced 0.5-2 m apart
    whose interpolated-vs-original elevations differ by small
    amounts — registering as huge grade percentages
    (e.g. 0.4 m / 0.5 m = 73.9 %) that are WORSE than the original
    violation we were trying to fix.

    Validation: a sub-polygon is only accepted when its own worst
    within-shape vertex pair is BETTER (smaller grade) than the
    parent's worst pair.  If the cut would produce a sub-polygon
    that's MORE violating, the original is kept and the cut is
    abandoned — prevents the iterative subdivision from making
    things worse.

    Returns the number of polygons that were subdivided
    (informational).
    """
    if not layout.shapes:
        return 0
    from shapely.geometry import LineString
    from shapely.ops import split as _shapely_split
    n_subdivided = 0
    new_shapes: List[BuiltShape] = []
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            new_shapes.append(s)
            continue
        if s.polygon is None or s.polygon.is_empty:
            new_shapes.append(s)
            continue
        try:
            ring = list(s.polygon.exterior.coords)
        except Exception:
            new_shapes.append(s)
            continue
        if ring and ring[0] == ring[-1]:
            ring = ring[:-1]
        n = len(ring)
        if n < 4:
            new_shapes.append(s)
            continue
        if s.altitude is not None:
            elevs = [float(s.altitude)] * n
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == n + 1:
                na = na[:-1]
            if len(na) != n:
                new_shapes.append(s)
                continue
            elevs = [float(e) for e in na]
        else:
            new_shapes.append(s)
            continue

        def _worst_grade(pts: List[Tuple[float, float]],
                         es: List[float]) -> Tuple[float,
                                                    Optional[Tuple[int, int]]]:
            """Return (worst_grade, worst_pair_idx) over all
            vertex pairs of ``pts`` within
            SUBDIVIDE_MAX_PAIR_DIST_M and > 0.5 m apart.  The
            distance floor is critical: < 0.5 m pairs are
            essentially the same vertex with rounding noise
            and produce huge spurious grades."""
            radius2 = SUBDIVIDE_MAX_PAIR_DIST_M ** 2
            min_d2 = 0.5 ** 2
            wg = 0.0
            wp = None
            m = len(pts)
            for a in range(m):
                xa, ya = pts[a]
                ea = es[a]
                for b in range(a + 1, m):
                    xb, yb = pts[b]
                    dx_ = xa - xb
                    dy_ = ya - yb
                    d2_ = dx_ * dx_ + dy_ * dy_
                    if d2_ < min_d2 or d2_ > radius2:
                        continue
                    d_ = math.sqrt(d2_)
                    de_ = abs(ea - es[b])
                    if de_ <= TAXI_MAX_GRADE * d_ + 0.10:
                        continue
                    g_ = de_ / d_
                    if g_ > wg:
                        wg = g_
                        wp = (a, b)
            return wg, wp

        worst_grade, worst_pair = _worst_grade(ring, elevs)
        if (worst_pair is None
                or worst_grade < SUBDIVIDE_VIOLATION_GRADE):
            new_shapes.append(s)
            continue
        # Build the perpendicular cut line through the midpoint.
        i, j = worst_pair
        pi = ring[i]
        pj = ring[j]
        mx = 0.5 * (pi[0] + pj[0])
        my = 0.5 * (pi[1] + pj[1])
        dx = pj[0] - pi[0]
        dy = pj[1] - pi[1]
        seg_len = math.hypot(dx, dy)
        if seg_len < 1e-6:
            new_shapes.append(s)
            continue
        # Perpendicular unit vector (rotate +90°).
        px = -dy / seg_len
        py = dx / seg_len
        bx_min, by_min, bx_max, by_max = s.polygon.bounds
        bbox_diag = math.hypot(bx_max - bx_min, by_max - by_min)
        L = max(bbox_diag * 2.0, 1000.0)
        cut = LineString([
            (mx - L * px, my - L * py),
            (mx + L * px, my + L * py),
        ])
        try:
            parts = _shapely_split(s.polygon, cut)
        except Exception:
            new_shapes.append(s)
            continue
        sub_polys: List[Polygon] = []
        try:
            for g in getattr(parts, "geoms", [parts]):
                if (g is None or g.is_empty
                        or g.geom_type != "Polygon"
                        or g.area < SUBDIVIDE_MIN_AREA_M2):
                    continue
                sub_polys.append(g)
        except Exception:
            new_shapes.append(s)
            continue
        if len(sub_polys) < 2:
            new_shapes.append(s)
            continue

        # Pre-compute snap-target ring vertices keyed by their
        # squared snap radius for O(n) lookup per sub-vertex.
        snap_r2 = SUBDIVIDE_SNAP_RADIUS_M ** 2

        def _snap_to_ring(qx: float, qy: float
                          ) -> Tuple[float, float, int]:
            """Return the closest ring vertex within snap radius
            and its index, or ``(qx, qy, -1)`` if no ring vertex
            is close enough.
            """
            best_k = -1
            best_d2 = snap_r2
            for k in range(n):
                rx, ry = ring[k]
                d2 = (rx - qx) ** 2 + (ry - qy) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    best_k = k
            if best_k >= 0:
                rx, ry = ring[best_k]
                return rx, ry, best_k
            return qx, qy, -1

        def _lookup_elev(qx: float, qy: float, hint_idx: int
                         ) -> float:
            """Look up elevation for sub-polygon vertex.  When
            ``hint_idx`` ≥ 0 (snapped to ring vertex) use the
            original elevation directly.  Otherwise interpolate
            along the closest original ring edge.
            """
            if hint_idx >= 0:
                return elevs[hint_idx]
            best_d2 = float("inf")
            best_e = elevs[0]
            for k in range(n):
                ax, ay = ring[k]
                bx, by = ring[(k + 1) % n]
                edx = bx - ax
                edy = by - ay
                seg2 = edx * edx + edy * edy
                if seg2 < 1e-9:
                    continue
                t = ((qx - ax) * edx + (qy - ay) * edy) / seg2
                if t < 0.0:
                    t = 0.0
                elif t > 1.0:
                    t = 1.0
                ccx = ax + t * edx
                ccy = ay + t * edy
                d2 = (qx - ccx) ** 2 + (qy - ccy) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    ea = elevs[k]
                    eb = elevs[(k + 1) % n]
                    best_e = ea + t * (eb - ea)
            return round(float(best_e), 1)

        # Build snapped sub-polygons + per-vertex elevations.
        # Acceptance criterion: each sub-polygon's own worst pair
        # must be MEASURABLY better (≥ 0.5 % grade improvement)
        # than the parent's worst.  This prevents cuts that
        # technically separate the worst pair but introduce new
        # cut-line-vertex pairs of similar magnitude (the bug
        # that made the relaxed version worse than the strict).
        validated_subs: List[Tuple[Polygon, List[Tuple[float, float]],
                                    List[float]]] = []
        cut_was_useful = True
        for sp in sub_polys:
            sub_ring_raw = list(sp.exterior.coords)
            if sub_ring_raw and sub_ring_raw[0] == sub_ring_raw[-1]:
                sub_ring_raw = sub_ring_raw[:-1]
            snapped_pts: List[Tuple[float, float]] = []
            snapped_hints: List[int] = []
            for (qx, qy) in sub_ring_raw:
                sx, sy, hint = _snap_to_ring(qx, qy)
                if (snapped_pts and abs(snapped_pts[-1][0] - sx) < 1e-9
                        and abs(snapped_pts[-1][1] - sy) < 1e-9):
                    continue  # consecutive dup after snap
                snapped_pts.append((sx, sy))
                snapped_hints.append(hint)
            while (len(snapped_pts) >= 2
                   and abs(snapped_pts[0][0]
                           - snapped_pts[-1][0]) < 1e-9
                   and abs(snapped_pts[0][1]
                           - snapped_pts[-1][1]) < 1e-9):
                snapped_pts.pop()
                snapped_hints.pop()
            if len(snapped_pts) < 3:
                cut_was_useful = False
                break
            try:
                snapped_poly = Polygon(snapped_pts)
                if not snapped_poly.is_valid:
                    snapped_poly = snapped_poly.buffer(0)
                if (snapped_poly.is_empty
                        or snapped_poly.geom_type != "Polygon"
                        or snapped_poly.area < SUBDIVIDE_MIN_AREA_M2):
                    cut_was_useful = False
                    break
            except Exception:
                cut_was_useful = False
                break
            sub_elevs = [_lookup_elev(qx, qy, h)
                         for (qx, qy), h
                         in zip(snapped_pts, snapped_hints)]
            sub_worst, _ = _worst_grade(snapped_pts, sub_elevs)
            if sub_worst >= worst_grade - 0.005:
                cut_was_useful = False
                break
            validated_subs.append(
                (snapped_poly, snapped_pts, sub_elevs))

        if not cut_was_useful or len(validated_subs) < 2:
            new_shapes.append(s)
            continue

        for sp, sub_pts, sub_elevs in validated_subs:
            sub_shape = BuiltShape(
                polygon=sp, role=ROLE_JUNCTION, ref=s.ref)
            elev_range = max(sub_elevs) - min(sub_elevs)
            if elev_range < 0.05:
                sub_shape.altitude = round(
                    sum(sub_elevs) / len(sub_elevs), 1)
            else:
                closed = list(sub_elevs) + [sub_elevs[0]]
                sub_shape.node_altitudes = closed
            new_shapes.append(sub_shape)
        n_subdivided += 1
    layout.shapes = new_shapes
    return n_subdivided


def _merge_sliver_junctions_into_neighbours(
        layout: "PavementLayout",
        icao: str = "",
        sliver_area_m2: float = 1000.0,
        sliver_ratio: float = 0.05,
        shared_vertex_tol_m: float = 0.5,
        ) -> int:
    """Merge small junction polygons into adjacent larger ones.

    A "sliver" is a junction polygon whose area is below
    ``sliver_area_m2`` AND whose ratio to a neighbour's area is
    below ``sliver_ratio``.  Two junctions are "adjacent" if they
    share at least 2 boundary vertices within
    ``shared_vertex_tol_m``.

    Common cause: ``_decompose_polygon_with_holes`` cuts a
    polygon-with-holes into simple pieces, occasionally carving
    off a tiny strip when the cut grazes the polygon's edge.
    The strip and main piece share a boundary segment (the cut
    line); the strip should be merged back.  Subdivision passes
    in ``_compute_elevations`` can produce similar slivers.

    Returns the number of slivers merged.
    """
    junction_idxs = [i for i, s in enumerate(layout.shapes)
                     if s.role == ROLE_JUNCTION
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if len(junction_idxs) < 2:
        return 0
    # Cache per-shape vertex sets in meter coords for fast tests.
    j_verts: Dict[int, List[Tuple[float, float]]] = {}
    for i in junction_idxs:
        try:
            coords = list(layout.shapes[i].polygon.exterior.coords)
        except Exception:
            j_verts[i] = []
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        j_verts[i] = coords
    tol2 = shared_vertex_tol_m * shared_vertex_tol_m
    merge_into: Dict[int, int] = {}
    for i in junction_idxs:
        ai = layout.shapes[i].polygon.area
        if ai >= sliver_area_m2:
            continue
        best_idx: Optional[int] = None
        best_area = 0.0
        for j in junction_idxs:
            if j == i:
                continue
            aj = layout.shapes[j].polygon.area
            if aj <= ai:
                continue
            if (ai / aj) > sliver_ratio:
                continue
            shared = 0
            for vx, vy in j_verts[i]:
                for ux, uy in j_verts[j]:
                    if (vx - ux) ** 2 + (vy - uy) ** 2 <= tol2:
                        shared += 1
                        break
                if shared >= 2:
                    break
            if shared >= 2 and aj > best_area:
                best_idx = j
                best_area = aj
        if best_idx is not None:
            merge_into[i] = best_idx
    if not merge_into:
        return 0
    for sliver_i, target_j in merge_into.items():
        try:
            target_shape = layout.shapes[target_j]
            sliver_shape = layout.shapes[sliver_i]
            merged = unary_union([
                target_shape.polygon,
                sliver_shape.polygon])
            if (merged.geom_type == "Polygon"
                    and not merged.is_empty):
                # Build a per-vertex elevation lookup from the
                # ORIGINAL target + sliver vertices, then re-derive
                # node_altitudes for the merged polygon by nearest-
                # neighbour match.  Per user 2026-04-29 (CYXY apron
                # regression): just setting ``node_altitudes = None``
                # leaves the merged polygon with NO elevation at all
                # (no downstream re-derivation step exists), which
                # makes X-Plane interpolate from neighbour shapes
                # and produces the "terrain all over the place"
                # apron the user saw.
                lookup: List[Tuple[float, float, float]] = []
                for src_shape in (target_shape, sliver_shape):
                    if not src_shape.node_altitudes:
                        # Sloped or flat alternatives.
                        if (src_shape.altitude_high is not None
                                and src_shape.altitude_low is not None):
                            avg = 0.5 * (
                                src_shape.altitude_high
                                + src_shape.altitude_low)
                        elif src_shape.altitude is not None:
                            avg = float(src_shape.altitude)
                        else:
                            continue
                        try:
                            sc = list(
                                src_shape.polygon.exterior.coords)
                        except Exception:
                            continue
                        if sc and sc[0] == sc[-1]:
                            sc = sc[:-1]
                        for sx, sy in sc:
                            lookup.append((sx, sy, float(avg)))
                        continue
                    src_alts = list(src_shape.node_altitudes)
                    try:
                        sc = list(src_shape.polygon.exterior.coords)
                    except Exception:
                        continue
                    if sc and sc[0] == sc[-1]:
                        sc = sc[:-1]
                    if (len(src_alts) == len(sc) + 1
                            and src_alts[0] == src_alts[-1]):
                        src_alts = src_alts[:-1]
                    for k, (sx, sy) in enumerate(sc):
                        if k >= len(src_alts):
                            break
                        lookup.append(
                            (sx, sy, float(src_alts[k])))
                if not lookup:
                    layout.shapes[target_j].polygon = merged
                    layout.shapes[target_j].node_altitudes = None
                    continue
                merged_coords = list(merged.exterior.coords)
                if (merged_coords
                        and merged_coords[0] == merged_coords[-1]):
                    merged_coords_open = merged_coords[:-1]
                else:
                    merged_coords_open = merged_coords
                new_alts: List[float] = []
                for mx, my in merged_coords_open:
                    best_d2 = float("inf")
                    best_alt = 0.0
                    for sx, sy, sa in lookup:
                        d2 = (mx - sx) ** 2 + (my - sy) ** 2
                        if d2 < best_d2:
                            best_d2 = d2
                            best_alt = sa
                    new_alts.append(round(best_alt, 1))
                # Closing-vertex repeat for the ring.
                if new_alts:
                    new_alts.append(new_alts[0])
                target_shape.polygon = merged
                target_shape.node_altitudes = new_alts
        except Exception:
            continue
    sliver_set = set(merge_into.keys())
    layout.shapes = [
        s for k, s in enumerate(layout.shapes)
        if k not in sliver_set]
    try:
        UI.vprint(1,
            f"  [pav-builder] {icao}: merged "
            f"{len(merge_into)} sliver junction(s) into "
            f"adjacent larger junctions.")
    except Exception:
        pass
    return len(merge_into)


def _report_within_shape_violations(
        layout: "PavementLayout", icao: str) -> None:
    """Layer 3: scan every emitted polygon for vertex-pair grade
    violations and emit a stderr WARN summary.

    Pair set: only pairs within ``WITHIN_SHAPE_VIOLATION_RADIUS_M``
    of each other (the same Triangle4XP-plausible-edge radius
    check_grade.py uses).  Far-pair grades are noisy false
    positives — Triangle4XP would interpose a Steiner point and
    never connect them directly.
    """
    if not layout.shapes:
        return
    radius = WITHIN_SHAPE_VIOLATION_RADIUS_M
    radius2 = radius * radius
    cos0 = math.cos(math.radians(layout.anchor[0]))
    n_viol = 0
    worst_pct = 0.0
    worst_info: Optional[Tuple[str, str, float, float, float, float]] = None
    for s in layout.shapes:
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        n = len(coords)
        if n < 3:
            continue
        if s.altitude is not None:
            elevs = [float(s.altitude)] * n
        elif s.node_altitudes is not None:
            na = list(s.node_altitudes)
            if len(na) == n + 1:
                na = na[:-1]
            if len(na) != n:
                continue
            elevs = [float(e) for e in na]
        elif (s.altitude_high is not None and s.altitude_low is not None
              and n == 4):
            elevs = [float(s.altitude_high), float(s.altitude_low),
                     float(s.altitude_low), float(s.altitude_high)]
        else:
            continue
        # Convert ring to local meter coords for grade check.
        coords_m = []
        for (lon, lat) in coords:
            # Some shapes' polygons store (x, y) in meters already
            # (the new pavement builder works in meters); others
            # might store (lon, lat).  Detect by magnitude.
            if abs(lon) > 180.0 or abs(lat) > 180.0:
                coords_m.append((lon, lat))  # already meters
            else:
                # lat/lon → meters
                x = math.radians(lon - layout.anchor[1]) * R_EARTH * cos0
                y = math.radians(lat - layout.anchor[0]) * R_EARTH
                coords_m.append((x, y))
        for i in range(n):
            xi, yi = coords_m[i]
            ei = elevs[i]
            for j in range(i + 1, n):
                xj, yj = coords_m[j]
                dx = xi - xj
                dy = yi - yj
                d2 = dx * dx + dy * dy
                if d2 > radius2 or d2 < 0.25:
                    continue
                d = math.sqrt(d2)
                de = abs(ei - elevs[j])
                if de <= TAXI_MAX_GRADE * d + 0.10:
                    continue
                pct = (de / d) * 100.0
                if pct > worst_pct:
                    worst_pct = pct
                    worst_info = (
                        s.role or "?", s.ref or "",
                        ei, elevs[j], d, de)
                n_viol += 1
    if n_viol > 0:
        try:
            import sys as _sys
            msg = (f"  [pav-builder] WARN: {icao}: {n_viol} within-shape "
                   f"grade violations (> {TAXI_MAX_GRADE * 100:.1f}%, "
                   f"checked pairs within {radius:.0f} m)")
            if worst_info is not None:
                role, ref, ea, eb, d, de = worst_info
                rstr = f"/{ref}" if ref else ""
                msg += (f"; worst {worst_pct:.1f}% on "
                        f"{role}{rstr} ({ea:.1f} → {eb:.1f}, "
                        f"d={d:.1f}m, de={de:.1f}m)")
            msg += "."
            UI.vprint(1, msg)
        except Exception:
            pass


WITHIN_SHAPE_VIOLATION_RADIUS_M = 60.0   # max pair distance for
                                          # within-shape violation
                                          # reporting; matches
                                          # check_grade.py's
                                          # WITHIN_SHAPE_MAX_PAIR_DIST_M.


SHARED_VERTEX_CLUSTER_TOL_M = 1.5


def _drop_overlap_against_fixed_shapes(
        layout: "PavementLayout",
        icao: str = "") -> None:
    """Enforce the no-overlap invariant on the layout.

    Walks every shape and, where it overlaps another shape, modifies
    or drops it so no two shapes overlap.  Order of priority (the
    LATER a role appears in this list, the more it "yields"):

    1. RUNWAY corners (CIFP-anchored, immutable footprint).
    2. TAXI rect (one per OSM centerline; rect dedup ran upstream).
    3. TERMINAL pad (OSM building).
    4. JUNCTION (residue — must clip to fit around all of the above).

    Strategy:

    * Drop duplicate TERMINAL polygons (a terminal entirely inside
      another terminal is a duplicate from OSM relation parsing).
    * Drop duplicate or heavily-overlapping RUNWAY segments
      (apron-merged runway segmentation can produce overlap with
      the original single-rect runway).
    * Clip each JUNCTION against every fixed shape (rect / runway /
      terminal) and against larger junctions.  Iterate up to 4
      passes so chained clips converge.

    Mutates ``layout.shapes`` in place.
    """
    from shapely.strtree import STRtree
    MIN_KEEP_AREA_M2 = 0.5
    NOISE_OVERLAP_M2 = 1.0   # ignore sub-1 m² overlaps as float noise

    def _valid_poly(p: Optional[Polygon]) -> Optional[Polygon]:
        if p is None or p.is_empty:
            return None
        if p.geom_type != "Polygon":
            return None
        if not p.is_valid:
            try:
                p = p.buffer(0)
            except Exception:
                return None
            if p.is_empty or p.geom_type != "Polygon":
                return None
        return p

    def _clip_keep_largest(p: Polygon, c: Polygon
                           ) -> Optional[Polygon]:
        """Return ``p.difference(c)``, picking the largest piece if
        the difference is a MultiPolygon.  Returns None if the
        result is empty / below MIN_KEEP_AREA_M2."""
        try:
            d = p.difference(c)
        except Exception:
            return p
        if d.is_empty:
            return None
        if d.geom_type == "Polygon":
            return d if d.area >= MIN_KEEP_AREA_M2 else None
        if d.geom_type == "MultiPolygon":
            pieces = [g for g in d.geoms
                      if g.geom_type == "Polygon"
                      and g.area >= MIN_KEEP_AREA_M2]
            if not pieces:
                return None
            pieces.sort(key=lambda g: -g.area)
            return pieces[0]
        return None

    n_dropped = 0
    n_clipped = 0
    DUPLICATE_FRAC = 0.80

    # ── Step 1: enforce same-role no-overlap.  Two shapes of the
    # same role (two terminals, two runway segments, two taxi
    # rects) must never overlap.  Three behaviours:
    #
    #   * If one shape is mostly inside the other (≥ DUPLICATE_FRAC
    #     of its area), drop it as a duplicate.
    #   * Otherwise, clip the smaller shape against the larger so
    #     the overlap region is removed from the smaller (the
    #     larger is "the more authoritative" footprint).
    #   * If the clip leaves no usable polygon, drop it.
    for role_set in [
            {ROLE_TERMINAL},
            {ROLE_RUNWAY},
            {ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
             ROLE_STUB, ROLE_CROSS_CONNECTOR}]:
        # Iterate to a fixed point in case clipping creates new
        # adjacencies that need further clipping.
        for _ in range(4):
            candidates: List[int] = [
                i for i, s in enumerate(layout.shapes)
                if s.role in role_set
                and _valid_poly(s.polygon) is not None]
            # Sort by area DESC so smaller shapes are clipped
            # against larger ones (we walk pairs (i, j) with i<j
            # and clip the SMALLER of the pair).
            candidates.sort(
                key=lambda i: -layout.shapes[i].polygon.area)
            any_change = False
            for ai in range(len(candidates)):
                i = candidates[ai]
                if layout.shapes[i].polygon is None:
                    continue
                pi = layout.shapes[i].polygon
                for bi in range(ai + 1, len(candidates)):
                    j = candidates[bi]
                    if layout.shapes[j].polygon is None:
                        continue
                    pj = layout.shapes[j].polygon
                    try:
                        if not pi.intersects(pj):
                            continue
                        inter = pi.intersection(pj)
                        if (inter.is_empty
                                or inter.area <= NOISE_OVERLAP_M2):
                            continue
                        # Duplicate test.
                        a_min = min(pi.area, pj.area)
                        if (a_min > 0
                                and inter.area / a_min
                                >= DUPLICATE_FRAC):
                            # j is the smaller (sorted desc) —
                            # drop it.
                            layout.shapes[j].polygon = None
                            n_dropped += 1
                            any_change = True
                            continue
                        # Partial overlap — clip j against i.
                        clipped = _clip_keep_largest(pj, pi)
                        if clipped is None:
                            layout.shapes[j].polygon = None
                            n_dropped += 1
                        else:
                            layout.shapes[j].polygon = clipped
                            n_clipped += 1
                        any_change = True
                    except Exception:
                        continue
            if not any_change:
                break
    layout.shapes = [s for s in layout.shapes
                     if s.polygon is not None]

    # ── Step 2: priority-ordered clip pass.  Each role yields to
    # everything LISTED ABOVE it in the ``priority`` list:
    #   * RUNWAY (CIFP-anchored) — never modified.
    #   * TERMINAL — yields to runway only.
    #   * TAXI rects — yield to runway + terminal.
    #   * JUNCTION (residue) — yields to everything.
    # Each shape is clipped against every higher-priority shape it
    # overlaps; the result keeps only the largest piece if the clip
    # produces multiple disjoint fragments.  Same-priority shapes
    # of the JUNCTION class additionally yield to LARGER junctions
    # so two junctions can't both claim the same residue area.
    priority: List[set] = [
        {ROLE_RUNWAY},
        {ROLE_TERMINAL},
        {ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
         ROLE_STUB, ROLE_CROSS_CONNECTOR},
        {ROLE_JUNCTION},
    ]
    for outer in range(4):
        any_change = False
        # For each tier (after the first), clip its shapes against
        # all higher-priority shapes.
        for tier_idx in range(1, len(priority)):
            tier_roles = priority[tier_idx]
            higher_polys: List[Polygon] = []
            for s in layout.shapes:
                if s.polygon is None:
                    continue
                role_tier = next(
                    (ti for ti, rs in enumerate(priority)
                     if s.role in rs), -1)
                if 0 <= role_tier < tier_idx:
                    p = _valid_poly(s.polygon)
                    if p is not None:
                        higher_polys.append(p)
            higher_tree = (STRtree(higher_polys)
                           if higher_polys else None)
            # Targets in this tier, sorted by area (largest first
            # — within the JUNCTION tier this lets smaller junctions
            # later be clipped against the already-finalised
            # larger ones).
            target_idx: List[int] = [
                i for i, s in enumerate(layout.shapes)
                if s.role in tier_roles
                and _valid_poly(s.polygon) is not None]
            target_idx.sort(
                key=lambda i: -layout.shapes[i].polygon.area)
            for k, i in enumerate(target_idx):
                tp = layout.shapes[i].polygon
                if tp is None:
                    continue
                new_p: Optional[Polygon] = tp
                # Clip against higher-priority shapes.
                if higher_tree is not None:
                    for hit in higher_tree.query(new_p):
                        fp = higher_polys[hit]
                        try:
                            if not new_p.intersects(fp):
                                continue
                            inter = new_p.intersection(fp)
                            if (inter.is_empty
                                    or inter.area
                                    <= NOISE_OVERLAP_M2):
                                continue
                            clipped = _clip_keep_largest(new_p, fp)
                            if clipped is None:
                                new_p = None
                                break
                            new_p = clipped
                            any_change = True
                            n_clipped += 1
                        except Exception:
                            continue
                if (new_p is not None
                        and tier_roles == {ROLE_JUNCTION}):
                    # Also clip against LARGER same-tier junctions.
                    for k2 in range(k):
                        i2 = target_idx[k2]
                        tp2 = layout.shapes[i2].polygon
                        if tp2 is None:
                            continue
                        try:
                            if not new_p.intersects(tp2):
                                continue
                            inter = new_p.intersection(tp2)
                            if (inter.is_empty
                                    or inter.area
                                    <= NOISE_OVERLAP_M2):
                                continue
                            clipped = _clip_keep_largest(new_p, tp2)
                            if clipped is None:
                                new_p = None
                                break
                            new_p = clipped
                            any_change = True
                            n_clipped += 1
                        except Exception:
                            continue
                if new_p is None:
                    layout.shapes[i].polygon = None
                    n_dropped += 1
                    continue
                if new_p is not tp:
                    layout.shapes[i].polygon = new_p
        if not any_change:
            break

    layout.shapes = [s for s in layout.shapes
                     if s.polygon is not None]

    if (n_clipped + n_dropped) > 0:
        try:
            UI.vprint(1,
                f"  [pav-builder] {icao}: overlap-clip pass — "
                f"{n_clipped} clip operation(s), "
                f"{n_dropped} shape(s) dropped.")
        except Exception:
            pass
