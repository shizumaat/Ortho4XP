"""End-to-end orchestration of the airport pavement builder.

Sequences the phase-1 (geometry + role) and phase-2 (elevation)
passes by calling out into the focused ``O4_Pavement_*`` modules:

* OSM tile + airport extraction.
* apt.dat selection (best-of-OSM vs custom-pack heuristics).
* Phase-1 construction → ``O4_Pavement_Rects``,
  ``O4_Pavement_Centerlines``, ``O4_Pavement_Stubs``,
  ``O4_Pavement_Junctions``, ``O4_Pavement_Terminals``.
* Phase-2 elevation → ``O4_Pavement_Elevation``.
* Feature emit → ``O4_Pavement_Boundary``,
  ``O4_Pavement_Groundside``, ``O4_Pavement_Bridges``.
* Output via ``PavementLayout.to_osm`` (in ``O4_Pavement_Layout``).

Public API:

    build_airport_pavement(icao, xplane_root, *, compute_elevations=True,
                            taxiway_data=None, tile_dem=None,
                            airport_boundary=None)

Backward-compat shim ``O4_Airport_Pavement_Builder`` re-exports
``build_airport_pavement`` (and a few helpers used by other
modules) so existing call sites keep working.
"""
from __future__ import annotations

import math
import os
import re
from typing import Dict, List, Optional, Sequence, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.ops import (
    linemerge, nearest_points, transform as shp_transform, unary_union)

# Narrow exception tuple for shapely / numeric-geometry failure
# modes + file I/O.  Programming errors propagate so they surface
# immediately rather than being silently masked at runtime.
_GEOM_EXC = (OSError, ValueError, TypeError, KeyError,
             IndexError, RuntimeError,
             GEOSException, TopologicalError)

import O4_File_Names as FNAMES
import O4_UI_Utils as UI

from . import apt_dat_reader as APR
from .pavement import strips as PS


# ──────────────────────────────────────────────────────────────────
# Constants (re-exported from O4_Pavement_Config + O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────
from .config import (
    MIN_SEGMENT_LEN_M,
    LOAD_DSF_PAVEMENT,
    RUNWAY_APRON_AREA_RATIO,
)
from .layout import (
    BuiltShape,
    PavementLayout,
    R_EARTH,
    ROLE_RUNWAY,
    ROLE_STUB,
    ROLE_TERMINAL,
    _airport_anchor,
    _projection,
)


# ──────────────────────────────────────────────────────────────────
# Data model (re-exported from O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Input loaders
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# OSM cache + apt.dat selection helpers
# ──────────────────────────────────────────────────────────────────
from .osm_load import (
    _load_osm_airports,
    _load_osm_big_roads,
    _pick_best_apt_dat_against_osm,
)
from . import finalize, junction_emit


# ──────────────────────────────────────────────────────────────────
# Meter-space projection (re-exported from O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Runway rects, crossings, shoulders (re-exported from
# O4_Pavement_Runways)
# ──────────────────────────────────────────────────────────────────
from .pavement.runways import (
    _detect_runway_shoulders,
    _runway_rect_m,
)


# ──────────────────────────────────────────────────────────────────
# Pavement union helpers (re-exported from pavement.union_helpers)
# ──────────────────────────────────────────────────────────────────
from .pavement.union_helpers import (
    _merge_near_touching,
)


# ──────────────────────────────────────────────────────────────────
# Stub residue / runway-bridge helpers
# (re-exported from pavement.stubs)
# ──────────────────────────────────────────────────────────────────


from .pavement.absorption import (
    _drop_primary_parallels_embedded_in_pavement,
)


# ──────────────────────────────────────────────────────────────────
# Top-level builder
# ──────────────────────────────────────────────────────────────────

def build_airport_pavement(icao: str, xplane_root: str,
                            *,
                            compute_elevations: bool = True,
                            taxiway_data=None,
                            tile_dem=None,
                            airport_boundary=None,
                            current_tile_lat=None,
                            current_tile_lon=None,
                            ) -> PavementLayout:
    """Build the complete role-classified layout for ``icao``.

    The layout is ready to compare against a target OSM via
    ``tools/compare_target.py``.

    When ``compute_elevations`` is True (default), a Phase-2
    elevation pass runs at the end:
      * Single runway shapes are replaced with per-100m segmented
        runway rects produced by the legacy CIFP+DEM generator.
      * Taxi rects get ``altitude_high``/``altitude_low`` tags
        from DEM sampling at axis endpoints, anchored to the
        adjacent runway segment when within 30 m.
      * Terminal pads get ``altitude`` from the DEM at centroid.
      * Junctions, buildings, aprons are left un-elevated (X-Plane
        triangulator interpolates them from neighbouring shared
        vertices).

    Optional integration parameters (all default to ``None`` so
    standalone callers — ``tools/build_target_osm.py``, the test
    harness — work unchanged):

      * ``taxiway_data``: per-airport list from
        ``osm_aeroway.extract_taxiway_info``.  Threaded through to
        Pipeline's centerline-union helper so reclassification /
        downstream consumers can see OSM taxiway centerlines that
        didn't survive into rects.
      * ``tile_dem``: pre-loaded ``O4_DEM_Utils.DEM`` for the
        containing tile (typically Ortho4XP's smoothed DEM after
        ``smooth_raster_over_airports``).  When supplied, the
        Phase-2 elevation solver and boundary-shape emit consume
        this DEM directly instead of calling
        ``_load_airport_dem`` per-airport — avoids redundant DEM
        loads in the tile-pipeline driver and keeps auto_patch's
        elevation field aligned with Ortho4XP's smoothed terrain.
      * ``airport_boundary``: optional pre-computed airport
        boundary (``dico_airports[icao]['boundary']`` from
        Ortho4XP's ``update_airport_boundaries``).  Currently
        reserved — the parameter is plumbed through but the
        boundary-shape emit still derives its outline from
        apt.dat row-130 until the source-of-truth question is
        decided.
    """
    apt_path = _pick_best_apt_dat_against_osm(xplane_root, icao)
    if apt_path is None:
        raise RuntimeError(f"No apt.dat found for {icao}")
    apt = APR.load_airport(apt_path, icao)
    if apt is None:
        raise RuntimeError(f"Could not load airport block for {icao}")

    anchor = _airport_anchor(apt)
    to_m = _projection(anchor)

    layout = PavementLayout(icao=icao, anchor=anchor,
                             apt_dat_path=apt_path)

    # Project the apt.dat row-130 airport boundary (in lat/lon) into
    # meter space so downstream emission has it ready when the
    # boundary shape is built.
    if apt.boundary is not None and not apt.boundary.is_empty:
        try:
            from shapely.ops import transform as _shp_transform
            layout.airport_boundary = _shp_transform(
                lambda lon, lat, z=None: to_m(lon, lat),
                apt.boundary)
        except _GEOM_EXC:
            layout.airport_boundary = None

    # ── Runways ──────────────────────────────────────────────────
    runway_polys: List[Polygon] = []
    rwy_bearings: List[float] = []
    rwy_centerlines: List[LineString] = []
    for r in apt.runways:
        rect = _runway_rect_m(r, to_m)
        if rect.is_empty:
            continue
        runway_polys.append(rect)
        ref = f"{r.desig_a}/{r.desig_b}"
        layout.shapes.append(BuiltShape(
            polygon=rect, role=ROLE_RUNWAY, ref=ref))
        ax, ay = to_m(r.lon_a, r.lat_a)
        bx, by = to_m(r.lon_b, r.lat_b)
        if math.hypot(bx - ax, by - ay) > 1.0:
            rwy_centerlines.append(LineString([(ax, ay), (bx, by)]))
            rwy_bearings.append(
                math.degrees(math.atan2(bx - ax, by - ay)) % 180.0)
    layout.runway_union = unary_union(runway_polys) if runway_polys else None

    # ── Pavement union in meter space ────────────────────────────
    pav_polys: List[Polygon] = []
    for pav in apt.pavements:
        if pav.polygon is None or pav.polygon.is_empty:
            continue
        pm = shp_transform(to_m, pav.polygon)
        if pm.is_empty:
            continue
        if pm.geom_type == "Polygon":
            pav_polys.append(pm)
        else:
            pav_polys.extend(g for g in getattr(pm, "geoms", [])
                             if g.geom_type == "Polygon")
    # Snapshot the apt.dat-only polygon list before DSF additions.
    # Terminal-pad selection prefers the SMALLEST containing
    # polygon, and DSF often ships small overlay-style polygons
    # over apt.dat pavement; without this snapshot a small DSF
    # overlay covering part of the apron will win over the larger
    # apt.dat terminal pavement and the resulting terminal pad
    # loses most of its area (SPJC terminal1 regressed from
    # 105 K m² → 35 K m² before this fix).
    apt_only_pav_polys: List[Polygon] = list(pav_polys)

    # ── Runway shoulder absorption ─────────────────────────────────
    # Long thin row-110 polygons parallel to a runway and touching
    # or overlapping it are runway shoulders (or, when wider than
    # the apt.dat row-100 width and centered on the runway, the
    # runway's own envelope polygon — apt.dat sometimes labels
    # these as "taxiways", e.g. HECA's "New Taxiway 1").  Fold them
    # into the runway: widen the runway emit to the union of
    # perpendicular extents, mutate the runway's apt.dat record so
    # downstream CIFP segmenting picks up the new width, and remove
    # the absorbed polygons from the pavement set so they don't
    # re-emit as junction polygons wrapping the runway.
    #
    # Asymmetric shoulders (one side only — common at gravel
    # crosswind runways like CYXY's 02/20) are handled by shifting
    # the runway centerline toward the shoulder midpoint while
    # widening to the union extent.  The CIFP threshold elevations
    # (anchored at the original runway thresholds) still apply
    # because the threshold lat/lon stays paired with the same
    # apt.dat designation; the perpendicular shift moves the
    # centerline by the shoulder offset (typically < 20 m, well
    # within DEM noise tolerance).
    absorbed_pav_indices: set = set()
    for ridx, r in enumerate(apt.runways):
        new_left, new_right, absorbed = _detect_runway_shoulders(
            r, to_m, pav_polys)
        new_width = new_right - new_left
        # Only widen when the new extent meaningfully exceeds the
        # apt.dat row-100 width (≥ 0.5 m on top of current).
        if new_width <= r.width_m + 0.5:
            continue
        offset = 0.5 * (new_left + new_right)
        # Shift centerline by ``offset`` perpendicular.  The
        # perpendicular vector in meter space is (nx, ny) =
        # (-uy, ux), where (ux, uy) is the runway's
        # along-axis unit vector.
        ax_m, ay_m = to_m(r.lon_a, r.lat_a)
        bx_m, by_m = to_m(r.lon_b, r.lat_b)
        udx = bx_m - ax_m
        udy = by_m - ay_m
        L_axis = math.hypot(udx, udy)
        if L_axis < 1.0:
            continue
        ux_m, uy_m = udx / L_axis, udy / L_axis
        nx_m, ny_m = -uy_m, ux_m
        # Apply offset back through the inverse of ``to_m``.  Our
        # to_m projection (see ``_projection``) is anchored at
        # ``layout.anchor`` = (lat0, lon0) and uses cos(lat0).
        lat0, lon0 = layout.anchor
        cos0 = math.cos(math.radians(lat0))
        d_lat_per_m = 1.0 / R_EARTH
        d_lon_per_m = 1.0 / (R_EARTH * cos0) if cos0 > 1e-9 else 0.0
        d_lat = math.degrees(ny_m * offset * d_lat_per_m)
        d_lon = math.degrees(nx_m * offset * d_lon_per_m)
        old_w = r.width_m
        old_lat_a, old_lon_a = r.lat_a, r.lon_a
        if abs(offset) > 0.05:
            r.lat_a = r.lat_a + d_lat
            r.lon_a = r.lon_a + d_lon
            r.lat_b = r.lat_b + d_lat
            r.lon_b = r.lon_b + d_lon
        r.width_m = new_width
        new_rect = _runway_rect_m(r, to_m)
        if new_rect.is_empty:
            r.lat_a, r.lon_a = old_lat_a, old_lon_a
            r.width_m = old_w
            continue
        runway_polys[ridx] = new_rect
        ref = f"{r.desig_a}/{r.desig_b}"
        for s in layout.shapes:
            if s.role == ROLE_RUNWAY and s.ref == ref:
                s.polygon = new_rect
                break
        absorbed_pav_indices.update(absorbed)
        try:
            UI.vprint(1,
                f"  [pav-builder] {icao}: widened runway "
                f"{r.desig_a}/{r.desig_b}: {old_w:.1f}m → "
                f"{r.width_m:.1f}m"
                + (f" (centerline shifted {offset:+.1f}m)"
                if abs(offset) > 0.5 else "")
                + f" — absorbed {len(absorbed)} shoulder polygon(s).")
        except _GEOM_EXC:
            pass

    if absorbed_pav_indices:
        # Filter both pav_polys and the apt-only snapshot.  Use
        # WKB-identity to filter apt_only_pav_polys (its indices
        # don't necessarily match pav_polys' if either was modified;
        # the snapshot was taken just above so they're identical at
        # this point, but using WKB is robust to future changes).
        absorbed_wkbs = {pav_polys[i].wkb for i in absorbed_pav_indices
                         if 0 <= i < len(pav_polys)}
        pav_polys = [p for i, p in enumerate(pav_polys)
                     if i not in absorbed_pav_indices]
        apt_only_pav_polys = [p for p in apt_only_pav_polys
                              if p.wkb not in absorbed_wkbs]
        layout.runway_union = (unary_union(runway_polys)
                                if runway_polys else None)

    # ── Pavement-runway intersection points (per user 2026-05-05) ──
    # Walk each apt.dat row-110 pavement polygon's exterior; collect
    # the vertices that sit within ``INTERSECTION_PROX_M`` of a
    # runway's 4-corner rect boundary AND project to a centerline
    # parameter strictly between 0 and 1 (not at the runway ends).
    # These t-values become segment seam corners during Phase 2's
    # runway segmentation, so chain corners (= the runway-union
    # exterior) align exactly with apt.dat-pavement boundaries.
    # Without this, the junction-widening pass has to bridge the
    # gap with boundary-trace waypoints — the alignment makes that
    # unnecessary.  Per user direction: dedup intersections within
    # 2 m centerline distance (a junction can span a 2 m gap
    # without needing a node there).
    #
    # Per user 2026-05-11: tolerance widened from 0.5 m to 3.0 m so
    # apt.dat pavement boundaries drawn slightly INSIDE the row-100
    # runway rect (1-2 m offsets are common — SPJC's row-110 stops
    # 1.75 m short of 16R/34L at the V1 throat) still register as
    # runway intersection points.  Without this, the segmenter
    # doesn't insert a seam corner at the chart-level pavement
    # transition, and the downstream junction-widening pass can't
    # share a vertex with the runway there — leaving a visible
    # 1-2 m sliver gap between every taxi-junction and the runway
    # boundary.  The runway segmenter projects each kept point onto
    # the centerline and emits the seam corner there; the actual
    # corner sits on the runway boundary (rect-perpendicular at
    # half-width) regardless of how far the source pavement vertex
    # was off-edge, so widening the tolerance just unlocks more
    # near-runway pavement landmarks as segmentation breakpoints
    # without distorting the segmenter's output geometry.
    INTERSECTION_PROX_M = 3.0
    # Dedup proportionally to PROX so multi-vertex clusters of a
    # single pavement transition (row-110 boundaries drawn with 3-4
    # vertices within a 3 m span at the runway edge) collapse to one
    # seam corner instead of fragmenting the runway segment.
    INTERSECTION_DEDUP_M = 5.0
    pav_runway_intersections: dict = {}
    for ridx, r in enumerate(apt.runways):
        if ridx >= len(runway_polys):
            continue
        rect = runway_polys[ridx]
        if rect is None or rect.is_empty:
            continue
        rect_boundary = rect.exterior
        cl_ax, cl_ay = to_m(r.lon_a, r.lat_a)
        cl_bx, cl_by = to_m(r.lon_b, r.lat_b)
        cl_dx = cl_bx - cl_ax
        cl_dy = cl_by - cl_ay
        cl_L2 = cl_dx * cl_dx + cl_dy * cl_dy
        if cl_L2 < 1.0:
            continue
        # Extend centerline endpoints by blast pads so ``t`` is
        # computed in the blast-extended frame — matching both the
        # rect_boundary (which includes blast pads, see
        # ``_runway_rect_m``) and the segmenter's own phys_end_a/b
        # parameterisation (which also absorbs blast pads).  Without
        # this, a row-110 vertex sitting in the blast-pad zone (e.g.
        # SPJC 16R V1 throat at 1.79 m off the runway boundary,
        # ax≈−5 m in row-100 frame) lands at t<0 and gets dropped
        # by ``end_skirt`` — the segmenter never sees it as a
        # candidate breakpoint, so the sloped end-segment (which
        # must stay 4-corner) can't split there to give the apron
        # junction a shared snap node.
        blast_a_m = r.blast_a_m or 0.0
        blast_b_m = r.blast_b_m or 0.0
        if blast_a_m > 0.0 or blast_b_m > 0.0:
            row100_dist = math.sqrt(cl_L2)
            ux = cl_dx / row100_dist
            uy = cl_dy / row100_dist
            cl_ax -= ux * blast_a_m
            cl_ay -= uy * blast_a_m
            cl_bx += ux * blast_b_m
            cl_by += uy * blast_b_m
            cl_dx = cl_bx - cl_ax
            cl_dy = cl_by - cl_ay
            cl_L2 = cl_dx * cl_dx + cl_dy * cl_dy
        phys_dist = math.sqrt(cl_L2)
        # Avoid the runway end zones — the segmenter handles those
        # via thresholds + physical-end anchors and we don't want
        # spurious end-zone seams.
        end_skirt_t = 5.0 / phys_dist
        intersections: List[Tuple[float, float]] = []
        for pav_poly in apt_only_pav_polys:
            try:
                ring = pav_poly.exterior
            except _GEOM_EXC:
                continue
            coords = list(ring.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            for px, py in coords:
                if rect_boundary.distance(Point(px, py)) > INTERSECTION_PROX_M:
                    continue
                t = ((px - cl_ax) * cl_dx
                     + (py - cl_ay) * cl_dy) / cl_L2
                if t <= end_skirt_t or t >= 1.0 - end_skirt_t:
                    continue
                intersections.append((t, px, py))
        if not intersections:
            continue
        # Sort by centerline t and dedup.  Per user 2026-05-11: dedup
        # by EUCLIDEAN distance between consecutive candidate points
        # rather than centerline-t alone.  A row-110 boundary that
        # approaches the runway with a slight angle puts multiple
        # close-together apt.dat vertices on the runway edge —
        # each at a slightly different axial position but only ~3 m
        # apart in real space.  Centerline-t dedup keeps them all
        # (their t values differ by ≥ dedup_t_gap); euclidean dedup
        # merges them into a single seam corner, which is what the
        # runway segmenter actually needs.  Without this, every
        # close-together row-110 vertex becomes a runway-segment
        # seam corner and the downstream junction polygon has to
        # wrap around all of them (the cluster of -349/-351/-352
        # corners on -10109's east edge that pinched -10182).
        intersections.sort(key=lambda x: x[0])
        dedup_m2 = INTERSECTION_DEDUP_M * INTERSECTION_DEDUP_M
        deduped: List[Tuple[float, float, float]] = []
        for t, px, py in intersections:
            if deduped:
                dpx = px - deduped[-1][1]
                dpy = py - deduped[-1][2]
                if dpx * dpx + dpy * dpy < dedup_m2:
                    continue
            deduped.append((t, px, py))
        # Convert intersection meter-coords back to lat/lon via the
        # layout's m_to_ll (the segmenter consumes lat/lon).  Store
        # under both designator orderings so the segmenter lookup
        # finds them regardless of which key it tries.
        ll_pts = [layout.m_to_ll(px, py) for _, px, py in deduped]
        for key in (
                (r.desig_a, r.desig_b),
                (r.desig_b, r.desig_a),
                ("RW" + r.desig_a.lstrip("RW"),
                 "RW" + r.desig_b.lstrip("RW")),
                ("RW" + r.desig_b.lstrip("RW"),
                 "RW" + r.desig_a.lstrip("RW"))):
            pav_runway_intersections[key] = list(ll_pts)
    layout._pav_runway_intersections = pav_runway_intersections

    # Add draped pavement polygons from every available DSF for
    # this airport.  Some scenery packs (e.g. CYXY Whitehorse) ship
    # pavement geometry as DSF draped polygons referencing
    # ``lib/airport/pavement/*.pol`` definitions, with little or no
    # apt.dat row-110 coverage; for those airports the DSF is the
    # primary pavement source and we need to admit it.  Other packs
    # (e.g. SPJC Custom Scenery) ship apt.dat row-110 pavement AND
    # add layered visual overlays on top via DSF — emitting those
    # overlays as pavement duplicates the apt.dat coverage and pulls
    # non-pavement decoration into the layout.
    #
    # Three-tier filtering:
    #   1. ``O4_DSF_Reader._is_pavement_def`` admits only X-Plane
    #      stock pavement library paths (``lib/airport/pavement/...``
    #      and ``lib/airport/ground/pavement/...``).  Third-party
    #      libraries are dropped at this stage.
    #   2. Distance gate: the DSF tile is 1° × 1° (~110 km a side),
    #      and a single tile covers many airports' pavement.  Any
    #      DSF polygon whose bbox lies more than
    #      ``DSF_AIRPORT_RADIUS_M`` (5 km) from THIS airport's
    #      runway-bbox is somebody else's pavement — drop it.
    #      Caught the SPJC regression where 9 junctions ended up
    #      ~20 km away at SPLP because the SPJC custom scenery's
    #      DSF tile contains both airports' pavement.
    #   3. Overlay check: each surviving DSF polygon is compared
    #      against the apt.dat pavement union built so far.  If the
    #      polygon mostly overlaps existing pavement (≥ 80 %
    #      inside), treat it as an overlay and drop it entirely —
    #      preserves the apt.dat geometry.  Only DSF polygons that
    #      contribute substantially NEW coverage are appended.
    DSF_OVERLAY_FRAC = 0.80
    DSF_AIRPORT_RADIUS_M = 5_000.0
    # Per user 2026-04-29: drop DSF pavement polygons whose area
    # exceeds the largest apt.dat pavement polygon by more than
    # DSF_MAX_AREA_VS_APT_DAT_RATIO×.  Real airport pavement
    # polygons (taxiways, aprons, runway aprons) are bounded in
    # scale by the largest features apt.dat already represents at
    # the airport.  A DSF polygon dramatically larger than apt.dat's
    # biggest is a coarse "ground tile" — pavement-textured
    # decorative geometry painted across the whole airport surface
    # rather than a real pavement feature.  Confirmed at HECA
    # (Tai Models scenery), where ``lib/airport/ground/pavement/
    # asphalt/patched.pol`` instances of 3.5 M m² (with only 22
    # vertices, perimeter ~8 km) and 1.25 M m² (36 verts) overlay
    # the entire airport, dwarfing apt.dat's largest pavement
    # polygon at 378 k m² and inflating the rect-detection's
    # half-width probes to 100 m+ across what should be a 30 m
    # taxi corridor.  Safe at SPJC / CYXY / SPLP: their largest
    # legitimate DSF pavement polygons are within 2.2× apt.dat's
    # largest, well under the 3× cap.
    DSF_MAX_AREA_VS_APT_DAT_RATIO = 3.0
    apt_pav_union: Optional[Polygon] = None
    apt_pav_largest_area: float = 0.0
    if pav_polys:
        try:
            apt_pav_union = unary_union(pav_polys)
        except _GEOM_EXC:
            apt_pav_union = None
        apt_pav_largest_area = max(
            (p.area for p in pav_polys), default=0.0)

    # ── OSM-aeroway-footprint vs apt.dat coverage check ───────────
    # Per user 2026-04-29: prioritize apt.dat as the pavement
    # source.  Only fall back to DSF when there's a meaningful
    # discrepancy between apt.dat and what OSM aeroway data tells
    # us the airport actually has.  The OSM aeroway tags
    # (aeroway=apron, taxiway, taxi_lane, stand) are the user-
    # mapped truth about the airport's pavement extent.  If
    # apt.dat already covers that extent, DSF additions are at
    # best decorative overlays and at worst inflated ground tiles
    # filling non-pavement areas (HECA, where DSF added pavement
    # between runways and taxiways).  When apt.dat is missing
    # significant OSM-known pavement, DSF is admitted only inside
    # the gap.
    #
    # Load OSM here (early) so the DSF loop below can use the
    # aeroway footprint to gate DSF additions.
    nodes, ways, relations = _load_osm_airports(
        xplane_root, icao, anchor[0], anchor[1])
    osm_aeroway_footprint = _build_osm_aeroway_footprint(
        nodes, ways, to_m)
    # Gap = OSM-known pavement that apt.dat doesn't cover.  When
    # this is a small fraction, apt.dat is sufficient; skip DSF
    # entirely.
    osm_gap: Optional[Polygon] = None
    DSF_OSM_GAP_BUFFER_M = 1.0          # widen gap by 1 m for
                                         # tile-alignment slop.
                                         # User 2026-04-29 (HECA R
                                         # absorption): 5 m was too
                                         # generous; DSF clipped
                                         # with a 5 m fringe
                                         # extends ~3 m past the
                                         # OSM-tagged taxi corridor
                                         # and trips the long-edge-
                                         # adjacent absorption probe
                                         # (which fires at 2 m
                                         # outside the rect edge).
                                         # 1 m fringe matches the
                                         # apt.dat / DSF tile
                                         # alignment precision
                                         # without spilling enough
                                         # to look like apron-
                                         # adjacency.
    # Coverage threshold: above this, the airport's apt.dat is
    # considered "comprehensive" — apt.dat captures most of what
    # OSM thinks the airport has, so DSF is restricted to filling
    # the small remaining gap (typically a couple of taxi corridors
    # apt.dat happens to miss).  Below this threshold, apt.dat is
    # sparse (e.g. CYXY where apt.dat has ~52% of OSM) and OSM is
    # also incomplete; DSF is the primary pavement source there
    # and we trust it broadly (drop overlays only).
    APT_COMPREHENSIVE_OSM_FRAC = 0.80
    apt_is_comprehensive = False
    if (osm_aeroway_footprint is not None
            and not osm_aeroway_footprint.is_empty
            and apt_pav_union is not None
            and not apt_pav_union.is_empty):
        try:
            apt_in_osm = apt_pav_union.intersection(
                osm_aeroway_footprint).area
            osm_area = osm_aeroway_footprint.area
            if osm_area > 1.0:
                apt_is_comprehensive = (
                    apt_in_osm / osm_area
                    >= APT_COMPREHENSIVE_OSM_FRAC)
            gap = osm_aeroway_footprint.difference(apt_pav_union)
            if not gap.is_empty:
                osm_gap = gap.buffer(DSF_OSM_GAP_BUFFER_M)
        except _GEOM_EXC:
            pass
    # Compute the airport's bounding box from runway corners +
    # apt.dat pavement.  DSF polygons farther than
    # DSF_AIRPORT_RADIUS_M from this bbox are not this airport's.
    apt_bbox_m: Optional[Tuple[float, float, float, float]] = None
    bbox_polys = list(runway_polys) + list(pav_polys)
    if bbox_polys:
        try:
            uni = unary_union(bbox_polys)
            if not uni.is_empty:
                bx_min, by_min, bx_max, by_max = uni.bounds
                apt_bbox_m = (bx_min - DSF_AIRPORT_RADIUS_M,
                               by_min - DSF_AIRPORT_RADIUS_M,
                               bx_max + DSF_AIRPORT_RADIUS_M,
                               by_max + DSF_AIRPORT_RADIUS_M)
        except _GEOM_EXC:
            apt_bbox_m = None
    try:
        if not LOAD_DSF_PAVEMENT:
            raise StopIteration  # skip the DSF block entirely
        from . import dsf_reader as _DSFR
        seen_dsf: set = set()
        all_apt_dats = APR.find_all_airport_apt_dats(xplane_root, icao)
        n_dsf_kept = 0
        n_dsf_dropped_overlay = 0
        n_dsf_dropped_far = 0
        n_dsf_dropped_oversized = 0
        n_dsf_dropped_outside_osm_gap = 0
        for ad in all_apt_dats:
            dsf = _DSFR.find_associated_dsf(ad, anchor[0], anchor[1])
            if dsf is None or dsf in seen_dsf:
                continue
            seen_dsf.add(dsf)
            for ring in _DSFR.read_dsf_pavements(dsf):
                if len(ring) < 3:
                    continue
                try:
                    poly_ll = Polygon([(lon, lat) for (lon, lat) in ring])
                    if not poly_ll.is_valid:
                        poly_ll = poly_ll.buffer(0)
                    if (poly_ll.is_empty
                            or poly_ll.geom_type != "Polygon"):
                        continue
                    pm = shp_transform(to_m, poly_ll)
                    if pm.is_empty or pm.geom_type != "Polygon":
                        continue
                    # Distance gate: skip polygons outside this
                    # airport's expanded bbox.
                    if apt_bbox_m is not None:
                        px_min, py_min, px_max, py_max = pm.bounds
                        if (px_max < apt_bbox_m[0]
                                or px_min > apt_bbox_m[2]
                                or py_max < apt_bbox_m[1]
                                or py_min > apt_bbox_m[3]):
                            n_dsf_dropped_far += 1
                            continue
                    # apt.dat-priority gate (user 2026-04-29):
                    # when apt.dat is comprehensive (covers ≥ 80%
                    # of the OSM-aeroway footprint), CLIP each DSF
                    # polygon to the buffered OSM-vs-apt.dat gap.
                    # An intersect-test alone is too lax: a wide
                    # DSF polygon that grazes the gap by 1 m² gets
                    # kept entirely, and its 50–80 m extension
                    # past the gap drags non-pavement coverage
                    # into pav_union (HECA F case where the DSF
                    # along F's corridor extends 80 m onto the
                    # adjacent ramp, inflating ``_natural_half_
                    # width`` to ~50 m and triggering the apron-
                    # interior corner check on every F segment).
                    # Clipping keeps only the part of the DSF
                    # polygon that's actually filling an OSM-
                    # tagged corridor apt.dat happens to miss.
                    # When apt.dat is SPARSE (< 80 % of OSM), we
                    # skip the clip — apt.dat alone is too thin
                    # and OSM is also incomplete, so DSF is the
                    # primary source and we trust it broadly
                    # (CYXY, SPLP).
                    if (apt_is_comprehensive
                            and osm_gap is not None
                            and not osm_gap.is_empty):
                        try:
                            clipped_pm = pm.intersection(osm_gap)
                            if clipped_pm.is_empty:
                                n_dsf_dropped_outside_osm_gap += 1
                                continue
                            # Take the largest Polygon piece if
                            # the clip produced a MultiPolygon —
                            # narrow slivers from a wide DSF
                            # polygon grazing the gap aren't
                            # useful pavement either.
                            if (clipped_pm.geom_type
                                    == "MultiPolygon"):
                                clipped_pm = max(
                                    clipped_pm.geoms,
                                    key=lambda g: g.area)
                            if (clipped_pm.geom_type != "Polygon"
                                    or clipped_pm.is_empty
                                    or clipped_pm.area < 5.0):
                                n_dsf_dropped_outside_osm_gap += 1
                                continue
                            pm = clipped_pm
                        except _GEOM_EXC:
                            pass
                    # Oversized-vs-apt.dat gate: a DSF polygon
                    # dramatically larger than the airport's
                    # biggest apt.dat pavement polygon is a coarse
                    # "ground tile" overlay, not real pavement —
                    # drop it.  Only meaningful when apt.dat has
                    # ANY pavement; airports with no apt.dat
                    # pavement (CYXY-style sparse data) are
                    # unaffected.
                    if (apt_pav_largest_area > 0
                            and pm.area
                            > (apt_pav_largest_area
                               * DSF_MAX_AREA_VS_APT_DAT_RATIO)):
                        n_dsf_dropped_oversized += 1
                        continue
                    # Overlay check: drop the polygon if most of its
                    # area lies inside the existing apt.dat pavement
                    # union (it's a decorative overlay rather than
                    # new pavement).
                    if apt_pav_union is not None:
                        try:
                            inter_area = pm.intersection(
                                apt_pav_union).area
                            if (pm.area > 0
                                    and inter_area / pm.area
                                    >= DSF_OVERLAY_FRAC):
                                n_dsf_dropped_overlay += 1
                                continue
                        except _GEOM_EXC:
                            pass
                    pav_polys.append(pm)
                    n_dsf_kept += 1
                except _GEOM_EXC:
                    continue
        if (n_dsf_kept or n_dsf_dropped_overlay
                or n_dsf_dropped_far or n_dsf_dropped_oversized
                or n_dsf_dropped_outside_osm_gap):
            try:
                import sys as _sys
                msg = (f"  [pav-builder] {icao}: DSF pavement: "
                       f"{n_dsf_kept} kept, "
                       f"{n_dsf_dropped_overlay} dropped as overlay, "
                       f"{n_dsf_dropped_far} dropped as off-airport")
                if n_dsf_dropped_oversized:
                    msg += (f", {n_dsf_dropped_oversized} dropped "
                            f"as oversized-vs-apt.dat")
                if n_dsf_dropped_outside_osm_gap:
                    msg += (f", {n_dsf_dropped_outside_osm_gap} "
                            f"dropped: outside OSM-aeroway gap")
                UI.vprint(1, msg + ".")
            except _GEOM_EXC:
                pass
    except StopIteration:
        # DSF read intentionally disabled.
        try:
            UI.vprint(1,
                f"  [pav-builder] {icao}: DSF pavement read "
                f"disabled (LOAD_DSF_PAVEMENT=False).")
        except _GEOM_EXC:
            pass
    except _GEOM_EXC:
        pass
    pav_union = unary_union(pav_polys) if pav_polys else None
    # Merge near-touching apt.dat polygons so the union is one big
    # coverage (with real holes only) — see ``_merge_near_touching``.
    pav_union = _merge_near_touching(pav_union)
    # Per user 2026-05-05: simplify pav_union FIRST so all
    # downstream consumers (rect snap, junction emit) see a clean
    # 1 m-resolution coverage polygon.  Apt.dat row-110 polygons
    # routinely contain over-resolved curves (sub-meter steps) and
    # 0.2 m doubled-vertex needles that would otherwise survive
    # into the residue.  Rect corners snap to this simplified
    # boundary, so subtracting rects from pav_union should align
    # perfectly.
    from .pavement.union_helpers import _simplify_pavement_polygon
    pav_union = _simplify_pavement_polygon(pav_union, tol=1.0)
    # Stash the pre-runway-subtraction pavement polygon list for
    # the apron-merged-runway detection in _compute_elevations.
    # A runway segment is "apron-merged" when the apt.dat polygon
    # CONTAINING it is much larger than the segment itself —
    # apron pavement enclosing a runway is far wider than the
    # runway, while a normal runway lies inside a runway-shaped
    # apt.dat polygon that's only marginally larger than itself.
    apron_candidates = list(pav_polys)  # apt.dat + DSF, pre-subtract
    if pav_union is not None and layout.runway_union is not None:
        # Per user 2026-04-28: where a runway passes through a much
        # larger apron polygon, the runway is "apron-merged" — the
        # apron physically covers the runway pavement and the
        # downstream runway-segment-chain processing will drop the
        # apron-merged segments.  Don't subtract those parts from
        # pav_union now: the apron junctions should cover them
        # naturally, with no runway-shaped void to fill later.
        #
        # Detection mirrors ``_compute_elevations``'s segment-level
        # check (line ~3469) but applied to the original runway
        # polygons: the runway/candidate intersection counts as
        # apron-merged when the candidate is ≥
        # RUNWAY_APRON_AREA_RATIO × the intersection area.  A small
        # taxiway-sized candidate doesn't qualify (intersection is
        # most of the candidate); only big apron polygons do.
        apron_merged_regions: List[Polygon] = []
        for r_poly in runway_polys:
            for cand in apron_candidates:
                try:
                    inter = r_poly.intersection(cand)
                    if inter.is_empty or inter.area < 1.0:
                        continue
                    if cand.area > inter.area * RUNWAY_APRON_AREA_RATIO:
                        # Take the intersection as the apron-merged
                        # region — extracted as Polygon parts only.
                        if inter.geom_type == "Polygon":
                            apron_merged_regions.append(inter)
                        elif hasattr(inter, "geoms"):
                            for g in inter.geoms:
                                if (g.geom_type == "Polygon"
                                        and not g.is_empty):
                                    apron_merged_regions.append(g)
                except _GEOM_EXC:
                    continue
        if apron_merged_regions:
            try:
                merged_union = unary_union(apron_merged_regions)
                effective_runway = layout.runway_union.difference(
                    merged_union)
            except _GEOM_EXC:
                effective_runway = layout.runway_union
        else:
            effective_runway = layout.runway_union
        # Two pav_union variants:
        #   * ``pav_union_for_rects`` — full runway subtraction.
        #     Used by ``_build_taxi_rects`` for centerline clipping
        #     and the apron-interior boundary check.  Keeps F-style
        #     rects from extending into apron-merged-runway regions
        #     and failing the corner-on-boundary check (regression
        #     observed at CYXY's North F when residue switched to
        #     effective-runway subtraction).
        #   * ``pav_union`` (mutated below) — effective_runway
        #     subtraction so the residue / apron junctions cover
        #     apron-merged regions naturally.
        pav_union_for_rects = pav_union.difference(layout.runway_union)
        pav_union = pav_union.difference(effective_runway)
        layout._effective_runway_union = effective_runway
        layout._pav_union_for_rects = pav_union_for_rects

    # Collect all apt.dat pavement vertices (pre-union, real apt.dat
    # coord set) + runway corners.  This is the authoritative vertex
    # set the target snapper uses; rect corners will snap to these
    # preferentially so output shares vertices with target.
    apt_pav_vertices: List[Tuple[float, float]] = []
    for _pp in pav_polys:
        if _pp.is_empty or _pp.geom_type != "Polygon":
            continue
        _ec = list(_pp.exterior.coords)
        if _ec and _ec[0] == _ec[-1]:
            _ec = _ec[:-1]
        apt_pav_vertices.extend(_ec)
        for _ring in _pp.interiors:
            _rc = list(_ring.coords)
            if _rc and _rc[0] == _rc[-1]:
                _rc = _rc[:-1]
            apt_pav_vertices.extend(_rc)
    for _rp in runway_polys:
        _rc = list(_rp.exterior.coords)
        if _rc and _rc[0] == _rc[-1]:
            _rc = _rc[:-1]
        apt_pav_vertices.extend(_rc)

    # ── Taxi centerlines (apt.dat primary, OSM fallback) ─────────
    # Per user 2026-05-12: apt.dat row 1201/1202 taxi-network is
    # the authoritative source for the taxi graph at airports that
    # have it.  Since the pavement footprint is also drawn from
    # apt.dat row-110 polygons, the taxi-network endpoints align
    # exactly with the pavement boundary — eliminating the OSM-vs-
    # apt.dat boundary mismatch where OSM centerlines clipped
    # against ``pav_union`` produced empty / too-short intersections
    # (CYXY's long E parallel, all of F and G — the user's
    # "taxiways turning into big junctions" report).  Fall back to
    # OSM only when the apt.dat block has no taxi-network at all
    # (some custom packs omit rows 1201/1202).
    apt_centerlines = APR.taxi_centerlines(
        apt, to_m, rwy_centerlines=rwy_centerlines)
    if apt_centerlines:
        osm_centerlines = apt_centerlines
        UI.vprint(1,
            f"  [pav-builder] {icao}: using {len(apt_centerlines)} "
            f"apt.dat taxi-network centerline(s) "
            f"({len(apt.taxi_nodes)} nodes, "
            f"{len(apt.taxi_edges)} edges).")
    else:
        osm_centerlines = _extract_osm_taxi_centerlines(
            nodes, ways, to_m, rwy_centerlines=rwy_centerlines)
        UI.vprint(1,
            f"  [pav-builder] {icao}: apt.dat has no taxi network; "
            f"using {len(osm_centerlines)} OSM aeroway-taxiway "
            f"centerline(s).")

    # ── Terminal groundside-pavement subtraction (user 2026-04-29):
    # remove curbside / drop-off / parking pavement from pav_union
    # before downstream rect / junction construction sees it.
    # Groundside pavement sits at a different elevation than the
    # building's airside apron, so allowing it to become an apron
    # junction grade-clamps the building to the wrong altitude.
    # Subtract a perpendicular outward strip from each terminal
    # building's groundside edges (classified by OSM aeroway /
    # highway adjacency + apt.dat-pavement connectivity).
    try:
        _osm_terminal_buildings = _extract_osm_terminals(
            nodes, ways, relations, to_m)
        # Per user 2026-04-30 (CYXY -10123 NW phantom groundside):
        # pass the FULL apt.dat pavement polygon list so the
        # classifier can BFS from runway-touching polys through
        # transitive touches.  Without this, edges of the
        # terminal next to apron pavement that's not directly
        # touching a runway (e.g. the apron extends NW past the
        # terminal) classified UNKNOWN and got mis-promoted to
        # groundside.
        _ground_zone = _terminal_groundside_zone(
            _osm_terminal_buildings, nodes, ways, to_m,
            apt_pavement_seeds=runway_polys,
            apt_pavement_polys=apt_only_pav_polys)
        if _ground_zone is not None and not _ground_zone.is_empty:
            # Capture the groundside-only visible pavement BEFORE
            # the subtraction below empties pav_union of it.  Per
            # user 2026-04-29: groundside pavement should remain
            # in the output but follow DEM (with a 0.1 m gap from
            # the terminal building) rather than being flattened
            # to airside-apron elevation.  The shapes captured
            # here are emitted later with per-vertex DEM altitudes;
            # the 0.1 m terminal gap is enforced inside the emit
            # function using the LAYOUT's terminal shapes (which
            # may differ slightly from the OSM source extracts due
            # to apt.dat row-110 / DSF residue absorption).
            try:
                _groundside_visible = pav_union.intersection(
                    _ground_zone)
                _gs_polys: List[Polygon] = []
                if _groundside_visible is not None:
                    if _groundside_visible.geom_type == "Polygon":
                        if (not _groundside_visible.is_empty
                                and _groundside_visible.area >= 5.0):
                            _gs_polys.append(_groundside_visible)
                    elif (_groundside_visible.geom_type
                            == "MultiPolygon"):
                        for _g in _groundside_visible.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 5.0):
                                _gs_polys.append(_g)
                # Stash on the layout so the elevation pass can
                # find them once DEM is loaded.
                layout._groundside_polys = _gs_polys
            except _GEOM_EXC:
                layout._groundside_polys = []
            try:
                pav_union = pav_union.difference(_ground_zone)
                if hasattr(layout, "_pav_union_for_rects"):
                    layout._pav_union_for_rects = (
                        layout._pav_union_for_rects.difference(
                            _ground_zone))
                # Also subtract from the granular polygon list so
                # downstream apron-merged-runway / rect-corner
                # detection sees a consistent pavement footprint.
                _new_pav_polys: List[Polygon] = []
                for _p in pav_polys:
                    try:
                        _q = _p.difference(_ground_zone)
                    except _GEOM_EXC:
                        _new_pav_polys.append(_p)
                        continue
                    if _q.is_empty:
                        continue
                    if _q.geom_type == "Polygon":
                        _new_pav_polys.append(_q)
                    elif _q.geom_type == "MultiPolygon":
                        for _g in _q.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 1.0):
                                _new_pav_polys.append(_g)
                pav_polys[:] = _new_pav_polys
                # Same for apt_only_pav_polys (used for terminal
                # containment + rect-corner snapping).
                _new_apt_only: List[Polygon] = []
                for _p in apt_only_pav_polys:
                    try:
                        _q = _p.difference(_ground_zone)
                    except _GEOM_EXC:
                        _new_apt_only.append(_p)
                        continue
                    if _q.is_empty:
                        continue
                    if _q.geom_type == "Polygon":
                        _new_apt_only.append(_q)
                    elif _q.geom_type == "MultiPolygon":
                        for _g in _q.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 1.0):
                                _new_apt_only.append(_g)
                apt_only_pav_polys[:] = _new_apt_only
                # Also subtract from apron_candidates — captured at
                # line 2519 BEFORE this subtract — so apron-junction
                # construction in _compute_elevations can't wrap
                # around the terminal into the groundside zone (per
                # user 2026-04-29 / way -10125 vs -10111: the apron
                # junction was extending past the terminal to the
                # upper-side roads, sharing an edge with the new
                # DEM-following groundside pavement and creating a
                # 7 m vertical cliff at CYXY).
                _new_apron_cand: List[Polygon] = []
                for _p in apron_candidates:
                    try:
                        _q = _p.difference(_ground_zone)
                    except _GEOM_EXC:
                        _new_apron_cand.append(_p)
                        continue
                    if _q.is_empty:
                        continue
                    if _q.geom_type == "Polygon":
                        _new_apron_cand.append(_q)
                    elif _q.geom_type == "MultiPolygon":
                        for _g in _q.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 1.0):
                                _new_apron_cand.append(_g)
                apron_candidates[:] = _new_apron_cand
                try:
                    UI.vprint(1,
                        f"  [pav-builder] {icao}: subtracted "
                        f"{_ground_zone.area:,.0f} m² of "
                        f"groundside pavement (terminal "
                        f"curbside / drop-off / parking).")
                except _GEOM_EXC:
                    pass
            except _GEOM_EXC:
                pass
    except _GEOM_EXC:
        pass

    # ── Per-ref OVERALL chord bearings (pre-split) ───────────────
    # Used by ``_classify_role`` to disambiguate diagonal-overall
    # taxis whose curving ends happen to align near-parallel to
    # the runway locally.  Without this, a B/C/E/G stub at SPJC —
    # which enters the runway at a shallow angle — gets a small
    # post-curve segment classified as PRIMARY_PARALLEL because
    # the segment's local bearing falls inside the 20° parallel
    # window, even though the OVERALL B/C/E/G chord is diagonal.
    # The parent's overall chord bearing is the right reference.
    ref_overall_bearings: Dict[str, float] = {}
    _ref_longest_len: Dict[str, float] = {}
    for _ls, _ref in osm_centerlines:
        if not _ref:
            continue
        if _ls.length <= _ref_longest_len.get(_ref, 0.0):
            continue
        _coords = list(_ls.coords)
        if len(_coords) < 2:
            continue
        _dx = _coords[-1][0] - _coords[0][0]
        _dy = _coords[-1][1] - _coords[0][1]
        if math.hypot(_dx, _dy) < 1.0:
            continue
        ref_overall_bearings[_ref] = (
            math.degrees(math.atan2(_dx, _dy)) % 180.0)
        _ref_longest_len[_ref] = _ls.length

    # ── Primary-parallel SPINE LINES (pre-split) ─────────────────
    # Per user 2026-04-27: identify each primary-parallel taxiway
    # (overall db < 20° to the nearest runway) and build a single
    # extended SPINE LINE running through it.  Diagonal stubs
    # (B/C/D/E/G at SPJC) should END where they cross this spine —
    # not where the diagonal's OSM polyline happens to terminate
    # (which can be inside the apron, past where the parallel
    # ought to "continue" through).  The spine line lets us
    # imagine the parallel's centerline continuing through gaps
    # / aprons in the OSM data, providing a stable inland-side
    # bound for diagonal-stub trimming.
    parallel_spines: List[LineString] = []
    parallel_corridors: List[Tuple[str, Polygon]] = []
    SPINE_EXTEND_M = 800.0  # extend each spine ±800 m past its
                              # OSM-fragment endpoints so it acts
                              # as a guide line through aprons.
    PARALLEL_CORRIDOR_HALF_WIDTH_M = 30.0
                              # half-width of the imagined
                              # primary-parallel pavement corridor.
                              # Used to TRIM diagonal stubs at the
                              # corridor's runway-facing edge — so
                              # B/C/D/E at SPJC stop where they
                              # enter A's pavement (real or
                              # imagined-via-apron) rather than
                              # extending deep into the apron.
                              # 30 m is wider than a typical taxi
                              # half-width (22 m) so the corridor
                              # forgives small OSM/apt.dat
                              # misalignment.
    if rwy_centerlines:
        # Bearing of the FIRST runway centerline; we use it to
        # decide whether a ref qualifies as a primary parallel
        # (db < 20°).  All SPJC runways are parallel so any
        # runway works as the reference.
        _r0 = rwy_centerlines[0]
        _rc = list(_r0.coords)
        _rdx = _rc[-1][0] - _rc[0][0]
        _rdy = _rc[-1][1] - _rc[0][1]
        if math.hypot(_rdx, _rdy) > 1e-6:
            _rwy_bearing = (
                math.degrees(math.atan2(_rdx, _rdy)) % 180.0)
            for _ref, _bearing in ref_overall_bearings.items():
                _db = abs(_bearing - _rwy_bearing)
                _db = min(_db, 180.0 - _db)
                if _db >= 20.0:
                    continue  # not a primary parallel
                # Find the longest OSM centerline for this ref;
                # use its endpoints to define the spine direction
                # and base position.
                best_ls: Optional[LineString] = None
                best_len = 0.0
                for _ls, _r2 in osm_centerlines:
                    if _r2 != _ref:
                        continue
                    if _ls.length > best_len:
                        best_len = _ls.length
                        best_ls = _ls
                if best_ls is None or best_len < 100.0:
                    continue
                _ec = list(best_ls.coords)
                _eax, _eay = _ec[0]
                _ebx, _eby = _ec[-1]
                _edx = _ebx - _eax
                _edy = _eby - _eay
                _emag = math.hypot(_edx, _edy)
                if _emag < 1e-6:
                    continue
                _ux = _edx / _emag
                _uy = _edy / _emag
                _start = (_eax - _ux * SPINE_EXTEND_M,
                          _eay - _uy * SPINE_EXTEND_M)
                _end = (_ebx + _ux * SPINE_EXTEND_M,
                        _eby + _uy * SPINE_EXTEND_M)
                try:
                    _spine = LineString([_start, _end])
                    parallel_spines.append(_spine)
                    # Build the corridor polygon: spine buffered to
                    # PARALLEL_CORRIDOR_HALF_WIDTH_M, square ends so
                    # the corridor's apron-facing extension stays
                    # rectangular.
                    _corridor = _spine.buffer(
                        PARALLEL_CORRIDOR_HALF_WIDTH_M,
                        cap_style=2, join_style=2)
                    if (not _corridor.is_empty
                            and _corridor.geom_type == "Polygon"):
                        parallel_corridors.append((_ref, _corridor))
                except _GEOM_EXC:
                    pass

    # ── Pavement source-of-truth (user 2026-04-28): apt.dat row-110
    # ∪ DSF pavement, period.  Earlier revisions augmented pav_union
    # with a 30 m-wide synthetic buffer around any OSM centerline
    # that didn't intersect apt.dat/DSF; the rationale was to keep
    # rect-extraction working at airports where the OSM taxiway
    # network is more complete than the apt.dat coverage.  That
    # workaround is dropped: OSM centerlines drive WHICH taxiways
    # exist (geometry, ref tag, role), but the actual pavement
    # surface comes from apt.dat ∪ DSF only.  Centerlines without
    # matching apt.dat/DSF coverage produce no rect — that's an
    # apt.dat data gap to be fixed at the source, not papered over
    # with a synthetic strip whose width arbitrarily differs from
    # the OSM-tagged taxi width.

    # ── Terminals: expand OSM building outlines to the containing
    # apt.dat pavement polygon (or buffer if no polygon contains).
    # The target terminal is the "pad" — apt.dat pavement up to the
    # apron boundary.  OSM aeroway=terminal gives the building
    # footprint; we use that as a seed.
    osm_terminal_polys = _extract_osm_terminals(
        nodes, ways, relations, to_m)
    # Min-spacing simplification: vertices closer than this to a
    # neighbour are redundant for the airport-scale render and only
    # serve to spawn sliver triangles in the eventual ear-clip.
    # Applied to terminals (curved building footprints often
    # inherit closely-spaced OSM vertices) and to the junction
    # boundaries below.
    MIN_VERTEX_SPACING_M = 2.0
    terminal_polys: List[Polygon] = []
    for otp in osm_terminal_polys:
        # Apt.dat-only candidates — DSF polygons (overlays, gap
        # fills) shouldn't compete for terminal-pad selection.
        pad = _terminal_pad_from_building(otp, apt_only_pav_polys)
        if pad is None:
            continue
        try:
            simp = pad.simplify(
                MIN_VERTEX_SPACING_M, preserve_topology=True)
            if (simp.geom_type == "Polygon"
                    and not simp.is_empty
                    and simp.area >= 100.0):
                pad = simp
        except _GEOM_EXC:
            pass
        terminal_polys.append(pad)
    terminal_union = (unary_union(terminal_polys)
                      if terminal_polys else None)
    for i, tp in enumerate(terminal_polys):
        layout.shapes.append(BuiltShape(
            polygon=tp, role=ROLE_TERMINAL, ref=f"terminal{i+1}"))

    # ── Identify junction node CLUSTERS ──────────────────────────
    # Per user 2026-05-12: when the taxi graph comes from apt.dat
    # row 1201/1202, the junction nodes ARE explicit in the data —
    # any taxi node referenced by edges of ≥ 2 distinct names, or
    # by a runway-cross edge, is a chart-level junction.  Use
    # those points so ``_split_centerlines_at_points`` trims the
    # apt.dat centerlines at the right places (otherwise long
    # parallel taxis terminate deep inside apron polygons and
    # ``_snap_corners_to_pavement`` rejects them as
    # apron-interior).  Fall back to OSM topology when apt.dat
    # has no taxi network.
    if apt_centerlines:
        junction_points = APR.taxi_junction_points(apt, to_m)
        # Also add geometric crossing points between centerlines of
        # different names — covers the case where a taxi node was
        # tagged with one name but the connecting edge has another.
        # Same behaviour ``_find_junction_points`` provides for OSM.
        for i in range(len(osm_centerlines)):
            ls1, ref1 = osm_centerlines[i]
            for j in range(i + 1, len(osm_centerlines)):
                ls2, ref2 = osm_centerlines[j]
                if ref1 and ref2 and ref1 == ref2:
                    continue
                try:
                    if not ls1.intersects(ls2):
                        continue
                    inter = ls1.intersection(ls2)
                except _GEOM_EXC:
                    continue
                if inter.is_empty:
                    continue
                if inter.geom_type == "Point":
                    junction_points.append((inter.x, inter.y))
                elif inter.geom_type == "MultiPoint":
                    for p in inter.geoms:
                        junction_points.append((p.x, p.y))
    else:
        # Any OSM node referenced by ≥2 taxi ways is a potential
        # junction point.  Nodes within JUNCTION_CLUSTER_DIST of
        # each other are merged into one cluster (target junctions
        # often span a whole multi-way intersection, not just a
        # single OSM node).
        junction_points = _find_junction_points(
            nodes, ways, to_m, osm_centerlines=osm_centerlines)

    # ── Diagonal-stub trim at primary-parallel SPINES ────────────
    # Per user 2026-04-27: a diagonal stub (B/C/D/E/G overall db
    # ∈ [20°, 45°)) should END where its centerline crosses the
    # nearest primary-parallel SPINE LINE — even when the spine
    # passes through an apron with no OSM coverage.  Without this
    # the stub's apron-side OSM endpoint can sit deep inside the
    # apron, producing an over-long rect that the surrounding
    # junction has to wrap around.  Adding the spine intersection
    # to ``junction_points`` lets the existing
    # ``_split_centerlines_at_points`` machinery clip the diagonal
    # at the right place using its standard junction-margin rule.
    if parallel_spines and rwy_centerlines:
        _rwy_first = rwy_centerlines[0]
        _rfc = list(_rwy_first.coords)
        _rfdx = _rfc[-1][0] - _rfc[0][0]
        _rfdy = _rfc[-1][1] - _rfc[0][1]
        _rfmag = math.hypot(_rfdx, _rfdy)
        if _rfmag > 1e-6:
            _rfb = (math.degrees(math.atan2(_rfdx, _rfdy))
                    % 180.0)
            for _ls, _ref in osm_centerlines:
                if not _ref:
                    continue
                _b = ref_overall_bearings.get(_ref)
                if _b is None:
                    continue
                _db = abs(_b - _rfb)
                _db = min(_db, 180.0 - _db)
                # Only diagonal stubs (db ∈ [20°, 45°)) need this.
                if not (20.0 <= _db < 45.0):
                    continue
                # Find the closest spine and intersect.
                best_pt = None
                best_d = float("inf")
                for _spine in parallel_spines:
                    try:
                        _x = _ls.intersection(_spine)
                    except _GEOM_EXC:
                        continue
                    if _x.is_empty:
                        # Try extending the centerline ends to reach
                        # the spine — handles diagonals whose OSM
                        # path stops short of the spine's location.
                        continue
                    if _x.geom_type == "Point":
                        _pt = (_x.x, _x.y)
                    elif _x.geom_type == "MultiPoint":
                        # Pick the intersection furthest from the
                        # closest runway centerline — that's the
                        # apron-side cut we want for trimming.
                        _pt = None
                        _far = -1.0
                        for _g in _x.geoms:
                            _gp = (_g.x, _g.y)
                            try:
                                _gd = min(
                                    Point(_gp).distance(_r)
                                    for _r in rwy_centerlines)
                            except _GEOM_EXC:
                                _gd = 0.0
                            if _gd > _far:
                                _far = _gd
                                _pt = _gp
                        if _pt is None:
                            continue
                    else:
                        continue
                    # Distance along centerline from runway end —
                    # used to pick the closest spine intersection.
                    try:
                        _proj = _ls.project(Point(_pt))
                    except _GEOM_EXC:
                        continue
                    _d = abs(_proj - _ls.length / 2)
                    if _d < best_d:
                        best_d = _d
                        best_pt = _pt
                if best_pt is not None:
                    junction_points.append(best_pt)

    # Trim runway-approaching centerlines at a BUFFERED runway
    # polygon so the last rect on a taxi that crosses or enters
    # the runway stops short of the runway-junction approach.
    # The runway boundary IS an intersection (rule 70 % of gap
    # between intersections), but the physical junction — where
    # the taxi widens into the runway apron — extends some
    # distance OUTSIDE the runway polygon too.  Pulling the
    # centerline end back from the approach widening keeps the
    # rect from "encroaching into intersections" (user
    # 2026-04-21).
    #
    # Two thresholds:
    #   * Perpendicular (perp_diff < 25°): pull back 30 m.  These
    #     taxis hit the runway head-on and the widening is large.
    #   * Diagonal (25° ≤ perp_diff < 70°): pull back 15 m.
    #     The widening is gentler at oblique angles (B/C/D/E/G at
    #     SPJC) but still significant.  Per user 2026-04-27: "B,
    #     C, D, E, G should be handled like V3 stub" — V3 is a
    #     near-perpendicular sub-ref and gets the 30 m buffer; the
    #     diagonals get a smaller buffer so they don't over-shrink
    #     while still getting the same kind of pull-back.
    RWY_JUNCTION_BUFFER_M = 30.0
    RWY_DIAG_BUFFER_M = 15.0
    PERP_TRIM_MAX_DEG = 25.0
    DIAG_TRIM_MAX_DEG = 70.0
    if (layout.runway_union is not None
            and not layout.runway_union.is_empty
            and rwy_centerlines):
        try:
            rwy_buffered = layout.runway_union.buffer(
                RWY_JUNCTION_BUFFER_M)
        except _GEOM_EXC:
            rwy_buffered = layout.runway_union
        # Nearest-runway bearing for angle check
        def _perp_diff_to_runway(ls: LineString) -> float:
            c = list(ls.coords)
            if len(c) < 2:
                return 90.0
            dx = c[-1][0] - c[0][0]
            dy = c[-1][1] - c[0][1]
            mag = math.hypot(dx, dy)
            if mag < 1e-6:
                return 90.0
            ax_bearing = math.degrees(
                math.atan2(dx, dy)) % 180.0
            mid = ls.interpolate(ls.length / 2)
            best = min(rwy_centerlines,
                       key=lambda r: mid.distance(r))
            rc = list(best.coords)
            rx = rc[-1][0] - rc[0][0]
            ry = rc[-1][1] - rc[0][1]
            rmag = math.hypot(rx, ry)
            if rmag < 1e-6:
                return 90.0
            rwy_bearing = math.degrees(
                math.atan2(rx, ry)) % 180.0
            delta = abs(ax_bearing - rwy_bearing)
            delta = min(delta, 180.0 - delta)
            return abs(delta - 90.0)

        trimmed_centerlines: List[Tuple[LineString, str]] = []
        for ls, ref in osm_centerlines:
            # Only PERPENDICULAR centerlines (perp_diff < 25°)
            # need this 30 m buffer pull-back.  Diagonals are
            # handled by the SEPARATE corridor-trim pass below
            # (which subtracts the imagined primary-parallel
            # corridor from the diagonal's apron-side end).
            # Stacking both gave a too-aggressive shortening at
            # SPJC (B 203 m → 92 m).
            pd = _perp_diff_to_runway(ls)
            if pd >= PERP_TRIM_MAX_DEG:
                trimmed_centerlines.append((ls, ref))
                continue
            _buf = rwy_buffered
            try:
                diff = ls.difference(_buf)
            except _GEOM_EXC:
                trimmed_centerlines.append((ls, ref))
                continue
            if diff.is_empty:
                trimmed_centerlines.append((ls, ref))
                continue
            if diff.geom_type == "LineString":
                if diff.length >= MIN_SEGMENT_LEN_M:
                    trimmed_centerlines.append((diff, ref))
                else:
                    trimmed_centerlines.append((ls, ref))
            elif diff.geom_type == "MultiLineString":
                longest = max(diff.geoms, key=lambda g: g.length)
                if longest.length >= MIN_SEGMENT_LEN_M:
                    trimmed_centerlines.append((longest, ref))
                else:
                    trimmed_centerlines.append((ls, ref))
            else:
                trimmed_centerlines.append((ls, ref))
        osm_centerlines = trimmed_centerlines

    # Trim PERPENDICULAR centerlines at a BUFFERED PARALLEL-
    # CENTERLINE polygon.  Per user (2026-04-22): at SPLP
    # (unrefed airport), perpendicular stubs come from long
    # multipurpose taxis whose gap is bounded by a POINT
    # crossing with the primary axis.  That crossing doesn't
    # account for the primary's WIDTH — the 70 % rect then
    # extends into the primary's rect zone, reading as
    # "off-center toward the taxiway".  A 15 m buffer around
    # each parallel-to-runway centerline (= primary half-width)
    # pulls the perpendicular stub's primary-side endpoint
    # out to the primary's physical edge.  SPJC doesn't need
    # this because each sub-ref stub has its own dedicated OSM
    # way — the geometry inherently excludes the primary
    # width.
    PARALLEL_BUFFER_M = 15.0
    PARALLEL_MIN_LEN_M = 200.0
    if rwy_centerlines:
        parallel_polys = []
        for ls, _ref in osm_centerlines:
            if ls.length < PARALLEL_MIN_LEN_M:
                continue
            # Check if this centerline is PARALLEL to runway
            # (perp_diff > 75°).
            c = list(ls.coords)
            if len(c) < 2:
                continue
            dx = c[-1][0] - c[0][0]
            dy = c[-1][1] - c[0][1]
            mag = math.hypot(dx, dy)
            if mag < 1e-6:
                continue
            ax_bearing = math.degrees(
                math.atan2(dx, dy)) % 180.0
            mid = ls.interpolate(ls.length / 2)
            best_r = min(rwy_centerlines,
                         key=lambda r: mid.distance(r))
            rc = list(best_r.coords)
            rx = rc[-1][0] - rc[0][0]
            ry = rc[-1][1] - rc[0][1]
            rmag = math.hypot(rx, ry)
            if rmag < 1e-6:
                continue
            rwy_bearing = math.degrees(
                math.atan2(rx, ry)) % 180.0
            delta = abs(ax_bearing - rwy_bearing)
            delta = min(delta, 180.0 - delta)
            perp_diff = abs(delta - 90.0)
            if perp_diff > 75.0:
                # Parallel to runway → create buffer polygon
                try:
                    parallel_polys.append(
                        ls.buffer(PARALLEL_BUFFER_M))
                except _GEOM_EXC:
                    pass
        if parallel_polys:
            try:
                parallel_union = unary_union(parallel_polys)
            except _GEOM_EXC:
                parallel_union = None
            if parallel_union is not None and not parallel_union.is_empty:
                def _perp_diff_to_rwy(ls):
                    c = list(ls.coords)
                    if len(c) < 2:
                        return 90.0
                    dx = c[-1][0] - c[0][0]
                    dy = c[-1][1] - c[0][1]
                    mag = math.hypot(dx, dy)
                    if mag < 1e-6:
                        return 90.0
                    ax_b = math.degrees(
                        math.atan2(dx, dy)) % 180.0
                    mid_p = ls.interpolate(ls.length / 2)
                    best = min(rwy_centerlines,
                               key=lambda r: mid_p.distance(r))
                    rc = list(best.coords)
                    rx = rc[-1][0] - rc[0][0]
                    ry = rc[-1][1] - rc[0][1]
                    rmag = math.hypot(rx, ry)
                    if rmag < 1e-6:
                        return 90.0
                    rw_b = math.degrees(
                        math.atan2(rx, ry)) % 180.0
                    dlt = abs(ax_b - rw_b)
                    dlt = min(dlt, 180.0 - dlt)
                    return abs(dlt - 90.0)

                trimmed_perp: List[Tuple[LineString, str]] = []
                for ls, ref in osm_centerlines:
                    # Trim perpendicular AND diagonal centerlines
                    # — both cross the parallel-parallel corridor
                    # and benefit from starting at the primary
                    # EDGE not axis.  Skip only near-parallel ones
                    # (perp_diff >= 75°) which are themselves
                    # primaries.
                    if _perp_diff_to_rwy(ls) >= 75.0:
                        trimmed_perp.append((ls, ref))
                        continue
                    try:
                        diff = ls.difference(parallel_union)
                    except _GEOM_EXC:
                        trimmed_perp.append((ls, ref))
                        continue
                    if diff.is_empty:
                        trimmed_perp.append((ls, ref))
                        continue
                    if diff.geom_type == "LineString":
                        if diff.length >= MIN_SEGMENT_LEN_M:
                            trimmed_perp.append((diff, ref))
                        else:
                            trimmed_perp.append((ls, ref))
                    elif diff.geom_type == "MultiLineString":
                        # Multiple pieces — keep ALL so each stub
                        # between parallels emits separately.
                        for g in diff.geoms:
                            if g.length >= MIN_SEGMENT_LEN_M:
                                trimmed_perp.append((g, ref))
                    else:
                        trimmed_perp.append((ls, ref))
                osm_centerlines = trimmed_perp

    # ── Diagonal-stub corridor trim ──────────────────────────────
    # Per user 2026-04-27: when a diagonal stub (B/C/D/E/G at SPJC)
    # connects directly to a large apron rather than crossing a
    # primary parallel's actual rect, the stub's apron-end ends up
    # floating in the apron (no nearby pavement edge to snap to).
    # The fix: subtract the IMAGINED PRIMARY PARALLEL CORRIDOR
    # (spine ± PARALLEL_CORRIDOR_HALF_WIDTH_M, extended through
    # gaps in the OSM coverage) from each diagonal stub's
    # centerline.  The diagonal then ends at the corridor's
    # runway-facing edge — exactly where the imagined primary
    # parallel's "near" pavement boundary sits.  Downstream rect-
    # build + corner-snap-to-pavement then puts the rect's apron-
    # end corners on the pavement edge between runway and apron.
    if parallel_corridors and ref_overall_bearings:
        _r0c = list(rwy_centerlines[0].coords)
        _r0db = math.degrees(math.atan2(
            _r0c[-1][0] - _r0c[0][0],
            _r0c[-1][1] - _r0c[0][1])) % 180.0
        corridor_trimmed: List[Tuple[LineString, str]] = []
        for ls, ref in osm_centerlines:
            if not ref:
                corridor_trimmed.append((ls, ref))
                continue
            _b = ref_overall_bearings.get(ref)
            if _b is None:
                corridor_trimmed.append((ls, ref))
                continue
            _db_local = abs(_b - _r0db)
            _db_local = min(_db_local, 180.0 - _db_local)
            # Only diagonal stubs (overall db ∈ [20°, 45°)) need
            # this corridor trim.  Sub-refs (V3 etc) would also
            # qualify and benefit, so include them.
            if not (20.0 <= _db_local < 45.0):
                corridor_trimmed.append((ls, ref))
                continue
            # Subtract every OTHER ref's corridor from this
            # centerline.  Skip the centerline's own corridor
            # (so e.g. an A1 stub doesn't subtract A's corridor
            # from itself if it was somehow classified as
            # diagonal).
            current = ls
            for c_ref, corridor in parallel_corridors:
                if c_ref == ref:
                    continue
                try:
                    diff = current.difference(corridor)
                except _GEOM_EXC:
                    continue
                if diff.is_empty:
                    continue
                if diff.geom_type == "LineString":
                    if diff.length >= MIN_SEGMENT_LEN_M:
                        current = diff
                elif diff.geom_type == "MultiLineString":
                    longest = max(diff.geoms,
                                  key=lambda g: g.length)
                    if longest.length >= MIN_SEGMENT_LEN_M:
                        current = longest
            corridor_trimmed.append((current, ref))
        osm_centerlines = corridor_trimmed

    # Per user rule (2026-04-18): "Implement splitting at cross ref
    # crossings" — split every centerline at each multi-ref
    # junction node along its path, so each straight section
    # between intersections emits as a single rect.
    # Split per user rules 1+4 (2026-04-20):
    #   Rule 1: stop rects at intersections; cluster close
    #           intersections into single junction region.
    #   Rule 4: cover only narrowest straight sections; widened
    #           pavement belongs to junctions.
    # Uses OSM multi-ref nodes as intersection anchors, gated by
    # gap + pav-width clustering at the midpoint.
    # Uniform pipeline per user (2026-04-20): intersections +
    # sharp curves define break points; 70% rect between
    # consecutive breaks.  No width analysis.
    # Per user 2026-05-05: use the SAME pav_union for rect building
    # and junction emit.  Previously rects used ``pav_union_for_rects``
    # (full runway subtraction) while junctions used ``pav_union``
    # (effective_runway subtraction, retaining apron-merged-runway
    # pavement).  The split was added to fix CYXY's F primary
    # parallel, but it caused SPJC rect corners to land on
    # pav_for_rects.boundary that doesn't exist on pav_union.boundary
    # — leaving thin tabs in the residue.  Single-source-of-truth
    # boundary lets corners snap consistently.
    osm_centerlines = _split_centerlines_at_points(
        osm_centerlines, junction_points, approach_tol_m=25.0,
        pav_union=pav_union, rwy_union=layout.runway_union,
        rwy_centerlines=rwy_centerlines)

    # ── Build taxi rects from centerlines ────────────────────────
    taxi_rects = _build_taxi_rects(
        osm_centerlines, pav_union, layout.runway_union,
        rwy_centerlines, apt_vertices=apt_pav_vertices,
        ref_overall_bearings=ref_overall_bearings)

    # Filter stubs by user's runway-connection rule: a stub rect is
    # kept only if its OSM centerline reaches a runway.  Stubs whose
    # pavement is only between apron and parallel (or apron-internal)
    # get dropped — that pavement folds into the apron or junction.
    # This removes OSM sub-refs like A1-A6, D1/D2, F1, M1-M3, R1/R2, N
    # that the user's target doesn't emit.
    # Look up raw OSM ways per ref (untrimmed).  A stub "connects to
    # a runway" iff one of its raw endpoints lies within
    # RUNWAY_ENDPOINT_DIST_M of the runway footprint.  The threshold
    # accounts for the junction polygon between the stub rect and
    # the runway: target stub-rect corners are 22–66 m from the
    # runway edge (measured from the actual target), so 80 m gives
    # a small margin for freehand drift.
    #
    # Per user 2026-05-11: SKIP this filter when the taxi graph
    # comes from apt.dat (rows 1201/1202).  apt.dat doesn't carry
    # OSM-noise sub-refs; every row-1202 edge is a canonical taxi
    # path.  CYXY's G runs apron-to-apron without ever touching
    # a runway — under OSM-only filtering it gets dropped here as
    # "apron-internal" even though apt.dat declares it a real
    # taxiway.  The downstream absorption pass
    # (``_drop_primary_parallels_embedded_in_pavement``) is the
    # authoritative place to decide whether an apron-running taxi
    # rect should be partially absorbed or kept intact.
    RUNWAY_ENDPOINT_DIST_M = 80.0
    raw_endpoints_by_ref: Dict[str, List[Tuple[float, float]]] = {}
    if not apt_centerlines:
        for wid, nds, tags in ways:
            if tags.get("aeroway") != "taxiway":
                continue
            ref = tags.get("ref", "")
            pts = []
            for n in nds:
                if n in nodes:
                    lat, lon = nodes[n]
                    pts.append(to_m(lon, lat))
            if len(pts) >= 2:
                raw_endpoints_by_ref.setdefault(ref, []).append(pts[0])
                raw_endpoints_by_ref.setdefault(ref, []).append(pts[-1])

    if not apt_centerlines and layout.runway_union is not None:
        rwy_boundary = layout.runway_union.boundary
        filtered: List[Tuple[Polygon, LineString, str, str]] = []
        for rect, axis, role, ref in taxi_rects:
            if role != ROLE_STUB:
                filtered.append((rect, axis, role, ref))
                continue
            # Reject stubs whose RECT touches (or sits inside) the
            # runway polygon.  SPLP has a short unrefed taxi at the
            # SW end whose rect at (-463,-1266) has one corner
            # exactly on the runway boundary — user wants "only a
            # single stub at the south end" so this one must go.
            try:
                if rect.distance(rwy_boundary) < 5.0 and (
                        rect.intersects(layout.runway_union)
                        or rect.distance(layout.runway_union) < 2.0):
                    continue
            except _GEOM_EXC:
                pass
            reaches_runway = False
            for (px, py) in raw_endpoints_by_ref.get(ref, []):
                if Point(px, py).distance(
                        layout.runway_union) <= RUNWAY_ENDPOINT_DIST_M:
                    reaches_runway = True
                    break
            if reaches_runway:
                filtered.append((rect, axis, role, ref))
        taxi_rects = filtered

    # ── Runway-end stubs for primary parallels ─────────────────────
    # Per user (2026-04-21): primary parallels whose OSM path
    # terminates INSIDE the runway polygon (i.e. the taxi merges
    # onto runway pavement) should emit a wider-than-normal STUB
    # at the transition — the "ramp" where A / F meet the runway
    # short-edge.  Target A stub (688,1562) L=79 W=73 and F stub
    # (2243,-1666) L=93 W=78 both sit ~100 m out from the runway
    # polygon boundary.  Add an extra stub rect at the path vertex
    # just OUTSIDE the runway along the path.
    extra_stubs = _emit_primary_parallel_runway_stubs(
        nodes, ways, to_m, layout.runway_union, pav_union,
        apt_pav_vertices, taxi_rects)
    taxi_rects.extend(extra_stubs)

    # ── Drop overlapping taxi rects ───────────────────────────────
    # Pavement layout invariant (user 2026-04-26): rects MUST NOT
    # overlap each other.  Junctions are constructed as the residue
    # of pavement minus rects/runways/terminals; if rects overlap,
    # the residue is wrong and junctions can't tile the gaps.
    #
    # Sources of rect-rect overlap from upstream rect extraction:
    #   - OSM centerlines too close to each other (a main taxi and a
    #     service road parallel to it; each produces a rect whose
    #     half-width buffer covers the other).
    #   - Apron-emit retries that try to add more rects to simplify
    #     a junction polygon — without an explicit overlap check the
    #     new rect can land on top of an existing one.
    #
    # Drop rule: walk all rect pairs; if their intersection exceeds
    # RECT_OVERLAP_NOISE_M2 (a tiny float-noise threshold so corners
    # that "kiss" but don't actually overlap aren't flagged), drop
    # the rect with the SHORTER axis — the longer rect is more
    # likely the real centerline.  Tie-break by larger area.
    # Iterate until no rect-pair overlap remains.
    RECT_OVERLAP_NOISE_M2 = 1.0
    # Drop only when the overlap is a SIGNIFICANT fraction of the
    # smaller rect's area; a diagonal stub joining a primary parallel
    # legitimately produces a small corner-overlap triangle (e.g. at
    # SPJC, C joining A produces a 24 m² corner overlap on a 10 051 m²
    # C rect — 0.2 %, not a duplicate).  5 % is the empirical
    # boundary: above that, the rects clearly cover the same area;
    # below that, it's a diagonal-stub corner kiss.
    RECT_OVERLAP_FRAC_TOL = 0.05
    if len(taxi_rects) >= 2:
        kept = list(taxi_rects)
        changed = True
        while changed:
            changed = False
            n = len(kept)
            drop_idx: set = set()
            for i in range(n):
                if i in drop_idx:
                    continue
                rect_i, axis_i, role_i, ref_i = kept[i]
                len_i = axis_i.length if axis_i is not None else 0.0
                for j in range(i + 1, n):
                    if j in drop_idx:
                        continue
                    rect_j, axis_j, role_j, ref_j = kept[j]
                    len_j = axis_j.length if axis_j is not None else 0.0
                    try:
                        if not rect_i.intersects(rect_j):
                            continue
                        inter = rect_i.intersection(rect_j)
                        if inter.is_empty:
                            continue
                        if inter.area <= RECT_OVERLAP_NOISE_M2:
                            continue
                        # Diagonal stubs that join a primary parallel
                        # (e.g. SPJC's B/C/E/G connecting to A/F/L/V)
                        # legitimately share a small corner triangle
                        # with their parent.  Drop only when the
                        # overlap is a SIGNIFICANT fraction of the
                        # smaller rect — not just a corner kiss.
                        smaller_area = min(rect_i.area, rect_j.area)
                        if (smaller_area > 0
                                and inter.area
                                / smaller_area
                                <= RECT_OVERLAP_FRAC_TOL):
                            continue
                        # Drop the shorter-axis rect; tie-break by
                        # smaller area.
                        if len_i < len_j or (
                                abs(len_i - len_j) < 0.5
                                and rect_i.area < rect_j.area):
                            drop_idx.add(i)
                            break
                        else:
                            drop_idx.add(j)
                    except _GEOM_EXC:
                        continue
            if drop_idx:
                kept = [k for idx, k in enumerate(kept)
                        if idx not in drop_idx]
                changed = True
        if len(kept) < len(taxi_rects):
            try:
                UI.vprint(1,
                    f"  [pav-builder] {icao}: dropped "
                    f"{len(taxi_rects) - len(kept)} taxi rect(s) "
                    f"that overlapped another rect.")
            except _GEOM_EXC:
                pass
            taxi_rects = kept

    # Per user 2026-04-27 invariant: a junction polygon must NEVER run
    # along a sloping rect's sloping edge.  When a primary_parallel rect
    # sits FULLY INSIDE the airport's pavement union (apt.dat row-110
    # ⊕ DSF ⊕ OSM-synthetic — i.e. the same union that becomes the
    # residue), the surrounding apron junction unavoidably wraps the
    # rect's sloping edges as it traces around the rect-shaped hole.
    # The fix is not to emit the rect at all: let the apron absorb
    # the rect's footprint and slope multi-directionally.  Primary
    # parallels that legitimately cross unpaved area are unaffected
    # (their sloping edges aren't inside pavement).
    #
    # Note: ``pav_union`` is now apt.dat ∪ DSF only (no OSM-synth);
    # that's the right denominator for the absorption check.
    # Including DSF is essential — at SPJC's SE apron, F's long
    # edges are 0 % / 24 % inside row-110 alone but ≈100 % inside
    # apt.dat ∪ DSF.
    taxi_rects = _drop_primary_parallels_embedded_in_pavement(
        taxi_rects, pav_union, runway_polys=runway_polys)

    # ── Detect bridge taxi rects from OSM (user 2026-04-29) ───────
    # Per OSM convention, bridge taxiways carry ``bridge=yes`` (or
    # any non-empty bridge=*) on their parent way.  Build a list of
    # raw OSM bridge-tagged taxiway LineStrings; a rect is a
    # "bridge rect" if its axis lies within a small lateral
    # tolerance of one of those LineStrings.  We use spatial
    # proximity rather than tag-pass-through because
    # ``_extract_osm_taxi_centerlines`` linemerges per-ref, which
    # loses the per-way bridge tag.
    bridge_lines: List[LineString] = []
    for _wid, _nds, _tags in ways:
        if _tags.get("aeroway") != "taxiway":
            continue
        b = _tags.get("bridge", "")
        if not b or b == "no":
            continue
        _pts = []
        for _n in _nds:
            if _n in nodes:
                _lat, _lon = nodes[_n]
                _pts.append(to_m(_lon, _lat))
        if len(_pts) < 2:
            continue
        try:
            _ls = LineString(_pts)
        except _GEOM_EXC:
            continue
        if _ls.is_empty or _ls.length < 5.0:
            continue
        bridge_lines.append(_ls)
    bridge_rect_indices: set = set()
    BRIDGE_AXIS_PROXIMITY_M = 5.0
    if bridge_lines:
        for ri, (_rect, _axis, _role, _ref) in enumerate(taxi_rects):
            if _axis is None or _axis.is_empty:
                continue
            for _bl in bridge_lines:
                # An axis is on a bridge if it lies within
                # BRIDGE_AXIS_PROXIMITY_M of the bridge LineString
                # along most of its length AND has overlap.
                try:
                    if _axis.distance(_bl) > BRIDGE_AXIS_PROXIMITY_M:
                        continue
                    inter = _axis.intersection(
                        _bl.buffer(BRIDGE_AXIS_PROXIMITY_M))
                    if (not inter.is_empty
                            and hasattr(inter, "length")
                            and inter.length
                            >= 0.5 * _axis.length):
                        bridge_rect_indices.add(ri)
                        break
                except _GEOM_EXC:
                    continue

    # ── Hole-aware sloping-edge snap (user 2026-05-04) ────────────
    # When a sloping rect's long edge runs near and parallel to an
    # apt.dat row-110 hole boundary, snap the rect so its sloping
    # edge LIES ON the hole boundary.  Without this, the rect's body
    # sits inside the hole's interior or crosses the hole's perimeter,
    # which then gets absorbed by ``pav_union.difference(rect)`` (the
    # "tunnel" effect): the apron polygon then spans across the hole
    # and ends up sharing boundary with the rect's sloping side.
    # Aligning the rect's sloping edge with the hole boundary makes
    # the surrounding junction wrap around the hole via the rect's
    # CROSS edges instead.
    taxi_rects = _snap_rect_sloping_edges_to_holes(taxi_rects, pav_union)

    # Emit taxi rects (already trimmed to narrow-width portion).
    emitted_taxi_rects: List[Polygon] = []
    for ri, (rect, axis, role, ref) in enumerate(taxi_rects):
        emitted_taxi_rects.append(rect)
        layout.shapes.append(BuiltShape(
            polygon=rect, role=role, ref=ref, source_axis=axis,
            is_bridge=(ri in bridge_rect_indices)))

    # ── Junction emission + pre-Phase-2 geometry finalize ────────
    junction_emit.emit_junctions_and_finalize(
        layout,
        pav_union=pav_union,
        emitted_taxi_rects=emitted_taxi_rects,
        terminal_union=terminal_union,
        taxi_rects=taxi_rects,
        icao=icao)


    # ── Phase-2 elevations + feature emit ────────────────────────
    if compute_elevations:
        finalize.run_phase2(
            layout, icao, xplane_root, apt,
            nodes=nodes, ways=ways, to_m=to_m,
            apron_candidates=apron_candidates,
            tile_dem=tile_dem,
            current_tile_lat=current_tile_lat,
            current_tile_lon=current_tile_lon)

        # Final Rule 2 enforcement (user 2026-05-01).  Triangulation
        # densification and Laplacian-solver vertex insertions can
        # leave a few junction vertices within SLOPING_EDGE_SNAP_M of a
        # sloping rect's sloping edge despite the in-densify guard.  A
        # final post-elevation snap clears these residual cases.
        # Aligned per-vertex altitudes are preserved by index, so the
        # snap is altitude-safe (a 10 m planar move at typical
        # taxi-grade ≤ 1.5 % shifts elevation by ≤ 15 cm — well
        # within the within-shape grade tolerance).
        from .junction_rules import (
            _align_rect_slope_to_axis,
            _enforce_runway_1to1_sharing,
            _snap_to_sloping_edge_corners,
            stitch_pavement_to_flat_runways,
            widen_junctions_to_runway_corners,
        )
        # Slope alignment runs FIRST post-elevation: rects whose
        # slope is purely perpendicular to source_axis become flat
        # (single altitude) so subsequent rules treat them as
        # multi-connection-allowed (per user 2026-05-02).
        _align_rect_slope_to_axis(layout)
        _snap_to_sloping_edge_corners(layout)
        _enforce_runway_1to1_sharing(layout)
        # Rule 1 v6 widening (user 2026-05-02): runs ONLY here,
        # post-elevation, after the runway is segmented.  Inserts
        # outboard runway corners as new junction vertices with
        # matching altitudes.  Pre-elevation widening was disabled
        # because the cascade with Rule 4 + segmentation re-run
        # over-grew junctions past the 4-node cap.
        widen_junctions_to_runway_corners(layout)
        # Stitch pavement to flat runway shapes (user 2026-05-09):
        # for blast pads / flat-interior runway segments, insert a
        # shared vertex on the runway boundary at the projection of
        # every adjacent pavement vertex that sits within edge
        # tolerance.  These new shared vertices become HARD anchors
        # at the runway altitude when the per-surface solver runs,
        # cutting cap-projection distance from the runway corners
        # (often 100s of m apart on long blast pads) down to tens
        # of metres — adjacent junctions / stubs lift toward the
        # runway elevation instead of stalling at terrain.
        stitch_pavement_to_flat_runways(layout)

        # Per user 2026-05-03: per-surface solver runs AS THE LAST
        # STEP of the pipeline, after every junction rule and
        # runway-corner insertion.  ``widen_junctions_to_runway_
        # corners`` inserts vertices with raw runway altitudes that
        # may violate the surrounding junction's per-axis grade
        # rule with respect to neighbouring (terrain-following)
        # vertices; the solver pass cap-projects them into
        # compliance.
        from .elevation import USE_PER_SURFACE_SOLVER
        if USE_PER_SURFACE_SOLVER and layout.anchor is not None:
            from .elevation import _load_airport_dem
            from .elevation_per_surface import solve as per_surface_solve
            dem = tile_dem if tile_dem is not None else _load_airport_dem(
                layout.anchor[0], layout.anchor[1])
            tile_lat = math.floor(layout.anchor[0])
            tile_lon = math.floor(layout.anchor[1])
            per_surface_solve(layout, icao,
                               dem=dem,
                               tile_lat=tile_lat, tile_lon=tile_lon)

            # Grade-based subdivide (``_subdivide_violating_junctions``,
            # threshold 2 %) — catches residual within-shape grade
            # violations the solver alone can't relax.  Iterates up
            # to 4 rounds, then re-runs the solver to integrate the
            # new geometry.
            from .junction_repair import _subdivide_violating_junctions
            n_grade = 0
            for _ in range(4):
                n = _subdivide_violating_junctions(layout)
                if n == 0:
                    break
                n_grade += n
            if n_grade > 0:
                per_surface_solve(layout, icao,
                                   dem=dem,
                                   tile_lat=tile_lat, tile_lon=tile_lon)

        # Stitch pavement to terminal pads (user 2026-05-04): make
        # the two share an identical vertex sequence on every shared
        # edge — pavement vertices near a terminal corner snap to it,
        # vertices in the edge interior get inserted into the
        # terminal polygon.  Eliminates the 4 sub-metre "step"
        # artefacts that survived the densify-skip guard.
        from .junction_rules import (
            stitch_pavement_polygons,
            stitch_pavement_to_terminals,
        )
        stitch_pavement_to_terminals(layout)
        # Adjacent junction polygons whose rings have parallel-but-
        # near-coincident edges should share OSM nids on every shared
        # boundary segment.  Inserts vertices into the other polygon's
        # ring at the projected point with z linearly interpolated
        # along the host edge.  Companion to
        # ``stitch_pavement_to_terminals`` for junction-junction
        # adjacency.  Runs before the corner-snap reconciliation so
        # ``_enforce_shared_vertex_altitudes`` can average shared-
        # bucket altitudes the stitch makes coincident.
        stitch_pavement_polygons(layout)

        # Final cross-shape reconciliation pass (user 2026-05-08):
        # snap each junction vertex whose bucket coincides with a
        # rect / runway / terminal corner to that authoritative
        # shape's altitude tag value, then average junction-to-
        # junction shared-bucket altitudes so neighbouring junctions
        # agree at their seam.  The per-surface solver writes one
        # elevation per bucket but the writeback then averages
        # terminal corners into a single ``altitude`` tag, losing
        # per-corner precision; subsequent geometry passes (overlap
        # clip, sliver merge, stitch) can also drift junction
        # vertices off the solver's value at shared buckets.
        # Without these passes, ``test_pavement_grade``'s cross-
        # shape and step checks find 0.2-7 m gaps at every shared
        # corner where the junction's per-vertex altitude disagrees
        # with the terminal's flat tag, the rect's altitude_high/
        # low, or another junction's value at the same point.
        from .elevation import (
            _snap_junction_altitudes_to_rect_corners,
            _enforce_shared_vertex_altitudes,
        )
        _snap_junction_altitudes_to_rect_corners(layout)
        _enforce_shared_vertex_altitudes(layout)
        # Re-run the rect-corner snap after the junction-pair
        # average, since averaging can pull a shared-with-rect
        # bucket away from the rect's tag value.
        _snap_junction_altitudes_to_rect_corners(layout)

        # The corner-snap above can introduce within-junction grade
        # violations: when one corner of a long junction sits on a
        # runway (snapped to z=77) while another corner is anchored
        # to lower-elevation pavement (z=70), the junction's ring
        # spans 7 m of elevation over ~10 m of distance — 60 %+
        # grade.  Re-run the grade-based subdivide loop to split
        # those polygons into shorter pieces with consistent
        # altitudes, then re-snap so the new sub-polygon corners
        # adopt their respective rect/runway anchor values.
        # SPLP junction-10053 is the canonical case (user 2026-05-08).
        from .junction_repair import _subdivide_violating_junctions
        n_post_snap = 0
        for _ in range(4):
            n = _subdivide_violating_junctions(layout)
            if n == 0:
                break
            n_post_snap += n
        if n_post_snap > 0:
            _snap_junction_altitudes_to_rect_corners(
                layout, interior_proximity_m=3.0)
            _enforce_shared_vertex_altitudes(layout)
            _snap_junction_altitudes_to_rect_corners(
                layout, interior_proximity_m=3.0)

        # Per user 2026-05-10: shapes cannot cross integer lat/lon
        # tile boundaries (X-Plane / Ortho4XP render each 1°x1° tile
        # separately).  Cut a 10 m gap along every tile boundary
        # line that passes through the airport's pavement footprint;
        # Ortho4XP + X-Plane stitch the seam at render time.
        from .tile_cut import cut_layout_at_tile_boundaries
        n_tile_delta = cut_layout_at_tile_boundaries(
            layout,
            current_tile_lat=current_tile_lat,
            current_tile_lon=current_tile_lon,
        )
        if n_tile_delta != 0:
            UI.vprint(1,
                f"  [pav-builder] {icao}: tile-boundary cut "
                f"adjusted shape count by {n_tile_delta:+d}.")

        # Final within-shape grade WARN reflects the absolute
        # final state — junction / apron / terminal Euclidean caps
        # post-final-solver.  Per user 2026-05-03 the WARN was
        # previously firing mid-pipeline with stale numbers.
        from .elevation import _report_within_shape_violations
        _report_within_shape_violations(layout, icao)

    return layout


# ══════════════════════════════════════════════════════════════════
# Phase-2: elevations
# ══════════════════════════════════════════════════════════════════


# ──────────────────────────────────────────────────────────────────
# DEM + CIFP + main elevation pipeline
# (re-exported from O4_Pavement_Junctions)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# DEM + CIFP + main elevation pipeline
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────
from .elevation import (
    _compute_elevations,
)


# ──────────────────────────────────────────────────────────────────
# Airport boundary shape (re-exported from O4_Pavement_Boundary)
# ──────────────────────────────────────────────────────────────────
from .boundary import _emit_airport_boundary_shape


# ──────────────────────────────────────────────────────────────────
# Groundside (curbside / drop-off) pavement
# (re-exported from O4_Pavement_Groundside)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Boundary→DEM bridge polygons (re-exported from O4_Pavement_Boundary)
# ──────────────────────────────────────────────────────────────────
from .boundary import _emit_boundary_dem_bridge


# ──────────────────────────────────────────────────────────────────
# Taxi/road bridges + tunnel portals + depressed-road segments
# (re-exported from O4_Pavement_Bridges; gated by EMIT_BRIDGES_AND_TUNNELS)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Per-shape elevation field + altitude reconciliation
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Elevation finalization (corner buckets, clamp, sliver, overlap)
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# OSM terminal pad extraction (re-exported from O4_Pavement_Terminals)
# ──────────────────────────────────────────────────────────────────
from .terminals import (
    _build_osm_aeroway_footprint,
    _extract_osm_terminals,
    _terminal_groundside_zone,
    _terminal_pad_from_building,
)


# ──────────────────────────────────────────────────────────────────
# Junction-polygon construction
# (re-exported from O4_Pavement_Junctions)
# ──────────────────────────────────────────────────────────────────
from .pavement.junctions import (
    _find_junction_points,
)


# Centerline-slice constants moved to O4_Pavement_Config (MIN_SEGMENT_LEN_M)
# and to O4_Pavement_Centerlines (RDP_SIMPLIFY_TOL_M, SIGNIFICANT_BEND_DEG,
# BEND_CLUSTER_M, GAP_BRIDGE_MAX_M) by slice 3e.
CLOSE_INTERSECTION_M = 200.0  # dead — kept until next cleanup pass
STUB_MAX_LEN_M = 250.0        # dead — kept until next cleanup pass


# ──────────────────────────────────────────────────────────────────
# Same-ref polyline bridging
# (re-exported from O4_Pavement_Centerlines)
# ──────────────────────────────────────────────────────────────────
from .pavement.centerlines import _bridge_same_ref_polylines


# ──────────────────────────────────────────────────────────────────
# Runway-end primary-parallel stub emission
# (re-exported from O4_Pavement_Stubs)
# ──────────────────────────────────────────────────────────────────
from .pavement.stubs import _emit_primary_parallel_runway_stubs


# ──────────────────────────────────────────────────────────────────
# OSM aeroway centerline extraction + splitting
# (re-exported from O4_Pavement_Centerlines)
# ──────────────────────────────────────────────────────────────────
from .pavement.centerlines import (
    _extract_osm_taxi_centerlines,
    _split_centerlines_at_points,
)


# ──────────────────────────────────────────────────────────────────
# Taxi rect construction (re-exported from O4_Pavement_Rects)
# ──────────────────────────────────────────────────────────────────
from .pavement.rects import (
    _build_taxi_rects,
    _classify_role,
    _snap_rect_sloping_edges_to_holes,
)
