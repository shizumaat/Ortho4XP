"""Runway-segment patch emission driven by CIFP threshold elevations.

For each paired runway, samples the DEM along the centerline at
``RUNWAY_SEGMENT_LENGTH`` intervals, anchors CIFP threshold
elevations at each end, then applies grade-limited smoothing
(``MAX_RUNWAY_GRADE``) plus FAA vertical-curve relaxation
(``MAX_RUNWAY_GRADE_CHANGE_PER_M``) to produce a chain of sloped
rectangles that follow the runway's natural contour.

This module owns the **runway-segment auto-patch emit** —
conceptually pavement geometry, not orchestration.  Lives outside
``O4_Auto_Patch`` (the tile-level driver) so the orchestrator
stays a thin caller.

Public API:
    generate_patch_osm
"""
from __future__ import annotations

from math import cos, pi, sqrt

import O4_UI_Utils as UI

from .runway_geometry import (
    DEFAULT_RUNWAY_WIDTH,
    extend_point,
    runway_corners,
)


# ──────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────
DEG_TO_M = 111120.0  # approximate meters per degree of latitude

# FAA AC 150/5300-13B grade limits for Approach Category C-E airports.
MAX_RUNWAY_GRADE = 0.015      # 1.5% max longitudinal grade for runways

# FAA vertical-curve rules: L >= 1000 ft x |delta-G| for runways
# (Design Group III+) means a 1% grade change requires a 305 m
# vertical curve, i.e. ~0.0033% grade change per metre of pavement.
MAX_RUNWAY_GRADE_CHANGE_PER_M = 1.0 / 30000.0

# Iteration cap for the runway-segment grade-relaxation passes
# inside ``generate_patch_osm`` (one for hard-cap, one for vertical
# curve).
GRADE_RELAX_ITERATIONS = 80
RUNWAY_MARGIN = 3.0  # meters added to each side of the runway
DEFAULT_CELL_SIZE = 2.0  # meters between interpolation points (KBNA finding)
DEFAULT_PROFILE = "spline"  # spline profile for natural terrain transitions
# How far beyond the physical runway end to extend as a flat apron.
OVERRUN_EXTENSION = 30.0
RUNWAY_SEGMENT_LENGTH = 100.0  # meters — length of each runway segment


__all__ = [
    "DEFAULT_CELL_SIZE",
    "DEFAULT_PROFILE",
    "DEG_TO_M",
    "GRADE_RELAX_ITERATIONS",
    "MAX_RUNWAY_GRADE",
    "MAX_RUNWAY_GRADE_CHANGE_PER_M",
    "OVERRUN_EXTENSION",
    "RUNWAY_MARGIN",
    "RUNWAY_SEGMENT_LENGTH",
    "generate_patch_osm",
]


def generate_patch_osm(icao, runway_pairs, runway_widths=None, tile=None,
                       apt_runways=None, extra_anchors=None,
                       pav_intersections=None):
    """Generate OSM XML content for segmented runway auto-patches.

    For each paired runway, samples the DEM along the centerline at
    RUNWAY_SEGMENT_LENGTH intervals, anchors CIFP threshold elevations
    at each end, then applies grade-limited smoothing (1.75% max) to
    produce a series of sloped rectangles that follow the runway's
    natural contour — dips, rises, and all.

    Each segment shares its endpoint elevations with its neighbours,
    so the rectangles join seamlessly end-to-end.  Flat overrun
    extensions are added beyond each physical runway end.

    Args:
        icao: Airport ICAO code (e.g. 'SPJC')
        runway_pairs: List from pair_runways()
        runway_widths: Dict of {designator: width_m} from apt.dat, or None.
        tile: Optional Tile object with .dem, .lat, .lon for DEM sampling.
        apt_runways: Optional dict of
            ``{designator: (lat, lon, width_m, displaced_m, blast_m)}``
            parsed from apt.dat row-100 records.  When provided,
            apt.dat is the **sole source of truth for runway
            footprint geometry** — lat/lon, width, displaced
            threshold, and blast-pad / stopway length all come from
            apt.dat for each matching designator.  CIFP is used
            only for threshold elevations.  This keeps the emitted
            runway rectangles pixel-aligned with what X-Plane
            renders, and avoids any geometry-source mismatch.
        extra_anchors: Optional dict of
            ``{(desig_a, desig_b): [(lat, lon, elev_m), ...]}`` —
            additional anchored elevation points along each runway.
            Each anchor is projected onto the runway centerline and
            inserted into the sample list as ANCHORED, so the
            envelope-clamp + grade-cap solver treats it as a hard
            constraint alongside the CIFP threshold anchors.  Used
            to inject cross-runway / taxi-crossing constraints —
            e.g. "where this runway is crossed by a taxi anchored at
            another runway's threshold, this runway must be at the
            other threshold's elevation (within taxi-grade × distance)".
        pav_intersections: Optional dict of
            ``{(desig_a, desig_b): [(lat, lon), ...]}`` — apt.dat
            pavement-polygon vertices that touch the runway boundary
            (computed by ``pipeline.py`` at runway-rect build time).
            Each point is projected onto the centerline and inserted
            as a NON-ANCHORED segment break — the seam corner sits
            there but the elevation comes from DEM/anchor-profile
            interpolation, not a hard constraint.  Aligning seam
            corners with apt.dat boundary intersections lets the
            junction-widening pass reach runway corners via the
            existing single-step chain walk, without needing
            boundary-trace waypoints (per user 2026-05-05).

    Returns:
        str: Complete OSM XML content for the patch file.
    """
    if runway_widths is None:
        runway_widths = {}
    if apt_runways is None:
        apt_runways = {}
    if extra_anchors is None:
        extra_anchors = {}
    node_id = -1
    way_id = -1
    nodes = []  # list of (id, lat, lon)
    ways = []  # list of (id, [node_ids], {tags})
    # Chain of emitted runway segments, captured for downstream
    # consumers that need the authoritative runway elevation at an
    # arbitrary (lat, lon).  Each entry:
    #   (lat_a, lon_a, elev_a, lat_b, lon_b, elev_b, width_m)
    # — exactly the same 7 values add_rect_patch receives when
    # emitting a runway segment, so the elevation model any caller
    # queries from this list matches the patch output byte-for-byte.
    runway_chain = []

    def add_node(lat, lon):
        nonlocal node_id
        nid = node_id
        node_id -= 2
        nodes.append((nid, lat, lon))
        return nid

    def add_way(node_ids, tags):
        nonlocal way_id
        wid = way_id
        way_id -= 2
        ways.append((wid, node_ids, tags))
        return wid

    def add_rect_patch(lat_a, lon_a, elev_a, lat_b, lon_b, elev_b,
                       width):
        """Helper: add a rectangular patch polygon.

        Uses altitude_high/altitude_low when elevations differ,
        otherwise flat altitude tag.  The high-elevation end MUST be
        the first point passed to runway_corners() because
        altitude_high applies to corners 0-1 (the first end).
        """
        if abs(elev_a - elev_b) >= 0.1:
            # Order so the high end is first (corners 0,1)
            if elev_a >= elev_b:
                corners = runway_corners(lat_a, lon_a, lat_b, lon_b, width)
                eh, el = elev_a, elev_b
            else:
                corners = runway_corners(lat_b, lon_b, lat_a, lon_a, width)
                eh, el = elev_b, elev_a
            if corners is None:
                return
            n0 = add_node(corners[0][0], corners[0][1])
            n1 = add_node(corners[1][0], corners[1][1])
            n2 = add_node(corners[2][0], corners[2][1])
            n3 = add_node(corners[3][0], corners[3][1])
            tags = {
                "altitude_high": "{:.1f}".format(eh),
                "altitude_low": "{:.1f}".format(el),
                "cell_size": str(int(DEFAULT_CELL_SIZE)),
                "profile": DEFAULT_PROFILE,
            }
        else:
            corners = runway_corners(lat_a, lon_a, lat_b, lon_b, width)
            if corners is None:
                return
            n0 = add_node(corners[0][0], corners[0][1])
            n1 = add_node(corners[1][0], corners[1][1])
            n2 = add_node(corners[2][0], corners[2][1])
            n3 = add_node(corners[3][0], corners[3][1])
            avg = round((elev_a + elev_b) / 2.0, 1)
            tags = {"altitude": "{:.1f}".format(avg)}
        add_way([n0, n1, n2, n3, n0], tags)

    def add_flat_multi_rect(samples_ll, elev, width):
        """Helper: add a multi-node FLAT runway polygon.

        ``samples_ll`` is the centerline sample list ordered start to
        end; the polygon is built by extending each sample
        perpendicular by ±width/2.  Used when consecutive sloped/flat
        segments collapse into a single flat run — intermediate
        samples are kept only at pav_intersection positions where
        adjacent junctions need to snap.

        Per user 2026-05-09: flat shapes use single ``altitude=`` tag
        and may carry an arbitrary number of corners.
        """
        if len(samples_ll) < 2:
            return
        lat_a, lon_a = samples_ll[0]
        lat_b, lon_b = samples_ll[-1]
        bcorners = runway_corners(lat_a, lon_a, lat_b, lon_b, width)
        if bcorners is None:
            return
        # ``runway_corners`` returns [a-left, b-left, b-right, a-right]
        # — perpendicular offset is (corner[0] - sample[0]).
        perp_dlat = bcorners[0][0] - lat_a
        perp_dlon = bcorners[0][1] - lon_a
        # Build ring: A→B on left side, then B→A on right side.
        ring: List[Tuple[float, float]] = []
        for s_lat, s_lon in samples_ll:
            ring.append((s_lat + perp_dlat, s_lon + perp_dlon))
        for s_lat, s_lon in reversed(samples_ll):
            ring.append((s_lat - perp_dlat, s_lon - perp_dlon))
        node_ids = [add_node(la, lo) for la, lo in ring]
        node_ids.append(node_ids[0])
        tags = {"altitude": "{:.1f}".format(round(float(elev), 1))}
        add_way(node_ids, tags)

    def _sample_dem(lat, lon):
        """Sample DEM elevation at a lat/lon, returning 0 on failure."""
        if tile is None or not hasattr(tile, "dem") or tile.dem is None:
            return None
        try:
            return tile.dem.alt((lon - tile.lon, lat - tile.lat))
        except (IndexError, ValueError, TypeError, ZeroDivisionError):
            return None

    # ── Auto cross-runway anchor pre-pass ──────────────────────
    # Per user 2026-04-28: project every paired runway's threshold
    # onto every OTHER paired runway's centerline.  When the
    # projection falls inside that other runway's length AND the
    # perpendicular distance is plausibly walkable by a taxi
    # (≤ ``MAX_CROSS_RUNWAY_LATERAL_M``), the projection becomes
    # an additional anchor on the receiver runway.
    #
    # Per user 2026-04-28 (refined): the anchor's elevation is
    # NOT pinned to the source threshold's elevation.  Instead it
    # follows DEM at the projection point, but is CLAMPED so a
    # connecting taxi at MAX_TAXI_GRADE (1.5 %) over the
    # perpendicular distance can still reach the source threshold
    # — i.e. the anchor lies in
    # ``[src_elev − perp × 0.015, src_elev + perp × 0.015]``.  If
    # DEM is in band, use DEM; if outside, clamp to the nearest
    # band edge.  This keeps the receiver runway as close to its
    # natural terrain as possible while still guaranteeing the
    # connecting taxi can be built.
    MAX_CROSS_RUNWAY_LATERAL_M = 300.0
    MAX_TAXI_GRADE_FOR_CROSS = 0.015  # 1.5 % FAA cap for taxiways
    auto_extra_anchors: dict = {}
    paired_list = [(da, dat_a, db, dat_b)
                   for da, dat_a, db, dat_b in runway_pairs
                   if db is not None and dat_b is not None]
    for ti, (da_t, dat_a_t, db_t, dat_b_t) in enumerate(paired_list):
        for src_desig, src_data in (
                (da_t, dat_a_t), (db_t, dat_b_t)):
            for ri, (da_r, dat_a_r, db_r, dat_b_r) in enumerate(
                    paired_list):
                if ri == ti:
                    continue
                mid_lat = 0.5 * (dat_a_r["lat"] + dat_b_r["lat"])
                cl_v = cos(mid_lat * pi / 180.0)
                if cl_v < 1e-6:
                    cl_v = 1e-6
                rdx = (dat_b_r["lon"] - dat_a_r["lon"]) * cl_v * DEG_TO_M
                rdy = (dat_b_r["lat"] - dat_a_r["lat"]) * DEG_TO_M
                rL2 = rdx * rdx + rdy * rdy
                if rL2 < 1.0:
                    continue
                vx = (src_data["lon"] - dat_a_r["lon"]) * cl_v * DEG_TO_M
                vy = (src_data["lat"] - dat_a_r["lat"]) * DEG_TO_M
                t = (vx * rdx + vy * rdy) / rL2
                if t <= 0.05 or t >= 0.95:
                    continue
                proj_x = t * rdx
                proj_y = t * rdy
                perp = sqrt((vx - proj_x) ** 2 + (vy - proj_y) ** 2)
                if perp > MAX_CROSS_RUNWAY_LATERAL_M:
                    continue
                p_lat = (dat_a_r["lat"]
                         + t * (dat_b_r["lat"] - dat_a_r["lat"]))
                p_lon = (dat_a_r["lon"]
                         + t * (dat_b_r["lon"] - dat_a_r["lon"]))
                # DEM-preferred elevation, clamped to taxi-grade
                # band from the source threshold.
                src_elev = src_data["elevation_m"]
                band = perp * MAX_TAXI_GRADE_FOR_CROSS
                lo_band = src_elev - band
                hi_band = src_elev + band
                dem_e = _sample_dem(p_lat, p_lon)
                if dem_e is None:
                    # No DEM available → midpoint of the band as a
                    # safe seed (equivalent to "as close to source
                    # threshold as the grade lets us").
                    anchor_e = src_elev
                elif dem_e < lo_band:
                    anchor_e = lo_band
                elif dem_e > hi_band:
                    anchor_e = hi_band
                else:
                    anchor_e = dem_e
                key = (da_r, db_r)
                auto_extra_anchors.setdefault(key, []).append(
                    (p_lat, p_lon, anchor_e))
    # Merge user-supplied extra_anchors on top of auto-detected
    # ones — user values take precedence (replace auto if same
    # exact lat/lon, else append).
    for k, v in (extra_anchors or {}).items():
        auto_extra_anchors.setdefault(k, []).extend(v)
    extra_anchors = auto_extra_anchors

    for desig_a, data_a, desig_b, data_b in runway_pairs:
        if desig_b is not None and data_b is not None:
            # ── Paired runway ────────────────────────────────────────────
            # apt.dat is the sole source of truth for runway
            # footprint geometry (lat/lon, width, displaced
            # thresholds, blast pads).  CIFP contributes ONLY the
            # threshold elevations used to seed the per-segment
            # elevation profile.
            elev_a = data_a["elevation_m"]
            elev_b = data_b["elevation_m"]

            apt_a = apt_runways.get(desig_a)
            apt_b = apt_runways.get(desig_b)
            have_apt_geom = apt_a is not None and apt_b is not None

            if have_apt_geom:
                lat_a, lon_a = apt_a[0], apt_a[1]
                lat_b, lon_b = apt_b[0], apt_b[1]
                rwy_width = apt_a[2]
                displaced_a = apt_a[3]
                displaced_b = apt_b[3]
                blast_a = apt_a[4]
                blast_b = apt_b[4]
            else:
                # Legacy fallback: no apt.dat geometry available,
                # use CIFP lat/lon + displaced, runway_widths dict
                # for width.  Blast pads are not in CIFP so fall
                # back to OVERRUN_EXTENSION.
                lat_a, lon_a = data_a["lat"], data_a["lon"]
                lat_b, lon_b = data_b["lat"], data_b["lon"]
                displaced_a = data_a["displaced_m"]
                displaced_b = data_b["displaced_m"]
                blast_a = OVERRUN_EXTENSION
                blast_b = OVERRUN_EXTENSION
                rwy_width = runway_widths.get(
                    desig_a,
                    runway_widths.get(desig_b, DEFAULT_RUNWAY_WIDTH),
                )
            patch_width = rwy_width + 2 * RUNWAY_MARGIN

            # cos(lat) for meter conversions
            mid_lat = (lat_a + lat_b) / 2.0
            cos_lat_v = cos(mid_lat * pi / 180.0)
            if cos_lat_v < 1e-6:
                cos_lat_v = 1e-6

            # ── Physical runway ends ─────────────────────────────────────
            # apt.dat row-100 lat/lon ARE the physical ends of the
            # runway surface (excluding blast pads).  CIFP thresholds
            # are at (lat_a/b) + displaced_a/b inward.
            #
            # Legacy CIFP fallback: lat_a/b in CIFP are at the
            # displaced threshold; extend outward by displaced_m to
            # approximate the physical end.
            if have_apt_geom:
                phys_end_a = (lat_a, lon_a)
                phys_end_b = (lat_b, lon_b)
            else:
                if displaced_a > 0:
                    phys_end_a = extend_point(
                        lat_b, lon_b, lat_a, lon_a, displaced_a)
                else:
                    phys_end_a = (lat_a, lon_a)
                if displaced_b > 0:
                    phys_end_b = extend_point(
                        lat_a, lon_a, lat_b, lon_b, displaced_b)
                else:
                    phys_end_b = (lat_b, lon_b)

            # Per user 2026-05-09: treat blast-pad / overrun
            # extensions as displaced-threshold continuations.
            # Extending phys_end_a/b outward by ``blast_a/b`` and
            # absorbing those distances into ``displaced_a/b`` makes
            # the existing chain cover the full physical extent
            # (blast pad → runway proper → blast pad) in one pass.
            # CIFP threshold elevations stay anchored at the
            # displaced-threshold positions (now interior to the
            # extended chain), and the blast-pad area gets DEM-
            # sampled + grade-limited just like the runway interior.
            # No separate flat-rect emit is needed for blast pads.
            if blast_a > 0.1:
                ext_a = extend_point(
                    phys_end_b[0], phys_end_b[1],
                    phys_end_a[0], phys_end_a[1],
                    blast_a,
                )
                phys_end_a = ext_a
                displaced_a += blast_a
                blast_a = 0.0
            if blast_b > 0.1:
                ext_b = extend_point(
                    phys_end_a[0], phys_end_a[1],
                    phys_end_b[0], phys_end_b[1],
                    blast_b,
                )
                phys_end_b = ext_b
                displaced_b += blast_b
                blast_b = 0.0

            # Full physical runway length (phys_end to phys_end).
            dx_phys = (phys_end_b[1] - phys_end_a[1]) * cos_lat_v * DEG_TO_M
            dy_phys = (phys_end_b[0] - phys_end_a[0]) * DEG_TO_M
            phys_dist = sqrt(dx_phys ** 2 + dy_phys ** 2)
            if phys_dist < 1.0:
                continue

            # Threshold-to-threshold distance (between the two
            # DISPLACED thresholds, where CIFP elevations are
            # anchored).  This is what "grade" is measured over.
            thresh_dist = phys_dist - displaced_a - displaced_b
            if thresh_dist < 1.0:
                thresh_dist = phys_dist  # degenerate, both disp=0
            grade = (elev_b - elev_a) / thresh_dist

            # Elevation at each physical end: shift the CIFP
            # threshold elevation by the grade over the displaced
            # distance (0 when disp=0 → elev_a/b unchanged).
            elev_phys_a = elev_a - grade * displaced_a
            elev_phys_b = elev_b + grade * displaced_b

            # ── Build segment sample points along centerline ─────────────
            # Always include: physical end A, threshold A, threshold B,
            # physical end B.  Add DEM sample points in between.
            n_segs = max(1, int(phys_dist / RUNWAY_SEGMENT_LENGTH))
            # Sample points as fraction of phys_dist (0 = end A, 1 = end B)
            fractions = [float(i) / n_segs for i in range(n_segs + 1)]

            # Ensure thresholds are in the list.  ``anchored_t``
            # tracks t-values that MUST stay (physical ends and
            # CIFP threshold positions) so the pav_intersection
            # dedup below can distinguish them from uniform seams.
            anchored_t: List[float] = [0.0, 1.0]
            if phys_dist > 0:
                t_a = (displaced_a / phys_dist) if displaced_a > 0 else 0.0
                t_b = 1.0 - (displaced_b / phys_dist) if displaced_b > 0 else 1.0
                for t in [t_a, t_b]:
                    if 0.0 < t < 1.0:
                        # Insert if not already near an existing fraction
                        if not any(abs(t - f) < 0.01 for f in fractions):
                            fractions.append(t)
                        if not any(abs(t - a) < 1e-6 for a in anchored_t):
                            anchored_t.append(t)
                fractions.sort()

            # ``pav_int_t_vals`` accumulates the t-values of every
            # pav_intersection that survives dedup into ``fractions``.
            # Used downstream to identify which sample positions are
            # legitimate junction-snap points so they get retained as
            # intermediate corners when consecutive flat segments are
            # consolidated into a single multi-node flat polygon
            # (user 2026-05-09).
            pav_int_t_vals: List[float] = []
            # Per user 2026-05-05: inject pav_intersection breakpoints
            # so segment seam corners align with apt.dat-pavement
            # boundary points where the apron / taxiway meets the
            # runway.  These are NOT elevation anchors — just
            # geometric segment seams.
            #
            # Dedup behaviour (user 2026-05-05 followup): when a
            # pav_intersection lands close to a NON-anchored
            # (uniform 100 m) seam, REPLACE the uniform seam with
            # the apt.dat position.  Otherwise the runway ends up
            # with 3-12 m sliver segments between a uniform seam
            # and an apt.dat intersection, and junctions only share
            # the apt.dat one — the uniform seam is dead weight.
            # Threshold 12 m of centerline distance: large enough
            # to consolidate the typical sliver, small enough to
            # leave genuinely independent intersections alone.
            # Anchored fractions (thresholds, physical ends) are
            # never replaced; pav_intersections within
            # ``anchor_dedup_m`` of an anchor are dropped.
            if pav_intersections and phys_dist > 0:
                pav_pts = []
                for pkey in (
                        (desig_a, desig_b),
                        (desig_b, desig_a),
                        ("RW" + desig_a.lstrip("RW"),
                         "RW" + desig_b.lstrip("RW")),
                        ("RW" + desig_b.lstrip("RW"),
                         "RW" + desig_a.lstrip("RW"))):
                    if pkey in pav_intersections:
                        pav_pts = pav_intersections[pkey]
                        break
                rL2 = dx_phys * dx_phys + dy_phys * dy_phys
                merge_t = 12.0 / phys_dist
                anchor_dedup_t = 2.0 / phys_dist
                for pp_lat, pp_lon in pav_pts:
                    px = (pp_lon - phys_end_a[1]) * cos_lat_v * DEG_TO_M
                    py = (pp_lat - phys_end_a[0]) * DEG_TO_M
                    if rL2 <= 0:
                        break
                    pt = (px * dx_phys + py * dy_phys) / rL2
                    if pt <= 0.001 or pt >= 0.999:
                        continue
                    # If close to an anchor (threshold/end), drop
                    # the apt.dat intersection (anchor wins).
                    if any(abs(pt - a) < anchor_dedup_t for a in anchored_t):
                        continue
                    # Otherwise, replace the closest non-anchored
                    # fraction within merge_t (= uniform-seam
                    # snap), or append if none in range.
                    closest_idx = None
                    closest_d = merge_t
                    for i, f in enumerate(fractions):
                        if any(abs(f - a) < 1e-6 for a in anchored_t):
                            continue
                        d = abs(pt - f)
                        if d < closest_d:
                            closest_d = d
                            closest_idx = i
                    if closest_idx is not None:
                        fractions[closest_idx] = pt
                    else:
                        fractions.append(pt)
                    pav_int_t_vals.append(pt)
                fractions.sort()

            # NOTE: tile-boundary cuts intentionally happen at the
            # end of the pipeline in ``tile_cut.py``, NOT here.
            # Building the chain with awareness of tile-gap edges
            # interferes with junction widening (the widen pass
            # promotes the new "gap edge" runway corners into
            # adjacent junctions, which then misalign with the cut
            # the post-process applies later).  Per user 2026-05-12:
            # keep the chain construction agnostic; cut everything
            # uniformly at the end.

            # For each sample point, compute lat/lon and seed elevation
            sample_pts = []  # [(lat, lon, seeded_elev, is_anchored), ...]
            for frac in fractions:
                s_lat = phys_end_a[0] + frac * (phys_end_b[0] - phys_end_a[0])
                s_lon = phys_end_a[1] + frac * (phys_end_b[1] - phys_end_a[1])
                dist_from_a = frac * phys_dist

                # Is this a CIFP-anchored point?
                is_anchor = False
                if displaced_a > 0:
                    if abs(dist_from_a - displaced_a) < 1.0:
                        sample_pts.append((s_lat, s_lon, elev_a, True))
                        continue
                else:
                    if frac < 0.001:
                        sample_pts.append((s_lat, s_lon, elev_phys_a, True))
                        continue

                if displaced_b > 0:
                    if abs(dist_from_a - (phys_dist - displaced_b)) < 1.0:
                        sample_pts.append((s_lat, s_lon, elev_b, True))
                        continue
                else:
                    if frac > 0.999:
                        sample_pts.append((s_lat, s_lon, elev_phys_b, True))
                        continue

                # Physical ends are anchored
                if frac < 0.001:
                    sample_pts.append((s_lat, s_lon, elev_phys_a, True))
                    continue
                if frac > 0.999:
                    sample_pts.append((s_lat, s_lon, elev_phys_b, True))
                    continue

                # Interior point: use DEM if available, else interpolate
                dem_val = _sample_dem(s_lat, s_lon)
                if dem_val is not None:
                    sample_pts.append((s_lat, s_lon, dem_val, False))
                else:
                    # Linear interpolation between thresholds
                    interp = elev_phys_a + frac * (elev_phys_b - elev_phys_a)
                    sample_pts.append((s_lat, s_lon, interp, False))

            # Per user 2026-04-28: cross-runway anchor injection.
            # ``extra_anchors`` maps a (desig_a, desig_b) key (with
            # both orderings) to a list of (lat, lon, elev) anchor
            # points along this runway that came from intersections
            # with other runways or taxiways at known elevations.
            # Each one becomes an ANCHORED sample so the envelope-
            # clamp + grade-cap solver treats it as a hard
            # constraint along with the CIFP threshold anchors.
            ra_key = None
            extra = []
            if extra_anchors:
                for k in (
                        (desig_a, desig_b),
                        (desig_b, desig_a),
                        ("RW" + desig_a.lstrip("RW"),
                         "RW" + desig_b.lstrip("RW")),
                        ("RW" + desig_b.lstrip("RW"),
                         "RW" + desig_a.lstrip("RW"))):
                    if k in extra_anchors:
                        extra = extra_anchors[k]
                        ra_key = k
                        break
            for a_lat, a_lon, a_elev in extra:
                # Project anchor lat/lon onto the runway centerline
                # parameter (frac in [0, 1] from phys_end_a to
                # phys_end_b).  Use simple lat/lon-as-Cartesian since
                # the runway is short.
                ax = (a_lon - phys_end_a[1]) * cos_lat_v * DEG_TO_M
                ay = (a_lat - phys_end_a[0]) * DEG_TO_M
                # rwy direction in meters
                rdx = dx_phys
                rdy = dy_phys
                rL2 = rdx * rdx + rdy * rdy
                if rL2 <= 0:
                    continue
                t = (ax * rdx + ay * rdy) / rL2
                if t <= 0.001 or t >= 0.999:
                    continue
                # Interpolate the lat/lon onto the centerline at t.
                p_lat = phys_end_a[0] + t * (
                    phys_end_b[0] - phys_end_a[0])
                p_lon = phys_end_a[1] + t * (
                    phys_end_b[1] - phys_end_a[1])
                # Insert into sample_pts in fraction-order.
                inserted = False
                for j in range(len(sample_pts)):
                    if sample_pts[j][:2] == (p_lat, p_lon):
                        # Already a sample at this exact point — upgrade
                        # to anchored with the constraint elevation.
                        sample_pts[j] = (p_lat, p_lon, a_elev, True)
                        inserted = True
                        break
                if inserted:
                    continue
                # Find sorted insertion position.
                # Recompute frac for each existing sample (they were
                # appended in fractions order).
                for j in range(len(sample_pts)):
                    s_la, s_lo = sample_pts[j][0], sample_pts[j][1]
                    s_ax = (s_lo - phys_end_a[1]) * cos_lat_v * DEG_TO_M
                    s_ay = (s_la - phys_end_a[0]) * DEG_TO_M
                    s_t = (s_ax * rdx + s_ay * rdy) / rL2
                    if s_t > t:
                        sample_pts.insert(
                            j, (p_lat, p_lon, a_elev, True))
                        # Also extend fractions list to keep ordering
                        # state consistent.
                        fractions.insert(j, t)
                        inserted = True
                        break
                if not inserted:
                    sample_pts.append((p_lat, p_lon, a_elev, True))
                    fractions.append(t)

            # ── Anchor-profile baseline + DEM blend ───────────────────
            # Per user 2026-04-28: build the runway profile by first
            # collecting ALL anchors (CIFP thresholds + physical
            # ends + cross-runway projections), constructing a smooth
            # FAA-compliant baseline that passes through them, then
            # blending in DEM up to a bounded deviation.  This
            # replaces the previous "DEM-seed → envelope clamp" path
            # that produced V-shaped kinks at hard anchors when DEM
            # was higher than the anchor (the rate-of-change solver
            # couldn't smooth the kink because anchored samples are
            # immutable).
            #
            # Order of operations:
            #   1. Collect all anchored samples as (frac, elev)
            #      pairs, sorted by frac.
            #   2. Validate anchor feasibility — warn if any
            #      adjacent-anchor pair grade exceeds
            #      MAX_RUNWAY_GRADE.
            #   3. Build the linear-interpolation anchor profile:
            #      profile(frac) returns the smooth baseline at any
            #      fraction along the runway.
            #   4. For each non-anchor sample, set
            #         elev[i] = clamp(dem_e, base_e ± DEM_BAND)
            #      where base_e = profile(fractions[i]) and DEM_BAND
            #      is the maximum DEM deviation we allow before
            #      tightening to the baseline.
            #   5. Existing _pass_hard_cap + _pass_rate_of_change run
            #      below as a final cleanup; with the baseline-
            #      driven seed they have very little work to do.
            elevs = [s[2] for s in sample_pts]
            anchored = [s[3] for s in sample_pts]
            n_samples = len(elevs)

            # cumulative distance along the centerline (used by
            # validity check + spans-to-each-anchor logging).
            cum_dist = [0.0]
            for i in range(1, n_samples):
                cum_dist.append(
                    cum_dist[-1]
                    + abs(fractions[i] - fractions[i - 1]) * phys_dist)

            # 1. Collect anchors.
            profile_anchors: list = [
                (fractions[i], elevs[i])
                for i in range(n_samples)
                if anchored[i]]
            profile_anchors.sort()

            # 2. Validate feasibility.
            for k in range(len(profile_anchors) - 1):
                f0, e0 = profile_anchors[k]
                f1, e1 = profile_anchors[k + 1]
                seg_d = abs(f1 - f0) * phys_dist
                if seg_d <= 0.5:
                    continue
                g = abs(e1 - e0) / seg_d
                if g > MAX_RUNWAY_GRADE + 1e-6:
                    UI.vprint(1,
                        f"  [auto-patch] {icao} runway "
                        f"{desig_a}/{desig_b}: anchor pair grade "
                        f"{g * 100:.2f}% > "
                        f"{MAX_RUNWAY_GRADE * 100:.1f}% between "
                        f"fractions {f0:.3f} and {f1:.3f} "
                        f"({seg_d:.0f} m apart, ΔE={e1 - e0:+.2f} m) "
                        f"— profile will be infeasible at FAA "
                        f"runway max grade.")

            # 3. Linear-interpolation anchor profile.
            def _anchor_profile(frac: float):
                if not profile_anchors:
                    return None
                if frac <= profile_anchors[0][0]:
                    return profile_anchors[0][1]
                if frac >= profile_anchors[-1][0]:
                    return profile_anchors[-1][1]
                for k in range(len(profile_anchors) - 1):
                    f0, e0 = profile_anchors[k]
                    f1, e1 = profile_anchors[k + 1]
                    if f0 <= frac <= f1:
                        if f1 - f0 < 1e-9:
                            return e0
                        t = (frac - f0) / (f1 - f0)
                        return e0 + t * (e1 - e0)
                return profile_anchors[-1][1]

            # 4. DEM blend with band size tied to FAA absorption capacity.
            # Per user 2026-04-28: "max DEM following" — let DEM
            # influence the profile up to whatever the FAA rate-of-
            # change rule can absorb in a parabolic vertical curve
            # connecting back to the anchor profile.
            #
            # The maximum deviation of a parabolic vertical curve
            # from its tangent line at distance d from the curve's
            # PVI (point of vertical intersection) is
            #     max_dev = 0.5 × K × d²
            # where K is the rate-of-change constant
            # (MAX_RUNWAY_GRADE_CHANGE_PER_M = 1/30,000 for runways).
            # Treating each anchor as a PVI, the band at sample i
            # is the tightest such limit imposed by ANY anchor.
            # Near anchors the band is small (forces the profile to
            # stay near the linear baseline); far from anchors it
            # widens enough that DEM character can come through.
            DEM_BAND_M_MAX = 5.0   # absolute cap (sanity ceiling).
            for i in range(n_samples):
                if anchored[i]:
                    continue
                base_e = _anchor_profile(fractions[i])
                if base_e is None:
                    continue
                # Distance to nearest anchor along centerline.
                nearest_d = float('inf')
                for j in range(n_samples):
                    if not anchored[j]:
                        continue
                    d = abs(cum_dist[i] - cum_dist[j])
                    if d < nearest_d:
                        nearest_d = d
                # FAA vertical-curve absorption capacity.
                fa_cap = (0.5 * MAX_RUNWAY_GRADE_CHANGE_PER_M
                           * nearest_d * nearest_d)
                band = min(DEM_BAND_M_MAX, fa_cap)
                dem_e = elevs[i]
                if dem_e > base_e + band:
                    elevs[i] = base_e + band
                elif dem_e < base_e - band:
                    elevs[i] = base_e - band
                else:
                    elevs[i] = dem_e

            # ── FAA-absorption envelope pre-clamp ────────────────
            # Per user 2026-04-28: enforce the parabolic vertical-
            # curve envelope explicitly.  At distance d from any
            # anchor, the maximum elevation deviation that's
            # achievable while respecting MAX_GC × distance grade-
            # change-rate is:
            #   max_dev = 0.5 × MAX_GC × d²       (for d ≤ L_VC)
            #   max_dev = 0.5 × MAX_GC × L_VC²
            #             + MAX_GRADE × (d − L_VC)  (for d > L_VC)
            # where L_VC = MAX_GRADE / MAX_GC (the curve length to
            # reach max grade).
            #
            # Threshold anchors with blast pads ALSO need to
            # satisfy this envelope assuming g_in=0 from the flat
            # blast pad.  Interior (cross-runway) anchors don't
            # have a defined g_in, so the envelope around them
            # uses just the linear MAX_GRADE × d bound.
            #
            # This pre-clamp guarantees the profile is FAA-
            # compliant by construction.  The downstream
            # _pass_hard_cap and _pass_rate_of_change converge in
            # 1-2 iterations (mostly to round 0.001 m residuals).
            L_VC = MAX_RUNWAY_GRADE / MAX_RUNWAY_GRADE_CHANGE_PER_M

            def _faa_max_dev(d: float) -> float:
                if d <= L_VC:
                    return 0.5 * MAX_RUNWAY_GRADE_CHANGE_PER_M * d * d
                return (0.5 * MAX_RUNWAY_GRADE_CHANGE_PER_M
                        * L_VC * L_VC
                        + MAX_RUNWAY_GRADE * (d - L_VC))

            for i in range(n_samples):
                if anchored[i]:
                    continue
                lo = float('-inf')
                hi = float('inf')
                for j in range(n_samples):
                    if not anchored[j]:
                        continue
                    d_ij = abs(cum_dist[i] - cum_dist[j])
                    # Boundary anchors with blast pads use the
                    # FAA absorption envelope (the blast pad's
                    # g_in=0 is part of the constraint).  Interior
                    # anchors use only the linear grade envelope.
                    is_boundary_with_blast = (
                        (j == 0 and blast_a > 0.1)
                        or (j == n_samples - 1 and blast_b > 0.1))
                    if is_boundary_with_blast:
                        cap = _faa_max_dev(d_ij)
                    else:
                        cap = MAX_RUNWAY_GRADE * d_ij
                    lo = max(lo, elevs[j] - cap)
                    hi = min(hi, elevs[j] + cap)
                if lo <= hi:
                    if elevs[i] > hi:
                        elevs[i] = hi
                    elif elevs[i] < lo:
                        elevs[i] = lo
                else:
                    # Infeasible: anchors are inconsistent with
                    # the FAA envelope.  Use baseline as fallback.
                    base_e = _anchor_profile(fractions[i])
                    if base_e is not None:
                        elevs[i] = base_e

            def _pass_hard_cap():
                for _it in range(GRADE_RELAX_ITERATIONS):
                    changed = False
                    for idx in range(len(elevs)):
                        if anchored[idx]:
                            continue
                        lo = float("-inf")
                        hi = float("inf")
                        for nidx in (idx - 1, idx + 1):
                            if nidx < 0 or nidx >= len(elevs):
                                continue
                            seg_dist = (
                                abs(fractions[nidx] - fractions[idx])
                                * phys_dist)
                            if seg_dist < 0.1:
                                continue
                            max_rise = seg_dist * MAX_RUNWAY_GRADE
                            lo = max(lo, elevs[nidx] - max_rise)
                            hi = min(hi, elevs[nidx] + max_rise)
                        if lo == float("-inf") and hi == float("inf"):
                            continue
                        if lo > hi:
                            # Infeasible — neighbours diverge more
                            # than the grade cap allows.  Snap to
                            # the midpoint so the constraint
                            # propagates.
                            new_e = (lo + hi) / 2.0
                        else:
                            new_e = min(max(elevs[idx], lo), hi)
                        if abs(new_e - elevs[idx]) > 0.001:
                            elevs[idx] = new_e
                            changed = True
                    if not changed:
                        return

            # Joint solver: alternate hard-cap and rate-of-change
            # passes until neither one changes anything.  Running
            # each pass only once can leave small mutual residuals
            # (the cap pass can push samples outside the rate-of-
            # change envelope, and vice versa); with the envelope
            # pre-clamp above this usually converges in 2–3 outer
            # iterations.
            _pass_hard_cap()

            # FAA vertical-curve rate-of-change pass for RUNWAYS.
            # The runway rule is L ≥ 305 m × |ΔG| (vs 30.5 m for
            # taxiways), which means grade can only change by
            # ≈ 0.0033 % per metre of pavement.  Over a 100 m
            # segment adjacent grades can therefore differ by at
            # most 0.33 %, so going from 0 % to the 1.5 % cap takes
            # at least ~450 m of runway.  This is what produces the
            # "long gentle slopes" characteristic of real runways.
            #
            # For each interior sample i (1..n-2):
            #   g_left  = (e[i]   - e[i-1]) / L_left
            #   g_right = (e[i+1] - e[i])   / L_right
            #   |g_right - g_left| must be
            #     ≤ MAX_RUNWAY_GRADE_CHANGE_PER_M × ((L_left+L_right)/2)
            def _seg_len(i):
                return abs(fractions[i + 1] - fractions[i]) * phys_dist

            def _pass_rate_of_change():
                """FAA rate-of-grade-change pass.

                Per user 2026-04-28 (refined): treats blast pads as
                virtual anchored samples (grade=0 contributing
                segments) at the start and end of the sample list.
                With this, the constraint at the threshold/blast-pad
                interface propagates inward over MULTIPLE samples
                instead of being one-sample-deep — so the runway
                forms a proper FAA vertical curve coming out of the
                threshold, not a sharp grade kink.

                Also enforces ΔG at INTERIOR anchored samples (e.g.
                cross-runway projection anchors) by adjusting
                the anchor's neighbours when the anchor's incoming
                vs outgoing grade differ by more than the FAA limit.
                """
                # Build extended arrays with virtual blast pad
                # samples at frac < 0 and frac > 1.  The virtual
                # samples are anchored at the same elevation as
                # the threshold, so the segment between them and
                # the threshold has grade 0 — exactly modelling a
                # flat blast pad.
                elevs_ext = list(elevs)
                anchored_ext = list(anchored)
                # Segment lengths: one entry per gap between
                # consecutive samples in elevs_ext.
                seg_lens: list = [_seg_len(k)
                                   for k in range(len(elevs) - 1)]
                start_offset = 0
                if blast_a > 0.1 and anchored[0]:
                    elevs_ext.insert(0, elevs[0])
                    anchored_ext.insert(0, True)
                    seg_lens.insert(0, blast_a)
                    start_offset = 1
                end_offset = 0
                if blast_b > 0.1 and anchored[-1]:
                    elevs_ext.append(elevs[-1])
                    anchored_ext.append(True)
                    seg_lens.append(blast_b)
                    end_offset = 1

                def _eseg(i):
                    return seg_lens[i] if 0 <= i < len(seg_lens) else 0.0

                any_change = False
                for _it in range(GRADE_RELAX_ITERATIONS):
                    changed = False
                    # Walk every sample (including anchored interior
                    # ones — they don't move themselves but their
                    # ΔG is enforced by adjusting neighbours).
                    for i in range(1, len(elevs_ext) - 1):
                        ll = _eseg(i - 1)
                        lr = _eseg(i)
                        if ll < 0.1 or lr < 0.1:
                            continue
                        g_left = (
                            (elevs_ext[i] - elevs_ext[i - 1]) / ll)
                        g_right = (
                            (elevs_ext[i + 1] - elevs_ext[i]) / lr)
                        max_dg = (MAX_RUNWAY_GRADE_CHANGE_PER_M
                                   * ((ll + lr) / 2.0))
                        dg = g_right - g_left
                        if abs(dg) <= max_dg:
                            continue
                        target_dg = (max_dg if dg > 0 else -max_dg)
                        excess = dg - target_dg
                        if not anchored_ext[i]:
                            # Standard case: move sample i to make
                            # ΔG = target_dg.
                            denom = 1.0 / lr + 1.0 / ll
                            new_e = (elevs_ext[i + 1] / lr
                                      + elevs_ext[i - 1] / ll
                                      - target_dg) / denom
                            if abs(new_e - elevs_ext[i]) > 0.001:
                                elevs_ext[i] = new_e
                                changed = True
                                any_change = True
                        else:
                            # Anchored sample: we cannot move it.
                            # Move its non-anchored neighbours to
                            # make ΔG = target_dg.
                            #
                            # We want g_left' = g_left + δ_l,
                            # g_right' = g_right - δ_r, such that
                            # new dg = dg − δ_l − δ_r = target_dg.
                            # → δ_l + δ_r = excess.
                            free_l = not anchored_ext[i - 1]
                            free_r = not anchored_ext[i + 1]
                            if not free_l and not free_r:
                                continue  # both neighbours anchored
                            if free_l and free_r:
                                delta_l = excess / 2.0
                                delta_r = excess / 2.0
                            elif free_l:
                                delta_l = excess
                                delta_r = 0.0
                            else:
                                delta_l = 0.0
                                delta_r = excess
                            if free_l and abs(delta_l) > 1e-9:
                                # increase g_left by delta_l
                                # → decrease elev[i-1] by delta_l*ll
                                new_lo = (
                                    elevs_ext[i - 1] - delta_l * ll)
                                if abs(new_lo
                                        - elevs_ext[i - 1]) > 0.001:
                                    elevs_ext[i - 1] = new_lo
                                    changed = True
                                    any_change = True
                            if free_r and abs(delta_r) > 1e-9:
                                # decrease g_right by delta_r
                                # → decrease elev[i+1] by delta_r*lr
                                new_hi = (
                                    elevs_ext[i + 1] - delta_r * lr)
                                if abs(new_hi
                                        - elevs_ext[i + 1]) > 0.001:
                                    elevs_ext[i + 1] = new_hi
                                    changed = True
                                    any_change = True
                    if not changed:
                        break

                # Copy interior values back to the real arrays
                # (skip the virtual blast pad samples).
                for j in range(len(elevs)):
                    elevs[j] = elevs_ext[j + start_offset]
                return any_change

            _pass_rate_of_change()

            # Outer joint convergence loop: keep alternating the two
            # passes until neither one moves anything.
            for _outer in range(8):
                prev_elevs = list(elevs)
                _pass_hard_cap()
                _pass_rate_of_change()
                if max(abs(a - b) for a, b in zip(elevs, prev_elevs)) < 0.005:
                    break

            # Identify which sample indices correspond to pav_inter-
            # section breakpoints (junction-snap points).  These are
            # the only intermediate vertices retained when consecutive
            # flat segments are consolidated.
            pav_int_indices = set()
            for k, frac in enumerate(fractions):
                for pt in pav_int_t_vals:
                    if abs(frac - pt) < 1e-6:
                        pav_int_indices.add(k)
                        break

            # ── Emit segmented rectangles ────────────────────────────────
            # Per user 2026-05-09: consolidate consecutive flat
            # segments into a single multi-node flat polygon.  Outer
            # corners at the run's start / end samples; intermediate
            # corners at any pav_intersection sample positions inside
            # the run (so adjacent junctions still get their snap
            # points).  Sloped pairs and single-segment flats keep
            # the legacy 4-corner emit.
            FLAT_TOL = 0.05
            n_samples = len(sample_pts)
            idx = 0
            while idx < n_samples - 1:
                end_idx = idx
                while (end_idx < n_samples - 1
                        and abs(elevs[end_idx + 1] - elevs[end_idx])
                        < FLAT_TOL):
                    end_idx += 1
                if end_idx > idx + 1:
                    # Per user 2026-05-12: keep ALL intermediate
                    # samples as corners (both pav_intersection
                    # breakpoints and uniform 100 m seams).  Pre-
                    # Task-1 the runway emitted as separate per-100m
                    # rects whose corners served as junction-snap
                    # points via ``widen_junctions_to_runway_corners``.
                    # Task 1 consolidates consecutive flat segments
                    # into one multi-node polygon — but dropping the
                    # uniform-seam samples broke that snap chain;
                    # junctions in the middle of the flat zone had
                    # no nearby runway corners to share.  Keeping
                    # every sample as a corner preserves the
                    # per-segment granularity inside a single
                    # polygon.
                    intermediate: List[Tuple[float, float, float, bool]] = [
                        sample_pts[k] for k in range(idx + 1, end_idx)
                    ]
                    flat_pts = [sample_pts[idx]] + intermediate \
                        + [sample_pts[end_idx]]
                    samples_ll = [(s[0], s[1]) for s in flat_pts]
                    flat_elev = float(elevs[idx])
                    add_flat_multi_rect(samples_ll, flat_elev, patch_width)
                    runway_chain.append((
                        "MULTI_FLAT", samples_ll, flat_elev, patch_width,
                    ))
                    idx = end_idx
                else:
                    s_a = sample_pts[idx]
                    s_b = sample_pts[idx + 1]
                    add_rect_patch(
                        s_a[0], s_a[1], elevs[idx],
                        s_b[0], s_b[1], elevs[idx + 1],
                        patch_width,
                    )
                    runway_chain.append((
                        s_a[0], s_a[1], elevs[idx],
                        s_b[0], s_b[1], elevs[idx + 1],
                        patch_width,
                    ))
                    idx += 1

            # ── Flat blast-pad / overrun rectangles beyond ends ─────────
            # Length comes from apt.dat row-100 blast_a/blast_b when
            # available; otherwise the legacy 30 m OVERRUN_EXTENSION
            # fallback is already assigned to blast_a/blast_b above.
            # Zero-length blast pads produce nothing.
            if blast_a > 0.1:
                ext_a = extend_point(
                    phys_end_b[0], phys_end_b[1],
                    phys_end_a[0], phys_end_a[1],
                    blast_a,
                )
                add_rect_patch(
                    ext_a[0], ext_a[1], elevs[0],
                    phys_end_a[0], phys_end_a[1], elevs[0],
                    patch_width,
                )
                runway_chain.append((
                    ext_a[0], ext_a[1], elevs[0],
                    phys_end_a[0], phys_end_a[1], elevs[0],
                    patch_width,
                ))
            if blast_b > 0.1:
                ext_b = extend_point(
                    phys_end_a[0], phys_end_a[1],
                    phys_end_b[0], phys_end_b[1],
                    blast_b,
                )
                add_rect_patch(
                    phys_end_b[0], phys_end_b[1], elevs[-1],
                    ext_b[0], ext_b[1], elevs[-1],
                    patch_width,
                )
                runway_chain.append((
                    phys_end_b[0], phys_end_b[1], elevs[-1],
                    ext_b[0], ext_b[1], elevs[-1],
                    patch_width,
                ))

        else:
            # ── Unpaired runway: flat patch at known elevation ───────────
            lat_a, lon_a = data_a["lat"], data_a["lon"]
            elev_a = data_a["elevation_m"]
            displaced_a = data_a.get("displaced_m", 0)
            rwy_width = runway_widths.get(desig_a, DEFAULT_RUNWAY_WIDTH)
            patch_width = rwy_width + 2 * RUNWAY_MARGIN
            ext_dist = max(displaced_a, 50.0)
            ext = extend_point(lat_a, lon_a, lat_a + 0.0001, lon_a, ext_dist)
            ext2 = extend_point(lat_a, lon_a, lat_a - 0.0001, lon_a, ext_dist)
            add_rect_patch(
                ext[0], ext[1], elev_a,
                ext2[0], ext2[1], elev_a,
                patch_width,
            )
            runway_chain.append((
                ext[0], ext[1], elev_a,
                ext2[0], ext2[1], elev_a,
                patch_width,
            ))

    # ── Assemble OSM XML ─────────────────────────────────────────────────────
    lines = [
        "<?xml version='1.0' encoding='UTF-8'?>",
        "<osm version='0.6' upload='false' generator='Ortho4XP_AutoPatch'>",
        "  <!-- Auto-generated runway patch for {} -->".format(icao),
        "  <!-- Source: CIFP/AIRAC threshold elevation data -->",
        "  <!-- This file is overwritten on each build. Manual edits will be lost. -->",
    ]
    for nid, lat, lon in nodes:
        lines.append(
            "  <node id='{}' action='modify' visible='true'"
            " lat='{:.11f}' lon='{:.11f}' />".format(nid, lat, lon)
        )
    for wid, nids, tags in ways:
        lines.append(
            "  <way id='{}' action='modify' visible='true'>".format(wid)
        )
        for nid in nids:
            lines.append("    <nd ref='{}' />".format(nid))
        for k, v in sorted(tags.items()):
            lines.append("    <tag k='{}' v='{}' />".format(k, v))
        lines.append("  </way>")
    lines.append("</osm>")
    return "\n".join(lines) + "\n", runway_chain
