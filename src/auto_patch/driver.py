"""Auto-generate runway slope patches from CIFP/AIRAC aeronautical data.

This module parses ARINC 424 (CIFP) data files to extract precise runway
threshold elevations and coordinates, then generates .patch.osm files that
provide accurate runway slope profiles. These auto-patches replace the
default polynomial-fit altitude model with authoritative aeronautical data.

Auto-generated patches are named {ICAO}_auto.patch.osm and are given lower
priority than user-provided manual patches.
"""

import os
import re
from math import cos, sin, pi, sqrt, floor, atan2, acos

from shapely import geometry as shp_geom
from shapely import ops as shp_ops
from shapely.errors import GEOSException, TopologicalError

# Driver harness tuple — covers expected runtime failure modes for a
# per-airport pass.  Specifically OMITS NameError / AttributeError /
# ImportError so typos and broken imports propagate immediately
# rather than being silently logged and skipped.
_DRIVER_EXC = (OSError, ValueError, TypeError, KeyError,
               IndexError, RuntimeError,
               GEOSException, TopologicalError)

import O4_UI_Utils as UI
import O4_File_Names as FNAMES
from .cifp_reader import (
    airport_in_tile,
    discover_cifp_airports,
    find_aptdat,
    parse_cifp_file,
    xplane_root_from_cifp_path,
)
from .pavement.runway_geometry import (
    DEFAULT_RUNWAY_WIDTH,
    extend_point,
    pair_runways,
    parse_aptdat_runway_widths,
    runway_corners,
)

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
FT_TO_M = 0.3048  # left over from the dead surface-patches code (slice 0)
# DEFAULT_RUNWAY_WIDTH imported from O4_Runway_Geometry above.

# FAA AC 150/5300-13B grade limits for Approach Category C-E airports.
MAX_TAXIWAY_GRADE = 0.015     # 1.5% max longitudinal grade for taxiways

# FAA vertical-curve rules — taxiway counterpart of the runway value.
MAX_TAXIWAY_GRADE_CHANGE_PER_M = 1.0 / 3000.0

DEFAULT_STEEPNESS = 2
# Maximum number of runway chunks for a single patch polygon
MAX_NODE_ID = -1  # will be decremented for each new node

# Debug toggle: when set to "1" via O4_DEBUG_TAXI_ONLY env var, the surface
# generator skips every non-taxi emission phase (apron, buildings, coverage
# fill, transition strips, junction triangles, boundary band, tunnel portals,
# drainage).  Runway segments still come out via generate_patch_osm.  This
# isolates the Phase C0 taxi rect output for visual debugging without the
# rest of the pipeline obscuring it.
DEBUG_TAXI_ONLY = os.environ.get("O4_DEBUG_TAXI_ONLY", "0") == "1"

# Phase C0 source toggle: when "1" (the default), Phase C0 emits taxi rects
# from OSM centerlines clipped against the apt.dat pavement union.
DEBUG_OSM_CENTERLINES = os.environ.get("O4_OSM_CENTERLINES", "1") == "1"

# New pavement model toggle.
DEBUG_NEW_MODEL = os.environ.get("O4_NEW_MODEL", "0") == "1"

# Strip-model toggle.
DEBUG_STRIP_MODEL = os.environ.get("O4_STRIP_MODEL", "0") == "1"

# True when ANY of the replacement pavement models is active.
DEBUG_REPLACE_LEGACY_PAVEMENT = DEBUG_NEW_MODEL or DEBUG_STRIP_MODEL


# ──────────────────────────────────────────────────────────────────────────────
# Runway-segment patch emission (re-exported from
# O4_Pavement_Runway_Segments)
# ──────────────────────────────────────────────────────────────────────────────
from .pavement.runway_segments import (
    DEFAULT_CELL_SIZE,
    DEFAULT_PROFILE,
    DEG_TO_M,
    GRADE_RELAX_ITERATIONS,
    MAX_RUNWAY_GRADE,
    MAX_RUNWAY_GRADE_CHANGE_PER_M,
    OVERRUN_EXTENSION,
    RUNWAY_MARGIN,
    RUNWAY_SEGMENT_LENGTH,
    generate_patch_osm,
)




# ──────────────────────────────────────────────────────────────────────────────
# Main Entry Point
# ──────────────────────────────────────────────────────────────────────────────
def generate_auto_patches(tile, cifp_path, taxiway_data=None,
                          building_data=None, dico_airports=None,
                          road_data=None, mode="ICAO"):
    """Generate auto-patch files for all CIFP airports within a tile.

    Scans the CIFP directory for airport data files, parses runway threshold
    data, and writes {ICAO}_auto.patch.osm files into the tile's Patches
    directory.

    Auto-patches cover the full airport surface as a single non-overlapping
    mesh when building data is available:
    1. Runway slope patches from CIFP threshold elevations
    2. Building flattening (flat altitude=N footprints)
    3. Grade-limited transition triangles (taxiways, aprons, surrounding area)

    When no building data is available, falls back to runway-only patches
    using altitude_high/altitude_low rectangles.

    Args:
        tile: Tile object with .lat, .lon, and .dem attributes.
        cifp_path: Path to the CIFP data directory.
        taxiway_data: Optional dict from extract_taxiway_info().
        building_data: Optional dict from extract_building_info().
        dico_airports: Optional dict with processed airport data (provides
                       apron geometry and boundaries).
        mode: "ICAO" (default) only patches airports with a 4-letter ICAO
              code; "All" patches every CIFP airport regardless of code
              format. ("None" is handled at the call site by skipping this
              function entirely.)

    Returns:
        list: ICAO codes of airports for which auto-patches were generated.
    """
    if not cifp_path or not os.path.isdir(cifp_path):
        UI.vprint(
            1,
            "   Auto-patch: CIFP directory not found at",
            cifp_path,
            ", skipping.",
        )
        return []

    if building_data is None:
        building_data = {}
    if taxiway_data is None:
        taxiway_data = {}

    tile_lat = int(floor(tile.lat))
    tile_lon = int(floor(tile.lon))
    patch_dir = FNAMES.patch_dir(tile_lat, tile_lon)

    # Discover which manual patches already exist
    manual_patches = set()
    if os.path.exists(patch_dir):
        for fname in os.listdir(patch_dir):
            if fname.endswith(".patch.osm") and "_auto.patch.osm" not in fname:
                # Extract probable ICAO code from filename
                base = fname[:-10]  # strip .patch.osm
                # The ICAO prefix is the part before any underscore, or the
                # whole base name if no underscore
                icao_prefix = base.split("_")[0].upper()
                manual_patches.add(icao_prefix)
            elif os.path.isdir(os.path.join(patch_dir, fname)):
                manual_patches.add(fname.upper())

    # Scan all CIFP airports
    cifp_airports = discover_cifp_airports(cifp_path)
    auto_patched = []

    for icao, filepath in sorted(cifp_airports.items()):
        # In ICAO mode, only patch airports with a real 4-letter ICAO code
        # (skip 3-letter FAA codes and alphanumeric local-use codes like "1A2")
        if mode == "ICAO" and not (len(icao) == 4 and icao.isalpha()):
            UI.vprint(
                2,
                "   Auto-patch: Skipping",
                icao,
                "(non-ICAO code, mode=ICAO).",
            )
            continue
        # Skip if a manual patch already covers this airport
        if icao in manual_patches:
            UI.vprint(
                2,
                "   Auto-patch: Skipping",
                icao,
                "(manual patch exists).",
            )
            continue

        # Parse runway data
        runways = parse_cifp_file(filepath)
        if not runways:
            continue

        # Check if any runway falls within this tile
        if not airport_in_tile(runways, tile_lat, tile_lon):
            continue

        # Pair runways and generate patch
        pairs = pair_runways(runways)
        if not pairs:
            continue

        # Look up actual runway widths from apt.dat
        runway_widths = {}
        aptdat_path = find_aptdat(cifp_path)
        if aptdat_path:
            runway_widths = parse_aptdat_runway_widths(aptdat_path, icao)
            if runway_widths:
                UI.vprint(
                    2,
                    "   Auto-patch: Got runway widths from apt.dat for",
                    icao,
                )

        # Build runway pair data for elevation interpolation (shared by
        # taxiway and building patch generation)
        rwy_pairs_for_elev = []
        for desig_a, data_a, desig_b, data_b in pairs:
            if data_b is not None:
                rwy_pairs_for_elev.append({
                    "data_a": data_a,
                    "data_b": data_b,
                    "desig_a": desig_a,
                    "desig_b": desig_b,
                })

        has_dem = hasattr(tile, "dem") and tile.dem is not None

        # Collect taxiway and building data for this airport
        airport_taxiways = (
            taxiway_data.get(icao)
            or taxiway_data.get(icao.upper())
            or taxiway_data.get(icao.lower())
        ) if taxiway_data else None
        airport_buildings = (
            building_data.get(icao)
            or building_data.get(icao.upper())
            or building_data.get(icao.lower())
        ) if building_data else None

        # Look up the airport's processed data (aprons, boundary, etc.)
        dico_apt_entry = {}
        if dico_airports:
            dico_apt_entry = (
                dico_airports.get(icao)
                or dico_airports.get(icao.upper())
                or dico_airports.get(icao.lower())
                or {}
            )

        # ── Generate the patch content ──────────────────────────────────
        # Use the new O4_Airport_Pavement_Builder pipeline.  It produces
        # segmented sloped runway rects, grade-compliant taxi rects,
        # non-overlapping junction polygons, and terminal pads, all in
        # one self-contained call from the same CIFP + apt.dat + OSM +
        # DEM inputs the legacy pipeline used.
        xp_root = xplane_root_from_cifp_path(cifp_path)
        if xp_root is None:
            UI.vprint(
                1, "   Auto-patch: Skipping", icao,
                "(cannot resolve X-Plane root from CIFP path).")
            continue
        try:
            from .pipeline import build_airport_pavement
            # Forward Ortho4XP-side per-airport data already
            # computed during the tile pipeline:
            #   * ``airport_taxiways`` from ``extract_taxiway_info``
            #     (above) — OSM centerlines, used by Pipeline's
            #     centerline-union helper.
            #   * ``tile.dem`` — Ortho4XP's post-smoothing DEM,
            #     reused by Phase-2 elevation + boundary emit so
            #     auto_patch reads the SAME smoothed DEM that
            #     drives flattening, instead of loading its own.
            #     User 2026-05-07: tested raw-DEM switch; reverted
            #     because it didn't fix the SPJC apron rough spot
            #     and risked re-introducing terrain artefacts the
            #     smoothing was added to remove.
            #   * ``dico_apt_entry['boundary']`` — currently
            #     reserved (parameter slot only); auto_patch's
            #     boundary emit still derives its outline from
            #     apt.dat row-130 until source-of-truth chosen.
            layout = build_airport_pavement(
                icao, xp_root,
                taxiway_data=airport_taxiways,
                tile_dem=getattr(tile, "dem", None),
                airport_boundary=dico_apt_entry.get("boundary")
                                  if dico_apt_entry else None,
                # Per user 2026-05-12: pass the CURRENT tile being
                # processed by Ortho4XP (not the airport's anchor
                # tile) so ``tile_cut`` drops the right shape
                # pieces.  Cross-tile airports (e.g. SPLP at
                # lon ~ -77.0) get processed multiple times — once
                # per tile they touch — and each pass should
                # produce only the in-tile shapes.
                current_tile_lat=tile_lat,
                current_tile_lon=tile_lon,
            )
        except _DRIVER_EXC as _e:
            UI.vprint(
                1, "   Auto-patch: Pavement builder failed for",
                icao, ":", str(_e))
            continue

        # Write the auto-patch file
        if not os.path.exists(patch_dir):
            os.makedirs(patch_dir)

        auto_patch_file = os.path.join(
            patch_dir, "{}_auto.patch.osm".format(icao)
        )
        try:
            layout.to_osm(auto_patch_file)
            # Classify shapes for the status line.
            from collections import Counter as _Counter
            counts = _Counter(s.role for s in layout.shapes)
            summary = " + ".join(
                "{} {}".format(n, r) for r, n in
                sorted(counts.items(), key=lambda x: -x[1]))
            UI.vprint(
                1, "   Auto-patch: Generated", icao,
                "(" + summary + ")")
            auto_patched.append(icao)
        except _DRIVER_EXC as e:
            UI.vprint(
                1,
                "   Auto-patch: Failed to write",
                auto_patch_file,
                ":",
                str(e),
            )

    if auto_patched:
        UI.vprint(
            0,
            "   Auto-patch: Generated patches for {} airports.".format(
                len(auto_patched)
            ),
        )
    else:
        UI.vprint(2, "   Auto-patch: No airports with CIFP data in this tile.")

    return auto_patched


