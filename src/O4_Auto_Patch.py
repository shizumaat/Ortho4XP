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

import O4_UI_Utils as UI
import O4_File_Names as FNAMES

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
FT_TO_M = 0.3048
DEG_TO_M = 111120.0  # approximate meters per degree of latitude
DEFAULT_RUNWAY_WIDTH = 45.0  # meters - typical for major runways
RUNWAY_MARGIN = 3.0  # meters added to each side of the runway
DEFAULT_CELL_SIZE = 2.0  # meters between interpolation points (KBNA finding)
DEFAULT_PROFILE = "spline"  # spline profile for natural terrain transitions (KBNA finding)
DEFAULT_STEEPNESS = 2
# How far beyond the physical runway end to extend as a flat apron (meters)
OVERRUN_EXTENSION = 30.0
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
# from OSM centerlines clipped against the apt.dat pavement union.  This is
# the centerline-driven model — every named OSM taxi (A, A1, B, E, F, V…)
# produces aligned rects regardless of which apt.dat polygon contains the
# paint.  Setting "0" reverts to the polygon-driven path that walks
# `apt_twy_rect_chains` from Phase A3 (the prior MRR/skeleton/decompose model).
DEBUG_OSM_CENTERLINES = os.environ.get("O4_OSM_CENTERLINES", "1") == "1"

# New pavement model toggle: when "1" (set via O4_NEW_MODEL=1), Phase C0/C1/C2
# and Phase D are replaced by a single pass that decomposes the pavement union
# into flat N-vertex polygons (one per connected pavement component, at a
# DEM-driven target elevation) plus 4-vertex sloped rectangles at runway-edge
# transitions (using _runway_elev_lookup for the runway-side elevation).
# Phase C3 (building pads) and Phase C4 (transition-strip detection between
# adjacent flat shapes) still run, so building-to-pavement cliffs get bridged
# automatically.  Phase E/F stay gated behind DEBUG_TAXI_ONLY.  Default OFF
# so existing tests keep passing; rollback is `unset O4_NEW_MODEL`.
DEBUG_NEW_MODEL = os.environ.get("O4_NEW_MODEL", "0") == "1"

# Strip-model toggle: when "1" (via O4_STRIP_MODEL=1), pavement is
# decomposed by O4_Pavement_Strips and emitted as a GEOMETRY
# PREVIEW with zero elevations — every way is flat at 0.0 m.  This
# lets you inspect the strip/junction/apron partitioning in JOSM
# through the real patch pipeline, without the elevation solver in
# place.  Mutually exclusive with O4_NEW_MODEL.  Default OFF.
DEBUG_STRIP_MODEL = os.environ.get("O4_STRIP_MODEL", "0") == "1"

# True when ANY of the replacement pavement models is active — the
# legacy Phase C0/C1/C2/C4/D/E1 steps are skipped in that case.
DEBUG_REPLACE_LEGACY_PAVEMENT = DEBUG_NEW_MODEL or DEBUG_STRIP_MODEL


# ──────────────────────────────────────────────────────────────────────────────
# CIFP Coordinate Parsing
# ──────────────────────────────────────────────────────────────────────────────
def parse_cifp_lat(s):
    """Parse ARINC 424 latitude: 'S12002744' → -12.00762222 degrees."""
    hem = s[0]
    deg = int(s[1:3])
    mins = int(s[3:5])
    secs = int(s[5:9]) / 100.0
    decimal = deg + mins / 60.0 + secs / 3600.0
    if hem in ("S", "s"):
        decimal = -decimal
    return decimal


def parse_cifp_lon(s):
    """Parse ARINC 424 longitude: 'W077071686' → -77.12135 degrees."""
    hem = s[0]
    deg = int(s[1:4])
    mins = int(s[4:6])
    secs = int(s[6:10]) / 100.0
    decimal = deg + mins / 60.0 + secs / 3600.0
    if hem in ("W", "w"):
        decimal = -decimal
    return decimal


# ──────────────────────────────────────────────────────────────────────────────
# CIFP File Parsing
# ──────────────────────────────────────────────────────────────────────────────
def parse_cifp_file(filepath):
    """Parse a CIFP .dat file and extract runway threshold data.

    CIFP RWY record format:
        RWY:RW16L,+0580,      ,00044, ,IJCH,3,   ;S12002744,W077071686,0000;

    Fields before semicolon (comma-separated):
        [0] designator  (RW16L)
        [1] mag heading  (+0580 = 058.0°)
        [2] (reserved)
        [3] threshold elevation in feet (00044 = 44 ft)
        [4] (reserved)
        [5] ILS identifier
        [6] ILS category
        [7] (reserved)

    After first semicolon (comma-separated):
        [0] latitude   (S12002744 = S12°00'27.44")
        [1] longitude  (W077071686 = W077°07'16.86")
        [2] displaced threshold distance in feet (0000)

    Returns:
        dict: {designator: {lat, lon, elevation_m, displaced_m}} e.g.
              {'RW16L': {'lat': -12.0076, 'lon': -77.1214, 'elevation_m': 13.41, ...}}
    """
    runways = {}
    try:
        with open(filepath, "r", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line.startswith("RWY:"):
                    continue

                parts = line[4:].split(";")
                if len(parts) < 2:
                    continue

                fields = parts[0].split(",")
                coord_fields = parts[1].split(",")
                if len(fields) < 4 or len(coord_fields) < 3:
                    continue

                designator = fields[0].strip()
                if not designator.startswith("RW"):
                    continue

                # Threshold elevation in feet → meters
                elev_str = fields[3].strip()
                if not elev_str or not elev_str.isdigit():
                    continue
                elevation_m = int(elev_str) * FT_TO_M

                # Threshold coordinates
                lat_str = coord_fields[0].strip()
                lon_str = coord_fields[1].strip()
                if not lat_str or not lon_str:
                    continue
                if len(lat_str) < 9 or len(lon_str) < 10:
                    continue

                try:
                    lat = parse_cifp_lat(lat_str)
                    lon = parse_cifp_lon(lon_str)
                except (ValueError, IndexError):
                    continue

                # Displaced threshold distance (feet → meters)
                displaced_str = coord_fields[2].strip().rstrip(";")
                displaced_m = (
                    int(displaced_str) * FT_TO_M
                    if displaced_str.isdigit()
                    else 0.0
                )

                runways[designator] = {
                    "lat": lat,
                    "lon": lon,
                    "elevation_m": elevation_m,
                    "displaced_m": displaced_m,
                }
    except Exception as e:
        UI.vprint(
            1,
            "   Warning: Could not parse CIFP file",
            filepath,
            ":",
            str(e),
        )
    return runways


# ──────────────────────────────────────────────────────────────────────────────
# Apt.dat Runway Width Parsing
# ──────────────────────────────────────────────────────────────────────────────
def parse_aptdat_runway_widths(aptdat_path, icao):
    """Extract runway widths from an X-Plane apt.dat file for a given airport.

    Apt.dat row code 100 format:
        100 <width_m> <surface> <shoulder> <smoothness> ... <rwy1_name> <lat> <lon> ...
                                                             <rwy2_name> <lat> <lon> ...

    Scans the (potentially large) apt.dat for the airport header (row 1 with
    matching ICAO), then reads its row-100 runways until the next airport.

    Args:
        aptdat_path: Path to apt.dat file.
        icao: ICAO code to search for (e.g. 'SPJC').

    Returns:
        dict: {designator: width_m} e.g. {'RW16R': 45.0, 'RW34L': 45.0}
              Empty dict if airport not found or on error.
    """
    widths = {}
    if not aptdat_path or not os.path.isfile(aptdat_path):
        return widths
    try:
        in_airport = False
        with open(aptdat_path, "r", errors="replace") as f:
            for line in f:
                fields = line.strip().split()
                if not fields:
                    continue
                row_code = fields[0]

                # Airport header: row code 1 (land), 16 (seaplane), 17 (heliport)
                if row_code in ("1", "16", "17"):
                    if in_airport:
                        break  # We've passed our airport, stop
                    # Check if this is our airport (ICAO is field 4)
                    if len(fields) >= 5 and fields[4].upper() == icao.upper():
                        in_airport = True
                    continue

                if not in_airport:
                    continue

                # Row code 100 = runway definition
                if row_code == "100" and len(fields) >= 9:
                    try:
                        width_m = float(fields[1])
                    except ValueError:
                        continue
                    # Runway end 1 name is field 8, end 2 name is further along
                    rwy1_name = fields[8]
                    # Runway end 2: the field index depends on the format.
                    # After rwy1's fields (name, lat, lon, displaced, overrun,
                    # markings, approach_lights, tdz, reil) = 9 fields starting
                    # at index 8, so rwy2 name is at index 17.
                    if len(fields) >= 18:
                        rwy2_name = fields[17]
                    else:
                        rwy2_name = None

                    # Normalize to our RWxx format
                    desig1 = "RW" + rwy1_name if not rwy1_name.startswith("RW") else rwy1_name
                    widths[desig1] = width_m
                    if rwy2_name:
                        desig2 = "RW" + rwy2_name if not rwy2_name.startswith("RW") else rwy2_name
                        widths[desig2] = width_m

    except Exception as e:
        UI.vprint(
            2,
            "   Auto-patch: Could not read apt.dat runway widths:",
            str(e),
        )
    return widths


def find_aptdat(cifp_path):
    """Attempt to locate apt.dat relative to the CIFP directory.

    X-Plane directory structure:
        <X-Plane>/Custom Data/CIFP/           ← cifp_path
        <X-Plane>/Resources/default scenery/default apt dat/Earth nav data/apt.dat
        <X-Plane>/Custom Scenery/Global Airports/Earth nav data/apt.dat

    Also checks for a custom apt.dat in Custom Data/.

    Returns:
        str or None: Path to apt.dat if found.
    """
    if not cifp_path:
        return None

    # cifp_path is typically <X-Plane>/Custom Data/CIFP
    custom_data = os.path.dirname(cifp_path)  # <X-Plane>/Custom Data
    xplane_root = os.path.dirname(custom_data)  # <X-Plane>

    candidates = [
        # Custom apt.dat (user overrides)
        os.path.join(custom_data, "apt.dat"),
        # Global Airports scenery pack
        os.path.join(
            xplane_root,
            "Custom Scenery",
            "Global Airports",
            "Earth nav data",
            "apt.dat",
        ),
        # Default apt dat
        os.path.join(
            xplane_root,
            "Resources",
            "default scenery",
            "default apt dat",
            "Earth nav data",
            "apt.dat",
        ),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def xplane_root_from_cifp_path(cifp_path):
    """Derive the X-Plane installation root from a CIFP directory path.

    cifp_path is typically ``<X-Plane>/Custom Data/CIFP`` (or the
    similar ``Resources/default data/CIFP`` location).  Walks up two
    directory levels to reach the X-Plane root.  Returns ``None`` if
    the path doesn't look right.
    """
    if not cifp_path:
        return None
    try:
        custom_data = os.path.dirname(os.path.normpath(cifp_path))
        root = os.path.dirname(custom_data)
        # Basic sanity check: the derived root should contain a
        # "Custom Scenery" or "Resources" directory.
        if (os.path.isdir(os.path.join(root, "Custom Scenery"))
                or os.path.isdir(os.path.join(root, "Resources"))):
            return root
    except Exception:
        pass
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Runway Pairing
# ──────────────────────────────────────────────────────────────────────────────
def get_reciprocal(designator):
    """Get the reciprocal runway designator. RW16L → RW34R, RW09 → RW27."""
    match = re.match(r"RW(\d{2})([LRC]?)", designator)
    if not match:
        return None
    num = int(match.group(1))
    suffix = match.group(2)
    recip_num = num + 18
    if recip_num > 36:
        recip_num -= 36
    recip_suffix = {"L": "R", "R": "L", "C": "C", "": ""}.get(suffix, "")
    return "RW{:02d}{}".format(recip_num, recip_suffix)


def pair_runways(runways):
    """Match runway thresholds into pairs.

    Returns list of tuples:
        (desig_a, data_a, desig_b, data_b)
    where a is the higher-numbered threshold (higher heading number) by
    convention, and b is the reciprocal. If unpaired, desig_b/data_b are None.
    """
    paired = set()
    pairs = []
    for desig in sorted(runways.keys()):
        if desig in paired:
            continue
        data = runways[desig]
        recip = get_reciprocal(desig)
        if recip and recip in runways:
            paired.add(desig)
            paired.add(recip)
            pairs.append((desig, data, recip, runways[recip]))
        else:
            pairs.append((desig, data, None, None))
    return pairs


# ──────────────────────────────────────────────────────────────────────────────
# Geometry Helpers
# ──────────────────────────────────────────────────────────────────────────────
def runway_corners(lat1, lon1, lat2, lon2, width_m):
    """Compute the 4 corners of a runway rectangle.

    The corners are ordered for the altitude_high/altitude_low patch convention:
        node0 (high-left)  → node1 (low-left) → node2 (low-right) → node3 (high-right)

    In include_patches():
        short_high = way[-2:] = [node3, node0]   (the "high" altitude side)
        short_low  = way[1:3] = [node1, node2]   (the "low" altitude side)

    The "high" side is at (lat1, lon1) and the "low" side at (lat2, lon2).

    Returns:
        list of 4 (lat, lon) tuples, or None if degenerate.
    """
    mid_lat = (lat1 + lat2) / 2.0
    cos_lat = cos(mid_lat * pi / 180.0)
    if cos_lat < 1e-6:
        cos_lat = 1e-6

    # Direction vector in meters
    dx_m = (lon2 - lon1) * cos_lat * DEG_TO_M
    dy_m = (lat2 - lat1) * DEG_TO_M
    length_m = sqrt(dx_m ** 2 + dy_m ** 2)
    if length_m < 1.0:
        return None

    # Perpendicular unit vector (90° clockwise rotation)
    perp_dx_m = -dy_m / length_m
    perp_dy_m = dx_m / length_m

    # Half-width offset in degrees
    half_w = width_m / 2.0
    perp_dlon = (perp_dx_m * half_w) / (cos_lat * DEG_TO_M)
    perp_dlat = (perp_dy_m * half_w) / DEG_TO_M

    # Corners: high-left, low-left, low-right, high-right
    c0 = (lat1 + perp_dlat, lon1 + perp_dlon)
    c1 = (lat2 + perp_dlat, lon2 + perp_dlon)
    c2 = (lat2 - perp_dlat, lon2 - perp_dlon)
    c3 = (lat1 - perp_dlat, lon1 - perp_dlon)
    return [c0, c1, c2, c3]


def extend_point(lat_from, lon_from, lat_to, lon_to, distance_m):
    """Extend a point beyond lat_to/lon_to by distance_m meters
    along the direction from lat_from/lon_from to lat_to/lon_to.

    Returns (lat, lon) of the extended point.
    """
    mid_lat = (lat_from + lat_to) / 2.0
    cos_lat = cos(mid_lat * pi / 180.0)
    if cos_lat < 1e-6:
        cos_lat = 1e-6

    dx_m = (lon_to - lon_from) * cos_lat * DEG_TO_M
    dy_m = (lat_to - lat_from) * DEG_TO_M
    length_m = sqrt(dx_m ** 2 + dy_m ** 2)
    if length_m < 1.0:
        return (lat_to, lon_to)

    # Unit direction vector
    ux = dx_m / length_m
    uy = dy_m / length_m

    ext_dlon = (ux * distance_m) / (cos_lat * DEG_TO_M)
    ext_dlat = (uy * distance_m) / DEG_TO_M
    return (lat_to + ext_dlat, lon_to + ext_dlon)


# Adaptive triangulation lives in O4_Surface_Mesh as
# ``adaptive_triangulate`` and is exercised by tests/test_surface_mesh.py.
# Phases C2 (apron) and D (junction) will be wired to call it in
# commits 5 and 6.  The deleted inline definition below this comment
# was a transitional copy.


# ──────────────────────────────────────────────────────────────────────────────
# OSM Patch Generation
# ──────────────────────────────────────────────────────────────────────────────
RUNWAY_SEGMENT_LENGTH = 100.0  # meters — length of each runway segment


def generate_patch_osm(icao, runway_pairs, runway_widths=None, tile=None,
                       apt_runways=None, extra_anchors=None):
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

    def _sample_dem(lat, lon):
        """Sample DEM elevation at a lat/lon, returning 0 on failure."""
        if tile is None or not hasattr(tile, "dem") or tile.dem is None:
            return None
        try:
            return tile.dem.alt((lon - tile.lon, lat - tile.lat))
        except Exception:
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

            # Ensure thresholds are in the list
            if phys_dist > 0:
                t_a = (displaced_a / phys_dist) if displaced_a > 0 else 0.0
                t_b = 1.0 - (displaced_b / phys_dist) if displaced_b > 0 else 1.0
                for t in [t_a, t_b]:
                    if 0.0 < t < 1.0:
                        # Insert if not already near an existing fraction
                        if not any(abs(t - f) < 0.01 for f in fractions):
                            fractions.append(t)
                fractions.sort()

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
                    try:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [auto-patch] {icao} runway "
                            f"{desig_a}/{desig_b}: anchor pair grade "
                            f"{g * 100:.2f}% > "
                            f"{MAX_RUNWAY_GRADE * 100:.1f}% between "
                            f"fractions {f0:.3f} and {f1:.3f} "
                            f"({seg_d:.0f} m apart, ΔE={e1 - e0:+.2f} m) "
                            f"— profile will be infeasible at FAA "
                            f"runway max grade.\n")
                    except Exception:
                        pass

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

            # ── Emit segmented rectangles ────────────────────────────────
            for idx in range(len(sample_pts) - 1):
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


# ──────────────────────────────────────────────────────────────────────────────
# CIFP Directory Discovery
# ──────────────────────────────────────────────────────────────────────────────
def discover_cifp_airports(cifp_path):
    """Scan a CIFP directory for airport .dat files.

    Returns:
        dict: {icao_code: filepath} for all discovered airports.
    """
    airports = {}
    if not cifp_path or not os.path.isdir(cifp_path):
        return airports
    for fname in os.listdir(cifp_path):
        if fname.lower().endswith(".dat"):
            icao = fname[:-4].upper()
            # Basic ICAO code validation: 2-4 alphanumeric characters
            if 2 <= len(icao) <= 4 and icao.replace("-", "").isalnum():
                airports[icao] = os.path.join(cifp_path, fname)
    return airports


def airport_in_tile(runways, tile_lat, tile_lon):
    """Check if any runway threshold falls within a 1°×1° tile.

    Args:
        runways: dict from parse_cifp_file()
        tile_lat: integer latitude of tile's SW corner
        tile_lon: integer longitude of tile's SW corner

    Returns:
        bool
    """
    for data in runways.values():
        lat = data["lat"]
        lon = data["lon"]
        if tile_lat <= lat < tile_lat + 1 and tile_lon <= lon < tile_lon + 1:
            return True
    return False


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
            from O4_Airport_Pavement_Builder import build_airport_pavement
            layout = build_airport_pavement(icao, xp_root)
        except Exception as _e:
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
        except Exception as e:
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


# ──────────────────────────────────────────────────────────────────────────────
# Taxiway Geometry Extraction and Patch Generation
# ──────────────────────────────────────────────────────────────────────────────
DEFAULT_TAXIWAY_WIDTH = 15.0  # meters - typical taxiway width
MAX_TAXIWAY_GRADE = 0.015     # 1.5% max grade for taxiways


def extract_taxiway_info(airport_layer, dico_airports, tile):
    """Extract taxiway centerline data from the OSM airport layer.

    Pulls taxiway way coordinates from the airport_layer's node/way dicts,
    keyed by the airport identifier used in dico_airports (typically ICAO).

    Args:
        airport_layer: OSM_layer object with dicosmn (nodes) and dicosmw (ways).
        dico_airports: dict from discover_airport_names / attach_surfaces, keyed
                       by airport identifier.
        tile: Tile object with .lat and .lon.

    Returns:
        dict: {airport_key: [{'centerline': [(lon, lat), ...],
                              'wayid': int,
                              'name': str}, ...]}
              Coordinates are ABSOLUTE (not tile-relative).
              ``name`` is the OSM ``ref`` tag (e.g. "A", "A1", "F",
              "V") if present, otherwise an empty string.  Used by
              the centerline-driven Phase C0 emission to group
              centerlines by named taxiway.
    """
    way_tags = getattr(airport_layer, "dicosmtags", {}).get("w", {})
    result = {}
    for airport in dico_airports:
        apt = dico_airports[airport]
        # taxiway is (MultiPolygon_area, wayid_list) after build_taxiway_areas
        if not isinstance(apt.get("taxiway"), tuple) or len(apt["taxiway"]) < 2:
            continue
        wayid_list = apt["taxiway"][1]
        if not wayid_list:
            continue
        taxiways = []
        for wayid in wayid_list:
            if wayid not in airport_layer.dicosmw:
                continue
            node_ids = airport_layer.dicosmw[wayid]
            coords = []
            for nid in node_ids:
                if nid in airport_layer.dicosmn:
                    coords.append(tuple(airport_layer.dicosmn[nid]))
            if len(coords) >= 2:
                tags = way_tags.get(wayid, {})
                taxiways.append({
                    "centerline": coords,
                    "wayid": wayid,
                    "name": tags.get("ref", "") or "",
                })
        if taxiways:
            result[airport] = taxiways
    return result


def extract_building_info(airport_layer, dico_airports, tile,
                          building_layer=None):
    """Extract building footprints within airport boundaries.

    Collects building polygons from three sources:
    1. Hangars already attached to airports (aeroway=hangar from dico_airports)
    2. Terminals from the airport_layer (aeroway=terminal ways)
    3. General buildings (building=*) from a separate building_layer, filtered
       to those intersecting airport boundaries.

    Args:
        airport_layer: OSM_layer with aeroway data (includes terminal ways).
        dico_airports: dict with processed airport data including 'hangar'
                       MultiPolygon and 'boundary' MultiPolygon.
        tile: Tile object with .lat and .lon.
        building_layer: Optional OSM_layer with building=* data. If None,
                        only hangars and terminals are extracted.

    Returns:
        dict: {airport_key: [{'footprint': [(lon, lat), ...], 'source': str}, ...]}
              Coordinates are ABSOLUTE (not tile-relative).
              Each footprint is a list of (lon, lat) ring coordinates.
    """
    result = {}
    for airport in dico_airports:
        apt = dico_airports[airport]
        buildings = []

        # Source 1: Hangars from airport_layer (already attached to airports)
        hangar_wayids = []
        # Before build_hangar_areas, apt["hangar"] is a list of wayids.
        # After build_hangar_areas, it's a MultiPolygon (tile-relative).
        # We need to handle both cases, but typically this runs after
        # build_hangar_areas, so we work with the MultiPolygon.
        hangar_geom = apt.get("hangar")
        if hangar_geom is not None and hasattr(hangar_geom, "geoms"):
            # It's a MultiPolygon in tile-relative coords - convert back
            for poly in hangar_geom.geoms:
                if poly.is_empty or not poly.is_valid or poly.area < 1e-10:
                    continue
                # Convert from tile-relative to absolute
                coords = [
                    (x + tile.lon, y + tile.lat)
                    for x, y in poly.exterior.coords
                ]
                buildings.append({
                    "footprint": coords,
                    "source": "hangar",
                })

        # Source 2: Terminals from airport_layer (aeroway=terminal)
        # Check both simple ways and multipolygon relations, since large
        # terminals are often mapped as relations in OSM.
        boundary = apt.get("boundary")
        if boundary is not None:
            from shapely import geometry as shp_geom

            def _try_add_terminal(coords, source_label):
                """Add terminal footprint if it intersects this airport."""
                if len(coords) < 4:
                    return
                try:
                    term_rel = shp_geom.Polygon([
                        (c[0] - tile.lon, c[1] - tile.lat)
                        for c in coords
                    ])
                    if (term_rel.is_valid and not term_rel.is_empty
                            and boundary.intersects(term_rel)):
                        buildings.append({
                            "footprint": coords,
                            "source": source_label,
                        })
                except Exception:
                    pass

            # 2a: Terminal ways
            for wayid in airport_layer.dicosmw:
                wtags = airport_layer.dicosmtags.get("w", {}).get(wayid, {})
                if wtags.get("aeroway") != "terminal":
                    continue
                node_ids = airport_layer.dicosmw[wayid]
                if len(node_ids) < 4 or node_ids[0] != node_ids[-1]:
                    continue
                coords = []
                for nid in node_ids:
                    if nid in airport_layer.dicosmn:
                        coords.append(tuple(airport_layer.dicosmn[nid]))
                _try_add_terminal(coords, "terminal")

            # 2b: Terminal relations (multipolygons)
            for relid in airport_layer.dicosmr:
                rtags = airport_layer.dicosmtags.get("r", {}).get(relid, {})
                if rtags.get("aeroway") != "terminal":
                    continue
                rel_data = airport_layer.dicosmr[relid]
                outer_rings = rel_data.get("outer", [])
                for ring_nids in outer_rings:
                    coords = []
                    for nid in ring_nids:
                        if nid in airport_layer.dicosmn:
                            coords.append(
                                tuple(airport_layer.dicosmn[nid])
                            )
                    _try_add_terminal(coords, "terminal")

        # Source 3: General buildings from building_layer
        if building_layer is not None and boundary is not None:
            try:
                for wayid in building_layer.dicosmw:
                    if wayid not in building_layer.dicosmtags.get("w", {}):
                        continue
                    wtags = building_layer.dicosmtags["w"][wayid]
                    if "building" not in wtags:
                        continue
                    node_ids = building_layer.dicosmw[wayid]
                    if len(node_ids) < 4:
                        continue
                    if node_ids[0] != node_ids[-1]:
                        continue
                    coords = []
                    for nid in node_ids:
                        if nid in building_layer.dicosmn:
                            coords.append(tuple(building_layer.dicosmn[nid]))
                    if len(coords) < 4:
                        continue
                    # Coords are absolute (lon, lat); boundary is
                    # tile-relative, so convert building to tile-relative
                    try:
                        from shapely import geometry as shp_geom
                        bldg_rel = shp_geom.Polygon([
                            (c[0] - tile.lon, c[1] - tile.lat)
                            for c in coords
                        ])
                        if (bldg_rel.is_valid
                                and not bldg_rel.is_empty
                                and boundary.intersects(bldg_rel)):
                            # Check it's not already covered by a hangar
                            # or terminal
                            is_dup = False
                            for existing in buildings:
                                try:
                                    existing_rel = shp_geom.Polygon([
                                        (c[0] - tile.lon, c[1] - tile.lat)
                                        for c in existing["footprint"]
                                    ])
                                    if (existing_rel.intersection(bldg_rel).area
                                            > 0.5 * bldg_rel.area):
                                        is_dup = True
                                        break
                                except Exception:
                                    pass
                            if not is_dup:
                                buildings.append({
                                    "footprint": coords,
                                    "source": "building",
                                })
                    except Exception:
                        pass
            except Exception:
                pass

        if buildings:
            result[airport] = buildings
            sources = {}
            for b in buildings:
                s = b.get("source", "unknown")
                sources[s] = sources.get(s, 0) + 1
            parts = ", ".join(
                "{} {}s".format(v, k) for k, v in sorted(sources.items())
            )
            UI.vprint(1, "   Auto-patch: {} buildings for {}: {}".format(
                len(buildings), airport, parts
            ))
    return result


# FAA AC 150/5300-13B grade limits for Approach Category C-E airports:
#   Runway longitudinal:  max 1.5%
#   Taxiway longitudinal: max 1.5%
#   Apron any direction:  max 1.0%
#   Rate of grade change: max 1% per 30m vertical curve
MAX_TAXIWAY_GRADE = 0.015     # 1.5% max longitudinal grade for taxiways
MAX_APRON_GRADE = 0.010       # 1.0% max grade in any direction for aprons
MAX_SURFACE_GRADE = 0.015     # 1.5% general transition grade limit
MAX_RUNWAY_GRADE = 0.015      # 1.5% max longitudinal grade for runways (Cat C-E)

# FAA AC 150/5300-13B vertical-curve rules:
#
#   Runways, Design Group III / C-III and up (every commercial
#   jetport, SPJC included):
#     L ≥ 1000 ft × |ΔG|   (L in feet, ΔG in %)
#     ≡ L ≥ 305 m × |ΔG|    (L in metres, ΔG in %)
#   A 1% grade change therefore requires a 305 m vertical curve,
#   giving a maximum grade change per metre of pavement of
#   ≈ 0.0000328.  This is what makes real runways look like long
#   gentle slopes with very gradual transitions.
#
#   Taxiways, same design group:
#     L ≥ 100 ft × |ΔG|  ≡  L ≥ 30.5 m × |ΔG|
#   i.e. 1% grade change per 30 m, ~10× steeper allowed than a
#   runway.
MAX_RUNWAY_GRADE_CHANGE_PER_M = 1.0 / 30000.0  # ≈ 0.0033% per metre
MAX_TAXIWAY_GRADE_CHANGE_PER_M = 1.0 / 3000.0   # ≈ 0.0333% per metre
# Back-compat alias: older code paths still reference this.
MAX_GRADE_CHANGE_PER_M = MAX_TAXIWAY_GRADE_CHANGE_PER_M

GRADE_RELAX_ITERATIONS = 80   # iterations for the elevation solver (increased)
TAXIWAY_BUFFER_WIDTH = 12.0   # meters half-width for taxiway surface area
APRON_BUFFER = 5.0            # meters extra buffer around aprons
DEM_VARIANCE_THRESHOLD = 0.5  # meters — only emit transition patches above this

BOUNDARY_SEG_LENGTH = 30.0    # meters — max edge for densification (matches FAA
                              #          rate-of-change rule: 1% per 30m)
TRANSITION_BUFFER_M = 50.0    # meters — buffer around buildings for transition zone
FLAT_TERRAIN_THRESHOLD = 2.0  # meters — max DEM range to consider "flat"
COMPLEX_APRON_ELEV_RANGE = 1.0  # meters — apron elev range triggering triangulation
MAX_TRANSITION_GAP = 5.0      # meters — detect transitions within this gap
TRANSITION_STRIP_MIN_WIDTH = 15.0  # meters — minimum transition strip width
# Road-related constants
DEFAULT_ROAD_WIDTH = 12.0     # meters — typical highway lane width × 2
MAX_HIGHWAY_GRADE = 0.06      # 6% max grade for highways (AASHTO)
ROAD_TERRAIN_INFLUENCE = 50.0 # meters — how far road elevation influences terrain
# Drainage constants (revised: sloped ditch model replaces flat-drop)
DRAINAGE_MAX_DEPTH = 1.5      # meters max — center of drainage ditch below edge
DRAINAGE_MIN_AREA = 200.0     # m² — minimum infield area to model drainage
DRAINAGE_EDGE_BUFFER = 2.0    # meters — drainage starts this far from pavement


def extract_road_info(dico_airports, tile, road_layer=None):
    """Extract road geometry near airports for terrain context.

    Roads (especially highways) near airports constrain terrain modeling:
    - A road along the airport boundary at a different elevation implies
      a retaining wall or grade transition.
    - Road centerlines provide known terrain elevations outside the airport.
    - Tunnels under the airport imply unrelated road elevation.

    Args:
        dico_airports: dict with processed airport data.
        tile: Tile object with .lat, .lon, .dem.
        road_layer: Optional OSM_layer with highway=* data.

    Returns:
        dict: {airport_key: [{'centerline': [(lon, lat), ...],
                              'highway_type': str, 'tunnel': bool,
                              'bridge': bool}, ...]}
    """
    result = {}
    if road_layer is None:
        return result

    for airport in dico_airports:
        apt = dico_airports[airport]
        boundary = apt.get("boundary")
        if boundary is None:
            continue

        roads = []
        try:
            from shapely import geometry as shp_geom
            # Buffer boundary to find nearby roads
            boundary_buf = boundary.buffer(
                ROAD_TERRAIN_INFLUENCE / DEG_TO_M)

            for wayid in road_layer.dicosmw:
                wtags = road_layer.dicosmtags.get("w", {}).get(wayid, {})
                hw_type = wtags.get("highway", "")
                if hw_type not in ("motorway", "trunk", "primary",
                                   "secondary", "tertiary",
                                   "motorway_link", "trunk_link",
                                   "primary_link"):
                    continue

                node_ids = road_layer.dicosmw[wayid]
                coords = []
                for nid in node_ids:
                    if nid in road_layer.dicosmn:
                        coords.append(tuple(road_layer.dicosmn[nid]))
                if len(coords) < 2:
                    continue

                # Check if road is near airport (tile-relative coords)
                road_rel = [(c[0] - tile.lon, c[1] - tile.lat)
                            for c in coords]
                try:
                    ls = shp_geom.LineString(road_rel)
                    if not boundary_buf.intersects(ls):
                        continue
                except Exception:
                    continue

                is_tunnel = wtags.get("tunnel") in ("yes", "building_passage")
                is_bridge = wtags.get("bridge") == "yes"
                roads.append({
                    "centerline": coords,
                    "highway_type": hw_type,
                    "tunnel": is_tunnel,
                    "bridge": is_bridge,
                })
        except Exception:
            pass

        if roads:
            result[airport] = roads
            UI.vprint(2, "   Auto-patch: {} roads near {}".format(
                len(roads), airport))
    return result


