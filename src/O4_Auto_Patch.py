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
# apt.dat → dico_apt_entry adapter
# ──────────────────────────────────────────────────────────────────────────────
# Phase A0.5 in generate_airport_surface_patches calls this to replace
# the OSM-derived shapes in dico_apt_entry with apt.dat-derived ones.
# The downstream phases (A1 runway clip poly, A3 apron polys, etc.)
# read from the dict, so this adapter returns a dict with the same
# keys in the same format — but with apt.dat geometry.
#
# Key design choice for commit 8 (STATUS.md option (a)): every apt.dat
# row-110 pavement goes into dico_apt_entry["apron"], and
# dico_apt_entry["taxiway"] is left empty.  Downstream Phase C2
# triangulates everything at the apron grade (1.0 %).  This loses
# the looser 1.5 % taxiway grade but keeps the pipeline simple; name-
# based classification is a follow-up commit.
#
# Coordinates: apt.dat polygons are in absolute lat/lon, but
# dico_apt_entry stores tile-relative shapes (x = lon - tile.lon,
# y = lat - tile.lat).  This adapter does the conversion per vertex.
def _dico_from_apt_dat(apt_data, tile, fallback_dico):
    """Build a dico_apt_entry-shaped dict from an apt.dat Airport.

    Args:
        apt_data: ``O4_Apt_Dat_Reader.Airport`` parsed from apt.dat.
        tile: the Tile we're rendering (for tile-relative coords).
        fallback_dico: the existing dico_apt_entry dict; used to
            inherit OSM-derived fields we don't have in apt.dat
            (hangar, repr_node, boundary fallback, etc.).

    Returns:
        A new dict that shadows ``fallback_dico`` but with
        apt.dat-derived runway, apron, taxiway, and boundary shapes.
        Always returns a dict even if apt.dat is partially
        populated — missing fields stay as whatever fallback_dico
        had.
    """
    from shapely import geometry as _shp_geom
    from shapely import ops as _shp_ops

    def _to_tile_relative(polygon):
        """Convert absolute lat/lon polygon to tile-relative
        (lon - tile.lon, lat - tile.lat) with holes preserved.
        """
        try:
            exterior = [
                (x - tile.lon, y - tile.lat)
                for x, y in polygon.exterior.coords
            ]
            holes = []
            for ring in polygon.interiors:
                holes.append(
                    [(x - tile.lon, y - tile.lat)
                     for x, y in ring.coords])
            p = _shp_geom.Polygon(exterior, holes if holes else None)
            if not p.is_valid:
                p = p.buffer(0)
            if p.is_empty:
                return None
            if not isinstance(p, _shp_geom.Polygon):
                if hasattr(p, "geoms"):
                    p = max(p.geoms, key=lambda g: g.area)
                else:
                    return None
            return p
        except Exception:
            return None

    new_dico = dict(fallback_dico) if fallback_dico else {}

    # ── 1. Aprons: every apt.dat pavement polygon ────────────────────
    apron_polys = []
    for pav in apt_data.pavements:
        p_tr = _to_tile_relative(pav.polygon)
        if p_tr is not None and not p_tr.is_empty:
            apron_polys.append(p_tr)
    if apron_polys:
        try:
            apron_union = _shp_ops.unary_union(apron_polys)
            # Store as a tuple (geom, wayid_list) to match the
            # format the legacy expects — see O4_Airport_Utils
            # build_apron_areas.  We have no OSM way IDs, so the
            # list is empty.
            new_dico["apron"] = (apron_union, [])
        except Exception:
            pass

    # ── 2. Taxiway: empty for commit 8 (all pavement is "apron") ────
    # Leaving this as an empty tuple signals the downstream Phase A2
    # loop to skip taxiway centerline / buffer processing entirely.
    new_dico["taxiway"] = (_shp_geom.Polygon(), [])

    # ── 3. Runway polygons: rectangles from apt.dat row 100 ─────────
    # Used by Phase A1 for clipping aprons and buildings off the
    # runway area.  Each Runway gives us two endpoints plus a width,
    # which we turn into a rectangle via runway_corners().
    rwy_polys = []
    for rwy in apt_data.runways:
        try:
            patch_width = rwy.width_m + 2 * RUNWAY_MARGIN
            corners = runway_corners(
                rwy.lat_a, rwy.lon_a,
                rwy.lat_b, rwy.lon_b,
                patch_width)
            if corners is None:
                continue
            # corners is [(lat, lon), ...] — flip to (x, y) and shift.
            poly = _shp_geom.Polygon([
                (lon - tile.lon, lat - tile.lat)
                for (lat, lon) in corners
            ])
            if not poly.is_valid:
                poly = poly.buffer(0)
            if not poly.is_empty and isinstance(poly, _shp_geom.Polygon):
                rwy_polys.append(poly)
        except Exception:
            pass
    if rwy_polys:
        try:
            rwy_union = _shp_ops.unary_union(rwy_polys)
            new_dico["runway"] = (rwy_union, [])
        except Exception:
            pass

    # ── 4. Boundary: prefer apt.dat row 130 when present ────────────
    if apt_data.boundary is not None and not apt_data.boundary.is_empty:
        bnd_tr = _to_tile_relative(apt_data.boundary)
        if bnd_tr is not None and not bnd_tr.is_empty:
            new_dico["boundary"] = bnd_tr

    return new_dico


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


def generate_patch_osm(icao, runway_pairs, runway_widths=None, tile=None):
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

    Returns:
        str: Complete OSM XML content for the patch file.
    """
    if runway_widths is None:
        runway_widths = {}
    node_id = -1
    way_id = -1
    nodes = []  # list of (id, lat, lon)
    ways = []  # list of (id, [node_ids], {tags})

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

    for desig_a, data_a, desig_b, data_b in runway_pairs:
        if desig_b is not None and data_b is not None:
            # ── Paired runway ────────────────────────────────────────────
            lat_a, lon_a = data_a["lat"], data_a["lon"]
            lat_b, lon_b = data_b["lat"], data_b["lon"]
            elev_a = data_a["elevation_m"]
            elev_b = data_b["elevation_m"]
            displaced_a = data_a["displaced_m"]
            displaced_b = data_b["displaced_m"]

            # Determine runway width for this pair (apt.dat → default)
            rwy_width = runway_widths.get(
                desig_a,
                runway_widths.get(desig_b, DEFAULT_RUNWAY_WIDTH),
            )
            patch_width = rwy_width + 2 * RUNWAY_MARGIN

            # Compute threshold-to-threshold distance for grade calculation
            mid_lat = (lat_a + lat_b) / 2.0
            cos_lat_v = cos(mid_lat * pi / 180.0)
            if cos_lat_v < 1e-6:
                cos_lat_v = 1e-6
            dx_m = (lon_b - lon_a) * cos_lat_v * DEG_TO_M
            dy_m = (lat_b - lat_a) * DEG_TO_M
            thresh_dist = sqrt(dx_m ** 2 + dy_m ** 2)
            if thresh_dist < 1.0:
                continue

            grade = (elev_b - elev_a) / thresh_dist

            # ── Physical runway ends ─────────────────────────────────────
            if displaced_a > 0:
                phys_end_a = extend_point(
                    lat_b, lon_b, lat_a, lon_a, displaced_a
                )
                elev_phys_a = elev_a - grade * displaced_a
            else:
                phys_end_a = (lat_a, lon_a)
                elev_phys_a = elev_a

            if displaced_b > 0:
                phys_end_b = extend_point(
                    lat_a, lon_a, lat_b, lon_b, displaced_b
                )
                elev_phys_b = elev_b + grade * displaced_b
            else:
                phys_end_b = (lat_b, lon_b)
                elev_phys_b = elev_b

            # Full physical runway length
            dx_phys = (phys_end_b[1] - phys_end_a[1]) * cos_lat_v * DEG_TO_M
            dy_phys = (phys_end_b[0] - phys_end_a[0]) * DEG_TO_M
            phys_dist = sqrt(dx_phys ** 2 + dy_phys ** 2)
            if phys_dist < 1.0:
                continue

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

            # ── Grade-limited smoothing ──────────────────────────────────
            # Two passes:
            #   1. Hard cap on per-segment longitudinal grade
            #      (MAX_RUNWAY_GRADE = 1.5 % per FAA AC 150/5300-13B
            #      and EASA CS-ADR-DSN for Code 3-4 / Cat C-E runways).
            #   2. Vertical-curve rate-of-change cap
            #      (MAX_GRADE_CHANGE_PER_M = 1/3000): the FAA rule
            #      "L_min ≥ 30 m per 1 % of grade change" means the
            #      change in grade between adjacent segments must be
            #      ≤ MAX_GRADE_CHANGE_PER_M × average_segment_length.
            #   3. Re-cap pass — the rate-of-change pass can push some
            #      segments past the absolute cap; tidy it up.
            #
            # CIFP-anchored samples (displaced thresholds + physical
            # ends) are immutable — they never move during relaxation.
            elevs = [s[2] for s in sample_pts]
            anchored = [s[3] for s in sample_pts]

            def _pass_hard_cap():
                for _it in range(GRADE_RELAX_ITERATIONS):
                    changed = False
                    for idx in range(len(elevs)):
                        if anchored[idx]:
                            continue
                        for nidx in (idx - 1, idx + 1):
                            if nidx < 0 or nidx >= len(elevs):
                                continue
                            seg_dist = (
                                abs(fractions[nidx] - fractions[idx])
                                * phys_dist)
                            if seg_dist < 0.1:
                                continue
                            max_rise = seg_dist * MAX_RUNWAY_GRADE
                            if elevs[idx] > elevs[nidx] + max_rise:
                                elevs[idx] = elevs[nidx] + max_rise
                                changed = True
                            elif elevs[idx] < elevs[nidx] - max_rise:
                                elevs[idx] = elevs[nidx] - max_rise
                                changed = True
                    if not changed:
                        return

            _pass_hard_cap()

            # FAA vertical-curve rate-of-change pass.  For each
            # interior sample i (1..n-2):
            #   g_left  = (e[i]   - e[i-1]) / L_left
            #   g_right = (e[i+1] - e[i])   / L_right
            #   |g_right - g_left| must be
            #     ≤ MAX_GRADE_CHANGE_PER_M × ((L_left + L_right) / 2)
            # If violated, push e[i] toward the value that exactly
            # meets the constraint.  Anchored samples cannot move.
            #
            # At the default 100 m segment length and 1.5 % grade cap,
            # this rule is automatically satisfied (max possible
            # |Δgrade| is 3 % × 100 m / 100 m = 0.03, vs allowed
            # 1/3000 × 100 = 0.0333), so this pass usually no-ops.
            # It exists so that any future tightening of the segment
            # length or grade cap stays compliant without code changes.
            def _seg_len(i):
                return abs(fractions[i + 1] - fractions[i]) * phys_dist

            for _it in range(GRADE_RELAX_ITERATIONS):
                changed = False
                for i in range(1, len(elevs) - 1):
                    if anchored[i]:
                        continue
                    ll = _seg_len(i - 1)
                    lr = _seg_len(i)
                    if ll < 0.1 or lr < 0.1:
                        continue
                    g_left = (elevs[i] - elevs[i - 1]) / ll
                    g_right = (elevs[i + 1] - elevs[i]) / lr
                    max_dg = MAX_GRADE_CHANGE_PER_M * ((ll + lr) / 2.0)
                    dg = g_right - g_left
                    if abs(dg) <= max_dg:
                        continue
                    # Solve elevs[i] s.t. dg = ±max_dg
                    target_dg = max_dg if dg > 0 else -max_dg
                    denom = 1.0 / lr + 1.0 / ll
                    new_e = (elevs[i + 1] / lr
                             + elevs[i - 1] / ll
                             - target_dg) / denom
                    if abs(new_e - elevs[i]) > 0.001:
                        elevs[i] = new_e
                        changed = True
                if not changed:
                    break

            # Re-apply the hard cap in case the rate-of-change pass
            # pushed any segment past the absolute grade limit.
            _pass_hard_cap()

            # ── Emit segmented rectangles ────────────────────────────────
            for idx in range(len(sample_pts) - 1):
                s_a = sample_pts[idx]
                s_b = sample_pts[idx + 1]
                add_rect_patch(
                    s_a[0], s_a[1], elevs[idx],
                    s_b[0], s_b[1], elevs[idx + 1],
                    patch_width,
                )

            # ── Flat overrun extensions beyond physical ends ─────────────
            ext_a = extend_point(
                phys_end_b[0], phys_end_b[1],
                phys_end_a[0], phys_end_a[1],
                OVERRUN_EXTENSION,
            )
            add_rect_patch(
                ext_a[0], ext_a[1], elevs[0],
                phys_end_a[0], phys_end_a[1], elevs[0],
                patch_width,
            )

            ext_b = extend_point(
                phys_end_a[0], phys_end_a[1],
                phys_end_b[0], phys_end_b[1],
                OVERRUN_EXTENSION,
            )
            add_rect_patch(
                phys_end_b[0], phys_end_b[1], elevs[-1],
                ext_b[0], ext_b[1], elevs[-1],
                patch_width,
            )

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
    return "\n".join(lines) + "\n"


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
                          road_data=None):
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
        # Start from the legacy per-segment runway patch; the surface
        # generator (Phase A-F) appends buildings, taxiways, aprons,
        # triangle wedges, and boundary band on top.
        osm_content = generate_patch_osm(
            icao, pairs, runway_widths=runway_widths, tile=tile
        )
        num_surface_patches = 0

        # Collect road data for this airport
        airport_roads = None
        if road_data:
            airport_roads = (
                road_data.get(icao)
                or road_data.get(icao.upper())
                or road_data.get(icao.lower())
            )

        # Surface generator requires building data to do anything useful
        # (taxiways alone don't trigger it).
        run_surface = has_dem and rwy_pairs_for_elev and airport_buildings
        if run_surface:
            (surface_lines, _) = generate_airport_surface_patches(
                icao, airport_taxiways or [], airport_buildings or [],
                rwy_pairs_for_elev, tile, dico_apt_entry,
                start_node_id=-10000,
                road_data=airport_roads,
                xplane_root=xplane_root_from_cifp_path(cifp_path),
            )
            if surface_lines:
                num_surface_patches = sum(
                    1 for l in surface_lines if "<way " in l
                )
                osm_content = osm_content.replace(
                    "</osm>",
                    "  <!-- Building/transition patches "
                    "({} triangles) -->\n".format(num_surface_patches)
                    + "\n".join(surface_lines) + "\n"
                    + "</osm>\n",
                )

        # Write the auto-patch file
        if not os.path.exists(patch_dir):
            os.makedirs(patch_dir)

        auto_patch_file = os.path.join(
            patch_dir, "{}_auto.patch.osm".format(icao)
        )
        try:
            with open(auto_patch_file, "w") as f:
                f.write(osm_content)
            parts = [icao, " ({} runway pairs".format(len(pairs))]
            if num_surface_patches:
                parts.append(
                    ", {} surface shapes".format(num_surface_patches)
                )
            parts.append(")")
            UI.vprint(1, "   Auto-patch: Generated", "".join(parts))
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
# CIFP Elevation Data for Taxiway Grade Anchoring
# ──────────────────────────────────────────────────────────────────────────────
def get_cifp_elevation_data(cifp_path, tile_lat, tile_lon):
    """Load CIFP runway elevation data for all airports within a tile.

    This provides authoritative threshold elevations that can be used to
    anchor taxiway elevation models instead of relying on noisy DEM data.

    For each airport with CIFP data, returns paired runway threshold
    coordinates and elevations in tile-relative coordinates (matching the
    coordinate system used by encode_runways_taxiways_and_aprons()).

    Args:
        cifp_path: Path to the CIFP data directory.
        tile_lat: Integer latitude of tile's SW corner.
        tile_lon: Integer longitude of tile's SW corner.

    Returns:
        dict: {icao: [{'start': (x, y), 'end': (x, y),
                        'elev_start': float, 'elev_end': float}, ...]}
              Coordinates are tile-relative (lon - tile_lon, lat - tile_lat).
              Elevations are in meters.
    """
    if not cifp_path or not os.path.isdir(cifp_path):
        return {}

    cifp_airports = discover_cifp_airports(cifp_path)
    result = {}

    for icao, filepath in cifp_airports.items():
        runways = parse_cifp_file(filepath)
        if not runways:
            continue
        if not airport_in_tile(runways, tile_lat, tile_lon):
            continue

        pairs = pair_runways(runways)
        if not pairs:
            continue

        runway_data = []
        for desig_a, data_a, desig_b, data_b in pairs:
            if desig_b is None or data_b is None:
                continue  # Need both thresholds for a slope

            # Tile-relative coordinates (lon, lat) → (x, y)
            start_x = data_a["lon"] - tile_lon
            start_y = data_a["lat"] - tile_lat
            end_x = data_b["lon"] - tile_lon
            end_y = data_b["lat"] - tile_lat

            runway_data.append({
                "start": (start_x, start_y),
                "end": (end_x, end_y),
                "elev_start": data_a["elevation_m"],
                "elev_end": data_b["elevation_m"],
                "desig_a": desig_a,
                "desig_b": desig_b,
            })

        if runway_data:
            result[icao] = runway_data

    return result


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
        dict: {airport_key: [{'centerline': [(lon, lat), ...], 'wayid': int}, ...]}
              Coordinates are ABSOLUTE (not tile-relative).
    """
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
                taxiways.append({
                    "centerline": coords,
                    "wayid": wayid,
                })
        if taxiways:
            result[airport] = taxiways
    return result


def compute_airport_bboxes(dico_airports, tile, buffer_m=1000.0):
    """Compute individual bounding boxes around each airport boundary.

    Returns one bbox per airport so that only buildings near each airport
    are downloaded, rather than filling the gaps between distant airports.

    Args:
        dico_airports: dict with processed airport data including 'boundary'.
        tile: Tile object with .lat and .lon.
        buffer_m: Buffer distance in meters around each airport boundary.

    Returns:
        list of (airport_key, (south, west, north, east)) tuples.
              Empty list if no airport boundaries are available.
    """
    results = []
    for airport in dico_airports:
        boundary = dico_airports[airport].get("boundary")
        if boundary is None or boundary.is_empty:
            continue
        bounds = boundary.bounds  # (minx, miny, maxx, maxy) tile-relative
        bmin_lon = bounds[0] + tile.lon
        bmin_lat = bounds[1] + tile.lat
        bmax_lon = bounds[2] + tile.lon
        bmax_lat = bounds[3] + tile.lat
        mid_lat = (bmin_lat + bmax_lat) / 2.0
        lat_buf = buffer_m / DEG_TO_M
        cos_lat = cos(mid_lat * pi / 180.0)
        if cos_lat < 1e-6:
            cos_lat = 1e-6
        lon_buf = buffer_m / (DEG_TO_M * cos_lat)
        results.append((
            airport,
            (bmin_lat - lat_buf, bmin_lon - lon_buf,
             bmax_lat + lat_buf, bmax_lon + lon_buf),
        ))
    return results


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


def _project_point_onto_runway(lat, lon, runway_pairs, max_dist=150.0):
    """Project a point onto the nearest CIFP runway and return interpolated elevation.

    Uses vector math to find the closest point on each runway centerline,
    then interpolates the CIFP threshold elevations at that position.

    Args:
        lat, lon: Point coordinates (absolute).
        runway_pairs: List of CIFP runway pair dicts with 'data_a' and 'data_b'
                      keys containing 'lat', 'lon', 'elevation_m'.
        max_dist: Maximum perpendicular distance in meters to consider a match.

    Returns:
        Interpolated elevation in meters, or None if no runway is within max_dist.
    """
    best_elev = None
    best_dist = max_dist

    for rwy in runway_pairs:
        r1_lat = rwy["data_a"]["lat"]
        r1_lon = rwy["data_a"]["lon"]
        r2_lat = rwy["data_b"]["lat"]
        r2_lon = rwy["data_b"]["lon"]
        r1_elev = rwy["data_a"]["elevation_m"]
        r2_elev = rwy["data_b"]["elevation_m"]

        mid_lat = (r1_lat + r2_lat) / 2.0
        cos_lat_factor = cos(mid_lat * pi / 180.0)

        bx = (r2_lon - r1_lon) * cos_lat_factor * DEG_TO_M
        by = (r2_lat - r1_lat) * DEG_TO_M
        px = (lon - r1_lon) * cos_lat_factor * DEG_TO_M
        py = (lat - r1_lat) * DEG_TO_M

        rwy_len_sq = bx * bx + by * by
        if rwy_len_sq < 1.0:
            continue
        t = (px * bx + py * by) / rwy_len_sq
        # Allow slight extension beyond runway ends (10% of length)
        t_clamped = max(-0.1, min(1.1, t))

        proj_x = t_clamped * bx
        proj_y = t_clamped * by
        dist = sqrt((px - proj_x) ** 2 + (py - proj_y) ** 2)

        if dist < best_dist:
            best_dist = dist
            t_elev = max(0.0, min(1.0, t))
            best_elev = r1_elev + t_elev * (r2_elev - r1_elev)

    return best_elev


def _compute_taxiway_elevations(centerline, runway_pairs, tile,
                                platform_anchors=None):
    """Compute elevation at each node of a taxiway centerline.

    The elevation model treats taxiways as *connectors* between different
    elevation zones.  The order of truth is:

    1. CIFP runway projections — highest confidence for nodes near runways.
    2. Platform anchors — building/terminal centroids sampled from DEM.
       Taxiway nodes near a platform inherit its elevation because the
       taxiway was built to reach that platform.
    3. DEM — fallback for mid-span nodes far from any anchor.
    4. Grade constraint — the whole chain is smoothed so no segment
       exceeds FAA max taxiway grade (1.5%).

    Args:
        centerline: list of (lon, lat) tuples - absolute coordinates.
        runway_pairs: list of CIFP runway pair dicts (from parse_cifp_file +
                      pair_runways) with lat/lon/elevation_m for each threshold.
        tile: Tile object with .dem, .lat, .lon.
        platform_anchors: optional list of (lat, lon, elevation) tuples for
                          building/terminal platforms that taxiways connect to.

    Returns:
        list of float: Elevation in meters at each centerline node.
    """
    n = len(centerline)
    elevations = [None] * n
    anchored = [False] * n

    # Step 1a: Anchor nodes near CIFP runways
    for i, (lon, lat) in enumerate(centerline):
        elev = _project_point_onto_runway(lat, lon, runway_pairs, max_dist=150.0)
        if elev is not None:
            elevations[i] = elev
            anchored[i] = True

    # Step 1b: Anchor nodes near building/terminal platforms.
    # A taxiway node within ~50m of a platform inherits the platform
    # elevation — the taxiway was graded to meet the building.
    # Only apply if the node wasn't already runway-anchored.
    PLATFORM_ANCHOR_DIST = 50.0  # meters
    if platform_anchors:
        for i, (lon, lat) in enumerate(centerline):
            if anchored[i]:
                continue
            mid_lat = lat
            cos_lat_factor = cos(mid_lat * pi / 180.0)
            best_dist = PLATFORM_ANCHOR_DIST
            best_elev = None
            for (p_lat, p_lon, p_elev) in platform_anchors:
                dx = (lon - p_lon) * cos_lat_factor * DEG_TO_M
                dy = (lat - p_lat) * DEG_TO_M
                dist = sqrt(dx * dx + dy * dy)
                if dist < best_dist:
                    best_dist = dist
                    best_elev = p_elev
            if best_elev is not None:
                elevations[i] = best_elev
                anchored[i] = True

    # Step 2: Fill remaining nodes from DEM
    for i, (lon, lat) in enumerate(centerline):
        if elevations[i] is None:
            x = lon - tile.lon
            y = lat - tile.lat
            try:
                elevations[i] = tile.dem.alt((x, y))
            except Exception:
                elevations[i] = 0.0

    # Step 3: Apply grade constraint — propagate from anchored points
    # outward.  Two passes (forward + backward) ensure the constraint
    # is enforced from both ends of any anchored segment.
    for i in range(1, n):
        if anchored[i]:
            continue
        lon0, lat0 = centerline[i - 1]
        lon1, lat1 = centerline[i]
        mid_lat = (lat0 + lat1) / 2.0
        cos_lat_factor = cos(mid_lat * pi / 180.0)
        seg_len = sqrt(
            ((lon1 - lon0) * cos_lat_factor * DEG_TO_M) ** 2
            + ((lat1 - lat0) * DEG_TO_M) ** 2
        )
        if seg_len < 0.1:
            continue
        max_rise = seg_len * MAX_TAXIWAY_GRADE
        if elevations[i] > elevations[i - 1] + max_rise:
            elevations[i] = elevations[i - 1] + max_rise
        elif elevations[i] < elevations[i - 1] - max_rise:
            elevations[i] = elevations[i - 1] - max_rise

    for i in range(n - 2, -1, -1):
        if anchored[i]:
            continue
        lon0, lat0 = centerline[i]
        lon1, lat1 = centerline[i + 1]
        mid_lat = (lat0 + lat1) / 2.0
        cos_lat_factor = cos(mid_lat * pi / 180.0)
        seg_len = sqrt(
            ((lon1 - lon0) * cos_lat_factor * DEG_TO_M) ** 2
            + ((lat1 - lat0) * DEG_TO_M) ** 2
        )
        if seg_len < 0.1:
            continue
        max_rise = seg_len * MAX_TAXIWAY_GRADE
        if elevations[i] > elevations[i + 1] + max_rise:
            elevations[i] = elevations[i + 1] + max_rise
        elif elevations[i] < elevations[i + 1] - max_rise:
            elevations[i] = elevations[i + 1] - max_rise

    return elevations


def _interpolate_elevation_along_centerline(point_m, centerline_m,
                                              elevations):
    """Interpolate elevation for a point by projecting onto a centerline.

    Projects the point onto the polyline and returns the linearly
    interpolated elevation at the projected position.

    Args:
        point_m: (x, y) tuple in meter-scaled coordinates.
        centerline_m: list of (x, y) tuples in meter-scaled coordinates.
        elevations: list of float elevations at each centerline node.

    Returns:
        float: Interpolated elevation.
    """
    px, py = point_m
    best_elev = elevations[0]
    best_dist_sq = float("inf")

    for i in range(len(centerline_m) - 1):
        ax, ay = centerline_m[i]
        bx, by = centerline_m[i + 1]
        dx, dy = bx - ax, by - ay
        seg_len_sq = dx * dx + dy * dy
        if seg_len_sq < 1e-6:
            continue
        t = ((px - ax) * dx + (py - ay) * dy) / seg_len_sq
        t = max(0.0, min(1.0, t))
        proj_x = ax + t * dx
        proj_y = ay + t * dy
        dist_sq = (px - proj_x) ** 2 + (py - proj_y) ** 2
        if dist_sq < best_dist_sq:
            best_dist_sq = dist_sq
            best_elev = elevations[i] + t * (elevations[i + 1] - elevations[i])

    return best_elev


def generate_taxiway_patches_osm(icao, taxiway_data, runway_pairs, tile,
                                  runway_widths=None, start_node_id=-10000):
    """DEPRECATED: Use generate_airport_surface_patches() instead.

    Legacy function that generates taxiway patches without platform anchors.
    For each taxiway, creates a single buffered polygon around the centerline
    with per-node elevations interpolated from the CIFP-anchored and
    grade-constrained centerline profile. Uses the 'node_altitudes' tag for
    per-vertex elevation specification.

    Args:
        icao: Airport ICAO code.
        taxiway_data: List of dicts with 'centerline' key from extract_taxiway_info.
        runway_pairs: List of runway pair dicts with threshold data.
        tile: Tile object with .dem, .lat, .lon.
        runway_widths: Optional dict of {designator: width_m} from apt.dat.
        start_node_id: Starting negative node ID (to avoid collisions with
                       runway patch nodes).

    Returns:
        tuple: (osm_lines, next_node_id) - list of OSM XML line strings and
               the next available node ID.
    """
    from shapely import geometry as shp_geom

    lines = []
    nid = start_node_id
    wid = start_node_id

    half_width = (DEFAULT_TAXIWAY_WIDTH + 2 * RUNWAY_MARGIN) / 2.0

    for twy in taxiway_data:
        centerline = twy["centerline"]
        if len(centerline) < 2:
            continue

        # Compute elevations at all centerline nodes
        elevations = _compute_taxiway_elevations(
            centerline, runway_pairs, tile
        )

        # Convert centerline to meter-scaled coords for accurate buffering
        # Use the midpoint latitude for the cosine correction
        mid_lat = sum(c[1] for c in centerline) / len(centerline)
        cos_lat_factor = cos(mid_lat * pi / 180.0)
        if cos_lat_factor < 1e-6:
            cos_lat_factor = 1e-6

        centerline_m = [
            (lon * cos_lat_factor * DEG_TO_M, lat * DEG_TO_M)
            for lon, lat in centerline
        ]

        # Build LineString and buffer to get the taxiway outline
        try:
            ls = shp_geom.LineString(centerline_m)
            buffered = ls.buffer(
                half_width,
                cap_style=2,   # flat caps at endpoints
                join_style=2,  # mitre joins at bends
                mitre_limit=2.0,
            )
        except Exception:
            continue

        if buffered.is_empty or not buffered.is_valid:
            continue

        # Extract exterior ring of the buffered polygon
        # (for MultiPolygon, take the largest)
        if buffered.geom_type == "MultiPolygon":
            buffered = max(buffered.geoms, key=lambda g: g.area)
        ring_m = list(buffered.exterior.coords)

        if len(ring_m) < 4:
            continue

        # Convert buffered polygon nodes back to lon/lat and compute
        # per-node elevation by projecting each node onto the centerline
        node_coords = []
        node_elevations = []
        for mx, my in ring_m:
            blon = mx / (cos_lat_factor * DEG_TO_M)
            blat = my / DEG_TO_M
            elev = _interpolate_elevation_along_centerline(
                (mx, my), centerline_m, elevations
            )
            node_coords.append((blat, blon))
            node_elevations.append(elev)

        # Write OSM nodes
        node_ids = []
        for (nlat, nlon) in node_coords:
            lines.append(
                "  <node id='{}' action='modify' visible='true'"
                " lat='{:.11f}' lon='{:.11f}' />".format(nid, nlat, nlon)
            )
            node_ids.append(nid)
            nid -= 1

        # Write way with per-node altitude tag
        alt_str = ",".join("{:.2f}".format(e) for e in node_elevations)
        lines.append(
            "  <way id='{}' action='modify' visible='true'>".format(wid)
        )
        for nd_id in node_ids:
            lines.append("    <nd ref='{}' />".format(nd_id))
        lines.append("    <tag k='node_altitudes' v='{}' />".format(alt_str))
        lines.append("  </way>")
        wid -= 1

    return (lines, nid)


# ──────────────────────────────────────────────────────────────────────────────
# Unified Airport Surface Generator
# ──────────────────────────────────────────────────────────────────────────────
# FAA AC 150/5300-13B grade limits for Approach Category C-E airports:
#   Runway longitudinal:  max 1.5%
#   Taxiway longitudinal: max 1.5%
#   Apron any direction:  max 1.0%
#   Rate of grade change: max 1% per 30m vertical curve
MAX_TAXIWAY_GRADE = 0.015     # 1.5% max longitudinal grade for taxiways
MAX_APRON_GRADE = 0.010       # 1.0% max grade in any direction for aprons
MAX_SURFACE_GRADE = 0.015     # 1.5% general transition grade limit
MAX_RUNWAY_GRADE = 0.015      # 1.5% max longitudinal grade for runways (Cat C-E)

# FAA AC 150/5300-13B vertical-curve rule for runways and taxiways:
#   L_min = 30 m × (Δgrade as % per 1%) = (Δgrade × 100) × 30 m
# i.e. a 1.0% grade change requires at least 30 m of vertical curve, so
# the maximum allowable change in grade per metre of pavement is 1/3000.
# Used as a second relaxation pass after the hard MAX_RUNWAY_GRADE cap
# (see generate_patch_osm).
MAX_GRADE_CHANGE_PER_M = 1.0 / 3000.0  # ≈ 0.0333% per metre

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
# Rate-of-grade-change: FAA AC 150/5300-13B requires max 1% change per 30m
# vertical curve.  A grade jump from +1.5% to -1.5% (3% change) needs 90m.
MAX_GRADE_CHANGE_PER_M = 0.01 / 30.0  # 1% per 30m ≈ 0.000333 per meter
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


def generate_airport_surface_patches(icao, taxiway_data, building_data,
                                      runway_pairs, tile, dico_apt_entry,
                                      start_node_id=-10000, road_data=None,
                                      xplane_root=None):
    """Generate efficient airport surface patches using JOSM-friendly shapes.

    Uses the simplest shape that achieves smooth, grade-limited slopes.
    All shapes use tags that are easy to inspect and edit in JOSM:
      - altitude          → flat polygon (single elevation)
      - altitude_high/low → sloped rectangle (slope in one direction)
      - node_altitudes    → triangle (3 vertices, compound slope)

    Shape hierarchy (prefer simpler over complex):
    1. Taxiways  → segmented sloped rectangles along each centerline,
                   exactly like runways.  Flat segments use altitude,
                   sloped segments use altitude_high/altitude_low.
    2. Aprons    → flat altitude=N polygons (most ramp areas are level).
    3. Buildings → flat altitude=N polygons.
    4. Triangles → ONLY where needed for compound slopes: complex
                   building transition zones on sloping terrain,
                   or junction areas between different surfaces.
    5. Blend rects → sloped rects from airport edge to surrounding DEM.
    6. Drainage  → triangulated ditches in infield areas, sloping from
                   pavement edge to center at natural drainage depth.

    The emitted-geometry accumulator (all_emitted_parts_m) ensures no
    phase emits shapes overlapping previously emitted ones.

    Args:
        icao: Airport ICAO code.
        taxiway_data: List of dicts with 'centerline' key.  Ignored
            when apt.dat data is available (see Phase A0.5).
        building_data: List of dicts with 'footprint' key.
        runway_pairs: List of CIFP runway pair dicts.
        tile: Tile object with .dem, .lat, .lon.
        dico_apt_entry: The dico_airports[airport] dict.
        start_node_id: Starting negative node ID.
        road_data: Optional list of road dicts from extract_road_info().
        xplane_root: X-Plane install root.  If set, Phase A0.5 searches
            ``<xplane_root>/Custom Scenery`` (per-airport packs first)
            and falls back to Global Airports / default scenery to
            locate the airport's apt.dat file.  When found, the
            apt.dat pavement polygons REPLACE the OSM-derived
            taxiway/apron/runway shapes in dico_apt_entry.  apt.dat
            polygons are disjoint by construction, eliminating the
            buffering / precision-drift artefacts that cause
            flat-flat and flat-triangle overlap.

    Returns:
        tuple: (osm_lines, next_node_id)
    """
    from shapely import geometry as shp_geom
    from shapely import ops as shp_ops

    # ════════════════════════════════════════════════════════════════
    # PHASE A0.5: Prefer apt.dat pavement data over OSM when available
    # ════════════════════════════════════════════════════════════════
    # apt.dat polygons are the authoritative pavement geometry for
    # what X-Plane actually renders as the ground texture.  When we
    # can find an apt.dat for this airport (per-airport Custom Scenery
    # pack first, then Global Airports, then default), we load it and
    # replace the OSM-derived shapes in dico_apt_entry.
    #
    # For commit 8 every apt.dat pavement goes into dico_apt_entry
    # ["apron"] and taxiway_data is cleared.  Name-based classification
    # (TWY / RAMP / etc.) is a follow-up commit.
    #
    # If no apt.dat is found, the OSM path runs unchanged.
    apt_dat_used = False
    if xplane_root:
        try:
            import O4_Apt_Dat_Reader as _APR
            aptdat_path = _APR.find_airport_apt_dat(xplane_root, icao)
            if aptdat_path:
                apt_data = _APR.load_airport(aptdat_path, icao)
                if apt_data is not None and apt_data.pavements:
                    dico_apt_entry = _dico_from_apt_dat(
                        apt_data, tile, dico_apt_entry)
                    taxiway_data = []    # apt.dat ⇒ everything is apron
                    apt_dat_used = True
                    UI.vprint(2,
                        "    {}: apt.dat loaded from {} — "
                        "{} pavements, {} runways".format(
                            icao, aptdat_path,
                            len(apt_data.pavements),
                            len(apt_data.runways)))
        except Exception as e:
            UI.vprint(2,
                "    {}: apt.dat load failed ({}), "
                "falling back to OSM".format(icao, e))
    if not apt_dat_used:
        UI.vprint(2, "    {}: using OSM aerodrome data".format(icao))

    ZONE_BUILDING = 1
    ZONE_APRON = 2
    ZONE_SURFACE = 3

    lines = []
    nid = start_node_id
    wid = start_node_id

    # ── Helpers ──────────────────────────────────────────────────────
    def _add_node(lat, lon):
        nonlocal nid
        lines.append(
            "  <node id='{}' action='modify' visible='true'"
            " lat='{:.11f}' lon='{:.11f}' />".format(nid, lat, lon)
        )
        cur = nid
        nid -= 1
        return cur

    def _add_way(node_ids, tags):
        nonlocal wid
        lines.append(
            "  <way id='{}' action='modify' visible='true'>".format(wid)
        )
        for nd_id in node_ids:
            lines.append("    <nd ref='{}' />".format(nd_id))
        for k, v in tags.items():
            lines.append("    <tag k='{}' v='{}' />".format(k, v))
        lines.append("  </way>")
        wid -= 1

    def _validate_ring(coords_ll):
        """Validate a coordinate ring.  Returns a list of valid Polygon
        coordinate rings, or an empty list if the geometry is degenerate.
        Repairs self-intersections by buffering with zero width.
        Input/output coords are (lon, lat) tuples."""
        if len(coords_ll) < 4:
            return []
        cl = list(coords_ll)
        if cl[0] != cl[-1]:
            cl.append(cl[0])
        try:
            p = shp_geom.Polygon(cl)
            if p.is_empty:
                return []
            if not p.is_valid:
                p = p.buffer(0)
            if p.is_empty:
                return []
            # buffer(0) can produce MultiPolygon
            if p.geom_type == "Polygon":
                return [list(p.exterior.coords)]
            elif hasattr(p, "geoms"):
                result = []
                for g in p.geoms:
                    if (g.geom_type == "Polygon"
                            and not g.is_empty
                            and g.area > 0):
                        result.append(list(g.exterior.coords))
                return result
        except Exception:
            pass
        return []

    def _emit_flat_poly(coords_ll, elev):
        """Emit a flat altitude=N polygon from (lon, lat) coords.
        Validates geometry before emission."""
        for ring in _validate_ring(coords_ll):
            if len(ring) < 4:
                continue
            node_ids = [_add_node(lat, lon) for (lon, lat) in ring]
            _add_way(node_ids, {"altitude": "{:.1f}".format(elev)})

    def _emit_sloped_rect(lat_a, lon_a, elev_a, lat_b, lon_b, elev_b,
                          width):
        """Emit a sloped or flat rectangle, exactly like a runway segment.

        Uses altitude_high/altitude_low when elevations differ, otherwise
        flat altitude.  High-elevation end is always at corners 0-1.
        """
        if abs(elev_a - elev_b) >= 0.1:
            if elev_a >= elev_b:
                corners = runway_corners(lat_a, lon_a, lat_b, lon_b, width)
                eh, el = elev_a, elev_b
            else:
                corners = runway_corners(lat_b, lon_b, lat_a, lon_a, width)
                eh, el = elev_b, elev_a
            if corners is None:
                return
            n0 = _add_node(corners[0][0], corners[0][1])
            n1 = _add_node(corners[1][0], corners[1][1])
            n2 = _add_node(corners[2][0], corners[2][1])
            n3 = _add_node(corners[3][0], corners[3][1])
            _add_way([n0, n1, n2, n3, n0], {
                "altitude_high": "{:.1f}".format(eh),
                "altitude_low": "{:.1f}".format(el),
                "cell_size": str(int(DEFAULT_CELL_SIZE)),
                "profile": DEFAULT_PROFILE,
            })
        else:
            corners = runway_corners(lat_a, lon_a, lat_b, lon_b, width)
            if corners is None:
                return
            n0 = _add_node(corners[0][0], corners[0][1])
            n1 = _add_node(corners[1][0], corners[1][1])
            n2 = _add_node(corners[2][0], corners[2][1])
            n3 = _add_node(corners[3][0], corners[3][1])
            _add_way([n0, n1, n2, n3, n0], {
                "altitude": "{:.1f}".format((elev_a + elev_b) / 2.0),
            })

    def _emit_triangle(lat0, lon0, e0, lat1, lon1, e1, lat2, lon2, e2):
        """Emit a single triangle with node_altitudes."""
        n0 = _add_node(lat0, lon0)
        n1 = _add_node(lat1, lon1)
        n2 = _add_node(lat2, lon2)
        alt_str = "{:.1f},{:.1f},{:.1f},{:.1f}".format(e0, e1, e2, e0)
        _add_way([n0, n1, n2, n0], {"node_altitudes": alt_str})

    # ── Mid-latitude for meter scaling ──────────────────────────────
    all_lats = []
    for bldg in (building_data or []):
        for lon, lat in bldg["footprint"]:
            all_lats.append(lat)
    for twy in (taxiway_data or []):
        for lon, lat in twy["centerline"]:
            all_lats.append(lat)
    for rp in (runway_pairs or []):
        for side in ("data_a", "data_b"):
            d = rp.get(side)
            if d:
                all_lats.append(d["lat"])
    if not all_lats:
        return (lines, nid)
    mid_lat = sum(all_lats) / len(all_lats)
    cos_lat = cos(mid_lat * pi / 180.0)
    if cos_lat < 1e-6:
        cos_lat = 1e-6

    def to_m(lon, lat):
        return (lon * cos_lat * DEG_TO_M, lat * DEG_TO_M)

    def to_ll(mx, my):
        return (mx / (cos_lat * DEG_TO_M), my / DEG_TO_M)

    # Full taxiway width (both sides of centerline)
    TWY_FULL_WIDTH = TAXIWAY_BUFFER_WIDTH * 2.0

    # ── Emitted-geometry accumulator ──────────────────────────────
    # Tracks ALL shapes emitted so far.  Each phase subtracts this
    # before emission to guarantee zero overlap across all phases.
    # Updated incrementally — after each batch of emissions, new
    # shapes are appended and union recomputed.
    all_emitted_parts_m = []  # list of Shapely geoms

    def _emitted_union():
        """Return the union of all emitted shapes so far.

        Defensively handles invalid geometries: if a straight
        unary_union fails (TopologicalError on a self-intersecting
        polygon, etc.), we run buffer(0) on each part to clean it
        up and retry.  Returning an empty Polygon on failure means
        Phase E and Phase F downstream would emit *everywhere* in
        the airport footprint, defeating the no-overlap invariant.
        """
        if not all_emitted_parts_m:
            return shp_geom.Polygon()
        try:
            return shp_ops.unary_union(all_emitted_parts_m)
        except Exception:
            pass
        # Retry with cleaned polygons.
        cleaned = []
        for p in all_emitted_parts_m:
            try:
                if p.is_empty:
                    continue
                if p.is_valid:
                    cleaned.append(p)
                else:
                    fixed = p.buffer(0)
                    if fixed.is_valid and not fixed.is_empty:
                        cleaned.append(fixed)
            except Exception:
                pass
        if not cleaned:
            return shp_geom.Polygon()
        try:
            return shp_ops.unary_union(cleaned)
        except Exception:
            # As a last resort, return a unary_union of just the
            # bounds of every part.  Coarse but never empty.
            try:
                from shapely.geometry import box
                bboxes = []
                for p in cleaned:
                    if not p.is_empty:
                        b = p.bounds
                        bboxes.append(box(*b))
                return shp_ops.unary_union(bboxes)
            except Exception:
                return shp_geom.Polygon()

    # ── Road geometry for terrain anchoring ───────────────────────
    road_lines_m = []    # [(LineString_m, highway_type, tunnel, bridge), ...]
    road_anchors_m = []  # [(x_m, y_m, elevation), ...] for IDW
    if road_data:
        for rd in road_data:
            cl = rd["centerline"]
            if len(cl) < 2:
                continue
            cl_m = [to_m(lon, lat) for lon, lat in cl]
            try:
                ls = shp_geom.LineString(cl_m)
                if not ls.is_empty:
                    road_lines_m.append((
                        ls, rd["highway_type"],
                        rd.get("tunnel", False),
                        rd.get("bridge", False)))
                    # Sample DEM at road nodes as terrain anchors
                    # (skip tunnels — they don't represent surface elevation)
                    if not rd.get("tunnel", False):
                        for lon, lat in cl:
                            try:
                                dem_e = tile.dem.alt(
                                    (lon - tile.lon, lat - tile.lat))
                                mx, my = to_m(lon, lat)
                                road_anchors_m.append((mx, my, dem_e))
                            except Exception:
                                pass
            except Exception:
                pass
    if road_lines_m:
        UI.vprint(2, "    {}: {} roads ({} surface anchors)".format(
            icao, len(road_lines_m), len(road_anchors_m)))

    # ════════════════════════════════════════════════════════════════
    # PHASE A: Build all geometries in meter space (no emission yet)
    # ════════════════════════════════════════════════════════════════

    # A0: Airport boundary from OSM data
    # The OSM aerodrome boundary is the authoritative perimeter of the
    # prepared surface area.  Inside it, terrain has been graded (cut or
    # filled) and is relatively smooth.  Outside it, natural terrain can
    # drop or rise steeply.  Used by Phase E (blend rects) and Phase F
    # (drainage) as the airport footprint.
    osm_boundary_m = shp_geom.Polygon()
    boundary_geom = dico_apt_entry.get("boundary")
    if boundary_geom is not None:
        try:
            # boundary can be Polygon, MultiPolygon, or tuple
            if isinstance(boundary_geom, tuple):
                boundary_geom = boundary_geom[0]
            bnd_parts = (boundary_geom.geoms
                         if hasattr(boundary_geom, "geoms")
                         else [boundary_geom])
            bnd_polys_m = []
            for bp in bnd_parts:
                if bp.is_empty or not bp.is_valid:
                    continue
                coords_m = [
                    to_m(x + tile.lon, y + tile.lat)
                    for x, y in bp.exterior.coords
                ]
                try:
                    p_m = shp_geom.Polygon(coords_m)
                    if p_m.is_valid and not p_m.is_empty:
                        bnd_polys_m.append(p_m)
                except Exception:
                    pass
            if bnd_polys_m:
                osm_boundary_m = shp_ops.unary_union(bnd_polys_m)
        except Exception:
            pass

    if not osm_boundary_m.is_empty:
        UI.vprint(2, "    {}: OSM airport boundary: {:.0f} m²".format(
            icao, osm_boundary_m.area))
    else:
        UI.vprint(2, "    {}: No OSM airport boundary available".format(icao))

    # A1: Runway polygons
    runway_polys_m = []
    runway_data_geom = dico_apt_entry.get("runway")
    if runway_data_geom and isinstance(runway_data_geom, (list, tuple)):
        rwy_geom = (runway_data_geom[0]
                     if len(runway_data_geom) > 0 else None)
        if rwy_geom is not None:
            rwy_parts = (rwy_geom.geoms if hasattr(rwy_geom, "geoms")
                         else [rwy_geom])
            for poly in rwy_parts:
                if poly.is_empty or not poly.is_valid:
                    continue
                coords_m = [
                    to_m(x + tile.lon, y + tile.lat)
                    for x, y in poly.exterior.coords
                ]
                try:
                    p = shp_geom.Polygon(coords_m)
                    if p.is_valid and not p.is_empty:
                        runway_polys_m.append(p)
                except Exception:
                    pass

    rwy_union_m = (shp_ops.unary_union(runway_polys_m)
                   if runway_polys_m else shp_geom.Polygon())

    # A2: Taxiway geometry — build buffer polygons and meter-space
    # centerlines.  Elevations are computed AFTER building platforms
    # are established (Step A2b below) so taxiways can anchor to both
    # CIFP runways AND building platforms.
    twy_buffers_m = []   # individual shapely polygons
    twy_centerlines = []  # [(centerline, elevations), ...] — filled in A2b
    twy_centerlines_m = []
    twy_raw_centerlines = []  # raw (lon, lat) lists for deferred elevation

    for twy in (taxiway_data or []):
        centerline = twy["centerline"]
        if len(centerline) < 2:
            twy_raw_centerlines.append(None)
            twy_centerlines_m.append(None)
            twy_buffers_m.append(shp_geom.Polygon())
            continue

        twy_raw_centerlines.append(centerline)

        cl_m = [to_m(lon, lat) for lon, lat in centerline]
        twy_centerlines_m.append(cl_m)

        try:
            ls = shp_geom.LineString(cl_m)
            buf = ls.buffer(
                TAXIWAY_BUFFER_WIDTH,
                cap_style=2, join_style=2, mitre_limit=2.0,
            )
            if buf.geom_type == "MultiPolygon":
                buf = max(buf.geoms, key=lambda g: g.area)
            if buf.is_valid and not buf.is_empty:
                twy_buffers_m.append(buf)
            else:
                twy_buffers_m.append(shp_geom.Polygon())
        except Exception:
            twy_buffers_m.append(shp_geom.Polygon())

    twy_union_m = (shp_ops.unary_union(
        [b for b in twy_buffers_m if not b.is_empty])
        if twy_buffers_m else shp_geom.Polygon())

    # (CIFP surface model and taxiway elevations are computed in
    #  Phase A5 below, after building platforms are established.)

    # A3: Apron polygons
    # When apt.dat polygons are the source (Phase A0.5 replaced
    # dico_apt_entry["apron"]), we take the polygons as-is: they're
    # already the authoritative shape X-Plane renders.  In
    # particular we do NOT apply APRON_BUFFER — that was an
    # ad-hoc outward inflation to patch over OSM geometry quality,
    # and it massively inflates apt.dat output because the pavements
    # are already clean.  We also preserve interior rings (holes)
    # so any carved-out regions in the apt.dat polygon survive into
    # meter space.
    #
    # For OSM-derived data the legacy behaviour is preserved: 5 m
    # outward buffer, exterior-only.
    apron_polys_m = []
    apron_data = dico_apt_entry.get("apron")
    if apron_data is not None:
        apron_geom = (apron_data[0] if isinstance(apron_data, tuple)
                      else apron_data)
        raw_polys = []
        if hasattr(apron_geom, "geoms"):
            raw_polys = list(apron_geom.geoms)
        elif hasattr(apron_geom, "exterior"):
            raw_polys = [apron_geom]
        for poly in raw_polys:
            if poly.is_empty or not poly.is_valid:
                continue
            try:
                exterior_m = [
                    to_m(x + tile.lon, y + tile.lat)
                    for x, y in poly.exterior.coords
                ]
                holes_m = []
                for ring in poly.interiors:
                    holes_m.append([
                        to_m(x + tile.lon, y + tile.lat)
                        for x, y in ring.coords
                    ])
                p_m = shp_geom.Polygon(
                    exterior_m, holes_m if holes_m else None)
                if not p_m.is_valid:
                    p_m = p_m.buffer(0)
                if p_m.is_empty:
                    continue
                if apt_dat_used:
                    # Trust the source — don't inflate or simplify.
                    apron_polys_m.append(p_m)
                else:
                    # Legacy OSM path: 5 m outward buffer.
                    apron_polys_m.append(p_m.buffer(APRON_BUFFER))
            except Exception:
                pass

    try:
        apron_union_m = shp_ops.unary_union(apron_polys_m) if apron_polys_m \
            else shp_geom.Polygon()
    except Exception:
        apron_union_m = shp_geom.Polygon()

    # A4: Building pads
    BLDG_PAD = 5.0
    BLDG_SIMPLIFY = 8.0
    building_polys_m = []

    for bldg in (building_data or []):
        fp = bldg["footprint"]
        if len(fp) < 4:
            continue
        coords_m = [to_m(lon, lat) for lon, lat in fp]
        try:
            p = shp_geom.Polygon(coords_m)
            if p.is_valid and not p.is_empty and p.area > 20.0:
                padded = p.buffer(BLDG_PAD, join_style=2, mitre_limit=2.0)
                simple = padded.simplify(BLDG_SIMPLIFY, preserve_topology=True)
                if simple.is_valid and not simple.is_empty:
                    building_polys_m.append(simple)
                else:
                    building_polys_m.append(padded)
        except Exception:
            pass

    bldg_union_m = (shp_ops.unary_union(building_polys_m)
                    if building_polys_m else shp_geom.Polygon())

    # A4b: Merge nearby buildings — absorb gate protuberances into
    # the main terminal shape.  At SPJC and similar airports, the
    # terminal has small protuberances for gates that create tiny
    # shapes overlapping the surrounding apron area.  By merging
    # them into the parent terminal, the result is one clean shape
    # with real concavities preserved.
    #
    # Critical invariants (do NOT regress):
    #   * Merge zone is measured from the *original* terminal footprint
    #     (term_buf is computed once before the inner loop) — never from
    #     the growing union, otherwise each absorption expands reach and
    #     the merge snowballs across the airfield.
    #   * The merged result is the *unary union* (concavities and all),
    #     NEVER the convex hull, which produces a giant blob that
    #     swallows aprons and unrelated buildings inside the hull.
    BLDG_MERGE_DIST = 8.0   # meters — gate protuberances within this
                            # radius of the ORIGINAL terminal footprint
                            # are absorbed.  Tuned for SPJC; smaller
                            # than the legacy 15 m default to prevent
                            # over-aggressive merging.
    if len(building_polys_m) > 1:
        bldg_areas = [(bp.area, i) for i, bp in enumerate(building_polys_m)
                       if not bp.is_empty]
        bldg_areas.sort(reverse=True)
        absorbed = set()
        for _, term_idx in bldg_areas:
            if term_idx in absorbed:
                continue
            original_terminal = building_polys_m[term_idx]
            if original_terminal.area < 500.0:
                continue
            # Reach measured from the ORIGINAL footprint, never the
            # growing union — single-pass, no snowball.
            term_buf = original_terminal.buffer(BLDG_MERGE_DIST)
            to_merge = [original_terminal]
            for _, other_idx in bldg_areas:
                if other_idx == term_idx or other_idx in absorbed:
                    continue
                other = building_polys_m[other_idx]
                if term_buf.intersects(other):
                    to_merge.append(other)
                    absorbed.add(other_idx)
            if len(to_merge) > 1:
                # unary_union keeps the actual outline (concavities
                # intact); simplify trims OSM micro-jitter.  No convex
                # hull anywhere — that would swallow apron area.
                merged = shp_ops.unary_union(to_merge)
                # Buildings within BLDG_MERGE_DIST of the original
                # terminal may not touch each other directly — the
                # union is then a MultiPolygon.  Bridge those sub-
                # BLDG_MERGE_DIST gaps with a morphological closing so
                # the result is one connected outline.  Buffer/unbuffer
                # by half the merge distance + small epsilon: any gap
                # of ≤ BLDG_MERGE_DIST gets filled, outer boundary
                # returns to its original shape (within shapely
                # precision).
                if hasattr(merged, "geoms"):
                    eps = BLDG_MERGE_DIST / 2.0 + 0.5
                    merged = merged.buffer(
                        eps, join_style=2, mitre_limit=2.0).buffer(
                        -eps, join_style=2, mitre_limit=2.0)
                # If still multi (very wide separation), take the
                # largest component containing the original terminal.
                if hasattr(merged, "geoms"):
                    merged = max(
                        merged.geoms, key=lambda g: g.area)
                # Smooth out small-scale features (< 10 m wide) via
                # a morphological opening: shrink by 5 m, grow by
                # 5 m.  The opening removes anything narrower than
                # 10 m — gate protuberances, service jetty bumps,
                # OSM tracing noise — and leaves the real terminal
                # outline.  Preserves concavities larger than 10 m
                # wide (real cut-outs stay intact).
                TERMINAL_OPEN_R = 5.0   # half of "feature size" 10 m
                try:
                    opened = merged.buffer(
                        -TERMINAL_OPEN_R,
                        join_style=2, mitre_limit=2.0)
                    if opened.is_empty:
                        # Building is so thin that the shrink fully
                        # consumed it — keep the original merged.
                        pass
                    else:
                        opened = opened.buffer(
                            TERMINAL_OPEN_R,
                            join_style=2, mitre_limit=2.0)
                        if hasattr(opened, "geoms"):
                            opened = max(
                                opened.geoms, key=lambda g: g.area)
                        if (not opened.is_empty
                                and hasattr(opened, "exterior")
                                and opened.area > 0.5 * merged.area):
                            # Accept the opened outline only if it
                            # kept at least half the original area
                            # (otherwise the building was mostly
                            # thin fingers and we'd lose it).
                            merged = opened
                except Exception:
                    pass
                merged = merged.simplify(
                    BLDG_SIMPLIFY, preserve_topology=True)
                if (merged.is_valid and not merged.is_empty
                        and hasattr(merged, "exterior")):
                    building_polys_m[term_idx] = merged

        if absorbed:
            new_polys = []
            for i, bp in enumerate(building_polys_m):
                if i not in absorbed:
                    new_polys.append(bp)
            building_polys_m = new_polys
            UI.vprint(2, "    {}: merged {} gate/building protuberances"
                      " into terminal shapes".format(icao, len(absorbed)))
            bldg_union_m = (shp_ops.unary_union(building_polys_m)
                            if building_polys_m else shp_geom.Polygon())

    # ════════════════════════════════════════════════════════════════
    # A5: Elevation Pipeline — anchors → platforms → connectors
    # ════════════════════════════════════════════════════════════════
    # The elevation model is built in this order:
    #
    #   1. CIFP RUNWAY THRESHOLDS — immutable truth from aeronautical data.
    #   2. BUILDING PLATFORMS — buildings sit on natural grade; DEM is most
    #      reliable at surveyed structure locations.  These become "platform
    #      anchors" at their centroid elevation.  Validated against nearest
    #      CIFP anchor with a sanity check.
    #   3. TAXIWAY ELEVATIONS — taxiways are connectors between zones.
    #      Nodes near runways anchor to CIFP; nodes near buildings anchor
    #      to the platform elevation; interior nodes use DEM then grade-
    #      constrain the whole chain at FAA 1.5% max.
    #   4. CIFP SURFACE MODEL — IDW interpolation from all anchors
    #      (thresholds + platforms + CIFP-projected taxiway nodes) for
    #      querying elevation at any airport point.
    #
    # This ensures taxiways correctly *negotiate* between the runway
    # elevation and the building/terminal elevation, rather than blindly
    # projecting the runway slope outward.

    # ── Step 1: CIFP runway threshold + physical end anchors ────────
    cifp_anchors_m = []  # [(x_m, y_m, elevation), ...]
    for rp in (runway_pairs or []):
        da = rp.get("data_a")
        db = rp.get("data_b")
        if not da or not db:
            continue
        for d in (da, db):
            if "lat" in d and "lon" in d and "elevation_m" in d:
                mx, my = to_m(d["lon"], d["lat"])
                cifp_anchors_m.append((mx, my, d["elevation_m"]))
        # Also add physical runway ends (past displaced thresholds)
        # so the surface model covers the full pavement extent.
        if "lat" in da and "lat" in db:
            lat_a, lon_a = da["lat"], da["lon"]
            lat_b, lon_b = db["lat"], db["lon"]
            elev_a, elev_b = da["elevation_m"], db["elevation_m"]
            mid_lat_rwy = (lat_a + lat_b) / 2.0
            cos_rwy = cos(mid_lat_rwy * pi / 180.0)
            if cos_rwy < 1e-6:
                cos_rwy = 1e-6
            dx_r = (lon_b - lon_a) * cos_rwy * DEG_TO_M
            dy_r = (lat_b - lat_a) * DEG_TO_M
            td = sqrt(dx_r ** 2 + dy_r ** 2)
            if td > 1.0:
                grade = (elev_b - elev_a) / td
                disp_a = da.get("displaced_m", 0)
                disp_b = db.get("displaced_m", 0)
                if disp_a > 0:
                    phys_a = extend_point(lat_b, lon_b, lat_a, lon_a, disp_a)
                    elev_pa = elev_a - grade * disp_a
                    mx, my = to_m(phys_a[1], phys_a[0])
                    cifp_anchors_m.append((mx, my, elev_pa))
                if disp_b > 0:
                    phys_b = extend_point(lat_a, lon_a, lat_b, lon_b, disp_b)
                    elev_pb = elev_b + grade * disp_b
                    mx, my = to_m(phys_b[1], phys_b[0])
                    cifp_anchors_m.append((mx, my, elev_pb))
    # Add road DEM anchors (non-tunnel surface roads near airport)
    cifp_anchors_m.extend(road_anchors_m)

    # ── Step 2: Building platform elevations ───────────────────────
    # Sample DEM at building centroid.  Buildings were constructed at
    # the natural grade — DEM is most reliable at surveyed structures.
    # Validate: if platform elevation differs from nearest CIFP anchor
    # by more than plausible grade would allow, prefer CIFP-derived.
    bldg_elevations = {}
    complex_bldg_indices = set()
    platform_anchors_ll = []  # (lat, lon, elevation) for taxiway anchoring

    for bi, bp in enumerate(building_polys_m):
        cx_m, cy_m = bp.centroid.x, bp.centroid.y
        clon, clat = to_ll(cx_m, cy_m)

        # Sample DEM at centroid and boundary
        dem_samples = []
        try:
            dem_samples.append(
                tile.dem.alt((clon - tile.lon, clat - tile.lat))
            )
        except Exception:
            pass
        for coord in bp.exterior.coords[:-1]:
            blon, blat = to_ll(coord[0], coord[1])
            try:
                dem_samples.append(
                    tile.dem.alt((blon - tile.lon, blat - tile.lat))
                )
            except Exception:
                pass

        if len(dem_samples) >= 2:
            dem_range = max(dem_samples) - min(dem_samples)
        else:
            dem_range = 0.0

        # Primary: DEM centroid elevation (building sits on natural grade)
        dem_centroid = dem_samples[0] if dem_samples else None

        # Validation: check against nearest CIFP anchor
        cifp_elev = _project_point_onto_runway(
            clat, clon, runway_pairs, max_dist=300.0
        )

        if dem_centroid is not None and cifp_elev is not None:
            # Find distance to nearest runway for grade validation
            min_rwy_dist = float("inf")
            for rp in (runway_pairs or []):
                for side in ("data_a", "data_b"):
                    d = rp.get(side)
                    if d:
                        rdx = (cx_m - to_m(d["lon"], d["lat"])[0])
                        rdy = (cy_m - to_m(d["lon"], d["lat"])[1])
                        rdist = sqrt(rdx * rdx + rdy * rdy)
                        if rdist < min_rwy_dist:
                            min_rwy_dist = rdist
            max_plausible_delta = min_rwy_dist * MAX_SURFACE_GRADE + 2.0
            if abs(dem_centroid - cifp_elev) > max_plausible_delta:
                # DEM is implausible — likely bad data (creek, artifact)
                # Use CIFP-derived elevation instead
                elev = cifp_elev
            else:
                # DEM is plausible — trust building sits on natural grade
                elev = dem_centroid
        elif dem_centroid is not None:
            elev = dem_centroid
        elif cifp_elev is not None:
            elev = cifp_elev
        else:
            elev = 0.0

        bldg_elevations[bi] = round(elev, 1)
        platform_anchors_ll.append((clat, clon, round(elev, 1)))
        if dem_range > FLAT_TERRAIN_THRESHOLD:
            complex_bldg_indices.add(bi)

        # Add building centroid as a CIFP surface anchor
        cifp_anchors_m.append((cx_m, cy_m, round(elev, 1)))

    UI.vprint(2, "    {}: {} buildings ({} simple, {} complex), "
              "{} platform anchors".format(
                  icao,
                  len(building_polys_m),
                  len(building_polys_m) - len(complex_bldg_indices),
                  len(complex_bldg_indices),
                  len(platform_anchors_ll)))

    # ── Step 3: Taxiway elevations ─────────────────────────────────
    # Now that building platforms are known, compute taxiway node
    # elevations.  Taxiways anchor to CIFP runways AND building
    # platforms, then grade-constrain the path between them.
    for twy_idx, raw_cl in enumerate(twy_raw_centerlines):
        if raw_cl is None:
            twy_centerlines.append(None)
            continue
        elevations = _compute_taxiway_elevations(
            raw_cl, runway_pairs, tile,
            platform_anchors=platform_anchors_ll
        )
        twy_centerlines.append((raw_cl, elevations))

    # ── Step 3b: Adjust building elevations for grade compliance ───
    # Now that taxiway elevations are grade-constrained, check if any
    # building's DEM elevation creates an impossible grade to/from its
    # nearest taxiway.  If the grade from building to nearest taxiway
    # node exceeds MAX_SURFACE_GRADE (1.5%), adjust the building
    # elevation toward the taxiway until the grade is achievable.
    # This is the "workable compromise" step — runways are fixed
    # (CIFP), taxiways are grade-constrained, buildings yield when
    # the geometry doesn't work.
    n_bldg_adjusted = 0
    for bi, bp in enumerate(building_polys_m):
        if bp.is_empty or bi not in bldg_elevations:
            continue
        bldg_elev = bldg_elevations[bi]
        cx_m, cy_m = bp.centroid.x, bp.centroid.y

        # Find nearest taxiway node and its grade-constrained elevation
        best_dist = float("inf")
        best_twy_elev = None
        for twy_idx, twy_cl in enumerate(twy_centerlines):
            if twy_cl is None:
                continue
            _, elevations = twy_cl
            cl_m = twy_centerlines_m[twy_idx]
            if cl_m is None:
                continue
            for i, (nx, ny) in enumerate(cl_m):
                if i >= len(elevations):
                    continue
                d = sqrt((nx - cx_m) ** 2 + (ny - cy_m) ** 2)
                if d < best_dist:
                    best_dist = d
                    best_twy_elev = elevations[i]

        if best_twy_elev is None or best_dist < 1.0:
            continue

        # Check grade from building to nearest taxiway
        delta = abs(bldg_elev - best_twy_elev)
        grade = delta / best_dist
        if grade > MAX_SURFACE_GRADE:
            # Grade exceeds limit — adjust building elevation
            # Move building toward taxiway elevation until grade = limit
            max_delta = best_dist * MAX_SURFACE_GRADE
            if bldg_elev > best_twy_elev:
                new_elev = best_twy_elev + max_delta
            else:
                new_elev = best_twy_elev - max_delta
            new_elev = round(new_elev, 1)
            bldg_elevations[bi] = new_elev
            # Update the corresponding anchor in cifp_anchors_m
            # (building centroid was added at line ~2298)
            for ai in range(len(cifp_anchors_m)):
                ax, ay, ae = cifp_anchors_m[ai]
                if (abs(ax - cx_m) < 0.1 and abs(ay - cy_m) < 0.1
                        and abs(ae - bldg_elev) < 0.2):
                    cifp_anchors_m[ai] = (ax, ay, new_elev)
                    break
            # Update platform_anchors_ll too
            clon, clat = to_ll(cx_m, cy_m)
            for pai in range(len(platform_anchors_ll)):
                pa_lat, pa_lon, pa_e = platform_anchors_ll[pai]
                if (abs(pa_lat - clat) < 0.00001
                        and abs(pa_lon - clon) < 0.00001):
                    platform_anchors_ll[pai] = (
                        pa_lat, pa_lon, new_elev)
                    break
            n_bldg_adjusted += 1

    if n_bldg_adjusted > 0:
        UI.vprint(1, "    {}: adjusted {} building elevations for "
                  "grade compliance".format(icao, n_bldg_adjusted))

    # ── Step 4: Full CIFP surface model ────────────────────────────
    # Add grade-constrained taxiway node elevations as anchors.
    # These are the actual negotiated elevations from Step 3 (not
    # runway projections), so the IDW interpolation reflects the
    # real surface: runways → taxiways → buildings, all connected
    # at grade-limited slopes.
    for twy_idx, twy_cl in enumerate(twy_centerlines):
        if twy_cl is None:
            continue
        centerline, elevations = twy_cl
        cl_m = twy_centerlines_m[twy_idx]
        if cl_m is None:
            continue
        for i in range(len(centerline)):
            if i < len(cl_m) and i < len(elevations):
                cifp_anchors_m.append(
                    (cl_m[i][0], cl_m[i][1], elevations[i]))

    UI.vprint(2, "    {}: {} total elevation anchors (CIFP + platforms + "
              "taxiway + road)".format(icao, len(cifp_anchors_m)))

    # Build spatial index for fast anchor lookup
    _anchor_points = [shp_geom.Point(ax, ay) for ax, ay, _ in cifp_anchors_m]
    _anchor_tree = None
    if _anchor_points:
        from shapely.strtree import STRtree
        _anchor_tree = STRtree(_anchor_points)

    def _cifp_surface_elevation(x_m, y_m, max_search=500.0):
        """Query the CIFP surface model at a point in meter coords.

        Uses inverse-distance-weighted interpolation from the nearest
        anchors (CIFP thresholds, building platforms, projected taxiway
        nodes, road DEM samples).  Spatial-indexed for performance.

        Returns: elevation in meters, or None if no data available.
        """
        if not cifp_anchors_m or _anchor_tree is None:
            return None

        # Spatial query: find anchors within max_search radius
        search_area = shp_geom.Point(x_m, y_m).buffer(max_search)
        try:
            idxs = _anchor_tree.query(search_area)
        except Exception:
            return None

        candidates = []
        for idx in idxs:
            i = int(idx)
            ax, ay, ae = cifp_anchors_m[i]
            dist = sqrt((x_m - ax) ** 2 + (y_m - ay) ** 2)
            if dist < max_search:
                candidates.append((dist, ae))

        if not candidates:
            return None

        candidates.sort(key=lambda c: c[0])

        # Very close to an anchor — use it directly
        if candidates[0][0] < 5.0:
            return candidates[0][1]

        # IDW from up to 8 nearest anchors
        top_n = candidates[:8]
        weight_sum = 0.0
        elev_sum = 0.0
        for dist, ae in top_n:
            w = 1.0 / (dist * dist)
            weight_sum += w
            elev_sum += w * ae

        if weight_sum < 1e-12:
            return None

        return elev_sum / weight_sum

    # ── Terminal paved zones ──────────────────────────────────────
    # Identify large contiguous paved areas around terminals.
    # These encompass the terminal building + surrounding aprons +
    # any taxiways running through.  They'll be treated as unified
    # surfaces with large triangles for compound slopes, rather than
    # hundreds of small independent shapes.
    terminal_zone_m = shp_geom.Polygon()
    terminal_zones_list = []  # individual terminal zone polygons
    try:
        large_bldgs = [bp for bp in building_polys_m
                       if not bp.is_empty and bp.area > 500.0]
        if large_bldgs:
            for lb in large_bldgs:
                # Buffer by 15m to create the terminal influence zone
                tz = lb.buffer(BLDG_MERGE_DIST)
                # Union with any aprons that intersect
                for ap in apron_polys_m:
                    if not ap.is_empty and tz.intersects(ap):
                        tz = tz.union(ap)
                # Union with any taxiway buffers that are mostly inside
                for tb in twy_buffers_m:
                    if tb.is_empty:
                        continue
                    try:
                        overlap = tz.intersection(tb)
                        if not overlap.is_empty and overlap.area > 0.5 * tb.area:
                            tz = tz.union(tb)
                    except Exception:
                        pass
                tz = tz.simplify(2.0, preserve_topology=True)
                if tz.is_valid and not tz.is_empty and tz.area > 1000.0:
                    terminal_zones_list.append(tz)
            if terminal_zones_list:
                terminal_zone_m = shp_ops.unary_union(terminal_zones_list)
                UI.vprint(2, "    {}: {} terminal paved zone(s), total {:.0f} m²"
                          .format(icao, len(terminal_zones_list),
                                  terminal_zone_m.area))
    except Exception:
        pass

    # ════════════════════════════════════════════════════════════════
    # PHASE B: Compute junction zones (areas needing triangles)
    # ════════════════════════════════════════════════════════════════
    # Junctions = any area where surfaces overlap:
    #   - Taxiway buffer ∩ runway polygon
    #   - Taxiway buffer ∩ another taxiway buffer
    # These zones get triangles; the overlapping shapes get clipped.

    junction_parts = []

    # B1: Taxiway–runway overlaps
    if not rwy_union_m.is_empty and not twy_union_m.is_empty:
        try:
            twy_rwy = twy_union_m.intersection(rwy_union_m)
            if not twy_rwy.is_empty:
                junction_parts.append(twy_rwy)
        except Exception:
            pass

    # B2: Taxiway–taxiway overlaps (pairwise)
    for i in range(len(twy_buffers_m)):
        if twy_buffers_m[i].is_empty:
            continue
        for j in range(i + 1, len(twy_buffers_m)):
            if twy_buffers_m[j].is_empty:
                continue
            try:
                overlap = twy_buffers_m[i].intersection(twy_buffers_m[j])
                if not overlap.is_empty and overlap.area > 1.0:
                    junction_parts.append(overlap)
            except Exception:
                pass

    junction_zone_raw = (shp_ops.unary_union(junction_parts)
                         if junction_parts else shp_geom.Polygon())
    # Buffer the junction zone slightly so it captures the full
    # transition area and doesn't leave tiny flat-poly scraps at
    # the taxiway–junction boundary.
    try:
        junction_zone_m = (junction_zone_raw.buffer(3.0)
                           if not junction_zone_raw.is_empty
                           else shp_geom.Polygon())
    except Exception:
        junction_zone_m = junction_zone_raw

    UI.vprint(2, "    {}: junction zone area {:.0f} m²".format(
        icao, junction_zone_m.area if not junction_zone_m.is_empty else 0))

    # (Building elevations and CIFP surface model already computed
    #  in Phase A5 above.)

    # ════════════════════════════════════════════════════════════════
    # PHASE C: Clip and emit non-overlapping shapes
    # ════════════════════════════════════════════════════════════════
    # Priority: runways (already emitted) > junctions > taxiways >
    #           aprons > buildings (buildings nest inside others).
    # Each shape subtracts higher-priority coverage before emission.
    #
    # emitted_flat_shapes tracks ALL flat shapes with elevations for
    # Phase C4 transition strip detection between adjacent shapes.

    emitted_flat_shapes = []  # [(poly_m, elevation), ...]

    # C1: Taxiway shapes — clip away runways + junction zones
    n_twy_flat = 0
    n_twy_sloped = 0
    sloped_twy_parts_m = []   # curve areas sent to Phase D for triangulation
    emitted_twy_quads_m = []  # meter-space footprints of emitted sloped quads
    emitted_flat_twy_m = []   # meter-space footprints of emitted flat twy polys

    # Merge parallel taxiways within 15m into single surfaces
    PARALLEL_TWY_MERGE_DIST = 15.0
    merged_twy_groups = {}  # ti -> group_id
    twy_group_id = 0
    for i in range(len(twy_buffers_m)):
        if twy_buffers_m[i].is_empty or i in merged_twy_groups:
            continue
        group = [i]
        for j in range(i + 1, len(twy_buffers_m)):
            if twy_buffers_m[j].is_empty or j in merged_twy_groups:
                continue
            try:
                dist = twy_buffers_m[i].distance(twy_buffers_m[j])
                if dist < PARALLEL_TWY_MERGE_DIST:
                    group.append(j)
            except Exception:
                pass
        if len(group) > 1:
            # Merge these into a single wider buffer
            merged = shp_ops.unary_union(
                [twy_buffers_m[g] for g in group])
            if merged.is_valid and not merged.is_empty:
                twy_buffers_m[group[0]] = merged
                for g in group[1:]:
                    twy_buffers_m[g] = shp_geom.Polygon()
                    merged_twy_groups[g] = group[0]
            twy_group_id += 1

    if merged_twy_groups:
        # Recompute twy_union_m
        twy_union_m = shp_ops.unary_union(
            [tb for tb in twy_buffers_m if not tb.is_empty]
        ) if any(not tb.is_empty for tb in twy_buffers_m) else shp_geom.Polygon()
        UI.vprint(2, "    {}: merged {} parallel taxiway groups".format(
            icao, len(set(merged_twy_groups.values()))))

    for ti, entry in enumerate(twy_centerlines):
        if entry is None:
            continue
        centerline, elevations = entry
        cl_m = twy_centerlines_m[ti]
        buf = twy_buffers_m[ti]
        if buf.is_empty:
            continue

        # Skip taxiways mostly inside terminal zones — they'll be
        # handled as part of the unified terminal paved surface
        if not terminal_zone_m.is_empty:
            try:
                tz_overlap = buf.intersection(terminal_zone_m)
                if (not tz_overlap.is_empty
                        and tz_overlap.area > 0.5 * buf.area):
                    continue  # skip — terminal zone handles this
            except Exception:
                pass

        # Clip this taxiway buffer: subtract runways + junctions +
        # other taxiway buffers (prevents cross-taxiway overlap).
        try:
            other_twy_m = twy_union_m.difference(buf)
        except Exception:
            other_twy_m = shp_geom.Polygon()

        clipped = buf
        try:
            for subtract in (rwy_union_m, junction_zone_m, other_twy_m):
                if not subtract.is_empty:
                    clipped = clipped.difference(subtract)
        except Exception:
            continue
        if clipped.is_empty:
            continue

        # Determine elevation range along the centerline
        e_min, e_max = min(elevations), max(elevations)
        elev_range = e_max - e_min

        if elev_range <= DEM_VARIANCE_THRESHOLD:
            # FLAT: emit the clipped buffer as a single flat polygon
            avg_e = round((e_min + e_max) / 2.0, 1)
            # Handle MultiPolygon from clipping
            polys = ([clipped] if clipped.geom_type == "Polygon"
                     else list(clipped.geoms)
                     if hasattr(clipped, "geoms") else [])
            for p in polys:
                if p.is_empty or not hasattr(p, "exterior"):
                    continue
                if p.area < 5.0:
                    continue    # skip tiny slivers from clipping
                # Small fragments near junction zones should go to
                # triangulation — emitting them flat creates elevation
                # cliffs against the junction's triangulated surface.
                MIN_FLAT_AREA = 200.0  # m² — smaller goes to Phase D
                near_junction = False
                if not junction_zone_m.is_empty and p.area < MIN_FLAT_AREA:
                    try:
                        near_junction = (p.distance(junction_zone_m) < 5.0)
                    except Exception:
                        pass
                if near_junction:
                    sloped_twy_parts_m.append(p)
                    continue
                ring_ll = [to_ll(mx, my) for mx, my in p.exterior.coords]
                _emit_flat_poly(ring_ll, avg_e)
                emitted_flat_twy_m.append(p)
                emitted_flat_shapes.append((p, avg_e))
                n_twy_flat += 1
        else:
            # SLOPING: build an offset polygon from the centerline
            # so that adjacent quads share exact edge vertices.
            # Zero same-taxiway overlap by construction.
            # other_twy_m is already computed above for clipping.
            #
            # At each interior vertex, compute the miter (bisector)
            # offset so left[i]/right[i] are shared between segment
            # i-1 and segment i.  At endpoints, use simple perpendicular.
            n_verts = len(cl_m)
            hw = TAXIWAY_BUFFER_WIDTH   # half-width in meters
            left_m = []
            right_m = []

            for vi in range(n_verts):
                if vi == 0:
                    # First vertex: perpendicular to first segment
                    dx = cl_m[1][0] - cl_m[0][0]
                    dy = cl_m[1][1] - cl_m[0][1]
                elif vi == n_verts - 1:
                    # Last vertex: perpendicular to last segment
                    dx = cl_m[-1][0] - cl_m[-2][0]
                    dy = cl_m[-1][1] - cl_m[-2][1]
                else:
                    # Interior vertex: miter join of two segments
                    dx1 = cl_m[vi][0] - cl_m[vi - 1][0]
                    dy1 = cl_m[vi][1] - cl_m[vi - 1][1]
                    dx2 = cl_m[vi + 1][0] - cl_m[vi][0]
                    dy2 = cl_m[vi + 1][1] - cl_m[vi][1]
                    len1 = sqrt(dx1 ** 2 + dy1 ** 2)
                    len2 = sqrt(dx2 ** 2 + dy2 ** 2)
                    if len1 < 0.01 and len2 < 0.01:
                        left_m.append(cl_m[vi])
                        right_m.append(cl_m[vi])
                        continue
                    if len1 < 0.01:
                        dx, dy = dx2, dy2
                    elif len2 < 0.01:
                        dx, dy = dx1, dy1
                    else:
                        # Unit normals pointing left (rotate 90° CCW)
                        nx1, ny1 = -dy1 / len1, dx1 / len1
                        nx2, ny2 = -dy2 / len2, dx2 / len2
                        # Miter direction = average of the two normals
                        mx = nx1 + nx2
                        my = ny1 + ny2
                        mlen = sqrt(mx ** 2 + my ** 2)
                        if mlen < 0.01:
                            # ~180° turn, fall back to first normal
                            mx, my = nx1, ny1
                            mlen = 1.0
                        mx /= mlen
                        my /= mlen
                        # Adjust for miter angle so perpendicular
                        # distance stays at half-width
                        cos_half = mx * nx1 + my * ny1
                        if cos_half < 0.3:
                            cos_half = 0.3   # cap extreme miter length
                        adj_hw = hw / cos_half
                        left_m.append((cl_m[vi][0] + mx * adj_hw,
                                       cl_m[vi][1] + my * adj_hw))
                        right_m.append((cl_m[vi][0] - mx * adj_hw,
                                        cl_m[vi][1] - my * adj_hw))
                        continue

                    # (fall-through for degenerate single-segment case)

                # Simple perpendicular for endpoints (or degenerate)
                seg_len = sqrt(dx ** 2 + dy ** 2)
                if seg_len < 0.01:
                    left_m.append(cl_m[vi])
                    right_m.append(cl_m[vi])
                    continue
                nx = -dy / seg_len   # left normal
                ny = dx / seg_len
                left_m.append((cl_m[vi][0] + nx * hw,
                               cl_m[vi][1] + ny * hw))
                right_m.append((cl_m[vi][0] - nx * hw,
                                cl_m[vi][1] - ny * hw))

            if len(left_m) != n_verts or len(right_m) != n_verts:
                continue

            # Emit per-segment quads.  Each quad is defined by the
            # shared offset points — adjacent quads share an edge
            # exactly, so there is zero overlap between siblings.
            for si in range(n_verts - 1):
                elev_a = elevations[si]
                elev_b = elevations[si + 1]

                # Quad corners in meter space:
                #   left[si] → left[si+1] → right[si+1] → right[si]
                quad_coords_m = [
                    left_m[si], left_m[si + 1],
                    right_m[si + 1], right_m[si],
                    left_m[si],   # close ring
                ]
                try:
                    quad_m = shp_geom.Polygon(quad_coords_m)
                    if not quad_m.is_valid or quad_m.is_empty:
                        continue
                except Exception:
                    continue

                # Check if this segment is at a significant turn with
                # elevation change — these need triangulation, not a
                # flat step that creates sharp grade breaks.
                if si > 0 and si < n_verts - 2:
                    prev_dx = cl_m[si][0] - cl_m[si - 1][0]
                    prev_dy = cl_m[si][1] - cl_m[si - 1][1]
                    curr_dx = cl_m[si + 1][0] - cl_m[si][0]
                    curr_dy = cl_m[si + 1][1] - cl_m[si][1]
                    prev_len = sqrt(prev_dx**2 + prev_dy**2)
                    curr_len = sqrt(curr_dx**2 + curr_dy**2)
                    if prev_len > 0.1 and curr_len > 0.1:
                        cos_angle = ((prev_dx * curr_dx + prev_dy * curr_dy)
                                     / (prev_len * curr_len))
                        cos_angle = max(-1.0, min(1.0, cos_angle))
                        turn_angle = abs(acos(cos_angle))
                        elev_diff = abs(elev_a - elev_b)
                        if turn_angle > 0.35 and elev_diff > 0.3:
                            # Sharp turn + elevation change → triangulate
                            sloped_twy_parts_m.append(quad_m)
                            continue

                # Check whether this quad overlaps runways, junction
                # zones, or other taxiways' buffers.  If it does, we
                # CANNOT emit it as-is.  Send the clipped portion to
                # Phase D instead.
                #
                # other_twy_m = union of all taxiway buffers except this
                # one.  Catches quads whose miter-join corners extend
                # into a neighboring taxiway's area.
                has_external_overlap = False
                try:
                    for zone in (rwy_union_m, junction_zone_m,
                                 other_twy_m):
                        if not zone.is_empty:
                            ix = quad_m.intersection(zone)
                            if ix.area > 0.1:
                                has_external_overlap = True
                                break
                except Exception:
                    has_external_overlap = True

                if has_external_overlap:
                    # Clip and send to Phase D triangulation
                    try:
                        clipped_q = quad_m
                        for zone in (rwy_union_m, junction_zone_m):
                            if not zone.is_empty:
                                clipped_q = clipped_q.difference(zone)
                        if not clipped_q.is_empty:
                            sloped_twy_parts_m.append(clipped_q)
                    except Exception:
                        pass
                    continue

                # No external overlap — emit directly.
                # Convert corners to lat/lon.  to_ll returns (lon, lat).
                ll_i_L = to_ll(left_m[si][0], left_m[si][1])
                ll_ip_L = to_ll(left_m[si + 1][0], left_m[si + 1][1])
                ll_ip_R = to_ll(right_m[si + 1][0], right_m[si + 1][1])
                ll_i_R = to_ll(right_m[si][0], right_m[si][1])

                lat_iL, lon_iL = ll_i_L[1], ll_i_L[0]
                lat_ipL, lon_ipL = ll_ip_L[1], ll_ip_L[0]
                lat_ipR, lon_ipR = ll_ip_R[1], ll_ip_R[0]
                lat_iR, lon_iR = ll_i_R[1], ll_i_R[0]

                if abs(elev_a - elev_b) >= 0.1:
                    # Sloped quad.  altitude_high/altitude_low needs
                    # corner order: high-left → low-left → low-right
                    # → high-right.
                    if elev_a >= elev_b:
                        eh, el = elev_a, elev_b
                        n0 = _add_node(lat_iL, lon_iL)
                        n1 = _add_node(lat_ipL, lon_ipL)
                        n2 = _add_node(lat_ipR, lon_ipR)
                        n3 = _add_node(lat_iR, lon_iR)
                    else:
                        eh, el = elev_b, elev_a
                        n0 = _add_node(lat_ipL, lon_ipL)
                        n1 = _add_node(lat_iL, lon_iL)
                        n2 = _add_node(lat_iR, lon_iR)
                        n3 = _add_node(lat_ipR, lon_ipR)
                    _add_way([n0, n1, n2, n3, n0], {
                        "altitude_high": "{:.1f}".format(eh),
                        "altitude_low": "{:.1f}".format(el),
                        "cell_size": str(int(DEFAULT_CELL_SIZE)),
                        "profile": DEFAULT_PROFILE,
                    })
                else:
                    # Flat segment
                    avg_e = round((elev_a + elev_b) / 2.0, 1)
                    n0 = _add_node(lat_iL, lon_iL)
                    n1 = _add_node(lat_ipL, lon_ipL)
                    n2 = _add_node(lat_ipR, lon_ipR)
                    n3 = _add_node(lat_iR, lon_iR)
                    _add_way([n0, n1, n2, n3, n0], {
                        "altitude": "{:.1f}".format(avg_e),
                    })

                emitted_twy_quads_m.append(quad_m)
                quad_avg_e = round((elev_a + elev_b) / 2.0, 1)
                emitted_flat_shapes.append((quad_m, quad_avg_e))
                n_twy_sloped += 1

    UI.vprint(2, "    {}: {} flat taxiway polys, {} sloped quads".format(
        icao, n_twy_flat, n_twy_sloped))

    # ── C2: Apron shapes ────────────────────────────────────────────
    # New apron emission strategy (commit 5):
    #
    #   1. Each OSM apron polygon is clipped against runways,
    #      taxiways, junction zones, and buildings (so it never
    #      overlaps any of those features).
    #   2. Anchor elevations are collected from neighbouring features
    #      using the elevation priority documented in STATUS.md:
    #        * building edges that touch the apron → bldg_elevations
    #        * taxiway centerline crossings        → twy_centerlines
    #      CIFP runway elevations are NOT used as apron anchors —
    #      they belong to the runways alone.
    #   3. The apron is triangulated by O4_Surface_Mesh.adaptive_-
    #      triangulate, which produces the smallest set of
    #      grade-compliant DEM-faithful triangles (typically 2-15
    #      per apron, vs the legacy's ~50).
    #   4. Triangles are emitted via _emit_triangle (node_altitudes
    #      tag, allowed for 3-vertex shapes only).
    #
    # The "complex vs flat" branch is gone: adaptive_triangulate
    # naturally yields just 2 triangles for a flat apron and adds
    # interior detail only where DEM departs from the plane.
    import O4_Surface_Mesh as SM

    def _dem_at(x_m, y_m):
        """DEM sampler closure for the current `to_ll`/tile."""
        try:
            lon, lat = to_ll(x_m, y_m)
            return float(tile.dem.alt((lon - tile.lon, lat - tile.lat)))
        except Exception:
            return None

    APRON_SIMPLIFY_M = 10.0      # input polygon simplification
                                 # tolerance — strips OSM vertex
                                 # noise so adaptive_triangulate
                                 # starts with a small seed set.
                                 # 10 m is well below the 1.0 %
                                 # apron grade rule's resolution
                                 # (a 10 m horizontal step at 1 %
                                 # is 10 cm vertical) so simplifying
                                 # at this scale loses no real
                                 # geometry that the grade rule cares
                                 # about.
    APRON_ANCHOR_GAP_M = 50.0    # minimum spacing between anchors
                                 # of the same kind, to keep the
                                 # triangulation seed sparse

    def _sparsify_anchors(anchors, min_gap):
        """Return a subset of `anchors` with no two within `min_gap`
        of each other.  Greedy: walk the input order, skip any
        candidate that's already covered.
        """
        kept = []
        gap2 = min_gap * min_gap
        for ax, ay, az in anchors:
            covered = False
            for kx, ky, _kz in kept:
                if (ax - kx) ** 2 + (ay - ky) ** 2 < gap2:
                    covered = True
                    break
            if not covered:
                kept.append((ax, ay, az))
        return kept

    def _apron_anchors(apron_poly):
        """Collect (x, y, z) anchor elevations from neighbouring
        features that touch this apron.  Three anchor sources:

        1. Building edges: sample points along the shared boundary
           between each touching building and the apron, at the
           building's pad elevation.

        2. **Flat terminal buffer** (new in commit 9): for each
           touching building, also sample points along an OUTWARD
           offset ring at TERMINAL_FLAT_BUFFER_M metres from the
           building's footprint, all at the building's elevation.
           This creates a band of "flat" anchors around each
           terminal so the apron triangulation is forced to stay at
           the terminal's elevation for ~10 m around the building
           before sloping away to the nearest taxiway.  Result: the
           ground around terminals is flat, then slopes gently
           (≤ 1.5 %) to the taxiway, matching how real airports are
           graded.

        3. Taxiway centerline crossings at the apron boundary,
           elev interpolated from the taxiway's routed elevation.

        Per STATUS.md elevation rule, CIFP runway elevations are
        NOT used here.  The result is sparsified to
        APRON_ANCHOR_GAP_M minimum spacing.
        """
        raw = []
        # Building edge anchors.  A building "touches" the apron if
        # its padded footprint comes within 6 m of the apron polygon
        # (the legacy A4b padding + simplification can leave a small
        # gap between the two outlines even when they share an edge
        # in the source data).
        TOUCH_M = 6.0
        TERMINAL_FLAT_BUFFER_M = 10.0   # flat zone half-width
        for bi, bp in enumerate(building_polys_m):
            if bp.is_empty or bi not in bldg_elevations:
                continue
            try:
                if bp.distance(apron_poly) > TOUCH_M:
                    continue
            except Exception:
                continue
            elev = float(bldg_elevations[bi])
            try:
                shared = bp.boundary.intersection(
                    apron_poly.buffer(TOUCH_M))
            except Exception:
                shared = None
            if shared is None or shared.is_empty:
                raw.append((bp.centroid.x, bp.centroid.y, elev))
            else:
                lines = []
                if hasattr(shared, "geoms"):
                    for g in shared.geoms:
                        if hasattr(g, "coords"):
                            lines.append(g)
                elif hasattr(shared, "coords"):
                    lines.append(shared)
                for ln in lines:
                    length = ln.length
                    if length < 1.0:
                        p = ln.interpolate(0.5, normalized=True)
                        raw.append((p.x, p.y, elev))
                        continue
                    # Two endpoints + midpoint.  Sparsifier will
                    # trim further if the building is very small.
                    for t in (0.0, 0.5, 1.0):
                        p = ln.interpolate(t, normalized=True)
                        raw.append((p.x, p.y, elev))
            # Flat terminal buffer ring: sample an outward offset
            # of the building's boundary every ~20 m, clipped to
            # the apron polygon.  These anchors force the apron
            # triangulation to be flat at the terminal elevation
            # for TERMINAL_FLAT_BUFFER_M around the building.
            try:
                ring_geom = bp.buffer(TERMINAL_FLAT_BUFFER_M).boundary
                clipped_ring = ring_geom.intersection(apron_poly)
            except Exception:
                clipped_ring = None
            if clipped_ring is not None and not clipped_ring.is_empty:
                ring_lines = []
                if hasattr(clipped_ring, "geoms"):
                    for g in clipped_ring.geoms:
                        if hasattr(g, "coords"):
                            ring_lines.append(g)
                elif hasattr(clipped_ring, "coords"):
                    ring_lines.append(clipped_ring)
                for rl in ring_lines:
                    total = rl.length
                    if total < 1.0:
                        p = rl.interpolate(0.5, normalized=True)
                        raw.append((p.x, p.y, elev))
                        continue
                    n = max(2, int(total / 20.0) + 1)
                    for k in range(n):
                        t = k / (n - 1) if n > 1 else 0.5
                        p = rl.interpolate(t, normalized=True)
                        raw.append((p.x, p.y, elev))
        # Taxiway crossing anchors.  Use a slightly buffered apron
        # so taxiways that stub into the edge also contribute.
        TWY_TOUCH_M = 5.0
        zone = apron_poly.buffer(TWY_TOUCH_M)
        for ti2, cl_m2 in enumerate(twy_centerlines_m):
            if cl_m2 is None:
                continue
            entry2 = twy_centerlines[ti2]
            if entry2 is None:
                continue
            cl_elevs = entry2[1]
            if not cl_elevs or len(cl_elevs) != len(cl_m2):
                continue
            try:
                ls2 = shp_geom.LineString(cl_m2)
                crossing = ls2.intersection(zone)
            except Exception:
                continue
            if crossing.is_empty:
                continue
            # One representative point per crossing component.
            comps = (list(crossing.geoms)
                     if hasattr(crossing, "geoms") else [crossing])
            for comp in comps:
                if not hasattr(comp, "coords"):
                    continue
                coords = list(comp.coords)
                if not coords:
                    continue
                # Use the midpoint of the crossing
                mid = coords[len(coords) // 2]
                px, py = mid[0], mid[1]
                ie = _interpolate_elevation_along_centerline(
                    (px, py), cl_m2, cl_elevs)
                raw.append((px, py, float(ie)))
        return _sparsify_anchors(raw, APRON_ANCHOR_GAP_M)

    n_apron = 0
    n_apron_tris = 0
    complex_apron_parts_m = []   # left in scope for Phase D / coverage

    for ai, p_m in enumerate(apron_polys_m):
        clipped = p_m
        try:
            for subtract in (rwy_union_m, twy_union_m,
                             junction_zone_m, bldg_union_m):
                if not subtract.is_empty:
                    clipped = clipped.difference(subtract)
        except Exception:
            continue
        if clipped.is_empty:
            continue

        pieces = ([clipped] if clipped.geom_type == "Polygon"
                  else list(clipped.geoms)
                  if hasattr(clipped, "geoms") else [])
        for piece in pieces:
            if (piece.is_empty or not hasattr(piece, "exterior")
                    or piece.area < 10.0):
                continue

            # Simplify the piece's outline before triangulation so
            # adaptive_triangulate starts with a small seed set.
            # OSM apron polygons commonly have 50-200 contour
            # vertices that don't carry meaningful shape information
            # at the patch-mesh scale.  10 m tolerance trims those
            # without losing the apron's recognisable outline.
            #
            # When apt.dat is the source, SKIP simplification: the
            # polygons are already clean, and a 10 m tolerance can
            # silently delete small interior rings (building-shaped
            # holes carved by the rwy/bldg subtraction a few lines
            # above) which then causes triangles to cover building
            # pads and produce flat-triangle overlap.
            if apt_dat_used:
                seed_poly = piece
            else:
                try:
                    seed_poly = piece.simplify(
                        APRON_SIMPLIFY_M, preserve_topology=True)
                    if (seed_poly.is_empty
                            or not hasattr(seed_poly, "exterior")):
                        seed_poly = piece
                except Exception:
                    seed_poly = piece

            anchors = _apron_anchors(seed_poly)
            # apt.dat pavements lump taxiways, aprons, and ramps into
            # one "pavement" category.  Per user direction we treat
            # them all as taxiway for now (1.5 % grade) rather than
            # the stricter 1.0 % apron grade — a finer name-based
            # classification is a follow-up commit.
            pavement_grade = (MAX_TAXIWAY_GRADE if apt_dat_used
                              else MAX_APRON_GRADE)
            try:
                tris = SM.adaptive_triangulate(
                    seed_poly, anchors, _dem_at,
                    max_grade=pavement_grade,
                    fidelity_tol=1.0,
                    max_extra_points=50)
            except Exception:
                tris = []

            if not tris:
                # Degenerate fallback: emit the piece as a flat
                # polygon at the average anchor elevation (or DEM
                # centroid if no anchors).
                if anchors:
                    avg_e = sum(a[2] for a in anchors) / len(anchors)
                else:
                    cx_m, cy_m = piece.centroid.x, piece.centroid.y
                    avg_e = _dem_at(cx_m, cy_m) or 0.0
                ring_ll = [to_ll(mx, my)
                           for mx, my in piece.exterior.coords]
                _emit_flat_poly(ring_ll, round(avg_e, 1))
                emitted_flat_shapes.append((piece, round(avg_e, 1)))
                n_apron += 1
                continue

            # Emit each triangle from the adaptive mesh.
            avg_e = (sum(v[2] for t in tris for v in t)
                     / (3.0 * len(tris)))
            for (v0, v1, v2) in tris:
                lon0, lat0 = to_ll(v0[0], v0[1])
                lon1, lat1 = to_ll(v1[0], v1[1])
                lon2, lat2 = to_ll(v2[0], v2[1])
                _emit_triangle(
                    lat0, lon0, round(v0[2], 1),
                    lat1, lon1, round(v1[2], 1),
                    lat2, lon2, round(v2[2], 1))
                n_apron_tris += 1
            # Record the apron piece for emitted_flat_shapes overlap
            # tracking — the actual emission was a triangle set, but
            # the overlap accumulator only needs a footprint and a
            # representative elevation.
            emitted_flat_shapes.append((piece, round(avg_e, 1)))
            complex_apron_parts_m.append(piece)
            n_apron += 1

    UI.vprint(2, "    {}: {} apron pieces emitted as {} triangles"
              " (adaptive mesh, no CIFP anchors)".format(
                  icao, n_apron, n_apron_tris))
    # Legacy variables that downstream code still reads:
    n_apron_flat = 0          # all aprons go through the unified
    n_apron_complex = n_apron #   adaptive path now

    # C3: Building flat pads — clip against runways + taxiways so pads
    # don't overlap pavement shapes.  Buildings CAN nest inside aprons
    # (aprons already have building-shaped holes from C2 subtraction).
    # Elevations are pre-computed before Phase C.

    for bi, bp in enumerate(building_polys_m):
        if bi not in bldg_elevations:
            continue
        # Clip building pad against runways + taxiways + junctions
        clipped_bp = bp
        try:
            for subtract in (rwy_union_m, twy_union_m, junction_zone_m):
                if not subtract.is_empty:
                    clipped_bp = clipped_bp.difference(subtract)
        except Exception:
            continue
        if clipped_bp.is_empty:
            continue
        # Update for Phase D use (complex building transitions)
        building_polys_m[bi] = clipped_bp
        # Handle MultiPolygon (clipping can split the footprint)
        pieces = ([clipped_bp] if clipped_bp.geom_type == "Polygon"
                  else list(clipped_bp.geoms)
                  if hasattr(clipped_bp, "geoms") else [])
        for piece in pieces:
            if piece.is_empty or not hasattr(piece, "exterior"):
                continue
            if piece.area < 10.0:
                continue    # skip tiny slivers from clipping
            fp_ll = [to_ll(mx, my) for mx, my in piece.exterior.coords]
            if len(fp_ll) < 4:
                continue
            _emit_flat_poly(fp_ll, bldg_elevations[bi])
            emitted_flat_shapes.append((piece, bldg_elevations[bi]))

    # Recompute bldg_union_m now that C3 has clipped building pads
    try:
        bldg_union_m = shp_ops.unary_union(
            [bp for bp in building_polys_m if not bp.is_empty]
        ) if building_polys_m else shp_geom.Polygon()
    except Exception:
        bldg_union_m = shp_geom.Polygon()

    # ── Coverage sweep: fill remaining uncovered paved areas ──────
    # Any area within runway/taxiway/apron unions that hasn't been
    # covered by the specific Phase C emissions gets a flat patch at
    # the CIFP surface elevation.  This ensures every bit of pavement
    # has an elevation shape, preventing DEM artifacts from distorting
    # the surface.
    n_coverage_fill = 0
    try:
        paved_area = shp_geom.Polygon()
        paved_parts = []
        if not rwy_union_m.is_empty:
            paved_parts.append(rwy_union_m)
        if not twy_union_m.is_empty:
            paved_parts.append(twy_union_m)
        if not apron_union_m.is_empty:
            paved_parts.append(apron_union_m)
        if paved_parts:
            paved_area = shp_ops.unary_union(paved_parts)

        if not paved_area.is_empty:
            # Compute what's already covered (or will be by Phase D)
            covered_parts = [rwy_union_m]  # runways always covered
            if not terminal_zone_m.is_empty:
                covered_parts.append(terminal_zone_m)
            if not junction_zone_m.is_empty:
                covered_parts.append(junction_zone_m)
            for s_poly, _ in emitted_flat_shapes:
                if not s_poly.is_empty:
                    covered_parts.append(s_poly)
            for q in emitted_twy_quads_m:
                if not q.is_empty:
                    covered_parts.append(q)
            covered = shp_ops.unary_union(
                [c for c in covered_parts if not c.is_empty])

            uncovered = paved_area.difference(covered.buffer(1.0))
            if not uncovered.is_empty:
                uncovered = uncovered.simplify(
                    1.0, preserve_topology=True)
                unc_polys = (list(uncovered.geoms)
                             if hasattr(uncovered, "geoms")
                             else [uncovered])
                for up in unc_polys:
                    if up.is_empty or up.area < 10.0:
                        continue
                    if not hasattr(up, "exterior"):
                        continue
                    cx_u, cy_u = up.centroid.x, up.centroid.y
                    u_elev = _cifp_surface_elevation(
                        cx_u, cy_u, max_search=500.0)
                    if u_elev is None:
                        ulon, ulat = to_ll(cx_u, cy_u)
                        try:
                            u_elev = tile.dem.alt(
                                (ulon - tile.lon, ulat - tile.lat))
                        except Exception:
                            continue
                    ring_ll = [to_ll(mx, my)
                               for mx, my in up.exterior.coords]
                    _emit_flat_poly(ring_ll, round(u_elev, 1))
                    emitted_flat_shapes.append(
                        (up, round(u_elev, 1)))
                    n_coverage_fill += 1
    except Exception:
        pass

    UI.vprint(2, "    {}: {} coverage fill patches".format(
        icao, n_coverage_fill))

    # ════════════════════════════════════════════════════════════════
    # PHASE C4: Detect elevation discontinuities between adjacent
    # flat shapes and create transition strips for Phase D.
    # ════════════════════════════════════════════════════════════════
    # Adjacent flat shapes at different elevations create vertical
    # cliffs that violate FAA grade limits.  For each pair of nearby
    # shapes with excessive grade, we create a transition strip that
    # Phase D will triangulate with grade-constrained relaxation.
    #
    # Example: terminal at 20.5m touching apron at 23.7m needs a
    # graded transition zone ≥ (3.2m / 0.015) ≈ 213m wide.
    # In practice, the complex apron detection in C2 handles most
    # cases (the apron itself becomes triangulated).  C4 catches
    # remaining flat-to-flat interfaces.

    transition_strip_parts = []

    if len(emitted_flat_shapes) >= 2:
        from shapely.strtree import STRtree
        fs_polys = [s[0] for s in emitted_flat_shapes]
        fs_elevs = [s[1] for s in emitted_flat_shapes]
        tree = STRtree(fs_polys)

        checked_pairs = set()
        for i, (poly_i, elev_i) in enumerate(emitted_flat_shapes):
            # Query shapes within MAX_TRANSITION_GAP of this one
            try:
                search_buf = poly_i.buffer(MAX_TRANSITION_GAP)
                candidates = tree.query(search_buf)
            except Exception:
                continue
            for j_idx in candidates:
                j = int(j_idx)
                if j <= i:
                    continue
                pair_key = (i, j) if i < j else (j, i)
                if pair_key in checked_pairs:
                    continue
                checked_pairs.add(pair_key)

                poly_j = fs_polys[j]
                elev_j = fs_elevs[j]
                delta_e = abs(elev_i - elev_j)
                if delta_e < 0.2:
                    continue   # negligible difference

                # Compute gap distance between the two shapes
                try:
                    gap_dist = poly_i.distance(poly_j)
                except Exception:
                    continue
                if gap_dist > MAX_TRANSITION_GAP:
                    continue

                # Required transition distance for grade compliance
                required_dist = delta_e / MAX_SURFACE_GRADE
                actual_dist = gap_dist + 1.0  # +1 for tolerance

                if actual_dist >= required_dist:
                    continue   # grade is already OK

                # Create transition strip: buffer the shared boundary
                # region by the needed transition width.
                strip_width = max(
                    TRANSITION_STRIP_MIN_WIDTH,
                    (required_dist - actual_dist) / 2.0
                )
                try:
                    # Find the shared boundary region
                    near_zone_i = poly_i.buffer(strip_width)
                    near_zone_j = poly_j.buffer(strip_width)
                    strip = near_zone_i.intersection(near_zone_j)
                    if strip.is_empty or strip.area < 5.0:
                        continue
                    # Subtract all existing emitted shapes
                    for s_poly, _ in emitted_flat_shapes:
                        if not s_poly.is_empty:
                            strip = strip.difference(s_poly)
                    # Also subtract runways, taxiways, junctions
                    for subtract in (rwy_union_m, twy_union_m,
                                     junction_zone_m):
                        if not subtract.is_empty:
                            strip = strip.difference(subtract)
                    if not strip.is_empty and strip.area > 5.0:
                        transition_strip_parts.append(strip)
                except Exception:
                    pass

    UI.vprint(2, "    {}: {} transition strips between flat shapes"
              .format(icao, len(transition_strip_parts)))

    # Recompute apron_union_m to include only flat aprons (complex
    # ones are going to Phase D, not subtracting from triangle zone)
    flat_apron_polys_m = [s[0] for s in emitted_flat_shapes
                          if any(s[0].equals(ap)
                                 for ap in apron_polys_m
                                 if not ap.is_empty)]
    # Keep the original apron_union_m for subtract operations — it
    # includes all aprons, both flat and complex.

    # ════════════════════════════════════════════════════════════════
    # PHASE D: Triangulate junction zones + sloping taxiway curves
    # ════════════════════════════════════════════════════════════════
    # Junctions are areas where two paved features physically overlap:
    #   * taxiway buffer ∩ runway polygon
    #   * taxiway buffer ∩ another taxiway buffer
    # Sloped taxiway curves (sent from Phase C1) are added too — their
    # per-segment rectangle emission would fan and overlap on curves,
    # so they get triangulated instead.
    #
    # Each connected component of the combined triangle zone is fed
    # through O4_Surface_Mesh.adaptive_triangulate, which produces the
    # smallest set of grade-compliant DEM-faithful triangles for that
    # component.  The legacy code Delaunay-triangulated a single
    # densified point cloud across the entire zone (1394 triangles at
    # SPJC); the new per-component approach typically yields 5-30
    # triangles per junction.
    #
    # Anchors come from:
    #   * taxiway centerline crossings → twy_centerlines elevations
    #   * runway-touch boundary points → linear CIFP interpolation
    #     along the touched runway centerline (this is the ONE place
    #     non-runway features legitimately read CIFP elevations: at
    #     the actual runway boundary, where the junction must meet
    #     the runway exactly)
    #   * adjacent emitted-flat-shape edges (taxiways, aprons,
    #     buildings already emitted) — anchors them to the
    #     surrounding mesh
    # Per STATUS.md elevation rule, NO global CIFP surface model
    # queries (no _cifp_surface_elevation) are used for non-touch
    # vertices.
    triangle_zone_parts = []

    if not junction_zone_m.is_empty:
        jz_clipped = junction_zone_m
        if not rwy_union_m.is_empty:
            try:
                jz_clipped = jz_clipped.difference(rwy_union_m)
            except Exception:
                pass
        if not jz_clipped.is_empty:
            triangle_zone_parts.append(jz_clipped)

    # Sloping taxiway buffers from Phase C1 (curves where the
    # per-segment rectangle emission would fan and overlap).
    for part in sloped_twy_parts_m:
        if not part.is_empty:
            triangle_zone_parts.append(part)

    UI.vprint(2, "    {}: {} complex apron parts, {} transition strips"
              " → triangle zone".format(
                  icao, len(complex_apron_parts_m),
                  len(transition_strip_parts)))

    triangle_zone = shp_geom.Polygon()
    if triangle_zone_parts:
        try:
            triangle_zone = shp_ops.unary_union(triangle_zone_parts)
        except Exception:
            triangle_zone = shp_geom.Polygon()

    # Junction polygons are NOT simplified.  Their boundaries are
    # SHARED edges with adjacent flat shapes (taxiway flats, apron
    # pieces, building pads) and ANY vertex shift moves the
    # triangulation off the shared edge, causing cross-category
    # overlap.  The polygon's vertex count is a small constant per
    # junction (most junctions have 4-12 vertices anyway), so
    # passing them all through to adaptive_triangulate is fine.
    JUNCTION_ANCHOR_GAP_M = 30.0

    def _runway_touch_elev(x_m, y_m):
        """Return the CIFP-interpolated runway elevation at (x, y)
        if the point lies on or near a runway boundary.  Used only
        for anchors at the actual runway touch line — never as a
        global query.
        """
        try:
            lon, lat = to_ll(x_m, y_m)
            return _project_point_onto_runway(
                lat, lon, runway_pairs, max_dist=10.0)
        except Exception:
            return None

    def _junction_anchors(zone_poly):
        """Collect (x, y, z) anchors for a single junction zone
        component.  Sources (in priority order):
          1. Taxiway centerlines that pass through or stub into it
          2. Runway touch points along the zone's intersection with
             a runway boundary
          3. Edges of already-emitted flat shapes (taxiway flats,
             apron pieces, building pads) within ~2 m of the zone
        """
        raw = []
        # 1. Taxiway crossings.
        TWY_REACH = TAXIWAY_BUFFER_WIDTH * 2
        zone_buf = zone_poly.buffer(TWY_REACH)
        for ti2, cl_m2 in enumerate(twy_centerlines_m):
            if cl_m2 is None:
                continue
            entry2 = twy_centerlines[ti2]
            if entry2 is None:
                continue
            cl_elevs = entry2[1]
            if not cl_elevs or len(cl_elevs) != len(cl_m2):
                continue
            try:
                ls2 = shp_geom.LineString(cl_m2)
                if ls2.distance(zone_poly) > TWY_REACH:
                    continue
                inside = ls2.intersection(zone_buf)
            except Exception:
                continue
            if inside.is_empty:
                continue
            comps = (list(inside.geoms)
                     if hasattr(inside, "geoms") else [inside])
            for comp in comps:
                if not hasattr(comp, "coords"):
                    continue
                coords = list(comp.coords)
                if not coords:
                    continue
                # Use both endpoints of each crossing segment so
                # adaptive_triangulate sees the full taxiway profile
                # through the junction.
                for c in (coords[0], coords[-1]):
                    px, py = c[0], c[1]
                    ie = _interpolate_elevation_along_centerline(
                        (px, py), cl_m2, cl_elevs)
                    raw.append((px, py, float(ie)))
        # 2. Runway touch points: sample along the zone's boundary
        # wherever it lies inside a runway polygon.
        if not rwy_union_m.is_empty:
            try:
                touch = zone_poly.boundary.intersection(rwy_union_m)
            except Exception:
                touch = None
            if touch is not None and not touch.is_empty:
                comps = (list(touch.geoms)
                         if hasattr(touch, "geoms") else [touch])
                for comp in comps:
                    if not hasattr(comp, "coords"):
                        continue
                    length = comp.length if hasattr(comp, "length") else 0.0
                    if length < 1.0:
                        for c in comp.coords:
                            elev = _runway_touch_elev(c[0], c[1])
                            if elev is not None:
                                raw.append((c[0], c[1], float(elev)))
                        continue
                    n = max(2, int(length / 30.0) + 1)
                    for k in range(n):
                        t = k / (n - 1) if n > 1 else 0.5
                        try:
                            p = comp.interpolate(t, normalized=True)
                            elev = _runway_touch_elev(p.x, p.y)
                            if elev is not None:
                                raw.append((p.x, p.y, float(elev)))
                        except Exception:
                            pass
        # 3. Adjacent emitted flat shape edges.  Dense anchors along
        # the shared boundary pin junction vertices to neighbour
        # elevations and prevent the triangulation from drifting
        # over the boundary into the flat shape.
        if emitted_flat_shapes:
            for sp, se in emitted_flat_shapes:
                if sp.is_empty:
                    continue
                try:
                    if sp.distance(zone_poly) > 2.0:
                        continue
                    shared = sp.boundary.intersection(zone_poly.buffer(2.0))
                except Exception:
                    continue
                if shared is None or shared.is_empty:
                    continue
                lines = []
                if hasattr(shared, "geoms"):
                    for g in shared.geoms:
                        if hasattr(g, "coords"):
                            lines.append(g)
                elif hasattr(shared, "coords"):
                    lines.append(shared)
                for ln in lines:
                    length = ln.length
                    if length < 1.0:
                        p = ln.interpolate(0.5, normalized=True)
                        raw.append((p.x, p.y, float(se)))
                        continue
                    for t in (0.0, 0.5, 1.0):
                        p = ln.interpolate(t, normalized=True)
                        raw.append((p.x, p.y, float(se)))
        return _sparsify_anchors(raw, JUNCTION_ANCHOR_GAP_M)

    total_tri = 0
    emitted_triangle_polys_m = []   # actual triangle polygons in
                                    # meter space, fed into
                                    # all_emitted_parts_m so Phase
                                    # E/F's _emitted_union() sees
                                    # the precise triangulated
                                    # geometry (not the abstract
                                    # triangle_zone, which can drift)
    # Pre-compute the unions that will be subtracted from each
    # component (also re-subtracted after per-component simplify so
    # the simplified outline can't drift over an adjacent painted
    # shape).
    all_flat_polys = [s[0] for s in emitted_flat_shapes
                      if not s[0].is_empty]
    flat_all_union = shp_geom.Polygon()
    quad_union = shp_geom.Polygon()
    try:
        if all_flat_polys:
            flat_all_union = shp_ops.unary_union(all_flat_polys)
        if emitted_twy_quads_m:
            quad_union = shp_ops.unary_union(emitted_twy_quads_m)
    except Exception:
        pass

    if not triangle_zone.is_empty:
        # Subtract emitted flat shapes so triangles don't overlap
        # what's already painted by Phase C.
        try:
            if not flat_all_union.is_empty:
                triangle_zone = triangle_zone.difference(flat_all_union)
            if not quad_union.is_empty:
                triangle_zone = triangle_zone.difference(quad_union)
        except Exception:
            pass

        if not triangle_zone.is_empty:
            zone_components = (list(triangle_zone.geoms)
                               if hasattr(triangle_zone, "geoms")
                               else [triangle_zone])
            for comp in zone_components:
                if (comp.is_empty or not hasattr(comp, "exterior")
                        or comp.area < 5.0):
                    continue
                seed_comp = comp
                # Re-subtraction can split the component into a
                # MultiPolygon.  Iterate over each sub-piece.
                sub_pieces = (list(seed_comp.geoms)
                              if hasattr(seed_comp, "geoms")
                              else [seed_comp])
                for sub in sub_pieces:
                    if (sub.is_empty
                            or not hasattr(sub, "exterior")
                            or sub.area < 5.0):
                        continue
                    anchors = _junction_anchors(sub)
                    try:
                        tris = SM.adaptive_triangulate(
                            sub, anchors, _dem_at,
                            max_grade=MAX_SURFACE_GRADE,
                            fidelity_tol=1.0,
                            max_extra_points=50)
                    except Exception:
                        tris = []
                    for (v0, v1, v2) in tris:
                        # Skip triangles whose elevations all match
                        # DEM within DEM_VARIANCE_THRESHOLD — those
                        # are "flat" pieces where the underlying
                        # DEM is already correct, no patch needed.
                        # Matches the legacy's variance filter and
                        # cuts the triangle count dramatically.
                        has_variance = False
                        for v in (v0, v1, v2):
                            dem_z = _dem_at(v[0], v[1])
                            if (dem_z is not None
                                    and abs(v[2] - dem_z)
                                    > DEM_VARIANCE_THRESHOLD):
                                has_variance = True
                                break
                        if not has_variance:
                            continue
                        lon0, lat0 = to_ll(v0[0], v0[1])
                        lon1, lat1 = to_ll(v1[0], v1[1])
                        lon2, lat2 = to_ll(v2[0], v2[1])
                        _emit_triangle(
                            lat0, lon0, round(v0[2], 1),
                            lat1, lon1, round(v1[2], 1),
                            lat2, lon2, round(v2[2], 1))
                        total_tri += 1
                        # Record the actual triangle for the
                        # downstream emit accumulator.
                        try:
                            tri_poly = shp_geom.Polygon(
                                [(v0[0], v0[1]),
                                 (v1[0], v1[1]),
                                 (v2[0], v2[1])])
                            if (tri_poly.is_valid
                                    and not tri_poly.is_empty):
                                emitted_triangle_polys_m.append(tri_poly)
                        except Exception:
                            pass
            if total_tri:
                UI.vprint(2, "    {}: {} junction triangles "
                          "(adaptive mesh)".format(icao, total_tri))

    n_twy_shapes = n_twy_flat + n_twy_sloped

    # ── Accumulate ALL emitted shapes so far (Phases C + D) ──────
    # Phase E (boundary band) and Phase F (drainage) call
    # _emitted_union() to subtract everything that's already been
    # painted.  We feed it the actual triangle polygons (not the
    # abstract triangle_zone) so the union is precise: drainage and
    # boundary bands can't accidentally land on top of a Phase D
    # triangle just because the zone calculation drifted.
    for s_poly, _ in emitted_flat_shapes:
        if not s_poly.is_empty:
            all_emitted_parts_m.append(s_poly)
    for q in emitted_twy_quads_m:
        if not q.is_empty:
            all_emitted_parts_m.append(q)
    for tri_poly in emitted_triangle_polys_m:
        if not tri_poly.is_empty:
            all_emitted_parts_m.append(tri_poly)
    if not rwy_union_m.is_empty:
        all_emitted_parts_m.append(rwy_union_m)

    UI.vprint(2, "    {}: accumulated {} emitted parts before Phase E"
              .format(icao, len(all_emitted_parts_m)))

    # ════════════════════════════════════════════════════════════════
    # PHASE E: Boundary band + tunnel portals
    # ════════════════════════════════════════════════════════════════
    # Simplified approach:
    #   1. A 15m-wide flat band following the airport boundary inward
    #      at the airport surface elevation.  This represents the
    #      perimeter access road that exists at every real airport.
    #      Ortho4XP smooths the transition to outside DEM automatically.
    #   2. Tunnel portal handling: where a road enters a tunnel under
    #      the airport, emit grade-down rects from airport surface to
    #      tunnel depth (~6m below grade).

    BOUNDARY_BAND_WIDTH = 15.0      # meters inward from boundary
    TUNNEL_GRADE_DOWN_LENGTH = 120.0  # meters to descend to tunnel depth
    TUNNEL_DEPTH_DEFAULT = 6.0      # meters below airport surface
    TUNNEL_RECT_WIDTH = 22.0        # meters — tunnel road rect width
    RETAINING_WALL_WIDTH = 6.0      # meters — wall strip width
    GRADE_STEP_LENGTH = 20.0        # meters per grade-down segment
    BAND_SEG_LENGTH = 40.0          # meters — boundary band segment length

    n_blend = 0
    blend_polys_m = []

    # Build the airport footprint for Phase E and F.
    if not osm_boundary_m.is_empty:
        airport_footprint = osm_boundary_m
    else:
        all_airport_parts = []
        if not rwy_union_m.is_empty:
            all_airport_parts.append(rwy_union_m)
        if not twy_union_m.is_empty:
            all_airport_parts.append(twy_union_m)
        if not apron_union_m.is_empty:
            all_airport_parts.append(apron_union_m)
        airport_footprint = (shp_ops.unary_union(all_airport_parts)
                             if all_airport_parts else shp_geom.Polygon())

    try:
        if not airport_footprint.is_empty and airport_footprint.area > 100:
            emitted_union = _emitted_union()

            # ── E1: Boundary band (perimeter access road) ──────────
            # Overlap-free band via polygon difference, with per-piece
            # elevation sampling from CIFP surface model. Each piece gets
            # the elevation at its centroid, which naturally follows the
            # airport slope.
            try:
                inner_ring = airport_footprint.buffer(
                    -BOUNDARY_BAND_WIDTH)
                if not inner_ring.is_empty:
                    band = airport_footprint.difference(inner_ring)
                    # Subtract already-emitted shapes to avoid overlaps
                    if not emitted_union.is_empty:
                        band = band.difference(emitted_union)
                    if not band.is_empty:
                        # Simplify lightly to reduce vertex count
                        band = band.simplify(
                            1.0, preserve_topology=True)
                        band_polys = (
                            list(band.geoms) if hasattr(band, "geoms")
                            else [band])
                        for bp in band_polys:
                            if (bp.is_empty or bp.area < 20.0
                                    or not hasattr(bp, "exterior")):
                                continue
                            # Sample elevation at centroid
                            cx_b, cy_b = bp.centroid.x, bp.centroid.y
                            b_elev = _cifp_surface_elevation(
                                cx_b, cy_b, max_search=500.0)
                            if b_elev is None:
                                blon, blat = to_ll(cx_b, cy_b)
                                try:
                                    b_elev = tile.dem.alt(
                                        (blon - tile.lon,
                                         blat - tile.lat))
                                except Exception:
                                    continue
                            ring_ll = [to_ll(mx, my)
                                       for mx, my in
                                       bp.exterior.coords]
                            _emit_flat_poly(
                                ring_ll, round(b_elev, 1))
                            blend_polys_m.append(bp)
                            n_blend += 1
            except Exception:
                pass

            # ── E2: Tunnel portals ─────────────────────────────────
            # Where a tunnel road crosses under or near the airport
            # boundary, emit grade-down rects and retaining wall flanks
            # at each portal (entrance/exit).
            if road_lines_m:
                from shapely.strtree import STRtree as _STRtree
                boundary_line = airport_footprint.boundary
                if boundary_line.geom_type == "MultiLineString":
                    boundary_line = max(
                        boundary_line.geoms, key=lambda g: g.length)

                for r_ls, r_hwy, r_tunnel, r_bridge in road_lines_m:
                    if not r_tunnel:
                        continue
                    # Check if this tunnel road is near the boundary
                    if boundary_line.distance(r_ls) > 100.0:
                        continue

                    # Find where the tunnel crosses the boundary
                    try:
                        xing = r_ls.intersection(boundary_line)
                    except Exception:
                        xing = None

                    # Also check endpoints — tunnel may start/end at
                    # boundary rather than crossing it
                    r_coords = list(r_ls.coords)
                    portal_points = []

                    if xing is not None and not xing.is_empty:
                        if hasattr(xing, "geoms"):
                            for g in xing.geoms:
                                if hasattr(g, "x"):
                                    portal_points.append((g.x, g.y))
                        elif hasattr(xing, "x"):
                            portal_points.append((xing.x, xing.y))

                    # If no crossing found, check endpoints near boundary
                    if not portal_points:
                        for rc in [r_coords[0], r_coords[-1]]:
                            pt = shp_geom.Point(rc)
                            if boundary_line.distance(pt) < 80.0:
                                portal_points.append(rc)

                    if not portal_points:
                        continue

                    for px, py in portal_points:
                        # Compute local road direction at this portal
                        try:
                            # Interpolate direction along road ±15m from portal
                            portal_pt = shp_geom.Point(px, py)
                            d_along = r_ls.project(portal_pt)
                            p1 = r_ls.interpolate(max(0, d_along - 15))
                            p2 = r_ls.interpolate(
                                min(r_ls.length, d_along + 15))
                            rdx = p2.x - p1.x
                            rdy = p2.y - p1.y
                            r_dir_len = sqrt(rdx * rdx + rdy * rdy)
                            if r_dir_len < 0.1:
                                continue
                            road_dir = (rdx / r_dir_len, rdy / r_dir_len)
                        except Exception:
                            continue

                        # Airport elevation at portal
                        apt_elev = _cifp_surface_elevation(
                            px, py, max_search=500.0)
                        if apt_elev is None:
                            plon, plat = to_ll(px, py)
                            try:
                                apt_elev = tile.dem.alt(
                                    (plon - tile.lon,
                                     plat - tile.lat))
                            except Exception:
                                continue

                        # Direction from portal outward (away from
                        # airport center)
                        apt_cx = airport_footprint.centroid.x
                        apt_cy = airport_footprint.centroid.y
                        away_x = px - apt_cx
                        away_y = py - apt_cy
                        away_len = sqrt(
                            away_x ** 2 + away_y ** 2)
                        if away_len < 0.1:
                            continue
                        out_dir = (away_x / away_len,
                                   away_y / away_len)

                        # Choose road direction that aligns with
                        # outward direction
                        dot_out = (road_dir[0] * out_dir[0]
                                   + road_dir[1] * out_dir[1])
                        if dot_out < 0:
                            eff_dir = (-road_dir[0], -road_dir[1])
                        else:
                            eff_dir = road_dir

                        # Log tunnel portal detection
                        UI.vprint(2, "      Tunnel portal at ({:.1f}, {:.1f})".format(px, py))

                        # Grade-down rects from portal outward
                        n_steps = max(2, int(
                            TUNNEL_GRADE_DOWN_LENGTH
                            / GRADE_STEP_LENGTH))
                        step_len = (TUNNEL_GRADE_DOWN_LENGTH
                                    / n_steps)

                        for si in range(n_steps):
                            frac_s = si / n_steps
                            frac_e = (si + 1) / n_steps
                            e_s = apt_elev - (
                                TUNNEL_DEPTH_DEFAULT * frac_s)
                            e_e = apt_elev - (
                                TUNNEL_DEPTH_DEFAULT * frac_e)

                            sx = px + eff_dir[0] * step_len * si
                            sy = py + eff_dir[1] * step_len * si
                            ex = px + eff_dir[0] * step_len * (si + 1)
                            ey = py + eff_dir[1] * step_len * (si + 1)

                            seg_mid = shp_geom.Point(
                                (sx + ex) / 2, (sy + ey) / 2)
                            # Use a small inward buffer for edge cases
                            apt_buffer = airport_footprint.buffer(-2.0) if not airport_footprint.is_empty else airport_footprint
                            if not apt_buffer.is_empty and apt_buffer.contains(seg_mid):
                                continue

                            slon, slat = to_ll(sx, sy)
                            elon, elat = to_ll(ex, ey)
                            try:
                                tc = runway_corners(
                                    slat, slon, elat, elon,
                                    TUNNEL_RECT_WIDTH)
                                if tc is None:
                                    continue
                                tp = shp_geom.Polygon(
                                    [to_m(c[1], c[0]) for c in tc]
                                    + [to_m(tc[0][1], tc[0][0])])
                                if not tp.is_valid:
                                    tp = tp.buffer(0)
                                if tp.is_empty:
                                    continue
                            except Exception:
                                continue

                            _emit_sloped_rect(
                                slat, slon, round(e_s, 1),
                                elat, elon, round(e_e, 1),
                                TUNNEL_RECT_WIDTH
                            )
                            blend_polys_m.append(tp)
                            n_blend += 1

                        # Retaining wall flanks
                        perp_nx = -eff_dir[1]
                        perp_ny = eff_dir[0]
                        for side in (1.0, -1.0):
                            wall_off = (TUNNEL_RECT_WIDTH / 2.0
                                        + RETAINING_WALL_WIDTH / 2.0)
                            w_sx = px + side * perp_nx * wall_off
                            w_sy = py + side * perp_ny * wall_off
                            w_ex = (px
                                    + eff_dir[0]
                                    * TUNNEL_GRADE_DOWN_LENGTH
                                    + side * perp_nx * wall_off)
                            w_ey = (py
                                    + eff_dir[1]
                                    * TUNNEL_GRADE_DOWN_LENGTH
                                    + side * perp_ny * wall_off)
                            wm = shp_geom.Point(
                                (w_sx + w_ex) / 2,
                                (w_sy + w_ey) / 2)
                            # Use a small inward buffer for edge cases
                            apt_buffer = airport_footprint.buffer(-2.0) if not airport_footprint.is_empty else airport_footprint
                            if not apt_buffer.is_empty and apt_buffer.contains(wm):
                                continue
                            wslon, wslat = to_ll(w_sx, w_sy)
                            welon, welat = to_ll(w_ex, w_ey)
                            try:
                                wc = runway_corners(
                                    wslat, wslon, welat, welon,
                                    RETAINING_WALL_WIDTH)
                                if wc is None:
                                    continue
                                wp = shp_geom.Polygon(
                                    [to_m(c[1], c[0]) for c in wc]
                                    + [to_m(wc[0][1], wc[0][0])])
                                if not wp.is_valid:
                                    wp = wp.buffer(0)
                                if wp.is_empty:
                                    continue
                            except Exception:
                                continue
                            _emit_sloped_rect(
                                wslat, wslon, round(apt_elev, 1),
                                welat, welon, round(apt_elev, 1),
                                RETAINING_WALL_WIDTH
                            )
                            blend_polys_m.append(wp)
                            n_blend += 1

    except Exception:
        pass

    # Add Phase E shapes to emitted accumulator
    all_emitted_parts_m.extend(blend_polys_m)

    # ════════════════════════════════════════════════════════════════
    # PHASE F: Drainage zones
    # ════════════════════════════════════════════════════════════════
    # Two approaches based on the elevation variance of surrounding
    # pavement:
    #
    # Type A — Flat infield (surrounding elevations differ < 0.5m):
    #   Emit a polygon matching the infield shape but inset to ~15%
    #   of its size, set 2m below surrounding pavement.  This creates
    #   a catch basin; Ortho4XP smooths the transition from the
    #   pavement edges down to the depressed polygon automatically.
    #
    # Type B — Sloped infield (surrounding elevations differ ≥ 0.5m):
    #   Emit a sloping rectangle ~1m wide and ~70% of the infield's
    #   length, oriented along the slope, creating a long linear
    #   ditch that drains downhill.

    DRAIN_INSET_FRACTION = 0.15     # Type A: 15% of area (by linear scale)
    DRAIN_FLAT_DEPTH = 2.0          # Type A: meters below surrounding
    DRAIN_SLOPE_THRESHOLD = 0.5     # meters — switch from Type A to B
    DRAIN_DITCH_WIDTH = 1.0         # Type B: ditch width in meters
    DRAIN_DITCH_LENGTH_FRAC = 0.70  # Type B: 70% of infield length

    n_drainage = 0
    try:
        if not airport_footprint.is_empty and airport_footprint.area > 100:
            full_emitted = _emitted_union()
            enclosure = (airport_footprint if not osm_boundary_m.is_empty
                         else airport_footprint.convex_hull)

            infield = enclosure
            if not full_emitted.is_empty:
                infield = enclosure.difference(
                    full_emitted.buffer(DRAINAGE_EDGE_BUFFER))

            if not infield.is_empty:
                infield_polys = (list(infield.geoms)
                                 if hasattr(infield, "geoms")
                                 else [infield])
                for ip in infield_polys:
                    if ip.is_empty or ip.area < DRAINAGE_MIN_AREA:
                        continue
                    if not hasattr(ip, "exterior"):
                        continue

                    # Only create drainage in areas fully enclosed by
                    # pavement (runways, taxiways, aprons).  Skip
                    # infields that border the airport perimeter or
                    # apron/terminal areas.
                    ip_boundary = ip.boundary
                    ip_boundary_len = ip_boundary.length
                    if ip_boundary_len < 10.0:
                        continue
                    # Drainage must be enclosed by any kind of
                    # pavement.  The legacy only counted runway and
                    # taxiway contact, which ruled out the interior
                    # grass islands of apron-rich airports (the
                    # usual drainage candidates).  With apt.dat as
                    # the source, taxiway_union is empty — all
                    # pavement is in the apron union — so we have
                    # to include apron contact too.
                    paved_contact = 0.0
                    try:
                        paved_parts = []
                        if not rwy_union_m.is_empty:
                            paved_parts.append(rwy_union_m)
                        if not twy_union_m.is_empty:
                            paved_parts.append(twy_union_m)
                        if not apron_union_m.is_empty:
                            paved_parts.append(apron_union_m)
                        if paved_parts:
                            paved_union = shp_ops.unary_union(paved_parts)
                            paved_buf = paved_union.buffer(
                                DRAINAGE_EDGE_BUFFER + 5.0)
                            contact = ip_boundary.intersection(paved_buf)
                            if not contact.is_empty:
                                paved_contact = (
                                    contact.length / ip_boundary_len)
                    except Exception:
                        pass
                    if paved_contact < 0.80:
                        continue  # not enclosed by pavement
                    # Also skip if touching buildings/terminals
                    try:
                        if not bldg_union_m.is_empty:
                            bldg_contact = ip.intersection(
                                bldg_union_m.buffer(5.0))
                            if (not bldg_contact.is_empty
                                    and bldg_contact.area > 0.05 * ip.area):
                                continue
                    except Exception:
                        pass

                    # Sample surrounding elevations from boundary
                    try:
                        bnd_dense = ip.exterior.segmentize(30.0)
                        bnd_coords = list(bnd_dense.coords)
                    except Exception:
                        continue
                    bnd_elevs = []
                    for bx, by in bnd_coords:
                        se = _cifp_surface_elevation(
                            bx, by, max_search=500.0)
                        if se is None:
                            slon, slat = to_ll(bx, by)
                            try:
                                se = tile.dem.alt(
                                    (slon - tile.lon, slat - tile.lat))
                            except Exception:
                                continue
                        bnd_elevs.append((bx, by, se))

                    if len(bnd_elevs) < 3:
                        continue

                    elevs_only = [e for _, _, e in bnd_elevs]
                    elev_range = max(elevs_only) - min(elevs_only)
                    avg_elev = sum(elevs_only) / len(elevs_only)

                    if elev_range < DRAIN_SLOPE_THRESHOLD:
                        # ── Type A: Flat infield — inset polygon ───
                        # Inset by shrinking toward centroid.
                        # Linear scale factor = sqrt(DRAIN_INSET_FRACTION)
                        # ≈ 0.387 → the polygon is ~39% of linear size
                        # → ~15% of area.
                        scale = sqrt(DRAIN_INSET_FRACTION)
                        try:
                            # Negative buffer to inset
                            # Compute inset distance from area/perimeter
                            inset_dist = (1.0 - scale) * sqrt(
                                ip.area / pi) * 0.5
                            inset_poly = ip.buffer(-inset_dist)
                            if inset_poly.is_empty:
                                continue
                            if hasattr(inset_poly, "geoms"):
                                inset_poly = max(
                                    inset_poly.geoms,
                                    key=lambda g: g.area)
                            if not hasattr(inset_poly, "exterior"):
                                continue
                            if inset_poly.area < 20.0:
                                continue
                        except Exception:
                            continue

                        # Emit the inset polygon at depressed elevation
                        low_elev = round(avg_elev - DRAIN_FLAT_DEPTH, 1)
                        try:
                            ring_ll = [to_ll(mx, my)
                                       for mx, my in
                                       inset_poly.exterior.coords]
                            _emit_flat_poly(ring_ll, low_elev)
                            n_drainage += 1
                        except Exception:
                            pass

                    else:
                        # ── Type B: Sloped infield — linear ditch ──
                        # Align ditch parallel to nearest runway or taxiway centerline.
                        # Fall back to highest→lowest if no centerlines found.

                        # Build reference lines from runways and taxiways
                        # Each entry: (LineString, unit_direction)
                        ref_lines = []
                        if runway_pairs:
                            for rp_d in runway_pairs:
                                da = rp_d.get("data_a")
                                db = rp_d.get("data_b")
                                if da and db:
                                    rwy_mx1, rwy_my1 = to_m(da["lon"], da["lat"])
                                    rwy_mx2, rwy_my2 = to_m(db["lon"], db["lat"])
                                    rdx = rwy_mx2 - rwy_mx1
                                    rdy = rwy_my2 - rwy_my1
                                    rlen = sqrt(rdx * rdx + rdy * rdy)
                                    if rlen > 1.0:
                                        ref_lines.append((
                                            shp_geom.LineString(
                                                [(rwy_mx1, rwy_my1),
                                                 (rwy_mx2, rwy_my2)]),
                                            (rdx / rlen, rdy / rlen)))

                        if twy_centerlines_m:
                            for tcl in twy_centerlines_m:
                                if tcl and len(tcl) >= 2:
                                    tdx = tcl[-1][0] - tcl[0][0]
                                    tdy = tcl[-1][1] - tcl[0][1]
                                    tlen = sqrt(tdx * tdx + tdy * tdy)
                                    if tlen > 1.0:
                                        ref_lines.append((
                                            shp_geom.LineString(tcl),
                                            (tdx / tlen, tdy / tlen)))

                        # Infield centroid
                        ic_x = ip.centroid.x
                        ic_y = ip.centroid.y

                        # Find best reference direction (nearest centerline)
                        best_dir = None
                        best_dist = float('inf')
                        if ref_lines:
                            centroid_pt = shp_geom.Point(ic_x, ic_y)
                            for line, udir in ref_lines:
                                d = line.distance(centroid_pt)
                                if d < best_dist:
                                    best_dist = d
                                    best_dir = udir

                        # Fallback: highest to lowest boundary point
                        if best_dir is None:
                            high_pt = max(bnd_elevs, key=lambda t: t[2])
                            low_pt = min(bnd_elevs, key=lambda t: t[2])
                            dx = low_pt[0] - high_pt[0]
                            dy = low_pt[1] - high_pt[1]
                            ditch_full_len = sqrt(dx * dx + dy * dy)
                            if ditch_full_len < 5.0:
                                continue
                            best_dir = (dx / ditch_full_len, dy / ditch_full_len)
                            high_e_val = high_pt[2]
                            low_e_val = low_pt[2]
                        else:
                            # Compute elevations in the best_dir direction
                            # Project boundary points onto this direction
                            proj_vals = []
                            for bx, by, be in bnd_elevs:
                                proj = (bx - ic_x) * best_dir[0] + (by - ic_y) * best_dir[1]
                                proj_vals.append((proj, be))
                            if proj_vals:
                                proj_vals.sort(key=lambda t: t[0])
                                high_e_val = proj_vals[-1][1]
                                low_e_val = proj_vals[0][1]
                            else:
                                high_e_val = avg_elev
                                low_e_val = avg_elev

                        # Use best_dir to position ditch
                        # Estimate infield length along the ditch direction
                        # by projecting all boundary coords onto best_dir
                        proj_extent = []
                        for bx, by, _ in bnd_elevs:
                            proj_extent.append(
                                (bx - ic_x) * best_dir[0]
                                + (by - ic_y) * best_dir[1])
                        if proj_extent:
                            ditch_full_len = max(proj_extent) - min(proj_extent)
                        else:
                            ditch_full_len = sqrt(ip.area)
                        if ditch_full_len < 5.0:
                            continue
                        ditch_len = ditch_full_len * DRAIN_DITCH_LENGTH_FRAC
                        half_len = ditch_len / 2.0
                        sx = ic_x - best_dir[0] * half_len
                        sy = ic_y - best_dir[1] * half_len
                        ex = ic_x + best_dir[0] * half_len
                        ey = ic_y + best_dir[1] * half_len

                        # Verify ditch is inside the infield
                        ditch_mid = shp_geom.Point(ic_x, ic_y)
                        if not ip.contains(ditch_mid):
                            continue

                        # Elevations: depressed 2m below surrounding
                        high_e = round(high_e_val - DRAIN_FLAT_DEPTH, 1)
                        low_e = round(low_e_val - DRAIN_FLAT_DEPTH, 1)

                        slon, slat = to_ll(sx, sy)
                        elon, elat = to_ll(ex, ey)
                        _emit_sloped_rect(
                            slat, slon, high_e,
                            elat, elon, low_e,
                            DRAIN_DITCH_WIDTH
                        )
                        n_drainage += 1

    except Exception:
        pass

    UI.vprint(1, "    {}: {} taxiway ({} flat, {} sloped)"
              " + {} apron ({} flat, {} complex)"
              " + {} building pads"
              " + {} junction triangles"
              " + {} boundary/road shapes + {} drainage low-points".format(
                  icao, n_twy_shapes, n_twy_flat, n_twy_sloped,
                  n_apron_flat + n_apron_complex,
                  n_apron_flat, n_apron_complex,
                  len(building_polys_m), total_tri,
                  n_blend, n_drainage))

    return (lines, nid)


# Legacy standalone functions kept for backward compatibility but
# deprecated in favor of the unified generate_airport_surface_patches()
# pipeline which uses the full "anchors → platforms → connectors" model.
# These functions use the old model (project onto nearest runway) and
# should NOT be called for new code paths.


def _compute_building_elevation(footprint, runway_pairs, tile):
    """DEPRECATED: Use the A5 elevation pipeline in generate_airport_surface_patches.

    Legacy function that computes building elevation by projecting onto the
    nearest CIFP runway.  Does not use platform anchors or surface model.
    """
    n = len(footprint)
    if n < 3:
        return 0.0
    if footprint[0] == footprint[-1]:
        n -= 1
    centroid_lon = sum(c[0] for c in footprint[:n]) / n
    centroid_lat = sum(c[1] for c in footprint[:n]) / n
    best_elev = _project_point_onto_runway(
        centroid_lat, centroid_lon, runway_pairs, max_dist=200.0
    )
    if best_elev is not None:
        return best_elev
    x = centroid_lon - tile.lon
    y = centroid_lat - tile.lat
    try:
        return tile.dem.alt((x, y))
    except Exception:
        return 0.0
