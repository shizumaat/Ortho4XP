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
    # The union into dico_apt_entry["apron"] is kept for backwards
    # compatibility (fallback consumers read only "apron").  For the
    # taxiway/apron split we ALSO stash the per-pavement objects with
    # their names in dico_apt_entry["_apt_pavements"] so Phase A3 can
    # classify each polygon individually before converting to meter
    # space.
    apron_polys = []
    pav_tuples = []   # list of (tile_relative_polygon, name)
    for pav in apt_data.pavements:
        p_tr = _to_tile_relative(pav.polygon)
        if p_tr is not None and not p_tr.is_empty:
            apron_polys.append(p_tr)
            pav_tuples.append((p_tr, pav.name))
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
    new_dico["_apt_pavements"] = pav_tuples

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
            # Extend the rectangle by the blast-pad / stopway length
            # at each end so the emitted footprint covers what
            # generate_patch_osm draws (runway surface + blast pads).
            end_a = (rwy.lat_a, rwy.lon_a)
            end_b = (rwy.lat_b, rwy.lon_b)
            if rwy.blast_a_m > 0.1:
                end_a = extend_point(
                    rwy.lat_b, rwy.lon_b, rwy.lat_a, rwy.lon_a,
                    rwy.blast_a_m)
            if rwy.blast_b_m > 0.1:
                end_b = extend_point(
                    rwy.lat_a, rwy.lon_a, rwy.lat_b, rwy.lon_b,
                    rwy.blast_b_m)
            corners = runway_corners(
                end_a[0], end_a[1],
                end_b[0], end_b[1],
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
    # an additional anchor on the other runway with the source
    # threshold's elevation.  This forces parallel runways with
    # offset thresholds (e.g. CYXY 14R/32L vs 14L/32R, where
    # 32R end at 701 m sits ~200 m off 14R/32L's interior) to
    # respect each other's hard CIFP elevations rather than
    # letting DEM-seeded interior samples drift to the upper
    # envelope and break grade for taxi crossings.
    MAX_CROSS_RUNWAY_LATERAL_M = 300.0
    auto_extra_anchors: dict = {}
    paired_list = [(da, dat_a, db, dat_b)
                   for da, dat_a, db, dat_b in runway_pairs
                   if db is not None and dat_b is not None]
    for ti, (da_t, dat_a_t, db_t, dat_b_t) in enumerate(paired_list):
        # Each pair has two thresholds; project each onto every
        # OTHER pair's centerline.
        for src_desig, src_data in (
                (da_t, dat_a_t), (db_t, dat_b_t)):
            for ri, (da_r, dat_a_r, db_r, dat_b_r) in enumerate(
                    paired_list):
                if ri == ti:
                    continue
                # Centerline from dat_a_r → dat_b_r in meters at
                # the pair's mid-latitude.
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
                # Skip if projection is outside the threshold span
                # (with a small inset so we don't double-anchor at
                # the receiver runway's own threshold).
                if t <= 0.05 or t >= 0.95:
                    continue
                # Perpendicular distance from src threshold to
                # receiver runway centerline.
                proj_x = t * rdx
                proj_y = t * rdy
                perp = sqrt((vx - proj_x) ** 2 + (vy - proj_y) ** 2)
                if perp > MAX_CROSS_RUNWAY_LATERAL_M:
                    continue
                # Build anchor at the projection lat/lon on the
                # receiver runway, with the SRC threshold's
                # elevation.  The receiver runway's solver will
                # then envelope-clamp itself around this anchor as
                # an additional hard constraint.
                p_lat = (dat_a_r["lat"]
                         + t * (dat_b_r["lat"] - dat_a_r["lat"]))
                p_lon = (dat_a_r["lon"]
                         + t * (dat_b_r["lon"] - dat_a_r["lon"]))
                key = (da_r, db_r)
                auto_extra_anchors.setdefault(key, []).append(
                    (p_lat, p_lon, src_data["elevation_m"]))
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

            # ── Wide-window smoothing of the DEM profile ─────────────
            # A real runway is a graded surface that approximates the
            # underlying terrain only on average — local DEM bumps
            # from buildings, vegetation or surface-model noise do
            # not belong on the runway profile.  Before the
            # grade-clamp pass we replace each interior DEM sample
            # with a moving average over a wide window so the final
            # elevation profile is a single long gentle slope rather
            # than a staircase of 1.5 % ramps in alternating
            # directions following every local DEM bump.
            elevs = [s[2] for s in sample_pts]
            anchored = [s[3] for s in sample_pts]
            n_samples = len(elevs)
            if n_samples >= 5:
                half_win = max(4, n_samples // 4)
                smoothed = list(elevs)
                for i in range(n_samples):
                    if anchored[i]:
                        continue
                    lo_w = max(0, i - half_win)
                    hi_w = min(n_samples, i + half_win + 1)
                    window = elevs[lo_w:hi_w]
                    smoothed[i] = sum(window) / len(window)
                elevs = smoothed

            # ── Envelope pre-clamp from every anchor ─────────────────
            # For each sample i compute the tightest allowed band
            # [lower, upper] imposed by every anchor at distance d,
            # using the rule |elev[i] - elev[anchor]| ≤ d × cap.
            # Clamping DEM to this envelope guarantees the profile
            # is anchor-consistent before any local smoothing runs,
            # and collapses the number of local-cap iterations to
            # near zero.  The remaining local pass only has to
            # reconcile adjacent DEM samples.
            cum_dist = [0.0]
            for i in range(1, n_samples):
                cum_dist.append(
                    cum_dist[-1]
                    + abs(fractions[i] - fractions[i - 1]) * phys_dist)
            for i in range(n_samples):
                if anchored[i]:
                    continue
                upper = float("inf")
                lower = float("-inf")
                for j in range(n_samples):
                    if not anchored[j]:
                        continue
                    d_ij = abs(cum_dist[i] - cum_dist[j])
                    upper = min(upper, elevs[j] + d_ij * MAX_RUNWAY_GRADE)
                    lower = max(lower, elevs[j] - d_ij * MAX_RUNWAY_GRADE)
                if upper < lower:
                    # Anchors are inconsistent (shouldn't happen with
                    # CIFP data but guard anyway) — use midpoint.
                    elevs[i] = (upper + lower) / 2.0
                elif elevs[i] > upper:
                    elevs[i] = upper
                elif elevs[i] < lower:
                    elevs[i] = lower

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
                any_change = False
                for _it in range(GRADE_RELAX_ITERATIONS):
                    changed = False

                    # Boundary constraint at the blast-pad / runway
                    # interface: the flat blast pad effectively has
                    # grade 0 on the outside of the anchor, so at
                    # sample[0] the grade change from 0 into the
                    # first interior segment must not exceed
                    # MAX_RUNWAY_GRADE_CHANGE_PER_M × L.  The anchor
                    # itself cannot move, so we enforce the rule by
                    # clamping elevs[1] toward elevs[0].  Same at
                    # the far end.
                    if (len(elevs) >= 2 and anchored[0]
                            and not anchored[1]):
                        lr0 = _seg_len(0)
                        max_dg0 = (MAX_RUNWAY_GRADE_CHANGE_PER_M
                                   * lr0)
                        max_delta0 = max_dg0 * lr0
                        target_hi = elevs[0] + max_delta0
                        target_lo = elevs[0] - max_delta0
                        if elevs[1] > target_hi:
                            elevs[1] = target_hi
                            changed = True
                            any_change = True
                        elif elevs[1] < target_lo:
                            elevs[1] = target_lo
                            changed = True
                            any_change = True
                    if (len(elevs) >= 2 and anchored[-1]
                            and not anchored[-2]):
                        llN = _seg_len(len(elevs) - 2)
                        max_dgN = (MAX_RUNWAY_GRADE_CHANGE_PER_M
                                   * llN)
                        max_deltaN = max_dgN * llN
                        target_hi = elevs[-1] + max_deltaN
                        target_lo = elevs[-1] - max_deltaN
                        if elevs[-2] > target_hi:
                            elevs[-2] = target_hi
                            changed = True
                            any_change = True
                        elif elevs[-2] < target_lo:
                            elevs[-2] = target_lo
                            changed = True
                            any_change = True

                    for i in range(1, len(elevs) - 1):
                        if anchored[i]:
                            continue
                        ll = _seg_len(i - 1)
                        lr = _seg_len(i)
                        if ll < 0.1 or lr < 0.1:
                            continue
                        g_left = (elevs[i] - elevs[i - 1]) / ll
                        g_right = (elevs[i + 1] - elevs[i]) / lr
                        max_dg = (MAX_RUNWAY_GRADE_CHANGE_PER_M
                                  * ((ll + lr) / 2.0))
                        dg = g_right - g_left
                        if abs(dg) <= max_dg:
                            continue
                        target_dg = max_dg if dg > 0 else -max_dg
                        denom = 1.0 / lr + 1.0 / ll
                        new_e = (elevs[i + 1] / lr
                                 + elevs[i - 1] / ll
                                 - target_dg) / denom
                        if abs(new_e - elevs[i]) > 0.001:
                            elevs[i] = new_e
                            changed = True
                            any_change = True
                    if not changed:
                        break
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


def _emit_pavement_strip_model(
    apt_twy_polys_m, apron_polys_m, rwy_union_m,
    to_ll, to_m, apt_data,
    _emit_flat_poly, _emit_sloped_rect, _emit_triangle,
    icao,
):
    """Zero-elevation preview of the apron-deformation model.

    Calls :func:`O4_Pavement_Strips.decompose_pavement` (which now
    returns a flat tuple of :class:`Shape`) and emits:

    * each taxi shape as one flat sloped rect along its axis at
      the shape's representative width;
    * each apron connected component as a flat N-gon of its
      polygon exterior.

    There are no junctions, no wedges, no triangle zones yet —
    that's the job of the elevation solver + triangulated emitter
    which will replace this preview.
    """
    import O4_Pavement_Strips as _PS
    import O4_Pavement_Classifier as _PC
    from shapely.ops import transform as _shp_transform, unary_union \
        as _unary_union

    PREVIEW_ELEV = 0.0
    # Any apt.dat pavement polygon whose area sits >this fraction
    # inside the runway union is assumed to BE the runway surface
    # and is dropped so we don't emit it as a huge apron under the
    # runway segments.
    RUNWAY_OVERLAP_DROP_FRAC = 0.5

    # Prefer the RAW apt.dat pavements over the Phase-A3 rectifiable
    # subset: my strip decomposition does its own mega-polygon
    # splitting and wants the original polygons.  Fall back to the
    # pre-decomposed polys when apt.dat isn't available.
    raw_taxi_m = []
    raw_apron_m = []
    if apt_data is not None and apt_data.pavements:
        def _ll_to_m(lon, lat, z=None):
            x, y = to_m(lon, lat)
            return (x, y) if z is None else (x, y, z)
        for pav in apt_data.pavements:
            if pav.polygon is None or pav.polygon.is_empty:
                continue
            try:
                pav_m = _shp_transform(_ll_to_m, pav.polygon)
            except Exception:
                continue
            if pav_m.is_empty:
                continue
            kind = _PC.classify_pavement_m(pav_m, pav.name or "")
            if kind.kind == "taxiway":
                raw_taxi_m.append(pav_m)
            else:
                raw_apron_m.append(pav_m)
    else:
        raw_taxi_m = list(apt_twy_polys_m or [])
        raw_apron_m = list(apron_polys_m or [])

    # Filter out polygons whose extent is mostly the runway surface,
    # then subtract the runway union from the rest.  The first
    # filter removes apt.dat polygons that ARE the runway (pavement
    # polygons matching the runway rect + shoulders).  The second
    # trims polygons that merely overlap a runway corner.
    def _runway_drop_and_trim(polys):
        out = []
        rwy_present = (rwy_union_m is not None
                       and not rwy_union_m.is_empty)
        for p in polys:
            if p is None or p.is_empty:
                continue
            if rwy_present:
                try:
                    inter_area = p.intersection(rwy_union_m).area
                except Exception:
                    inter_area = 0.0
                if p.area > 0 and inter_area / p.area \
                        >= RUNWAY_OVERLAP_DROP_FRAC:
                    continue   # polygon IS (mostly) a runway surface
                try:
                    q = p.difference(rwy_union_m)
                except Exception:
                    q = p
            else:
                q = p
            if q is None or q.is_empty:
                continue
            if q.geom_type == "Polygon":
                out.append(q)
            elif hasattr(q, "geoms"):
                for g in q.geoms:
                    if g.geom_type == "Polygon" and not g.is_empty:
                        out.append(g)
        return out

    safe_taxi = _runway_drop_and_trim(raw_taxi_m)
    safe_apron = _runway_drop_and_trim(raw_apron_m)

    # Compute runway bearings in meter-space so the trunk extractor
    # can prefer "parallel to runway" pairings at cross-junctions.
    # This is what makes V, L, A come out as single long rects
    # instead of zig-zags through their cross-connectors.
    rwy_bearings: list = []
    if apt_data is not None:
        for r in apt_data.runways:
            try:
                ax, ay = to_m(r.lon_a, r.lat_a)
                bx, by = to_m(r.lon_b, r.lat_b)
                dx, dy = bx - ax, by - ay
                import math as _math
                if _math.hypot(dx, dy) > 1.0:
                    # Compass bearing: 0 = +Y (north), 90 = +X (east).
                    bearing = _math.degrees(_math.atan2(dx, dy))
                    rwy_bearings.append(bearing % 180.0)
            except Exception:
                pass

    shapes = _PS.decompose_pavement(
        safe_taxi, safe_apron,
        preferred_bearings=rwy_bearings or None,
    )

    # Role classification → category per shape.
    from shapely.geometry import LineString as _LineString
    adjacencies = _PS.build_adjacency_graph(shapes)
    runway_cls: list = []
    if apt_data is not None:
        for r in apt_data.runways:
            try:
                ax, ay = to_m(r.lon_a, r.lat_a)
                bx, by = to_m(r.lon_b, r.lat_b)
                if ((bx - ax) ** 2 + (by - ay) ** 2) > 1.0:
                    runway_cls.append(_LineString([(ax, ay), (bx, by)]))
            except Exception:
                pass
    roles = _PS.classify_shape_roles(shapes, adjacencies, runway_cls)

    # ── Build taxi rectangles ────────────────────────────────────
    # Every taxi rect is a clean rectangle from (trimmed axis +
    # width).  10 m is trimmed off each short end regardless of
    # role — V, L, A, every stub, every cross-connector.  Long
    # sides of the rect touch adjacent apron without any buffer;
    # that's enforced by the apron carve below.
    #
    # Parallels (primary/secondary) subdivide into ~100 m chunks
    # so the elevation solver can apply a per-segment sloped
    # profile; stubs and cross-connectors emit as one rect.
    SEGMENT_LEN_M = 100.0
    GAP_M = 10.0
    MIN_TRIMMED_LEN_M = 10.0
    from shapely.ops import substring as _substring
    from shapely.geometry import Polygon as _Polygon

    def _subdivide_axis(axis, seg_len):
        total = axis.length
        if total <= seg_len:
            cc = list(axis.coords)
            return [(cc[0], cc[-1])] if len(cc) >= 2 else []
        n = max(1, int(round(total / seg_len)))
        step = total / n
        pts = []
        for i in range(n + 1):
            p = axis.interpolate(i * step)
            pts.append((p.x, p.y))
        return list(zip(pts[:-1], pts[1:]))

    def _trim_axis(axis, trim_m):
        total = axis.length
        if total <= 2.0 * trim_m + MIN_TRIMMED_LEN_M:
            return None
        try:
            return _substring(axis, trim_m, total - trim_m)
        except Exception:
            return axis

    def _rect_between(p1, p2, width):
        """Axis-aligned rectangle (oriented to the p1→p2 direction)
        of the given width, centered on the line from p1 to p2.
        Returns a Polygon in meter space."""
        import math as _math
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        mag = _math.hypot(dx, dy)
        if mag < 1e-9:
            return None
        ux, uy = dx / mag, dy / mag
        px, py = -uy, ux   # perpendicular (left)
        half = width / 2.0
        p1L = (p1[0] + px * half, p1[1] + py * half)
        p2L = (p2[0] + px * half, p2[1] + py * half)
        p1R = (p1[0] - px * half, p1[1] - py * half)
        p2R = (p2[0] - px * half, p2[1] - py * half)
        return _Polygon([p1L, p2L, p2R, p1R])

    # Sort shapes by axis length descending — longest taxi (V)
    # claims its area first, so shorter overlapping taxi candidates
    # get trimmed.
    taxi_shape_indices = sorted(
        (i for i, s in enumerate(shapes)
         if s.kind == "taxi" and s.axis is not None),
        key=lambda i: -shapes[i].axis.length,
    )

    taxi_rect_records = []  # (segments[], width, role)
    taxi_rect_polys = []    # meter-space Polygons used to carve apron
    claimed_union = None
    n_trimmed_too_short = 0
    n_overlapping_skipped = 0

    for i in taxi_shape_indices:
        s = shapes[i]
        role = roles[i]
        trimmed = _trim_axis(s.axis, GAP_M)
        if trimmed is None:
            n_trimmed_too_short += 1
            continue
        if role in (_PS.ROLE_PRIMARY_PARALLEL,
                    _PS.ROLE_SECONDARY_PARALLEL):
            segments = _subdivide_axis(trimmed, SEGMENT_LEN_M)
        else:
            cc = list(trimmed.coords)
            segments = [(cc[j], cc[j + 1])
                        for j in range(len(cc) - 1)]

        # Check overlap with already-claimed taxi area.  If a newly-
        # built rect overlaps more than 20 % with prior claims, skip
        # it entirely — it's a redundant skeleton branch through
        # already-covered pavement.
        candidate_rects = []
        skip_this = False
        for (p1, p2) in segments:
            rect = _rect_between(p1, p2, s.width_m)
            if rect is None or rect.is_empty:
                continue
            if claimed_union is not None and not claimed_union.is_empty:
                try:
                    overlap_area = rect.intersection(claimed_union).area
                except Exception:
                    overlap_area = 0.0
                if rect.area > 0 and overlap_area / rect.area > 0.2:
                    skip_this = True
                    break
            candidate_rects.append(((p1, p2), rect))
        if skip_this or not candidate_rects:
            n_overlapping_skipped += 1
            continue

        # Commit the rects.
        kept_segments = []
        for (seg, rect) in candidate_rects:
            kept_segments.append(seg)
            taxi_rect_polys.append(rect)
        taxi_rect_records.append((kept_segments, s.width_m, role))

        # Update claimed union.
        try:
            new_union = _unary_union([r for _, r in candidate_rects])
            if claimed_union is None:
                claimed_union = new_union
            else:
                claimed_union = _unary_union(
                    [claimed_union, new_union])
        except Exception:
            pass

    n_taxi_rects = sum(len(segs) for segs, _, _ in taxi_rect_records)

    # ── Emit taxi rects ─────────────────────────────────────────
    for (segments, width, role) in taxi_rect_records:
        for (p1, p2) in segments:
            lon1, lat1 = to_ll(p1[0], p1[1])
            lon2, lat2 = to_ll(p2[0], p2[1])
            _emit_sloped_rect(
                lat1, lon1, PREVIEW_ELEV,
                lat2, lon2, PREVIEW_ELEV,
                width)

    # ── Apron carving ────────────────────────────────────────────
    # Apron footprint = full pavement polygon MINUS the union of
    # emitted taxi rects (no additional buffer — taxi rect long
    # sides share a boundary with apron, short ends leave the
    # 10 m gap that's already baked into the trimmed axis).  Every
    # input apron polygon gets trimmed against the taxi rect
    # union so no apron shape overlaps a taxi rect.
    taxi_rect_union = (_unary_union(taxi_rect_polys)
                       if taxi_rect_polys else None)

    n_apron = 0
    for s in shapes:
        if s.kind != "apron":
            continue
        poly = s.polygon
        if poly is None or poly.is_empty:
            continue
        if taxi_rect_union is not None:
            try:
                poly = poly.difference(taxi_rect_union)
            except Exception:
                pass
        if poly is None or poly.is_empty:
            continue
        if poly.geom_type == "Polygon":
            apron_parts = [poly]
        elif hasattr(poly, "geoms"):
            apron_parts = [g for g in poly.geoms
                           if g.geom_type == "Polygon" and not g.is_empty]
        else:
            apron_parts = []
        for p in apron_parts:
            if p.area < 50.0:
                continue   # sliver from carving
            coords_m = list(p.exterior.coords)
            coords_ll = [to_ll(x, y) for (x, y) in coords_m]
            _emit_flat_poly(coords_ll, PREVIEW_ELEV)
            n_apron += 1

    try:
        from O4_UI_Utils import vprint as _vprint
    except Exception:
        _vprint = None
    if _vprint:
        by_role: dict = {}
        for r in roles:
            by_role[r] = by_role.get(r, 0) + 1
        _vprint(1, "[{}] strip-model preview: "
                "{} taxi rects from {} shapes, {} apron polys "
                "({} stubs/crosses too short after 10m trim) — "
                "roles: {}".format(
                    icao, n_taxi_rects,
                    sum(1 for s in shapes if s.kind == "taxi"),
                    n_apron, n_trimmed_too_short,
                    ", ".join(f"{k}={v}" for k, v in sorted(by_role.items()))))


def _emit_pavement_new_model(
    apt_twy_polys_m, apron_polys_m, rwy_union_m, bldg_union_m,
    building_polys_m, bldg_elevations,
    _dem_at, _runway_elev_lookup, _project_point_onto_runway_fn,
    runway_pairs, to_ll, to_m,
    _emit_flat_poly, _emit_sloped_rect, _emit_triangle,
    emitted_flat_shapes, emitted_taxi_rects_m, emitted_twy_quads_m,
    all_emitted_parts_m, icao,
):
    """Minimum-shape pavement decomposition via bottom-up merging.

    Pipeline:
      1. PREPARE — union apt.dat pavement, subtract runways,
         morphological-close + simplify.
      2. SEED — Delaunay of boundary vertices + anchor points only
         (no densification, no interior grid).  DEM elevation at
         each vertex, grade-clamped, rounded to 0.5 m.
      3. MERGE — greedily merge adjacent coplanar triangles.
      4. EMIT — flat N-gon / sloped 4-vertex rect / 3-vertex triangle.
         Shape type decided by elevation pattern:
           - z_range < 0.5 m → flat N-gon
           - slope primarily one direction → sloped 4-vertex rect
           - > 0.5 m change in multiple directions → triangles
      5. VALIDATE — check equal-altitudes-at-joins for all adjacent
         pairs; report violations.
    """
    import math as _math
    from shapely.geometry import Polygon as _Poly, Point as _Pt, MultiPoint
    from shapely import ops as _shp_ops
    from shapely.ops import triangulate as _triangulate
    from shapely.strtree import STRtree

    CLOSE_RADIUS = 25.0
    SIMPLIFY_TOL = 15.0
    MAX_GRADE = 0.015
    MERGE_PLANE_TOL = 0.5    # m — max residual for merge
    MIN_AREA = 200.0
    ANCHOR_SPACING = 30.0
    ELEV_ROUND = 0.5          # m — round all elevations to this step
    JOIN_TOL = 0.5            # m — max allowed mismatch at joins

    def _round_z(z):
        return round(z / ELEV_ROUND) * ELEV_ROUND

    # ── Step 1: PREPARE ──────────────────────────────────────────────
    pavement_sources = [p for p in (apt_twy_polys_m + apron_polys_m)
                        if p is not None and not p.is_empty]
    if not pavement_sources:
        return
    try:
        pav_raw = _shp_ops.unary_union(pavement_sources)
    except Exception:
        return
    if not rwy_union_m.is_empty:
        try:
            pav_raw = pav_raw.difference(rwy_union_m)
        except Exception:
            pass
    try:
        pav_closed = pav_raw.buffer(CLOSE_RADIUS).buffer(-CLOSE_RADIUS)
    except Exception:
        pav_closed = pav_raw
    try:
        pav_simple = pav_closed.simplify(SIMPLIFY_TOL,
                                         preserve_topology=True)
    except Exception:
        pav_simple = pav_closed

    components = []
    if pav_simple.is_empty:
        return
    if hasattr(pav_simple, "geoms"):
        for g in pav_simple.geoms:
            if hasattr(g, "exterior") and g.area >= MIN_AREA:
                components.append(g)
    elif hasattr(pav_simple, "exterior"):
        if pav_simple.area >= MIN_AREA:
            components.append(pav_simple)
    if not components:
        return

    UI.vprint(2, "    {}: new model — {} components, {:.0f} m²".format(
        icao, len(components), sum(c.area for c in components)))

    # ── Helpers ──────────────────────────────────────────────────────
    def _rwy_elev_at(x, y):
        try:
            v = _runway_elev_lookup(x, y)
            if v is not None:
                return float(v)
        except Exception:
            pass
        try:
            lo, la = to_ll(x, y)
            v = _project_point_onto_runway_fn(
                la, lo, runway_pairs, max_dist=30.0)
            if v is not None:
                return float(v)
        except Exception:
            pass
        return None

    def _fit_plane(pts):
        if len(pts) < 3:
            if pts:
                return (0.0, 0.0, sum(p[2] for p in pts) / len(pts))
            return None
        cx = sum(p[0] for p in pts) / len(pts)
        cy = sum(p[1] for p in pts) / len(pts)
        ATA = [[0.0]*3 for _ in range(3)]
        ATz = [0.0]*3
        for (x, y, z) in pts:
            row = [x - cx, y - cy, 1.0]
            for i in range(3):
                for j in range(3):
                    ATA[i][j] += row[i] * row[j]
                ATz[i] += row[i] * z
        def det(m):
            return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
                   -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                   +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
        D = det(ATA)
        if abs(D) < 1e-12:
            return (0.0, 0.0, sum(p[2] for p in pts) / len(pts))
        result = []
        for k in range(3):
            Mk = [list(r) for r in ATA]
            for i in range(3):
                Mk[i][k] = ATz[i]
            result.append(det(Mk) / D)
        a, b, c0 = result
        c = c0 - a * cx - b * cy
        if _math.hypot(a, b) > 0.5:
            return (0.0, 0.0, sum(p[2] for p in pts) / len(pts))
        return (a, b, c)

    def _pz(plane, x, y):
        return plane[0]*x + plane[1]*y + plane[2]

    def _pgrade(plane):
        return _math.hypot(plane[0], plane[1])

    # ── Steps 2–4 per component ──────────────────────────────────────
    total_seed = 0
    total_merged = 0
    total_emitted = 0
    # Track all emitted shapes with their vertex elevations for validation.
    emitted_for_validation = []  # list of (poly_m, [(x,y,z), ...])

    for comp in components:
        # ── 2a. Collect seed points: boundary verts + anchors only ────
        seed_pts = {}  # (x,y) → z
        anchored_keys = set()

        # Boundary vertices (exterior + holes) — already simplified,
        # no densification needed.
        for x, y in comp.exterior.coords[:-1]:
            z = _dem_at(x, y)
            seed_pts[(round(x, 2), round(y, 2))] = _round_z(z if z else 0)
        for hole in comp.interiors:
            for x, y in hole.coords[:-1]:
                z = _dem_at(x, y)
                seed_pts[(round(x, 2), round(y, 2))] = _round_z(z if z else 0)

        # Runway-boundary anchors
        if not rwy_union_m.is_empty:
            try:
                touch = comp.boundary.intersection(rwy_union_m.buffer(1.0))
            except Exception:
                touch = None
            if touch is not None and not touch.is_empty:
                segs = (list(touch.geoms) if hasattr(touch, "geoms")
                        else [touch])
                for seg in segs:
                    if not hasattr(seg, "length") or seg.length < 1:
                        continue
                    n = max(2, int(seg.length / ANCHOR_SPACING) + 1)
                    for k in range(n):
                        t = k / (n - 1) if n > 1 else 0.5
                        try:
                            p = seg.interpolate(t, normalized=True)
                        except Exception:
                            continue
                        rz = _rwy_elev_at(p.x, p.y)
                        if rz is not None:
                            key = (round(p.x, 2), round(p.y, 2))
                            seed_pts[key] = _round_z(rz)
                            anchored_keys.add(key)

        # Building-boundary anchors
        if building_polys_m and bldg_elevations:
            for bi, bp in enumerate(building_polys_m):
                if bp.is_empty or bi not in bldg_elevations:
                    continue
                try:
                    if bp.distance(comp) > 2.0:
                        continue
                    shared = bp.boundary.intersection(comp.buffer(2.0))
                except Exception:
                    continue
                if shared.is_empty:
                    continue
                pad_z = _round_z(float(bldg_elevations[bi]))
                segs = (list(shared.geoms) if hasattr(shared, "geoms")
                        else [shared])
                for seg in segs:
                    if not hasattr(seg, "coords"):
                        continue
                    for cx, cy in seg.coords:
                        key = (round(cx, 2), round(cy, 2))
                        seed_pts[key] = pad_z
                        anchored_keys.add(key)

        pts_list = list(seed_pts.keys())
        if len(pts_list) < 3:
            continue

        # ── 2b. Grade-clamp (FAA smoothing, bounded) ────────────────
        # Build edges from the Delaunay before we use elevations.
        raw_tris = _triangulate(MultiPoint(pts_list))
        inside_tris = [t for t in raw_tris
                       if comp.contains(t.representative_point())]
        if not inside_tris:
            continue

        edges = set()
        for tri in inside_tris:
            vs = list(tri.exterior.coords)[:-1]
            for k in range(3):
                e = tuple(sorted([
                    (round(vs[k][0], 2), round(vs[k][1], 2)),
                    (round(vs[(k+1)%3][0], 2), round(vs[(k+1)%3][1], 2))
                ]))
                edges.add(e)

        MAX_CLAMP_ADJUST = 2.0
        original_z = dict(seed_pts)
        for _iter in range(10):
            moved = False
            for (v1, v2) in edges:
                z1 = seed_pts.get(v1, 0.0)
                z2 = seed_pts.get(v2, 0.0)
                dist = _math.hypot(v1[0]-v2[0], v1[1]-v2[1])
                if dist < 1.0:
                    continue
                max_dz = dist * MAX_GRADE
                if abs(z1 - z2) <= max_dz:
                    continue
                a1 = v1 in anchored_keys
                a2 = v2 in anchored_keys
                if a1 and a2:
                    continue
                o1 = original_z.get(v1, z1)
                o2 = original_z.get(v2, z2)
                if a1:
                    new_z2 = z1 - (1 if z1 > z2 else -1) * max_dz
                    new_z2 = max(o2 - MAX_CLAMP_ADJUST,
                                 min(o2 + MAX_CLAMP_ADJUST, new_z2))
                    if abs(new_z2 - z2) > 0.01:
                        seed_pts[v2] = _round_z(new_z2)
                        moved = True
                elif a2:
                    new_z1 = z2 - (1 if z2 > z1 else -1) * max_dz
                    new_z1 = max(o1 - MAX_CLAMP_ADJUST,
                                 min(o1 + MAX_CLAMP_ADJUST, new_z1))
                    if abs(new_z1 - z1) > 0.01:
                        seed_pts[v1] = _round_z(new_z1)
                        moved = True
                else:
                    excess = abs(z1 - z2) - max_dz
                    adj = excess / 2.0
                    sign = 1 if z1 > z2 else -1
                    n1 = z1 - sign * adj
                    n2 = z2 + sign * adj
                    n1 = max(o1 - MAX_CLAMP_ADJUST,
                             min(o1 + MAX_CLAMP_ADJUST, n1))
                    n2 = max(o2 - MAX_CLAMP_ADJUST,
                             min(o2 + MAX_CLAMP_ADJUST, n2))
                    if abs(n1-z1) > 0.01 or abs(n2-z2) > 0.01:
                        seed_pts[v1] = _round_z(n1)
                        seed_pts[v2] = _round_z(n2)
                        moved = True
            if not moved:
                break

        total_seed += len(inside_tris)

        # ── 2c. Build regions with clamped, rounded elevations ───────
        regions = []
        for tri in inside_tris:
            verts = list(tri.exterior.coords)[:-1]
            pts_3d = []
            for vx, vy in verts:
                key = (round(vx, 2), round(vy, 2))
                z = seed_pts.get(key, 0.0)
                pts_3d.append((vx, vy, float(z)))
            plane = _fit_plane(pts_3d)
            if plane is None:
                continue
            regions.append({'poly': tri, 'pts': pts_3d, 'plane': plane})

        # ── 3. MERGE ─────────────────────────────────────────────────
        while True:
            if len(regions) < 2:
                break
            polys_list = [r['poly'] for r in regions]
            tree = STRtree(polys_list)
            used = [False] * len(regions)
            new_regions = []
            any_merged = False
            for i in range(len(regions)):
                if used[i]:
                    continue
                cur = dict(regions[i])
                cur['pts'] = list(cur['pts'])
                try:
                    cands = tree.query(cur['poly'].buffer(0.5))
                except Exception:
                    cands = []
                for j_idx in cands:
                    j = int(j_idx)
                    if j <= i or used[j]:
                        continue
                    other = regions[j]
                    try:
                        if not cur['poly'].intersects(
                                other['poly'].buffer(0.1)):
                            continue
                    except Exception:
                        continue
                    try:
                        merged_poly = cur['poly'].union(other['poly'])
                        if (not hasattr(merged_poly, "exterior")
                                or merged_poly.is_empty):
                            continue
                    except Exception:
                        continue
                    merged_pts = cur['pts'] + list(other['pts'])
                    merged_plane = _fit_plane(merged_pts)
                    if merged_plane is None:
                        continue
                    if _pgrade(merged_plane) > MAX_GRADE:
                        continue
                    max_res = max(
                        abs(_pz(merged_plane, x, y) - z)
                        for (x, y, z) in merged_pts)
                    if max_res > MERGE_PLANE_TOL:
                        continue
                    cur['poly'] = merged_poly
                    cur['pts'] = merged_pts
                    cur['plane'] = merged_plane
                    used[j] = True
                    any_merged = True
                used[i] = True
                new_regions.append(cur)
            regions = new_regions
            if not any_merged:
                break

        total_merged += len(regions)

        # ── 3.5 PRE-EMIT FIX: insert transition rects at all cliffs ─
        # Scan adjacent region pairs for elevation mismatch.  Where
        # the boundary z differs by > ELEV_ROUND, build a sloped
        # transition rect (length = diff / MAX_GRADE) and subtract
        # it from both regions.  Computed all-at-once before any
        # emission so there's no cascade.
        transition_rects = []  # (rect_poly, z_high, z_low, seg_len)
        for i in range(len(regions)):
            ri = regions[i]
            if ri['poly'].is_empty:
                continue
            for j in range(i + 1, len(regions)):
                rj = regions[j]
                if rj['poly'].is_empty:
                    continue
                try:
                    if not ri['poly'].intersects(rj['poly'].buffer(0.5)):
                        continue
                    shared = ri['poly'].boundary.intersection(
                        rj['poly'].boundary.buffer(1.0))
                except Exception:
                    continue
                if shared.is_empty:
                    continue
                segs = (list(shared.geoms) if hasattr(shared, 'geoms')
                        else [shared])
                best = max(segs,
                           key=lambda s: s.length
                           if hasattr(s, 'length') else 0,
                           default=None)
                if (best is None or not hasattr(best, 'length')
                        or best.length < 5.0):
                    continue
                try:
                    mid = best.interpolate(0.5, normalized=True)
                except Exception:
                    continue
                zi = _pz(ri['plane'], mid.x, mid.y)
                zj = _pz(rj['plane'], mid.x, mid.y)
                diff = abs(zi - zj)
                if diff <= ELEV_ROUND:
                    continue
                # Build transition rect
                rect_len = min(diff / MAX_GRADE, 500.0)
                seg_coords = list(best.coords)
                if len(seg_coords) < 2:
                    continue
                sx = seg_coords[-1][0] - seg_coords[0][0]
                sy = seg_coords[-1][1] - seg_coords[0][1]
                seg_len = _math.hypot(sx, sy)
                if seg_len < 5.0:
                    continue
                ux, uy = sx / seg_len, sy / seg_len
                px, py = -uy, ux
                # Orient toward higher-elev region
                test = _Pt(mid.x + px * 5, mid.y + py * 5)
                higher = ri['poly'] if zi >= zj else rj['poly']
                if not higher.contains(test):
                    px, py = uy, -ux
                z_high = _round_z(max(zi, zj))
                z_low = _round_z(min(zi, zj))
                half = rect_len / 2.0
                c0 = (seg_coords[0][0]+px*half, seg_coords[0][1]+py*half)
                c1 = (seg_coords[-1][0]+px*half, seg_coords[-1][1]+py*half)
                c2 = (seg_coords[-1][0]-px*half, seg_coords[-1][1]-py*half)
                c3 = (seg_coords[0][0]-px*half, seg_coords[0][1]-py*half)
                try:
                    rect_poly = _Poly([c0, c1, c2, c3, c0])
                    if not rect_poly.is_valid or rect_poly.is_empty:
                        continue
                except Exception:
                    continue
                transition_rects.append(
                    (rect_poly, z_high, z_low, seg_len))

        # Subtract all transition rects from merged regions
        if transition_rects:
            all_rects = _shp_ops.unary_union([tr[0] for tr in transition_rects])
            for r in regions:
                try:
                    r['poly'] = r['poly'].difference(all_rects.buffer(0.1))
                except Exception:
                    pass

        # Emit transition rects
        n_trans = 0
        for (rect_poly, z_high, z_low, seg_len) in transition_rects:
            c = list(rect_poly.exterior.coords)
            mh_x = (c[0][0]+c[1][0])/2; mh_y = (c[0][1]+c[1][1])/2
            ml_x = (c[2][0]+c[3][0])/2; ml_y = (c[2][1]+c[3][1])/2
            mh_lo, mh_la = to_ll(mh_x, mh_y)
            ml_lo, ml_la = to_ll(ml_x, ml_y)
            try:
                _emit_sloped_rect(mh_la, mh_lo, z_high,
                                  ml_la, ml_lo, z_low, seg_len)
            except Exception:
                continue
            all_emitted_parts_m.append(rect_poly)
            total_emitted += 1
            n_trans += 1

        if n_trans:
            UI.vprint(2, "    {}: inserted {} transition rects"
                      .format(icao, n_trans))

        # ── 4. EMIT with shape-type decision ─────────────────────────
        for r in regions:
            poly = r['poly']
            plane = r['plane']
            if poly.is_empty or poly.area < MIN_AREA:
                continue

            # Simplify to remove collinear Delaunay artifacts
            try:
                poly = poly.simplify(5.0, preserve_topology=True)
                if poly.is_empty or not hasattr(poly, "exterior"):
                    continue
            except Exception:
                pass

            verts = list(poly.exterior.coords)[:-1]
            n_v = len(verts)
            zs = [_round_z(_pz(plane, x, y)) for x, y in verts]
            z_range = max(zs) - min(zs)

            # Track vertex elevations for validation
            vert_elevs = [(verts[k][0], verts[k][1], zs[k])
                          for k in range(n_v)]

            if z_range < ELEV_ROUND:
                # FLAT — all vertices within one rounding step
                avg_z = _round_z(sum(zs) / len(zs))
                ring_ll = [to_ll(x, y)
                           for x, y in poly.exterior.coords]
                try:
                    _emit_flat_poly(ring_ll, avg_z)
                except Exception:
                    continue
                emitted_flat_shapes.append((poly, avg_z))
                all_emitted_parts_m.append(poly)
                emitted_for_validation.append(
                    (poly, [(v[0], v[1], avg_z) for v in verts]))
                total_emitted += 1

            elif n_v == 3:
                # TRIANGLE — 3 vertices with per-vertex elevations
                lo0, la0 = to_ll(*verts[0])
                lo1, la1 = to_ll(*verts[1])
                lo2, la2 = to_ll(*verts[2])
                _emit_triangle(la0, lo0, zs[0],
                               la1, lo1, zs[1],
                               la2, lo2, zs[2])
                all_emitted_parts_m.append(poly)
                emitted_for_validation.append((poly, vert_elevs))
                total_emitted += 1

            elif n_v == 4:
                # 4-VERTEX — decide: sloped rect vs flat
                # Check if slope is primarily one direction
                zs4 = list(zs)
                best_s = 0
                best_m = -1e18
                for s in range(4):
                    m = (zs4[s] + zs4[(s+1) % 4]) / 2.0
                    if m > best_m:
                        best_m = m
                        best_s = s
                z_hi = _round_z(
                    (zs4[best_s] + zs4[(best_s+1) % 4]) / 2.0)
                z_lo = _round_z(
                    (zs4[(best_s+2) % 4] + zs4[(best_s+3) % 4]) / 2.0)
                if abs(z_hi - z_lo) < ELEV_ROUND:
                    # Actually flat
                    avg_z = _round_z(sum(zs4) / 4.0)
                    ring_ll = [to_ll(x, y)
                               for x, y in poly.exterior.coords]
                    try:
                        _emit_flat_poly(ring_ll, avg_z)
                    except Exception:
                        continue
                    emitted_flat_shapes.append((poly, avg_z))
                    all_emitted_parts_m.append(poly)
                    emitted_for_validation.append(
                        (poly, [(v[0], v[1], avg_z) for v in verts]))
                else:
                    # Sloped rect
                    hi_a = verts[best_s]
                    hi_b = verts[(best_s + 1) % 4]
                    mid_hi = ((hi_a[0]+hi_b[0])/2,
                              (hi_a[1]+hi_b[1])/2)
                    lo_a = verts[(best_s + 2) % 4]
                    lo_b = verts[(best_s + 3) % 4]
                    mid_lo = ((lo_a[0]+lo_b[0])/2,
                              (lo_a[1]+lo_b[1])/2)
                    seg_len = _math.hypot(
                        hi_a[0]-hi_b[0], hi_a[1]-hi_b[1])
                    mh_lo, mh_la = to_ll(*mid_hi)
                    ml_lo, ml_la = to_ll(*mid_lo)
                    try:
                        _emit_sloped_rect(
                            mh_la, mh_lo, z_hi,
                            ml_la, ml_lo, z_lo, seg_len)
                    except Exception:
                        continue
                    all_emitted_parts_m.append(poly)
                    emitted_for_validation.append(
                        (poly, vert_elevs))
                total_emitted += 1

            else:
                # >4 VERTICES, sloped — fan-triangulate from centroid
                # so each triangle preserves the plane's per-vertex
                # elevations.  This avoids the "flat at mean" cliff
                # factory.  For an n-vertex polygon this emits n
                # triangles, each sharing an edge with the next.
                cx = sum(v[0] for v in verts) / n_v
                cy = sum(v[1] for v in verts) / n_v
                cz = _round_z(_pz(plane, cx, cy))
                for k in range(n_v):
                    v0 = verts[k]
                    v1 = verts[(k + 1) % n_v]
                    z0 = zs[k]
                    z1 = zs[(k + 1) % n_v]
                    lo0, la0 = to_ll(*v0)
                    lo1, la1 = to_ll(*v1)
                    loc, lac = to_ll(cx, cy)
                    try:
                        _emit_triangle(la0, lo0, z0,
                                       la1, lo1, z1,
                                       lac, loc, cz)
                    except Exception:
                        continue
                    try:
                        tri_poly = _Poly([v0, v1, (cx, cy)])
                        if tri_poly.is_valid and not tri_poly.is_empty:
                            all_emitted_parts_m.append(tri_poly)
                            emitted_for_validation.append(
                                (tri_poly,
                                 [(v0[0],v0[1],z0),
                                  (v1[0],v1[1],z1),
                                  (cx,cy,cz)]))
                    except Exception:
                        pass
                    total_emitted += 1

    # ── 5. VALIDATE: equal altitudes at joins ────────────────────────
    n_violations = 0
    max_violation = 0.0
    if len(emitted_for_validation) >= 2:
        val_polys = [ev[0] for ev in emitted_for_validation]
        val_tree = STRtree(val_polys)
        for i, (pi, vi) in enumerate(emitted_for_validation):
            try:
                cands = val_tree.query(pi.buffer(1.0))
            except Exception:
                continue
            for j_idx in cands:
                j = int(j_idx)
                if j <= i:
                    continue
                pj, vj = emitted_for_validation[j]
                try:
                    if not pi.intersects(pj.buffer(0.5)):
                        continue
                    shared = pi.boundary.intersection(pj.buffer(1.0))
                except Exception:
                    continue
                if shared.is_empty:
                    continue
                # Sample points along the shared boundary and
                # compare each shape's elevation there.
                test_pts = []
                if hasattr(shared, 'length') and shared.length > 0:
                    for t in (0.0, 0.25, 0.5, 0.75, 1.0):
                        try:
                            p = shared.interpolate(t, normalized=True)
                            test_pts.append((p.x, p.y))
                        except Exception:
                            pass
                elif hasattr(shared, 'geoms'):
                    for g in shared.geoms:
                        if hasattr(g, 'coords'):
                            for c in g.coords:
                                test_pts.append((c[0], c[1]))

                for tx, ty in test_pts:
                    # Elevation from shape i
                    zi = _interp_shape_z(vi, tx, ty)
                    zj = _interp_shape_z(vj, tx, ty)
                    if zi is not None and zj is not None:
                        diff = abs(zi - zj)
                        if diff > JOIN_TOL:
                            n_violations += 1
                            max_violation = max(max_violation, diff)
                            break  # one violation per pair is enough

    UI.vprint(1,
        "    {}: new model — {} seed → {} merged → {} emitted"
        " | join violations: {} (max {:.1f}m)".format(
            icao, total_seed, total_merged, total_emitted,
            n_violations, max_violation))


def _interp_shape_z(vert_elevs, x, y):
    """Interpolate elevation at (x, y) from a shape's vertex elevations.
    Uses inverse-distance weighting from the shape's vertices.
    Returns None if no vertices."""
    if not vert_elevs:
        return None
    import math
    total_w = 0.0
    total_z = 0.0
    for vx, vy, vz in vert_elevs:
        d = math.hypot(x - vx, y - vy)
        if d < 0.1:
            return vz
        w = 1.0 / (d * d)
        total_w += w
        total_z += w * vz
    if total_w < 1e-12:
        return None
    return total_z / total_w


def generate_airport_surface_patches(icao, taxiway_data, building_data,
                                      runway_pairs, tile, dico_apt_entry,
                                      start_node_id=-10000, road_data=None,
                                      xplane_root=None,
                                      runway_segment_chain=None):
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

    def _dem_at(x_m, y_m):
        """DEM sampler closure for the current ``to_ll`` / tile.
        Returns ``None`` on lookup failure so callers can fall back
        to a nearest-anchor or interpolated estimate.
        """
        try:
            lon, lat = to_ll(x_m, y_m)
            return float(tile.dem.alt((lon - tile.lon, lat - tile.lat)))
        except Exception:
            return None

    # Full taxiway width (both sides of centerline)
    TWY_FULL_WIDTH = TAXIWAY_BUFFER_WIDTH * 2.0

    # ── Emitted-geometry accumulator ──────────────────────────────
    # Tracks ALL shapes emitted so far.  Each phase subtracts this
    # before emission to guarantee zero overlap across all phases.
    # Updated incrementally — after each batch of emissions, new
    # shapes are appended and union recomputed.
    all_emitted_parts_m = []  # list of Shapely geoms

    # Taxiway sloped-quad accumulator — hoisted to the top because
    # Phase C0 (apt.dat taxiway rect-chain, runs from A3) writes to
    # it and the coverage-fill pass (Phase C2 tail) reads from it to
    # avoid re-emitting flat patches over already-emitted taxi rects.
    emitted_twy_quads_m = []

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
    # Keep a pristine copy BEFORE the 0.5 m safety inflate below.
    # The unbuffered version is used by Phase A3's runway-taxiway
    # intersection detection: we only want to anchor taxi
    # centerlines at points that ACTUALLY CROSS into the runway
    # rectangle, not at points merely within 0.5 m of the edge
    # (which would wrongly anchor a taxi running parallel to a
    # runway).
    rwy_union_raw_m = rwy_union_m
    # Sub-meter safety inflate: generate_patch_osm emits per-segment
    # rectangles whose corners may drift from the whole-rectangle
    # corners by a fraction of a meter due to float precision in the
    # runway_corners perpendicular offset.  A 0.5 m outward buffer on
    # rwy_union_m (used only for subtraction from aprons, never
    # emitted) absorbs those slivers without visibly affecting any
    # apron edge.
    if not rwy_union_m.is_empty:
        rwy_union_m = rwy_union_m.buffer(0.5)

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
    # When apt.dat is used, we also run each pavement through
    # O4_Pavement_Classifier and split the polygons into two lists:
    # taxiway-classified pieces go into apt_twy_polys_m and apron-
    # classified pieces go into apron_polys_m.  Task 1 in the
    # refinement plan does classification only — both lists are
    # unioned for downstream consumers so the emitted output is
    # bit-identical.  Task 2 consumes apt_twy_polys_m through a
    # dedicated rect-chain taxiway emission path.
    #
    # For OSM-derived data the legacy behaviour is preserved: 5 m
    # outward buffer, exterior-only, everything into apron_polys_m.
    apron_polys_m = []
    apt_twy_polys_m = []   # populated only when apt_dat_used

    # apt.dat pavement simplification tolerance.  Per user direction:
    # "we can simplify taxiway shapes — let's say anything less than
    # a meter — we don't need to follow detailed curves; too many
    # shapes hurts performance."  1 m is well below the resolution
    # of the 1 %-grade apron rule and the 1.5 %-grade taxiway rule
    # (a 1 m horizontal step at 1 % is 1 cm vertical), so losing
    # sub-metre polygon detail costs nothing in the emitted mesh
    # while significantly improving rect-fit of gently curved
    # strips.
    APT_DAT_SIMPLIFY_M = 1.0

    def _poly_to_meter_space(poly_tile_rel):
        """Convert a tile-relative (lon-tile.lon, lat-tile.lat)
        shapely Polygon to meter space, preserving interior rings.
        Applies APT_DAT_SIMPLIFY_M tolerance after conversion.
        Returns None if the converted polygon is degenerate.
        """
        try:
            exterior_m = [
                to_m(x + tile.lon, y + tile.lat)
                for x, y in poly_tile_rel.exterior.coords
            ]
            holes_m = []
            for ring in poly_tile_rel.interiors:
                holes_m.append([
                    to_m(x + tile.lon, y + tile.lat)
                    for x, y in ring.coords
                ])
            p_m = shp_geom.Polygon(
                exterior_m, holes_m if holes_m else None)
            if not p_m.is_valid:
                p_m = p_m.buffer(0)
            if p_m.is_empty or not hasattr(p_m, "exterior"):
                return None
            # Simplify sub-metre detail.  preserve_topology=True keeps
            # interior rings from collapsing and avoids self-intersect.
            try:
                simple = p_m.simplify(
                    APT_DAT_SIMPLIFY_M, preserve_topology=True)
                if (simple.is_valid and not simple.is_empty
                        and hasattr(simple, "exterior")):
                    p_m = simple
            except Exception:
                pass
            return p_m
        except Exception:
            return None

    if apt_dat_used:
        # Per-pavement path: classify each polygon individually so
        # its name (from apt.dat row 110 header) can steer the
        # decision.  The raw per-pavement list lives in
        # dico_apt_entry["_apt_pavements"] (stashed by
        # _dico_from_apt_dat).
        #
        # apt.dat pavements are NOT disjoint in general: custom
        # scenery packs commonly layer a giant "Base Ramp" underlay
        # beneath named sub-ramps, and taxiway mega-polygons that
        # concatenate many real taxiways into one row-110 entry
        # routinely overlap those ramps.  To guarantee both
        # (a) complete pavement coverage and (b) no cross-feature
        # overlap, we run this two-pass pipeline:
        #
        #   Pass 1 — classify + probe: every polygon is classified
        #   apron/taxiway.  Each TAXIWAY-classified polygon is
        #   then probed by O4_Taxiway_Rects.build_taxiway_rects
        #   with the same tunables Phase C0 will use.  If the
        #   probe succeeds (fit ratio ≥ 0.80, strip-like shape),
        #   the polygon is a "rectifiable taxiway" and is held
        #   for C0 emission as sloping rects.  If it fails, it is
        #   demoted to the apron pool — mega multi-taxiway blobs
        #   and chunky non-strip shapes land here.
        #
        #   Pass 2 — union-and-diff: everything in the apron pool
        #   is unioned to eliminate self-overlap ("Base Ramp" ∪
        #   "Main Ramp" etc. collapse into one disjoint geometry).
        #   The successful taxiway regions are subtracted from the
        #   apron union so the Phase C0 rects and the C2
        #   triangulation cover complementary ground with no
        #   overlap.  Fallback taxiways never participate in the
        #   subtract so their full area stays in the apron.
        try:
            import O4_Pavement_Classifier as _PC
            import O4_Taxiway_Rects as _TR
            import O4_Taxiway_Decompose as _TD
            import O4_Taxiway_Skeleton as _TS
        except Exception:
            _PC = None
            _TR = None
            _TD = None
            _TS = None
        apt_pavements = dico_apt_entry.get("_apt_pavements") or []
        classify_counts = {"taxiway": 0, "apron": 0}
        classify_reasons = {}
        rectifiable_twy_m = []    # survived build_taxiway_rects probe
        rectifiable_twy_rects = []  # cached rect chains so C0 does
                                    # not re-run the probe
        pool_apron_m = []         # apron-classified OR fallback twy

        # ────────────────────────────────────────────────────────────
        # Runway elevation lookup closure.
        #
        # User rule: "the goal is a perfectly smooth transition
        # from runway to taxiway".  "Perfectly smooth" requires
        # using the EXACT same elevations as the runway segments
        # emitted by generate_patch_osm — NOT the raw CIFP linear
        # interpolation, which can differ by up to ~1 m on curved
        # runways because generate_patch_osm applies wide-window
        # DEM smoothing and a FAA vertical-curve solver to its
        # profile.
        #
        # We get byte-identical matching by accepting the emitted
        # runway chain (from generate_patch_osm) as a parameter
        # and projecting each taxi endpoint onto the nearest
        # segment of that chain.
        #
        # The chain format: a list of
        #   (lat_a, lon_a, elev_a, lat_b, lon_b, elev_b, width_m)
        # tuples — one per emitted runway rectangle.  We convert
        # each to meter space here and build a spatial index once
        # per airport.
        _runway_segs_m = []   # list of (a_pt, b_pt, elev_a, elev_b,
                              #          length, unit_x, unit_y, width)
        if runway_segment_chain:
            for (la, lo_a, ea, lb, lob, eb, w) in runway_segment_chain:
                try:
                    ax, ay = to_m(lo_a, la)
                    bx, by = to_m(lob, lb)
                except Exception:
                    continue
                dx = bx - ax
                dy = by - ay
                L = sqrt(dx * dx + dy * dy)
                if L < 1e-3:
                    continue
                ux = dx / L
                uy = dy / L
                _runway_segs_m.append((
                    (ax, ay), (bx, by), float(ea), float(eb),
                    L, ux, uy, float(w)))

        def _runway_elev_lookup(x_m, y_m):
            """Return the elevation of the nearest point on the
            emitted runway chain to (x_m, y_m).  None if the chain
            is empty or every runway segment is beyond RADIUS_M.
            The caller is responsible for deciding WHETHER to ask
            (e.g. "is this taxi endpoint near a runway?") — this
            function just interpolates once the decision is made.
            """
            RADIUS_M = 100.0
            if not _runway_segs_m:
                return None
            best_elev = None
            best_dist = 1e18
            for (pa, pb, ea, eb, L, ux, uy, w) in _runway_segs_m:
                # Project onto the runway segment's centerline.
                px = x_m - pa[0]
                py = y_m - pa[1]
                t = px * ux + py * uy
                t_clamped = max(0.0, min(L, t))
                cx = pa[0] + ux * t_clamped
                cy = pa[1] + uy * t_clamped
                d = sqrt((x_m - cx) ** 2 + (y_m - cy) ** 2)
                if d > RADIUS_M:
                    continue
                if d < best_dist:
                    best_dist = d
                    frac = t_clamped / L if L > 0 else 0.0
                    best_elev = ea + frac * (eb - ea)
            return best_elev

        def _try_rectify(strip_poly):
            """Return a rect chain if the strip is rectifiable, else
            None.  Factored so the decomposition branch can re-test
            each branch after morphological splitting.  The MRR-
            aligned builder gets the same runway_polygon +
            runway_elev_lookup as the centerline builder so its
            endpoints are anchored to the runway the same way."""
            if _TR is None:
                return None
            try:
                return _TR.build_taxiway_rects(
                    strip_poly, _dem_at,
                    max_grade=MAX_TAXIWAY_GRADE,
                    runway_polygon=rwy_union_raw_m,
                    runway_elev_lookup=_runway_elev_lookup)
            except Exception:
                return None

        n_twy_decomposed = 0      # mega-polys that needed decomposition
        n_twy_strip_branches = 0  # branches emitted as rects after decomp
        # Global rect-union across ALL taxiway source polygons so a
        # rect chain built from polygon A cannot overlap a rect
        # already emitted from polygon B (apt.dat packs commonly
        # layer overlapping pavements — e.g. "Base Ramp" under
        # "Taxiway Aux B, C, D, E, F, G").
        global_twy_rect_union = shp_geom.Polygon()
        for (poly_tr, pav_name) in apt_pavements:
            if poly_tr is None or poly_tr.is_empty or not poly_tr.is_valid:
                continue
            p_m = _poly_to_meter_space(poly_tr)
            if p_m is None:
                continue
            if _PC is None:
                pool_apron_m.append(p_m)
                continue
            kind = _PC.classify_pavement_m(p_m, pav_name or "")
            classify_counts[kind.kind] = classify_counts.get(
                kind.kind, 0) + 1
            classify_reasons[kind.reason] = classify_reasons.get(
                kind.reason, 0) + 1
            if kind.kind != "taxiway" or _TR is None:
                pool_apron_m.append(p_m)
                continue

            # Fast path: try the polygon whole.
            rects = _try_rectify(p_m)
            if rects:
                rectifiable_twy_m.append(p_m)
                rectifiable_twy_rects.append(rects)
                continue

            # Fast path failed.  Attempt recursive morphological
            # decomposition — a "Taxiway V / U / Q / R / L / M"-style
            # mega-polygon splits into strip branches per taxiway
            # plus junction hubs at the letters' crossings.  Each
            # branch that still fails the rect probe (e.g. an L-bend
            # within one taxiway) is re-decomposed, and so on, until
            # everything is either rectifiable, too small to bother
            # with, or recursion-depth exhausted.
            if _TD is None:
                pool_apron_m.append(p_m)
                continue

            def _try_skeleton(poly):
                """Extract Voronoi centerlines and build a rect
                chain along each one.  Centerlines are processed
                longest-first, and each new chain is clipped
                against the running union of previously-emitted
                chains so sibling branches at a Y-junction don't
                overlap each other — a new rect is DROPPED if its
                area already overlaps an earlier rect by more
                than 5 m².

                Returns (emitted_pairs, residue) where
                emitted_pairs is a list of (centerline, chain)
                and residue is the portion of ``poly`` not
                covered by any chain (for triangulation).
                """
                nonlocal global_twy_rect_union
                if _TS is None or _TR is None:
                    return [], poly
                try:
                    centerlines = _TS.extract_centerlines(poly)
                except Exception:
                    return [], poly
                if not centerlines:
                    return [], poly

                emitted_pairs = []
                # Running union of meter-space rect polygons — seeded
                # with the global union so a rect from this polygon
                # can't collide with an earlier polygon's emission.
                emitted_union = global_twy_rect_union
                OVERLAP_EPS_M2 = 1.0

                for ln in centerlines:
                    try:
                        rc = _TR.build_rects_along_centerline(
                            ln, poly, _dem_at,
                            max_grade=MAX_TAXIWAY_GRADE,
                            seg_length=50.0,
                            # Loose elevation fidelity: user
                            # directive is "25-30 rects for all
                            # taxiways at SPJC", so a single
                            # straight skeleton segment should
                            # usually emit ONE rect even with
                            # typical DEM noise.  3 m tolerance
                            # is still well inside the 1.5 %
                            # grade budget at the 200 m+ segment
                            # lengths we now produce.
                            fidelity_tol=3.0,
                            runway_polygon=rwy_union_raw_m,
                            runway_elev_lookup=_runway_elev_lookup)
                    except Exception:
                        rc = None
                    if not rc:
                        continue

                    # Keep only rects that don't collide with any
                    # earlier emission — both in previous chains
                    # AND earlier rects of THIS chain.  Update the
                    # running union after every accepted rect so
                    # consecutive rects in one chain also see each
                    # other at curved joins.
                    kept_chain = []
                    for r in rc:
                        try:
                            rp = shp_geom.Polygon(r.corners_m)
                            if not rp.is_valid or rp.is_empty:
                                continue
                        except Exception:
                            continue
                        try:
                            ov = (rp.intersection(
                                emitted_union).area
                                if not emitted_union.is_empty else 0.0)
                        except Exception:
                            ov = 0.0
                        if ov > OVERLAP_EPS_M2:
                            continue
                        kept_chain.append(r)
                        try:
                            emitted_union = (
                                emitted_union.union(rp)
                                if not emitted_union.is_empty else rp)
                        except Exception:
                            pass

                    if not kept_chain:
                        continue
                    emitted_pairs.append((ln, kept_chain))

                if not emitted_pairs:
                    return [], poly

                try:
                    residue = poly.difference(emitted_union.buffer(0.5))
                except Exception:
                    residue = poly
                # Publish the updated union to the outer scope so
                # subsequent polygons see this polygon's rects.
                global_twy_rect_union = emitted_union
                return emitted_pairs, residue

            def _recursive_rectify(poly, depth):
                """Walk decomposition recursively, with a final
                fall-through to Voronoi-skeleton centerline
                emission for polygons that don't decompose.
                Appends rectifiable chains to rectifiable_twy_m /
                rectifiable_twy_rects and any non-rectifiable
                residue to pool_apron_m.
                """
                nonlocal n_twy_strip_branches
                if poly is None or poly.is_empty or poly.area < 50.0:
                    if poly is not None and not poly.is_empty:
                        pool_apron_m.append(poly)
                    return

                # 1. Try rectifying the whole polygon as a single
                #    MRR-aligned rect chain (cheap fast path).
                rects = _try_rectify(poly)
                if rects:
                    rectifiable_twy_m.append(poly)
                    rectifiable_twy_rects.append(rects)
                    n_twy_strip_branches += 1
                    return

                # 2. Voronoi-skeleton centerline chains.  Handles
                #    curved or branched polygons — walks the
                #    medial axis and emits rects along each path,
                #    which is what works for the common
                #    "mega-polygon of 6 taxiways concatenated
                #    through a hub" shape.  Tried BEFORE the
                #    morphological decomposition because the
                #    decomposition only succeeds when the hub is
                #    dramatically wider than the strips, whereas
                #    skeleton emission works at all width ratios.
                emitted_pairs, residue = _try_skeleton(poly)
                if emitted_pairs:
                    for strip_ln, chain in emitted_pairs:
                        # Build a tracking polygon for the chain so
                        # the downstream dedup works.
                        chain_shapes = [shp_geom.Polygon(r.corners_m)
                                        for r in chain
                                        if r is not None]
                        if chain_shapes:
                            chain_union = shp_ops.unary_union(chain_shapes)
                            if not chain_union.is_empty:
                                rectifiable_twy_m.append(chain_union)
                                rectifiable_twy_rects.append(chain)
                                n_twy_strip_branches += 1
                    # Residue (bays/slivers not covered by any
                    # centerline chain) goes to the apron path.
                    if residue is not None and not residue.is_empty:
                        if hasattr(residue, "geoms"):
                            for g in residue.geoms:
                                if (hasattr(g, "exterior")
                                        and g.area >= 50.0):
                                    pool_apron_m.append(g)
                        elif hasattr(residue, "exterior"):
                            if residue.area >= 50.0:
                                pool_apron_m.append(residue)
                    return

                # 4. Nothing worked — triangulate the whole
                #    polygon via the apron path.
                pool_apron_m.append(poly)

            n_twy_decomposed += 1
            _recursive_rectify(p_m, depth=4)
        UI.vprint(2,
            "    apt.dat pavement classification: "
            "{} taxiway ({} rectifiable; {} mega-polys decomposed "
            "into {} strip branches), {} apron  (reasons: {})".format(
                classify_counts.get("taxiway", 0),
                len(rectifiable_twy_m),
                n_twy_decomposed,
                n_twy_strip_branches,
                classify_counts.get("apron", 0),
                ", ".join("{}={}".format(k, v)
                          for k, v in sorted(classify_reasons.items()))))

        # Union the rectifiable twy polygons and the apron pool
        # separately.  Subtract the twy union from the apron union
        # so the two sets are disjoint.  Break each union into its
        # component polygons.
        def _union_components(polys):
            if not polys:
                return []
            try:
                u = shp_ops.unary_union(polys)
            except Exception:
                return polys
            if u.is_empty:
                return []
            if hasattr(u, "geoms"):
                return [g for g in u.geoms
                        if hasattr(g, "exterior") and not g.is_empty]
            return [u] if hasattr(u, "exterior") else []

        twy_union_all = (shp_ops.unary_union(rectifiable_twy_m)
                         if rectifiable_twy_m else shp_geom.Polygon())

        apron_union_all = (shp_ops.unary_union(pool_apron_m)
                           if pool_apron_m else shp_geom.Polygon())
        if not twy_union_all.is_empty and not apron_union_all.is_empty:
            try:
                apron_union_all = apron_union_all.difference(
                    twy_union_all)
            except Exception:
                pass
        if apron_union_all.is_empty:
            apron_components = []
        elif hasattr(apron_union_all, "geoms"):
            apron_components = [g for g in apron_union_all.geoms
                                if hasattr(g, "exterior")
                                and not g.is_empty]
        else:
            apron_components = ([apron_union_all]
                                if hasattr(apron_union_all, "exterior")
                                else [])

        # apt_twy_polys_m holds the rectifiable polygons (each
        # associated with a pre-computed rect chain kept in a
        # parallel list of lists).  Phase C0 iterates these in
        # step and just emits the cached rects — no re-probe.
        apt_twy_polys_m = rectifiable_twy_m
        apt_twy_rect_chains = rectifiable_twy_rects
        apron_polys_m = apron_components

        UI.vprint(2,
            "    apt.dat pavement deduplication: "
            "rectifiable taxiway {:.0f} m² in {} components, "
            "apron {:.0f} m² in {} components, "
            "total {:.0f} m² (no self-overlap)".format(
                sum(p.area for p in apt_twy_polys_m),
                len(apt_twy_polys_m),
                sum(p.area for p in apron_polys_m),
                len(apron_polys_m),
                sum(p.area for p in apt_twy_polys_m)
                + sum(p.area for p in apron_polys_m)))
    else:
        apt_twy_rect_chains = []
        # Legacy OSM path: one unified apron MultiPolygon, 5 m buffer.
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
                p_m = _poly_to_meter_space(poly)
                if p_m is None:
                    continue
                apron_polys_m.append(p_m.buffer(APRON_BUFFER))

    # Phase C0: emit classified taxiway polygons as chains of
    # sloping rectangles.
    #
    # Each apt_twy_polys_m entry is a rectifiable strip pavement
    # polygon that already passed the build_taxiway_rects fit-
    # ratio gate during Phase A3's probe pass.  The pre-computed
    # rect chain is in apt_twy_rect_chains at the matching index,
    # so we just emit each segment via _emit_sloped_rect and track
    # the 4 meter-space corners in emitted_taxi_rects_m for
    # downstream subtraction.
    emitted_taxi_rects_m = []
    # Counts feed into the final progress line; initialised here so
    # the summary works whether or not Phase C0 runs.
    n_apt_twy_rects = 0        # sloped rect-chain segments (C0)
    n_apt_twy_flat_rects = 0   # flat rect-chain segments (C0)
    # DEBUG_NEW_MODEL replaces the apt.dat rect-chain emission with a
    # unified pavement decomposition run later; skip the legacy emit
    # loop entirely so it doesn't populate twy_union_m with outdated
    # rects the new model would have to work around.
    if (apt_twy_polys_m and apt_twy_rect_chains
            and not DEBUG_REPLACE_LEGACY_PAVEMENT):
        twy_rect_count = 0
        twy_flat_count = 0
        for twy_poly, rects in zip(apt_twy_polys_m, apt_twy_rect_chains):
            if not rects:
                continue
            for r in rects:
                lon_hi, lat_hi = to_ll(r.center_high[0],
                                       r.center_high[1])
                lon_lo, lat_lo = to_ll(r.center_low[0],
                                       r.center_low[1])
                _emit_sloped_rect(
                    lat_hi, lon_hi, r.elev_high,
                    lat_lo, lon_lo, r.elev_low,
                    r.width_m)
                twy_rect_count += 1
                if r.is_flat:
                    twy_flat_count += 1
                try:
                    rp = shp_geom.Polygon(r.corners_m)
                    if rp.is_valid and not rp.is_empty:
                        emitted_taxi_rects_m.append(rp)
                except Exception:
                    pass
        twy_fallback_count = 0  # dedup pass already reclassified them
        n_apt_twy_rects = twy_rect_count - twy_flat_count
        n_apt_twy_flat_rects = twy_flat_count
        if emitted_taxi_rects_m:
            all_emitted_parts_m.extend(emitted_taxi_rects_m)
            # Feed the rects into emitted_twy_quads_m so the later
            # coverage-fill pass sees them as "already covered" and
            # skips re-emitting flat patches for the same area.
            emitted_twy_quads_m.extend(emitted_taxi_rects_m)
            # Union the emitted taxi rects into twy_union_m so the
            # existing C2 apron subtraction (`clipped = clipped.
            # difference(twy_union_m)` at line ~4005) carves them
            # out cleanly.  Without this the apron triangulation
            # would produce triangles overlapping the rect chain.
            # Buffer by 0.5 m to swallow sub-meter slivers from
            # corner rounding in runway_corners precision.
            try:
                taxi_union = shp_ops.unary_union(emitted_taxi_rects_m)
                if not taxi_union.is_empty:
                    if twy_union_m.is_empty:
                        twy_union_m = taxi_union.buffer(0.5)
                    else:
                        twy_union_m = shp_ops.unary_union(
                            [twy_union_m, taxi_union.buffer(0.5)])
            except Exception:
                pass
        UI.vprint(2,
            "    apt.dat taxiway rect-chain: {} rects ({} flat) "
            "over {} rectifiable polygons".format(
                twy_rect_count, twy_flat_count, len(apt_twy_polys_m)))

    try:
        apron_union_m = shp_ops.unary_union(apron_polys_m) if apron_polys_m \
            else shp_geom.Polygon()
    except Exception:
        apron_union_m = shp_geom.Polygon()

    # A4: Building pads
    # Commit 10: BLDG_PAD is 0 — use the EXACT footprint.  The
    # legacy's 5 m outward buffer was a workaround for OSM tracing
    # imprecision and to glue adjacent buildings together, but it
    # directly caused both overlap classes:
    #   * flat∩flat: adjacent building footprints padded by 5 m each
    #     developed 10 m of overlap space.
    #   * flat∩triangle: the apron hole was shaped to the padded
    #     union, so apron triangles could extend up to 5 m into the
    #     actual wall line.
    # With precise cut-outs, the apron subtraction exactly equals
    # the emitted pad and overlap drops to zero by construction.
    # Adjacent-building merging is now handled by the morphological
    # closing in Phase A4b, which bridges gaps up to BLDG_MERGE_DIST
    # wide without inflating the outer boundary.
    BLDG_PAD = 0.0
    BLDG_SIMPLIFY = 2.0   # was 8 m — tighter now that we're not
                          # hiding a 5 m pad.  2 m still kills OSM
                          # cm-scale tracing jitter.
    building_polys_m = []

    for bldg in (building_data or []):
        fp = bldg["footprint"]
        if len(fp) < 4:
            continue
        coords_m = [to_m(lon, lat) for lon, lat in fp]
        try:
            p = shp_geom.Polygon(coords_m)
            if not p.is_valid:
                p = p.buffer(0)
            if p.is_empty or p.area <= 20.0:
                continue
            # Simplify only (no outward buffer).
            simple = p.simplify(BLDG_SIMPLIFY, preserve_topology=True)
            if simple.is_valid and not simple.is_empty:
                building_polys_m.append(simple)
            else:
                building_polys_m.append(p)
        except Exception:
            pass

    # Source OSM data can include buildings whose footprints overlap
    # each other (common when a terminal and a jetway are both tagged
    # building=* in OSM and the jetway's outline clips into the
    # terminal).  Without padding to hide it, this produces directly
    # overlapping building pads.  Resolve by unioning overlapping
    # pairs: sort largest-first, and subtract each already-processed
    # building from the next.  Result: every pad is interior-disjoint
    # from every other pad.
    if len(building_polys_m) > 1:
        building_polys_m.sort(key=lambda p: p.area, reverse=True)
        cleaned = []
        for bp in building_polys_m:
            cur = bp
            for prev in cleaned:
                if prev.is_empty or cur.is_empty:
                    continue
                try:
                    if prev.intersects(cur):
                        cur = cur.difference(prev)
                        if not cur.is_valid:
                            cur = cur.buffer(0)
                except Exception:
                    pass
                if cur.is_empty:
                    break
            if cur.is_empty or cur.area < 20.0:
                continue
            # Difference can produce a MultiPolygon — keep the
            # largest piece (a building split in two by a larger
            # neighbour is unusual; the smaller fragment is
            # typically noise).
            if hasattr(cur, "geoms"):
                cur = max(cur.geoms, key=lambda g: g.area)
            if (not cur.is_empty
                    and hasattr(cur, "exterior")
                    and cur.area >= 20.0):
                cleaned.append(cur)
        building_polys_m = cleaned

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

    # ════════════════════════════════════════════════════════════════
    # NEW MODEL: pavement decomposition
    # ════════════════════════════════════════════════════════════════
    # When O4_NEW_MODEL=1, replace the legacy Phase C0 (apt.dat rect-
    # chain) + C1 (OSM taxi buffer) + C2 (apron triangulation) + D
    # (junction triangulation) with one unified pass that emits:
    #
    #   - One flat N-vertex polygon per connected pavement component
    #     at a DEM-driven target elevation (the exterior ring only —
    #     interior holes for buildings are filled and then re-carved
    #     by Phase C3 emitting building pads at their own elevations;
    #     Phase C4 then detects any flat-to-pad elevation mismatch
    #     and auto-creates transition strips).
    #   - 4-vertex sloped rectangles at runway-touch boundaries,
    #     using _runway_elev_lookup (the emitted chain, not raw CIFP
    #     interp) for the runway-side elevation.  Strip depth sized
    #     by grade budget up to 1.5% per user allowance.
    #
    # Runs BEFORE Phase C3 so building pads can carve into the flat
    # fills cleanly.
    if DEBUG_STRIP_MODEL:
        _emit_pavement_strip_model(
            apt_twy_polys_m, apron_polys_m, rwy_union_m,
            to_ll, to_m, apt_data if apt_dat_used else None,
            _emit_flat_poly, _emit_sloped_rect, _emit_triangle,
            icao)
    elif DEBUG_NEW_MODEL:
        _emit_pavement_new_model(
            apt_twy_polys_m, apron_polys_m, rwy_union_m, bldg_union_m,
            building_polys_m, bldg_elevations,
            _dem_at, _runway_elev_lookup, _project_point_onto_runway,
            runway_pairs, to_ll, to_m,
            _emit_flat_poly, _emit_sloped_rect, _emit_triangle,
            emitted_flat_shapes, emitted_taxi_rects_m, emitted_twy_quads_m,
            all_emitted_parts_m, icao)
        # Update twy_union_m to include the new strips so downstream
        # Phase C3 building-pad clipping avoids them.
        try:
            if emitted_twy_quads_m:
                new_twy_union = shp_ops.unary_union(emitted_twy_quads_m)
                if twy_union_m.is_empty:
                    twy_union_m = new_twy_union.buffer(0.5)
                else:
                    twy_union_m = shp_ops.unary_union(
                        [twy_union_m, new_twy_union.buffer(0.5)])
        except Exception:
            pass

    # C1: Taxiway shapes — clip away runways + junction zones
    n_twy_flat = 0
    n_twy_sloped = 0
    sloped_twy_parts_m = []   # curve areas sent to Phase D for triangulation
    emitted_flat_twy_m = []   # meter-space footprints of emitted flat twy polys

    # Merge parallel taxiways within 15m into single surfaces
    PARALLEL_TWY_MERGE_DIST = 15.0
    merged_twy_groups = {}  # ti -> group_id
    twy_group_id = 0
    _twy_buffer_iter = ([] if DEBUG_REPLACE_LEGACY_PAVEMENT
                        else list(range(len(twy_buffers_m))))
    for i in _twy_buffer_iter:
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
            # NOTE: commit 9 added a TERMINAL_FLAT_BUFFER_M flat ring
            # around each building at the pad elevation, to keep the
            # ground flat for ~10 m around each terminal.  At airports
            # with densely-packed buildings at varying pad elevations
            # (SPJC has 313 buildings spanning ~40 m of pad elevation
            # across ~2 km of pavement) the overlapping flat rings
            # create incompatible anchor constraints — two adjacent
            # buildings' flat rings demand two different elevations
            # at the same apron point, and the grade clamp cannot
            # resolve it.  Per the user's rule "aprons share an edge
            # with the building pad, THEN slope gently to taxiways",
            # the flat ring is unnecessary: vertex coincidence at the
            # actual building edge (handled by the shared-boundary
            # block above) is enough, and the apron can slope freely
            # from there at ≤ MAX_APRON_GRADE.  Ring code removed.
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
    emitted_apron_triangles_m = []   # actual Phase C2 triangles (for
                                     # Phase F drainage detection to see
                                     # via the emit accumulator)

    # DEBUG_TAXI_ONLY: skip apron triangulation entirely.  Apron polygons
    # still feed taxiway classification upstream, but no Phase C2 apron
    # triangles are emitted so the patch contains only taxi rects.
    # DEBUG_NEW_MODEL: the new unified decomposition already handles
    # apron pavement as part of the connected-component flat fills;
    # skip the triangulation path so we don't double-emit.
    _apron_polys_iter = ([] if (DEBUG_TAXI_ONLY
                                or DEBUG_REPLACE_LEGACY_PAVEMENT)
                         else apron_polys_m)
    for ai, p_m in enumerate(_apron_polys_iter):
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
            # Taxiway-classified pavement is already handled by the
            # Phase C0 rect-chain emission — anything reaching this
            # apron triangulation loop is apron-classified and must
            # respect the stricter 1.0 % apron grade rule.
            pavement_grade = MAX_APRON_GRADE
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

            # Emit each triangle from the adaptive mesh.  Skip any
            # triangle whose geometry actually overlaps a building
            # pad (real intersection area, not just centroid — the
            # user rightly pointed out that edge-crossing triangles
            # can straddle a building while their centroid is
            # outside it).  A triangle that merely shares an edge
            # with a building boundary has zero intersection area
            # and is kept; a triangle that reaches any interior
            # slice of a building pad is dropped.
            #
            # Each surviving triangle is also added to
            # emitted_apron_triangles_m so Phase F's drainage
            # detection sees the precise triangulated shape via
            # _emitted_union() instead of just the apron-piece
            # outline (which has building-shaped holes that
            # drainage could land inside).
            avg_e = (sum(v[2] for t in tris for v in t)
                     / (3.0 * len(tris)))
            OVERLAP_EPS = 0.5   # m² — ignore vanishing boundary slop
            # Sliver filters: Delaunay on complex polygons routinely
            # produces a handful of near-colinear "sliver" triangles
            # whose shortest edge is 1-3 m.  adaptive_triangulate's
            # fit_plane returns None for colinear triples so the
            # grade clamp silently skips them, and the sliver is
            # emitted with raw DEM/anchor elevations — which then
            # look like 50-100 % grade in the final patch because
            # the sliver's short edge amplifies any ΔZ into an
            # absurd edge gradient.  Drop sliver triangles here and
            # rely on the per-piece coverage-fill below to paint
            # their area as flat residue at the piece's average
            # elevation.
            SLIVER_MIN_AREA_M2 = 5.0
            SLIVER_MIN_EDGE_M = 2.0
            # Plane gradient filter: drop any triangle whose plane
            # slope exceeds 2× MAX_APRON_GRADE.  These are slivers
            # whose plane fit is dominated by a short edge with a
            # large Δz — the grade clamp cannot fix them because
            # pulling corners toward the centroid just produces a
            # different degenerate plane.  Dropped triangles become
            # holes in the triangulation that the per-piece coverage
            # fill below paints as flat residue at the piece's
            # average elevation, which preserves coverage and is a
            # much better approximation than a sliver with a 300 %
            # plane gradient.
            MAX_PLANE_GRADIENT = 2.0 * MAX_APRON_GRADE
            for (v0, v1, v2) in tris:
                tri_poly = None
                try:
                    tri_poly = shp_geom.Polygon([
                        (v0[0], v0[1]),
                        (v1[0], v1[1]),
                        (v2[0], v2[1])])
                    if (not tri_poly.is_valid
                            or tri_poly.is_empty):
                        continue
                    if tri_poly.area < SLIVER_MIN_AREA_M2:
                        continue
                    # Min edge length check — catches long-thin
                    # slivers whose area is just above threshold.
                    e01 = sqrt((v0[0] - v1[0]) ** 2
                               + (v0[1] - v1[1]) ** 2)
                    e12 = sqrt((v1[0] - v2[0]) ** 2
                               + (v1[1] - v2[1]) ** 2)
                    e20 = sqrt((v2[0] - v0[0]) ** 2
                               + (v2[1] - v0[1]) ** 2)
                    if min(e01, e12, e20) < SLIVER_MIN_EDGE_M:
                        continue
                    # Plane gradient check via inline fit_plane.
                    v1x = v1[0] - v0[0]; v1y = v1[1] - v0[1]
                    v1z = v1[2] - v0[2]
                    v2x = v2[0] - v0[0]; v2y = v2[1] - v0[1]
                    v2z = v2[2] - v0[2]
                    nz = v1x * v2y - v1y * v2x
                    if abs(nz) > 1e-9:
                        nx = v1y * v2z - v1z * v2y
                        ny = v1z * v2x - v1x * v2z
                        grad = sqrt((nx * nx + ny * ny)) / abs(nz)
                        if grad > MAX_PLANE_GRADIENT:
                            continue
                    if (not bldg_union_m.is_empty
                            and tri_poly.intersection(
                                bldg_union_m).area > OVERLAP_EPS):
                        continue
                except Exception:
                    continue
                lon0, lat0 = to_ll(v0[0], v0[1])
                lon1, lat1 = to_ll(v1[0], v1[1])
                lon2, lat2 = to_ll(v2[0], v2[1])
                _emit_triangle(
                    lat0, lon0, round(v0[2], 1),
                    lat1, lon1, round(v1[2], 1),
                    lat2, lon2, round(v2[2], 1))
                n_apron_tris += 1
                emitted_apron_triangles_m.append(tri_poly)

            # (The per-piece flat-residue coverage fill that lived
            # here was a bandaid: it emitted flat polygons at the
            # piece MEAN elevation over concave bays that Delaunay
            # missed, producing visible multi-metre elevation steps
            # at the flat-to-triangle edges.  The correct fix is
            # upstream — decompose multi-taxiway mega-polygons so
            # they go through the rect-chain path, and let the
            # global coverage-fill pass at the end of C2 paint any
            # genuinely missing slivers at a locally-sampled
            # elevation.  See commit log for details.)

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

    # DEBUG_TAXI_ONLY: skip building pad emission.
    _bldg_iter = ([] if DEBUG_TAXI_ONLY
                  else list(enumerate(building_polys_m)))
    for bi, bp in _bldg_iter:
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
    # DEBUG_TAXI_ONLY / DEBUG_NEW_MODEL: skip coverage-fill sweep.
    try:
        if DEBUG_TAXI_ONLY or DEBUG_REPLACE_LEGACY_PAVEMENT:
            raise RuntimeError("skip coverage fill")
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
            # Compute what's already covered (or will be by Phase D).
            # IMPORTANT: we track the ACTUAL emitted apron triangles
            # (emitted_apron_triangles_m), not the apron piece
            # outlines in emitted_flat_shapes.  adaptive_triangulate
            # can leave tiny concave-bay gaps that the piece outline
            # papers over; using the triangle geometry directly lets
            # the coverage-fill pass detect and patch those gaps.
            # terminal_zone_m is a classification mask (used by the
            # legacy taxiway path to identify "terminal-adjacent"
            # aprons); it is NOT emitted as a shape and must NOT be
            # treated as covered, or the coverage-fill pass blindly
            # skips 1 M+ m² of pavement.
            covered_parts = [rwy_union_m]  # runways always covered
            if not junction_zone_m.is_empty:
                covered_parts.append(junction_zone_m)
            # Non-apron emitted flats (flat taxiway polys, transition
            # strip substitutes, etc.) still track via piece outlines.
            # Apron pieces are skipped here and accounted for via
            # emitted_apron_triangles_m below so concave-bay gaps
            # become candidate coverage fills.
            apron_piece_set = {id(p) for p in complex_apron_parts_m}
            for s_poly, _ in emitted_flat_shapes:
                if s_poly.is_empty:
                    continue
                if id(s_poly) in apron_piece_set:
                    continue
                covered_parts.append(s_poly)
            if emitted_apron_triangles_m:
                covered_parts.extend(emitted_apron_triangles_m)
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

    # DEBUG_TAXI_ONLY / DEBUG_NEW_MODEL: skip transition strips.
    if (not DEBUG_TAXI_ONLY and not DEBUG_REPLACE_LEGACY_PAVEMENT
            and len(emitted_flat_shapes) >= 2):
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
                    # Subtract existing emitted flat shapes that
                    # actually touch the strip's bbox.  The previous
                    # version scanned ALL emitted_flat_shapes
                    # (~800 polys at SPJC) for every pair, producing
                    # 260 000 difference() calls — 4.4 s of the
                    # pre-optimisation 45 s total.  The STRtree
                    # query cuts that to the 0-20 neighbours that
                    # actually overlap the strip.
                    try:
                        nbrs = tree.query(strip)
                    except Exception:
                        nbrs = []
                    for n_idx in nbrs:
                        n = int(n_idx)
                        if n == i or n == j:
                            continue
                        s_poly = fs_polys[n]
                        if s_poly.is_empty:
                            continue
                        if not s_poly.intersects(strip):
                            continue
                        strip = strip.difference(s_poly)
                        if strip.is_empty:
                            break
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

    # DEBUG_TAXI_ONLY: skip junction triangulation entirely.
    if not DEBUG_TAXI_ONLY and not junction_zone_m.is_empty:
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

    if DEBUG_TAXI_ONLY or DEBUG_REPLACE_LEGACY_PAVEMENT:
        triangle_zone = shp_geom.Polygon()
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
    for tri_poly in emitted_apron_triangles_m:
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
        if DEBUG_TAXI_ONLY or DEBUG_REPLACE_LEGACY_PAVEMENT:
            raise RuntimeError("skip Phase E")
        if not airport_footprint.is_empty and airport_footprint.area > 100:
            emitted_union = _emitted_union()

            # ── Portal exclusion zones ────────────────────────────────
            # Pre-scan tunnel roads to find portal points so Phase E1
            # can carve a gap in the boundary band around each portal,
            # leaving room for the Phase E2 grade-down/retaining-wall
            # rects that start ON the boundary and dip into the band
            # for a few metres.
            PORTAL_EXCLUSION_RADIUS = 35.0
            portal_exclusion = shp_geom.Polygon()
            if road_lines_m:
                boundary_line_pre = airport_footprint.boundary
                if boundary_line_pre.geom_type == "MultiLineString":
                    boundary_line_pre = max(
                        boundary_line_pre.geoms,
                        key=lambda g: g.length)
                portal_pts = []
                for r_ls, _h, r_tun, _b in road_lines_m:
                    if not r_tun:
                        continue
                    if boundary_line_pre.distance(r_ls) > 100.0:
                        continue
                    try:
                        xing = r_ls.intersection(boundary_line_pre)
                    except Exception:
                        xing = None
                    if xing is not None and not xing.is_empty:
                        if hasattr(xing, "geoms"):
                            for g in xing.geoms:
                                if hasattr(g, "x"):
                                    portal_pts.append((g.x, g.y))
                        elif hasattr(xing, "x"):
                            portal_pts.append((xing.x, xing.y))
                    r_coords = list(r_ls.coords)
                    for rc in (r_coords[0], r_coords[-1]):
                        pt = shp_geom.Point(rc)
                        if boundary_line_pre.distance(pt) < 80.0:
                            portal_pts.append(rc)
                if portal_pts:
                    circles = [
                        shp_geom.Point(pp[0], pp[1]).buffer(
                            PORTAL_EXCLUSION_RADIUS)
                        for pp in portal_pts]
                    portal_exclusion = shp_ops.unary_union(circles)

            # ── E1: Sloped boundary band (segmented) ───────────────
            # Walk the airport boundary in BAND_SEG_LENGTH steps and
            # emit one rect per segment.  The rect is band_width deep
            # (inward from the boundary) and spans seg_length along
            # the boundary.  Elevations at the two short edges come
            # from the CIFP surface model (which projects the runway
            # slope outward via IDW over all runway/building/taxiway
            # anchors); the rect is emitted sloped when the
            # endpoint elevations differ by ≥ 0.1 m, otherwise flat.
            #
            # This replaces the old "one flat polygon per connected
            # band piece" approach, which painted huge areas at a
            # single elevation and didn't track the airport slope.
            # Each emitted rect is clipped against emitted_union +
            # portal_exclusion; if the clipped result has ≥ 95 % of
            # the original area we emit the full sloped rect, else
            # we emit just the clipped geometry as a flat polygon at
            # the averaged elevation (handles partial overlaps with
            # buildings, runway overruns and tunnel portal exclusion
            # zones without leaving gaps).
            try:
                if emitted_union.is_empty:
                    subtract_base = portal_exclusion
                elif portal_exclusion.is_empty:
                    subtract_base = emitted_union.buffer(1.0)
                else:
                    subtract_base = shp_ops.unary_union([
                        emitted_union.buffer(1.0),
                        portal_exclusion])
                # Rects already emitted in this loop — we subtract
                # them from each new rect so adjacent segments at
                # polygon corners don't overlap each other.
                e1_emitted_union = shp_geom.Polygon()

                boundary_line_e1 = airport_footprint.boundary
                e1_lines = (
                    list(boundary_line_e1.geoms)
                    if boundary_line_e1.geom_type == "MultiLineString"
                    else [boundary_line_e1])

                for e1_line in e1_lines:
                    if e1_line.length < BAND_SEG_LENGTH:
                        continue
                    n_segs = max(1, int(
                        e1_line.length / BAND_SEG_LENGTH))
                    step = e1_line.length / n_segs
                    for i in range(n_segs):
                        s_pt = e1_line.interpolate(i * step)
                        e_pt = e1_line.interpolate((i + 1) * step)
                        sx_b, sy_b = s_pt.x, s_pt.y
                        ex_b, ey_b = e_pt.x, e_pt.y
                        dxs = ex_b - sx_b
                        dys = ey_b - sy_b
                        slen = sqrt(dxs * dxs + dys * dys)
                        if slen < 1.0:
                            continue

                        # Inward perpendicular: left-90 of along,
                        # probed against the footprint in case the
                        # ring is CW-oriented.
                        along_b = (dxs / slen, dys / slen)
                        perp_b = (-along_b[1], along_b[0])
                        mid_x = (sx_b + ex_b) / 2.0
                        mid_y = (sy_b + ey_b) / 2.0
                        probe_pt = shp_geom.Point(
                            mid_x + perp_b[0] * 2.0,
                            mid_y + perp_b[1] * 2.0)
                        if not airport_footprint.contains(probe_pt):
                            perp_b = (-perp_b[0], -perp_b[1])

                        # Centerline offset inward by band_width/2
                        # so the runway_corners-built rect spans
                        # [boundary, boundary + band_width inward].
                        half_bw = BOUNDARY_BAND_WIDTH / 2.0
                        cs_x = sx_b + perp_b[0] * half_bw
                        cs_y = sy_b + perp_b[1] * half_bw
                        ce_x = ex_b + perp_b[0] * half_bw
                        ce_y = ey_b + perp_b[1] * half_bw

                        # Sample elevation at each short-edge center
                        # (= boundary point at start / end of the
                        # segment).  CIFP surface model first; DEM
                        # fallback if out of range.
                        e_s = _cifp_surface_elevation(
                            sx_b, sy_b, max_search=500.0)
                        if e_s is None:
                            slon_b, slat_b = to_ll(sx_b, sy_b)
                            try:
                                e_s = tile.dem.alt(
                                    (slon_b - tile.lon,
                                     slat_b - tile.lat))
                            except Exception:
                                continue
                        e_e = _cifp_surface_elevation(
                            ex_b, ey_b, max_search=500.0)
                        if e_e is None:
                            elon_b, elat_b = to_ll(ex_b, ey_b)
                            try:
                                e_e = tile.dem.alt(
                                    (elon_b - tile.lon,
                                     elat_b - tile.lat))
                            except Exception:
                                continue
                        if e_s is None or e_e is None:
                            continue

                        # Build the rect polygon via runway_corners
                        # for the clipping step.
                        cs_lon, cs_lat = to_ll(cs_x, cs_y)
                        ce_lon, ce_lat = to_ll(ce_x, ce_y)
                        try:
                            rc_corners = runway_corners(
                                cs_lat, cs_lon, ce_lat, ce_lon,
                                BOUNDARY_BAND_WIDTH)
                            if rc_corners is None:
                                continue
                            rp = shp_geom.Polygon([
                                to_m(c[1], c[0]) for c in rc_corners])
                            if not rp.is_valid:
                                rp = rp.buffer(0)
                            if rp.is_empty or rp.area < 5.0:
                                continue
                        except Exception:
                            continue

                        full_area = rp.area
                        # Clip against prior phase emitted shapes +
                        # portal exclusion (static `subtract_base`)
                        # followed by rects already emitted earlier
                        # in THIS walk (`e1_emitted_union`).  The
                        # previous version unioned both into one
                        # polygon per iteration, which made a single
                        # line of the Phase E1 loop account for 12 s
                        # of the 32 s runtime at SPJC because the
                        # static half was re-unioned 400 times.  Two
                        # sequential differences are equivalent to
                        # one difference-of-union, and the static
                        # polygon is pre-computed once above.
                        try:
                            rp_clean = rp
                            if not subtract_base.is_empty:
                                rp_clean = rp_clean.difference(subtract_base)
                            if (not rp_clean.is_empty
                                    and not e1_emitted_union.is_empty):
                                rp_clean = rp_clean.difference(
                                    e1_emitted_union)
                        except Exception:
                            rp_clean = rp
                        if (rp_clean.is_empty
                                or rp_clean.area < 10.0):
                            continue

                        clean_ratio = rp_clean.area / full_area
                        # Full-sloped-rect branch only when the
                        # overlap is BOTH <0.1% of rect area AND
                        # <0.25 m² absolute — otherwise a 0.5 m²
                        # sliver survives the ratio test on a 601
                        # m² rect and shows up in the overlap
                        # audit.
                        if (clean_ratio >= 0.999
                                and full_area - rp_clean.area < 0.25):
                            # Emit full sloped rect.
                            _emit_sloped_rect(
                                cs_lat, cs_lon, round(e_s, 1),
                                ce_lat, ce_lon, round(e_e, 1),
                                BOUNDARY_BAND_WIDTH)
                            blend_polys_m.append(rp)
                            n_blend += 1
                            try:
                                e1_emitted_union = shp_ops.unary_union(
                                    [e1_emitted_union, rp])
                            except Exception:
                                pass
                        else:
                            # Partial overlap — emit the clipped
                            # geometry as a flat polygon at the
                            # averaged elevation.  This is one of
                            # the "small flat polygons for flat
                            # areas" connecting adjacent sloped
                            # rects where the band is cut short by
                            # a building edge, runway apron, or a
                            # previously-emitted E1 rect at a
                            # polygon corner.
                            avg_e = round((e_s + e_e) / 2.0, 1)
                            pieces = (
                                [rp_clean]
                                if rp_clean.geom_type == "Polygon"
                                else list(rp_clean.geoms))
                            for piece in pieces:
                                if (piece.is_empty
                                        or piece.area < 10.0
                                        or not hasattr(
                                            piece, "exterior")):
                                    continue
                                ring_ll = [
                                    to_ll(mx, my)
                                    for mx, my in
                                    piece.exterior.coords]
                                _emit_flat_poly(ring_ll, avg_e)
                                blend_polys_m.append(piece)
                                n_blend += 1
                                try:
                                    e1_emitted_union = (
                                        shp_ops.unary_union(
                                            [e1_emitted_union, piece]))
                                except Exception:
                                    pass
            except Exception:
                pass

            # ── E2: Tunnel portals ─────────────────────────────────
            # Model: the tunnel way (OSM tunnel=yes) has a start node
            # OUTSIDE the airport where the road is at surface
            # elevation, and it crosses the airport boundary at the
            # portal, where the road is already at tunnel depth (~6 m
            # below airport grade).  The ramp is the outside-airport
            # sub-segment of the tunnel way.  We emit:
            #
            #   1. A single sloped rect covering the whole ramp,
            #      HIGH (surface DEM) at the OSM tunnel node, LOW
            #      (apt_elev − DEPTH) at the airport boundary.
            #   2. Retaining walls flat at apt_elev along both sides
            #      of the ramp plus an end-cap wrapping around the
            #      LOW (portal) end, with a small gap between wall
            #      and ramp so they are not touching.
            #
            # Divided-highway tunnels are two OSM ways — we cluster
            # ramps by portal location and emit one wider combined
            # rect covering all carriageways per cluster.
            TUNNEL_WALL_GAP = 0.5       # m: gap between wall and ramp
            CARRIAGEWAY_BASE_WIDTH = TUNNEL_RECT_WIDTH  # per carriageway
            if road_lines_m:
                boundary_line = airport_footprint.boundary
                if boundary_line.geom_type == "MultiLineString":
                    boundary_line = max(
                        boundary_line.geoms, key=lambda g: g.length)

                # ── Collect ramps: outside sub-LineStrings of each
                # tunnel way.  Each ramp's portal end = the endpoint
                # on/near the boundary; outside end = the OSM tunnel
                # node (far from boundary).
                ramps = []  # (portal_xy, outside_xy, length, r_hwy)
                for r_ls, r_hwy, r_tunnel, _ in road_lines_m:
                    if not r_tunnel:
                        continue
                    if boundary_line.distance(r_ls) > 100.0:
                        continue
                    try:
                        outside = r_ls.difference(airport_footprint)
                    except Exception:
                        outside = None
                    if outside is None or outside.is_empty:
                        continue
                    pieces = (list(outside.geoms)
                              if hasattr(outside, "geoms")
                              else [outside])
                    for piece in pieces:
                        if (piece.is_empty
                                or piece.geom_type != "LineString"
                                or piece.length < 5.0):
                            continue
                        coords = list(piece.coords)
                        p0 = coords[0]
                        p1 = coords[-1]
                        d0 = boundary_line.distance(shp_geom.Point(p0))
                        d1 = boundary_line.distance(shp_geom.Point(p1))
                        # Must touch the boundary at one end
                        if min(d0, d1) > 5.0:
                            continue
                        if d0 <= d1:
                            portal_xy = p0
                            outside_xy = p1
                        else:
                            portal_xy = p1
                            outside_xy = p0
                        ramps.append(
                            (portal_xy, outside_xy,
                             piece.length, r_hwy))

                # ── Cluster ramps by portal proximity.  Two
                # carriageways of a divided highway produce two
                # ramps with portals a few metres apart; they get
                # merged into one emission.
                PORTAL_CLUSTER_DIST = 40.0
                clusters = []
                for ramp in ramps:
                    pxy = ramp[0]
                    placed = False
                    for cl in clusters:
                        cref = cl[0][0]
                        if sqrt((cref[0] - pxy[0]) ** 2
                                + (cref[1] - pxy[1]) ** 2
                                ) < PORTAL_CLUSTER_DIST:
                            cl.append(ramp)
                            placed = True
                            break
                    if not placed:
                        clusters.append([ramp])

                for cluster in clusters:
                    # Centroid portal and outside points
                    cpx = sum(r[0][0] for r in cluster) / len(cluster)
                    cpy = sum(r[0][1] for r in cluster) / len(cluster)
                    cox = sum(r[1][0] for r in cluster) / len(cluster)
                    coy = sum(r[1][1] for r in cluster) / len(cluster)

                    # Ramp direction: from outside to portal
                    dxr = cpx - cox
                    dyr = cpy - coy
                    dlen = sqrt(dxr * dxr + dyr * dyr)
                    if dlen < 5.0:
                        continue
                    ramp_dir = (dxr / dlen, dyr / dlen)
                    perp = (-ramp_dir[1], ramp_dir[0])

                    # Combined width: base carriageway + the
                    # perpendicular spread of all member portals
                    # (so a divided highway with two carriageways
                    # gets a rect wide enough to cover both).
                    projs = [
                        ((r[0][0] - cpx) * perp[0]
                         + (r[0][1] - cpy) * perp[1])
                        for r in cluster
                    ]
                    span = max(projs) - min(projs) if projs else 0.0
                    ramp_width = CARRIAGEWAY_BASE_WIDTH + span

                    # Airport elevation at the portal
                    apt_elev = _cifp_surface_elevation(
                        cpx, cpy, max_search=500.0)
                    if apt_elev is None:
                        plon_p, plat_p = to_ll(cpx, cpy)
                        try:
                            apt_elev = tile.dem.alt(
                                (plon_p - tile.lon,
                                 plat_p - tile.lat))
                        except Exception:
                            continue
                    if apt_elev is None:
                        continue

                    # Surface elevation at the outside (OSM tunnel
                    # node) end — DEM at that point.
                    olon, olat = to_ll(cox, coy)
                    try:
                        out_elev = tile.dem.alt(
                            (olon - tile.lon, olat - tile.lat))
                    except Exception:
                        out_elev = apt_elev
                    if out_elev is None:
                        out_elev = apt_elev

                    # HIGH end (outside, surface DEM), LOW end
                    # (portal, apt_elev − DEPTH).
                    elev_high = out_elev
                    elev_low = apt_elev - TUNNEL_DEPTH_DEFAULT

                    # Log the portal
                    _plon, _plat = to_ll(cpx, cpy)
                    UI.vprint(
                        2,
                        "      {}: tunnel portal at ({:.5f}, {:.5f}) "
                        "hwy={} len={:.0f}m width={:.0f}m "
                        "elev {:.1f}→{:.1f}m".format(
                            icao, _plat, _plon, cluster[0][3],
                            dlen, ramp_width,
                            elev_high, elev_low))

                    # ── Ramp rect: sloped from outside to portal ──
                    # _emit_sloped_rect expects the HIGH end to be
                    # passed first (altitude_high applies to its
                    # first two corners).
                    oslon, oslat = to_ll(cox, coy)
                    pslon, pslat = to_ll(cpx, cpy)
                    try:
                        rc = runway_corners(
                            oslat, oslon, pslat, pslon, ramp_width)
                        if rc is None:
                            continue
                        rp = shp_geom.Polygon(
                            [to_m(c[1], c[0]) for c in rc])
                        if not rp.is_valid:
                            rp = rp.buffer(0)
                        if rp.is_empty:
                            continue
                    except Exception:
                        continue
                    _emit_sloped_rect(
                        oslat, oslon, round(elev_high, 1),
                        pslat, pslon, round(elev_low, 1),
                        ramp_width)
                    blend_polys_m.append(rp)
                    n_blend += 1

                    # ── Retaining walls: U-shape around the LOW end
                    # with a small gap to the ramp on every edge.
                    #
                    # The ramp rect sits in a local frame where
                    # "along" = ramp_dir (from outside HIGH to portal
                    # LOW) and "across" = perp.  Walls:
                    #   * two sides along the ramp, offset outward
                    #     by ramp_width/2 + GAP + wall_width/2.
                    #     Each side runs from the outside end to
                    #     the portal end PLUS a short overhang so
                    #     it meets the end-cap without a gap at the
                    #     corner.
                    #   * one end-cap across the portal end, offset
                    #     beyond the portal by GAP + wall_width/2
                    #     (away from the ramp).
                    wall_off = (ramp_width / 2.0
                                + TUNNEL_WALL_GAP
                                + RETAINING_WALL_WIDTH / 2.0)
                    # End-cap runs across the ramp LOW end, offset
                    # past the portal into the airport by
                    # (TUNNEL_WALL_GAP + RETAINING_WALL_WIDTH/2).
                    cap_back = (TUNNEL_WALL_GAP
                                + RETAINING_WALL_WIDTH / 2.0)
                    # Cap center, one wall-width past the portal
                    ccx = cpx + ramp_dir[0] * cap_back
                    ccy = cpy + ramp_dir[1] * cap_back
                    # Cap length = ramp width + both side-wall gaps,
                    # so the cap fits BETWEEN the two side walls
                    # without overlapping them at the corners.
                    cap_len = (ramp_width
                               + 2 * TUNNEL_WALL_GAP)
                    cap_a_x = ccx + perp[0] * cap_len / 2.0
                    cap_a_y = ccy + perp[1] * cap_len / 2.0
                    cap_b_x = ccx - perp[0] * cap_len / 2.0
                    cap_b_y = ccy - perp[1] * cap_len / 2.0
                    cap_rect_corners = None
                    try:
                        alon, alat = to_ll(cap_a_x, cap_a_y)
                        blon, blat = to_ll(cap_b_x, cap_b_y)
                        cap_rect_corners = runway_corners(
                            alat, alon, blat, blon,
                            RETAINING_WALL_WIDTH)
                    except Exception:
                        cap_rect_corners = None

                    side_wall_rects = []
                    for side in (1.0, -1.0):
                        # Side wall start (outside/high end): at the
                        # outside node projected out by wall_off.
                        s_sx = cox + side * perp[0] * wall_off
                        s_sy = coy + side * perp[1] * wall_off
                        # Side wall end: at the portal + cap_back
                        # so it meets the cap cleanly.
                        s_ex = (cpx + ramp_dir[0] * cap_back
                                + side * perp[0] * wall_off)
                        s_ey = (cpy + ramp_dir[1] * cap_back
                                + side * perp[1] * wall_off)
                        try:
                            wslon, wslat = to_ll(s_sx, s_sy)
                            welon, welat = to_ll(s_ex, s_ey)
                            wc = runway_corners(
                                wslat, wslon, welat, welon,
                                RETAINING_WALL_WIDTH)
                            if wc is None:
                                continue
                            wp = shp_geom.Polygon(
                                [to_m(c[1], c[0]) for c in wc])
                            if not wp.is_valid:
                                wp = wp.buffer(0)
                            if wp.is_empty:
                                continue
                        except Exception:
                            continue
                        side_wall_rects.append(
                            (wslat, wslon, welat, welon, wp))

                    # Emit each side wall flat at apt_elev
                    for wslat, wslon, welat, welon, wp in side_wall_rects:
                        _emit_flat_poly(
                            [(c[1], c[0])
                             for c in runway_corners(
                                 wslat, wslon, welat, welon,
                                 RETAINING_WALL_WIDTH)],
                            round(apt_elev, 1))
                        blend_polys_m.append(wp)
                        n_blend += 1

                    # Emit the end-cap (flat, wrapped around LOW end)
                    if cap_rect_corners is not None:
                        try:
                            cp = shp_geom.Polygon(
                                [to_m(c[1], c[0])
                                 for c in cap_rect_corners])
                            if not cp.is_valid:
                                cp = cp.buffer(0)
                            if not cp.is_empty:
                                _emit_flat_poly(
                                    [(c[1], c[0])
                                     for c in cap_rect_corners],
                                    round(apt_elev, 1))
                                blend_polys_m.append(cp)
                                n_blend += 1
                        except Exception:
                            pass

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
        if DEBUG_TAXI_ONLY or DEBUG_REPLACE_LEGACY_PAVEMENT:
            raise RuntimeError("skip Phase F")
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

                        # Build the ditch rect polygon in meter
                        # space and reject it if any part extends
                        # past the infield (into adjacent pavement
                        # triangles/pads).  Strict containment —
                        # guarantees the ditch doesn't overlap any
                        # already-emitted geometry.
                        try:
                            perp_x = -best_dir[1]
                            perp_y = best_dir[0]
                            hw = DRAIN_DITCH_WIDTH / 2.0
                            ditch_poly = shp_geom.Polygon([
                                (sx + perp_x * hw, sy + perp_y * hw),
                                (ex + perp_x * hw, ey + perp_y * hw),
                                (ex - perp_x * hw, ey - perp_y * hw),
                                (sx - perp_x * hw, sy - perp_y * hw),
                            ])
                            # Allow 0.5 m² boundary slop; anything
                            # larger means the ditch crosses the
                            # infield boundary.
                            outside = ditch_poly.difference(ip).area
                            if outside > 0.5:
                                continue
                        except Exception:
                            pass

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

    # Final summary counts.  Phase C0 (apt.dat taxiway rect-chain)
    # and Phase C1 (legacy OSM-centerline taxiway path) both feed
    # the same two buckets: flat and sloped taxiway quads.
    total_twy_flat = n_twy_flat + n_apt_twy_flat_rects
    total_twy_sloped = n_twy_sloped + n_apt_twy_rects
    total_twy_shapes = total_twy_flat + total_twy_sloped
    UI.vprint(1, "    {}: {} taxiway ({} flat, {} sloped)"
              " + {} apron ({} flat, {} complex)"
              " + {} building pads"
              " + {} junction triangles"
              " + {} boundary/road shapes + {} drainage low-points".format(
                  icao, total_twy_shapes, total_twy_flat, total_twy_sloped,
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
