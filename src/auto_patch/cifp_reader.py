"""CIFP (Coded Instrument Flight Procedure) data reader.

Parses ARINC 424 (CIFP) airport data files to extract precise runway
threshold elevations and coordinates, plus utilities for locating CIFP
files within an X-Plane installation and discovering which airports
fall in a given tile.

Public API:
    parse_cifp_lat(s)             — ARINC 424 latitude string -> float
    parse_cifp_lon(s)             — ARINC 424 longitude string -> float
    parse_cifp_file(filepath)     — full CIFP RWY-record extraction
    find_aptdat(cifp_path)        — locate the apt.dat alongside a CIFP file
    xplane_root_from_cifp_path(p) — derive X-Plane install root from CIFP dir
    discover_cifp_airports(dir)   — enumerate all airports in a CIFP directory
    airport_in_tile(rwys, lat, lon) — does any runway fall in this 1°×1° tile?

Used by:
    * O4_Auto_Patch.generate_auto_patches (the tile-level driver)
    * O4_Airport_Pavement_Builder (Phase-2 elevation anchoring)
    * tests/conftest.py (test airport discovery)
"""
import os
import re

from math import floor

import O4_UI_Utils as UI

__all__ = [
    "parse_cifp_lat",
    "parse_cifp_lon",
    "parse_cifp_file",
    "find_aptdat",
    "xplane_root_from_cifp_path",
    "discover_cifp_airports",
    "airport_in_tile",
]


FT_TO_M = 0.3048


# ──────────────────────────────────────────────────────────────────────
# CIFP Coordinate Parsing
# ──────────────────────────────────────────────────────────────────────
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



# ──────────────────────────────────────────────────────────────────────
# CIFP Runway Data
# ──────────────────────────────────────────────────────────────────────
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



# ──────────────────────────────────────────────────────────────────────
# CIFP File / X-Plane Install Layout Helpers
# ──────────────────────────────────────────────────────────────────────
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



# ──────────────────────────────────────────────────────────────────────
# Tile-level Airport Discovery
# ──────────────────────────────────────────────────────────────────────
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


