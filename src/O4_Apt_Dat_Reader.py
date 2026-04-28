"""Parser for X-Plane ``apt.dat`` airport data files.

Loads runway, pavement (taxiway / apron / ramp) and boundary geometry
for a single airport.  This is the *patch-mesh authoritative* source
for airport pavement shapes — apt.dat polygons are exactly what the
X-Plane simulator renders as the pavement texture, so elevation
patches generated against these polygons align perfectly with the
ground texture (no visible seams).

Compared to OSM aerodrome data:

* apt.dat polygons are **disjoint by construction** — no overlapping
  shapes to subtract, no precision-drift artefacts.
* Taxiways are stored as **outline polygons**, not centerlines that
  we have to buffer to a guessed width.
* Curved taxiway shoulders use **Bezier control points** (rows 112,
  114) which we sample into polygon vertices.
* Each pavement carries its **surface type**, **roughness** and
  **orientation** so downstream code can apply per-surface grade
  rules.

Compared to CIFP, apt.dat does NOT have per-runway-threshold
elevations — runway row 100 only stores lat/lon, width, and
displaced-threshold offsets.  CIFP is still the source of truth for
runway elevations.

The parser is read-only and side-effect-free: given a path to an
``apt.dat`` file and an ICAO code, it returns an :class:`Airport`
object containing the parsed geometry.  Use :func:`find_airport_apt_dat`
to locate the right ``apt.dat`` for a given airport, preferring a
per-airport Custom Scenery pack over the global one.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Iterator, List, Optional, Tuple

from shapely.geometry import Polygon
from shapely.ops import unary_union


# ──────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────
FT_TO_M = 0.3048

# Default number of straight-line segments to sample each Bezier curve
# into.  4 produces a visibly smooth corner without exploding the
# vertex count.  Tunable via the load_airport(..., bezier_segments=N)
# parameter.
DEFAULT_BEZIER_SEGMENTS = 4

# Row type codes (X-Plane apt.dat 1100 / 1200 spec).
ROW_AIRPORT_HEADER = 1
ROW_RUNWAY = 100
ROW_HELIPAD = 102
ROW_PAVEMENT_HEADER = 110
ROW_NODE = 111
ROW_NODE_BEZIER = 112
ROW_CLOSE = 113
ROW_CLOSE_BEZIER = 114
ROW_BOUNDARY_HEADER = 130


# ──────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────
@dataclass
class Runway:
    """One paired runway (apt.dat row 100)."""
    desig_a: str
    desig_b: str
    lat_a: float
    lon_a: float
    lat_b: float
    lon_b: float
    width_m: float
    surface_code: int
    displaced_a_m: float
    displaced_b_m: float
    blast_a_m: float = 0.0   # blast pad / overrun length beyond end a
    blast_b_m: float = 0.0   # blast pad / overrun length beyond end b


@dataclass
class Pavement:
    """One pavement polygon (apt.dat row 110)."""
    polygon: Polygon            # exterior + interior holes
    surface_code: int
    roughness: float            # 0.0 (smooth) … 1.0 (rough)
    orientation: float          # texture rotation in degrees from N
    name: str = ""              # pavement label, e.g. "TWY A", "RAMP 1"


@dataclass
class Airport:
    """Parsed airport geometry from one apt.dat block."""
    icao: str
    name: str
    reference_elev_ft: int      # row 1 elevation, in feet (0 if absent)
    runways: List[Runway] = field(default_factory=list)
    pavements: List[Pavement] = field(default_factory=list)
    boundary: Optional[Polygon] = None
    source_path: str = ""

    @property
    def reference_elev_m(self) -> float:
        return self.reference_elev_ft * FT_TO_M


# ──────────────────────────────────────────────────────────────────────
# Public API: locating the right apt.dat
# ──────────────────────────────────────────────────────────────────────
def find_airport_apt_dat(xplane_root: str, icao: str) -> Optional[str]:
    """Locate the most-specific ``apt.dat`` containing the given ICAO.

    Search priority:

    1. Per-airport packs in ``<X-Plane>/Custom Scenery/<pack>/Earth nav
       data/apt.dat``.  Any pack whose apt.dat starts an airport block
       for the ICAO wins.  This is what the user almost always wants:
       a custom-built scenery for that specific airport.
    2. ``<X-Plane>/Custom Scenery/Global Airports/Earth nav data/apt.dat``.
       The community-curated global file shipped with X-Plane.
    3. ``<X-Plane>/Resources/default scenery/default apt dat/Earth nav
       data/apt.dat``.  Laminar's stock fallback.

    Returns the path to the chosen apt.dat, or ``None`` if no apt.dat
    on the search path contains a header for the ICAO.

    Notes:
        * The check is "does the file contain a row 1 line whose ICAO
          field matches?"  We don't actually parse the airport — that
          would be wasteful when scanning many packs.
        * The "Global Airports" pack is itself a Custom Scenery
          directory; we explicitly defer it to step 2 so per-airport
          packs win.
    """
    if not xplane_root or not os.path.isdir(xplane_root):
        return None

    icao = icao.strip().upper()
    if not icao:
        return None

    custom_scenery = os.path.join(xplane_root, "Custom Scenery")
    # X-Plane 11 layout:
    global_pack_v11 = os.path.join(
        custom_scenery, "Global Airports", "Earth nav data", "apt.dat")
    # X-Plane 12 layout (shipped pack moved to Global Scenery):
    global_pack_v12 = os.path.join(
        xplane_root, "Global Scenery", "Global Airports",
        "Earth nav data", "apt.dat")
    default_pack = os.path.join(
        xplane_root, "Resources", "default scenery",
        "default apt dat", "Earth nav data", "apt.dat")

    # Two-pass search: prefer files that contain proper row-110
    # pavement polygons (our pipeline needs those), then fall back to
    # any file that just contains the airport.  This handles e.g. the
    # KBNA Custom Scenery pack which uses row-120 linear features but
    # no row-110 pavements — the Global pack is the right source for
    # pavement geometry there.
    custom_packs: List[str] = []
    if os.path.isdir(custom_scenery):
        for entry in sorted(os.listdir(custom_scenery)):
            if entry == "Global Airports":
                continue
            pack_apt = os.path.join(
                custom_scenery, entry, "Earth nav data", "apt.dat")
            if os.path.isfile(pack_apt) and _file_has_airport(pack_apt, icao):
                custom_packs.append(pack_apt)

    candidates: List[str] = list(custom_packs)
    for cand in (global_pack_v11, global_pack_v12):
        if os.path.isfile(cand) and _file_has_airport(cand, icao):
            candidates.append(cand)
    if os.path.isfile(default_pack) and _file_has_airport(default_pack, icao):
        candidates.append(default_pack)

    # First pass: prefer the most-specific source that ALSO has pavement.
    for cand in candidates:
        if _file_has_airport_with_pavement(cand, icao):
            return cand
    # Second pass: any file with the airport header (lets the rest of
    # the pipeline at least parse runways even if no pavements exist).
    if candidates:
        return candidates[0]
    return None


# ──────────────────────────────────────────────────────────────────────
# Public API: parsing an airport block
# ──────────────────────────────────────────────────────────────────────
def load_airport(
    aptdat_path: str,
    icao: str,
    bezier_segments: int = DEFAULT_BEZIER_SEGMENTS,
) -> Optional[Airport]:
    """Parse the airport block for ``icao`` out of ``aptdat_path``.

    Args:
        aptdat_path: filesystem path to an apt.dat file.
        icao: 4-letter airport code (case-insensitive).
        bezier_segments: how many straight-line segments to subdivide
            each Bezier curve into.  4 is a good default for taxiway
            corners; raise it for tight curves.

    Returns:
        An :class:`Airport` object, or ``None`` if the airport block
        could not be found in the file.
    """
    if not aptdat_path or not os.path.isfile(aptdat_path):
        return None

    block = _read_airport_block(aptdat_path, icao)
    if block is None:
        return None

    header = block[0]
    # Row 1 format: ``1 elevation_ft tower_height beacon_type ICAO airport_name``
    # The name is everything after the ICAO and may contain spaces.
    parts = header.split(maxsplit=5)
    try:
        ref_elev = int(parts[1])
    except (IndexError, ValueError):
        ref_elev = 0
    name = parts[5] if len(parts) > 5 else ""

    airport = Airport(
        icao=icao.upper(),
        name=name.strip(),
        reference_elev_ft=ref_elev,
        source_path=aptdat_path,
    )

    pavement_rows: List[List[str]] = []
    boundary_rows: List[List[str]] = []
    in_pavement = False
    in_boundary = False

    def flush_pavement():
        if pavement_rows:
            pav = _parse_pavement(pavement_rows, bezier_segments)
            if pav is not None:
                airport.pavements.append(pav)
            pavement_rows.clear()

    def flush_boundary():
        if boundary_rows:
            poly = _parse_boundary(boundary_rows, bezier_segments)
            if poly is not None:
                # If we already have a boundary, union with the new one.
                if airport.boundary is None:
                    airport.boundary = poly
                else:
                    try:
                        merged = unary_union([airport.boundary, poly])
                        if isinstance(merged, Polygon):
                            airport.boundary = merged
                    except Exception:
                        pass
            boundary_rows.clear()

    for line in block[1:]:
        stripped = line.strip()
        if not stripped:
            continue
        toks = stripped.split()
        try:
            row_type = int(toks[0])
        except ValueError:
            continue

        # Pavement / boundary blocks accumulate consecutive node rows.
        if row_type == ROW_PAVEMENT_HEADER:
            flush_pavement()
            flush_boundary()
            in_pavement = True
            in_boundary = False
            pavement_rows.append(toks)
            continue
        if row_type == ROW_BOUNDARY_HEADER:
            flush_pavement()
            flush_boundary()
            in_pavement = False
            in_boundary = True
            boundary_rows.append(toks)
            continue
        if row_type in (ROW_NODE, ROW_NODE_BEZIER,
                        ROW_CLOSE, ROW_CLOSE_BEZIER):
            if in_pavement:
                pavement_rows.append(toks)
            elif in_boundary:
                boundary_rows.append(toks)
            continue

        # Anything else terminates the current pavement / boundary
        # block (and contributes its own data).
        flush_pavement()
        flush_boundary()
        in_pavement = False
        in_boundary = False

        if row_type == ROW_RUNWAY:
            rwy = _parse_runway(toks)
            if rwy is not None:
                airport.runways.append(rwy)

    # Final flush in case the block ends mid-pavement.
    flush_pavement()
    flush_boundary()

    return airport


# ──────────────────────────────────────────────────────────────────────
# Internal: file scanning
# ──────────────────────────────────────────────────────────────────────
def find_all_airport_apt_dats(xplane_root: str,
                              icao: str) -> List[str]:
    """Return EVERY apt.dat path under ``xplane_root`` that contains
    a row-1 header for ``icao`` (any pack — Custom Scenery,
    Global Airports, default).

    Different packs commonly carry different geometry for the same
    airport: a community pack might add row-110 pavement that the
    Global pack lacks, AND a custom DSF with draped polygons that
    neither has.  Callers that want the union of all available
    pavement geometry walk this list.
    """
    if not xplane_root or not os.path.isdir(xplane_root):
        return []
    icao = icao.strip().upper()
    if not icao:
        return []
    out: List[str] = []
    custom_scenery = os.path.join(xplane_root, "Custom Scenery")
    if os.path.isdir(custom_scenery):
        for entry in sorted(os.listdir(custom_scenery)):
            pack_apt = os.path.join(
                custom_scenery, entry, "Earth nav data", "apt.dat")
            if (os.path.isfile(pack_apt)
                    and _file_has_airport(pack_apt, icao)):
                out.append(pack_apt)
    global_v11 = os.path.join(
        xplane_root, "Custom Scenery", "Global Airports",
        "Earth nav data", "apt.dat")
    global_v12 = os.path.join(
        xplane_root, "Global Scenery", "Global Airports",
        "Earth nav data", "apt.dat")
    for cand in (global_v11, global_v12):
        if (os.path.isfile(cand)
                and _file_has_airport(cand, icao)
                and cand not in out):
            out.append(cand)
    default = os.path.join(
        xplane_root, "Resources", "default scenery",
        "default apt dat", "Earth nav data", "apt.dat")
    if (os.path.isfile(default)
            and _file_has_airport(default, icao)
            and default not in out):
        out.append(default)
    return out


# Process-wide cache for apt.dat header scans.  Keyed by
# (path, mtime_ns, size) so a stale entry is invalidated automatically
# if the file is rewritten during the same process run.  Populated
# lazily by :func:`_index_apt_dat` on first access; subsequent
# `_file_has_airport` / `_file_has_airport_with_pavement` calls become
# O(1) dict lookups.
#
# Why this matters: the auto-patch pipeline calls
# ``build_airport_pavement(icao, ...)`` once per airport in a tile, and
# each call invokes ``find_airport_apt_dat`` and
# ``find_all_airport_apt_dats``.  Before the cache, every airport
# invocation re-scanned every apt.dat file in the X-Plane install
# (Custom Scenery + Global + default), which on a typical setup with
# ~25 airports per tile pulled ~10 GB through the line scanner per
# tile build.  After the cache, each apt.dat is scanned exactly once
# per process.
_APT_DAT_INDEX_CACHE: dict = {}


def _index_apt_dat(aptdat_path: str) -> Tuple[frozenset, frozenset]:
    """Return ``(icaos_present, icaos_with_pavement)`` for the file.

    Both sets are uppercase ICAO codes.  An entry in
    ``icaos_with_pavement`` means the airport block has at least one
    row 110 (pavement header).  Result is cached process-wide; if the
    file is rewritten (mtime / size changes) the cache entry is
    invalidated and the file is rescanned.
    """
    try:
        st = os.stat(aptdat_path)
    except OSError:
        return frozenset(), frozenset()
    key = (aptdat_path, st.st_mtime_ns, st.st_size)
    cached = _APT_DAT_INDEX_CACHE.get(key)
    if cached is not None:
        return cached
    # Drop any stale entry for this path (different mtime/size).
    for k in [k for k in _APT_DAT_INDEX_CACHE if k[0] == aptdat_path]:
        _APT_DAT_INDEX_CACHE.pop(k, None)

    icaos = set()
    with_pavement = set()
    current: Optional[str] = None
    saw_pavement_in_current = False
    try:
        with open(aptdat_path, "r", encoding="utf-8",
                  errors="replace") as f:
            for line in f:
                stripped = line.lstrip()
                if stripped.startswith("1 ") or stripped.startswith("1\t"):
                    parts = stripped.split()
                    if len(parts) >= 5 and parts[0] == "1":
                        # Close out the previous airport block.
                        if current is not None and saw_pavement_in_current:
                            with_pavement.add(current)
                        current = parts[4].upper()
                        saw_pavement_in_current = False
                        icaos.add(current)
                        continue
                if (current is not None
                        and not saw_pavement_in_current
                        and (stripped.startswith("110 ")
                             or stripped.startswith("110\t"))):
                    saw_pavement_in_current = True
        # Close out the last block at EOF.
        if current is not None and saw_pavement_in_current:
            with_pavement.add(current)
    except Exception:
        # Cache an empty result so we don't re-attempt every call.
        result = (frozenset(), frozenset())
        _APT_DAT_INDEX_CACHE[key] = result
        return result

    result = (frozenset(icaos), frozenset(with_pavement))
    _APT_DAT_INDEX_CACHE[key] = result
    return result


def _file_has_airport(aptdat_path: str, icao: str) -> bool:
    """Return True if `aptdat_path` contains a row 1 header for ICAO.

    Backed by :func:`_index_apt_dat`'s process-wide cache; the file is
    fully scanned at most once per (path, mtime, size).
    """
    icaos, _ = _index_apt_dat(aptdat_path)
    return icao.upper() in icaos


def _file_has_airport_with_pavement(aptdat_path: str, icao: str) -> bool:
    """Return True if `aptdat_path` contains a row 1 header for ICAO
    AND the airport block has at least one row 110 (pavement header).

    Some Custom Scenery packs (e.g. KBNA) replace pavement polygons
    with linear-feature markup (row 120 + 111 nodes), leaving the
    airport block with 0 row-110 records.  Our pavement pipeline
    needs row-110 polygons to compute the residue/junction set, so
    such packs are unusable for pavement geometry — we fall back to
    the Global apt.dat which does have proper row-110 pavements.

    Backed by :func:`_index_apt_dat`'s process-wide cache.
    """
    _, with_pavement = _index_apt_dat(aptdat_path)
    return icao.upper() in with_pavement


def _read_airport_block(aptdat_path: str, icao: str) -> Optional[List[str]]:
    """Return all lines from the row-1 header for `icao` up to (but
    not including) the next row-1 header.  None if not found.
    """
    icao = icao.upper()
    block: List[str] = []
    in_block = False
    try:
        with open(aptdat_path, "r", encoding="utf-8",
                  errors="replace") as f:
            for line in f:
                if line.startswith("1 ") or line.startswith("1\t"):
                    parts = line.split()
                    if len(parts) >= 5 and parts[0] == "1":
                        if in_block:
                            # Reached the next airport.
                            return block
                        if parts[4].upper() == icao:
                            in_block = True
                            block.append(line)
                            continue
                if in_block:
                    block.append(line)
    except Exception:
        return None
    return block if in_block else None


# ──────────────────────────────────────────────────────────────────────
# Internal: row parsers
# ──────────────────────────────────────────────────────────────────────
def _parse_runway(toks: List[str]) -> Optional[Runway]:
    """Parse an apt.dat row 100 into a Runway.  Format:

    ``100 width surface shoulder smoothness centerline edge_lights distance_signs
         <end_a:9> <end_b:9>``

    Each end-of-runway block is 9 tokens:
    ``desig lat lon displaced blastpad markings approach_lights tdz_lights reil``

    ``blastpad`` (index 4 within the end block) is the length in
    metres of the blast pad / stopway / overrun surface beyond the
    threshold on that end.
    """
    if len(toks) < 25:
        return None
    try:
        width_m = float(toks[1])
        surface_code = int(toks[2])
        # toks[3] = shoulder, toks[4] = smoothness, toks[5] = centerline,
        # toks[6] = edge_lights, toks[7] = distance_signs
        end_a = toks[8:17]   # 9 fields
        end_b = toks[17:26]
        desig_a = end_a[0]
        lat_a = float(end_a[1])
        lon_a = float(end_a[2])
        displaced_a_m = float(end_a[3])
        blast_a_m = float(end_a[4])
        desig_b = end_b[0]
        lat_b = float(end_b[1])
        lon_b = float(end_b[2])
        displaced_b_m = float(end_b[3])
        blast_b_m = float(end_b[4])
    except (ValueError, IndexError):
        return None

    return Runway(
        desig_a=desig_a, desig_b=desig_b,
        lat_a=lat_a, lon_a=lon_a, lat_b=lat_b, lon_b=lon_b,
        width_m=width_m, surface_code=surface_code,
        displaced_a_m=displaced_a_m, displaced_b_m=displaced_b_m,
        blast_a_m=blast_a_m, blast_b_m=blast_b_m,
    )


def _parse_pavement(rows: List[List[str]],
                    bezier_segments: int) -> Optional[Pavement]:
    """Parse a row-110 header + node rows into a Pavement.

    The first contour (terminated by 113/114) is the exterior; any
    subsequent contours within the same pavement are interior holes.
    """
    if not rows:
        return None
    header = rows[0]
    try:
        surface_code = int(header[1])
        roughness = float(header[2])
        orientation = float(header[3])
    except (IndexError, ValueError):
        return None
    name = " ".join(header[4:]) if len(header) > 4 else ""

    contours = _split_contours(rows[1:])
    if not contours:
        return None

    rings = []
    for contour in contours:
        ring = _interpolate_contour(contour, bezier_segments)
        if len(ring) >= 3:
            # Close the ring explicitly.
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            rings.append(ring)
    if not rings:
        return None

    try:
        polygon = Polygon(rings[0], rings[1:] if len(rings) > 1 else None)
    except Exception:
        return None
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    if polygon.is_empty:
        return None
    # buffer(0) on a self-intersecting source polygon can return a
    # MultiPolygon (the cleaned result split into disjoint pieces).
    # Take the largest component — at SPJC's custom apt.dat about 8
    # of 51 pavements take this path; the dropped slivers are tiny
    # geometric artefacts of the source data, not real pavement.
    if not isinstance(polygon, Polygon):
        if hasattr(polygon, "geoms"):
            try:
                polygon = max(polygon.geoms, key=lambda g: g.area)
            except Exception:
                return None
        else:
            return None
        if polygon.is_empty or not isinstance(polygon, Polygon):
            return None

    return Pavement(
        polygon=polygon,
        surface_code=surface_code,
        roughness=roughness,
        orientation=orientation,
        name=name.strip(),
    )


def _parse_boundary(rows: List[List[str]],
                    bezier_segments: int) -> Optional[Polygon]:
    """Parse a row-130 header + node rows into a boundary Polygon.

    Boundaries follow the same node row format as pavements.  We
    treat the first contour as the exterior and any extra contours
    as interior holes (rare for boundaries but allowed by the spec).
    """
    if not rows:
        return None
    contours = _split_contours(rows[1:])
    if not contours:
        return None
    rings = []
    for contour in contours:
        ring = _interpolate_contour(contour, bezier_segments)
        if len(ring) >= 3:
            if ring[0] != ring[-1]:
                ring.append(ring[0])
            rings.append(ring)
    if not rings:
        return None
    try:
        poly = Polygon(rings[0], rings[1:] if len(rings) > 1 else None)
    except Exception:
        return None
    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty or not isinstance(poly, Polygon):
        return None
    return poly


def _split_contours(node_rows: List[List[str]]) -> List[List[List[str]]]:
    """Walk a list of 111/112/113/114 rows and split into contours.

    A contour starts at the first row after the header (or after the
    previous contour's closing row) and ends at the next 113 / 114
    closing row.  Each returned contour is a list of node rows
    INCLUDING its closing 113/114 row.
    """
    contours: List[List[List[str]]] = []
    current: List[List[str]] = []
    for row in node_rows:
        if not row:
            continue
        try:
            row_type = int(row[0])
        except ValueError:
            continue
        if row_type not in (ROW_NODE, ROW_NODE_BEZIER,
                            ROW_CLOSE, ROW_CLOSE_BEZIER):
            continue
        current.append(row)
        if row_type in (ROW_CLOSE, ROW_CLOSE_BEZIER):
            contours.append(current)
            current = []
    # Drop a trailing un-closed contour.
    return contours


# ──────────────────────────────────────────────────────────────────────
# Internal: Bezier interpolation
# ──────────────────────────────────────────────────────────────────────
def _node_xy(row: List[str]) -> Tuple[float, float]:
    """Return (lon, lat) for a node row (we use lon-first internally
    so shapely Polygons get the standard (x, y) order)."""
    return (float(row[2]), float(row[1]))


def _node_ctrl(row: List[str]) -> Optional[Tuple[float, float]]:
    """Return the Bezier control point for a 112/114 node, or None
    for a plain 111/113 node.
    """
    try:
        rt = int(row[0])
    except ValueError:
        return None
    if rt not in (ROW_NODE_BEZIER, ROW_CLOSE_BEZIER):
        return None
    try:
        return (float(row[4]), float(row[3]))
    except (IndexError, ValueError):
        return None


def _quadratic_bezier(p0, p1, p2, n_segments):
    """Sample a quadratic Bezier from p0 → p2 with control p1.
    Returns a list of (n_segments + 1) points starting at p0 and
    ending at p2.
    """
    pts = []
    for i in range(n_segments + 1):
        t = i / n_segments
        omt = 1.0 - t
        x = omt * omt * p0[0] + 2 * omt * t * p1[0] + t * t * p2[0]
        y = omt * omt * p0[1] + 2 * omt * t * p1[1] + t * t * p2[1]
        pts.append((x, y))
    return pts


def _cubic_bezier(p0, p1, p2, p3, n_segments):
    """Sample a cubic Bezier from p0 → p3 with controls p1, p2.
    Returns a list of (n_segments + 1) points.
    """
    pts = []
    for i in range(n_segments + 1):
        t = i / n_segments
        omt = 1.0 - t
        b0 = omt * omt * omt
        b1 = 3 * omt * omt * t
        b2 = 3 * omt * t * t
        b3 = t * t * t
        x = b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0]
        y = b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1]
        pts.append((x, y))
    return pts


def _mirror(point, anchor):
    """Reflect `point` through `anchor`."""
    return (2 * anchor[0] - point[0], 2 * anchor[1] - point[1])


def _interpolate_contour(contour: List[List[str]],
                         bezier_segments: int) -> List[Tuple[float, float]]:
    """Convert a contour (list of 111/112 rows ending in 113/114)
    into a flat list of (x, y) polygon vertices, sampling Bezier
    curves into straight-line segments.

    Bezier convention used here (matches X-Plane apt.dat 1100 spec):

    For each consecutive pair of nodes A → B:

    * If neither carries a control point: straight line A → B.
    * If only A has a control point (a 112 followed by a 111/113):
      quadratic Bezier from A through ctrl_a to B.
    * If only B has a control point (a 111 followed by a 112/114):
      quadratic Bezier from A through (mirror of ctrl_b across B) to B.
    * If both A and B carry control points: cubic Bezier from A,
      with controls ctrl_a and (mirror of ctrl_b across B), to B.
    """
    n = len(contour)
    if n < 2:
        return []

    # Build a closed ring of nodes (the closing 113/114 brings us
    # back to the first vertex; segment "last → first" closes the
    # contour).
    ring_nodes = list(contour)

    out: List[Tuple[float, float]] = []
    for i in range(n):
        a_row = ring_nodes[i]
        b_row = ring_nodes[(i + 1) % n]
        a_xy = _node_xy(a_row)
        b_xy = _node_xy(b_row)
        a_ctrl = _node_ctrl(a_row)
        b_ctrl = _node_ctrl(b_row)

        # Append A only (B will be appended by the next iteration).
        if not out or out[-1] != a_xy:
            out.append(a_xy)

        if a_ctrl is None and b_ctrl is None:
            # Straight line — nothing to interpolate, B will be added
            # next iteration.
            continue

        if a_ctrl is not None and b_ctrl is None:
            curve = _quadratic_bezier(a_xy, a_ctrl, b_xy, bezier_segments)
        elif a_ctrl is None and b_ctrl is not None:
            mirrored = _mirror(b_ctrl, b_xy)
            curve = _quadratic_bezier(a_xy, mirrored, b_xy, bezier_segments)
        else:
            mirrored = _mirror(b_ctrl, b_xy)
            curve = _cubic_bezier(a_xy, a_ctrl, mirrored, b_xy,
                                  bezier_segments)
        # Drop the first point (= a_xy, already in out) and the last
        # (= b_xy, will be appended next iteration).  Append only the
        # interior curve samples.
        for pt in curve[1:-1]:
            if not out or out[-1] != pt:
                out.append(pt)
    return out


# ──────────────────────────────────────────────────────────────────────
# Aggregate helpers
# ──────────────────────────────────────────────────────────────────────
def airport_pavement_summary(airport: Airport) -> str:
    """Short multi-line summary of an Airport's parsed contents.
    Useful for diagnostic logging during integration.
    """
    lines = [
        "Airport {}: {}".format(airport.icao, airport.name),
        "  ref elev: {} ft ({:.1f} m)".format(
            airport.reference_elev_ft, airport.reference_elev_m),
        "  source:   {}".format(airport.source_path),
        "  runways:  {}".format(len(airport.runways)),
        "  pavements: {}".format(len(airport.pavements)),
        "  boundary: {}".format(
            "yes ({:.0f} m² in lat-lon space)".format(
                airport.boundary.area * 12_345_679_000.0)
            if airport.boundary is not None else "no"),
    ]
    return "\n".join(lines)
