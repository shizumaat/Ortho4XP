"""Build the full role-classified pavement layout for an airport.

The target of this module is to auto-generate the shape + role
layout that a human would hand-craft in JOSM for an airport, within
**0.5 m symmetric vertex tolerance** against a reference target OSM
file.  See ``tests/fixtures/SPJC_target.osm`` and ``SPLP_target.osm``
for the authoritative schema.

Phase 1 (this file) emits **geometry + role only**; no elevation
data is produced.  Phase 2 will layer elevations on top using the
same shape topology.

Output roles:

* ``runway`` — 4-vertex rect from apt.dat row 100.
* ``primary_parallel`` — long taxi parallel to the runway, close by.
* ``secondary_parallel`` — parallel taxi, farther from runway or
  shorter than primary threshold.
* ``stub`` — short connector to a primary/secondary parallel.
* ``cross_connector`` — perpendicular between two parallels.
* ``apron`` — concave pavement area (terminal apron, engine test).
* ``terminal`` — building pad polygon.
* ``junction`` — explicit polygon filling the gap between two
  adjacent rects, sharing their corner vertices exactly.

Usage::

    from O4_Airport_Pavement_Builder import build_airport_pavement
    layout = build_airport_pavement("SPJC", "/path/to/X-Plane 12")
    layout.to_osm("/tmp/SPJC_auto.osm")

The layout's shapes all live in meter space anchored at
``layout.anchor = (lat0, lon0)``.
"""
from __future__ import annotations

import bz2
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from shapely.geometry import LineString, MultiLineString, Point, Polygon
from shapely.ops import (
    linemerge, nearest_points, transform as shp_transform, unary_union)

import O4_Apt_Dat_Reader as APR
import O4_Pavement_Classifier as PC
import O4_Pavement_Strips as PS


# ──────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────
R_EARTH = 6_378_137.0
SHARED_VERTEX_TOL_M = 0.5    # snap vertices closer than this together
RUNWAY_INSIDE_APRON_FRAC = 0.95  # if ≥95% of a runway segment lies
                                  # inside an apt.dat/DSF apron
                                  # polygon (and that polygon is
                                  # much larger than the segment —
                                  # see RUNWAY_APRON_AREA_RATIO),
                                  # treat the segment as apron-
                                  # merged and drop the separate rect
RUNWAY_APRON_AREA_RATIO = 3.0     # the containing apt.dat/DSF
                                   # polygon must be ≥ this ratio
                                   # times the segment area to count
                                   # as an apron.  Normal runway
                                   # segments sit in runway-shaped
                                   # polygons whose area is only
                                   # marginally larger than the
                                   # segment itself.

ROLE_RUNWAY = "runway"
ROLE_PRIMARY_PARALLEL = PS.ROLE_PRIMARY_PARALLEL
ROLE_SECONDARY_PARALLEL = PS.ROLE_SECONDARY_PARALLEL
ROLE_STUB = PS.ROLE_STUB
ROLE_CROSS_CONNECTOR = PS.ROLE_CROSS_CONNECTOR
ROLE_APRON = PS.ROLE_APRON
ROLE_TERMINAL = "terminal"
ROLE_JUNCTION = "junction"
ROLE_BOUNDARY = "boundary"

# TEMP 2026-04-20: when False, the builder only emits rects +
# runways + terminals + aprons, suppressing all junction polygons.
# User requested this while iterating on rect correctness.
EMIT_JUNCTIONS = True
# TEMP 2026-04-21: when False, aprons are also suppressed so we can
# focus exclusively on getting taxiway rects right.  Junctions and
# aprons are treated interchangeably for now.
EMIT_APRONS = False

AEROWAY_FOR_ROLE = {
    ROLE_RUNWAY: "runway",
    ROLE_PRIMARY_PARALLEL: "taxiway",
    ROLE_SECONDARY_PARALLEL: "taxiway",
    ROLE_STUB: "taxiway",
    ROLE_CROSS_CONNECTOR: "taxiway",
    ROLE_JUNCTION: "taxiway",
    ROLE_APRON: "apron",
    ROLE_TERMINAL: "building",
    ROLE_BOUNDARY: "aerodrome",
}


# ──────────────────────────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────────────────────────

@dataclass
class BuiltShape:
    """A single emitted shape: polygon + classification tags.

    Polygons live in meter space anchored at the layout's origin.
    ``ref`` is optional (runway designator, taxi ref from OSM, or
    generated label).  ``source_axis`` is kept on taxi rects for
    elevation sampling along their axis.

    Phase-2 elevation: exactly one of these options is set at any
    time:
      * all four None (no elevation yet)
      * only ``altitude`` set (flat polygon at that elevation, m)
      * ``altitude_high`` + ``altitude_low`` set (linearly sloped
        between the two parallel edges).  Rects use this.
      * ``node_altitudes`` set (per-vertex elevation list, one
        value per ring vertex INCLUDING the closing repeat —
        i.e. len(node_altitudes) == len(closed_ring_nids)).
        Used for triangulated junction polygons that slope in
        more than one direction.
    Runway segments carry altitude_high/low per the legacy patch
    convention.
    """
    polygon: Polygon
    role: str
    ref: str = ""
    source_axis: Optional[LineString] = None
    altitude: Optional[float] = None
    altitude_high: Optional[float] = None
    altitude_low: Optional[float] = None
    node_altitudes: Optional[List[float]] = None


@dataclass
class PavementLayout:
    icao: str
    anchor: Tuple[float, float]          # (lat0, lon0)
    shapes: List[BuiltShape] = field(default_factory=list)
    # ancillary:
    airport_boundary: Optional[Polygon] = None
    runway_union: Optional[Polygon] = None

    # ---- coordinate helpers ------------------------------------------
    def m_to_ll(self, x: float, y: float) -> Tuple[float, float]:
        lat0, lon0 = self.anchor
        cos0 = math.cos(math.radians(lat0))
        lon = lon0 + math.degrees(x / (R_EARTH * cos0))
        lat = lat0 + math.degrees(y / R_EARTH)
        return lat, lon

    def ll_to_m(self, lat: float, lon: float) -> Tuple[float, float]:
        lat0, lon0 = self.anchor
        cos0 = math.cos(math.radians(lat0))
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return x, y

    # ---- serialization -----------------------------------------------
    def to_osm(self, path: str) -> None:
        """Emit to a JOSM-readable OSM file with shared node IDs.

        Vertices within ``SHARED_VERTEX_TOL_M`` are assigned the same
        node id, matching the target-OSM convention.  Polygons with
        interior rings (holes, typical for junction polygons that
        wrap around rects) are emitted as OSM multipolygon
        relations — each ring becomes a closed way, and a relation
        ties them together with role=outer / role=inner.
        """
        bucket_size = SHARED_VERTEX_TOL_M
        node_key_to_id: Dict[Tuple[int, int], int] = {}
        node_id_to_ll: Dict[int, Tuple[float, float]] = {}
        next_nid = [-1]

        def _intern(x: float, y: float) -> int:
            kx = int(round(x / bucket_size))
            ky = int(round(y / bucket_size))
            key = (kx, ky)
            if key in node_key_to_id:
                return node_key_to_id[key]
            nid = next_nid[0]
            next_nid[0] -= 1
            node_key_to_id[key] = nid
            # Use the ACTUAL coordinates of the first vertex that
            # landed in this bucket (not the bucket center).  Bucket
            # centers can lie up to bucket_size/2 away from the real
            # vertex, which introduces sub-metre overlaps between
            # adjacent shapes after OSM round-trip.
            node_id_to_ll[nid] = self.m_to_ll(x, y)
            return nid

        def _ring_to_nids(ring_coords, ring_elevs=None):
            """Build a closed-ring nid list from coords.

            Returns ``(nids, elevs_or_None)``.  ``elevs_or_None`` is
            an aligned per-vertex elevation list when ``ring_elevs``
            is provided (used for ``node_altitudes`` polygons);
            otherwise None.  Both share the same dedup logic so the
            element-count invariant survives.

            Defensive against upstream polygon-build bugs:

            * Drops consecutive duplicate nids (two ring vertices
              colliding in the SHARED_VERTEX_TOL_M bucket — would
              produce a zero-length edge that crashes downstream
              meshers).
            * Drops non-consecutive duplicate nids (a ring revisits
              the same node — figure-8 / self-touching polygon —
              keeping only the first occurrence).
            """
            coords = list(ring_coords)
            elevs = list(ring_elevs) if ring_elevs is not None else None
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
                if elevs is not None and len(elevs) > 1 and elevs[0] == elevs[-1]:
                    elevs = elevs[:-1]
            if len(coords) < 3:
                return None, None
            nids = [_intern(x, y) for (x, y) in coords]
            # Dedup any duplicate nid (consecutive OR not).
            seen: set = set()
            deduped_nids: List[int] = []
            deduped_elevs: List[float] = []
            for k, nid in enumerate(nids):
                if nid in seen:
                    continue
                seen.add(nid)
                deduped_nids.append(nid)
                if elevs is not None and k < len(elevs):
                    deduped_elevs.append(elevs[k])
            if len(deduped_nids) < 3:
                return None, None
            deduped_nids.append(deduped_nids[0])
            if elevs is not None:
                deduped_elevs.append(deduped_elevs[0])
                return deduped_nids, deduped_elevs
            return deduped_nids, None

        # Emit one simple way per shape (exterior ring only, with
        # all tags on that way).  Interior rings — which appear
        # on junction polygons that wrap around rect-shaped holes
        # — are dropped for X-Plane patch compatibility: the
        # Ortho4XP patch parser ([O4_Vector_Map.include_patches])
        # iterates ways only, so tags on an OSM multipolygon
        # relation never reach the outer way.  A junction ring
        # emitted without its holes will slightly overlap the
        # rects that used to punch those holes — X-Plane
        # triangulator handles the overlap by seed-region
        # processing; the rects' altitude_high/low tags prevail
        # where they cover.
        way_blocks: List[Tuple[int, List[int], Dict[str, str]]] = []
        next_wid = [-10001]
        for s in self.shapes:
            # Validate the polygon's geometry before emission.
            # Upstream pipeline stages (decomposition, seam-point
            # injection, shared-vertex enforcement) can occasionally
            # produce a self-touching ring that's geometrically
            # invalid; X-Plane's mesh builder crashes on these.
            poly = s.polygon
            if poly is None or poly.is_empty:
                continue
            if not poly.is_valid:
                try:
                    repaired = poly.buffer(0)
                    if (repaired.is_empty
                            or repaired.geom_type
                            not in ("Polygon", "MultiPolygon")):
                        continue
                    if repaired.geom_type == "MultiPolygon":
                        repaired = max(repaired.geoms,
                                       key=lambda g: g.area)
                    if (repaired.is_empty
                            or repaired.geom_type != "Polygon"):
                        continue
                    # node_altitudes from the original ring no longer
                    # aligns with the repaired ring; degrade to a
                    # flat polygon at the mean of the original
                    # vertex elevations to preserve emission.
                    if s.node_altitudes:
                        valid_elevs = [
                            e for e in s.node_altitudes[:-1]]
                        if valid_elevs:
                            s.altitude = round(
                                sum(valid_elevs) / len(valid_elevs),
                                1)
                        s.node_altitudes = None
                    poly = repaired
                except Exception:
                    continue
            # Pass node_altitudes alongside ring coords so dedup of
            # duplicate nids drops the matching elevations too,
            # keeping the per-vertex count invariant.
            ext_nids, ext_elevs = _ring_to_nids(
                poly.exterior.coords,
                s.node_altitudes)
            if ext_nids is None:
                continue
            # Final validity check: rebuild the polygon from the
            # POST-DEDUP lat/lon coords AT THE PRECISION THE OSM
            # FILE WILL CONTAIN (.11f, ≈ 1 mm at the equator).
            # Polygons that are valid at full float precision can
            # become spike-vertex-on-non-adjacent-edge invalid
            # after this truncation; X-Plane's mesh builder
            # crashes on those.  Drop the whole shape rather than
            # ship a polygon X-Plane can't handle.
            try:
                latlon_ring = []
                for nid in ext_nids[:-1]:
                    lat, lon = node_id_to_ll[nid]
                    latlon_ring.append(
                        (float(f"{lat:.11f}"),
                         float(f"{lon:.11f}")))
                check_poly = Polygon(
                    [(lon, lat) for lat, lon in latlon_ring])
                if not check_poly.is_valid:
                    try:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] WARN: dropping "
                            f"invalid polygon (role={s.role}, "
                            f"nids={len(ext_nids) - 1}): "
                            f"X-Plane mesh builder would crash.\n")
                    except Exception:
                        pass
                    continue
                # Sliver-corner safety net: if any interior angle is
                # below SLIVER_ANGLE_THRESHOLD_DEG, drop the polygon.
                # The source-fix _drop_sliver_corners normally
                # eliminates these for junctions, but this catches
                # any path (rect emission, decomposition fragment,
                # repaired-by-buffer(0) polygon) that could
                # reintroduce a needle tip Triangle4XP can't handle.
                ring_m = [self.ll_to_m(lat, lon)
                          for (lat, lon) in latlon_ring]
                cos_thresh = math.cos(
                    math.radians(SLIVER_ANGLE_THRESHOLD_DEG))
                m = len(ring_m)
                worst_ang = None
                for vi in range(m):
                    ax, ay = ring_m[(vi - 1) % m]
                    bx, by = ring_m[vi]
                    cx, cy = ring_m[(vi + 1) % m]
                    v1x, v1y = ax - bx, ay - by
                    v2x, v2y = cx - bx, cy - by
                    n1 = math.hypot(v1x, v1y)
                    n2 = math.hypot(v2x, v2y)
                    if n1 < 1e-9 or n2 < 1e-9:
                        continue
                    cos = (v1x * v2x + v1y * v2y) / (n1 * n2)
                    if cos > cos_thresh:
                        worst_ang = math.degrees(
                            math.acos(max(-1.0, min(1.0, cos))))
                        break
                if worst_ang is not None:
                    try:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] WARN: dropping "
                            f"sliver-corner polygon (role={s.role}, "
                            f"nids={len(ext_nids) - 1}, "
                            f"min angle {worst_ang:.2f}°): "
                            f"X-Plane mesh builder would crash.\n")
                    except Exception:
                        pass
                    continue
            except Exception:
                continue
            tags = {
                "aeroway": AEROWAY_FOR_ROLE.get(s.role, "taxiway"),
                "role": s.role,
            }
            if s.ref:
                tags["ref"] = s.ref
            # Phase-2 elevation tags.  Sloped rects also carry
            # cell_size + profile so X-Plane uses spline
            # interpolation between the high and low short
            # edges, matching the legacy auto-patch format.
            if s.altitude_high is not None and s.altitude_low is not None:
                tags["altitude_high"] = f"{s.altitude_high:.1f}"
                tags["altitude_low"] = f"{s.altitude_low:.1f}"
                tags["cell_size"] = "2"
                tags["profile"] = "spline"
            elif (ext_elevs is not None
                  and len(ext_elevs) == len(ext_nids)):
                # Per-vertex elevation: comma-separated list, one
                # value per ring nid (including the closing repeat).
                # X-Plane mesh builder triangulates the polygon and
                # interpolates linearly between vertex elevations.
                tags["node_altitudes"] = ",".join(
                    f"{e:.1f}" for e in ext_elevs)
            elif s.altitude is not None:
                tags["altitude"] = f"{s.altitude:.1f}"
            way_blocks.append((next_wid[0], ext_nids, tags))
            next_wid[0] -= 1
        rel_blocks: List[Tuple[int, List[Tuple[int, str]],
                               Dict[str, str]]] = []

        lines = [
            "<?xml version='1.0' encoding='UTF-8'?>",
            "<osm version='0.6' upload='false' generator='O4_Airport_Pavement_Builder'>",
        ]
        for nid, (lat, lon) in sorted(node_id_to_ll.items(), reverse=True):
            lines.append(
                f"  <node id='{nid}' action='modify' visible='true' "
                f"lat='{lat:.11f}' lon='{lon:.11f}' />"
            )
        for wid, nids, tags in way_blocks:
            lines.append(
                f"  <way id='{wid}' action='modify' visible='true'>"
            )
            for nid in nids:
                lines.append(f"    <nd ref='{nid}' />")
            for k, v in sorted(tags.items()):
                lines.append(f"    <tag k='{k}' v='{v}' />")
            lines.append("  </way>")
        for rid, members, tags in rel_blocks:
            lines.append(
                f"  <relation id='{rid}' action='modify' visible='true'>"
            )
            for mwid, role in members:
                lines.append(
                    f"    <member type='way' ref='{mwid}' role='{role}' />"
                )
            for k, v in sorted(tags.items()):
                lines.append(f"    <tag k='{k}' v='{v}' />")
            lines.append("  </relation>")
        lines.append("</osm>")
        Path(path).write_text("\n".join(lines) + "\n")


# ──────────────────────────────────────────────────────────────────
# Input loaders
# ──────────────────────────────────────────────────────────────────

def _load_osm_tile(path: str) -> Tuple[Dict[str, Tuple[float, float]],
                                       List[Tuple[str, List[str], Dict[str, str]]],
                                       List[Tuple[str, List[str], Dict[str, str]]]]:
    """Parse an Ortho4XP-cached OSM tile (.osm.bz2 or .osm).

    Returns (nodes, ways, relations) where:
      * nodes: {id: (lat, lon)}
      * ways:  [(id, [nd_ref, ...], {tag: val})]
      * relations: [(id, [member_way_ref, ...], {tag: val})]
        (only outer-role way members are included)
    """
    if path.endswith(".bz2"):
        with bz2.open(path, "rt") as f:
            txt = f.read()
    else:
        txt = Path(path).read_text()

    node_re = re.compile(
        r"""<node\s+id=["'](-?\d+)["'][^>]*?lat=["']([^"']+)["']\s+lon=["']([^"']+)["']"""
    )
    way_re = re.compile(r"""<way[^>]*?id=["'](-?\d+)["'][^>]*>(.*?)</way>""", re.S)
    rel_re = re.compile(
        r"""<relation[^>]*?id=["'](-?\d+)["'][^>]*>(.*?)</relation>""", re.S)
    nd_re = re.compile(r"""<nd\s+ref=["'](-?\d+)["']""")
    tag_re = re.compile(r"""<tag\s+k=["']([^"']+)["']\s+v=["']([^"']+)["']""")
    outer_member_re = re.compile(
        r"""<member\s+type=["']way["']\s+ref=["'](-?\d+)["']\s+role=["']outer["']""")

    nodes: Dict[str, Tuple[float, float]] = {}
    for m in node_re.finditer(txt):
        try:
            nodes[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        except ValueError:
            continue

    ways = []
    for m in way_re.finditer(txt):
        wid = m.group(1)
        body = m.group(2)
        nds = nd_re.findall(body)
        tags = dict(tag_re.findall(body))
        ways.append((wid, nds, tags))

    relations = []
    for m in rel_re.finditer(txt):
        rid = m.group(1)
        body = m.group(2)
        outer = outer_member_re.findall(body)
        tags = dict(tag_re.findall(body))
        relations.append((rid, outer, tags))
    return nodes, ways, relations


def _load_osm_airports(xplane_root: str, icao: str,
                       apt_lat: float, apt_lon: float,
                       radius_deg: float = 0.05
                       ) -> Tuple[Dict[str, Tuple[float, float]],
                                  List[Tuple[str, List[str], Dict[str, str]]],
                                  List[Tuple[str, List[str], Dict[str, str]]]]:
    """Load the airports-layer OSM cache covering the given lat/lon.

    Returns nodes + ways filtered to a bbox around the airport.
    """
    def _fmt(v, pad):
        sign = "+" if v >= 0 else "-"
        return f"{sign}{abs(v):0{pad}d}"

    def _tile_path(lat_tile: int, lon_tile: int) -> str:
        # Ortho4XP's tile layout: e.g. -20-080/-13-078
        lat_group = (lat_tile // 10) * 10
        lon_group = (lon_tile // 10) * 10
        return os.path.join(
            "OSM_data",
            f"{_fmt(lat_group, 2)}{_fmt(lon_group, 3)}",
            f"{_fmt(lat_tile, 2)}{_fmt(lon_tile, 3)}",
            f"{_fmt(lat_tile, 2)}{_fmt(lon_tile, 3)}"
            "_airports.osm.bz2",
        )

    # An airport near a tile boundary may be cached in an adjacent
    # tile — try the natural tile + all 8 neighbours and merge.
    base_lat = int(math.floor(apt_lat))
    base_lon = int(math.floor(apt_lon))
    # If the natural tile's airport-OSM cache is missing, download it
    # via the same Overpass query Ortho4XP's main pipeline uses.  This
    # makes build_airport_pavement self-sufficient when called outside
    # the full Ortho4XP run (e.g. for standalone testing of new
    # airports like CYXY whose OSM tile hasn't been pre-cached).
    natural_path = _tile_path(base_lat, base_lon)
    if not os.path.isfile(natural_path):
        try:
            import O4_OSM_Utils as _OSM
            os.makedirs(os.path.dirname(natural_path), exist_ok=True)
            layer = _OSM.OSM_layer()
            queries = [('node["aeroway"]', 'way["aeroway"]',
                        'rel["aeroway"]')]
            ok = _OSM.OSM_queries_to_OSM_layer(
                queries, layer, base_lat, base_lon,
                tags_of_interest=["all"],
                cached_suffix="airports")
            if not ok:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] WARN: Overpass download failed "
                    f"for airports tile +{base_lat}{base_lon:+04d}; "
                    f"junctions/rects will be empty.\n")
        except Exception as exc:
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] WARN: airport OSM download error: "
                f"{exc}\n")
    # Per user 2026-04-29: OSM Overpass exports use locally-
    # generated NEGATIVE IDs that aren't globally unique.  At
    # HECA (and any airport with adjacent OSM tiles), tile A's
    # node ``-1`` and tile B's node ``-1`` typically refer to
    # entirely different geographic coordinates.  The previous
    # ``nodes.update(n2)`` merge across the 9-tile grid silently
    # overwrites tile A's coordinates with tile B's whenever
    # there's an ID collision — leaving every way that
    # references the colliding ID with the WRONG coordinates.
    #
    # At HECA: 2762 node IDs collide between +30+030 and +30+031
    # (every node ID in +30+030 also appears in +30+031, with
    # different coordinates).  The post-merge HECA centerlines
    # spanned 154 km and the 14 explicit ``aeroway=terminal``
    # ways all had their centroids dragged outside the airport
    # bbox by the overwritten coordinates.
    #
    # Fix: namespace each tile's node IDs by prefixing with the
    # tile coordinate.  Way node-id references are rewritten to
    # use the same prefix at load time, so each way only ever
    # resolves to its own tile's nodes.  After this no ID
    # collision is possible — and the previous "centroid in
    # bbox" filter alone is sufficient (no need for the cross-
    # tile-span guard since corrupted ways no longer exist).
    nodes: Dict[str, Tuple[float, float]] = {}
    ways: List[Tuple[str, List[str], Dict[str, str]]] = []
    relations: List[Tuple[str, List[str], Dict[str, str]]] = []
    seen_paths = set()
    for dlat in (0, -1, 1):
        for dlon in (0, -1, 1):
            tile_lat_n = base_lat + dlat
            tile_lon_n = base_lon + dlon
            osm_path = _tile_path(tile_lat_n, tile_lon_n)
            if osm_path in seen_paths or not os.path.isfile(osm_path):
                continue
            seen_paths.add(osm_path)
            n2, w2, r2 = _load_osm_tile(osm_path)
            # Namespace this tile's node IDs.
            tile_prefix = f"t{tile_lat_n:+03d}{tile_lon_n:+04d}:"
            for nid, coord in n2.items():
                nodes[tile_prefix + nid] = coord
            for wid, nds, tags in w2:
                ways.append(
                    (tile_prefix + wid,
                     [tile_prefix + n for n in nds],
                     tags))
            for rid, outer_ids, tags in r2:
                relations.append(
                    (tile_prefix + rid,
                     [tile_prefix + w for w in outer_ids],
                     tags))
    if not nodes:
        return {}, [], []

    # Filter ways: keep only those whose node centroid lies
    # within ``radius_deg`` of the airport AND whose own bbox
    # spans no more than ``MAX_WAY_SPAN_DEG``.  With per-tile
    # node namespacing above, ID collisions can no longer
    # pollute a way's resolved coordinates — but the span check
    # is cheap and provides a defensive guard against any other
    # data quality issues (e.g. cross-tile ways that genuinely
    # span more than an airport's worth of geography).
    MAX_WAY_SPAN_DEG = 0.1  # ~11 km; airport ways are at most a few km.

    def _in_box(lat, lon):
        return (abs(lat - apt_lat) <= radius_deg and
                abs(lon - apt_lon) <= radius_deg)

    def _way_passes_filters(nds):
        pts = [nodes[n] for n in nds if n in nodes]
        if not pts:
            return False
        lats = [p[0] for p in pts]
        lons = [p[1] for p in pts]
        clat = sum(lats) / len(pts)
        clon = sum(lons) / len(pts)
        if not _in_box(clat, clon):
            return False
        if (max(lats) - min(lats) > MAX_WAY_SPAN_DEG
                or max(lons) - min(lons) > MAX_WAY_SPAN_DEG):
            return False
        return True

    kept_ways = []
    way_by_id: Dict[str, Tuple[str, List[str], Dict[str, str]]] = {}
    for wid, nds, tags in ways:
        way_by_id[wid] = (wid, nds, tags)
        if _way_passes_filters(nds):
            kept_ways.append((wid, nds, tags))
    # Relations: keep if ANY member way passes the filters.
    kept_rels = []
    for rid, outer_ids, tags in relations:
        for wid in outer_ids:
            if wid not in way_by_id:
                continue
            _, nds, _ = way_by_id[wid]
            if _way_passes_filters(nds):
                kept_rels.append((rid, outer_ids, tags))
                break
    return nodes, kept_ways, kept_rels


# ──────────────────────────────────────────────────────────────────
# Meter-space projection
# ──────────────────────────────────────────────────────────────────

def _projection(anchor: Tuple[float, float]):
    lat0, lon0 = anchor
    cos0 = math.cos(math.radians(lat0))

    def to_m(lon: float, lat: float, z=None):
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return (x, y) if z is None else (x, y, z)

    return to_m


def _airport_anchor(apt: APR.Airport) -> Tuple[float, float]:
    if apt.runways:
        r = apt.runways[0]
        return ((r.lat_a + r.lat_b) / 2.0,
                (r.lon_a + r.lon_b) / 2.0)
    if apt.boundary:
        c = apt.boundary.centroid
        return (c.y, c.x)
    return (0.0, 0.0)


# ──────────────────────────────────────────────────────────────────
# Runway rects
# ──────────────────────────────────────────────────────────────────

def _runway_rect_m(runway, to_m) -> Polygon:
    """4-vertex runway rect spanning end-to-end including blast pads.

    Targets are drawn with blast pads included (matches how the
    runway paint looks in satellite imagery).
    """
    ax, ay = to_m(runway.lon_a, runway.lat_a)
    bx, by = to_m(runway.lon_b, runway.lat_b)
    dx, dy = bx - ax, by - ay
    mag = math.hypot(dx, dy)
    if mag < 1e-9:
        return Polygon()
    ux, uy = dx / mag, dy / mag
    a_extra = runway.blast_a_m or 0.0
    b_extra = runway.blast_b_m or 0.0
    ax2 = ax - ux * a_extra
    ay2 = ay - uy * a_extra
    bx2 = bx + ux * b_extra
    by2 = by + uy * b_extra
    px, py = -uy, ux
    half = runway.width_m / 2.0
    return Polygon([
        (ax2 + px * half, ay2 + py * half),
        (bx2 + px * half, by2 + py * half),
        (bx2 - px * half, by2 - py * half),
        (ax2 - px * half, ay2 - py * half),
    ])


def _sample_runway_segment_elev(
        shape: "BuiltShape",
        x: float,
        y: float) -> "Optional[float]":
    """Linearly interpolate a runway segment's elevation at point
    ``(x, y)`` from its ``altitude_high`` / ``altitude_low``
    (sloped 4-corner rect) or ``altitude`` (flat).

    Returns None when the segment carries no elevation info.

    Convention used by :func:`_compute_elevations` segment emit
    (line ~2742): corners ``[0, 3]`` are the HIGH-elevation short
    edge, corners ``[1, 2]`` are the LOW-elevation short edge.
    Point's projection onto the high → low axis gives the
    interpolation parameter ``t`` (clipped to [0, 1] outside the
    segment's extent).
    """
    if shape.altitude is not None:
        return float(shape.altitude)
    if (shape.altitude_high is None
            or shape.altitude_low is None
            or shape.polygon is None):
        return None
    try:
        coords = list(shape.polygon.exterior.coords)
    except Exception:
        return None
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) < 4:
        return 0.5 * (float(shape.altitude_high)
                      + float(shape.altitude_low))
    high_mid_x = 0.5 * (coords[0][0] + coords[3][0])
    high_mid_y = 0.5 * (coords[0][1] + coords[3][1])
    low_mid_x = 0.5 * (coords[1][0] + coords[2][0])
    low_mid_y = 0.5 * (coords[1][1] + coords[2][1])
    ax = low_mid_x - high_mid_x
    ay = low_mid_y - high_mid_y
    L2 = ax * ax + ay * ay
    if L2 < 1e-6:
        return 0.5 * (float(shape.altitude_high)
                      + float(shape.altitude_low))
    t = ((x - high_mid_x) * ax + (y - high_mid_y) * ay) / L2
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    return (float(shape.altitude_high) + t
            * (float(shape.altitude_low)
               - float(shape.altitude_high)))


def _snap_polygon_vertices_to_rect_corners(
        poly: "Polygon",
        sloping_rect_polys: "List[Polygon]",
        snap_tol_m: float = 5.0,
        ) -> "Polygon":
    """Snap every vertex of ``poly`` that lies within ``snap_tol_m``
    of any sloping-rect corner to that corner.

    Per user 2026-04-28 invariant: a sloping rect (runway, primary/
    secondary parallel, stub, cross-connector) can only share a
    *corner* with an adjacent junction polygon — never a point on
    one of its four edges.  Edge-interior coincidence breaks the
    rect's straight-line slope by injecting an extra elevation
    constraint at a non-corner location.

    The runway-crossing-junction emit (``_resolve_runway_crossings``)
    builds a junction polygon from the union of crossing runway
    segments.  Shapely's ``unary_union`` produces vertices at every
    boundary intersection point; some of those points land 2-5 m
    *along* a surviving rect's long edge instead of *at* the rect's
    corner because the dropped (in-crossing) and surviving (out-of-
    crossing) runway segments don't perfectly tile.  Snapping near-
    corner vertices fixes the immediate violation without distorting
    the union polygon's overall footprint.

    Consecutive duplicate vertices produced by the snap are deduped.
    Returns the input polygon unchanged if snapping would leave
    fewer than 3 distinct vertices or produce an invalid polygon.
    """
    try:
        coords = list(poly.exterior.coords)
    except Exception:
        return poly
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) < 3:
        return poly

    corners: List[Tuple[float, float]] = []
    for r in sloping_rect_polys:
        if r is None or r.is_empty:
            continue
        try:
            rc = list(r.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        corners.extend((float(x), float(y)) for x, y in rc)
    if not corners:
        return poly

    snap_tol2 = snap_tol_m * snap_tol_m
    snapped: List[Tuple[float, float]] = []
    for vx, vy in coords:
        best_corner: Optional[Tuple[float, float]] = None
        best_d2 = snap_tol2
        for cx, cy in corners:
            d2 = (vx - cx) ** 2 + (vy - cy) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_corner = (cx, cy)
        if best_corner is not None:
            snapped.append(best_corner)
        else:
            snapped.append((float(vx), float(vy)))

    # Dedupe consecutive duplicates (within 1 cm).
    deduped: List[Tuple[float, float]] = []
    for v in snapped:
        if (not deduped
                or (v[0] - deduped[-1][0]) ** 2
                + (v[1] - deduped[-1][1]) ** 2 > 1e-4):
            deduped.append(v)
    if (len(deduped) > 1
            and (deduped[0][0] - deduped[-1][0]) ** 2
            + (deduped[0][1] - deduped[-1][1]) ** 2 < 1e-4):
        deduped.pop()
    if len(deduped) < 3:
        return poly

    try:
        new_poly = Polygon(deduped + [deduped[0]])
        if not new_poly.is_valid:
            new_poly = new_poly.buffer(0)
        if (new_poly.is_empty
                or new_poly.geom_type != "Polygon"):
            return poly
        return new_poly
    except Exception:
        return poly


def _resolve_runway_crossings(
        layout: "PavementLayout",
        min_overlap_m2: float = 20.0,
        proximity_buffer_m: float = 2.0,
        ) -> int:
    """When two runway segments cross significantly, replace BOTH
    with a single junction polygon covering their union, with
    per-vertex altitudes interpolated from the source segments.

    Per user 2026-04-27: the existing overlap-clip pass clips one
    runway segment against another when they cross, producing a
    non-rectangular shape (5+ vertices) that still carries
    ``altitude_high`` / ``altitude_low`` tags — but X-Plane's patch
    format only renders 4-corner rects with H/L tagging, so the
    extra vertex breaks the slope rendering.  Crossings need
    multi-directional sloping which only a junction polygon
    (with per-vertex altitudes and triangulated rendering) can
    express.

    Detection uses an STRtree + union-find so transitively-
    overlapping segment groups are resolved together (e.g. CYXY's
    crosswind 02/20 crosses BOTH 14R/32L and 14L/32R; all three
    segment groups merge into one junction at the triple crossing).

    Returns the number of crossing groups resolved.
    """
    rwy_indices = [i for i, s in enumerate(layout.shapes)
                    if s.role == ROLE_RUNWAY
                    and s.polygon is not None
                    and not s.polygon.is_empty]
    if len(rwy_indices) < 2:
        return 0
    rwy_polys = [layout.shapes[i].polygon for i in rwy_indices]
    from shapely.strtree import STRtree
    try:
        tree = STRtree(rwy_polys)
    except Exception:
        return 0

    # Union-find for transitive grouping.
    n = len(rwy_indices)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union_uf(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for ai in range(n):
        pa = rwy_polys[ai]
        try:
            cands = tree.query(pa)
        except Exception:
            continue
        for ci in cands:
            bi = int(ci)
            if bi <= ai:
                continue
            pb = rwy_polys[bi]
            try:
                inter = pa.intersection(pb)
                if inter.is_empty or inter.area < min_overlap_m2:
                    continue
                union_uf(ai, bi)
            except Exception:
                continue

    # Group by root; ignore singleton groups.
    groups: Dict[int, List[int]] = {}
    for i in range(n):
        r = find(i)
        groups.setdefault(r, []).append(i)

    # Pre-compute sloping-rect corner snap targets for the corner-
    # alignment pass below.  Sloping rects = runway + parallel +
    # stub + cross-connector; junction vertices can only land on
    # their corners, never on their edges.
    sloping_roles_for_snap = (ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                              ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                              ROLE_CROSS_CONNECTOR)
    all_sloping_indices = [i for i, s in enumerate(layout.shapes)
                            if s.role in sloping_roles_for_snap
                            and s.polygon is not None
                            and not s.polygon.is_empty]

    drop_set: set = set()
    new_shapes: List[BuiltShape] = []
    n_resolved = 0
    for members in groups.values():
        if len(members) <= 1:
            continue
        seg_shapes = [layout.shapes[rwy_indices[m]] for m in members]
        seg_polys = [s.polygon for s in seg_shapes]
        member_set = {rwy_indices[m] for m in members}
        try:
            union_poly = unary_union(seg_polys)
            if not union_poly.is_valid:
                union_poly = union_poly.buffer(0)
        except Exception:
            continue
        if union_poly.is_empty:
            continue
        if union_poly.geom_type != "Polygon":
            polys = [g for g in getattr(union_poly, "geoms", [])
                      if g.geom_type == "Polygon"]
            polys.sort(key=lambda p: -p.area)
            if not polys:
                continue
            union_poly = polys[0]
        # Snap any union-polygon vertex that lands within 5 m of a
        # surviving sloping-rect corner to that corner.  Without
        # this, ``unary_union``'s intersection points can sit a few
        # metres along a surviving rect's edge — violating the
        # corner-only-junction-vertex invariant.
        other_sloping_polys = [
            layout.shapes[i].polygon
            for i in all_sloping_indices
            if i not in member_set]
        union_poly = _snap_polygon_vertices_to_rect_corners(
            union_poly, other_sloping_polys, snap_tol_m=5.0)
        try:
            coords = list(union_poly.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) < 3:
            continue

        # Per-vertex altitudes: average over segments whose buffered
        # polygon contains the vertex (i.e. the segments physically
        # adjacent to that corner).  Falls back to all-segment
        # average if no segment contains the vertex (defensive — a
        # corner introduced by Shapely's union may sit ε outside
        # every input).
        seg_buffers = []
        for s in seg_shapes:
            try:
                seg_buffers.append(
                    s.polygon.buffer(proximity_buffer_m))
            except Exception:
                seg_buffers.append(None)
        ring_alts: List[Optional[float]] = []
        for (cx, cy) in coords:
            pt = Point(cx, cy)
            samples: List[float] = []
            for buf, s in zip(seg_buffers, seg_shapes):
                if buf is None:
                    continue
                try:
                    if buf.contains(pt):
                        e = _sample_runway_segment_elev(s, cx, cy)
                        if e is not None:
                            samples.append(e)
                except Exception:
                    continue
            if not samples:
                # Fallback: all segments.
                for s in seg_shapes:
                    e = _sample_runway_segment_elev(s, cx, cy)
                    if e is not None:
                        samples.append(e)
            if samples:
                ring_alts.append(round(
                    sum(samples) / len(samples), 1))
            else:
                ring_alts.append(None)
        if any(a is None for a in ring_alts):
            continue
        # node_altitudes spans the closed ring.
        closed_alts: List[float] = list(ring_alts) + [ring_alts[0]]
        ref_combined = "+".join(
            s.ref for s in seg_shapes if s.ref)
        new_shape = BuiltShape(
            polygon=union_poly,
            role=ROLE_JUNCTION,
            ref=ref_combined,
            node_altitudes=closed_alts)
        new_shapes.append(new_shape)
        for m in members:
            drop_set.add(rwy_indices[m])
        n_resolved += 1

    if drop_set:
        layout.shapes = [s for i, s in enumerate(layout.shapes)
                          if i not in drop_set]
        layout.shapes.extend(new_shapes)
    return n_resolved


def _detect_runway_shoulders(
        runway,
        to_m,
        pav_polys: "List[Polygon]",
        max_lat_gap_m: float = 1.0,
        min_self_inside_frac: float = 0.80,
        min_runway_overlap_frac: float = 0.40,
        min_axial_aspect: float = 2.5,
        min_strip_length_m: float = 50.0,
        ) -> "Tuple[float, float, List[int]]":
    """Scan ``pav_polys`` for polygons that are runway pavement
    (shoulders or the runway's own envelope polygon often labelled
    as a "taxiway" by apt.dat) and return the perpendicular extent
    those polygons + the runway itself jointly occupy.

    Per user 2026-04-27: long thin pavement polygons parallel to a
    runway and touching it (within 1 m of the runway edge) are
    SHOULDERS — fold them into the runway emit so the runway
    pavement reflects actual paved area, not just the apt.dat row-100
    designation.  At HECA, runway 05R/23L's apt.dat row-100 width is
    60 m but a row-110 polygon "New Taxiway 1" sits centered on the
    runway at 75.6 m wide (the runway-with-shoulders envelope).
    Without this absorption, the 7.8 m-per-side shoulder strip
    emits as junctions wrapping the runway.  At CYXY's short
    crosswind 02/20 (22.9 m wide), narrow strips "North of 02" and
    "South of 02" act as shoulders extending the runway pavement
    asymmetrically.

    Returns ``(new_left, new_right, absorbed_indices)``:
      * ``new_left``  — most-negative perpendicular offset from the
        runway centerline that emitted runway pavement should reach.
      * ``new_right`` — most-positive perpendicular offset.
      * ``absorbed_indices`` — indices into ``pav_polys`` of the
        polygons folded into the runway.  Caller should remove them
        from ``pav_polys`` (and ``apt_only_pav_polys`` if applicable)
        so they don't re-emit as separate junction shapes.

    When no shoulders are found, returns ``(-runway_half, +runway_half,
    [])``.

    Detection rules per polygon:
      1. Long axis aligned with runway (axial extent ≥
         ``min_axial_aspect`` × perpendicular extent, axial extent ≥
         ``min_strip_length_m``).
      2. Polygon mostly inside the runway's longitudinal extent:
         ``polygon's u-extent inside [0, L] ≥ min_self_inside_frac
         × polygon's total u-extent``.
      3. Polygon covers a significant fraction of the runway's
         length: ``u-overlap with [0, L] ≥ min_runway_overlap_frac
         × runway length L``.  Filters out tiny end-only blast pads
         (e.g. CYXY's '32L' polygon, 79 m at one end of a 2899 m
         runway) while admitting partial-length shoulders that
         cover roughly half the runway (e.g. CYXY's "North of 02"
         covering 47 % of the 02/20 runway).
      4. Perpendicular interval overlaps or touches
         ``[-runway_half, +runway_half]`` within ``max_lat_gap_m``.
      5. Side-extension limit: the polygon's extension PAST the
         runway edge on each side is less than the runway width.
         Filters out large aprons (e.g. CYXY's "Apron 1 and E",
         329 m wide, that nibbles the runway edge by 4 m) while
         admitting both wide envelopes (HECA's runway-with-
         shoulders polygon, 7.8 m past each edge) and narrow
         single-side shoulders (CYXY's "North of 02", 18 m past
         the north edge).
    """
    ax, ay = to_m(runway.lon_a, runway.lat_a)
    bx, by = to_m(runway.lon_b, runway.lat_b)
    udx = bx - ax
    udy = by - ay
    L = math.hypot(udx, udy)
    runway_half = runway.width_m / 2.0
    if L < 1.0:
        return (-runway_half, runway_half, [])
    ux, uy = udx / L, udy / L
    nx, ny = -uy, ux  # perpendicular (rotated 90° CCW from u)

    runway_width = 2.0 * runway_half
    new_left = -runway_half
    new_right = runway_half
    absorbed: List[int] = []
    for idx, pav in enumerate(pav_polys):
        if pav is None or pav.is_empty:
            continue
        try:
            coords = list(pav.exterior.coords)
        except Exception:
            continue
        if not coords:
            continue
        u_proj = [(cx - ax) * ux + (cy - ay) * uy for cx, cy in coords]
        n_proj = [(cx - ax) * nx + (cy - ay) * ny for cx, cy in coords]
        u_min, u_max = min(u_proj), max(u_proj)
        n_min, n_max = min(n_proj), max(n_proj)
        u_extent = u_max - u_min
        n_extent = n_max - n_min
        if u_extent < min_strip_length_m:
            continue
        if u_extent < min_axial_aspect * n_extent:
            continue
        # Polygon must be mostly inside the runway's longitudinal
        # extent.
        u_overlap = min(u_max, L) - max(u_min, 0.0)
        if u_overlap < min_self_inside_frac * u_extent:
            continue
        # Polygon must cover a significant fraction of the runway
        # length (filters out end-only blast pads / stopways).
        if u_overlap < min_runway_overlap_frac * L:
            continue
        # Perpendicular adjacency / overlap.
        if n_max < -runway_half - max_lat_gap_m:
            continue
        if n_min > runway_half + max_lat_gap_m:
            continue
        # Side-extension limit: the polygon's extension PAST the
        # runway edge on each side must be less than the runway
        # width.  Catches both:
        #   - HECA's "New Taxiway 1" envelope (8.7 m past each edge,
        #     well under 60 m runway width — accepted).
        #   - CYXY's "North of 02" shoulder (18 m past the north
        #     edge, under 22.9 m runway width — accepted).
        # Rejects:
        #   - CYXY's "Apron 1 and E" (extends 325 m past the east
        #     edge of 14R/32L, way over the 45.7 m runway width).
        north_ext = max(0.0, n_max - runway_half)
        south_ext = max(0.0, -runway_half - n_min)
        if north_ext > runway_width or south_ext > runway_width:
            continue
        absorbed.append(idx)
        if n_min < new_left:
            new_left = n_min
        if n_max > new_right:
            new_right = n_max
    return (new_left, new_right, absorbed)


# ──────────────────────────────────────────────────────────────────
# Pavement union helpers
# ──────────────────────────────────────────────────────────────────

# Tolerance for bridging numerical / sub-meter gaps between near-
# touching apt.dat polygons.  apt.dat at busy airports often stores
# adjacent apron areas as separate row-110 polygons whose shared
# edges have sub-millimeter-to-cm floating-point differences (e.g.
# at SPJC's SE apron a 6.5 m thin "strip" appears between two large
# apron polygons because their shared boundary y-values differ by
# 5 cm).  unary_union doesn't merge these because they don't
# overlap — but for our purposes they ARE one continuous coverage.
# Closing 0.5 m gaps via buffer-shrink merges them while preserving
# real holes (typically meters-wide non-pavement islands).
PAVEMENT_BRIDGE_GAP_M = 0.1


def _merge_near_touching(geom: Optional[Polygon],
                         eps: float = PAVEMENT_BRIDGE_GAP_M
                         ) -> Optional[Polygon]:
    """Force-merge near-touching components of ``geom`` (a possibly
    MultiPolygon) by a buffer-then-shrink.  Returns the same kind
    of geometry (Polygon if single, MultiPolygon if truly disjoint).

    Per user 2026-04-27: we want a single pavement union with holes
    only — apt.dat polygon boundaries shouldn't survive into the
    output.  This helper bridges sub-meter precision gaps between
    apt.dat polygons that ``unary_union`` leaves disjoint.
    """
    if geom is None or geom.is_empty:
        return geom
    try:
        merged = geom.buffer(eps, join_style=2).buffer(
            -eps, join_style=2)
        if merged.is_empty:
            return geom
        if merged.geom_type not in ("Polygon", "MultiPolygon"):
            return geom
        return merged
    except Exception:
        return geom


def _clip_residue_at_stub_long_edges(
        residue: "Polygon",
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        outer_buffer_m: float = 30.0,
        ) -> "Polygon":
    """Subtract a thin strip just OUTSIDE each STUB rect's long
    edges from the residue.  Per user 2026-04-27 invariant: a
    junction polygon must never run along a sloping rect's long
    edge — the rect's long edges are TAXI BOUNDARIES, not junction
    boundaries.  When apt.dat pavement bulges past a stub's long
    edge between its two short corners, the bulge becomes a
    "wrap" on the adjacent junction.  Removing the bulge from the
    residue forces the junction to stop at the stub's short-edge
    corners.

    The strip extends ``outer_buffer_m`` past each long edge —
    enough to swallow typical apt.dat curvature noise (the SPJC F
    south-edge bulge is 21 m).  Stubs follow
    ``_rect_from_axis_extended``'s corner convention: corners
    [0,1] form long-edge "side1", [2,3] form long-edge "side2".
    """
    if residue is None or residue.is_empty:
        return residue
    for rect, axis, role, ref in taxi_rects:
        if role != ROLE_STUB:
            continue
        try:
            rc = list(rect.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Long edges per ``_rect_from_axis_extended`` convention.
        long_edges = [(rc[0], rc[1]), (rc[2], rc[3])]
        for (e0, e1) in long_edges:
            ex = e1[0] - e0[0]
            ey = e1[1] - e0[1]
            mag = math.hypot(ex, ey)
            if mag < 0.5:
                continue
            ux, uy = ex / mag, ey / mag
            # Outward normal (away from the rect centroid).
            nx, ny = -uy, ux
            cx_r, cy_r = rect.centroid.x, rect.centroid.y
            mid_x = 0.5 * (e0[0] + e1[0])
            mid_y = 0.5 * (e0[1] + e1[1])
            if (cx_r - mid_x) * nx + (cy_r - mid_y) * ny > 0:
                nx, ny = -nx, -ny
            # Build a thin rectangle outside the long edge.
            o0 = (e0[0] + nx * outer_buffer_m,
                  e0[1] + ny * outer_buffer_m)
            o1 = (e1[0] + nx * outer_buffer_m,
                  e1[1] + ny * outer_buffer_m)
            try:
                strip = Polygon([e0, e1, o1, o0])
                if strip.is_valid and not strip.is_empty:
                    residue = residue.difference(strip)
            except Exception:
                continue
    return residue


def _drop_primary_parallels_embedded_in_pavement(
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        apt_pav_union: "Optional[Polygon]",
        runway_polys: "Optional[List[Polygon]]" = None,
        adjacency_frac: float = 0.10,
        proximity_m: float = 1.0,
        embed_frac: float = 0.95,         # legacy param, unused
        long_edge_buffer_m: float = 5.0,  # legacy param, unused
        ) -> "List[Tuple[Polygon, LineString, str, str]]":
    """Drop ``primary_parallel`` rects whose long edges sit entirely
    inside apt.dat row-110 pavement.

    Per user 2026-04-27 invariant: a junction polygon must NEVER run
    along a sloping rect's long edge — X-Plane's elevation engine
    over-constrains the junction's spline to F's long-edge altitudes
    and produces visible elevation glitches.

    There are two ways the invariant gets violated for a primary
    parallel:

    1. apt.dat pavement bulges past the long edge between the rect's
       short corners.  The residue then runs along the long edge to
       reach the bulge.  ``_clip_residue_at_stub_long_edges`` already
       handles this for STUBs by carving the bulge away.
    2. The rect sits FULLY INSIDE a paved area (apron, ramp, big
       terminal area).  apt.dat covers the rect's footprint AND
       everything around it; the residue is the apron, and it
       inevitably traces the rect's long edges as it wraps around
       the rect-shaped hole.

    Case 2 doesn't have a "carve away the bulge" fix — there's no
    bulge, the entire surrounding pavement is legitimate apron.  The
    correct behaviour is to NOT emit the rect at all: let the apron
    junction cover the rect's footprint and slope multi-
    directionally.  The primary-parallel grade rule isn't worth
    enforcing for a taxi lane that's just one of many possible
    paths through a paved apron — the apron's own elevation
    propagation produces a perfectly serviceable surface.

    SPJC's F is the canonical example.  F is a 121×60 m primary
    parallel sitting in the middle of the SE apron; row-110 pavement
    covers F + everything around it.  Dropping F lets the SE apron
    (junction -10132) absorb F's footprint instead of wrapping
    around F's long edges.

    Implementation: a primary-parallel rect is considered embedded
    when AT LEAST ONE of its two long edges is ≥ ``embed_frac``
    inside ``apt_pav_union``, AND the strip just outside that edge
    (``long_edge_buffer_m`` wide) is also covered.  Why "at least
    one" instead of "both": a primary parallel can sit on the edge
    of an apron, with one long edge inside the apron (where the
    junction wraps it — the violation we're fixing) and the other
    long edge along the apron's outer pavement boundary (where
    there's no junction at all, so no wrap to worry about).  SPJC's
    F is exactly this: south long edge fully inside the SE apron
    (100 %), north long edge on the apron's outer boundary (27 %).
    Dropping F is correct in both cases — the apron absorbs the
    embedded side, and the boundary side becomes the apron's new
    edge.  Only when NEITHER long edge is embedded does the rect
    actually carry the pavement (a primary parallel through grass)
    — those we must keep.
    """
    if apt_pav_union is None or apt_pav_union.is_empty:
        return taxi_rects
    # Per user 2026-04-28: a sloping rect cannot have a junction or
    # apron polygon running alongside its long edge.  The slope
    # along the rect's long edge is uniform (altitude_high at one
    # short edge, altitude_low at the other); any adjacent junction-
    # class pavement would have to match that slope along the seam,
    # which produces visible elevation glitches in X-Plane when the
    # junction's natural DEM slope differs from the rect's straight-
    # line slope.  Absorb the rect into the junction so the whole
    # region becomes one polygon with per-vertex node_altitudes
    # capturing the natural slope.
    #
    # "Junction-class pavement" = pav_union − runway_polys − OTHER
    # taxi rects.  This excludes:
    #   - Runway adjacency (parallel taxis legitimately run alongside
    #     a runway; the runway has its own slope already).
    #   - Other taxi-rect adjacency (two adjacent taxi rects each
    #     own their slope along their own axes).
    junction_pav = apt_pav_union
    if runway_polys:
        try:
            for r in runway_polys:
                if r is not None and not r.is_empty:
                    junction_pav = junction_pav.difference(r)
        except Exception:
            junction_pav = apt_pav_union
    # Subtract every taxi-rect's footprint (including the rect being
    # tested — that's CORRECT, because we want to know whether the
    # JUNCTION polygon will be adjacent after the rect is emitted
    # and the residue is computed).  After subtraction, the apron
    # polygon containing the rect has a rect-shaped hole; the hole's
    # boundary IS within 1 m of the rect's long edge.
    try:
        all_taxi_polys = [
            r for (r, _ax, _role, _ref) in taxi_rects
            if r is not None and not r.is_empty]
        if all_taxi_polys:
            taxi_union = unary_union(all_taxi_polys)
            if not taxi_union.is_empty:
                junction_pav = junction_pav.difference(taxi_union)
    except Exception:
        pass
    if junction_pav.is_empty:
        return taxi_rects
    # All sloping-rect roles are subject to the absorption rule.
    sloping_rect_roles = {ROLE_PRIMARY_PARALLEL,
                          ROLE_SECONDARY_PARALLEL,
                          ROLE_STUB, ROLE_CROSS_CONNECTOR}

    # Per user 2026-04-28 (final): absorb wherever EITHER long edge
    # has junction-class pavement running alongside it.  Sloping
    # rects cannot share a long edge with an apron/junction polygon
    # — the rect's straight-line slope along the long edge has to
    # match the junction's natural DEM slope along the seam, which
    # it generally won't.  Whether the apron sits on one side or
    # both, the violation is the same.  Apply PARTIAL absorption —
    # only the axial range where adjacency holds gets absorbed; the
    # rest of the rect survives as shorter rect(s).
    #
    # Probe semantics: at every 5 m axis step, compute a point
    # ``OUTER_PROBE_M`` METRES OUTSIDE each long edge and ask
    # whether that point lies *directly inside* junction-pav (no
    # buffer).  A real apron extends many meters past the rect's
    # long edge, so a point 2 m beyond the edge will hit it.  The
    # 1 m-and-buffer formulation we replaced was sensitive to
    # sub-1 m slivers caused by apt.dat / DSF polygons that render
    # ~0.5 m wider than the OSM-tagged taxi width — those slivers
    # are polygon-imprecision noise, NOT real apron adjacency, and
    # spuriously absorbed every rect at CYXY.
    #
    # Walk each rect's axis at 5 m steps.  At each step, sample the
    # 2 m-outside-left and 2 m-outside-right points.  A step is
    # "adjacent" when EITHER outside point is in junction-pav.  Find
    # contiguous adjacent runs ≥ 10 % of axial length; keep non-
    # adjacent intervals (≥ 30 m fragments) as new rects.
    #
    # Use cases:
    #   - F primary parallel at CYXY: south 30 % is INSIDE the apron
    #     (apron extends past both long edges) → adjacent →
    #     absorbed.  North 70 % extends out of the apron — only
    #     polygon-imprecision slivers within 0.5 m of the long edge,
    #     no substantive apron at 2 m → not adjacent → kept,
    #     extending from apron edge to the centerline bend.
    #   - E primary parallel at CYXY: apron runs alongside ONE long
    #     side for ~65 % of length, with several meters of apron
    #     past the edge → adjacent → absorbed.  Non-adjacent
    #     fragments (if any ≥ 30 m) survive.
    #   - F embedded in SPJC apron (both sides 100 %): fully
    #     absorbed.
    SAMPLE_STEP_M = 5.0
    MIN_KEPT_M = 30.0
    OUTER_PROBE_M = 2.0

    kept: List[Tuple[Polygon, LineString, str, str]] = []
    abs_refs: List[str] = []
    n_full = 0
    n_split = 0
    n_clipped = 0
    for entry in taxi_rects:
        rect, axis, role, ref = entry
        if role not in sloping_rect_roles:
            kept.append(entry)
            continue
        try:
            rc = list(rect.exterior.coords)
        except Exception:
            kept.append(entry)
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            kept.append(entry)
            continue
        # Axis from short-edge midpoints (corners 0,3 → one short
        # edge; 1,2 → other) per ``_rect_from_axis_extended``
        # convention.
        a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                  0.5 * (rc[0][1] + rc[3][1]))
        b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                  0.5 * (rc[1][1] + rc[2][1]))
        ax = b_mid[0] - a_mid[0]
        ay = b_mid[1] - a_mid[1]
        L = math.hypot(ax, ay)
        if L < 30.0:
            kept.append(entry)
            continue
        ux = ax / L
        uy = ay / L
        nx, ny = -uy, ux
        half_w = math.hypot(rc[0][0] - a_mid[0], rc[0][1] - a_mid[1])
        if half_w < 1.0:
            kept.append(entry)
            continue

        n_steps = max(2, int(L / SAMPLE_STEP_M) + 1)
        either_adj = [False] * n_steps
        outer = half_w + OUTER_PROBE_M
        for i in range(n_steps):
            u = min(L, i * SAMPLE_STEP_M)
            cx = a_mid[0] + u * ux
            cy = a_mid[1] + u * uy
            try:
                left_pt = Point(cx + nx * outer,
                                 cy + ny * outer)
                right_pt = Point(cx - nx * outer,
                                  cy - ny * outer)
                either_adj[i] = (
                    bool(junction_pav.contains(left_pt))
                    or bool(junction_pav.contains(right_pt)))
            except Exception:
                continue

        # Find contiguous "either adjacent" runs ≥ 10 % of axis.
        min_run_steps = max(1, int(adjacency_frac * n_steps))
        absorbed_intervals: List[Tuple[float, float]] = []
        i = 0
        while i < n_steps:
            if not either_adj[i]:
                i += 1
                continue
            j = i
            while j < n_steps and either_adj[j]:
                j += 1
            if (j - i) >= min_run_steps:
                u_start = i * SAMPLE_STEP_M
                u_end = min(L, j * SAMPLE_STEP_M)
                absorbed_intervals.append((u_start, u_end))
            i = j

        if not absorbed_intervals:
            kept.append(entry)
            continue

        # Compute kept intervals = [0, L] minus absorbed.
        kept_intervals: List[Tuple[float, float]] = []
        prev_end = 0.0
        for s, e in absorbed_intervals:
            if s > prev_end:
                kept_intervals.append((prev_end, s))
            prev_end = e
        if L > prev_end:
            kept_intervals.append((prev_end, L))
        kept_intervals = [(s, e) for s, e in kept_intervals
                          if (e - s) >= MIN_KEPT_M]

        if not kept_intervals:
            n_full += 1
            abs_refs.append(ref or "?")
            continue

        # No-op: kept covers (almost) the full rect.
        if (len(kept_intervals) == 1
                and kept_intervals[0][0] <= SAMPLE_STEP_M
                and kept_intervals[0][1] >= L - SAMPLE_STEP_M):
            kept.append(entry)
            continue

        # Build new rects from kept intervals.  Drop any kept
        # fragment that's apron-interior — i.e. ≥ 2 of its 4 corners
        # are off the pavement boundary.  These are tiny rects
        # floating inside an apron polygon, with no real corridor
        # geometry around them; emitting them just creates a rect-
        # shaped hole in the apron's residue (and a thin diagonal
        # junction polygon connecting the hole to the apron's
        # outer boundary — see CYXY -10005/-10129 regression).
        # Mirrors ``_build_taxi_rects``'s apron-interior check.
        APRON_INTERIOR_TOL_M = 2.0
        try:
            pav_boundary = apt_pav_union.boundary
        except Exception:
            pav_boundary = None
        n_dropped_interior = 0
        new_rects: List[Tuple[Polygon, LineString, str, str]] = []
        for u_lo, u_hi in kept_intervals:
            new_a_mid = (a_mid[0] + u_lo * ux,
                         a_mid[1] + u_lo * uy)
            new_b_mid = (a_mid[0] + u_hi * ux,
                         a_mid[1] + u_hi * uy)
            new_corners = [
                (new_a_mid[0] + nx * half_w,
                 new_a_mid[1] + ny * half_w),
                (new_b_mid[0] + nx * half_w,
                 new_b_mid[1] + ny * half_w),
                (new_b_mid[0] - nx * half_w,
                 new_b_mid[1] - ny * half_w),
                (new_a_mid[0] - nx * half_w,
                 new_a_mid[1] - ny * half_w),
            ]
            try:
                new_rect = Polygon(new_corners)
                if not new_rect.is_valid:
                    new_rect = new_rect.buffer(0)
                if (new_rect.is_empty
                        or new_rect.geom_type != "Polygon"):
                    continue
                # Apron-interior check on the kept fragment.  A
                # kept fragment naturally has 2 corners at the
                # absorbed/kept split (interior to the original
                # rect's footprint, slightly off-boundary by ~1-3 m
                # due to apt.dat boundary imprecision at corridor
                # narrowings).  The discriminator is "all 4 corners
                # off-boundary" — captures truly apron-floating
                # fragments (CYXY -10005: all 4 corners 2.5-23 m
                # off) without dropping normal split-end fragments
                # (CYXY -10004: 3 corners 0-2.9 m off but at most
                # one corner > 3 m off).
                if pav_boundary is not None:
                    n_off = sum(
                        1 for (cx, cy) in new_corners
                        if Point(cx, cy).distance(pav_boundary)
                        > APRON_INTERIOR_TOL_M)
                    max_off = max(
                        (Point(cx, cy).distance(pav_boundary)
                         for (cx, cy) in new_corners),
                        default=0.0)
                    # Drop if all 4 corners off boundary AND at
                    # least one is > 5 m off (indicating real
                    # apron-floating, not boundary imprecision).
                    if n_off == 4 and max_off > 5.0:
                        n_dropped_interior += 1
                        continue
                new_axis = LineString([new_a_mid, new_b_mid])
                new_rects.append((new_rect, new_axis, role, ref))
            except Exception:
                continue
        kept.extend(new_rects)
        if not new_rects:
            n_full += 1
        elif len(new_rects) >= 2:
            n_split += 1
        else:
            n_clipped += 1
        abs_refs.append(
            f"{ref or '?'}"
            f"{'/int=' + str(n_dropped_interior) if n_dropped_interior else ''}")

    if abs_refs:
        try:
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] long-edge-adjacent absorption: "
                f"{n_full} dropped, {n_split} split, "
                f"{n_clipped} clipped (refs: "
                f"{', '.join(abs_refs)}).\n")
        except Exception:
            pass
    return kept


def _split_primary_parallels_at_pavement_boundary(
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        pav_union: "Optional[Polygon]",
        step_m: float = 10.0,
        min_dropped_m: float = 100.0,
        min_kept_m: float = 50.0,
        ) -> "List[Tuple[Polygon, LineString, str, str]]":
    """Clip primary_parallel rects whose long edge has a contiguous
    embedded sub-range at one end of the rect's axis.

    Per user 2026-04-27: when a taxiway is bounded by apron pavement
    along most of one side, that bounded section should be ABSORBED
    INTO THE APRON (covered by it, removed from the taxi rect).
    The rule fires for FULL embedding via
    ``_drop_primary_parallels_embedded_in_pavement``; this helper
    handles PARTIAL embedding where only a sub-range of the rect's
    axis is bounded.  The unbounded portion stays as a (shorter)
    sloping rect; the absorbed portion's footprint reappears in
    the residue → junction polygon → apron extends to cover it.

    Algorithm:
        1. Walk the rect's axis at ``step_m``-metre steps.  At each
           step ``u``, sample the left and right long-edge midpoints
           at axial offset ``u`` and check whether either falls
           inside ``pav_union``.
        2. Mark each step as "embedded" or not.
        3. Find embedded contiguous prefixes / suffixes ≥
           ``min_dropped_m`` long.
        4. If clipping the embedded prefix and/or suffix leaves a
           kept range ≥ ``min_kept_m``, replace the rect with the
           clipped version.

    The clipped rect uses the same axis direction; its corners are
    at the new axial bounds projected to the original perpendicular
    half-width.

    At CYXY, taxiway E (469 m, NW-SE oriented) has its NW half
    bounded by SW apron pavement on its west side.  This helper
    clips E's NW portion (~290 m) and keeps the SE portion
    (~180 m), leaving the apron polygon to cover E's NW footprint.
    """
    if pav_union is None or pav_union.is_empty:
        return taxi_rects
    out: List[Tuple[Polygon, LineString, str, str]] = []
    n_clipped = 0
    for entry in taxi_rects:
        rect, axis, role, ref = entry
        if role != ROLE_PRIMARY_PARALLEL:
            out.append(entry)
            continue
        try:
            rc = list(rect.exterior.coords)
        except Exception:
            out.append(entry)
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            out.append(entry)
            continue
        # Axis from short-edge-A midpoint (corners 0,3) to short-
        # edge-B midpoint (corners 1,2).  Per the convention used by
        # ``_rect_from_axis_extended`` / runway-segment emit, corners
        # [0, 1] form one long edge and [2, 3] form the other.
        a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                  0.5 * (rc[0][1] + rc[3][1]))
        b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                  0.5 * (rc[1][1] + rc[2][1]))
        ax = b_mid[0] - a_mid[0]
        ay = b_mid[1] - a_mid[1]
        L = math.hypot(ax, ay)
        if L < 1.0:
            out.append(entry)
            continue
        ux, uy = ax / L, ay / L
        nx, ny = -uy, ux
        half_w = math.hypot(rc[0][0] - a_mid[0], rc[0][1] - a_mid[1])
        if half_w < 1.0:
            out.append(entry)
            continue
        n_steps = max(2, int(math.floor(L / step_m)) + 1)
        embedded: List[bool] = []
        try:
            for i in range(n_steps):
                u = min(L, i * step_m)
                cx = a_mid[0] + u * ux
                cy = a_mid[1] + u * uy
                left = Point(cx + nx * half_w, cy + ny * half_w)
                right = Point(cx - nx * half_w, cy - ny * half_w)
                emb = (pav_union.contains(left)
                       or pav_union.contains(right))
                embedded.append(emb)
        except Exception:
            out.append(entry)
            continue
        # Find embedded prefix length (in steps).
        pfx = 0
        while pfx < n_steps and embedded[pfx]:
            pfx += 1
        # Find embedded suffix length (in steps).
        sfx = 0
        while sfx < n_steps - pfx and embedded[n_steps - 1 - sfx]:
            sfx += 1
        # Convert to metres.  step_m is the axis step; the embedded
        # prefix ends at the LAST embedded step's u-coordinate, so
        # the prefix length in metres is ``pfx_steps × step_m``.
        pfx_m = pfx * step_m
        sfx_m = sfx * step_m
        # Allowable clip: prefix or suffix must be ≥ min_dropped_m,
        # and the kept range must be ≥ min_kept_m.
        clip_pfx = pfx_m >= min_dropped_m
        clip_sfx = sfx_m >= min_dropped_m
        if not (clip_pfx or clip_sfx):
            out.append(entry)
            continue
        u_lo = pfx_m if clip_pfx else 0.0
        u_hi = (L - sfx_m) if clip_sfx else L
        if u_hi - u_lo < min_kept_m:
            out.append(entry)
            continue
        # Build clipped rect: 4 corners at u_lo and u_hi axial
        # offsets, ±half_w perpendicular.
        new_a_mid = (a_mid[0] + u_lo * ux, a_mid[1] + u_lo * uy)
        new_b_mid = (a_mid[0] + u_hi * ux, a_mid[1] + u_hi * uy)
        new_corners = [
            (new_a_mid[0] + nx * half_w, new_a_mid[1] + ny * half_w),
            (new_b_mid[0] + nx * half_w, new_b_mid[1] + ny * half_w),
            (new_b_mid[0] - nx * half_w, new_b_mid[1] - ny * half_w),
            (new_a_mid[0] - nx * half_w, new_a_mid[1] - ny * half_w),
        ]
        try:
            new_rect = Polygon(new_corners)
            if not new_rect.is_valid:
                new_rect = new_rect.buffer(0)
            if (new_rect.is_empty
                    or new_rect.geom_type != "Polygon"):
                out.append(entry)
                continue
            new_axis = LineString([new_a_mid, new_b_mid])
            out.append((new_rect, new_axis, role, ref))
            n_clipped += 1
            try:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] clipped primary_parallel "
                    f"{ref!r}: {L:.0f}m → {(u_hi - u_lo):.0f}m "
                    f"(dropped "
                    + ("prefix " if clip_pfx else "")
                    + ("suffix " if clip_sfx else "")
                    + "embedded in pavement).\n")
            except Exception:
                pass
        except Exception:
            out.append(entry)
            continue
    return out


def _add_stub_to_runway_bridges(
        residue: "Polygon",
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        runway_union: Optional["Polygon"],
        max_bridge_m: float = 100.0,
        ) -> "Polygon":
    """For each STUB rect whose runway-facing short edge has a
    pavement GAP to the runway, ADD a synthetic quadrilateral from
    the stub's runway-facing edge straight to the runway boundary.

    Per user 2026-04-27: airport apt.dat data is sometimes
    inconsistent with OSM (construction-era discrepancies); when
    we see a stub that should connect to a runway but the apt.dat
    pavement doesn't bridge the gap, project a straight line over
    so the connecting junction has continuous coverage instead of
    a "weird shape with a big gap".

    The runway-facing short edge is identified as whichever of the
    stub's two short edges (corners [0,3] or [1,2]) has its
    midpoint nearest the runway boundary.  If that midpoint sits
    > ``max_bridge_m`` from the runway, we don't bridge (probably
    not a runway-side stub).
    """
    if residue is None or residue.is_empty:
        return residue
    if runway_union is None or runway_union.is_empty:
        return residue
    rwy_boundary = runway_union.boundary
    additions: List[Polygon] = []
    for rect, axis, role, ref in taxi_rects:
        if role != ROLE_STUB:
            continue
        try:
            rc = list(rect.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Two short edges per ``_rect_from_axis_extended``
        # convention: [3,0] and [1,2].
        short_edges = [(rc[3], rc[0]), (rc[1], rc[2])]
        # Pick the one closer to the runway.
        best_edge = None
        best_d = max_bridge_m
        for (e0, e1) in short_edges:
            mid = Point(0.5 * (e0[0] + e1[0]),
                        0.5 * (e0[1] + e1[1]))
            try:
                d = mid.distance(rwy_boundary)
            except Exception:
                continue
            if d < best_d:
                best_d = d
                best_edge = (e0, e1)
        if best_edge is None or best_d <= 1.0:
            # No runway-side short edge in range, OR the stub
            # already touches the runway — nothing to bridge.
            continue
        # Project each short-edge endpoint to the nearest runway-
        # boundary point.
        from shapely.ops import nearest_points
        e0, e1 = best_edge
        try:
            n0, _ = nearest_points(rwy_boundary, Point(e0))
            n1, _ = nearest_points(rwy_boundary, Point(e1))
        except Exception:
            continue
        # Build the bridge quadrilateral: stub edge → runway edge.
        # Order: e0, e1, n1, n0 so the bridge closes properly.
        try:
            bridge = Polygon([e0, e1,
                              (n1.x, n1.y), (n0.x, n0.y)])
            if not bridge.is_valid:
                bridge = bridge.buffer(0)
            if (bridge.is_empty
                    or bridge.geom_type != "Polygon"):
                continue
            # Don't overlap with the rect itself or the runway.
            bridge = bridge.difference(rect)
            if (not bridge.is_empty
                    and bridge.geom_type == "Polygon"):
                bridge = bridge.difference(runway_union)
            if (not bridge.is_empty
                    and bridge.geom_type == "Polygon"
                    and bridge.area >= 1.0):
                additions.append(bridge)
        except Exception:
            continue
    if additions:
        try:
            residue = unary_union([residue] + additions)
        except Exception:
            pass
    return residue


# ──────────────────────────────────────────────────────────────────
# Top-level builder
# ──────────────────────────────────────────────────────────────────

def build_airport_pavement(icao: str, xplane_root: str,
                            compute_elevations: bool = True
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
    """
    apt_path = APR.find_airport_apt_dat(xplane_root, icao)
    if apt_path is None:
        raise RuntimeError(f"No apt.dat found for {icao}")
    apt = APR.load_airport(apt_path, icao)
    if apt is None:
        raise RuntimeError(f"Could not load airport block for {icao}")

    anchor = _airport_anchor(apt)
    to_m = _projection(anchor)

    layout = PavementLayout(icao=icao, anchor=anchor)

    # Project the apt.dat row-130 airport boundary (in lat/lon) into
    # meter space so downstream emission has it ready when the
    # boundary shape is built.
    if apt.boundary is not None and not apt.boundary.is_empty:
        try:
            from shapely.ops import transform as _shp_transform
            layout.airport_boundary = _shp_transform(
                lambda lon, lat, z=None: to_m(lon, lat),
                apt.boundary)
        except Exception:
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
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] {icao}: widened runway "
                f"{r.desig_a}/{r.desig_b}: {old_w:.1f}m → "
                f"{r.width_m:.1f}m"
                + (f" (centerline shifted {offset:+.1f}m)"
                   if abs(offset) > 0.5 else "")
                + f" — absorbed {len(absorbed)} shoulder polygon(s).\n"
            )
        except Exception:
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
    apt_pav_union: Optional[Polygon] = None
    if pav_polys:
        try:
            apt_pav_union = unary_union(pav_polys)
        except Exception:
            apt_pav_union = None
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
        except Exception:
            apt_bbox_m = None
    try:
        import O4_DSF_Reader as _DSFR
        seen_dsf: set = set()
        all_apt_dats = APR.find_all_airport_apt_dats(xplane_root, icao)
        n_dsf_kept = 0
        n_dsf_dropped_overlay = 0
        n_dsf_dropped_far = 0
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
                        except Exception:
                            pass
                    pav_polys.append(pm)
                    n_dsf_kept += 1
                except Exception:
                    continue
        if (n_dsf_kept or n_dsf_dropped_overlay
                or n_dsf_dropped_far):
            try:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] {icao}: DSF pavement: "
                    f"{n_dsf_kept} kept, "
                    f"{n_dsf_dropped_overlay} dropped as overlay, "
                    f"{n_dsf_dropped_far} dropped as off-airport.\n")
            except Exception:
                pass
    except Exception:
        pass
    pav_union = unary_union(pav_polys) if pav_polys else None
    # Merge near-touching apt.dat polygons so the union is one big
    # coverage (with real holes only) — see ``_merge_near_touching``.
    pav_union = _merge_near_touching(pav_union)
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
                except Exception:
                    continue
        if apron_merged_regions:
            try:
                merged_union = unary_union(apron_merged_regions)
                effective_runway = layout.runway_union.difference(
                    merged_union)
            except Exception:
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

    # ── Load OSM centerlines + relations ─────────────────────────
    nodes, ways, relations = _load_osm_airports(
        xplane_root, icao, anchor[0], anchor[1])
    osm_centerlines = _extract_osm_taxi_centerlines(
        nodes, ways, to_m, rwy_centerlines=rwy_centerlines)

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
                except Exception:
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
        except Exception:
            pass
        terminal_polys.append(pad)
    terminal_union = (unary_union(terminal_polys)
                      if terminal_polys else None)
    for i, tp in enumerate(terminal_polys):
        layout.shapes.append(BuiltShape(
            polygon=tp, role=ROLE_TERMINAL, ref=f"terminal{i+1}"))

    # ── Identify junction node CLUSTERS from OSM topology ───────
    # Any OSM node referenced by ≥2 taxi ways is a potential junction
    # point.  Nodes within JUNCTION_CLUSTER_DIST of each other are
    # merged into one cluster (target junctions often span a whole
    # multi-way intersection, not just a single OSM node).
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
                    except Exception:
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
                            except Exception:
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
                    except Exception:
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
        except Exception:
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
            except Exception:
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
                except Exception:
                    pass
        if parallel_polys:
            try:
                parallel_union = unary_union(parallel_polys)
            except Exception:
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
                    except Exception:
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
                except Exception:
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
    # Use ``pav_union_for_rects`` (full runway subtraction, no
    # apron-merged-runway carve-out) for centerline splitting + rect
    # building.  The apron-merge exclusion was added to ``pav_union``
    # to keep apron junctions seamless across apron-merged runway
    # ends, but it makes the pavement boundary expand outward into
    # those former runway regions — primary parallels at the airport
    # NW (e.g. CYXY's F) then have their natural corridor end (where
    # the pavement narrows back to taxi width) lost in the expanded
    # pav_union, and the apron-interior corner check rejects them.
    _pav_for_rects = getattr(layout, "_pav_union_for_rects", pav_union)
    osm_centerlines = _split_centerlines_at_points(
        osm_centerlines, junction_points, approach_tol_m=25.0,
        pav_union=_pav_for_rects, rwy_union=layout.runway_union,
        rwy_centerlines=rwy_centerlines)

    # ── Build taxi rects from centerlines ────────────────────────
    taxi_rects = _build_taxi_rects(
        osm_centerlines, _pav_for_rects, layout.runway_union,
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
    RUNWAY_ENDPOINT_DIST_M = 80.0
    raw_endpoints_by_ref: Dict[str, List[Tuple[float, float]]] = {}
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

    if layout.runway_union is not None:
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
            except Exception:
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
                    except Exception:
                        continue
            if drop_idx:
                kept = [k for idx, k in enumerate(kept)
                        if idx not in drop_idx]
                changed = True
        if len(kept) < len(taxi_rects):
            try:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] {icao}: dropped "
                    f"{len(taxi_rects) - len(kept)} taxi rect(s) "
                    f"that overlapped another rect.\n")
            except Exception:
                pass
            taxi_rects = kept

    # Per user 2026-04-27 invariant: a junction polygon must NEVER run
    # along a sloping rect's long edge.  When a primary_parallel rect
    # sits FULLY INSIDE the airport's pavement union (apt.dat row-110
    # ⊕ DSF ⊕ OSM-synthetic — i.e. the same union that becomes the
    # residue), the surrounding apron junction unavoidably wraps the
    # rect's long edges as it traces around the rect-shaped hole.
    # The fix is not to emit the rect at all: let the apron absorb
    # the rect's footprint and slope multi-directionally.  Primary
    # parallels that legitimately cross unpaved area are unaffected
    # (their long edges aren't inside pavement).
    #
    # Note: ``pav_union`` is now apt.dat ∪ DSF only (no OSM-synth);
    # that's the right denominator for the absorption check.
    # Including DSF is essential — at SPJC's SE apron, F's long
    # edges are 0 % / 24 % inside row-110 alone but ≈100 % inside
    # apt.dat ∪ DSF.
    taxi_rects = _drop_primary_parallels_embedded_in_pavement(
        taxi_rects, pav_union, runway_polys=runway_polys)

    # Emit taxi rects (already trimmed to narrow-width portion).
    emitted_taxi_rects: List[Polygon] = []
    for rect, axis, role, ref in taxi_rects:
        emitted_taxi_rects.append(rect)
        layout.shapes.append(BuiltShape(
            polygon=rect, role=role, ref=ref, source_axis=axis))

    # ── Junction emission (user 2026-04-23): junctions and
    # aprons are treated identically going forward — both will
    # triangulate and slope multi-directionally (unlike rects,
    # which slope only along their axis).  No distinction needed;
    # emit every connected non-rect non-terminal pavement region
    # as a SINGLE junction polygon, with rect + terminal corners
    # injected as exact shared boundary vertices.  Seamless
    # coverage follows by construction.
    MIN_JUNCTION_AREA_M2 = 50.0   # drop only sliver noise from
                                   # rect-snap inexactness
    SIMPLIFY_TOL_M = 2.0          # min-spacing simplification —
                                   # drops apt.dat-curve vertices
                                   # closer than this to their
                                   # neighbours.  Larger ⇒ fewer
                                   # triangles, fewer slivers; rect /
                                   # terminal seam corners survive
                                   # because they're sharp 90°
                                   # turns that DP keeps.
    RECT_CORNER_TOL_M = 2.0       # corner-to-boundary injection range

    taxi_rect_union = (unary_union(emitted_taxi_rects)
                       if emitted_taxi_rects else None)

    # ── Junction polygons: every pavement region not covered by a
    # rect or a terminal.  User 2026-04-23: simplest polygon that
    # covers the remaining pavement and connects to all rects /
    # terminals; each one will triangulate + slope
    # multi-directionally at elevation time.
    if pav_union is not None:
        residue = pav_union
        if taxi_rect_union is not None:
            residue = residue.difference(taxi_rect_union)
        if terminal_union is not None and not terminal_union.is_empty:
            residue = residue.difference(terminal_union)
        # Defensive: subtract runway even though pav_union already
        # had it removed — floating-point boundary artifacts can
        # leave sub-meter residue slivers overlapping runway.  Use
        # the EFFECTIVE runway union (i.e. with apron-merged
        # regions excluded) so the residue covers parts of runways
        # that pass through aprons.
        _eff_rwy = getattr(layout, "_effective_runway_union",
                           layout.runway_union)
        if _eff_rwy is not None and not _eff_rwy.is_empty:
            residue = residue.difference(_eff_rwy)

        # Per user 2026-04-27 invariant: NO polygon along the long
        # edge of a sloping rect.  Even if apt.dat has pavement
        # extending past a stub's long edge (because the boundary
        # bulges between the stub's two short edges), the residue
        # polygon there must NOT become a junction — it would wrap
        # around the stub's long edge.  Subtract a thin strip just
        # OUTSIDE each stub's long edges from the residue so the
        # resulting junctions stop at the stub's short-edge corners.
        try:
            residue = _clip_residue_at_stub_long_edges(
                residue, taxi_rects)
        except Exception:
            pass

        # Per user 2026-04-27 exception: when a stub's runway-facing
        # short edge has a pavement GAP to the runway (apt.dat
        # construction-era discrepancy), project a synthetic
        # quadrilateral from the stub's short edge straight to the
        # runway boundary so the connecting junction has continuous
        # coverage.
        try:
            residue = _add_stub_to_runway_bridges(
                residue, taxi_rects, layout.runway_union)
        except Exception:
            pass

        # Merge near-touching residue parts (post-subtraction
        # MultiPolygon split is often a numerical artifact of
        # apt.dat polygon-boundary precision rather than a real
        # disconnect — see ``_merge_near_touching``).
        residue = _merge_near_touching(residue)
        parts = ([residue] if residue.geom_type == "Polygon"
                 else list(getattr(residue, "geoms", [])))

        # Collect rect corners and terminal corners — both get
        # injected as shared boundary vertices so every junction
        # seams cleanly to its rect/terminal neighbours.
        seam_points: List[Tuple[float, float]] = []
        for rect, axis, role, ref in taxi_rects:
            rc = list(rect.exterior.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            seam_points.extend(rc)
        for shape in layout.shapes:
            if shape.role != ROLE_TERMINAL:
                continue
            tc = list(shape.polygon.exterior.coords)
            if tc and tc[0] == tc[-1]:
                tc = tc[:-1]
            seam_points.extend(tc)
        # Also the runway corners + axis-aligned runway vertices
        # so stub-to-runway edges are shared.
        if layout.runway_union is not None and not layout.runway_union.is_empty:
            ru = layout.runway_union
            geoms = [ru] if ru.geom_type == "Polygon" else list(
                getattr(ru, "geoms", []))
            for g in geoms:
                if g.geom_type != "Polygon":
                    continue
                rc = list(g.exterior.coords)
                if rc and rc[0] == rc[-1]:
                    rc = rc[:-1]
                seam_points.extend(rc)

        # Dedup seam points that coincide within 0.5 m.
        uniq: List[Tuple[float, float]] = []
        for c in seam_points:
            if not any(math.hypot(c[0] - u[0], c[1] - u[1]) < 0.5
                       for u in uniq):
                uniq.append(c)
        seam_points = uniq

        for part in parts:
            if part.geom_type != "Polygon":
                continue
            if part.area < MIN_JUNCTION_AREA_M2:
                continue
            # Inject seam points as exact boundary vertices.
            part = _insert_points_on_boundary(
                part, seam_points, tol=RECT_CORNER_TOL_M)
            # Light simplification of the remaining apt.dat
            # sub-vertex noise.  Tol kept small so injected seam
            # vertices aren't dropped as near-colinear.
            try:
                simp = part.simplify(SIMPLIFY_TOL_M,
                                     preserve_topology=True)
            except Exception:
                simp = part
            if (simp.is_empty
                    or simp.geom_type != "Polygon"
                    or simp.area < MIN_JUNCTION_AREA_M2):
                continue
            # Re-inject seam points in case simplify dropped any
            # (Douglas-Peucker will remove a vertex on an almost-
            # straight segment even when we need it as a seam).
            simp = _insert_points_on_boundary(
                simp, seam_points, tol=RECT_CORNER_TOL_M)
            if EMIT_JUNCTIONS:
                # If the residue polygon wraps a fully-enclosed
                # rect (holes in the polygon), decompose it into
                # multiple simple polygons that join AROUND the
                # rect instead.  Cleaner human-editable output and
                # avoids hole-splicing artefacts in triangulation.
                pieces = _decompose_polygon_with_holes(
                    simp, min_area_m2=MIN_JUNCTION_AREA_M2)
                for piece in pieces:
                    # Re-inject seam points: cut lines added by
                    # the decomposition split may have introduced
                    # new boundary vertices that DON'T align with
                    # any rect/terminal corner; we still need
                    # corners that lie on this piece's boundary.
                    piece = _insert_points_on_boundary(
                        piece, seam_points,
                        tol=RECT_CORNER_TOL_M)
                    if (piece.is_empty
                            or piece.geom_type != "Polygon"
                            or piece.area < MIN_JUNCTION_AREA_M2):
                        continue
                    layout.shapes.append(BuiltShape(
                        polygon=piece, role=ROLE_JUNCTION))

    # ── Runway-taxiway shared-vertex sync ──
    # Stubs that widen into the runway apron (V1-style) need the
    # runway polygon to carry the projection of their outer
    # corners as vertices, so the junction that wraps around the
    # stub has exact vertex coincidence with the runway.
    rwy_vertex_inserts: List[Tuple[float, float]] = []
    if layout.runway_union is not None and not layout.runway_union.is_empty:
        rwy_boundary = layout.runway_union.boundary
        for rect, axis, role, ref in taxi_rects:
            coords_ax = list(axis.coords)
            if len(coords_ax) < 2:
                continue
            for ax_pt in (coords_ax[0], coords_ax[-1]):
                d_rwy = Point(ax_pt).distance(rwy_boundary)
                if d_rwy > 60.0:
                    continue
                pairs = _rect_end_corners(rect, axis)
                if len(pairs) < 2:
                    continue
                end_idx = 0 if ax_pt == coords_ax[0] else 1
                c1, c2 = pairs[end_idx]
                try:
                    rp1 = nearest_points(rwy_boundary, Point(c1))[0]
                    rp2 = nearest_points(rwy_boundary, Point(c2))[0]
                except Exception:
                    continue
                rwy_vertex_inserts.append((rp1.x, rp1.y))
                rwy_vertex_inserts.append((rp2.x, rp2.y))

    if rwy_vertex_inserts:
        for shape in layout.shapes:
            if shape.role != ROLE_RUNWAY:
                continue
            shape.polygon = _insert_points_on_boundary(
                shape.polygon, rwy_vertex_inserts, tol=2.0)

    # ── Global shared-vertex enforcement (user rule 16) ─────────
    # Cluster all emitted-shape vertices within SHARED_VERTEX_TOL_M
    # and replace each with the cluster centroid.  This guarantees
    # adjacent shapes have EXACT vertex coincidence, which target
    # files enforce and compare_target's v_tgt metric measures.
    _enforce_shared_vertices(layout, tol=SHARED_VERTEX_CLUSTER_TOL_M)

    # Pavement layout invariant (user 2026-04-26): NO shape may
    # overlap another, period.  Runs AFTER shared-vertex collapse
    # because the collapse step can shift boundaries by up to
    # SHARED_VERTEX_CLUSTER_TOL_M / 2 and create tiny new overlaps.
    # Re-run shared-vertex collapse afterwards because clipping
    # introduces new vertices at intersection points that may sit
    # within the cluster tol of an existing vertex on the
    # higher-priority shape's edge.
    _drop_overlap_against_fixed_shapes(layout, icao=icao)
    _enforce_shared_vertices(layout, tol=SHARED_VERTEX_CLUSTER_TOL_M)

    # Validate the invariant: every vertex of every shape must either
    # be unique (distance > 2 × tol to any other shape's vertex) OR
    # exactly equal to a vertex on an adjacent shape.  No "close but
    # not equal" drift is permitted.
    _validate_shared_vertex_invariant(layout,
                                      tol=SHARED_VERTEX_CLUSTER_TOL_M)

    # ── Phase-2 elevations ──────────────────────────────────────
    if compute_elevations:
        _compute_elevations(
            layout, icao, xplane_root, apt,
            osm_nodes=nodes, osm_ways=ways, to_m=to_m,
            apron_candidates_m=apron_candidates)
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
        # Final WARN summary — emitted after every elevation pass
        # has run so the count reflects what the OSM emitter will
        # actually write to disk.  Earlier reports (mid-pipeline)
        # over-counted because they ran before the smoother /
        # shared-vertex agree / rect-corner snap chain converged.
        _report_within_shape_violations(layout, icao)
        # Per user 2026-04-28: emit a 5 m-wide ribbon polygon
        # tracing the airport boundary (apt.dat row-130) with
        # per-vertex altitudes clamped to ≤ 3 % grade from the
        # nearest runway within 400 m, falling back to DEM
        # beyond.  Provides the elevation transition between
        # airport pavement and surrounding terrain.
        try:
            _lat0, _lon0 = layout.anchor
            _tile_lat = int(math.floor(_lat0))
            _tile_lon = int(math.floor(_lon0))
            _dem = _load_airport_dem(_lat0, _lon0)
            n_b = _emit_airport_boundary_shape(
                layout, _dem, _tile_lat, _tile_lon)
            if n_b:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] emitted "
                    f"{n_b} airport-boundary shape piece(s).\n")
            # Then emit DEM-bridge polygons inside the boundary
            # wherever the clamped boundary altitude differs from
            # raw DEM by > 5 m (per user 2026-04-28).
            try:
                n_br = _emit_boundary_dem_bridge(
                    layout, _dem, _tile_lat, _tile_lon)
                if n_br:
                    import sys as _sys
                    _sys.stderr.write(
                        f"  [pav-builder] emitted "
                        f"{n_br} boundary→DEM bridge "
                        f"polygon(s).\n")
            except Exception:
                pass
        except Exception:
            pass

    return layout


# ══════════════════════════════════════════════════════════════════
# Phase-2: elevations
# ══════════════════════════════════════════════════════════════════

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
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] {fname} missing; downloading "
                    f"DEM tile via Ortho4XP elevation provider...\n")
            except Exception:
                pass
            dem = _DEM.DEM(tile_lat, tile_lon)
    except Exception as exc:
        try:
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] WARN: DEM load/download failed for "
                f"{fname}: {exc}\n")
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
            import O4_Auto_Patch as _AP
            cifp_runways = _AP.parse_cifp_file(cifp_path)
            if cifp_runways:
                pairs = _AP.pair_runways(cifp_runways)

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
                    continue
                kept_shapes.append(sh)
                if sh.polygon is not None and not sh.polygon.is_empty:
                    kept_polys.append(sh.polygon)
            if n_dropped:
                try:
                    import sys as _sys
                    _sys.stderr.write(
                        f"  [pav-builder] {icao}: dropped "
                        f"{n_dropped} runway segment(s) "
                        f"apron-merged.\n")
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
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] {icao}: resolved "
                    f"{n_crossings} runway crossing(s) into "
                    f"junction polygon(s).\n")
            except Exception:
                pass
            new_runway_polys = [
                s.polygon for s in layout.shapes
                if s.role == ROLE_RUNWAY
                and s.polygon is not None
                and not s.polygon.is_empty]

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
            _sys.stderr.write(
                f"  [pav-builder] terminal({shape.ref or '?'}) "
                f"altitude {new_alt} m (max grade-compliant from "
                f"runway corners within "
                f"{TERMINAL_NEIGHBOUR_RADIUS_M:.0f} m, "
                f"DEM-median ceiling = {dem_str}).\n")
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


def _push_junction_vertices_off_taxi_rect_edges(
        layout: "PavementLayout",
        edge_tol_m: float = 0.5,
        corner_tol_m: float = 2.0,
        edge_gap_m: float = 1.0,
        ) -> int:
    """Per the user 2026-04-24 invariant: a junction polygon may
    share a vertex with a taxi rect ONLY at one of the rect's 4
    corners.  A junction vertex landing on the INTERIOR of a rect
    edge would split that edge at render time and break the rect's
    altitude_high/altitude_low slope convention.

    Per user 2026-04-29: junctions and taxiway/stub/runway rects
    should be handled the SAME WAY as runway-runway crossings —
    the junction should connect to either the LOW or HIGH short
    edge of a rect (i.e. its corners), never to the rect's edge
    interior.

    Two-stage policy applied to every junction ring vertex:

      Stage 1 — collapse redundant edge-interior vertices.
        If a vertex lies on the interior of a rect edge AND its
        ring-adjacent neighbours are both at corners of the SAME
        rect (one each, on the two ends of THAT edge), the vertex
        is redundant.  The junction's ring already connects the
        two corners; the intermediate vertex just splits a single
        rect edge into two pieces.  Drop the vertex — the junction
        edge then runs corner-to-corner along the rect's short
        edge (HIGH-side or LOW-side) cleanly.
      Stage 2 — corner snap / push for survivors.
        For each remaining vertex:
          * Within ``corner_tol_m`` of a rect corner ⇒ snap to
            that corner.
          * Within ``edge_tol_m`` of a rect edge interior ⇒ push
            ``edge_gap_m`` perpendicular outside the rect.
          * Otherwise ⇒ leave alone.

    Geometric only: doesn't touch elevations.  Runway corners are
    treated identically to taxi rect corners so junction vertices
    snap there too — the user's "same logic for all runway
    intersections" requirement.

    Returns the number of junction polygons modified.
    """
    rect_roles = {
        ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
        ROLE_STUB, ROLE_CROSS_CONNECTOR, ROLE_RUNWAY}
    rects: List[Tuple[Polygon, List[Tuple[float, float]]]] = []
    for s in layout.shapes:
        if s.role not in rect_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
        except Exception:
            continue
        if len(coords) != 4:
            continue
        rects.append((s.polygon, coords))
    if not rects:
        return 0

    corner_tol2 = corner_tol_m * corner_tol_m

    def _on_edge_between_corners(
            x: float, y: float,
            corners: List[Tuple[float, float]],
            ) -> Optional[int]:
        """Return the edge index (0-3) the point lies on (within
        ``edge_tol_m`` of an edge interior, t ∈ (ε, 1-ε)), or
        None if the point isn't on any rect edge interior."""
        for i in range(4):
            ax, ay = corners[i]
            bx, by = corners[(i + 1) % 4]
            dx = bx - ax
            dy = by - ay
            seg_len_sq = dx * dx + dy * dy
            if seg_len_sq <= 0.01:
                continue
            t = ((x - ax) * dx + (y - ay) * dy) / seg_len_sq
            if t <= 0.001 or t >= 0.999:
                continue
            cx_proj = ax + t * dx
            cy_proj = ay + t * dy
            d_sq = ((x - cx_proj) ** 2
                    + (y - cy_proj) ** 2)
            if d_sq <= edge_tol_m * edge_tol_m:
                return i
        return None

    def _at_corner_index(
            x: float, y: float,
            corners: List[Tuple[float, float]],
            ) -> Optional[int]:
        """Return the corner index (0-3) the point lies at
        (within ``corner_tol_m``), or None."""
        for ci, (cx, cy) in enumerate(corners):
            if (x - cx) ** 2 + (y - cy) ** 2 <= corner_tol2:
                return ci
        return None

    def _push_off(
            x: float, y: float,
            rect_poly: Polygon,
            corners: List[Tuple[float, float]],
            edge_idx: int,
            ) -> Tuple[float, float]:
        """Push the point ``edge_gap_m`` perpendicular to the
        rect edge ``edge_idx``, toward the OUTSIDE of the rect."""
        ax, ay = corners[edge_idx]
        bx, by = corners[(edge_idx + 1) % 4]
        dx = bx - ax
        dy = by - ay
        seg_len = math.sqrt(dx * dx + dy * dy)
        if seg_len <= 0.01:
            return (x, y)
        t = ((x - ax) * dx + (y - ay) * dy) / (seg_len * seg_len)
        cx_proj = ax + t * dx
        cy_proj = ay + t * dy
        perp_x = -dy / seg_len
        perp_y = dx / seg_len
        test_x = cx_proj + perp_x * 0.1
        test_y = cy_proj + perp_y * 0.1
        if rect_poly.contains(Point(test_x, test_y)):
            perp_x = -perp_x
            perp_y = -perp_y
        return (cx_proj + perp_x * edge_gap_m,
                cy_proj + perp_y * edge_gap_m)

    n_modified = 0
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        try:
            ring = list(shape.polygon.exterior.coords)
        except Exception:
            continue
        # Drop closing repeat for ring traversal.
        if ring and ring[0] == ring[-1]:
            ring_open = ring[:-1]
        else:
            ring_open = ring
        n_v = len(ring_open)
        if n_v < 3:
            continue

        # Stage 1: collapse redundant edge-interior vertices.
        # A vertex v is redundant if its ring-prev and ring-next
        # neighbours are at the two corners of the same rect edge
        # AND v itself lies on that edge's interior.
        keep_mask = [True] * n_v
        for i in range(n_v):
            vx, vy = ring_open[i]
            px, py = ring_open[(i - 1) % n_v]
            nx, ny = ring_open[(i + 1) % n_v]
            for rect_poly, corners in rects:
                p_corner = _at_corner_index(px, py, corners)
                n_corner = _at_corner_index(nx, ny, corners)
                if p_corner is None or n_corner is None:
                    continue
                # The two neighbour corners must be adjacent
                # corners (i.e. share an edge).  Adjacent corner
                # pairs: (0,1), (1,2), (2,3), (3,0).
                diff = abs(p_corner - n_corner)
                if diff != 1 and diff != 3:
                    continue
                # The edge between them is index = min(...) if
                # adjacent, but for the wrap (0,3 / 3,0) it's
                # edge 3.  Just identify by the corner pair.
                edge_idx = (
                    min(p_corner, n_corner)
                    if diff == 1 else 3)
                v_edge = _on_edge_between_corners(vx, vy, corners)
                if v_edge != edge_idx:
                    continue
                # Vertex v sits on the edge between two corners
                # that are already in the ring as neighbours.
                # Drop it — the junction's ring will go
                # corner-to-corner along the rect edge.
                keep_mask[i] = False
                break
        # Survivor index list lets us drop the matching entries from
        # ``shape.node_altitudes`` after the rebuild — the dropped
        # vertices' altitudes are no longer needed but the kept
        # vertices' altitudes ARE still valid (the solver's elevation
        # field depends on the polygon being valid AND on per-vertex
        # values; nuking them all forces re-derivation from scratch
        # and reintroduces grade violations the solver already fixed).
        survivor_idx = [i for i in range(n_v) if keep_mask[i]]
        ring_after_collapse = [ring_open[i] for i in survivor_idx]
        n_collapsed = n_v - len(ring_after_collapse)

        # Stage 2: corner-snap / edge-push the survivors.
        new_ring: List[Tuple[float, float]] = []
        n_snapped = 0
        n_pushed = 0
        for vx, vy in ring_after_collapse:
            target = (vx, vy)
            for rect_poly, corners in rects:
                ci = _at_corner_index(vx, vy, corners)
                if ci is not None:
                    target = corners[ci]
                    if target != (vx, vy):
                        n_snapped += 1
                    break
                ei = _on_edge_between_corners(vx, vy, corners)
                if ei is not None:
                    target = _push_off(
                        vx, vy, rect_poly, corners, ei)
                    if target != (vx, vy):
                        n_pushed += 1
                    break
            new_ring.append(target)

        if n_collapsed == 0 and n_snapped == 0 and n_pushed == 0:
            continue

        # Re-close the ring and rebuild the polygon.
        if new_ring and new_ring[0] != new_ring[-1]:
            new_ring_closed = new_ring + [new_ring[0]]
        else:
            new_ring_closed = new_ring
        try:
            new_poly = Polygon(new_ring_closed,
                                list(shape.polygon.interiors))
            buffer_repaired = False
            if not new_poly.is_valid:
                # The original ring may already be self-intersecting
                # — buffer(0) can return a MultiPolygon with one
                # main piece + tiny artifacts.  Take the largest
                # Polygon piece so the rect-edge fix still applies.
                fixed = new_poly.buffer(0)
                if fixed.geom_type == "Polygon":
                    new_poly = fixed
                    buffer_repaired = True
                elif (fixed.geom_type == "MultiPolygon"
                        and not fixed.is_empty):
                    new_poly = max(
                        fixed.geoms, key=lambda g: g.area)
                    buffer_repaired = True
                else:
                    new_poly = None
            if (new_poly is not None
                    and new_poly.geom_type == "Polygon"
                    and not new_poly.is_empty):
                shape.polygon = new_poly
                n_new = len(list(new_poly.exterior.coords)) - 1
                # Preserve per-vertex altitudes where we can.  When
                # buffer(0) restructured the ring (vertex order /
                # count not predictable from the input), fall back
                # to None and let downstream re-derive.  When Stage
                # 1 collapsed K vertices but Stage 2 only snapped/
                # pushed in place, the survivor index list maps the
                # new ring to the original altitudes 1:1.
                if shape.node_altitudes is None:
                    pass
                elif buffer_repaired:
                    if len(shape.node_altitudes) - 1 != n_new:
                        shape.node_altitudes = None
                elif n_new == len(survivor_idx):
                    # Closing-vertex repeat: keep one trailing slot.
                    src = shape.node_altitudes
                    src_open = (
                        src[:-1]
                        if len(src) - 1 == n_v
                        else src[:n_v])
                    if len(src_open) == n_v:
                        new_alts = [src_open[i] for i in survivor_idx]
                        shape.node_altitudes = (
                            new_alts + [new_alts[0]]
                            if new_alts else None)
                    else:
                        if len(shape.node_altitudes) - 1 != n_new:
                            shape.node_altitudes = None
                else:
                    if len(shape.node_altitudes) - 1 != n_new:
                        shape.node_altitudes = None
                n_modified += 1
        except Exception:
            pass
    return n_modified


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
        import sys as _sys
        _sys.stderr.write(
            f"  [pav-builder] {icao}: unified Laplacian solver "
            f"converged in {it + 1}/{max_iters} iters "
            f"({elapsed:.2f} s); applied to "
            f"{n_terms} terminal(s), {n_rects} rect(s), "
            f"{n_junctions} junction(s).\n")
    except Exception:
        pass


def _emit_airport_boundary_shape(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        strip_half_width_m: float = 2.5,
        runway_clamp_radius_m: float = 400.0,
        runway_clamp_grade: float = 0.03,
        densify_step_m: float = 25.0,
        ) -> int:
    """Emit a node_altitudes polygon tracing the airport boundary
    (apt.dat row-130) at ``2 × strip_half_width_m`` width.

    Per user 2026-04-28: an airport-perimeter "ribbon" with
    controlled per-vertex altitudes provides the elevation
    transition between the airport's pavement and the surrounding
    DEM.  Vertices within ``runway_clamp_radius_m`` of any runway
    are clamped to the runway elevation ± ``runway_clamp_grade``
    × distance (default 3 % grade); vertices beyond the radius
    follow the DEM directly.

    Implementation:
      1. Buffer the boundary's exterior LineString by
         ``strip_half_width_m`` to produce a closed strip polygon.
         The strip naturally has an interior ring (the airport
         interior shrunk inward by the buffer).
      2. Decompose the holed strip into simple non-holed pieces
         via ``_decompose_polygon_with_holes`` so X-Plane's patch
         parser (which drops interior rings) renders the strip
         correctly.
      3. For each piece, densify boundary segments to
         ``densify_step_m`` so per-vertex altitude clamping
         resolves at a useful spatial frequency.
      4. Compute per-vertex altitudes against the runway-distance
         rule + DEM.
      5. Append each piece as a ``ROLE_BOUNDARY`` BuiltShape.

    Returns the number of boundary shape pieces emitted.
    """
    if layout.airport_boundary is None or layout.airport_boundary.is_empty:
        return 0
    from shapely.geometry import LineString as _LS, Polygon as _Polygon
    from shapely.geometry import Point as _Point
    from shapely.ops import nearest_points as _nearest_points

    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH
    def m_to_ll(x: float, y: float) -> Tuple[float, float]:
        lat = lat0 + math.degrees(y / R)
        lon = lon0 + math.degrees(x / (R * cos0))
        return lat, lon

    # Pre-collect runway polygons + their elevation samplers for
    # the per-vertex distance / clamp lookup.
    runway_shapes: List[BuiltShape] = [
        s for s in layout.shapes
        if s.role == ROLE_RUNWAY
        and s.polygon is not None
        and not s.polygon.is_empty]
    if not runway_shapes:
        return 0

    def _runway_clamped_alt(x: float, y: float) -> Optional[float]:
        """Return DEM at (x, y) clamped to ``[runway_e - g·d,
        runway_e + g·d]`` when within ``runway_clamp_radius_m`` of
        any runway, else raw DEM, else None."""
        try:
            lat, lon = m_to_ll(x, y)
            dem_e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            dem_e = None
        # Find nearest runway and its elevation at the nearest point.
        best_d = float('inf')
        best_e = None
        pt = _Point(x, y)
        for s in runway_shapes:
            try:
                d = s.polygon.distance(pt)
            except Exception:
                continue
            if d >= best_d:
                continue
            try:
                if d == 0.0:
                    np_x, np_y = x, y
                else:
                    np = _nearest_points(s.polygon, pt)[0]
                    np_x, np_y = np.x, np.y
                e = _sample_runway_segment_elev(s, np_x, np_y)
            except Exception:
                e = None
            if e is None:
                continue
            best_d = d
            best_e = e
        if best_e is None:
            return dem_e
        if best_d > runway_clamp_radius_m:
            return dem_e
        band = best_d * runway_clamp_grade
        lo = best_e - band
        hi = best_e + band
        if dem_e is None:
            return 0.5 * (lo + hi)
        if dem_e < lo:
            return lo
        if dem_e > hi:
            return hi
        return dem_e

    def _densify_ring(coords: List[Tuple[float, float]]
                      ) -> List[Tuple[float, float]]:
        """Insert intermediate points so consecutive vertices are
        ≤ ``densify_step_m`` apart.  Closes the ring at the end."""
        if not coords:
            return coords
        if coords[0] == coords[-1]:
            coords = coords[:-1]
        out: List[Tuple[float, float]] = []
        n = len(coords)
        for i in range(n):
            a = coords[i]
            b = coords[(i + 1) % n]
            out.append(a)
            d = math.hypot(b[0] - a[0], b[1] - a[1])
            if d > densify_step_m:
                steps = max(1, int(d / densify_step_m))
                for k in range(1, steps):
                    t = k / steps
                    out.append((a[0] + t * (b[0] - a[0]),
                                a[1] + t * (b[1] - a[1])))
        out.append(out[0])
        return out

    # 1. Boundary line → 5 m strip polygon (with inner ring if the
    #    airport is large enough that 2.5 m × 2 < interior radius).
    boundary_geom = layout.airport_boundary
    if boundary_geom.geom_type == "Polygon":
        ext_rings = [boundary_geom.exterior]
    elif boundary_geom.geom_type == "MultiPolygon":
        ext_rings = [g.exterior for g in boundary_geom.geoms]
    else:
        return 0
    # Build a union of all existing pavement shapes — the boundary
    # ribbon is meant to control elevations OUTSIDE the pavement
    # (grass / approach lights / ramp).  Subtract pavement from the
    # strip so the boundary doesn't overlap runway / taxi rects /
    # apron junctions etc. (would otherwise fail the no-self-
    # overlap regression test and double-up elevation tags at the
    # airport perimeter).
    pavement_polys = [
        s.polygon for s in layout.shapes
        if s.polygon is not None
        and not s.polygon.is_empty
        and s.role != ROLE_BOUNDARY]
    pavement_union: Optional[Polygon] = None
    if pavement_polys:
        try:
            pavement_union = unary_union(pavement_polys)
        except Exception:
            pavement_union = None
    n_emitted = 0
    for ring in ext_rings:
        ring_coords = list(ring.coords)
        try:
            line = _LS(ring_coords)
            strip = line.buffer(strip_half_width_m,
                                cap_style=2, join_style=2)
            if not strip.is_valid:
                strip = strip.buffer(0)
        except Exception:
            continue
        if strip.is_empty:
            continue
        if pavement_union is not None and not pavement_union.is_empty:
            try:
                strip = strip.difference(pavement_union)
            except Exception:
                pass
            if strip.is_empty:
                continue
            if not strip.is_valid:
                strip = strip.buffer(0)
        if strip.geom_type == "MultiPolygon":
            strip_polys = list(strip.geoms)
        elif strip.geom_type == "Polygon":
            strip_polys = [strip]
        else:
            continue
        # 2. Decompose any holes.
        all_pieces: List[Polygon] = []
        for sp in strip_polys:
            try:
                pieces = _decompose_polygon_with_holes(
                    sp, min_area_m2=10.0, max_depth=8)
            except Exception:
                pieces = [sp]
            for p in pieces:
                if p.is_empty or p.geom_type != "Polygon":
                    continue
                all_pieces.append(p)
        # 3-5. Densify, compute altitudes, emit.
        for piece in all_pieces:
            try:
                exterior = list(piece.exterior.coords)
            except Exception:
                continue
            dense = _densify_ring(exterior)
            if len(dense) < 4:
                continue
            try:
                new_poly = _Polygon(dense)
                if not new_poly.is_valid:
                    new_poly = new_poly.buffer(0)
                if (new_poly.is_empty
                        or new_poly.geom_type != "Polygon"):
                    continue
            except Exception:
                continue
            # Re-extract the (post-buffer-cleanup) exterior so the
            # node_altitudes count matches polygon.exterior.coords.
            new_coords = list(new_poly.exterior.coords)
            alts: List[float] = []
            for (cx, cy) in new_coords:
                e = _runway_clamped_alt(cx, cy)
                if e is None:
                    e = 0.0
                alts.append(round(float(e), 1))
            layout.shapes.append(BuiltShape(
                polygon=new_poly,
                role=ROLE_BOUNDARY,
                ref="airport_boundary",
                node_altitudes=alts))
            n_emitted += 1
    return n_emitted


def _emit_boundary_dem_bridge(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        gap_threshold_m: float = 5.0,
        bridge_depth_m: float = 100.0,
        densify_step_m: float = 25.0,
        runway_clamp_radius_m: float = 400.0,
        runway_clamp_grade: float = 0.03,
        ) -> int:
    """Emit a wider "bridge" polygon INSIDE the airport boundary
    where the boundary's clamped altitude differs from the raw DEM
    by more than ``gap_threshold_m``.

    Per user 2026-04-28: when the boundary ribbon is forced (by the
    runway-distance clamp at ≤ 3 % grade) to a value that disagrees
    with the natural terrain DEM by > 5 m, X-Plane renders a
    valley/cliff between the 5 m boundary ribbon and the surrounding
    terrain inside the airport perimeter.  The bridge polygon is a
    larger transition strip whose OUTER edge sits on the airport
    perimeter at the boundary's clamped altitude and whose INNER
    edge sits ``bridge_depth_m`` further inside the airport at the
    raw DEM altitude.  Per-vertex altitudes interpolate linearly
    between the two edges, giving X-Plane a gradual surface to
    descend / ascend over instead of a single hard step.

    OUTSIDE the airport boundary X-Plane keeps falling directly to
    DEM (no bridge needed there) — the user explicitly scoped this
    feature to the interior side only.

    Implementation:
      1. Densify the airport-boundary line to ≤ ``densify_step_m``.
      2. For each densified vertex, sample raw DEM and the
         runway-clamped altitude (same rule as the 5 m ribbon).
         Mark vertex if |gap| > ``gap_threshold_m``.
      3. Group consecutive marked vertices into "bridge runs"
         (with a 1-vertex slack so isolated unmarked vertices in
         the middle of a long gap don't split the run).
      4. For each run, build an inward-offset polygon
         (``bridge_depth_m`` inward from the boundary line) and
         clip it against any existing pavement / boundary ribbon.
      5. Emit per-vertex altitudes: outer edge = clamped, inner
         edge = DEM, with shape vertices on the boundary side
         tagged ``clamped`` and inner-edge vertices tagged DEM.
    """
    if (layout.airport_boundary is None
            or layout.airport_boundary.is_empty):
        return 0
    from shapely.geometry import LineString as _LS, Point as _Point
    from shapely.geometry import Polygon as _Polygon

    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def m_to_ll(x: float, y: float) -> Tuple[float, float]:
        lat = lat0 + math.degrees(y / R)
        lon = lon0 + math.degrees(x / (R * cos0))
        return lat, lon

    runway_shapes: List[BuiltShape] = [
        s for s in layout.shapes
        if s.role == ROLE_RUNWAY
        and s.polygon is not None
        and not s.polygon.is_empty]
    if not runway_shapes:
        return 0

    def _clamped_alt(x: float, y: float) -> Optional[float]:
        try:
            lat, lon = m_to_ll(x, y)
            dem_e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            dem_e = None
        best_d = float('inf')
        best_e = None
        from shapely.ops import nearest_points as _np
        pt = _Point(x, y)
        for s in runway_shapes:
            try:
                d = s.polygon.distance(pt)
            except Exception:
                continue
            if d >= best_d:
                continue
            try:
                if d == 0.0:
                    np_x, np_y = x, y
                else:
                    np = _np(s.polygon, pt)[0]
                    np_x, np_y = np.x, np.y
                e = _sample_runway_segment_elev(s, np_x, np_y)
            except Exception:
                e = None
            if e is None:
                continue
            best_d = d
            best_e = e
        if best_e is None:
            return dem_e
        if best_d > runway_clamp_radius_m:
            return dem_e
        band = best_d * runway_clamp_grade
        lo = best_e - band
        hi = best_e + band
        if dem_e is None:
            return 0.5 * (lo + hi)
        if dem_e < lo:
            return lo
        if dem_e > hi:
            return hi
        return dem_e

    def _dem_alt(x: float, y: float) -> Optional[float]:
        try:
            lat, lon = m_to_ll(x, y)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            return None

    boundary_geom = layout.airport_boundary
    if boundary_geom.geom_type == "Polygon":
        rings = [boundary_geom]
    elif boundary_geom.geom_type == "MultiPolygon":
        rings = list(boundary_geom.geoms)
    else:
        return 0

    # Compose existing PAVEMENT union (excluding ROLE_BOUNDARY
    # shapes — the just-emitted 5 m ribbon's centerline IS the
    # boundary line, so the ribbon would reject every boundary
    # vertex from the pre-filter below).  The bridge is meant to
    # avoid overlapping real pavement (runways / taxis / aprons /
    # terminals); it's placed alongside the ribbon, not on top of
    # other pavement.
    pavement_polys = [
        s.polygon for s in layout.shapes
        if s.role != ROLE_BOUNDARY
        and s.polygon is not None
        and not s.polygon.is_empty]
    pavement_union: Optional[Polygon] = None
    if pavement_polys:
        try:
            pavement_union = unary_union(pavement_polys)
        except Exception:
            pavement_union = None
    # Separately track the boundary ribbon — its centerline matches
    # the boundary line, so the bridge polygon overlaps the ribbon
    # in its inner 2.5 m by construction.  The bridge must be
    # trimmed against the ribbon to satisfy the no-self-overlap
    # geometry test.
    ribbon_polys = [
        s.polygon for s in layout.shapes
        if s.role == ROLE_BOUNDARY
        and s.ref == "airport_boundary"
        and s.polygon is not None
        and not s.polygon.is_empty]
    ribbon_union: Optional[Polygon] = None
    if ribbon_polys:
        try:
            ribbon_union = unary_union(ribbon_polys)
        except Exception:
            ribbon_union = None

    # Pre-collect pavement EDGE points with altitudes — used for
    # nearest-pavement lookup when assigning per-vertex altitudes
    # to the bridge polygon.  Per user 2026-04-28: bridge vertices
    # adjacent to pavement must match the pavement's altitude (not
    # raw DEM) so the bridge actually FILLS the gap between
    # boundary and pavement instead of creating its own valley.
    pav_edge_pts: List[Tuple[float, float, float]] = []
    for s in layout.shapes:
        if s.role == ROLE_BOUNDARY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        # Per-vertex altitudes for junctions / boundary; rect tags
        # for sloping rects.
        if s.role in (ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                       ROLE_CROSS_CONNECTOR):
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
                pav_edge_pts.append((float(x), float(y), float(a)))
        elif s.node_altitudes:
            for (x, y), a in zip(coords,
                                  s.node_altitudes[:len(coords)]):
                pav_edge_pts.append((float(x), float(y), float(a)))
        elif s.altitude is not None:
            for x, y in coords:
                pav_edge_pts.append((float(x), float(y),
                                     float(s.altitude)))

    def _nearest_pav_alt(x: float, y: float,
                         max_d_m: float = 500.0
                         ) -> Optional[Tuple[float, float]]:
        """Return ``(alt, distance_m)`` for the nearest pavement
        edge point within ``max_d_m`` of ``(x, y)``, or None when
        no pavement is in range."""
        best_d2 = max_d_m * max_d_m
        best_alt: Optional[float] = None
        for px, py, pa in pav_edge_pts:
            d2 = (x - px) * (x - px) + (y - py) * (y - py)
            if d2 < best_d2:
                best_d2 = d2
                best_alt = pa
        if best_alt is None:
            return None
        return (best_alt, math.sqrt(best_d2))

    n_emitted = 0
    for boundary_poly in rings:
        try:
            ext_coords = list(boundary_poly.exterior.coords)
        except Exception:
            continue
        if len(ext_coords) < 4:
            continue
        # Densify the boundary line.
        if ext_coords[0] == ext_coords[-1]:
            ext_coords = ext_coords[:-1]
        dense: List[Tuple[float, float]] = []
        n = len(ext_coords)
        for i in range(n):
            ax, ay = ext_coords[i]
            bx, by = ext_coords[(i + 1) % n]
            dense.append((ax, ay))
            d = math.hypot(bx - ax, by - ay)
            if d > densify_step_m:
                steps = max(1, int(d / densify_step_m))
                for k in range(1, steps):
                    t = k / steps
                    dense.append((ax + t * (bx - ax),
                                  ay + t * (by - ay)))
        if len(dense) < 4:
            continue
        # Per-vertex clamped + DEM + gap.
        per_vert: List[Tuple[float, float, float, float]] = []
        for x, y in dense:
            ca = _clamped_alt(x, y)
            da = _dem_alt(x, y)
            if ca is None or da is None:
                per_vert.append((x, y, float('nan'), float('nan')))
                continue
            per_vert.append((x, y, float(ca), float(da)))

        # A vertex is "needs-bridge" only if (a) the gap exceeds
        # the threshold AND (b) the boundary line at that vertex
        # is NOT already inside pavement (a runway / taxi rect
        # extending to the perimeter doesn't need a transition —
        # pavement is right there).  Pre-filtering on (b) avoids
        # building bridge polygons that overlap pavement; the
        # subsequent pavement-difference would otherwise leave
        # vertices stranded on sloping-rect edges (test
        # ``test_no_vertex_on_sloping_rect_edge``).
        from shapely.geometry import Point as _P2
        marked = []
        for i, v in enumerate(per_vert):
            if math.isnan(v[2]) or math.isnan(v[3]):
                continue
            if abs(v[2] - v[3]) <= gap_threshold_m:
                continue
            if pavement_union is not None and not pavement_union.is_empty:
                try:
                    if pavement_union.distance(
                            _P2(v[0], v[1])) < 5.0:
                        continue
                except Exception:
                    pass
            marked.append((i, v))
        if not marked:
            continue
        # Group consecutive marked vertices into runs (treat the
        # boundary as cyclic; allow 1-vertex unmarked slack).
        marked_idx = sorted(set(m[0] for m in marked))
        N = len(per_vert)
        runs: List[List[int]] = []
        if marked_idx:
            cur = [marked_idx[0]]
            for idx in marked_idx[1:]:
                # Distance along ring, accounting for wrap.
                gap_idx = idx - cur[-1]
                if gap_idx <= 2:
                    cur.append(idx)
                else:
                    runs.append(cur)
                    cur = [idx]
            runs.append(cur)
            # Wrap merge: last run end-of-ring + first run
            # start-of-ring close ⇒ merge.
            if len(runs) >= 2:
                tail = runs[-1][-1]
                head = runs[0][0]
                if (N - tail) + head <= 2:
                    runs[0] = runs[-1] + runs[0]
                    runs.pop()

        for run in runs:
            if len(run) < 2:
                continue
            # Outer edge: the boundary line vertices for the run,
            # in order.
            outer_pts = [(per_vert[i][0], per_vert[i][1])
                         for i in run]
            if len(outer_pts) < 2:
                continue
            try:
                outer_line = _LS(outer_pts)
            except Exception:
                continue
            if outer_line.is_empty or outer_line.length < 1.0:
                continue
            # Build inner offset on whichever side is INSIDE the
            # airport boundary polygon.
            inner_line = None
            for side in ("left", "right"):
                try:
                    off = outer_line.parallel_offset(
                        bridge_depth_m, side=side, join_style=2)
                except Exception:
                    off = None
                if off is None or off.is_empty:
                    continue
                # Probe a midpoint of the offset to test
                # containment in the airport boundary.
                try:
                    mid = off.interpolate(0.5, normalized=True)
                    if boundary_poly.contains(mid):
                        inner_line = off
                        break
                except Exception:
                    continue
            if inner_line is None or inner_line.is_empty:
                continue
            # Build the bridge polygon: outer line + reversed
            # inner line.  parallel_offset on the LEFT side returns
            # a line in REVERSED order; on the RIGHT side it's in
            # the same order.  Either way, we walk the outer
            # forward then close along the inner.  Determine
            # winding by trying both and keeping the valid one.
            in_coords = list(inner_line.coords)
            ring1 = list(outer_pts) + list(reversed(in_coords))
            ring2 = list(outer_pts) + list(in_coords)
            bridge_poly: Optional[Polygon] = None
            for cand in (ring1, ring2):
                try:
                    p = _Polygon(cand)
                    if not p.is_valid:
                        p = p.buffer(0)
                    if (not p.is_empty
                            and p.geom_type == "Polygon"
                            and p.area > 100.0):
                        bridge_poly = p
                        break
                except Exception:
                    continue
            if bridge_poly is None:
                continue
            # Per user 2026-04-29: bridge polygons must stop 5 m
            # short of any runway — never connect directly to
            # runway pavement.  The bridge is a transition strip
            # between the boundary ribbon and natural terrain;
            # forcing a runway corner / edge into its outline
            # would re-introduce sloping-rect-edge vertex
            # violations and create grade conflicts at the
            # runway interface.  Subtract a 5 m-buffered runway
            # union so the bridge keeps a clean gap.
            try:
                runway_union = unary_union(
                    [s.polygon for s in runway_shapes])
                if runway_union is not None and not runway_union.is_empty:
                    bridge_poly = bridge_poly.difference(
                        runway_union.buffer(5.0))
            except Exception:
                pass
            if (bridge_poly.is_empty
                    or bridge_poly.geom_type
                    not in ("Polygon", "MultiPolygon")):
                continue
            # Subtract NON-SLOPING pavement (junctions, terminals)
            # and the boundary ribbon from the bridge.  These can
            # overlap the bridge by hundreds of square metres
            # without sharing vertices with sloping rects, so
            # subtracting them is safe and necessary for the
            # no-self-overlap test.
            non_sloping_pav_polys: List[Polygon] = [
                s.polygon for s in layout.shapes
                if s.role in (ROLE_JUNCTION, ROLE_TERMINAL)
                and s.polygon is not None
                and not s.polygon.is_empty]
            for sub_geom in non_sloping_pav_polys:
                try:
                    bridge_poly = bridge_poly.difference(sub_geom)
                except Exception:
                    pass
                if bridge_poly.is_empty:
                    break
            if (bridge_poly.is_empty
                    or bridge_poly.geom_type
                    not in ("Polygon", "MultiPolygon")):
                continue
            if not bridge_poly.is_valid:
                bridge_poly = bridge_poly.buffer(0)
            if bridge_poly.geom_type == "MultiPolygon":
                parts = sorted(bridge_poly.geoms,
                               key=lambda g: -g.area)
                bridge_poly = parts[0] if parts else None
            if (bridge_poly is None
                    or bridge_poly.is_empty
                    or bridge_poly.geom_type != "Polygon"
                    or bridge_poly.area < 100.0):
                continue
            # Subtract the boundary ribbon so the bridge starts at
            # the ribbon's INNER edge instead of overlapping the
            # ribbon's inner half.
            if (ribbon_union is not None
                    and not ribbon_union.is_empty):
                try:
                    bridge_poly = bridge_poly.difference(ribbon_union)
                except Exception:
                    pass
                if bridge_poly.is_empty:
                    continue
                if not bridge_poly.is_valid:
                    bridge_poly = bridge_poly.buffer(0)
                if bridge_poly.geom_type == "MultiPolygon":
                    parts = sorted(bridge_poly.geoms,
                                   key=lambda g: -g.area)
                    bridge_poly = parts[0] if parts else None
                if (bridge_poly is None
                        or bridge_poly.is_empty
                        or bridge_poly.geom_type != "Polygon"
                        or bridge_poly.area < 100.0):
                    continue
            # If the bridge polygon overlaps any sloping rect, the
            # bridge run extended too close to pavement despite
            # pre-filtering — trim against the rect union with a
            # 0.1 m safety buffer (rather than risk creating
            # vertices on a sloping rect's edge interior).
            sloping_rect_polys: List[Polygon] = [
                s.polygon for s in layout.shapes
                if s.role in (
                    ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                    ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                    ROLE_CROSS_CONNECTOR)
                and s.polygon is not None
                and not s.polygon.is_empty]
            overlaps_rect = False
            for r in sloping_rect_polys:
                try:
                    if (bridge_poly.intersection(r).area > 1.0):
                        overlaps_rect = True
                        break
                except Exception:
                    continue
            if overlaps_rect:
                # Trim the bridge against the sloping-rect union
                # using buffered-shrink to avoid creating
                # edge-interior vertices, then snap any near-corner
                # vertex to its nearest sloping-rect corner.
                try:
                    rect_union = unary_union(sloping_rect_polys)
                    # Buffer the rect union by > the
                    # ``EDGE_PROX_M`` test tolerance (0.5 m) so any
                    # vertex from the difference operation lands
                    # outside that proximity band.
                    bridge_poly = bridge_poly.difference(
                        rect_union.buffer(1.0))
                except Exception:
                    bridge_poly = None
                if (bridge_poly is None
                        or bridge_poly.is_empty):
                    continue
                if bridge_poly.geom_type == "MultiPolygon":
                    parts = sorted(bridge_poly.geoms,
                                   key=lambda g: -g.area)
                    bridge_poly = parts[0] if parts else None
                if (bridge_poly is None
                        or bridge_poly.geom_type != "Polygon"
                        or bridge_poly.area < 100.0):
                    continue
                try:
                    bridge_poly = (
                        _snap_polygon_vertices_to_rect_corners(
                            bridge_poly,
                            sloping_rect_polys,
                            snap_tol_m=5.0))
                except Exception:
                    pass
                if (bridge_poly is None
                        or bridge_poly.is_empty
                        or bridge_poly.geom_type != "Polygon"
                        or bridge_poly.area < 100.0):
                    continue
            # Per-vertex altitudes.  Per user 2026-04-28:
            # the bridge is meant to FILL the gap between the
            # boundary (at clamped altitude) and the nearest
            # pavement (at the pavement's emitted altitude).  So
            # each vertex gets a distance-weighted blend of those
            # two values:
            #
            #   alt(V) = (w_o · clamped_at_outer + w_p · pav_alt)
            #            / (w_o + w_p)
            #
            # with weights w_o = 1 / max(d_outer, ε),
            # w_p = 1 / max(d_pav, ε) — i.e. inverse-distance
            # interpolation.  Vertices on the outer edge land at
            # clamped; vertices touching pavement land at the
            # pavement altitude; interior vertices smoothly
            # interpolate.  This eliminates the previous
            # "DEM-everywhere" inner edge that sat 12–35 m below
            # the surrounding pavement at CYXY and was the cause
            # of the persistent valley the user reported.
            outer_set = set((round(x, 1), round(y, 1))
                            for x, y in outer_pts)
            new_coords = list(bridge_poly.exterior.coords)
            alts: List[float] = []
            EPS_M = 0.5
            for cx, cy in new_coords:
                # Nearest clamped (outer-edge) point and its alt.
                best_o_alt: Optional[float] = None
                best_o_d = float('inf')
                for idx in run:
                    ox, oy, ca, _da = per_vert[idx]
                    if math.isnan(ca):
                        continue
                    d = math.hypot(cx - ox, cy - oy)
                    if d < best_o_d:
                        best_o_d = d
                        best_o_alt = ca
                # Nearest pavement edge point and its alt.
                pav_hit = _nearest_pav_alt(cx, cy)
                key = (round(cx, 1), round(cy, 1))
                # Quick exits: vertex sits exactly on outer edge or
                # on a pavement edge.
                if key in outer_set and best_o_alt is not None:
                    alts.append(round(float(best_o_alt), 1))
                    continue
                if pav_hit is not None and pav_hit[1] < EPS_M:
                    alts.append(round(float(pav_hit[0]), 1))
                    continue
                # Distance-weighted blend.
                if best_o_alt is None and pav_hit is None:
                    # No reference — fall back to DEM, then 0.
                    d = _dem_alt(cx, cy) or _clamped_alt(cx, cy)
                    alts.append(round(float(d or 0.0), 1))
                    continue
                if pav_hit is None:
                    alts.append(round(float(best_o_alt), 1))
                    continue
                if best_o_alt is None:
                    alts.append(round(float(pav_hit[0]), 1))
                    continue
                d_o = max(best_o_d, EPS_M)
                d_p = max(pav_hit[1], EPS_M)
                w_o = 1.0 / d_o
                w_p = 1.0 / d_p
                blended = ((w_o * best_o_alt + w_p * pav_hit[0])
                            / (w_o + w_p))
                alts.append(round(float(blended), 1))
            layout.shapes.append(BuiltShape(
                polygon=bridge_poly,
                role=ROLE_BOUNDARY,
                ref="boundary_dem_bridge",
                node_altitudes=alts))
            n_emitted += 1
    return n_emitted


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
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] terminal({s.ref or '?'}) "
                    f"altitude {old_alt} → {new_alt} m "
                    f"(median of {len(closest)} closest hard "
                    f"anchors within {sample_radius_m:.0f} m; "
                    f"replaces DEM-median).\n")
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
        layout: "PavementLayout") -> int:
    """For every junction polygon, snap any vertex whose bucket
    coincides with a runway / sloping-rect corner to that rect's
    corresponding altitude tag value (``altitude_high`` for HIGH
    corners 0,3; ``altitude_low`` for LOW corners 1,2; ``altitude``
    for flat shapes).

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
        elif s.altitude is not None:
            # Flat shapes (terminal pads, pre-elevation rects).
            # Any number of vertices; all share a single altitude.
            for (cx, cy) in coords:
                b = _corner_elevation_bucket(cx, cy)
                rwy_corner_alt.setdefault(b, float(s.altitude))
    if not rwy_corner_alt:
        return 0
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
            b = _corner_elevation_bucket(cx, cy)
            target_e = rwy_corner_alt.get(b)
            if target_e is None:
                continue
            if abs(s.node_altitudes[i] - target_e) < 0.05:
                continue
            s.node_altitudes[i] = round(target_e, 1)
            n_changed += 1
        # Maintain closed-ring invariant: last == first.
        if (s.node_altitudes
                and len(s.node_altitudes) >= 2
                and s.node_altitudes[0] != s.node_altitudes[-1]):
            s.node_altitudes[-1] = s.node_altitudes[0]
    return n_changed


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


def _decompose_polygon_with_holes(polygon: Polygon,
                                  min_area_m2: float = 50.0,
                                  max_depth: int = 8
                                  ) -> List[Polygon]:
    """Return a list of simple (no-hole) polygons that tile the
    same area as ``polygon``.  Cuts through each hole's centroid
    along a direction chosen to follow the hole's natural axes (or
    the centroid-spread of multiple holes), then recurses on each
    side.

    Per user 2026-04-28: previous implementation always cut
    horizontally, which produced thin (~2–12 m) strips between
    parallel taxi-rect holes whose y-centroids happened to be
    close.  Those strips were arbitrary geometric artifacts —
    their cut edges didn't follow any real pavement feature —
    and required a downstream merge band-aid to suppress.  The
    new policy:

      * Multiple holes ⇒ cut PERPENDICULAR to the centroid-spread
        direction.  When holes are spread along x (typical
        parallel-taxiway apron), the cut runs vertically BETWEEN
        the holes rather than horizontally between two close
        cut lines.
      * Single hole ⇒ cut along the hole's MRR long axis.  The
        cut is collinear with the rect's natural orientation and
        produces two pieces straddling the rect rather than
        slicing it across.

    Both rules align cut lines with real geometric features
    (rect axes, hole-cluster alignment) instead of an arbitrary
    horizontal direction.  The thin-piece merge fallback is kept
    at a small threshold purely for floating-point sliver clean-up.
    """
    from shapely.ops import split as _shp_split
    from shapely.geometry import LineString as _LS

    if (polygon.is_empty or polygon.geom_type != "Polygon"):
        return []
    if not polygon.interiors:
        return [polygon]
    if max_depth <= 0:
        # Recursion guard — emit the exterior with holes dropped
        # rather than retry forever.  Should never trigger for
        # realistic airport geometry.
        return [Polygon(polygon.exterior.coords)]
    # Pick the largest remaining hole — it'll be on one side of
    # the cut after slicing.
    interiors = list(polygon.interiors)
    interiors.sort(key=lambda h: -Polygon(h).area)
    hole = interiors[0]
    cx = float(hole.centroid.x)
    cy = float(hole.centroid.y)
    minx, miny, maxx, maxy = polygon.bounds

    # Decide cut direction.  The cut is a line through the hole
    # centroid; ``angle_rad`` is the angle of the cut LINE (not
    # the perpendicular), measured from +x.  Default = horizontal.
    angle_rad = 0.0
    if len(interiors) >= 2:
        # Multi-hole: cut perpendicular to the centroid-spread
        # axis so the cut SEPARATES the holes rather than slicing
        # parallel to their alignment.  E.g. holes spread along
        # +x → cut vertically (angle = π/2) between them.
        cs = [h.centroid for h in interiors]
        x_spread = max(c.x for c in cs) - min(c.x for c in cs)
        y_spread = max(c.y for c in cs) - min(c.y for c in cs)
        # We want the cut LINE to be perpendicular to the spread
        # direction.  spread_along_x ⇒ cut vertical (π/2);
        # spread_along_y ⇒ cut horizontal (0).
        if x_spread > y_spread:
            angle_rad = math.pi / 2.0  # vertical
        else:
            angle_rad = 0.0            # horizontal
    else:
        # Single hole: cut along its MRR long-axis direction.
        try:
            mrr = hole.minimum_rotated_rectangle
            if (mrr is not None and not mrr.is_empty
                    and mrr.geom_type == "Polygon"):
                mc = list(mrr.exterior.coords)
                if len(mc) >= 5:
                    sides = []
                    for i in range(4):
                        ax, ay = mc[i]
                        bx, by = mc[i + 1]
                        sides.append(
                            (math.hypot(bx - ax, by - ay),
                             math.atan2(by - ay, bx - ax)))
                    sides.sort(reverse=True)
                    angle_rad = sides[0][1]
        except Exception:
            angle_rad = 0.0

    # Build a cut line through (cx, cy) at angle_rad, extended well
    # past the polygon bounds on both sides.
    span = max(maxx - minx, maxy - miny) + 2.0
    dx = math.cos(angle_rad)
    dy = math.sin(angle_rad)
    cut = _LS([(cx - span * dx, cy - span * dy),
               (cx + span * dx, cy + span * dy)])
    try:
        result = _shp_split(polygon, cut)
    except Exception:
        # Fallback: emit the exterior with holes dropped.  Should
        # not occur for valid simple geometries.
        return [Polygon(polygon.exterior.coords)]
    pieces: List[Polygon] = []
    geoms = (list(getattr(result, "geoms", []))
             if result.geom_type != "Polygon" else [result])
    for g in geoms:
        if g.geom_type != "Polygon" or g.is_empty:
            continue
        if g.area < min_area_m2:
            continue
        pieces.extend(_decompose_polygon_with_holes(
            g, min_area_m2=min_area_m2, max_depth=max_depth - 1))
    # Sliver clean-up: smart-cut alignment eliminates the most
    # egregious wide-band strips (5 m × 67 m, 10 m × 119 m) that
    # appeared with horizontal-only cuts, but recursive splits can
    # still leave narrow corner slivers when a hole's MRR axis
    # nearly parallels a polygon edge.  Merge anything thinner
    # than 12 m into its largest-shared-boundary neighbour so the
    # apron stays one continuous polygon and we don't get
    # cliff-rendering on thin corners.
    MIN_PIECE_THICKNESS_M = 12.0
    pieces = _merge_thin_decomposed_pieces(
        pieces, min_thickness_m=MIN_PIECE_THICKNESS_M)
    return pieces


def _polygon_min_thickness(poly: "Polygon") -> float:
    """Approximate minimum thickness of a polygon: half-width of
    the rotated minimum bounding rectangle.  Fast computation: try
    the polygon's minimum-rotated-rectangle and return the shorter
    side length."""
    try:
        mrr = poly.minimum_rotated_rectangle
        if mrr.is_empty or mrr.geom_type != "Polygon":
            return 0.0
        coords = list(mrr.exterior.coords)
        if len(coords) < 5:
            return 0.0
        sides = []
        for i in range(4):
            ax, ay = coords[i]
            bx, by = coords[i + 1]
            sides.append(math.hypot(bx - ax, by - ay))
        return min(sides)
    except Exception:
        return 0.0


def _merge_thin_decomposed_pieces(
        pieces: "List[Polygon]",
        min_thickness_m: float = 4.0,
        max_iters: int = 50,
        ) -> "List[Polygon]":
    """Merge any piece in ``pieces`` whose minimum-rotated-rectangle
    thickness is less than ``min_thickness_m`` into the neighbouring
    piece sharing the longest boundary.  Used by
    ``_decompose_polygon_with_holes`` to suppress 2 m-thick horizontal
    strips that arise when multiple holes have close y-centroids.
    Returns a possibly-shorter list with thin strips absorbed.
    """
    if not pieces:
        return pieces
    work = list(pieces)
    for _it in range(max_iters):
        # Find the thinnest piece below threshold.
        thin_idx: Optional[int] = None
        thin_thick = float('inf')
        for i, p in enumerate(work):
            if p is None or p.is_empty:
                continue
            t = _polygon_min_thickness(p)
            if t < min_thickness_m and t < thin_thick:
                thin_idx = i
                thin_thick = t
        if thin_idx is None:
            break
        thin = work[thin_idx]
        # Find the neighbour with the longest shared boundary.
        best_j: Optional[int] = None
        best_share = 0.0
        thin_boundary = thin.boundary
        for j, p in enumerate(work):
            if j == thin_idx or p is None or p.is_empty:
                continue
            try:
                shared = thin_boundary.intersection(
                    p.boundary).length
            except Exception:
                shared = 0.0
            if shared > best_share:
                best_share = shared
                best_j = j
        if best_j is None or best_share <= 0.0:
            # No neighbour to merge with; drop the thin piece so
            # it doesn't render as a cliff.
            work[thin_idx] = None
            continue
        try:
            merged = unary_union([thin, work[best_j]])
            if merged.is_empty:
                work[thin_idx] = None
                continue
            if merged.geom_type == "MultiPolygon":
                # Pick the largest piece — the union didn't fully
                # bridge.  Drop the thin one.
                work[thin_idx] = None
                continue
            if merged.geom_type != "Polygon":
                work[thin_idx] = None
                continue
            work[best_j] = merged
            work[thin_idx] = None
        except Exception:
            work[thin_idx] = None
    return [p for p in work if p is not None and not p.is_empty]


# ── Hole splicing ────────────────────────────────────────────────
#
# Defensive fallback for any polygon with holes that slips through
# decomposition (e.g. shapely.ops.split failed).  Splices each
# hole into the exterior via a zero-width bridge so ear-clipping
# can operate on a single ring; sliver triangles along the bridge
# are filtered downstream by the area threshold.


def _splice_holes(polygon: Polygon) -> List[Tuple[float, float]]:
    """Return the vertex list of the spliced single-ring polygon
    (without closing repeat).  Holes are inserted one at a time by
    finding the closest exterior vertex / hole vertex pair and
    threading the hole into the exterior at that bridge.
    """
    ext = list(polygon.exterior.coords)
    if ext and ext[0] == ext[-1]:
        ext = ext[:-1]
    holes_list: List[List[Tuple[float, float]]] = []
    for h in polygon.interiors:
        h_coords = list(h.coords)
        if h_coords and h_coords[0] == h_coords[-1]:
            h_coords = h_coords[:-1]
        if len(h_coords) >= 3:
            holes_list.append(h_coords)
    if not holes_list:
        return ext
    # Process holes from largest to smallest so big holes get the
    # "best" bridge slots; small holes thread into the still-clean
    # remainder.
    holes_list.sort(key=lambda h: -_polygon_area(h))
    ring = list(ext)
    for hole in holes_list:
        ring = _splice_one_hole(ring, hole)
    return ring


def _polygon_area(coords: Sequence[Tuple[float, float]]) -> float:
    s = 0.0
    n = len(coords)
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) * 0.5


def _splice_one_hole(ring: List[Tuple[float, float]],
                     hole: List[Tuple[float, float]]
                     ) -> List[Tuple[float, float]]:
    """Find the closest (ring_vertex, hole_vertex) pair and splice
    the hole into the ring at that bridge.  Hole is walked in its
    native (CW relative to a CCW exterior) direction so the spliced
    ring stays simple."""
    best_i = best_j = 0
    best_d2 = float("inf")
    for i, (rx, ry) in enumerate(ring):
        for j, (hx, hy) in enumerate(hole):
            d2 = (rx - hx) * (rx - hx) + (ry - hy) * (ry - hy)
            if d2 < best_d2:
                best_d2 = d2
                best_i, best_j = i, j
    # Spliced ring:
    #   ring[0..best_i] + hole[best_j..end] + hole[0..best_j]
    #   + ring[best_i..end]
    # The boundary touches ring[best_i] and hole[best_j] twice —
    # this is the bridge corridor.
    spliced: List[Tuple[float, float]] = []
    spliced.extend(ring[: best_i + 1])
    spliced.extend(hole[best_j:])
    spliced.extend(hole[: best_j + 1])
    spliced.extend(ring[best_i:])
    return spliced


# ── Ear-clip triangulation ───────────────────────────────────────
#
# Pure-Python ear-clipping for simple polygons (no holes).  Returns
# index triples (i,j,k) into the input vertex list.  For an N-vertex
# simple polygon, exactly N-2 triangles are produced — the proven
# minimum count when no Steiner points are added.


MAX_BOUNDARY_EDGE_M = 30.0  # max length of any junction-polygon
                             # boundary segment.  Long segments
                             # (e.g. cuts from hole-decomposition)
                             # leave Triangle4XP's quality refinement
                             # unable to fit small triangles between
                             # the cut endpoints; interior Steiners
                             # span the full distance and produce
                             # steep local gradients.  Densifying
                             # the boundary with linearly-interpolated
                             # midpoints gives Triangle4XP closer
                             # boundary anchors to connect to.


SHARED_NEIGHBOUR_EDGE_TOL_M = 5.0  # max distance from a midpoint
                                    # to a neighbour-shape edge to
                                    # treat the midpoint as on a
                                    # shared boundary; matches
                                    # check_grade's edge-step
                                    # search radius


def _densify_long_boundary_edges(
    ring: List[Tuple[float, float]],
    vert_elev: List[float],
    neighbour_edges: List[Tuple[float, float, float, float,
                                 float, float]],
) -> Tuple[List[Tuple[float, float]], List[float]]:
    """Insert interpolated midpoints along long ring segments.

    For each candidate midpoint, check distance to the nearest
    rect/runway/terminal edge.  If within
    ``SHARED_NEIGHBOUR_EDGE_TOL_M``, the midpoint sits on a
    shared boundary; use the NEIGHBOUR'S edge-interpolated
    elevation (matches whatever the neighbour renders at that
    point — including segmented-runway piecewise profiles).
    Otherwise use linear interpolation between the segment's
    endpoint elevations.

    This guarantees no shared-boundary step is introduced and
    handles both the simple "junction-on-rect-edge" case (linear
    interp matches) and the "junction-on-segmented-runway-edge"
    case (use the runway segment's interp).
    """
    n = len(ring)
    if n < 3 or len(vert_elev) != n:
        return ring, vert_elev

    def _interp_at(mx: float, my: float, fallback: float) -> float:
        """Return neighbour-edge interpolation at (mx, my) if a
        neighbour edge passes within tolerance; else fallback."""
        best_d2 = SHARED_NEIGHBOUR_EDGE_TOL_M * SHARED_NEIGHBOUR_EDGE_TOL_M
        best_e: Optional[float] = None
        for ax, ay, bx, by, ea, eb in neighbour_edges:
            dx = bx - ax
            dy = by - ay
            seg2 = dx * dx + dy * dy
            if seg2 < 0.04:
                continue
            t = ((mx - ax) * dx + (my - ay) * dy) / seg2
            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0
            cx = ax + t * dx
            cy = ay + t * dy
            d2 = (mx - cx) * (mx - cx) + (my - cy) * (my - cy)
            if d2 < best_d2:
                best_d2 = d2
                best_e = ea + t * (eb - ea)
        return best_e if best_e is not None else fallback

    def _point_on_neighbour(mx: float, my: float) -> bool:
        """True if point (mx, my) lies within SHARED_VERTEX_TOL_M of
        some neighbour edge's interior.

        A densification midpoint that lands this close to a
        rect/runway/terminal edge will be snapped onto that edge by
        X-Plane's vector_map encroachment check, creating a T-
        junction that breaks the neighbour rect's 4-corner slope
        rendering.  Skip these midpoints — the junction edge stays
        un-densified at that point, which is the right call because
        Triangle4XP doesn't need an interior-anchor near a shared
        boundary anyway (the neighbour rect's own subdivision
        provides the anchors).
        """
        tol2 = SHARED_VERTEX_TOL_M * SHARED_VERTEX_TOL_M
        for nax, nay, nbx, nby, _, _ in neighbour_edges:
            ndx = nbx - nax
            ndy = nby - nay
            seg2 = ndx * ndx + ndy * ndy
            if seg2 < 0.04:
                continue
            t = ((mx - nax) * ndx + (my - nay) * ndy) / seg2
            if t < 0.0 or t > 1.0:
                continue
            cx = nax + t * ndx
            cy = nay + t * ndy
            d2 = (mx - cx) * (mx - cx) + (my - cy) * (my - cy)
            if d2 < tol2:
                return True
        return False

    # Pre-compute the bucket-key of every existing ring vertex so we
    # can reject midpoints that would collide with one via
    # ``to_osm``'s SHARED_VERTEX_TOL_M intern.  A midpoint that
    # collides would emit the same OSM nid as a non-adjacent ring
    # vertex, producing a polygon that visits the same node twice
    # — duplicate-consecutive or self-intersection (figure-8) at
    # OSM-write time.  Either crashes X-Plane's mesh builder.
    existing_buckets = {
        _corner_elevation_bucket(x, y) for (x, y) in ring}
    new_ring: List[Tuple[float, float]] = []
    new_elev: List[float] = []
    for i in range(n):
        a = ring[i]
        b = ring[(i + 1) % n]
        ea = vert_elev[i]
        eb = vert_elev[(i + 1) % n]
        new_ring.append(a)
        new_elev.append(ea)
        d = math.hypot(b[0] - a[0], b[1] - a[1])
        if d <= MAX_BOUNDARY_EDGE_M:
            continue
        n_subs = int(math.ceil(d / MAX_BOUNDARY_EDGE_M))
        for k in range(1, n_subs):
            t = k / n_subs
            mx = a[0] + t * (b[0] - a[0])
            my = a[1] + t * (b[1] - a[1])
            mb = _corner_elevation_bucket(mx, my)
            if mb in existing_buckets:
                continue  # would collide with an existing ring nid
            # Skip midpoints that would land near a neighbour edge
            # — X-Plane's vector_map would snap them onto the
            # neighbour, creating a T-junction that breaks the
            # neighbour rect's 4-corner slope rendering.
            if _point_on_neighbour(mx, my):
                continue
            linear_me = ea + t * (eb - ea)
            me = _interp_at(mx, my, linear_me)
            new_ring.append((mx, my))
            new_elev.append(me)
            existing_buckets.add(mb)
    return new_ring, new_elev


COLINEAR_DROP_M = 3.0  # max perpendicular distance to neighbours
                        # for a non-anchor vertex to be removed.
                        # Bumped from 0.5 m → 3.0 m (user 2026-04-25)
                        # to thin apt.dat boundary curves more
                        # aggressively — eliminates most ear-clip
                        # sliver triangles whose 3 anchored vertices
                        # are nearly colinear.  Preserves rect /
                        # runway / terminal corners (sharp 90° turns
                        # are well above this threshold by definition)
                        # and any vertex shared between junctions.


SPIKE_VERTEX_TOL_M = 0.005  # max perpendicular distance from a
                             # ring vertex to a non-adjacent edge
                             # of the same polygon for the vertex
                             # to count as a "spike" — the ring
                             # ventured out and returned to (or
                             # very near) itself.  Sub-mm spikes
                             # are valid in shapely's eyes but
                             # become hard self-intersections
                             # after .11f OSM-format truncation,
                             # which would crash X-Plane.


def _drop_spike_vertices(
    ring: List[Tuple[float, float]],
) -> List[Tuple[float, float]]:
    """Drop any ring vertex that lies within ``SPIKE_VERTEX_TOL_M``
    of a NON-adjacent edge of the same ring.

    Such a vertex represents a degenerate "stick-out-and-return"
    in the polygon boundary — the ring went away from a straight
    section and returned right onto it, leaving a near-zero-area
    lobe that's a self-touch / self-intersection in any reasonable
    coordinate precision.

    Iterates to a fixed point (dropping one spike can expose
    another).  Capped at 8 passes against pathological inputs.
    """
    if len(ring) < 4:
        return ring
    for _ in range(8):
        n = len(ring)
        if n < 4:
            break
        keep = [True] * n
        for i in range(n):
            vx, vy = ring[i]
            for j in range(n):
                if abs(i - j) <= 1 or (i == 0 and j == n - 1) or (j == 0 and i == n - 1):
                    continue
                ax, ay = ring[j]
                bx, by = ring[(j + 1) % n]
                dx = bx - ax
                dy = by - ay
                seg2 = dx * dx + dy * dy
                if seg2 < 1e-6:
                    continue
                t = ((vx - ax) * dx + (vy - ay) * dy) / seg2
                if t < 0.0 or t > 1.0:
                    continue
                cx = ax + t * dx
                cy = ay + t * dy
                d2 = (vx - cx) * (vx - cx) + (vy - cy) * (vy - cy)
                if d2 < SPIKE_VERTEX_TOL_M * SPIKE_VERTEX_TOL_M:
                    keep[i] = False
                    break
        new_ring = [r for r, k in zip(ring, keep) if k]
        if len(new_ring) == n:
            break
        ring = new_ring
    return ring


SLIVER_ANGLE_THRESHOLD_DEG = 2.0  # interior angles below this count
                                   # as needle-tip slivers.  Source:
                                   # residue construction can leave
                                   # thin wedges where rect/terminal
                                   # edges meet the apt.dat boundary at
                                   # near-collinear angles.  The
                                   # polygon is shapely-valid but a
                                   # sub-2° corner forces Triangle4XP
                                   # to emit a near-degenerate triangle
                                   # there — crashes X-Plane's mesh
                                   # builder.  Caught at junction-
                                   # emission time by _drop_sliver_corners
                                   # (drops just the tip vertex,
                                   # preserves the rest of the polygon)
                                   # and again by a to_osm safety net
                                   # (drops the whole shape if any
                                   # slipped through, e.g. via buffer(0)
                                   # repair or post-densification).


def _drop_sliver_corners(
    ring: List[Tuple[float, float]],
) -> List[Tuple[float, float]]:
    """Drop ring vertices whose interior angle is below
    ``SLIVER_ANGLE_THRESHOLD_DEG``.

    A sliver corner is a needle-tip vertex: the polygon comes in
    along one edge, makes a near-180° fold, and goes back out
    almost on top of the incoming edge, leaving a thin wedge.
    Source: residue construction (apt.dat pav − rects − terminals)
    leaves wedges where two rect/terminal edges meet the boundary at
    nearly-collinear angles.  Shapely calls these polygons valid
    (the two long edges are parallel-but-not-equal, no self-
    intersection), but the polygon's needle-tip corner forces
    Triangle4XP downstream to emit at least one corner triangle
    with that interior angle — for sub-2° tips that triangle is
    near-degenerate (NaN normal) and crashes X-Plane's mesh builder.

    Dropping the tip vertex collapses the wedge into a single edge
    between the two flanking vertices.  Coverage cost: the area
    between the tip and the truncation chord (typically < 50 m²).

    Iterates to a fixed point — dropping one tip can expose
    another.  Capped at 8 passes.
    """
    if len(ring) < 4:
        return ring
    cos_thresh = math.cos(math.radians(SLIVER_ANGLE_THRESHOLD_DEG))
    for _ in range(8):
        n = len(ring)
        if n < 4:
            break
        keep = [True] * n
        for i in range(n):
            ax, ay = ring[(i - 1) % n]
            bx, by = ring[i]
            cx, cy = ring[(i + 1) % n]
            v1x, v1y = ax - bx, ay - by
            v2x, v2y = cx - bx, cy - by
            n1 = math.hypot(v1x, v1y)
            n2 = math.hypot(v2x, v2y)
            if n1 < 1e-9 or n2 < 1e-9:
                continue
            cos = (v1x * v2x + v1y * v2y) / (n1 * n2)
            if cos > cos_thresh:
                # Angle = acos(cos) is below threshold.
                keep[i] = False
        new_ring = [r for r, k in zip(ring, keep) if k]
        if len(new_ring) == n:
            break
        ring = new_ring
    return ring


def _drop_colinear_boundary_vertices(
    ring: List[Tuple[float, float]],
    corner_elev: Dict[Tuple[int, int], float],
    shared_junction_buckets: set,
) -> List[Tuple[float, float]]:
    """Remove vertices whose perpendicular distance to the line
    through their immediate neighbours is below COLINEAR_DROP_M,
    EXCEPT vertices that are shared corners (rect/runway/terminal
    or cross-junction shared buckets) — those carry topological
    meaning and must not be dropped.

    Eliminates the "ear-clip can only emit a sliver here" geometry
    that produces visible step artefacts on long thin apron strips.
    """
    if len(ring) < 4:
        return ring
    # Iterate to fixed point: dropping one vertex may make a
    # neighbour droppable too.
    for _ in range(8):
        n = len(ring)
        if n < 4:
            break
        keep = [True] * n
        for i in range(n):
            bucket = _corner_elevation_bucket(*ring[i])
            if bucket in corner_elev or bucket in shared_junction_buckets:
                continue  # anchor — keep no matter what
            ax, ay = ring[(i - 1) % n]
            bx, by = ring[(i + 1) % n]
            cx, cy = ring[i]
            # Perpendicular distance from C to line AB.
            dx = bx - ax
            dy = by - ay
            seg_len = math.hypot(dx, dy)
            if seg_len < 0.1:
                continue
            # Cross-product / line-length = perpendicular distance.
            perp = abs((cx - ax) * dy - (cy - ay) * dx) / seg_len
            if perp < COLINEAR_DROP_M:
                keep[i] = False
        new_ring = [c for c, k in zip(ring, keep) if k]
        if len(new_ring) == n:
            break
        ring = new_ring
    return ring


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


def _smooth_polygon_grid(
        polygon: Polygon,
        hard_anchors: List[Tuple[float, float, float]],
        dem,
        tile_lat: int,
        tile_lon: int,
        layout_anchor: Tuple[float, float],
        grid_step_m: float = ELEVATION_GRID_STEP_M,
        max_iters: int = ELEVATION_SMOOTH_MAX_ITERS,
        convergence_tol_m: float = ELEVATION_SMOOTH_CONVERGE_M):
    """Smooth a 2D elevation field within ``polygon``.

    ``hard_anchors`` is a list of ``(x, y, elev)`` triples whose
    elevation must be preserved exactly (rect / runway / terminal
    corners; cross-junction shared vertices already locked from a
    previous outer iteration).

    Returns ``(sampler, sample)`` where ``sampler(x, y)`` returns
    the bilinearly-interpolated elevation, or ``None`` if the
    polygon was too small for grid construction.  When no anchors
    are reachable for the polygon, falls back to the global
    elevation graph / DEM at each query point — same path
    ``_vertex_elev_anchored`` step 5 takes today.
    """
    import numpy as np

    # Degenerate polygon — return a graph/DEM sampler.
    if polygon is None or polygon.is_empty:
        return None
    minx, miny, maxx, maxy = polygon.bounds
    span_x = maxx - minx
    span_y = maxy - miny
    if span_x <= 0 or span_y <= 0:
        return None

    # Build grid bbox with a 1-cell margin so boundary vertices land
    # comfortably inside the active region.
    margin = grid_step_m
    minx -= margin
    miny -= margin
    maxx += margin
    maxy += margin
    nx = max(2, int(math.ceil((maxx - minx) / grid_step_m)) + 1)
    ny = max(2, int(math.ceil((maxy - miny) / grid_step_m)) + 1)
    # Cap grid size — should never trigger for realistic airport
    # polygons but prevents pathological all-airport polygons from
    # eating all RAM.
    if nx * ny > 250_000:
        return None

    # Cell-center coords.
    xs = minx + np.arange(nx) * grid_step_m
    ys = miny + np.arange(ny) * grid_step_m

    # ── INSIDE mask ─────────────────────────────────────────────
    # A cell is INSIDE if its center sits inside the polygon OR
    # within ``grid_step_m`` of the polygon boundary (so cells
    # straddling the boundary stay active).  Bilinear sampling at
    # boundary vertices needs the surrounding cells active.
    from shapely.prepared import prep
    prepped = prep(polygon)
    boundary = polygon.boundary
    inside = np.zeros((nx, ny), dtype=bool)
    for i in range(nx):
        for j in range(ny):
            pt = Point(float(xs[i]), float(ys[j]))
            if prepped.contains(pt):
                inside[i, j] = True
            else:
                try:
                    if boundary.distance(pt) < grid_step_m:
                        inside[i, j] = True
                except Exception:
                    pass

    if not inside.any():
        return None

    # ── Pin anchors ─────────────────────────────────────────────
    # Each hard anchor pins its NEAREST INSIDE cell to its elevation.
    # Multiple anchors hitting the same cell are averaged (rare —
    # implies two rect corners within one grid cell, which means
    # they should be the same shared-vertex bucket already).
    pinned = np.zeros((nx, ny), dtype=bool)
    pinned_elev = np.zeros((nx, ny), dtype=float)
    pinned_count = np.zeros((nx, ny), dtype=int)
    eff_anchors: List[Tuple[float, float, float]] = []
    for (ax, ay, az) in hard_anchors:
        i = int(round((ax - minx) / grid_step_m))
        j = int(round((ay - miny) / grid_step_m))
        if 0 <= i < nx and 0 <= j < ny:
            # Snap to the nearest INSIDE cell within a 1-cell radius
            # (anchors at the boundary may map to an INACTIVE cell).
            if not inside[i, j]:
                snapped = False
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        ii, jj = i + di, j + dj
                        if (0 <= ii < nx and 0 <= jj < ny
                                and inside[ii, jj]):
                            i, j = ii, jj
                            snapped = True
                            break
                    if snapped:
                        break
                if not snapped:
                    continue
            if pinned[i, j]:
                pinned_elev[i, j] = (
                    pinned_elev[i, j] * pinned_count[i, j] + az
                ) / (pinned_count[i, j] + 1)
                pinned_count[i, j] += 1
            else:
                pinned[i, j] = True
                pinned_elev[i, j] = az
                pinned_count[i, j] = 1
            eff_anchors.append((float(xs[i]), float(ys[j]),
                                float(pinned_elev[i, j])))

    # ── Per-cell feasibility band ───────────────────────────────
    INF = float("inf")
    lo = np.full((nx, ny), -INF, dtype=float)
    hi = np.full((nx, ny), INF, dtype=float)
    if eff_anchors:
        XX, YY = np.meshgrid(xs, ys, indexing="ij")
        for (ax, ay, az) in eff_anchors:
            dist = np.hypot(XX - ax, YY - ay)
            band = dist * TAXI_MAX_GRADE
            lo = np.maximum(lo, az - band)
            hi = np.minimum(hi, az + band)

    # ── Initialize cells ────────────────────────────────────────
    # Pinned: their anchor elevation.  Free: prefer the elevation
    # graph's already-grade-compliant network value (smoothed at
    # 1.5 % across the centerline graph), falling back to DEM
    # only when the cell is too far from the network for graph
    # sampling to return.  Using the graph first preserves grade
    # compliance for free cells beyond any anchor cone's reach —
    # raw DEM-init lets real-terrain bumps leak into the polygon's
    # interior, producing wild boundary samples.
    elev = pinned_elev.copy()
    cos0 = math.cos(math.radians(layout_anchor[0]))
    lat0, lon0 = layout_anchor
    for i in range(nx):
        for j in range(ny):
            if not inside[i, j] or pinned[i, j]:
                continue
            d = None
            if dem is not None:
                lat = lat0 + math.degrees(ys[j] / R_EARTH)
                lon = (lon0 + math.degrees(
                    xs[i] / (R_EARTH * cos0)))
                d = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            l = lo[i, j]
            h = hi[i, j]
            if d is None:
                if l == -INF and h == INF:
                    elev[i, j] = 0.0
                elif l == -INF:
                    elev[i, j] = h
                elif h == INF:
                    elev[i, j] = l
                else:
                    elev[i, j] = 0.5 * (l + h)
            else:
                if l > h:
                    elev[i, j] = 0.5 * (l + h)
                else:
                    elev[i, j] = max(l, min(h, d))

    # ── Iterate Laplacian + grade cap ───────────────────────────
    step_cap = grid_step_m * TAXI_MAX_GRADE
    inside_arr = inside
    pinned_arr = pinned
    free_inside = inside_arr & ~pinned_arr
    damping = 0.5

    for _ in range(max_iters):
        # Laplacian — average of 4 cardinal INSIDE neighbours.
        e_pad = np.pad(elev, 1, mode="edge")
        in_pad = np.pad(inside_arr, 1, mode="constant",
                        constant_values=False)
        n_e = e_pad[2:, 1:-1]
        n_w = e_pad[:-2, 1:-1]
        n_n = e_pad[1:-1, 2:]
        n_s = e_pad[1:-1, :-2]
        m_e = in_pad[2:, 1:-1]
        m_w = in_pad[:-2, 1:-1]
        m_n = in_pad[1:-1, 2:]
        m_s = in_pad[1:-1, :-2]
        cnt = (m_e.astype(np.int8)
               + m_w.astype(np.int8)
               + m_n.astype(np.int8)
               + m_s.astype(np.int8))
        sum_n = (n_e * m_e + n_w * m_w + n_n * m_n + n_s * m_s)
        with np.errstate(divide="ignore", invalid="ignore"):
            mean_n = np.where(cnt > 0,
                              sum_n / np.maximum(cnt, 1),
                              elev)
        target = elev + damping * (mean_n - elev)
        new_elev = np.where(free_inside, target, elev)
        # Clip free cells back into feasibility bands.
        new_elev = np.where(
            free_inside,
            np.maximum(np.minimum(new_elev, hi), lo),
            new_elev)

        # Edge grade cap on each axis-aligned cell pair.
        for axis in (0, 1):
            if axis == 0:
                a = new_elev[:-1, :]
                b = new_elev[1:, :]
                ia = inside_arr[:-1, :]
                ib = inside_arr[1:, :]
                pa = pinned_arr[:-1, :]
                pb = pinned_arr[1:, :]
            else:
                a = new_elev[:, :-1]
                b = new_elev[:, 1:]
                ia = inside_arr[:, :-1]
                ib = inside_arr[:, 1:]
                pa = pinned_arr[:, :-1]
                pb = pinned_arr[:, 1:]
            both_inside = ia & ib
            diff = a - b
            need = both_inside & (np.abs(diff) > step_cap)
            sign = np.sign(diff)
            excess = np.maximum(np.abs(diff) - step_cap, 0.0)
            both_free = (~pa) & (~pb)
            a_pinned = pa & (~pb)
            b_pinned = pb & (~pa)
            half = 0.5 * excess * sign
            corr_a_both = -half
            corr_b_both = half
            corr_b_apf = excess * sign
            corr_a_bpf = -excess * sign
            apply_both = need & both_free
            apply_b_apf = need & a_pinned
            apply_a_bpf = need & b_pinned
            if axis == 0:
                new_elev[:-1, :] = np.where(
                    apply_both, a + corr_a_both, new_elev[:-1, :])
                new_elev[1:, :] = np.where(
                    apply_both, b + corr_b_both, new_elev[1:, :])
                new_elev[1:, :] = np.where(
                    apply_b_apf, b + corr_b_apf, new_elev[1:, :])
                new_elev[:-1, :] = np.where(
                    apply_a_bpf, a + corr_a_bpf, new_elev[:-1, :])
            else:
                new_elev[:, :-1] = np.where(
                    apply_both, a + corr_a_both, new_elev[:, :-1])
                new_elev[:, 1:] = np.where(
                    apply_both, b + corr_b_both, new_elev[:, 1:])
                new_elev[:, 1:] = np.where(
                    apply_b_apf, b + corr_b_apf, new_elev[:, 1:])
                new_elev[:, :-1] = np.where(
                    apply_a_bpf, a + corr_a_bpf, new_elev[:, :-1])
            # Re-pin: pinned cells must never have their elevation
            # mutated by the cap pass.
            new_elev = np.where(pinned_arr, pinned_elev, new_elev)

        max_change = float(np.max(np.abs(new_elev - elev)))
        elev = new_elev
        if max_change < convergence_tol_m:
            break

    # ── Bilinear sampler closure ────────────────────────────────
    def sampler(x: float, y: float) -> Optional[float]:
        fi = (x - minx) / grid_step_m
        fj = (y - miny) / grid_step_m
        i0 = int(math.floor(fi))
        j0 = int(math.floor(fj))
        if i0 < 0:
            i0 = 0
        if j0 < 0:
            j0 = 0
        if i0 > nx - 2:
            i0 = nx - 2
        if j0 > ny - 2:
            j0 = ny - 2
        u = fi - i0
        v = fj - j0
        if u < 0:
            u = 0.0
        elif u > 1:
            u = 1.0
        if v < 0:
            v = 0.0
        elif v > 1:
            v = 1.0
        m00 = inside[i0, j0]
        m10 = inside[i0 + 1, j0]
        m01 = inside[i0, j0 + 1]
        m11 = inside[i0 + 1, j0 + 1]
        if m00 and m10 and m01 and m11:
            c00 = elev[i0, j0]
            c10 = elev[i0 + 1, j0]
            c01 = elev[i0, j0 + 1]
            c11 = elev[i0 + 1, j0 + 1]
            e0 = c00 * (1 - u) + c10 * u
            e1 = c01 * (1 - u) + c11 * u
            return float(e0 * (1 - v) + e1 * v)
        # Fallback: nearest INSIDE cell within a 2-cell radius.
        best_d2 = float("inf")
        best_e: Optional[float] = None
        for di in range(-2, 4):
            for dj in range(-2, 4):
                ii = i0 + di
                jj = j0 + dj
                if (0 <= ii < nx and 0 <= jj < ny
                        and inside[ii, jj]):
                    d2 = (xs[ii] - x) ** 2 + (ys[jj] - y) ** 2
                    if d2 < best_d2:
                        best_d2 = d2
                        best_e = float(elev[ii, jj])
        return best_e

    return sampler


# ── Junction triangulation pass ──────────────────────────────────


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


def _triangulate_junctions(
    layout: "PavementLayout",
    dem,
    tile_lat: int,
    tile_lon: int,
    m_to_ll,
) -> int:
    """Replace each junction shape with its ear-clip triangulation,
    setting per-triangle ``node_altitudes`` from the corner-elevation
    map (with neighbour-corner, elevation-graph, and DEM fallbacks).

    Returns the number of triangle shapes produced (informational).
    """
    corner_elev = _corner_elev_map(layout)

    # Cross-junction shared-vertex anchoring.  When two adjacent
    # junction polygons share a boundary point (e.g. either side of
    # a decomposition cut, or both edges of a former hole), they
    # MUST agree on its elevation.  Independent per-junction
    # smoothing can otherwise drift them apart, producing a
    # vertical step at the shared edge.
    #
    # Approach: count how many junction polygons reference each
    # SHARED_VERTEX_TOL_M bucket.  Any bucket touched by ≥ 2
    # junctions is "shared"; we anchor those vertices to a
    # deterministic value (the first junction's graph-sampled or
    # corner-derived elevation) so all junctions read the same
    # value when looking up that bucket.
    junction_bucket_count: Dict[Tuple[int, int], int] = {}
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        try:
            ring = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if ring and ring[0] == ring[-1]:
            ring = ring[:-1]
        seen: set = set()
        for (x, y) in ring:
            key = _corner_elevation_bucket(x, y)
            if key in seen:
                continue
            seen.add(key)
            junction_bucket_count[key] = (
                junction_bucket_count.get(key, 0) + 1)
    shared_junction_buckets = {
        k for k, c in junction_bucket_count.items() if c >= 2}
    # Computed lazily and cached: the first junction to look up a
    # shared-bucket vertex computes its elevation; subsequent
    # junctions read the same value.
    # Two-tier cache so we can distinguish where a shared-bucket
    # value came from: rect-corner / near-corner / edge-interp
    # (HARD — same physical pavement element pinned the value, must
    # be preserved) versus graph or DEM (SOFT — sampled from a 1.5 %
    # network-distance-compliant model, but not 2D-Euclidean
    # compliant against other graph samples that may be far on the
    # network but close in 2D).  Junction-local smoothing must be
    # free to move SOFT values into 2D grade compliance; the cross-
    # junction iteration below averages SOFT shared-bucket values
    # across every junction that touches the bucket so junctions
    # still agree on a single shared elevation.
    shared_junction_elev: Dict[Tuple[int, int], float] = {}      # HARD
    shared_junction_elev_soft: Dict[Tuple[int, int], float] = {}  # SOFT

    # Flat list of (cx, cy, elev) for nearest-corner search beyond
    # the 0.5 m bucket.  Catches boundary vertices that landed close
    # to but not exactly on a neighbour corner (e.g. apt.dat
    # boundary-trace points that abut a terminal pad edge).
    rect_like_roles = {ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL,
                       ROLE_STUB, ROLE_CROSS_CONNECTOR}
    NEAR_CORNER_M = 6.0  # search radius for off-bucket corner matches
    NEAR_EDGE_M = 5.0    # search radius for rect-edge interpolation
    near_corner_list: List[Tuple[float, float, float]] = []
    # Each ``neighbour_edges`` entry is one polygon edge of a non-
    # junction shape, paired with the elevations at its two ends:
    # ``(ax, ay, bx, by, e_a, e_b)``.  Junction vertices that fall
    # close to such an edge are anchored to the linearly-interpolated
    # elevation along the edge — guaranteeing the junction triangle
    # meets a sloped rect's long edge at the same height the rect is
    # rendering at, rather than at the centerline-graph value (which
    # is offset by the rect's half-width).
    neighbour_edges: List[Tuple[float, float, float, float,
                                float, float]] = []
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
        # Per-corner elevations.
        if (s.role in rect_like_roles
                and s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            elevs = [s.altitude_high, s.altitude_low,
                     s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            elevs = [float(s.altitude)] * len(coords)
        else:
            continue
        for (cx, cy), e in zip(coords, elevs):
            near_corner_list.append((cx, cy, float(e)))
        # Edge list (closed ring).
        m = len(coords)
        for i in range(m):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % m]
            ea = float(elevs[i])
            eb = float(elevs[(i + 1) % m])
            neighbour_edges.append((ax, ay, bx, by, ea, eb))

    def _edge_interp_elev(x: float, y: float,
                          max_dist: float = NEAR_EDGE_M
                          ) -> Optional[float]:
        """Return the elevation interpolated along the closest
        non-junction polygon edge within ``max_dist`` of (x, y),
        or None if no edge is within range.

        Projects the query point onto each edge's segment, clamps
        to [0, 1], and linearly interpolates between the edge's
        two endpoint elevations.  This matches the elevation X-Plane
        renders for a sloped rect at any point along its long edge
        — so junction triangles abutting a sloped rect at this
        point will meet it without a step.
        """
        best_e: Optional[float] = None
        best_d2 = max_dist * max_dist
        for ax, ay, bx, by, ea, eb in neighbour_edges:
            dx = bx - ax
            dy = by - ay
            seg_len2 = dx * dx + dy * dy
            if seg_len2 < 0.04:  # < 0.2 m segment, skip
                continue
            t = ((x - ax) * dx + (y - ay) * dy) / seg_len2
            if t < 0.0:
                t = 0.0
            elif t > 1.0:
                t = 1.0
            cx = ax + t * dx
            cy = ay + t * dy
            d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy)
            if d2 < best_d2:
                best_d2 = d2
                best_e = ea + t * (eb - ea)
        return best_e

    def _vertex_elev_anchored(x: float, y: float
                              ) -> Tuple[Optional[float], bool]:
        """Return ``(elev, is_anchor)``.  ``is_anchor`` is True only
        for HARD anchor sources: rect / runway / terminal corner
        elevations (whether matched by exact bucket or off-bucket
        near-corner / edge-interp search), or a value already cached
        from such a source on a prior shared-bucket lookup.  Returns
        is_anchor=False for graph- or DEM-sampled values so the
        boundary smoother and Mode B grid are free to pull them into
        2D Euclidean grade compliance with the real anchors —
        graph-derived values respect 1.5 % over **network distance**
        (centerline graph), which can yield arbitrarily large 2D
        steps between graph nodes that are far on the network but
        close in 2D.  Cross-junction agreement on shared SOFT
        buckets is restored by the cross-junction averaging pass in
        the outer iteration below."""
        bucket = _corner_elevation_bucket(x, y)
        # 1. Bucket lookup — exact shared-vertex match against
        # rect / runway / terminal corners.  HARD.
        e = corner_elev.get(bucket)
        if e is not None:
            return e, True
        # 2a. Cross-junction shared bucket previously locked from a
        # HARD source (rect corner / near-corner / edge-interp).
        e_shared = shared_junction_elev.get(bucket)
        if e_shared is not None:
            return e_shared, True
        # 2b. Cross-junction shared bucket previously sampled from a
        # SOFT source (graph / DEM).  Reuse the value so this call
        # returns deterministically consistent with the first
        # junction that touched the bucket — but mark it FREE so
        # smoothing can still move it.
        e_shared_soft = shared_junction_elev_soft.get(bucket)
        if e_shared_soft is not None:
            return e_shared_soft, False
        # 3. Wider linear search for off-bucket near-corner matches.
        # Treated as HARD: the value comes from a real rect / runway /
        # terminal corner physically nearby.
        best_e: Optional[float] = None
        best_d2 = NEAR_CORNER_M * NEAR_CORNER_M
        for cx, cy, ce in near_corner_list:
            d2 = (cx - x) * (cx - x) + (cy - y) * (cy - y)
            if d2 < best_d2:
                best_d2 = d2
                best_e = ce
        if best_e is not None:
            if bucket in shared_junction_buckets:
                shared_junction_elev[bucket] = best_e
            return best_e, True
        # 4. Rect-edge interpolation — pulls junction vertices
        # pushed 1m off a long edge (or any boundary-trace vertex
        # within NEAR_EDGE_M of a rect/runway/terminal edge) onto
        # the rect's slope at the projected position.  HARD.
        e_edge = _edge_interp_elev(x, y)
        if e_edge is not None:
            if bucket in shared_junction_buckets:
                shared_junction_elev[bucket] = e_edge
            return e_edge, True
        # 5. Free sample from graph / DEM.  SOFT — caller's smoother
        # is free to move this value.  Cache shared-bucket values
        # in the SOFT cache so subsequent junctions see the same
        # initial value.
        e_free: Optional[float] = None
        if dem is not None:
            lat, lon = m_to_ll(x, y)
            e_free = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        if e_free is None:
            return None, False
        if bucket in shared_junction_buckets:
            shared_junction_elev_soft[bucket] = e_free
        return e_free, False

    def _smooth_junction_boundary(
        ring: List[Tuple[float, float]],
        elev: List[float],
        is_anchor: List[bool],
    ) -> List[float]:
        """Bounds-propagation smoothing: each anchor vertex
        constrains every other vertex's elevation to lie within
        ``anchor ± dist × TAXI_MAX_GRADE`` (straight-line distance
        through the junction).  Non-anchor vertices clip to the
        intersection of these bands; anchor vertices stay put.

        Then a Laplacian pass smooths non-anchor elevations toward
        their neighbours' average, re-clipped to the bands so
        anchor compliance is preserved.
        """
        n = len(ring)
        e = list(elev)
        if n < 2:
            return e
        INF = float("inf")
        lo = [-INF] * n
        hi = [INF] * n
        # Anchors lock themselves and contribute bands to others.
        for i in range(n):
            if is_anchor[i]:
                lo[i] = e[i]
                hi[i] = e[i]
        anchor_idx = [i for i in range(n) if is_anchor[i]]
        for i in range(n):
            if is_anchor[i]:
                continue
            xi, yi = ring[i]
            for ai in anchor_idx:
                xa, ya = ring[ai]
                d = math.hypot(xi - xa, yi - ya)
                band = d * TAXI_MAX_GRADE
                lo_i = e[ai] - band
                hi_i = e[ai] + band
                if lo_i > lo[i]:
                    lo[i] = lo_i
                if hi_i < hi[i]:
                    hi[i] = hi_i
        # Initial clip into feasibility intervals.
        for i in range(n):
            if is_anchor[i]:
                continue
            if lo[i] > hi[i]:
                # Conflicting anchors — fall back to midpoint of
                # the conflicting bounds.  Will pull this vertex
                # closer to the closest anchor in practice.
                e[i] = 0.5 * (lo[i] + hi[i])
            else:
                if e[i] < lo[i]:
                    e[i] = lo[i]
                if e[i] > hi[i]:
                    e[i] = hi[i]
        # Laplacian-style smoothing pass between non-anchor
        # neighbours, clipped to bands each iteration so the band
        # constraint stays satisfied.  Convergence in ~20 iters
        # for our junction sizes.
        damping = 0.4
        for _ in range(20):
            new_e = list(e)
            max_change = 0.0
            for i in range(n):
                if is_anchor[i]:
                    continue
                # Mean of all OTHER non-anchor + anchor neighbours,
                # weighted by inverse distance (closer pulls more).
                xi, yi = ring[i]
                num = 0.0
                den = 0.0
                for j in range(n):
                    if j == i:
                        continue
                    xj, yj = ring[j]
                    d = math.hypot(xi - xj, yi - yj)
                    if d < 0.5:
                        continue
                    w = 1.0 / d
                    num += e[j] * w
                    den += w
                if den <= 0:
                    continue
                mean = num / den
                target = e[i] + (mean - e[i]) * damping
                if target < lo[i]:
                    target = lo[i]
                if target > hi[i]:
                    target = hi[i]
                if abs(target - e[i]) > max_change:
                    max_change = abs(target - e[i])
                new_e[i] = target
            e = new_e
            if max_change < 1e-3:
                break
        return e

    # ── Per-junction ring cleanup (independent of elevation) ────
    # Build the cleaned ring + Polygon once per junction; reuse
    # across every cross-junction iteration so we don't redo
    # geometric cleanup at each pass.
    junction_cleaned: List[Tuple["BuiltShape",
                                  List[Tuple[float, float]],
                                  Polygon]] = []
    junction_dropped_shapes: List["BuiltShape"] = []
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        try:
            ring = _splice_holes(shape.polygon)
        except Exception:
            try:
                ring = list(shape.polygon.exterior.coords)
                if ring and ring[0] == ring[-1]:
                    ring = ring[:-1]
            except Exception:
                # Couldn't extract a ring at all — preserve the
                # original shape unchanged downstream.
                junction_dropped_shapes.append(shape)
                continue
        if len(ring) < 3:
            continue
        ring = _drop_colinear_boundary_vertices(
            ring, corner_elev, shared_junction_buckets)
        if len(ring) < 3:
            continue
        ring = _drop_spike_vertices(ring)
        if len(ring) < 3:
            continue
        ring = _drop_sliver_corners(ring)
        if len(ring) < 3:
            continue
        try:
            poly = Polygon(ring)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if (poly.is_empty or poly.geom_type != "Polygon"
                    or poly.area < 0.5):
                continue
        except Exception:
            continue
        junction_cleaned.append((shape, ring, poly))

    # ── Single-pass per-junction smoothing + shared-vertex
    # reconciliation.  Per-junction smoothing alone moves SOFT
    # (graph-derived) boundary vertices independently in each
    # junction that touches the same shared bucket, breaking the
    # shared-vertex invariant.  After all junctions have smoothed,
    # we walk every shared SOFT bucket and overwrite each junction's
    # value with the cross-junction average; this restores
    # consistency at shared vertices.
    #
    # Per the elevation field plan (2026-04-26), Step 3's iterative
    # variant of this — average then lock then re-smooth — was
    # tried and made within-shape grade worse, because forcing the
    # average at a shared vertex pulls it away from each junction's
    # locally-optimal smooth solution and constrains the next
    # smoothing round.  A single average pass keeps the per-
    # junction smoothing's local optimum and only sacrifices a
    # small grade residual at shared vertices.
    iter_results: List[Tuple["BuiltShape",
                             List[Tuple[float, float]],
                             List[float]]] = []
    for (shape, ring, poly) in junction_cleaned:
        ev_pairs = [_vertex_elev_anchored(x, y) for (x, y) in ring]
        vert_elev_raw: List[Optional[float]] = [p[0] for p in ev_pairs]
        is_anchor_list: List[bool] = [p[1] for p in ev_pairs]
        known = [e for e in vert_elev_raw if e is not None]
        fallback = sum(known) / len(known) if known else 0.0
        vert_elev = [e if e is not None else fallback
                     for e in vert_elev_raw]
        # Mode B (gated): replace free-vertex elevations with
        # bilinear samples from a 2D smoothed field.  Hard anchors
        # keep their corner-derived value.
        if USE_PER_POLYGON_ELEVATION_FIELD:
            try:
                _anchors_b: List[
                    Tuple[float, float, float]] = [
                        (rx, ry, e)
                        for (rx, ry), e, isa in zip(
                            ring, vert_elev, is_anchor_list)
                        if isa]
                sampler_b = _smooth_polygon_grid(
                    poly, _anchors_b, dem,
                    tile_lat, tile_lon, layout.anchor)
                if sampler_b is not None:
                    for k, (rx, ry) in enumerate(ring):
                        if is_anchor_list[k]:
                            continue
                        e_b = sampler_b(rx, ry)
                        if e_b is not None:
                            vert_elev[k] = float(e_b)
            except Exception:
                pass
        # Junction-local 1D boundary smoothing — pull free vertices
        # into anchor-band compliance with the polygon's hard
        # anchors.
        vert_elev = _smooth_junction_boundary(
            ring, vert_elev, is_anchor_list)
        iter_results.append((shape, ring, vert_elev))

    # Cross-junction shared-vertex reconciliation.  Per-junction
    # smoothing moves SOFT (graph-derived) shared boundary vertices
    # independently — at most shared buckets the moves agree across
    # all junctions to within a small tolerance (junctions have
    # similar local boundaries near a shared vertex), but at some
    # buckets the per-junction smoothing diverges and we'd ship a
    # visible cross-shape cliff.  For the divergent buckets, revert
    # to the cached SOFT graph value (which both junctions would
    # have read consistently anyway — same as baseline behaviour
    # for shared buckets).  Where smoothing converged consistently
    # across junctions, keep the smoothed value (this is the
    # within-shape-grade win that the SOFT classification
    # delivers).
    # Cross-junction reconciliation at shared SOFT buckets.  Each
    # junction's per-junction smoother moved its copy of the
    # shared vertex toward its own neighbour mean, so two
    # junctions touching the same bucket can disagree on the
    # final elevation.  For buckets where the disagreement
    # exceeds SHARED_AGREE_TOL_M we average across all junctions
    # and overwrite — restoring the shared-vertex invariant at
    # the cost of a small within-shape grade residual local to
    # that vertex.  For buckets where every junction's smoother
    # already converged to within tolerance, we leave the
    # per-junction values alone — preserving the SOFT
    # classification's within-shape-grade win.  The tolerance is
    # generous (0.10 m, 1 % grade across a 10 m shared edge) so
    # most shared SOFT buckets pass without averaging.  Constant
    # is module-level so both legacy and Mode B paths can use it.
    shared_smooth_samples: Dict[
        Tuple[int, int], List[float]] = {}
    for (shape, ring, vert_elev) in iter_results:
        seen: set = set()
        for (vx, vy), e in zip(ring, vert_elev):
            bucket = _corner_elevation_bucket(vx, vy)
            if bucket in seen:
                continue
            seen.add(bucket)
            if bucket in corner_elev:
                continue
            if bucket not in shared_junction_buckets:
                continue
            if bucket not in shared_junction_elev_soft:
                continue
            shared_smooth_samples.setdefault(
                bucket, []).append(e)
    shared_lock: Dict[Tuple[int, int], float] = {}
    for bucket, vals in shared_smooth_samples.items():
        if len(vals) >= 2 and (
                max(vals) - min(vals) > SHARED_AGREE_TOL_M):
            shared_lock[bucket] = sum(vals) / len(vals)
    if shared_lock:
        for (shape, ring, vert_elev) in iter_results:
            for k, (vx, vy) in enumerate(ring):
                bucket = _corner_elevation_bucket(vx, vy)
                if bucket in shared_lock:
                    vert_elev[k] = shared_lock[bucket]

    new_shapes: List[BuiltShape] = []
    triangle_count = 0
    grade_violations = 0
    # Carry every non-junction shape through unchanged, plus any
    # junction shape that failed pre-loop cleanup.
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            new_shapes.append(shape)
    for shape in junction_dropped_shapes:
        new_shapes.append(shape)
    # Emit each cleaned junction with its post-iteration vert_elev
    # passed through the existing densify + FLAT/PLANAR/COMPOUND
    # classifier.
    for (shape, ring, vert_elev) in iter_results:
        # Densify long boundary edges (user 2026-04-25): for any
        # ring segment longer than MAX_BOUNDARY_EDGE_M, insert
        # interpolated midpoints.  Each midpoint's elevation is
        # the linear interpolation between the edge's endpoints,
        # so the rendered surface along the edge is unchanged ON
        # THE EDGE — adjacent shapes that share this edge (rect
        # boundaries) interpolate to the same value at the midpoint
        # so no step is introduced.  The benefit is INTERIOR
        # triangulation: Triangle4XP's quality refinement can
        # connect interior Steiners to the new closer boundary
        # vertices, producing smaller triangles and gentler
        # gradients within the polygon (especially helpful for
        # long cut edges from hole-decomposition).
        densified_ring, densified_elev = _densify_long_boundary_edges(
            ring, vert_elev, neighbour_edges)
        # Validate: a densification midpoint can occasionally land
        # on a non-adjacent ring edge (concave polygons with
        # near-touches), turning a valid polygon into a self-
        # touching one.  X-Plane crashes on those.  Revert to the
        # pre-densification ring if the densified polygon fails
        # validity.
        if len(densified_ring) > len(ring):
            try:
                test_poly = Polygon(densified_ring)
                if test_poly.is_valid:
                    ring, vert_elev = densified_ring, densified_elev
            except Exception:
                pass

        # ── Surface-complexity classification ───────────────────
        # User 2026-04-25: only triangulate where the surface has
        # a compound slope.  Flat or planar regions can stay as a
        # single polygon — fewer shapes, cleaner OSM output, less
        # work for X-Plane's mesh builder.
        #
        #   * FLAT       — vertex elevations vary by < 0.2 % of the
        #                  polygon's bbox extent (with a floor of
        #                  ``FLAT_ABS_FLOOR_M``).  Emit one polygon
        #                  with a single ``altitude`` tag.
        #   * PLANAR     — every vertex sits within
        #                  ``PLANAR_RESIDUAL_M`` of the best-fit
        #                  plane through them.  Emit one polygon
        #                  with ``node_altitudes`` (X-Plane
        #                  triangulates internally; since the
        #                  surface is planar the result is identical
        #                  to our pre-triangulated mesh).
        #   * COMPOUND   — triangulate as before.
        FLAT_GRADE_THRESHOLD = 0.001      # 0.1 % (so 0.05 m at 50 m
                                           # bbox) — tightened to keep
                                           # the FLAT mean within
                                           # check_grade's
                                           # SHARED_NID_TOLERANCE_M
                                           # (0.15 m) of every anchor.
                                           # Plateau snapping makes
                                           # most flat regions
                                           # already fit this; what
                                           # doesn't falls to PLANAR
                                           # which emits per-vertex
                                           # elevations matching each
                                           # anchor exactly.
        FLAT_ABS_FLOOR_M = 0.05           # 5 cm — half of one stored
                                           # decimal, so the rounded
                                           # mean can't differ from
                                           # any vertex by more than
                                           # rounding noise.
        PLANAR_RESIDUAL_M = 0.30          # 30 cm — vertices
                                           # within this of best-fit
                                           # plane render virtually
                                           # identical to a
                                           # triangulated mesh
        # Polygon bbox extent.
        xs = [p[0] for p in ring]
        ys = [p[1] for p in ring]
        extent = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
        flat_thresh = max(FLAT_ABS_FLOOR_M,
                          FLAT_GRADE_THRESHOLD * extent)
        elev_range = max(vert_elev) - min(vert_elev)
        if elev_range < flat_thresh:
            # FLAT — single polygon, single altitude tag.
            try:
                flat_poly = Polygon(ring)
                if not flat_poly.is_valid:
                    flat_poly = flat_poly.buffer(0)
                if (flat_poly.geom_type == "Polygon"
                        and not flat_poly.is_empty
                        and flat_poly.area >= 0.5):
                    new_shape = BuiltShape(
                        polygon=flat_poly,
                        role=ROLE_JUNCTION,
                        ref=shape.ref)
                    new_shape.altitude = round(
                        sum(vert_elev) / len(vert_elev), 1)
                    new_shapes.append(new_shape)
                    continue
            except Exception:
                pass  # fall through to triangulation
        # Best-fit plane: solve ax + by + c = z via 3×3 normal eqns.
        # Classify as PLANAR only when the fit's residuals are tight
        # AND the plane's slope itself is grade-compliant.  A polygon
        # whose vertices lie on a 5 % plane has tiny residuals but
        # X-Plane would render a 5 % cross-grade across the polygon
        # — well above TAXI_MAX_GRADE.  Demote those to COMPOUND so
        # _smooth_junction_boundary's anchor-band clamping pulls the
        # free vertices into compliance instead.
        plane = _planar_fit(ring, vert_elev)
        residuals = plane[3] if plane is not None else None
        if plane is not None:
            slope_mag = math.hypot(plane[0], plane[1])
        else:
            slope_mag = float("inf")
        if (residuals is not None
                and max(residuals) < PLANAR_RESIDUAL_M
                and slope_mag <= TAXI_MAX_GRADE):
            # PLANAR — single polygon, per-vertex altitudes.
            try:
                planar_poly = Polygon(ring)
                if not planar_poly.is_valid:
                    planar_poly = planar_poly.buffer(0)
                if (planar_poly.geom_type == "Polygon"
                        and not planar_poly.is_empty
                        and planar_poly.area >= 0.5):
                    # node_altitudes spans the closed ring; match
                    # the order Polygon(ring).exterior.coords
                    # produced (which is ring + closing first).
                    closed_ring = list(planar_poly.exterior.coords)
                    closed_elev = [
                        round(float(_match_elev(rx, ry, ring,
                                                vert_elev)), 1)
                        for (rx, ry) in closed_ring]
                    new_shape = BuiltShape(
                        polygon=planar_poly,
                        role=ROLE_JUNCTION,
                        ref=shape.ref)
                    new_shape.node_altitudes = closed_elev
                    new_shapes.append(new_shape)
                    continue
            except Exception:
                pass  # fall through to triangulation

        # ── Compound-slope: emit as a single polygon with per-
        # vertex node_altitudes; defer triangulation to Triangle4XP.
        #
        # Earlier versions ear-clipped or Delaunay-triangulated
        # the junction here, emitting one patch.osm way per output
        # triangle.  Each triangle was geometrically determined
        # (3 vertices), giving Triangle4XP nothing to refine —
        # any sliver / steep-plane triangle from our triangulator
        # made it into the rendered mesh untouched.
        #
        # Submitting the WHOLE polygon as a single constraint
        # boundary lets Triangle4XP's quality-refinement pass
        # ([O4_Mesh_Utils.py:670] flag ``-pq``) insert interior
        # Steiners with min-angle ≥ 20°.  Steiner elevations are
        # bilinear-interpolated from the polygon's boundary
        # ``node_altitudes`` at mesh time — exactly how the legacy
        # Ortho4XP pavement smoothing avoided cliffs without any
        # explicit grade enforcement.
        try:
            poly_for_emit = Polygon(ring)
            if not poly_for_emit.is_valid:
                poly_for_emit = poly_for_emit.buffer(0)
            if (poly_for_emit.is_empty
                    or poly_for_emit.geom_type != "Polygon"
                    or poly_for_emit.area < 0.5):
                continue
        except Exception:
            continue
        # node_altitudes spans the closed ring (one value per
        # vertex including the closing-repeat).  Build it from the
        # smoothed vert_elev in shapely's emitted ring order.
        closed = list(poly_for_emit.exterior.coords)
        elev_for_ring: List[float] = []
        for (rx, ry) in closed:
            best_e = vert_elev[0]
            best_d2 = float("inf")
            for (tx, ty), te in zip(ring, vert_elev):
                d2 = (rx - tx) * (rx - tx) + (ry - ty) * (ry - ty)
                if d2 < best_d2:
                    best_d2 = d2
                    best_e = te
            elev_for_ring.append(round(float(best_e), 1))
        new_shape = BuiltShape(
            polygon=poly_for_emit,
            role=ROLE_JUNCTION,
            ref=shape.ref)
        if (max(elev_for_ring) - min(elev_for_ring)) < 0.05:
            # Pre-classifier missed; emit flat.
            new_shape.altitude = round(
                sum(elev_for_ring[:-1]) / len(elev_for_ring[:-1]),
                1)
        else:
            new_shape.node_altitudes = elev_for_ring
        new_shapes.append(new_shape)
        triangle_count += 1

    layout.shapes = new_shapes
    if grade_violations:
        # Surfaced via stderr so the user sees it during the test
        # tool run; not a hard failure.
        try:
            import sys
            sys.stderr.write(
                f"  [pav-builder] WARN: {grade_violations} junction "
                f"triangle(s) exceed {TAXI_MAX_GRADE * 100:.1f}% grade.\n")
        except Exception:
            pass
    return triangle_count


# Layer 2: free-vertex clamping with neighbour-boundary lookup ─────
#
# After _triangulate_junctions sets every junction's per-vertex
# elevation, walk each junction's boundary and tighten any FREE
# vertex (one whose elevation is NOT pinned to a rect/runway/
# terminal corner or to another junction's shared vertex) so that
# it satisfies grade compliance with every nearby shape boundary
# point — not just the same-polygon anchored vertices considered
# during the first smoothing pass.
#
# Why a separate pass:
#   The first-pass smoothing (_smooth_junction_boundary) only
#   considers anchors WITHIN the same polygon.  Adjacent junction
#   polygons that are 1.5–5 m apart but don't share node IDs
#   (HECA's "adjacent-but-not-shared" pattern) can have their
#   nearby boundary vertices drift to incompatible elevations,
#   producing the visible sunken-area / elevated-plateau cliffs
#   the user reported at SPJC and the 159% within-shape grade
#   violations at HECA.
#
# What this pass does NOT change:
#   * Rect / runway / terminal altitudes — those are derived from
#     the elevation graph and the 1.5 % grade rule along the taxi
#     network; this pass operates only on junction polygons.
#   * Anchored junction vertices (shared with another shape via
#     the OSM nid).  Moving them would break the shared-vertex
#     invariant.
#
# Where it operates:
#   For each junction polygon's free vertex V at coord (x, y) with
#   current elevation e_v, scan every other shape's boundary
#   edges within NEIGHBOUR_CLAMP_RADIUS_M.  Each nearby boundary
#   point E at distance d contributes a feasibility band
#   ``[E_elev − d × TAXI_MAX_GRADE, E_elev + d × TAXI_MAX_GRADE]``.
#   Intersect these bands and clip e_v.  When the intersection is
#   empty (anchors disagree), fall back to the midpoint — Layer 1
#   will eventually reconcile the conflicting anchors.
NEIGHBOUR_CLAMP_RADIUS_M = 5.0


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
            msg += ".\n"
            _sys.stderr.write(msg)
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
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] {icao}: overlap-clip pass — "
                f"{n_clipped} clip operation(s), "
                f"{n_dropped} shape(s) dropped.\n")
        except Exception:
            pass


def _enforce_shared_vertices(layout: "PavementLayout",
                             tol: float = 1.5) -> None:
    """Collapse all emitted-shape vertices that lie within ``tol``
    of each other to a single canonical point (the cluster mean),
    then rewrite each shape's polygon with those canonical vertices.

    Implements rule 16 (exact shared vertices between adjacent
    shapes).  Must run AFTER all shapes are emitted.
    """
    # Gather every vertex with a (shape_idx, is_interior, ring_idx,
    # vert_idx) handle so we can rewrite them in place.
    handles: List[Tuple[int, int, int, int, Tuple[float, float]]] = []
    for si, shape in enumerate(layout.shapes):
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        ext = list(poly.exterior.coords)
        if ext and ext[0] == ext[-1]:
            ext = ext[:-1]
        for vi, v in enumerate(ext):
            handles.append((si, 0, 0, vi, (v[0], v[1])))
        for ri, ring in enumerate(poly.interiors):
            rc = list(ring.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            for vi, v in enumerate(rc):
                handles.append((si, 1, ri, vi, (v[0], v[1])))
    if not handles:
        return

    from collections import defaultdict

    # Union-find with O(n²) pair scan.  n is typically 200-2000
    # across both airports, well within millisecond range, and the
    # simpler code eliminates any spatial-index off-by-one bugs.
    n = len(handles)
    parent = list(range(n))

    def _find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def _union(a, b):
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[ra] = rb

    coords_only = [h[4] for h in handles]
    for i in range(n):
        ix, iy = coords_only[i]
        for j in range(i + 1, n):
            jx, jy = coords_only[j]
            dx = ix - jx
            dy = iy - jy
            if dx > tol or dx < -tol or dy > tol or dy < -tol:
                continue
            if math.hypot(dx, dy) <= tol:
                _union(i, j)

    # Compute cluster centroids (mean of member coords).
    cluster_members: Dict[int, List[int]] = defaultdict(list)
    for i in range(len(handles)):
        cluster_members[_find(i)].append(i)
    canonical: Dict[int, Tuple[float, float]] = {}
    for root, members in cluster_members.items():
        sx = sum(handles[m][4][0] for m in members) / len(members)
        sy = sum(handles[m][4][1] for m in members) / len(members)
        canonical[root] = (sx, sy)

    # Rewrite each shape's rings with the canonical coords.
    new_coords_by_shape: Dict[int, Dict[Tuple[int, int, int],
                                        Tuple[float, float]]] = defaultdict(dict)
    for i, h in enumerate(handles):
        si, is_int, ri, vi, _orig = h
        new_coords_by_shape[si][(is_int, ri, vi)] = canonical[_find(i)]

    for si, shape in enumerate(layout.shapes):
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        if si not in new_coords_by_shape:
            continue
        # Rebuild exterior.
        ext = list(poly.exterior.coords)
        if ext and ext[0] == ext[-1]:
            ext = ext[:-1]
        new_ext = [new_coords_by_shape[si].get((0, 0, vi), ext[vi])
                   for vi in range(len(ext))]
        # Drop consecutive duplicates that arose from clustering.
        dedup_ext: List[Tuple[float, float]] = []
        for c in new_ext:
            if not dedup_ext or math.hypot(
                    c[0] - dedup_ext[-1][0],
                    c[1] - dedup_ext[-1][1]) > 0.05:
                dedup_ext.append(c)
        if (len(dedup_ext) >= 2
                and math.hypot(dedup_ext[0][0] - dedup_ext[-1][0],
                               dedup_ext[0][1] - dedup_ext[-1][1]) < 0.05):
            dedup_ext = dedup_ext[:-1]
        if len(dedup_ext) < 3:
            # Shape collapsed to a degenerate sliver after cluster
            # rewrite (e.g. a thin grid-decomposition sliver whose
            # vertices got pulled together).  Empty the polygon so
            # the un-clustered original isn't left behind to
            # violate the shared-vertex invariant.  Empty polygons
            # are skipped by the validator and ``to_osm``.
            shape.polygon = Polygon()
            continue
        # Rebuild interiors.
        new_interiors: List[List[Tuple[float, float]]] = []
        for ri, ring in enumerate(poly.interiors):
            rc = list(ring.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            new_ring = [new_coords_by_shape[si].get((1, ri, vi), rc[vi])
                        for vi in range(len(rc))]
            dedup_ring: List[Tuple[float, float]] = []
            for c in new_ring:
                if not dedup_ring or math.hypot(
                        c[0] - dedup_ring[-1][0],
                        c[1] - dedup_ring[-1][1]) > 0.05:
                    dedup_ring.append(c)
            if len(dedup_ring) >= 3:
                new_interiors.append(dedup_ring)
        try:
            new_poly = Polygon(dedup_ext, new_interiors)
            # Always apply the rewrite — the invariant (shared
            # vertices with adjacent shapes) is our priority.  If
            # the rewrite makes the polygon self-intersect, accept
            # it: the OSM/JOSM representation stores the vertex
            # list as-is; downstream tools can validate separately.
            if (new_poly.geom_type == "Polygon"
                    and not new_poly.is_empty):
                shape.polygon = new_poly
        except Exception:
            pass


def _validate_shared_vertex_invariant(layout: "PavementLayout",
                                      tol: float = 1.5) -> None:
    """Assert that every pair of shape vertices is EITHER exactly
    equal (< 0.01 m after clustering) OR > ``tol`` apart.  A
    "close but not equal" pair violates rule 16 and signals a
    clustering bug.  Raises RuntimeError on violation.
    """
    verts: List[Tuple[int, Tuple[float, float]]] = []
    for si, shape in enumerate(layout.shapes):
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        ext = list(poly.exterior.coords)
        if ext and ext[0] == ext[-1]:
            ext = ext[:-1]
        for v in ext:
            verts.append((si, (v[0], v[1])))
    if len(verts) < 2:
        return
    # Grid-bucket check: every pair within tol must be within 0.01.
    from collections import defaultdict
    cell = tol * 2.0
    buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, (_, (x, y)) in enumerate(verts):
        buckets[(int(x // cell), int(y // cell))].append(i)
    for (gx, gy), idxs in buckets.items():
        neigh: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                neigh.extend(buckets.get((gx + dx, gy + dy), []))
        for a in idxs:
            ax, ay = verts[a][1]
            for b in neigh:
                if b <= a:
                    continue
                bx, by = verts[b][1]
                d = math.hypot(ax - bx, ay - by)
                if 0.01 < d <= tol:
                    raise RuntimeError(
                        f"Shared-vertex invariant violated: shape {verts[a][0]}"
                        f" @ ({ax:.3f},{ay:.3f}) and shape {verts[b][0]}"
                        f" @ ({bx:.3f},{by:.3f}) are {d:.3f} m apart"
                        f" (within tol={tol} m but not exactly equal).")


# ──────────────────────────────────────────────────────────────────
# Centerline-based taxi rect builder
# ──────────────────────────────────────────────────────────────────

JUNCTION_CLUSTER_DIST_M = 40.0  # merge junction nodes within this distance
                                # — between 80 (too coarse, merged
                                # distinct crossings) and 25 (too
                                # fine, created spurious crossings at
                                # every sub-way endpoint)
JUNCTION_RADIUS_SCALE = 1.5     # disc radius = local_half_width × this


def _terminal_pad_from_building(
    building: Polygon,
    pav_polys: List[Polygon],
) -> Optional[Polygon]:
    """Expand an OSM building outline into the terminal pad.

    The pad is the apt.dat pavement polygon that **contains** the
    building (the pavement area dedicated to the terminal).  If no
    pavement polygon contains the building centroid, fall back to
    a buffered version of the building.
    """
    if building.is_empty:
        return None
    ctr = building.centroid
    best = None
    best_area = -1.0
    for pav in pav_polys:
        if pav.contains(ctr):
            # Prefer the SMALLEST containing polygon (most specific).
            if best is None or pav.area < best_area:
                best = pav
                best_area = pav.area
    if best is not None:
        return best
    # Fallback: buffered building
    try:
        buf = building.buffer(20.0)
    except Exception:
        return building
    return buf if buf.geom_type == "Polygon" else None


def _extract_osm_terminals(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    relations: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
) -> List[Polygon]:
    """Extract terminal-class building polygons (ways OR
    multipolygon relations with outer rings) in meter space.

    Per user 2026-04-28: also recognize ``aeroway=hangar`` and
    ``aeroway=tower`` because OSM mappers at some airports
    (notably HECA Cairo) tag passenger-terminal-class buildings
    as hangars rather than terminals.  Functionally these are
    identical for the pavement-grading pipeline: flat, fixed-
    altitude structures on apron pavement that the surrounding
    apron should grade to.

    Guard against false positives: ``aeroway=hangar`` /
    ``aeroway=tower`` are only used when the airport has NO
    ``aeroway=terminal`` items.  At airports where mappers DID
    use ``aeroway=terminal`` (e.g. SPJC), the explicit terminals
    are authoritative and the hangar/tower buildings are likely
    actual hangars / towers that overlap pavement and would
    cause overlap-clip to malform sloping rects.

    All accepted categories are emitted as ROLE_TERMINAL.
    """
    # Detect whether this airport uses explicit aeroway=terminal.
    has_explicit_terminal = any(
        tags.get("aeroway") == "terminal"
        for _wid, _nds, tags in ways)
    has_explicit_terminal = has_explicit_terminal or any(
        tags.get("aeroway") == "terminal"
        for _rid, _wids, tags in relations)
    if has_explicit_terminal:
        terminal_aeroway_tags = {"terminal"}
    else:
        terminal_aeroway_tags = {"terminal", "hangar", "tower"}
    TERMINAL_AEROWAY_TAGS = terminal_aeroway_tags
    out: List[Polygon] = []
    way_by_id = {wid: (nds, tags) for wid, nds, tags in ways}

    def _ring_polygon(nds: List[str]) -> Optional[Polygon]:
        pts = []
        for n in nds:
            if n in nodes:
                lat, lon = nodes[n]
                pts.append(to_m(lon, lat))
        if len(pts) < 3:
            return None
        try:
            p = Polygon(pts).buffer(0)
        except Exception:
            return None
        if p.is_empty:
            return None
        if p.geom_type == "MultiPolygon":
            p = max(p.geoms, key=lambda g: g.area)
        return p if p.geom_type == "Polygon" else None

    # Way terminals
    for wid, nds, tags in ways:
        if tags.get("aeroway") not in TERMINAL_AEROWAY_TAGS:
            continue
        p = _ring_polygon(nds)
        if p is not None and p.area >= 100.0:
            out.append(p)

    # Relation terminals — emit ONE simplified polygon per
    # aeroway=terminal relation that captures the full building
    # extent without per-jet-bridge fine detail.
    #
    # Previously the code took only the largest connected component
    # of the union of all outer rings.  That works when the largest
    # piece dominates (e.g. SPJC rel -222: 28,883 m² largest covers
    # the bulk of the terminal-between-runways).  It FAILS when the
    # building is fragmented into many comparable pieces (e.g. SPJC
    # rel -221: largest 6,557 m² is only 49 % of the total
    # 13,500 m² building footprint — missing the satellite concourses).
    #
    # Strategy (user 2026-04-24): take the convex hull of the union
    # of all "significant" components (≥ ``MIN_TERMINAL_COMPONENT_M2``).
    # If the largest component is already > LARGEST_DOMINATES_FRAC of
    # the total significant area, use it as-is (preserves the
    # well-shaped output for dominant-piece terminals).  Otherwise
    # use the convex hull (captures full multi-piece extent in a
    # single simplified polygon).
    MIN_TERMINAL_COMPONENT_M2 = 500.0
    LARGEST_DOMINATES_FRAC = 0.7
    for rid, outer_wids, tags in relations:
        if tags.get("aeroway") not in TERMINAL_AEROWAY_TAGS:
            continue
        rings = []
        for wid in outer_wids:
            if wid not in way_by_id:
                continue
            nds, _ = way_by_id[wid]
            p = _ring_polygon(nds)
            if p is not None:
                rings.append(p)
        if not rings:
            continue
        merged = unary_union(rings).buffer(0)
        # Collect significant components.
        if merged.geom_type == "Polygon":
            components = [merged] if merged.area >= 100.0 else []
        elif merged.geom_type == "MultiPolygon":
            components = [g for g in merged.geoms
                          if g.geom_type == "Polygon"
                          and g.area >= MIN_TERMINAL_COMPONENT_M2]
        else:
            continue
        if not components:
            continue
        components.sort(key=lambda g: -g.area)
        total = sum(g.area for g in components)
        # Single dominant component → use directly.
        if components[0].area / total >= LARGEST_DOMINATES_FRAC:
            out.append(components[0])
            continue
        # Multi-piece terminal → emit the convex hull as a single
        # simplified polygon spanning the full footprint.
        try:
            hull = unary_union(components).convex_hull
        except Exception:
            out.append(components[0])
            continue
        if hull.geom_type == "Polygon" and hull.area >= 100.0:
            out.append(hull)
        else:
            out.append(components[0])
    return out


def _find_junction_points(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    osm_centerlines: Optional[List[Tuple[LineString, str]]] = None,
) -> List[Tuple[float, float]]:
    """Identify junction POINTS — OSM nodes shared by ≥ 2 DIFFERENT refs.

    Per user rule: pure bends within one taxi (two same-ref ways
    meeting at a node) are NOT junctions — they emit as adjacent
    same-role rects sharing a vertex line.  Only nodes where
    multiple distinct refs (or an unrefed way + a refed way) meet
    are junction candidates.

    Candidates within ``JUNCTION_CLUSTER_DIST_M`` of each other
    collapse into one cluster; the cluster centroid is the
    junction point.
    """
    from collections import defaultdict
    refs_at_node: Dict[str, set] = defaultdict(set)
    refed_taxi_nodes: set = set()
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        ref = tags.get("ref", "")
        if not ref:
            continue
        refed_taxi_nodes.update(nds)
        for n in nds:
            refs_at_node[n].add(ref)

    # Second pass: unrefed taxi ways act as CONNECTORS between
    # refed taxis at SPJC (e.g. 6 short unrefed ways bridge V to U,
    # marking the chart-level V↔U intersection points that split V
    # into multiple rects).  Only contribute a connector node when
    # it's also on a refed taxi way — this filters pure apron-area
    # markings (whose nodes touch only other unrefed ways).  Length
    # cap filters the long apron-boundary unrefed ways.
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        if tags.get("ref", ""):
            continue
        # Compute unrefed way length
        path_len = 0.0
        prev = None
        for n in nds:
            if n not in nodes:
                continue
            lat, lon = nodes[n]
            cur = to_m(lon, lat)
            if prev is not None:
                path_len += math.hypot(cur[0] - prev[0], cur[1] - prev[1])
            prev = cur
        if path_len > 200.0 or path_len < 20.0:
            continue
        for n in nds:
            if n in refed_taxi_nodes:
                refs_at_node[n].add("_conn")

    candidates: List[Tuple[float, float]] = []
    for nid, refs in refs_at_node.items():
        if len(refs) < 2:
            continue
        if nid not in nodes:
            continue
        lat, lon = nodes[nid]
        candidates.append(to_m(lon, lat))

    # ALSO add geometric crossing points between different-ref
    # centerlines (helps SPLP where few OSM nodes are shared).
    if osm_centerlines:
        for i in range(len(osm_centerlines)):
            ls1, ref1 = osm_centerlines[i]
            for j in range(i+1, len(osm_centerlines)):
                ls2, ref2 = osm_centerlines[j]
                if ref1 and ref2 and ref1 == ref2:
                    continue
                if not ls1.intersects(ls2):
                    continue
                try:
                    inter = ls1.intersection(ls2)
                except Exception:
                    continue
                if inter.is_empty:
                    continue
                if inter.geom_type == "Point":
                    candidates.append((inter.x, inter.y))
                elif inter.geom_type == "MultiPoint":
                    for p in inter.geoms:
                        candidates.append((p.x, p.y))
                elif inter.geom_type == "LineString":
                    candidates.append(inter.centroid.coords[0])

    # Cluster within JUNCTION_CLUSTER_DIST_M (greedy single-link)
    clusters: List[List[Tuple[float, float]]] = []
    for pt in candidates:
        placed = False
        for cl in clusters:
            if any(math.hypot(pt[0]-q[0], pt[1]-q[1]) <= JUNCTION_CLUSTER_DIST_M
                   for q in cl):
                cl.append(pt)
                placed = True
                break
        if not placed:
            clusters.append([pt])

    return [(sum(p[0] for p in cl)/len(cl),
             sum(p[1] for p in cl)/len(cl)) for cl in clusters]


def _build_junctions_from_rect_endpoints(
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    merge_dist: float,
    pav_union: Optional[Polygon],
    terminal_union: Optional[Polygon] = None,
) -> List[Polygon]:
    """Build junctions from rect endpoint clusters (user's approach).

    Algorithm:
      1. Collect each rect's 2 axis endpoints + 2 corner vertices at
         each end (total: 2 endpoints × 2 corners = 4 corners per rect).
      2. Cluster axis endpoints by single-link within ``merge_dist``.
      3. For each cluster of ≥ 2 endpoints:
         * If all endpoints share the same ref → same-taxi bend (no
           junction emitted; same-ref rects connect via shared vertices
           handled elsewhere).
         * Else → emit a junction polygon whose vertices are the 2
           corner vertices of each participating rect at the cluster.

    The polygon vertices are ordered angularly around the cluster
    centroid, giving a star polygon that wraps through each rect's
    corner pair.
    """
    if not taxi_rects:
        return []

    # Endpoint records: (rect_idx, end_index, axis_pt, corner_pair, ref)
    endpoints = []
    for i, (rect, axis, role, ref) in enumerate(taxi_rects):
        pairs = _rect_end_corners(rect, axis)
        if len(pairs) < 2:
            continue
        coords = list(axis.coords)
        endpoints.append((i, 0, coords[0], pairs[0], ref))
        endpoints.append((i, 1, coords[-1], pairs[1], ref))

    # Single-link cluster by axis-endpoint proximity
    n = len(endpoints)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i+1, n):
            ax = endpoints[i][2]
            bx = endpoints[j][2]
            if math.hypot(ax[0]-bx[0], ax[1]-bx[1]) <= merge_dist:
                union(i, j)

    clusters: Dict[int, List[int]] = {}
    for i in range(n):
        r = find(i)
        clusters.setdefault(r, []).append(i)

    out: List[Polygon] = []
    for cl in clusters.values():
        if len(cl) < 2:
            continue
        # Unique refs in cluster
        refs = {endpoints[i][4] for i in cl}
        # Unique rect ids (cluster can contain multiple ends of same rect)
        unique_rects = {endpoints[i][0] for i in cl}
        # Pure same-ref bend (only one ref AND only one rect pairs bends) — skip junction
        if len(refs) == 1 and len(unique_rects) <= 2:
            continue
        # Collect all corner vertices
        all_corners = []
        for i in cl:
            c1, c2 = endpoints[i][3]
            all_corners.append(c1)
            all_corners.append(c2)
        if len(all_corners) < 3:
            continue
        # Deduplicate near-identical corners
        uniq: List[Tuple[float, float]] = []
        for c in all_corners:
            if not any(math.hypot(c[0]-u[0], c[1]-u[1]) < 0.1 for u in uniq):
                uniq.append(c)
        if len(uniq) < 3:
            continue
        cx = sum(p[0] for p in uniq) / len(uniq)
        cy = sum(p[1] for p in uniq) / len(uniq)
        ordered = sorted(uniq, key=lambda p: math.atan2(p[1]-cy, p[0]-cx))
        try:
            poly = Polygon(ordered).buffer(0)
        except Exception:
            continue
        if poly.is_empty:
            continue
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon" or poly.area < 100.0:
            continue
        # Don't let junction bleed into a terminal
        if terminal_union is not None:
            try:
                poly = poly.difference(terminal_union)
            except Exception:
                pass
            if poly.is_empty:
                continue
            if poly.geom_type == "MultiPolygon":
                poly = max(poly.geoms, key=lambda g: g.area)
            if poly.geom_type != "Polygon":
                continue
        out.append(poly)
    return out


def _rect_end_corners(rect: Polygon, axis: LineString
                      ) -> List[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Return the 2 pairs of corners at the rect's 2 short ends.

    For a 4-corner rect built as [p1+h*perp, p2+h*perp, p2-h*perp,
    p1-h*perp], the "start" end has corners [idx 0, idx 3] and
    the "end" end has corners [idx 1, idx 2].  This pairing
    matters because junction polygons are built from the pair at
    whichever end of the rect meets the junction centroid.
    """
    coords = list(rect.exterior.coords)
    if len(coords) < 5:
        return []
    # idx 0 = p1+perp, idx 1 = p2+perp, idx 2 = p2-perp, idx 3 = p1-perp
    start_pair = (coords[0], coords[3])   # corners at axis start
    end_pair = (coords[1], coords[2])     # corners at axis end
    return [start_pair, end_pair]


def _build_junction_constructive(
    cluster_centroid: Tuple[float, float],
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    pav_union,
    terminal_union,
    max_corner_dist_m: float = 80.0,
    local_disc_radius_m: float = 120.0,
    max_arc_vertices: int = 4,
) -> Optional[Polygon]:
    """Constructive junction-polygon build per the user's
    authoritative shape rule
    (memory: feedback_shape_rules):

      1. One vertex per incoming rect corner (outer corners at the
         rect end that meets this junction).
      2. Between consecutive corners belonging to DIFFERENT rects,
         trace apt.dat pavement vertices that lie on the local pav
         boundary arc between them.
      3. Between consecutive corners of the SAME rect, connect
         directly (that's the rect's short edge — the junction's
         inner face on that side).

    Arcs are bounded by clipping pav_union to a disc of radius
    ``local_disc_radius_m`` around the cluster centroid; this keeps
    the boundary walk local and avoids wrap-around issues at large
    multi-component pavements.

    Returns ``None`` if fewer than 2 rect ends cluster here or
    construction fails.
    """
    if pav_union is None or pav_union.is_empty:
        return None
    cx, cy = cluster_centroid
    cp = Point(cx, cy)
    if terminal_union is not None and not terminal_union.is_empty:
        try:
            if terminal_union.contains(cp):
                return None
        except Exception:
            pass

    # Gather rect ends near the cluster centroid.
    rect_ends: List[Tuple[int, Tuple[float, float], Tuple[float, float]]] = []
    for i, (rect, axis, role, ref) in enumerate(taxi_rects):
        coords_ax = list(axis.coords)
        if len(coords_ax) < 2:
            continue
        ax_start = coords_ax[0]
        ax_end = coords_ax[-1]
        d_start = math.hypot(ax_start[0] - cx, ax_start[1] - cy)
        d_end = math.hypot(ax_end[0] - cx, ax_end[1] - cy)
        if min(d_start, d_end) > max_corner_dist_m:
            continue
        pairs = _rect_end_corners(rect, axis)
        if len(pairs) < 2:
            continue
        idx = 0 if d_start <= d_end else 1
        c1, c2 = pairs[idx]
        rect_ends.append((i, c1, c2))

    if len(rect_ends) < 2:
        return None

    # Local pavement: intersect pav with a disc around the cluster.
    try:
        disc = cp.buffer(local_disc_radius_m)
        local_pav = pav_union.intersection(disc)
    except Exception:
        return None
    if local_pav.is_empty:
        return None
    if local_pav.geom_type == "MultiPolygon":
        # Take the component containing (or closest to) the centroid.
        best = None
        best_d = float('inf')
        for g in local_pav.geoms:
            if g.geom_type != "Polygon":
                continue
            d = g.distance(cp)
            if d < best_d:
                best_d = d
                best = g
        if best is None:
            return None
        local_pav = best
    if local_pav.geom_type != "Polygon":
        return None

    # Projection helper: use the exterior ring as a line for param.
    ext_coords = list(local_pav.exterior.coords)
    if len(ext_coords) < 4:
        return None
    # LinearRing is closed (first == last); use as LineString for project()
    ext_ls = LineString(ext_coords)
    ext_length = ext_ls.length
    if ext_length <= 0:
        return None

    # Project each rect outer corner onto the exterior ring.
    corner_data: List[Tuple[float, Tuple[float, float], int]] = []
    for rect_idx, c1, c2 in rect_ends:
        for c in (c1, c2):
            try:
                param = ext_ls.project(Point(c))
            except Exception:
                continue
            proj_pt = ext_ls.interpolate(param)
            # If the projection is far (corner inside pav interior,
            # not on boundary — e.g. the local disc cut through a
            # rect short edge), use the original corner position.
            if proj_pt.distance(Point(c)) > 20.0:
                continue
            corner_data.append((param, c, rect_idx))

    # Need at least 3 corners for a polygon.
    if len(corner_data) < 3:
        return None

    # Sort corners by boundary param.  This orders them around the
    # local pav exterior ring; same-rect corner pairs will typically
    # be adjacent (a rect's 2 outer corners land close together on
    # the boundary).
    corner_data.sort(key=lambda t: t[0])

    # Collect all exterior-ring vertices with their params for arc walks.
    ext_verts = ext_coords[:-1] if ext_coords[0] == ext_coords[-1] else ext_coords
    ext_vert_params: List[Tuple[float, Tuple[float, float]]] = []
    acc = 0.0
    for i, v in enumerate(ext_verts):
        if i > 0:
            acc += math.hypot(v[0] - ext_verts[i-1][0],
                              v[1] - ext_verts[i-1][1])
        ext_vert_params.append((acc, (v[0], v[1])))

    # Build polygon by walking sorted corners.  Between consecutive
    # corners of different rects, insert all exterior-ring vertices
    # whose params lie between the 2 corner params (forward direction,
    # with wrap-around from last to first).
    n = len(corner_data)
    poly_coords: List[Tuple[float, float]] = []
    for i in range(n):
        cur_param, cur_xy, cur_rect = corner_data[i]
        nxt_param, nxt_xy, nxt_rect = corner_data[(i + 1) % n]
        poly_coords.append(cur_xy)
        if cur_rect == nxt_rect:
            # Same-rect: direct connection (rect short edge), no arc.
            continue
        # Different rect: walk exterior ring from cur_param to nxt_param
        # in increasing-param direction (wrap at end).
        if i == n - 1 or nxt_param < cur_param:
            arc_verts = (
                [v for (vp, v) in ext_vert_params if vp > cur_param] +
                [v for (vp, v) in ext_vert_params if vp < nxt_param])
        else:
            arc_verts = [v for (vp, v) in ext_vert_params
                         if cur_param < vp < nxt_param]
        # Per user rule (2026-04-18): "use a maximum of 4 points
        # between each" rect corner pair.  Sub-sample evenly.
        if len(arc_verts) > max_arc_vertices:
            step = len(arc_verts) / max_arc_vertices
            arc_verts = [arc_verts[int(k * step)]
                         for k in range(max_arc_vertices)]
        for v in arc_verts:
            poly_coords.append(v)

    if len(poly_coords) < 3:
        return None
    try:
        poly = Polygon(poly_coords).buffer(0)
    except Exception:
        return None
    if poly.is_empty:
        return None
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    if poly.geom_type != "Polygon" or poly.area < 80.0:
        return None

    # Clip to the local pav (in case rect short-edges extend past it)
    try:
        poly = poly.intersection(local_pav)
    except Exception:
        pass
    if poly.is_empty:
        return None
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)
    if poly.geom_type != "Polygon":
        return None

    # Exclude terminal overlap
    if terminal_union is not None and not terminal_union.is_empty:
        try:
            poly = poly.difference(terminal_union)
        except Exception:
            pass
        if poly.is_empty or poly.geom_type not in ("Polygon", "MultiPolygon"):
            return None
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon":
            return None

    return poly


def _build_junction_polys_from_corners(
    junction_points: List[Tuple[float, float]],
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    pav: Optional[Polygon],
    terminal_union: Optional[Polygon] = None,
    max_corner_dist_m: float = 100.0,
) -> List[Polygon]:
    """Build each junction polygon from the CORNER VERTICES of the
    adjacent rects, per the user's rule:

        "One vertex for each corner vertex of the rects it's joining."

    For each cluster centroid:
      1. Gather rects whose axis start/end is within
         ``max_corner_dist_m`` of the centroid.
      2. For each gathered rect, take the 2 corner vertices at the
         NEARER end (start or end).
      3. Order all corner vertices angularly around the centroid.
      4. Emit as the junction polygon (simple polygon through these
         vertices).

    If fewer than 2 rect ends cluster here, no junction polygon is
    emitted (would be degenerate).

    Junction polys that overlap a terminal aren't emitted; terminal
    boundary is adjacent-direct (per user rule "aprons can join
    directly to taxiways without a junction").
    """
    if pav is None or not junction_points:
        return []

    polys: List[Polygon] = []
    for (cx, cy) in junction_points:
        jpt = Point(cx, cy)
        # Avoid emitting into a terminal
        if terminal_union is not None and terminal_union.contains(jpt):
            continue

        corner_pts: List[Tuple[float, float]] = []
        for rect, axis, role, ref in taxi_rects:
            coords_ax = list(axis.coords)
            if len(coords_ax) < 2:
                continue
            ax_start = coords_ax[0]
            ax_end = coords_ax[-1]
            d_start = math.hypot(ax_start[0]-cx, ax_start[1]-cy)
            d_end = math.hypot(ax_end[0]-cx, ax_end[1]-cy)
            near_d = min(d_start, d_end)
            if near_d > max_corner_dist_m:
                continue
            pairs = _rect_end_corners(rect, axis)
            if not pairs:
                continue
            # Take the pair at the closer axis end
            idx = 0 if d_start <= d_end else 1
            corner_pts.extend(pairs[idx])

        if len(corner_pts) < 3:
            continue
        # Order corners angularly around centroid
        ordered = sorted(corner_pts,
                         key=lambda p: math.atan2(p[1]-cy, p[0]-cx))
        try:
            poly = Polygon(ordered).buffer(0)
        except Exception:
            continue
        if poly.is_empty:
            continue
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon":
            continue
        if poly.area < 80.0:
            continue
        # Clip to pavement so junction doesn't escape the pavement footprint
        poly = poly.intersection(pav)
        if poly.is_empty:
            continue
        if poly.geom_type == "MultiPolygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        if poly.geom_type != "Polygon":
            continue
        polys.append(poly)
    return polys


RDP_SIMPLIFY_TOL_M = 1.0      # RDP tolerance after ref-merge
MIN_SEGMENT_LEN_M = 15.0      # drop segments shorter than this
SIGNIFICANT_BEND_DEG = 5.0    # only split parallels at bends this sharp
BEND_CLUSTER_M = 100.0        # cluster consecutive bends within this
                              # distance.  100m because a taxi's
                              # direction change at an intersection
                              # area (e.g. V bends slightly around
                              # V2's junction) may span 70-90m with
                              # 2 small bends marking the start and
                              # end of the transition.
CLOSE_INTERSECTION_M = 200.0  # intersections within this distance
                              # on one centerline may merge into a
                              # single junction region (rects stop
                              # short; interval is junction, no rect).
                              # The 60-200m zone is gated by a pav-
                              # width midpoint check (combined
                              # junction only if midpoint is wider
                              # than narrow).  Measured parallel
                              # spacings: Q-R=95m, L-M=47m, V-U=74m.
                              # OSM-fragmented V2/V3 junctions span
                              # up to 180m on V's axis.
                              # (RDP keeps small wobbles; target subdivides
                              # only at chart-level direction changes)
GAP_BRIDGE_MAX_M = 120.0       # bridge same-ref polyline gaps up to this
STUB_MAX_LEN_M = 250.0         # polylines <= this emit as one rect


def _bridge_same_ref_polylines(lines: List[LineString]
                               ) -> List[LineString]:
    """Greedily connect endpoints of same-ref polylines within
    ``GAP_BRIDGE_MAX_M`` by concatenation.  Produces fewer, longer
    polylines covering the ref's full extent.
    """
    if len(lines) < 2:
        return lines

    remaining = list(lines)
    merged_lines: List[LineString] = []
    while remaining:
        cur = remaining.pop(0)
        while True:
            cur_coords = list(cur.coords)
            cur_start = cur_coords[0]
            cur_end = cur_coords[-1]
            best_idx = -1
            best_d = GAP_BRIDGE_MAX_M
            best_order = None  # "append_end", "append_start", "append_end_rev", "append_start_rev"
            for i, other in enumerate(remaining):
                oc = list(other.coords)
                o_start, o_end = oc[0], oc[-1]
                for order, pair in (
                    ("append_end", (cur_end, o_start)),
                    ("append_end_rev", (cur_end, o_end)),
                    ("append_start", (cur_start, o_end)),
                    ("append_start_rev", (cur_start, o_start)),
                ):
                    d = math.hypot(pair[0][0]-pair[1][0],
                                   pair[0][1]-pair[1][1])
                    if d < best_d:
                        best_d = d
                        best_idx = i
                        best_order = order
            if best_idx < 0:
                merged_lines.append(cur)
                break
            other = remaining.pop(best_idx)
            oc = list(other.coords)
            if best_order == "append_end":
                cur = LineString(cur_coords + oc)
            elif best_order == "append_end_rev":
                cur = LineString(cur_coords + oc[::-1])
            elif best_order == "append_start":
                cur = LineString(oc + cur_coords)
            elif best_order == "append_start_rev":
                cur = LineString(oc[::-1] + cur_coords)
    return merged_lines


def _emit_primary_parallel_runway_stubs(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    runway_union: Optional[Polygon],
    pav_union: Optional[Polygon],
    apt_vertices: Optional[List[Tuple[float, float]]],
    existing_taxi_rects: List[Tuple[Polygon, LineString, str, str]],
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Emit an extra STUB rect at each primary parallel OSM path
    endpoint that terminates INSIDE the runway polygon.

    A / F / L OSM primary taxis at SPJC extend their polyline
    onto the runway pavement itself (the endpoint vertex sits
    inside the runway polygon).  The target hand-drawn OSM has
    a short wide STUB at that transition — the RAMP where the
    taxi meets the runway short-edge.  Detection:

      1. Merge each primary-parallel ref's OSM ways
         (linemerge + gap-bridge, same as the main extraction).
      2. Check each merged polyline's 2 endpoints.
         If an endpoint is INSIDE the runway polygon (or within
         10 m of it), the taxi terminates on runway pavement.
      3. Walk along the polyline from that endpoint toward the
         interior, until the path vertex distance to runway
         boundary exceeds ``STUB_EXIT_D_M`` (~80 m).  The
         "exit" vertex is where the path leaves the runway-
         apron ramp and enters normal taxi corridor.
      4. Emit a STUB rect CENTERED on the exit vertex along the
         local path direction, length ``STUB_LEN_M`` (~80 m).
         The rect width is the perpendicular pav half-width
         at the exit point × 2 (full pav width — the ramp is
         wider than a normal taxi).

    Returns a list of extra (rect, axis, role, ref) tuples to
    append to the main ``taxi_rects`` list.
    """
    if runway_union is None or runway_union.is_empty:
        return []
    if pav_union is None or pav_union.is_empty:
        return []

    # SPJC A/F have loop ramps at runway ends — target stubs sit at
    # the APEX of the loop where the ramp meets the normal taxi
    # corridor (d_rwy ≈ 100 m at a path vertex).  Vertex-based
    # exit at threshold 80 m lands on that apex vertex.  SPLP has
    # no loops — the primary curves smoothly into the runway and
    # target stubs sit MID-CURVE BETWEEN vertices at d_rwy ≈ 75 m.
    # Vertex-based exit over-shoots or under-shoots; interpolate
    # to the exact target d-value.  Separate thresholds:
    STUB_EXIT_D_M = 80.0
    STUB_INTERP_TARGET_D_UNREFED = 75.0
    STUB_LEN_M = 80.0       # target A/F stubs are 79–93 m
    ENDPOINT_INSIDE_TOL_M = 10.0  # allow near-boundary endpoints
    OUTSIDE_NEAR_RWY_M = 135.0  # NE-end endpoint within 135 m of
                                # runway but not inside (SPLP main
                                # taxi ends at 127 m from runway).
                                # SPJC L's internal endpoints sit
                                # at 148-150 m so 135 m excludes
                                # them while catching SPLP's
                                # runway-facing taxi end.
    UNREFED_MIN_LEN_M = 800.0   # unrefed ways that act as long
                                # primary parallels (SPLP main taxi
                                # is 2640 m unrefed)

    # Gather per-ref OSM lines (for parallel refs like A/F/L at
    # SPJC), PLUS unrefed taxi ways (for SPLP whose primary taxis
    # are all unrefed).  The unrefed ways are merged together and
    # only very long (>UNREFED_MIN_LEN_M) merged polylines qualify
    # as "primary parallels" for stub emission.
    by_ref: Dict[str, List[LineString]] = {}
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        ref = tags.get("ref", "")
        # Sub-refs (letter+digit, e.g. V1, L3) are short connector
        # spurs at every airport — never the long parallel taxi
        # we're hunting for here.  All other refs (or no ref) are
        # candidates; the length filter at the end (UNREFED_MIN_LEN_M
        # for unrefed; the by-ref endpoint test for everything else)
        # keeps only the long ones that touch the runway.
        if ref and any(c.isdigit() for c in ref):
            continue
        pts = []
        for n in nds:
            if n in nodes:
                lat, lon = nodes[n]
                pts.append(to_m(lon, lat))
        if len(pts) >= 2:
            try:
                by_ref.setdefault(ref, []).append(LineString(pts))
            except Exception:
                pass

    # Pre-compute existing rect union for overlap detection
    existing_rects_union = None
    if existing_taxi_rects:
        try:
            existing_rects_union = unary_union(
                [r for r, _, _, _ in existing_taxi_rects])
        except Exception:
            existing_rects_union = None

    new_stubs: List[Tuple[Polygon, LineString, str, str]] = []
    rwy_boundary = runway_union.boundary
    emitted_centers: List[Tuple[float, float]] = []
    DEDUP_DIST_M = 50.0  # de-dup stub centers within 50 m
    for ref, lines in by_ref.items():
        # Process INDIVIDUAL OSM ways (not merged).  Shared
        # endpoints between ways (e.g. SPLP -696729 and -696733
        # both end at (226,1182) which is an internal junction
        # to the runway-apron ramp) would become internal
        # vertices after linemerge — and thus hidden from
        # endpoint checking.  Per-way processing exposes each
        # OSM endpoint; the ``emitted_centers`` dedup step
        # coalesces identical runway-facing endpoints.
        if ref == "":
            # Unrefed airports: each individual way must be long
            # enough to represent a primary-parallel taxi.
            # SPLP main taxi (-696729) is 2640 m and -696733 is
            # 1456 m; both exceed the 800 m threshold.
            processed_lines = [ml for ml in lines
                               if ml.length >= UNREFED_MIN_LEN_M]
        else:
            # Refed parallels: each way is already part of a
            # named taxi.  Process all of them.
            processed_lines = list(lines)

        for ml in processed_lines:
            coords = list(ml.coords)
            if len(coords) < 2:
                continue
            # Check both endpoints for runway-terminating condition
            for end_idx in (0, -1):
                ep_pt = Point(coords[end_idx])
                d_ep = ep_pt.distance(rwy_boundary)
                endpoint_inside = (
                    runway_union.contains(ep_pt)
                    or d_ep <= ENDPOINT_INSIDE_TOL_M
                )
                endpoint_outside_near = (
                    not endpoint_inside
                    and d_ep <= OUTSIDE_NEAR_RWY_M
                )
                if not (endpoint_inside or endpoint_outside_near):
                    continue

                if endpoint_inside:
                    # Walk from the endpoint toward the interior
                    # until d_rwy > STUB_EXIT_D_M.  The "exit"
                    # vertex sits just outside the runway-apron
                    # ramp.  For REFED taxis (SPJC A/F) use the
                    # vertex directly as stub center — target
                    # happens to sit at a vertex (apex of loop
                    # ramp).  For UNREFED (SPLP curving primary)
                    # INTERPOLATE back to ``STUB_INTERP_TARGET_D_UNREFED``
                    # because the target sits BETWEEN vertices
                    # in a smooth curve.
                    step = 1 if end_idx == 0 else -1
                    exit_idx = None
                    prev_i = None
                    prev_d = 0.0
                    i = end_idx if end_idx >= 0 else len(coords) - 1
                    while 0 <= i < len(coords):
                        d = Point(coords[i]).distance(rwy_boundary)
                        inside = runway_union.contains(
                            Point(coords[i]))
                        if not inside and d > STUB_EXIT_D_M:
                            exit_idx = i
                            exit_d_val = d
                            break
                        prev_i = i
                        prev_d = d
                        i += step
                    if exit_idx is None:
                        continue
                    # Compute stub center (cx, cy)
                    if (not ref and prev_i is not None
                            and exit_d_val > prev_d
                            and prev_d < STUB_INTERP_TARGET_D_UNREFED
                            < exit_d_val):
                        frac = ((STUB_INTERP_TARGET_D_UNREFED
                                 - prev_d)
                                / (exit_d_val - prev_d))
                        interp_cx = (coords[prev_i][0]
                                     + frac * (coords[exit_idx][0]
                                               - coords[prev_i][0]))
                        interp_cy = (coords[prev_i][1]
                                     + frac * (coords[exit_idx][1]
                                               - coords[prev_i][1]))
                    else:
                        interp_cx = coords[exit_idx][0]
                        interp_cy = coords[exit_idx][1]
                else:
                    # Endpoint is OUTSIDE the runway but within
                    # OUTSIDE_NEAR_RWY_M.  The taxi curves to
                    # runway at this end but doesn't enter runway
                    # pavement (SPLP NE end at 127 m).  Use the
                    # endpoint itself as the stub center — target
                    # rect sits at the ramp where the taxi
                    # approaches runway.
                    exit_idx = 0 if end_idx == 0 else len(coords) - 1
                    interp_cx = coords[exit_idx][0]
                    interp_cy = coords[exit_idx][1]

                # Center stub on interpolated center (interp_cx,
                # interp_cy), length STUB_LEN_M along local path
                # direction (from prev to next of exit vertex).
                prev_idx = max(0, exit_idx - 1)
                next_idx = min(len(coords) - 1, exit_idx + 1)
                dx = coords[next_idx][0] - coords[prev_idx][0]
                dy = coords[next_idx][1] - coords[prev_idx][1]
                mag = math.hypot(dx, dy)
                if mag < 1e-6:
                    continue
                ux, uy = dx / mag, dy / mag
                cx, cy = interp_cx, interp_cy
                # First-pass axis at default length, used only to
                # probe the local pavement width.
                ax_start = (cx - ux * STUB_LEN_M / 2,
                            cy - uy * STUB_LEN_M / 2)
                ax_end = (cx + ux * STUB_LEN_M / 2,
                          cy + uy * STUB_LEN_M / 2)
                try:
                    probe_axis = LineString([ax_start, ax_end])
                except Exception:
                    continue
                _nat, _p90, narrow = _natural_half_width(
                    probe_axis, pav_union)
                if narrow < 3.5 or narrow > 50.0:
                    continue
                width = 2.0 * narrow
                # Width-based stub length (user 2026-04-27 spec):
                # the stub should be roughly square so its long
                # edges sit on the apron-narrowing pavement
                # boundary, corners snap there, and surrounding
                # junctions only connect at the short edges instead
                # of wrapping around the long edges.  Cap range
                # 50..80 m so the stub never collapses to a sliver
                # nor extends beyond the natural runway-apron
                # transition length.  (Diagonal stubs would use a
                # tighter cap; A/F/L at SPJC are perpendicular by
                # construction so the same formula applies.)
                target_len = max(50.0, min(STUB_LEN_M, width + 5.0))
                # ----- L-style pull-back -----
                # Distinguish L-style "primary parallel curving
                # into a stub" (narrow connector pavement) from
                # A/F-style "loop ramp" (wide rwy-end pavement).
                # The signal: pavement WIDTH at exit_idx.
                # • Wide  (≥ NARROW_PAV_M): loop ramp — keep stub
                #   centred AT exit_idx (the apex of the loop is
                #   exactly where the user wants the rect).
                # • Narrow (< NARROW_PAV_M): smooth curve from a
                #   primary parallel into the runway — the user
                #   2026-04-27 spec calls for the diagonal rule:
                #   pull the rect centre BACK along the path
                #   toward the runway by 0.35 × gap (gap = path
                #   length from runway-boundary crossing to
                #   exit_idx).  Matches the diagonal-rule 35 %
                #   retention used by ``_split_centerlines_at_points``
                #   for V3-style diagonal stubs.  At SPJC this
                #   lands the L stub at d_rwy ≈ 65 m (matching
                #   user's target hand-edit from d_rwy ≈ 103 m
                #   the loop-ramp rule gave).
                NARROW_PAV_M = 60.0
                PULL_BACK_FRAC = 0.35
                if (endpoint_inside and ref
                        and width < NARROW_PAV_M):
                    pull_path: List[Tuple[float, float]] = []
                    # Find runway-boundary crossing (last in-rwy
                    # vertex → first out-of-rwy vertex; intersect
                    # the connecting segment with the runway
                    # boundary).
                    s = 1 if end_idx == 0 else -1
                    k = (end_idx if end_idx >= 0
                         else len(coords) - 1)
                    last_in = None
                    while 0 <= k < len(coords):
                        ptk = Point(coords[k])
                        dk = ptk.distance(rwy_boundary)
                        if (runway_union.contains(ptk)
                                or dk <= ENDPOINT_INSIDE_TOL_M):
                            last_in = k
                            k += s
                        else:
                            break
                    if last_in is not None and 0 <= k < len(coords):
                        cross_pt = coords[k]
                        try:
                            seg = LineString(
                                [coords[last_in], coords[k]])
                            cd = seg.difference(runway_union)
                            if (not cd.is_empty
                                    and cd.geom_type == "LineString"):
                                cc = list(cd.coords)
                                d0 = math.hypot(
                                    cc[0][0] - coords[last_in][0],
                                    cc[0][1] - coords[last_in][1])
                                d1 = math.hypot(
                                    cc[-1][0] - coords[last_in][0],
                                    cc[-1][1] - coords[last_in][1])
                                cross_pt = (cc[0] if d0 < d1
                                            else cc[-1])
                        except Exception:
                            pass
                        pull_path.append(
                            (cross_pt[0], cross_pt[1]))
                        # Walk from there to exit_idx (inclusive).
                        m = k
                        while True:
                            pull_path.append(
                                (coords[m][0], coords[m][1]))
                            if m == exit_idx:
                                break
                            m += s
                            if not (0 <= m < len(coords)):
                                break
                    if len(pull_path) >= 2:
                        try:
                            gap_curve = LineString(pull_path)
                            gap = gap_curve.length
                        except Exception:
                            gap = 0.0
                        if gap > 30.0:
                            # New centre at (1 - PULL_BACK_FRAC)
                            # along the path FROM the runway side,
                            # i.e. PULL_BACK_FRAC * gap inland of
                            # the boundary crossing.
                            new_along = (1.0 - PULL_BACK_FRAC) * gap
                            cpt = gap_curve.interpolate(new_along)
                            cx, cy = cpt.x, cpt.y
                            # Local tangent at the new centre.
                            eps_ = max(1.0, gap * 0.02)
                            ta = max(0.0, new_along - eps_)
                            tb = min(gap, new_along + eps_)
                            pa = gap_curve.interpolate(ta)
                            pb = gap_curve.interpolate(tb)
                            tdx = pb.x - pa.x
                            tdy = pb.y - pa.y
                            tmag = math.hypot(tdx, tdy)
                            if tmag > 1e-6:
                                ux, uy = tdx / tmag, tdy / tmag
                ax_start = (cx - ux * target_len / 2,
                            cy - uy * target_len / 2)
                ax_end = (cx + ux * target_len / 2,
                          cy + uy * target_len / 2)
                try:
                    stub_axis = LineString([ax_start, ax_end])
                except Exception:
                    continue
                rect = _rect_from_axis_extended(
                    stub_axis, width, pav_union,
                    apt_vertices=apt_vertices)
                if rect is None or rect.is_empty:
                    continue
                if not rect.is_valid:
                    try:
                        rect = rect.buffer(0)
                    except Exception:
                        continue
                    if (rect.is_empty
                            or rect.geom_type != "Polygon"):
                        continue
                # Per user 2026-04-27: F-style runway-end stubs sit
                # in the wide runway-end ramp; the half-width probe
                # (capped at RAY_CAP_M = 40 m) under-sizes the rect
                # so the ramp extends past the rect's long edges and
                # the surrounding junction wraps around them.  After
                # the rect is built and snapped, ray-cast each corner
                # outward perpendicular to the axis until it hits
                # the apt.dat pavement boundary — turning the rect
                # into a trapezoid that covers the FULL ramp width
                # at each end independently.
                rect = _extend_rect_corners_perpendicular(
                    rect, stub_axis, pav_union)
                if rect is None or rect.is_empty:
                    continue
                # Skip if the stub would overlap an existing rect
                # or a same-ref rect's 30 m buffer (duplicates the
                # L SE stub that the main pipeline already emits).
                skip = False
                if existing_rects_union is not None:
                    try:
                        overlap = rect.intersection(
                            existing_rects_union).area
                        if overlap > rect.area * 0.2:
                            skip = True
                    except Exception:
                        pass
                if not skip and ref:
                    # Same-ref near-duplicate guard: applies only
                    # when ref is non-empty (avoid dropping all
                    # unrefed-airport stubs since every existing
                    # rect has ref="" too).
                    for er, _, _, eref in existing_taxi_rects:
                        if eref != ref:
                            continue
                        try:
                            if er.buffer(30.0).intersects(rect):
                                skip = True
                                break
                        except Exception:
                            pass
                if skip:
                    continue
                # Dedup by stub CENTER position — for unrefed
                # airports, two ways can share the same runway-
                # facing endpoint (e.g. SPLP -696729 and -696733
                # both end at (226,1182)) and produce near-
                # identical stubs.
                c = rect.centroid
                dup = False
                for (ex, ey) in emitted_centers:
                    if math.hypot(c.x - ex, c.y - ey) < DEDUP_DIST_M:
                        dup = True
                        break
                if dup:
                    continue
                emitted_centers.append((c.x, c.y))
                new_stubs.append(
                    (rect, stub_axis, ROLE_STUB, ref))
    return new_stubs


def _extract_osm_taxi_centerlines(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    rwy_centerlines: Optional[List[LineString]] = None,
) -> List[Tuple[LineString, str]]:
    """Extract one polyline segment per (ref, straight-run).

    Algorithm:
      1. Gather all OSM taxi ways per ref.
      2. linemerge() the per-ref ways into the fewest possible
         contiguous polylines.  OSM splits one physical taxi
         across many <way> rows at each node — merging restores
         the single polyline per strip.
      3. RDP-simplify each merged polyline at
         ``RDP_SIMPLIFY_TOL_M``.  Minor GPS wobble collapses;
         genuine bends remain as internal vertices.
      4. Split the simplified polyline at each remaining vertex
         into individual straight-run segments.
      5. Drop segments shorter than ``MIN_SEGMENT_LEN_M`` and
         (at airports with any refs) drop unrefed segments.

    Per user convention, stubs and cross-connectors typically merge
    down to a single polyline with a single straight run (emitted
    as 1 rect).  Long parallel taxis that physically bend (like L
    at SPJC) retain bend vertices and emit multiple rects that
    share corner vertices at the bend.
    """
    by_ref: Dict[str, List[LineString]] = {}
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        ref = tags.get("ref", "")
        pts = []
        for n in nds:
            if n in nodes:
                lat, lon = nodes[n]
                pts.append(to_m(lon, lat))
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        by_ref.setdefault(ref, []).append(ls)

    out: List[Tuple[LineString, str]] = []
    for ref, lines in by_ref.items():
        # Stage 1: contiguous-endpoint linemerge.
        if len(lines) > 1:
            try:
                merged = linemerge(MultiLineString(lines))
            except Exception:
                merged = None
            if merged is None or merged.is_empty:
                merged_lines = lines
            elif merged.geom_type == "LineString":
                merged_lines = [merged]
            else:
                merged_lines = list(merged.geoms)
        else:
            merged_lines = lines

        # Stage 2: gap bridging.  Bridge gaps for PRIMARY refs
        # (long continuous taxis) where OSM fragments across
        # intersections.  For SUB-REFS (letter+digit like V2, L3)
        # each OSM way is typically a separate short stub;
        # bridging their gaps creates fake segments through
        # non-pavement and confuses downstream width-profile
        # narrow-corridor detection.
        is_sub_ref = ref and any(c.isdigit() for c in ref)
        if ref and len(merged_lines) > 1 and not is_sub_ref:
            merged_lines = _bridge_same_ref_polylines(merged_lines)

        for ls in merged_lines:
            try:
                simp = ls.simplify(RDP_SIMPLIFY_TOL_M,
                                   preserve_topology=False)
            except Exception:
                continue
            scoords = list(simp.coords)
            if len(scoords) < 2:
                continue
            # Geometrically-straight-enough centerlines emit as ONE
            # rect rather than being bend-split.  This catches
            # continuous diagonal taxis at any airport (e.g. SPJC's
            # B/C/E/G, CYXY's E parallel) where the simplified
            # polyline has small bends that would otherwise get
            # bend-split into too-short fragments.  Sub-refs are
            # excluded because they're typically already short
            # connector spurs that benefit from bend-splitting at
            # their natural curve points.
            #
            # The chord/path-only test is INSUFFICIENT for taxis like
            # SPJC's L that have a long mostly-straight middle plus
            # tight curves at the ends (chord/path = 0.955 even though
            # L bends 23° at one end and 14-22° at the other).  A
            # post-pass STRAIGHT-ENOUGH check requires chord/path
            # close to 1 AND every interior bend below
            # ``MAX_INTERIOR_BEND_DEG``.  L's max bend is 23° → fails;
            # B/C/E/G/CYXY-E's max wobble is ~5° → still passes.
            MAX_INTERIOR_BEND_DEG = 15.0
            has_digit = bool(ref) and any(c.isdigit() for c in ref)
            if not has_digit:
                path_len = ls.length
                sc = list(simp.coords)
                if len(sc) >= 2 and path_len > 1e-6:
                    chord = math.hypot(sc[-1][0] - sc[0][0],
                                       sc[-1][1] - sc[0][1])
                    chord_ratio = chord / path_len
                    max_interior_bend = 0.0
                    for k in range(1, len(sc) - 1):
                        ax, ay = sc[k - 1]
                        bx, by = sc[k]
                        cx, cy = sc[k + 1]
                        v1x, v1y = bx - ax, by - ay
                        v2x, v2y = cx - bx, cy - by
                        m1 = math.hypot(v1x, v1y)
                        m2 = math.hypot(v2x, v2y)
                        if m1 < 1e-6 or m2 < 1e-6:
                            continue
                        d = (v1x * v2x + v1y * v2y) / (m1 * m2)
                        if d > 1.0:
                            d = 1.0
                        elif d < -1.0:
                            d = -1.0
                        ang = math.degrees(math.acos(d))
                        if ang > max_interior_bend:
                            max_interior_bend = ang
                    if (chord_ratio > 0.95
                            and max_interior_bend
                            < MAX_INTERIOR_BEND_DEG):
                        out.append((simp, ref))
                        continue
            # SHORT UNREFED runway-connecting stubs: SPLP has short
            # curvy unrefed taxis (e.g. way -696731, 165 m chord
            # 144 m) that link runway to apron/primary.  Target
            # emits a single rect in the middle of each.  Bend-
            # splitting fragments them into pieces too small to
            # survive the 40 m floor in `_split_centerlines_at_points`.
            # Emit atomically when: ref="" (unrefed) AND path < 300 m
            # AND one endpoint is inside the runway polygon.
            if (not ref
                    and ls.length < 300.0
                    and rwy_centerlines):
                try:
                    sc = list(simp.coords)
                    if len(sc) >= 2:
                        ep0 = Point(sc[0])
                        ep1 = Point(sc[-1])
                        ep0_near = any(
                            ep0.distance(r) < 30.0 for r in rwy_centerlines)
                        ep1_near = any(
                            ep1.distance(r) < 30.0 for r in rwy_centerlines)
                        if ep0_near or ep1_near:
                            out.append((simp, ref))
                            continue
                except Exception:
                    pass
            # All refs (including sub-refs) split at significant
            # bends.  Per user (2026-04-20 refined): intersections
            # + sharp curves define rect break points; there's no
            # reason sub-refs should be exempt from curve detection.
            is_parallel = True  # unified: all refs use bend-split
            if is_parallel:
                # Split at INTERNAL bends with angle change ≥
                # SIGNIFICANT_BEND_DEG, but cluster consecutive
                # bends within BEND_CLUSTER_M together.  A curve
                # (many tiny bends adding up to a big turn) counts
                # as ONE break point at its midpoint — matching
                # how the target treats a curve as a single logical
                # transition between rects.
                candidate_bends: List[int] = []
                for i in range(1, len(scoords) - 1):
                    a = scoords[i - 1]
                    b = scoords[i]
                    c = scoords[i + 1]
                    v1 = (b[0] - a[0], b[1] - a[1])
                    v2 = (c[0] - b[0], c[1] - b[1])
                    m1 = math.hypot(*v1)
                    m2 = math.hypot(*v2)
                    if m1 < 1e-6 or m2 < 1e-6:
                        continue
                    dot = (v1[0] * v2[0] + v1[1] * v2[1]) / (m1 * m2)
                    dot = max(-1.0, min(1.0, dot))
                    angle_change = math.degrees(math.acos(dot))
                    if angle_change >= SIGNIFICANT_BEND_DEG:
                        candidate_bends.append(i)
                # Cluster consecutive bends within BEND_CLUSTER_M of
                # each other.  Per user rule (2026-04-20): when a
                # primary taxi curves at the runway, emit the straight
                # portions as rects and leave the curve itself as
                # junction territory (no rect emitted for the curve).
                # → Each cluster yields TWO break indices
                #   (cluster start, cluster end) with the curve
                #   interval between them skipped.
                # A single isolated bend (cluster of 1) gives only
                # ONE break at itself.
                # Sub-refs (letter+digit: V3, V5, L1, …) are usually
                # short stub taxis whose straight portion is much
                # less than the primary's junction-bend-transition
                # span.  Using the primary's 100 m cluster distance
                # swallows a true 90°-corner stub's straight middle
                # run (e.g. V5 indices [13..15] are a 97 m straight
                # between two tight curves).  For sub-refs we cluster
                # bends far more conservatively so the straight run
                # between two curves can survive as its own segment.
                cluster_m = BEND_CLUSTER_M
                if ref and any(c.isdigit() for c in ref):
                    cluster_m = 30.0
                clusters: List[List[int]] = []
                for bi in candidate_bends:
                    if clusters and (scoords[bi][0] - scoords[clusters[-1][-1]][0])**2 + \
                            (scoords[bi][1] - scoords[clusters[-1][-1]][1])**2 \
                            <= cluster_m * cluster_m:
                        clusters[-1].append(bi)
                    else:
                        clusters.append([bi])
                # Build an ordered list of (break_index, kind) where
                # kind='point' (single bend) or 'interval_start' /
                # 'interval_end' (curve boundaries).
                events: List[Tuple[int, str]] = [(0, "point")]
                for cl in clusters:
                    if len(cl) == 1:
                        events.append((cl[0], "point"))
                    else:
                        # Only treat curve as junction interval if
                        # NEAR A RUNWAY (per user rule 3: "primary
                        # taxiway curves and intersects the runway"
                        # → straight rect, curve = junction, perp =
                        # stub).  Curves in the middle of the
                        # airport (e.g. A's gentle bend) stay as
                        # single break points.
                        near_rwy = False
                        if rwy_centerlines:
                            cluster_mid = scoords[cl[len(cl) // 2]]
                            cp = Point(cluster_mid)
                            for r in rwy_centerlines:
                                if cp.distance(r) < 200.0:
                                    near_rwy = True
                                    break
                        if near_rwy:
                            events.append((cl[0], "interval_start"))
                            events.append((cl[-1], "interval_end"))
                        else:
                            events.append((cl[len(cl) // 2], "point"))
                events.append((len(scoords) - 1, "point"))
                events.sort()
                # Walk events pair-wise; skip intervals between
                # interval_start and interval_end (that's the curve).
                for k in range(len(events) - 1):
                    i0, k0 = events[k]
                    i1, k1 = events[k + 1]
                    # Skip the curve interval itself.
                    if k0 == "interval_start" and k1 == "interval_end":
                        continue
                    if i0 == i1:
                        continue
                    try:
                        seg = LineString(scoords[i0:i1 + 1])
                    except Exception:
                        continue
                    if seg.is_empty or seg.length < MIN_SEGMENT_LEN_M:
                        continue
                    out.append((seg, ref))

    # Drop unrefed centerlines AT AIRPORTS THAT HAVE ANY REFED
    # CENTERLINES (per user 2026-04-27).  At SPJC etc. the OSM data
    # has comprehensive refs on real taxiways; unrefed lines that
    # remain are typically apron decorations / vehicle paths /
    # painted markings that, if extracted into rects, get inserted
    # INSIDE apron polygons — junction polygons then wrap around
    # them and produce visible elevation ridges where the rect's
    # short edges meet the junction at slightly different heights.
    # The user's directive: don't insert rects inside junction /
    # apron polygons.
    #
    # Previously this was an "if any_ref drop unrefed" filter in
    # this same place; an in-flight Session 8 change removed it on
    # the theory that downstream geometric-overlap dedup would
    # catch spurious unrefed sub-segments.  But spurious unrefed
    # apron lines DON'T overlap any refed rect (they sit inside an
    # apron region, not along a real taxi corridor) so the dedup
    # never fires for them, and HEAD-clean's clean baseline (47
    # rects, all refed at SPJC) regressed to 120 rects (65 of them
    # unrefed) inside apron areas.  Restoring the filter here.
    #
    # At airports with NO refed centerlines (CYXY where every OSM
    # taxi is unrefed) the filter is a no-op — every centerline is
    # kept.
    any_ref = any(r for (_, r) in out)
    if any_ref:
        out = [(ls, r) for (ls, r) in out if r]
    return out


def _insert_points_on_ring(
    ring_coords: List[Tuple[float, float]],
    pts: List[Tuple[float, float]],
    tol: float,
) -> List[Tuple[float, float]]:
    """Insert each point in ``pts`` as a vertex at its projected
    position on the closed ring (list of coords, first == last),
    if within ``tol``.  Returns the new ring coords (closed).
    Pure helper so both exterior and interior rings are handled
    uniformly."""
    if not pts or len(ring_coords) < 4:
        return ring_coords
    ring = LineString(ring_coords)
    inserts: List[Tuple[float, Tuple[float, float]]] = []
    for (x, y) in pts:
        p = Point(x, y)
        if p.distance(ring) > tol:
            continue
        try:
            param = ring.project(p)
            proj = ring.interpolate(param)
        except Exception:
            continue
        inserts.append((param, (proj.x, proj.y)))
    if not inserts:
        return ring_coords
    inserts.sort()
    coords = list(ring_coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    new_coords: List[Tuple[float, float]] = []
    cur_param = 0.0
    insert_i = 0
    for i in range(len(coords)):
        new_coords.append(coords[i])
        next_i = (i + 1) % len(coords)
        seg_len = math.hypot(coords[next_i][0] - coords[i][0],
                             coords[next_i][1] - coords[i][1])
        seg_end = cur_param + seg_len
        while (insert_i < len(inserts)
               and inserts[insert_i][0] < seg_end):
            new_coords.append(inserts[insert_i][1])
            insert_i += 1
        cur_param = seg_end
    new_coords.append(new_coords[0])
    return new_coords


def _insert_points_on_boundary(
    poly: Polygon,
    pts: List[Tuple[float, float]],
    tol: float = 2.0,
) -> Polygon:
    """Insert each point in ``pts`` as a vertex on the polygon's
    boundary (exterior + all interior rings) at its projected
    position, if within ``tol`` of the boundary.  Interior rings
    are preserved — critical when the polygon represents a
    pavement residue with rect-shaped holes.  Used to seam
    junction polygons with their neighbouring rect / terminal
    corners."""
    if not pts:
        return poly
    try:
        ext = list(poly.exterior.coords)
        new_ext = _insert_points_on_ring(ext, pts, tol)
        new_ints = []
        for ring in poly.interiors:
            ri = list(ring.coords)
            new_ints.append(_insert_points_on_ring(ri, pts, tol))
        new_poly = Polygon(new_ext, new_ints)
        if not new_poly.is_valid:
            new_poly = new_poly.buffer(0)
        if (new_poly.geom_type == "Polygon"
                and new_poly.is_valid and not new_poly.is_empty):
            return new_poly
    except Exception:
        pass
    return poly


def _split_by_width_profile(
    centerlines: List[Tuple[LineString, str]],
    pav_union: Polygon,
    probe_step_m: float = 5.0,
    wide_factor: float = 1.20,
    min_rect_len_m: float = 30.0,
) -> List[Tuple[LineString, str]]:
    """Split each centerline into NARROW-CORRIDOR intervals per
    user rule 4 (2026-04-20): rects cover only the narrowest
    straight sections; any widening (around intersections or
    terminal aprons) is junction territory and is skipped.

    For each line:
      1. Probe pav half-width at ``probe_step_m`` intervals.
      2. narrow_hw = 10th-percentile probe (robust narrow baseline).
      3. Interval flag per probe: NARROW if hw ≤ wide_factor × narrow_hw.
      4. Emit contiguous NARROW intervals ≥ min_rect_len_m as rects.
    """
    if not centerlines or pav_union is None or pav_union.is_empty:
        return centerlines
    from shapely.ops import substring
    pav_boundary = pav_union.boundary

    result: List[Tuple[LineString, str]] = []
    for ls, ref in centerlines:
        if ls.length < min_rect_len_m:
            result.append((ls, ref))
            continue
        n_probes = max(10, int(ls.length / probe_step_m))
        # sample (param, hw) pairs along the line
        samples: List[Tuple[float, float]] = []
        for i in range(n_probes + 1):
            t = i / n_probes * ls.length
            pt = ls.interpolate(t)
            hw = pt.distance(pav_boundary) if pav_union.contains(pt) else 0.0
            samples.append((t, hw))
        # narrow_hw = 10th-percentile of positive hw values
        hws = sorted(h for _, h in samples if h > 0)
        if not hws:
            result.append((ls, ref))
            continue
        narrow_hw = hws[max(1, len(hws) // 10)]
        wide_thresh = wide_factor * narrow_hw
        # Flag each sample narrow or wide
        is_narrow = [h > 0 and h <= wide_thresh for (_, h) in samples]
        # Find contiguous narrow intervals
        intervals: List[Tuple[float, float]] = []
        i = 0
        while i < len(samples):
            if not is_narrow[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(samples) and is_narrow[j + 1]:
                j += 1
            start_t = samples[i][0]
            end_t = samples[j][0]
            if end_t - start_t >= min_rect_len_m:
                intervals.append((start_t, end_t))
            i = j + 1
        if not intervals:
            # Whole line is "wide" — probably an apron traverse.
            # Drop it; user's target wouldn't emit a rect here.
            continue
        for (s, e) in intervals:
            try:
                seg = substring(ls, s, e)
            except Exception:
                continue
            if (seg.geom_type == "LineString"
                    and not seg.is_empty
                    and seg.length >= min_rect_len_m):
                result.append((seg, ref))
    return result


def _sub_ref_narrow_corridor(
    centerlines: List[Tuple[LineString, str]],
    pav_union: Polygon,
    probe_step_m: float = 4.0,
    wide_factor: float = 1.30,
    narrow_margin_frac: float = 0.15,
) -> List[Tuple[LineString, str]]:
    """For each sub-ref (ref like V1/V3/A1/L3 — letter+digit),
    replace its centerline(s) with the 70% middle slice of the
    LONGEST narrow-corridor interval.

    Algorithm:
      1. Perpendicular ray-cast half-width probes every
         ``probe_step_m`` along the centerline.
      2. narrow_hw = min probe (≥ 3.5 m floor).
      3. Flag each probe narrow/wide by ``wide_factor × narrow_hw``.
      4. Longest contiguous narrow interval → emit 70% middle.

    Also de-dupes multiple OSM ways with the same sub-ref by
    keeping the one whose selected-slice is longest and narrowest.
    """
    if not centerlines or pav_union is None or pav_union.is_empty:
        return centerlines
    from shapely.ops import substring

    RAY_CAP_M = 40.0
    RAY_STEP_M = 0.5

    def _perp_hw(line: LineString, t: float) -> float:
        dt = min(2.0, line.length * 0.05)
        t0 = max(0.0, t - dt)
        t1 = min(line.length, t + dt)
        a = line.interpolate(t0)
        b = line.interpolate(t1)
        tx, ty = b.x - a.x, b.y - a.y
        mag = math.hypot(tx, ty)
        if mag < 1e-6:
            return 0.0
        ux, uy = tx / mag, ty / mag
        nx, ny = -uy, ux
        pt = line.interpolate(t)
        ox, oy = pt.x, pt.y
        best = RAY_CAP_M
        for sign in (-1, 1):
            d = 0.0
            while d <= RAY_CAP_M:
                qx = ox + sign * nx * d
                qy = oy + sign * ny * d
                if not pav_union.contains(Point(qx, qy)):
                    if d < best:
                        best = d
                    break
                d += RAY_STEP_M
        return best

    # ICAO Code E taxi: 23m wide = 11.5m half-width.
    # Allow up to 16m half-width as "still on the taxi strip";
    # beyond that we're in a widening (intersection or apron).
    NARROW_TAXI_HW_M = 16.0

    def _narrow_slice(ls: LineString) -> Optional[Tuple[LineString, float]]:
        """Return (slice, avg_hw_in_narrow) for the 70% middle of
        the longest narrow-corridor interval along ls.  Uses a
        FIXED narrow-width threshold (NARROW_TAXI_HW_M) based on
        ICAO standards rather than per-line percentiles, which
        are unreliable on short or highly-curved sub-refs."""
        if ls.length < MIN_SEGMENT_LEN_M:
            return None
        n = max(10, int(ls.length / probe_step_m))
        samples: List[Tuple[float, float]] = []
        for i in range(n + 1):
            t = i / n * ls.length
            hw = _perp_hw(ls, t)
            samples.append((t, hw))
        if not any(h > 0 for _, h in samples):
            return None
        is_narrow = [3.5 <= h <= NARROW_TAXI_HW_M for (_, h) in samples]
        # Longest contiguous narrow interval.
        best_i, best_j = -1, -1
        i = 0
        while i < len(samples):
            if not is_narrow[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(samples) and is_narrow[j + 1]:
                j += 1
            if j - i > best_j - best_i:
                best_i, best_j = i, j
            i = j + 1
        if best_i < 0:
            return None
        a_t = samples[best_i][0]
        b_t = samples[best_j][0]
        interval_len = b_t - a_t
        if interval_len < MIN_SEGMENT_LEN_M:
            return None
        margin = narrow_margin_frac * interval_len
        s_t = a_t + margin
        e_t = b_t - margin
        if e_t - s_t < MIN_SEGMENT_LEN_M:
            return None
        try:
            seg = substring(ls, s_t, e_t)
        except Exception:
            return None
        if (seg.geom_type != "LineString"
                or seg.is_empty
                or seg.length < MIN_SEGMENT_LEN_M):
            return None
        # Average hw within the narrow interval (for dedup scoring).
        narrow_hws = [h for (t, h) in samples
                      if best_i <= samples.index((t, h)) <= best_j
                      if 3.5 <= h <= NARROW_TAXI_HW_M]
        avg_hw = sum(narrow_hws) / len(narrow_hws) if narrow_hws else 0.0
        return seg, avg_hw

    # Group sub-ref lines; keep all other lines as-is.
    from collections import defaultdict
    sub_ref_lines: Dict[str, List[LineString]] = defaultdict(list)
    result: List[Tuple[LineString, str]] = []
    for ls, ref in centerlines:
        if ref and any(c.isdigit() for c in ref):
            sub_ref_lines[ref].append(ls)
        else:
            result.append((ls, ref))

    # For each sub-ref, pick the best slice.
    for ref, lines in sub_ref_lines.items():
        slices: List[Tuple[LineString, float]] = []
        for l in lines:
            r = _narrow_slice(l)
            if r is not None:
                slices.append(r)
        if not slices:
            # Fallback: keep the longest raw polyline.
            lines_sorted = sorted(lines, key=lambda l: -l.length)
            if lines_sorted and lines_sorted[0].length >= MIN_SEGMENT_LEN_M:
                result.append((lines_sorted[0], ref))
            continue
        # Prefer the LONGEST slice — the physical taxi corridor
        # typically has the longest continuous narrow interval.
        slices.sort(key=lambda sh: -sh[0].length)
        best = slices[0][0]
        result.append((best, ref))

    return result


def _split_centerlines_at_points(
    centerlines: List[Tuple[LineString, str]],
    split_points: List[Tuple[float, float]],
    approach_tol_m: float = 25.0,
    endpoint_guard_m: float = 5.0,
    pav_union: Optional[Polygon] = None,
    rwy_union: Optional[Polygon] = None,
    rwy_centerlines: Optional[List[LineString]] = None,
) -> List[Tuple[LineString, str]]:
    """Split each centerline at intersection points; emit between-
    break rects (15 % margin normally, 30 % for non-perpendicular
    taxis).

    Per user (2026-04-20 refined + 2026-04-21): rects are defined
    by intersections + sharp curves.  Two consecutive cut params
    merge into ONE junction region if the pavement at their
    midpoint is WIDER than the centerline's own narrow half-width
    (factor 1.2) — i.e. the pavement is widening in between
    (intersection widening).  Otherwise they remain separate
    junctions with a rect emitted between them.  This replaces
    the previous fixed ``CLOSE_INTERSECTION_M`` distance
    threshold, which was too coarse: 200 m was needed for
    OSM-fragmented V3 on V but merged real Q/R 3-rect splits too.

    Non-perpendicular taxis (45° stubs like V3) use a GAP_MARGIN_FRAC
    of 0.30 instead of 0.15 because the intersection point on the
    primary and on the runway sit farther down the taxi's own axis
    — without the larger margin the rect overlaps both junctions.
    """
    if not centerlines:
        return centerlines
    from shapely.ops import substring

    pav_for_probe = pav_union
    if pav_for_probe is not None and rwy_union is not None:
        try:
            pav_for_probe = pav_for_probe.union(rwy_union)
        except Exception:
            pass

    # Per user 2026-04-28: at sub-segment endpoints that are SHARED
    # with another sub-segment endpoint (i.e. the centerline was
    # bend-split there in ``_extract_osm_taxi_centerlines``), the
    # natural break is the bend itself.  A small junction polygon
    # at the bend is unavoidable (rect axes are straight; bent
    # centerlines need a junction at every angle change).  But the
    # default 15-30% margin on the segment overshoots the bend by
    # tens of metres, leaving a long uncovered corridor that the
    # downstream junction balloons into.  At bend-shared endpoints
    # use a small fixed margin (``BEND_ENDPOINT_MARGIN_M``) instead
    # of the percentage so the rect extends right up to the bend.
    BEND_ENDPOINT_MARGIN_M = 5.0
    BEND_SHARED_TOL_M = 25.0
    bend_share_tol2 = BEND_SHARED_TOL_M * BEND_SHARED_TOL_M
    centerline_endpoints: List[Tuple[Tuple[float, float],
                                     Tuple[float, float]]] = []
    for ls, _ref in centerlines:
        try:
            cs = list(ls.coords)
            centerline_endpoints.append((cs[0], cs[-1]))
        except Exception:
            centerline_endpoints.append(((0.0, 0.0), (0.0, 0.0)))

    def _is_bend_shared(idx: int, endpoint: Tuple[float, float]) -> bool:
        """True iff ``endpoint`` of centerline ``idx`` lies within
        ``BEND_SHARED_TOL_M`` of any other centerline's endpoint —
        signalling that the two centerlines were bend-split apart
        from one continuous OSM way at that point."""
        ex, ey = endpoint
        for j, (s2, e2) in enumerate(centerline_endpoints):
            if j == idx:
                continue
            for px, py in (s2, e2):
                if (ex - px) ** 2 + (ey - py) ** 2 <= bend_share_tol2:
                    return True
        return False

    def _avg_perp_halfwidth(ls: LineString, t: float) -> float:
        """(left+right)/2 perpendicular half-width at axis param t,
        so widening detection is comparable to narrow_hw (which is
        also derived from (left+right)/2 per-probe averages)."""
        if pav_for_probe is None or pav_for_probe.is_empty:
            return 0.0
        RAY_CAP_M = 40.0
        RAY_STEP_M = 0.5
        dt = min(2.0, ls.length * 0.05)
        a = ls.interpolate(max(0.0, t - dt))
        b = ls.interpolate(min(ls.length, t + dt))
        tx, ty = b.x - a.x, b.y - a.y
        mag = math.hypot(tx, ty)
        if mag < 1e-6:
            return 0.0
        ux, uy = tx / mag, ty / mag
        nx, ny = -uy, ux
        pt = ls.interpolate(t)
        sides: List[float] = []
        for sign in (-1, 1):
            side = RAY_CAP_M
            d = 0.0
            while d <= RAY_CAP_M:
                qx = pt.x + sign * nx * d
                qy = pt.y + sign * ny * d
                if not pav_for_probe.contains(Point(qx, qy)):
                    side = d
                    break
                d += RAY_STEP_M
            sides.append(side)
        return sum(sides) / 2.0 if sides else 0.0

    def _rect_margin_frac_for(ls: LineString, ref: str) -> float:
        # Stubs / cross-connectors oriented > 30° off perpendicular
        # to the nearest runway get a larger margin because the
        # intersection points (primary and runway) sit farther along
        # the taxi's axis due to the oblique crossing.
        if not rwy_centerlines:
            return 0.15
        c = list(ls.coords)
        if len(c) < 2:
            return 0.15
        dx = c[-1][0] - c[0][0]
        dy = c[-1][1] - c[0][1]
        mag = math.hypot(dx, dy)
        if mag < 1e-6:
            return 0.15
        axis_bearing = math.degrees(math.atan2(dx, dy)) % 180.0
        # Nearest runway to centerline mid
        mid = ls.interpolate(ls.length / 2)
        best_r = None
        best_d = float("inf")
        for r in rwy_centerlines:
            d = mid.distance(r)
            if d < best_d:
                best_d = d
                best_r = r
        if best_r is None:
            return 0.15
        rc = list(best_r.coords)
        if len(rc) < 2:
            return 0.15
        rx = rc[-1][0] - rc[0][0]
        ry = rc[-1][1] - rc[0][1]
        rmag = math.hypot(rx, ry)
        if rmag < 1e-6:
            return 0.15
        rwy_bearing = math.degrees(math.atan2(rx, ry)) % 180.0
        delta = abs(axis_bearing - rwy_bearing)
        delta = min(delta, 180.0 - delta)
        # Perpendicular = 90°.  Taxis that connect to the runway
        # AT AN ANGLE (not parallel, not perpendicular) get the
        # diagonal-stub treatment: 35 % rect length biased 25 %
        # of the gap toward the runway-facing endpoint.  Primary
        # parallels (delta ≈ 0°, perp_diff ≈ 90°) stay at 15 %;
        # perpendicular cross-connectors (delta ≈ 90°,
        # perp_diff ≈ 0°) stay at 15 %.  B/C/E/G at SPJC measure
        # perp_diff ≈ 69° — just outside the old 25-65 window —
        # so widen to 20-75 to cover the full "at an angle"
        # band while still excluding pure parallels and
        # perpendiculars.
        #
        # Per user 2026-04-27: this used to be gated on length
        # < 250 m (long diagonals fell back to 15 %).  But long
        # diagonals like SPJC's B/C/E (575 m / 572 m / 593 m at
        # ~21° to runway) ARE diagonal stubs in the same sense as
        # short ones — they should get the same 35 % retention so
        # the rect doesn't extend deep into the adjacent apron.
        # Removed the length gate; the angle alone classifies.
        perp_diff = abs(delta - 90.0)
        if 20.0 < perp_diff < 75.0:
            # Per user 2026-04-28: the 30 % diagonal-stub margin is
            # appropriate ONLY when at least one endpoint sits AT a
            # runway boundary — i.e. the segment IS the stub between
            # a parallel taxi and the runway.  When NEITHER endpoint
            # is near a runway, the segment is a non-stub diagonal
            # connector (e.g. CYXY E nodes 2-4 transitioning between
            # the south-of-apron parallel section and the apron-
            # internal parallel section, perp_diff ≈ 60° but both
            # ends far from any runway).  Such segments shouldn't
            # lose 60 % of their length to junction-margin trim;
            # they're not bordered by junctions on both sides.
            STUB_ENDPOINT_RUNWAY_M = 50.0
            ep0 = Point(c[0])
            ep1 = Point(c[-1])
            ep0_near = any(
                ep0.distance(r) <= STUB_ENDPOINT_RUNWAY_M
                for r in rwy_centerlines)
            ep1_near = any(
                ep1.distance(r) <= STUB_ENDPOINT_RUNWAY_M
                for r in rwy_centerlines)
            if ep0_near or ep1_near:
                return 0.30
            # Neither endpoint near a runway — treat as a long
            # diagonal connector with the parallel-style 15 % margin.
        # Short unrefed parallel-to-runway rect (perp_diff >= 75°,
        # length < 150 m) sitting between two diagonal stubs on
        # SPLP's south chain — apply 30 % margin each side so the
        # resulting primary rect is half its default length and
        # doesn't overlap the adjacent diagonals.
        if not ref and ls.length < 150.0 and perp_diff >= 75.0:
            return 0.30
        return 0.15

    result: List[Tuple[LineString, str]] = []
    for ls_idx, (ls, ref) in enumerate(centerlines):
        gap_margin_frac = _rect_margin_frac_for(ls, ref)
        # Detect whether the centerline's start / end is a bend-shared
        # endpoint (continuation with a neighbouring sub-segment via
        # bend).  Used below to clamp the margin at those endpoints.
        try:
            _ls_cs = list(ls.coords)
            _start_endpoint = _ls_cs[0]
            _end_endpoint = _ls_cs[-1]
        except Exception:
            _start_endpoint = (0.0, 0.0)
            _end_endpoint = (0.0, 0.0)
        start_is_bend = _is_bend_shared(ls_idx, _start_endpoint)
        end_is_bend = _is_bend_shared(ls_idx, _end_endpoint)
        # Estimate centerline's narrow half-width for midpoint check.
        if pav_for_probe is not None and not pav_for_probe.is_empty:
            _nat, _p90, narrow_hw = _natural_half_width(ls, pav_for_probe)
        else:
            narrow_hw = 0.0

        # Collect cut params for intersections that lie on this line.
        cut_params: List[float] = []
        for (sx, sy) in split_points or ():
            sp = Point(sx, sy)
            if ls.distance(sp) > approach_tol_m:
                continue
            try:
                param = ls.project(sp)
            except Exception:
                continue
            if param < endpoint_guard_m:
                continue
            if param > ls.length - endpoint_guard_m:
                continue
            cut_params.append(param)

        cut_params.sort()
        # Pav-width midpoint cluster: two consecutive cut_params
        # merge when the midpoint half-width > narrow_hw × 1.2.
        # Always merge when they're within 25 m (same-crossing
        # multi-node noise).  Never merge past 400 m apart.
        clusters: List[List[float]] = []
        WIDEN_FACTOR = 1.2
        MIN_ALWAYS_MERGE = 25.0
        MAX_CLUSTER_SPAN_M = 400.0
        for p in cut_params:
            if not clusters:
                clusters.append([p])
                continue
            prev = clusters[-1][-1]
            gap = p - prev
            merge = False
            if gap <= MIN_ALWAYS_MERGE:
                merge = True
            elif gap > MAX_CLUSTER_SPAN_M:
                merge = False
            elif narrow_hw > 0:
                try:
                    mid_hw = _avg_perp_halfwidth(ls, (prev + p) / 2.0)
                    if mid_hw > narrow_hw * WIDEN_FACTOR:
                        merge = True
                except Exception:
                    pass
            if merge:
                clusters[-1].append(p)
            else:
                clusters.append([p])

        breaks: List[float] = [0.0]
        for cl in clusters:
            breaks.append(cl[0])
            breaks.append(cl[-1])
        breaks.append(ls.length)

        # Enumerate candidate segments and identify which are the
        # first/last ones that would actually emit.  Cross-connector
        # taxis (perpendicular to the runway, terminating into wider
        # parallel taxis) use 30 % margin each side on first/last
        # emitted segments because the parallel-taxi widening zone
        # extends into the cross-connector's axis.  Detected by
        # geometry: bearing within 20° of perpendicular to nearest
        # runway, axis midpoint > 250 m from runway centerline.
        n_breaks = len(breaks)
        # Collect only segments that will actually emit (post-margin
        # length >= 40 m under ANY margin we might apply, so first/
        # last indexing is stable).  The 40 m floor matches the
        # post-emit filter below.  Use smallest possible retained
        # fraction (35 % for diagonal, 40 % for cross-connector
        # ends) to test.
        candidates: List[Tuple[float, float]] = []
        for i in range(0, n_breaks - 1, 2):
            p0, p1 = breaks[i], breaks[i + 1]
            gap = p1 - p0
            if gap < MIN_SEGMENT_LEN_M:
                continue
            # Will this segment emit under any plausible margin?
            min_retained_frac = 0.35
            if gap * min_retained_frac < 40.0 and gap * (1 - 2 * gap_margin_frac) < 40.0:
                # Won't emit — skip so it doesn't shift end-indexing.
                continue
            candidates.append((p0, p1))
        # Cross-connector detection: full centerline is perpendicular
        # (within 20°) to nearest runway AND its midpoint is > 250 m
        # from any runway centerline (i.e. it's a connector BETWEEN
        # parallels, not a runway-touching stub).
        is_cross = False
        if rwy_centerlines:
            try:
                lc = list(ls.coords)
                if len(lc) >= 2:
                    ldx = lc[-1][0] - lc[0][0]
                    ldy = lc[-1][1] - lc[0][1]
                    lmag = math.hypot(ldx, ldy)
                    if lmag > 1e-6:
                        l_bearing = math.degrees(
                            math.atan2(ldx, ldy)) % 180.0
                        mid = ls.interpolate(ls.length / 2.0)
                        rmid = min(rwy_centerlines,
                                   key=lambda r: mid.distance(r))
                        rcc = list(rmid.coords)
                        if len(rcc) >= 2:
                            rdx2 = rcc[-1][0] - rcc[0][0]
                            rdy2 = rcc[-1][1] - rcc[0][1]
                            rmag2 = math.hypot(rdx2, rdy2)
                            if rmag2 > 1e-6:
                                r_bearing = math.degrees(
                                    math.atan2(rdx2, rdy2)) % 180.0
                                d = abs(l_bearing - r_bearing)
                                d = min(d, 180.0 - d)
                                d_perp = abs(d - 90.0)
                                d_to_rwy = mid.distance(rmid)
                                if d_perp <= 20.0 and d_to_rwy > 250.0:
                                    is_cross = True
            except Exception:
                pass
        for idx, (p0, p1) in enumerate(candidates):
            gap = p1 - p0
            is_end_seg = (idx == 0 or idx == len(candidates) - 1)
            # Per-segment parallel check: a SHORT SLICE of a long
            # unrefed curving taxi can be parallel-to-runway even
            # when the full centerline's orientation is diagonal.
            # Classify such a slice for the "short primary between
            # diagonals" shrinkage (30 % margin each side, no bias)
            # — fixes SPLP's (-500,-1054) which sits between two
            # diagonal stubs and would otherwise overlap them.
            is_short_parallel_slice = False
            if (not ref and gap < 150.0 and rwy_centerlines):
                try:
                    seg_a = ls.interpolate(p0)
                    seg_b = ls.interpolate(p1)
                    sdx = seg_b.x - seg_a.x
                    sdy = seg_b.y - seg_a.y
                    smag = math.hypot(sdx, sdy)
                    if smag > 1e-6:
                        seg_bearing = math.degrees(
                            math.atan2(sdx, sdy)) % 180.0
                        rseg_best = min(
                            rwy_centerlines,
                            key=lambda r: ls.interpolate(
                                (p0 + p1) / 2.0).distance(r))
                        rseg = list(rseg_best.coords)
                        rdx = rseg[-1][0] - rseg[0][0]
                        rdy = rseg[-1][1] - rseg[0][1]
                        rmag2 = math.hypot(rdx, rdy)
                        if rmag2 > 1e-6:
                            seg_rwy_bearing = math.degrees(
                                math.atan2(rdx, rdy)) % 180.0
                            seg_delta = abs(
                                seg_bearing - seg_rwy_bearing)
                            seg_delta = min(
                                seg_delta, 180.0 - seg_delta)
                            seg_perp = abs(seg_delta - 90.0)
                            if seg_perp >= 75.0:
                                is_short_parallel_slice = True
                except Exception:
                    pass
            if is_short_parallel_slice:
                # Shrink to half length, no bias (it's parallel,
                # not diagonal).  Using 22 % margin each side
                # (56 % retained) keeps the piece above the 40 m
                # post-margin floor even for short 60-70 m slices
                # of SPLP's main taxi.
                m_start = 0.22 * gap
                m_end = 0.22 * gap
            elif is_cross and is_end_seg:
                m_start = 0.30 * gap
                m_end = 0.30 * gap
            elif gap_margin_frac >= 0.25:
                # Non-perpendicular diagonal stub (V3-like): 32.5 %
                # margin each side (35 % retained), biased 25 % of
                # the gap toward the runway-facing axis endpoint.
                retained = 0.35 * gap
                remaining_margin = gap - retained
                # Bias: shift the rect center 20 % of the gap toward
                # the endpoint nearer the runway.  User (2026-04-21):
                # adjusted 5 % back from the original 25 % because
                # biasing too close to the runway made diagonal stubs
                # encroach on the runway-ramp widening zone.
                bias = 0.20 * gap
                # Which endpoint is closer to a runway?
                if rwy_centerlines:
                    ep0 = ls.interpolate(p0)
                    ep1 = ls.interpolate(p1)
                    d0 = min(ep0.distance(r) for r in rwy_centerlines)
                    d1 = min(ep1.distance(r) for r in rwy_centerlines)
                    runway_is_p0_side = d0 < d1
                else:
                    runway_is_p0_side = True
                if runway_is_p0_side:
                    m_start = remaining_margin / 2.0 - bias
                    m_end = remaining_margin / 2.0 + bias
                else:
                    m_start = remaining_margin / 2.0 + bias
                    m_end = remaining_margin / 2.0 - bias
                # Clamp margins to be non-negative
                m_start = max(0.0, m_start)
                m_end = max(0.0, m_end)
                # Recompute retained so p0+m_start..p1-m_end fits
                if gap - m_start - m_end < MIN_SEGMENT_LEN_M:
                    continue
            else:
                m_start = gap_margin_frac * gap
                m_end = gap_margin_frac * gap
            # Per user 2026-04-28: at bend-shared centerline endpoints,
            # the rect should extend right up to the bend (only the
            # tiny natural triangular junction at the angle change is
            # unavoidable).  Override the percentage margin with a
            # small fixed value when the corresponding endpoint of
            # this segment touches a bend-shared end of the
            # centerline.
            #
            # BUT: cap the extension at the point where the corridor
            # widens past 1.3 × narrow_hw — past that the rect would
            # extend deep into an apron, fail the apron-interior
            # check in ``_build_taxi_rects`` (≥ 2 corners off-
            # boundary), and never be emitted (e.g. CYXY's North F
            # bend-extends 57 m into an apron).  Walk inward from
            # the centerline's end probing the half-width; stop
            # where the corridor is back to within 1.3 × narrow_hw.
            CORRIDOR_WIDTH_FACTOR = 1.3
            def _bend_margin_at(end_param: float, sign: int) -> float:
                """``end_param`` = 0 (start) or ls.length (end);
                ``sign`` = +1 (walk forward into the line) or -1
                (walk backward).  Returns a margin in metres at
                least ``BEND_ENDPOINT_MARGIN_M`` and at most the
                point where the corridor narrows back to
                ``CORRIDOR_WIDTH_FACTOR × narrow_hw``."""
                base = BEND_ENDPOINT_MARGIN_M
                if narrow_hw <= 0:
                    return base
                target_hw = narrow_hw * CORRIDOR_WIDTH_FACTOR
                STEP = 5.0
                MAX = max(base, gap / 2.0)
                u = base
                while u <= MAX:
                    t = end_param + sign * u
                    if t < 0 or t > ls.length:
                        break
                    try:
                        hw_here = _avg_perp_halfwidth(ls, t)
                    except Exception:
                        hw_here = 0.0
                    if 0 < hw_here <= target_hw:
                        return u
                    u += STEP
                # Corridor never narrowed to ≤ target_hw within
                # half the gap — fall back to the percentage margin
                # so the rect doesn't extend into apron territory.
                return float('inf')
            if start_is_bend and abs(p0) < 0.5:
                bm = _bend_margin_at(0.0, +1)
                if bm != float('inf'):
                    m_start = min(m_start, bm)
            if end_is_bend and abs(p1 - ls.length) < 0.5:
                bm = _bend_margin_at(ls.length, -1)
                if bm != float('inf'):
                    m_end = min(m_end, bm)
            rect_p0 = p0 + m_start
            rect_p1 = p1 - m_end
            if rect_p1 - rect_p0 < MIN_SEGMENT_LEN_M:
                continue
            try:
                piece = substring(ls, rect_p0, rect_p1)
            except Exception:
                continue
            # Drop short between-junction fragments (< 40 m).  Target
            # cross_connector smallest = 52 m, Q smallest = 59 m,
            # so 40 m post-margin is a safe floor that still drops
            # spurious junction-approach tails (e.g. R switchback
            # tails at the V/Q/R triple junction).
            if (piece.geom_type == "LineString"
                    and not piece.is_empty
                    and piece.length >= 40.0):
                result.append((piece, ref))
    return result


def _build_taxi_rects(
    centerlines: List[Tuple[LineString, str]],
    pav_union: Optional[Polygon],
    rwy_union: Optional[Polygon],
    rwy_centerlines: List[LineString],
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
    ref_overall_bearings: Optional[Dict[str, float]] = None,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Convert each usable centerline into a 4-vertex rect.

    For each centerline we:
      1. Probe the half-width (distance to apt.dat pavement boundary)
         at many points along the axis.
      2. Determine the ``natural half-width`` of the strip as the
         median of those probes.
      3. TRIM axis endpoints inward until the probe there is
         ≤ 1.3 × natural_half_width — this is the user's rule of
         "pull rects back to the narrowest part of the taxiway."
         Everything past the trim is widening territory, reserved
         for junction polygons.
      4. Emit a 4-vertex rect over the trimmed axis with width =
         2 × natural_half_width.  The rect's 4 corners sit at the
         trimmed axis endpoints ± perpendicular half-width.

    Returns list of (rect, clipped_axis, role, ref).
    """
    if pav_union is None:
        return []

    pav_non_rwy = pav_union
    if rwy_union is not None:
        pav_non_rwy = pav_non_rwy.difference(rwy_union)

    centerlines = sorted(centerlines, key=lambda x: -x[0].length)

    emitted: List[Tuple[Polygon, LineString, str, str]] = []
    emitted_union: Optional[Polygon] = None

    for axis, ref in centerlines:
        # Clip to the full pavement (including runway).  The rect may
        # overlap runway slightly at stubs that reach the runway edge;
        # we prefer a stub that correctly reaches the runway over
        # clipping it at the runway boundary and losing most of it.
        try:
            clipped = axis.intersection(pav_union)
        except Exception:
            continue
        if clipped.is_empty:
            continue
        if clipped.geom_type == "MultiLineString":
            clipped = max(clipped.geoms, key=lambda g: g.length)
        if clipped.geom_type != "LineString":
            continue
        if clipped.length < 20.0:
            continue

        # Probe half-widths along axis.  ``narrow_hw`` (p10) is the
        # strip's narrowest portion; we use it both as the trim
        # baseline AND as the rect's half-width per the user's
        # authoritative rule.  The factor 1.15 means trim ends
        # where pavement grows more than 15% wider than the rect —
        # that extra width is junction territory.
        _natural_hw, _max_hw, narrow_hw = _natural_half_width(
            clipped, pav_non_rwy)
        if narrow_hw < 3.5 or narrow_hw > 40.0:
            continue

        # Per user (2026-04-20): trimming is now handled by the
        # upstream _split_centerlines_at_points which emits each
        # rect axis at 70% of the distance between intersections
        # (15% margin on each junction-facing end).  Skip the
        # width-based trim here — the axis is already cut to the
        # rect's intended length.
        trimmed = clipped
        trim_narrow_hw = narrow_hw

        # Dedup against trimmed axis
        if emitted_union is not None and not emitted_union.is_empty:
            try:
                inside_len = trimmed.intersection(emitted_union).length
                if inside_len / trimmed.length > 0.7:
                    continue
            except Exception:
                pass

        width = 2.0 * trim_narrow_hw
        rect = _rect_from_axis_extended(trimmed, width, pav_non_rwy,
                                        apt_vertices=apt_vertices)
        if rect is None or rect.is_empty:
            continue
        # Skip invalid rects (self-intersecting after snap).
        if not rect.is_valid:
            try:
                rect = rect.buffer(0)
            except Exception:
                continue
            if (rect.is_empty or rect.geom_type != "Polygon"
                    or not rect.is_valid):
                continue

        # ── Apron-interior rect rejection (user 2026-04-27) ────────
        # A rect whose 4 corners aren't on (or very near) the
        # pavement boundary is sitting INSIDE an apron — the
        # surrounding pavement wraps around it, downstream junction
        # construction has to wrap a junction around it too, and the
        # junction's elevation has to bridge the rect's slope on
        # both long edges → visible elevation ridges in JOSM and at
        # render time.
        #
        # Real taxi rects have their 4 corners at pavement-boundary
        # points (the intersections where the taxi corridor meets
        # the adjacent apron / parallel / runway).  If the rect's
        # corners are well INSIDE the pavement, the centerline runs
        # through an apron and shouldn't emit a separate rect — the
        # apron stays as one continuous junction.
        BOUNDARY_TOL_M = 2.0
        rect_coords = list(rect.exterior.coords)
        if rect_coords and rect_coords[0] == rect_coords[-1]:
            rect_coords = rect_coords[:-1]
        boundary = pav_union.boundary
        n_off_boundary = sum(
            1 for (cx, cy) in rect_coords
            if Point(cx, cy).distance(boundary) > BOUNDARY_TOL_M)
        if n_off_boundary >= 2:
            # ≥ 2 corners away from any pavement edge — the rect
            # sits inside an apron.  Skip it; the apron pavement
            # stays as residue → junction.
            continue

        role = _classify_role(trimmed, width, rwy_centerlines,
                               rwy_union, ref=ref,
                               ref_overall_bearings=ref_overall_bearings)
        emitted.append((rect, trimmed, role, ref))
        try:
            emitted_union = (unary_union([emitted_union, rect])
                             if emitted_union is not None else rect)
        except Exception:
            # Self-intersection of accumulated union — skip update.
            pass

    # Post-classify: secondary passes to fix roles based on neighbour
    # topology (stubs touch parallels, cross_connector touches 2 parallels)
    _refine_roles(emitted, rwy_centerlines)

    # Stub-ref dedup: OSM often has multiple disjoint ways with the
    # same sub-ref label (e.g. V2 has 3 separate OSM pieces, each
    # contributing a rect).  Target emits ONE rect per stub ref.
    # Keep only the LONGEST rect per stub ref.  Applies to both
    # sub-refs (letter+digit like V1, L3) AND letter-only stubs
    # (B, C, E, G) which are single-stub taxis that must not
    # fragment across internal bends.
    def _should_dedup(ref_str: str, role_str: str) -> bool:
        if not ref_str:
            return False
        if any(c.isdigit() for c in ref_str):
            return True
        # Letter-only: dedup when classified as stub (B, C, E, G, D).
        if role_str == ROLE_STUB:
            return True
        return False
    # Per-ref, group same-ref rects into geometrically-overlapping
    # clusters; within each cluster, keep only the longest.  Non-
    # overlapping rects of the same ref (e.g. multiple sub-segments
    # of one OSM way that bend-split into rects covering distinct
    # parts of the corridor) coexist — only OSM fragmentation that
    # produces actual duplicates gets deduped.
    OVERLAP_PROX_M = 5.0
    by_ref: Dict[str, List[int]] = {}
    for i, (_r, _a, role, ref) in enumerate(emitted):
        if not _should_dedup(ref, role):
            continue
        by_ref.setdefault(ref, []).append(i)
    drop: set = set()
    for ref, idxs in by_ref.items():
        # Build overlap clusters.
        n = len(idxs)
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a in range(n):
            ra = emitted[idxs[a]][0]
            for b in range(a + 1, n):
                rb = emitted[idxs[b]][0]
                try:
                    overlap = ra.intersection(rb).area
                    proximate = ra.distance(rb) < OVERLAP_PROX_M
                except Exception:
                    overlap, proximate = 0.0, False
                if overlap > 1.0 or proximate:
                    pa, pb = find(a), find(b)
                    if pa != pb:
                        parent[pa] = pb
        clusters: Dict[int, List[int]] = {}
        for k in range(n):
            clusters.setdefault(find(k), []).append(idxs[k])
        # Within each cluster, keep only the longest axis.
        for members in clusters.values():
            if len(members) <= 1:
                continue
            members.sort(key=lambda m: -emitted[m][1].length)
            for m in members[1:]:
                drop.add(m)
    keep: List[Tuple[Polygon, LineString, str, str]] = [
        item for i, item in enumerate(emitted) if i not in drop]
    return keep


def _merge_collinear_rects_principled(
    emitted: List[Tuple[Polygon, LineString, str, str]],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
    angle_tol_deg: float = 4.0,
    gap_tol_m: float = 12.0,
    width_uniformity_tol: float = 1.10,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Merge adjacent same-ref rects whose joining point shows
    NO widening — the pavement runs straight at uniform narrow
    width through the joint.  This is the only case where "single
    rect between junctions on straight sections" applies.
    """
    if not emitted or pav is None or pav.is_empty:
        return emitted
    boundary = pav.boundary
    changed = True
    work = list(emitted)
    while changed:
        changed = False
        for i in range(len(work)):
            for j in range(i + 1, len(work)):
                ri, ai, roli, refi = work[i]
                rj, aj, rolj, refj = work[j]
                if refi != refj or refi == "":
                    continue  # refless airports handled separately
                if roli != rolj:
                    continue
                # Bearings (mod 180°).
                def _bearing(a):
                    c = list(a.coords)
                    return math.degrees(
                        math.atan2(c[-1][0] - c[0][0],
                                   c[-1][1] - c[0][1])) % 180.0
                bi = _bearing(ai)
                bj = _bearing(aj)
                db = abs(bi - bj)
                db = min(db, 180.0 - db)
                if db > angle_tol_deg:
                    continue
                # Closest endpoints & far endpoints.
                coords_i = list(ai.coords)
                coords_j = list(aj.coords)
                pairs = [
                    (coords_i[0], coords_j[0], 0, 0),
                    (coords_i[0], coords_j[-1], 0, 1),
                    (coords_i[-1], coords_j[0], 1, 0),
                    (coords_i[-1], coords_j[-1], 1, 1),
                ]
                best = min(pairs, key=lambda p: math.hypot(
                    p[0][0] - p[1][0], p[0][1] - p[1][1]))
                endp_i, endp_j, ei, ej = best
                gap = math.hypot(endp_i[0] - endp_j[0],
                                 endp_i[1] - endp_j[1])
                if gap > gap_tol_m:
                    continue
                # Joining-region pavement half-width: probe at the
                # midpoint of the two touching endpoints.
                mid = Point((endp_i[0] + endp_j[0]) / 2,
                            (endp_i[1] + endp_j[1]) / 2)
                hw_joint = mid.distance(boundary) if pav.contains(mid) else 0
                if hw_joint <= 0:
                    continue
                # Each rect's own half-width (MRR short side / 2).
                def _rect_hw(p):
                    mrr = p.minimum_rotated_rectangle
                    c = list(mrr.exterior.coords)
                    if len(c) < 5:
                        return 0.0
                    s1 = math.hypot(c[1][0] - c[0][0], c[1][1] - c[0][1])
                    s2 = math.hypot(c[2][0] - c[1][0], c[2][1] - c[1][1])
                    return min(s1, s2) / 2.0
                hwi = _rect_hw(ri)
                hwj = _rect_hw(rj)
                if hwi <= 0 or hwj <= 0:
                    continue
                # No widening: joint hw is within uniformity_tol of
                # each rect's own hw (equivalently, joint hw ≤
                # max(hwi, hwj) × uniformity_tol AND rects have
                # similar widths).
                max_hw = max(hwi, hwj)
                if hw_joint > max_hw * width_uniformity_tol:
                    continue
                ratio_ij = max(hwi, hwj) / min(hwi, hwj)
                if ratio_ij > width_uniformity_tol:
                    continue
                # All checks pass — merge.
                far_i = coords_i[0] if ei == 1 else coords_i[-1]
                far_j = coords_j[0] if ej == 1 else coords_j[-1]
                try:
                    merged_axis = LineString([far_i, far_j])
                except Exception:
                    continue
                merged_rect = _rect_from_axis_extended(
                    merged_axis, 2.0 * (hwi + hwj) / 2.0, pav,
                    apt_vertices=apt_vertices)
                if merged_rect is None or merged_rect.is_empty:
                    continue
                work[i] = (merged_rect, merged_axis, roli, refi)
                del work[j]
                changed = True
                break
            if changed:
                break
    return work


def _merge_collinear_rects(
    emitted: List[Tuple[Polygon, LineString, str, str]],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
    angle_tol_deg: float = 4.0,
    gap_tol_m: float = 8.0,
    width_ratio_tol: float = 1.15,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Merge adjacent same-ref rects whose axes are nearly
    collinear and whose meeting point sits in the narrow corridor
    (no widening).  Produces one rect per straight pavement
    section between real junctions.
    """
    if not emitted:
        return emitted
    # Group by ref and role so we only merge within a ref's rects.
    merged = True
    work = list(emitted)
    while merged:
        merged = False
        n = len(work)
        for i in range(n):
            for j in range(i + 1, n):
                ri, ai, roli, refi = work[i]
                rj, aj, rolj, refj = work[j]
                # Only merge same ref + compatible role.
                if refi != refj or roli != rolj:
                    continue
                if refi == "" or refj == "":
                    continue  # refless: don't auto-merge (SPLP)
                # Axis bearings (mod 180°).
                def _bearing(a):
                    c = list(a.coords)
                    return math.degrees(
                        math.atan2(c[-1][0] - c[0][0],
                                   c[-1][1] - c[0][1])) % 180.0
                bi = _bearing(ai)
                bj = _bearing(aj)
                db = abs(bi - bj)
                db = min(db, 180.0 - db)
                if db > angle_tol_deg:
                    continue
                # Closest endpoints between axes.
                coords_i = list(ai.coords)
                coords_j = list(aj.coords)
                pairs = [
                    (coords_i[0], coords_j[0], 0, 0),
                    (coords_i[0], coords_j[-1], 0, 1),
                    (coords_i[-1], coords_j[0], 1, 0),
                    (coords_i[-1], coords_j[-1], 1, 1),
                ]
                best = min(pairs, key=lambda p: math.hypot(
                    p[0][0] - p[1][0], p[0][1] - p[1][1]))
                endp_i, endp_j, ei, ej = best
                gap = math.hypot(endp_i[0] - endp_j[0],
                                 endp_i[1] - endp_j[1])
                if gap > gap_tol_m:
                    continue
                # Width similarity: compare rect widths (from their
                # polygons' minimum-rotated-rect short-side).
                def _rect_width(p):
                    mrr = p.minimum_rotated_rectangle
                    coords = list(mrr.exterior.coords)
                    if len(coords) < 5:
                        return 0.0
                    s1 = math.hypot(coords[1][0] - coords[0][0],
                                    coords[1][1] - coords[0][1])
                    s2 = math.hypot(coords[2][0] - coords[1][0],
                                    coords[2][1] - coords[1][1])
                    return min(s1, s2)
                wi = _rect_width(ri)
                wj = _rect_width(rj)
                if wi < 1.0 or wj < 1.0:
                    continue
                if max(wi, wj) / min(wi, wj) > width_ratio_tol:
                    continue
                # Build merged axis: take the 2 FAR endpoints.
                far_i = coords_i[0] if ei == 1 else coords_i[-1]
                far_j = coords_j[0] if ej == 1 else coords_j[-1]
                try:
                    merged_axis = LineString([far_i, far_j])
                except Exception:
                    continue
                # Build merged rect from merged axis at average width.
                avg_width = (wi + wj) / 2.0
                merged_rect = _rect_from_axis_extended(
                    merged_axis, avg_width, pav,
                    apt_vertices=apt_vertices)
                if merged_rect is None or merged_rect.is_empty:
                    continue
                # Replace i with merged, drop j.
                work[i] = (merged_rect, merged_axis, roli, refi)
                del work[j]
                merged = True
                break
            if merged:
                break
    return work


def _natural_half_width(axis: LineString, pav: Polygon,
                        n_probes: int = 15) -> Tuple[float, float, float]:
    """Return (natural_hw, max_hw, narrow_hw) LOCAL half-width probes
    along the axis.

    Uses PERPENDICULAR RAY CAST (not distance-to-boundary) so the
    probe measures the taxi's own local half-width on EACH side
    rather than the distance to some faraway edge.  Per user
    rule 4 (2026-04-20): the rect half-width should be the
    NARROWEST pavement width (the taxi's own strip width), not
    an inflated value from adjacent aprons or runway clearance.

    For each probe point:
      * cast a ray perpendicular LEFT from axis; find where ray
        first exits the pavement polygon.
      * cast a ray perpendicular RIGHT from axis similarly.
      * half-width at this probe = min(left, right), capped at
        RAY_CAP_M to avoid saturating across an apron.
    """
    RAY_CAP_M = 40.0
    RAY_STEP_M = 0.5
    if axis.length < 1e-3:
        return 0.0, 0.0, 0.0

    def _perpendicular_half_at(t: float) -> float:
        """Cast perpendicular rays left/right at axis param t.

        Returns the AVERAGE of the two sides — (left + right) / 2 —
        so the width reflects the full pavement strip centered on
        the pavement (not a narrow corridor seen from an off-center
        axis).  Corner snapping downstream pulls the 4 rect corners
        onto the pav boundary, centering the rect on the actual
        pavement regardless of the axis's offset.
        """
        # Local tangent: use points slightly before/after t.
        dt = min(2.0, axis.length * 0.05)
        t0 = max(0.0, t - dt)
        t1 = min(axis.length, t + dt)
        a = axis.interpolate(t0)
        b = axis.interpolate(t1)
        tx, ty = b.x - a.x, b.y - a.y
        mag = math.hypot(tx, ty)
        if mag < 1e-6:
            return 0.0
        ux, uy = tx / mag, ty / mag
        nx, ny = -uy, ux  # left-perp
        pt = axis.interpolate(t)
        ox, oy = pt.x, pt.y
        sides: List[float] = []
        for sign in (-1, 1):
            side = RAY_CAP_M
            d = 0.0
            while d <= RAY_CAP_M:
                qx = ox + sign * nx * d
                qy = oy + sign * ny * d
                if not pav.contains(Point(qx, qy)):
                    side = d
                    break
                d += RAY_STEP_M
            sides.append(side)
        return sum(sides) / 2.0 if sides else RAY_CAP_M

    dists: List[float] = []
    for k in range(n_probes):
        t = (k + 1) / (n_probes + 1) * axis.length
        hw = _perpendicular_half_at(t)
        if hw > 0.1:
            dists.append(hw)
    if not dists:
        return 0.0, 0.0, 0.0
    dists.sort()
    median = dists[len(dists) // 2]
    p90_idx = max(0, int(len(dists) * 0.9) - 1)
    p90 = dists[p90_idx] if p90_idx < len(dists) else dists[-1]
    # ``narrow`` = the MIN half-width probe with a floor to guard
    # against grazing a building corner (skip probes < 3.5 m as
    # noise).  Per user rule 4: rect width = the ACTUAL narrowest
    # section, so the rect fits snugly along the taxi's own narrow
    # corridor, leaving widened areas to junctions.
    filtered = [d for d in dists if d >= 3.5]
    narrow = filtered[0] if filtered else dists[0]
    return median, p90, narrow


def _trim_to_narrow(axis: LineString, pav: Polygon, natural_hw: float,
                    widen_factor: float = 1.3) -> Optional[LineString]:
    """Trim the axis inward from each end until the PERPENDICULAR
    half-width at the endpoint drops below ``widen_factor × natural_hw``.

    Uses the same perpendicular ray-cast probing as
    ``_natural_half_width`` so the trim threshold is applied to
    the taxi's LOCAL half-width (not distance-to-boundary).

    Trim step is 2 m.  We never trim more than 50 % of the axis
    length.
    """
    RAY_CAP_M = 40.0
    RAY_STEP_M = 0.5
    total_len = axis.length
    thresh = natural_hw * widen_factor
    step = 2.0
    max_trim = total_len * 0.45

    def _perp_hw(t: float) -> float:
        dt = min(2.0, total_len * 0.05)
        t0 = max(0.0, t - dt)
        t1 = min(total_len, t + dt)
        a = axis.interpolate(t0)
        b = axis.interpolate(t1)
        tx, ty = b.x - a.x, b.y - a.y
        mag = math.hypot(tx, ty)
        if mag < 1e-6:
            return 0.0
        ux, uy = tx / mag, ty / mag
        nx, ny = -uy, ux
        pt = axis.interpolate(t)
        ox, oy = pt.x, pt.y
        best = RAY_CAP_M
        for sign in (-1, 1):
            d = 0.0
            while d <= RAY_CAP_M:
                qx = ox + sign * nx * d
                qy = oy + sign * ny * d
                if not pav.contains(Point(qx, qy)):
                    if d < best:
                        best = d
                    break
                d += RAY_STEP_M
        return best

    trim_a = 0.0
    while trim_a < max_trim:
        if _perp_hw(trim_a) <= thresh:
            break
        trim_a += step

    trim_b = total_len
    min_b = total_len - max_trim
    while trim_b > min_b:
        if _perp_hw(trim_b) <= thresh:
            break
        trim_b -= step

    if trim_b - trim_a < 10.0:
        return None
    from shapely.ops import substring
    return substring(axis, trim_a, trim_b)


def _probe_axis_width(axis: LineString, pav: Polygon,
                     n_probes: int = 9) -> float:
    """Return 2× the MEDIAN distance-to-boundary along the axis.

    At a large connected pavement (SPJC, where apron + taxi + runway
    are all one blob), ray-casting perpendicular overshoots into the
    apron.  The distance-to-boundary gives the local narrow-corridor
    half-width — which for a taxi centered in its strip is the
    half-strip-width.  Using the median over several probe points
    is robust to both (a) endpoints that sit at wider apron junctions
    and (b) narrow bottlenecks from adjacent building edges.
    """
    if axis.length < 1e-3:
        return 0.0
    boundary = pav.boundary
    dists = []
    for k in range(n_probes):
        t = (k + 1) / (n_probes + 1)
        pt = axis.interpolate(t, normalized=True)
        if not pav.contains(pt):
            continue
        d = pt.distance(boundary)
        if d > 0.1:
            dists.append(d)
    if not dists:
        return 0.0
    dists.sort()
    # Median, then 2× for full width
    return dists[len(dists) // 2] * 2.0


def _extend_rect_corners_perpendicular(
        rect: Polygon, axis: LineString,
        pav: Polygon, max_dist: float = 80.0,
        ) -> Polygon:
    """Extend each rect corner OUTWARD perpendicular to the rect's
    AXIS until the apt.dat pavement boundary (or ``max_dist``).

    Used for runway-end stubs (e.g. SPJC's F) where the half-width
    probe is capped at ``RAY_CAP_M = 40 m`` in
    ``_natural_half_width``, so a stub sitting in a wide runway-end
    ramp ends up under-sized.  The pavement extends past the
    rect's long edges, and the surrounding junction wraps around
    them (wrap-around = polygon along long edge of sloping rect,
    forbidden by the user's invariant).  Per user 2026-04-27: the
    rect should cover the FULL pavement width at each end —
    turning into a trapezoid where each corner sits exactly on the
    pavement boundary independently.

    The returned polygon has the same 4 corners in the same order
    (so X-Plane's altitude_high/low convention is preserved); each
    corner is just shifted outward to its respective pavement edge.
    """
    coords = list(rect.exterior.coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return rect
    a = list(axis.coords)
    if len(a) < 2:
        return rect
    p1, p2 = a[0], a[-1]
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    mag = math.hypot(dx, dy)
    if mag < 1e-6:
        return rect
    ux = dx / mag
    uy = dy / mag
    px = -uy
    py = ux
    STEP = 0.5

    def _ray_extent(ax: float, ay: float,
                    dir_x: float, dir_y: float,
                    base_dist: float) -> float:
        """From axis point (ax, ay) heading (dir_x, dir_y), find
        the FURTHEST distance d such that (ax+d*dir,ay+d*dir) is
        inside ``pav``.  Start at base_dist and walk OUTWARD only
        (we never shrink past base_dist; the rect keeps at least
        its nominal half-width)."""
        # If even base_dist is OUTSIDE pav, the original rect
        # corner is already outside; keep base.
        if not pav.contains(
                Point(ax + dir_x * base_dist,
                      ay + dir_y * base_dist)):
            return base_dist
        d = base_dist
        while d < max_dist:
            d_test = d + STEP
            if not pav.contains(
                    Point(ax + dir_x * d_test,
                          ay + dir_y * d_test)):
                return d
            d = d_test
        return d

    new_corners: List[Tuple[float, float]] = []
    for cx, cy in coords:
        # For each corner: project onto axis to determine which
        # endpoint (p1 or p2) it belongs to and which perp side.
        vx = cx - p1[0]
        vy = cy - p1[1]
        proj = vx * ux + vy * uy
        # axis endpoint nearer this corner
        ax_pt = p1 if proj < mag * 0.5 else p2
        # perpendicular signed distance from axis
        rel_x = cx - ax_pt[0]
        rel_y = cy - ax_pt[1]
        perp_signed = rel_x * px + rel_y * py
        if perp_signed >= 0:
            dir_x, dir_y = px, py
        else:
            dir_x, dir_y = -px, -py
        base = abs(perp_signed)
        if base < 1.0:
            new_corners.append((cx, cy))
            continue
        d = _ray_extent(ax_pt[0], ax_pt[1], dir_x, dir_y, base)
        new_corners.append(
            (ax_pt[0] + dir_x * d, ax_pt[1] + dir_y * d))

    try:
        new_rect = Polygon(new_corners)
        if new_rect.is_valid and not new_rect.is_empty:
            return new_rect
    except Exception:
        pass
    return rect


def _rect_from_axis_extended(axis: LineString, width: float,
                            pav: Polygon,
                            apt_vertices: Optional[
                                List[Tuple[float, float]]] = None,
                            ) -> Optional[Polygon]:
    """Build a rect around the axis at its first-to-last direction.

    The 4 corners are placed at axis endpoints ± perpendicular half-
    width, then each corner is snapped FIRST to the nearest apt.dat
    pavement vertex within ``VERTEX_SNAP_RADIUS_M``, ELSE to the
    nearest pavement edge point within ``EDGE_SNAP_RADIUS_M``.
    This matches the snapped target convention where every non-
    runway vertex sits on an apt.dat pavement vertex.

    Asymmetric-snap trim: when the two snapped end-widths differ
    by more than ``ASYM_WIDTH_TOL_M``, the rect has extended into
    a widening pavement area (apron or junction) on the wider end.
    Per user (2026-04-21): "if we're getting an asymmetric rect,
    most likely it's too long and needs to be shortened a bit so
    it's not pulled into a junction."  Retry once with the axis
    trimmed by ``ASYM_TRIM_FRAC`` of its length on the wider end.
    """
    from shapely.ops import substring

    # Symmetry tolerances: a proper rectangle has equal long sides
    # and equal short sides.  Use a RATIO check so small rects
    # aren't over-trimmed — a 5 m width delta on a 30 m-wide stub
    # (17 %) reads as near-symmetric, while a 5 m delta on a
    # 22 m-wide cross-connector (23 %) reads as a trapezoid.  User
    # (2026-04-21): "I don't see why it should trim the second
    # diagonal at the south end, which is already quite short and
    # appears symmetrical" — that stub had width Δ = 5 m / max
    # 29 m = 17 %, below threshold.
    ASYM_WIDTH_RATIO_TOL = 0.20   # width Δ / max_width > 20 % → trim
    ASYM_LENGTH_RATIO_TOL = 0.10  # length Δ / max_length > 10 % → trim
    # Per-iteration axis shrink per user (2026-04-21): trim BOTH
    # ends by 2.5 % of length each (5 % total) so the rect STAYS
    # CENTERED on its axis as it shrinks — trimming only the
    # wider end shifts the rect toward the narrower side and
    # leaves the problematic (wider) end snapped to the same
    # widening zone on the next iteration.
    ASYM_TRIM_EACH_END = 0.025     # 2.5 % off each end per iteration
    MAX_ASYM_RETRIES = 15          # 15 * 5 % = up to 75 % shrink

    cur_axis = axis
    for attempt in range(MAX_ASYM_RETRIES + 1):
        coords = list(cur_axis.coords)
        if len(coords) < 2:
            return None
        p1 = coords[0]
        p2 = coords[-1]
        dx, dy = p2[0] - p1[0], p2[1] - p1[1]
        mag = math.hypot(dx, dy)
        if mag < 1e-6:
            return None
        ux, uy = dx / mag, dy / mag
        px, py = -uy, ux
        half = width / 2.0
        corners = [
            (p1[0] + px * half, p1[1] + py * half),   # 0: end1 side1
            (p2[0] + px * half, p2[1] + py * half),   # 1: end2 side1
            (p2[0] - px * half, p2[1] - py * half),   # 2: end2 side2
            (p1[0] - px * half, p1[1] - py * half),   # 3: end1 side2
        ]
        snapped = _snap_corners_to_pavement(
            corners, pav, apt_vertices)
        # Reject degenerate rects where snap collapsed two corners
        # onto the same apt.dat vertex.
        degenerate = False
        for i in range(4):
            for j in range(i + 1, 4):
                if math.hypot(snapped[i][0] - snapped[j][0],
                              snapped[i][1] - snapped[j][1]) < 1.0:
                    degenerate = True
                    break
            if degenerate:
                break
        if degenerate:
            return None

        # Symmetry check: equal widths (end1 vs end2) AND equal
        # lengths (side1 vs side2).
        w_end1 = math.hypot(snapped[0][0] - snapped[3][0],
                            snapped[0][1] - snapped[3][1])
        w_end2 = math.hypot(snapped[1][0] - snapped[2][0],
                            snapped[1][1] - snapped[2][1])
        l_side1 = math.hypot(snapped[1][0] - snapped[0][0],
                             snapped[1][1] - snapped[0][1])
        l_side2 = math.hypot(snapped[2][0] - snapped[3][0],
                             snapped[2][1] - snapped[3][1])
        width_asym = abs(w_end1 - w_end2)
        length_asym = abs(l_side1 - l_side2)
        max_w = max(w_end1, w_end2)
        max_l = max(l_side1, l_side2)
        width_ratio = (width_asym / max_w) if max_w > 1e-6 else 0.0
        length_ratio = (length_asym / max_l) if max_l > 1e-6 else 0.0
        symmetric = (width_ratio <= ASYM_WIDTH_RATIO_TOL
                     and length_ratio <= ASYM_LENGTH_RATIO_TOL)
        if symmetric or attempt == MAX_ASYM_RETRIES:
            return Polygon(snapped)

        # Asymmetric — trim BOTH ends and retry, keeping the rect
        # centered on its axis (per user 2026-04-21 feedback).
        trim_each = ASYM_TRIM_EACH_END * cur_axis.length
        new_start = trim_each
        new_end = cur_axis.length - trim_each
        if new_end - new_start < MIN_SEGMENT_LEN_M:
            # Axis would become too short; accept current asymmetric
            # rect rather than discarding.
            return Polygon(snapped)
        try:
            cur_axis = substring(cur_axis, new_start, new_end)
        except Exception:
            return Polygon(snapped)
    return None


VERTEX_SNAP_RADIUS_M = 8.0  # first-choice: snap to real apt.dat vertex
EDGE_SNAP_RADIUS_M = 15.0   # fallback: nearest point on pav boundary


def _snap_corners_to_pavement(
    corners: List[Tuple[float, float]],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]] = None,
) -> List[Tuple[float, float]]:
    """Two-stage corner snap per the user's rule (2026-04-18):

    1. First try nearest apt.dat pavement VERTEX within
       ``VERTEX_SNAP_RADIUS_M``.  Apt.dat vertices are the
       authoritative coordinate set the target also snaps to, so
       snapping rect corners to them produces exact shared-vertex
       alignment.
    2. If no apt.dat vertex within range, fall back to the nearest
       POINT on the pav boundary within ``EDGE_SNAP_RADIUS_M``.
    3. If neither within range, leave the corner unsnapped.

    GUARD: a vertex-snap that would collapse two corners onto the
    SAME apt.dat vertex (producing a degenerate rect) is rejected
    — that corner falls through to edge snap instead.  Rects with
    two coincident corners violate rule 7 and are rejected
    downstream anyway, so we'd rather keep the 4 distinct corners.
    """
    boundary = pav.boundary
    # Pre-filter apt.dat vertices to ONLY those that actually sit
    # on the pav_union boundary.  When two apt.dat pavement
    # polygons overlap (common at aprons / terminal pads), a
    # corner vertex of one polygon ends up in the interior of the
    # union.  Snapping a rect corner to such an interior vertex
    # violates rule 7 ("corners ALWAYS on pavement boundary") and
    # produces a rect that floats inside the pavement, leaving a
    # sliver that a junction wraps around.  Discovered at SPJC V1
    # (2026-04-23): c0 snapped to an apt.dat vertex 5.88 m inside
    # the union, leaving the junction to wrap around V1's short
    # side.
    #
    # Per user 2026-04-27: tolerance is STRICT (1 cm).  Snap means
    # EXACTLY on the boundary, not "close to it" — a 0.36 m offset
    # at SPJC G's NW corner caused the surrounding junction to
    # add two extra nodes wrapping around the offset.  Vertices
    # that aren't on the boundary fall through to Stage 2 (edge
    # snap), which projects directly onto the boundary line.
    BOUNDARY_TOL_M = 0.01
    boundary_verts: Optional[List[Tuple[float, float]]] = None
    if apt_vertices:
        boundary_verts = [
            v for v in apt_vertices
            if Point(v[0], v[1]).distance(boundary) <= BOUNDARY_TOL_M
        ]
    # Stage 1: pick nearest apt.dat boundary vertex per corner.
    candidates: List[Optional[Tuple[float, float]]] = []
    for (cx, cy) in corners:
        best_v = None
        best_d = VERTEX_SNAP_RADIUS_M
        if boundary_verts:
            for (vx, vy) in boundary_verts:
                d = math.hypot(cx - vx, cy - vy)
                if d < best_d:
                    best_v = (vx, vy)
                    best_d = d
        candidates.append(best_v)

    # Collision handling: when two vertex-snap candidates collide
    # (would coincide), keep the nearer one on the vertex.  The
    # other corner uses the ORIGINAL pre-snap coordinate (NOT
    # edge-snap, which could pull back to the same region).  This
    # preserves a valid 4-corner rect while still placing the
    # kept corner exactly on an apt.dat vertex.
    # Use proximity (not identity) — two candidates within 1 m
    # would also collapse the rect corner.  Collapse collision:
    # drop the farther-from-original candidate, use pre-snap.
    use_original: List[bool] = [False] * len(candidates)
    COLLISION_TOL = 1.0
    for i in range(len(candidates)):
        if candidates[i] is None:
            continue
        for j in range(i + 1, len(candidates)):
            if candidates[j] is None:
                continue
            d_cand = math.hypot(candidates[i][0] - candidates[j][0],
                                candidates[i][1] - candidates[j][1])
            if d_cand <= COLLISION_TOL:
                di = math.hypot(corners[i][0] - candidates[i][0],
                                corners[i][1] - candidates[i][1])
                dj = math.hypot(corners[j][0] - candidates[j][0],
                                corners[j][1] - candidates[j][1])
                if di <= dj:
                    candidates[j] = None
                    use_original[j] = True
                else:
                    candidates[i] = None
                    use_original[i] = True
                    # candidates[i] is now None; subsequent j
                    # iterations would dereference it.  Break and
                    # let the outer i loop advance.
                    break

    snapped: List[Tuple[float, float]] = []
    for i, (cx, cy) in enumerate(corners):
        if candidates[i] is not None:
            snapped.append(candidates[i])
            continue
        if use_original[i]:
            # Collision fallback — keep pre-snap coord to avoid
            # coincident corners.
            snapped.append((cx, cy))
            continue
        # Stage 2: nearest pav edge point.
        p = Point(cx, cy)
        near, _ = nearest_points(boundary, p)
        if p.distance(near) <= EDGE_SNAP_RADIUS_M:
            snapped.append((near.x, near.y))
        else:
            snapped.append((cx, cy))

    # Final coincidence sweep: vertex-snap AND edge-snap can BOTH
    # pull two corners onto the same area (e.g. edge-snap falls
    # through to the same boundary vertex a vertex-snap picked
    # from the other corner).  If any two final coords are within
    # 1 m, revert the farther-from-original to its pre-snap coord.
    FINAL_COLLISION_TOL = 1.0
    for i in range(len(snapped)):
        for j in range(i + 1, len(snapped)):
            d = math.hypot(snapped[i][0] - snapped[j][0],
                           snapped[i][1] - snapped[j][1])
            if d <= FINAL_COLLISION_TOL:
                di = math.hypot(corners[i][0] - snapped[i][0],
                                corners[i][1] - snapped[i][1])
                dj = math.hypot(corners[j][0] - snapped[j][0],
                                corners[j][1] - snapped[j][1])
                if di <= dj:
                    snapped[j] = corners[j]
                else:
                    snapped[i] = corners[i]

    # Symmetry check: rect corners come in as [end1_side1, end2_side1,
    # end2_side2, end1_side2] per _rect_from_axis_extended.  The two
    # widths (end1_side1 → end1_side2 and end2_side1 → end2_side2)
    # should be equal for a proper rectangle.  When snap pulls
    # corners to apt.dat vertices at asymmetric offsets (one rect
    # end near a wider pavement apron than the other), the widths
    # diverge and the shape reads as a trapezoid.  Per user
    # (2026-04-21): "the perpendicular stub is asymmetrical and
    # looks like it's coming into the space of the primary
    # parallel."  If the snapped widths differ by more than
    # ``ASYM_WIDTH_TOL_M``, revert ONE corner on each side to
    # its pre-snap coord so the rect stays symmetric.  We keep
    # the NARROWER side's snap (matching the tighter pavement) and
    # revert the wider side's corners to the pre-snap perpendicular
    # offset.
    return snapped


def _cap_rect_length_to_width(
    taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    rwy_centerlines: List[LineString],
    pav: Polygon,
    apt_vertices: Optional[List[Tuple[float, float]]],
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Cap each rect's length-to-width ratio per the user's
    2026-04-27 spec: rects should be roughly square (length ≈
    width) so the long edges sit on the pavement-narrowing
    boundary, corners snap there, and surrounding junctions
    connect only at the short edges (never wrap around long
    edges).

    Cap depends on the bearing-to-nearest-runway:

      * Parallel  (db < 20°)  — NO CAP (long parallel taxis are
        legitimate, often running 500 m+ along the runway).
      * Diagonal  (20° ≤ db < 45°) — length ≤ 1.0 × width
        (truly square; matches the tighter 30 % margin used for
        diagonal stubs in ``_rect_margin_frac_for``).
      * Perpendicular (db ≥ 45°) — length ≤ 1.3 × width (small
        excess so the rect can extend slightly past the apron's
        narrow corridor without forcing surrounding junctions to
        wrap).

    Shrinks the axis symmetrically (same amount from both ends) so
    the rect's centre stays put, then re-runs
    ``_rect_from_axis_extended`` so corners re-snap to apt.dat
    pavement boundary on the new axis.
    """
    PERP_CAP_RATIO = 1.3
    DIAG_CAP_RATIO = 1.0
    from shapely.ops import substring
    out: List[Tuple[Polygon, LineString, str, str]] = []
    for rect, axis, role, ref in taxi_rects:
        try:
            coords = list(rect.exterior.coords)
        except Exception:
            out.append((rect, axis, role, ref))
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            out.append((rect, axis, role, ref))
            continue
        edge_lens = [
            math.hypot(coords[(i + 1) % 4][0] - coords[i][0],
                       coords[(i + 1) % 4][1] - coords[i][1])
            for i in range(4)]
        length = max(edge_lens)
        width = min(edge_lens)
        if width < 1.0 or length < 1.0:
            out.append((rect, axis, role, ref))
            continue
        db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
        if db is None:
            out.append((rect, axis, role, ref))
            continue
        if db < 20.0:
            # Parallel — no cap.
            out.append((rect, axis, role, ref))
            continue
        cap_ratio = (DIAG_CAP_RATIO if db < 45.0
                     else PERP_CAP_RATIO)
        max_len = cap_ratio * width
        if length <= max_len + 0.5:
            out.append((rect, axis, role, ref))
            continue
        axis_len = axis.length
        new_axis_len = max_len
        margin = (axis_len - new_axis_len) / 2.0
        if margin <= 0:
            out.append((rect, axis, role, ref))
            continue
        try:
            new_axis = substring(
                axis, margin, axis_len - margin)
        except Exception:
            out.append((rect, axis, role, ref))
            continue
        if (new_axis.is_empty
                or new_axis.geom_type != "LineString"
                or new_axis.length < 5.0):
            out.append((rect, axis, role, ref))
            continue
        new_rect = _rect_from_axis_extended(
            new_axis, width, pav, apt_vertices=apt_vertices)
        if (new_rect is None or new_rect.is_empty
                or new_rect.geom_type != "Polygon"
                or not new_rect.is_valid):
            out.append((rect, axis, role, ref))
            continue
        out.append((new_rect, new_axis, role, ref))
    return out


def _classify_role(axis: LineString, width: float,
                   rwy_centerlines: List[LineString],
                   rwy_union: Optional[Polygon],
                   ref: str = "",
                   ref_overall_bearings: Optional[Dict[str, float]]
                   = None) -> str:
    """Classify a taxi rect by axis geometry alone.

    The role is determined entirely by:
      * bearing to nearest runway (parallel within 20°,
        perpendicular beyond 45°),
      * straight-line distance from axis midpoint to nearest runway
        centerline (close → primary, far → secondary or cross),
      * axis length (must clear minimum length per role).

    Ref letters are NOT used as a classifier — they vary wildly
    between airports (SPJC's parallel taxis are A/F/L/V; CYXY's
    is E; KBNA uses different letters again; many CYXY taxis have
    no ref at all).  The only ref-pattern rules retained are:

      * Sub-ref (any digit in the label, e.g. V1, L3, A1) →
        always STUB.  Sub-ref tagging is universal: a digit
        suffix means a short connector spur regardless of airport.
      * Diagonal parent ref → always STUB.  When the parent OSM
        way's chord bearing is itself diagonal (db ≥ 20°), every
        segment of that ref is part of a diagonal stub even if
        a curving end-segment happens to align near-parallel
        locally.  At SPJC, B/C/E/G enter the runway at shallow
        angles; without this check, the post-curve sub-segment
        (B's 85 m piece, db_local = 18°) gets misclassified as
        PRIMARY_PARALLEL even though there's no actual B parallel
        taxiway.

    Roles:
      * PRIMARY_PARALLEL  — db < 20°, length ≥ 50 m, < 400 m from runway
      * SECONDARY_PARALLEL — db < 20°, length ≥ 50 m, ≥ 400 m from runway
      * CROSS_CONNECTOR   — db > 45°, length ≥ 80 m, > 250 m from runway
      * STUB              — everything else (short, runway-adjacent perp, etc.)
    """
    if ref and any(c.isdigit() for c in ref):
        # Sub-refs are always stubs.
        return ROLE_STUB

    db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
    if db is None:
        return ROLE_STUB

    # Diagonal-parent check: if this rect's REF has an overall-
    # DIAGONAL parent OSM way (parent db_overall ∈ [20°, 45°)),
    # force STUB regardless of the local segment bearing.  See
    # the ``Diagonal parent ref`` rule above.  Excludes
    # perpendicular parents (db ≥ 45°, e.g. cross-connectors
    # like Q/R at SPJC) which legitimately classify as
    # CROSS_CONNECTOR via the local-axis check below.
    if (ref and ref_overall_bearings
            and ref in ref_overall_bearings
            and rwy_centerlines):
        try:
            _rwy = min(rwy_centerlines,
                       key=lambda r: axis.distance(r))
            _rc = list(_rwy.coords)
            _rdx = _rc[-1][0] - _rc[0][0]
            _rdy = _rc[-1][1] - _rc[0][1]
            if math.hypot(_rdx, _rdy) > 1e-6:
                _rwy_bearing = (math.degrees(
                    math.atan2(_rdx, _rdy)) % 180.0)
                _ref_db = abs(ref_overall_bearings[ref]
                              - _rwy_bearing)
                _ref_db = min(_ref_db, 180.0 - _ref_db)
                if 20.0 <= _ref_db < 45.0:
                    return ROLE_STUB
        except Exception:
            pass
    try:
        mid = axis.interpolate(0.5, normalized=True)
        dist_rwy = min(mid.distance(r) for r in rwy_centerlines)
    except Exception:
        dist_rwy = 1e6
    length = axis.length
    if db < 20.0:
        if length >= 50.0:
            if dist_rwy < 400.0:
                return ROLE_PRIMARY_PARALLEL
            return ROLE_SECONDARY_PARALLEL
        return ROLE_STUB
    if db > 45.0 and length >= 80.0:
        if dist_rwy > 250.0:
            return ROLE_CROSS_CONNECTOR
        return ROLE_STUB
    return ROLE_STUB


def _axis_to_nearest_rwy_db(axis: LineString,
                            rwy_centerlines: List[LineString]
                            ) -> Optional[float]:
    """Return the bearing difference from ``axis`` to the nearest
    runway centerline, modulo 180°."""
    if not rwy_centerlines:
        return None
    rwy = min(rwy_centerlines, key=lambda r: axis.distance(r))

    def _bearing(ls):
        c = list(ls.coords)
        return math.degrees(math.atan2(c[-1][0] - c[0][0],
                                       c[-1][1] - c[0][1])) % 180.0
    db = abs(_bearing(axis) - _bearing(rwy))
    return min(db, 180.0 - db)


def _refine_roles(emitted, rwy_centerlines):
    """Post-classify: demote the stub-A / stub-F segment (the short
    runway-connector within a parallel ref's polyline) from
    primary_parallel to stub.

    Rule: for each parallel ref, find SEGMENTS that are
    significantly perpendicular (>= 40° off runway) AND short (< 150 m).
    Demote to stub.  Multiple per ref allowed.
    """
    if not rwy_centerlines or not emitted:
        return
    for i, (rect, axis, role, ref) in enumerate(emitted):
        if role != ROLE_PRIMARY_PARALLEL:
            continue
        db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
        if db is None:
            continue
        if db >= 40.0 and axis.length < 150.0:
            emitted[i] = (rect, axis, ROLE_STUB, ref)
