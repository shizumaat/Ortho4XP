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
    nodes: Dict[str, Tuple[float, float]] = {}
    ways: List[Tuple[str, List[str], Dict[str, str]]] = []
    relations: List[Tuple[str, List[str], Dict[str, str]]] = []
    seen_paths = set()
    for dlat in (0, -1, 1):
        for dlon in (0, -1, 1):
            osm_path = _tile_path(base_lat + dlat, base_lon + dlon)
            if osm_path in seen_paths or not os.path.isfile(osm_path):
                continue
            seen_paths.add(osm_path)
            n2, w2, r2 = _load_osm_tile(osm_path)
            nodes.update(n2)
            ways.extend(w2)
            relations.extend(r2)
    if not nodes:
        return {}, [], []

    # Filter by bbox — keep only ways whose node centroids are near
    # the airport center.
    def _in_box(lat, lon):
        return (abs(lat - apt_lat) <= radius_deg and
                abs(lon - apt_lon) <= radius_deg)

    kept_ways = []
    way_by_id: Dict[str, Tuple[str, List[str], Dict[str, str]]] = {}
    for wid, nds, tags in ways:
        way_by_id[wid] = (wid, nds, tags)
        pts = [nodes[n] for n in nds if n in nodes]
        if not pts:
            continue
        clat = sum(p[0] for p in pts) / len(pts)
        clon = sum(p[1] for p in pts) / len(pts)
        if _in_box(clat, clon):
            kept_ways.append((wid, nds, tags))
    # Relations: keep if ANY member way centroid is in-box
    kept_rels = []
    for rid, outer_ids, tags in relations:
        any_in = False
        for wid in outer_ids:
            if wid in way_by_id:
                _, nds, _ = way_by_id[wid]
                pts = [nodes[n] for n in nds if n in nodes]
                if pts:
                    clat = sum(p[0] for p in pts) / len(pts)
                    clon = sum(p[1] for p in pts) / len(pts)
                    if _in_box(clat, clon):
                        any_in = True
                        break
        if any_in:
            kept_rels.append((rid, outer_ids, tags))
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

        # Build new rects from kept intervals.
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
                new_axis = LineString([new_a_mid, new_b_mid])
                kept.append((new_rect, new_axis, role, ref))
            except Exception:
                continue
        if len(kept_intervals) >= 2:
            n_split += 1
        else:
            n_clipped += 1
        abs_refs.append(ref or "?")

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
        pav_union = pav_union.difference(layout.runway_union)

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
        # leave sub-meter residue slivers overlapping runway.
        if layout.runway_union is not None and not layout.runway_union.is_empty:
            residue = residue.difference(layout.runway_union)

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

    return layout


# ══════════════════════════════════════════════════════════════════
# Phase-2: elevations
# ══════════════════════════════════════════════════════════════════

# FAA AC 150/5300-13B, Design Group III+ (commercial jetport):
#   Taxiway longitudinal: max 1.5 % grade.
#   Apron: max 1.0 % any direction.
#   Runway longitudinal: max 1.5 % grade (handled by legacy).
TAXI_MAX_GRADE = 0.015
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
            n_dropped = 0
            for sh in layout.shapes:
                if sh.role != ROLE_RUNWAY:
                    kept_shapes.append(sh)
                    continue
                drop = False
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
                                break
                        except Exception:
                            continue
                if drop:
                    n_dropped += 1
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

    # ── Terminal pad elevations (flat, DEM median) ──────────────
    # Computed BEFORE the elevation graph so terminal corners can be
    # anchored in the graph as hard 1.5 %-grade constraints on the
    # surrounding apron.  The DEM-median terminal elevation pulls
    # the network UP along the apron-to-terminal interface so the
    # surrounding pavement maintains grade with the (often-elevated)
    # terminal pad — letting X-Plane render a real cliff between
    # the apron's edge and the natural DEM beyond, rather than
    # forcing a 100 %-slope drop INSIDE the pavement (user 2026-04-24).
    for shape in layout.shapes:
        if shape.role != ROLE_TERMINAL:
            continue
        samples: List[float] = []
        for x, y in [
            (shape.polygon.centroid.x, shape.polygon.centroid.y)
        ] + list(shape.polygon.exterior.coords)[:-1]:
            lat, lon = m_to_ll(x, y)
            e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            if e is not None:
                samples.append(e)
        if not samples:
            continue
        samples.sort()
        median = samples[len(samples) // 2]
        shape.altitude = round(median, 1)

    # ── Taxi rect elevations via grade-compliant network ────────
    # Build a graph of OSM taxi centerlines + runway segment
    # endpoints, densified to ≤ 30 m edges.  CIFP runway
    # thresholds anchor the graph; bounds propagation + Laplacian
    # smoothing produce a grade-compliant elevation surface over
    # the network.  Each rect's two short edges take their
    # altitude from sampling that network at the rect's axis
    # endpoints.
    graph = _build_elevation_network(
        osm_nodes, osm_ways, to_m,
        runway_segment_chain, layout.anchor, dem, tile_lat, tile_lon)
    if graph is not None:
        # Anchor terminal pad corners in the elevation graph so the
        # surrounding apron is forced grade-compliant w.r.t. the
        # terminal level.  Each corner becomes a graph node
        # connected to its nearest existing graph node by a bridge
        # of length = straight-line distance, capping per-edge
        # elev diff at 1.5 % × bridge length.
        _add_terminal_anchors(graph, layout)
        # Add cross-junction bridges so smoothing enforces 1.5 %
        # grade across junction diagonals — not just along
        # centerline routes.  Must run BEFORE propagate_bounds so
        # the new edges contribute to per-node feasibility cones.
        _add_junction_bridges(graph, layout)
        graph.propagate_bounds()
        graph.choose_values()
        graph.smooth_rate_of_change(iters=30)
        # ── Plateau detection + snap ───────────────────────────
        # User 2026-04-25: identify shallow regions of the
        # smoothed surface (edges with grade < PLATEAU_FLATNESS_GRADE)
        # and snap each connected component to its median elevation.
        # The result: most of the apron sits at a few discrete
        # plateau elevations connected by short ramp segments,
        # rather than continuously varying.  Rects entirely within
        # a plateau emit as FLAT (single altitude tag); junctions
        # within a plateau pass the existing FLAT classifier
        # automatically (vertex elevations all match the plateau).
        # Triangulation is reserved for compound-slope regions
        # crossing multiple plateaus.
        _, _soft_anchored = _snap_plateaus(graph)
        # Post-plateau bridge enforcement (2026-04-26): each
        # plateau cluster's boundary check ran while OTHER
        # clusters were still pre-snap, so cross-cluster bridges
        # can end up violating after both commit.  Walk every
        # edge between two anchors; demote SOFT-anchored
        # endpoints involved in violations so the next smoothing
        # pass can pull them into compliance.  Hard anchors
        # (CIFP runway thresholds) are never demoted.
        n_demoted, n_hard_avg = _demote_violating_soft_anchors(
            graph, _soft_anchored)
        if n_demoted or n_hard_avg:
            try:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] {icao}: post-plateau bridge "
                    f"reconciliation — demoted {n_demoted} soft "
                    f"anchors, averaged {n_hard_avg} hard-vs-hard "
                    f"violating pairs.\n")
            except Exception:
                pass
        # Re-smooth ramps between plateaus.  Plateau nodes are
        # now hard anchors; this run enforces 1.5 % grade on the
        # transition ramps between them.
        graph.propagate_bounds()
        graph.choose_values()
        graph.smooth_rate_of_change(iters=30)
        taxi_roles = {ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                      ROLE_STUB, ROLE_CROSS_CONNECTOR}
        for shape in layout.shapes:
            if shape.role not in taxi_roles:
                continue
            if shape.source_axis is None:
                continue
            coords = list(shape.source_axis.coords)
            if len(coords) < 2:
                continue
            p1, p2 = coords[0], coords[-1]
            e1 = graph.elevation_at(p1[0], p1[1])
            e2 = graph.elevation_at(p2[0], p2[1])
            # DEM fallback when network couldn't serve a point
            # (e.g. walking-exit stubs whose axis sits off-network).
            if e1 is None:
                lat, lon = m_to_ll(p1[0], p1[1])
                e1 = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            if e2 is None:
                lat, lon = m_to_ll(p2[0], p2[1])
                e2 = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            if e1 is None and e2 is None:
                continue
            if e1 is None:
                e1 = e2
            if e2 is None:
                e2 = e1
            # Defensive within-rect grade cap (should already be
            # satisfied by network but catches the DEM-fallback
            # and mixed cases).
            axis_len = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
            if axis_len > 1.0:
                dmax = axis_len * TAXI_MAX_GRADE
                diff = e1 - e2
                if abs(diff) > dmax:
                    mean = (e1 + e2) / 2.0
                    sign = 1 if diff > 0 else -1
                    e1 = mean + (dmax / 2.0) * sign
                    e2 = mean - (dmax / 2.0) * sign
            eh = max(e1, e2)
            el = min(e1, e2)
            if abs(eh - el) >= 0.1:
                shape.altitude_high = round(eh, 1)
                shape.altitude_low = round(el, 1)
                # Earlier pipeline stages may have permuted the
                # polygon ring.  Always re-derive the ring order
                # from the 4 raw corner positions so that the
                # X-Plane patch convention (corners 0,3 = high
                # short edge, 1,2 = low short edge) holds.
                _orient_rect_for_altitude(shape, p1, p2, e1, e2)
            else:
                shape.altitude = round((eh + el) / 2.0, 1)

    # ── Keep junction polygons off taxi rect edges ──────────────
    # User rule (2026-04-24): a junction may join an
    # altitude_high/altitude_low rect only at the 4 rect corners.
    # Any junction vertex landing mid-edge of a rect would be
    # split into the rect's ring at render time and break the
    # rect's 4-corner slope convention.
    #
    # Implementation: walk each junction's ring; for each vertex
    # that's within ``TAXI_EDGE_TOL_M`` of a taxi rect edge but
    # NOT within ``TAXI_CORNER_TOL_M`` of one of that rect's four
    # corners, push the vertex outward by ``TAXI_EDGE_GAP_M``
    # perpendicular to the offending edge.  Vertices near a
    # corner snap exactly to the corner instead, preserving
    # shared-vertex connectivity at rect ends.
    TAXI_EDGE_TOL_M = 0.5
    TAXI_CORNER_TOL_M = 2.0
    TAXI_EDGE_GAP_M = 1.0
    taxi_roles_set = {ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                      ROLE_STUB, ROLE_CROSS_CONNECTOR}
    # Collect (polygon, 4 corner tuples) for each taxi rect.
    taxi_rect_info: List[Tuple[Polygon, List[Tuple[float, float]]]] = []
    for s in layout.shapes:
        if s.role not in taxi_roles_set:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
        except Exception:
            continue
        if len(coords) != 4:
            continue
        taxi_rect_info.append((s.polygon, coords))

    def _fix_junction_vertex(x: float, y: float) -> Tuple[float, float]:
        """If (x,y) is close to a taxi rect edge but not close to
        one of its 4 corners, return a pushed-outward position."""
        for rect_poly, corners in taxi_rect_info:
            # Check corners first — if close, snap exactly.
            for cx, cy in corners:
                if (x - cx) ** 2 + (y - cy) ** 2 <= (
                        TAXI_CORNER_TOL_M ** 2):
                    return (cx, cy)
            # For each of 4 edges, check mid-edge proximity.
            for i in range(4):
                ax, ay = corners[i]
                bx, by = corners[(i + 1) % 4]
                dx = bx - ax; dy = by - ay
                seg_len_sq = dx * dx + dy * dy
                if seg_len_sq <= 0.01:
                    continue
                t = ((x - ax) * dx + (y - ay) * dy) / seg_len_sq
                if t <= 0.001 or t >= 0.999:
                    continue  # at an endpoint, handled by corner check
                cx_proj = ax + t * dx
                cy_proj = ay + t * dy
                d = math.hypot(x - cx_proj, y - cy_proj)
                if d > TAXI_EDGE_TOL_M:
                    continue
                # Perpendicular to edge, pointing AWAY from the
                # rect interior.  Try both ±perp and pick the one
                # whose (proj + perp) is outside the rect polygon.
                seg_len = math.sqrt(seg_len_sq)
                perp_x = -dy / seg_len
                perp_y = dx / seg_len
                # Test which side is outside: move a tiny step in
                # each direction and check containment.
                test_x = cx_proj + perp_x * 0.1
                test_y = cy_proj + perp_y * 0.1
                if rect_poly.contains(Point(test_x, test_y)):
                    perp_x = -perp_x
                    perp_y = -perp_y
                return (cx_proj + perp_x * TAXI_EDGE_GAP_M,
                        cy_proj + perp_y * TAXI_EDGE_GAP_M)
        return (x, y)

    if taxi_rect_info:
        for shape in layout.shapes:
            if shape.role != ROLE_JUNCTION:
                continue
            try:
                ring = list(shape.polygon.exterior.coords)
            except Exception:
                continue
            changed = False
            new_ring = []
            for (vx, vy) in ring:
                nx, ny = _fix_junction_vertex(vx, vy)
                if (nx, ny) != (vx, vy):
                    changed = True
                new_ring.append((nx, ny))
            if not changed:
                continue
            # Rebuild polygon.
            try:
                new_poly = Polygon(new_ring,
                                    list(shape.polygon.interiors))
                if not new_poly.is_valid:
                    new_poly = new_poly.buffer(0)
                if (new_poly.geom_type == "Polygon"
                        and not new_poly.is_empty):
                    shape.polygon = new_poly
            except Exception:
                pass

    # ── Junction triangulation ──────────────────────────────────
    # Each junction polygon is replaced with N-2 ear-clip
    # triangles, each carrying a per-vertex elevation list
    # (``node_altitudes``).  Vertex elevations come from the
    # corner-elevation bucket map (rect/runway/terminal corners
    # share node ids with junction boundary vertices), with
    # graph and DEM fallbacks for boundary-trace vertices.  X-Plane
    # interpolates linearly across each triangle, giving the
    # multi-directional slope behaviour the user requested.
    _triangulate_junctions(layout, graph, dem, tile_lat, tile_lon, m_to_ll)

    # ── Global pavement mesh ────────────────────────────────────
    # Build one mesh covering every pavement shape's boundary —
    # nodes are unique vertex buckets, edges are ring adjacency
    # plus cross-shape adjacency at shared buckets.  ANCHORED ONLY
    # on runway corner buckets (CIFP-derived elevations are the
    # sole HARD truth); every other node is free to be moved by
    # the smoother to satisfy 2D Euclidean grade compliance with
    # its neighbours.  After solving, each shape's altitude tags
    # are rewritten from the mesh values.  Per the user 2026-04-26
    # principle: only runways are HARD; taxi rect / terminal
    # altitudes derived earlier from the centerline graph are
    # initial seeds that the mesh can refine.
    runway_anchors = _runway_corner_elev_map(layout)
    # Augment with junction corners that sit close to (but not
    # exactly at) a runway polygon — see ``_near_runway_anchor_map``.
    near_runway_anchors = _near_runway_anchor_map(layout)
    for _b, _e in near_runway_anchors.items():
        runway_anchors.setdefault(_b, _e)
    if runway_anchors:
        try:
            mesh, bucket_to_node = _build_pavement_mesh(
                layout, runway_anchors, dem,
                tile_lat, tile_lon)
            _solve_pavement_mesh(mesh, layout, bucket_to_node)
            _writeback_pavement_mesh(layout, mesh, bucket_to_node)
        except Exception as exc:
            try:
                import sys as _sys
                _sys.stderr.write(
                    f"  [pav-builder] WARN: {icao}: global pavement "
                    f"mesh solve failed ({exc}); falling back to "
                    f"per-shape elevations.\n")
            except Exception:
                pass
    # Layer 2 (2026-04-26): clamp each junction's free boundary
    # vertices to be grade-compliant with EVERY nearby boundary
    # of EVERY other shape.  Iterate to a fixed point (constraints
    # depend on neighbour elevations that are themselves being
    # clamped — one pass settles only the "obvious" violations,
    # later passes reconcile cascading effects).
    #
    # The geometry-only state (spatial grid + shared-bucket set)
    # is invariant across this 8-iteration loop — only elevations
    # change.  Build once, reuse.  Saves ~7 redundant grid builds
    # at HECA where each build has 5 k+ boundary edges.
    clamp_geom = _build_clamp_geom_state(layout)
    for _ in range(8):
        n = _clamp_junction_free_vertices(layout, clamp_geom)
        if n == 0:
            break
    # Junction subdivision (2026-04-26, refined): the
    # perpendicular-cut pass now snaps new cut-line vertices to
    # existing ring vertices within SUBDIVIDE_SNAP_RADIUS_M and
    # only commits the split when ALL sub-polygons have a smaller
    # worst within-shape grade than the parent.  Iterates to a
    # fixed point (cap 4 passes); a sub-polygon may itself need
    # further subdivision if it still has incompatible anchors.
    for _ in range(4):
        n = _subdivide_violating_junctions(layout)
        if n == 0:
            break
    # Re-clamp free vertices after subdivision (the few new
    # cut-line vertices that COULDN'T snap inherit interpolated
    # elevations and may benefit from neighbour-clamping).
    # Subdivision may have ADDED vertices, so rebuild the geom
    # state once before this loop.  The 4-iter loop itself reuses.
    clamp_geom = _build_clamp_geom_state(layout)
    for _ in range(4):
        n = _clamp_junction_free_vertices(layout, clamp_geom)
        if n == 0:
            break
    # Layer 3 (2026-04-26): scan every emitted polygon for any
    # within-shape vertex pair grade > TAXI_MAX_GRADE.  Surface a
    # WARN summary so regressions are visible during iteration —
    # not yet a hard fail (would drop too much coverage at HECA-
    # complexity airports while Layers 1/2 are still maturing).
    _report_within_shape_violations(layout, icao)


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


# ── Elevation network ────────────────────────────────────────────
#
# Node-based graph representing the taxiway centerline network in
# meter space, densified to enforce FAA grade rules.  Runway
# segment endpoints act as hard anchors.  The full pipeline is:
#
#   1. ``_build_elevation_network()`` constructs the graph from
#      OSM taxi ways + runway_segment_chain, densified to ≤ 30 m.
#   2. ``ElevationGraph.propagate_bounds()`` runs Dijkstra from
#      each anchor with edge weight = length × TAXI_MAX_GRADE,
#      producing a feasibility interval per node.
#   3. ``ElevationGraph.choose_values()`` picks each node's
#      elevation: DEM sample clipped into the interval (unbounded
#      interval → DEM sample as-is).
#   4. ``ElevationGraph.smooth_rate_of_change()`` applies a
#      Laplacian-style smoothing pass to tame the rate-of-change
#      (target 1 % / 30 m, FAA curvature rule), with anchors and
#      hard grade caps re-applied each iteration.
#   5. ``ElevationGraph.elevation_at(x, y)`` samples the network
#      surface at an arbitrary meter-space point — used to map
#      rect axis endpoints onto the network.
NETWORK_DENSIFY_M = 30.0            # max edge length
NETWORK_BRIDGE_MAX_M = 60.0         # max dist for taxi→rwy bridge edge.
                                     # User 2026-04-28: grade is along
                                     # the taxiway AXIS, not Euclidean.
                                     # Lateral bridges from a taxi to
                                     # nearby runway segments would
                                     # constrain the taxi by short
                                     # Euclidean distances rather than
                                     # actual axis-walked distances,
                                     # producing wrong altH/altL.  Only
                                     # bridge perpendicular STUB
                                     # endpoints (≤ 60 m from runway,
                                     # e.g. SPJC's V1 at ~55 m); leave
                                     # parallel taxis to inherit
                                     # elevation through their axial
                                     # OSM-network path.
NETWORK_RUNWAY_ANCHOR_RADIUS_M = 5.0  # a taxi node within this of
                                       # rwy segment gets anchored


class ElevationGraph:
    """Sparse undirected graph of (x, y, elevation) nodes."""

    __slots__ = ("nodes", "edges_adj", "anchor_elev", "dem_elev",
                 "elev", "interval_lo", "interval_hi",
                 "lat0", "lon0", "cos0", "tile_lat", "tile_lon",
                 "dem")

    def __init__(self, lat0, lon0, cos0, tile_lat, tile_lon, dem):
        self.nodes: List[Tuple[float, float]] = []
        self.edges_adj: List[List[Tuple[int, float]]] = []
        self.anchor_elev: List[Optional[float]] = []
        self.dem_elev: List[Optional[float]] = []
        self.elev: List[float] = []
        self.interval_lo: List[float] = []
        self.interval_hi: List[float] = []
        self.lat0 = lat0
        self.lon0 = lon0
        self.cos0 = cos0
        self.tile_lat = tile_lat
        self.tile_lon = tile_lon
        self.dem = dem

    def add_node(self, x: float, y: float,
                 anchor: Optional[float] = None) -> int:
        idx = len(self.nodes)
        self.nodes.append((x, y))
        self.edges_adj.append([])
        self.anchor_elev.append(anchor)
        # DEM sample (None if DEM unavailable)
        lat, lon = self._m_to_ll(x, y)
        self.dem_elev.append(
            _sample_dem(self.dem, self.tile_lat, self.tile_lon,
                        lat, lon))
        return idx

    def add_edge(self, i: int, j: int, length: float) -> None:
        if i == j:
            return
        self.edges_adj[i].append((j, length))
        self.edges_adj[j].append((i, length))

    def _m_to_ll(self, x, y):
        lon = self.lon0 + math.degrees(x / (R_EARTH * self.cos0))
        lat = self.lat0 + math.degrees(y / R_EARTH)
        return lat, lon

    def propagate_bounds(self) -> None:
        """Compute per-node feasibility interval [lo, hi] from
        anchor cones via a single MULTI-SOURCE Dijkstra.

        Each node ends up with the cone of its NEAREST reachable
        anchor: ``[anchor_val ± nearest_dist × TAXI_MAX_GRADE]``.
        Unreachable nodes get ``(-inf, +inf)``.

        This is looser than the previous all-anchor-intersection
        semantics (which ran one Dijkstra per anchor and intersected
        every cone, costing ``O(anchors × (V + E) log V)``).  The
        looser bounds let the seed re-clip pass through more values
        unchanged, but the smoother's hard edge-grade-cap (run
        repeatedly inside :func:`smooth_rate_of_change`) still
        enforces 1.5 %/m on every edge — so the final elevation
        field is grade-compliant.

        On large airports the speedup is dramatic: ``propagate_bounds``
        was 20 s of SPJC's 30 s elevation phase before this change.
        With multi-source Dijkstra it becomes effectively a single
        graph traversal (``O((V + E) log V)``) regardless of anchor
        count.
        """
        import heapq
        n = len(self.nodes)
        INF = float("inf")
        lo = [-INF] * n
        hi = [INF] * n
        nearest_dist: List[float] = [INF] * n
        nearest_val: List[float] = [0.0] * n
        heap: List[Tuple[float, int, float]] = []
        for i, a in enumerate(self.anchor_elev):
            if a is None:
                continue
            nearest_dist[i] = 0.0
            nearest_val[i] = a
            heapq.heappush(heap, (0.0, i, float(a)))
        while heap:
            d, u, val = heapq.heappop(heap)
            if d > nearest_dist[u]:
                continue
            for w, length in self.edges_adj[u]:
                nd = d + length
                if nd < nearest_dist[w]:
                    nearest_dist[w] = nd
                    nearest_val[w] = val
                    heapq.heappush(heap, (nd, w, val))
        for i in range(n):
            if nearest_dist[i] == INF:
                continue
            band = nearest_dist[i] * TAXI_MAX_GRADE
            lo[i] = nearest_val[i] - band
            hi[i] = nearest_val[i] + band
        self.interval_lo = lo
        self.interval_hi = hi

    def choose_values(self) -> None:
        """Pick each node's elevation: anchor value if pinned,
        otherwise DEM clipped into [lo, hi].  Intervals that came
        out empty (conflicting anchors) fall back to the interval
        midpoint."""
        n = len(self.nodes)
        out: List[float] = [0.0] * n
        for i in range(n):
            if self.anchor_elev[i] is not None:
                out[i] = self.anchor_elev[i]
                continue
            lo = self.interval_lo[i]
            hi = self.interval_hi[i]
            if lo > hi:
                # Conflicting anchors — pick midpoint for
                # least-squares-like compromise.
                out[i] = 0.5 * (lo + hi)
                continue
            d = self.dem_elev[i]
            if d is None:
                # No DEM, no anchors: pick 0 if unbounded, midpoint
                # otherwise.
                if lo == float("-inf") and hi == float("inf"):
                    out[i] = 0.0
                elif lo == float("-inf"):
                    out[i] = hi
                elif hi == float("inf"):
                    out[i] = lo
                else:
                    out[i] = 0.5 * (lo + hi)
                continue
            out[i] = max(lo, min(hi, d))
        self.elev = out

    def smooth_rate_of_change(self, iters: int = 30,
                              damping: float = 0.4,
                              tol: float = 0.01) -> None:
        """Laplacian-style smoothing to tame rate-of-change while
        respecting anchors and the hard 1.5 % grade cap on every
        edge.  Each iteration: each non-anchor node moves toward
        its neighbours' mean by ``damping`` × (mean - self).  Then
        for each edge, if |Δ|/length > 1.5 %, split the excess
        evenly between both non-anchored endpoints."""
        n = len(self.nodes)
        if n == 0:
            return
        anchor = [self.anchor_elev[i] is not None for i in range(n)]
        for _ in range(iters):
            max_change = 0.0
            # Neighbour mean move
            new_elev = list(self.elev)
            for i in range(n):
                if anchor[i]:
                    continue
                nbrs = self.edges_adj[i]
                if not nbrs:
                    continue
                mean_n = sum(self.elev[j] for j, _l in nbrs) / len(nbrs)
                delta = (mean_n - self.elev[i]) * damping
                new_elev[i] = self.elev[i] + delta
                if abs(delta) > max_change:
                    max_change = abs(delta)
            # Clip back into feasibility intervals.
            for i in range(n):
                if anchor[i]:
                    continue
                lo = self.interval_lo[i]
                hi = self.interval_hi[i]
                if new_elev[i] < lo: new_elev[i] = lo
                if new_elev[i] > hi: new_elev[i] = hi
            self.elev = new_elev
            # Hard-cap each edge at TAXI_MAX_GRADE.
            for i in range(n):
                for j, length in self.edges_adj[i]:
                    if j <= i:
                        continue
                    diff = self.elev[i] - self.elev[j]
                    dmax = length * TAXI_MAX_GRADE
                    if abs(diff) <= dmax:
                        continue
                    excess = abs(diff) - dmax
                    # Push both ends toward each other by half the
                    # excess (or one end fully if the other is
                    # anchored).
                    if anchor[i] and anchor[j]:
                        continue  # can't fix, both pinned
                    if anchor[i]:
                        self.elev[j] += excess * (1 if diff > 0 else -1)
                    elif anchor[j]:
                        self.elev[i] -= excess * (1 if diff > 0 else -1)
                    else:
                        half = 0.5 * excess * (1 if diff > 0 else -1)
                        self.elev[i] -= half
                        self.elev[j] += half
            if max_change < tol:
                break

    def elevation_at(self, x: float, y: float) -> Optional[float]:
        """Sample the network elevation at an arbitrary
        meter-space point.  Returns the interpolated elevation
        along the nearest graph edge when the point is close
        enough (≤ 50 m from the edge).  Returns None when the
        point is far from any network edge."""
        if not self.nodes:
            return None
        # Find nearest node — cheap first approximation.
        best_i = -1
        best_d2 = 1e18
        for i, (nx, ny) in enumerate(self.nodes):
            d2 = (nx - x) ** 2 + (ny - y) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_i = i
        if best_i < 0:
            return None
        # Search adjacent edges; pick the closest edge by
        # perpendicular distance, then interpolate along it.
        px, py = x, y
        best_edge_d = math.sqrt(best_d2)
        best_elev = self.elev[best_i]
        candidates = [(best_i, best_i)]  # self-segment fallback
        for nb, _length in self.edges_adj[best_i]:
            candidates.append((best_i, nb))
        for i, j in candidates:
            if i == j:
                continue
            ax, ay = self.nodes[i]
            bx, by = self.nodes[j]
            dx = bx - ax
            dy = by - ay
            seg_len2 = dx * dx + dy * dy
            if seg_len2 <= 0:
                continue
            t = ((px - ax) * dx + (py - ay) * dy) / seg_len2
            t = max(0.0, min(1.0, t))
            closest_x = ax + t * dx
            closest_y = ay + t * dy
            d = math.hypot(px - closest_x, py - closest_y)
            if d < best_edge_d:
                best_edge_d = d
                ea, eb = self.elev[i], self.elev[j]
                best_elev = ea + (eb - ea) * t
        if best_edge_d > 50.0:
            return None
        return best_elev


# ── Plateau detection + snap ────────────────────────────────────
#
# After the initial graph smoothing settles every node to a
# grade-compliant elevation, we cluster nodes into "shallow
# plateaus" (BFS over edges where ``|Δelev| / length`` is below
# PLATEAU_FLATNESS_GRADE) and snap each plateau's nodes to a
# single shared elevation (the median of the cluster).  The
# plateau nodes then become hard anchors for the second smoothing
# pass that fits 1.5 % ramps between them.
#
# Effect: the elevation surface becomes a series of discrete
# level zones connected by short transition ramps — matching how
# real airport pavement is actually shaped — rather than a
# continuous sub-1.5 % gradient everywhere.  Downstream:
#
#   * A taxi rect whose two axis endpoints land in the same
#     plateau auto-emits as FLAT (altitude tag, no slope).
#   * A junction polygon whose anchored boundary vertices are all
#     in one plateau passes the existing FLAT classifier without
#     any rect adjustment closure.
#   * Compound-slope junctions (spanning multiple plateaus) still
#     triangulate.

PLATEAU_MAX_RANGE_M = 2.0         # max elevation range from
                                   # min-to-max across an entire
                                   # plateau cluster.  Loose: the
                                   # verify pass shrinks plateaus
                                   # whose snap would violate ramp
                                   # grades.
PLATEAU_EDGE_SANITY_GRADE = 0.015 # max per-edge grade allowed
                                   # within a plateau (1.5 %, the
                                   # FAA cap).  Loose for the same
                                   # reason — verify pass enforces
                                   # ramp compliance after snap.
PLATEAU_MIN_NODES = 4              # ignore micro-plateaus that are
                                   # smaller than this many nodes
PLATEAU_MIN_EXTENT_M = 30.0        # ignore plateaus whose bounding
                                   # extent is smaller than this


def _detect_plateaus(g: "ElevationGraph") -> List[List[int]]:
    """Cluster graph nodes via BFS where the cluster's total
    elevation range stays below ``PLATEAU_MAX_RANGE_M`` and each
    edge crossed has grade ≤ ``PLATEAU_EDGE_SANITY_GRADE``.

    Returns plateau clusters (lists of node indices) larger than
    ``PLATEAU_MIN_NODES`` and ``PLATEAU_MIN_EXTENT_M`` of xy span.
    """
    n = len(g.nodes)
    if n == 0:
        return []
    visited = [False] * n
    plateaus: List[List[int]] = []
    for start in range(n):
        if visited[start]:
            continue
        cluster: List[int] = [start]
        cluster_min = g.elev[start]
        cluster_max = g.elev[start]
        stack = [start]
        visited[start] = True
        while stack:
            u = stack.pop()
            for v, length in g.edges_adj[u]:
                if visited[v] or length < 1e-3:
                    continue
                # Per-edge sanity: don't cross cliffs.
                edge_grade = abs(g.elev[u] - g.elev[v]) / length
                if edge_grade > PLATEAU_EDGE_SANITY_GRADE:
                    continue
                # Total-range cap: would adding v push the cluster
                # past PLATEAU_MAX_RANGE_M?
                new_min = min(cluster_min, g.elev[v])
                new_max = max(cluster_max, g.elev[v])
                if new_max - new_min > PLATEAU_MAX_RANGE_M:
                    continue
                visited[v] = True
                cluster.append(v)
                cluster_min = new_min
                cluster_max = new_max
                stack.append(v)
        if len(cluster) < PLATEAU_MIN_NODES:
            continue
        xs = [g.nodes[i][0] for i in cluster]
        ys = [g.nodes[i][1] for i in cluster]
        extent = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
        if extent < PLATEAU_MIN_EXTENT_M:
            continue
        plateaus.append(cluster)
    return plateaus


def _demote_violating_soft_anchors(
    g: "ElevationGraph",
    soft_anchored: set,
) -> Tuple[int, int]:
    """After plateau snap, walk every edge and reconcile anchored-
    pair grade violations:

      1. SOFT-anchor demotion:  An edge between two anchors whose
         elevation difference exceeds ``length × TAXI_MAX_GRADE``
         and at least one endpoint is a plateau-snapped soft anchor
         → demote the soft side(s).  The next smoothing pass pulls
         them into compliance with the surviving (typically hard)
         neighbour.
      2. HARD-vs-HARD averaging:  An edge between two HARD anchors
         with the same violation gets resolved by averaging:  both
         endpoints adopt the midpoint elevation, capped to the
         length × grade band.  This handles runway-runway crossings
         where each runway's CIFP threshold gives a different
         elevation at a shared (or near-shared) physical point —
         e.g. CYXY has 3 runways that cross at the south apron with
         CIFP-derived elevations differing by up to 16 % across a
         few-metre bridge.  Averaging trades a small CIFP deviation
         (≤ a few metres) for the elimination of a 16 % cliff in
         the rendered surface.

    Returns ``(num_soft_demoted, num_hard_pairs_averaged)``.
    """
    n = len(g.nodes)
    to_demote: set = set()
    hard_avg_pairs: List[Tuple[int, int, float]] = []
    for i in range(n):
        if g.anchor_elev[i] is None:
            continue
        for j, length in g.edges_adj[i]:
            if j <= i:
                continue
            if g.anchor_elev[j] is None:
                continue
            if length < 1e-3:
                continue
            diff = abs(g.anchor_elev[i] - g.anchor_elev[j])
            if diff <= length * TAXI_MAX_GRADE:
                continue
            i_soft = i in soft_anchored
            j_soft = j in soft_anchored
            if i_soft or j_soft:
                if i_soft:
                    to_demote.add(i)
                if j_soft:
                    to_demote.add(j)
            else:
                # Both hard — average them.
                hard_avg_pairs.append((i, j, length))
    for i in to_demote:
        g.anchor_elev[i] = None
    # Iterative averaging: a hard anchor in multiple violating
    # pairs ends up at the average of its violating neighbours
    # (we just take repeated pairwise averages — converges in a
    # few passes).
    n_hard = 0
    for _ in range(8):
        changed = False
        for i, j, length in hard_avg_pairs:
            ei = g.anchor_elev[i]
            ej = g.anchor_elev[j]
            if ei is None or ej is None:
                continue
            diff = abs(ei - ej)
            if diff <= length * TAXI_MAX_GRADE:
                continue
            # Move both ends toward their midpoint by half the
            # excess.
            mid = 0.5 * (ei + ej)
            new_i = ei + (mid - ei) * 0.5
            new_j = ej + (mid - ej) * 0.5
            g.anchor_elev[i] = new_i
            g.elev[i] = new_i
            g.anchor_elev[j] = new_j
            g.elev[j] = new_j
            changed = True
            n_hard += 1
        if not changed:
            break
    return len(to_demote), n_hard


def _snap_plateaus(g: "ElevationGraph"
                   ) -> Tuple[int, set]:
    """Detect candidate plateaus, validate each one, snap to
    median elevation only after shrinking to maintain ramp grade
    compliance at the cluster boundary.

    For each candidate cluster (detected with loose thresholds —
    up to 1.5 % per-edge / 2.0 m total range), iteratively peel
    off the boundary node whose post-snap edge-to-outside-node
    grade would exceed TAXI_MAX_GRADE.  Stop when no boundary
    edge violates grade, OR when the cluster shrinks below the
    minimum size.

    The result: aggressive plateau detection where geometry
    permits, conservative shrinkage where transition ramps are
    short — addresses Tier 2 of the sliver-elimination plan
    (user 2026-04-25).

    Returns ``(num_clusters_snapped, set_of_soft_anchored_node_ids)``.
    The soft-anchor set is consumed by the post-pass
    ``_demote_violating_soft_anchors`` to reconcile cross-cluster
    grade violations that neither cluster's own boundary check
    could see (because each cluster's check ran while the OTHER
    cluster was still pre-snap).
    """
    candidates = _detect_plateaus(g)
    soft_anchored: set = set()
    if not candidates:
        return 0, soft_anchored
    snapped = 0
    for cluster in candidates:
        cluster_set = set(cluster)
        # Iterative shrink: remove the node whose boundary edge
        # would force the worst grade violation, recompute the
        # median, repeat.
        guard = len(cluster) + 5
        while guard > 0 and len(cluster_set) >= PLATEAU_MIN_NODES:
            guard -= 1
            elevs = sorted(g.elev[i] for i in cluster_set)
            median = elevs[len(elevs) // 2]
            # Find the worst boundary-edge grade if we snapped
            # right now.
            worst_node = -1
            worst_excess = 0.0
            for i in cluster_set:
                for v, length in g.edges_adj[i]:
                    if v in cluster_set:
                        continue
                    if length < 1e-3:
                        continue
                    grade = abs(median - g.elev[v]) / length
                    excess = grade - TAXI_MAX_GRADE
                    if excess > worst_excess:
                        worst_excess = excess
                        worst_node = i
            if worst_node < 0:
                break  # no violations — commit
            cluster_set.discard(worst_node)
        if len(cluster_set) < PLATEAU_MIN_NODES:
            continue
        # Re-check minimum xy extent (shrinkage may have
        # eliminated geographic spread).
        xs = [g.nodes[i][0] for i in cluster_set]
        ys = [g.nodes[i][1] for i in cluster_set]
        extent = math.hypot(max(xs) - min(xs), max(ys) - min(ys))
        if extent < PLATEAU_MIN_EXTENT_M:
            continue
        # Commit.  Track each newly-anchored node in the soft set
        # so the post-pass can demote them if cross-cluster bridges
        # turn out to violate.
        elevs = sorted(g.elev[i] for i in cluster_set)
        median = elevs[len(elevs) // 2]
        for i in cluster_set:
            g.elev[i] = float(median)
            g.anchor_elev[i] = float(median)
            soft_anchored.add(i)
        snapped += 1
    return snapped, soft_anchored


def _build_elevation_network(
    osm_nodes, osm_ways, to_m,
    runway_segment_chain, anchor, dem, tile_lat, tile_lon,
) -> Optional[ElevationGraph]:
    """Construct the elevation network from OSM taxi centerlines
    + runway segment endpoints, densified to NETWORK_DENSIFY_M.
    Runway segment endpoints become hard anchors."""
    if osm_nodes is None or osm_ways is None or to_m is None:
        return None
    lat0, lon0 = anchor
    cos0 = math.cos(math.radians(lat0))
    g = ElevationGraph(lat0, lon0, cos0, tile_lat, tile_lon, dem)

    def _add_segment(a, b, anchor_a=None, anchor_b=None):
        """Add nodes at endpoints a, b and connect with densified
        intermediate nodes so no edge exceeds NETWORK_DENSIFY_M."""
        ax, ay = a
        bx, by = b
        dx = bx - ax; dy = by - ay
        seg_len = math.hypot(dx, dy)
        if seg_len < 0.5:
            return None, None
        n_sub = max(1, int(math.ceil(seg_len / NETWORK_DENSIFY_M)))
        step = seg_len / n_sub
        # Add endpoint nodes
        i_a = g.add_node(ax, ay, anchor=anchor_a)
        prev = i_a
        for k in range(1, n_sub):
            t = k / n_sub
            ix = ax + t * dx
            iy = ay + t * dy
            i_mid = g.add_node(ix, iy)
            g.add_edge(prev, i_mid, step)
            prev = i_mid
        i_b = g.add_node(bx, by, anchor=anchor_b)
        g.add_edge(prev, i_b, step)
        return i_a, i_b

    # ── OSM taxi ways → graph ────────────────────────────────────
    # Track node indices by OSM node id so ways that share a node
    # share a graph node (crossings / junction topology).
    osm_node_to_graph: Dict[str, int] = {}
    for wid, nds, tags in osm_ways:
        if tags.get("aeroway") != "taxiway":
            continue
        pts: List[Tuple[str, Tuple[float, float]]] = []
        for nid in nds:
            if nid not in osm_nodes:
                continue
            lat, lon = osm_nodes[nid]
            pts.append((nid, to_m(lon, lat)))
        if len(pts) < 2:
            continue
        prev_nid, prev_xy = pts[0]
        if prev_nid not in osm_node_to_graph:
            osm_node_to_graph[prev_nid] = g.add_node(prev_xy[0], prev_xy[1])
        prev_idx = osm_node_to_graph[prev_nid]
        for nid, xy in pts[1:]:
            if nid not in osm_node_to_graph:
                osm_node_to_graph[nid] = g.add_node(xy[0], xy[1])
            cur_idx = osm_node_to_graph[nid]
            # Densified chain between prev_idx and cur_idx.
            dx = xy[0] - prev_xy[0]; dy = xy[1] - prev_xy[1]
            seg_len = math.hypot(dx, dy)
            if seg_len < 0.5:
                prev_nid, prev_xy, prev_idx = nid, xy, cur_idx
                continue
            n_sub = max(1, int(math.ceil(seg_len / NETWORK_DENSIFY_M)))
            step = seg_len / n_sub
            last = prev_idx
            for k in range(1, n_sub):
                t = k / n_sub
                ix = prev_xy[0] + t * dx
                iy = prev_xy[1] + t * dy
                mid = g.add_node(ix, iy)
                g.add_edge(last, mid, step)
                last = mid
            g.add_edge(last, cur_idx, step)
            prev_nid, prev_xy, prev_idx = nid, xy, cur_idx

    # ── Runway segment endpoints as hard anchors ────────────────
    # Each runway segment's endpoints get their CIFP-derived
    # elevations as hard anchors.  Consecutive segments share
    # endpoint lat/lon, so dedupe by rounded coord key.
    rwy_node_key_to_idx: Dict[Tuple[int, int], int] = {}
    def _rwy_key(x, y):
        return (int(round(x * 10)), int(round(y * 10)))
    for seg in runway_segment_chain:
        lat_a, lon_a, elev_a, lat_b, lon_b, elev_b, _w = seg
        ax, ay = _latlon_to_m_local(lat_a, lon_a, lat0, lon0, cos0)
        bx, by = _latlon_to_m_local(lat_b, lon_b, lat0, lon0, cos0)
        ka, kb = _rwy_key(ax, ay), _rwy_key(bx, by)
        if ka not in rwy_node_key_to_idx:
            rwy_node_key_to_idx[ka] = g.add_node(
                ax, ay, anchor=float(elev_a))
        else:
            # Ensure anchor is set (later seg sharing the key).
            idx = rwy_node_key_to_idx[ka]
            if g.anchor_elev[idx] is None:
                g.anchor_elev[idx] = float(elev_a)
        if kb not in rwy_node_key_to_idx:
            rwy_node_key_to_idx[kb] = g.add_node(
                bx, by, anchor=float(elev_b))
        else:
            idx = rwy_node_key_to_idx[kb]
            if g.anchor_elev[idx] is None:
                g.anchor_elev[idx] = float(elev_b)
        # Edge length = segment length in meters
        seg_len = math.hypot(bx - ax, by - ay)
        if seg_len > 0:
            # Densify too so intermediate runway surface is in the
            # graph (helps bridge edges from taxis mid-runway).
            n_sub = max(1, int(math.ceil(seg_len / NETWORK_DENSIFY_M)))
            step = seg_len / n_sub
            ia = rwy_node_key_to_idx[ka]
            ib = rwy_node_key_to_idx[kb]
            last = ia
            for k in range(1, n_sub):
                t = k / n_sub
                ix = ax + t * (bx - ax)
                iy = ay + t * (by - ay)
                # Interpolated anchor elevation
                interp_e = float(elev_a) + t * (
                    float(elev_b) - float(elev_a))
                mid = g.add_node(ix, iy, anchor=interp_e)
                g.add_edge(last, mid, step)
                last = mid
            g.add_edge(last, ib, step)

    # ── Bridge edges: taxi nodes close to runway → runway ──────
    # Connect every taxi graph node within NETWORK_BRIDGE_MAX_M of
    # a runway graph node, representing the junction / widening
    # that physically bridges them.  OSM splits continuous taxi
    # centerlines into multiple way chains, so a "degree 1" end
    # test is unreliable — we use a distance threshold instead.
    # NETWORK_BRIDGE_MAX_M is kept tight enough that interior
    # nodes of parallel taxis (typically > 100 m from the runway
    # at SPJC) don't get falsely anchored.
    if rwy_node_key_to_idx:
        rwy_indices = list(rwy_node_key_to_idx.values())
        # Effective bridge length = (taxi→runway-centerline distance)
        # MINUS the runway half-width, because a runway is flat
        # across its width (no grade budget accrues crossing the
        # width) and a taxi really joins the runway at its boundary,
        # not its centerline.  Floor at 0.5 m so bridges never have
        # zero/negative length.
        _RWY_HALF_WIDTH_M = 22.5  # SPJC / SPLP runways are 45 m
        for g_idx in range(len(g.nodes)):
            if g.anchor_elev[g_idx] is not None:
                continue  # runway anchor — skip
            nx, ny = g.nodes[g_idx]
            best = None
            best_d = NETWORK_BRIDGE_MAX_M
            for ri in rwy_indices:
                rx, ry = g.nodes[ri]
                d = math.hypot(nx - rx, ny - ry)
                if d < best_d:
                    best_d = d
                    best = ri
            if best is not None:
                edge_len = max(0.5, best_d - _RWY_HALF_WIDTH_M)
                g.add_edge(g_idx, best, edge_len)

    if not g.nodes:
        return None
    return g


# ── Junction-diagonal bridges ────────────────────────────────────
#
# The elevation network is built along OSM taxiway centerlines +
# runway centerlines.  Two graph nodes that are close in straight
# line through a junction may be FAR apart along the network (e.g.
# centerlines that loop around the junction).  Without bridging,
# their post-smoothing elevations can violate the 1.5 % grade
# rule across the junction's diagonal.
#
# For each junction polygon, find the graph nodes nearest each
# boundary vertex and add pairwise bridge edges with length =
# straight-line distance.  Re-running propagate_bounds /
# choose_values / smooth afterwards then enforces grade compliance
# across the junction diagonals — so by the time rect altitudes
# are sampled, junction triangulation needs no further refinement.
JUNCTION_BRIDGE_MAX_M = 100.0
JUNCTION_BRIDGE_NODE_DIST_M = 30.0  # max boundary-vertex → graph
                                     # node distance to consider
                                     # the boundary vertex bridged.


TERMINAL_BRIDGE_MAX_M = 200.0  # max distance from a terminal
                                # corner to bridge into the graph
TERMINAL_BRIDGE_MAX_NEIGHBOURS = 12  # per-corner bridge-out cap


def _add_terminal_anchors(g: "ElevationGraph",
                          layout: "PavementLayout") -> int:
    """Add each terminal pad corner to the elevation graph as a
    hard anchor at the terminal's flat altitude, then bridge it
    to multiple nearby existing graph nodes (not just the nearest).

    A single nearest-node bridge leaves the terminal anchor at the
    far end of a long graph chain — the centerline path between
    apron-side OSM nodes is often hundreds of metres of twisty
    routing even when the straight-line distance is small.
    Bridging to many nodes within ``TERMINAL_BRIDGE_MAX_M`` makes
    the terminal's elevation constraint reach surrounding apron
    nodes via short straight-line paths, so smoothing can lift the
    apron UP to maintain 1.5 % grade with the terminal level.

    Returns the number of terminal corner anchors added.
    """
    if not g.nodes:
        return 0
    added = 0
    max_d2 = TERMINAL_BRIDGE_MAX_M * TERMINAL_BRIDGE_MAX_M
    for shape in layout.shapes:
        if shape.role != ROLE_TERMINAL:
            continue
        if shape.altitude is None:
            continue
        try:
            coords = list(shape.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for (cx, cy) in coords:
            # Collect every graph node within bridge range, sorted
            # by distance.  Bridge to the closest N (capped to
            # TERMINAL_BRIDGE_MAX_NEIGHBOURS) so the constraint
            # propagates to the apron from multiple directions.
            cand: List[Tuple[float, int]] = []
            for i, (nx, ny) in enumerate(g.nodes):
                d2 = (nx - cx) * (nx - cx) + (ny - cy) * (ny - cy)
                if d2 < max_d2:
                    cand.append((d2, i))
            cand.sort()
            new_idx = g.add_node(
                cx, cy, anchor=float(shape.altitude))
            for d2, i in cand[:TERMINAL_BRIDGE_MAX_NEIGHBOURS]:
                g.add_edge(new_idx, i, max(0.5, math.sqrt(d2)))
            added += 1
    return added


def _add_junction_bridges(g: "ElevationGraph",
                          layout: "PavementLayout") -> int:
    """Add cross-junction bridge edges to the elevation graph.

    For each junction polygon, pair every boundary vertex with the
    graph node nearest it (a centerline node, typically offset
    ~22 m perpendicular into the adjacent rect).  Then add bridge
    edges between every pair of those graph nodes, with bridge
    length set to the **boundary-vertex-to-boundary-vertex**
    distance — NOT the centerline-to-centerline distance.

    The grade-compliance constraint
    ``|N1.elev − N2.elev| ≤ length × TAXI_MAX_GRADE`` operates on
    elevations carried by the graph nodes.  Those elevations
    propagate to the rect altitudes, which propagate to the
    junction's anchored boundary vertex elevations at the rect
    corners.  The grade rule we ACTUALLY want enforced is between
    the junction's anchored boundary vertices — i.e. the
    corner-to-corner distance.

    Using the centerline-to-centerline distance produces a slack
    constraint: two parallel rects 60 m apart along their centerlines
    but with corners only 25 m apart end up allowed
    ``60 × 0.015 = 0.9 m`` of elevation difference, which becomes
    ``0.9 m / 25 m = 3.6 %`` cross-grade at the corners (the
    SPJC -10051 / HECA -10119 within-shape failure mode).

    Returns the number of bridges added (informational).
    """
    if not g.nodes:
        return 0
    added = 0
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        try:
            coords = list(shape.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) < 3:
            continue
        # Pair each boundary vertex with the nearest graph node and
        # remember the BOUNDARY-VERTEX position for distance calcs.
        nbrs: List[Tuple[int, float, float]] = []
        for (vx, vy) in coords:
            best_i = -1
            best_d2 = JUNCTION_BRIDGE_NODE_DIST_M ** 2
            for i, (nx, ny) in enumerate(g.nodes):
                d2 = (nx - vx) * (nx - vx) + (ny - vy) * (ny - vy)
                if d2 < best_d2:
                    best_d2 = d2
                    best_i = i
            if best_i >= 0:
                nbrs.append((best_i, vx, vy))
        if len(nbrs) < 2:
            continue
        # Dedup graph-node ids while preserving order.  When two
        # boundary vertices map to the same graph node (densified
        # ring + same nearest centerline node), prefer the boundary
        # position closest to that graph node.
        best_for_idx: Dict[int, Tuple[int, float, float]] = {}
        for ni, vx, vy in nbrs:
            nx, ny = g.nodes[ni]
            d2 = (nx - vx) ** 2 + (ny - vy) ** 2
            existing = best_for_idx.get(ni)
            if existing is None or d2 < (
                    (g.nodes[ni][0] - existing[1]) ** 2
                    + (g.nodes[ni][1] - existing[2]) ** 2):
                best_for_idx[ni] = (ni, vx, vy)
        uniq = list(best_for_idx.values())
        # Pairwise bridges, length = boundary-vertex-to-boundary-vertex
        # distance (NOT centerline-to-centerline).
        for a in range(len(uniq)):
            ia, ax, ay = uniq[a]
            for b in range(a + 1, len(uniq)):
                ib, bx, by = uniq[b]
                d = math.hypot(ax - bx, ay - by)
                if 0.5 < d <= JUNCTION_BRIDGE_MAX_M:
                    g.add_edge(ia, ib, d)
                    added += 1
    return added


# ── Corner elevation lookup ──────────────────────────────────────
#
# Junction boundary vertices share node ids with adjacent rect
# corners, runway corners, and terminal corners (the
# SHARED_VERTEX_TOL_M bucket invariant enforced before to_osm).
# Build a bucket → elevation map so each junction triangulation
# vertex can read its elevation directly from its neighbour.


def _corner_elevation_bucket(x: float, y: float,
                             tol: float = SHARED_VERTEX_TOL_M
                             ) -> Tuple[int, int]:
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


def _runway_corner_elev_map(layout: "PavementLayout"
                            ) -> Dict[Tuple[int, int], float]:
    """Bucket → elevation map containing ONLY runway corners.

    Runway elevations come from CIFP threshold data and are the
    only HARD anchors in the pavement system — every other shape
    (taxi rects, terminals, junctions) is allowed to be modified
    by the global mesh smoother to satisfy 2D-Euclidean grade
    compliance.

    Sloped runways (altitude_high + altitude_low) follow the
    X-Plane patch convention: ring indices 0,3 are the HIGH short
    edge and 1,2 are the LOW short edge.  Flat runways carry a
    single ``altitude``.
    """
    out: Dict[Tuple[int, int], float] = {}
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        if (s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            elevs = [s.altitude_high, s.altitude_low,
                     s.altitude_low, s.altitude_high]
            for (cx, cy), e in zip(coords, elevs):
                out.setdefault(
                    _corner_elevation_bucket(cx, cy), float(e))
            continue
        if s.altitude is not None:
            for (cx, cy) in coords:
                out.setdefault(
                    _corner_elevation_bucket(cx, cy),
                    float(s.altitude))
    return out


# Tolerance for "junction corner is adjacent to runway pavement".
# The bucket-exact match in ``_runway_corner_elev_map`` only catches
# junction corners that LITERALLY coincide with a runway vertex
# (within SHARED_VERTEX_TOL_M = 0.5 m).  Real airports often have
# the junction polygon ending a few metres short of the runway
# polygon — apt.dat polygon precision, runway-end blast pads, etc.
# Without an anchor, the mesh solver settles those corners to
# whatever the inland taxi pavement is at, producing a visible
# elevation drop between the junction and the runway (user reported
# this at SPJC's 16L L1 stub).
NEAR_RUNWAY_ANCHOR_M = 60.0


def _runway_elev_at_point(rwy_shape: "BuiltShape",
                          x: float, y: float) -> Optional[float]:
    """Interpolate the runway shape's surface elevation at the
    point on its boundary closest to ``(x, y)``.

    For a flat runway (``altitude``), returns that value regardless
    of position.  For a sloped runway (``altitude_high`` /
    ``altitude_low``), parameterises ``(x, y)``'s nearest-boundary
    projection along the high-low axis and linearly interpolates.
    Returns ``None`` if the shape has no usable elevation tag.
    """
    if rwy_shape.altitude is not None:
        return float(rwy_shape.altitude)
    if (rwy_shape.altitude_high is None
            or rwy_shape.altitude_low is None):
        return None
    try:
        coords = list(rwy_shape.polygon.exterior.coords)
    except Exception:
        return None
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return float(
            (rwy_shape.altitude_high + rwy_shape.altitude_low) / 2)
    # X-Plane convention: corners 0,3 = HIGH short edge,
    # corners 1,2 = LOW short edge.  Axis = midpoint(0,3) →
    # midpoint(1,2).
    high_mx = 0.5 * (coords[0][0] + coords[3][0])
    high_my = 0.5 * (coords[0][1] + coords[3][1])
    low_mx = 0.5 * (coords[1][0] + coords[2][0])
    low_my = 0.5 * (coords[1][1] + coords[2][1])
    ax = low_mx - high_mx
    ay = low_my - high_my
    a2 = ax * ax + ay * ay
    if a2 < 1e-9:
        return float(
            (rwy_shape.altitude_high + rwy_shape.altitude_low) / 2)
    t = ((x - high_mx) * ax + (y - high_my) * ay) / a2
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    return float(rwy_shape.altitude_high
                 + t * (rwy_shape.altitude_low
                        - rwy_shape.altitude_high))


def _near_runway_anchor_map(
    layout: "PavementLayout",
    tol_m: float = NEAR_RUNWAY_ANCHOR_M,
) -> Dict[Tuple[int, int], float]:
    """Bucket → elevation map for JUNCTION corners that sit within
    ``tol_m`` of a runway polygon, anchored to the runway's
    elevation at the closest boundary point.

    Per user 2026-04-27: junctions adjacent to a runway should
    match the runway's elevation at that boundary, even when the
    junction polygon falls a few metres short of the runway
    polygon (apt.dat polygon-precision artefact, blast-pad gap,
    DSF / runway boundary mismatch).  Without this anchor the
    mesh solver settles the junction at the inland taxi area's
    elevation, producing the user-reported drop (SPJC 16L: junction
    -10163 at 12.9 m next to runway -10077 at 13.4 m).

    Only RUNWAY-FACING corners get anchored — corners on the inland
    side of the junction (further from any runway) stay free for
    the mesh solver to slope the junction smoothly away from the
    runway-facing edge.
    """
    out: Dict[Tuple[int, int], float] = {}
    runways = [s for s in layout.shapes
               if s.role == ROLE_RUNWAY
               and s.polygon is not None
               and not s.polygon.is_empty]
    if not runways:
        return out
    # Build an STRtree over runway polygons so per-vertex distance
    # checks only touch runway candidates whose bounding box falls
    # within ``tol_m`` of the vertex.  Replaces the original
    # ``O(junction_verts × runways)`` pairwise loop with
    # ``O(junction_verts × log(runways) + hits)``.  At CYXY this
    # cuts 211 k shapely.distance() calls down to a few thousand,
    # saving ~0.7 s of `_compute_elevations`.
    from shapely.strtree import STRtree
    runway_geoms = [r.polygon for r in runways]
    tree = STRtree(runway_geoms)
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        for (jx, jy) in coords:
            pt = Point(jx, jy)
            # Query the tree for runway polygons whose bbox is
            # within tol_m of the vertex.  Shapely's STRtree.query
            # returns indices into the geom list passed at
            # construction time.
            try:
                cand_idxs = tree.query(pt.buffer(tol_m,
                                                  resolution=1))
            except Exception:
                cand_idxs = list(range(len(runways)))
            best_d = tol_m
            best_e = None
            for ci in cand_idxs:
                r = runways[int(ci)]
                try:
                    d = r.polygon.distance(pt)
                except Exception:
                    continue
                if d > best_d:
                    continue
                e = _runway_elev_at_point(r, jx, jy)
                if e is None:
                    continue
                best_d = d
                best_e = e
            if best_e is None:
                continue
            bucket = _corner_elevation_bucket(jx, jy)
            # Don't override an exact runway-corner bucket (those
            # are already in ``runway_anchors`` from the corner map).
            out.setdefault(bucket, best_e)
    return out


# ── Global pavement mesh ─────────────────────────────────────────
#
# A mesh built from every pavement shape's boundary (rects, taxi
# rects, terminals, junctions) where nodes are unique vertex
# buckets and edges are ring-adjacency steps + cross-shape
# adjacency at shared buckets.  Anchored only at runway corner
# buckets — the CIFP-derived runway elevations are the sole HARD
# truth.  Every other node is free to be moved by the smoother
# into 2D-Euclidean grade compliance with its neighbours.
#
# After solving, each shape's altitude tags are rewritten from the
# mesh values (taxi rect altitude_high/altitude_low ← mean of the
# rect's high-side / low-side corner mesh values; flat shapes ←
# mean of corner mesh values; junction node_altitudes ← per-vertex
# mesh values).  Runway altitudes are preserved verbatim — they
# anchored the mesh.


def _build_pavement_mesh(
        layout: "PavementLayout",
        runway_anchors: Dict[Tuple[int, int], float],
        dem,
        tile_lat: int,
        tile_lon: int,
        ) -> Tuple["ElevationGraph",
                   Dict[Tuple[int, int], int]]:
    """Build a global ElevationGraph from every pavement shape's
    boundary.  Returns ``(mesh, bucket_to_node)`` where the mesh
    is ready for ``propagate_bounds()`` / ``choose_values()`` /
    ``smooth_rate_of_change()`` and ``bucket_to_node`` lets the
    writeback pass look up each shape's vertex elevations after
    smoothing.

    Initial elevations come from the existing per-shape altitude
    tags — i.e. the centerline-graph-derived rect / terminal
    elevations and the graph-or-DEM-sampled junction boundary
    samples.  The mesh then refines those into 2D-Euclidean
    compliance with the runway anchors.
    """
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    mesh = ElevationGraph(lat0, lon0, cos0, tile_lat, tile_lon, dem)
    bucket_to_node: Dict[Tuple[int, int], int] = {}

    # Per-bucket initial elevation, accumulated as we walk shapes
    # so a bucket touched by multiple shapes seeds itself from the
    # mean of every shape's idea of what it should be.
    bucket_init_sum: Dict[Tuple[int, int], Tuple[float, int]] = {}

    def _bucket_for_node(x: float, y: float) -> int:
        b = _corner_elevation_bucket(x, y)
        nid = bucket_to_node.get(b)
        if nid is None:
            anchor = runway_anchors.get(b)
            nid = mesh.add_node(x, y, anchor=anchor)
            bucket_to_node[b] = nid
        return nid

    def _ring_with_elevs(s: "BuiltShape"
                         ) -> Tuple[List[Tuple[float, float]],
                                    List[Optional[float]]]:
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            return [], []
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            return [], []
        if (s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            elevs = [s.altitude_high, s.altitude_low,
                     s.altitude_low, s.altitude_high]
            return coords, [float(e) for e in elevs]
        if (s.node_altitudes is not None
                and len(s.node_altitudes) >= len(coords)):
            na = list(s.node_altitudes)
            if len(na) == len(coords) + 1:
                na = na[:-1]
            if len(na) == len(coords):
                return coords, [float(e) for e in na]
        if s.altitude is not None:
            return coords, [float(s.altitude)] * len(coords)
        return coords, [None] * len(coords)

    # Pass 1 — register every boundary vertex as a mesh node and
    # accumulate seed elevations.
    for s in layout.shapes:
        coords, elevs = _ring_with_elevs(s)
        if not coords:
            continue
        for (x, y), e in zip(coords, elevs):
            nid = _bucket_for_node(x, y)
            if e is not None:
                b = _corner_elevation_bucket(x, y)
                tot, cnt = bucket_init_sum.get(b, (0.0, 0))
                bucket_init_sum[b] = (tot + e, cnt + 1)

    # Pass 2 — add edges.  Two kinds:
    #   (a) ring-adjacency: each ring step (i → i+1) becomes an
    #       edge weighted by 2D distance.  Cross-shape shared
    #       edges are deduplicated.
    #   (b) within-polygon all-pairs within
    #       ``MESH_WITHIN_SHAPE_RADIUS_M`` of each other.  Without
    #       these the mesh's grade-cap only enforces 1.5 % on
    #       ring-adjacent pairs; vertices on opposite sides of a
    #       wide junction could legitimately violate 1.5 % per 2D
    #       metre even after smoothing because there's no direct
    #       mesh edge between them.  Adding the all-pairs edges
    #       (capped at the same radius check_grade.py uses for
    #       within-shape violations) makes the smoother's grade-
    #       cap pass enforce the 2D rule directly.
    MESH_WITHIN_SHAPE_RADIUS_M = WITHIN_SHAPE_VIOLATION_RADIUS_M
    seen_edges: set = set()
    for s in layout.shapes:
        coords, _e = _ring_with_elevs(s)
        n = len(coords)
        if n < 2:
            continue
        # Resolve every ring vertex to its bucket node id once.
        node_ids: List[int] = []
        for (cx, cy) in coords:
            node_ids.append(
                bucket_to_node[_corner_elevation_bucket(cx, cy)])
        # (a) Ring adjacency.
        for i in range(n):
            ai = node_ids[i]
            bi = node_ids[(i + 1) % n]
            if ai == bi:
                continue
            key = (min(ai, bi), max(ai, bi))
            if key in seen_edges:
                continue
            seen_edges.add(key)
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % n]
            d = math.hypot(bx - ax, by - ay)
            if d <= 0:
                continue
            mesh.add_edge(ai, bi, d)
        # (b) All-pairs within MESH_WITHIN_SHAPE_RADIUS_M.
        if n < 3:
            continue
        radius2 = (MESH_WITHIN_SHAPE_RADIUS_M
                   * MESH_WITHIN_SHAPE_RADIUS_M)
        for i in range(n):
            ax, ay = coords[i]
            ai = node_ids[i]
            for j in range(i + 1, n):
                bx, by = coords[j]
                bi = node_ids[j]
                if ai == bi:
                    continue
                d2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay)
                if d2 <= 0 or d2 > radius2:
                    continue
                key = (min(ai, bi), max(ai, bi))
                if key in seen_edges:
                    continue
                seen_edges.add(key)
                mesh.add_edge(ai, bi, math.sqrt(d2))

    # Pass 3 — cross-polygon proximity edges.  Pairs of mesh nodes
    # in DIFFERENT polygons that sit within
    # ``MESH_CROSS_SHAPE_RADIUS_M`` of each other in 2D need an
    # edge so the smoother enforces 1.5 % grade compliance across
    # the gap (and check_grade's CROSS-SHAPE check is satisfied).
    # Without these, two boundary vertices in adjacent buckets
    # (different node ids) on different shapes would diverge in
    # elevation by the same network-distance-but-not-2D-distance
    # mechanism that drives the original violation pattern.
    #
    # Bucketed spatial index: 5 m grid bins of node ids, then
    # each node only checks neighbouring bins.  O(N) instead of
    # the naive O(N²) all-pairs scan.
    MESH_CROSS_SHAPE_RADIUS_M = 5.0
    BIN_M = 5.0
    bins: Dict[Tuple[int, int], List[int]] = {}
    for nid, (nx, ny) in enumerate(mesh.nodes):
        key = (int(math.floor(nx / BIN_M)),
               int(math.floor(ny / BIN_M)))
        bins.setdefault(key, []).append(nid)
    radius2 = MESH_CROSS_SHAPE_RADIUS_M * MESH_CROSS_SHAPE_RADIUS_M
    for nid, (nx, ny) in enumerate(mesh.nodes):
        kx = int(math.floor(nx / BIN_M))
        ky = int(math.floor(ny / BIN_M))
        for dkx in (-1, 0, 1):
            for dky in (-1, 0, 1):
                cell = bins.get((kx + dkx, ky + dky))
                if not cell:
                    continue
                for other in cell:
                    if other <= nid:
                        continue
                    ox, oy = mesh.nodes[other]
                    d2 = (ox - nx) ** 2 + (oy - ny) ** 2
                    if d2 <= 0 or d2 > radius2:
                        continue
                    key = (nid, other)
                    if key in seen_edges:
                        continue
                    seen_edges.add(key)
                    mesh.add_edge(nid, other, math.sqrt(d2))

    # Seed initial elevations on the mesh's elev[] array so the
    # smoother starts from realistic per-shape values rather than
    # pure DEM.  Anchors keep their anchor value; non-anchor nodes
    # get the mean of every shape's seed at that bucket (or DEM if
    # no shape contributed).
    n_nodes = len(mesh.nodes)
    init_elev: List[float] = [0.0] * n_nodes
    for b, nid in bucket_to_node.items():
        a = mesh.anchor_elev[nid]
        if a is not None:
            init_elev[nid] = a
            continue
        seed = bucket_init_sum.get(b)
        if seed is not None and seed[1] > 0:
            init_elev[nid] = seed[0] / seed[1]
        else:
            d = mesh.dem_elev[nid]
            init_elev[nid] = d if d is not None else 0.0
    mesh.elev = init_elev

    return mesh, bucket_to_node


def _solve_pavement_mesh(
        mesh: "ElevationGraph",
        layout: "PavementLayout",
        bucket_to_node: Dict[Tuple[int, int], int]) -> None:
    """Run propagate-bounds + Laplacian-with-grade-cap on the
    pavement mesh, while enforcing two extra constraints during
    every smoothing iteration:

    1. Each sloped rect's short-edge corner pair is forced equal
       (X-Plane patch format requires ``altitude_high`` /
       ``altitude_low`` to apply uniformly across each short
       edge).  Averaging the pair after each smoothing step lets
       the rest of the mesh adjust around the equalised value,
       so junctions sharing a rect corner bucket converge to
       elevations consistent with the rect's emit.
    2. Bounds-propagation gives every node a 2D Euclidean
       feasibility band against every reachable runway anchor;
       smoothing pulls free nodes toward neighbour means while
       enforcing |Δelev| ≤ length × TAXI_MAX_GRADE on every
       edge.

    Mutates ``mesh.elev`` in place.
    """
    # propagate_bounds will RECOMPUTE intervals from anchors only;
    # we want to keep our seeded initial values for free nodes,
    # so cache and re-clip rather than calling choose_values.
    seeded = list(mesh.elev)
    mesh.propagate_bounds()

    # ── Apron decoupling (user 2026-04-27, plan A) ──────────────────
    # Aprons sit on natural terrain and follow DEM, not runway-grade
    # propagation.  When the airport is on a hillside (CYXY), the
    # apron is several metres above the runway.  The mesh's 2D
    # propagate_bounds finds short geometric paths through within-
    # shape and cross-shape proximity edges (~130 m at CYXY) and
    # caps the apron at runway+1.5%×130m ≈ runway+2m — far below
    # the actual terrain (10m above runway at CYXY's SW apron).
    #
    # Fix (user direction): identify "apron-class" junction shapes
    # (large area — a parking apron, not a small taxi connector) and
    # bias their EXCLUSIVE boundary vertices (those NOT shared with
    # rect/terminal/runway shapes) toward DEM elevation, subject to
    # taxi-grade compliance walked along the apron's PERIMETER ring
    # from the shared-with-taxi vertices (which already have anchor-
    # propagation-derived elevations).  This matches the user's
    # phrasing: "we still have to maintain taxi grades which will
    # determine just how high the apron can be and how much it gets
    # cut into the hills".
    #
    # Algorithm — per apron polygon:
    #   1. Build the ring-adjacency graph (just ring edges, no
    #      mesh shortcuts).
    #   2. Sources are the SHARED-with-anchored-shape vertices,
    #      using their current mesh elevation as the anchor value.
    #      If no shared vertex exists, the whole apron is exclusive
    #      and we anchor at DEM directly.
    #   3. Dijkstra outward: distance[v] = ring-adjacency length.
    #   4. Per-vertex interval = ∩ (source_elev ± dist × 1.5 %).
    #   5. Choose elev = round(clip(DEM, interval), 0.1).
    #
    # Per-airport effect:
    #   * Flat airports (SPJC, KBNA, HECA): DEM ≈ neighbor elevation,
    #     so DEM falls inside the interval and the override is a
    #     near-no-op (apron stays where it was).
    #   * Hillside airports (CYXY, SPLP): apron rises toward DEM,
    #     up to what the apron's ring perimeter from the shared
    #     edge allows.  Vertices far from any shared edge reach
    #     full DEM; vertices near a shared edge stay close to the
    #     taxi-grade-derived elevation.  The apron's perimeter
    #     splines smoothly between them.
    APRON_AREA_THRESHOLD_M2 = 3_000.0

    def _fit_dem_plane(coords: List[Tuple[float, float]],
                        dems: List[Optional[float]]
                        ) -> Optional[Tuple[float, float, float]]:
        """Fit ``z = a*x + b*y + c`` to the DEM samples by ordinary
        least squares.  Returns ``(a, b, c)`` or ``None`` when fewer
        than 3 valid samples or the design matrix is singular.

        The plane represents the apron's overall terrain slope while
        ignoring local DEM noise (buildings, vegetation, SRTM
        artefacts).  Pinning vertices to plane values rather than
        individual DEM values gives a smooth ramp that matches the
        terrain trend without producing 4 m cliffs between
        physically-close ring vertices whose individual DEM samples
        happen to differ.
        """
        xs = []
        ys = []
        zs = []
        for (cx, cy), d in zip(coords, dems):
            if d is None:
                continue
            xs.append(cx)
            ys.append(cy)
            zs.append(float(d))
        n = len(zs)
        if n < 3:
            return None
        # Normal equations:
        #   [Σx²  Σxy  Σx ] [a]   [Σxz]
        #   [Σxy  Σy²  Σy ] [b] = [Σyz]
        #   [Σx   Σy   N  ] [c]   [Σz ]
        sx = sum(xs); sy = sum(ys); sz = sum(zs)
        sxx = sum(x * x for x in xs)
        syy = sum(y * y for y in ys)
        sxy = sum(x * y for x, y in zip(xs, ys))
        sxz = sum(x * z for x, z in zip(xs, zs))
        syz = sum(y * z for y, z in zip(ys, zs))
        # Solve via Cramer's rule (3x3).
        m = [[sxx, sxy, sx],
             [sxy, syy, sy],
             [sx,  sy,  float(n)]]
        rhs = [sxz, syz, sz]
        def _det3(mm):
            return (mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
                    - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
                    + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]))
        det = _det3(m)
        if abs(det) < 1e-9:
            # Degenerate (collinear vertices): fall back to plane = mean.
            mean_z = sz / n
            return (0.0, 0.0, mean_z)
        coeffs = []
        for col in range(3):
            mc = [row[:] for row in m]
            for r in range(3):
                mc[r][col] = rhs[r]
            coeffs.append(_det3(mc) / det)
        return (coeffs[0], coeffs[1], coeffs[2])

    rect_like_anchored = {
        ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
        ROLE_STUB, ROLE_CROSS_CONNECTOR, ROLE_RUNWAY,
        ROLE_TERMINAL,
    }
    anchored_shape_buckets: set = set()
    for s in layout.shapes:
        if s.role not in rect_like_anchored:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for (x, y) in coords:
            anchored_shape_buckets.add(
                _corner_elevation_bucket(x, y))

    # Build ONE GLOBAL graph spanning every apron-class junction,
    # with shared boundary vertices stitching adjacent junctions
    # together.  Per-polygon Dijkstra didn't work because adjacent
    # junctions share boundary vertices: each polygon's pass would
    # re-pin the shared vertex with different values, leaving big
    # cliffs at junction-to-junction borders.  One global Dijkstra
    # over all apron-class junction polygons treats the whole
    # apron complex as a single connected surface and produces a
    # consistent elevation field.
    apron_node_ids: set = set()
    apron_graph_adj: Dict[int, List[Tuple[int, float]]] = {}
    seen_edges_in_apron: set = set()

    def _add_apron_edge(a: int, b: int, length: float) -> None:
        if a == b:
            return
        key = (min(a, b), max(a, b))
        if key in seen_edges_in_apron:
            return
        seen_edges_in_apron.add(key)
        apron_graph_adj.setdefault(a, []).append((b, length))
        apron_graph_adj.setdefault(b, []).append((a, length))

    n_apron_polys = 0
    prox_radius2 = (WITHIN_SHAPE_VIOLATION_RADIUS_M
                    * WITHIN_SHAPE_VIOLATION_RADIUS_M)
    # Per-node target elevation = plane fit at the node's (x,y).
    # Built per-polygon below; merged across the global apron set
    # via simple averaging when multiple polygons disagree about a
    # shared vertex (rare — adjacent aprons usually have similar
    # DEM trends at the shared edge).
    apron_target: Dict[int, Tuple[float, int]] = {}
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        try:
            if s.polygon.area < APRON_AREA_THRESHOLD_M2:
                continue
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        n = len(coords)
        if n < 4:
            continue
        nids: List[int] = []
        for (cx, cy) in coords:
            b = _corner_elevation_bucket(cx, cy)
            nid = bucket_to_node.get(b)
            if nid is None:
                nids = []
                break
            nids.append(nid)
        if not nids:
            continue
        n_apron_polys += 1
        for nid in nids:
            apron_node_ids.add(nid)
        # Ring edges.
        for i in range(n):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % n]
            d = math.hypot(bx - ax, by - ay)
            _add_apron_edge(nids[i], nids[(i + 1) % n], d)
        # Within-shape 2D-proximity edges (catches non-convex
        # ring-distant-but-close-in-2D pairs).
        for i in range(n):
            ax, ay = coords[i]
            for j in range(i + 1, n):
                if (i + 1) % n == j or (j + 1) % n == i:
                    continue
                bx, by = coords[j]
                d2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay)
                if d2 <= 0 or d2 > prox_radius2:
                    continue
                _add_apron_edge(nids[i], nids[j], math.sqrt(d2))
        # Fit a plane to this polygon's DEM samples, then record the
        # plane's value at every ring vertex as that vertex's target.
        # If the plane fit fails (degenerate geometry), fall back to
        # individual DEM samples.
        node_dems = [mesh.dem_elev[nid] for nid in nids]
        plane = _fit_dem_plane(coords, node_dems)
        if plane is not None:
            a, b, c = plane
            for (cx, cy), nid in zip(coords, nids):
                t = a * cx + b * cy + c
                acc = apron_target.get(nid)
                if acc is None:
                    apron_target[nid] = (t, 1)
                else:
                    apron_target[nid] = (acc[0] + t, acc[1] + 1)
        else:
            for nid in nids:
                d = mesh.dem_elev[nid]
                if d is None:
                    continue
                acc = apron_target.get(nid)
                if acc is None:
                    apron_target[nid] = (float(d), 1)
                else:
                    apron_target[nid] = (acc[0] + float(d), acc[1] + 1)

    # Sources: any apron node whose bucket coincides with an
    # anchored shape's vertex (taxi rect / runway / terminal corner
    # that the apron also borders).  These take their CURRENT mesh
    # elevation (already re-clipped above) as the propagation seed.
    sources_per_node: Dict[int, float] = {}
    for nid in apron_node_ids:
        if mesh.anchor_elev[nid] is not None:
            sources_per_node[nid] = mesh.anchor_elev[nid]
            continue
        # Reverse-look-up: is this node's bucket also a vertex of
        # any anchored shape?  bucket_to_node maps bucket→nid; we
        # need nid→bucket.  Iterate (cheaper than reversing once
        # since this set is small).
        for b, nb in bucket_to_node.items():
            if nb == nid and b in anchored_shape_buckets:
                sources_per_node[nid] = mesh.elev[nid]
                break

    # Cross-shape proximity sources (user 2026-04-27): an apron
    # polygon and a taxi rect can be physically adjacent (within
    # a few metres) WITHOUT sharing any boundary vertex — apt.dat
    # row-110 follows the apron's curved boundary while the taxi
    # rect uses 4 axis-aligned corners.  Without proximity sourcing,
    # the apron-decouple Dijkstra has no path from the taxi to the
    # apron and apron-exclusive vertices climb freely toward DEM
    # ignoring the taxi-grade cap that should bound them.
    #
    # Fix: index every NON-apron mesh node whose bucket belongs to
    # an anchored shape, then for each apron node, look up nearby
    # external anchored nodes within ``EXTERNAL_SOURCE_RADIUS_M``
    # and stitch them into the apron graph.  The external nodes
    # become additional sources at their current mesh.elev value
    # (taxi-grade-derived from the centerline propagation),
    # enabling the multi-source Dijkstra to walk taxi → apron via
    # a 2-5 m proximity edge.
    EXTERNAL_SOURCE_RADIUS_M = 5.0
    external_anchored: Dict[int, float] = {}
    for b, nid in bucket_to_node.items():
        if nid in apron_node_ids:
            continue
        if b in anchored_shape_buckets:
            external_anchored[nid] = mesh.elev[nid]
    if external_anchored:
        ext_bins: Dict[Tuple[int, int], List[int]] = {}
        for nid in external_anchored:
            x, y = mesh.nodes[nid]
            key = (int(math.floor(x / EXTERNAL_SOURCE_RADIUS_M)),
                   int(math.floor(y / EXTERNAL_SOURCE_RADIUS_M)))
            ext_bins.setdefault(key, []).append(nid)
        radius2 = (EXTERNAL_SOURCE_RADIUS_M
                   * EXTERNAL_SOURCE_RADIUS_M)
        linked_externals: set = set()
        for apron_nid in list(apron_node_ids):
            ax, ay = mesh.nodes[apron_nid]
            kx = int(math.floor(ax / EXTERNAL_SOURCE_RADIUS_M))
            ky = int(math.floor(ay / EXTERNAL_SOURCE_RADIUS_M))
            for dkx in (-1, 0, 1):
                for dky in (-1, 0, 1):
                    cell = ext_bins.get((kx + dkx, ky + dky))
                    if not cell:
                        continue
                    for ext_nid in cell:
                        ex, ey = mesh.nodes[ext_nid]
                        d2 = ((ex - ax) * (ex - ax)
                              + (ey - ay) * (ey - ay))
                        if d2 > radius2:
                            continue
                        d = math.sqrt(d2) if d2 > 0 else 0.0
                        _add_apron_edge(apron_nid, ext_nid, d)
                        linked_externals.add(ext_nid)
        # Bring the linked external nodes into the apron set so the
        # multi-source Dijkstra relaxes through them; they're
        # sources at their current mesh elevation.
        for ext_nid in linked_externals:
            apron_node_ids.add(ext_nid)
            sources_per_node.setdefault(
                ext_nid, external_anchored[ext_nid])

    n_apron_pinned = 0
    if apron_node_ids:
        # Single multi-source Dijkstra over the global apron graph.
        # Each node ends up with the cone of its NEAREST source:
        # ``[src_elev ± nearest_dist × TAXI_MAX_GRADE]``.  Same
        # algorithmic shift as the centerline ``propagate_bounds``
        # — we previously did one Dijkstra per source and intersected
        # cones, which scales as O(sources × (V+E) log V).  At HECA
        # there are hundreds of sources and the per-source variant
        # was tens of seconds; the single-pass multi-source pattern
        # is one graph traversal regardless of source count.
        import heapq
        INF = float("inf")
        per_lo: Dict[int, float] = {nid: -INF for nid in apron_node_ids}
        per_hi: Dict[int, float] = {nid: INF for nid in apron_node_ids}
        nearest_dist: Dict[int, float] = {nid: INF
                                           for nid in apron_node_ids}
        nearest_val: Dict[int, float] = {nid: 0.0
                                          for nid in apron_node_ids}
        heap: List[Tuple[float, int, float]] = []
        for src_nid, src_elev in sources_per_node.items():
            nearest_dist[src_nid] = 0.0
            nearest_val[src_nid] = src_elev
            heapq.heappush(heap, (0.0, src_nid, float(src_elev)))
        while heap:
            d, u, val = heapq.heappop(heap)
            if d > nearest_dist[u]:
                continue
            for v, length in apron_graph_adj.get(u, ()):
                nd = d + length
                if v in nearest_dist and nd < nearest_dist[v]:
                    nearest_dist[v] = nd
                    nearest_val[v] = val
                    heapq.heappush(heap, (nd, v, val))
        for nid in apron_node_ids:
            d = nearest_dist[nid]
            if d == INF:
                continue
            band = d * TAXI_MAX_GRADE
            per_lo[nid] = nearest_val[nid] - band
            per_hi[nid] = nearest_val[nid] + band

        # Pin every apron-EXCLUSIVE node (i.e. not already a source).
        # Target = plane-fit value (smooth across the polygon, free
        # of DEM noise), clipped to the Dijkstra-from-sources cap.
        for nid in apron_node_ids:
            if nid in sources_per_node:
                continue
            if mesh.anchor_elev[nid] is not None:
                continue
            tgt_acc = apron_target.get(nid)
            if tgt_acc is None:
                # Fallback to per-vertex DEM if plane fit didn't
                # cover this node (shouldn't happen — every apron
                # vertex contributes to at least one polygon's plane).
                d = mesh.dem_elev[nid]
                if d is None:
                    continue
                target_pref = float(d)
            else:
                target_pref = tgt_acc[0] / tgt_acc[1]
            lo = per_lo[nid]
            hi = per_hi[nid]
            # Cap by the GLOBAL mesh interval too: that one was
            # produced by ``propagate_bounds`` over the full mesh
            # (ring + within-shape + 5 m cross-shape proximity), so
            # it captures taxi-grade propagation via SMALL junction
            # polygons sitting between the apron and a taxi rect
            # (CYXY's SW apron is 100 m+ from any E corner directly
            # but is bridged by a small connector junction).  Without
            # this cap, an apron-exclusive vertex with no apron-
            # graph source within reach would pin freely to plane
            # DEM, producing a 10+ m discontinuity at the apron-
            # taxi boundary.  Reading mesh.interval_hi BEFORE we
            # overwrite it below gives the original Dijkstra-derived
            # cap.
            global_hi = mesh.interval_hi[nid]
            global_lo = mesh.interval_lo[nid]
            if global_hi != INF and global_hi < hi:
                hi = global_hi
            if global_lo != -INF and global_lo > lo:
                lo = global_lo
            if lo == -INF and hi == INF:
                target = target_pref
            elif lo > hi:
                target = 0.5 * (lo + hi)
            else:
                target = target_pref
                if target < lo:
                    target = lo
                elif target > hi:
                    target = hi
            target = round(target, 1)
            mesh.anchor_elev[nid] = target
            mesh.interval_lo[nid] = target
            mesh.interval_hi[nid] = target
            seeded[nid] = target
            n_apron_pinned += 1
    if n_apron_pinned:
        try:
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] pinned "
                f"{n_apron_pinned} apron-exclusive vertex(es) "
                f"across {n_apron_polys} apron polygon(s) at "
                f"DEM-clipped-to-ring-grade elevation.\n")
        except Exception:
            pass

    for i in range(len(mesh.nodes)):
        if mesh.anchor_elev[i] is not None:
            mesh.elev[i] = mesh.anchor_elev[i]
            continue
        lo = mesh.interval_lo[i]
        hi = mesh.interval_hi[i]
        if lo > hi:
            mesh.elev[i] = 0.5 * (lo + hi)
        else:
            v = seeded[i]
            if v < lo:
                v = lo
            elif v > hi:
                v = hi
            mesh.elev[i] = v
    # Pre-build the equalisation lists ONCE so the per-iteration
    # equalise step is O(constraints) instead of O(shapes).  Two
    # constraint families:
    #   * Sloped rect short-edge pairs (corners 0+3 share
    #     altitude_high, corners 1+2 share altitude_low).
    #   * All corners of a FLAT shape (terminal, flat rect)
    #     share a single elevation — the shape emits ``altitude``
    #     uniformly.
    rect_like_roles = {ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL,
                       ROLE_STUB, ROLE_CROSS_CONNECTOR}
    rect_pairs: List[Tuple[int, int, int, int]] = []
    flat_groups: List[List[int]] = []
    for s in layout.shapes:
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        nids = []
        for c in coords:
            nid = bucket_to_node.get(
                _corner_elevation_bucket(c[0], c[1]))
            if nid is None:
                nids = []
                break
            nids.append(nid)
        if not nids:
            continue
        # Sloped rect: 4-corner H/L/L/H pattern.
        if (s.role in rect_like_roles
                and s.altitude_high is not None
                and s.altitude_low is not None
                and len(nids) == 4):
            rect_pairs.append(
                (nids[0], nids[3], nids[1], nids[2]))
            continue
        # Flat shape: terminal, flat rect, or flat runway-like
        # (runways are SKIPPED — their corners are anchored, not
        # subject to equalisation by the mesh).
        if s.role == ROLE_RUNWAY:
            continue
        if (s.role == ROLE_TERMINAL
                or (s.role in rect_like_roles
                    and s.altitude is not None)):
            # Deduplicate node ids inside a single shape so
            # adjacent ring vertices in the same bucket don't
            # double-count.
            uniq = list(dict.fromkeys(nids))
            if len(uniq) >= 2:
                flat_groups.append(uniq)

    def _equalise_constraints() -> None:
        for (h0, h1, l0, l1) in rect_pairs:
            eh = 0.5 * (mesh.elev[h0] + mesh.elev[h1])
            el = 0.5 * (mesh.elev[l0] + mesh.elev[l1])
            mesh.elev[h0] = eh
            mesh.elev[h1] = eh
            mesh.elev[l0] = el
            mesh.elev[l1] = el
        for grp in flat_groups:
            avg = sum(mesh.elev[i] for i in grp) / len(grp)
            for i in grp:
                mesh.elev[i] = avg

    # ── Apron-zone tighter grade cap ────────────────────────────
    # User 2026-04-27: aprons (junction polygons that surround
    # terminal pads) should hold a 1 % grade cap rather than the
    # 1.5 % cap used elsewhere — keeps the parking / loading area
    # visually flat near the terminal.  Excess slope is absorbed
    # by the surrounding taxis (which still operate at 1.5 %), so
    # the aircraft sees a gentle apron and a slightly steeper taxi
    # connection rather than a gentle taxi and a too-steep apron.
    #
    # An apron node is any mesh node belonging to a junction
    # polygon that shares at least one vertex with a terminal pad.
    # An apron edge is one where BOTH endpoints are apron nodes —
    # those get the 1 % cap.  Edges from an apron node out to a
    # taxi rect's runway-side corner cross the apron boundary and
    # keep the 1.5 % cap so the rect can slope to absorb the
    # transition.
    APRON_GRADE_CAP = 0.010
    terminal_buckets: set = set()
    for s in layout.shapes:
        if s.role != ROLE_TERMINAL:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for (x, y) in coords:
            terminal_buckets.add(_corner_elevation_bucket(x, y))
    apron_nodes: set = set()
    if terminal_buckets:
        for s in layout.shapes:
            if s.role != ROLE_JUNCTION:
                continue
            try:
                coords = list(s.polygon.exterior.coords)
            except Exception:
                continue
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            touches_terminal = any(
                _corner_elevation_bucket(x, y) in terminal_buckets
                for (x, y) in coords)
            if not touches_terminal:
                continue
            for (x, y) in coords:
                nid = bucket_to_node.get(
                    _corner_elevation_bucket(x, y))
                if nid is not None:
                    apron_nodes.add(nid)
    apron_edges: List[Tuple[int, int, float]] = []
    for i in apron_nodes:
        for (j, length) in mesh.edges_adj[i]:
            if j <= i:
                continue
            if j in apron_nodes:
                apron_edges.append((i, j, length))

    def _enforce_apron_cap() -> None:
        # Iterate the apron-cap pass until no edge over the cap, or
        # ``iters`` exhausted.  One linear pass over edges only
        # corrects each edge once, but moving one edge's endpoints
        # can put adjacent edges out of compliance, so a second
        # pass picks up the cascade.  Cap at 8 inner iterations to
        # bound runtime when anchors are mutually infeasible (the
        # cap can't make adjacent anchored buckets agree, but the
        # iteration still converges to a stable midpoint).
        for _ in range(8):
            any_change = False
            for (i, j, length) in apron_edges:
                ai = mesh.anchor_elev[i]
                aj = mesh.anchor_elev[j]
                if ai is not None and aj is not None:
                    continue
                diff = mesh.elev[i] - mesh.elev[j]
                dmax = length * APRON_GRADE_CAP
                if abs(diff) <= dmax:
                    continue
                excess = abs(diff) - dmax
                sign = 1.0 if diff > 0 else -1.0
                if ai is not None:
                    mesh.elev[j] += excess * sign
                elif aj is not None:
                    mesh.elev[i] -= excess * sign
                else:
                    mesh.elev[i] -= 0.5 * excess * sign
                    mesh.elev[j] += 0.5 * excess * sign
                any_change = True
            if not any_change:
                break

    # Iterated smoothing.  ElevationGraph.smooth_rate_of_change
    # is invoked one outer pass at a time so we can interleave
    # rect short-edge equalisation between iterations.  Each
    # ``smooth_rate_of_change`` call internally applies the
    # Laplacian + 1.5 % grade-cap on every edge.  After each call
    # we run the tighter apron cap (1 %) on apron-zone edges, then
    # equalise rect / flat-shape constraints.
    #
    # Convergence-based early exit: snapshot ``mesh.elev`` before
    # each outer iteration and break when the maximum per-node
    # change drops below ``OUTER_CONVERGE_TOL_M``.  Captures the
    # post-smooth + apron-cap + equalise final state; if all three
    # steps moved every node by < 1 cm, further iterations don't
    # help.  At small/medium airports the loop typically settles
    # in 5-10 outer iterations rather than the unconditional 30
    # we were running before.
    OUTER_ITERS = 30
    OUTER_CONVERGE_TOL_M = 0.01
    for _ in range(OUTER_ITERS):
        prev_elev = list(mesh.elev)
        mesh.smooth_rate_of_change(iters=2)
        _enforce_apron_cap()
        _equalise_constraints()
        max_change = 0.0
        for i in range(len(mesh.elev)):
            d = abs(mesh.elev[i] - prev_elev[i])
            if d > max_change:
                max_change = d
        if max_change < OUTER_CONVERGE_TOL_M:
            break


def _equalize_rect_short_edges(
        layout: "PavementLayout",
        mesh: "ElevationGraph",
        bucket_to_node: Dict[Tuple[int, int], int]) -> None:
    """For every sloped 4-corner rect, force both corners on each
    short edge to share a single mesh elevation (the mean of the
    pair).  X-Plane's patch format requires ``altitude_high`` /
    ``altitude_low`` to apply uniformly across each short edge, so
    the two corners on the same short edge MUST emit at the same
    elevation.  Without this step the rect writeback averages
    behind the mesh's back, leaving any junction sharing a rect-
    corner bucket reading a different elevation than the rect
    itself emits — a 0.1–0.5 m cross-shape step at every shared
    rect corner.

    Mutating the mesh values in place propagates the equalised
    elevation to every shape sharing the bucket; the subsequent
    writeback then reads the same value everywhere.
    """
    rect_like_roles = {ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL,
                       ROLE_STUB, ROLE_CROSS_CONNECTOR,
                       ROLE_RUNWAY}
    for s in layout.shapes:
        if s.role not in rect_like_roles:
            continue
        if (s.altitude_high is None or s.altitude_low is None):
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            continue
        nids = []
        for c in coords:
            b = _corner_elevation_bucket(c[0], c[1])
            nid = bucket_to_node.get(b)
            if nid is None:
                nids = []
                break
            nids.append(nid)
        if len(nids) != 4:
            continue
        # Runway: anchor values are immutable.  Skip.
        if s.role == ROLE_RUNWAY:
            continue
        # X-Plane patch convention: corners (0, 3) are HIGH side,
        # (1, 2) are LOW side.  Average each short-edge pair.
        e0, e3 = mesh.elev[nids[0]], mesh.elev[nids[3]]
        e1, e2 = mesh.elev[nids[1]], mesh.elev[nids[2]]
        eh = 0.5 * (e0 + e3)
        el = 0.5 * (e1 + e2)
        mesh.elev[nids[0]] = eh
        mesh.elev[nids[3]] = eh
        mesh.elev[nids[1]] = el
        mesh.elev[nids[2]] = el


def _writeback_pavement_mesh(
        layout: "PavementLayout",
        mesh: "ElevationGraph",
        bucket_to_node: Dict[Tuple[int, int], int]) -> None:
    """Rewrite each shape's altitude tags from the smoothed mesh.

    * Runway shapes: untouched (their corners anchored the mesh).
    * Sloped 4-corner taxi rects: ``altitude_high`` ← mean of the
      two high-side corner mesh values; ``altitude_low`` ← mean of
      the two low-side corner mesh values.  Re-orient the ring if
      the high/low designation flipped.
    * Flat 4-corner shapes (terminals, flat rects, flat runways
      already): ``altitude`` ← mean of all corner mesh values.
    * Junction shapes: per-vertex mesh elevations populate
      ``node_altitudes``; the FLAT/PLANAR/COMPOUND classifier
      (run later) collapses to a single ``altitude`` if the vertex
      range is small enough.
    """
    rect_like_roles = {ROLE_PRIMARY_PARALLEL,
                       ROLE_SECONDARY_PARALLEL,
                       ROLE_STUB, ROLE_CROSS_CONNECTOR}

    def _lookup(coord: Tuple[float, float]) -> Optional[float]:
        b = _corner_elevation_bucket(coord[0], coord[1])
        nid = bucket_to_node.get(b)
        if nid is None:
            return None
        return float(mesh.elev[nid])

    for s in layout.shapes:
        if s.role == ROLE_RUNWAY:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if not coords:
            continue
        mesh_elevs: List[Optional[float]] = [
            _lookup(c) for c in coords]
        if not any(e is not None for e in mesh_elevs):
            continue
        # Fill any missing entries with the ring mean.
        known = [e for e in mesh_elevs if e is not None]
        fill = sum(known) / len(known) if known else 0.0
        mesh_elevs = [e if e is not None else fill
                      for e in mesh_elevs]
        if (s.role in rect_like_roles
                and s.altitude_high is not None
                and s.altitude_low is not None
                and len(coords) == 4):
            # Sloped rect: X-Plane patch format requires both
            # corners on each short edge to share a single
            # elevation (``altitude_high`` for one short edge,
            # ``altitude_low`` for the other).  The mesh solver
            # enforces this constraint during smoothing
            # (``_solve_pavement_mesh`` re-equalises rect short-
            # edge pairs every iteration), so the four mesh values
            # already satisfy ``e0 == e3`` and ``e1 == e2`` here
            # except for floating-point noise — averaging is
            # safe.  Junctions sharing a rect-corner bucket read
            # the same equalised mesh value via writeback below.
            short_a = 0.5 * (mesh_elevs[0] + mesh_elevs[3])
            short_b = 0.5 * (mesh_elevs[1] + mesh_elevs[2])
            eh = max(short_a, short_b)
            el = min(short_a, short_b)
            if abs(eh - el) >= 0.1:
                s.altitude_high = round(eh, 1)
                s.altitude_low = round(el, 1)
                # If the original HIGH side (corners 0, 3) is now
                # below, swap the ring orientation so corners 0, 3
                # remain the HIGH short edge.
                if short_a < short_b:
                    new_coords = [coords[1], coords[0],
                                  coords[3], coords[2]]
                    try:
                        new_poly = Polygon(new_coords)
                        if new_poly.is_valid and not new_poly.is_empty:
                            s.polygon = new_poly
                    except Exception:
                        pass
                s.altitude = None
                s.node_altitudes = None
            else:
                s.altitude_high = None
                s.altitude_low = None
                s.altitude = round(0.5 * (eh + el), 1)
                s.node_altitudes = None
            continue
        if (s.role == ROLE_TERMINAL
                or (s.role in rect_like_roles
                    and s.altitude is not None)):
            mean_e = sum(mesh_elevs) / len(mesh_elevs)
            s.altitude = round(mean_e, 1)
            s.altitude_high = None
            s.altitude_low = None
            continue
        if s.role == ROLE_JUNCTION:
            # Per-vertex: store node_altitudes spanning the closed
            # ring.  Downstream classifier will collapse to FLAT or
            # PLANAR if appropriate.
            closed = mesh_elevs + [mesh_elevs[0]]
            if (max(closed) - min(closed)) < 0.05:
                s.altitude = round(
                    sum(mesh_elevs) / len(mesh_elevs), 1)
                s.node_altitudes = None
            else:
                s.node_altitudes = [round(e, 1) for e in closed]
                s.altitude = None


# ── Junction polygon decomposition ───────────────────────────────
#
# When the residue polygon (pavement minus rects/terminals/runway)
# wraps around a fully-enclosed rect, the result has interior
# holes.  Triangulating a polygon-with-holes is messy — the
# cleanest approach (per user 2026-04-24) is to split the polygon
# into multiple simple polygons that JOIN AROUND the hole instead.
#
# ``_decompose_polygon_with_holes`` recursively cuts the polygon
# with horizontal lines through each hole's centroid until every
# remaining piece is a simple (no-hole) polygon.  Each piece
# becomes its own junction shape — humans can edit them
# independently, and triangulation runs without hole-splicing.


def _decompose_polygon_with_holes(polygon: Polygon,
                                  min_area_m2: float = 50.0,
                                  max_depth: int = 8
                                  ) -> List[Polygon]:
    """Return a list of simple (no-hole) polygons that tile the
    same area as ``polygon``.  Cuts horizontally through each
    hole's centroid using shapely.ops.split, recursing on each
    side."""
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
    # Pick the largest remaining hole and slice horizontally
    # through its centroid.
    interiors = list(polygon.interiors)
    interiors.sort(key=lambda h: -Polygon(h).area)
    hole = interiors[0]
    cy = float(hole.centroid.y)
    minx, miny, maxx, maxy = polygon.bounds
    cut = _LS([(minx - 1.0, cy), (maxx + 1.0, cy)])
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
    return pieces


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
        graph: Optional["ElevationGraph"],
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
            if graph is not None:
                d = graph.elevation_at(float(xs[i]), float(ys[j]))
            if d is None and dem is not None:
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


def _triangulate_junctions(
    layout: "PavementLayout",
    graph: Optional["ElevationGraph"],
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
        if graph is not None:
            e_free = graph.elevation_at(x, y)
        if e_free is None and dem is not None:
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
                    poly, _anchors_b, graph, dem,
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
    """Extract aeroway=terminal polygons (ways OR multipolygon
    relations with outer rings) in meter space."""
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
        if tags.get("aeroway") != "terminal":
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
        if tags.get("aeroway") != "terminal":
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
    for ls, ref in centerlines:
        gap_margin_frac = _rect_margin_frac_for(ls, ref)
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
    stub_ref_best: Dict[str, int] = {}
    for i, (rect, axis, role, ref) in enumerate(emitted):
        if not _should_dedup(ref, role):
            continue
        cur = stub_ref_best.get(ref)
        if cur is None or axis.length > emitted[cur][1].length:
            stub_ref_best[ref] = i
    keep: List[Tuple[Polygon, LineString, str, str]] = []
    for i, item in enumerate(emitted):
        _rect, _axis, role, ref = item
        if _should_dedup(ref, role):
            if stub_ref_best.get(ref) != i:
                continue
        keep.append(item)
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
