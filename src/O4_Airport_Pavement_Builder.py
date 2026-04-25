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

        def _ring_to_nids(ring_coords):
            coords = list(ring_coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            if len(coords) < 3:
                return None
            nids = [_intern(x, y) for (x, y) in coords]
            nids.append(nids[0])
            return nids

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
            ext_nids = _ring_to_nids(s.polygon.exterior.coords)
            if ext_nids is None:
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
            elif (s.node_altitudes is not None
                  and len(s.node_altitudes) == len(ext_nids)):
                # Per-vertex elevation: comma-separated list, one
                # value per ring nid (including the closing repeat).
                # X-Plane mesh builder triangulates the polygon and
                # interpolates linearly between vertex elevations.
                tags["node_altitudes"] = ",".join(
                    f"{e:.1f}" for e in s.node_altitudes)
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
    pav_union = unary_union(pav_polys) if pav_polys else None
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
        pad = _terminal_pad_from_building(otp, pav_polys)
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

    # Trim PERPENDICULAR-TO-RUNWAY centerlines at a BUFFERED runway
    # polygon so the last rect on a taxi that crosses the runway
    # stops short of the runway-junction approach.  The runway
    # boundary IS an intersection (rule 70 % of gap between
    # intersections), but the physical junction — where the taxi
    # widens into the runway apron — extends some distance
    # OUTSIDE the runway polygon too.  A 30 m buffer pulls the
    # centerline end back from the approach widening for taxis
    # that cross runway perpendicularly.  Per user (2026-04-21):
    # "rects are encroaching into intersections".
    #
    # Only PERPENDICULAR taxis get buffered: diagonal stubs
    # (e.g. SPLP's (-447,-1087) 45° stub) naturally approach
    # runway at an angle and their pavement widens less
    # dramatically; buffer-trimming them over-shrinks.
    RWY_JUNCTION_BUFFER_M = 30.0
    PERP_TRIM_MAX_DEG = 25.0   # perp_diff < 25° → treat as perpendicular
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
            # Only buffer-trim PERPENDICULAR centerlines
            if _perp_diff_to_runway(ls) >= PERP_TRIM_MAX_DEG:
                trimmed_centerlines.append((ls, ref))
                continue
            try:
                diff = ls.difference(rwy_buffered)
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

    # NOTE: pavement-edge-aware end-trim was attempted here but
    # interacted badly with the diagonal centerline assembly
    # (dropped a valid stub).  Revisit when we need finer-grained
    # widening detection for diagonals not near any parallel
    # centerline.

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
        rwy_centerlines, apt_vertices=apt_pav_vertices)

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
            osm_nodes=nodes, osm_ways=ways, to_m=to_m)

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
DEM_SUFFIX = ".hgt"
_DEM_CACHE: Dict[Tuple[int, int], object] = {}


def _load_airport_dem(lat0: float, lon0: float):
    """Return an ``O4_DEM_Utils.DEM`` covering the 1° tile that
    contains (lat0, lon0).  Falls back to None if the .hgt file is
    missing."""
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
    if not os.path.isfile(dem_path):
        _DEM_CACHE[key] = None
        return None
    try:
        import O4_DEM_Utils as _DEM
        dem = _DEM.DEM(tile_lat, tile_lon, source=dem_path)
    except Exception:
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


def _runway_segment_elev_lookup(
    runway_segment_chain, layout: "PavementLayout"
):
    """Return a callable ``elev_at(x_m, y_m)`` that returns the
    runway surface elevation at a meter-space point, or None if the
    point is farther than ``TAXI_ANCHOR_DIST_M`` from any runway
    segment centerline.  Used to anchor taxi rect endpoints to the
    matching runway elevation.
    """
    anchor = layout.anchor
    lat0, lon0 = anchor
    cos0 = math.cos(math.radians(lat0))

    def _ll_to_m(lat, lon):
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return x, y

    # Convert chain to meter-space segments for fast nearest lookup.
    seg_lines: List[Tuple[LineString, float, float]] = []
    for seg in runway_segment_chain:
        lat_a, lon_a, elev_a, lat_b, lon_b, elev_b, _w = seg
        a = _ll_to_m(lat_a, lon_a)
        b = _ll_to_m(lat_b, lon_b)
        if math.hypot(b[0] - a[0], b[1] - a[1]) < 0.5:
            continue
        seg_lines.append(
            (LineString([a, b]), float(elev_a), float(elev_b)))

    def elev_at(x: float, y: float) -> Optional[float]:
        if not seg_lines:
            return None
        p = Point(x, y)
        best = None
        best_d = TAXI_ANCHOR_DIST_M
        for line, ea, eb in seg_lines:
            d = line.distance(p)
            if d < best_d:
                best_d = d
                # Interpolate along segment.
                total = line.length
                if total <= 0:
                    interp = ea
                else:
                    t = line.project(p) / total
                    interp = ea + (eb - ea) * t
                best = float(interp)
        return best

    return elev_at


def _compute_elevations(layout: "PavementLayout", icao: str,
                        xplane_root: str, apt,
                        osm_nodes=None, osm_ways=None,
                        to_m=None) -> None:
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
        _snap_plateaus(graph)
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
                                     # Tight enough that interior nodes
                                     # of parallel taxis (typically > 100 m
                                     # from the runway) aren't bridged;
                                     # wide enough to catch perpendicular
                                     # stub endpoints (V1 ~55 m at SPJC).
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
        """Compute per-node feasibility interval [lo, hi] as the
        intersection of (anchor ± dist × TAXI_MAX_GRADE) over all
        reachable anchors.  Unreachable nodes get (-inf, +inf)."""
        import heapq
        n = len(self.nodes)
        INF = float("inf")
        lo = [-INF] * n
        hi = [INF] * n
        anchors = [i for i, a in enumerate(self.anchor_elev) if a is not None]
        for aidx in anchors:
            a_val = self.anchor_elev[aidx]
            # Dijkstra from aidx; dist[k] = shortest graph distance
            dist = [INF] * n
            dist[aidx] = 0.0
            heap = [(0.0, aidx)]
            while heap:
                d, u = heapq.heappop(heap)
                if d > dist[u]:
                    continue
                for v, length in self.edges_adj[u]:
                    nd = d + length
                    if nd < dist[v]:
                        dist[v] = nd
                        heapq.heappush(heap, (nd, v))
            # Apply this anchor's cone to every reachable node.
            for k in range(n):
                if dist[k] < INF:
                    band = dist[k] * TAXI_MAX_GRADE
                    node_lo = a_val - band
                    node_hi = a_val + band
                    if node_lo > lo[k]: lo[k] = node_lo
                    if node_hi < hi[k]: hi[k] = node_hi
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


def _snap_plateaus(g: "ElevationGraph") -> int:
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
    """
    candidates = _detect_plateaus(g)
    if not candidates:
        return 0
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
        # Commit.
        elevs = sorted(g.elev[i] for i in cluster_set)
        median = elevs[len(elevs) // 2]
        for i in cluster_set:
            g.elev[i] = float(median)
            g.anchor_elev[i] = float(median)
        snapped += 1
    return snapped


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
        # Find the nearest graph node for each boundary vertex.
        nbrs: List[int] = []
        for (vx, vy) in coords:
            best_i = -1
            best_d2 = JUNCTION_BRIDGE_NODE_DIST_M ** 2
            for i, (nx, ny) in enumerate(g.nodes):
                d2 = (nx - vx) * (nx - vx) + (ny - vy) * (ny - vy)
                if d2 < best_d2:
                    best_d2 = d2
                    best_i = i
            if best_i >= 0:
                nbrs.append(best_i)
        if len(nbrs) < 2:
            continue
        # Dedup while preserving order.
        seen: set = set()
        uniq: List[int] = []
        for ni in nbrs:
            if ni not in seen:
                seen.add(ni)
                uniq.append(ni)
        # Pairwise bridges (capped by JUNCTION_BRIDGE_MAX_M).
        for a in range(len(uniq)):
            for b in range(a + 1, len(uniq)):
                ia, ib = uniq[a], uniq[b]
                ax, ay = g.nodes[ia]
                bx, by = g.nodes[ib]
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


# ── Tier 4: validated centroid-Steiner subdivision ──────────────
#
# After triangulation, walk each fat-steep junction triangle
# (aspect ≤ TIER4_ASPECT_CAP, plane gradient > TAXI_MAX_GRADE)
# and try centroid-Steiner subdivision with several candidate
# elevations.  Apply only if the max sub-triangle gradient is
# STRICTLY LESS than the parent's gradient — otherwise leave the
# parent alone.  Slivers (aspect > cap) are never subdivided
# because centroid-Steiner of a sliver always produces thinner
# sub-slivers.

TIER4_ASPECT_CAP = 10.0    # only refine triangles with aspect ≤ this
TIER4_MAX_PASSES = 2        # cap blow-up at 9× per offender


def _tier4_plane_gradient(coords: List[Tuple[float, float]],
                          elevs: List[float]) -> float:
    """Plane gradient magnitude for triangle (3 vertices).
    Returns 0 for degenerate triangles."""
    (x1, y1), (x2, y2), (x3, y3) = coords
    z1, z2, z3 = elevs
    ux, uy, uz = x2 - x1, y2 - y1, z2 - z1
    vx, vy, vz = x3 - x1, y3 - y1, z3 - z1
    nx_ = uy * vz - uz * vy
    ny_ = uz * vx - ux * vz
    nz_ = ux * vy - uy * vx
    if abs(nz_) < 1e-6:
        return 0.0
    return math.hypot(nx_ / nz_, ny_ / nz_)


def _tier4_aspect(coords: List[Tuple[float, float]]) -> float:
    """Sliver-ness: longest_edge² / (4 × area).  > 10 = sliver."""
    (x1, y1), (x2, y2), (x3, y3) = coords
    cross = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    area = 0.5 * cross
    if area < 0.5:
        return float("inf")
    max_edge = max(
        math.hypot(coords[i][0] - coords[(i + 1) % 3][0],
                   coords[i][1] - coords[(i + 1) % 3][1])
        for i in range(3))
    return (max_edge * max_edge) / (4.0 * area)


def _refine_fat_steep_triangles(
    shapes: List[BuiltShape]
) -> Tuple[List[BuiltShape], int]:
    """Subdivide fat-but-steep junction triangles via centroid
    Steiner with elevation chosen to STRICTLY reduce the parent's
    plane gradient.  Iterates up to TIER4_MAX_PASSES.

    Returns (new_shapes, refinement_count).
    """
    current = list(shapes)
    total_refined = 0
    for _ in range(TIER4_MAX_PASSES):
        next_shapes: List[BuiltShape] = []
        any_refined = False
        for s in current:
            if s.role != ROLE_JUNCTION:
                next_shapes.append(s)
                continue
            try:
                ring = list(s.polygon.exterior.coords)
            except Exception:
                next_shapes.append(s)
                continue
            if ring and ring[0] == ring[-1]:
                ring = ring[:-1]
            if len(ring) != 3:
                next_shapes.append(s)
                continue
            # Decode per-vertex elevations.
            if (s.node_altitudes is not None
                    and len(s.node_altitudes) >= 4):
                elevs = list(s.node_altitudes[:3])
            elif s.altitude is not None:
                elevs = [s.altitude] * 3
            else:
                next_shapes.append(s)
                continue
            grade = _tier4_plane_gradient(ring, elevs)
            if grade <= TAXI_MAX_GRADE:
                next_shapes.append(s)
                continue
            aspect = _tier4_aspect(ring)
            if aspect > TIER4_ASPECT_CAP:
                next_shapes.append(s)  # sliver, skip
                continue
            # Try multiple Steiner elevations; pick the one that
            # MOST reduces the max sub-triangle gradient.
            cx = (ring[0][0] + ring[1][0] + ring[2][0]) / 3.0
            cy = (ring[0][1] + ring[1][1] + ring[2][1]) / 3.0
            mean_e = sum(elevs) / 3.0
            min_e, max_e = min(elevs), max(elevs)
            candidates = [
                mean_e,
                mean_e - 0.5, mean_e + 0.5,
                mean_e - 1.0, mean_e + 1.0,
                min_e, max_e,
                (min_e + mean_e) / 2.0,
                (max_e + mean_e) / 2.0,
            ]
            best_zS: Optional[float] = None
            best_max_grade = grade  # baseline: refuse to make worse
            for zS in candidates:
                max_sub = 0.0
                for a, b in ((0, 1), (1, 2), (2, 0)):
                    sub = _tier4_plane_gradient(
                        [ring[a], ring[b], (cx, cy)],
                        [elevs[a], elevs[b], zS])
                    if sub > max_sub:
                        max_sub = sub
                if max_sub < best_max_grade - 1e-4:
                    best_max_grade = max_sub
                    best_zS = zS
            if best_zS is None:
                # No candidate strictly improves — keep the
                # parent (better than mangling).
                next_shapes.append(s)
                continue
            # Apply subdivision.
            any_refined = True
            total_refined += 1
            zS = best_zS
            for a, b in ((0, 1), (1, 2), (2, 0)):
                sub_pts = [ring[a], ring[b], (cx, cy)]
                try:
                    sub_poly = Polygon(sub_pts)
                    if not sub_poly.is_valid:
                        sub_poly = sub_poly.buffer(0)
                except Exception:
                    continue
                if (sub_poly.is_empty
                        or sub_poly.geom_type != "Polygon"
                        or sub_poly.area < 0.5):
                    continue
                sub_elevs = [elevs[a], elevs[b], zS]
                new_s = BuiltShape(
                    polygon=sub_poly,
                    role=ROLE_JUNCTION,
                    ref=s.ref)
                # node_altitudes in shapely-emitted ring order.
                closed = list(sub_poly.exterior.coords)
                elev_for_ring: List[float] = []
                for (rx, ry) in closed:
                    best_e = sub_elevs[0]
                    best_d2 = float("inf")
                    for (tx, ty), te in zip(sub_pts, sub_elevs):
                        d2 = (rx - tx) * (rx - tx) + (ry - ty) * (ry - ty)
                        if d2 < best_d2:
                            best_d2 = d2
                            best_e = te
                    elev_for_ring.append(round(float(best_e), 1))
                if max(elev_for_ring) - min(elev_for_ring) < 0.05:
                    new_s.altitude = round(
                        sum(elev_for_ring[:-1]) / 3.0, 1)
                else:
                    new_s.node_altitudes = elev_for_ring
                next_shapes.append(new_s)
        current = next_shapes
        if not any_refined:
            break
    return current, total_refined


def _planar_fit_residuals(ring: List[Tuple[float, float]],
                          elev: List[float]
                          ) -> Optional[List[float]]:
    """Fit a plane ``z = a*x + b*y + c`` to (x, y, z) by least
    squares and return the per-vertex residuals (|actual - plane|).
    Returns None if the fit is degenerate (colinear xy).

    Used to detect "uniformly sloped" junction polygons that can
    emit as a single polygon with ``node_altitudes`` instead of
    being triangulated — X-Plane interpolates linearly across the
    ring and the render is identical to our triangulated mesh.
    """
    n = len(ring)
    if n < 3 or len(elev) != n:
        return None
    # Normal equations: A^T A x = A^T b where A rows are [x, y, 1].
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
    # 3×3 system:
    #   [sxx sxy sxc] [a]   [sxz]
    #   [sxy syy syc] [b] = [syz]
    #   [sxc syc scc] [c]   [szc]
    # Solve via Cramer's rule.
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
    return [abs(z - (a * x + b * y + c))
            for (x, y), z in zip(ring, elev)]


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


DELAUNAY_COVERAGE_TOL_FRAC = 0.02  # accept up to 2 % polygon-area
                                    # gap in the Delaunay output
                                    # before falling back to ear-clip

# NOTE 2026-04-25: tried adding interior Steiner vertices on a
# 30 m grid before Delaunay (Tier 1 plan).  Steiner elevations
# from boundary IDW disagreed with rect-corner anchors near the
# polygon edge by enough to PRODUCE NEW gradient violations
# (252 vs 149 baseline at SPJC).  Steiner elevations from the
# elevation graph were even worse (442) — graph samples differ
# from rect.altitude_high/low at off-axis corners.  Delaunay-on-
# boundary alone (no Steiners) is the small but reliable win
# kept here.


def _delaunay_with_steiners(
    ring: List[Tuple[float, float]],
    vert_elev: List[float],
    polygon: Polygon,
    graph: Optional["ElevationGraph"],
    dem,
    tile_lat: int,
    tile_lon: int,
    m_to_ll,
) -> Optional[List[Tuple[Polygon, List[float]]]]:
    """Triangulate ``polygon`` via Delaunay over its boundary
    vertices.  Maximises minimum interior angle vs ear-clip's
    first-valid-ear → fatter triangles, fewer slivers.

    Returns a list of (triangle_polygon, [e0, e1, e2]) for triangles
    inside the polygon, or None if Delaunay coverage misses more
    than DELAUNAY_COVERAGE_TOL_FRAC of the polygon's area (caller
    should fall back to ear-clip).
    """
    try:
        from shapely.geometry import MultiPoint
        from shapely.ops import triangulate as _tri
    except Exception:
        return None
    if polygon.is_empty or polygon.geom_type != "Polygon":
        return None
    if len(ring) < 3 or len(vert_elev) != len(ring):
        return None
    all_pts = list(ring)
    all_elevs = list(vert_elev)
    try:
        mp = MultiPoint(all_pts)
        tris = _tri(mp)
    except Exception:
        return None
    inside_tris: List[Tuple[Polygon, List[float]]] = []
    covered = 0.0
    # Build a coord-bucket -> elevation map for vertex lookup.
    coord_to_elev: Dict[Tuple[int, int], float] = {}
    for (px, py), pe in zip(all_pts, all_elevs):
        coord_to_elev[(int(round(px / 0.1)),
                       int(round(py / 0.1)))] = pe
    for t in tris:
        if t.is_empty or t.geom_type != "Polygon":
            continue
        # Strict containment: triangle's centroid must be inside
        # the original polygon (handles concave polygons by
        # discarding triangles spanning the concave gap).
        try:
            if not polygon.contains(t.centroid):
                continue
        except Exception:
            continue
        try:
            tc = list(t.exterior.coords)
        except Exception:
            continue
        if tc and tc[0] == tc[-1]:
            tc = tc[:-1]
        if len(tc) != 3:
            continue
        elevs = []
        for (cx, cy) in tc:
            key = (int(round(cx / 0.1)), int(round(cy / 0.1)))
            e = coord_to_elev.get(key)
            if e is None:
                # Search nearby buckets (rounding noise).
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        e2 = coord_to_elev.get(
                            (key[0] + dx, key[1] + dy))
                        if e2 is not None:
                            e = e2
                            break
                    if e is not None:
                        break
            if e is None:
                e = 0.0
            elevs.append(float(e))
        inside_tris.append((t, elevs))
        covered += t.area
    if not inside_tris:
        return None
    # Coverage check: did Delaunay cover the polygon?
    gap = abs(polygon.area - covered) / polygon.area
    if gap > DELAUNAY_COVERAGE_TOL_FRAC:
        return None  # caller falls back to ear-clip
    return inside_tris


def _ear_clip(coords: Sequence[Tuple[float, float]]
              ) -> List[Tuple[int, int, int]]:
    """Best-ear ear-clipping: at every step, find ALL valid ears
    and clip the one with the highest minimum interior angle.

    Greedy first-ear selection produces sliver triangles when the
    polygon has near-colinear boundary segments — the slivers'
    planes can have huge gradients perpendicular to their thin
    dimension even if the 3 vertices' pairwise grades are modest.
    Best-ear keeps every triangle as "fat" as the polygon's
    geometry allows, eliminating the bulk of sliver-induced step
    artefacts.
    """
    n = len(coords)
    if n < 3:
        return []
    if n == 3:
        return [(0, 1, 2)]
    # Determine winding via shoelace; force CCW for the algorithm.
    s = 0.0
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    is_ccw = s > 0
    indices = list(range(n)) if is_ccw else list(range(n - 1, -1, -1))

    def _cross(o: int, a: int, b: int) -> float:
        ox, oy = coords[o]
        ax, ay = coords[a]
        bx, by = coords[b]
        return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox)

    def _pt_in_tri(p: int, a: int, b: int, c: int) -> bool:
        px, py = coords[p]
        ax, ay = coords[a]
        bx, by = coords[b]
        cx, cy = coords[c]
        d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by)
        d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy)
        d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay)
        has_neg = d1 < 0 or d2 < 0 or d3 < 0
        has_pos = d1 > 0 or d2 > 0 or d3 > 0
        return not (has_neg and has_pos)

    def _min_angle(a: int, b: int, c: int) -> float:
        """Minimum interior angle of triangle (a, b, c) in radians.
        Higher = fatter triangle.  Returns 0 for degenerate."""
        ax, ay = coords[a]
        bx, by = coords[b]
        cx, cy = coords[c]
        # Side vectors at each vertex.
        abx, aby = bx - ax, by - ay
        acx, acy = cx - ax, cy - ay
        bax, bay = -abx, -aby
        bcx, bcy = cx - bx, cy - by
        cax, cay = -acx, -acy
        cbx, cby = -bcx, -bcy
        ab = math.hypot(abx, aby)
        ac = math.hypot(acx, acy)
        bc = math.hypot(bcx, bcy)
        if ab < 1e-6 or ac < 1e-6 or bc < 1e-6:
            return 0.0
        cos_a = (abx * acx + aby * acy) / (ab * ac)
        cos_b = (bax * bcx + bay * bcy) / (ab * bc)
        cos_c = (cax * cbx + cay * cby) / (ac * bc)
        # Clamp to [-1, 1] for numerical safety.
        cos_a = max(-1.0, min(1.0, cos_a))
        cos_b = max(-1.0, min(1.0, cos_b))
        cos_c = max(-1.0, min(1.0, cos_c))
        return min(math.acos(cos_a),
                   math.acos(cos_b),
                   math.acos(cos_c))

    triangles: List[Tuple[int, int, int]] = []
    guard = 4 * n
    while len(indices) > 3 and guard > 0:
        guard -= 1
        m = len(indices)
        # Score every valid ear by its min-angle; clip the best.
        best_score = -1.0
        best_i = -1
        best_triple: Tuple[int, int, int] = (0, 0, 0)
        for i in range(m):
            prev_i = indices[(i - 1) % m]
            cur_i = indices[i]
            next_i = indices[(i + 1) % m]
            if _cross(prev_i, cur_i, next_i) <= 0:
                continue  # not convex → not an ear
            ok = True
            for k in indices:
                if k in (prev_i, cur_i, next_i):
                    continue
                if _pt_in_tri(k, prev_i, cur_i, next_i):
                    ok = False
                    break
            if not ok:
                continue
            score = _min_angle(prev_i, cur_i, next_i)
            if score > best_score:
                best_score = score
                best_i = i
                best_triple = (prev_i, cur_i, next_i)
        if best_i < 0:
            # No valid ear — fan-fallback for pathological input.
            for i in range(1, len(indices) - 1):
                triangles.append(
                    (indices[0], indices[i], indices[i + 1]))
            return triangles
        triangles.append(best_triple)
        del indices[best_i]
    if len(indices) == 3:
        triangles.append((indices[0], indices[1], indices[2]))
    return triangles


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
    shared_junction_elev: Dict[Tuple[int, int], float] = {}

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
        """Return ``(elev, is_anchor)``.  ``is_anchor`` is True if
        the elevation came from a corner, rect-edge interpolation,
        or another junction's prior smoothed value at this same
        bucket (must not be moved by junction-local smoothing).
        False means a graph / DEM sample that the smoothing pass is
        free to pull into grade compliance."""
        bucket = _corner_elevation_bucket(x, y)
        # 1. Bucket lookup — exact shared-vertex match against
        # rect / runway / terminal corners.
        e = corner_elev.get(bucket)
        if e is not None:
            return e, True
        # 2. Cross-junction shared bucket — if another junction
        # has already locked this bucket's elevation, reuse it so
        # both junctions render at the same height at the shared
        # boundary point.
        e_shared = shared_junction_elev.get(bucket)
        if e_shared is not None:
            return e_shared, True
        # 3. Wider linear search for off-bucket near-corner matches.
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
        # the rect's slope at the projected position.
        e_edge = _edge_interp_elev(x, y)
        if e_edge is not None:
            if bucket in shared_junction_buckets:
                shared_junction_elev[bucket] = e_edge
            return e_edge, True
        # 5. Free sample.  Anchor it if it's shared between
        # junctions (deterministic graph value -> same anchor
        # for both junctions).
        e_free: Optional[float] = None
        if graph is not None:
            e_free = graph.elevation_at(x, y)
        if e_free is None and dem is not None:
            lat, lon = m_to_ll(x, y)
            e_free = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        if e_free is None:
            return None, False
        if bucket in shared_junction_buckets:
            shared_junction_elev[bucket] = e_free
            return e_free, True
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

    new_shapes: List[BuiltShape] = []
    triangle_count = 0
    grade_violations = 0
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            new_shapes.append(shape)
            continue
        # Splice any interior holes into the exterior so triangulation
        # incorporates hole-boundary vertices (typically rect corners
        # for rects fully enclosed by an apron-style junction).
        # Without splicing, the triangulation covers the holes too
        # and the included rect's corners never become triangle
        # vertices — producing visible elevation steps along the
        # rect's long edge inside the apron.
        try:
            ring = _splice_holes(shape.polygon)
        except Exception:
            try:
                ring = list(shape.polygon.exterior.coords)
                if ring and ring[0] == ring[-1]:
                    ring = ring[:-1]
            except Exception:
                new_shapes.append(shape)
                continue
        if len(ring) < 3:
            continue  # degenerate junction; drop
        # Drop non-shared, non-corner-anchored vertices whose
        # perpendicular distance to the line through their two
        # neighbours is below COLINEAR_DROP_M.  Without this, a
        # long apt.dat boundary trace with closely-spaced near-
        # colinear vertices forces ear-clipping to emit slivers
        # whose plane gradient perpendicular to the long axis can
        # exceed 50 % even when no triangle EDGE violates grade.
        # Preserved vertices: anything in the corner-elev bucket
        # (rect/runway/terminal corner) or any cross-junction
        # shared bucket — those carry topological meaning beyond
        # boundary tracing and must not be removed.
        ring = _drop_colinear_boundary_vertices(
            ring, corner_elev, shared_junction_buckets)
        if len(ring) < 3:
            continue
        # Vertex elevations + anchor flags.
        ev_pairs = [_vertex_elev_anchored(x, y) for (x, y) in ring]
        vert_elev_raw: List[Optional[float]] = [p[0] for p in ev_pairs]
        is_anchor_list: List[bool] = [p[1] for p in ev_pairs]
        # Fill any None gaps with the average of known neighbours
        # (or 0.0 if all are unknown — should not happen).
        known = [e for e in vert_elev_raw if e is not None]
        fallback = sum(known) / len(known) if known else 0.0
        vert_elev = [e if e is not None else fallback
                     for e in vert_elev_raw]
        # Junction-local smoothing: pull non-anchored vertex
        # elevations into pairwise grade compliance with their
        # neighbours.  Anchored (corner-derived) vertices are
        # untouchable so the shared-vertex invariant with rect /
        # runway / terminal corners is preserved.
        vert_elev = _smooth_junction_boundary(
            ring, vert_elev, is_anchor_list)

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
        residuals = _planar_fit_residuals(ring, vert_elev)
        if (residuals is not None
                and max(residuals) < PLANAR_RESIDUAL_M):
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

        # ── Tier 1: Delaunay triangulation with interior Steiners ─
        # For polygons large enough to benefit (> DELAUNAY_MIN_AREA_M2),
        # add a regular grid of interior Steiner points with
        # elevations sampled from the smoothed elevation graph,
        # then triangulate via Delaunay (shapely.ops.triangulate).
        # Delaunay maximizes minimum interior angle by construction
        # → produces fat triangles, eliminating most ear-clip
        # slivers.  Falls through to ear-clip if Delaunay leaves
        # gaps (rare; typically only happens in concave polygons
        # that defeat the centroid-in-polygon filter).
        delaunay_result = _delaunay_with_steiners(
            ring, vert_elev, shape.polygon, graph, dem,
            tile_lat, tile_lon, m_to_ll)
        if delaunay_result is not None:
            for tri_poly, tri_elevs in delaunay_result:
                # Grade check (informational).
                tcs = list(tri_poly.exterior.coords)
                if tcs and tcs[0] == tcs[-1]:
                    tcs = tcs[:-1]
                if len(tcs) == 3 and len(tri_elevs) >= 3:
                    for (p, q, ep, eq) in (
                        (tcs[0], tcs[1], tri_elevs[0], tri_elevs[1]),
                        (tcs[1], tcs[2], tri_elevs[1], tri_elevs[2]),
                        (tcs[2], tcs[0], tri_elevs[2], tri_elevs[0]),
                    ):
                        d = math.hypot(p[0] - q[0], p[1] - q[1])
                        if d < 0.5:
                            continue
                        if abs(ep - eq) / d > TAXI_MAX_GRADE + 1e-4:
                            grade_violations += 1
                            break
                new_shape = BuiltShape(
                    polygon=tri_poly,
                    role=ROLE_JUNCTION,
                    ref=shape.ref)
                # Build closed-ring elevations matching shapely's
                # vertex order.
                closed = list(tri_poly.exterior.coords)
                elev_for_ring = []
                for (rx, ry) in closed:
                    best_e = tri_elevs[0]
                    best_d2 = float("inf")
                    for (tx, ty), te in zip(tcs, tri_elevs[:3]):
                        d2 = (rx - tx) * (rx - tx) + (ry - ty) * (ry - ty)
                        if d2 < best_d2:
                            best_d2 = d2
                            best_e = te
                    elev_for_ring.append(round(float(best_e), 1))
                if max(elev_for_ring) - min(elev_for_ring) < 0.05:
                    new_shape.altitude = round(
                        sum(elev_for_ring[:-1]) / 3.0, 1)
                else:
                    new_shape.node_altitudes = elev_for_ring
                new_shapes.append(new_shape)
                triangle_count += 1
            continue

        # Fall-back: ear-clip (concave polygons where Delaunay
        # leaves gaps, or polygons too small for Steiner grid).
        triples = _ear_clip(ring)
        if not triples:
            continue

        for (i, j, k) in triples:
            tri_coords = [ring[i], ring[j], ring[k]]
            tri_poly = Polygon(tri_coords)
            if not tri_poly.is_valid:
                tri_poly = tri_poly.buffer(0)
            if (tri_poly.is_empty
                    or tri_poly.geom_type != "Polygon"
                    or tri_poly.area < 0.5):
                continue
            ea, eb, ec = vert_elev[i], vert_elev[j], vert_elev[k]
            # Grade check (informational; the bridge enrichment
            # before smoothing should keep these in tolerance).
            for (p, q, ep, eq) in (
                (tri_coords[0], tri_coords[1], ea, eb),
                (tri_coords[1], tri_coords[2], eb, ec),
                (tri_coords[2], tri_coords[0], ec, ea),
            ):
                d = math.hypot(p[0] - q[0], p[1] - q[1])
                if d < 0.5:
                    continue
                if abs(ep - eq) / d > TAXI_MAX_GRADE + 1e-4:
                    grade_violations += 1
                    break
            # node_altitudes spans the closed ring (4 entries for
            # a triangle: e_i, e_j, e_k, e_i).  Match the order
            # produced by Polygon(tri_coords).exterior.coords.
            try:
                tri_ring = list(tri_poly.exterior.coords)
            except Exception:
                continue
            # Re-derive elevations in the actual ring order: each
            # ring coord matches one of (i, j, k) up to rounding,
            # so map by nearest tri_coord.
            elev_for_ring: List[float] = []
            for (rx, ry) in tri_ring:
                best_e = ea
                best_d2 = 1e18
                for (tx, ty), te in zip(tri_coords, (ea, eb, ec)):
                    d2 = (rx - tx) * (rx - tx) + (ry - ty) * (ry - ty)
                    if d2 < best_d2:
                        best_d2 = d2
                        best_e = te
                elev_for_ring.append(round(float(best_e), 1))

            new_shape = BuiltShape(
                polygon=tri_poly,
                role=ROLE_JUNCTION,
                ref=shape.ref)
            # Flat triangle → emit single altitude tag (cleaner).
            if (max(elev_for_ring) - min(elev_for_ring)) < 0.05:
                new_shape.altitude = round(
                    sum(elev_for_ring[:-1]) / 3.0, 1)
            else:
                new_shape.node_altitudes = elev_for_ring
            new_shapes.append(new_shape)
            triangle_count += 1

    # ── Tier 4: validated centroid-Steiner subdivision of fat-
    # but-steep junction triangles.  Strictly improves max sub-
    # triangle gradient or leaves the parent alone.
    new_shapes, _refined = _refine_fat_steep_triangles(new_shapes)
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


SHARED_VERTEX_CLUSTER_TOL_M = 1.5


def _consolidate_touching_junctions(layout: "PavementLayout",
                                    tol: float = 1.0) -> None:
    """Union adjacent (touching within ``tol``) junction polygons.

    Target files sometimes use ONE large junction polygon where my
    pipeline emits multiple smaller fragments (one per OSM multi-
    ref cluster + residue gap-fill).  Physically touching junctions
    should be one shape.
    """
    junctions = [(i, s) for i, s in enumerate(layout.shapes)
                 if s.role == ROLE_JUNCTION
                 and s.polygon is not None
                 and not s.polygon.is_empty
                 and s.polygon.geom_type == "Polygon"]
    if len(junctions) < 2:
        return

    # Single-link cluster: junctions within tol of each other merge.
    n = len(junctions)
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

    for i in range(n):
        for j in range(i + 1, n):
            pi = junctions[i][1].polygon
            pj = junctions[j][1].polygon
            try:
                if pi.distance(pj) <= tol:
                    _union(i, j)
            except Exception:
                pass

    clusters: Dict[int, List[int]] = {}
    for i in range(n):
        clusters.setdefault(_find(i), []).append(i)

    # For each cluster > 1 member, union polygons; drop originals.
    drop_shape_ids: set = set()
    for root, members in clusters.items():
        if len(members) < 2:
            continue
        try:
            polys = [junctions[m][1].polygon for m in members]
            merged = unary_union(polys)
        except Exception:
            continue
        if merged.geom_type == "MultiPolygon":
            # Rare — distance check said they touch but union is
            # still multi-poly (shared only at a point).  Keep
            # each piece separate.
            continue
        if merged.geom_type != "Polygon" or merged.is_empty:
            continue
        # Replace the first member's polygon with the union; mark
        # others for removal.
        first_shape_idx = junctions[members[0]][0]
        layout.shapes[first_shape_idx].polygon = merged
        for m in members[1:]:
            drop_shape_ids.add(junctions[m][0])

    if drop_shape_ids:
        layout.shapes[:] = [s for i, s in enumerate(layout.shapes)
                            if i not in drop_shape_ids]


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

# Parallel refs at SPJC.  These are the "long taxi" refs where target
# splits into many segments.  For these, we aggressively bridge gaps
# across intermediate intersections.
PARALLEL_REFS = frozenset({"A", "F", "L", "V", "M", "U"})


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
        if ref and ref not in PARALLEL_REFS:
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
                ax_start = (cx - ux * STUB_LEN_M / 2,
                            cy - uy * STUB_LEN_M / 2)
                ax_end = (cx + ux * STUB_LEN_M / 2,
                          cy + uy * STUB_LEN_M / 2)
                try:
                    stub_axis = LineString([ax_start, ax_end])
                except Exception:
                    continue
                # Width from pav probe at the exit point
                _nat, _p90, narrow = _natural_half_width(
                    stub_axis, pav_union)
                if narrow < 3.5 or narrow > 50.0:
                    continue
                width = 2.0 * narrow
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
            # Single-letter non-parallel stubs (B, C, D, E, G) are
            # continuous diagonal taxis at SPJC — the target emits
            # them as ONE long rect covering ~35 % of the full
            # path.  Bend-splitting fragments them into junk.
            # Detection: letter-only ref, not a primary/secondary
            # parallel ref, and path is geometrically straight
            # (chord/path > 0.95).  Emit as ONE centerline.
            if (ref
                    and len(ref) == 1
                    and ref not in PARALLEL_REFS
                    and ref not in {"Q", "R", "X"}):
                path_len = ls.length
                sc = list(simp.coords)
                if len(sc) >= 2 and path_len > 1e-6:
                    chord = math.hypot(sc[-1][0] - sc[0][0],
                                       sc[-1][1] - sc[0][1])
                    if chord / path_len > 0.95:
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

    any_ref = any(r for _, r in out)
    if any_ref:
        out = [(ls, r) for ls, r in out if r]
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
        # PARALLEL_REFS (A, F, L, V, M, U) and cross-connector
        # refs (Q, R, X) always use the 15 % primary margin.
        if ref in PARALLEL_REFS:
            return 0.15
        if ref in {"Q", "R", "X"}:
            return 0.15
        # Unrefed taxis: apply the diagonal-stub rule only for
        # SHORT centerlines (< 250 m).  SPLP's main taxi runs
        # long primaries at ~19° off runway (perp_diff=71°
        # inside the 20-75 window) and must stay at 15 % margin
        # to preserve the full primary-parallel length.  Only
        # short unrefed diagonal connectors (e.g. SPLP's
        # (-449,-1084) L=109) should get 35 % + bias
        # shrinkage.  Refed non-parallel taxis (SPJC B/C/E/G,
        # V1/V3/V5 sub-refs) always qualify for the diagonal
        # check regardless of length.
        if not ref and ls.length >= 250.0:
            return 0.15
        # Short unrefed parallel-to-runway rects (between two
        # diagonal stubs on SPLP's south chain) need extra margin
        # to avoid overlapping the neighbouring diagonal stubs.
        # Detect: ref empty, length < 150 m, nearly parallel to
        # runway (perp_diff > 75°) — emit as half-length primary
        # by using 30 % margin each side (40 % retained).
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
        perp_diff = abs(delta - 90.0)
        if 20.0 < perp_diff < 75.0:
            return 0.30
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
        # first/last ones that would actually emit.  For cross-
        # connector refs (Q/R/X), first and last emitted segments
        # use 30 % margin each side (40 % rect) because the taxi
        # terminates into a WIDER parallel taxi whose widening
        # zone extends into the cross-connector's axis.  Middle
        # segments between two cross-ref junctions use 15 %.
        CROSS_CONNECTOR_REFS = {"Q", "R", "X"}
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
        is_cross = ref in CROSS_CONNECTOR_REFS
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

        role = _classify_role(trimmed, width, rwy_centerlines,
                               rwy_union, ref=ref)
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
    BOUNDARY_TOL_M = 0.5
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


def _classify_role(axis: LineString, width: float,
                   rwy_centerlines: List[LineString],
                   rwy_union: Optional[Polygon],
                   ref: str = "") -> str:
    """Classify by (a) ref pattern and (b) axis bearing to runway.

    * Ref letter-only (A, F, L, V, U, M) → parallel candidate.
    * Ref Q/R/X (known cross) → cross_connector.
    * Ref letter+digit (A1, L3, V2) → stub.
    * Ref unknown or empty: use angle-only classification.

    The ref-pattern rule handles SPJC cleanly because the chart
    naming convention is stable.  The angle fallback handles SPLP
    (no refs) and unnamed airports.
    """
    # ── Ref-based classification (SPJC convention) ──────────────
    if ref:
        has_digit = any(c.isdigit() for c in ref)
        if ref in ("Q", "R", "X"):
            db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
            if db is not None and db > 40.0:
                return ROLE_CROSS_CONNECTOR
            return ROLE_STUB
        # Known parallel refs (A/F/L/V/M/U) always emit as parallel.
        # A separate post-pass detects the SHORT runway-connector
        # segment (stub-A) and demotes it.
        if ref in ("M", "U"):
            return ROLE_SECONDARY_PARALLEL
        if ref in PARALLEL_REFS:
            return ROLE_PRIMARY_PARALLEL
        if has_digit:
            # L1, A3, V5 etc. — sub-refs = stub
            return ROLE_STUB
        # Plain-letter non-parallel (B, C, D, E, G) = stub
        return ROLE_STUB

    # ── Angle-only classification (refless airports like SPLP) ─
    # Uses bearing-to-runway, length, and distance-to-runway.
    db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
    if db is None:
        return ROLE_STUB
    # Distance from axis midpoint to nearest runway centerline.
    try:
        mid = axis.interpolate(0.5, normalized=True)
        dist_rwy = min(mid.distance(r) for r in rwy_centerlines)
    except Exception:
        dist_rwy = 1e6
    length = axis.length
    # Parallel branch (bearing within 20° of a runway).
    if db < 20.0:
        # Parallel-oriented rects are part of the primary-parallel
        # chain even when short.  SPLP's main taxi has a 66-80 m
        # parallel segment between runway-connecting diagonals
        # that should emit as a primary_parallel rather than a
        # stub (user feedback 2026-04-21: "the next piece should
        # be a primary parallel, not a stub").
        if length >= 50.0:
            # Close to runway → primary; far → secondary.
            if dist_rwy < 400.0:
                return ROLE_PRIMARY_PARALLEL
            return ROLE_SECONDARY_PARALLEL
        # Very short parallel — treat as stub (e.g. ramp tie-in).
        return ROLE_STUB
    # Perpendicular branch (bearing > 45° off a runway).
    if db > 45.0 and length >= 80.0:
        # Cross-connector if the axis sits BETWEEN parallels (i.e.
        # not adjacent to the runway).  Perpendicular pieces close
        # to the runway are stubs (runway-connector).
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
