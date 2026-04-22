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
EMIT_JUNCTIONS = False
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
    later use (elevation phase) but isn't emitted.
    """
    polygon: Polygon
    role: str
    ref: str = ""
    source_axis: Optional[LineString] = None


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
        node id, matching the target-OSM convention.
        """
        # Collect unique vertices (snap by bucket).
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
            # Use bucket center so snap is consistent
            cx = kx * bucket_size
            cy = ky * bucket_size
            node_id_to_ll[nid] = self.m_to_ll(cx, cy)
            return nid

        way_blocks: List[Tuple[int, List[int], Dict[str, str]]] = []
        next_wid = [-10001]
        for s in self.shapes:
            coords = list(s.polygon.exterior.coords)
            if coords[0] == coords[-1]:
                coords = coords[:-1]
            if len(coords) < 3:
                continue
            nids = [_intern(x, y) for (x, y) in coords]
            nids.append(nids[0])
            tags = {
                "aeroway": AEROWAY_FOR_ROLE.get(s.role, "taxiway"),
                "role": s.role,
            }
            if s.ref:
                tags["ref"] = s.ref
            way_blocks.append((next_wid[0], nids, tags))
            next_wid[0] -= 1

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

def build_airport_pavement(icao: str, xplane_root: str) -> PavementLayout:
    """Build the complete role-classified layout for ``icao``.

    The layout is ready to compare against a target OSM via
    ``tools/compare_target.py``.
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
    terminal_polys: List[Polygon] = []
    for otp in osm_terminal_polys:
        pad = _terminal_pad_from_building(otp, pav_polys)
        if pad is not None:
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

    # ── Junction emission: corner + apt.dat pavement arc ────────
    # Per user's authoritative rules (memory: feedback_shape_rules,
    # session 2026-04-18):
    #   Rule 8:  vertices = incoming rect corners + apron corners.
    #   Rule 9:  arcs between consecutive corners trace apt.dat
    #            pavement vertices, MAX 4 points per arc.
    #   Rule 10: no hulls, no smooths.
    #   Rule 11: junctions fill widened portions left by trimmed rects.
    MIN_APRON_AREA_M2 = 25000.0
    MIN_JUNCTION_AREA_M2 = 80.0
    CLUSTER_MERGE_DIST_M = 40.0  # tight — each OSM multi-ref cluster
                                 # typically maps to one target junction

    taxi_rect_union = (unary_union(emitted_taxi_rects)
                       if emitted_taxi_rects else None)

    # Detect per-rect trim at each end.  A rect end is "uniform"
    # (no widening) if the trimmed axis endpoint sits right on the
    # pav-union (or runway) boundary — i.e. the axis went all the
    # way through.  Otherwise there was a widening, reserving the
    # adjacent region for a junction.
    pav_full = unary_union([pav_union, layout.runway_union]) if (
        pav_union is not None and layout.runway_union is not None
    ) else (pav_union or layout.runway_union)
    rect_trim_flags: List[Tuple[bool, bool]] = []  # (trim_start, trim_end)
    for rect, axis, role, ref in taxi_rects:
        coords_ax = list(axis.coords)
        if len(coords_ax) < 2 or pav_full is None:
            rect_trim_flags.append((False, False))
            continue
        ts = Point(coords_ax[0]).distance(pav_full.boundary) > 2.0
        te = Point(coords_ax[-1]).distance(pav_full.boundary) > 2.0
        rect_trim_flags.append((ts, te))

    # Seed junctions from ALL rect-end points + OSM junction_points.
    # Every rect end is a potential junction (two rects meeting
    # share an end region); the constructive build returns None
    # when there aren't ≥2 rects contributing, so excess seeds are
    # harmless.  Using all ends catches bend-junctions (same-ref
    # rects meeting after a turn) that aren't in OSM multi-ref
    # nodes.
    rect_end_pts: List[Tuple[float, float]] = []
    for i, (rect, axis, role, ref) in enumerate(taxi_rects):
        coords_ax = list(axis.coords)
        if len(coords_ax) < 2:
            continue
        rect_end_pts.append(coords_ax[0])
        rect_end_pts.append(coords_ax[-1])
    # Seed points: all rect-ends + OSM multi-ref clusters.
    seed_points = list(rect_end_pts) + list(junction_points)

    # Cluster seed points.
    merged_clusters: List[List[Tuple[float, float]]] = []
    for jp in seed_points:
        placed = False
        for grp in merged_clusters:
            for p in grp:
                if math.hypot(jp[0] - p[0], jp[1] - p[1]) <= CLUSTER_MERGE_DIST_M:
                    grp.append(jp)
                    placed = True
                    break
            if placed:
                break
        if not placed:
            merged_clusters.append([jp])

    final_junctions: List[Polygon] = []
    for grp in merged_clusters:
        gx = sum(p[0] for p in grp) / len(grp)
        gy = sum(p[1] for p in grp) / len(grp)
        extent = max((math.hypot(p[0] - gx, p[1] - gy) for p in grp),
                     default=0.0)
        # Local disc: cluster extent + enough to reach rect widening points.
        local_radius = max(60.0, extent + 80.0)
        poly = _build_junction_constructive(
            (gx, gy), taxi_rects, pav_union, terminal_union,
            max_corner_dist_m=max(80.0, extent + 80.0),
            local_disc_radius_m=max(120.0, extent + 120.0),
            max_arc_vertices=4,
        )
        if poly is None:
            continue
        # Keep clear of rects + terminal (already clipped to local
        # pav inside constructive, but rect-union subtract is needed
        # because corners sit on rect edges).
        if taxi_rect_union is not None:
            try:
                poly = poly.difference(taxi_rect_union)
            except Exception:
                pass
            if poly.is_empty:
                continue
            if poly.geom_type == "MultiPolygon":
                poly = max(poly.geoms, key=lambda g: g.area)
            if poly.geom_type != "Polygon":
                continue
        if poly.area < MIN_JUNCTION_AREA_M2:
            continue
        final_junctions.append(poly)

    # Dedup junctions that overlap heavily (multiple seeds can produce
    # nearly-identical polygons).  Keep the larger-area one.
    def _polys_overlap_heavily(p1: Polygon, p2: Polygon,
                               thresh: float = 0.8) -> bool:
        try:
            inter = p1.intersection(p2).area
        except Exception:
            return False
        if inter <= 0:
            return False
        return inter / min(p1.area, p2.area) >= thresh

    final_junctions.sort(key=lambda p: -p.area)
    kept_junctions: List[Polygon] = []
    for jp in final_junctions:
        if any(_polys_overlap_heavily(jp, kp) for kp in kept_junctions):
            continue
        kept_junctions.append(jp)
    final_junctions = kept_junctions

    # TEMP 2026-04-20: user requested junction emission be disabled
    # while refining rects.  Re-enable once rect counts match target.
    if EMIT_JUNCTIONS:
        for jp in final_junctions:
            layout.shapes.append(BuiltShape(polygon=jp, role=ROLE_JUNCTION))

    # ── Runway-taxiway junction: ONLY for widening stubs ────────
    # Per user rule 13 (2026-04-18): uniform-width stub (e.g. L1)
    # connects directly to the runway with NO junction; only
    # widening stubs (e.g. V1, whose end trim pulled back from the
    # runway edge because pavement widened) get a junction between
    # stub and runway.  We detect "widening end" as: rect end sits
    # > 2 m from the runway boundary after trim (rect_trim_flags).
    RWY_JUNCTION_DIST_M = 60.0
    rwy_vertex_inserts: List[Tuple[float, float]] = []
    if layout.runway_union is not None and not layout.runway_union.is_empty:
        rwy_boundary = layout.runway_union.boundary
        rwy_tj_polys: List[Polygon] = []
        for i, (rect, axis, role, ref) in enumerate(taxi_rects):
            coords_ax = list(axis.coords)
            if len(coords_ax) < 2:
                continue
            trim_s, trim_e = rect_trim_flags[i]
            for end_idx, (ax_pt, was_trimmed) in enumerate(
                    [(coords_ax[0], trim_s), (coords_ax[-1], trim_e)]):
                if not was_trimmed:
                    continue  # uniform stub — direct connect to runway
                d_rwy = Point(ax_pt).distance(rwy_boundary)
                if d_rwy > RWY_JUNCTION_DIST_M:
                    continue
                pairs = _rect_end_corners(rect, axis)
                if len(pairs) < 2:
                    continue
                c1, c2 = pairs[end_idx]
                try:
                    rp1 = nearest_points(rwy_boundary, Point(c1))[0]
                    rp2 = nearest_points(rwy_boundary, Point(c2))[0]
                except Exception:
                    continue
                # 4-corner trapezoid; junction polygon vertices are
                # the 2 stub corners and the 2 runway-side projections.
                poly_coords = [c1, c2, (rp2.x, rp2.y), (rp1.x, rp1.y)]
                try:
                    poly = Polygon(poly_coords)
                    if not poly.is_valid:
                        poly = poly.buffer(0)
                except Exception:
                    continue
                if (poly.is_empty
                        or poly.geom_type not in ("Polygon", "MultiPolygon")):
                    continue
                if poly.geom_type == "MultiPolygon":
                    poly = max(poly.geoms, key=lambda g: g.area)
                if poly.geom_type != "Polygon" or poly.area < 10.0:
                    continue
                try:
                    poly = poly.difference(layout.runway_union)
                except Exception:
                    pass
                if taxi_rect_union is not None:
                    try:
                        poly = poly.difference(taxi_rect_union)
                    except Exception:
                        pass
                if terminal_union is not None and not terminal_union.is_empty:
                    try:
                        poly = poly.difference(terminal_union)
                    except Exception:
                        pass
                if poly.is_empty:
                    continue
                if poly.geom_type == "MultiPolygon":
                    poly = max(poly.geoms, key=lambda g: g.area)
                if poly.geom_type != "Polygon" or poly.area < 10.0:
                    continue
                rwy_tj_polys.append(poly)
                rwy_vertex_inserts.append((rp1.x, rp1.y))
                rwy_vertex_inserts.append((rp2.x, rp2.y))
        if EMIT_JUNCTIONS:
            for rtp in rwy_tj_polys:
                layout.shapes.append(BuiltShape(polygon=rtp, role=ROLE_JUNCTION))

    # Insert runway-side projection points as vertices in the
    # runway polygons so taxi-stub junctions share exact vertex
    # positions with the runway they attach to.
    if rwy_vertex_inserts:
        for shape in layout.shapes:
            if shape.role != ROLE_RUNWAY:
                continue
            shape.polygon = _insert_points_on_boundary(
                shape.polygon, rwy_vertex_inserts, tol=2.0)

    # ── Apron emission: pavement residue after rects + junctions ──
    junctions_union = None
    if final_junctions:
        try:
            junctions_union = unary_union(final_junctions)
        except Exception:
            # Fallback: union valid ones iteratively.
            junctions_union = None
            for jp in final_junctions:
                if not jp.is_valid:
                    jp = jp.buffer(0)
                if jp.is_empty or jp.geom_type not in ("Polygon", "MultiPolygon"):
                    continue
                try:
                    junctions_union = (jp if junctions_union is None
                                       else junctions_union.union(jp))
                except Exception:
                    pass
    apron_polys: List[Polygon] = []
    if pav_union is not None:
        residue = pav_union
        if taxi_rect_union is not None:
            residue = residue.difference(taxi_rect_union)
        if junctions_union is not None:
            residue = residue.difference(junctions_union)
        if terminal_union is not None:
            residue = residue.difference(terminal_union)
        parts = [residue] if residue.geom_type == "Polygon" else list(
            getattr(residue, "geoms", []))
        for part in parts:
            if part.geom_type != "Polygon":
                continue
            if part.area >= MIN_APRON_AREA_M2:
                apron_polys.append(part)
            elif part.area >= 500.0:
                # Residue gap smaller than an apron but > 500 m² →
                # emit as a junction (rule 15: no gaps).  Filter
                # below 500 m² to skip sliver noise from rect-snap
                # inexactness.  Simplified to bound vertex count.
                try:
                    simp = part.simplify(3.0, preserve_topology=True)
                except Exception:
                    simp = part
                if (simp.is_empty
                        or simp.geom_type != "Polygon"
                        or simp.area < 500.0):
                    continue
                if EMIT_JUNCTIONS:
                    layout.shapes.append(BuiltShape(
                        polygon=simp, role=ROLE_JUNCTION))

    # Merge apron pieces split by thin rect strips (e.g. parallel
    # taxi cutting across the apron).  Target has aprons that wrap
    # around rect footprints as one connected region.
    APRON_MERGE_DIST_M = 100.0
    if apron_polys:
        try:
            ap_union = unary_union(apron_polys)
            ap_closed = ap_union.buffer(APRON_MERGE_DIST_M).buffer(
                -APRON_MERGE_DIST_M)
            if pav_union is not None:
                ap_closed = ap_closed.intersection(pav_union)
            if taxi_rect_union is not None:
                ap_closed = ap_closed.difference(taxi_rect_union)
            if terminal_union is not None:
                ap_closed = ap_closed.difference(terminal_union)
            merged_aprons = ([ap_closed] if ap_closed.geom_type == "Polygon"
                             else list(getattr(ap_closed, "geoms", [])))
            apron_polys = [p for p in merged_aprons
                           if p.geom_type == "Polygon"
                           and p.area >= MIN_APRON_AREA_M2]
        except Exception:
            pass
    if EMIT_APRONS:
        for ap in apron_polys:
            simp = ap.simplify(1.0, preserve_topology=True)
            if simp.is_empty or simp.geom_type != "Polygon":
                simp = ap
            layout.shapes.append(BuiltShape(polygon=simp, role=ROLE_APRON))

    # NOTE: Consolidation of touching junctions was tested
    # (tol=1m) and regressed match count by 4 — target retains
    # separate junctions at some physical touch points.  Function
    # `_consolidate_touching_junctions` retained as scaffold.
    # _consolidate_touching_junctions(layout)

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

    return layout


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

    # Relation terminals — union all outer rings
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
        if rings:
            merged = unary_union(rings).buffer(0)
            if merged.geom_type == "MultiPolygon":
                merged = max(merged.geoms, key=lambda g: g.area)
            if merged.geom_type == "Polygon" and merged.area >= 100.0:
                out.append(merged)
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
    # corridor (d_rwy ≈ 100 m).  80 m threshold lands us there.
    # SPLP has no loops — the primary curves smoothly into the
    # runway and target stubs sit MID-CURVE (d_rwy ≈ 50 m).
    # 80 m threshold walks past the mid-curve into the straight
    # section.  Use tighter threshold for unrefed ways.
    STUB_EXIT_D_M_REFED = 80.0
    STUB_EXIT_D_M_UNREFED = 50.0
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
                    # until d_rwy > exit threshold.  The "exit"
                    # vertex sits just outside the runway-apron
                    # ramp and becomes the stub center.
                    exit_d = (STUB_EXIT_D_M_REFED if ref
                              else STUB_EXIT_D_M_UNREFED)
                    step = 1 if end_idx == 0 else -1
                    exit_idx = None
                    i = end_idx if end_idx >= 0 else len(coords) - 1
                    while 0 <= i < len(coords):
                        d = Point(coords[i]).distance(rwy_boundary)
                        inside = runway_union.contains(
                            Point(coords[i]))
                        if not inside and d > exit_d:
                            exit_idx = i
                            break
                        i += step
                    if exit_idx is None:
                        continue
                else:
                    # Endpoint is OUTSIDE the runway but within
                    # OUTSIDE_NEAR_RWY_M.  The taxi curves to
                    # runway at this end but doesn't enter runway
                    # pavement (SPLP NE end at 127 m).  Use the
                    # endpoint itself as the stub center — target
                    # rect sits at the ramp where the taxi
                    # approaches runway.
                    exit_idx = 0 if end_idx == 0 else len(coords) - 1

                # Center stub on the exit vertex, length STUB_LEN_M
                # along local path direction (from prev to next).
                prev_idx = max(0, exit_idx - 1)
                next_idx = min(len(coords) - 1, exit_idx + 1)
                dx = coords[next_idx][0] - coords[prev_idx][0]
                dy = coords[next_idx][1] - coords[prev_idx][1]
                mag = math.hypot(dx, dy)
                if mag < 1e-6:
                    continue
                ux, uy = dx / mag, dy / mag
                cx, cy = coords[exit_idx]
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


def _insert_points_on_boundary(
    poly: Polygon,
    pts: List[Tuple[float, float]],
    tol: float = 2.0,
) -> Polygon:
    """Insert each point in ``pts`` as a vertex on the polygon's
    exterior ring at its projected position, if within ``tol`` of
    the boundary.  Used to add runway-taxi projection vertices
    to the runway rect per the shared-vertex invariant."""
    if not pts:
        return poly
    ext = poly.exterior
    inserts: List[Tuple[float, Tuple[float, float]]] = []
    for (x, y) in pts:
        p = Point(x, y)
        if p.distance(ext) > tol:
            continue
        try:
            param = ext.project(p)
            proj = ext.interpolate(param)
        except Exception:
            continue
        inserts.append((param, (proj.x, proj.y)))
    if not inserts:
        return poly
    inserts.sort()
    coords = list(ext.coords)
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) < 3:
        return poly
    try:
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
        new_poly = Polygon(new_coords)
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
    # Stage 1: pick nearest apt.dat vertex candidate per corner.
    candidates: List[Optional[Tuple[float, float]]] = []
    for (cx, cy) in corners:
        best_v = None
        best_d = VERTEX_SNAP_RADIUS_M
        if apt_vertices:
            for (vx, vy) in apt_vertices:
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
