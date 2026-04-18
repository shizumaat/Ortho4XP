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

    # ── Load OSM centerlines + relations ─────────────────────────
    nodes, ways, relations = _load_osm_airports(
        xplane_root, icao, anchor[0], anchor[1])
    osm_centerlines = _extract_osm_taxi_centerlines(nodes, ways, to_m)

    # ── Terminals from OSM (way or relation aeroway=terminal) ────
    terminal_polys = _extract_osm_terminals(nodes, ways, relations, to_m)
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
    junction_points = _find_junction_points(nodes, ways, to_m)

    # ── Build taxi rects from centerlines ────────────────────────
    taxi_rects = _build_taxi_rects(
        osm_centerlines, pav_union, layout.runway_union, rwy_centerlines)

    # Emit taxi rects (already trimmed to narrow-width portion).
    emitted_taxi_rects: List[Polygon] = []
    for rect, axis, role, ref in taxi_rects:
        emitted_taxi_rects.append(rect)
        layout.shapes.append(BuiltShape(
            polygon=rect, role=role, ref=ref, source_axis=axis))

    # ── Junction polys from rect corner vertices ───────────────
    # For each cluster where >= 2 rect ends meet, junction polygon
    # has one vertex per rect corner (user's rule).  These are
    # emitted FIRST; residue fills the rest.  No overlap.
    corner_junctions = _build_junction_polys_from_corners(
        junction_points, taxi_rects, pav_union,
        terminal_union=terminal_union)
    corner_junction_union = (unary_union(corner_junctions)
                             if corner_junctions else None)
    for jp in corner_junctions:
        layout.shapes.append(BuiltShape(polygon=jp, role=ROLE_JUNCTION))

    # ── Residue (fills areas not covered by rects/terminals/corner-jns)
    # Classify each connected component by area:
    #   - >= MIN_APRON_AREA  → apron (terminal apron / engine test)
    #   - MIN_JUNCTION .. MIN_APRON  → junction
    # Residue polygon vertices landing within SNAP_TO_CORNER_M of any
    # rect corner are snapped to that corner — junctions share
    # vertices exactly with adjacent rects.
    MIN_APRON_AREA_M2 = 25000.0
    MIN_JUNCTION_AREA_M2 = 80.0
    SNAP_TO_CORNER_M = 1.5
    taxi_rect_union = (unary_union(emitted_taxi_rects)
                       if emitted_taxi_rects else None)
    rect_corners = []
    for rect in emitted_taxi_rects:
        rc = list(rect.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        rect_corners.extend(rc)

    def _snap_to_rect_corners(poly: Polygon) -> Polygon:
        coords = list(poly.exterior.coords)
        out = []
        for (x, y) in coords:
            best = None
            best_d = SNAP_TO_CORNER_M
            for (rx, ry) in rect_corners:
                d = math.hypot(x - rx, y - ry)
                if d < best_d:
                    best = (rx, ry)
                    best_d = d
            out.append(best if best else (x, y))
        try:
            return Polygon(out).buffer(0)
        except Exception:
            return poly

    if pav_union is not None:
        residue = pav_union
        if taxi_rect_union is not None:
            residue = residue.difference(taxi_rect_union)
        if terminal_union is not None:
            residue = residue.difference(terminal_union)
        if corner_junction_union is not None:
            residue = residue.difference(corner_junction_union)
        parts = [residue] if residue.geom_type == "Polygon" else list(
            getattr(residue, "geoms", []))
        for part in parts:
            if part.geom_type != "Polygon":
                continue
            if part.area < MIN_JUNCTION_AREA_M2:
                continue
            role = ROLE_APRON if part.area >= MIN_APRON_AREA_M2 else ROLE_JUNCTION
            simp = part.simplify(1.0, preserve_topology=True)
            if simp.is_empty or simp.geom_type != "Polygon":
                simp = part
            snapped = _snap_to_rect_corners(simp)
            if snapped.is_empty or snapped.geom_type != "Polygon":
                snapped = simp
            layout.shapes.append(BuiltShape(polygon=snapped, role=role))

    return layout


# ──────────────────────────────────────────────────────────────────
# Centerline-based taxi rect builder
# ──────────────────────────────────────────────────────────────────

JUNCTION_CLUSTER_DIST_M = 80.0  # merge junction nodes within this distance
JUNCTION_RADIUS_SCALE = 1.5     # disc radius = local_half_width × this


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
) -> List[Tuple[float, float]]:
    """Identify junction POINTS (clusters of shared OSM nodes).

    Any OSM node referenced by ≥ 2 taxi ways is a candidate junction
    point.  Candidates within ``JUNCTION_CLUSTER_DIST_M`` of each
    other merge into one cluster; the cluster centroid is the
    junction point.
    """
    from collections import defaultdict
    in_ways: Dict[str, int] = defaultdict(int)
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        for n in nds:
            in_ways[n] += 1

    candidates: List[Tuple[float, float]] = []
    for nid, count in in_ways.items():
        if count < 2:
            continue
        if nid not in nodes:
            continue
        lat, lon = nodes[nid]
        candidates.append(to_m(lon, lat))

    # Cluster within JUNCTION_CLUSTER_DIST_M (greedy single-link)
    clusters: List[List[Tuple[float, float]]] = []
    for pt in candidates:
        placed = False
        for cl in clusters:
            # Check distance to any member
            if any(math.hypot(pt[0]-q[0], pt[1]-q[1]) <= JUNCTION_CLUSTER_DIST_M
                   for q in cl):
                cl.append(pt)
                placed = True
                break
        if not placed:
            clusters.append([pt])

    # Return centroids
    return [(sum(p[0] for p in cl)/len(cl),
             sum(p[1] for p in cl)/len(cl)) for cl in clusters]


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


def _extract_osm_taxi_centerlines(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
) -> List[Tuple[LineString, str]]:
    """Return list of (linestring_m, ref_tag) for aeroway=taxiway ways.

    We split each OSM way at its internal vertices so the emitter
    gets ONE SEGMENT PER STRAIGHT RUN between consecutive topology
    nodes.  This matches how targets are drawn: the user puts one
    rect per straight taxi segment and fills every bend/junction
    with a separate junction polygon.  Merging same-ref ways (e.g.
    A1–A6) would go the opposite direction and produce 400+ m
    rects that span multiple target segments.

    An OSM way that's already one straight segment yields one
    centerline.  An OSM way with N internal vertices yields N
    centerlines.
    """
    out = []
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
        total_len = sum(
            math.hypot(pts[i+1][0]-pts[i][0], pts[i+1][1]-pts[i][1])
            for i in range(len(pts)-1))
        # Unrefed OSM taxiways at SPJC are parking-position / gate
        # access paths that the target rolls into apron shapes.
        # Drop them when the airport has refs (SPJC).  At airports
        # without any refs (SPLP), we have to use unrefed centerlines
        # as the only geometric signal — keep them when long enough.
        # Caller-level heuristic below: if ANY refed way exists in
        # the batch, unrefed ways are suppressed.
        if not ref:
            # Defer filtering to post-pass so we know whether the
            # airport has any refs at all.
            pass

        # RDP-simplify with 1 m tolerance — almost lossless, keeps
        # every bend vertex.  The user's target convention is to
        # emit simple straight rect segments and fill bends / curves
        # / intersections with enlarged junction polygons.  So we
        # need MAXIMAL segmentation: any vertex where the centerline
        # changes direction breaks into a new segment.
        try:
            simp = LineString(pts).simplify(1.0, preserve_topology=False)
        except Exception:
            continue
        scoords = list(simp.coords)
        if len(scoords) < 2:
            continue
        # Split the simplified polyline at each remaining vertex.
        for i in range(len(scoords) - 1):
            try:
                seg = LineString([scoords[i], scoords[i+1]])
            except Exception:
                continue
            if seg.is_empty or seg.length < 15.0:
                continue
            out.append((seg, ref))

    # Post-pass: if the airport has any refs, drop unrefed
    # centerlines (they're parking/gate paths the target skips).
    # Keep unrefed centerlines at airports with NO refs (SPLP),
    # since they're the only geometric seed available.
    any_ref = any(r for _, r in out)
    if any_ref:
        out = [(ls, r) for ls, r in out if r]
    return out


def _build_taxi_rects(
    centerlines: List[Tuple[LineString, str]],
    pav_union: Optional[Polygon],
    rwy_union: Optional[Polygon],
    rwy_centerlines: List[LineString],
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
        try:
            clipped = axis.intersection(pav_non_rwy)
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

        # Sample half-widths along clipped axis to find natural width
        natural_hw = _natural_half_width(clipped, pav_non_rwy)
        if natural_hw < 3.5 or natural_hw > 40.0:
            continue

        # Trim endpoints inward where half-width exceeds widening threshold.
        # Aggressive widen_factor = 1.05 pulls the rect back to where the
        # apt.dat pavement is essentially at its natural width — the
        # trimmed region becomes junction territory.
        trimmed = _trim_to_narrow(clipped, pav_non_rwy, natural_hw,
                                 widen_factor=1.05)
        if trimmed is None or trimmed.length < 15.0:
            continue

        # Dedup (against trimmed axis so shorter dupes don't fall through)
        if emitted_union is not None and not emitted_union.is_empty:
            try:
                inside_len = trimmed.intersection(emitted_union).length
                if inside_len / trimmed.length > 0.7:
                    continue
            except Exception:
                pass

        width = 2.0 * natural_hw
        rect = _rect_from_axis_extended(trimmed, width, pav_non_rwy)
        if rect is None or rect.is_empty:
            continue

        role = _classify_role(trimmed, width, rwy_centerlines,
                               rwy_union, ref=ref)
        emitted.append((rect, trimmed, role, ref))
        emitted_union = (unary_union([emitted_union, rect])
                         if emitted_union is not None else rect)

    # Post-classify: secondary passes to fix roles based on neighbour
    # topology (stubs touch parallels, cross_connector touches 2 parallels)
    _refine_roles(emitted, rwy_centerlines)
    return emitted


def _natural_half_width(axis: LineString, pav: Polygon,
                        n_probes: int = 11) -> float:
    """Return the MEDIAN distance-to-boundary along the axis.

    This is the natural half-width of the strip — where the apt.dat
    pavement is narrowest along the centerline.  Endpoints that
    widen into junctions skew the mean, so median is robust.
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
    return dists[len(dists) // 2]


def _trim_to_narrow(axis: LineString, pav: Polygon, natural_hw: float,
                    widen_factor: float = 1.3) -> Optional[LineString]:
    """Trim the axis inward from each end until the half-width at
    the endpoint drops below ``widen_factor × natural_hw``.

    The ``widen_factor`` threshold is where we declare the pavement
    is "widening into a junction" — anything past that point is
    not part of the rect, it's junction territory.

    Trim step is 2 m.  We never trim more than 50 % of the axis
    length; if the trim would eat too much, fall back to the
    original axis (rect emits as-is, which may overlap a junction
    that later subtracts from it).
    """
    boundary = pav.boundary
    total_len = axis.length
    thresh = natural_hw * widen_factor
    step = 2.0
    max_trim = total_len * 0.45

    def _hw_at(t: float) -> float:
        pt = axis.interpolate(t)
        if not pav.contains(pt):
            return 0.0
        return pt.distance(boundary)

    trim_a = 0.0
    while trim_a < max_trim:
        if _hw_at(trim_a) <= thresh:
            break
        trim_a += step

    trim_b = total_len
    min_b = total_len - max_trim
    while trim_b > min_b:
        if _hw_at(trim_b) <= thresh:
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
                            pav: Polygon) -> Optional[Polygon]:
    """Build a rect around the axis at its first-to-last direction.

    The 4 corners are placed at axis endpoints ± perpendicular half-
    width, then each corner is snapped to the nearest apt.dat
    pavement boundary point within SNAP_RADIUS meters.  This matches
    the snapped target convention where every non-runway vertex sits
    on a pavement edge.

    If no pavement edge is within SNAP_RADIUS, the corner stays at
    its un-snapped position.
    """
    coords = list(axis.coords)
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
        (p1[0] + px * half, p1[1] + py * half),
        (p2[0] + px * half, p2[1] + py * half),
        (p2[0] - px * half, p2[1] - py * half),
        (p1[0] - px * half, p1[1] - py * half),
    ]
    snapped = _snap_corners_to_pavement(corners, pav)
    return Polygon(snapped)


SNAP_RADIUS_M = 5.0  # output rect corners snapped within this range


def _snap_corners_to_pavement(
    corners: List[Tuple[float, float]],
    pav: Polygon,
) -> List[Tuple[float, float]]:
    """Snap each corner to its nearest point on the pavement boundary,
    but only if the nearest point is within ``SNAP_RADIUS_M``."""
    boundary = pav.boundary
    snapped = []
    for (cx, cy) in corners:
        p = Point(cx, cy)
        near, _ = nearest_points(boundary, p)
        if p.distance(near) <= SNAP_RADIUS_M:
            snapped.append((near.x, near.y))
        else:
            snapped.append((cx, cy))
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
            # Q / R have multiple segments at varied bearings; the
            # combined role is cross_connector as long as ANY segment
            # is reasonably perpendicular.  Relax to Δ > 40°.
            db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
            if db is not None and db > 40.0:
                return ROLE_CROSS_CONNECTOR
            return ROLE_STUB
        if has_digit:
            # L1, A3, V5 etc. — always stubs in the SPJC target
            return ROLE_STUB
        # Plain-letter refs come in two flavours:
        #   a. Named parallel taxis (A, F, L, V) — primary_parallel.
        #   b. Apron-traversing connectors (B, C, D, E, G) — stubs
        #      even though their bearing may be near-parallel.
        # Distinguisher: (a) runs NEAR a runway centerline
        # (< PRIMARY_CORRIDOR_DIST_M = 500 m) AND is parallel
        # (Δ < 15°).  (b) sits deeper in the apron.
        db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
        if db is not None and db < 15.0 and rwy_centerlines:
            mid = axis.interpolate(0.5, normalized=True)
            dmin = min(mid.distance(r) for r in rwy_centerlines)
            if dmin < 500.0:
                if ref in ("M", "U"):
                    return ROLE_SECONDARY_PARALLEL
                return ROLE_PRIMARY_PARALLEL
        # Parallel-but-far or non-parallel plain-letter refs are
        # apron connectors = stubs.
        return ROLE_STUB

    # ── Angle-only fallback (SPLP, unnamed airports) ────────────
    db = _axis_to_nearest_rwy_db(axis, rwy_centerlines)
    if db is None:
        return ROLE_STUB
    if db < 15.0:
        # parallel
        return ROLE_PRIMARY_PARALLEL
    if db > 60.0 and axis.length >= 100.0:
        return ROLE_CROSS_CONNECTOR
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
    """Second pass: currently a no-op.

    The earlier heuristic (cross_connector must touch >= 2 parallels)
    produced false negatives at SPJC where Q / R rects are corner-
    snapped and don't share an exact metric boundary with the primary
    rects they visually connect.  We trust the ref-based classifier
    instead.
    """
    pass
