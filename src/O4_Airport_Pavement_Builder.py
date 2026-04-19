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
    osm_centerlines = _extract_osm_taxi_centerlines(nodes, ways, to_m)


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

    # Per user rule (2026-04-18): "Implement splitting at cross ref
    # crossings" — split every centerline at each multi-ref
    # junction node along its path, so each straight section
    # between intersections emits as a single rect.
    # NOTE: widening-aware split (skipping OSM multi-ref nodes
    # where pav runs straight through) was tested with factor 1.05
    # and 1.15; both regressed match count by 14-17 at SPJC because
    # the target subdivides at some nodes that show only 1-4%
    # widening.  Reverting to unconditional split; the merge-
    # collinear-rects approach (Phase B2) remains open — it should
    # happen AFTER rect building (merge adjacent same-ref rects
    # whose axis angle differs by < 2° and whose joining edge sits
    # on rect interior, not at a junction polygon).
    osm_centerlines = _split_centerlines_at_points(
        osm_centerlines, junction_points, pav_union=None,
        approach_tol_m=25.0)

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
        filtered: List[Tuple[Polygon, LineString, str, str]] = []
        for rect, axis, role, ref in taxi_rects:
            if role != ROLE_STUB:
                filtered.append((rect, axis, role, ref))
                continue
            reaches_runway = False
            for (px, py) in raw_endpoints_by_ref.get(ref, []):
                if Point(px, py).distance(
                        layout.runway_union) <= RUNWAY_ENDPOINT_DIST_M:
                    reaches_runway = True
                    break
            if reaches_runway:
                filtered.append((rect, axis, role, ref))
        taxi_rects = filtered

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
    for ap in apron_polys:
        simp = ap.simplify(1.0, preserve_topology=True)
        if simp.is_empty or simp.geom_type != "Polygon":
            simp = ap
        layout.shapes.append(BuiltShape(polygon=simp, role=ROLE_APRON))

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

JUNCTION_CLUSTER_DIST_M = 80.0  # merge junction nodes within this distance
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
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
        ref = tags.get("ref", "") or "<unrefed>"
        for n in nds:
            refs_at_node[n].add(ref)

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
                pts = []
                if inter.geom_type == "Point":
                    pts.append((inter.x, inter.y))
                elif inter.geom_type == "MultiPoint":
                    pts.extend((p.x, p.y) for p in inter.geoms)
                elif inter.geom_type == "LineString":
                    pts.append(inter.centroid.coords[0])
                for pt in pts:
                    candidates.append(pt)

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


def _extract_osm_taxi_centerlines(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
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

        # Stage 2: gap bridging.  For NON-parallel refs (stubs,
        # cross-connectors, sub-refs), bridge gaps up to
        # GAP_BRIDGE_MAX_M.  For parallel refs (A/F/L/V/M/U), DON'T
        # bridge — user keeps their parallel polylines separate at
        # intermediate intersections.
        if ref and len(merged_lines) > 1 and ref not in PARALLEL_REFS:
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
            # Parallel refs (A/F/L/V/M/U) and cross-connector refs
            # (Q/R/X) split at bends — user draws multi-piece long
            # taxiways.  Everything else (stubs, sub-refs) emits as
            # ONE rect regardless of internal bends or length.
            # Unrefed airports: treat long polylines as parallels.
            SPLIT_AT_BENDS_REFS = PARALLEL_REFS | {"Q", "R", "X"}
            is_parallel = (
                ref in SPLIT_AT_BENDS_REFS or
                (ref == "" and simp.length > 500.0)
            )
            if is_parallel:
                for i in range(len(scoords) - 1):
                    try:
                        seg = LineString([scoords[i], scoords[i+1]])
                    except Exception:
                        continue
                    if seg.is_empty or seg.length < MIN_SEGMENT_LEN_M:
                        continue
                    out.append((seg, ref))
            else:
                # For stubs, use the ORIGINAL (un-simplified) polyline
                # so curves don't get chord-cut.  Pass the full line
                # in; downstream trim + rect-build works with it.
                if ls.length >= MIN_SEGMENT_LEN_M:
                    out.append((ls, ref))

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


def _split_centerlines_at_points(
    centerlines: List[Tuple[LineString, str]],
    split_points: List[Tuple[float, float]],
    pav_union: Optional[Polygon] = None,
    approach_tol_m: float = 25.0,
    endpoint_guard_m: float = 5.0,
    widening_factor: float = 1.15,
    adjacent_probe_m: float = 15.0,
) -> List[Tuple[LineString, str]]:
    """Split each centerline at every widening intersection.

    Per user rule (2026-04-18): cut at every curve or intersection,
    and "single rect between junctions on straight sections" —
    i.e. don't split where the pavement is STRAIGHT through the
    OSM multi-ref node.  A node counts as a widening intersection
    only if the pavement is materially wider there than along the
    adjacent centerline sections.

    For each candidate split point close to a centerline, we probe
    the pavement half-width at the split point and at a reference
    point ``adjacent_probe_m`` away along the line.  We only split
    if probe_at_split > ``widening_factor`` × probe_adjacent.
    """
    if not centerlines or not split_points:
        return centerlines
    from shapely.ops import substring

    pav_boundary = (pav_union.boundary
                    if pav_union is not None and not pav_union.is_empty
                    else None)

    def _hw_at_point(pt: Point) -> float:
        if pav_boundary is None:
            return 0.0
        try:
            return pt.distance(pav_boundary)
        except Exception:
            return 0.0

    result: List[Tuple[LineString, str]] = []
    for ls, ref in centerlines:
        # Collect split-parameters along the line for each split
        # point within ``approach_tol_m`` that shows actual
        # pavement widening.
        cut_params: List[float] = []
        for (sx, sy) in split_points:
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
            # Widening check: pav half-width at split point vs at
            # narrow reference points on BOTH sides of the split.
            # We split only if the split is wider than the
            # narrowest adjacent corridor — a true widening (i.e.
            # the corridor is narrow nearby, so this IS the junction).
            if pav_boundary is not None:
                split_pt = ls.interpolate(param)
                # Probe both sides, 10-30 m away, take narrow probe.
                side_probes: List[float] = []
                for step in (10.0, 20.0, 30.0):
                    for sign in (-1, 1):
                        t = param + sign * step
                        if t < endpoint_guard_m or t > ls.length - endpoint_guard_m:
                            continue
                        side_probes.append(_hw_at_point(ls.interpolate(t)))
                if side_probes:
                    hw_adj = min(side_probes)
                    hw_split = _hw_at_point(split_pt)
                    if hw_adj > 0 and hw_split < widening_factor * hw_adj:
                        # No widening — pavement runs straight through
                        # this OSM node.  Don't split.
                        continue
            cut_params.append(param)
        if not cut_params:
            result.append((ls, ref))
            continue
        # Dedupe near-duplicate params (within 2 m along the line).
        cut_params.sort()
        deduped: List[float] = []
        for p in cut_params:
            if not deduped or p - deduped[-1] > 2.0:
                deduped.append(p)
        # Build segments via substring between 0, cuts..., line.length
        breaks = [0.0] + deduped + [ls.length]
        for i in range(len(breaks) - 1):
            p0, p1 = breaks[i], breaks[i+1]
            if p1 - p0 < MIN_SEGMENT_LEN_M:
                continue
            try:
                piece = substring(ls, p0, p1)
            except Exception:
                continue
            if (piece.geom_type == "LineString"
                    and not piece.is_empty
                    and piece.length >= MIN_SEGMENT_LEN_M):
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

        # Trim endpoints inward where half-width exceeds widening threshold.
        # Per user rule (2026-04-18): "If pavement widens 1% cut the
        # rect."  Baseline = narrow_hw (p10), factor 1.01 → any
        # pavement more than 1% wider than the rect's own half-width
        # is already junction territory.
        trimmed = _trim_to_narrow(clipped, pav_non_rwy, narrow_hw,
                                 widen_factor=1.01)
        if trimmed is None or trimmed.length < 15.0:
            continue

        # Re-probe along the TRIMMED axis.  Per user rule: rect
        # half-width = narrowest probe within the covered segment;
        # anything wider belongs to the adjacent junction.
        _, _, trim_narrow_hw = _natural_half_width(trimmed, pav_non_rwy)
        if trim_narrow_hw < 3.5:
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

        role = _classify_role(trimmed, width, rwy_centerlines,
                               rwy_union, ref=ref)
        emitted.append((rect, trimmed, role, ref))
        emitted_union = (unary_union([emitted_union, rect])
                         if emitted_union is not None else rect)

    # Post-classify: secondary passes to fix roles based on neighbour
    # topology (stubs touch parallels, cross_connector touches 2 parallels)
    _refine_roles(emitted, rwy_centerlines)

    # Phase B2 (2026-04-19 principled) tested with width-uniformity
    # check and still regressed SPJC match count by 2 — target
    # subdivides at some joint that meets the uniformity test.
    # Leaving the function available as scaffolding; disabled.
    # emitted = _merge_collinear_rects_principled(
    #     emitted, pav_non_rwy, apt_vertices=apt_vertices)
    return emitted


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
                        n_probes: int = 11) -> Tuple[float, float, float]:
    """Return (natural_hw, max_hw, narrow_hw) distance-to-boundary along axis.

    * ``natural_hw`` = 25-th-percentile probe distance.  Used by the
      trim algorithm as the widening-detector baseline.
    * ``max_hw`` = 90-th-percentile probe distance.  Retained for
      legacy callers that still want the "widest section" width.
    * ``narrow_hw`` = min(p10, 2nd-lowest) probe distance.  Per the
      user's authoritative rule (memory: feedback_shape_rules),
      the rect half-width is the NARROWEST pavement width along the
      covered segment — the widened portions near each intersection
      are junction territory, not rect territory.
    """
    if axis.length < 1e-3:
        return 0.0, 0.0, 0.0
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
        return 0.0, 0.0, 0.0
    dists.sort()
    median = dists[len(dists) // 2]
    p90_idx = max(0, int(len(dists) * 0.9) - 1)
    p90 = dists[p90_idx] if p90_idx < len(dists) else dists[-1]
    # Narrow = 10th-percentile (or second-lowest for small n).  Pulling
    # the absolute min sometimes grabs a pinch from a neighbouring
    # building; p10 is a robust "narrowest straight section" proxy.
    p10_idx = max(0, int(len(dists) * 0.1))
    narrow = dists[p10_idx] if p10_idx < len(dists) else dists[0]
    return median, p90, narrow


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
    snapped = _snap_corners_to_pavement(corners, pav, apt_vertices)
    # Reject degenerate rects where snap collapsed two corners onto
    # the same apt.dat vertex (produces a zero-area or triangle-
    # shaped polygon that breaks rule 7 "corners on pav boundary"
    # by having 2 coincident corners).
    for i in range(4):
        for j in range(i + 1, 4):
            if math.hypot(snapped[i][0] - snapped[j][0],
                          snapped[i][1] - snapped[j][1]) < 1.0:
                return None
    return Polygon(snapped)


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
        if length >= 80.0:
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
