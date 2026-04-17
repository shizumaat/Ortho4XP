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
from shapely.ops import linemerge, transform as shp_transform, unary_union

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
                                       List[Tuple[str, List[str], Dict[str, str]]]]:
    """Parse an Ortho4XP-cached OSM tile (.osm.bz2 or .osm).

    Returns (nodes, ways) where:
      * nodes: {id: (lat, lon)}
      * ways:  [(id, [nd_ref, ...], {tag: val})]
    """
    if path.endswith(".bz2"):
        with bz2.open(path, "rt") as f:
            txt = f.read()
    else:
        txt = Path(path).read_text()

    # Quotes may be ' or "
    node_re = re.compile(
        r"""<node\s+id=["'](-?\d+)["'][^>]*?lat=["']([^"']+)["']\s+lon=["']([^"']+)["']"""
    )
    way_re = re.compile(r"""<way[^>]*?id=["'](-?\d+)["'][^>]*>(.*?)</way>""", re.S)
    nd_re = re.compile(r"""<nd\s+ref=["'](-?\d+)["']""")
    tag_re = re.compile(r"""<tag\s+k=["']([^"']+)["']\s+v=["']([^"']+)["']""")

    nodes: Dict[str, Tuple[float, float]] = {}
    for m in node_re.finditer(txt):
        try:
            nodes[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        except ValueError:
            continue

    ways: List[Tuple[str, List[str], Dict[str, str]]] = []
    for m in way_re.finditer(txt):
        wid = m.group(1)
        body = m.group(2)
        nds = nd_re.findall(body)
        tags = dict(tag_re.findall(body))
        ways.append((wid, nds, tags))
    return nodes, ways


def _load_osm_airports(xplane_root: str, icao: str,
                       apt_lat: float, apt_lon: float,
                       radius_deg: float = 0.05
                       ) -> Tuple[Dict[str, Tuple[float, float]],
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
    seen_paths = set()
    for dlat in (0, -1, 1):
        for dlon in (0, -1, 1):
            osm_path = _tile_path(base_lat + dlat, base_lon + dlon)
            if osm_path in seen_paths or not os.path.isfile(osm_path):
                continue
            seen_paths.add(osm_path)
            n2, w2 = _load_osm_tile(osm_path)
            nodes.update(n2)
            ways.extend(w2)
    if not nodes:
        return {}, []

    # Filter by bbox — keep only ways whose node centroids are near
    # the airport center.
    def _in_box(lat, lon):
        return (abs(lat - apt_lat) <= radius_deg and
                abs(lon - apt_lon) <= radius_deg)

    kept_ways = []
    for wid, nds, tags in ways:
        pts = [nodes[n] for n in nds if n in nodes]
        if not pts:
            continue
        clat = sum(p[0] for p in pts) / len(pts)
        clon = sum(p[1] for p in pts) / len(pts)
        if _in_box(clat, clon):
            kept_ways.append((wid, nds, tags))
    return nodes, kept_ways


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

    # ── Load OSM centerlines ─────────────────────────────────────
    nodes, ways = _load_osm_airports(
        xplane_root, icao, anchor[0], anchor[1])
    osm_centerlines = _extract_osm_taxi_centerlines(nodes, ways, to_m)

    # ── Build taxi rects from centerlines ────────────────────────
    taxi_rects = _build_taxi_rects(
        osm_centerlines, pav_union, layout.runway_union, rwy_centerlines)

    emitted_taxi_rects: List[Polygon] = []
    for rect, axis, role, ref in taxi_rects:
        emitted_taxi_rects.append(rect)
        layout.shapes.append(BuiltShape(
            polygon=rect, role=role, ref=ref, source_axis=axis))

    # ── Aprons: pavement minus taxi rects, one per CC ────────────
    taxi_rect_union = (unary_union(emitted_taxi_rects)
                       if emitted_taxi_rects else None)
    if pav_union is not None:
        apron = pav_union
        if taxi_rect_union is not None:
            apron = apron.difference(taxi_rect_union)
        parts = [apron] if apron.geom_type == "Polygon" else list(
            getattr(apron, "geoms", []))
        for part in parts:
            if part.geom_type != "Polygon":
                continue
            if part.area < 500.0:
                continue
            layout.shapes.append(BuiltShape(
                polygon=part, role=ROLE_APRON))

    return layout


# ──────────────────────────────────────────────────────────────────
# Centerline-based taxi rect builder
# ──────────────────────────────────────────────────────────────────

def _extract_osm_taxi_centerlines(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
) -> List[Tuple[LineString, str]]:
    """Return list of (linestring_m, ref_tag) for aeroway=taxiway ways."""
    out = []
    for wid, nds, tags in ways:
        if tags.get("aeroway") != "taxiway":
            continue
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
        out.append((ls, tags.get("ref", "")))
    return out


def _build_taxi_rects(
    centerlines: List[Tuple[LineString, str]],
    pav_union: Optional[Polygon],
    rwy_union: Optional[Polygon],
    rwy_centerlines: List[LineString],
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Convert each usable centerline into a 4-vertex rect.

    Returns list of (rect, clipped_axis, role, ref).
    Overlapping rects (from parallel parking-position centerlines, or
    duplicate refs) are deduplicated: a rect whose centerline lies
    mostly inside an earlier-emitted rect is dropped.
    """
    if pav_union is None:
        return []

    # Subtract runway from pavement so centerlines near runway stop at edge
    pav_non_rwy = pav_union
    if rwy_union is not None:
        pav_non_rwy = pav_non_rwy.difference(rwy_union)

    # Sort by length descending — longer wins when refs duplicate
    centerlines = sorted(centerlines, key=lambda x: -x[0].length)

    emitted: List[Tuple[Polygon, LineString, str, str]] = []
    emitted_union: Optional[Polygon] = None

    for axis, ref in centerlines:
        # Clip axis to pavement (outside runway).
        try:
            clipped = axis.intersection(pav_non_rwy)
        except Exception:
            continue
        if clipped.is_empty:
            continue
        # If MultiLine, take longest
        if clipped.geom_type == "MultiLineString":
            clipped = max(clipped.geoms, key=lambda g: g.length)
        if clipped.geom_type != "LineString":
            continue
        if clipped.length < 30.0:
            continue

        # Width probe along clipped axis
        width = _probe_axis_width(clipped, pav_non_rwy)
        if width < 7.0 or width > 80.0:
            continue

        # Dedup: skip if axis lies >70 % inside already-emitted union
        if emitted_union is not None and not emitted_union.is_empty:
            try:
                inside_len = clipped.intersection(emitted_union).length
                if inside_len / clipped.length > 0.7:
                    continue
            except Exception:
                pass

        # Rect: endpoints extended so the rect touches pavement edges
        # at both ends (avoids mid-pavement square edges).
        rect = _rect_from_axis_extended(clipped, width, pav_non_rwy)
        if rect is None or rect.is_empty:
            continue

        role = _classify_role(clipped, width, rwy_centerlines, rwy_union)
        emitted.append((rect, clipped, role, ref))
        emitted_union = (unary_union([emitted_union, rect])
                         if emitted_union is not None else rect)

    # Post-classify: secondary passes to fix roles based on neighbour
    # topology (stubs touch parallels, cross_connector touches 2 parallels)
    _refine_roles(emitted, rwy_centerlines)
    return emitted


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
    """Build a rect around axis with endpoints extended to pavement edge.

    Extension makes the rect reach both ends of the paved strip so
    shared corners with neighbours (runway, other taxis) land on
    apt.dat pavement boundaries.
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

    # Extend endpoints outward along axis until we leave the pavement
    # (then back off 1 m).  Cap extension so we don't shoot far.
    def _extend(x, y, dir_x, dir_y, max_d=40.0):
        d = 0.0
        last_inside = 0.0
        while d < max_d:
            d += 0.5
            qx = x + dir_x * d
            qy = y + dir_y * d
            if pav.contains(Point(qx, qy)):
                last_inside = d
            else:
                break
        return (x + dir_x * last_inside, y + dir_y * last_inside)

    p1e = _extend(p1[0], p1[1], -ux, -uy)
    p2e = _extend(p2[0], p2[1], ux, uy)

    px, py = -uy, ux
    half = width / 2.0
    return Polygon([
        (p1e[0] + px * half, p1e[1] + py * half),
        (p2e[0] + px * half, p2e[1] + py * half),
        (p2e[0] - px * half, p2e[1] - py * half),
        (p1e[0] - px * half, p1e[1] - py * half),
    ])


def _classify_role(axis: LineString, width: float,
                   rwy_centerlines: List[LineString],
                   rwy_union: Optional[Polygon]) -> str:
    """Baseline role classification.

    * parallel to runway, within 300 m, length ≥ 500 m → primary_parallel
    * parallel to runway, otherwise → secondary_parallel
    * perpendicular to runway → cross_connector
    * short (length < 300 m) + not cross → stub
    """
    if not rwy_centerlines:
        return ROLE_STUB

    rwy = rwy_centerlines[0]
    for r in rwy_centerlines:
        if axis.distance(r) < rwy.distance(axis):
            rwy = r

    # Bearing difference modulo 180
    def _bearing(ls):
        c = list(ls.coords)
        return math.degrees(math.atan2(c[-1][0] - c[0][0],
                                       c[-1][1] - c[0][1])) % 180.0
    db = abs(_bearing(axis) - _bearing(rwy))
    db = min(db, 180.0 - db)

    # Closest distance from axis midpoint to runway
    mid = axis.interpolate(0.5, normalized=True)
    dist_to_rwy = mid.distance(rwy)

    if db < 20.0:  # parallel
        if axis.length >= 500.0 and dist_to_rwy < 300.0:
            return ROLE_PRIMARY_PARALLEL
        return ROLE_SECONDARY_PARALLEL
    if 70.0 < db <= 90.0:  # perpendicular
        if axis.length >= 150.0:
            return ROLE_CROSS_CONNECTOR
    return ROLE_STUB


def _refine_roles(emitted, rwy_centerlines):
    """Second pass: cross_connector must touch ≥ 2 parallels, else demote
    to stub.  stubs touching no parallel get demoted to 'nowhere' and
    removed — actually we leave them as stub for now and let the
    operator decide from the output.
    """
    # Find parallels
    parallels = [(i, s) for i, s in enumerate(emitted)
                 if s[2] in (ROLE_PRIMARY_PARALLEL,
                             ROLE_SECONDARY_PARALLEL)]
    for i, (rect, axis, role, ref) in enumerate(emitted):
        if role != ROLE_CROSS_CONNECTOR:
            continue
        # Does it touch ≥ 2 parallel rects?
        touching = 0
        for j, ps in parallels:
            if j == i:
                continue
            if rect.intersects(ps[0]):
                touching += 1
        if touching < 2:
            emitted[i] = (rect, axis, ROLE_STUB, ref)
