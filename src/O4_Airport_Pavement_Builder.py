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
# Constants (re-exported from O4_Pavement_Config + O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Config import (
    JUNCTION_CLUSTER_DIST_M,
    MIN_SEGMENT_LEN_M,
    SLIVER_ANGLE_THRESHOLD_DEG,
    EMIT_APRONS,
    EMIT_BRIDGES_AND_TUNNELS,
    EMIT_JUNCTIONS,
    LOAD_DSF_PAVEMENT,
    RUNWAY_APRON_AREA_RATIO,
    RUNWAY_INSIDE_APRON_FRAC,
)
from O4_Pavement_Layout import (
    AEROWAY_FOR_ROLE,
    BuiltShape,
    PavementLayout,
    R_EARTH,
    ROLE_APRON,
    ROLE_BOUNDARY,
    ROLE_CROSS_CONNECTOR,
    ROLE_GROUNDSIDE_PAVEMENT,
    ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RETAINING_WALL,
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
    ROLE_TUNNEL_RAMP,
    SHARED_VERTEX_TOL_M,
    _airport_anchor,
    _projection,
)
from O4_Pavement_Vertices import (
    _drop_spike_vertices,
    _enforce_shared_vertices,
    _push_junction_vertices_off_taxi_rect_edges,
    _snap_polygon_vertices_to_rect_corners,
    _validate_shared_vertex_invariant,
)



# ──────────────────────────────────────────────────────────────────
# Data model (re-exported from O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────


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


def _osm_tile_path(lat_tile: int, lon_tile: int,
                    cached_suffix: str = "airports") -> str:
    """Ortho4XP's per-tile cached OSM filename layout."""
    def _fmt(v, pad):
        sign = "+" if v >= 0 else "-"
        return f"{sign}{abs(v):0{pad}d}"
    lat_group = (lat_tile // 10) * 10
    lon_group = (lon_tile // 10) * 10
    return os.path.join(
        "OSM_data",
        f"{_fmt(lat_group, 2)}{_fmt(lon_group, 3)}",
        f"{_fmt(lat_tile, 2)}{_fmt(lon_tile, 3)}",
        f"{_fmt(lat_tile, 2)}{_fmt(lon_tile, 3)}"
        f"_{cached_suffix}.osm.bz2",
    )


def _load_osm_airports(xplane_root: str, icao: str,
                       apt_lat: float, apt_lon: float,
                       radius_deg: float = 0.05
                       ) -> Tuple[Dict[str, Tuple[float, float]],
                                  List[Tuple[str, List[str], Dict[str, str]]],
                                  List[Tuple[str, List[str], Dict[str, str]]]]:
    """Load the airports-layer OSM cache covering the given lat/lon.

    Returns nodes + ways filtered to a bbox around the airport.
    """
    def _tile_path(lat_tile: int, lon_tile: int) -> str:
        return _osm_tile_path(lat_tile, lon_tile, "airports")

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


def _score_apt_dat_against_osm(
        apt_path: str,
        icao: str,
        nodes: Dict[str, Tuple[float, float]],
        ways: List[Tuple[str, List[str], Dict[str, str]]],
        taxi_buffer_m: float = 5.0,
        ) -> Tuple[float, float]:
    """Score how well ``apt_path`` covers OSM-known features at
    this airport.

    Returns ``(apron_coverage, taxi_coverage)`` where:

      * ``apron_coverage`` = (area of OSM ``aeroway=apron`` polygons
        that lies inside the apt.dat row-110 pavement union) /
        (total area of OSM apron polygons).  1.0 = perfect, 0.0 =
        no apt.dat coverage of OSM-known apron areas.

      * ``taxi_coverage`` = (length of OSM ``aeroway=taxiway``
        centerline LineStrings within ``taxi_buffer_m`` of any
        apt.dat pavement) / (total length of OSM taxi
        centerlines).  Buffer matches typical taxiway half-width.
        1.0 = every OSM-known taxiway centerline has apt.dat
        pavement underneath it.

    Returns (1.0, 1.0) if no OSM apron / taxi data exists (the
    apt.dat passes by default).  Returns (0.0, 0.0) if the apt.dat
    can't be loaded.

    Per user 2026-04-30: when a custom-scenery apt.dat scores
    below 0.7 on either metric, it's "missing a lot vs OSM" and
    we fall back to the global apt.dat instead.
    """
    try:
        apt = APR.load_airport(apt_path, icao)
    except Exception:
        return (0.0, 0.0)
    if apt is None:
        return (0.0, 0.0)
    if not apt.pavements and not apt.runways:
        return (0.0, 0.0)
    # Anchor + projection (use first runway threshold).
    if not apt.runways:
        return (1.0, 1.0)
    r0 = apt.runways[0]
    lat0, lon0 = r0.lat_a, r0.lon_a
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)
    # Build apt.dat pavement union in meter space.  apt.dat
    # polygons are stored in lat/lon — project to meters for
    # consistent area / distance math.
    from shapely.ops import transform as shp_transform
    pav_polys_m: List[Polygon] = []
    for pav in apt.pavements:
        if pav.polygon is None or pav.polygon.is_empty:
            continue
        try:
            pm = shp_transform(
                lambda lon, lat, z=None:
                    to_m(lon, lat) if z is None
                    else (*to_m(lon, lat), z),
                pav.polygon)
            if pm.is_empty:
                continue
            if not pm.is_valid:
                pm = pm.buffer(0)
            if pm.geom_type == "Polygon" and not pm.is_empty:
                pav_polys_m.append(pm)
            elif pm.geom_type == "MultiPolygon":
                for g in pm.geoms:
                    if (g.geom_type == "Polygon"
                            and not g.is_empty):
                        pav_polys_m.append(g)
        except Exception:
            continue
    # Also include runway rects so a taxi centerline that ends on
    # the runway counts as covered.
    for rwy in apt.runways:
        try:
            rect = _runway_rect_m(rwy, to_m)
            if rect is not None and not rect.is_empty:
                pav_polys_m.append(rect)
        except Exception:
            continue
    if not pav_polys_m:
        return (0.0, 0.0)
    try:
        pav_union = unary_union(pav_polys_m)
    except Exception:
        return (0.0, 0.0)
    if pav_union.is_empty:
        return (0.0, 0.0)
    pav_buf = pav_union.buffer(taxi_buffer_m)
    # Apron coverage: OSM apron polygons (closed ways tagged
    # aeroway=apron) → fraction inside pav_union.
    apron_total = 0.0
    apron_inside = 0.0
    taxi_total = 0.0
    taxi_inside = 0.0
    for wid, nrefs, tags in ways:
        ay = tags.get("aeroway", "")
        if ay == "apron":
            pts = []
            for n in nrefs:
                if n in nodes:
                    la, lo = nodes[n]
                    pts.append(to_m(lo, la))
            if (len(pts) >= 4
                    and abs(pts[0][0] - pts[-1][0]) < 0.5
                    and abs(pts[0][1] - pts[-1][1]) < 0.5):
                try:
                    poly = Polygon(pts)
                    if not poly.is_valid:
                        poly = poly.buffer(0)
                    if (poly.is_empty
                            or poly.geom_type != "Polygon"):
                        continue
                    apron_total += poly.area
                    inter = poly.intersection(pav_union)
                    if not inter.is_empty:
                        apron_inside += inter.area
                except Exception:
                    continue
        elif ay == "taxiway":
            # Open ways (centerlines).  Closed taxi polygons are
            # rare in OSM; treat as line either way.
            pts = []
            for n in nrefs:
                if n in nodes:
                    la, lo = nodes[n]
                    pts.append(to_m(lo, la))
            if len(pts) < 2:
                continue
            try:
                from shapely.geometry import LineString as _LS
                ls = _LS(pts)
                if ls.is_empty or ls.length < 1.0:
                    continue
                taxi_total += ls.length
                inter = ls.intersection(pav_buf)
                if not inter.is_empty:
                    if hasattr(inter, "length"):
                        taxi_inside += inter.length
                    elif hasattr(inter, "geoms"):
                        for g in inter.geoms:
                            if hasattr(g, "length"):
                                taxi_inside += g.length
            except Exception:
                continue
    apron_cov = (apron_inside / apron_total
                  if apron_total > 0 else 1.0)
    taxi_cov = (taxi_inside / taxi_total
                 if taxi_total > 0 else 1.0)
    return (apron_cov, taxi_cov)


def _pick_best_apt_dat_against_osm(
        xplane_root: str,
        icao: str,
        apron_threshold: float = 0.7,
        taxi_threshold: float = 0.7,
        ) -> Optional[str]:
    """Find the best apt.dat for ``icao``, falling back from a
    sparse custom-scenery pack to the global apt.dat when the
    custom one is missing too much OSM-known geometry.

    Per user 2026-04-30: when a custom apt.dat is missing a lot
    of pavement vs OSM, the X-Plane global definition often has
    more complete row-110 polygons.  This selector evaluates
    every candidate apt.dat (custom packs in priority order,
    then global, then default) by computing apron + taxi
    coverage against OSM, and picks the first candidate whose
    coverage clears both thresholds.  Falls back to the
    legacy first-with-pavement selection if no candidate
    qualifies.

    Thresholds:
      ``apron_threshold`` = 0.7 — apron polygons should be 70 %
        covered by apt.dat row-110 pavement.
      ``taxi_threshold`` = 0.7 — OSM taxi centerlines should be
        70 % covered by apt.dat row-110 pavement (buffered 5 m).
    """
    candidates = APR.find_all_airport_apt_dats(xplane_root, icao)
    if not candidates:
        return APR.find_airport_apt_dat(xplane_root, icao)
    if len(candidates) == 1:
        return candidates[0]
    # Load OSM data once (use first candidate's airport as the
    # anchor — OSM tile loading needs a lat/lon hint).
    anchor_apt = None
    for cand in candidates:
        try:
            anchor_apt = APR.load_airport(cand, icao)
            if (anchor_apt is not None
                    and anchor_apt.runways):
                break
        except Exception:
            continue
    if anchor_apt is None or not anchor_apt.runways:
        return APR.find_airport_apt_dat(xplane_root, icao)
    r0 = anchor_apt.runways[0]
    try:
        nodes_o, ways_o, _ = _load_osm_airports(
            xplane_root, icao, r0.lat_a, r0.lon_a)
    except Exception:
        return APR.find_airport_apt_dat(xplane_root, icao)
    if not ways_o:
        return APR.find_airport_apt_dat(xplane_root, icao)
    scores: List[Tuple[str, float, float]] = []
    for cand in candidates:
        ac, tc = _score_apt_dat_against_osm(
            cand, icao, nodes_o, ways_o)
        scores.append((cand, ac, tc))
    # Walk in priority order; pick the first that clears both
    # thresholds.  Log every candidate's score.
    chosen: Optional[str] = None
    for cand, ac, tc in scores:
        passed = ac >= apron_threshold and tc >= taxi_threshold
        if passed and chosen is None:
            chosen = cand
        try:
            import sys as _sys
            label = "PICK" if (passed and cand == chosen) else (
                "ok" if passed else "skip")
            _sys.stderr.write(
                f"  [pav-builder] {icao}: apt.dat candidate "
                f"[{label}] apron={ac:.0%} taxi={tc:.0%}  "
                f"{cand}\n")
        except Exception:
            pass
    if chosen is not None:
        return chosen
    # Nothing qualified — pick whichever has the highest
    # combined coverage to avoid emitting nothing.
    if scores:
        scores.sort(key=lambda s: -(s[1] + s[2]))
        try:
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] {icao}: no apt.dat met "
                f"thresholds (apron≥{apron_threshold:.0%}, "
                f"taxi≥{taxi_threshold:.0%}); falling back to "
                f"highest combined coverage: "
                f"apron={scores[0][1]:.0%} taxi={scores[0][2]:.0%}.\n")
        except Exception:
            pass
        return scores[0][0]
    return APR.find_airport_apt_dat(xplane_root, icao)


def _load_osm_big_roads(apt_lat: float, apt_lon: float,
                        radius_deg: float = 0.05
                        ) -> Tuple[Dict[str, Tuple[float, float]],
                                   List[Tuple[str, List[str], Dict[str, str]]]]:
    """Load the ``big_roads`` OSM cache (motorway / trunk / primary /
    secondary / railway ways with ``tunnel`` and ``bridge`` tag
    annotations).  Same multi-tile namespace + bbox-filter logic as
    ``_load_osm_airports``.

    Returns ``(nodes, ways)``.  Relations are not used by the road
    pipeline.  Returns empty containers when no cache exists at
    this tile (tile builds without road data — e.g. SPLP, where
    big_roads.osm.bz2 was never generated — silently skip tunnel
    emission).
    """
    base_lat = int(math.floor(apt_lat))
    base_lon = int(math.floor(apt_lon))
    nodes: Dict[str, Tuple[float, float]] = {}
    ways: List[Tuple[str, List[str], Dict[str, str]]] = []
    seen_paths = set()
    for dlat in (0, -1, 1):
        for dlon in (0, -1, 1):
            tile_lat_n = base_lat + dlat
            tile_lon_n = base_lon + dlon
            osm_path = _osm_tile_path(
                tile_lat_n, tile_lon_n, "big_roads")
            if osm_path in seen_paths or not os.path.isfile(osm_path):
                continue
            seen_paths.add(osm_path)
            n2, w2, _r2 = _load_osm_tile(osm_path)
            tile_prefix = (
                f"r{tile_lat_n:+03d}{tile_lon_n:+04d}:")
            for nid, coord in n2.items():
                nodes[tile_prefix + nid] = coord
            for wid, nds, tags in w2:
                ways.append(
                    (tile_prefix + wid,
                     [tile_prefix + n for n in nds],
                     tags))
    if not nodes:
        return {}, []

    def _in_box(lat, lon):
        return (abs(lat - apt_lat) <= radius_deg
                and abs(lon - apt_lon) <= radius_deg)

    kept = []
    for wid, nds, tags in ways:
        pts = [nodes[n] for n in nds if n in nodes]
        if not pts:
            continue
        clat = sum(p[0] for p in pts) / len(pts)
        clon = sum(p[1] for p in pts) / len(pts)
        # A road may cross the airport bbox without its centroid
        # landing inside — keep the way if EITHER (centroid in box)
        # OR (any vertex within ``radius_deg``).
        if _in_box(clat, clon):
            kept.append((wid, nds, tags))
            continue
        for lat, lon in pts:
            if _in_box(lat, lon):
                kept.append((wid, nds, tags))
                break
    return nodes, kept


# ──────────────────────────────────────────────────────────────────
# Meter-space projection (re-exported from O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────



# ──────────────────────────────────────────────────────────────────
# Runway rects, crossings, shoulders (re-exported from
# O4_Pavement_Runways)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Runways import (
    _detect_runway_shoulders,
    _insert_runway_chain_bridges,
    _resolve_runway_crossings,
    _runway_rect_m,
    _sample_runway_segment_elev,
)


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


from O4_Pavement_Absorption import (
    _drop_primary_parallels_embedded_in_pavement,
    _split_primary_parallels_at_pavement_boundary,
)



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
    apt_path = _pick_best_apt_dat_against_osm(xplane_root, icao)
    if apt_path is None:
        raise RuntimeError(f"No apt.dat found for {icao}")
    apt = APR.load_airport(apt_path, icao)
    if apt is None:
        raise RuntimeError(f"Could not load airport block for {icao}")

    anchor = _airport_anchor(apt)
    to_m = _projection(anchor)

    layout = PavementLayout(icao=icao, anchor=anchor,
                             apt_dat_path=apt_path)

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
    # Per user 2026-04-29: drop DSF pavement polygons whose area
    # exceeds the largest apt.dat pavement polygon by more than
    # DSF_MAX_AREA_VS_APT_DAT_RATIO×.  Real airport pavement
    # polygons (taxiways, aprons, runway aprons) are bounded in
    # scale by the largest features apt.dat already represents at
    # the airport.  A DSF polygon dramatically larger than apt.dat's
    # biggest is a coarse "ground tile" — pavement-textured
    # decorative geometry painted across the whole airport surface
    # rather than a real pavement feature.  Confirmed at HECA
    # (Tai Models scenery), where ``lib/airport/ground/pavement/
    # asphalt/patched.pol`` instances of 3.5 M m² (with only 22
    # vertices, perimeter ~8 km) and 1.25 M m² (36 verts) overlay
    # the entire airport, dwarfing apt.dat's largest pavement
    # polygon at 378 k m² and inflating the rect-detection's
    # half-width probes to 100 m+ across what should be a 30 m
    # taxi corridor.  Safe at SPJC / CYXY / SPLP: their largest
    # legitimate DSF pavement polygons are within 2.2× apt.dat's
    # largest, well under the 3× cap.
    DSF_MAX_AREA_VS_APT_DAT_RATIO = 3.0
    apt_pav_union: Optional[Polygon] = None
    apt_pav_largest_area: float = 0.0
    if pav_polys:
        try:
            apt_pav_union = unary_union(pav_polys)
        except Exception:
            apt_pav_union = None
        apt_pav_largest_area = max(
            (p.area for p in pav_polys), default=0.0)

    # ── OSM-aeroway-footprint vs apt.dat coverage check ───────────
    # Per user 2026-04-29: prioritize apt.dat as the pavement
    # source.  Only fall back to DSF when there's a meaningful
    # discrepancy between apt.dat and what OSM aeroway data tells
    # us the airport actually has.  The OSM aeroway tags
    # (aeroway=apron, taxiway, taxi_lane, stand) are the user-
    # mapped truth about the airport's pavement extent.  If
    # apt.dat already covers that extent, DSF additions are at
    # best decorative overlays and at worst inflated ground tiles
    # filling non-pavement areas (HECA, where DSF added pavement
    # between runways and taxiways).  When apt.dat is missing
    # significant OSM-known pavement, DSF is admitted only inside
    # the gap.
    #
    # Load OSM here (early) so the DSF loop below can use the
    # aeroway footprint to gate DSF additions.
    nodes, ways, relations = _load_osm_airports(
        xplane_root, icao, anchor[0], anchor[1])
    osm_aeroway_footprint = _build_osm_aeroway_footprint(
        nodes, ways, to_m)
    # Gap = OSM-known pavement that apt.dat doesn't cover.  When
    # this is a small fraction, apt.dat is sufficient; skip DSF
    # entirely.
    osm_gap: Optional[Polygon] = None
    DSF_OSM_GAP_BUFFER_M = 1.0          # widen gap by 1 m for
                                         # tile-alignment slop.
                                         # User 2026-04-29 (HECA R
                                         # absorption): 5 m was too
                                         # generous; DSF clipped
                                         # with a 5 m fringe
                                         # extends ~3 m past the
                                         # OSM-tagged taxi corridor
                                         # and trips the long-edge-
                                         # adjacent absorption probe
                                         # (which fires at 2 m
                                         # outside the rect edge).
                                         # 1 m fringe matches the
                                         # apt.dat / DSF tile
                                         # alignment precision
                                         # without spilling enough
                                         # to look like apron-
                                         # adjacency.
    # Coverage threshold: above this, the airport's apt.dat is
    # considered "comprehensive" — apt.dat captures most of what
    # OSM thinks the airport has, so DSF is restricted to filling
    # the small remaining gap (typically a couple of taxi corridors
    # apt.dat happens to miss).  Below this threshold, apt.dat is
    # sparse (e.g. CYXY where apt.dat has ~52% of OSM) and OSM is
    # also incomplete; DSF is the primary pavement source there
    # and we trust it broadly (drop overlays only).
    APT_COMPREHENSIVE_OSM_FRAC = 0.80
    apt_is_comprehensive = False
    if (osm_aeroway_footprint is not None
            and not osm_aeroway_footprint.is_empty
            and apt_pav_union is not None
            and not apt_pav_union.is_empty):
        try:
            apt_in_osm = apt_pav_union.intersection(
                osm_aeroway_footprint).area
            osm_area = osm_aeroway_footprint.area
            if osm_area > 1.0:
                apt_is_comprehensive = (
                    apt_in_osm / osm_area
                    >= APT_COMPREHENSIVE_OSM_FRAC)
            gap = osm_aeroway_footprint.difference(apt_pav_union)
            if not gap.is_empty:
                osm_gap = gap.buffer(DSF_OSM_GAP_BUFFER_M)
        except Exception:
            pass
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
        if not LOAD_DSF_PAVEMENT:
            raise StopIteration  # skip the DSF block entirely
        import O4_DSF_Reader as _DSFR
        seen_dsf: set = set()
        all_apt_dats = APR.find_all_airport_apt_dats(xplane_root, icao)
        n_dsf_kept = 0
        n_dsf_dropped_overlay = 0
        n_dsf_dropped_far = 0
        n_dsf_dropped_oversized = 0
        n_dsf_dropped_outside_osm_gap = 0
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
                    # apt.dat-priority gate (user 2026-04-29):
                    # when apt.dat is comprehensive (covers ≥ 80%
                    # of the OSM-aeroway footprint), CLIP each DSF
                    # polygon to the buffered OSM-vs-apt.dat gap.
                    # An intersect-test alone is too lax: a wide
                    # DSF polygon that grazes the gap by 1 m² gets
                    # kept entirely, and its 50–80 m extension
                    # past the gap drags non-pavement coverage
                    # into pav_union (HECA F case where the DSF
                    # along F's corridor extends 80 m onto the
                    # adjacent ramp, inflating ``_natural_half_
                    # width`` to ~50 m and triggering the apron-
                    # interior corner check on every F segment).
                    # Clipping keeps only the part of the DSF
                    # polygon that's actually filling an OSM-
                    # tagged corridor apt.dat happens to miss.
                    # When apt.dat is SPARSE (< 80 % of OSM), we
                    # skip the clip — apt.dat alone is too thin
                    # and OSM is also incomplete, so DSF is the
                    # primary source and we trust it broadly
                    # (CYXY, SPLP).
                    if (apt_is_comprehensive
                            and osm_gap is not None
                            and not osm_gap.is_empty):
                        try:
                            clipped_pm = pm.intersection(osm_gap)
                            if clipped_pm.is_empty:
                                n_dsf_dropped_outside_osm_gap += 1
                                continue
                            # Take the largest Polygon piece if
                            # the clip produced a MultiPolygon —
                            # narrow slivers from a wide DSF
                            # polygon grazing the gap aren't
                            # useful pavement either.
                            if (clipped_pm.geom_type
                                    == "MultiPolygon"):
                                clipped_pm = max(
                                    clipped_pm.geoms,
                                    key=lambda g: g.area)
                            if (clipped_pm.geom_type != "Polygon"
                                    or clipped_pm.is_empty
                                    or clipped_pm.area < 5.0):
                                n_dsf_dropped_outside_osm_gap += 1
                                continue
                            pm = clipped_pm
                        except Exception:
                            pass
                    # Oversized-vs-apt.dat gate: a DSF polygon
                    # dramatically larger than the airport's
                    # biggest apt.dat pavement polygon is a coarse
                    # "ground tile" overlay, not real pavement —
                    # drop it.  Only meaningful when apt.dat has
                    # ANY pavement; airports with no apt.dat
                    # pavement (CYXY-style sparse data) are
                    # unaffected.
                    if (apt_pav_largest_area > 0
                            and pm.area
                            > (apt_pav_largest_area
                               * DSF_MAX_AREA_VS_APT_DAT_RATIO)):
                        n_dsf_dropped_oversized += 1
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
                or n_dsf_dropped_far or n_dsf_dropped_oversized
                or n_dsf_dropped_outside_osm_gap):
            try:
                import sys as _sys
                msg = (f"  [pav-builder] {icao}: DSF pavement: "
                       f"{n_dsf_kept} kept, "
                       f"{n_dsf_dropped_overlay} dropped as overlay, "
                       f"{n_dsf_dropped_far} dropped as off-airport")
                if n_dsf_dropped_oversized:
                    msg += (f", {n_dsf_dropped_oversized} dropped "
                            f"as oversized-vs-apt.dat")
                if n_dsf_dropped_outside_osm_gap:
                    msg += (f", {n_dsf_dropped_outside_osm_gap} "
                            f"dropped: outside OSM-aeroway gap")
                _sys.stderr.write(msg + ".\n")
            except Exception:
                pass
    except StopIteration:
        # DSF read intentionally disabled.
        try:
            import sys as _sys
            _sys.stderr.write(
                f"  [pav-builder] {icao}: DSF pavement read "
                f"disabled (LOAD_DSF_PAVEMENT=False).\n")
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

    # ── OSM centerlines from already-loaded OSM ──────────────────
    # ``nodes`` / ``ways`` / ``relations`` were loaded earlier so
    # the DSF-loading loop could compare apt.dat coverage against
    # the OSM-aeroway footprint.
    osm_centerlines = _extract_osm_taxi_centerlines(
        nodes, ways, to_m, rwy_centerlines=rwy_centerlines)

    # ── Terminal groundside-pavement subtraction (user 2026-04-29):
    # remove curbside / drop-off / parking pavement from pav_union
    # before downstream rect / junction construction sees it.
    # Groundside pavement sits at a different elevation than the
    # building's airside apron, so allowing it to become an apron
    # junction grade-clamps the building to the wrong altitude.
    # Subtract a perpendicular outward strip from each terminal
    # building's groundside edges (classified by OSM aeroway /
    # highway adjacency + apt.dat-pavement connectivity).
    try:
        _osm_terminal_buildings = _extract_osm_terminals(
            nodes, ways, relations, to_m)
        # Per user 2026-04-30 (CYXY -10123 NW phantom groundside):
        # pass the FULL apt.dat pavement polygon list so the
        # classifier can BFS from runway-touching polys through
        # transitive touches.  Without this, edges of the
        # terminal next to apron pavement that's not directly
        # touching a runway (e.g. the apron extends NW past the
        # terminal) classified UNKNOWN and got mis-promoted to
        # groundside.
        _ground_zone = _terminal_groundside_zone(
            _osm_terminal_buildings, nodes, ways, to_m,
            apt_pavement_seeds=runway_polys,
            apt_pavement_polys=apt_only_pav_polys)
        if _ground_zone is not None and not _ground_zone.is_empty:
            # Capture the groundside-only visible pavement BEFORE
            # the subtraction below empties pav_union of it.  Per
            # user 2026-04-29: groundside pavement should remain
            # in the output but follow DEM (with a 0.1 m gap from
            # the terminal building) rather than being flattened
            # to airside-apron elevation.  The shapes captured
            # here are emitted later with per-vertex DEM altitudes;
            # the 0.1 m terminal gap is enforced inside the emit
            # function using the LAYOUT's terminal shapes (which
            # may differ slightly from the OSM source extracts due
            # to apt.dat row-110 / DSF residue absorption).
            try:
                _groundside_visible = pav_union.intersection(
                    _ground_zone)
                _gs_polys: List[Polygon] = []
                if _groundside_visible is not None:
                    if _groundside_visible.geom_type == "Polygon":
                        if (not _groundside_visible.is_empty
                                and _groundside_visible.area >= 5.0):
                            _gs_polys.append(_groundside_visible)
                    elif (_groundside_visible.geom_type
                            == "MultiPolygon"):
                        for _g in _groundside_visible.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 5.0):
                                _gs_polys.append(_g)
                # Stash on the layout so the elevation pass can
                # find them once DEM is loaded.
                layout._groundside_polys = _gs_polys
            except Exception:
                layout._groundside_polys = []
            try:
                pav_union = pav_union.difference(_ground_zone)
                if hasattr(layout, "_pav_union_for_rects"):
                    layout._pav_union_for_rects = (
                        layout._pav_union_for_rects.difference(
                            _ground_zone))
                # Also subtract from the granular polygon list so
                # downstream apron-merged-runway / rect-corner
                # detection sees a consistent pavement footprint.
                _new_pav_polys: List[Polygon] = []
                for _p in pav_polys:
                    try:
                        _q = _p.difference(_ground_zone)
                    except Exception:
                        _new_pav_polys.append(_p)
                        continue
                    if _q.is_empty:
                        continue
                    if _q.geom_type == "Polygon":
                        _new_pav_polys.append(_q)
                    elif _q.geom_type == "MultiPolygon":
                        for _g in _q.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 1.0):
                                _new_pav_polys.append(_g)
                pav_polys[:] = _new_pav_polys
                # Same for apt_only_pav_polys (used for terminal
                # containment + rect-corner snapping).
                _new_apt_only: List[Polygon] = []
                for _p in apt_only_pav_polys:
                    try:
                        _q = _p.difference(_ground_zone)
                    except Exception:
                        _new_apt_only.append(_p)
                        continue
                    if _q.is_empty:
                        continue
                    if _q.geom_type == "Polygon":
                        _new_apt_only.append(_q)
                    elif _q.geom_type == "MultiPolygon":
                        for _g in _q.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 1.0):
                                _new_apt_only.append(_g)
                apt_only_pav_polys[:] = _new_apt_only
                # Also subtract from apron_candidates — captured at
                # line 2519 BEFORE this subtract — so apron-junction
                # construction in _compute_elevations can't wrap
                # around the terminal into the groundside zone (per
                # user 2026-04-29 / way -10125 vs -10111: the apron
                # junction was extending past the terminal to the
                # upper-side roads, sharing an edge with the new
                # DEM-following groundside pavement and creating a
                # 7 m vertical cliff at CYXY).
                _new_apron_cand: List[Polygon] = []
                for _p in apron_candidates:
                    try:
                        _q = _p.difference(_ground_zone)
                    except Exception:
                        _new_apron_cand.append(_p)
                        continue
                    if _q.is_empty:
                        continue
                    if _q.geom_type == "Polygon":
                        _new_apron_cand.append(_q)
                    elif _q.geom_type == "MultiPolygon":
                        for _g in _q.geoms:
                            if (_g.geom_type == "Polygon"
                                    and not _g.is_empty
                                    and _g.area >= 1.0):
                                _new_apron_cand.append(_g)
                apron_candidates[:] = _new_apron_cand
                try:
                    import sys as _sys
                    _sys.stderr.write(
                        f"  [pav-builder] {icao}: subtracted "
                        f"{_ground_zone.area:,.0f} m² of "
                        f"groundside pavement (terminal "
                        f"curbside / drop-off / parking).\n")
                except Exception:
                    pass
            except Exception:
                pass
    except Exception:
        pass

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

    # ── Detect bridge taxi rects from OSM (user 2026-04-29) ───────
    # Per OSM convention, bridge taxiways carry ``bridge=yes`` (or
    # any non-empty bridge=*) on their parent way.  Build a list of
    # raw OSM bridge-tagged taxiway LineStrings; a rect is a
    # "bridge rect" if its axis lies within a small lateral
    # tolerance of one of those LineStrings.  We use spatial
    # proximity rather than tag-pass-through because
    # ``_extract_osm_taxi_centerlines`` linemerges per-ref, which
    # loses the per-way bridge tag.
    bridge_lines: List[LineString] = []
    for _wid, _nds, _tags in ways:
        if _tags.get("aeroway") != "taxiway":
            continue
        b = _tags.get("bridge", "")
        if not b or b == "no":
            continue
        _pts = []
        for _n in _nds:
            if _n in nodes:
                _lat, _lon = nodes[_n]
                _pts.append(to_m(_lon, _lat))
        if len(_pts) < 2:
            continue
        try:
            _ls = LineString(_pts)
        except Exception:
            continue
        if _ls.is_empty or _ls.length < 5.0:
            continue
        bridge_lines.append(_ls)
    bridge_rect_indices: set = set()
    BRIDGE_AXIS_PROXIMITY_M = 5.0
    if bridge_lines:
        for ri, (_rect, _axis, _role, _ref) in enumerate(taxi_rects):
            if _axis is None or _axis.is_empty:
                continue
            for _bl in bridge_lines:
                # An axis is on a bridge if it lies within
                # BRIDGE_AXIS_PROXIMITY_M of the bridge LineString
                # along most of its length AND has overlap.
                try:
                    if _axis.distance(_bl) > BRIDGE_AXIS_PROXIMITY_M:
                        continue
                    inter = _axis.intersection(
                        _bl.buffer(BRIDGE_AXIS_PROXIMITY_M))
                    if (not inter.is_empty
                            and hasattr(inter, "length")
                            and inter.length
                            >= 0.5 * _axis.length):
                        bridge_rect_indices.add(ri)
                        break
                except Exception:
                    continue

    # Emit taxi rects (already trimmed to narrow-width portion).
    emitted_taxi_rects: List[Polygon] = []
    for ri, (rect, axis, role, ref) in enumerate(taxi_rects):
        emitted_taxi_rects.append(rect)
        layout.shapes.append(BuiltShape(
            polygon=rect, role=role, ref=ref, source_axis=axis,
            is_bridge=(ri in bridge_rect_indices)))

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
        # Per user 2026-04-29: merge small junction slivers into
        # adjacent larger junctions.  Polygon-with-holes
        # decomposition + post-elevation subdivisions can carve
        # off small (< 1000 m²) pieces sharing a boundary segment
        # with a much larger neighbour (HECA -10244 = 485 m²
        # adjacent to -10243 = 30 k m²).  Merge them back so
        # JOSM doesn't show two near-duplicate polygons.
        _merge_sliver_junctions_into_neighbours(layout, icao=icao)
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
            # Per user 2026-04-29: re-emit groundside pavement
            # captured before the airside-apron subtraction, with
            # per-vertex DEM altitudes and a 0.1 m gap from the
            # terminal building.  Lets curbside / drop-off /
            # parking pavement render at local terrain elevation
            # (CYXY: terminal cut into hill, ~4 m higher than the
            # airside apron).
            try:
                n_gs = _emit_groundside_pavement_dem(
                    layout, _dem, _tile_lat, _tile_lon)
                if n_gs:
                    import sys as _sys
                    _sys.stderr.write(
                        f"  [pav-builder] emitted "
                        f"{n_gs} groundside pavement "
                        f"polygon(s) with DEM altitudes.\n")
            except Exception:
                pass
            # Per user 2026-04-29 (CYXY -10111 / -10115): drop
            # junction polygons that are connected ONLY to non-
            # airside pavement (groundside polygons or each
            # other), with no shared vertices on any airside
            # rect/terminal/runway.  These slivers got assigned
            # the airside-flat altitude during the Laplacian
            # solver but visually they're not airside — they
            # share an edge with the DEM-following groundside
            # at a different elevation, creating "valleys" that
            # X-Plane renders as cliffs.  Eliminating them lets
            # the boundary ribbon and the groundside polygon
            # define the surface in that area.
            try:
                n_orph = _drop_groundside_orphan_junctions(layout)
                if n_orph:
                    import sys as _sys
                    _sys.stderr.write(
                        f"  [pav-builder] dropped {n_orph} "
                        f"junction(s) sharing vertices with "
                        f"groundside pavement.\n")
            except Exception:
                pass
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
            # TODO(bridges): re-enable bridge / tunnel emission
            # after core pavement geometry is stable + refactored.
            # See ``EMIT_BRIDGES_AND_TUNNELS`` at module top for
            # the full to-do list.  Each gated call must carve its
            # footprint out of overlapping airside / groundside
            # pavement before emitting, or
            # ``test_no_self_overlap`` will fail.
            if EMIT_BRIDGES_AND_TUNNELS:
                # Per user 2026-04-29 (KPHX Sky Harbor Blvd): when
                # a public road passes under any airport bridge,
                # its ENTIRE inside-airport stretch must be at
                # apt_elev − 8 m, not just at the bridge crossings.
                # Run BEFORE _emit_tunnel_portals so the latter can
                # skip OSM ways already depressed here.
                _depressed_way_ids: set = set()
                try:
                    n_dep, _depressed_way_ids = (
                        _emit_through_airport_depressed_roads(
                            layout, _dem, _tile_lat, _tile_lon,
                            xplane_root=xplane_root, icao=icao))
                    if n_dep:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] emitted "
                            f"{n_dep} through-airport depressed "
                            f"road segment(s).\n")
                except Exception:
                    _depressed_way_ids = set()
                # Per user 2026-04-29: re-enable tunnel-portal
                # emission.  For each big-roads tunnel crossing the
                # airport boundary, emit a sloped ramp + 3
                # retaining walls that transition the road from
                # outside-DEM elevation down to airport-elevation
                # − 6 m at the portal.  Subtracts tunnel zones from
                # the boundary ribbon and DEM-bridge polygons so
                # they don't overlap.  OSM way ids handled by the
                # through-airport depressed-road emit are excluded
                # so we don't double-emit on Sky-Harbor-style
                # multi-bridge crossings.
                try:
                    n_tun = _emit_tunnel_portals(
                        layout, _dem, _tile_lat, _tile_lon,
                        excluded_way_ids=_depressed_way_ids)
                    if n_tun:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] emitted "
                            f"{n_tun} tunnel-portal cluster(s) "
                            f"(ramp + walls along approach).\n")
                except Exception:
                    pass
                # Per user 2026-04-29: emit retaining walls along
                # taxi bridges (KBNA Taxiway A, KPHX taxis over
                # Sky Harbor Blvd) and road-following approach
                # shapes that descend from outside-DEM down to
                # apt_elev − 8 m under the bridge.  When the
                # scenery pack already includes 3D bridge OBJs
                # (KBNA), the user wants the road to cut straight
                # through and let the OBJ be the bridge — skip our
                # walls and emit the under-bridge flat polygon.
                # When it doesn't (KPHX), keep the terrain flat
                # for the taxiway and only ramp the road up to the
                # bridge edge — emit walls + skip under-bridge.
                try:
                    _scn_bridge = _scenery_has_bridge_objects(layout)
                except Exception:
                    _scn_bridge = False
                try:
                    n_brg = _emit_taxi_bridges(
                        layout, _dem, _tile_lat, _tile_lon,
                        scenery_has_bridge_objects=_scn_bridge)
                    if n_brg:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] emitted "
                            f"{n_brg} taxi-bridge wall pair(s).\n")
                    elif _scn_bridge:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] {icao}: scenery has "
                            f"3D bridge OBJ(s); skipping wall "
                            f"emission.\n")
                except Exception:
                    pass
                try:
                    n_app = _emit_underpass_road_approaches(
                        layout, _dem, _tile_lat, _tile_lon,
                        scenery_has_bridge_objects=_scn_bridge)
                    if n_app:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] emitted underpass-"
                            f"road approaches for {n_app} "
                            f"surface(s)"
                            f"{' (cut through under bridge OBJ)' if _scn_bridge else ' (ramp up to bridge edge)'}.\n")
                except Exception:
                    pass
        except Exception:
            pass

    return layout


# ══════════════════════════════════════════════════════════════════
# Phase-2: elevations
# ══════════════════════════════════════════════════════════════════


# ──────────────────────────────────────────────────────────────────
# DEM + CIFP + main elevation pipeline
# (re-exported from O4_Pavement_Junctions)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Junctions import (
    _decompose_polygon_with_holes,
    _densify_long_boundary_edges,
    _drop_colinear_boundary_vertices,
    _drop_sliver_corners,
    _merge_thin_decomposed_pieces,
    _polygon_area,
    _polygon_min_thickness,
    _splice_holes,
    _splice_one_hole,
)


# ──────────────────────────────────────────────────────────────────
# DEM + CIFP + main elevation pipeline
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Elevation import (
    APRON_MAX_GRADE,
    DEM_SUFFIX,
    ELEVATION_GRID_STEP_M,
    ELEVATION_SMOOTH_CONVERGE_M,
    ELEVATION_SMOOTH_MAX_ITERS,
    SHARED_VERTEX_CLUSTER_TOL_M,
    TAXI_ANCHOR_DIST_M,
    TAXI_MAX_GRADE,
    _apply_geometric_finalization,
    _compute_elevations,
    _find_cifp_path,
    _load_airport_dem,
    _resample_node_altitudes_nn,
    _sample_dem,
    _solve_pavement_elevations_unified,
)


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


def _emit_groundside_pavement_dem(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        densify_step_m: float = 15.0,
        terminal_gap_m: float = 0.1,
        ) -> int:
    """Emit each saved groundside pavement polygon as a DEM-following
    shape with per-vertex altitudes.

    Per user 2026-04-29: pavement that wraps around the GROUNDSIDE
    of a terminal building (curbside, drop-off, parking) sits at a
    different elevation than the airside apron — at CYXY the
    terminal is cut into the hill so the airside apron is several
    metres lower than the road frontage.  The earlier subtraction
    pass (see ``_terminal_groundside_zone``) keeps the airside
    pavement clean of these strips, but they still belong in the
    output: they should render at local DEM elevation and should
    NOT touch the terminal building footprint (a 0.1 m gap is
    already applied during capture).

    Implementation:
      1. Iterate ``layout._groundside_polys`` (captured during
         ``build_airport_pavement`` immediately before the
         groundside subtraction).
      2. Densify each polygon's exterior to ``densify_step_m`` so
         per-vertex altitudes resolve at the same spatial
         frequency as the boundary ribbon (15 m step → typical
         curbside has 5–10 vertices per side).
      3. Sample DEM at every vertex; emit as ``BuiltShape`` with
         role ``ROLE_GROUNDSIDE_PAVEMENT``, ``node_altitudes`` set,
         and ``altitude``/``altitude_high``/``altitude_low`` left
         None so the OSM emitter writes per-vertex altitude tags.

    Returns the number of polygons emitted.
    """
    polys = list(getattr(layout, "_groundside_polys", []) or [])
    if not polys:
        return 0
    # Build a buffered union of every emitted terminal shape — we
    # subtract this from each groundside polygon so the result
    # leaves a ``terminal_gap_m`` clearance to every actual
    # terminal polygon in the final layout.  Using layout shapes
    # (not OSM source) handles cases where apt.dat row-110 /
    # DSF residue absorption produced a slightly different ring.
    _term_buf = None
    try:
        _t_polys = [s.polygon for s in layout.shapes
                    if s.role == ROLE_TERMINAL
                    and s.polygon is not None
                    and not s.polygon.is_empty]
        if _t_polys:
            _term_buf = unary_union(
                [tp.buffer(terminal_gap_m) for tp in _t_polys])
            if _term_buf.is_empty:
                _term_buf = None
    except Exception:
        _term_buf = None
    # Also subtract every other pavement-bearing layout shape so
    # the groundside pavement never overlaps a rect / junction /
    # apron / runway / terminal / wall / ramp.  The boundary
    # ribbon is excluded — by design it traces over everything.
    NON_OVERLAP_ROLES = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
        ROLE_STUB, ROLE_CROSS_CONNECTOR, ROLE_APRON, ROLE_JUNCTION,
        ROLE_TUNNEL_RAMP, ROLE_RETAINING_WALL,
    }
    _other_buf = None
    try:
        _other_polys = [s.polygon for s in layout.shapes
                        if s.role in NON_OVERLAP_ROLES
                        and s.polygon is not None
                        and not s.polygon.is_empty]
        if _other_polys:
            _other_buf = unary_union(_other_polys)
            if _other_buf.is_empty:
                _other_buf = None
    except Exception:
        _other_buf = None
    cuts = []
    if _term_buf is not None:
        cuts.append(_term_buf)
    if _other_buf is not None:
        cuts.append(_other_buf)
    if cuts:
        try:
            cut_union = unary_union(cuts) if len(cuts) > 1 else cuts[0]
        except Exception:
            cut_union = None
        if cut_union is not None and not cut_union.is_empty:
            clipped: List[Polygon] = []
            for p in polys:
                try:
                    q = p.difference(cut_union)
                except Exception:
                    continue
                if q is None or q.is_empty:
                    continue
                if q.geom_type == "Polygon":
                    if q.area >= 5.0:
                        clipped.append(q)
                elif q.geom_type == "MultiPolygon":
                    for g in q.geoms:
                        if (g.geom_type == "Polygon"
                                and not g.is_empty
                                and g.area >= 5.0):
                            clipped.append(g)
            polys = clipped
    if not polys:
        return 0
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))

    def _dem_at(x: float, y: float) -> Optional[float]:
        try:
            lat, lon = _m_to_ll(x, y)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            return None
    n_emitted = 0
    for p in polys:
        if p is None or p.is_empty or p.geom_type != "Polygon":
            continue
        try:
            ring = list(p.exterior.coords)
        except Exception:
            continue
        if not ring:
            continue
        # Drop trailing repeat to operate on open ring; we'll
        # re-close at the end.
        if ring[0] == ring[-1]:
            ring = ring[:-1]
        if len(ring) < 3:
            continue
        # Densify so per-vertex altitudes resolve well across long
        # straight edges.
        densified: List[Tuple[float, float]] = []
        n_r = len(ring)
        for i in range(n_r):
            ax, ay = ring[i]
            bx, by = ring[(i + 1) % n_r]
            densified.append((ax, ay))
            edge_len = math.hypot(bx - ax, by - ay)
            if edge_len <= densify_step_m:
                continue
            n_intermediate = int(edge_len // densify_step_m)
            for k in range(1, n_intermediate + 1):
                t = (k * densify_step_m) / edge_len
                if t >= 1.0:
                    break
                densified.append((ax + (bx - ax) * t,
                                  ay + (by - ay) * t))
        if len(densified) < 3:
            continue
        # Sample DEM at every densified vertex.  Fall back to the
        # nearest neighbour with a valid sample if a single point
        # lands outside the DEM tile.
        alts: List[Optional[float]] = []
        for x, y in densified:
            alts.append(_dem_at(x, y))
        if all(a is None for a in alts):
            continue
        for k, a in enumerate(alts):
            if a is not None:
                continue
            # Walk outward from k looking for the closest valid
            # sample.
            found: Optional[float] = None
            for off in range(1, len(alts)):
                left = (k - off) % len(alts)
                right = (k + off) % len(alts)
                if alts[left] is not None:
                    found = alts[left]
                    break
                if alts[right] is not None:
                    found = alts[right]
                    break
            alts[k] = found if found is not None else 0.0
        # Rebuild the polygon from densified coords (so it matches
        # the node_altitudes list 1-for-1) and append the closing
        # repeat.
        try:
            new_poly = Polygon(densified)
            if not new_poly.is_valid:
                new_poly = new_poly.buffer(0)
            if (new_poly.geom_type != "Polygon"
                    or new_poly.is_empty):
                continue
        except Exception:
            continue
        # The buffer(0) cleanup may rebuild the ring; re-extract
        # coords and re-sample DEM if the vertex count changed.
        rebuilt = list(new_poly.exterior.coords)
        if rebuilt and rebuilt[0] == rebuilt[-1]:
            rebuilt = rebuilt[:-1]
        if len(rebuilt) != len(densified):
            alts = []
            for x, y in rebuilt:
                a = _dem_at(x, y)
                alts.append(round(float(a), 1) if a is not None
                            else 0.0)
        else:
            alts = [round(float(a), 1) for a in alts]
        # Closing repeat for the OSM emitter's convention.
        node_alts = alts + [alts[0]]
        layout.shapes.append(BuiltShape(
            polygon=new_poly,
            role=ROLE_GROUNDSIDE_PAVEMENT,
            ref="groundside",
            node_altitudes=node_alts))
        n_emitted += 1
    return n_emitted


def _drop_groundside_orphan_junctions(
        layout: "PavementLayout",
        vertex_match_tol_m: float = 0.5,
        ) -> int:
    """Drop junction polygons that connect ONLY to groundside
    pavement (no path through shared vertices to any airside
    rect / runway / terminal).

    Per user 2026-04-29 (CYXY -10111 + -10115): the rect/
    junction tessellator can leave small junction polygons
    sitting next to a groundside polygon when the apt.dat
    row-110 union has a thin strip outside the groundside-
    zone subtraction's perpendicular extent.  Those junctions
    get assigned the airside-flat altitude during the unified
    Laplacian solver (because they were classified as airside
    junctions even though they don't actually touch any
    airside pavement), then they share an edge with the DEM-
    following groundside polygon at a 7 m altitude mismatch —
    X-Plane renders that as a cliff.

    Detection rule:
      A junction is dropped if BOTH:
        1. It shares ≥1 vertex with a
           ``ROLE_GROUNDSIDE_PAVEMENT`` polygon.
        2. It does NOT share any vertex with an airside seed
           shape (runway / primary_parallel / secondary_parallel
           / stub / cross_connector / terminal), whether
           directly or transitively through other junction
           polygons (BFS over junction-junction shared vertices).

      Rationale: a junction between the groundside polygon and
      surrounding terrain is the cliff-creating sliver — drop
      it.  But a junction that legitimately connects an apron
      to a runway or terminal must be kept even if its outer
      edge happens to touch a groundside polygon, otherwise we
      tear a hole in the airside surface (SPJC primary_parallels
      U and M had their short edges become disconnected when
      the simpler rule dropped their connecting junctions).

    Returns the number of junctions dropped.
    """
    AIRSIDE_SEED_ROLES = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL, ROLE_STUB,
        ROLE_CROSS_CONNECTOR, ROLE_TERMINAL,
    }
    bucket_size = vertex_match_tol_m

    def _verts_buckets(s: "BuiltShape") -> List[Tuple[int, int]]:
        if s.polygon is None or s.polygon.is_empty:
            return []
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            return []
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        out = []
        for x, y in coords:
            out.append((int(round(x / bucket_size)),
                        int(round(y / bucket_size))))
        return out
    # Index every junction's vertex buckets (1-bucket halo so
    # near-misses still match neighbours).
    junction_idxs = [i for i, s in enumerate(layout.shapes)
                      if s.role == ROLE_JUNCTION
                      and s.polygon is not None
                      and not s.polygon.is_empty]
    if not junction_idxs:
        return 0
    junction_buckets: Dict[int, set] = {}
    bucket_to_jidx: Dict[Tuple[int, int], List[int]] = {}
    for ji in junction_idxs:
        bs = _verts_buckets(layout.shapes[ji])
        halo: set = set()
        for bx, by in bs:
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    halo.add((bx + dx, by + dy))
        junction_buckets[ji] = halo
        for b in bs:
            bucket_to_jidx.setdefault(b, []).append(ji)
    # Airside seeds (rect / runway / terminal vertex buckets).
    seed_buckets: set = set()
    for s in layout.shapes:
        if s.role not in AIRSIDE_SEED_ROLES:
            continue
        for b in _verts_buckets(s):
            seed_buckets.add(b)
    # Build airside connectivity component over junctions: BFS
    # starting from junctions that share a bucket with any
    # airside seed, propagating through junction-junction
    # shared buckets.
    airside_set: set = set()
    for ji in junction_idxs:
        if junction_buckets[ji] & seed_buckets:
            airside_set.add(ji)
    queue = list(airside_set)
    while queue:
        ji = queue.pop()
        for b in junction_buckets[ji]:
            for kj in bucket_to_jidx.get(b, []):
                if kj in airside_set:
                    continue
                airside_set.add(kj)
                queue.append(kj)
    # Groundside vertex buckets (1-bucket halo).
    gs_buckets: set = set()
    for s in layout.shapes:
        if s.role != ROLE_GROUNDSIDE_PAVEMENT:
            continue
        for bx, by in _verts_buckets(s):
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    gs_buckets.add((bx + dx, by + dy))
    if not gs_buckets:
        return 0
    drop_set: set = set()
    for ji in junction_idxs:
        if ji in airside_set:
            continue
        if junction_buckets[ji] & gs_buckets:
            drop_set.add(ji)
    if not drop_set:
        return 0
    layout.shapes = [s for i, s in enumerate(layout.shapes)
                     if i not in drop_set]
    return len(drop_set)


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


def _emit_tunnel_portals(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        tunnel_depth_m: float = 8.0,
        max_ramp_grade: float = 0.04,
        ramp_min_length_m: float = 200.0,
        arm_max_length_m: float = 500.0,
        carriageway_width_m: float = 22.0,
        retaining_wall_width_m: float = 1.0,
        wall_gap_m: float = 0.5,
        portal_cluster_dist_m: float = 40.0,
        boundary_clearance_m: float = 1.0,
        excluded_way_ids: Optional[set] = None,
        ) -> int:
    """For each tunnel portal (each end of an OSM ``aeroway=*``
    ``tunnel=yes|building_passage`` way), emit the visible road-
    surface structure that transitions outside-DEM elevation down
    to ``apt_elev − tunnel_depth_m`` at the portal:

      1. A flat ``ROLE_RETAINING_WALL`` CAP polygon AT the portal
         node (perpendicular to road direction at portal).  The
         cap's centre line is the portal node — its width spans
         the road (carriageway width + wall_gap on each side),
         its thickness is ``retaining_wall_width_m`` (1 m).
      2. Two flat ``ROLE_RETAINING_WALL`` ARM polygons reaching
         OUTWARD from the cap along the surface roadway, on each
         side of the carriageway.  The arms follow the OSM
         surface road's polyline — multi-segment when the road
         curves.  Arm length adapts to terrain: we first walk up
         to ``arm_max_length_m`` (default 500 m), sample DEM at
         the far end, then truncate the walk so the grade from
         ``apt_elev − tunnel_depth_m`` (portal) up to that DEM
         height never exceeds ``max_ramp_grade``.  Floor at
         ``ramp_min_length_m`` (default 200 m) so a flat road
         still gets a substantial visible approach.
      3. A chain of sloped ``ROLE_TUNNEL_RAMP`` polygons matching
         the surface-road segments under the arms.  Each segment's
         elevation interpolates linearly from
         ``apt_elev − tunnel_depth_m`` at the portal to outside
         DEM at the far end of the walk, in proportion to its
         cumulative distance from the portal.

    Per user 2026-04-29 (SPJC review): the previous geometry put
    the cap PAST the portal INSIDE the airport, with arms going
    AWAY from the tunnel only as far as the OSM way extended
    outside the airport.  That made SPJC's SW tunnel "too short
    and inside the tunnel" because the OSM way ended right at the
    boundary.  Walking the SURFACE road OUTWARD from the OSM
    portal node correctly lands the cap at the tunnel mouth and
    the arms on the highway approach.

    Two-carriageway tunnels (divided highways) cluster by portal-
    node proximity (within ``portal_cluster_dist_m``).  Each
    carriageway in a cluster gets its own arm pair; the caps form
    a perpendicular line across all member portals.

    Boundary coordination: subtract the tunnel-polygon union
    (buffered by ``boundary_clearance_m``, default 1 m which
    exceeds the OSM-emit vertex bucket size of 0.5 m so boundary
    nodes never collapse onto wall nodes) from every
    ``ROLE_BOUNDARY`` shape.

    Returns the number of tunnel PORTALS emitted (each contributing
    1 cap + 2 arm walls + a ramp chain).
    """
    # Load big-roads OSM cache for this tile.
    nodes_r, ways_r = _load_osm_big_roads(
        layout.anchor[0], layout.anchor[1])
    if not ways_r:
        return 0
    # Project nodes to meter space.
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)

    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))
    nodes_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_r.items():
        nodes_m[nid] = _to_m(lon, lat)
    HW_TUNNEL_TYPES = {
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "motorway_link", "trunk_link",
        "primary_link", "residential", "service",
    }
    TUNNEL_VALUES = {"yes", "building_passage"}
    # Build node-to-way and way-by-id indices for surface-road walking.
    way_by_id: Dict[str, Tuple[List[str], Dict[str, str]]] = {}
    node_to_ways: Dict[str, List[str]] = {}
    for wid, nrefs, tags in ways_r:
        way_by_id[wid] = (nrefs, tags)
        for n in nrefs:
            node_to_ways.setdefault(n, []).append(wid)
    # We walk a generous maximum, then truncate per-portal based
    # on actual DEM at the far end so the resulting grade never
    # exceeds ``max_ramp_grade``.  The truncated length is at
    # least ``ramp_min_length_m`` so even flat-road portals get
    # a recognisable approach.  The planning grade is reduced by
    # 0.005 to leave headroom for the 0.1 m altitude rounding —
    # without it, short segments (e.g. 9 m) could round up to
    # ~4.4 % when the design grade is exactly 4 %.
    grade_safety_margin = 0.005
    plan_grade = max(max_ramp_grade - grade_safety_margin, 1e-3)
    arm_walk_max_m = max(arm_max_length_m,
                         ramp_min_length_m,
                         tunnel_depth_m / plan_grade)
    # Helper: airport surface elevation at (cx, cy).  Use the
    # boundary-ribbon ``node_altitudes`` (CIFP-anchored, grade-
    # clamped) when a vertex is nearby, else fall back to DEM.
    def _airport_elevation_at(cx: float, cy: float) -> Optional[float]:
        best_d = float('inf')
        best_alt: Optional[float] = None
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                continue
            if s.ref != "airport_boundary":
                continue
            if not s.node_altitudes:
                continue
            try:
                rcoords = list(s.polygon.exterior.coords)
            except Exception:
                continue
            if rcoords and rcoords[0] == rcoords[-1]:
                rcoords = rcoords[:-1]
            for k, (vx, vy) in enumerate(rcoords):
                if k >= len(s.node_altitudes):
                    break
                d = math.hypot(vx - cx, vy - cy)
                if d < best_d:
                    best_d = d
                    best_alt = s.node_altitudes[k]
        if best_alt is not None and best_d <= 200.0:
            return float(best_alt)
        # Fall back to DEM at the point.
        try:
            lat, lon = _m_to_ll(cx, cy)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            return None
    # Helper: from a portal node, walk the connecting non-tunnel
    # surface road OUTWARD for ``length_m`` metres.  Returns the
    # walked path as a list of (x, y) points starting at the
    # portal, or None if no valid surface way connects.
    def _walk_surface(portal_nid: str,
                      tunnel_wid: str,
                      length_m: float
                      ) -> Optional[List[Tuple[float, float]]]:
        if portal_nid not in nodes_m:
            return None
        # Look at every way that shares this node.  Prefer surface
        # roads (highway tag, not tunnel-tagged); skip the tunnel
        # way itself.  When multiple candidates exist (e.g., a
        # roundabout), pick the first valid one.
        for other_wid in node_to_ways.get(portal_nid, []):
            if other_wid == tunnel_wid:
                continue
            if other_wid not in way_by_id:
                continue
            o_nrefs, o_tags = way_by_id[other_wid]
            if o_tags.get("tunnel") in TUNNEL_VALUES:
                continue
            if o_tags.get("highway") not in HW_TUNNEL_TYPES:
                continue
            try:
                idx = o_nrefs.index(portal_nid)
            except ValueError:
                continue
            # Decide which direction to walk (away from tunnel).
            forward = o_nrefs[idx:]
            backward = list(reversed(o_nrefs[:idx + 1]))
            # If portal is at one end, walk the other way.  If
            # mid-way, walk the longer side.
            walk_refs: List[str]
            if idx == 0:
                walk_refs = forward
            elif idx == len(o_nrefs) - 1:
                walk_refs = backward
            else:
                fl = sum(
                    math.hypot(
                        nodes_m[forward[i + 1]][0]
                        - nodes_m[forward[i]][0],
                        nodes_m[forward[i + 1]][1]
                        - nodes_m[forward[i]][1])
                    for i in range(len(forward) - 1)
                    if (forward[i] in nodes_m
                        and forward[i + 1] in nodes_m))
                bl = sum(
                    math.hypot(
                        nodes_m[backward[i + 1]][0]
                        - nodes_m[backward[i]][0],
                        nodes_m[backward[i + 1]][1]
                        - nodes_m[backward[i]][1])
                    for i in range(len(backward) - 1)
                    if (backward[i] in nodes_m
                        and backward[i + 1] in nodes_m))
                walk_refs = forward if fl >= bl else backward
            # Collect walk points up to length_m.
            pts: List[Tuple[float, float]] = []
            cum = 0.0
            for ni, n in enumerate(walk_refs):
                if n not in nodes_m:
                    break
                p = nodes_m[n]
                if pts:
                    seg_len = math.hypot(
                        p[0] - pts[-1][0], p[1] - pts[-1][1])
                    if cum + seg_len >= length_m:
                        # Truncate the last segment to hit length_m
                        if seg_len > 0:
                            t = (length_m - cum) / seg_len
                            tx = pts[-1][0] + t * (p[0] - pts[-1][0])
                            ty = pts[-1][1] + t * (p[1] - pts[-1][1])
                            pts.append((tx, ty))
                        cum = length_m
                        break
                    cum += seg_len
                pts.append(p)
            if len(pts) >= 2 and cum > 5.0:
                return pts
            # Surface way too short — could chain to next way at
            # the far end, but for simplicity we accept short walks
            # if we got at least 5 m.
            if len(pts) >= 2:
                return pts
        return None
    # Collect portal data: (portal_node_id, tunnel_wid, walk_pts,
    # hw_type, apt_elev_at_portal, dem_at_far_end).
    portal_data: List[Tuple[str, str, List[Tuple[float, float]],
                              str, float, float]] = []
    excluded = excluded_way_ids or set()
    for tw_id, t_nrefs, t_tags in ways_r:
        if t_tags.get("tunnel") not in TUNNEL_VALUES:
            continue
        hw = t_tags.get("highway")
        if hw not in HW_TUNNEL_TYPES:
            continue
        if len(t_nrefs) < 2:
            continue
        # Skip OSM way IDs already handled by the through-
        # airport depressed-road emit (which produces a single
        # uniform depression instead of per-bridge ramps).
        if tw_id in excluded:
            continue
        for portal_idx in (0, len(t_nrefs) - 1):
            portal_nid = t_nrefs[portal_idx]
            if portal_nid not in nodes_m:
                continue
            walk = _walk_surface(portal_nid, tw_id, arm_walk_max_m)
            if walk is None or len(walk) < 2:
                continue
            # Merge very short consecutive segments so altitude
            # rounding to 0.1 m can't push the per-segment grade
            # above ``max_ramp_grade``.  At a 0.05 m worst-case
            # round-up, a 15 m segment rounds to ≤ 0.33 %/error.
            min_segment_m = 15.0
            merged: List[Tuple[float, float]] = [walk[0]]
            for k in range(1, len(walk)):
                d = math.hypot(walk[k][0] - merged[-1][0],
                               walk[k][1] - merged[-1][1])
                if d < min_segment_m and k != len(walk) - 1:
                    continue
                merged.append(walk[k])
            walk = merged
            if len(walk) < 2:
                continue
            portal_xy = walk[0]
            apt_elev = _airport_elevation_at(*portal_xy)
            if apt_elev is None:
                continue
            elev_low = apt_elev - tunnel_depth_m
            # Find the truncation point that keeps grade ≤
            # max_ramp_grade.  Walk along the polyline summing
            # cumulative distance; sample DEM at each vertex; the
            # required-length-so-far is (DEM - elev_low) /
            # max_ramp_grade.  Stop walking once the actual walk
            # length matches or exceeds that requirement (i.e. the
            # grade from portal to here is ≤ max_ramp_grade).
            cum = 0.0
            kept_pts: List[Tuple[float, float]] = [walk[0]]
            grade_ok_at: float = 0.0  # cum dist where grade is OK
            for i in range(1, len(walk)):
                seg_len = math.hypot(
                    walk[i][0] - walk[i - 1][0],
                    walk[i][1] - walk[i - 1][1])
                cum += seg_len
                kept_pts.append(walk[i])
                try:
                    plat, plon = _m_to_ll(*walk[i])
                    dem_h = _sample_dem(
                        dem, tile_lat, tile_lon, plat, plon)
                except Exception:
                    dem_h = None
                if dem_h is None:
                    continue
                drop = float(dem_h) - elev_low
                req = (drop / plan_grade if drop > 0 else 0.0)
                if cum >= req and cum >= ramp_min_length_m:
                    grade_ok_at = cum
                    break
                grade_ok_at = cum
            # Truncate walk to the OK length (or to last vertex
            # reached if we never satisfied the grade — best
            # effort, the resulting ramp will be the gentlest
            # achievable on the available roadway).
            walk = kept_pts
            far_xy = walk[-1]
            try:
                far_lat, far_lon = _m_to_ll(*far_xy)
                far_dem = _sample_dem(
                    dem, tile_lat, tile_lon, far_lat, far_lon)
            except Exception:
                far_dem = None
            if far_dem is None:
                far_dem = apt_elev
            # If the resulting ramp would still violate grade
            # (DEM very high, road too short), cap far_dem at
            # elev_low + grade*length so we still emit something
            # sensible.  The visible top edge then sits LOWER
            # than DEM (subtle terrain dip) — preferable to a
            # spike or a mid-tunnel cap.
            max_drop = plan_grade * grade_ok_at
            if (far_dem - elev_low) > max_drop:
                far_dem = elev_low + max_drop
            portal_data.append(
                (portal_nid, tw_id, walk, hw,
                 float(apt_elev), float(far_dem)))
    if not portal_data:
        return 0
    # Cluster portals by node-coord proximity (divided highways).
    clusters: List[List[int]] = []
    used: set = set()
    for i in range(len(portal_data)):
        if i in used:
            continue
        nid_i = portal_data[i][0]
        cl = [i]
        used.add(i)
        pi = nodes_m[nid_i]
        for j in range(i + 1, len(portal_data)):
            if j in used:
                continue
            pj = nodes_m[portal_data[j][0]]
            if (math.hypot(pi[0] - pj[0], pi[1] - pj[1])
                    < portal_cluster_dist_m):
                cl.append(j)
                used.add(j)
        clusters.append(cl)
    # Per-cluster: build cap + arm walls + ramp chain.
    exclusion_zones: List[Polygon] = []
    n_emitted = 0
    half_carriage = 0.5 * carriageway_width_m
    half_wall_w = retaining_wall_width_m / 2.0
    for cl in clusters:
        # All portals in cluster share approximately the same
        # location.  Use the first portal's walk as the canonical
        # arm path; combine widths for divided highways.
        head = portal_data[cl[0]]
        portal_nid, _wid_unused, walk_pts, _hw, apt_elev, far_dem = head
        if len(walk_pts) < 2:
            continue
        elev_low = apt_elev - tunnel_depth_m
        elev_high = far_dem
        # Compute the walk's cumulative distance for elevation
        # interpolation.
        cum_dists = [0.0]
        for i in range(1, len(walk_pts)):
            cum_dists.append(cum_dists[-1] + math.hypot(
                walk_pts[i][0] - walk_pts[i - 1][0],
                walk_pts[i][1] - walk_pts[i - 1][1]))
        total_walk = cum_dists[-1]
        if total_walk < 5.0:
            continue
        # Cluster spread for combined width: project each cluster
        # member's portal node onto the perpendicular at the head
        # portal.
        first_seg = (walk_pts[1][0] - walk_pts[0][0],
                     walk_pts[1][1] - walk_pts[0][1])
        first_len = math.hypot(*first_seg)
        if first_len < 0.1:
            continue
        first_dir = (first_seg[0] / first_len,
                     first_seg[1] / first_len)
        first_perp = (-first_dir[1], first_dir[0])
        spans = []
        for k in cl:
            ni = portal_data[k][0]
            if ni not in nodes_m:
                continue
            p = nodes_m[ni]
            spans.append(
                (p[0] - walk_pts[0][0]) * first_perp[0]
                + (p[1] - walk_pts[0][1]) * first_perp[1])
        cluster_span = max(spans) - min(spans) if spans else 0.0
        combined_half = half_carriage + 0.5 * cluster_span

        def _build_wall_segment(p_a: Tuple[float, float],
                                 p_b: Tuple[float, float],
                                 perp_off: float
                                 ) -> Optional[Polygon]:
            """4-corner wall polygon parallel to segment ``a-b``,
            offset by ``perp_off`` from the segment's centre line,
            ``retaining_wall_width_m`` thick."""
            seg = (p_b[0] - p_a[0], p_b[1] - p_a[1])
            slen = math.hypot(*seg)
            if slen < 0.1:
                return None
            ux, uy = seg[0] / slen, seg[1] / slen
            nx, ny = -uy, ux
            # Sign of perp_off picks which side.  Inner edge of
            # wall is at perp_off, outer edge at perp_off ± width.
            inner = perp_off
            outer = (perp_off + half_wall_w * 2.0
                     if perp_off >= 0
                     else perp_off - half_wall_w * 2.0)
            corners = [
                (p_a[0] + nx * inner, p_a[1] + ny * inner),
                (p_b[0] + nx * inner, p_b[1] + ny * inner),
                (p_b[0] + nx * outer, p_b[1] + ny * outer),
                (p_a[0] + nx * outer, p_a[1] + ny * outer),
            ]
            try:
                p = Polygon(corners)
                if not p.is_valid:
                    p = p.buffer(0)
                if p.geom_type == "Polygon" and not p.is_empty:
                    return p
            except Exception:
                return None
            return None
        # 1) Cap wall AT the portal node, perpendicular to the
        #    first segment.  The cap's centre line passes through
        #    the portal; its width spans the carriageway + 2 ×
        #    wall_gap; its thickness is retaining_wall_width_m.
        cap_half_len = combined_half + wall_gap_m
        portal_xy = walk_pts[0]
        # Cap polygon: centred at portal, perpendicular to first
        # segment direction, thickness facing INTO the tunnel
        # (opposite first_dir).  We put the cap's outer face
        # right at the portal node so the back of the wall is at
        # OSM's tunnel-portal point.
        c0 = (portal_xy[0] + first_perp[0] * cap_half_len,
              portal_xy[1] + first_perp[1] * cap_half_len)
        c1 = (portal_xy[0] - first_perp[0] * cap_half_len,
              portal_xy[1] - first_perp[1] * cap_half_len)
        # Move cap thickness INTO the tunnel direction (negative
        # first_dir) — cap occupies the strip from portal back
        # by retaining_wall_width_m.
        c0_back = (c0[0] - first_dir[0] * retaining_wall_width_m,
                   c0[1] - first_dir[1] * retaining_wall_width_m)
        c1_back = (c1[0] - first_dir[0] * retaining_wall_width_m,
                   c1[1] - first_dir[1] * retaining_wall_width_m)
        try:
            cap_poly = Polygon([c0, c1, c1_back, c0_back])
            if not cap_poly.is_valid:
                cap_poly = cap_poly.buffer(0)
            if (cap_poly.geom_type == "Polygon"
                    and not cap_poly.is_empty):
                layout.shapes.append(BuiltShape(
                    polygon=cap_poly,
                    role=ROLE_RETAINING_WALL,
                    ref="tunnel_cap",
                    altitude=round(apt_elev, 1)))
                exclusion_zones.append(cap_poly)
        except Exception:
            pass
        # 2) Arm walls + 3) Ramp polygons — one per walk segment.
        # Pre-compute per-vertex offset corners using the bisector
        # of adjacent segments at interior bends.  This makes
        # consecutive segments share their boundary vertices
        # exactly — no overlap, no gap.
        arm_off = combined_half + wall_gap_m + half_wall_w
        n_w = len(walk_pts)
        # Per-vertex perpendicular direction (unit vector).
        verts_perp: List[Tuple[float, float]] = []
        verts_scale: List[float] = []  # extension scale (1/cos(θ/2))
        for i in range(n_w):
            if i == 0:
                s = (walk_pts[1][0] - walk_pts[0][0],
                     walk_pts[1][1] - walk_pts[0][1])
                sl = math.hypot(*s)
                verts_perp.append((-s[1] / sl, s[0] / sl))
                verts_scale.append(1.0)
            elif i == n_w - 1:
                s = (walk_pts[i][0] - walk_pts[i - 1][0],
                     walk_pts[i][1] - walk_pts[i - 1][1])
                sl = math.hypot(*s)
                verts_perp.append((-s[1] / sl, s[0] / sl))
                verts_scale.append(1.0)
            else:
                s1 = (walk_pts[i][0] - walk_pts[i - 1][0],
                      walk_pts[i][1] - walk_pts[i - 1][1])
                s2 = (walk_pts[i + 1][0] - walk_pts[i][0],
                      walk_pts[i + 1][1] - walk_pts[i][1])
                l1 = math.hypot(*s1)
                l2 = math.hypot(*s2)
                u1 = (s1[0] / l1, s1[1] / l1)
                u2 = (s2[0] / l2, s2[1] / l2)
                avg = ((u1[0] + u2[0]) / 2.0,
                       (u1[1] + u2[1]) / 2.0)
                al = math.hypot(*avg)
                if al < 1e-6:
                    # Near-180° doubleback — use first segment's
                    # perpendicular and scale 1.
                    verts_perp.append((-u1[1], u1[0]))
                    verts_scale.append(1.0)
                    continue
                tangent = (avg[0] / al, avg[1] / al)
                perp = (-tangent[1], tangent[0])
                # Adjust offset to compensate for bend angle.
                # cos(θ/2) ≈ sqrt((1 + u1·u2) / 2).
                dot = u1[0] * u2[0] + u1[1] * u2[1]
                cos_half = max(0.1, math.sqrt(
                    max(0.0, (1.0 + dot) / 2.0)))
                verts_perp.append(perp)
                verts_scale.append(1.0 / cos_half)

        def _vertex_offset(idx: int, off: float
                           ) -> Tuple[float, float]:
            px, py = walk_pts[idx]
            nx, ny = verts_perp[idx]
            scaled = off * verts_scale[idx]
            return (px + nx * scaled, py + ny * scaled)

        for i in range(n_w - 1):
            p_a = walk_pts[i]
            p_b = walk_pts[i + 1]
            d_a = cum_dists[i]
            d_b = cum_dists[i + 1]
            seg_len = d_b - d_a
            if seg_len < 0.5:
                continue
            frac_a = d_a / total_walk if total_walk > 0 else 0.0
            frac_b = d_b / total_walk if total_walk > 0 else 0.0
            e_a = (1 - frac_a) * elev_low + frac_a * elev_high
            e_b = (1 - frac_b) * elev_low + frac_b * elev_high
            # Arm walls (one per side).  Inner edge at +/- (arm_off
            # − half_wall_w); outer edge at +/- (arm_off +
            # half_wall_w).  Using the per-vertex bisector so
            # adjacent segments share their join.
            for sign in (+1, -1):
                inner = sign * (arm_off - half_wall_w)
                outer = sign * (arm_off + half_wall_w)
                ai = _vertex_offset(i, inner)
                bi = _vertex_offset(i + 1, inner)
                bo = _vertex_offset(i + 1, outer)
                ao = _vertex_offset(i, outer)
                try:
                    wp = Polygon([ai, bi, bo, ao])
                    if not wp.is_valid:
                        wp = wp.buffer(0)
                    if (wp.geom_type == "Polygon"
                            and not wp.is_empty
                            and wp.area > 0.5):
                        layout.shapes.append(BuiltShape(
                            polygon=wp,
                            role=ROLE_RETAINING_WALL,
                            ref="tunnel_wall",
                            altitude=round(apt_elev, 1)))
                        exclusion_zones.append(wp)
                except Exception:
                    continue
            # Ramp polygon (single segment, sloped).  Corners
            # share with adjacent segments via verts_perp.
            ra = _vertex_offset(i, +combined_half)
            rb = _vertex_offset(i + 1, +combined_half)
            rc = _vertex_offset(i + 1, -combined_half)
            rd = _vertex_offset(i, -combined_half)
            # Rect convention (see _sample_runway_segment_elev):
            # corners [0, 3] are the HIGH-elevation short edge
            # (across the road at the high end), corners [1, 2]
            # are the LOW-elevation short edge.  ra/rb sit on
            # the +side, rd/rc on the -side; ra/rd are at walk[i]
            # and rb/rc at walk[i+1].  Order corners so the slope
            # axis runs ALONG the road (i ↔ i+1), not across it.
            if e_b >= e_a:
                # walk[i+1] = HIGH end → corners 0,3 at i+1
                ramp_corners = [rb, ra, rd, rc]
                eh, el = e_b, e_a
            else:
                # walk[i] = HIGH end → corners 0,3 at i
                ramp_corners = [ra, rb, rc, rd]
                eh, el = e_a, e_b
            try:
                rp = Polygon(ramp_corners)
                if not rp.is_valid:
                    rp = rp.buffer(0)
                if (rp.geom_type == "Polygon"
                        and not rp.is_empty
                        and rp.area > 0.5):
                    if abs(eh - el) >= 0.1:
                        layout.shapes.append(BuiltShape(
                            polygon=rp,
                            role=ROLE_TUNNEL_RAMP,
                            ref="tunnel_ramp",
                            altitude_high=round(eh, 1),
                            altitude_low=round(el, 1)))
                    else:
                        layout.shapes.append(BuiltShape(
                            polygon=rp,
                            role=ROLE_TUNNEL_RAMP,
                            ref="tunnel_ramp",
                            altitude=round(
                                0.5 * (eh + el), 1)))
                    exclusion_zones.append(rp)
            except Exception:
                pass
        n_emitted += 1
    # Boundary coordination: clip every ROLE_BOUNDARY shape so
    # it doesn't overlap the actual tunnel-polygon footprint.
    # The boundary ribbon is built by line-buffering the apt.dat
    # row-130 line, which extends ~2.5 m to either side of the
    # boundary line — so when a tunnel ramp crosses the boundary,
    # the ribbon overlaps both the inside-airport (LOW-end) and
    # outside-airport (HIGH-end) parts of the ramp.  Subtracting
    # the actual ramp+walls union (buffered by 0.5 m for a small
    # visible gap) handles both cases without carving exclusion
    # discs outside the tunnel's footprint — the boundary still
    # traces the rest of the perimeter unchanged.
    if exclusion_zones:
        try:
            tunnel_union = unary_union(exclusion_zones)
        except Exception:
            tunnel_union = None
        if tunnel_union is None or tunnel_union.is_empty:
            return n_emitted
        excl_union = tunnel_union.buffer(boundary_clearance_m)
        kept_shapes: List[BuiltShape] = []
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                kept_shapes.append(s)
                continue
            # Capture the old ring + altitudes BEFORE the clip so
            # we can NN-resample altitudes for the new ring.  Per
            # user 2026-04-29 (HECA crash investigation): leaving
            # the boundary ribbon with no altitudes after a clip
            # made X-Plane crash on load — the 1066-vertex polygon
            # without elevation guidance was unrenderable.
            try:
                _old_ring = list(s.polygon.exterior.coords)
            except Exception:
                _old_ring = []
            if _old_ring and _old_ring[0] == _old_ring[-1]:
                _old_ring = _old_ring[:-1]
            _old_alts = (list(s.node_altitudes)
                          if s.node_altitudes else None)
            try:
                new_poly = s.polygon.difference(excl_union)
            except Exception:
                kept_shapes.append(s)
                continue
            if new_poly.is_empty:
                continue
            if new_poly.geom_type == "Polygon":
                s.polygon = new_poly
                resampled = _resample_node_altitudes_nn(
                    new_poly, _old_ring, _old_alts)
                if resampled is not None:
                    s.node_altitudes = resampled
                # else: keep existing s.altitude / node_altitudes
                # (possibly mismatched but better than nothing).
                kept_shapes.append(s)
            elif new_poly.geom_type == "MultiPolygon":
                # Split the boundary shape into the resulting pieces.
                for g in new_poly.geoms:
                    if (g.geom_type != "Polygon"
                            or g.is_empty
                            or g.area < 5.0):
                        continue
                    resampled = _resample_node_altitudes_nn(
                        g, _old_ring, _old_alts)
                    kept_shapes.append(BuiltShape(
                        polygon=g,
                        role=s.role,
                        ref=s.ref,
                        altitude=(s.altitude if resampled is None
                                  else None),
                        node_altitudes=resampled))
        layout.shapes = kept_shapes
    return n_emitted


def _scenery_has_bridge_objects(
        layout: "PavementLayout",
        bridge_proximity_m: float = 30.0,
        ) -> bool:
    """Return True when the X-Plane scenery pack containing
    ``layout.apt_dat_path`` places at least one taxi-bridge OBJ
    near a bridge taxi rect.

    Detection logic (per user 2026-04-29):
      1. Find the DSF for the scenery pack via
         ``O4_DSF_Reader.find_associated_dsf``.
      2. Convert it to text (cached alongside the .dsf as
         ``.dsf.text``) using DSFTool when not already cached.
      3. Walk OBJECT_DEF lines.  Mark a def as a "bridge def" when
         its path matches ``bridge|elevated|viaduct|overpass``
         AND does NOT match ``sign|signage|trafficsign|wall|
         truss|crane`` — KPHX has 3 ``lib/g10/roadsigns/SignBridge
         *.obj`` defs that are road sign gantries, NOT taxi
         bridges, and the exclude regex filters them out.
      4. Walk OBJECT placement lines.  When a placement uses a
         "bridge def" AND its lat/lon lies within
         ``bridge_proximity_m`` of any taxi rect with
         ``is_bridge=True``, return True.

    Confirmed signal at the test set:
      KBNA  — 23 OBJECT_DEFs match ``Objects/KBNA Bridges/...``
              (KBNA_Bridge_Taxiway-L_p1..p6, KBNA_Crossing_Bridge,
              elevated_edge_twy_B, ...) — many placements within
              the bridge taxi rect → True.
      KPHX  — only ``lib/g10/roadsigns/SignBridge*`` defs which
              the exclude regex drops → False.
    """
    if not layout.apt_dat_path:
        return False
    bridge_rects = [s.polygon for s in layout.shapes
                     if getattr(s, "is_bridge", False)
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not bridge_rects:
        return False
    try:
        import O4_DSF_Reader as _DSFR
    except Exception:
        return False
    dsf_path = _DSFR.find_associated_dsf(
        layout.apt_dat_path,
        layout.anchor[0], layout.anchor[1])
    if dsf_path is None or not os.path.isfile(dsf_path):
        return False
    text_path = dsf_path + ".text"
    needs_convert = (
        not os.path.isfile(text_path)
        or os.path.getmtime(text_path) < os.path.getmtime(dsf_path))
    if needs_convert:
        tool = _DSFR._dsftool_path()
        if tool is None:
            return False
        try:
            import subprocess as _sp
            _sp.run(
                [tool, "--dsf2text", dsf_path, text_path],
                check=True, capture_output=True, timeout=120)
        except Exception:
            return False
    BRIDGE_RE = re.compile(
        r"(?i)bridge|elevated|viaduct|overpass")
    EXCLUDE_RE = re.compile(
        r"(?i)sign|signage|trafficsign|truss|crane")
    bridge_def_idx: set = set()
    object_def_count = 0
    placements: List[Tuple[int, float, float]] = []
    try:
        with open(text_path, "r", encoding="utf-8",
                  errors="replace") as f:
            for line in f:
                if line.startswith("OBJECT_DEF"):
                    parts = line.strip().split(maxsplit=1)
                    path = parts[1] if len(parts) > 1 else ""
                    if (BRIDGE_RE.search(path)
                            and not EXCLUDE_RE.search(path)):
                        bridge_def_idx.add(object_def_count)
                    object_def_count += 1
                elif line.startswith("OBJECT "):
                    tok = line.split()
                    if len(tok) >= 4:
                        try:
                            idx = int(tok[1])
                            lon = float(tok[2])
                            lat = float(tok[3])
                            placements.append((idx, lon, lat))
                        except ValueError:
                            continue
    except Exception:
        return False
    if not bridge_def_idx or not placements:
        return False
    # Project bridge rects to lat/lon for proximity check.  The
    # rects' polygons are in meter space anchored at the layout —
    # convert each placement to meters and test against the
    # buffered rect union.
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _to_m(lon_v: float, lat_v: float) -> Tuple[float, float]:
        return (math.radians(lon_v - lon0) * R * cos0,
                math.radians(lat_v - lat0) * R)
    try:
        bridge_buf = unary_union(bridge_rects).buffer(
            bridge_proximity_m)
    except Exception:
        return False
    for idx, lon_v, lat_v in placements:
        if idx not in bridge_def_idx:
            continue
        x, y = _to_m(lon_v, lat_v)
        if bridge_buf.contains(Point(x, y)):
            return True
    return False


def _emit_taxi_bridges(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        retaining_wall_width_m: float = 1.0,
        wall_gap_m: float = 0.5,
        boundary_clearance_m: float = 1.0,
        scenery_has_bridge_objects: bool = False,
        ) -> int:
    """For each taxi rect marked ``is_bridge=True``, optionally
    emit two flat retaining walls along its long edges at the
    rect's average elevation (the bridge deck altitude).

    Per user 2026-04-29 (KBNA vs KPHX): when the X-Plane scenery
    pack ALREADY contains a 3D taxi-bridge OBJ (detected by
    ``_scenery_has_bridge_objects``, e.g. KBNA's
    ``Objects/KBNA Bridges/KBNA_Bridge_Taxiway-L_*.obj``), the
    scenery's own bridge model is the visible structure — we
    skip the wall emission entirely so we don't double up.  When
    the scenery has NO bridge OBJ (KPHX), our emitted walls
    provide the only visible side-face structure.

    No end-cap walls — the bridge's short edges connect to
    adjacent taxis or junctions at apt_elev (the deck continues
    onto the surrounding airport surface), so an end-cap would
    visually block that join.

    Boundary coordination: when a bridge rect lies inside the
    airport boundary (typical at SPJC / KBNA / KPHX), the walls
    are also inside.  The inside-airport portion of (rect ∪
    walls) gets buffered by ``boundary_clearance_m`` (0.5 m) and
    subtracted from each ``ROLE_BOUNDARY`` shape — same pattern
    as ``_emit_tunnel_portals``.

    Returns the number of bridge rects whose walls were emitted
    (0 when the scenery already has bridge OBJs).
    """
    bridge_shapes = [s for s in layout.shapes
                     if getattr(s, "is_bridge", False)
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not bridge_shapes:
        return 0
    if scenery_has_bridge_objects:
        # The scenery's own 3D bridge OBJs are the visible
        # structure.  Emit nothing here — the deck itself
        # already exists as the taxi rect, and the surrounding
        # mesh is handled by the road-approach helper.
        return 0
    n_emitted = 0
    exclusion_zones: List[Polygon] = []
    for s in bridge_shapes:
        rc = list(s.polygon.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Per ``_rect_from_axis_extended`` corner convention:
        # corners 0,3 = one short edge, corners 1,2 = other.
        # Long edges: corners (0,1) and (2,3).
        # Compute the deck elevation for the wall: average of
        # altitude_high/altitude_low on a sloped rect, or
        # altitude on a flat rect.  Walls match the deck so the
        # join is seamless.
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            deck_elev = 0.5 * (s.altitude_high + s.altitude_low)
        elif s.altitude is not None:
            deck_elev = s.altitude
        else:
            # No elevation set yet — skip.  The post-elevation
            # pass calls _emit_taxi_bridges after rect altitudes
            # are filled in by the unified solver.
            continue
        for (a, b) in ((0, 1), (2, 3)):
            ax, ay = rc[a]
            bx, by = rc[b]
            ex = bx - ax
            ey = by - ay
            elen = math.hypot(ex, ey)
            if elen < 1.0:
                continue
            ux = ex / elen
            uy = ey / elen
            # Outward normal — flip if the test point is INSIDE
            # the rect.
            n_x = -uy
            n_y = ux
            mid_x = 0.5 * (ax + bx)
            mid_y = 0.5 * (ay + by)
            if s.polygon.contains(
                    Point(mid_x + n_x * 0.1,
                          mid_y + n_y * 0.1)):
                n_x = -n_x
                n_y = -n_y
            # Wall sits ``wall_gap_m`` outboard of the rect's
            # long edge, ``retaining_wall_width_m`` thick.
            inner = wall_gap_m
            outer = wall_gap_m + retaining_wall_width_m
            wc = [
                (ax + n_x * inner, ay + n_y * inner),
                (bx + n_x * inner, by + n_y * inner),
                (bx + n_x * outer, by + n_y * outer),
                (ax + n_x * outer, ay + n_y * outer),
            ]
            try:
                wall_poly = Polygon(wc)
                if not wall_poly.is_valid:
                    wall_poly = wall_poly.buffer(0)
                if (wall_poly.geom_type == "Polygon"
                        and not wall_poly.is_empty):
                    layout.shapes.append(BuiltShape(
                        polygon=wall_poly,
                        role=ROLE_RETAINING_WALL,
                        ref="bridge_wall",
                        altitude=round(float(deck_elev), 1)))
                    exclusion_zones.append(wall_poly)
            except Exception:
                continue
        # Track the bridge rect itself so the boundary subtraction
        # below also clears the deck area.
        exclusion_zones.append(s.polygon)
        n_emitted += 1
    # Boundary coordination: subtract the actual (rect ∪ walls)
    # footprint, buffered by 0.5 m, from each ROLE_BOUNDARY shape.
    # Same pattern as ``_emit_tunnel_portals``.
    if exclusion_zones:
        try:
            bridge_union = unary_union(exclusion_zones)
        except Exception:
            bridge_union = None
        if bridge_union is not None and not bridge_union.is_empty:
            try:
                excl = bridge_union.buffer(boundary_clearance_m)
                kept_shapes: List[BuiltShape] = []
                for s in layout.shapes:
                    if s.role != ROLE_BOUNDARY:
                        kept_shapes.append(s)
                        continue
                    try:
                        _old_ring = list(s.polygon.exterior.coords)
                    except Exception:
                        _old_ring = []
                    if (_old_ring
                            and _old_ring[0] == _old_ring[-1]):
                        _old_ring = _old_ring[:-1]
                    _old_alts = (list(s.node_altitudes)
                                  if s.node_altitudes else None)
                    try:
                        new_poly = s.polygon.difference(excl)
                    except Exception:
                        kept_shapes.append(s)
                        continue
                    if new_poly.is_empty:
                        continue
                    if new_poly.geom_type == "Polygon":
                        s.polygon = new_poly
                        resampled = _resample_node_altitudes_nn(
                            new_poly, _old_ring, _old_alts)
                        if resampled is not None:
                            s.node_altitudes = resampled
                        kept_shapes.append(s)
                    elif new_poly.geom_type == "MultiPolygon":
                        for g in new_poly.geoms:
                            if (g.geom_type != "Polygon"
                                    or g.is_empty
                                    or g.area < 5.0):
                                continue
                            resampled = _resample_node_altitudes_nn(
                                g, _old_ring, _old_alts)
                            kept_shapes.append(BuiltShape(
                                polygon=g, role=s.role,
                                ref=s.ref,
                                altitude=(s.altitude
                                          if resampled is None
                                          else None),
                                node_altitudes=resampled))
                layout.shapes = kept_shapes
            except Exception:
                pass
    return n_emitted


def _emit_underpass_road_approaches(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        clearance_depth_m: float = 8.0,
        approach_length_m: float = 80.0,
        road_width_m: float = 22.0,
        ramp_step_m: float = 20.0,
        scenery_has_bridge_objects: bool = False,
        ) -> int:
    """For each underpass case (taxi BRIDGE rect or road TUNNEL
    portal), emit a chain of sloped road-following polygons that
    transition the road surface from outside-DEM elevation down
    to ``apt_elev − clearance_depth_m`` near the underpass and
    back up to DEM after.

    Per user 2026-04-29 (KBNA / KPHX taxi-bridge case): without
    these polygons, the patch's mesh under a bridge sits at
    apt_elev (interpolated from the surrounding airport pavement),
    and OSM-tagged roads rendered at DEM elevation by Ortho4XP's
    road layer get buried.  Emitting a chain of road-following
    polygons forces the mesh under the road to step down to a
    height low enough for the road to clear the bridge underside
    (default 8 m below airport surface), then ramp back up to DEM
    away from the airport.

    Algorithm per underpass surface (bridge rect or tunnel portal
    region):

      1. Find OSM road LineStrings that cross the surface's
         footprint (or pass within 5 m of its long edges, for
         tunnel portals where the road just barely touches).
      2. For each crossing road, find the exterior segments
         immediately approaching and departing the surface.
      3. Walk each approach segment in ``ramp_step_m`` (20 m)
         increments, emitting a sloped 4-corner rect per step
         where elevation interpolates between DEM (start) and
         ``apt_elev − clearance_depth_m`` (end).
      4. Inside the underpass surface itself, emit a flat road-
         following polygon at ``apt_elev − clearance_depth_m``.

    Width of every emitted road polygon: ``road_width_m``
    (22 m by default; matches tunnel-ramp width).

    Returns the number of UNDERPASS surfaces processed.
    """
    # Collect underpass surfaces.
    bridge_shapes = [s for s in layout.shapes
                     if getattr(s, "is_bridge", False)
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not bridge_shapes:
        return 0
    nodes_r, ways_r = _load_osm_big_roads(
        layout.anchor[0], layout.anchor[1])
    if not ways_r:
        return 0
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH
    def _to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)
    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))
    nodes_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_r.items():
        nodes_m[nid] = _to_m(lon, lat)
    HW_TYPES = {
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "motorway_link", "trunk_link",
        "primary_link", "residential", "service",
    }
    # Collect candidate road LineStrings (not bridge or tunnel
    # tagged — those are special cases).  A road that PASSES
    # UNDER a taxi bridge typically isn't tagged bridge=yes
    # itself; only the airport surface is.  But a ROAD that
    # itself bridges over something else IS tagged bridge=yes;
    # we skip those because we don't want to emit road shapes
    # for road-on-road bridges.
    road_lines: List[LineString] = []
    for _wid, nrefs, tags in ways_r:
        if tags.get("highway") not in HW_TYPES:
            continue
        if tags.get("bridge") and tags.get("bridge") != "no":
            continue
        if (tags.get("tunnel")
                and tags.get("tunnel") != "no"):
            continue
        pts = [nodes_m[n] for n in nrefs if n in nodes_m]
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        road_lines.append(ls)
    if not road_lines:
        return 0
    n_processed = 0
    for s in bridge_shapes:
        # Bridge deck elevation.
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            deck_elev = 0.5 * (s.altitude_high + s.altitude_low)
        elif s.altitude is not None:
            deck_elev = s.altitude
        else:
            continue
        low_elev = float(deck_elev) - clearance_depth_m
        # Find road LineStrings that cross the bridge footprint.
        for road_ls in road_lines:
            try:
                inside = road_ls.intersection(s.polygon)
            except Exception:
                continue
            if inside.is_empty:
                continue
            # Pick the longest contiguous piece if multiple.
            if hasattr(inside, "geoms"):
                cand = [g for g in inside.geoms
                        if g.geom_type == "LineString"
                        and g.length > 1.0]
                if not cand:
                    continue
                inside = max(cand, key=lambda g: g.length)
            elif inside.geom_type != "LineString":
                continue
            # Inside-bridge flat polygon at low_elev — only emit
            # when the scenery has a 3D bridge OBJ (KBNA case).
            # In that case the road cuts straight through under
            # the bridge model.  Without a bridge OBJ (KPHX
            # case) the taxi rect itself acts as a flat plate
            # and the road approaches stop at the bridge edge.
            if scenery_has_bridge_objects:
                try:
                    inside_buf = inside.buffer(
                        road_width_m / 2.0,
                        cap_style=2, join_style=2)
                    if (inside_buf.geom_type == "Polygon"
                            and not inside_buf.is_empty):
                        layout.shapes.append(BuiltShape(
                            polygon=inside_buf,
                            role=ROLE_TUNNEL_RAMP,
                            ref="bridge_underpass",
                            altitude=round(low_elev, 1)))
                except Exception:
                    pass
            # Approach + departure ramp chains.  Find the parts
            # of the road OUTSIDE the bridge polygon, then walk
            # each in ramp_step_m steps emitting one sloped rect
            # per step.
            try:
                outside = road_ls.difference(s.polygon)
            except Exception:
                outside = None
            if outside is None or outside.is_empty:
                continue
            outside_pieces = (
                list(outside.geoms)
                if hasattr(outside, "geoms")
                else [outside])
            for piece in outside_pieces:
                if (piece.is_empty
                        or piece.geom_type != "LineString"):
                    continue
                # Decide which end of the piece TOUCHES the bridge.
                pcoords = list(piece.coords)
                if len(pcoords) < 2:
                    continue
                p_start = Point(pcoords[0])
                p_end = Point(pcoords[-1])
                d_start = s.polygon.distance(p_start)
                d_end = s.polygon.distance(p_end)
                if d_start <= d_end:
                    # piece runs FROM bridge edge OUT to far end —
                    # walk it as approach (low → high).
                    # Cap to approach_length_m.
                    walk_len = min(piece.length, approach_length_m)
                    bridge_end_pt = piece.interpolate(0.0)
                    far_end_pt = piece.interpolate(walk_len)
                    walk = LineString([(bridge_end_pt.x,
                                          bridge_end_pt.y),
                                        (far_end_pt.x,
                                          far_end_pt.y)])
                    # Use the actual sub-LineString shape for
                    # better following; but for simplicity use the
                    # straight sub-segment for ramp emission.
                    walk = LineString(pcoords[:])
                    # Trim walk to first ``walk_len`` metres.
                else:
                    walk_len = min(piece.length, approach_length_m)
                    walk = LineString(list(reversed(pcoords)))
                # Step through walk_len in ramp_step_m steps.
                # Each step is a sloped rect of (low_elev → DEM).
                u_prev = 0.0
                while u_prev < walk_len - 1.0:
                    u_next = min(walk_len, u_prev + ramp_step_m)
                    p0 = walk.interpolate(u_prev)
                    p1 = walk.interpolate(u_next)
                    seg_len = math.hypot(p1.x - p0.x, p1.y - p0.y)
                    if seg_len < 1.0:
                        break
                    # Tangent + perpendicular.
                    tx = (p1.x - p0.x) / seg_len
                    ty = (p1.y - p0.y) / seg_len
                    nx = -ty
                    ny = tx
                    half_w = road_width_m / 2.0
                    # Elevation interpolation: u_prev → low_elev
                    # at the bridge, DEM at the far end.
                    frac0 = u_prev / walk_len
                    frac1 = u_next / walk_len
                    try:
                        lat0_p, lon0_p = _m_to_ll(p0.x, p0.y)
                        lat1_p, lon1_p = _m_to_ll(p1.x, p1.y)
                        dem0 = _sample_dem(
                            dem, tile_lat, tile_lon,
                            lat0_p, lon0_p)
                        dem1 = _sample_dem(
                            dem, tile_lat, tile_lon,
                            lat1_p, lon1_p)
                    except Exception:
                        dem0 = dem1 = None
                    if dem0 is None or dem1 is None:
                        break
                    e0 = (1.0 - frac0) * low_elev + frac0 * dem0
                    e1 = (1.0 - frac1) * low_elev + frac1 * dem1
                    corners = [
                        (p0.x + nx * half_w,
                         p0.y + ny * half_w),
                        (p1.x + nx * half_w,
                         p1.y + ny * half_w),
                        (p1.x - nx * half_w,
                         p1.y - ny * half_w),
                        (p0.x - nx * half_w,
                         p0.y - ny * half_w),
                    ]
                    try:
                        seg_poly = Polygon(corners)
                        if not seg_poly.is_valid:
                            seg_poly = seg_poly.buffer(0)
                        if (seg_poly.geom_type == "Polygon"
                                and not seg_poly.is_empty):
                            # corners 0,3 = "high" end vs corners
                            # 1,2 = "low" end depends on which
                            # end is closer to the bridge.  e0
                            # (start = bridge side) is lower.
                            if abs(e0 - e1) >= 0.1:
                                layout.shapes.append(BuiltShape(
                                    polygon=seg_poly,
                                    role=ROLE_TUNNEL_RAMP,
                                    ref="bridge_approach",
                                    altitude_high=round(
                                        max(e0, e1), 1),
                                    altitude_low=round(
                                        min(e0, e1), 1)))
                            else:
                                layout.shapes.append(BuiltShape(
                                    polygon=seg_poly,
                                    role=ROLE_TUNNEL_RAMP,
                                    ref="bridge_approach",
                                    altitude=round(
                                        0.5 * (e0 + e1), 1)))
                    except Exception:
                        pass
                    u_prev = u_next
        n_processed += 1
    return n_processed


def _emit_through_airport_depressed_roads(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        xplane_root: str,
        icao: str,
        depression_depth_m: float = 8.0,
        max_ramp_grade: float = 0.04,
        ramp_min_length_m: float = 200.0,
        arm_max_length_m: float = 500.0,
        road_width_m: float = 22.0,
        retaining_wall_width_m: float = 1.0,
        wall_gap_m: float = 0.5,
        boundary_clearance_m: float = 1.0,
        ) -> Tuple[int, set]:
    """For each public road that ENTERS the airport boundary
    AND passes under a tagged ``aeroway=*, bridge=yes`` way (or
    is connected via OSM-graph node-sharing to a road that does)
    inside the airport, emit:

      1. A flat road-following polygon along the entire inside-
         boundary stretch at ``apt_elev − depression_depth_m``.
         Subsequent OSM road rendering (Ortho4XP's own road
         layer) sits on top of this flat plate.
      2. At each boundary entry/exit point, a ramp polygon
         OUTSIDE the airport that climbs from the depressed
         level back up to the local DEM, capped at
         ``max_ramp_grade`` (default 4 %).

    Per user 2026-04-29 (KPHX Sky Harbor Blvd): when a road
    passes under multiple airport bridges, the road MUST be
    depressed for its entire inside-airport stretch — not just
    at the bridges.  Two bridges spanning the same continuous
    road imply the road can never come back up between them.
    The general solution: any road that enters the airport
    boundary and crosses any aeroway=bridge inside is treated
    this way, plus any road CONNECTED to it via shared OSM
    nodes within the boundary (so on/off ramps inside the
    airport stay coherent with the main road).

    OSM tags consulted:
      * ``aeroway=*, bridge=yes|viaduct`` — bridge LineStrings
        (the airport surface above the depression).
      * ``highway=*`` (motorway, trunk, primary, secondary,
        tertiary, residential, service + their *_link forms) —
        the road network candidates.
      * ``bridge=yes|viaduct`` on a highway — that road is
        ITSELF a bridge over something else (skip; we're
        looking for the road UNDERNEATH).
      * ``tunnel=yes|building_passage`` on a highway — INCLUDED
        as seeds (the airport-bridge case typically tags the
        under-bridge road segment as building_passage in OSM).
        ``_emit_tunnel_portals`` is told to skip every OSM way
        depressed here so we don't double-emit.

    Boundary coordination: the union of every emitted polygon
    is buffered by ``boundary_clearance_m`` and subtracted from
    each ``ROLE_BOUNDARY`` shape (same pattern as the tunnel
    and taxi-bridge emitters), with NN-resampling of per-vertex
    altitudes so the boundary ribbon retains its altitude tags
    after the clip.

    Returns ``(n_emitted, depressed_way_ids)``.  The way-id
    set is intended for ``_emit_tunnel_portals`` to skip — it
    contains every OSM highway way (raw road network) handled
    by this pass, so the explicit-tunnel emitter doesn't
    double-process the same building_passage segments.
    """
    if (layout.airport_boundary is None
            or layout.airport_boundary.is_empty):
        return (0, set())
    boundary = layout.airport_boundary
    # Build a slightly contracted boundary for inside-vs-outside
    # tests so a road point exactly ON the boundary doesn't bounce
    # between true / false on numeric jitter.
    try:
        boundary_strict = boundary.buffer(-0.5)
        if boundary_strict.is_empty:
            boundary_strict = boundary
    except Exception:
        boundary_strict = boundary

    # Load OSM airport-layer tile (for aeroway=bridge LineStrings).
    try:
        nodes_a, ways_a, _ = _load_osm_airports(
            xplane_root, icao,
            layout.anchor[0], layout.anchor[1])
    except Exception:
        return (0, set())
    if not ways_a:
        return (0, set())
    # Load OSM big_roads (for highway ways).
    nodes_r, ways_r = _load_osm_big_roads(
        layout.anchor[0], layout.anchor[1])
    if not ways_r:
        return (0, set())

    grade_safety_margin = 0.005
    plan_grade = max(max_ramp_grade - grade_safety_margin, 1e-3)
    arm_walk_max_m = max(arm_max_length_m, ramp_min_length_m,
                         depression_depth_m / plan_grade)

    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)

    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))

    def _airport_elevation_at(cx: float, cy: float) -> Optional[float]:
        # Reuse pattern from _emit_tunnel_portals: prefer the
        # boundary-ribbon's per-vertex altitude near the point,
        # fall back to DEM.
        best_d = float('inf')
        best_alt: Optional[float] = None
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                continue
            if s.ref != "airport_boundary":
                continue
            if not s.node_altitudes:
                continue
            try:
                rcoords = list(s.polygon.exterior.coords)
            except Exception:
                continue
            if rcoords and rcoords[0] == rcoords[-1]:
                rcoords = rcoords[:-1]
            for k, (vx, vy) in enumerate(rcoords):
                if k >= len(s.node_altitudes):
                    break
                d = math.hypot(vx - cx, vy - cy)
                if d < best_d:
                    best_d = d
                    best_alt = s.node_altitudes[k]
        if best_alt is not None and best_d <= 400.0:
            return float(best_alt)
        try:
            lat, lon = _m_to_ll(cx, cy)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            return None

    # ── Bridge LineStrings (airport-layer OSM) ─────────────────
    bridge_lines: List[LineString] = []
    nodes_a_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_a.items():
        nodes_a_m[nid] = _to_m(lon, lat)
    for wid, nrefs, tags in ways_a:
        if not tags.get("aeroway"):
            continue
        if tags.get("bridge", "") not in ("yes", "viaduct"):
            continue
        pts = [nodes_a_m[n] for n in nrefs if n in nodes_a_m]
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        bridge_lines.append(ls)
    if not bridge_lines:
        return (0, set())

    # ── Highway candidates (big_roads OSM) ─────────────────────
    HW_TYPES = {
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "motorway_link", "trunk_link",
        "primary_link", "secondary_link", "tertiary_link",
        "residential", "service", "unclassified",
    }
    nodes_r_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_r.items():
        nodes_r_m[nid] = _to_m(lon, lat)
    way_data: List[Tuple[str, LineString, List[str]]] = []
    for wid, nrefs, tags in ways_r:
        if tags.get("highway") not in HW_TYPES:
            continue
        if tags.get("bridge", "") in ("yes", "viaduct"):
            # The road IS the bridge, not what's under — skip.
            continue
        # tunnel=building_passage tagging is INCLUDED.  At KPHX,
        # the under-bridge road segments use this tag; they're
        # exactly the seeds we want.
        pts = [nodes_r_m[n] for n in nrefs if n in nodes_r_m]
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        way_data.append((wid, ls, list(nrefs)))
    if not way_data:
        return (0, set())

    # ── Seed: ways whose inside-boundary section crosses a bridge ──
    BRIDGE_PROXIMITY_M = 5.0
    seed_depressed: set = set()
    inside_geom_by_wid: Dict[str, "BaseGeometry"] = {}
    for wid, ls, _nrefs in way_data:
        try:
            inside = ls.intersection(boundary)
        except Exception:
            continue
        if inside.is_empty:
            continue
        inside_geom_by_wid[wid] = inside
        # Iterate inside segments and check bridge proximity.
        segs = []
        if inside.geom_type == "LineString":
            segs = [inside]
        elif inside.geom_type == "MultiLineString":
            segs = list(inside.geoms)
        for seg in segs:
            if seg.is_empty or seg.length < 1.0:
                continue
            for bls in bridge_lines:
                try:
                    if seg.distance(bls) < BRIDGE_PROXIMITY_M:
                        seed_depressed.add(wid)
                        break
                except Exception:
                    continue
            if wid in seed_depressed:
                break
    if not seed_depressed:
        return (0, set())

    # ── BFS over OSM-graph node-sharing INSIDE the boundary ────
    # On-ramps/off-ramps inside the airport are separate OSM
    # ways; if they connect to a depressed seed at any node
    # INSIDE the boundary, they must be depressed too (otherwise
    # the seed and the connecting way disagree on altitude at
    # their shared node and X-Plane renders a cliff).
    node_to_ways: Dict[str, List[str]] = {}
    way_lookup: Dict[str, Tuple[LineString, List[str]]] = {}
    for wid, ls, nrefs in way_data:
        way_lookup[wid] = (ls, nrefs)
        for n in nrefs:
            node_to_ways.setdefault(n, []).append(wid)
    depressed_set: set = set(seed_depressed)
    queue: List[str] = list(seed_depressed)
    while queue:
        wid = queue.pop()
        ls, nrefs = way_lookup[wid]
        for n in nrefs:
            n_xy = nodes_r_m.get(n)
            if n_xy is None:
                continue
            try:
                if not boundary_strict.contains(Point(n_xy)):
                    continue
            except Exception:
                continue
            for other_wid in node_to_ways.get(n, []):
                if other_wid in depressed_set:
                    continue
                # Only propagate if the other way also has an
                # inside-boundary portion (otherwise it's just a
                # surface road glancing the boundary node).
                _o_ls, _o_nrefs = way_lookup[other_wid]
                try:
                    if _o_ls.intersection(boundary).is_empty:
                        continue
                except Exception:
                    continue
                depressed_set.add(other_wid)
                queue.append(other_wid)

    # ── Emit one set of polygons per depressed way ─────────────
    n_emitted = 0
    exclusion_zones: List[Polygon] = []
    half_w = road_width_m / 2.0

    def _smooth_walk(pts: List[Tuple[float, float]],
                      min_segment_m: float = 15.0
                      ) -> List[Tuple[float, float]]:
        """Drop near-colinear / closely-spaced intermediate
        vertices so altitude rounding can't push per-segment
        grade above the design limit."""
        if len(pts) < 3:
            return list(pts)
        merged: List[Tuple[float, float]] = [pts[0]]
        for k in range(1, len(pts)):
            d = math.hypot(pts[k][0] - merged[-1][0],
                           pts[k][1] - merged[-1][1])
            if d < min_segment_m and k != len(pts) - 1:
                continue
            merged.append(pts[k])
        return merged

    for wid in sorted(depressed_set):
        ls, nrefs = way_lookup[wid]
        # 1) Inside-boundary flat plate(s).
        try:
            inside = ls.intersection(boundary)
        except Exception:
            continue
        if inside.is_empty:
            continue
        inside_segs = []
        if inside.geom_type == "LineString":
            inside_segs = [inside]
        elif inside.geom_type == "MultiLineString":
            inside_segs = list(inside.geoms)
        for seg in inside_segs:
            if seg.is_empty or seg.length < 5.0:
                continue
            ctr = seg.centroid
            apt_elev = _airport_elevation_at(ctr.x, ctr.y)
            if apt_elev is None:
                continue
            elev_low = apt_elev - depression_depth_m
            try:
                flat_poly = seg.buffer(
                    half_w, cap_style=2, join_style=2)
                if not flat_poly.is_valid:
                    flat_poly = flat_poly.buffer(0)
            except Exception:
                continue
            if (flat_poly.is_empty
                    or flat_poly.geom_type != "Polygon"):
                continue
            layout.shapes.append(BuiltShape(
                polygon=flat_poly,
                role=ROLE_TUNNEL_RAMP,
                ref="depressed_road",
                altitude=round(elev_low, 1)))
            exclusion_zones.append(flat_poly)
            n_emitted += 1

        # 2) Outside-boundary ramp(s) — one per side of the
        #    boundary the way crosses.  Take the OSM polyline
        #    OUTSIDE the boundary; for each outside piece, walk
        #    OUTWARD from the boundary edge, truncate where the
        #    grade requirement is satisfied, then emit ramp
        #    polygons with bisector vertex sharing at bends.
        try:
            outside = ls.difference(boundary)
        except Exception:
            outside = None
        if outside is None or outside.is_empty:
            continue
        outside_pieces = ([outside]
                          if outside.geom_type == "LineString"
                          else (list(outside.geoms)
                                if outside.geom_type == "MultiLineString"
                                else []))
        for piece in outside_pieces:
            if piece.is_empty or piece.length < 5.0:
                continue
            coords = list(piece.coords)
            # Order coords so coords[0] is at the boundary,
            # coords[-1] is far outside.  The boundary edge is
            # the closest of the two endpoints to the boundary
            # exterior (distance 0 vs distance > 0).
            try:
                d_start = boundary.exterior.distance(
                    Point(coords[0]))
            except Exception:
                d_start = float('inf')
            try:
                d_end = boundary.exterior.distance(
                    Point(coords[-1]))
            except Exception:
                d_end = float('inf')
            if d_end < d_start:
                coords = list(reversed(coords))
            walk = _smooth_walk(coords)
            if len(walk) < 2:
                continue
            # Truncate the walk to a length whose grade keeps
            # ≤ plan_grade given DEM at the far end.
            apt_elev = _airport_elevation_at(*walk[0])
            if apt_elev is None:
                continue
            elev_low = apt_elev - depression_depth_m
            cum = 0.0
            kept_pts: List[Tuple[float, float]] = [walk[0]]
            grade_ok_at: float = 0.0
            for i in range(1, len(walk)):
                seg_len = math.hypot(
                    walk[i][0] - walk[i - 1][0],
                    walk[i][1] - walk[i - 1][1])
                cum += seg_len
                kept_pts.append(walk[i])
                if cum > arm_walk_max_m:
                    grade_ok_at = cum
                    break
                try:
                    plat, plon = _m_to_ll(*walk[i])
                    dem_h = _sample_dem(
                        dem, tile_lat, tile_lon, plat, plon)
                except Exception:
                    dem_h = None
                if dem_h is None:
                    continue
                drop = float(dem_h) - elev_low
                req = (drop / plan_grade if drop > 0 else 0.0)
                if cum >= req and cum >= ramp_min_length_m:
                    grade_ok_at = cum
                    break
                grade_ok_at = cum
            walk = kept_pts
            if len(walk) < 2:
                continue
            far_xy = walk[-1]
            try:
                far_lat, far_lon = _m_to_ll(*far_xy)
                far_dem = _sample_dem(
                    dem, tile_lat, tile_lon, far_lat, far_lon)
            except Exception:
                far_dem = None
            if far_dem is None:
                far_dem = apt_elev
            max_drop = plan_grade * grade_ok_at
            if (far_dem - elev_low) > max_drop:
                far_dem = elev_low + max_drop
            elev_high = far_dem
            cum_dists = [0.0]
            for i in range(1, len(walk)):
                cum_dists.append(cum_dists[-1] + math.hypot(
                    walk[i][0] - walk[i - 1][0],
                    walk[i][1] - walk[i - 1][1]))
            total_walk = cum_dists[-1]
            if total_walk < 5.0:
                continue
            # Per-vertex bisector perpendicular for shared corners
            # at bends (same pattern as _emit_tunnel_portals).
            n_w = len(walk)
            verts_perp: List[Tuple[float, float]] = []
            verts_scale: List[float] = []
            for i in range(n_w):
                if i == 0:
                    s = (walk[1][0] - walk[0][0],
                         walk[1][1] - walk[0][1])
                    sl = math.hypot(*s) or 1e-6
                    verts_perp.append((-s[1] / sl, s[0] / sl))
                    verts_scale.append(1.0)
                elif i == n_w - 1:
                    s = (walk[i][0] - walk[i - 1][0],
                         walk[i][1] - walk[i - 1][1])
                    sl = math.hypot(*s) or 1e-6
                    verts_perp.append((-s[1] / sl, s[0] / sl))
                    verts_scale.append(1.0)
                else:
                    s1 = (walk[i][0] - walk[i - 1][0],
                          walk[i][1] - walk[i - 1][1])
                    s2 = (walk[i + 1][0] - walk[i][0],
                          walk[i + 1][1] - walk[i][1])
                    l1 = math.hypot(*s1) or 1e-6
                    l2 = math.hypot(*s2) or 1e-6
                    u1 = (s1[0] / l1, s1[1] / l1)
                    u2 = (s2[0] / l2, s2[1] / l2)
                    avg = ((u1[0] + u2[0]) / 2.0,
                           (u1[1] + u2[1]) / 2.0)
                    al = math.hypot(*avg)
                    if al < 1e-6:
                        verts_perp.append((-u1[1], u1[0]))
                        verts_scale.append(1.0)
                        continue
                    tangent = (avg[0] / al, avg[1] / al)
                    perp = (-tangent[1], tangent[0])
                    dot = u1[0] * u2[0] + u1[1] * u2[1]
                    cos_half = max(0.1, math.sqrt(
                        max(0.0, (1.0 + dot) / 2.0)))
                    verts_perp.append(perp)
                    verts_scale.append(1.0 / cos_half)

            def _vertex_offset(idx: int, off: float
                               ) -> Tuple[float, float]:
                px, py = walk[idx]
                nx, ny = verts_perp[idx]
                scaled = off * verts_scale[idx]
                return (px + nx * scaled, py + ny * scaled)

            for i in range(n_w - 1):
                d_a = cum_dists[i]
                d_b = cum_dists[i + 1]
                seg_len = d_b - d_a
                if seg_len < 0.5:
                    continue
                frac_a = d_a / total_walk if total_walk > 0 else 0.0
                frac_b = d_b / total_walk if total_walk > 0 else 0.0
                e_a = (1 - frac_a) * elev_low + frac_a * elev_high
                e_b = (1 - frac_b) * elev_low + frac_b * elev_high
                ra = _vertex_offset(i, +half_w)
                rb = _vertex_offset(i + 1, +half_w)
                rc = _vertex_offset(i + 1, -half_w)
                rd = _vertex_offset(i, -half_w)
                # Corners [0, 3] = HIGH end, [1, 2] = LOW end.
                if e_b >= e_a:
                    ramp_corners = [rb, ra, rd, rc]
                    eh, el = e_b, e_a
                else:
                    ramp_corners = [ra, rb, rc, rd]
                    eh, el = e_a, e_b
                try:
                    rp = Polygon(ramp_corners)
                    if not rp.is_valid:
                        rp = rp.buffer(0)
                    if (rp.geom_type != "Polygon"
                            or rp.is_empty
                            or rp.area < 0.5):
                        continue
                except Exception:
                    continue
                if abs(eh - el) >= 0.1:
                    layout.shapes.append(BuiltShape(
                        polygon=rp,
                        role=ROLE_TUNNEL_RAMP,
                        ref="depressed_approach",
                        altitude_high=round(eh, 1),
                        altitude_low=round(el, 1)))
                else:
                    layout.shapes.append(BuiltShape(
                        polygon=rp,
                        role=ROLE_TUNNEL_RAMP,
                        ref="depressed_approach",
                        altitude=round(0.5 * (eh + el), 1)))
                exclusion_zones.append(rp)

    # ── Boundary coordination ─────────────────────────────────
    if exclusion_zones:
        try:
            depressed_union = unary_union(exclusion_zones)
        except Exception:
            depressed_union = None
        if (depressed_union is None
                or depressed_union.is_empty):
            return (n_emitted, depressed_set)
        excl_union = depressed_union.buffer(boundary_clearance_m)
        kept_shapes: List[BuiltShape] = []
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                kept_shapes.append(s)
                continue
            try:
                _old_ring = list(s.polygon.exterior.coords)
            except Exception:
                _old_ring = []
            if _old_ring and _old_ring[0] == _old_ring[-1]:
                _old_ring = _old_ring[:-1]
            _old_alts = (list(s.node_altitudes)
                          if s.node_altitudes else None)
            try:
                new_poly = s.polygon.difference(excl_union)
            except Exception:
                kept_shapes.append(s)
                continue
            if new_poly.is_empty:
                continue
            if new_poly.geom_type == "Polygon":
                s.polygon = new_poly
                resampled = _resample_node_altitudes_nn(
                    new_poly, _old_ring, _old_alts)
                if resampled is not None:
                    s.node_altitudes = resampled
                kept_shapes.append(s)
            elif new_poly.geom_type == "MultiPolygon":
                for g in new_poly.geoms:
                    if (g.geom_type != "Polygon"
                            or g.is_empty
                            or g.area < 5.0):
                        continue
                    resampled = _resample_node_altitudes_nn(
                        g, _old_ring, _old_alts)
                    kept_shapes.append(BuiltShape(
                        polygon=g,
                        role=s.role,
                        ref=s.ref,
                        altitude=(s.altitude
                                  if resampled is None
                                  else None),
                        node_altitudes=resampled))
        layout.shapes = kept_shapes
    return (n_emitted, depressed_set)



# ──────────────────────────────────────────────────────────────────
# Per-shape elevation field + altitude reconciliation
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Elevation import (
    SHARED_AGREE_TOL_M,
    USE_PER_POLYGON_ELEVATION_FIELD,
    _enforce_shared_vertex_altitudes,
    _latlon_to_m_local,
    _match_elev,
    _orient_rect_for_altitude,
    _planar_fit,
    _planar_fit_residuals,
    _re_emit_apron_merged_runway_segments,
    _rederive_terminal_altitude_from_apron_neighbours,
    _smooth_polygon_grid,
    _smooth_within_junction_adjacent_pair_grade,
    _snap_junction_altitudes_to_rect_corners,
)



# ──────────────────────────────────────────────────────────────────
# Elevation finalization (corner buckets, clamp, sliver, overlap)
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Elevation import (
    SUBDIVIDE_MAX_PAIR_DIST_M,
    SUBDIVIDE_SNAP_RADIUS_M,
    _build_clamp_geom_state,
    _clamp_junction_free_vertices,
    _corner_elev_map,
    _corner_elevation_bucket,
    _drop_overlap_against_fixed_shapes,
    _merge_sliver_junctions_into_neighbours,
    _report_within_shape_violations,
    _subdivide_violating_junctions,
    _triangulate_junctions,
)







# ──────────────────────────────────────────────────────────────────
# OSM terminal pad extraction (re-exported from O4_Pavement_Terminals)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Terminals import (
    _build_osm_aeroway_footprint,
    _extract_osm_terminals,
    _terminal_groundside_zone,
    _terminal_pad_from_building,
)



# ──────────────────────────────────────────────────────────────────
# Junction-polygon construction
# (re-exported from O4_Pavement_Junctions)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Junctions import (
    _build_junction_constructive,
    _build_junction_polys_from_corners,
    _build_junctions_from_rect_endpoints,
    _find_junction_points,
    _rect_end_corners,
)





# Centerline-slice constants moved to O4_Pavement_Config (MIN_SEGMENT_LEN_M)
# and to O4_Pavement_Centerlines (RDP_SIMPLIFY_TOL_M, SIGNIFICANT_BEND_DEG,
# BEND_CLUSTER_M, GAP_BRIDGE_MAX_M) by slice 3e.
CLOSE_INTERSECTION_M = 200.0  # dead — kept until next cleanup pass
STUB_MAX_LEN_M = 250.0        # dead — kept until next cleanup pass





# ──────────────────────────────────────────────────────────────────
# Same-ref polyline bridging
# (re-exported from O4_Pavement_Centerlines)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Centerlines import _bridge_same_ref_polylines



# ──────────────────────────────────────────────────────────────────
# Runway-end primary-parallel stub emission
# (re-exported from O4_Pavement_Stubs)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Stubs import _emit_primary_parallel_runway_stubs





# ──────────────────────────────────────────────────────────────────
# OSM aeroway centerline extraction + splitting
# (re-exported from O4_Pavement_Centerlines)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Centerlines import (
    _extract_osm_taxi_centerlines,
    _insert_points_on_boundary,
    _insert_points_on_ring,
    _split_by_width_profile,
    _split_centerlines_at_points,
    _sub_ref_narrow_corridor,
)




# ──────────────────────────────────────────────────────────────────
# Taxi rect construction (re-exported from O4_Pavement_Rects)
# ──────────────────────────────────────────────────────────────────
from O4_Pavement_Rects import (
    EDGE_SNAP_RADIUS_M,
    VERTEX_SNAP_RADIUS_M,
    _axis_to_nearest_rwy_db,
    _build_taxi_rects,
    _cap_rect_length_to_width,
    _classify_role,
    _extend_rect_corners_perpendicular,
    _merge_collinear_rects,
    _merge_collinear_rects_principled,
    _natural_half_width,
    _probe_axis_width,
    _rect_from_axis_extended,
    _refine_roles,
    _snap_corners_to_pavement,
    _trim_to_narrow,
)
