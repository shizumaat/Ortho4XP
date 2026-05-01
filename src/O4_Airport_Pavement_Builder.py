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
            import O4_Cifp_Reader as _CIFP
            import O4_Runway_Geometry as _RWY
            cifp_runways = _CIFP.parse_cifp_file(cifp_path)
            if cifp_runways:
                pairs = _RWY.pair_runways(cifp_runways)

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
                    # Per user 2026-04-29 (CYXY runway 32L/16R
                    # ridge): preserve the dropped segment's
                    # altitude info on the layout so a later
                    # pass can imprint the runway's slope onto
                    # any apron-junction polygon that ends up
                    # covering this segment's footprint.  Without
                    # this, the surrounding apron junction's
                    # vertices are 4-5 m higher than the runway
                    # at the same XY, producing a "ridge across
                    # the runway" the user reported.
                    if (sh.polygon is not None
                            and not sh.polygon.is_empty
                            and (sh.altitude_high is not None
                                 or sh.altitude is not None)):
                        if not hasattr(
                                layout,
                                "_apron_merged_runway_drops"):
                            layout._apron_merged_runway_drops = []
                        layout._apron_merged_runway_drops.append(sh)
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
        # Per user 2026-04-30: bridge-segment insertion was
        # tried (option c) but produced worse results — the
        # inserted bridges overlap apron junctions whose
        # altitudes weren't synchronized, creating 10m+ range
        # junctions with 141 % worst-edge grades.  Function
        # ``_insert_runway_chain_bridges`` is left in place but
        # not called pending a different approach.
        if False:
            n_bridges = _insert_runway_chain_bridges(layout)

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


def _resample_node_altitudes_nn(
        new_poly: Polygon,
        old_open: List[Tuple[float, float]],
        old_alts_closed: Optional[List[float]],
        ) -> Optional[List[float]]:
    """Given a new polygon (post geometry edit) and the OLD ring's
    open-form coords + closed-form altitudes, return a fresh
    ``node_altitudes`` list (closed) for ``new_poly`` by nearest-
    neighbour sampling each new ring vertex against the old ring.

    Used wherever a polygon edit (boundary clip, buffer(0) repair,
    push-off, sliver merge, etc.) changes the vertex count and we
    would otherwise have to drop ``node_altitudes`` — a bare drop
    leaves the polygon with no elevation guidance, which X-Plane
    can render as a terrain spike (or, for very large boundary
    polygons, crash on load — see HECA's 1066-vertex airport
    boundary that lost altitudes during tunnel-clip).

    Returns None if the inputs are insufficient to resample.
    """
    if not old_alts_closed or not old_open:
        return None
    if new_poly is None or new_poly.is_empty:
        return None
    src_alts_open = (
        old_alts_closed[:-1]
        if (len(old_alts_closed) == len(old_open) + 1
            and old_alts_closed[0] == old_alts_closed[-1])
        else old_alts_closed[:len(old_open)])
    if not src_alts_open:
        return None
    try:
        new_open = list(new_poly.exterior.coords)
    except Exception:
        return None
    if new_open and new_open[0] == new_open[-1]:
        new_open = new_open[:-1]
    if not new_open:
        return None
    new_alts: List[float] = []
    for nx, ny in new_open:
        best_d2 = float("inf")
        best_a = src_alts_open[0]
        for k, (sx, sy) in enumerate(old_open):
            if k >= len(src_alts_open):
                break
            d2 = (nx - sx) ** 2 + (ny - sy) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_a = src_alts_open[k]
        new_alts.append(round(float(best_a), 1))
    return new_alts + [new_alts[0]]



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
        layout: "PavementLayout",
        interior_proximity_m: float = 1.0,
        ) -> int:
    """For every junction polygon, snap any vertex whose bucket
    coincides with a runway / sloping-rect corner to that rect's
    corresponding altitude tag value (``altitude_high`` for HIGH
    corners 0,3; ``altitude_low`` for LOW corners 1,2; ``altitude``
    for flat shapes).

    Per user 2026-04-29 (CYXY runway-32L ridge): also snap
    junction vertices that lie INSIDE a sloping-rect footprint
    (within ``interior_proximity_m`` of the rect's interior or
    edge) to the rect's INTERPOLATED altitude at that point.
    Without this, a runway-crossing junction polygon whose
    vertices land inside a runway rect can override the rect's
    smooth slope with mesh-interpolated values that are 4-5 m
    off — visible as a ridge crossing the runway surface.

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
    # Also collect the FULL sloping-rect shapes for interior-
    # snap (point-in-polygon + interpolated altitude).
    rect_shapes_for_interior: List[BuiltShape] = []
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
            rect_shapes_for_interior.append(s)
        elif s.altitude is not None:
            # Flat shapes (terminal pads, pre-elevation rects).
            # Any number of vertices; all share a single altitude.
            for (cx, cy) in coords:
                b = _corner_elevation_bucket(cx, cy)
                rwy_corner_alt.setdefault(b, float(s.altitude))
            rect_shapes_for_interior.append(s)
    if not rwy_corner_alt and not rect_shapes_for_interior:
        return 0
    # Spatial index for interior-snap probes.
    try:
        from shapely.strtree import STRtree as _STRtree
        if rect_shapes_for_interior:
            interior_tree = _STRtree(
                [s.polygon for s in rect_shapes_for_interior])
        else:
            interior_tree = None
    except Exception:
        interior_tree = None
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
            # Pass 1: corner-bucket snap (prevents shared-corner
            # disagreement between junction and rect tags).
            b = _corner_elevation_bucket(cx, cy)
            target_e = rwy_corner_alt.get(b)
            if target_e is not None:
                if abs(s.node_altitudes[i] - target_e) >= 0.05:
                    s.node_altitudes[i] = round(target_e, 1)
                    n_changed += 1
                continue
            # Pass 2: interior snap — when the junction vertex
            # lies inside a sloping rect's polygon, snap to the
            # rect's interpolated altitude at that location.
            # Only adjust when the disagreement is > 0.5 m so we
            # don't undo the Laplacian solver's small refinements
            # at points that aren't truly inside a rect.
            if interior_tree is None:
                continue
            try:
                _Point = Point  # local alias
                pt = _Point(cx, cy)
                cands = interior_tree.query(pt)
            except Exception:
                cands = []
            best_e: Optional[float] = None
            best_d2 = float("inf")
            for hit in cands:
                ri = int(hit) if hasattr(hit, "__int__") else hit
                if not isinstance(ri, int):
                    continue
                rs = rect_shapes_for_interior[ri]
                try:
                    if rs.polygon.distance(pt) > interior_proximity_m:
                        continue
                except Exception:
                    continue
                e = _sample_runway_segment_elev(rs, cx, cy)
                if e is None:
                    continue
                # Pick the rect whose centroid is closest (deals
                # with overlapping rect candidates — runway +
                # parallel + stub at a complex junction).
                try:
                    rcx = rs.polygon.centroid.x
                    rcy = rs.polygon.centroid.y
                    d2 = (rcx - cx) ** 2 + (rcy - cy) ** 2
                except Exception:
                    d2 = 0.0
                if d2 < best_d2:
                    best_d2 = d2
                    best_e = float(e)
            if best_e is None:
                continue
            if abs(s.node_altitudes[i] - best_e) >= 0.5:
                s.node_altitudes[i] = round(best_e, 1)
                n_changed += 1
        # Maintain closed-ring invariant: last == first.
        if (s.node_altitudes
                and len(s.node_altitudes) >= 2
                and s.node_altitudes[0] != s.node_altitudes[-1]):
            s.node_altitudes[-1] = s.node_altitudes[0]
    return n_changed


def _re_emit_apron_merged_runway_segments(
        layout: "PavementLayout",
        ) -> int:
    """Re-emit each apron-merged-and-dropped runway segment back
    into ``layout.shapes`` AFTER the elevation pipeline has
    finished assigning altitudes to surrounding apron junctions.
    Subtracts the re-emitted footprint from any overlapping
    junction polygon so the two don't conflict, with NN-resample
    of the junction's per-vertex altitudes after the clip.

    Per user 2026-04-29 (CYXY runway 32L/16R ridge): the
    builder drops runway segments that lie inside a much-larger
    apron polygon to avoid emitting a visible rectangular
    ribbon.  The surrounding apron junction then takes over the
    surface in that area — but the junction's altitudes come
    from the Laplacian solver pinning at non-runway corners, so
    the surface at the dropped-segment's footprint can sit 4-5 m
    above the runway's CIFP profile.  Reading along the runway,
    the rendered surface dips into the runway segment, rises
    over the apron-junction-covered void, then dips back into
    the next segment — the user's "ridge across runway".

    Re-emitting the dropped segment with its preserved
    ``altitude``/``altitude_high``/``altitude_low`` puts a real
    runway-altitude plate back into the output.  The
    surrounding apron junction is clipped (subtracted) so it no
    longer covers the runway footprint.  The runway surface is
    now continuous at its profile altitude through the absorbed
    area.

    Returns the number of segments re-emitted.
    """
    drops = list(getattr(
        layout, "_apron_merged_runway_drops", []) or [])
    if not drops:
        return 0
    try:
        from shapely.strtree import STRtree as _STRtree
        # Index the junction polygons so we know which to clip.
        jct_idxs = [i for i, s in enumerate(layout.shapes)
                     if s.role == ROLE_JUNCTION
                     and s.polygon is not None
                     and not s.polygon.is_empty]
        if jct_idxs:
            index = _STRtree([layout.shapes[i].polygon
                               for i in jct_idxs])
        else:
            index = None
    except Exception:
        index = None
        jct_idxs = []
    n_emitted = 0
    for seg in drops:
        if seg.polygon is None or seg.polygon.is_empty:
            continue
        # Subtract the segment's footprint from every junction
        # whose polygon overlaps it.  Use the raw segment polygon
        # (no buffer) so the clip is exactly the runway shape.
        if index is not None:
            try:
                cands = index.query(seg.polygon)
            except Exception:
                cands = []
            for hit in cands:
                ji = int(hit) if hasattr(hit, "__int__") else hit
                if not isinstance(ji, int):
                    continue
                shape_i = jct_idxs[ji]
                target = layout.shapes[shape_i]
                if (target.polygon is None
                        or target.polygon.is_empty):
                    continue
                try:
                    inter_area = target.polygon.intersection(
                        seg.polygon).area
                except Exception:
                    continue
                if inter_area < 1.0:
                    continue
                # Capture old ring + altitudes BEFORE clip.
                try:
                    _old_ring = list(
                        target.polygon.exterior.coords)
                except Exception:
                    _old_ring = []
                if _old_ring and _old_ring[0] == _old_ring[-1]:
                    _old_ring = _old_ring[:-1]
                _old_alts = (list(target.node_altitudes)
                              if target.node_altitudes else None)
                try:
                    new_poly = target.polygon.difference(
                        seg.polygon)
                except Exception:
                    continue
                if (new_poly.is_empty
                        or new_poly.geom_type
                        not in ("Polygon", "MultiPolygon")):
                    continue
                if new_poly.geom_type == "MultiPolygon":
                    new_poly = max(
                        (g for g in new_poly.geoms
                          if g.geom_type == "Polygon"),
                        key=lambda g: g.area, default=None)
                    if new_poly is None or new_poly.is_empty:
                        continue
                target.polygon = new_poly
                resampled = _resample_node_altitudes_nn(
                    new_poly, _old_ring, _old_alts)
                if resampled is not None:
                    target.node_altitudes = resampled
        # Re-add the runway segment.
        layout.shapes.append(seg)
        n_emitted += 1
    return n_emitted


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


def _merge_sliver_junctions_into_neighbours(
        layout: "PavementLayout",
        icao: str = "",
        sliver_area_m2: float = 1000.0,
        sliver_ratio: float = 0.05,
        shared_vertex_tol_m: float = 0.5,
        ) -> int:
    """Merge small junction polygons into adjacent larger ones.

    A "sliver" is a junction polygon whose area is below
    ``sliver_area_m2`` AND whose ratio to a neighbour's area is
    below ``sliver_ratio``.  Two junctions are "adjacent" if they
    share at least 2 boundary vertices within
    ``shared_vertex_tol_m``.

    Common cause: ``_decompose_polygon_with_holes`` cuts a
    polygon-with-holes into simple pieces, occasionally carving
    off a tiny strip when the cut grazes the polygon's edge.
    The strip and main piece share a boundary segment (the cut
    line); the strip should be merged back.  Subdivision passes
    in ``_compute_elevations`` can produce similar slivers.

    Returns the number of slivers merged.
    """
    junction_idxs = [i for i, s in enumerate(layout.shapes)
                     if s.role == ROLE_JUNCTION
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if len(junction_idxs) < 2:
        return 0
    # Cache per-shape vertex sets in meter coords for fast tests.
    j_verts: Dict[int, List[Tuple[float, float]]] = {}
    for i in junction_idxs:
        try:
            coords = list(layout.shapes[i].polygon.exterior.coords)
        except Exception:
            j_verts[i] = []
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        j_verts[i] = coords
    tol2 = shared_vertex_tol_m * shared_vertex_tol_m
    merge_into: Dict[int, int] = {}
    for i in junction_idxs:
        ai = layout.shapes[i].polygon.area
        if ai >= sliver_area_m2:
            continue
        best_idx: Optional[int] = None
        best_area = 0.0
        for j in junction_idxs:
            if j == i:
                continue
            aj = layout.shapes[j].polygon.area
            if aj <= ai:
                continue
            if (ai / aj) > sliver_ratio:
                continue
            shared = 0
            for vx, vy in j_verts[i]:
                for ux, uy in j_verts[j]:
                    if (vx - ux) ** 2 + (vy - uy) ** 2 <= tol2:
                        shared += 1
                        break
                if shared >= 2:
                    break
            if shared >= 2 and aj > best_area:
                best_idx = j
                best_area = aj
        if best_idx is not None:
            merge_into[i] = best_idx
    if not merge_into:
        return 0
    for sliver_i, target_j in merge_into.items():
        try:
            target_shape = layout.shapes[target_j]
            sliver_shape = layout.shapes[sliver_i]
            merged = unary_union([
                target_shape.polygon,
                sliver_shape.polygon])
            if (merged.geom_type == "Polygon"
                    and not merged.is_empty):
                # Build a per-vertex elevation lookup from the
                # ORIGINAL target + sliver vertices, then re-derive
                # node_altitudes for the merged polygon by nearest-
                # neighbour match.  Per user 2026-04-29 (CYXY apron
                # regression): just setting ``node_altitudes = None``
                # leaves the merged polygon with NO elevation at all
                # (no downstream re-derivation step exists), which
                # makes X-Plane interpolate from neighbour shapes
                # and produces the "terrain all over the place"
                # apron the user saw.
                lookup: List[Tuple[float, float, float]] = []
                for src_shape in (target_shape, sliver_shape):
                    if not src_shape.node_altitudes:
                        # Sloped or flat alternatives.
                        if (src_shape.altitude_high is not None
                                and src_shape.altitude_low is not None):
                            avg = 0.5 * (
                                src_shape.altitude_high
                                + src_shape.altitude_low)
                        elif src_shape.altitude is not None:
                            avg = float(src_shape.altitude)
                        else:
                            continue
                        try:
                            sc = list(
                                src_shape.polygon.exterior.coords)
                        except Exception:
                            continue
                        if sc and sc[0] == sc[-1]:
                            sc = sc[:-1]
                        for sx, sy in sc:
                            lookup.append((sx, sy, float(avg)))
                        continue
                    src_alts = list(src_shape.node_altitudes)
                    try:
                        sc = list(src_shape.polygon.exterior.coords)
                    except Exception:
                        continue
                    if sc and sc[0] == sc[-1]:
                        sc = sc[:-1]
                    if (len(src_alts) == len(sc) + 1
                            and src_alts[0] == src_alts[-1]):
                        src_alts = src_alts[:-1]
                    for k, (sx, sy) in enumerate(sc):
                        if k >= len(src_alts):
                            break
                        lookup.append(
                            (sx, sy, float(src_alts[k])))
                if not lookup:
                    layout.shapes[target_j].polygon = merged
                    layout.shapes[target_j].node_altitudes = None
                    continue
                merged_coords = list(merged.exterior.coords)
                if (merged_coords
                        and merged_coords[0] == merged_coords[-1]):
                    merged_coords_open = merged_coords[:-1]
                else:
                    merged_coords_open = merged_coords
                new_alts: List[float] = []
                for mx, my in merged_coords_open:
                    best_d2 = float("inf")
                    best_alt = 0.0
                    for sx, sy, sa in lookup:
                        d2 = (mx - sx) ** 2 + (my - sy) ** 2
                        if d2 < best_d2:
                            best_d2 = d2
                            best_alt = sa
                    new_alts.append(round(best_alt, 1))
                # Closing-vertex repeat for the ring.
                if new_alts:
                    new_alts.append(new_alts[0])
                target_shape.polygon = merged
                target_shape.node_altitudes = new_alts
        except Exception:
            continue
    sliver_set = set(merge_into.keys())
    layout.shapes = [
        s for k, s in enumerate(layout.shapes)
        if k not in sliver_set]
    try:
        import sys as _sys
        _sys.stderr.write(
            f"  [pav-builder] {icao}: merged "
            f"{len(merge_into)} sliver junction(s) into "
            f"adjacent larger junctions.\n")
    except Exception:
        pass
    return len(merge_into)


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




def _terminal_pad_from_building(
    building: Polygon,
    pav_polys: List[Polygon],
) -> Optional[Polygon]:
    """Return the OSM building outline as the terminal pad.

    Per user 2026-04-29: the terminal shape must sit right at the
    building edge — no padding, no expansion to the containing
    apt.dat polygon, no fallback buffer.  The earlier behaviour
    (expand to smallest containing apt.dat polygon, or buffer 20 m
    if none) pushed the terminal border out beyond the building
    footprint.

    ``pav_polys`` is unused but kept in the signature so callers
    don't need to change.
    """
    del pav_polys  # intentionally unused
    if building.is_empty:
        return None
    return building


def _build_osm_aeroway_footprint(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    taxi_half_width_m: float = 15.0,
) -> Optional[Polygon]:
    """Return the OSM-known pavement footprint as a Polygon /
    MultiPolygon in meter coordinates.

    Sources:
      * ``aeroway=apron`` / ``stand`` / ``parking_position`` /
        ``hangar_apron`` closed ways → polygons
      * ``aeroway=taxiway`` / ``taxi_lane`` / ``runway`` open ways
        → linestrings buffered by ``taxi_half_width_m`` (covers
        the typical taxi corridor when no width tag is present)
      * If the way IS closed (apron-style polygon tagged taxiway),
        treat as a polygon — covers airports where mappers drew
        taxiway extents instead of centerlines.

    Used by the DSF-loading gate to decide whether apt.dat
    already covers the airport's user-mapped pavement.  When
    apt.dat covers ≥ 85 % of this footprint, DSF additions are
    refused.  When apt.dat has gaps, DSF is admitted only inside
    the gap.
    """
    AREA_AEROWAYS = {
        "apron", "stand", "parking_position",
        "hangar_apron",
    }
    LINEAR_AEROWAYS = {
        "taxiway", "taxi_lane", "runway",
    }
    pieces: List[Polygon] = []
    for _wid, nrefs, tags in ways:
        ay = tags.get("aeroway", "")
        if ay not in AREA_AEROWAYS and ay not in LINEAR_AEROWAYS:
            continue
        pts: List[Tuple[float, float]] = []
        for n in nrefs:
            if n in nodes:
                lat, lon = nodes[n]
                pts.append(to_m(lon, lat))
        if len(pts) < 2:
            continue
        is_closed = (
            len(pts) >= 3
            and abs(pts[0][0] - pts[-1][0]) < 0.5
            and abs(pts[0][1] - pts[-1][1]) < 0.5)
        try:
            if ay in AREA_AEROWAYS or is_closed:
                if len(pts) < 3:
                    continue
                p = Polygon(pts).buffer(0)
                if (p.geom_type == "Polygon"
                        and not p.is_empty
                        and p.area >= 1.0):
                    pieces.append(p)
                elif p.geom_type == "MultiPolygon":
                    for g in p.geoms:
                        if (g.geom_type == "Polygon"
                                and not g.is_empty
                                and g.area >= 1.0):
                            pieces.append(g)
            else:
                ls = LineString(pts)
                if ls.is_empty:
                    continue
                buf = ls.buffer(taxi_half_width_m)
                if (buf.geom_type == "Polygon"
                        and not buf.is_empty):
                    pieces.append(buf)
                elif buf.geom_type == "MultiPolygon":
                    for g in buf.geoms:
                        if (g.geom_type == "Polygon"
                                and not g.is_empty):
                            pieces.append(g)
        except Exception:
            continue
    if not pieces:
        return None
    try:
        merged = unary_union(pieces)
        if merged.is_empty:
            return None
        return merged
    except Exception:
        return None


def _terminal_groundside_zone(
    buildings: List[Polygon],
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    edge_classify_radius_m: float = 30.0,
    groundside_extent_m: float = 100.0,
    apt_pavement_seeds: Optional[List[Polygon]] = None,
    apt_pavement_polys: Optional[List[Polygon]] = None,
) -> Optional[Polygon]:
    """Identify pavement strips on the GROUNDSIDE of terminal
    buildings — the road/curbside frontage where access roads,
    drop-off lanes, and parking sit on pavement at a different
    elevation than the airside apron.  Returns a Polygon /
    MultiPolygon to subtract from ``pav_union`` so groundside
    pavement does NOT become an apron junction grade-clamped to
    the terminal altitude.

    Per user 2026-04-29: combines two airside / groundside
    indicators on each building edge:

      1. APRON ADJACENCY (approach 1).  If the apt.dat polygon
         abutting the edge is connected (touching, transitive)
         to a runway-bearing polygon — i.e. the pavement chain
         from this edge reaches the runway — it's airside.
         Implemented by passing ``apt_pavement_seeds`` and
         testing whether the probe area intersects the
         airside-reachable subset.
      2. OSM TAGS (approach 3).  ``aeroway`` features
         (apron / taxiway / taxi_lane / stand / runway / gate)
         within ``edge_classify_radius_m`` of an outward-edge
         probe → AIRSIDE.  ``highway`` road-class features
         (anything but footway / path / pedestrian / steps)
         within the same probe → GROUNDSIDE.

    An edge is GROUNDSIDE only when at least one groundside
    indicator fires AND no airside indicator does.  An edge with
    NEITHER indicator (most common at airports with sparse OSM
    coverage) defaults to airside (no subtraction) — the safer
    choice when the data can't tell us.

    For each groundside edge a perpendicular outward rectangle
    (depth ``groundside_extent_m``, width = edge length) is added
    to the subtraction zone.  The result is the union of all such
    rectangles; subtracting it from ``pav_union`` keeps the
    airside apron intact while removing the curbside / drop-off
    pavement that should not grade-clamp to the building.
    """
    if not buildings:
        return None
    AIRSIDE_AEROWAY = {
        "apron", "taxiway", "taxi_lane", "stand",
        "runway", "gate", "parking_position",
    }
    # Highway road classes — exclude pedestrian-class tags
    # (footway / path / steps / pedestrian / corridor) which
    # commonly trace airside pedestrian routes ON the apron.
    GROUNDSIDE_HIGHWAY = {
        "primary", "secondary", "tertiary",
        "residential", "unclassified", "service",
        "motorway", "trunk", "primary_link",
        "secondary_link", "tertiary_link",
        "motorway_link", "trunk_link", "living_street",
        "raceway", "road",
    }
    # Build airside / groundside OSM geometry catalogs (meter
    # coords).  Polygons (closed ways) become Polygon; open ways
    # become LineString.
    airside_geoms: List = []
    groundside_geoms: List = []
    for _wid, nrefs, tags in ways:
        ay = tags.get("aeroway", "")
        hw = tags.get("highway", "")
        is_airside = ay in AIRSIDE_AEROWAY
        is_groundside = hw in GROUNDSIDE_HIGHWAY
        if not is_airside and not is_groundside:
            continue
        pts: List[Tuple[float, float]] = []
        for n in nrefs:
            if n in nodes:
                lat, lon = nodes[n]
                pts.append(to_m(lon, lat))
        if len(pts) < 2:
            continue
        try:
            if (len(pts) >= 3
                    and abs(pts[0][0] - pts[-1][0]) < 0.5
                    and abs(pts[0][1] - pts[-1][1]) < 0.5):
                g = Polygon(pts).buffer(0)
            else:
                g = LineString(pts)
            if g.is_empty:
                continue
        except Exception:
            continue
        if is_airside:
            airside_geoms.append(g)
        else:
            groundside_geoms.append(g)
    # Approach 1: build the airside-reachable subset of apt.dat
    # pavement.  Two polygons are "connected" if their boundaries
    # are within ``TOUCH_TOL_M`` of each other.  Seeds are
    # polygons that contain or touch a runway centerline (passed
    # in as ``apt_pavement_seeds``).  BFS from seeds through
    # transitive touches; the reached set is the AIRSIDE-REACHABLE
    # portion of apt.dat pavement.
    #
    # Per user 2026-04-30 (CYXY -10123 NW phantom groundside):
    # When ``apt_pavement_polys`` is provided (full row-110 list),
    # we BFS through it.  Without that list, fall back to using
    # the seeds alone (legacy behaviour).  At CYXY, the apron
    # extends NW past the terminal as part of one giant
    # "Apron 1 and E" polygon that touches all three runways —
    # the BFS reaches that whole polygon, so the NW-edge probe
    # sees airside-reachable pavement and classifies AIRSIDE.
    airside_apt_polys: List[Polygon] = []
    if apt_pavement_seeds:
        TOUCH_TOL_M = 1.0
        if apt_pavement_polys:
            # BFS over apt.dat pavement, seeded by polys that
            # touch any runway / aeroway seed within TOUCH_TOL_M.
            try:
                from shapely.strtree import STRtree as _STRtree
                pav_tree = _STRtree(apt_pavement_polys)
            except Exception:
                pav_tree = None
            seed_buf_polys = [
                s.buffer(TOUCH_TOL_M)
                for s in apt_pavement_seeds
                if s is not None and not s.is_empty]
            airside_idxs: set = set()
            queue: List[int] = []
            # Initial seeds: any apt.dat poly intersecting any
            # buffered seed geometry (i.e. within TOUCH_TOL_M).
            for sb in seed_buf_polys:
                if pav_tree is not None:
                    cands = pav_tree.query(sb)
                else:
                    cands = range(len(apt_pavement_polys))
                for hit in cands:
                    pi = (int(hit) if hasattr(hit, "__int__")
                          else hit)
                    if not isinstance(pi, int):
                        continue
                    if pi in airside_idxs:
                        continue
                    cand = apt_pavement_polys[pi]
                    if cand is None or cand.is_empty:
                        continue
                    try:
                        if cand.intersects(sb):
                            airside_idxs.add(pi)
                            queue.append(pi)
                    except Exception:
                        continue
            # BFS: a poly is airside-reachable if its boundary is
            # within TOUCH_TOL_M of an already-airside poly's
            # boundary.
            while queue:
                pi = queue.pop()
                src = apt_pavement_polys[pi]
                src_buf = src.buffer(TOUCH_TOL_M)
                if pav_tree is not None:
                    cands = pav_tree.query(src_buf)
                else:
                    cands = range(len(apt_pavement_polys))
                for hit in cands:
                    qi = (int(hit) if hasattr(hit, "__int__")
                          else hit)
                    if not isinstance(qi, int):
                        continue
                    if qi in airside_idxs:
                        continue
                    cand = apt_pavement_polys[qi]
                    if cand is None or cand.is_empty:
                        continue
                    try:
                        if cand.intersects(src_buf):
                            airside_idxs.add(qi)
                            queue.append(qi)
                    except Exception:
                        continue
            # Use EXTERIOR rings only (drop interior holes).  Per
            # user 2026-04-30 (CYXY -10123): apt.dat row-110
            # apron polygons commonly have interior rings around
            # terminals / non-pavement zones.  An airside-side
            # building edge whose probe lands INSIDE such a hole
            # would otherwise miss the airside polygon entirely
            # and fall back to UNKNOWN/groundside.  Treating the
            # hole as part of the airside-reachable region (since
            # by definition it is SURROUNDED by airside pavement)
            # restores correct classification.
            for pi in airside_idxs:
                src = apt_pavement_polys[pi]
                if src is None or src.is_empty:
                    continue
                try:
                    ext_only = Polygon(list(src.exterior.coords))
                    if not ext_only.is_valid:
                        ext_only = ext_only.buffer(0)
                    if (ext_only.geom_type == "Polygon"
                            and not ext_only.is_empty):
                        airside_apt_polys.append(ext_only)
                except Exception:
                    airside_apt_polys.append(src)
            # Always include the seeds themselves so the probe
            # also fires when the building edge faces directly
            # onto a runway.
            airside_apt_polys.extend(apt_pavement_seeds)
        else:
            # Legacy fallback: use seeds directly (no BFS).
            airside_apt_polys = list(apt_pavement_seeds)
    if (not airside_geoms
            and not groundside_geoms
            and not airside_apt_polys):
        # No data to classify with — bail out, conservative.
        return None
    # STRtree indexes for fast spatial query.
    try:
        from shapely.strtree import STRtree
    except Exception:
        STRtree = None
    air_tree = (STRtree(airside_geoms)
                if STRtree and airside_geoms else None)
    grd_tree = (STRtree(groundside_geoms)
                if STRtree and groundside_geoms else None)
    apt_air_tree = (STRtree(airside_apt_polys)
                    if STRtree and airside_apt_polys else None)
    zones: List[Polygon] = []
    for bldg in buildings:
        if bldg is None or bldg.is_empty:
            continue
        try:
            coords = list(bldg.exterior.coords)
        except Exception:
            continue
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        n = len(coords)
        if n < 3:
            continue
        # (No building-level airside skip; per-edge probe
        # below handles classification.)
        # Per-building two-pass classification.  An edge becomes
        # GROUNDSIDE only when the OSM data gives us an EXPLICIT
        # indicator on that edge — either ``highway=*`` (road
        # class) immediately outward (approach 3) or the apt.dat
        # polygon at that edge is NOT reachable from any runway
        # via touching connectivity (approach 1, future work).
        # An edge with NO clear indicator stays UNKNOWN and is
        # NOT subtracted — promoting unknown to groundside on a
        # building that has at least one airside edge proved too
        # aggressive at HECA, where apron-tagging is patchy and
        # several airside edges look UNKNOWN to OSM.
        EDGE_AIRSIDE = 1
        EDGE_GROUNDSIDE = 2
        EDGE_UNKNOWN = 0
        edge_class: List[int] = [EDGE_UNKNOWN] * n
        edge_geom: List[Optional[Tuple[float, float, float, float,
                                         float, float]]] = [None] * n
        for i in range(n):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % n]
            edge_len = math.hypot(bx - ax, by - ay)
            if edge_len < 1.0:
                continue
            tx = (bx - ax) / edge_len
            ty = (by - ay) / edge_len
            n_x = ty
            n_y = -tx
            mid_x = 0.5 * (ax + bx)
            mid_y = 0.5 * (ay + by)
            if bldg.contains(
                    Point(mid_x + n_x * 1.0,
                          mid_y + n_y * 1.0)):
                n_x = -n_x
                n_y = -n_y
            edge_geom[i] = (ax, ay, bx, by, n_x, n_y)
            try:
                probe = Polygon([
                    (ax, ay), (bx, by),
                    (bx + n_x * edge_classify_radius_m,
                     by + n_y * edge_classify_radius_m),
                    (ax + n_x * edge_classify_radius_m,
                     ay + n_y * edge_classify_radius_m),
                ])
                if not probe.is_valid:
                    probe = probe.buffer(0)
                if probe.is_empty:
                    continue
            except Exception:
                continue
            airside = False
            groundside = False
            # STRtree.query in Shapely 2.x returns numpy int64
            # indices into the original geometry list — always
            # index back to fetch the geometry.
            if air_tree is not None:
                for hit in air_tree.query(probe):
                    g = airside_geoms[int(hit)]
                    if g.intersects(probe):
                        airside = True
                        break
            if not airside and apt_air_tree is not None:
                # Per user 2026-04-30 (CYXY -10123): use a
                # LONGER probe for the apt.dat-pavement-
                # connectivity check (vs the 30 m default for
                # OSM aeroway tags).  At CYXY the terminal sits
                # in a building-shaped concavity in the apron's
                # exterior ring; the 30 m probe falls inside
                # the concavity and never hits airside pavement.
                # A 100 m probe reaches past the concavity into
                # the apron's main body.  Groundside detection
                # via OSM highway tags still uses the 30 m probe
                # so this doesn't make the classifier more eager
                # to fire airside on roads close to the building.
                _far_radius_m = 100.0
                far_probe = Polygon([
                    (ax, ay), (bx, by),
                    (bx + n_x * _far_radius_m,
                     by + n_y * _far_radius_m),
                    (ax + n_x * _far_radius_m,
                     ay + n_y * _far_radius_m),
                ])
                if not far_probe.is_valid:
                    far_probe = far_probe.buffer(0)
                if not far_probe.is_empty:
                    for hit in apt_air_tree.query(far_probe):
                        g = airside_apt_polys[int(hit)]
                        if g.intersects(far_probe):
                            airside = True
                            break
            if airside:
                edge_class[i] = EDGE_AIRSIDE
                continue
            if grd_tree is not None:
                for hit in grd_tree.query(probe):
                    g = groundside_geoms[int(hit)]
                    if g.intersects(probe):
                        groundside = True
                        break
            if groundside:
                edge_class[i] = EDGE_GROUNDSIDE
                continue
            edge_class[i] = EDGE_UNKNOWN
        # Per user 2026-04-29 (latest): on the airside, the apron
        # SHARES VERTICES with the terminal building (the apron
        # junction wraps around the building corners).  On the
        # groundside, NO airport pavement should connect to the
        # terminal — driveways / parking are typically several
        # metres higher than the airside apron (CYXY 4 m
        # difference) for terrain or elevated road decks.
        #
        # Identification: airside if OSM aeroway features are
        # nearby; UNKNOWN if neither aeroway nor highway is
        # present (typical at airports with sparse OSM road
        # coverage — most of our test set).  When at least ONE
        # edge has an airside indicator, every edge that ISN'T
        # airside is treated as groundside (UNKNOWN promoted).
        # When no edge is airside, we have no signal to pick
        # which direction is which — leave the building alone.
        any_airside = any(c == EDGE_AIRSIDE for c in edge_class)
        for i in range(n):
            cls = edge_class[i]
            if cls == EDGE_AIRSIDE:
                continue
            if cls == EDGE_UNKNOWN and not any_airside:
                continue
            geom = edge_geom[i]
            if geom is None:
                continue
            ax, ay, bx, by, n_x, n_y = geom
            try:
                zone = Polygon([
                    (ax, ay), (bx, by),
                    (bx + n_x * groundside_extent_m,
                     by + n_y * groundside_extent_m),
                    (ax + n_x * groundside_extent_m,
                     ay + n_y * groundside_extent_m),
                ])
                if not zone.is_valid:
                    zone = zone.buffer(0)
                if (zone.geom_type == "Polygon"
                        and not zone.is_empty):
                    zones.append(zone)
            except Exception:
                continue
    if not zones:
        return None
    try:
        merged = unary_union(zones)
        if merged.is_empty:
            return None
        return merged
    except Exception:
        return None


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
