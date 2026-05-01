"""End-to-end orchestration of the airport pavement builder.

Sequences the phase-1 (geometry + role) and phase-2 (elevation)
passes by calling out into the focused ``O4_Pavement_*`` modules:

* OSM tile + airport extraction.
* apt.dat selection (best-of-OSM vs custom-pack heuristics).
* Phase-1 construction → ``O4_Pavement_Rects``,
  ``O4_Pavement_Centerlines``, ``O4_Pavement_Stubs``,
  ``O4_Pavement_Junctions``, ``O4_Pavement_Terminals``.
* Phase-2 elevation → ``O4_Pavement_Elevation``.
* Feature emit → ``O4_Pavement_Boundary``,
  ``O4_Pavement_Groundside``, ``O4_Pavement_Bridges``.
* Output via ``PavementLayout.to_osm`` (in ``O4_Pavement_Layout``).

Public API:

    build_airport_pavement(icao, xplane_root, *, compute_elevations=True)

Backward-compat shim ``O4_Airport_Pavement_Builder`` re-exports
``build_airport_pavement`` (and a few helpers used by other
modules) so existing call sites keep working.
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

from . import apt_dat_reader as APR
from .pavement import classifier as PC
from .pavement import strips as PS


# ──────────────────────────────────────────────────────────────────
# Constants (re-exported from O4_Pavement_Config + O4_Pavement_Layout)
# ──────────────────────────────────────────────────────────────────
from .config import (
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
from .layout import (
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
from .pavement.vertices import (
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
from .pavement.runways import (
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


from .pavement.absorption import (
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
        from . import dsf_reader as _DSFR
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
from .pavement.junctions import (
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
from .elevation import (
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



# ──────────────────────────────────────────────────────────────────
# Airport boundary shape (re-exported from O4_Pavement_Boundary)
# ──────────────────────────────────────────────────────────────────
from .boundary import _emit_airport_boundary_shape



# ──────────────────────────────────────────────────────────────────
# Groundside (curbside / drop-off) pavement
# (re-exported from O4_Pavement_Groundside)
# ──────────────────────────────────────────────────────────────────
from .groundside import (
    _drop_groundside_orphan_junctions,
    _emit_groundside_pavement_dem,
)



# ──────────────────────────────────────────────────────────────────
# Boundary→DEM bridge polygons (re-exported from O4_Pavement_Boundary)
# ──────────────────────────────────────────────────────────────────
from .boundary import _emit_boundary_dem_bridge



# ──────────────────────────────────────────────────────────────────
# Taxi/road bridges + tunnel portals + depressed-road segments
# (re-exported from O4_Pavement_Bridges; gated by EMIT_BRIDGES_AND_TUNNELS)
# ──────────────────────────────────────────────────────────────────
from .bridges import (
    _emit_taxi_bridges,
    _emit_through_airport_depressed_roads,
    _emit_tunnel_portals,
    _emit_underpass_road_approaches,
    _scenery_has_bridge_objects,
)





# ──────────────────────────────────────────────────────────────────
# Per-shape elevation field + altitude reconciliation
# (re-exported from O4_Pavement_Elevation)
# ──────────────────────────────────────────────────────────────────
from .elevation import (
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
from .elevation import (
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
from .terminals import (
    _build_osm_aeroway_footprint,
    _extract_osm_terminals,
    _terminal_groundside_zone,
    _terminal_pad_from_building,
)



# ──────────────────────────────────────────────────────────────────
# Junction-polygon construction
# (re-exported from O4_Pavement_Junctions)
# ──────────────────────────────────────────────────────────────────
from .pavement.junctions import (
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
from .pavement.centerlines import _bridge_same_ref_polylines



# ──────────────────────────────────────────────────────────────────
# Runway-end primary-parallel stub emission
# (re-exported from O4_Pavement_Stubs)
# ──────────────────────────────────────────────────────────────────
from .pavement.stubs import _emit_primary_parallel_runway_stubs





# ──────────────────────────────────────────────────────────────────
# OSM aeroway centerline extraction + splitting
# (re-exported from O4_Pavement_Centerlines)
# ──────────────────────────────────────────────────────────────────
from .pavement.centerlines import (
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
from .pavement.rects import (
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
