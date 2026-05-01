"""OSM cache + apt.dat selection helpers used by the pipeline.

These functions sit between the orchestrator (``driver``) and the
per-airport build (``pipeline.build_airport_pavement``):

* ``_load_osm_tile``: parse one Ortho4XP-cached OSM tile via
  ``O4_OSM_Utils.OSM_layer.update_dicosm``.  Returns the tile's
  nodes, ways, and outer-relation members in auto_patch's tuple
  shape.
* ``_load_osm_airports``: load + namespace-merge the 9 adjacent
  airport-layer tiles around a given airport coord.  Auto-downloads
  the natural tile on first use if not cached.
* ``_score_apt_dat_against_osm``: rate how well an apt.dat covers
  the OSM-known apron / taxiway features at a given airport.
* ``_pick_best_apt_dat_against_osm``: choose the highest-scoring
  apt.dat candidate (custom-pack vs global) for an airport.
* ``_load_osm_big_roads``: load the highway-layer OSM tile for the
  airport's bbox (used by groundside / boundary / bridges).

All public names are leading-underscored to keep them
package-internal; ``pipeline.build_airport_pavement`` imports them
via ``from .osm_load import ...``.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

from shapely.geometry import LineString, MultiPolygon, Polygon
from shapely.ops import unary_union

import O4_File_Names as FNAMES
import O4_UI_Utils as UI

from . import apt_dat_reader as APR
from .layout import R_EARTH
from .pavement.runways import _runway_rect_m


__all__ = [
    "_load_osm_airports",
    "_load_osm_big_roads",
    "_load_osm_tile",
    "_pick_best_apt_dat_against_osm",
    "_score_apt_dat_against_osm",
]


def _load_osm_tile(path: str) -> Tuple[Dict[str, Tuple[float, float]],
                                       List[Tuple[str, List[str], Dict[str, str]]],
                                       List[Tuple[str, List[str], Dict[str, str]]]]:
    """Parse an Ortho4XP-cached OSM tile (.osm.bz2 or .osm).

    Delegates to ``O4_OSM_Utils.OSM_layer.update_dicosm`` for the
    actual XML parsing (handles bz2 + plain, encoding edge cases,
    self-closing tags, escaped quotes, etc.) then adapts the layer's
    data structures into the auto_patch tuple shape:

      * nodes: ``{id: (lat, lon)}``  — note auto_patch uses ``(lat, lon)``;
        OSM_layer stores ``(lon, lat)``, so coords are flipped here.
      * ways:  ``[(id, [nd_ref, ...], {tag: val})]``
      * relations: ``[(id, [member_way_ref, ...], {tag: val})]``
        — only outer-role way members are included.

    All IDs are stringified at the boundary because every downstream
    auto_patch caller treats node/way IDs as opaque strings (and
    builds tile-prefixed keys like ``"t+30+031:-42"`` for the
    cross-tile-merge step in ``_load_osm_airports``).
    """
    import O4_OSM_Utils as OSM
    layer = OSM.OSM_layer()
    if not layer.update_dicosm(path):
        return {}, [], []

    nodes: Dict[str, Tuple[float, float]] = {
        str(nid): (lat, lon)
        for nid, (lon, lat) in layer.dicosmn.items()
    }
    way_tags = layer.dicosmtags.get("w", {})
    ways: List[Tuple[str, List[str], Dict[str, str]]] = [
        (str(wid), [str(nid) for nid in nds], way_tags.get(wid, {}))
        for wid, nds in layer.dicosmw.items()
    ]
    rel_tags = layer.dicosmtags.get("r", {})
    relations: List[Tuple[str, List[str], Dict[str, str]]] = []
    for rid, role_dict in layer.dicosmrorig.items():
        outer = role_dict.get("outer", []) if isinstance(role_dict, dict) else []
        relations.append(
            (str(rid), [str(wid) for wid in outer], rel_tags.get(rid, {})))
    return nodes, ways, relations


# `_osm_tile_path` removed — was a hand-rolled reimplementation of
# `O4_File_Names.osm_cached(lat, lon, suffix)`.  Call sites now use
# `FNAMES.osm_cached` directly.


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
        return FNAMES.osm_cached(lat_tile, lon_tile, "airports")

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
                UI.vprint(1,
                    f"  [pav-builder] WARN: Overpass download failed "
                    f"for airports tile +{base_lat}{base_lon:+04d}; "
                    f"junctions/rects will be empty.")
        except Exception as exc:
            UI.vprint(1,
                f"  [pav-builder] WARN: airport OSM download error: "
                f"{exc}")
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
            UI.vprint(1,
                f"  [pav-builder] {icao}: apt.dat candidate "
                f"[{label}] apron={ac:.0%} taxi={tc:.0%}  "
                f"{cand}")
        except Exception:
            pass
    if chosen is not None:
        return chosen
    # Nothing qualified — pick whichever has the highest
    # combined coverage to avoid emitting nothing.
    if scores:
        scores.sort(key=lambda s: -(s[1] + s[2]))
        try:
            UI.vprint(1,
                f"  [pav-builder] {icao}: no apt.dat met "
                f"thresholds (apron≥{apron_threshold:.0%}, "
                f"taxi≥{taxi_threshold:.0%}); falling back to "
                f"highest combined coverage: "
                f"apron={scores[0][1]:.0%} taxi={scores[0][2]:.0%}.")
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
            osm_path = FNAMES.osm_cached(
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
