"""OSM aeroway extractors: taxiway centerlines, terminal/hangar
buildings, big-roads near airports.

Higher-level extractors built on top of ``O4_OSM_Utils.OSM_layer``,
producing aviation-specific dicts keyed by airport identifier.
Consumed by ``O4_Vector_Map`` during the tile-build vector phase.

Public API:
    extract_taxiway_info(airport_layer, dico_airports, tile)
                                     — taxiway centerlines per airport
    extract_building_info(airport_layer, dico_airports, tile, ...)
                                     — terminal / hangar footprints per airport
    extract_road_info(dico_airports, tile, road_layer=None)
                                     — big-roads near airport boundaries
"""
import O4_UI_Utils as UI

__all__ = [
    "extract_taxiway_info",
    "extract_building_info",
    "extract_road_info",
]


DEFAULT_TAXIWAY_WIDTH = 15.0  # meters - typical taxiway width
MAX_TAXIWAY_GRADE = 0.015     # 1.5% max grade for taxiways


def extract_taxiway_info(airport_layer, dico_airports, tile):
    """Extract taxiway centerline data from the OSM airport layer.

    Pulls taxiway way coordinates from the airport_layer's node/way dicts,
    keyed by the airport identifier used in dico_airports (typically ICAO).

    Args:
        airport_layer: OSM_layer object with dicosmn (nodes) and dicosmw (ways).
        dico_airports: dict from discover_airport_names / attach_surfaces, keyed
                       by airport identifier.
        tile: Tile object with .lat and .lon.

    Returns:
        dict: {airport_key: [{'centerline': [(lon, lat), ...],
                              'wayid': int,
                              'name': str}, ...]}
              Coordinates are ABSOLUTE (not tile-relative).
              ``name`` is the OSM ``ref`` tag (e.g. "A", "A1", "F",
              "V") if present, otherwise an empty string.  Used by
              the centerline-driven Phase C0 emission to group
              centerlines by named taxiway.
    """
    way_tags = getattr(airport_layer, "dicosmtags", {}).get("w", {})
    result = {}
    for airport in dico_airports:
        apt = dico_airports[airport]
        # taxiway is (MultiPolygon_area, wayid_list) after build_taxiway_areas
        if not isinstance(apt.get("taxiway"), tuple) or len(apt["taxiway"]) < 2:
            continue
        wayid_list = apt["taxiway"][1]
        if not wayid_list:
            continue
        taxiways = []
        for wayid in wayid_list:
            if wayid not in airport_layer.dicosmw:
                continue
            node_ids = airport_layer.dicosmw[wayid]
            coords = []
            for nid in node_ids:
                if nid in airport_layer.dicosmn:
                    coords.append(tuple(airport_layer.dicosmn[nid]))
            if len(coords) >= 2:
                tags = way_tags.get(wayid, {})
                taxiways.append({
                    "centerline": coords,
                    "wayid": wayid,
                    "name": tags.get("ref", "") or "",
                })
        if taxiways:
            result[airport] = taxiways
    return result


def extract_building_info(airport_layer, dico_airports, tile,
                          building_layer=None):
    """Extract building footprints within airport boundaries.

    Collects building polygons from three sources:
    1. Hangars already attached to airports (aeroway=hangar from dico_airports)
    2. Terminals from the airport_layer (aeroway=terminal ways)
    3. General buildings (building=*) from a separate building_layer, filtered
       to those intersecting airport boundaries.

    Args:
        airport_layer: OSM_layer with aeroway data (includes terminal ways).
        dico_airports: dict with processed airport data including 'hangar'
                       MultiPolygon and 'boundary' MultiPolygon.
        tile: Tile object with .lat and .lon.
        building_layer: Optional OSM_layer with building=* data. If None,
                        only hangars and terminals are extracted.

    Returns:
        dict: {airport_key: [{'footprint': [(lon, lat), ...], 'source': str}, ...]}
              Coordinates are ABSOLUTE (not tile-relative).
              Each footprint is a list of (lon, lat) ring coordinates.
    """
    result = {}
    for airport in dico_airports:
        apt = dico_airports[airport]
        buildings = []

        # Source 1: Hangars from airport_layer (already attached to airports)
        hangar_wayids = []
        # Before build_hangar_areas, apt["hangar"] is a list of wayids.
        # After build_hangar_areas, it's a MultiPolygon (tile-relative).
        # We need to handle both cases, but typically this runs after
        # build_hangar_areas, so we work with the MultiPolygon.
        hangar_geom = apt.get("hangar")
        if hangar_geom is not None and hasattr(hangar_geom, "geoms"):
            # It's a MultiPolygon in tile-relative coords - convert back
            for poly in hangar_geom.geoms:
                if poly.is_empty or not poly.is_valid or poly.area < 1e-10:
                    continue
                # Convert from tile-relative to absolute
                coords = [
                    (x + tile.lon, y + tile.lat)
                    for x, y in poly.exterior.coords
                ]
                buildings.append({
                    "footprint": coords,
                    "source": "hangar",
                })

        # Source 2: Terminals from airport_layer (aeroway=terminal)
        # Check both simple ways and multipolygon relations, since large
        # terminals are often mapped as relations in OSM.
        boundary = apt.get("boundary")
        if boundary is not None:
            from shapely import geometry as shp_geom

            def _try_add_terminal(coords, source_label):
                """Add terminal footprint if it intersects this airport."""
                if len(coords) < 4:
                    return
                try:
                    term_rel = shp_geom.Polygon([
                        (c[0] - tile.lon, c[1] - tile.lat)
                        for c in coords
                    ])
                    if (term_rel.is_valid and not term_rel.is_empty
                            and boundary.intersects(term_rel)):
                        buildings.append({
                            "footprint": coords,
                            "source": source_label,
                        })
                except Exception:
                    pass

            # 2a: Terminal ways
            for wayid in airport_layer.dicosmw:
                wtags = airport_layer.dicosmtags.get("w", {}).get(wayid, {})
                if wtags.get("aeroway") != "terminal":
                    continue
                node_ids = airport_layer.dicosmw[wayid]
                if len(node_ids) < 4 or node_ids[0] != node_ids[-1]:
                    continue
                coords = []
                for nid in node_ids:
                    if nid in airport_layer.dicosmn:
                        coords.append(tuple(airport_layer.dicosmn[nid]))
                _try_add_terminal(coords, "terminal")

            # 2b: Terminal relations (multipolygons)
            for relid in airport_layer.dicosmr:
                rtags = airport_layer.dicosmtags.get("r", {}).get(relid, {})
                if rtags.get("aeroway") != "terminal":
                    continue
                rel_data = airport_layer.dicosmr[relid]
                outer_rings = rel_data.get("outer", [])
                for ring_nids in outer_rings:
                    coords = []
                    for nid in ring_nids:
                        if nid in airport_layer.dicosmn:
                            coords.append(
                                tuple(airport_layer.dicosmn[nid])
                            )
                    _try_add_terminal(coords, "terminal")

        # Source 3: General buildings from building_layer
        if building_layer is not None and boundary is not None:
            try:
                for wayid in building_layer.dicosmw:
                    if wayid not in building_layer.dicosmtags.get("w", {}):
                        continue
                    wtags = building_layer.dicosmtags["w"][wayid]
                    if "building" not in wtags:
                        continue
                    node_ids = building_layer.dicosmw[wayid]
                    if len(node_ids) < 4:
                        continue
                    if node_ids[0] != node_ids[-1]:
                        continue
                    coords = []
                    for nid in node_ids:
                        if nid in building_layer.dicosmn:
                            coords.append(tuple(building_layer.dicosmn[nid]))
                    if len(coords) < 4:
                        continue
                    # Coords are absolute (lon, lat); boundary is
                    # tile-relative, so convert building to tile-relative
                    try:
                        from shapely import geometry as shp_geom
                        bldg_rel = shp_geom.Polygon([
                            (c[0] - tile.lon, c[1] - tile.lat)
                            for c in coords
                        ])
                        if (bldg_rel.is_valid
                                and not bldg_rel.is_empty
                                and boundary.intersects(bldg_rel)):
                            # Check it's not already covered by a hangar
                            # or terminal
                            is_dup = False
                            for existing in buildings:
                                try:
                                    existing_rel = shp_geom.Polygon([
                                        (c[0] - tile.lon, c[1] - tile.lat)
                                        for c in existing["footprint"]
                                    ])
                                    if (existing_rel.intersection(bldg_rel).area
                                            > 0.5 * bldg_rel.area):
                                        is_dup = True
                                        break
                                except Exception:
                                    pass
                            if not is_dup:
                                buildings.append({
                                    "footprint": coords,
                                    "source": "building",
                                })
                    except Exception:
                        pass
            except Exception:
                pass

        if buildings:
            result[airport] = buildings
            sources = {}
            for b in buildings:
                s = b.get("source", "unknown")
                sources[s] = sources.get(s, 0) + 1
            parts = ", ".join(
                "{} {}s".format(v, k) for k, v in sorted(sources.items())
            )
            UI.vprint(1, "   Auto-patch: {} buildings for {}: {}".format(
                len(buildings), airport, parts
            ))
    return result


# How far a road's elevation should influence neighbouring terrain.
ROAD_TERRAIN_INFLUENCE = 50.0  # meters

# Approximate meters per degree of latitude (good enough for the
# bounding-box buffer math here).
DEG_TO_M = 111120.0


def extract_road_info(dico_airports, tile, road_layer=None):
    """Extract road geometry near airports for terrain context.

    Roads (especially highways) near airports constrain terrain modeling:
    - A road along the airport boundary at a different elevation implies
      a retaining wall or grade transition.
    - Road centerlines provide known terrain elevations outside the airport.
    - Tunnels under the airport imply unrelated road elevation.

    Args:
        dico_airports: dict with processed airport data.
        tile: Tile object with .lat, .lon, .dem.
        road_layer: Optional OSM_layer with highway=* data.

    Returns:
        dict: {airport_key: [{'centerline': [(lon, lat), ...],
                              'highway_type': str, 'tunnel': bool,
                              'bridge': bool}, ...]}
    """
    result = {}
    if road_layer is None:
        return result

    for airport in dico_airports:
        apt = dico_airports[airport]
        boundary = apt.get("boundary")
        if boundary is None:
            continue

        roads = []
        try:
            from shapely import geometry as shp_geom
            # Buffer boundary to find nearby roads
            boundary_buf = boundary.buffer(
                ROAD_TERRAIN_INFLUENCE / DEG_TO_M)

            for wayid in road_layer.dicosmw:
                wtags = road_layer.dicosmtags.get("w", {}).get(wayid, {})
                hw_type = wtags.get("highway", "")
                if hw_type not in ("motorway", "trunk", "primary",
                                   "secondary", "tertiary",
                                   "motorway_link", "trunk_link",
                                   "primary_link"):
                    continue

                node_ids = road_layer.dicosmw[wayid]
                coords = []
                for nid in node_ids:
                    if nid in road_layer.dicosmn:
                        coords.append(tuple(road_layer.dicosmn[nid]))
                if len(coords) < 2:
                    continue

                # Check if road is near airport (tile-relative coords)
                road_rel = [(c[0] - tile.lon, c[1] - tile.lat)
                            for c in coords]
                try:
                    ls = shp_geom.LineString(road_rel)
                    if not boundary_buf.intersects(ls):
                        continue
                except Exception:
                    continue

                is_tunnel = wtags.get("tunnel") in ("yes", "building_passage")
                is_bridge = wtags.get("bridge") == "yes"
                roads.append({
                    "centerline": coords,
                    "highway_type": hw_type,
                    "tunnel": is_tunnel,
                    "bridge": is_bridge,
                })
        except Exception:
            pass

        if roads:
            result[airport] = roads
            UI.vprint(2, "   Auto-patch: {} roads near {}".format(
                len(roads), airport))
    return result
