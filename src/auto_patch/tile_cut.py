"""Cut shapes along integer lat/lon tile boundaries.

X-Plane / Ortho4XP renders each 1° × 1° lat-lon tile as a separate
DSF file.  A single shape that spans two tiles is awkwardly bisected
at the seam, producing visual artifacts.  Per user 2026-05-10: for
each integer lat or lon line passing through the airport's pavement
footprint, build a buffered line (``half_width_m`` each side, so a
10 m strip by default) and subtract it from every shape.  Shapes
that split into multiple pieces are replaced with separate
``BuiltShape`` entries; sloped 4-corner rects convert to per-vertex
``node_altitudes`` (cut pieces are non-rectangular and the legacy
``[H, L, L, H]`` 4-corner convention no longer applies).

The mechanism mirrors how the pavement builder clips runway shapes
out of ``pav_union`` — just at a different geometric target.  Ortho4XP
and X-Plane stitch the resulting tile seams together at render time.

Public API:
    cut_layout_at_tile_boundaries
"""
from __future__ import annotations

import copy
import math
from typing import Callable, List, Optional

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from .layout import BuiltShape, PavementLayout, R_EARTH


# Same narrow exception set used in ``boundary.py`` — covers real
# shapely degeneracy without masking programming errors.
_GEOM_EXC = (ValueError, TypeError, GEOSException,
             TopologicalError, IndexError)


__all__ = ["cut_layout_at_tile_boundaries"]


def cut_layout_at_tile_boundaries(
        layout: PavementLayout,
        half_width_m: float = 5.0,
        min_piece_area_m2: float = 1.0,
        current_tile_lat: Optional[int] = None,
        current_tile_lon: Optional[int] = None) -> int:
    """Cut every shape crossing an integer lat or lon tile boundary,
    leaving a ``2 * half_width_m`` wide gap (default 10 m).

    Then DROP every shape piece whose representative point falls
    outside the current tile.  The neighbour-tile auto_patch run
    will generate the patch covering its portion.

    ``current_tile_lat`` / ``current_tile_lon`` identify the tile
    being processed by the driver — these can differ from the
    airport's anchor tile when a cross-tile airport is being
    processed during a NEIGHBOUR-tile build (e.g. Ortho4XP
    generates tile -13/-78, which includes SPLP because the
    airport extends into it, but SPLP's anchor is in -13/-77).
    When None, fall back to ``floor(layout.anchor)`` (the
    airport-anchor tile) — backward-compatible default for tests
    and direct ``build_airport_pavement`` calls.

    Mutates ``layout.shapes`` in place.  Returns the net change in
    shape count (positive when shapes split, negative when slivers
    fall below ``min_piece_area_m2`` and get dropped).
    """
    if not layout.shapes or layout.anchor is None:
        return 0
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))

    polys = [s.polygon for s in layout.shapes
             if s.polygon is not None and not s.polygon.is_empty]
    if not polys:
        return 0
    try:
        union = unary_union(polys)
    except _GEOM_EXC:
        return 0
    minx, miny, maxx, maxy = union.bounds

    # Footprint bounds in lat/lon.
    min_lat = lat0 + math.degrees(miny / R_EARTH)
    max_lat = lat0 + math.degrees(maxy / R_EARTH)
    min_lon = lon0 + math.degrees(minx / (R_EARTH * cos0))
    max_lon = lon0 + math.degrees(maxx / (R_EARTH * cos0))

    # Integer lat lines strictly inside the airport's lat range.
    cut_lines: List[LineString] = []
    for lat_int in range(
            int(math.ceil(min_lat)), int(math.floor(max_lat)) + 1):
        if min_lat < lat_int < max_lat:
            y_int = math.radians(lat_int - lat0) * R_EARTH
            cut_lines.append(LineString([
                (minx - 100.0, y_int), (maxx + 100.0, y_int)]))
    # Integer lon lines strictly inside the airport's lon range.
    for lon_int in range(
            int(math.ceil(min_lon)), int(math.floor(max_lon)) + 1):
        if min_lon < lon_int < max_lon:
            x_int = math.radians(lon_int - lon0) * R_EARTH * cos0
            cut_lines.append(LineString([
                (x_int, miny - 100.0), (x_int, maxy + 100.0)]))

    if not cut_lines:
        return 0

    try:
        cut_polys = [line.buffer(half_width_m, cap_style=2)
                     for line in cut_lines]
        cut_union = unary_union(cut_polys)
    except _GEOM_EXC:
        return 0

    # Cut the airport-boundary polygon itself.  Downstream
    # consumers (Ortho4XP encoder / boundary ribbon emit / DEM
    # bridge) see the cut MultiPolygon and naturally produce
    # per-tile geometry instead of spanning-tile artifacts.
    if (layout.airport_boundary is not None
            and not layout.airport_boundary.is_empty):
        try:
            ab_cut = layout.airport_boundary.difference(cut_union)
        except _GEOM_EXC:
            ab_cut = None
        if ab_cut is not None and not ab_cut.is_empty:
            if ab_cut.geom_type in ("Polygon", "MultiPolygon"):
                layout.airport_boundary = ab_cut

    # The CURRENT tile (the one this auto_patch run is generating)
    # is the airport-anchor tile.  Per user 2026-05-12: after the
    # cut, drop any shape (or shape piece) that's not inside the
    # current tile — when the neighbour tile is processed in its own
    # auto_patch run, IT generates the patch covering its portion of
    # the airport.  Without this drop, ``_runway_clamped_alt`` etc.
    # would try to sample the neighbour-tile DEM (which isn't
    # loaded) and substitute 0 m, producing altitude-0 boundary
    # rects in X-Plane.
    cur_tile_lat = (current_tile_lat if current_tile_lat is not None
                    else int(math.floor(lat0)))
    cur_tile_lon = (current_tile_lon if current_tile_lon is not None
                    else int(math.floor(lon0)))

    def _in_current_tile(poly: Polygon) -> bool:
        try:
            c = poly.representative_point()
        except _GEOM_EXC:
            try:
                c = poly.centroid
            except _GEOM_EXC:
                return True  # fail open
        lat = lat0 + math.degrees(c.y / R_EARTH)
        lon = lon0 + math.degrees(c.x / (R_EARTH * cos0))
        return (cur_tile_lat <= lat < cur_tile_lat + 1
                and cur_tile_lon <= lon < cur_tile_lon + 1)

    # Also clip ``layout.airport_boundary`` to the current tile.
    if (layout.airport_boundary is not None
            and not layout.airport_boundary.is_empty):
        ab = layout.airport_boundary
        if ab.geom_type == "MultiPolygon":
            kept = [g for g in ab.geoms
                    if g.geom_type == "Polygon" and not g.is_empty
                    and _in_current_tile(g)]
            if not kept:
                layout.airport_boundary = None
            elif len(kept) == 1:
                layout.airport_boundary = kept[0]
            else:
                from shapely.geometry import MultiPolygon
                layout.airport_boundary = MultiPolygon(kept)
        elif ab.geom_type == "Polygon":
            if not _in_current_tile(ab):
                layout.airport_boundary = None

    n_before = len(layout.shapes)
    new_shapes: List[BuiltShape] = []
    for s in layout.shapes:
        if s.polygon is None or s.polygon.is_empty:
            new_shapes.append(s)
            continue
        try:
            if not s.polygon.intersects(cut_union):
                # No cut — keep iff the shape is in the current tile.
                if _in_current_tile(s.polygon):
                    new_shapes.append(s)
                continue
            diff = s.polygon.difference(cut_union)
        except _GEOM_EXC:
            new_shapes.append(s)
            continue
        if diff.is_empty:
            # Source polygon entirely inside the cut buffer.  No
            # pavement pieces but a bridge will be emitted below.
            pieces: List[Polygon] = []
        elif diff.geom_type == "Polygon":
            pieces = [diff]
        elif diff.geom_type == "MultiPolygon":
            pieces = [g for g in diff.geoms
                      if g.geom_type == "Polygon" and not g.is_empty]
        else:
            # Unexpected result (e.g. non-empty GeometryCollection);
            # keep the original shape and skip the cut + bridge.
            if _in_current_tile(s.polygon):
                new_shapes.append(s)
            continue
        pieces = [p for p in pieces if p.area >= min_piece_area_m2]
        pieces = [p for p in pieces if _in_current_tile(p)]

        slope_sampler = _make_slope_sampler(s)
        for piece in pieces:
            new_s = _build_piece_shape(s, piece, slope_sampler)
            if new_s is not None:
                new_shapes.append(new_s)
    layout.shapes = new_shapes
    return len(layout.shapes) - n_before


def _make_slope_sampler(
        s: BuiltShape) -> Optional[Callable[[float, float], float]]:
    """Build a closure that samples a sloped 4-corner rect's
    elevation at any (x, y) by projecting onto the high-mid → low-mid
    axis.  Returns None when ``s`` isn't a 4-corner sloped rect.
    """
    if s.altitude_high is None or s.altitude_low is None:
        return None
    if s.polygon is None or s.polygon.is_empty:
        return None
    try:
        coords = list(s.polygon.exterior.coords)
    except _GEOM_EXC:
        return None
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return None
    high_mid_x = 0.5 * (coords[0][0] + coords[3][0])
    high_mid_y = 0.5 * (coords[0][1] + coords[3][1])
    low_mid_x = 0.5 * (coords[1][0] + coords[2][0])
    low_mid_y = 0.5 * (coords[1][1] + coords[2][1])
    ax = low_mid_x - high_mid_x
    ay = low_mid_y - high_mid_y
    L2 = ax * ax + ay * ay
    H = float(s.altitude_high)
    L = float(s.altitude_low)
    if L2 < 1e-6:
        avg = 0.5 * (H + L)
        return lambda x, y: avg

    def sample(x: float, y: float) -> float:
        t = ((x - high_mid_x) * ax + (y - high_mid_y) * ay) / L2
        if t < 0.0:
            t = 0.0
        elif t > 1.0:
            t = 1.0
        return H + t * (L - H)
    return sample


def _build_piece_shape(
        orig: BuiltShape,
        piece: Polygon,
        slope_sampler: Optional[Callable[[float, float], float]],
) -> Optional[BuiltShape]:
    """Construct a BuiltShape for one cut piece, copying tags from
    ``orig`` and resampling altitudes for the new polygon vertices.

    * Flat shape (``altitude`` set, ``altitude_high`` None):
      keeps ``altitude`` unchanged; ``node_altitudes`` cleared.
    * Sloped 4-corner rect: convert to ``node_altitudes`` by
      projecting each new vertex onto the original H→L axis.  The
      polygon's vertex count typically differs from 4 post-cut so
      the legacy [H, L, L, H] convention no longer applies.
    * Per-vertex ``node_altitudes``: resample via nearest-neighbour
      against the original ring.
    """
    new_s = copy.copy(orig)
    new_s.polygon = piece

    # Flat with single altitude — corner count irrelevant.
    if orig.altitude is not None and orig.altitude_high is None:
        new_s.altitude = orig.altitude
        new_s.altitude_high = None
        new_s.altitude_low = None
        new_s.node_altitudes = None
        return new_s

    # Sloped 4-corner rect → per-vertex node_altitudes.
    if slope_sampler is not None:
        try:
            coords = list(piece.exterior.coords)
        except _GEOM_EXC:
            return None
        alts = [round(float(slope_sampler(x, y)), 1)
                for (x, y) in coords]
        new_s.node_altitudes = alts
        new_s.altitude = None
        new_s.altitude_high = None
        new_s.altitude_low = None
        return new_s

    # Per-vertex node_altitudes → resample via NN.
    if orig.node_altitudes:
        try:
            old_coords = list(orig.polygon.exterior.coords)
        except _GEOM_EXC:
            return None
        if old_coords and old_coords[0] == old_coords[-1]:
            old_open = old_coords[:-1]
        else:
            old_open = old_coords
        from .elevation import _resample_node_altitudes_nn
        new_alts = _resample_node_altitudes_nn(
            piece, old_open, orig.node_altitudes)
        if new_alts is not None:
            new_s.node_altitudes = new_alts
        return new_s

    # No elevation data — keep as-is with the new polygon.
    return new_s
