"""Insert seam vertices at integer lat/lon tile-boundary lines.

For each pavement shape whose polygon's exterior boundary crosses an
integer latitude or longitude line within the airport footprint, this
module:

  1. Inserts new ring vertices at the crossing points (deterministic
     from polygon geometry alone, independent of which tile is being
     built — both tile builds produce the same vertices).
  2. Converts sloped 4-corner rects to ``node_altitudes`` representation
     so each vertex (including seam crossings) carries its own altitude.
  3. Records the seam-vertex bucket keys in
     ``layout._seam_anchor_keys`` for the Phase-2 elevation solver to
     HARD-anchor against ``dem.alt_strict``.

Cross-tile parity: each tile build runs over the same pavement
geometry, finds the same cut lines, and inserts vertices at the same
(x, y) positions.  Phase-2 then samples the same SRTM pixel at each
seam vertex (SRTM .hgt overlap row), so all tile builds compute the
same altitude.  ``cut_layout_at_tile_boundaries`` then keeps only the
shape pieces falling in the current tile.
"""
from __future__ import annotations

import math
from typing import List, Optional, Set, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, Polygon
from shapely.ops import unary_union

from .layout import (
    BuiltShape, PavementLayout, R_EARTH, SHARED_VERTEX_TOL_M,
    ROLE_APRON, ROLE_BOUNDARY, ROLE_CROSS_CONNECTOR, ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL, ROLE_RUNWAY, ROLE_SECONDARY_PARALLEL,
    ROLE_STUB, ROLE_TERMINAL, ROLE_TUNNEL_RAMP, ROLE_RETAINING_WALL,
    ROLE_GROUNDSIDE_PAVEMENT,
)

__all__ = ["split_pavement_at_seams", "apply_seam_dem_anchors"]

_GEOM_EXC = (ValueError, TypeError, GEOSException,
             TopologicalError, IndexError)

# Shape roles whose polygons participate in seam-splitting.
# Tile-cut bridges are intentionally excluded — they're emitted later
# (in tile_cut.py) for backwards-compat with bridge-based seam pinning;
# once the new seam-anchor pass proves out we can drop bridges.
_SEAM_SPLIT_ROLES = {
    ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
    ROLE_STUB, ROLE_CROSS_CONNECTOR, ROLE_APRON, ROLE_TERMINAL,
    ROLE_JUNCTION, ROLE_BOUNDARY, ROLE_TUNNEL_RAMP, ROLE_RETAINING_WALL,
    ROLE_GROUNDSIDE_PAVEMENT,
}

# Sub-meter tolerance for skipping insertions at existing vertices.
_EDGE_T_TOL = 1e-4


def _bucket_key(x: float, y: float) -> Tuple[int, int]:
    """Bucket key matching ``elevation._corner_elevation_bucket`` so
    Phase-2 can look up seam anchors directly against the solver's
    vertex graph."""
    s = 1.0 / SHARED_VERTEX_TOL_M  # 2.0
    return (int(round(x * s)), int(round(y * s)))


def split_pavement_at_seams(layout: PavementLayout) -> int:
    """Insert seam vertices and convert sloped rects to ``node_altitudes``.

    Records seam-vertex bucket keys on ``layout._seam_anchor_keys``.

    Returns the net change in shape count (always 0 today — this pass
    modifies existing shapes in place rather than splitting them).
    """
    layout._seam_anchor_keys = set()  # type: ignore[attr-defined]

    if not layout.shapes or layout.anchor is None:
        return 0
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))

    # Build footprint of all pavement shapes to identify which
    # integer lines pass through the airport.
    pav_polys = [s.polygon for s in layout.shapes
                 if s.polygon is not None and not s.polygon.is_empty]
    if not pav_polys:
        return 0
    try:
        pav_union = unary_union(pav_polys)
    except _GEOM_EXC:
        return 0
    minx, miny, maxx, maxy = pav_union.bounds
    min_lat = lat0 + math.degrees(miny / R_EARTH)
    max_lat = lat0 + math.degrees(maxy / R_EARTH)
    min_lon = lon0 + math.degrees(minx / (R_EARTH * cos0))
    max_lon = lon0 + math.degrees(maxx / (R_EARTH * cos0))

    cut_lines: List[LineString] = []
    for lat_int in range(int(math.ceil(min_lat)),
                          int(math.floor(max_lat)) + 1):
        if min_lat < lat_int < max_lat:
            y_int = math.radians(lat_int - lat0) * R_EARTH
            cut_lines.append(LineString([
                (minx - 100.0, y_int), (maxx + 100.0, y_int)]))
    for lon_int in range(int(math.ceil(min_lon)),
                          int(math.floor(max_lon)) + 1):
        if min_lon < lon_int < max_lon:
            x_int = math.radians(lon_int - lon0) * R_EARTH * cos0
            cut_lines.append(LineString([
                (x_int, miny - 100.0), (x_int, maxy + 100.0)]))
    if not cut_lines:
        return 0

    anchor_keys: Set[Tuple[int, int]] = set()
    for i, shape in enumerate(layout.shapes):
        if shape.role not in _SEAM_SPLIT_ROLES:
            continue
        if shape.polygon is None or shape.polygon.is_empty:
            continue
        # Quick reject: skip if no boundary crossing.
        try:
            if not shape.polygon.boundary.intersects(
                    unary_union(cut_lines)):
                continue
        except _GEOM_EXC:
            continue
        new_shape = _insert_seam_vertices(shape, cut_lines, anchor_keys)
        if new_shape is not None:
            layout.shapes[i] = new_shape

    # Per user 2026-05-13: when ANY sub-rect of a runway has been
    # seam-converted to node_altitudes, ALL sub-rects of that same
    # runway need node_altitudes too — otherwise altitude_high/low's
    # planar-surface assumption forces averaging across adjacent
    # sub-rects' shared corners that no longer agree (the
    # seam-crossing sub-rect has its H corner pinned to DEM while
    # the neighbour sub-rect has its L corner at CIFP).  Convert the
    # entire runway chain to node_altitudes so each corner carries
    # its own altitude through the solver.
    seam_runway_refs: Set[str] = set()
    for shape in layout.shapes:
        if (shape.role == ROLE_RUNWAY
                and shape.node_altitudes
                and shape.ref):
            seam_runway_refs.add(shape.ref)
    if seam_runway_refs:
        for shape in layout.shapes:
            if shape.role != ROLE_RUNWAY:
                continue
            if shape.ref not in seam_runway_refs:
                continue
            if shape.node_altitudes:
                continue
            if shape.polygon is None or shape.polygon.is_empty:
                continue
            ring = list(shape.polygon.exterior.coords)
            if ring and ring[0] == ring[-1]:
                ring = ring[:-1]
            if (len(ring) == 4
                    and shape.altitude_high is not None
                    and shape.altitude_low is not None):
                alts = [
                    float(shape.altitude_high),
                    float(shape.altitude_low),
                    float(shape.altitude_low),
                    float(shape.altitude_high),
                ]
                shape.node_altitudes = alts + [alts[0]]
                shape.altitude_high = None
                shape.altitude_low = None
            elif shape.altitude is not None:
                alts = [float(shape.altitude)] * len(ring)
                shape.node_altitudes = alts + [alts[0]]
                shape.altitude = None

    layout._seam_anchor_keys = anchor_keys  # type: ignore[attr-defined]
    return 0


def _insert_seam_vertices(
        shape: BuiltShape,
        cut_lines: List[LineString],
        anchor_keys: Set[Tuple[int, int]]) -> Optional[BuiltShape]:
    """Insert intersection points of cut_lines with the shape's
    exterior ring, return a new BuiltShape with ``node_altitudes`` set.

    The new vertices' altitudes are placeholders (0.0); Phase-2 sets
    them to ``dem.alt_strict``.  Existing-vertex altitudes are
    preserved (sloped rect → unpacked via [H, L, L, H] convention;
    flat → broadcast; per-vertex → carried through).
    """
    poly = shape.polygon
    ring = list(poly.exterior.coords)
    if ring and ring[0] == ring[-1]:
        ring = ring[:-1]
    n_orig = len(ring)
    if n_orig < 3:
        return None

    # Determine the original per-vertex altitudes.
    if shape.node_altitudes:
        old_alts = list(shape.node_altitudes[:n_orig])
        if len(old_alts) < n_orig:
            old_alts += [old_alts[-1]] * (n_orig - len(old_alts))
    elif (shape.altitude_high is not None
            and shape.altitude_low is not None
            and n_orig == 4):
        old_alts = [
            float(shape.altitude_high),
            float(shape.altitude_low),
            float(shape.altitude_low),
            float(shape.altitude_high),
        ]
    elif shape.altitude is not None:
        old_alts = [float(shape.altitude)] * n_orig
    else:
        # No elevation data yet — Phase 2 will fill all altitudes.
        old_alts = [0.0] * n_orig

    # Walk each edge, find intersections with each cut line, insert in
    # parametric order.  Track which inserted vertices are seam-anchored
    # AND which existing vertices sit on a seam.
    new_ring: List[Tuple[float, float]] = []
    new_alts: List[float] = []
    inserted_idxs: List[int] = []
    existing_on_seam: List[int] = []  # indices in new_ring of original
                                       # ring vertices that lie on a seam

    for i in range(n_orig):
        p1 = ring[i]
        p2 = ring[(i + 1) % n_orig]
        new_ring.append(p1)
        new_alts.append(old_alts[i])
        edge = LineString([p1, p2])
        edge_len = edge.length
        if edge_len < 1e-6:
            continue
        ips: List[Tuple[float, Tuple[float, float]]] = []
        for cl in cut_lines:
            try:
                inter = edge.intersection(cl)
            except _GEOM_EXC:
                continue
            if inter.is_empty:
                continue
            if inter.geom_type == "Point":
                pt = (inter.x, inter.y)
            else:
                continue
            dx = pt[0] - p1[0]
            dy = pt[1] - p1[1]
            t = math.hypot(dx, dy) / edge_len
            if t < _EDGE_T_TOL:
                # Cut line passes through p1 (current vertex).
                anchor_keys.add(_bucket_key(p1[0], p1[1]))
                # Record it for the shape-level conversion check.
                existing_on_seam.append(len(new_ring) - 1)
                continue
            if t > 1.0 - _EDGE_T_TOL:
                # Cut line passes through p2 (handled next iteration).
                anchor_keys.add(_bucket_key(p2[0], p2[1]))
                continue
            ips.append((t, pt))
        ips.sort(key=lambda x: x[0])
        a1 = old_alts[i]
        a2 = old_alts[(i + 1) % n_orig]
        for t, pt in ips:
            interp_alt = a1 + t * (a2 - a1)
            inserted_idxs.append(len(new_ring))
            new_ring.append(pt)
            new_alts.append(interp_alt)
            anchor_keys.add(_bucket_key(pt[0], pt[1]))

    if not inserted_idxs and not existing_on_seam:
        return None

    # Build the new BuiltShape.  Switch to node_altitudes representation
    # since the original [H, L, L, H] or flat scheme no longer applies
    # cleanly with N+ vertices.
    try:
        new_poly = Polygon(new_ring)
        if not new_poly.is_valid:
            new_poly = new_poly.buffer(0)
            if (new_poly.geom_type != "Polygon"
                    or new_poly.is_empty):
                return None
    except _GEOM_EXC:
        return None
    # node_altitudes carries the CLOSING repeat per layout convention.
    closed_alts = new_alts + [new_alts[0]]
    # If the source was a sloped 4-corner rect without an explicit
    # source_axis (typical of runway shapes built from CIFP), derive
    # one from the H→L pair so downstream Stage A regrade can project
    # vertices onto the runway centerline.  For non-rect shapes the
    # original source_axis (if any) is carried through.
    derived_axis = shape.source_axis
    if (derived_axis is None
            and shape.altitude_high is not None
            and shape.altitude_low is not None
            and n_orig == 4):
        c0, c1, c2, c3 = ring[0], ring[1], ring[2], ring[3]
        h_mid = (0.5 * (c0[0] + c3[0]), 0.5 * (c0[1] + c3[1]))
        l_mid = (0.5 * (c1[0] + c2[0]), 0.5 * (c1[1] + c2[1]))
        if (h_mid[0] - l_mid[0]) ** 2 + (h_mid[1] - l_mid[1]) ** 2 > 1e-6:
            derived_axis = LineString([h_mid, l_mid])

    new_shape = BuiltShape(
        polygon=new_poly,
        role=shape.role,
        ref=shape.ref,
        source_axis=derived_axis,
        altitude=None,
        altitude_high=None,
        altitude_low=None,
        node_altitudes=closed_alts,
        is_bridge=shape.is_bridge,
    )
    return new_shape


def apply_seam_dem_anchors(
    layout: PavementLayout,
    dem,
    tile_lat: int,
    tile_lon: int,
) -> int:
    """Sample DEM at every seam vertex and overwrite the placeholder
    interp altitude in ``node_altitudes`` with the DEM value.

    Uses ``dem.alt_strict`` so the sample is the raw HGT pixel
    (preserve_boundary keeps this at the .hgt overlap value, identical
    in both tiles' DEMs).  Both tiles' builds sample the same lat/lon
    point and get the same value.

    Must be called AFTER ``split_pavement_at_seams`` (which populates
    ``layout._seam_anchor_keys``) and BEFORE the elevation solver.

    Returns the number of vertices updated.
    """
    anchor_keys = getattr(layout, "_seam_anchor_keys", None)
    if not anchor_keys:
        return 0
    if dem is None:
        return 0
    nodata = getattr(dem, "nodata", -32768)
    n_updated = 0
    for shape in layout.shapes:
        if not shape.node_altitudes:
            continue
        if shape.polygon is None or shape.polygon.is_empty:
            continue
        ring = list(shape.polygon.exterior.coords)
        if ring and ring[0] == ring[-1]:
            ring = ring[:-1]
        alts = list(shape.node_altitudes[:len(ring)])
        changed = False
        for i, (x, y) in enumerate(ring):
            if _bucket_key(x, y) not in anchor_keys:
                continue
            lat, lon = layout.m_to_ll(x, y)
            try:
                v = float(dem.alt_strict(
                    (lon - tile_lon, lat - tile_lat)))
            except _GEOM_EXC:
                continue
            if v != v or v == nodata:  # NaN or no-data
                continue
            alts[i] = round(v, 1)
            changed = True
            n_updated += 1
        if changed:
            shape.node_altitudes = alts + [alts[0]]
    return n_updated
