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

# Sloping taxi rect roles — per user 2026-05-19 these stay 4-corner
# through seam-crossings (split-at-seam, not insert-vertex-at-seam).
# Runways are intentionally NOT in this set: the runway seam pipeline
# converts to node_altitudes to preserve per-vertex precision at
# DEM-noisy tile boundaries (user 2026-05-13).
_TAXI_RECT_ROLES = {
    ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
    ROLE_STUB, ROLE_CROSS_CONNECTOR,
}

# Sub-meter tolerance for skipping insertions at existing vertices.
_EDGE_T_TOL = 1e-4


def _bucket_key(x: float, y: float) -> Tuple[int, int]:
    """Bucket key matching ``elevation._corner_elevation_bucket`` so
    Phase-2 can look up seam anchors directly against the solver's
    vertex graph."""
    s = 1.0 / SHARED_VERTEX_TOL_M  # 2.0
    return (int(round(x * s)), int(round(y * s)))


def _split_ring_at_seam(ring, seam_line):
    """Split a 4-corner ring at one seam line into 2 4-corner rings.

    Returns ``[ring_a, ring_b]`` when the seam cleanly intersects 2
    non-adjacent edges, ``None`` otherwise (no crossing, single
    crossing, or seam clips a corner — fall back to insert-in-place
    in those cases since a clean 4-corner split isn't available).
    """
    n = len(ring)
    if n != 4:
        return None
    intersections: List[Tuple[int, float, Tuple[float, float]]] = []
    for i in range(n):
        ax, ay = ring[i]
        bx, by = ring[(i + 1) % n]
        edge = LineString([(ax, ay), (bx, by)])
        try:
            inter = edge.intersection(seam_line)
        except _GEOM_EXC:
            continue
        if inter.is_empty or inter.geom_type != "Point":
            continue
        dx = bx - ax
        dy = by - ay
        L2 = dx * dx + dy * dy
        if L2 < 1e-9:
            continue
        t = ((inter.x - ax) * dx + (inter.y - ay) * dy) / L2
        if t <= _EDGE_T_TOL or t >= 1.0 - _EDGE_T_TOL:
            continue
        intersections.append((i, t, (inter.x, inter.y)))
    if len(intersections) != 2:
        return None
    intersections.sort(key=lambda r: r[0])
    (ea, _ta, pa), (eb, _tb, pb) = intersections
    # The seam must cross 2 NON-ADJACENT edges (opposite sides of
    # the rect) for a clean 2 × 4-corner split.  Adjacent-edge
    # crossings clip a corner off and produce a triangle + pentagon
    # — not the canonical 4-corner form.
    if (eb - ea) % n != 2:
        return None
    # Build sub-rect A: ring[0..ea], pa, pb, ring[eb+1..n-1]
    ring_a: List[Tuple[float, float]] = []
    for i in range(ea + 1):
        ring_a.append(ring[i])
    ring_a.append(pa)
    ring_a.append(pb)
    for i in range(eb + 1, n):
        ring_a.append(ring[i])
    # Build sub-rect B: pa, ring[ea+1..eb], pb
    ring_b: List[Tuple[float, float]] = [pa]
    for i in range(ea + 1, eb + 1):
        ring_b.append(ring[i])
    ring_b.append(pb)
    if len(ring_a) != 4 or len(ring_b) != 4:
        return None
    return [ring_a, ring_b]


def _split_taxi_rect_at_seams(
        shape: BuiltShape,
        cut_lines: List[LineString],
        anchor_keys: Set[Tuple[int, int]],
        layout: PavementLayout,
) -> Optional[List[BuiltShape]]:
    """Replace a 4-corner taxi rect with sub-rects produced by
    splitting at each seam crossing, preserving 4-corner geometry.

    Per user 2026-05-19: a taxi rect's slope rendering depends on
    the canonical 4-corner [HI-LEFT, LO-LEFT, LO-RIGHT, HI-RIGHT]
    ring convention (or equivalent CCW rotation), and EVERY
    downstream pass that operates on rects assumes 4 corners
    (absorption, junction-rule tests, ``_collect_junction_axes``,
    sloping-edge identification).  Inserting seam vertices into a
    sloping edge breaks that assumption: even a vertex collinear
    with its neighbours produces 5- or 6-corner rings that the
    rest of the pipeline rejects.

    Splitting the rect at each seam produces N + 1 sub-rects, each
    still 4-corner.  Sub-rect altitudes are intentionally left
    unset — the elevation solver's first pass fills them after
    HARD-anchoring the seam corners (recorded in ``anchor_keys``)
    to ``dem.alt_strict``.

    Returns ``None`` if any seam crossing doesn't admit a clean
    2 × 4-corner split (e.g. a seam clips a single corner); the
    caller then falls back to the legacy insert-vertices path so
    no shape is dropped.
    """
    if shape.polygon is None or shape.polygon.is_empty:
        return None
    ring = list(shape.polygon.exterior.coords)
    if ring and ring[0] == ring[-1]:
        ring = ring[:-1]
    if len(ring) != 4:
        return None

    rings: List[List[Tuple[float, float]]] = [
        [(float(x), float(y)) for x, y in ring]]
    new_corner_pts: List[Tuple[float, float]] = []
    original_pt_keys = {_bucket_key(x, y) for x, y in ring}
    for seam_line in cut_lines:
        next_rings: List[List[Tuple[float, float]]] = []
        for r in rings:
            poly_r = Polygon(r)
            try:
                if not poly_r.boundary.intersects(seam_line):
                    next_rings.append(r)
                    continue
            except _GEOM_EXC:
                next_rings.append(r)
                continue
            split = _split_ring_at_seam(r, seam_line)
            if split is None:
                # Couldn't split cleanly — abort and let caller fall
                # back to insert-in-place.
                return None
            next_rings.extend(split)
            for sub in split:
                for x, y in sub:
                    if _bucket_key(x, y) in original_pt_keys:
                        continue
                    new_corner_pts.append((x, y))
        rings = next_rings

    if len(rings) < 2:
        return None

    # Route every new corner through the canonical-point registry
    # (so adjacent shapes see the same coords) and record bucket
    # keys for the solver's HARD-anchor pass.
    registry = getattr(layout, "canonical_points", None)

    def _canon(x: float, y: float) -> Tuple[float, float]:
        if registry is None:
            return (x, y)
        return registry.get_or_add(float(x), float(y))

    for x, y in new_corner_pts:
        cx, cy = _canon(x, y)
        anchor_keys.add(_bucket_key(cx, cy))

    out: List[BuiltShape] = []
    import copy
    for r in rings:
        canon_ring = [_canon(x, y) for x, y in r]
        # Drop degenerate rings (sub-tol vertices).
        canon_ring_dedup: List[Tuple[float, float]] = []
        for pt in canon_ring:
            if not canon_ring_dedup or canon_ring_dedup[-1] != pt:
                canon_ring_dedup.append(pt)
        if (len(canon_ring_dedup) >= 2
                and canon_ring_dedup[0] == canon_ring_dedup[-1]):
            canon_ring_dedup = canon_ring_dedup[:-1]
        if len(canon_ring_dedup) != 4:
            return None
        new_s = copy.copy(shape)
        new_s.polygon = Polygon(canon_ring_dedup + [canon_ring_dedup[0]])
        new_s.altitude = None
        new_s.altitude_high = None
        new_s.altitude_low = None
        new_s.node_altitudes = None
        # source_axis inherited via copy.copy — both sub-rects share
        # the parent's axis direction, which is what
        # ``_canonicalise_rect`` needs at writeback time.
        out.append(new_s)
    return out


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
    # Per user 2026-05-19: don't add vertices to a sloping taxi rect
    # — every downstream pass (absorption, junction-rule tests,
    # _collect_junction_axes) assumes a canonical 4-corner ring and
    # breaks when extra vertices appear on a sloping edge.  Instead,
    # for taxi rect roles, SPLIT the rect at each seam into 2
    # 4-corner sub-rects so each sub-rect remains 4-corner.  Altitude
    # handling: leave altitudes unset on the sub-rects; the elevation
    # solver fills them after seeding seam corners from DEM (the
    # seam keys recorded here drive ``_seed_elevations``' HARD-anchor
    # pass).
    new_shapes_extra: List[BuiltShape] = []
    indices_to_drop: List[int] = []
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
        if shape.role in _TAXI_RECT_ROLES:
            sub_rects = _split_taxi_rect_at_seams(
                shape, cut_lines, anchor_keys, layout)
            if sub_rects is not None and len(sub_rects) >= 2:
                indices_to_drop.append(i)
                new_shapes_extra.extend(sub_rects)
                continue
            # Fall through to insert-in-place if split didn't apply.
        new_shape = _insert_seam_vertices(shape, cut_lines, anchor_keys)
        if new_shape is not None:
            layout.shapes[i] = new_shape
    if indices_to_drop:
        keep_set = set(range(len(layout.shapes))) - set(indices_to_drop)
        layout.shapes = [layout.shapes[i] for i in sorted(keep_set)]
        layout.shapes.extend(new_shapes_extra)

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
