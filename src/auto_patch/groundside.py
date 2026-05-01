"""Groundside (curbside / drop-off / parking) pavement emit.

Pulls roads tagged as airport-access or service highway out of the
OSM cache, lifts them off the DEM, and emits matching curbside
ribbon polygons.  Then prunes orphan junction polygons that are
fully contained inside (or touch only) the groundside ribbon —
those are road-island fragments that don't belong with airside
pavement.

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _emit_groundside_pavement_dem
    _drop_groundside_orphan_junctions
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Set, Tuple

from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon
from shapely.ops import linemerge, nearest_points, unary_union

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
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
    ROLE_RETAINING_WALL,
    SHARED_VERTEX_TOL_M,
)
from .pavement.vertices import _snap_polygon_vertices_to_rect_corners
from .elevation import _sample_dem, _resample_node_altitudes_nn


__all__ = [
    "_drop_groundside_orphan_junctions",
    "_emit_groundside_pavement_dem",
]


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
    from .pipeline import _load_osm_big_roads
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


