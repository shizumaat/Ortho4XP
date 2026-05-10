"""Airport boundary ribbon + boundary→DEM bridge polygons.

Two emitters:

* ``_emit_airport_boundary_shape`` — closed-ring airport boundary
  shape that wraps the union of all emitted airside pavement,
  used by downstream consumers to clip terrain underlay.
* ``_emit_boundary_dem_bridge`` — wedge polygons that bridge the
  airport boundary to the surrounding DEM where the airport sits
  noticeably above or below the natural terrain (CYXY's plateau,
  HECA's berm), preventing visual cliffs.

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _emit_airport_boundary_shape
    _emit_boundary_dem_bridge
"""
from __future__ import annotations

import math
import os
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
from .pavement.junctions import _decompose_polygon_with_holes
from .elevation import _resample_node_altitudes_nn, _sample_dem


__all__ = [
    "_emit_airport_boundary_shape",
    "_emit_boundary_dem_bridge",
]


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
    from .pipeline import _load_osm_airports, _load_osm_big_roads
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
    from .pipeline import _load_osm_airports, _load_osm_big_roads
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
            if (s.altitude_high is not None
                    and s.altitude_low is not None):
                # Sloped 4-corner rect — strict convention.
                if len(coords) != 4:
                    continue
                per = [s.altitude_high, s.altitude_low,
                       s.altitude_low, s.altitude_high]
            elif s.altitude is not None:
                # Flat shape: any number of corners (multi-node flat
                # runway shapes from the segmenter).
                if len(coords) < 4:
                    continue
                per = [float(s.altitude)] * len(coords)
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


