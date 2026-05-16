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

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon
from shapely.ops import linemerge, nearest_points, unary_union

# Narrow exception tuple for shapely / geometry ops that signal
# degenerate input rather than a programming error.  Replaces the
# blanket ``except Exception`` blocks that previously silently
# swallowed ``NameError`` from a missing import (user 2026-05-10
# — boundary runway-elevation clamp had been broken since the
# slice-5 refactor because the import dance masked a NameError).
#
# Programming errors (``NameError``, ``ImportError``,
# ``AttributeError`` from typos / ``None``-leaks) intentionally
# propagate so they surface immediately during testing rather than
# being silently masked at runtime.  Real shapely degeneracy
# surfaces as ``GEOSException`` / ``TopologicalError`` /
# ``ValueError``; out-of-bounds DEM indexing surfaces as
# ``IndexError``.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)


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
from .pavement.runways import _sample_runway_segment_elev
from .elevation import _resample_node_altitudes_nn, _sample_dem


__all__ = [
    "_emit_airport_boundary_shape",
    "_emit_boundary_dem_bridge",
    "_clip_boundary_bridges_against_pavement",
]


def _clip_boundary_bridges_against_pavement(
        layout: "PavementLayout",
        min_area_m2: float = 25.0) -> int:
    """Post-process: re-subtract pavement (junction/terminal/rect/
    runway) from every ``boundary_dem_bridge`` shape.

    ``_emit_boundary_dem_bridge`` subtracts the junctions/terminals
    present AT EMIT TIME, but downstream passes (per_surface_solve,
    subdivide_violating_junctions, stitch_pavement_polygons,
    _split_sloped_rects_at_violations) reshape pavement polygons —
    a junction may merge with a neighbour, a subdivide may grow
    a junction across the bridge boundary, etc.  Any such growth
    creates a stale overlap because the bridge was clipped against
    the bridge's emit-time snapshot.

    This pass runs LAST, just before tile_cut, against the final
    pavement geometry.  Per user 2026-05-13 (CYXY way -10483
    overlap report): zero tolerance for bridge↔pavement overlap.

    Returns the number of bridge shapes modified (clipped or dropped).
    """
    bridges = [s for s in layout.shapes
               if s.role == ROLE_BOUNDARY
               and s.ref == "boundary_dem_bridge"
               and s.polygon is not None
               and not s.polygon.is_empty]
    if not bridges:
        return 0

    # Roles that bridges must NOT overlap.  We exclude other
    # boundary shapes (ribbon + DEM bridges) because they share
    # vertices by design at the airport perimeter.
    NON_BRIDGE_PAVEMENT = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
        ROLE_STUB, ROLE_CROSS_CONNECTOR,
        ROLE_JUNCTION, ROLE_TERMINAL, ROLE_APRON,
    }
    obstacles = [s for s in layout.shapes
                 if s.role in NON_BRIDGE_PAVEMENT
                 and s.polygon is not None
                 and not s.polygon.is_empty]
    if not obstacles:
        return 0

    n_modified = 0
    new_shapes: List[BuiltShape] = []
    for s in layout.shapes:
        if (s.role != ROLE_BOUNDARY
                or s.ref != "boundary_dem_bridge"
                or s.polygon is None
                or s.polygon.is_empty):
            new_shapes.append(s)
            continue
        bridge_poly = s.polygon
        old_alts = s.node_altitudes
        old_open = list(bridge_poly.exterior.coords)
        if old_open and old_open[0] == old_open[-1]:
            old_open = old_open[:-1]
        modified = False
        for obs in obstacles:
            try:
                if not bridge_poly.intersects(obs.polygon):
                    continue
                inter_area = bridge_poly.intersection(obs.polygon).area
                if inter_area <= 0.0:
                    continue
                bridge_poly = bridge_poly.difference(obs.polygon)
                modified = True
            except _GEOM_EXC:
                continue
            if bridge_poly.is_empty:
                break
            if bridge_poly.geom_type not in (
                    "Polygon", "MultiPolygon",
                    "GeometryCollection"):
                bridge_poly = None
                break
        if bridge_poly is None or bridge_poly.is_empty:
            n_modified += 1
            continue
        if not modified:
            new_shapes.append(s)
            continue
        # Extract Polygon members.  difference() can yield Polygon,
        # MultiPolygon, or GeometryCollection (when subtracted
        # boundaries touch at points/edges).
        if bridge_poly.geom_type == "Polygon":
            pieces = [bridge_poly]
        elif bridge_poly.geom_type == "MultiPolygon":
            pieces = list(bridge_poly.geoms)
        elif bridge_poly.geom_type == "GeometryCollection":
            pieces = [g for g in bridge_poly.geoms
                      if g.geom_type == "Polygon"]
        else:
            pieces = []
        pieces = [p for p in pieces
                  if p.is_valid and not p.is_empty
                  and p.area >= min_area_m2]
        if not pieces:
            n_modified += 1
            continue
        # Keep the largest piece (consistent with emit-time logic).
        pieces.sort(key=lambda g: -g.area)
        keep = pieces[0]
        new_s = BuiltShape(
            polygon=keep,
            role=s.role,
            ref=s.ref,
            source_axis=s.source_axis,
            altitude=s.altitude,
            altitude_high=s.altitude_high,
            altitude_low=s.altitude_low,
            node_altitudes=None,
            is_bridge=s.is_bridge,
        )
        # Resample node_altitudes via edge interpolation against the
        # ORIGINAL bridge ring's per-vertex altitudes.
        if old_alts is not None and old_open:
            new_alts = _resample_node_altitudes_nn(
                keep, old_open, old_alts)
            if new_alts is not None:
                new_s.node_altitudes = new_alts
        new_shapes.append(new_s)
        n_modified += 1

    layout.shapes = new_shapes
    return n_modified


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
        """Return DEM at (x, y) clamped UP toward the nearest runway
        when within ``runway_clamp_radius_m`` and DEM dips below
        ``runway_e - g·d``, else raw DEM, else None.

        Per user 2026-05-11: the clamp is ASYMMETRIC.  We only ever
        pull the boundary UP toward the runway (the original
        "graded up to runway elevation" rule for low terrain near
        the runway).  We never pull the boundary DOWN — if the
        surrounding terrain is higher than the runway-band, the
        boundary follows DEM so Ortho4XP's
        ``smooth_raster_over_airports`` doesn't drag the rendered
        terrain down into a 20 m canyon around the airport
        perimeter.
        """
        try:
            lat, lon = m_to_ll(x, y)
            dem_e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except _GEOM_EXC:
            dem_e = None
        # Find nearest runway and its elevation at the nearest point.
        best_d = float('inf')
        best_e = None
        pt = _Point(x, y)
        for s in runway_shapes:
            try:
                d = s.polygon.distance(pt)
            except _GEOM_EXC:
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
            except _GEOM_EXC:
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
        if dem_e is None:
            # No DEM available — fall back to the floor (the
            # closest the boundary can be to the runway at this
            # distance without violating the grade cap).
            return lo
        # Asymmetric clamp: only pull UP toward runway.  If DEM is
        # below the floor, lift it; otherwise follow DEM.
        if dem_e < lo:
            return lo
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

    # Per user 2026-05-12: emit the boundary as a CHAIN OF 4-corner
    # rectangles (one per densified boundary segment) instead of a
    # single buffered strip polygon.  Each rect is either flat
    # (single ``altitude=`` tag) or sloped (``altitude_high`` /
    # ``altitude_low`` with the [high, low, low, high] corner
    # convention), so debug tools like JOSM can read the altitude
    # profile along the perimeter directly off each rect's tags.
    boundary_geom = layout.airport_boundary
    if boundary_geom.geom_type == "Polygon":
        ext_rings = [boundary_geom.exterior]
    elif boundary_geom.geom_type == "MultiPolygon":
        ext_rings = [g.exterior for g in boundary_geom.geoms]
    else:
        return 0
    pavement_polys = [
        s.polygon for s in layout.shapes
        if s.polygon is not None
        and not s.polygon.is_empty
        and s.role != ROLE_BOUNDARY]
    pavement_union: Optional[Polygon] = None
    if pavement_polys:
        try:
            pavement_union = unary_union(pavement_polys)
        except _GEOM_EXC:
            pavement_union = None

    def _rect_for_segment(
            p0: Tuple[float, float],
            p1: Tuple[float, float],
            alt0: float, alt1: float,
            ) -> Optional[Tuple[Polygon, Optional[float], float]]:
        """Build a 4-corner rect spanning the boundary segment
        p0 → p1 with strip half-width.  Convention: corners 0, 3 at
        the HIGH-altitude end, corners 1, 2 at the LOW end (matches
        runway segment emit).  Returns
        ``(polygon, altitude_high, altitude_low)`` with
        ``altitude_high=None`` for flat segments (``altitude_low``
        carries the single flat value in that case).
        """
        # Order so p0 is the HIGH end (alt0 >= alt1).
        if abs(alt0 - alt1) < 0.1:
            eh: Optional[float] = None
            el = round((alt0 + alt1) / 2.0, 1)
        elif alt0 >= alt1:
            eh = round(alt0, 1)
            el = round(alt1, 1)
        else:
            p0, p1 = p1, p0
            alt0, alt1 = alt1, alt0
            eh = round(alt0, 1)
            el = round(alt1, 1)
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        L = math.hypot(dx, dy)
        if L < 0.5:
            return None
        # Perpendicular unit vector × half-width.  Sign matches
        # ``pavement.runway_geometry.runway_corners`` so adjacent
        # rects don't accidentally flip ring orientation.
        perp_x = -dy / L * strip_half_width_m
        perp_y = dx / L * strip_half_width_m
        corners = [
            (p0[0] + perp_x, p0[1] + perp_y),  # 0 high-left
            (p1[0] + perp_x, p1[1] + perp_y),  # 1 low-left
            (p1[0] - perp_x, p1[1] - perp_y),  # 2 low-right
            (p0[0] - perp_x, p0[1] - perp_y),  # 3 high-right
        ]
        try:
            poly = _Polygon(corners)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if poly.is_empty or poly.geom_type != "Polygon":
                return None
        except _GEOM_EXC:
            return None
        return poly, eh, el

    n_emitted = 0
    for ring in ext_rings:
        ring_coords = list(ring.coords)
        if ring_coords and ring_coords[0] == ring_coords[-1]:
            ring_coords = ring_coords[:-1]
        if len(ring_coords) < 3:
            continue
        # Densify to ``densify_step_m`` along the ring.  The closing
        # duplicate is added back at the end of ``_densify_ring``.
        dense = _densify_ring(ring_coords)
        if len(dense) < 4:
            continue
        # Walk consecutive pairs; emit a rect per pair.
        n_pairs = len(dense) - 1
        for i in range(n_pairs):
            p0 = dense[i]
            p1 = dense[i + 1]
            a0 = _runway_clamped_alt(p0[0], p0[1])
            a1 = _runway_clamped_alt(p1[0], p1[1])
            if a0 is None:
                a0 = 0.0
            if a1 is None:
                a1 = 0.0
            built = _rect_for_segment(p0, p1, float(a0), float(a1))
            if built is None:
                continue
            poly, eh, el = built
            # Skip rects entirely buried inside pavement — they
            # would just shadow runway / taxi / apron geometry and
            # fail the no-self-overlap test.  Partial overlaps are
            # OK; X-Plane resolves at render time and the rect
            # still labels its segment.
            if (pavement_union is not None
                    and not pavement_union.is_empty):
                try:
                    if pavement_union.contains(poly):
                        continue
                    # If pavement covers >80 % of the rect, skip too
                    # — keeps the chain coherent with what's
                    # actually visible.
                    inter = pavement_union.intersection(poly)
                    if (not inter.is_empty
                            and inter.area > 0.8 * poly.area):
                        continue
                    # Otherwise trim against pavement; if the
                    # trimmed result is still a Polygon, replace.
                    trimmed = poly.difference(pavement_union)
                    if (not trimmed.is_empty
                            and trimmed.geom_type == "Polygon"):
                        poly = trimmed
                except _GEOM_EXC:
                    pass
            shape = BuiltShape(
                polygon=poly,
                role=ROLE_BOUNDARY,
                ref="airport_boundary",
            )
            if eh is None:
                shape.altitude = el
            else:
                shape.altitude_high = eh
                shape.altitude_low = el
            layout.shapes.append(shape)
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
        except _GEOM_EXC:
            dem_e = None
        best_d = float('inf')
        best_e = None
        from shapely.ops import nearest_points as _np
        pt = _Point(x, y)
        for s in runway_shapes:
            try:
                d = s.polygon.distance(pt)
            except _GEOM_EXC:
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
            except _GEOM_EXC:
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
        if dem_e is None:
            return lo
        # Asymmetric clamp (user 2026-05-11): only pull UP toward
        # runway when DEM is below the floor; never pull DOWN.
        # See ``_runway_clamped_alt`` in ``_emit_airport_boundary_shape``
        # for the full rationale.
        if dem_e < lo:
            return lo
        return dem_e

    def _dem_alt(x: float, y: float) -> Optional[float]:
        try:
            lat, lon = m_to_ll(x, y)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except _GEOM_EXC:
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
        except _GEOM_EXC:
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
        except _GEOM_EXC:
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
        except _GEOM_EXC:
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
        except _GEOM_EXC:
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
                except _GEOM_EXC:
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

        # ── User's sequential-walk algorithm (2026-05-16) ─────────
        # 1. Walk boundary; mark vertices within ``runway_clamp_
        #    radius_m`` of any runway (B19's clamp radius is the
        #    upstream condition that creates altitude gaps).
        # 2. Each maximal contiguous "marked" stretch is a bridge
        #    run: the bridge's outer edge walks those vertices.
        # 3. For each run, snap from the run's end-vertex across to
        #    the nearest pavement_union outer-ring vertex; then walk
        #    pavement_union BACK toward the run start, collecting
        #    canonical pavement vertices along the way.  Close the
        #    polygon by snapping from the last pavement-walk vertex
        #    to the run start.  By construction every bridge node
        #    is either a boundary node (outer side) or a pavement-
        #    union outer-ring node (inner side); no synthesised
        #    intersection vertices; no overlap because both walks
        #    are monotonic along their respective polygon
        #    perimeters.
        from shapely.ops import nearest_points

        # Pre-build pavement_union outer ring (canonical nodes).
        # Include runways so the bridge inner edge wraps around them
        # (touching, not overlapping).
        pav_for_inner = [
            s.polygon for s in layout.shapes
            if s.role != ROLE_BOUNDARY
            and s.polygon is not None
            and not s.polygon.is_empty]
        pav_union_local: Optional[Polygon] = None
        pav_ring_coords: List[Tuple[float, float]] = []
        pav_ring_line: Optional[LineString] = None
        if pav_for_inner:
            try:
                pav_union_local = unary_union(pav_for_inner)
                if pav_union_local.geom_type == "Polygon":
                    rc = list(pav_union_local.exterior.coords)
                    if rc and rc[0] == rc[-1]:
                        rc = rc[:-1]
                    pav_ring_coords = [(float(x), float(y))
                                        for (x, y) in rc]
                elif pav_union_local.geom_type == "MultiPolygon":
                    # Pick the largest component — bridges face the
                    # main pavement mass; small disconnected pieces
                    # aren't bridged.
                    largest = max(
                        pav_union_local.geoms,
                        key=lambda g: g.area)
                    rc = list(largest.exterior.coords)
                    if rc and rc[0] == rc[-1]:
                        rc = rc[:-1]
                    pav_ring_coords = [(float(x), float(y))
                                        for (x, y) in rc]
            except _GEOM_EXC:
                pav_ring_coords = []
        if len(pav_ring_coords) >= 3:
            try:
                pav_ring_line = _LS(pav_ring_coords + [pav_ring_coords[0]])
            except _GEOM_EXC:
                pav_ring_line = None

        # Altitude lookup for pav_ring nodes (round to 0.1 m).
        pav_alt_lookup: Dict[Tuple[int, int], float] = {}
        for (px, py, pa) in pav_edge_pts:
            k = (int(round(px * 10)), int(round(py * 10)))
            pav_alt_lookup[k] = float(pa)
        def _pav_alt(px: float, py: float) -> float:
            k = (int(round(px * 10)), int(round(py * 10)))
            if k in pav_alt_lookup:
                return pav_alt_lookup[k]
            # Synthesised by unary_union — fall back to DEM
            sd = _dem_alt(px, py)
            return float(sd) if sd is not None else 0.0

        for run in runs:
            if len(run) < 2:
                continue
            # Outer side of the bridge sits at the airport_boundary
            # ribbon's INNER edge — offset inward by
            # ``strip_half_width_m`` (2.5m) from the boundary line.
            # This places the bridge's outer vertices on the same
            # locus as the ribbon's interior-side nodes, eliminating
            # the 2.5m ribbon overlap that walking the raw boundary
            # would produce.
            STRIP_HALF_WIDTH_M = 2.5
            raw_outer_pts: List[Tuple[float, float]] = []
            raw_outer_alts: List[float] = []
            for ii_in_run, i_dense in enumerate(run):
                raw_outer_pts.append(
                    (per_vert[i_dense][0], per_vert[i_dense][1]))
                raw_outer_alts.append(
                    round(float(per_vert[i_dense][2]), 1))
            # Compute inward-perpendicular offset per vertex from
            # the local boundary tangent (average of the two
            # adjacent segments).  The two perpendiculars are
            # disambiguated by ``boundary_poly.contains()`` on a
            # short probe.
            outer_pts = []
            outer_alts = list(raw_outer_alts)
            n_raw = len(raw_outer_pts)
            for k, (bx, by) in enumerate(raw_outer_pts):
                if 0 < k < n_raw - 1:
                    prev_pt = raw_outer_pts[k - 1]
                    next_pt = raw_outer_pts[k + 1]
                elif k > 0:
                    prev_pt = raw_outer_pts[k - 1]
                    next_pt = (bx, by)
                else:
                    prev_pt = (bx, by)
                    next_pt = (raw_outer_pts[k + 1]
                                if n_raw > 1 else (bx, by))
                tx = next_pt[0] - prev_pt[0]
                ty = next_pt[1] - prev_pt[1]
                tmag = math.hypot(tx, ty)
                if tmag < 1e-6:
                    outer_pts.append((bx, by))
                    continue
                ux = tx / tmag
                uy = ty / tmag
                probe = _Point(bx + (-uy) * 0.5, by + ux * 0.5)
                if boundary_poly.contains(probe):
                    perp_x = -uy
                    perp_y = ux
                else:
                    perp_x = uy
                    perp_y = -ux
                outer_pts.append((bx + perp_x * STRIP_HALF_WIDTH_M,
                                   by + perp_y * STRIP_HALF_WIDTH_M))
            if len(outer_pts) < 2:
                continue

            # When no pavement is available, fall back to a 100m-
            # inward synthesised inner edge.
            if pav_ring_line is None or not pav_ring_coords:
                # Build an inward-perpendicular polyline at
                # ``bridge_depth_m`` for the inner edge.
                # (Same logic as previous v1 fallback.)
                ctr = boundary_poly.centroid
                inner_pts: List[Tuple[float, float]] = []
                inner_alts: List[float] = []
                for (bx, by) in outer_pts:
                    perp_x = ctr.x - bx
                    perp_y = ctr.y - by
                    pmag = math.hypot(perp_x, perp_y)
                    if pmag < 1e-6:
                        inner_pts.append((bx, by))
                        inner_alts.append(0.0)
                        continue
                    perp_x /= pmag
                    perp_y /= pmag
                    sx = bx + perp_x * bridge_depth_m
                    sy = by + perp_y * bridge_depth_m
                    sd = _dem_alt(sx, sy)
                    inner_pts.append((sx, sy))
                    inner_alts.append(
                        round(float(sd), 1)
                        if sd is not None else 0.0)
                ring_pts = list(outer_pts) + list(reversed(inner_pts))
                ring_alts = list(outer_alts) + list(reversed(inner_alts))
            else:
                # ── Step 1: snap run endpoints onto pav_union ring ──
                start_b = outer_pts[0]
                end_b = outer_pts[-1]
                # Snap to nearest pav_union outer-ring VERTEX
                # (not just nearest point — the snap target must be
                # canonical).
                def _nearest_pav_vertex(x: float, y: float
                                          ) -> Tuple[int, float]:
                    best_i = -1
                    best_d = float('inf')
                    for ii, (px, py) in enumerate(pav_ring_coords):
                        d = math.hypot(px - x, py - y)
                        if d < best_d:
                            best_d = d
                            best_i = ii
                    return best_i, best_d
                start_i, start_d = _nearest_pav_vertex(*start_b)
                end_i, end_d = _nearest_pav_vertex(*end_b)
                # Reject runs with no nearby pavement on either end.
                if (start_i < 0 or end_i < 0
                        or start_d > bridge_depth_m * 2
                        or end_d > bridge_depth_m * 2):
                    # Fall back to synthesised inner edge.
                    ctr = boundary_poly.centroid
                    inner_pts = []
                    inner_alts = []
                    for (bx, by) in outer_pts:
                        perp_x = ctr.x - bx
                        perp_y = ctr.y - by
                        pmag = math.hypot(perp_x, perp_y)
                        if pmag < 1e-6:
                            inner_pts.append((bx, by))
                            inner_alts.append(0.0)
                            continue
                        perp_x /= pmag
                        perp_y /= pmag
                        sx = bx + perp_x * bridge_depth_m
                        sy = by + perp_y * bridge_depth_m
                        sd = _dem_alt(sx, sy)
                        inner_pts.append((sx, sy))
                        inner_alts.append(
                            round(float(sd), 1)
                            if sd is not None else 0.0)
                    ring_pts = list(outer_pts) + list(reversed(inner_pts))
                    ring_alts = list(outer_alts) + list(reversed(inner_alts))
                else:
                    # ── Step 2: walk pav_ring from end_i back to start_i ──
                    # Pick the direction whose initial step from
                    # end_i is CLOSER to start_b than the opposite
                    # direction's initial step (so we walk "back
                    # toward the run start").
                    n_ring = len(pav_ring_coords)
                    fwd_first = pav_ring_coords[(end_i + 1) % n_ring]
                    bwd_first = pav_ring_coords[(end_i - 1) % n_ring]
                    d_fwd = math.hypot(fwd_first[0] - start_b[0],
                                        fwd_first[1] - start_b[1])
                    d_bwd = math.hypot(bwd_first[0] - start_b[0],
                                        bwd_first[1] - start_b[1])
                    step = +1 if d_fwd < d_bwd else -1
                    # Walk pav_union from end_i back toward start;
                    # STOP when current vertex is > 400 m from
                    # start_b (the user's closure rule).  This
                    # bounds the inner walk so a bridge run on one
                    # side of the airport doesn't trace around to
                    # the opposite side.
                    CLOSURE_DIST_M = runway_clamp_radius_m  # 400m
                    inner_pts = []
                    inner_alts = []
                    idx = end_i
                    visited = 0
                    while visited < n_ring:
                        px, py = pav_ring_coords[idx]
                        d_to_start = math.hypot(px - start_b[0],
                                                 py - start_b[1])
                        if d_to_start > CLOSURE_DIST_M:
                            break
                        inner_pts.append((float(px), float(py)))
                        inner_alts.append(round(_pav_alt(px, py), 1))
                        if idx == start_i:
                            break
                        idx = (idx + step) % n_ring
                        visited += 1
                    if len(inner_pts) < 2:
                        continue
                    # Close: last inner vertex → start_b via the
                    # implicit short edge of the polygon ring.
                    ring_pts = list(outer_pts) + list(inner_pts)
                    ring_alts = list(outer_alts) + list(inner_alts)

            if len(ring_pts) < 4:
                continue
            try:
                bridge_poly = _Polygon(ring_pts)
                if not bridge_poly.is_valid:
                    fixed = bridge_poly.buffer(0)
                    if (fixed.is_empty
                            or fixed.geom_type != "Polygon"):
                        continue
                    fc = list(fixed.exterior.coords)
                    if fc and fc[0] == fc[-1]:
                        fc = fc[:-1]
                    if len(fc) != len(ring_pts):
                        continue
                    bridge_poly = fixed
                if bridge_poly.is_empty:
                    continue
                if bridge_poly.geom_type != "Polygon":
                    continue
                if bridge_poly.area < 100.0:
                    continue
            except _GEOM_EXC:
                continue

            # Final cleanup: subtract NON-RUNWAY pavement +
            # ribbon to trim small residual overlaps from
            # per-vertex / per-segment perpendicular mismatch.
            # Runways are excluded (they need the 1m sloping-rect
            # buffer per B12 to avoid creating bridge vertices on
            # runway long edges, which the
            # ``test_no_vertex_on_sloping_rect_edge`` invariant
            # catches).  Runway-vs-bridge overlap is small after
            # the canonical-node construction and stays under the
            # overlap-baseline cap on its own.
            cleanup_subs: List[Polygon] = []
            non_runway_pav = [
                s.polygon for s in layout.shapes
                if s.role not in (ROLE_BOUNDARY, ROLE_RUNWAY)
                and s.polygon is not None
                and not s.polygon.is_empty]
            if non_runway_pav:
                try:
                    nr_union = unary_union(non_runway_pav)
                    if nr_union is not None and not nr_union.is_empty:
                        cleanup_subs.append(nr_union)
                except _GEOM_EXC:
                    pass
            if ribbon_union is not None and not ribbon_union.is_empty:
                cleanup_subs.append(ribbon_union)
            # Sloping-rect (runway / parallel / stub / cross-conn)
            # union buffered by 1m so any intersection vertices the
            # subtraction creates land OUTSIDE the
            # ``EDGE_PROX_M`` test tolerance.
            sloping_polys = [
                s.polygon for s in layout.shapes
                if s.role in (ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL,
                              ROLE_SECONDARY_PARALLEL, ROLE_STUB,
                              ROLE_CROSS_CONNECTOR)
                and s.polygon is not None
                and not s.polygon.is_empty]
            if sloping_polys:
                try:
                    sr_union = unary_union(sloping_polys)
                    if sr_union is not None and not sr_union.is_empty:
                        cleanup_subs.append(sr_union.buffer(1.0))
                except _GEOM_EXC:
                    pass
            for sub in cleanup_subs:
                try:
                    trimmed = bridge_poly.difference(sub)
                except _GEOM_EXC:
                    trimmed = None
                if trimmed is None or trimmed.is_empty:
                    continue
                if trimmed.geom_type == "Polygon":
                    bridge_poly = trimmed
                elif trimmed.geom_type == "MultiPolygon":
                    parts = sorted(trimmed.geoms,
                                   key=lambda g: -g.area)
                    if parts and parts[0].area >= 100.0:
                        bridge_poly = parts[0]
                    else:
                        bridge_poly = None
                        break
                else:
                    bridge_poly = None
                    break
            if bridge_poly is None or bridge_poly.is_empty:
                continue
            if bridge_poly.area < 100.0:
                continue
            # Resample altitudes for the (possibly reshaped) ring.
            new_coords = list(bridge_poly.exterior.coords)
            if new_coords and new_coords[0] == new_coords[-1]:
                new_coords_open = new_coords[:-1]
            else:
                new_coords_open = new_coords
            if len(new_coords_open) < 3:
                continue
            canon_alt: Dict[Tuple[int, int], float] = {}
            for (cx, cy), ca in zip(
                    list(outer_pts) + list(inner_pts),
                    list(outer_alts) + list(inner_alts)):
                ck = (int(round(cx * 10)), int(round(cy * 10)))
                canon_alt[ck] = ca
            ring_alts = []
            for (cx, cy) in new_coords_open:
                ck = (int(round(cx * 10)),
                      int(round(cy * 10)))
                ca = canon_alt.get(ck)
                if ca is None:
                    sd = _dem_alt(cx, cy)
                    ca = (round(float(sd), 1)
                          if sd is not None else 0.0)
                ring_alts.append(ca)

            node_alts = list(ring_alts) + [ring_alts[0]]
            layout.shapes.append(BuiltShape(
                polygon=bridge_poly,
                role=ROLE_BOUNDARY,
                ref="boundary_dem_bridge",
                node_altitudes=node_alts,
            ))
            n_emitted += 1

    return n_emitted


