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
            except _GEOM_EXC:
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
                except _GEOM_EXC:
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
                except _GEOM_EXC:
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
                except _GEOM_EXC:
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
            except _GEOM_EXC:
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
                except _GEOM_EXC:
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
                except _GEOM_EXC:
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
                except _GEOM_EXC:
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
                except _GEOM_EXC:
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
                except _GEOM_EXC:
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


