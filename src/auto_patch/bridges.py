"""Taxi/road bridges + tunnel portals + depressed-road segments.

Five emitters:

* ``_emit_tunnel_portals`` — short ramp polygons at tunnel-OSM
  entrance / exit so the road dives under airport pavement.
* ``_scenery_has_bridge_objects`` — DSF-pavement check that gates
  bridge emission (no point emitting bridge polygons if the
  underlying scenery already provides bridge meshes).
* ``_emit_taxi_bridges`` — taxi-over-road bridges.
* ``_emit_underpass_road_approaches`` — road approaches sloping
  down toward an underpass.
* ``_emit_through_airport_depressed_roads`` — road segments
  depressed below airport surface where they cut through.

**All five are gated off via ``EMIT_BRIDGES_AND_TUNNELS=False`` in
``O4_Pavement_Config``** pending a future re-enable pass — they
emit ``tunnel_ramp`` / ``retaining_wall`` polygons that overlap
junction / groundside_pavement geometry on the post-refactor
output (KPHX 75 K m² overlap, caught by ``test_no_self_overlap``).
The module structure is preserved so re-engaging is a single flag
flip in Config plus per-emitter overlap-clip work.

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _emit_tunnel_portals
    _scenery_has_bridge_objects
    _emit_taxi_bridges
    _emit_underpass_road_approaches
    _emit_through_airport_depressed_roads
"""
from __future__ import annotations

import math
import os
import re
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
    ROLE_TUNNEL_RAMP,
    SHARED_VERTEX_TOL_M,
)
from .pavement.vertices import _snap_polygon_vertices_to_rect_corners
from .pavement.runways import _sample_runway_segment_elev
from .elevation import _sample_dem


__all__ = [
    "_emit_taxi_bridges",
    "_emit_through_airport_depressed_roads",
    "_emit_tunnel_portals",
    "_emit_underpass_road_approaches",
    "_scenery_has_bridge_objects",
]


# Per-OSM-highway-type carriageway width (user 2026-05-03).
# Was a single 22 m default, which made every tunnel look like a
# 6-lane motorway.  Real-world widths vary by classification; the
# numbers below match typical FAA-relevant standards (single
# carriageway including shoulders).
HIGHWAY_CARRIAGEWAY_WIDTH_M = {
    "motorway":         24.0,  # 6+ lanes per direction in some places
    "motorway_link":     8.0,
    "trunk":            22.0,
    "trunk_link":        8.0,
    "primary":          18.0,
    "primary_link":      7.0,
    "secondary":        11.0,  # ~half of trunk per user 2026-05-03
    "secondary_link":    7.0,
    "tertiary":          9.0,
    "tertiary_link":     6.0,
    "residential":       7.0,
    "service":           6.0,
}


def _carriageway_width_for(highway_type: Optional[str],
                            default_m: float) -> float:
    """Return the carriageway width in metres for an OSM highway
    type, falling back to ``default_m`` for unknown types.
    """
    if highway_type is None:
        return default_m
    return HIGHWAY_CARRIAGEWAY_WIDTH_M.get(highway_type, default_m)


def _emit_tunnel_portals(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        tunnel_depth_m: float = 8.0,
        max_ramp_grade: float = 0.04,
        ramp_min_length_m: float = 200.0,
        arm_max_length_m: float = 500.0,
        carriageway_width_m: float = 22.0,
        retaining_wall_width_m: float = 1.0,
        # ``wall_gap_m`` must exceed the OSM emit's vertex bucket
        # size (SHARED_VERTEX_TOL_M = 0.5 m) so the ramp's road-edge
        # corners and the wall's inner corners don't hash to the
        # same node id.  At the portal end the ramp altitude is
        # apt_elev − tunnel_depth; the wall altitude is apt_elev.
        # Sharing the vertex would emit one node with two altitudes,
        # rendering as a vertical glitch (user 2026-05-03).
        wall_gap_m: float = 0.6,
        portal_cluster_dist_m: float = 15.0,
        boundary_clearance_m: float = 0.5,
        excluded_way_ids: Optional[set] = None,
        ) -> int:
    """For each tunnel portal (each end of an OSM ``aeroway=*``
    ``tunnel=yes|building_passage`` way), emit the visible road-
    surface structure that transitions outside-DEM elevation down
    to ``apt_elev − tunnel_depth_m`` at the portal:

      1. A flat ``ROLE_RETAINING_WALL`` CAP polygon AT the portal
         node (perpendicular to road direction at portal).  The
         cap's centre line is the portal node — its width spans
         the road (carriageway width + wall_gap on each side),
         its thickness is ``retaining_wall_width_m`` (1 m).
      2. Two flat ``ROLE_RETAINING_WALL`` ARM polygons reaching
         OUTWARD from the cap along the surface roadway, on each
         side of the carriageway.  The arms follow the OSM
         surface road's polyline — multi-segment when the road
         curves.  Arm length adapts to terrain: we first walk up
         to ``arm_max_length_m`` (default 500 m), sample DEM at
         the far end, then truncate the walk so the grade from
         ``apt_elev − tunnel_depth_m`` (portal) up to that DEM
         height never exceeds ``max_ramp_grade``.  Floor at
         ``ramp_min_length_m`` (default 200 m) so a flat road
         still gets a substantial visible approach.
      3. A chain of sloped ``ROLE_TUNNEL_RAMP`` polygons matching
         the surface-road segments under the arms.  Each segment's
         elevation interpolates linearly from
         ``apt_elev − tunnel_depth_m`` at the portal to outside
         DEM at the far end of the walk, in proportion to its
         cumulative distance from the portal.

    Per user 2026-04-29 (SPJC review): the previous geometry put
    the cap PAST the portal INSIDE the airport, with arms going
    AWAY from the tunnel only as far as the OSM way extended
    outside the airport.  That made SPJC's SW tunnel "too short
    and inside the tunnel" because the OSM way ended right at the
    boundary.  Walking the SURFACE road OUTWARD from the OSM
    portal node correctly lands the cap at the tunnel mouth and
    the arms on the highway approach.

    Two-carriageway tunnels (divided highways) cluster by portal-
    node proximity (within ``portal_cluster_dist_m``).  Each
    carriageway in a cluster gets its own arm pair; the caps form
    a perpendicular line across all member portals.

    Boundary coordination: subtract the tunnel-polygon union
    (buffered by ``boundary_clearance_m``, default 1 m which
    exceeds the OSM-emit vertex bucket size of 0.5 m so boundary
    nodes never collapse onto wall nodes) from every
    ``ROLE_BOUNDARY`` shape.

    Returns the number of tunnel PORTALS emitted (each contributing
    1 cap + 2 arm walls + a ramp chain).
    """
    from .pipeline import _load_osm_airports, _load_osm_big_roads
    # Load big-roads OSM cache for this tile.
    nodes_r, ways_r = _load_osm_big_roads(
        layout.anchor[0], layout.anchor[1])
    if not ways_r:
        return 0
    # Project nodes to meter space.
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)

    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))
    nodes_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_r.items():
        nodes_m[nid] = _to_m(lon, lat)
    HW_TUNNEL_TYPES = {
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "motorway_link", "trunk_link",
        "primary_link", "residential", "service",
    }
    TUNNEL_VALUES = {"yes", "building_passage"}
    # Build node-to-way and way-by-id indices for surface-road walking.
    way_by_id: Dict[str, Tuple[List[str], Dict[str, str]]] = {}
    node_to_ways: Dict[str, List[str]] = {}
    for wid, nrefs, tags in ways_r:
        way_by_id[wid] = (nrefs, tags)
        for n in nrefs:
            node_to_ways.setdefault(n, []).append(wid)
    # We walk a generous maximum, then truncate per-portal based
    # on actual DEM at the far end so the resulting grade never
    # exceeds ``max_ramp_grade``.  The truncated length is at
    # least ``ramp_min_length_m`` so even flat-road portals get
    # a recognisable approach.  The planning grade is reduced by
    # 0.005 to leave headroom for the 0.1 m altitude rounding —
    # without it, short segments (e.g. 9 m) could round up to
    # ~4.4 % when the design grade is exactly 4 %.
    grade_safety_margin = 0.005
    plan_grade = max(max_ramp_grade - grade_safety_margin, 1e-3)
    arm_walk_max_m = max(arm_max_length_m,
                         ramp_min_length_m,
                         tunnel_depth_m / plan_grade)
    # Helper: airport surface elevation at (cx, cy).  Use the
    # boundary-ribbon ``node_altitudes`` (CIFP-anchored, grade-
    # clamped) when a vertex is nearby, else fall back to DEM.
    def _airport_elevation_at(cx: float, cy: float) -> Optional[float]:
        best_d = float('inf')
        best_alt: Optional[float] = None
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                continue
            if s.ref != "airport_boundary":
                continue
            if not s.node_altitudes:
                continue
            try:
                rcoords = list(s.polygon.exterior.coords)
            except Exception:
                continue
            if rcoords and rcoords[0] == rcoords[-1]:
                rcoords = rcoords[:-1]
            for k, (vx, vy) in enumerate(rcoords):
                if k >= len(s.node_altitudes):
                    break
                d = math.hypot(vx - cx, vy - cy)
                if d < best_d:
                    best_d = d
                    best_alt = s.node_altitudes[k]
        if best_alt is not None and best_d <= 200.0:
            return float(best_alt)
        # Fall back to DEM at the point.
        try:
            lat, lon = _m_to_ll(cx, cy)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            return None
    # Helper: from a portal node, walk the connecting non-tunnel
    # surface road OUTWARD for ``length_m`` metres.  Returns the
    # walked path as a list of (x, y) points starting at the
    # portal, or None if no valid surface way connects.
    def _walk_surface(portal_nid: str,
                      tunnel_wid: str,
                      length_m: float
                      ) -> Optional[List[Tuple[float, float]]]:
        if portal_nid not in nodes_m:
            return None
        # Look at every way that shares this node.  Prefer surface
        # roads (highway tag, not tunnel-tagged); skip the tunnel
        # way itself.  When multiple candidates exist (e.g., a
        # roundabout), pick the first valid one.
        for other_wid in node_to_ways.get(portal_nid, []):
            if other_wid == tunnel_wid:
                continue
            if other_wid not in way_by_id:
                continue
            o_nrefs, o_tags = way_by_id[other_wid]
            if o_tags.get("tunnel") in TUNNEL_VALUES:
                continue
            if o_tags.get("highway") not in HW_TUNNEL_TYPES:
                continue
            try:
                idx = o_nrefs.index(portal_nid)
            except ValueError:
                continue
            # Decide which direction to walk (away from tunnel).
            forward = o_nrefs[idx:]
            backward = list(reversed(o_nrefs[:idx + 1]))
            # If portal is at one end, walk the other way.  If
            # mid-way, walk the longer side.
            walk_refs: List[str]
            if idx == 0:
                walk_refs = forward
            elif idx == len(o_nrefs) - 1:
                walk_refs = backward
            else:
                fl = sum(
                    math.hypot(
                        nodes_m[forward[i + 1]][0]
                        - nodes_m[forward[i]][0],
                        nodes_m[forward[i + 1]][1]
                        - nodes_m[forward[i]][1])
                    for i in range(len(forward) - 1)
                    if (forward[i] in nodes_m
                        and forward[i + 1] in nodes_m))
                bl = sum(
                    math.hypot(
                        nodes_m[backward[i + 1]][0]
                        - nodes_m[backward[i]][0],
                        nodes_m[backward[i + 1]][1]
                        - nodes_m[backward[i]][1])
                    for i in range(len(backward) - 1)
                    if (backward[i] in nodes_m
                        and backward[i + 1] in nodes_m))
                walk_refs = forward if fl >= bl else backward
            # Collect walk points up to length_m.
            pts: List[Tuple[float, float]] = []
            cum = 0.0
            for ni, n in enumerate(walk_refs):
                if n not in nodes_m:
                    break
                p = nodes_m[n]
                if pts:
                    seg_len = math.hypot(
                        p[0] - pts[-1][0], p[1] - pts[-1][1])
                    if cum + seg_len >= length_m:
                        # Truncate the last segment to hit length_m
                        if seg_len > 0:
                            t = (length_m - cum) / seg_len
                            tx = pts[-1][0] + t * (p[0] - pts[-1][0])
                            ty = pts[-1][1] + t * (p[1] - pts[-1][1])
                            pts.append((tx, ty))
                        cum = length_m
                        break
                    cum += seg_len
                pts.append(p)
            if len(pts) >= 2 and cum > 5.0:
                return pts
            # Surface way too short — could chain to next way at
            # the far end, but for simplicity we accept short walks
            # if we got at least 5 m.
            if len(pts) >= 2:
                return pts
        return None
    # Collect portal data: (portal_node_id, tunnel_wid, walk_pts,
    # hw_type, apt_elev_at_portal, dem_at_far_end).
    portal_data: List[Tuple[str, str, List[Tuple[float, float]],
                              str, float, float]] = []
    excluded = excluded_way_ids or set()
    for tw_id, t_nrefs, t_tags in ways_r:
        if t_tags.get("tunnel") not in TUNNEL_VALUES:
            continue
        hw = t_tags.get("highway")
        if hw not in HW_TUNNEL_TYPES:
            continue
        if len(t_nrefs) < 2:
            continue
        # Skip OSM way IDs already handled by the through-
        # airport depressed-road emit (which produces a single
        # uniform depression instead of per-bridge ramps).
        if tw_id in excluded:
            continue
        for portal_idx in (0, len(t_nrefs) - 1):
            portal_nid = t_nrefs[portal_idx]
            if portal_nid not in nodes_m:
                continue
            walk = _walk_surface(portal_nid, tw_id, arm_walk_max_m)
            if walk is None or len(walk) < 2:
                continue
            # Merge very short consecutive segments so altitude
            # rounding to 0.1 m can't push the per-segment grade
            # above ``max_ramp_grade``.  At a 0.05 m worst-case
            # round-up, a 15 m segment rounds to ≤ 0.33 %/error.
            min_segment_m = 15.0
            merged: List[Tuple[float, float]] = [walk[0]]
            for k in range(1, len(walk)):
                d = math.hypot(walk[k][0] - merged[-1][0],
                               walk[k][1] - merged[-1][1])
                if d < min_segment_m and k != len(walk) - 1:
                    continue
                merged.append(walk[k])
            walk = merged
            if len(walk) < 2:
                continue
            # Densify long segments so the visible ramp tracks the
            # road with multiple sloped pieces — user 2026-05-03
            # ("SW tunnel only has one ramp segment, should be
            # multiple following the road up to DEM elevation").
            # Sparse OSM ways often have ~150-200 m gaps between
            # nodes; without densification a 200 m approach renders
            # as a single straight ramp.  Target ~50 m segments.
            target_seg_m = 50.0
            densified: List[Tuple[float, float]] = [walk[0]]
            for k in range(1, len(walk)):
                px, py = densified[-1]
                qx, qy = walk[k]
                d = math.hypot(qx - px, qy - py)
                n_sub = max(1, int(round(d / target_seg_m)))
                for s in range(1, n_sub + 1):
                    t = s / n_sub
                    densified.append(
                        (px + t * (qx - px), py + t * (qy - py)))
            walk = densified
            portal_xy = walk[0]
            apt_elev = _airport_elevation_at(*portal_xy)
            if apt_elev is None:
                continue
            elev_low = apt_elev - tunnel_depth_m
            # Find the truncation point that keeps grade ≤
            # max_ramp_grade.  Walk along the polyline summing
            # cumulative distance; sample DEM at each vertex; the
            # required-length-so-far is (DEM - elev_low) /
            # max_ramp_grade.  Stop walking once the actual walk
            # length matches or exceeds that requirement (i.e. the
            # grade from portal to here is ≤ max_ramp_grade).
            cum = 0.0
            kept_pts: List[Tuple[float, float]] = [walk[0]]
            grade_ok_at: float = 0.0  # cum dist where grade is OK
            for i in range(1, len(walk)):
                seg_len = math.hypot(
                    walk[i][0] - walk[i - 1][0],
                    walk[i][1] - walk[i - 1][1])
                cum += seg_len
                kept_pts.append(walk[i])
                try:
                    plat, plon = _m_to_ll(*walk[i])
                    dem_h = _sample_dem(
                        dem, tile_lat, tile_lon, plat, plon)
                except Exception:
                    dem_h = None
                if dem_h is None:
                    continue
                drop = float(dem_h) - elev_low
                req = (drop / plan_grade if drop > 0 else 0.0)
                if cum >= req and cum >= ramp_min_length_m:
                    grade_ok_at = cum
                    break
                grade_ok_at = cum
            # Truncate walk to the OK length (or to last vertex
            # reached if we never satisfied the grade — best
            # effort, the resulting ramp will be the gentlest
            # achievable on the available roadway).
            walk = kept_pts
            far_xy = walk[-1]
            try:
                far_lat, far_lon = _m_to_ll(*far_xy)
                far_dem = _sample_dem(
                    dem, tile_lat, tile_lon, far_lat, far_lon)
            except Exception:
                far_dem = None
            if far_dem is None:
                far_dem = apt_elev
            # If the resulting ramp would still violate grade
            # (DEM very high, road too short), cap far_dem at
            # elev_low + grade*length so we still emit something
            # sensible.  The visible top edge then sits LOWER
            # than DEM (subtle terrain dip) — preferable to a
            # spike or a mid-tunnel cap.
            max_drop = plan_grade * grade_ok_at
            if (far_dem - elev_low) > max_drop:
                far_dem = elev_low + max_drop
            portal_data.append(
                (portal_nid, tw_id, walk, hw,
                 float(apt_elev), float(far_dem)))
    if not portal_data:
        return 0
    # Cluster portals by node-coord proximity (divided highways).
    clusters: List[List[int]] = []
    used: set = set()
    for i in range(len(portal_data)):
        if i in used:
            continue
        nid_i = portal_data[i][0]
        cl = [i]
        used.add(i)
        pi = nodes_m[nid_i]
        for j in range(i + 1, len(portal_data)):
            if j in used:
                continue
            pj = nodes_m[portal_data[j][0]]
            if (math.hypot(pi[0] - pj[0], pi[1] - pj[1])
                    < portal_cluster_dist_m):
                cl.append(j)
                used.add(j)
        clusters.append(cl)
    # Per-cluster: build cap + arm walls + ramp chain.
    exclusion_zones: List[Polygon] = []
    n_emitted = 0
    half_wall_w = retaining_wall_width_m / 2.0
    for cl in clusters:
        # All portals in cluster share approximately the same
        # location.  Use the first portal's walk as the canonical
        # arm path; combine widths for divided highways.
        head = portal_data[cl[0]]
        portal_nid, _wid_unused, walk_pts, hw_type, apt_elev, far_dem = head
        if len(walk_pts) < 2:
            continue
        # Per-OSM-highway-type carriageway width (user 2026-05-03):
        # was a fixed 22 m default; now varies by classification
        # so secondary tunnels are ~half the width of trunk tunnels.
        carriage_w = _carriageway_width_for(
            hw_type, carriageway_width_m)
        half_carriage = 0.5 * carriage_w
        elev_low = apt_elev - tunnel_depth_m
        elev_high = far_dem
        # Compute the walk's cumulative distance for elevation
        # interpolation.
        cum_dists = [0.0]
        for i in range(1, len(walk_pts)):
            cum_dists.append(cum_dists[-1] + math.hypot(
                walk_pts[i][0] - walk_pts[i - 1][0],
                walk_pts[i][1] - walk_pts[i - 1][1]))
        total_walk = cum_dists[-1]
        if total_walk < 5.0:
            continue
        # Cluster spread for combined width: project each cluster
        # member's portal node onto the perpendicular at the head
        # portal.  The cap (and only the cap) is centred on the
        # cluster centroid, not on the head portal — user 2026-05-03
        # ("trunk highway tunnels not centered on OSM ways, offset
        # with one edge on one of the ways").  The arms still follow
        # the head walk because we don't yet emit per-carriageway
        # arms; centring at least the cap puts it symmetric across
        # both carriageways of a divided highway.
        first_seg = (walk_pts[1][0] - walk_pts[0][0],
                     walk_pts[1][1] - walk_pts[0][1])
        first_len = math.hypot(*first_seg)
        if first_len < 0.1:
            continue
        first_dir = (first_seg[0] / first_len,
                     first_seg[1] / first_len)
        first_perp = (-first_dir[1], first_dir[0])
        spans = []
        for k in cl:
            ni = portal_data[k][0]
            if ni not in nodes_m:
                continue
            p = nodes_m[ni]
            spans.append(
                (p[0] - walk_pts[0][0]) * first_perp[0]
                + (p[1] - walk_pts[0][1]) * first_perp[1])
        cluster_span = max(spans) - min(spans) if spans else 0.0
        cluster_perp_offset = (
            (max(spans) + min(spans)) / 2.0 if spans else 0.0)
        combined_half = half_carriage + 0.5 * cluster_span

        def _build_wall_segment(p_a: Tuple[float, float],
                                 p_b: Tuple[float, float],
                                 perp_off: float
                                 ) -> Optional[Polygon]:
            """4-corner wall polygon parallel to segment ``a-b``,
            offset by ``perp_off`` from the segment's centre line,
            ``retaining_wall_width_m`` thick."""
            seg = (p_b[0] - p_a[0], p_b[1] - p_a[1])
            slen = math.hypot(*seg)
            if slen < 0.1:
                return None
            ux, uy = seg[0] / slen, seg[1] / slen
            nx, ny = -uy, ux
            # Sign of perp_off picks which side.  Inner edge of
            # wall is at perp_off, outer edge at perp_off ± width.
            inner = perp_off
            outer = (perp_off + half_wall_w * 2.0
                     if perp_off >= 0
                     else perp_off - half_wall_w * 2.0)
            corners = [
                (p_a[0] + nx * inner, p_a[1] + ny * inner),
                (p_b[0] + nx * inner, p_b[1] + ny * inner),
                (p_b[0] + nx * outer, p_b[1] + ny * outer),
                (p_a[0] + nx * outer, p_a[1] + ny * outer),
            ]
            try:
                p = Polygon(corners)
                if not p.is_valid:
                    p = p.buffer(0)
                if p.geom_type == "Polygon" and not p.is_empty:
                    return p
            except Exception:
                return None
            return None
        # 1) Cap wall AT the portal cluster's centroid, perpendicular
        #    to the first segment.  The cap's centre line passes
        #    through the cluster centroid (so divided-highway
        #    tunnels are centered between the carriageways, user
        #    2026-05-03), its width spans the combined carriageways
        #    + 2 × wall_gap, its thickness is
        #    retaining_wall_width_m.
        cap_half_len = combined_half + wall_gap_m
        portal_xy = walk_pts[0]
        cap_centre = (
            portal_xy[0] + first_perp[0] * cluster_perp_offset,
            portal_xy[1] + first_perp[1] * cluster_perp_offset)
        c0 = (cap_centre[0] + first_perp[0] * cap_half_len,
              cap_centre[1] + first_perp[1] * cap_half_len)
        c1 = (cap_centre[0] - first_perp[0] * cap_half_len,
              cap_centre[1] - first_perp[1] * cap_half_len)
        # Move cap thickness INTO the tunnel direction (negative
        # first_dir) — cap occupies the strip from portal back
        # by retaining_wall_width_m.
        c0_back = (c0[0] - first_dir[0] * retaining_wall_width_m,
                   c0[1] - first_dir[1] * retaining_wall_width_m)
        c1_back = (c1[0] - first_dir[0] * retaining_wall_width_m,
                   c1[1] - first_dir[1] * retaining_wall_width_m)
        try:
            cap_poly = Polygon([c0, c1, c1_back, c0_back])
            if not cap_poly.is_valid:
                cap_poly = cap_poly.buffer(0)
            if (cap_poly.geom_type == "Polygon"
                    and not cap_poly.is_empty):
                layout.shapes.append(BuiltShape(
                    polygon=cap_poly,
                    role=ROLE_RETAINING_WALL,
                    ref="tunnel_cap",
                    altitude=round(apt_elev, 1)))
                exclusion_zones.append(cap_poly)
        except Exception:
            pass
        # 2) Arm walls + 3) Ramp polygons — one per walk segment.
        # Pre-compute per-vertex offset corners using the bisector
        # of adjacent segments at interior bends.  This makes
        # consecutive segments share their boundary vertices
        # exactly — no overlap, no gap.
        arm_off = combined_half + wall_gap_m + half_wall_w
        n_w = len(walk_pts)
        # Per-vertex perpendicular direction (unit vector).
        verts_perp: List[Tuple[float, float]] = []
        verts_scale: List[float] = []  # extension scale (1/cos(θ/2))
        for i in range(n_w):
            if i == 0:
                s = (walk_pts[1][0] - walk_pts[0][0],
                     walk_pts[1][1] - walk_pts[0][1])
                sl = math.hypot(*s)
                verts_perp.append((-s[1] / sl, s[0] / sl))
                verts_scale.append(1.0)
            elif i == n_w - 1:
                s = (walk_pts[i][0] - walk_pts[i - 1][0],
                     walk_pts[i][1] - walk_pts[i - 1][1])
                sl = math.hypot(*s)
                verts_perp.append((-s[1] / sl, s[0] / sl))
                verts_scale.append(1.0)
            else:
                s1 = (walk_pts[i][0] - walk_pts[i - 1][0],
                      walk_pts[i][1] - walk_pts[i - 1][1])
                s2 = (walk_pts[i + 1][0] - walk_pts[i][0],
                      walk_pts[i + 1][1] - walk_pts[i][1])
                l1 = math.hypot(*s1)
                l2 = math.hypot(*s2)
                u1 = (s1[0] / l1, s1[1] / l1)
                u2 = (s2[0] / l2, s2[1] / l2)
                avg = ((u1[0] + u2[0]) / 2.0,
                       (u1[1] + u2[1]) / 2.0)
                al = math.hypot(*avg)
                if al < 1e-6:
                    # Near-180° doubleback — use first segment's
                    # perpendicular and scale 1.
                    verts_perp.append((-u1[1], u1[0]))
                    verts_scale.append(1.0)
                    continue
                tangent = (avg[0] / al, avg[1] / al)
                perp = (-tangent[1], tangent[0])
                # Adjust offset to compensate for bend angle.
                # cos(θ/2) ≈ sqrt((1 + u1·u2) / 2).
                dot = u1[0] * u2[0] + u1[1] * u2[1]
                cos_half = max(0.1, math.sqrt(
                    max(0.0, (1.0 + dot) / 2.0)))
                verts_perp.append(perp)
                verts_scale.append(1.0 / cos_half)

        def _vertex_offset(idx: int, off: float
                           ) -> Tuple[float, float]:
            px, py = walk_pts[idx]
            nx, ny = verts_perp[idx]
            scaled = off * verts_scale[idx]
            return (px + nx * scaled, py + ny * scaled)

        for i in range(n_w - 1):
            p_a = walk_pts[i]
            p_b = walk_pts[i + 1]
            d_a = cum_dists[i]
            d_b = cum_dists[i + 1]
            seg_len = d_b - d_a
            if seg_len < 0.5:
                continue
            frac_a = d_a / total_walk if total_walk > 0 else 0.0
            frac_b = d_b / total_walk if total_walk > 0 else 0.0
            e_a = (1 - frac_a) * elev_low + frac_a * elev_high
            e_b = (1 - frac_b) * elev_low + frac_b * elev_high
            # Arm walls (one per side).  Inner edge at +/- (arm_off
            # − half_wall_w); outer edge at +/- (arm_off +
            # half_wall_w).  Using the per-vertex bisector so
            # adjacent segments share their join.
            for sign in (+1, -1):
                inner = sign * (arm_off - half_wall_w)
                outer = sign * (arm_off + half_wall_w)
                ai = _vertex_offset(i, inner)
                bi = _vertex_offset(i + 1, inner)
                bo = _vertex_offset(i + 1, outer)
                ao = _vertex_offset(i, outer)
                try:
                    wp = Polygon([ai, bi, bo, ao])
                    if not wp.is_valid:
                        wp = wp.buffer(0)
                    if (wp.geom_type == "Polygon"
                            and not wp.is_empty
                            and wp.area > 0.5):
                        layout.shapes.append(BuiltShape(
                            polygon=wp,
                            role=ROLE_RETAINING_WALL,
                            ref="tunnel_wall",
                            altitude=round(apt_elev, 1)))
                        exclusion_zones.append(wp)
                except Exception:
                    continue
            # Ramp polygon (single segment, sloped).  Corners
            # share with adjacent segments via verts_perp.
            ra = _vertex_offset(i, +combined_half)
            rb = _vertex_offset(i + 1, +combined_half)
            rc = _vertex_offset(i + 1, -combined_half)
            rd = _vertex_offset(i, -combined_half)
            # Rect convention (see _sample_runway_segment_elev):
            # corners [0, 3] are the HIGH-elevation short edge
            # (across the road at the high end), corners [1, 2]
            # are the LOW-elevation short edge.  ra/rb sit on
            # the +side, rd/rc on the -side; ra/rd are at walk[i]
            # and rb/rc at walk[i+1].  Order corners so the slope
            # axis runs ALONG the road (i ↔ i+1), not across it.
            if e_b >= e_a:
                # walk[i+1] = HIGH end → corners 0,3 at i+1
                ramp_corners = [rb, ra, rd, rc]
                eh, el = e_b, e_a
            else:
                # walk[i] = HIGH end → corners 0,3 at i
                ramp_corners = [ra, rb, rc, rd]
                eh, el = e_a, e_b
            try:
                rp = Polygon(ramp_corners)
                if not rp.is_valid:
                    rp = rp.buffer(0)
                if (rp.geom_type == "Polygon"
                        and not rp.is_empty
                        and rp.area > 0.5):
                    if abs(eh - el) >= 0.1:
                        layout.shapes.append(BuiltShape(
                            polygon=rp,
                            role=ROLE_TUNNEL_RAMP,
                            ref="tunnel_ramp",
                            altitude_high=round(eh, 1),
                            altitude_low=round(el, 1)))
                    else:
                        layout.shapes.append(BuiltShape(
                            polygon=rp,
                            role=ROLE_TUNNEL_RAMP,
                            ref="tunnel_ramp",
                            altitude=round(
                                0.5 * (eh + el), 1)))
                    exclusion_zones.append(rp)
            except Exception:
                pass
        n_emitted += 1
    # Boundary coordination: clip every ROLE_BOUNDARY shape so
    # it doesn't overlap the actual tunnel-polygon footprint.
    # The boundary ribbon is built by line-buffering the apt.dat
    # row-130 line, which extends ~2.5 m to either side of the
    # boundary line — so when a tunnel ramp crosses the boundary,
    # the ribbon overlaps both the inside-airport (LOW-end) and
    # outside-airport (HIGH-end) parts of the ramp.  Subtracting
    # the actual ramp+walls union (buffered by 0.5 m for a small
    # visible gap) handles both cases without carving exclusion
    # discs outside the tunnel's footprint — the boundary still
    # traces the rest of the perimeter unchanged.
    if exclusion_zones:
        try:
            tunnel_union = unary_union(exclusion_zones)
        except Exception:
            tunnel_union = None
        if tunnel_union is None or tunnel_union.is_empty:
            return n_emitted
        excl_union = tunnel_union.buffer(boundary_clearance_m)
        kept_shapes: List[BuiltShape] = []
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                kept_shapes.append(s)
                continue
            # Capture the old ring + altitudes BEFORE the clip so
            # we can NN-resample altitudes for the new ring.  Per
            # user 2026-04-29 (HECA crash investigation): leaving
            # the boundary ribbon with no altitudes after a clip
            # made X-Plane crash on load — the 1066-vertex polygon
            # without elevation guidance was unrenderable.
            try:
                _old_ring = list(s.polygon.exterior.coords)
            except Exception:
                _old_ring = []
            if _old_ring and _old_ring[0] == _old_ring[-1]:
                _old_ring = _old_ring[:-1]
            _old_alts = (list(s.node_altitudes)
                          if s.node_altitudes else None)
            try:
                new_poly = s.polygon.difference(excl_union)
            except Exception:
                kept_shapes.append(s)
                continue
            if new_poly.is_empty:
                continue
            if new_poly.geom_type == "Polygon":
                s.polygon = new_poly
                resampled = _resample_node_altitudes_nn(
                    new_poly, _old_ring, _old_alts)
                if resampled is not None:
                    s.node_altitudes = resampled
                # else: keep existing s.altitude / node_altitudes
                # (possibly mismatched but better than nothing).
                kept_shapes.append(s)
            elif new_poly.geom_type == "MultiPolygon":
                # Split the boundary shape into the resulting pieces.
                for g in new_poly.geoms:
                    if (g.geom_type != "Polygon"
                            or g.is_empty
                            or g.area < 5.0):
                        continue
                    resampled = _resample_node_altitudes_nn(
                        g, _old_ring, _old_alts)
                    kept_shapes.append(BuiltShape(
                        polygon=g,
                        role=s.role,
                        ref=s.ref,
                        altitude=(s.altitude if resampled is None
                                  else None),
                        node_altitudes=resampled))
        layout.shapes = kept_shapes
    return n_emitted


def _scenery_has_bridge_objects(
        layout: "PavementLayout",
        bridge_proximity_m: float = 30.0,
        ) -> bool:
    """Return True when the X-Plane scenery pack containing
    ``layout.apt_dat_path`` places at least one taxi-bridge OBJ
    near a bridge taxi rect.

    Detection logic (per user 2026-04-29):
      1. Find the DSF for the scenery pack via
         ``O4_DSF_Reader.find_associated_dsf``.
      2. Convert it to text (cached alongside the .dsf as
         ``.dsf.text``) using DSFTool when not already cached.
      3. Walk OBJECT_DEF lines.  Mark a def as a "bridge def" when
         its path matches ``bridge|elevated|viaduct|overpass``
         AND does NOT match ``sign|signage|trafficsign|wall|
         truss|crane`` — KPHX has 3 ``lib/g10/roadsigns/SignBridge
         *.obj`` defs that are road sign gantries, NOT taxi
         bridges, and the exclude regex filters them out.
      4. Walk OBJECT placement lines.  When a placement uses a
         "bridge def" AND its lat/lon lies within
         ``bridge_proximity_m`` of any taxi rect with
         ``is_bridge=True``, return True.

    Confirmed signal at the test set:
      KBNA  — 23 OBJECT_DEFs match ``Objects/KBNA Bridges/...``
              (KBNA_Bridge_Taxiway-L_p1..p6, KBNA_Crossing_Bridge,
              elevated_edge_twy_B, ...) — many placements within
              the bridge taxi rect → True.
      KPHX  — only ``lib/g10/roadsigns/SignBridge*`` defs which
              the exclude regex drops → False.
    """
    if not layout.apt_dat_path:
        return False
    bridge_rects = [s.polygon for s in layout.shapes
                     if getattr(s, "is_bridge", False)
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not bridge_rects:
        return False
    try:
        from . import dsf_reader as _DSFR
    except Exception:
        return False
    dsf_path = _DSFR.find_associated_dsf(
        layout.apt_dat_path,
        layout.anchor[0], layout.anchor[1])
    if dsf_path is None or not os.path.isfile(dsf_path):
        return False
    text_path = dsf_path + ".text"
    needs_convert = (
        not os.path.isfile(text_path)
        or os.path.getmtime(text_path) < os.path.getmtime(dsf_path))
    if needs_convert:
        tool = _DSFR._dsftool_path()
        if tool is None:
            return False
        try:
            import subprocess as _sp
            _sp.run(
                [tool, "--dsf2text", dsf_path, text_path],
                check=True, capture_output=True, timeout=120)
        except Exception:
            return False
    BRIDGE_RE = re.compile(
        r"(?i)bridge|elevated|viaduct|overpass")
    EXCLUDE_RE = re.compile(
        r"(?i)sign|signage|trafficsign|truss|crane")
    bridge_def_idx: set = set()
    object_def_count = 0
    placements: List[Tuple[int, float, float]] = []
    try:
        with open(text_path, "r", encoding="utf-8",
                  errors="replace") as f:
            for line in f:
                if line.startswith("OBJECT_DEF"):
                    parts = line.strip().split(maxsplit=1)
                    path = parts[1] if len(parts) > 1 else ""
                    if (BRIDGE_RE.search(path)
                            and not EXCLUDE_RE.search(path)):
                        bridge_def_idx.add(object_def_count)
                    object_def_count += 1
                elif line.startswith("OBJECT "):
                    tok = line.split()
                    if len(tok) >= 4:
                        try:
                            idx = int(tok[1])
                            lon = float(tok[2])
                            lat = float(tok[3])
                            placements.append((idx, lon, lat))
                        except ValueError:
                            continue
    except Exception:
        return False
    if not bridge_def_idx or not placements:
        return False
    # Project bridge rects to lat/lon for proximity check.  The
    # rects' polygons are in meter space anchored at the layout —
    # convert each placement to meters and test against the
    # buffered rect union.
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _to_m(lon_v: float, lat_v: float) -> Tuple[float, float]:
        return (math.radians(lon_v - lon0) * R * cos0,
                math.radians(lat_v - lat0) * R)
    try:
        bridge_buf = unary_union(bridge_rects).buffer(
            bridge_proximity_m)
    except Exception:
        return False
    for idx, lon_v, lat_v in placements:
        if idx not in bridge_def_idx:
            continue
        x, y = _to_m(lon_v, lat_v)
        if bridge_buf.contains(Point(x, y)):
            return True
    return False


def _emit_taxi_bridges(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        retaining_wall_width_m: float = 1.0,
        wall_gap_m: float = 0.5,
        boundary_clearance_m: float = 1.0,
        scenery_has_bridge_objects: bool = False,
        ) -> int:
    """For each taxi rect marked ``is_bridge=True``, optionally
    emit two flat retaining walls along its long edges at the
    rect's average elevation (the bridge deck altitude).

    Per user 2026-04-29 (KBNA vs KPHX): when the X-Plane scenery
    pack ALREADY contains a 3D taxi-bridge OBJ (detected by
    ``_scenery_has_bridge_objects``, e.g. KBNA's
    ``Objects/KBNA Bridges/KBNA_Bridge_Taxiway-L_*.obj``), the
    scenery's own bridge model is the visible structure — we
    skip the wall emission entirely so we don't double up.  When
    the scenery has NO bridge OBJ (KPHX), our emitted walls
    provide the only visible side-face structure.

    No end-cap walls — the bridge's short edges connect to
    adjacent taxis or junctions at apt_elev (the deck continues
    onto the surrounding airport surface), so an end-cap would
    visually block that join.

    Boundary coordination: when a bridge rect lies inside the
    airport boundary (typical at SPJC / KBNA / KPHX), the walls
    are also inside.  The inside-airport portion of (rect ∪
    walls) gets buffered by ``boundary_clearance_m`` (0.5 m) and
    subtracted from each ``ROLE_BOUNDARY`` shape — same pattern
    as ``_emit_tunnel_portals``.

    Returns the number of bridge rects whose walls were emitted
    (0 when the scenery already has bridge OBJs).
    """
    from .pipeline import _load_osm_airports, _load_osm_big_roads
    bridge_shapes = [s for s in layout.shapes
                     if getattr(s, "is_bridge", False)
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not bridge_shapes:
        return 0
    if scenery_has_bridge_objects:
        # The scenery's own 3D bridge OBJs are the visible
        # structure.  Emit nothing here — the deck itself
        # already exists as the taxi rect, and the surrounding
        # mesh is handled by the road-approach helper.
        return 0
    n_emitted = 0
    exclusion_zones: List[Polygon] = []
    for s in bridge_shapes:
        rc = list(s.polygon.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Per ``_rect_from_axis_extended`` corner convention:
        # corners 0,3 = one short edge, corners 1,2 = other.
        # Long edges: corners (0,1) and (2,3).
        # Compute the deck elevation for the wall: average of
        # altitude_high/altitude_low on a sloped rect, or
        # altitude on a flat rect.  Walls match the deck so the
        # join is seamless.
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            deck_elev = 0.5 * (s.altitude_high + s.altitude_low)
        elif s.altitude is not None:
            deck_elev = s.altitude
        else:
            # No elevation set yet — skip.  The post-elevation
            # pass calls _emit_taxi_bridges after rect altitudes
            # are filled in by the unified solver.
            continue
        for (a, b) in ((0, 1), (2, 3)):
            ax, ay = rc[a]
            bx, by = rc[b]
            ex = bx - ax
            ey = by - ay
            elen = math.hypot(ex, ey)
            if elen < 1.0:
                continue
            ux = ex / elen
            uy = ey / elen
            # Outward normal — flip if the test point is INSIDE
            # the rect.
            n_x = -uy
            n_y = ux
            mid_x = 0.5 * (ax + bx)
            mid_y = 0.5 * (ay + by)
            if s.polygon.contains(
                    Point(mid_x + n_x * 0.1,
                          mid_y + n_y * 0.1)):
                n_x = -n_x
                n_y = -n_y
            # Wall sits ``wall_gap_m`` outboard of the rect's
            # long edge, ``retaining_wall_width_m`` thick.
            inner = wall_gap_m
            outer = wall_gap_m + retaining_wall_width_m
            wc = [
                (ax + n_x * inner, ay + n_y * inner),
                (bx + n_x * inner, by + n_y * inner),
                (bx + n_x * outer, by + n_y * outer),
                (ax + n_x * outer, ay + n_y * outer),
            ]
            try:
                wall_poly = Polygon(wc)
                if not wall_poly.is_valid:
                    wall_poly = wall_poly.buffer(0)
                if (wall_poly.geom_type == "Polygon"
                        and not wall_poly.is_empty):
                    layout.shapes.append(BuiltShape(
                        polygon=wall_poly,
                        role=ROLE_RETAINING_WALL,
                        ref="bridge_wall",
                        altitude=round(float(deck_elev), 1)))
                    exclusion_zones.append(wall_poly)
            except Exception:
                continue
        # Track the bridge rect itself so the boundary subtraction
        # below also clears the deck area.
        exclusion_zones.append(s.polygon)
        n_emitted += 1
    # Boundary coordination: subtract the actual (rect ∪ walls)
    # footprint, buffered by 0.5 m, from each ROLE_BOUNDARY shape.
    # Same pattern as ``_emit_tunnel_portals``.
    if exclusion_zones:
        try:
            bridge_union = unary_union(exclusion_zones)
        except Exception:
            bridge_union = None
        if bridge_union is not None and not bridge_union.is_empty:
            try:
                excl = bridge_union.buffer(boundary_clearance_m)
                kept_shapes: List[BuiltShape] = []
                for s in layout.shapes:
                    if s.role != ROLE_BOUNDARY:
                        kept_shapes.append(s)
                        continue
                    try:
                        _old_ring = list(s.polygon.exterior.coords)
                    except Exception:
                        _old_ring = []
                    if (_old_ring
                            and _old_ring[0] == _old_ring[-1]):
                        _old_ring = _old_ring[:-1]
                    _old_alts = (list(s.node_altitudes)
                                  if s.node_altitudes else None)
                    try:
                        new_poly = s.polygon.difference(excl)
                    except Exception:
                        kept_shapes.append(s)
                        continue
                    if new_poly.is_empty:
                        continue
                    if new_poly.geom_type == "Polygon":
                        s.polygon = new_poly
                        resampled = _resample_node_altitudes_nn(
                            new_poly, _old_ring, _old_alts)
                        if resampled is not None:
                            s.node_altitudes = resampled
                        kept_shapes.append(s)
                    elif new_poly.geom_type == "MultiPolygon":
                        for g in new_poly.geoms:
                            if (g.geom_type != "Polygon"
                                    or g.is_empty
                                    or g.area < 5.0):
                                continue
                            resampled = _resample_node_altitudes_nn(
                                g, _old_ring, _old_alts)
                            kept_shapes.append(BuiltShape(
                                polygon=g, role=s.role,
                                ref=s.ref,
                                altitude=(s.altitude
                                          if resampled is None
                                          else None),
                                node_altitudes=resampled))
                layout.shapes = kept_shapes
            except Exception:
                pass
    return n_emitted


def _emit_underpass_road_approaches(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        clearance_depth_m: float = 8.0,
        approach_length_m: float = 80.0,
        road_width_m: float = 22.0,
        ramp_step_m: float = 20.0,
        scenery_has_bridge_objects: bool = False,
        ) -> int:
    """For each underpass case (taxi BRIDGE rect or road TUNNEL
    portal), emit a chain of sloped road-following polygons that
    transition the road surface from outside-DEM elevation down
    to ``apt_elev − clearance_depth_m`` near the underpass and
    back up to DEM after.

    Per user 2026-04-29 (KBNA / KPHX taxi-bridge case): without
    these polygons, the patch's mesh under a bridge sits at
    apt_elev (interpolated from the surrounding airport pavement),
    and OSM-tagged roads rendered at DEM elevation by Ortho4XP's
    road layer get buried.  Emitting a chain of road-following
    polygons forces the mesh under the road to step down to a
    height low enough for the road to clear the bridge underside
    (default 8 m below airport surface), then ramp back up to DEM
    away from the airport.

    Algorithm per underpass surface (bridge rect or tunnel portal
    region):

      1. Find OSM road LineStrings that cross the surface's
         footprint (or pass within 5 m of its long edges, for
         tunnel portals where the road just barely touches).
      2. For each crossing road, find the exterior segments
         immediately approaching and departing the surface.
      3. Walk each approach segment in ``ramp_step_m`` (20 m)
         increments, emitting a sloped 4-corner rect per step
         where elevation interpolates between DEM (start) and
         ``apt_elev − clearance_depth_m`` (end).
      4. Inside the underpass surface itself, emit a flat road-
         following polygon at ``apt_elev − clearance_depth_m``.

    Width of every emitted road polygon: ``road_width_m``
    (22 m by default; matches tunnel-ramp width).

    Returns the number of UNDERPASS surfaces processed.
    """
    from .pipeline import _load_osm_airports, _load_osm_big_roads
    # Collect underpass surfaces.
    bridge_shapes = [s for s in layout.shapes
                     if getattr(s, "is_bridge", False)
                     and s.polygon is not None
                     and not s.polygon.is_empty]
    if not bridge_shapes:
        return 0
    nodes_r, ways_r = _load_osm_big_roads(
        layout.anchor[0], layout.anchor[1])
    if not ways_r:
        return 0
    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH
    def _to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)
    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))
    nodes_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_r.items():
        nodes_m[nid] = _to_m(lon, lat)
    HW_TYPES = {
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "motorway_link", "trunk_link",
        "primary_link", "residential", "service",
    }
    # Collect candidate road LineStrings (not bridge or tunnel
    # tagged — those are special cases).  A road that PASSES
    # UNDER a taxi bridge typically isn't tagged bridge=yes
    # itself; only the airport surface is.  But a ROAD that
    # itself bridges over something else IS tagged bridge=yes;
    # we skip those because we don't want to emit road shapes
    # for road-on-road bridges.
    road_lines: List[LineString] = []
    for _wid, nrefs, tags in ways_r:
        if tags.get("highway") not in HW_TYPES:
            continue
        if tags.get("bridge") and tags.get("bridge") != "no":
            continue
        if (tags.get("tunnel")
                and tags.get("tunnel") != "no"):
            continue
        pts = [nodes_m[n] for n in nrefs if n in nodes_m]
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        road_lines.append(ls)
    if not road_lines:
        return 0
    n_processed = 0
    for s in bridge_shapes:
        # Bridge deck elevation.
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            deck_elev = 0.5 * (s.altitude_high + s.altitude_low)
        elif s.altitude is not None:
            deck_elev = s.altitude
        else:
            continue
        low_elev = float(deck_elev) - clearance_depth_m
        # Find road LineStrings that cross the bridge footprint.
        for road_ls in road_lines:
            try:
                inside = road_ls.intersection(s.polygon)
            except Exception:
                continue
            if inside.is_empty:
                continue
            # Pick the longest contiguous piece if multiple.
            if hasattr(inside, "geoms"):
                cand = [g for g in inside.geoms
                        if g.geom_type == "LineString"
                        and g.length > 1.0]
                if not cand:
                    continue
                inside = max(cand, key=lambda g: g.length)
            elif inside.geom_type != "LineString":
                continue
            # Inside-bridge flat polygon at low_elev — only emit
            # when the scenery has a 3D bridge OBJ (KBNA case).
            # In that case the road cuts straight through under
            # the bridge model.  Without a bridge OBJ (KPHX
            # case) the taxi rect itself acts as a flat plate
            # and the road approaches stop at the bridge edge.
            if scenery_has_bridge_objects:
                try:
                    inside_buf = inside.buffer(
                        road_width_m / 2.0,
                        cap_style=2, join_style=2)
                    if (inside_buf.geom_type == "Polygon"
                            and not inside_buf.is_empty):
                        layout.shapes.append(BuiltShape(
                            polygon=inside_buf,
                            role=ROLE_TUNNEL_RAMP,
                            ref="bridge_underpass",
                            altitude=round(low_elev, 1)))
                except Exception:
                    pass
            # Approach + departure ramp chains.  Find the parts
            # of the road OUTSIDE the bridge polygon, then walk
            # each in ramp_step_m steps emitting one sloped rect
            # per step.
            try:
                outside = road_ls.difference(s.polygon)
            except Exception:
                outside = None
            if outside is None or outside.is_empty:
                continue
            outside_pieces = (
                list(outside.geoms)
                if hasattr(outside, "geoms")
                else [outside])
            for piece in outside_pieces:
                if (piece.is_empty
                        or piece.geom_type != "LineString"):
                    continue
                # Decide which end of the piece TOUCHES the bridge.
                pcoords = list(piece.coords)
                if len(pcoords) < 2:
                    continue
                p_start = Point(pcoords[0])
                p_end = Point(pcoords[-1])
                d_start = s.polygon.distance(p_start)
                d_end = s.polygon.distance(p_end)
                if d_start <= d_end:
                    # piece runs FROM bridge edge OUT to far end —
                    # walk it as approach (low → high).
                    # Cap to approach_length_m.
                    walk_len = min(piece.length, approach_length_m)
                    bridge_end_pt = piece.interpolate(0.0)
                    far_end_pt = piece.interpolate(walk_len)
                    walk = LineString([(bridge_end_pt.x,
                                          bridge_end_pt.y),
                                        (far_end_pt.x,
                                          far_end_pt.y)])
                    # Use the actual sub-LineString shape for
                    # better following; but for simplicity use the
                    # straight sub-segment for ramp emission.
                    walk = LineString(pcoords[:])
                    # Trim walk to first ``walk_len`` metres.
                else:
                    walk_len = min(piece.length, approach_length_m)
                    walk = LineString(list(reversed(pcoords)))
                # Step through walk_len in ramp_step_m steps.
                # Each step is a sloped rect of (low_elev → DEM).
                u_prev = 0.0
                while u_prev < walk_len - 1.0:
                    u_next = min(walk_len, u_prev + ramp_step_m)
                    p0 = walk.interpolate(u_prev)
                    p1 = walk.interpolate(u_next)
                    seg_len = math.hypot(p1.x - p0.x, p1.y - p0.y)
                    if seg_len < 1.0:
                        break
                    # Tangent + perpendicular.
                    tx = (p1.x - p0.x) / seg_len
                    ty = (p1.y - p0.y) / seg_len
                    nx = -ty
                    ny = tx
                    half_w = road_width_m / 2.0
                    # Elevation interpolation: u_prev → low_elev
                    # at the bridge, DEM at the far end.
                    frac0 = u_prev / walk_len
                    frac1 = u_next / walk_len
                    try:
                        lat0_p, lon0_p = _m_to_ll(p0.x, p0.y)
                        lat1_p, lon1_p = _m_to_ll(p1.x, p1.y)
                        dem0 = _sample_dem(
                            dem, tile_lat, tile_lon,
                            lat0_p, lon0_p)
                        dem1 = _sample_dem(
                            dem, tile_lat, tile_lon,
                            lat1_p, lon1_p)
                    except Exception:
                        dem0 = dem1 = None
                    if dem0 is None or dem1 is None:
                        break
                    e0 = (1.0 - frac0) * low_elev + frac0 * dem0
                    e1 = (1.0 - frac1) * low_elev + frac1 * dem1
                    corners = [
                        (p0.x + nx * half_w,
                         p0.y + ny * half_w),
                        (p1.x + nx * half_w,
                         p1.y + ny * half_w),
                        (p1.x - nx * half_w,
                         p1.y - ny * half_w),
                        (p0.x - nx * half_w,
                         p0.y - ny * half_w),
                    ]
                    try:
                        seg_poly = Polygon(corners)
                        if not seg_poly.is_valid:
                            seg_poly = seg_poly.buffer(0)
                        if (seg_poly.geom_type == "Polygon"
                                and not seg_poly.is_empty):
                            # corners 0,3 = "high" end vs corners
                            # 1,2 = "low" end depends on which
                            # end is closer to the bridge.  e0
                            # (start = bridge side) is lower.
                            if abs(e0 - e1) >= 0.1:
                                layout.shapes.append(BuiltShape(
                                    polygon=seg_poly,
                                    role=ROLE_TUNNEL_RAMP,
                                    ref="bridge_approach",
                                    altitude_high=round(
                                        max(e0, e1), 1),
                                    altitude_low=round(
                                        min(e0, e1), 1)))
                            else:
                                layout.shapes.append(BuiltShape(
                                    polygon=seg_poly,
                                    role=ROLE_TUNNEL_RAMP,
                                    ref="bridge_approach",
                                    altitude=round(
                                        0.5 * (e0 + e1), 1)))
                    except Exception:
                        pass
                    u_prev = u_next
        n_processed += 1
    return n_processed


def _emit_through_airport_depressed_roads(
        layout: "PavementLayout",
        dem,
        tile_lat: int,
        tile_lon: int,
        xplane_root: str,
        icao: str,
        depression_depth_m: float = 8.0,
        max_ramp_grade: float = 0.04,
        ramp_min_length_m: float = 200.0,
        arm_max_length_m: float = 500.0,
        road_width_m: float = 22.0,
        retaining_wall_width_m: float = 1.0,
        wall_gap_m: float = 0.5,
        boundary_clearance_m: float = 1.0,
        ) -> Tuple[int, set]:
    """For each public road that ENTERS the airport boundary
    AND passes under a tagged ``aeroway=*, bridge=yes`` way (or
    is connected via OSM-graph node-sharing to a road that does)
    inside the airport, emit:

      1. A flat road-following polygon along the entire inside-
         boundary stretch at ``apt_elev − depression_depth_m``.
         Subsequent OSM road rendering (Ortho4XP's own road
         layer) sits on top of this flat plate.
      2. At each boundary entry/exit point, a ramp polygon
         OUTSIDE the airport that climbs from the depressed
         level back up to the local DEM, capped at
         ``max_ramp_grade`` (default 4 %).

    Per user 2026-04-29 (KPHX Sky Harbor Blvd): when a road
    passes under multiple airport bridges, the road MUST be
    depressed for its entire inside-airport stretch — not just
    at the bridges.  Two bridges spanning the same continuous
    road imply the road can never come back up between them.
    The general solution: any road that enters the airport
    boundary and crosses any aeroway=bridge inside is treated
    this way, plus any road CONNECTED to it via shared OSM
    nodes within the boundary (so on/off ramps inside the
    airport stay coherent with the main road).

    OSM tags consulted:
      * ``aeroway=*, bridge=yes|viaduct`` — bridge LineStrings
        (the airport surface above the depression).
      * ``highway=*`` (motorway, trunk, primary, secondary,
        tertiary, residential, service + their *_link forms) —
        the road network candidates.
      * ``bridge=yes|viaduct`` on a highway — that road is
        ITSELF a bridge over something else (skip; we're
        looking for the road UNDERNEATH).
      * ``tunnel=yes|building_passage`` on a highway — INCLUDED
        as seeds (the airport-bridge case typically tags the
        under-bridge road segment as building_passage in OSM).
        ``_emit_tunnel_portals`` is told to skip every OSM way
        depressed here so we don't double-emit.

    Boundary coordination: the union of every emitted polygon
    is buffered by ``boundary_clearance_m`` and subtracted from
    each ``ROLE_BOUNDARY`` shape (same pattern as the tunnel
    and taxi-bridge emitters), with NN-resampling of per-vertex
    altitudes so the boundary ribbon retains its altitude tags
    after the clip.

    Returns ``(n_emitted, depressed_way_ids)``.  The way-id
    set is intended for ``_emit_tunnel_portals`` to skip — it
    contains every OSM highway way (raw road network) handled
    by this pass, so the explicit-tunnel emitter doesn't
    double-process the same building_passage segments.
    """
    from .pipeline import _load_osm_airports, _load_osm_big_roads
    if (layout.airport_boundary is None
            or layout.airport_boundary.is_empty):
        return (0, set())
    boundary = layout.airport_boundary
    # Build a slightly contracted boundary for inside-vs-outside
    # tests so a road point exactly ON the boundary doesn't bounce
    # between true / false on numeric jitter.
    try:
        boundary_strict = boundary.buffer(-0.5)
        if boundary_strict.is_empty:
            boundary_strict = boundary
    except Exception:
        boundary_strict = boundary

    # Load OSM airport-layer tile (for aeroway=bridge LineStrings).
    try:
        nodes_a, ways_a, _ = _load_osm_airports(
            xplane_root, icao,
            layout.anchor[0], layout.anchor[1])
    except Exception:
        return (0, set())
    if not ways_a:
        return (0, set())
    # Load OSM big_roads (for highway ways).
    nodes_r, ways_r = _load_osm_big_roads(
        layout.anchor[0], layout.anchor[1])
    if not ways_r:
        return (0, set())

    grade_safety_margin = 0.005
    plan_grade = max(max_ramp_grade - grade_safety_margin, 1e-3)
    arm_walk_max_m = max(arm_max_length_m, ramp_min_length_m,
                         depression_depth_m / plan_grade)

    lat0, lon0 = layout.anchor
    cos0 = math.cos(math.radians(lat0))
    R = R_EARTH

    def _to_m(lon: float, lat: float) -> Tuple[float, float]:
        return (math.radians(lon - lon0) * R * cos0,
                math.radians(lat - lat0) * R)

    def _m_to_ll(x: float, y: float) -> Tuple[float, float]:
        return (lat0 + math.degrees(y / R),
                lon0 + math.degrees(x / (R * cos0)))

    def _airport_elevation_at(cx: float, cy: float) -> Optional[float]:
        # Reuse pattern from _emit_tunnel_portals: prefer the
        # boundary-ribbon's per-vertex altitude near the point,
        # fall back to DEM.
        best_d = float('inf')
        best_alt: Optional[float] = None
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                continue
            if s.ref != "airport_boundary":
                continue
            if not s.node_altitudes:
                continue
            try:
                rcoords = list(s.polygon.exterior.coords)
            except Exception:
                continue
            if rcoords and rcoords[0] == rcoords[-1]:
                rcoords = rcoords[:-1]
            for k, (vx, vy) in enumerate(rcoords):
                if k >= len(s.node_altitudes):
                    break
                d = math.hypot(vx - cx, vy - cy)
                if d < best_d:
                    best_d = d
                    best_alt = s.node_altitudes[k]
        if best_alt is not None and best_d <= 400.0:
            return float(best_alt)
        try:
            lat, lon = _m_to_ll(cx, cy)
            return _sample_dem(dem, tile_lat, tile_lon, lat, lon)
        except Exception:
            return None

    # ── Bridge LineStrings (airport-layer OSM) ─────────────────
    bridge_lines: List[LineString] = []
    nodes_a_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_a.items():
        nodes_a_m[nid] = _to_m(lon, lat)
    for wid, nrefs, tags in ways_a:
        if not tags.get("aeroway"):
            continue
        if tags.get("bridge", "") not in ("yes", "viaduct"):
            continue
        pts = [nodes_a_m[n] for n in nrefs if n in nodes_a_m]
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        bridge_lines.append(ls)
    if not bridge_lines:
        return (0, set())

    # ── Highway candidates (big_roads OSM) ─────────────────────
    HW_TYPES = {
        "motorway", "trunk", "primary", "secondary",
        "tertiary", "motorway_link", "trunk_link",
        "primary_link", "secondary_link", "tertiary_link",
        "residential", "service", "unclassified",
    }
    nodes_r_m: Dict[str, Tuple[float, float]] = {}
    for nid, (lat, lon) in nodes_r.items():
        nodes_r_m[nid] = _to_m(lon, lat)
    way_data: List[Tuple[str, LineString, List[str]]] = []
    for wid, nrefs, tags in ways_r:
        if tags.get("highway") not in HW_TYPES:
            continue
        if tags.get("bridge", "") in ("yes", "viaduct"):
            # The road IS the bridge, not what's under — skip.
            continue
        # tunnel=building_passage tagging is INCLUDED.  At KPHX,
        # the under-bridge road segments use this tag; they're
        # exactly the seeds we want.
        pts = [nodes_r_m[n] for n in nrefs if n in nodes_r_m]
        if len(pts) < 2:
            continue
        try:
            ls = LineString(pts)
        except Exception:
            continue
        if ls.is_empty or ls.length < 5.0:
            continue
        way_data.append((wid, ls, list(nrefs)))
    if not way_data:
        return (0, set())

    # ── Seed: ways whose inside-boundary section crosses a bridge ──
    BRIDGE_PROXIMITY_M = 5.0
    seed_depressed: set = set()
    inside_geom_by_wid: Dict[str, "BaseGeometry"] = {}
    for wid, ls, _nrefs in way_data:
        try:
            inside = ls.intersection(boundary)
        except Exception:
            continue
        if inside.is_empty:
            continue
        inside_geom_by_wid[wid] = inside
        # Iterate inside segments and check bridge proximity.
        segs = []
        if inside.geom_type == "LineString":
            segs = [inside]
        elif inside.geom_type == "MultiLineString":
            segs = list(inside.geoms)
        for seg in segs:
            if seg.is_empty or seg.length < 1.0:
                continue
            for bls in bridge_lines:
                try:
                    if seg.distance(bls) < BRIDGE_PROXIMITY_M:
                        seed_depressed.add(wid)
                        break
                except Exception:
                    continue
            if wid in seed_depressed:
                break
    if not seed_depressed:
        return (0, set())

    # ── BFS over OSM-graph node-sharing INSIDE the boundary ────
    # On-ramps/off-ramps inside the airport are separate OSM
    # ways; if they connect to a depressed seed at any node
    # INSIDE the boundary, they must be depressed too (otherwise
    # the seed and the connecting way disagree on altitude at
    # their shared node and X-Plane renders a cliff).
    node_to_ways: Dict[str, List[str]] = {}
    way_lookup: Dict[str, Tuple[LineString, List[str]]] = {}
    for wid, ls, nrefs in way_data:
        way_lookup[wid] = (ls, nrefs)
        for n in nrefs:
            node_to_ways.setdefault(n, []).append(wid)
    depressed_set: set = set(seed_depressed)
    queue: List[str] = list(seed_depressed)
    while queue:
        wid = queue.pop()
        ls, nrefs = way_lookup[wid]
        for n in nrefs:
            n_xy = nodes_r_m.get(n)
            if n_xy is None:
                continue
            try:
                if not boundary_strict.contains(Point(n_xy)):
                    continue
            except Exception:
                continue
            for other_wid in node_to_ways.get(n, []):
                if other_wid in depressed_set:
                    continue
                # Only propagate if the other way also has an
                # inside-boundary portion (otherwise it's just a
                # surface road glancing the boundary node).
                _o_ls, _o_nrefs = way_lookup[other_wid]
                try:
                    if _o_ls.intersection(boundary).is_empty:
                        continue
                except Exception:
                    continue
                depressed_set.add(other_wid)
                queue.append(other_wid)

    # ── Emit one set of polygons per depressed way ─────────────
    n_emitted = 0
    exclusion_zones: List[Polygon] = []
    half_w = road_width_m / 2.0

    def _smooth_walk(pts: List[Tuple[float, float]],
                      min_segment_m: float = 15.0
                      ) -> List[Tuple[float, float]]:
        """Drop near-colinear / closely-spaced intermediate
        vertices so altitude rounding can't push per-segment
        grade above the design limit."""
        if len(pts) < 3:
            return list(pts)
        merged: List[Tuple[float, float]] = [pts[0]]
        for k in range(1, len(pts)):
            d = math.hypot(pts[k][0] - merged[-1][0],
                           pts[k][1] - merged[-1][1])
            if d < min_segment_m and k != len(pts) - 1:
                continue
            merged.append(pts[k])
        return merged

    for wid in sorted(depressed_set):
        ls, nrefs = way_lookup[wid]
        # 1) Inside-boundary flat plate(s).
        try:
            inside = ls.intersection(boundary)
        except Exception:
            continue
        if inside.is_empty:
            continue
        inside_segs = []
        if inside.geom_type == "LineString":
            inside_segs = [inside]
        elif inside.geom_type == "MultiLineString":
            inside_segs = list(inside.geoms)
        for seg in inside_segs:
            if seg.is_empty or seg.length < 5.0:
                continue
            ctr = seg.centroid
            apt_elev = _airport_elevation_at(ctr.x, ctr.y)
            if apt_elev is None:
                continue
            elev_low = apt_elev - depression_depth_m
            try:
                flat_poly = seg.buffer(
                    half_w, cap_style=2, join_style=2)
                if not flat_poly.is_valid:
                    flat_poly = flat_poly.buffer(0)
            except Exception:
                continue
            if (flat_poly.is_empty
                    or flat_poly.geom_type != "Polygon"):
                continue
            layout.shapes.append(BuiltShape(
                polygon=flat_poly,
                role=ROLE_TUNNEL_RAMP,
                ref="depressed_road",
                altitude=round(elev_low, 1)))
            exclusion_zones.append(flat_poly)
            n_emitted += 1

        # 2) Outside-boundary ramp(s) — one per side of the
        #    boundary the way crosses.  Take the OSM polyline
        #    OUTSIDE the boundary; for each outside piece, walk
        #    OUTWARD from the boundary edge, truncate where the
        #    grade requirement is satisfied, then emit ramp
        #    polygons with bisector vertex sharing at bends.
        try:
            outside = ls.difference(boundary)
        except Exception:
            outside = None
        if outside is None or outside.is_empty:
            continue
        outside_pieces = ([outside]
                          if outside.geom_type == "LineString"
                          else (list(outside.geoms)
                                if outside.geom_type == "MultiLineString"
                                else []))
        for piece in outside_pieces:
            if piece.is_empty or piece.length < 5.0:
                continue
            coords = list(piece.coords)
            # Order coords so coords[0] is at the boundary,
            # coords[-1] is far outside.  The boundary edge is
            # the closest of the two endpoints to the boundary
            # exterior (distance 0 vs distance > 0).
            try:
                d_start = boundary.exterior.distance(
                    Point(coords[0]))
            except Exception:
                d_start = float('inf')
            try:
                d_end = boundary.exterior.distance(
                    Point(coords[-1]))
            except Exception:
                d_end = float('inf')
            if d_end < d_start:
                coords = list(reversed(coords))
            walk = _smooth_walk(coords)
            if len(walk) < 2:
                continue
            # Truncate the walk to a length whose grade keeps
            # ≤ plan_grade given DEM at the far end.
            apt_elev = _airport_elevation_at(*walk[0])
            if apt_elev is None:
                continue
            elev_low = apt_elev - depression_depth_m
            cum = 0.0
            kept_pts: List[Tuple[float, float]] = [walk[0]]
            grade_ok_at: float = 0.0
            for i in range(1, len(walk)):
                seg_len = math.hypot(
                    walk[i][0] - walk[i - 1][0],
                    walk[i][1] - walk[i - 1][1])
                cum += seg_len
                kept_pts.append(walk[i])
                if cum > arm_walk_max_m:
                    grade_ok_at = cum
                    break
                try:
                    plat, plon = _m_to_ll(*walk[i])
                    dem_h = _sample_dem(
                        dem, tile_lat, tile_lon, plat, plon)
                except Exception:
                    dem_h = None
                if dem_h is None:
                    continue
                drop = float(dem_h) - elev_low
                req = (drop / plan_grade if drop > 0 else 0.0)
                if cum >= req and cum >= ramp_min_length_m:
                    grade_ok_at = cum
                    break
                grade_ok_at = cum
            walk = kept_pts
            if len(walk) < 2:
                continue
            far_xy = walk[-1]
            try:
                far_lat, far_lon = _m_to_ll(*far_xy)
                far_dem = _sample_dem(
                    dem, tile_lat, tile_lon, far_lat, far_lon)
            except Exception:
                far_dem = None
            if far_dem is None:
                far_dem = apt_elev
            max_drop = plan_grade * grade_ok_at
            if (far_dem - elev_low) > max_drop:
                far_dem = elev_low + max_drop
            elev_high = far_dem
            cum_dists = [0.0]
            for i in range(1, len(walk)):
                cum_dists.append(cum_dists[-1] + math.hypot(
                    walk[i][0] - walk[i - 1][0],
                    walk[i][1] - walk[i - 1][1]))
            total_walk = cum_dists[-1]
            if total_walk < 5.0:
                continue
            # Per-vertex bisector perpendicular for shared corners
            # at bends (same pattern as _emit_tunnel_portals).
            n_w = len(walk)
            verts_perp: List[Tuple[float, float]] = []
            verts_scale: List[float] = []
            for i in range(n_w):
                if i == 0:
                    s = (walk[1][0] - walk[0][0],
                         walk[1][1] - walk[0][1])
                    sl = math.hypot(*s) or 1e-6
                    verts_perp.append((-s[1] / sl, s[0] / sl))
                    verts_scale.append(1.0)
                elif i == n_w - 1:
                    s = (walk[i][0] - walk[i - 1][0],
                         walk[i][1] - walk[i - 1][1])
                    sl = math.hypot(*s) or 1e-6
                    verts_perp.append((-s[1] / sl, s[0] / sl))
                    verts_scale.append(1.0)
                else:
                    s1 = (walk[i][0] - walk[i - 1][0],
                          walk[i][1] - walk[i - 1][1])
                    s2 = (walk[i + 1][0] - walk[i][0],
                          walk[i + 1][1] - walk[i][1])
                    l1 = math.hypot(*s1) or 1e-6
                    l2 = math.hypot(*s2) or 1e-6
                    u1 = (s1[0] / l1, s1[1] / l1)
                    u2 = (s2[0] / l2, s2[1] / l2)
                    avg = ((u1[0] + u2[0]) / 2.0,
                           (u1[1] + u2[1]) / 2.0)
                    al = math.hypot(*avg)
                    if al < 1e-6:
                        verts_perp.append((-u1[1], u1[0]))
                        verts_scale.append(1.0)
                        continue
                    tangent = (avg[0] / al, avg[1] / al)
                    perp = (-tangent[1], tangent[0])
                    dot = u1[0] * u2[0] + u1[1] * u2[1]
                    cos_half = max(0.1, math.sqrt(
                        max(0.0, (1.0 + dot) / 2.0)))
                    verts_perp.append(perp)
                    verts_scale.append(1.0 / cos_half)

            def _vertex_offset(idx: int, off: float
                               ) -> Tuple[float, float]:
                px, py = walk[idx]
                nx, ny = verts_perp[idx]
                scaled = off * verts_scale[idx]
                return (px + nx * scaled, py + ny * scaled)

            for i in range(n_w - 1):
                d_a = cum_dists[i]
                d_b = cum_dists[i + 1]
                seg_len = d_b - d_a
                if seg_len < 0.5:
                    continue
                frac_a = d_a / total_walk if total_walk > 0 else 0.0
                frac_b = d_b / total_walk if total_walk > 0 else 0.0
                e_a = (1 - frac_a) * elev_low + frac_a * elev_high
                e_b = (1 - frac_b) * elev_low + frac_b * elev_high
                ra = _vertex_offset(i, +half_w)
                rb = _vertex_offset(i + 1, +half_w)
                rc = _vertex_offset(i + 1, -half_w)
                rd = _vertex_offset(i, -half_w)
                # Corners [0, 3] = HIGH end, [1, 2] = LOW end.
                if e_b >= e_a:
                    ramp_corners = [rb, ra, rd, rc]
                    eh, el = e_b, e_a
                else:
                    ramp_corners = [ra, rb, rc, rd]
                    eh, el = e_a, e_b
                try:
                    rp = Polygon(ramp_corners)
                    if not rp.is_valid:
                        rp = rp.buffer(0)
                    if (rp.geom_type != "Polygon"
                            or rp.is_empty
                            or rp.area < 0.5):
                        continue
                except Exception:
                    continue
                if abs(eh - el) >= 0.1:
                    layout.shapes.append(BuiltShape(
                        polygon=rp,
                        role=ROLE_TUNNEL_RAMP,
                        ref="depressed_approach",
                        altitude_high=round(eh, 1),
                        altitude_low=round(el, 1)))
                else:
                    layout.shapes.append(BuiltShape(
                        polygon=rp,
                        role=ROLE_TUNNEL_RAMP,
                        ref="depressed_approach",
                        altitude=round(0.5 * (eh + el), 1)))
                exclusion_zones.append(rp)

    # ── Boundary coordination ─────────────────────────────────
    if exclusion_zones:
        try:
            depressed_union = unary_union(exclusion_zones)
        except Exception:
            depressed_union = None
        if (depressed_union is None
                or depressed_union.is_empty):
            return (n_emitted, depressed_set)
        excl_union = depressed_union.buffer(boundary_clearance_m)
        kept_shapes: List[BuiltShape] = []
        for s in layout.shapes:
            if s.role != ROLE_BOUNDARY:
                kept_shapes.append(s)
                continue
            try:
                _old_ring = list(s.polygon.exterior.coords)
            except Exception:
                _old_ring = []
            if _old_ring and _old_ring[0] == _old_ring[-1]:
                _old_ring = _old_ring[:-1]
            _old_alts = (list(s.node_altitudes)
                          if s.node_altitudes else None)
            try:
                new_poly = s.polygon.difference(excl_union)
            except Exception:
                kept_shapes.append(s)
                continue
            if new_poly.is_empty:
                continue
            if new_poly.geom_type == "Polygon":
                s.polygon = new_poly
                resampled = _resample_node_altitudes_nn(
                    new_poly, _old_ring, _old_alts)
                if resampled is not None:
                    s.node_altitudes = resampled
                kept_shapes.append(s)
            elif new_poly.geom_type == "MultiPolygon":
                for g in new_poly.geoms:
                    if (g.geom_type != "Polygon"
                            or g.is_empty
                            or g.area < 5.0):
                        continue
                    resampled = _resample_node_altitudes_nn(
                        g, _old_ring, _old_alts)
                    kept_shapes.append(BuiltShape(
                        polygon=g,
                        role=s.role,
                        ref=s.ref,
                        altitude=(s.altitude
                                  if resampled is None
                                  else None),
                        node_altitudes=resampled))
        layout.shapes = kept_shapes
    return (n_emitted, depressed_set)
