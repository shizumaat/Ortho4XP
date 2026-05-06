"""Shared-vertex enforcement and snap helpers.

Geometric-invariant utilities used after every shape-build pass to
enforce the "Shared vertices exact between adjacent shapes" rule
(memory: feedback_shape_rules.md, Coverage invariants).  These
functions run as a finalisation layer that walks every emitted
PavementLayout shape and snaps near-coincident vertices to the
same coordinate, drops spike vertices that would self-intersect
under .11f OSM precision, and validates that adjacent shapes share
exact corners.

Public API (all kept with leading-underscore names for backward
compatibility with internal callers in O4_Airport_Pavement_Builder):
    _snap_polygon_vertices_to_rect_corners
    _push_junction_vertices_off_taxi_rect_edges
    _drop_spike_vertices
    _enforce_shared_vertices
    _validate_shared_vertex_invariant
"""
from __future__ import annotations

import math
import sys
from typing import Dict, List, Optional, Set, Tuple

from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree

from ..layout import (
    BuiltShape,
    PavementLayout,
    ROLE_BOUNDARY,
    ROLE_CROSS_CONNECTOR,
    ROLE_GROUNDSIDE_PAVEMENT,
    ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
    SHARED_VERTEX_TOL_M,
)
from ..config import (
    JUNCTION_CLUSTER_DIST_M,
    SLIVER_ANGLE_THRESHOLD_DEG,
)



# Max perpendicular distance from a non-adjacent ring edge below
# which a vertex is treated as a "spike" and dropped.  The .11f
# OSM-emit truncation can collapse a sub-mm spike onto its neighbour
# edge and produce a self-intersecting polygon, which crashes
# X-Plane's mesh builder.  5 mm is well below any visible precision
# but well above float-precision noise.
SPIKE_VERTEX_TOL_M = 0.005


__all__ = [
    "SPIKE_VERTEX_TOL_M",
    "_drop_spike_vertices",
    "_enforce_shared_vertices",
    "_push_junction_vertices_off_taxi_rect_edges",
    "_snap_polygon_vertices_to_rect_corners",
    "_validate_shared_vertex_invariant",
]




def _snap_polygon_vertices_to_rect_corners(
        poly: "Polygon",
        sloping_rect_polys: "List[Polygon]",
        snap_tol_m: float = 5.0,
        ) -> "Polygon":
    """Snap every vertex of ``poly`` that lies within ``snap_tol_m``
    of any sloping-rect corner to that corner.

    Per user 2026-04-28 invariant: a sloping rect (runway, primary/
    secondary parallel, stub, cross-connector) can only share a
    *corner* with an adjacent junction polygon — never a point on
    one of its four edges.  Edge-interior coincidence breaks the
    rect's straight-line slope by injecting an extra elevation
    constraint at a non-corner location.

    The runway-crossing-junction emit (``_resolve_runway_crossings``)
    builds a junction polygon from the union of crossing runway
    segments.  Shapely's ``unary_union`` produces vertices at every
    boundary intersection point; some of those points land 2-5 m
    *along* a surviving rect's long edge instead of *at* the rect's
    corner because the dropped (in-crossing) and surviving (out-of-
    crossing) runway segments don't perfectly tile.  Snapping near-
    corner vertices fixes the immediate violation without distorting
    the union polygon's overall footprint.

    Consecutive duplicate vertices produced by the snap are deduped.
    Returns the input polygon unchanged if snapping would leave
    fewer than 3 distinct vertices or produce an invalid polygon.
    """
    try:
        coords = list(poly.exterior.coords)
    except Exception:
        return poly
    if coords and coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) < 3:
        return poly

    corners: List[Tuple[float, float]] = []
    for r in sloping_rect_polys:
        if r is None or r.is_empty:
            continue
        try:
            rc = list(r.exterior.coords)
        except Exception:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        corners.extend((float(x), float(y)) for x, y in rc)
    if not corners:
        return poly

    snap_tol2 = snap_tol_m * snap_tol_m
    snapped: List[Tuple[float, float]] = []
    for vx, vy in coords:
        best_corner: Optional[Tuple[float, float]] = None
        best_d2 = snap_tol2
        for cx, cy in corners:
            d2 = (vx - cx) ** 2 + (vy - cy) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_corner = (cx, cy)
        if best_corner is not None:
            snapped.append(best_corner)
        else:
            snapped.append((float(vx), float(vy)))

    # Dedupe consecutive duplicates (within 1 cm).
    deduped: List[Tuple[float, float]] = []
    for v in snapped:
        if (not deduped
                or (v[0] - deduped[-1][0]) ** 2
                + (v[1] - deduped[-1][1]) ** 2 > 1e-4):
            deduped.append(v)
    if (len(deduped) > 1
            and (deduped[0][0] - deduped[-1][0]) ** 2
            + (deduped[0][1] - deduped[-1][1]) ** 2 < 1e-4):
        deduped.pop()
    if len(deduped) < 3:
        return poly

    try:
        new_poly = Polygon(deduped + [deduped[0]])
        if not new_poly.is_valid:
            new_poly = new_poly.buffer(0)
        if (new_poly.is_empty
                or new_poly.geom_type != "Polygon"):
            return poly
        return new_poly
    except Exception:
        return poly


def _push_junction_vertices_off_taxi_rect_edges(
        layout: "PavementLayout",
        edge_tol_m: float = 0.5,
        corner_tol_m: float = 2.0,
        edge_gap_m: float = 1.0,
        ) -> int:
    """Per the user 2026-04-24 invariant: a junction polygon may
    share a vertex with a taxi rect ONLY at one of the rect's 4
    corners.  A junction vertex landing on the INTERIOR of a rect
    edge would split that edge at render time and break the rect's
    altitude_high/altitude_low slope convention.

    Per user 2026-04-29: junctions and taxiway/stub/runway rects
    should be handled the SAME WAY as runway-runway crossings —
    the junction should connect to either the LOW or HIGH short
    edge of a rect (i.e. its corners), never to the rect's edge
    interior.

    Two-stage policy applied to every junction ring vertex:

      Stage 1 — collapse redundant edge-interior vertices.
        If a vertex lies on the interior of a rect edge AND its
        ring-adjacent neighbours are both at corners of the SAME
        rect (one each, on the two ends of THAT edge), the vertex
        is redundant.  The junction's ring already connects the
        two corners; the intermediate vertex just splits a single
        rect edge into two pieces.  Drop the vertex — the junction
        edge then runs corner-to-corner along the rect's short
        edge (HIGH-side or LOW-side) cleanly.
      Stage 2 — corner snap / push for survivors.
        For each remaining vertex:
          * Within ``corner_tol_m`` of a rect corner ⇒ snap to
            that corner.
          * Within ``edge_tol_m`` of a rect edge interior ⇒ push
            ``edge_gap_m`` perpendicular outside the rect.
          * Otherwise ⇒ leave alone.

    Geometric only: doesn't touch elevations.  Runway corners are
    treated identically to taxi rect corners so junction vertices
    snap there too — the user's "same logic for all runway
    intersections" requirement.

    Returns the number of junction polygons modified.
    """
    rect_roles = {
        ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
        ROLE_STUB, ROLE_CROSS_CONNECTOR, ROLE_RUNWAY}
    rects: List[Tuple[Polygon, List[Tuple[float, float]]]] = []
    for s in layout.shapes:
        if s.role not in rect_roles:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
        except Exception:
            continue
        if len(coords) != 4:
            continue
        rects.append((s.polygon, coords))
    if not rects:
        return 0

    corner_tol2 = corner_tol_m * corner_tol_m

    def _on_edge_between_corners(
            x: float, y: float,
            corners: List[Tuple[float, float]],
            ) -> Optional[int]:
        """Return the edge index (0-3) the point lies on (within
        ``edge_tol_m`` of an edge interior, t ∈ (ε, 1-ε)), or
        None if the point isn't on any rect edge interior."""
        for i in range(4):
            ax, ay = corners[i]
            bx, by = corners[(i + 1) % 4]
            dx = bx - ax
            dy = by - ay
            seg_len_sq = dx * dx + dy * dy
            if seg_len_sq <= 0.01:
                continue
            t = ((x - ax) * dx + (y - ay) * dy) / seg_len_sq
            if t <= 0.001 or t >= 0.999:
                continue
            cx_proj = ax + t * dx
            cy_proj = ay + t * dy
            d_sq = ((x - cx_proj) ** 2
                    + (y - cy_proj) ** 2)
            if d_sq <= edge_tol_m * edge_tol_m:
                return i
        return None

    def _at_corner_index(
            x: float, y: float,
            corners: List[Tuple[float, float]],
            ) -> Optional[int]:
        """Return the corner index (0-3) the point lies at
        (within ``corner_tol_m``), or None."""
        for ci, (cx, cy) in enumerate(corners):
            if (x - cx) ** 2 + (y - cy) ** 2 <= corner_tol2:
                return ci
        return None

    def _push_off(
            x: float, y: float,
            rect_poly: Polygon,
            corners: List[Tuple[float, float]],
            edge_idx: int,
            ) -> Tuple[float, float]:
        """Push the point ``edge_gap_m`` perpendicular to the
        rect edge ``edge_idx``, toward the OUTSIDE of the rect."""
        ax, ay = corners[edge_idx]
        bx, by = corners[(edge_idx + 1) % 4]
        dx = bx - ax
        dy = by - ay
        seg_len = math.sqrt(dx * dx + dy * dy)
        if seg_len <= 0.01:
            return (x, y)
        t = ((x - ax) * dx + (y - ay) * dy) / (seg_len * seg_len)
        cx_proj = ax + t * dx
        cy_proj = ay + t * dy
        perp_x = -dy / seg_len
        perp_y = dx / seg_len
        test_x = cx_proj + perp_x * 0.1
        test_y = cy_proj + perp_y * 0.1
        if rect_poly.contains(Point(test_x, test_y)):
            perp_x = -perp_x
            perp_y = -perp_y
        return (cx_proj + perp_x * edge_gap_m,
                cy_proj + perp_y * edge_gap_m)

    n_modified = 0
    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        try:
            ring = list(shape.polygon.exterior.coords)
        except Exception:
            continue
        # Drop closing repeat for ring traversal.
        if ring and ring[0] == ring[-1]:
            ring_open = ring[:-1]
        else:
            ring_open = ring
        n_v = len(ring_open)
        if n_v < 3:
            continue

        # Stage 1: collapse redundant edge-interior vertices.
        # A vertex v is redundant if its ring-prev and ring-next
        # neighbours are at the two corners of the same rect edge
        # AND v itself lies on that edge's interior.
        keep_mask = [True] * n_v
        for i in range(n_v):
            vx, vy = ring_open[i]
            px, py = ring_open[(i - 1) % n_v]
            nx, ny = ring_open[(i + 1) % n_v]
            for rect_poly, corners in rects:
                p_corner = _at_corner_index(px, py, corners)
                n_corner = _at_corner_index(nx, ny, corners)
                if p_corner is None or n_corner is None:
                    continue
                # The two neighbour corners must be adjacent
                # corners (i.e. share an edge).  Adjacent corner
                # pairs: (0,1), (1,2), (2,3), (3,0).
                diff = abs(p_corner - n_corner)
                if diff != 1 and diff != 3:
                    continue
                # The edge between them is index = min(...) if
                # adjacent, but for the wrap (0,3 / 3,0) it's
                # edge 3.  Just identify by the corner pair.
                edge_idx = (
                    min(p_corner, n_corner)
                    if diff == 1 else 3)
                v_edge = _on_edge_between_corners(vx, vy, corners)
                if v_edge != edge_idx:
                    continue
                # Vertex v sits on the edge between two corners
                # that are already in the ring as neighbours.
                # Drop it — the junction's ring will go
                # corner-to-corner along the rect edge.
                keep_mask[i] = False
                break
        # Survivor index list lets us drop the matching entries from
        # ``shape.node_altitudes`` after the rebuild — the dropped
        # vertices' altitudes are no longer needed but the kept
        # vertices' altitudes ARE still valid (the solver's elevation
        # field depends on the polygon being valid AND on per-vertex
        # values; nuking them all forces re-derivation from scratch
        # and reintroduces grade violations the solver already fixed).
        survivor_idx = [i for i in range(n_v) if keep_mask[i]]
        ring_after_collapse = [ring_open[i] for i in survivor_idx]
        n_collapsed = n_v - len(ring_after_collapse)

        # Stage 2: corner-snap / edge-push the survivors.
        new_ring: List[Tuple[float, float]] = []
        n_snapped = 0
        n_pushed = 0
        for vx, vy in ring_after_collapse:
            target = (vx, vy)
            for rect_poly, corners in rects:
                ci = _at_corner_index(vx, vy, corners)
                if ci is not None:
                    target = corners[ci]
                    if target != (vx, vy):
                        n_snapped += 1
                    break
                ei = _on_edge_between_corners(vx, vy, corners)
                if ei is not None:
                    target = _push_off(
                        vx, vy, rect_poly, corners, ei)
                    if target != (vx, vy):
                        n_pushed += 1
                    break
            new_ring.append(target)

        if n_collapsed == 0 and n_snapped == 0 and n_pushed == 0:
            continue

        # Re-close the ring and rebuild the polygon.
        if new_ring and new_ring[0] != new_ring[-1]:
            new_ring_closed = new_ring + [new_ring[0]]
        else:
            new_ring_closed = new_ring
        try:
            new_poly = Polygon(new_ring_closed,
                                list(shape.polygon.interiors))
            buffer_repaired = False
            if not new_poly.is_valid:
                # The original ring may already be self-intersecting
                # — buffer(0) can return a MultiPolygon with one
                # main piece + tiny artifacts.  Take the largest
                # Polygon piece so the rect-edge fix still applies.
                fixed = new_poly.buffer(0)
                if fixed.geom_type == "Polygon":
                    new_poly = fixed
                    buffer_repaired = True
                elif (fixed.geom_type == "MultiPolygon"
                        and not fixed.is_empty):
                    new_poly = max(
                        fixed.geoms, key=lambda g: g.area)
                    buffer_repaired = True
                else:
                    new_poly = None
            if (new_poly is not None
                    and new_poly.geom_type == "Polygon"
                    and not new_poly.is_empty):
                # Capture the OLD ring + per-vertex altitudes before
                # we overwrite the polygon — needed for the nearest-
                # neighbour fallback below when the new ring's
                # vertex count differs.
                _old_alts = (list(shape.node_altitudes)
                              if shape.node_altitudes else None)
                _old_open = list(ring_open)  # captured pre-rebuild
                shape.polygon = new_poly
                n_new = len(list(new_poly.exterior.coords)) - 1
                # Preserve per-vertex altitudes where we can.  When
                # Stage 1 collapsed K vertices but Stage 2 only
                # snapped/pushed in place, the survivor index list
                # maps the new ring to the original altitudes 1:1.
                # When buffer(0) restructured the ring or vertex
                # counts otherwise don't match, fall back to a
                # NEAREST-NEIGHBOUR resampling against the old
                # ring so the polygon retains its elevation field.
                # Per user 2026-04-29 (CYXY apron regression):
                # just setting ``node_altitudes = None`` left
                # CYXY's 174 k m² apron junction with no altitude
                # at all — every X-Plane interpolation neighbour
                # was at a different elevation, producing the
                # "terrain all over the place" the user saw.
                # Nearest-neighbour resample of the new ring against
                # the captured OLD ring vertices.  Per user
                # 2026-04-29 (CYXY apron regression): the prior
                # behaviour of setting ``node_altitudes = None``
                # whenever the vertex count changed left the giant
                # 174 k m² apron junction with no altitude at all,
                # producing the "terrain all over the place" the
                # user saw.  NN resampling preserves the elevation
                # field through Stage-1 collapse, Stage-2 push, AND
                # buffer(0) MultiPolygon repair.
                if not _old_alts or not _old_open:
                    n_modified += 1
                    continue
                src_alts_open = (
                    _old_alts[:-1]
                    if (len(_old_alts) == len(_old_open) + 1
                        and _old_alts[0] == _old_alts[-1])
                    else _old_alts[:len(_old_open)])
                if not src_alts_open:
                    n_modified += 1
                    continue
                new_open = list(new_poly.exterior.coords)
                if new_open and new_open[0] == new_open[-1]:
                    new_open = new_open[:-1]
                new_alts: List[float] = []
                for nx, ny in new_open:
                    best_d2 = float("inf")
                    best_a = src_alts_open[0]
                    for k, (sx, sy) in enumerate(_old_open):
                        if k >= len(src_alts_open):
                            break
                        d2 = (nx - sx) ** 2 + (ny - sy) ** 2
                        if d2 < best_d2:
                            best_d2 = d2
                            best_a = src_alts_open[k]
                    new_alts.append(round(float(best_a), 1))
                if new_alts:
                    shape.node_altitudes = (
                        new_alts + [new_alts[0]])
                n_modified += 1
        except Exception:
            pass
    return n_modified


def _drop_spike_vertices(
    ring: List[Tuple[float, float]],
) -> List[Tuple[float, float]]:
    """Drop any ring vertex that lies within ``SPIKE_VERTEX_TOL_M``
    of a NON-adjacent edge of the same ring.

    Such a vertex represents a degenerate "stick-out-and-return"
    in the polygon boundary — the ring went away from a straight
    section and returned right onto it, leaving a near-zero-area
    lobe that's a self-touch / self-intersection in any reasonable
    coordinate precision.

    Iterates to a fixed point (dropping one spike can expose
    another).  Capped at 8 passes against pathological inputs.
    """
    if len(ring) < 4:
        return ring
    for _ in range(8):
        n = len(ring)
        if n < 4:
            break
        keep = [True] * n
        for i in range(n):
            vx, vy = ring[i]
            for j in range(n):
                if abs(i - j) <= 1 or (i == 0 and j == n - 1) or (j == 0 and i == n - 1):
                    continue
                ax, ay = ring[j]
                bx, by = ring[(j + 1) % n]
                dx = bx - ax
                dy = by - ay
                seg2 = dx * dx + dy * dy
                if seg2 < 1e-6:
                    continue
                t = ((vx - ax) * dx + (vy - ay) * dy) / seg2
                if t < 0.0 or t > 1.0:
                    continue
                cx = ax + t * dx
                cy = ay + t * dy
                d2 = (vx - cx) * (vx - cx) + (vy - cy) * (vy - cy)
                if d2 < SPIKE_VERTEX_TOL_M * SPIKE_VERTEX_TOL_M:
                    keep[i] = False
                    break
        new_ring = [r for r, k in zip(ring, keep) if k]
        if len(new_ring) == n:
            break
        ring = new_ring
    return ring




def _enforce_shared_vertices(layout: "PavementLayout",
                             tol: float = 1.5) -> None:
    """Collapse all emitted-shape vertices that lie within ``tol``
    of each other to a single canonical point (the cluster mean),
    then rewrite each shape's polygon with those canonical vertices.

    Implements rule 16 (exact shared vertices between adjacent
    shapes).  Must run AFTER all shapes are emitted.
    """
    # Gather every vertex with a (shape_idx, is_interior, ring_idx,
    # vert_idx) handle so we can rewrite them in place.
    handles: List[Tuple[int, int, int, int, Tuple[float, float]]] = []
    for si, shape in enumerate(layout.shapes):
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        ext = list(poly.exterior.coords)
        if ext and ext[0] == ext[-1]:
            ext = ext[:-1]
        for vi, v in enumerate(ext):
            handles.append((si, 0, 0, vi, (v[0], v[1])))
        for ri, ring in enumerate(poly.interiors):
            rc = list(ring.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            for vi, v in enumerate(rc):
                handles.append((si, 1, ri, vi, (v[0], v[1])))
    if not handles:
        return

    from collections import defaultdict

    # Union-find with O(n²) pair scan.  n is typically 200-2000
    # across both airports, well within millisecond range, and the
    # simpler code eliminates any spatial-index off-by-one bugs.
    n = len(handles)
    parent = list(range(n))

    def _find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def _union(a, b):
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[ra] = rb

    coords_only = [h[4] for h in handles]
    for i in range(n):
        ix, iy = coords_only[i]
        for j in range(i + 1, n):
            jx, jy = coords_only[j]
            dx = ix - jx
            dy = iy - jy
            if dx > tol or dx < -tol or dy > tol or dy < -tol:
                continue
            if math.hypot(dx, dy) <= tol:
                _union(i, j)

    # Compute cluster centroids (mean of member coords).
    cluster_members: Dict[int, List[int]] = defaultdict(list)
    for i in range(len(handles)):
        cluster_members[_find(i)].append(i)
    canonical: Dict[int, Tuple[float, float]] = {}
    for root, members in cluster_members.items():
        sx = sum(handles[m][4][0] for m in members) / len(members)
        sy = sum(handles[m][4][1] for m in members) / len(members)
        canonical[root] = (sx, sy)

    # Rewrite each shape's rings with the canonical coords.
    new_coords_by_shape: Dict[int, Dict[Tuple[int, int, int],
                                        Tuple[float, float]]] = defaultdict(dict)
    for i, h in enumerate(handles):
        si, is_int, ri, vi, _orig = h
        new_coords_by_shape[si][(is_int, ri, vi)] = canonical[_find(i)]

    for si, shape in enumerate(layout.shapes):
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        if si not in new_coords_by_shape:
            continue
        # Rebuild exterior.
        ext = list(poly.exterior.coords)
        if ext and ext[0] == ext[-1]:
            ext = ext[:-1]
        new_ext = [new_coords_by_shape[si].get((0, 0, vi), ext[vi])
                   for vi in range(len(ext))]
        # Drop consecutive duplicates that arose from clustering.
        dedup_ext: List[Tuple[float, float]] = []
        for c in new_ext:
            if not dedup_ext or math.hypot(
                    c[0] - dedup_ext[-1][0],
                    c[1] - dedup_ext[-1][1]) > 0.05:
                dedup_ext.append(c)
        if (len(dedup_ext) >= 2
                and math.hypot(dedup_ext[0][0] - dedup_ext[-1][0],
                               dedup_ext[0][1] - dedup_ext[-1][1]) < 0.05):
            dedup_ext = dedup_ext[:-1]
        if len(dedup_ext) < 3:
            # Shape collapsed to a degenerate sliver after cluster
            # rewrite (e.g. a thin grid-decomposition sliver whose
            # vertices got pulled together).  Empty the polygon so
            # the un-clustered original isn't left behind to
            # violate the shared-vertex invariant.  Empty polygons
            # are skipped by the validator and ``to_osm``.
            shape.polygon = Polygon()
            continue
        # Rebuild interiors.
        new_interiors: List[List[Tuple[float, float]]] = []
        for ri, ring in enumerate(poly.interiors):
            rc = list(ring.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            new_ring = [new_coords_by_shape[si].get((1, ri, vi), rc[vi])
                        for vi in range(len(rc))]
            dedup_ring: List[Tuple[float, float]] = []
            for c in new_ring:
                if not dedup_ring or math.hypot(
                        c[0] - dedup_ring[-1][0],
                        c[1] - dedup_ring[-1][1]) > 0.05:
                    dedup_ring.append(c)
            if len(dedup_ring) >= 3:
                new_interiors.append(dedup_ring)
        try:
            new_poly = Polygon(dedup_ext, new_interiors)
            # The cluster-rewrite step can create a self-touching
            # ring when two NON-adjacent ring vertices end up at
            # the same canonical cluster point — the consecutive-
            # dedup above only catches adjacent duplicates.
            # ``_drop_overlap_against_fixed_shapes`` runs buffer(0)
            # on each shape and DROPS the whole shape if the result
            # is a MultiPolygon.  Recover here by buffer(0) → keep
            # largest piece, before the clip sees it (per user
            # 2026-05-04 south-of-terminal1 130K-m² apron drop).
            if (new_poly.geom_type == "Polygon"
                    and not new_poly.is_empty
                    and not new_poly.is_valid):
                try:
                    fixed = new_poly.buffer(0)
                    if not fixed.is_empty:
                        if fixed.geom_type == "MultiPolygon":
                            fixed = max(fixed.geoms,
                                        key=lambda g: g.area)
                        if (fixed.geom_type == "Polygon"
                                and fixed.is_valid
                                and not fixed.is_empty):
                            new_poly = fixed
                except Exception:
                    pass
            if (new_poly.geom_type == "Polygon"
                    and not new_poly.is_empty):
                shape.polygon = new_poly
        except Exception:
            pass


def _validate_shared_vertex_invariant(layout: "PavementLayout",
                                      tol: float = 1.5) -> None:
    """Assert that every pair of shape vertices is EITHER exactly
    equal (< 0.01 m after clustering) OR > ``tol`` apart.  A
    "close but not equal" pair violates rule 16 and signals a
    clustering bug.  Raises RuntimeError on violation.
    """
    verts: List[Tuple[int, Tuple[float, float]]] = []
    for si, shape in enumerate(layout.shapes):
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        ext = list(poly.exterior.coords)
        if ext and ext[0] == ext[-1]:
            ext = ext[:-1]
        for v in ext:
            verts.append((si, (v[0], v[1])))
    if len(verts) < 2:
        return
    # Grid-bucket check: every pair within tol must be within 0.01.
    from collections import defaultdict
    cell = tol * 2.0
    buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, (_, (x, y)) in enumerate(verts):
        buckets[(int(x // cell), int(y // cell))].append(i)
    for (gx, gy), idxs in buckets.items():
        neigh: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                neigh.extend(buckets.get((gx + dx, gy + dy), []))
        for a in idxs:
            ax, ay = verts[a][1]
            for b in neigh:
                if b <= a:
                    continue
                bx, by = verts[b][1]
                d = math.hypot(ax - bx, ay - by)
                if 0.01 < d <= tol:
                    raise RuntimeError(
                        f"Shared-vertex invariant violated: shape {verts[a][0]}"
                        f" @ ({ax:.3f},{ay:.3f}) and shape {verts[b][0]}"
                        f" @ ({bx:.3f},{by:.3f}) are {d:.3f} m apart"
                        f" (within tol={tol} m but not exactly equal).")


# ──────────────────────────────────────────────────────────────────
# Centerline-based taxi rect builder
# ──────────────────────────────────────────────────────────────────

                                # — between 80 (too coarse, merged
                                # distinct crossings) and 25 (too
                                # fine, created spurious crossings at
                                # every sub-way endpoint)
JUNCTION_RADIUS_SCALE = 1.5     # disc radius = local_half_width × this
