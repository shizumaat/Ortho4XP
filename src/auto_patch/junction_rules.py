"""Junction-refinement rule passes (user 2026-05-01).

Plan: ``/Users/noah/.claude/plans/kind-meandering-sifakis.md``.

Four rules apply as post-emission passes inside
``junction_emit.emit_junctions_and_finalize``:

* **Rule 1** — junction-runway 1:1 vertex sharing.
* **Rule 2** — snap junction vertices within ``LONG_EDGE_SNAP_M`` of
  a sloping rect's long edge to the rect's nearest short-end corner.
* **Rule 3** — for every junction, every edge that isn't a
  pavement-boundary arc and isn't a shared anchor edge must run
  parallel or perpendicular to the longest runway axis.
* **Rule 4** — split junctions at narrow necks; absorb the smaller
  piece into a neighbour or emit standalone.

Each pass mutates ``layout.shapes`` in place.

Implementation phase order (per plan): Rule 2 → Rule 1 → Rule 4 →
Rule 3.  Stubs raise ``NotImplementedError`` until landed.
"""
from __future__ import annotations

import math
from typing import List, Optional, Sequence, Tuple

from shapely.geometry import LineString, Point, Polygon

from .config import (
    AXIS_ALIGN_TOL_DEG,
    LONG_EDGE_SNAP_M,
    NECK_ABSOLUTE_M,
    NECK_ABSORB_FRAC,
    NECK_RELATIVE,
    RUNWAY_BOUNDARY_TOL_M,
)
from .layout import (
    BuiltShape,
    PavementLayout,
    ROLE_JUNCTION,
    ROLE_RUNWAY,
    SHARED_VERTEX_TOL_M,
)


__all__ = [
    "apply_junction_rules",
    "longest_runway_axis_deg",
]


SLOPING_RECT_ROLES = (
    "primary_parallel",
    "secondary_parallel",
    "stub",
    "cross_connector",
)


def apply_junction_rules(layout: PavementLayout) -> None:
    """Run the rule passes in implementation-phase order.  Mutates
    ``layout.shapes`` in place.  Skips passes whose helpers haven't
    landed yet.
    """
    runway_axis_deg = longest_runway_axis_deg(layout)

    # Phase 1 (landed): Rule 2 — long-edge corner snap.
    # Use the FINAL rect polygons from layout.shapes — earlier passes
    # (overlap clip, shared-vertex centroid collapse) may have shifted
    # the original ``taxi_rects`` polygons; reading from layout
    # guarantees we snap against the geometry the test sees.
    _snap_to_long_edge_corners(layout)

    # Phase 2 (landed): Rule 1 — junction-runway 1:1 sharing.
    _enforce_runway_1to1_sharing(layout)

    # Phase 3 (landed): Rule 4 — split narrow necks.
    _split_narrow_necks(layout, runway_axis_deg)

    # Phase 4 (TODO): Rule 3 — axis-align cut lines.
    # _axis_align_non_pavement_borders(layout, runway_axis_deg)


# ── Runway-axis helper ───────────────────────────────────────────


def longest_runway_axis_deg(layout: PavementLayout) -> Optional[float]:
    """Return the bearing (degrees mod 180) of the longest runway
    polygon's MRR long axis.  ``None`` if no runway shape is present.
    0° = +Y (north), 90° = +X (east), per the rest of the codebase
    (see ``pavement/strips.py::_linestring_bearing_axis``).
    """
    longest_poly: Optional[Polygon] = None
    longest_len = 0.0
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        p = s.polygon
        if p is None or p.is_empty or p.geom_type != "Polygon":
            continue
        # MRR long-side gives the runway's main axis even when its
        # polygon is a long thin segment.
        try:
            mrr = p.minimum_rotated_rectangle
        except Exception:
            continue
        if mrr.is_empty or mrr.geom_type != "Polygon":
            continue
        coords = list(mrr.exterior.coords)
        if len(coords) < 5:
            continue
        sides = []
        for i in range(4):
            ax, ay = coords[i]
            bx, by = coords[i + 1]
            sides.append((math.hypot(bx - ax, by - ay), (ax, ay), (bx, by)))
        sides.sort(reverse=True)
        long_len, (ax, ay), (bx, by) = sides[0]
        if long_len > longest_len:
            longest_len = long_len
            longest_poly = p
            longest_axis = (ax, ay, bx, by)
    if longest_poly is None:
        return None
    ax, ay, bx, by = longest_axis
    dx = bx - ax
    dy = by - ay
    if math.hypot(dx, dy) < 1e-6:
        return 0.0
    return math.degrees(math.atan2(dx, dy)) % 180.0


# ── Rule 2: long-edge corner snap ────────────────────────────────


def _rect_long_edges(
    rect: Polygon,
) -> List[Tuple[Tuple[float, float], Tuple[float, float],
                Tuple[float, float], Tuple[float, float]]]:
    """Return the rect's two long edges as
    ``[(p1, p2, corner_a, corner_b), ...]`` where (p1, p2) is the
    long edge endpoints and (corner_a, corner_b) are the same two
    points (corners that bound this long edge).

    Per the rect-build convention used in
    ``_clip_residue_at_stub_long_edges`` (`pavement/stubs.py:495`),
    a 4-corner rect's exterior coords are ordered so that
    ``[(coords[0], coords[1]), (coords[2], coords[3])]`` are the
    long edges (parallel to the axis) and the two short ends are
    ``[(coords[1], coords[2]), (coords[3], coords[0])]``.
    """
    coords = list(rect.exterior.coords)
    if not coords:
        return []
    if coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return []
    return [
        (coords[0], coords[1], coords[0], coords[1]),
        (coords[2], coords[3], coords[2], coords[3]),
    ]


def _point_segment_distance(
    px: float, py: float,
    ax: float, ay: float, bx: float, by: float,
) -> Tuple[float, float, float]:
    """Distance from point (px, py) to segment (a, b).  Returns
    ``(distance, foot_x, foot_y)`` where (foot_x, foot_y) is the
    closest point on the segment.
    """
    dx = bx - ax
    dy = by - ay
    seg2 = dx * dx + dy * dy
    if seg2 < 1e-9:
        return math.hypot(px - ax, py - ay), ax, ay
    t = ((px - ax) * dx + (py - ay) * dy) / seg2
    if t < 0.0:
        t = 0.0
    elif t > 1.0:
        t = 1.0
    fx = ax + t * dx
    fy = ay + t * dy
    return math.hypot(px - fx, py - fy), fx, fy


def _point_perp_dist_within_segment(
    px: float, py: float,
    ax: float, ay: float, bx: float, by: float,
) -> Optional[float]:
    """Perpendicular distance from (px, py) to the line through (a, b),
    BUT only when the foot of perpendicular falls strictly within the
    segment (0 < t < 1).  Returns ``None`` if the projection lies at
    or past either endpoint.

    Per user 2026-05-01 clarification: Rule 2's 10 m exclusion is
    perpendicular-to-axis only, and only along the rect's axial
    extent.  A junction vertex that reaches toward the short-end
    corner sits beyond the long edge's endpoint and is NOT flagged
    even though its straight-line distance to the long edge is
    small — it's connecting at the short side, not running along
    the long side.
    """
    dx = bx - ax
    dy = by - ay
    seg2 = dx * dx + dy * dy
    if seg2 < 1e-9:
        return None
    t = ((px - ax) * dx + (py - ay) * dy) / seg2
    if t <= 0.0 or t >= 1.0:
        return None
    fx = ax + t * dx
    fy = ay + t * dy
    return math.hypot(px - fx, py - fy)


def _snap_to_long_edge_corners(layout: PavementLayout) -> None:
    """Rule 2: snap each junction vertex within ``LONG_EDGE_SNAP_M``
    of any sloping-rect long edge to the nearest of that long edge's
    two endpoints (which are rect corners).  After snapping, dedupe
    consecutive identical vertices.  When ``node_altitudes`` is set
    on a junction shape, drop the altitude entry alongside the
    vertex it corresponds to so the per-vertex altitude list stays
    aligned with the polygon's ring length.
    """
    # Pre-compute long edges + endpoint corners per rect.  Read from
    # ``layout.shapes`` so we snap against the FINAL rect polygons
    # (after overlap clip + shared-vertex collapse).
    long_edges: List[Tuple[float, float, float, float,
                           Tuple[float, float], Tuple[float, float]]] = []
    for shape in layout.shapes:
        if shape.role not in SLOPING_RECT_ROLES:
            continue
        rect = shape.polygon
        if rect is None or rect.is_empty or rect.geom_type != "Polygon":
            continue
        for ax, ay, bx, by in (
            (e[0][0], e[0][1], e[1][0], e[1][1])
            for e in _rect_long_edges(rect)
        ):
            long_edges.append((ax, ay, bx, by, (ax, ay), (bx, by)))
    if not long_edges:
        return

    snap_tol = LONG_EDGE_SNAP_M
    corner_tol = SHARED_VERTEX_TOL_M

    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        coords = list(poly.exterior.coords)
        if not coords:
            continue
        # node_altitudes (when set) spans the CLOSED ring and is
        # exactly one entry per vertex including the closing repeat.
        node_alts = shape.node_altitudes
        if node_alts is not None and len(node_alts) != len(coords):
            # Length mismatch — skip altitude tracking for this shape
            # (something earlier didn't keep them in sync; safer to
            # leave the snap untouched).
            node_alts = None
        had_close = (coords[0] == coords[-1])
        if had_close:
            coords = coords[:-1]
            if node_alts is not None:
                node_alts = list(node_alts[:-1])
        # snapped[i] = (new_xy, original_alt_or_None).
        snapped: List[Tuple[Tuple[float, float], Optional[float]]] = []
        changed = False
        for i, (vx, vy) in enumerate(coords):
            best_corner: Optional[Tuple[float, float]] = None
            best_dist = snap_tol
            for ax, ay, bx, by, c1, c2 in long_edges:
                # Perpendicular distance, only within the long edge's
                # axial extent — vertices reaching toward the short
                # end (projection past either corner) are allowed.
                d = _point_perp_dist_within_segment(
                    vx, vy, ax, ay, bx, by)
                if d is None or d >= best_dist:
                    continue
                d1 = math.hypot(vx - c1[0], vy - c1[1])
                d2 = math.hypot(vx - c2[0], vy - c2[1])
                if d1 <= corner_tol or d2 <= corner_tol:
                    continue
                candidate = c1 if d1 <= d2 else c2
                best_dist = d
                best_corner = candidate
            alt = node_alts[i] if node_alts is not None else None
            if best_corner is not None:
                snapped.append((best_corner, alt))
                changed = True
            else:
                snapped.append(((vx, vy), alt))

        if not changed:
            continue
        # Dedupe consecutive identical vertices, dropping the matching
        # altitude entry.
        deduped: List[Tuple[Tuple[float, float], Optional[float]]] = []
        for entry in snapped:
            (cx, cy), _alt = entry
            if deduped:
                (px, py), _ = deduped[-1]
                if math.hypot(cx - px, cy - py) <= corner_tol:
                    continue
            deduped.append(entry)
        if len(deduped) < 3:
            continue
        new_pts = [e[0] for e in deduped]
        try:
            new_poly = Polygon(new_pts).buffer(0)
        except Exception:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon":
            continue
        shape.polygon = new_poly
        if node_alts is not None:
            new_alts = [e[1] if e[1] is not None else 0.0 for e in deduped]
            shape.node_altitudes = new_alts + [new_alts[0]]


# ── Rule 1: junction-runway 1:1 vertex sharing ───────────────────


def _enforce_runway_1to1_sharing(layout: PavementLayout) -> None:
    """Rule 1 (user 2026-05-01): when a junction polygon shares an
    edge with a runway, its runway-side vertex sequence must match
    the runway's vertex sequence on that span exactly — no extras,
    no near-misses.

    Algorithm:
      1. For each (junction, runway) pair where the junction touches
         the runway boundary at one or more contiguous spans:
      2. For each span: drop junction vertices that lie within
         ``RUNWAY_BOUNDARY_TOL_M`` of the runway boundary but DON'T
         coincide with any runway vertex (within ``SHARED_VERTEX_TOL_M``).
      3. Snap the surviving runway-near junction vertices to the
         exact runway vertex they're closest to.

    The "widen out to next runway node" behaviour described by the
    user is achieved as a side-effect: when the junction's run of
    runway-near vertices is replaced by exact runway vertex matches,
    the OUTBOARD endpoints of the run get snapped to the nearest
    runway vertex, which is by definition outboard of where the
    junction's edge used to land.

    Per-vertex altitudes are preserved by index (a snapped vertex
    keeps its existing altitude entry; dropped vertices' altitudes
    are removed in lockstep).
    """
    runway_shapes = [
        s for s in layout.shapes
        if s.role == ROLE_RUNWAY
        and s.polygon is not None
        and not s.polygon.is_empty
        and s.polygon.geom_type == "Polygon"
    ]
    if not runway_shapes:
        return

    # Build runway boundary segments PAIRED with their two endpoints
    # (which are runway vertices).  Snap targets are the segment's
    # endpoints — never a runway vertex on a different segment.
    rwy_segs: List[Tuple[float, float, float, float,
                         Tuple[float, float], Tuple[float, float]]] = []
    for s in runway_shapes:
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        m = len(coords)
        for i in range(m):
            ax, ay = coords[i]
            bx, by = coords[(i + 1) % m]
            rwy_segs.append(
                (float(ax), float(ay), float(bx), float(by),
                 (float(ax), float(ay)), (float(bx), float(by))))

    if not rwy_segs:
        return

    boundary_tol = RUNWAY_BOUNDARY_TOL_M
    vertex_tol = SHARED_VERTEX_TOL_M

    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        coords = list(poly.exterior.coords)
        if not coords:
            continue
        node_alts = shape.node_altitudes
        if node_alts is not None and len(node_alts) != len(coords):
            node_alts = None
        had_close = (coords[0] == coords[-1])
        if had_close:
            coords = coords[:-1]
            if node_alts is not None:
                node_alts = list(node_alts[:-1])

        new_pts: List[Tuple[float, float]] = []
        new_alts: Optional[List[float]] = (
            [] if node_alts is not None else None)
        changed = False
        for i, (vx, vy) in enumerate(coords):
            # Find the runway segment closest to this vertex.
            best_seg = None
            best_seg_d = float("inf")
            for ax, ay, bx, by, c1, c2 in rwy_segs:
                d, _, _ = _point_segment_distance(vx, vy, ax, ay, bx, by)
                if d < best_seg_d:
                    best_seg_d = d
                    best_seg = (c1, c2)
                    if best_seg_d <= 1e-6:
                        break
            if best_seg_d > boundary_tol or best_seg is None:
                # Not on a runway boundary; keep as-is.
                new_pts.append((vx, vy))
                if new_alts is not None:
                    new_alts.append(node_alts[i])
                continue

            # Snap to the nearer endpoint of the closest segment.
            c1, c2 = best_seg
            d1 = math.hypot(vx - c1[0], vy - c1[1])
            d2 = math.hypot(vx - c2[0], vy - c2[1])
            target = c1 if d1 <= d2 else c2
            target_d = min(d1, d2)
            if target_d <= vertex_tol:
                # Already at this corner.
                new_pts.append(target)
                if new_alts is not None:
                    new_alts.append(node_alts[i])
                if target != (vx, vy):
                    changed = True
            else:
                new_pts.append(target)
                if new_alts is not None:
                    new_alts.append(node_alts[i])
                changed = True

        if not changed:
            continue
        # Dedupe consecutive identical vertices (snap may collapse
        # adjacent vertices onto the same runway node).
        deduped_pts: List[Tuple[float, float]] = []
        deduped_alts: Optional[List[float]] = (
            [] if new_alts is not None else None)
        for j, (cx, cy) in enumerate(new_pts):
            if deduped_pts:
                px, py = deduped_pts[-1]
                if math.hypot(cx - px, cy - py) <= vertex_tol:
                    continue
            deduped_pts.append((cx, cy))
            if deduped_alts is not None:
                deduped_alts.append(new_alts[j])
        if len(deduped_pts) < 3:
            continue
        try:
            new_poly = Polygon(deduped_pts).buffer(0)
        except Exception:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon":
            continue
        shape.polygon = new_poly
        if deduped_alts is not None:
            shape.node_altitudes = deduped_alts + [deduped_alts[0]]


# ── Rule 4: split narrow necks ───────────────────────────────────


def _polygon_neck_metrics(
    poly: Polygon,
) -> Tuple[float, float, Tuple[float, float, float, float]]:
    """Return ``(min_thickness_m, mrr_long_m, (long_a_x, long_a_y,
    long_b_x, long_b_y))`` for the polygon's minimum-rotated-rectangle.

    ``min_thickness_m`` is the shorter MRR side; ``mrr_long_m`` the
    longer; the tuple gives the endpoints of one long-side segment of
    the MRR (used to set cut direction perpendicular to the MRR's
    long axis when a neck split fires).
    """
    try:
        mrr = poly.minimum_rotated_rectangle
    except Exception:
        return 0.0, 0.0, (0.0, 0.0, 0.0, 0.0)
    if mrr.is_empty or mrr.geom_type != "Polygon":
        return 0.0, 0.0, (0.0, 0.0, 0.0, 0.0)
    coords = list(mrr.exterior.coords)
    if len(coords) < 5:
        return 0.0, 0.0, (0.0, 0.0, 0.0, 0.0)
    sides: List[Tuple[float, Tuple[float, float], Tuple[float, float]]] = []
    for i in range(4):
        ax, ay = coords[i]
        bx, by = coords[i + 1]
        sides.append((math.hypot(bx - ax, by - ay), (ax, ay), (bx, by)))
    sides.sort(key=lambda x: x[0])
    short_side, _, _ = sides[0]
    long_side, la, lb = sides[-1]
    return (short_side, long_side,
            (la[0], la[1], lb[0], lb[1]))


def _split_narrow_necks(
    layout: PavementLayout,
    runway_axis_deg: Optional[float],
) -> None:
    """Rule 4 (user 2026-05-01): when a junction polygon has a narrow
    neck (MRR short-side < ``NECK_ABSOLUTE_M`` OR MRR short/long
    ratio < ``NECK_RELATIVE``), split it at the neck.  Cut direction
    is the runway axis (or perpendicular to MRR long axis if no
    runway axis is available); cut location is the MRR long-axis
    midpoint.

    Disposition (v1, simplified):
      * Both pieces with area ≥ ``MIN_JUNCTION_AREA_M2`` → keep both
        as separate junctions.
      * One piece below threshold → discard.
      * Both below threshold → original polygon kept (no productive
        split possible).

    Absorption into adjacent rect / junction (per the plan's
    ``NECK_ABSORB_FRAC`` clause) is NOT implemented in v1 — keeping
    the simple symmetric split first to avoid disturbing existing
    geometry; we can add absorption when a real case demands it.
    """
    from .config import NECK_ABSOLUTE_M, NECK_RELATIVE
    # Mirror ``junction_emit.MIN_JUNCTION_AREA_M2`` (kept local there
    # for legacy; harmonise once both files reference one constant).
    MIN_JUNCTION_AREA_M2 = 50.0
    new_shapes: List[BuiltShape] = []
    drop_indices: List[int] = []
    for idx, shape in enumerate(layout.shapes):
        if shape.role != ROLE_JUNCTION:
            continue
        poly = shape.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        short_m, long_m, _ = _polygon_neck_metrics(poly)
        if long_m <= 0.0:
            continue
        if short_m >= NECK_ABSOLUTE_M and (short_m / long_m) >= NECK_RELATIVE:
            continue
        # Compute cut: line perpendicular to runway axis (or MRR long
        # axis as fallback) passing through the polygon's centroid.
        cx = poly.centroid.x
        cy = poly.centroid.y
        if runway_axis_deg is not None:
            # Cut runs PERPENDICULAR to the runway axis (so the two
            # pieces sit on either side of the runway-aligned cut).
            cut_axis_rad = math.radians((runway_axis_deg + 90.0) % 180.0)
        else:
            short_a, long_a, mrr_long = _polygon_neck_metrics(poly)
            ax, ay, bx, by = mrr_long
            cut_axis_rad = math.atan2(by - ay, bx - ax)
        span = max(poly.bounds[2] - poly.bounds[0],
                   poly.bounds[3] - poly.bounds[1]) + 10.0
        ux = math.cos(cut_axis_rad)
        uy = math.sin(cut_axis_rad)
        cut = LineString([(cx - span * ux, cy - span * uy),
                          (cx + span * ux, cy + span * uy)])
        try:
            from shapely.ops import split as _shp_split
            result = _shp_split(poly, cut)
        except Exception:
            continue
        pieces: List[Polygon] = []
        if result.geom_type == "Polygon":
            pieces.append(result)
        else:
            for g in getattr(result, "geoms", []):
                if g.geom_type == "Polygon" and not g.is_empty:
                    pieces.append(g)
        if len(pieces) < 2:
            continue
        big_pieces = [p for p in pieces if p.area >= MIN_JUNCTION_AREA_M2]
        if len(big_pieces) < 2:
            continue
        # Replace original with the first piece; queue the rest as
        # new shapes.  Drop node_altitudes (per-vertex altitudes don't
        # transfer through a polygon split).
        shape.polygon = big_pieces[0]
        shape.node_altitudes = None
        for p in big_pieces[1:]:
            new_shapes.append(BuiltShape(polygon=p, role=ROLE_JUNCTION))
    if new_shapes:
        layout.shapes.extend(new_shapes)
