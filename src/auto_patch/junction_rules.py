"""Junction-refinement rule passes (user 2026-05-01).

Plan: ``/Users/noah/.claude/plans/kind-meandering-sifakis.md``.

Four rules apply as post-emission passes inside
``junction_emit.emit_junctions_and_finalize``:

* **Rule 1** — junction-runway 1:1 vertex sharing.
* **Rule 2** — snap junction vertices within ``SLOPING_EDGE_SNAP_M`` of
  a sloping rect's SLOPING edge to the rect's nearest cross-edge corner.
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

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

# Narrow exception tuple for shapely / numeric-geometry failure
# modes.  Programming errors propagate so they surface immediately.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)

from .config import (
    AXIS_ALIGN_TOL_DEG,
    SLIVER_ANGLE_THRESHOLD_DEG,
    SLOPING_EDGE_SNAP_M,
    NECK_ABSOLUTE_M,
    NECK_ABSORB_FRAC,
    NECK_RELATIVE,
    RUNWAY_ADJACENCY_TOL_M,
    RUNWAY_BOUNDARY_TOL_M,
)
from .layout import (
    BuiltShape,
    PavementLayout,
    ROLE_APRON,
    ROLE_JUNCTION,
    ROLE_RUNWAY,
    ROLE_TERMINAL,
    SHARED_VERTEX_TOL_M,
)


__all__ = [
    "apply_junction_rules",
    "longest_runway_axis_deg",
    "stitch_pavement_polygons",
    "stitch_pavement_to_flat_runways",
    "stitch_pavement_to_terminals",
    "widen_junctions_to_runway_corners",
]


SLOPING_RECT_ROLES = (
    "primary_parallel",
    "secondary_parallel",
    "stub",
    "cross_connector",
)


SLOPING_RECT_FLAT_THRESHOLD_M = 0.05


def _align_rect_slope_to_axis(layout: PavementLayout) -> None:
    """Per user 2026-05-02: a sloping rect with a negligible high/low
    altitude delta should be converted to FLAT (single altitude).

    Convention-free: we compare the SCALAR ``altitude_high`` and
    ``altitude_low`` directly.  We do NOT try to determine which
    physical corner is high vs low — that depends on the
    ``[altitude_high, altitude_low, altitude_low, altitude_high]``
    polygon-vertex convention, which can be rotated by post-emit
    overlap-clip / shared-vertex collapse and is therefore unsafe
    to rely on here.

    Result:
      * If |altitude_high − altitude_low| < threshold → flatten
        (set ``altitude`` to the average; clear high/low).
      * Otherwise leave the rect alone — its slope direction is
        whatever the polygon vertex order encodes; re-aligning that
        to source_axis requires a polygon-reorder pass and is
        deferred.

    Flat rects are exempt from sloping-rect connection rules
    (Rule 2 etc.) and may receive junction connections on any side.
    """
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        if s.altitude_high is None or s.altitude_low is None:
            continue  # already flat
        if abs(s.altitude_high - s.altitude_low) >= SLOPING_RECT_FLAT_THRESHOLD_M:
            continue
        s.altitude = (s.altitude_high + s.altitude_low) / 2.0
        s.altitude_high = None
        s.altitude_low = None


def apply_junction_rules(layout: PavementLayout) -> None:
    """Run the rule passes in implementation-phase order.  Mutates
    ``layout.shapes`` in place.  Skips passes whose helpers haven't
    landed yet.
    """
    runway_axis_deg = longest_runway_axis_deg(layout)

    # Phase 0 (landed): align rect slopes to source_axis.  Rects
    # whose slope is purely perpendicular to source_axis become
    # FLAT (single altitude).  Sloping-rect connection rules
    # (Rule 2) skip flat rects per user 2026-05-02.
    _align_rect_slope_to_axis(layout)

    # Phase 1 (landed): Rule 2 — sloping-edge corner snap.
    # Use the FINAL rect polygons from layout.shapes — earlier passes
    # (overlap clip, shared-vertex centroid collapse) may have shifted
    # the original ``taxi_rects`` polygons; reading from layout
    # guarantees we snap against the geometry the test sees.
    _snap_to_sloping_edge_corners(layout)
    # Phase 1b (user 2026-05-16): also snap "almost-at-the-corner"
    # junction vertices that landed on a sloping rect's FLAT (cross)
    # edge.  Conservative tolerance — only pulls vertices that are
    # within 2 m perpendicular AND within 10 m of one corner — so
    # legitimate wrap-points further from corners are unaffected.
    _snap_junction_vertices_to_rect_flat_edge_corners(layout)

    # Phase 2 (landed): Rule 1 — junction-runway 1:1 sharing.
    _enforce_runway_1to1_sharing(layout)

    # Phase 3 (landed): Rule 4 — split narrow necks.
    _split_narrow_necks(layout, runway_axis_deg)


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
        except _GEOM_EXC:
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


# ── Rule 2: sloping-edge corner snap ─────────────────────────────


def _rect_sloping_edges(
    rect: Polygon,
    source_axis: Optional[LineString] = None,
) -> List[Tuple[Tuple[float, float], Tuple[float, float],
                Tuple[float, float], Tuple[float, float]]]:
    """Return the rect's two SLOPING edges — the edges parallel to
    its source_axis (where altitude varies linearly).  These are
    the edges junctions must NOT have nodes along (other than at
    the corners), because that breaks the rect's straight-line
    slope rendering.  Per user 2026-05-02 clarification:
    "long" vs "short" was misleading — what matters is sloping vs
    flat.  A rect can be wider than long and still have its slope
    along the short axis.

    Detection:
      * If ``source_axis`` provided: compute the absolute dot
        product between each edge direction and the axis direction;
        the 2 edges with the highest dot are the most parallel =
        sloping edges.
      * Fallback (no axis): use the 2 longest edges by length —
        works for typical sloping rects where the long dimension
        is the slope direction.

    Return format: ``[(p1, p2, corner_a, corner_b), ...]``.
    """
    coords = list(rect.exterior.coords)
    if not coords:
        return []
    if coords[0] == coords[-1]:
        coords = coords[:-1]
    if len(coords) != 4:
        return []
    edges = [(coords[i], coords[(i + 1) % 4]) for i in range(4)]
    if source_axis is not None and not source_axis.is_empty:
        ax_pts = list(source_axis.coords)
        if len(ax_pts) >= 2:
            axdx = ax_pts[-1][0] - ax_pts[0][0]
            axdy = ax_pts[-1][1] - ax_pts[0][1]
            axlen = math.hypot(axdx, axdy)
            if axlen >= 1e-6:
                aux, auy = axdx / axlen, axdy / axlen
                dots = []
                for a, b in edges:
                    ex, ey = b[0] - a[0], b[1] - a[1]
                    elen = math.hypot(ex, ey)
                    if elen < 1e-6:
                        dots.append(0.0)
                        continue
                    dots.append(abs(ex * aux + ey * auy) / elen)
                # The 2 edges with HIGHEST absolute dot are most
                # parallel to the axis = sloping.
                sloping_idx = sorted(
                    range(4), key=lambda i: -dots[i])[:2]
                return [(edges[i][0], edges[i][1],
                         edges[i][0], edges[i][1])
                        for i in sloping_idx]
    # Fallback: pick the 2 longest edges (typical heuristic).
    lengths = [math.hypot(b[0] - a[0], b[1] - a[1])
               for a, b in edges]
    long_idx = sorted(range(4), key=lambda i: -lengths[i])[:2]
    return [(edges[i][0], edges[i][1], edges[i][0], edges[i][1])
            for i in long_idx]


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
    corner sits beyond the sloping edge's endpoint and is NOT flagged
    even though its straight-line distance to the sloping edge is
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


def _snap_junction_vertices_to_rect_flat_edge_corners(
    layout: PavementLayout,
    perp_tol_m: float = 2.0,
    corner_max_m: float = 10.0,
) -> None:
    """Snap each junction vertex that lies within ``perp_tol_m``
    perpendicular of a sloping rect's FLAT (cross / short) edge AND
    within ``corner_max_m`` of one of that edge's two corners to
    the nearer corner.

    Per user 2026-05-16: sloping rects share their flat edge 1:1 —
    only the two corners are legal shared vertices.  A junction
    vertex sitting on the flat-edge interior (e.g. the CYXY apron
    boundary vertex (-348.79, 75.05) that's 1 m perpendicular and
    5.5 m from A2's corner (-344.85, 71.16)) violates that invariant
    and produces an elevation step at the boundary.

    Conservative tolerance: only snap "almost-at-the-corner" cases
    (perpendicular ≤ 2 m, corner ≤ 10 m).  Vertices that are
    perpendicular-close but corner-far are legitimate wrap-points
    where the junction polygon walks AROUND the rect via its flat
    edge — leave those alone (per the comment in
    ``_snap_to_sloping_edge_corners`` reverting commit 5a50d00's
    indiscriminate cross-edge snapping).
    """
    flat_edges: List[Tuple[float, float, float, float]] = []
    rect_corner_set: set = set()
    bucket = SHARED_VERTEX_TOL_M
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        if (s.altitude_high is None or s.altitude_low is None):
            continue
        rect = s.polygon
        if rect is None or rect.is_empty \
                or rect.geom_type != "Polygon":
            continue
        rc = list(rect.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Flat edges: corners (1,2) and (3,0) per
        # _rect_from_axis_extended convention.
        for ci_a, ci_b in ((1, 2), (3, 0)):
            ax, ay = rc[ci_a]
            bx, by = rc[ci_b]
            flat_edges.append((float(ax), float(ay),
                                float(bx), float(by)))
        for cx, cy in rc:
            rect_corner_set.add(
                (round(cx / bucket), round(cy / bucket)))
    if not flat_edges:
        return

    for shape in layout.shapes:
        if shape.role != ROLE_JUNCTION:
            continue
        poly = shape.polygon
        if poly is None or poly.is_empty \
                or poly.geom_type != "Polygon":
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
        if len(coords) < 3:
            continue
        new_coords: List[Tuple[float, float]] = []
        new_alts: Optional[List[float]] = (
            [] if node_alts is not None else None)
        changed = False
        for i, (vx, vy) in enumerate(coords):
            # If the vertex is already at a rect corner, keep as is.
            key = (round(vx / bucket), round(vy / bucket))
            if key in rect_corner_set:
                new_coords.append((vx, vy))
                if new_alts is not None:
                    new_alts.append(node_alts[i])
                continue
            best_corner: Optional[Tuple[float, float]] = None
            best_dc = corner_max_m
            for ax, ay, bx, by in flat_edges:
                d_perp = _point_perp_dist_within_segment(
                    vx, vy, ax, ay, bx, by)
                if d_perp is None or d_perp > perp_tol_m:
                    continue
                d_a = math.hypot(vx - ax, vy - ay)
                d_b = math.hypot(vx - bx, vy - by)
                near, dc = ((ax, ay), d_a) if d_a < d_b else (
                    (bx, by), d_b)
                if dc < best_dc:
                    best_dc = dc
                    best_corner = near
            if best_corner is not None:
                new_coords.append(best_corner)
                changed = True
                if new_alts is not None:
                    new_alts.append(node_alts[i])
            else:
                new_coords.append((vx, vy))
                if new_alts is not None:
                    new_alts.append(node_alts[i])
        if not changed:
            continue
        # Dedup consecutive identical vertices.
        deduped: List[Tuple[float, float]] = []
        deduped_alts: Optional[List[float]] = (
            [] if new_alts is not None else None)
        for j, (cx, cy) in enumerate(new_coords):
            if deduped:
                px, py = deduped[-1]
                if math.hypot(cx - px, cy - py) <= bucket:
                    continue
            deduped.append((cx, cy))
            if deduped_alts is not None:
                deduped_alts.append(new_alts[j])
        if len(deduped) < 3:
            continue
        try:
            new_poly = Polygon(deduped).buffer(0)
        except _GEOM_EXC:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon" or not new_poly.is_valid:
            continue
        # Reject if area collapsed significantly (snap created a
        # degenerate / self-intersecting shape).
        if new_poly.area < 0.5 * poly.area:
            continue
        shape.polygon = new_poly
        if deduped_alts is not None:
            shape.node_altitudes = deduped_alts + [deduped_alts[0]]


def _snap_to_sloping_edge_corners(layout: PavementLayout) -> None:
    """Rule 2: snap each junction vertex within ``SLOPING_EDGE_SNAP_M``
    of a sloping-rect's SLOPING edge to the nearest of that edge's
    two endpoints (which are rect corners).  After snapping, dedup
    consecutive identical vertices.

    Per user 2026-05-04 (regression analysis): the rule operates on
    SLOPING edges only, NOT cross edges.  Snapping on cross edges
    pulls junction vertices that legitimately sit BETWEEN a rect's
    cross-edge corners — and which the polygon walks AROUND the rect
    via — onto one of the corners, distorting the polygon's wrap
    path.  At ``best-elevation-model`` (commit 187d2cd) total sloping-
    rect-vs-junction overlap was 0 m²; the all-4-edges variant in
    commit 5a50d00 introduced the V3 stub overlap, which was then
    band-aided with a "drop intervening verts" rule that produced
    self-intersecting polygons → buffer(0) splits → slivers along
    sloping edges.  Both changes reverted here.

    Flat-rect gate: skip rects whose altitude_high/altitude_low
    aren't set yet (which means they're either flat or per-surface-
    solver hasn't run; in both cases the corner-only sharing rule
    is unnecessary).  Re-runs after the post-elevation phase pick up
    rects whose altitudes are now set.

    When ``node_altitudes`` is set on a junction shape, drop the
    altitude entry alongside the vertex it corresponds to so the
    per-vertex altitude list stays aligned with the polygon's ring
    length.
    """
    # Pre-compute sloping edges per rect.  Read from layout.shapes so
    # we snap against the FINAL rect polygons (after overlap clip +
    # shared-vertex collapse).  The rect_idx + corner_idx tuple is
    # carried for downstream adjacency-aware processing (e.g. dedup-
    # consecutive that crosses through identical-corner snap targets).
    rect_corners_per_rect: List[List[Tuple[float, float]]] = []
    rect_edges: List[Tuple[float, float, float, float, int, int, int]] = []
    for shape in layout.shapes:
        if shape.role not in SLOPING_RECT_ROLES:
            continue
        # Flat-rect exemption (restored from best-elevation-model).
        # Without an altitude pair, the rect either hasn't been
        # elevated yet or is intentionally flat — in both cases
        # junctions can connect anywhere on its boundary without
        # breaking the linear-slope rendering invariant.
        if (shape.altitude_high is None
                or shape.altitude_low is None):
            continue
        rect = shape.polygon
        if rect is None or rect.is_empty or rect.geom_type != "Polygon":
            continue
        rc = list(rect.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        rect_idx = len(rect_corners_per_rect)
        rect_corners_per_rect.append(
            [(float(c[0]), float(c[1])) for c in rc])
        # Sloping edges only — corners (0,1) and (2,3) per the
        # _rect_from_axis_extended convention.
        for ci_a, ci_b in ((0, 1), (2, 3)):
            ax, ay = rc[ci_a]
            bx, by = rc[ci_b]
            rect_edges.append(
                (float(ax), float(ay), float(bx), float(by),
                 rect_idx, ci_a, ci_b))
    if not rect_edges:
        return

    snap_tol = SLOPING_EDGE_SNAP_M
    corner_tol = SHARED_VERTEX_TOL_M

    def _corner_id_for(vx: float, vy: float
                        ) -> Optional[Tuple[int, int]]:
        """If (vx, vy) coincides with a rect corner (within
        ``corner_tol``), return ``(rect_idx, corner_idx)``; else None.
        Used for adjacency detection on already-snapped vertices."""
        for r_idx, corners in enumerate(rect_corners_per_rect):
            for c_idx, (cx, cy) in enumerate(corners):
                if math.hypot(vx - cx, vy - cy) <= corner_tol:
                    return (r_idx, c_idx)
        return None

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
        # snapped[i] = (new_xy, original_alt_or_None, corner_id_or_None).
        snapped: List[Tuple[Tuple[float, float], Optional[float],
                              Optional[Tuple[int, int]]]] = []
        changed = False
        for i, (vx, vy) in enumerate(coords):
            best_corner: Optional[Tuple[float, float]] = None
            best_corner_id: Optional[Tuple[int, int]] = None
            best_dist = snap_tol
            for ax, ay, bx, by, r_idx, ci_a, ci_b in rect_edges:
                # Perpendicular distance, only within the edge's
                # axial extent — vertices reaching toward a corner
                # past either endpoint are allowed (those reach to
                # an adjacent rect's edge legitimately).
                d = _point_perp_dist_within_segment(
                    vx, vy, ax, ay, bx, by)
                if d is None or d >= best_dist:
                    continue
                c1 = (ax, ay)
                c2 = (bx, by)
                d1 = math.hypot(vx - c1[0], vy - c1[1])
                d2 = math.hypot(vx - c2[0], vy - c2[1])
                if d1 <= corner_tol or d2 <= corner_tol:
                    continue
                if d1 <= d2:
                    candidate = c1
                    candidate_id = (r_idx, ci_a)
                else:
                    candidate = c2
                    candidate_id = (r_idx, ci_b)
                best_dist = d
                best_corner = candidate
                best_corner_id = candidate_id
            alt = node_alts[i] if node_alts is not None else None
            if best_corner is not None:
                snapped.append((best_corner, alt, best_corner_id))
                changed = True
            else:
                # Pre-existing corner coincidence (e.g. originally on
                # a corner): tag with corner id so adjacency detection
                # below sees it.
                snapped.append(((vx, vy), alt,
                                _corner_id_for(vx, vy)))

        if not changed:
            # Even if no NEW snap, we may still need to drop
            # intervening vertices between pre-existing corner
            # coincidences (pre-snap geometry already had two corners
            # in the polygon).  Continue if no corner ids found.
            if not any(e[2] is not None for e in snapped):
                continue
        # Dedupe consecutive identical vertices, dropping the matching
        # altitude entry.
        deduped: List[Tuple[Tuple[float, float], Optional[float],
                              Optional[Tuple[int, int]]]] = []
        for entry in snapped:
            (cx, cy), _alt, _cid = entry
            if deduped:
                (px, py), _, _ = deduped[-1]
                if math.hypot(cx - px, cy - py) <= corner_tol:
                    continue
            deduped.append(entry)
        # Note: a "drop intervening corner-snapped vertices" rule was
        # added in commit 5a50d00 to band-aid a V3 stub overlap that
        # the all-4-edges snap (also from 5a50d00) had introduced.
        # Both reverted here per regression analysis.  Sloping-edges-
        # only snap doesn't produce the V3 overlap, so the band-aid
        # isn't needed.
        if len(deduped) < 3:
            continue
        new_pts = [e[0] for e in deduped]
        try:
            new_poly = Polygon(new_pts).buffer(0)
        except _GEOM_EXC:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon":
            continue
        shape.polygon = new_poly
        if node_alts is not None:
            new_alts = [e[1] for e in deduped]
            shape.node_altitudes = new_alts + [new_alts[0]]


# ── Rule 1: junction-runway 1:1 vertex sharing ───────────────────


def _build_runway_union_chain(
    runway_shapes: Sequence[BuiltShape],
) -> Tuple[List[Tuple[float, float]],
           dict]:
    """Walk the runway-union boundary (one continuous loop per
    contiguous runway component) and return:
      * ``chain`` — list of corner positions in walk order.  For
        a multi-component layout (e.g. parallel runways) all
        components are concatenated; segment-internal seams are
        skipped because the union eliminates them.
      * ``corner_index`` — dict mapping bucketed (round to ~0.5 m)
        corner key → its position in ``chain``.
    """
    if not runway_shapes:
        return [], {}
    polys = [s.polygon for s in runway_shapes
             if s.polygon is not None and not s.polygon.is_empty]
    if not polys:
        return [], {}
    try:
        union = unary_union(polys)
    except _GEOM_EXC:
        return [], {}
    components: List[Polygon] = []
    if union.geom_type == "Polygon":
        components.append(union)
    else:
        for g in getattr(union, "geoms", []):
            if g.geom_type == "Polygon" and not g.is_empty:
                components.append(g)
    chain: List[Tuple[float, float]] = []
    corner_index: dict = {}
    for poly in components:
        coords = list(poly.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for c in coords:
            key = (round(c[0] * 2.0), round(c[1] * 2.0))
            corner_index[key] = len(chain)
            chain.append((float(c[0]), float(c[1])))
    return chain, corner_index


def _build_runway_corner_altitudes(
    runway_shapes: Sequence[BuiltShape],
) -> dict:
    """Map bucketed runway-corner key → altitude.  Per the convention
    in ``triangulation.py:165-170``: ``coords[0]`` and ``coords[3]``
    are at the ``altitude_high`` end; ``coords[1]`` and ``coords[2]``
    are at the ``altitude_low`` end.  Seam corners assigned twice
    (once from each adjacent segment) — the values should match
    because adjacent segments slope continuously through the seam.
    """
    out: dict = {}
    for s in runway_shapes:
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) < 4:
            continue
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            # Sloped 4-corner rect — strict convention.
            if len(coords) != 4:
                continue
            corner_alts = (
                (coords[0], float(s.altitude_high)),
                (coords[1], float(s.altitude_low)),
                (coords[2], float(s.altitude_low)),
                (coords[3], float(s.altitude_high)),
            )
        elif s.altitude is not None:
            # Flat shape: any number of corners (multi-node flat
            # runways carry intermediate snap corners along the long
            # sides).
            a = float(s.altitude)
            corner_alts = tuple((c, a) for c in coords)
        else:
            continue
        for c, alt in corner_alts:
            key = (round(c[0] * 2.0), round(c[1] * 2.0))
            # Take the FIRST encounter; matching with the second
            # encounter is asserted via the elevation pipeline's
            # continuity check.
            if key not in out:
                out[key] = alt
    return out


def widen_junctions_to_runway_corners(
    layout: PavementLayout,
) -> None:
    """Public entrypoint: widen each junction's runway-shared
    vertices by inserting the immediately-adjacent corners along
    the runway-union boundary.  Called POST-ELEVATION (after the
    runway is segmented and altitudes are committed) so:
      * the chain has all the segment seam corners
      * we can pull each new junction vertex's altitude from the
        runway shape's ``altitude_high`` / ``altitude_low``

    Per-junction cap of 4 runway-shared nodes (user 2026-05-02
    spec: 2, 3, or 4 nodes).
    """
    runway_shapes = [s for s in layout.shapes
                     if s.role == ROLE_RUNWAY
                     and s.polygon is not None
                     and not s.polygon.is_empty
                     and s.polygon.geom_type == "Polygon"]
    if not runway_shapes:
        return
    chain, corner_index = _build_runway_union_chain(runway_shapes)
    corner_alt = _build_runway_corner_altitudes(runway_shapes)
    _widen_runway_shared_corners(layout, chain, corner_index, corner_alt)


def _widen_runway_shared_corners(
    layout: PavementLayout,
    chain: Sequence[Tuple[float, float]],
    corner_index: dict,
    corner_alt: dict,
) -> None:
    # Pre-compute runway union for overlap rejection (per user
    # 2026-05-02: junctions and runways must never overlap; if a
    # widening insertion would extend the polygon body across the
    # runway boundary, revert the insertion).
    runway_polys = [s.polygon for s in layout.shapes
                    if s.role == ROLE_RUNWAY
                    and s.polygon is not None
                    and not s.polygon.is_empty]
    try:
        runway_union = unary_union(runway_polys) if runway_polys else None
    except _GEOM_EXC:
        runway_union = None
    pav_union = getattr(layout, "_apt_pav_union", None)
    return _do_widen(
        layout, chain, corner_index, corner_alt,
        runway_union, pav_union)


def _do_widen(
    layout: PavementLayout,
    chain: Sequence[Tuple[float, float]],
    corner_index: dict,
    corner_alt: dict,
    runway_union,
    pav_union=None,
) -> None:
    """Rule 1 v6 widening (user 2026-05-02): for each junction with
    at least one runway-shared vertex, insert the immediately-
    adjacent runway corners (one on each side, walking the runway
    union boundary) as new junction vertices.  The polygon may
    grow a thin arm extending along the runway — this is acceptable
    per the user direction.

    New vertices' altitudes are taken from the runway corner's
    altitude, providing smooth elevation continuity.
    """
    if not chain or len(chain) < 2:
        return
    n_chain = len(chain)
    vertex_tol = SHARED_VERTEX_TOL_M

    def _key(p):
        return (round(p[0] * 2.0), round(p[1] * 2.0))

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
        n = len(coords)
        if n < 3:
            continue

        # Snapshot the pre-widen polygon body for the body-
        # proximity cap (per user 2026-05-05 followup, implemented
        # 2026-05-16).  Multi-step walking + bumped cap let the
        # queue extend the runway-shared edge by 4-6 chain corners;
        # without a geometric cap, polygons whose body is small but
        # adjacent to a long runway-side stretch end up sharing 6+
        # corners that span 200+ m of runway, much further than the
        # polygon's natural extent.  Reject any chain neighbor whose
        # distance to the original body exceeds
        # ``WIDEN_BODY_PROX_M``.
        #
        # Per user 2026-05-16 (SPLP regression): with the runway
        # pre-cut at every apt.dat-pavement projection (commit
        # 15d3a89), the body's anchors already sit at the natural
        # extent of the junction along the runway.  Walking one chain
        # step further drops the corner ~30 m past the body's
        # runway-side edge — i.e. > WIDEN_BODY_PROX_M from the body,
        # so the prox cap correctly rejects it.  The runway-end-wrap
        # case (SPJC J-10131) is preserved because the wrap body's
        # boundary stays close to the chain corners along the wrap,
        # so the cap still permits the walk.
        original_body = poly
        WIDEN_BODY_PROX_M = 5.0

        # Identify runway-shared vertices in the polygon.
        # shared_in_poly: list of (poly_idx, chain_corner_position)
        existing_keys = set(_key(v) for v in coords)
        shared_in_poly: List[Tuple[int, Tuple[float, float]]] = []
        for i, v in enumerate(coords):
            k = _key(v)
            if k in corner_index:
                shared_in_poly.append((i, chain[corner_index[k]]))

        if not shared_in_poly:
            continue

        # User 2026-05-02 spec: 2-4 runway-shared nodes per junction
        # in the typical case.  Bumped to 5 per user 2026-05-02
        # follow-up for cases where 5 connections are needed (e.g.
        # diagonal stubs converging plus an internal runway-seam
        # corner falling within the joining region).  Bumped again
        # to 7 per user 2026-05-05 followup: with the apt.dat-aware
        # segmenter (commit 15d3a89), the chain has more corners
        # adjacent to each junction's runway-side edge, and the
        # 2-round multi-step walker can naturally reach 6+ chain
        # corners total (2 originals + 2 round-0 inserts + 2
        # round-1 inserts).  At SPJC J-10131 / way -10138 the cap-5
        # cut off the round-1 walk from A[8] to B[9], leaving the
        # Northern corner of taxiway C uncovered.
        max_total_shared = 7
        current_shared_count = len(shared_in_poly)
        max_inserts = max(0, max_total_shared - current_shared_count)
        if max_inserts == 0:
            continue

        # Per user 2026-05-04: validate each insertion individually
        # and commit one at a time.  The previous batch validation
        # rejected ALL insertions whenever ANY one of them caused a
        # geometric problem — at SPJC's 34R end the south-side chain
        # neighbor wraps around the runway end (1600 m² overlap), so
        # the legitimate north-side widening was getting thrown out
        # alongside it.
        current_coords = list(coords)
        current_alts: Optional[List[float]] = (
            list(node_alts) if node_alts is not None else None)
        n_committed = 0

        def _attempt_insert(insert_at, neighbor, alt):
            nonlocal current_coords, current_alts, existing_keys
            # Refuse to widen when the runway-corner altitude is
            # unknown but the polygon tracks per-vertex altitudes —
            # fabricating 0.0 would silently introduce a sea-level
            # vertex into a junction the solver then HARD-anchors.
            if current_alts is not None and alt is None:
                return False
            trial = list(current_coords)
            trial.insert(insert_at, neighbor)
            if len(trial) < 3:
                return False
            try:
                trial_poly = Polygon(trial).buffer(0)
            except _GEOM_EXC:
                return False
            if trial_poly.is_empty:
                return False
            if trial_poly.geom_type == "MultiPolygon":
                trial_poly = max(trial_poly.geoms, key=lambda g: g.area)
            if trial_poly.geom_type != "Polygon":
                return False
            if not trial_poly.is_valid or not trial_poly.is_simple:
                return False
            if trial_poly.area < 0.5 * poly.area:
                return False
            if runway_union is not None and not runway_union.is_empty:
                try:
                    ovl = trial_poly.intersection(runway_union).area
                except _GEOM_EXC:
                    ovl = 0.0
                if ovl > 1.0:
                    return False
            for other in layout.shapes:
                if other is shape or other.role != ROLE_JUNCTION:
                    continue
                if other.polygon is None or other.polygon.is_empty:
                    continue
                try:
                    if trial_poly.intersection(
                            other.polygon).area > 1.0:
                        return False
                except _GEOM_EXC:
                    pass
            # Sliver-corner pre-emption (per user 2026-05-05): the
            # OSM emitter drops any polygon with an interior angle
            # below ``SLIVER_ANGLE_THRESHOLD_DEG`` (≈ 2°) at emit
            # time.  Reject the insert here so the caller has the
            # option to skip the chain step entirely instead of
            # silently losing the whole junction downstream.
            try:
                ring_pts = list(trial_poly.exterior.coords)
                if ring_pts and ring_pts[0] == ring_pts[-1]:
                    ring_pts = ring_pts[:-1]
                m_ring = len(ring_pts)
                if m_ring >= 3:
                    sliver_cos = math.cos(
                        math.radians(SLIVER_ANGLE_THRESHOLD_DEG))
                    for vi in range(m_ring):
                        ax, ay = ring_pts[(vi - 1) % m_ring]
                        bx, by = ring_pts[vi]
                        cx, cy = ring_pts[(vi + 1) % m_ring]
                        v1x, v1y = ax - bx, ay - by
                        v2x, v2y = cx - bx, cy - by
                        n1 = math.hypot(v1x, v1y)
                        n2 = math.hypot(v2x, v2y)
                        if n1 < 1e-9 or n2 < 1e-9:
                            continue
                        cos = (v1x * v2x + v1y * v2y) / (n1 * n2)
                        if cos > sliver_cos:
                            return False
            except _GEOM_EXC:
                pass
            # Commit
            current_coords = trial
            if current_alts is not None:
                current_alts.insert(insert_at, alt)
            existing_keys.add(_key(neighbor))
            return True

        # Per user 2026-05-05 (after apt.dat-aware segmenter
        # landed): walk one chain step from each runway-shared
        # corner, INCLUDING newly-inserted ones, capped at 2
        # rounds.  With the segmenter now placing seam corners at
        # apt.dat-pavement intersections, the chain is denser; a
        # single-step walk from each ORIGINAL anchor reaches only
        # the immediately-adjacent (often very close) seam corners,
        # which doesn't extend the polygon to the runway corners
        # the user actually wants shared.  Walking one step from
        # each newly-inserted corner reaches one chain step
        # further, which corresponds to the natural extent of the
        # polygon's runway-shared edge.
        # Multi-step walking gate (per user 2026-05-05 followup):
        # round 0 always runs (single-step from each original
        # anchor).  Round 1+ (multi-step) runs only when the
        # polygon has a non-anchor vertex that's substantially
        # off the runway boundary — geometrically, this is the
        # runway-end-wrap case (J-10131 at SPJC) where the
        # polygon body extends well off the runway and we need
        # to wrap further along the runway to reach the corner
        # at the body's far edge.  When all non-anchor verts are
        # already adjacent to runway corners (the typical
        # apt.dat-aware-segmenter case at SPJC -10137), the
        # polygon is naturally aligned and round-1 walks would
        # extend the runway-shared edge past the body's natural
        # extent.
        # Distance from each non-anchor vertex to the nearest
        # runway-shared anchor in the same polygon: if all non-
        # anchors are within ``NON_ANCHOR_NEAR_ANCHOR_M`` of an
        # anchor, the polygon is tight against the runway with no
        # interior reach (the apt.dat-aware-segmenter case where
        # the polygon naturally aligns with runway corners — no
        # widening needed past round 0).  Otherwise the polygon
        # has a body that extends off the runway, and round 1
        # multi-step walking lets us wrap further along the
        # runway to reach the corner aligned with the body's far
        # extent (the runway-end-wrap case at SPJC J-10131).
        NON_ANCHOR_NEAR_ANCHOR_M = 70.0
        far_non_anchor_exists = False
        if shared_in_poly:
            anchor_pts = [c for _, c in shared_in_poly]
            for v in coords:
                if any(math.hypot(v[0] - a[0], v[1] - a[1]) < 1e-3
                       for a in anchor_pts):
                    continue
                d = min(math.hypot(v[0] - a[0], v[1] - a[1])
                        for a in anchor_pts)
                if d > NON_ANCHOR_NEAR_ANCHOR_M:
                    far_non_anchor_exists = True
                    break
        WIDEN_MAX_ROUNDS = 2 if far_non_anchor_exists else 1
        # Skip widening entirely when the polygon's anchors are
        # chain-adjacent (a contiguous run of chain corners) AND
        # no non-anchor reaches off the runway: the polygon is
        # already naturally aligned with runway corners and any
        # round-0 walk would extend its runway-shared edge past
        # its body extent (the SPJC -10137 case).
        if not far_non_anchor_exists and len(shared_in_poly) >= 2:
            anchor_cis = sorted(
                corner_index[_key(c)] for _, c in shared_in_poly
                if _key(c) in corner_index)
            if anchor_cis:
                contiguous = True
                for i in range(len(anchor_cis) - 1):
                    gap = (anchor_cis[i + 1] - anchor_cis[i]) % n_chain
                    if gap != 1:
                        contiguous = False
                        break
                if contiguous:
                    continue  # skip widening for this junction
        round_processed_keys = set(_key(c) for _, c in shared_in_poly)
        # Track which corners were newly inserted by widening, so
        # the post-widen interior-vert prune below knows which
        # adjacent non-anchor verts to consider redundant.
        newly_inserted_keys: set = set()
        widen_queue: List[Tuple[Tuple[float, float], int]] = [
            (c, 0) for _, c in shared_in_poly]

        while widen_queue and n_committed < max_inserts:
            corner, rnd = widen_queue.pop(0)
            ci = corner_index.get(_key(corner))
            if ci is None:
                continue
            chain_neighbors = (
                chain[(ci - 1) % n_chain],
                chain[(ci + 1) % n_chain],
            )
            # Find current poly_idx of this corner (insertions before
            # it shift its index forward by 1 each).
            try:
                poly_idx = next(
                    i for i, v in enumerate(current_coords)
                    if _key(v) == _key(corner))
            except StopIteration:
                continue
            n_cur = len(current_coords)
            prev_v = current_coords[(poly_idx - 1) % n_cur]
            next_v = current_coords[(poly_idx + 1) % n_cur]

            for neighbor in chain_neighbors:
                if n_committed >= max_inserts:
                    break
                if _key(neighbor) in existing_keys:
                    continue
                # Body-prox cap: the pre-widen polygon body already
                # marks the junction's natural extent along the
                # runway because the segmenter pre-cuts the runway
                # at every adjacent-pavement projection.  Skip chain
                # neighbors more than WIDEN_BODY_PROX_M from the
                # pre-widen body — those sit beyond the pavement the
                # junction actually borders.  Runway-end-wrap bodies
                # (SPJC J-10131) keep their chain corners close to
                # the wrapped boundary so the cap still permits the
                # walk.
                try:
                    if (original_body.distance(Point(neighbor))
                            > WIDEN_BODY_PROX_M):
                        continue
                except _GEOM_EXC:
                    pass
                # Decide insertion side: BEFORE poly_idx or AFTER.
                # Pick the side whose angle is closer to the
                # neighbor's bearing from the corner.
                ax_n, ay_n = (neighbor[0] - corner[0],
                              neighbor[1] - corner[1])
                ax_p, ay_p = (prev_v[0] - corner[0],
                              prev_v[1] - corner[1])
                ax_x, ay_x = (next_v[0] - corner[0],
                              next_v[1] - corner[1])
                a_n = math.atan2(ay_n, ax_n)
                a_p = math.atan2(ay_p, ax_p)
                a_x = math.atan2(ay_x, ax_x)

                def _ang_diff(a, b):
                    d = abs(a - b) % (2 * math.pi)
                    return min(d, 2 * math.pi - d)

                # Pick the side that produces the LESS-acute polygon
                # corner at the inserted neighbor (per user 2026-05-05
                # followup).  When both flanks are roughly opposite
                # the neighbor's bearing (multi-step walks where the
                # polygon's runway-shared edge sits between two near-
                # 180° flanks), the bearing-closeness tie-breaker
                # picks one side arbitrarily and the U-turn check
                # rejects the resulting near-spike — even though the
                # OTHER side would have produced an acceptable
                # corner.  Compute cos_turn for both sides; pick
                # whichever is higher (= more obtuse, less spike).
                def _side_cos(insert_at_b, flank_b, side_b):
                    if side_b == "before":
                        e1 = (neighbor[0] - flank_b[0],
                              neighbor[1] - flank_b[1])
                        e2 = (corner[0] - neighbor[0],
                              corner[1] - neighbor[1])
                    else:
                        e1 = (neighbor[0] - corner[0],
                              neighbor[1] - corner[1])
                        e2 = (flank_b[0] - neighbor[0],
                              flank_b[1] - neighbor[1])
                    m1 = math.hypot(*e1)
                    m2 = math.hypot(*e2)
                    if m1 < 1e-6 or m2 < 1e-6:
                        return None, e1, e2
                    return ((e1[0]*e2[0] + e1[1]*e2[1]) / (m1*m2),
                            e1, e2)
                cos_before, _, _ = _side_cos(poly_idx, prev_v, "before")
                cos_after, _, _ = _side_cos(poly_idx + 1, next_v,
                                            "after")
                if (cos_after is not None
                        and (cos_before is None
                             or cos_after > cos_before)):
                    insert_at = poly_idx + 1
                    flank_v = next_v
                    cos_turn = cos_after
                else:
                    insert_at = poly_idx
                    flank_v = prev_v
                    cos_turn = cos_before
                if cos_turn is not None:
                    # Per user 2026-05-04: only reject TRUE U-turns
                    # (cos < -0.99, > 172°).  The -0.95 threshold
                    # rejected legitimate widenings at SPJC's 34R
                    # west-side junction (-10130) where the angle
                    # was 162° (cos = -0.954) — a sharp but valid
                    # widening arm.  The per-insertion validity +
                    # overlap guards in ``_attempt_insert`` already
                    # catch the geometrically degenerate cases.
                    if cos_turn < -0.99:
                        continue

                alt = corner_alt.get(_key(neighbor))
                if _attempt_insert(insert_at, neighbor, alt):
                    n_committed += 1
                    newly_inserted_keys.add(_key(neighbor))
                    # Push the newly-inserted neighbor for one more
                    # round of walking, unless we've already hit the
                    # round cap.  ``round_processed_keys`` guards
                    # against re-walking the same corner twice.
                    if (rnd + 1 < WIDEN_MAX_ROUNDS
                            and _key(neighbor) not in round_processed_keys):
                        round_processed_keys.add(_key(neighbor))
                        widen_queue.append((neighbor, rnd + 1))
                    # poly_idx may shift if we inserted before it.
                    n_cur = len(current_coords)
                    try:
                        poly_idx = next(
                            i for i, v in enumerate(current_coords)
                            if _key(v) == _key(corner))
                    except StopIteration:
                        break
                    prev_v = current_coords[(poly_idx - 1) % n_cur]
                    next_v = current_coords[(poly_idx + 1) % n_cur]

        if n_committed == 0:
            continue

        # Targeted interior-vert prune (per user 2026-05-05
        # followup): when widening inserts a chain corner, the
        # polygon vertex IMMEDIATELY ADJACENT to that corner in
        # walk-order may be a leftover interior vertex from the
        # pre-widen polygon's path through pavement interior.
        # If that adjacent vert is NOT itself an anchor (chain
        # corner / rect corner / shape boundary point) and sits
        # > INTERIOR_PRUNE_M from any apt.dat-pavement /
        # runway-union boundary, drop it — the new arm to the
        # chain corner supersedes the old interior detour.
        # (Example at SPJC -10138: A9_INT vertex sits 14 m inside
        # pav_union after widen reaches B[9] — drop it so the
        # polygon walks A[8] → B[9] → A[0] cleanly.)
        INTERIOR_PRUNE_M = 5.0
        on_pav_boundary_obj = None
        if pav_union is not None and not pav_union.is_empty:
            try:
                on_pav_combined = (unary_union([pav_union, runway_union])
                                   if runway_union is not None
                                   and not runway_union.is_empty
                                   else pav_union)
                on_pav_boundary_obj = on_pav_combined.boundary
            except _GEOM_EXC:
                on_pav_boundary_obj = None
        if on_pav_boundary_obj is not None and newly_inserted_keys:
            pruned_coords: List[Tuple[float, float]] = []
            pruned_alts: Optional[List[float]] = (
                [] if current_alts is not None else None)
            n_cur_pre = len(current_coords)
            for i, v in enumerate(current_coords):
                k = _key(v)
                # Always keep chain corners and verts widening
                # just inserted.
                if k in corner_index or k in newly_inserted_keys:
                    pruned_coords.append(v)
                    if pruned_alts is not None:
                        pruned_alts.append(current_alts[i])
                    continue
                # Check if THIS vert is adjacent in walk order to a
                # newly-inserted chain corner.
                prev_k = _key(current_coords[(i - 1) % n_cur_pre])
                next_k = _key(current_coords[(i + 1) % n_cur_pre])
                adj_to_new = (prev_k in newly_inserted_keys
                              or next_k in newly_inserted_keys)
                if not adj_to_new:
                    pruned_coords.append(v)
                    if pruned_alts is not None:
                        pruned_alts.append(current_alts[i])
                    continue
                # Interior check: drop only if the vert sits >
                # INTERIOR_PRUNE_M from any pav_union/runway
                # boundary.
                try:
                    d = on_pav_boundary_obj.distance(Point(v))
                except _GEOM_EXC:
                    d = 0.0
                if d > INTERIOR_PRUNE_M:
                    continue  # drop
                pruned_coords.append(v)
                if pruned_alts is not None:
                    pruned_alts.append(current_alts[i])
            if (len(pruned_coords) >= 3
                    and len(pruned_coords) < len(current_coords)):
                # Validate the pruned polygon before committing.
                try:
                    test_poly = Polygon(pruned_coords).buffer(0)
                except _GEOM_EXC:
                    test_poly = None
                if (test_poly is not None
                        and not test_poly.is_empty
                        and test_poly.geom_type == "Polygon"
                        and test_poly.is_valid
                        and test_poly.is_simple
                        and test_poly.area >= 0.5 * poly.area):
                    current_coords = pruned_coords
                    if pruned_alts is not None:
                        current_alts = pruned_alts

        # Final polygon from the running coords (already validated
        # piecewise; guaranteed to be a single valid Polygon).
        try:
            new_poly = Polygon(current_coords).buffer(0)
        except _GEOM_EXC:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon":
            continue
        shape.polygon = new_poly
        if current_alts is not None:
            shape.node_altitudes = current_alts + [current_alts[0]]


def _enforce_runway_1to1_sharing(layout: PavementLayout) -> None:
    """Rule 1 v4 — surgical runway-edge rewrite (user 2026-05-02).

    For each junction polygon, find each contiguous run of vertices
    within ``RUNWAY_ADJACENCY_TOL_M`` of the runway boundary and
    REPLACE the run with a clean sequence of runway corners ordered
    to match the polygon's walk direction.

    Per-vertex snap rule (user 2026-05-02):
      * Each runway-near vertex's target = nearer endpoint of its
        nearest runway edge.
      * If that target collides with an existing junction vertex
        (snap would shrink), use the other endpoint of the same
        edge instead.

    Order rule:
      * Sort the unique snap targets by their progression from the
        polygon's vertex BEFORE the run to the polygon's vertex
        AFTER the run.
      * Drop targets that would force the polygon to backtrack
        (would create a degenerate spike) — the resulting
        junction-runway interface is then bounded by the polygon
        body's natural extent.

    Per-vertex altitudes drop alongside replaced vertices; new
    runway-corner vertices have no per-vertex altitude (the
    elevation pipeline interpolates from neighbours).
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

    # Build runway segment edges, each paired with its 2 endpoints
    # (= runway corners).  Snap targets are always one of these
    # endpoints — never a runway vertex from a different segment.
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

    adjacency_tol = RUNWAY_ADJACENCY_TOL_M
    vertex_tol = SHARED_VERTEX_TOL_M

    # Per user 2026-05-04: collect every sloping-rect corner.  A
    # junction vertex that already coincides with a rect corner must
    # NOT be replaced by a runway corner — the rect corner is a
    # legitimate 1:1 share that the rect-corner snap pass already
    # established.  Without this guard, the runway-snap run would
    # absorb V3-stub corner -84 (which sat 20 m from the runway) and
    # the polygon edge from V3 corner -85 to the new runway corner
    # would cut across V3, producing a 1012 m² overlap.
    rect_corner_buckets: set = set()
    bucket_size = SHARED_VERTEX_TOL_M
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        if s.polygon is None or s.polygon.is_empty \
                or s.polygon.geom_type != "Polygon":
            continue
        rc = list(s.polygon.exterior.coords)
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        for cx, cy in rc:
            rect_corner_buckets.add(
                (round(cx / bucket_size), round(cy / bucket_size)))
    # A snap "would shrink" when the candidate target is within
    # this distance of any OTHER existing junction vertex
    # (collision = unintended dedupe).
    shrink_collision_tol = SHARED_VERTEX_TOL_M * 2.0

    # Global "claimed runway corners" set (Rule 1 v5, user 2026-05-02):
    # accumulates runway corner positions that previously-processed
    # junctions have snapped to.  Used by the cross-junction shrink
    # check — if a corner is claimed, the next junction's snap to
    # that same corner is allowed (sharing) but if that snap would
    # collide with a NON-runway vertex of either junction, fall back
    # to the other endpoint of the same edge.  Also: PROCESSING
    # ORDER matters; we sort junctions by their nearest runway
    # boundary distance ASCENDING so junctions touching the runway
    # most directly snap first.
    claimed_corners: List[Tuple[float, float]] = []

    # Order junctions by min distance to runway boundary (closest
    # first) so confident snaps commit before borderline cases.
    def _min_d_to_runway(s):
        if s.role != ROLE_JUNCTION or s.polygon is None:
            return float("inf")
        c = list(s.polygon.exterior.coords)
        if c and c[0] == c[-1]:
            c = c[:-1]
        if not c:
            return float("inf")
        m = float("inf")
        for vx, vy in c:
            for ax, ay, bx, by, _, _ in rwy_segs:
                d, _, _ = _point_segment_distance(vx, vy, ax, ay, bx, by)
                if d < m:
                    m = d
                    if m <= 1e-6:
                        return m
        return m

    junction_indices_ordered = sorted(
        (i for i, s in enumerate(layout.shapes)
         if s.role == ROLE_JUNCTION and s.polygon is not None
         and not s.polygon.is_empty
         and s.polygon.geom_type == "Polygon"),
        key=lambda i: _min_d_to_runway(layout.shapes[i]))

    for shape_idx in junction_indices_ordered:
        shape = layout.shapes[shape_idx]
        poly = shape.polygon
        coords = list(poly.exterior.coords)
        node_alts = shape.node_altitudes
        if node_alts is not None and len(node_alts) != len(coords):
            node_alts = None
        had_close = (coords[0] == coords[-1])
        if had_close:
            coords = coords[:-1]
            if node_alts is not None:
                node_alts = list(node_alts[:-1])
        n = len(coords)
        if n < 3:
            continue

        # Pass 1: classify each vertex (runway-adjacent? closest segment?).
        # A vertex already at a sloping-rect corner stays put — the
        # rect-corner share is a legitimate anchor.
        nearest_seg: List[
            Optional[Tuple[Tuple[float, float], Tuple[float, float]]]
        ] = [None] * n
        for i, (vx, vy) in enumerate(coords):
            bk = (round(vx / bucket_size), round(vy / bucket_size))
            if bk in rect_corner_buckets:
                continue
            best_d = adjacency_tol
            best_endpoints = None
            for ax, ay, bx, by, c1, c2 in rwy_segs:
                d, _, _ = _point_segment_distance(vx, vy, ax, ay, bx, by)
                if d < best_d:
                    best_d = d
                    best_endpoints = (c1, c2)
                    if best_d <= 1e-6:
                        break
            if best_endpoints is not None:
                nearest_seg[i] = best_endpoints

        if not any(s is not None for s in nearest_seg):
            continue

        # Pass 2: find contiguous runs of runway-adjacent vertices
        # in circular order.
        runs = _find_circular_runs(
            [s is not None for s in nearest_seg], n)
        if not runs:
            continue

        # Pass 3: for each run, surgically replace it with the
        # ordered runway-corner sequence (Option B).  Pass the
        # GLOBAL claimed-corners list so we don't drop a target
        # already used by an earlier-processed junction (sharing
        # is allowed; the polygons just touch at that vertex).
        new_polygon = _rewrite_runway_runs(
            coords, node_alts, runs, nearest_seg,
            shrink_collision_tol,
            claimed_corners=claimed_corners)
        if new_polygon is None:
            continue
        new_pts, new_alts_out = new_polygon
        if len(new_pts) < 3:
            continue
        try:
            new_poly = Polygon(new_pts).buffer(0)
        except _GEOM_EXC:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon":
            continue
        if not new_poly.is_valid or not new_poly.is_simple:
            continue
        if new_poly.area < 0.5 * poly.area:
            continue
        shape.polygon = new_poly
        if new_alts_out is not None:
            shape.node_altitudes = new_alts_out + [new_alts_out[0]]
        for v in new_pts:
            on_rwy = False
            for ax, ay, bx, by, _, _ in rwy_segs:
                d, _, _ = _point_segment_distance(
                    v[0], v[1], ax, ay, bx, by)
                if d <= vertex_tol:
                    on_rwy = True
                    break
            if on_rwy:
                claimed_corners.append(v)

    # Phase 2 (v6 widening, DISABLED): widening via insertion at
    # outboard runway corners over-grows because the segmented
    # runway has dense corners and insertions cascade across multiple
    # passes.  Needs cleaner per-junction integration with Rule 4
    # split + post-elevation re-segmentation before re-enabling.
    # See user 2026-05-02 thread.


def _rewrite_runway_runs(
    coords: List[Tuple[float, float]],
    node_alts: Optional[List[float]],
    runs: List[List[int]],
    nearest_seg: List[
        Optional[Tuple[Tuple[float, float], Tuple[float, float]]]],
    shrink_collision_tol: float,
    claimed_corners: Optional[List[Tuple[float, float]]] = None,
) -> Optional[Tuple[List[Tuple[float, float]], Optional[List[float]]]]:
    """Surgically replace each runway-near vertex run with the
    ordered sequence of runway corners.  Returns the new (open)
    coord list and (optional) per-vertex altitude list, or None if
    no rewrite is possible.

    Per Option B: the run vertices are DROPPED entirely; runway
    corners are INSERTED in walk order (BEFORE_RUN → AFTER_RUN)
    such that the polygon's body vertices outside the run are
    preserved exactly.
    """
    n = len(coords)
    if not runs:
        return None

    # Build per-run replacement.  Returns list of (run_set,
    # replacement_corners, replacement_alts).
    replacements: List[
        Tuple[set, List[Tuple[float, float]],
              Optional[List[float]]]
    ] = []
    for run_indices in runs:
        # 1) Compute snap target per vertex (shrink-fallback).
        seen_targets: List[Tuple[float, float]] = []
        for idx in run_indices:
            seg_endpoints = nearest_seg[idx]
            if seg_endpoints is None:
                continue
            c1, c2 = seg_endpoints
            vx, vy = coords[idx]
            d1c = math.hypot(vx - c1[0], vy - c1[1])
            d2c = math.hypot(vx - c2[0], vy - c2[1])
            primary, alt = (c1, c2) if d1c <= d2c else (c2, c1)
            # Check collision with junction vertices NOT in this run.
            def _collides(target, excluded_set=set(run_indices)):
                for k in range(n):
                    if k in excluded_set:
                        continue
                    px, py = coords[k]
                    if math.hypot(target[0] - px,
                                  target[1] - py) <= shrink_collision_tol:
                        return True
                return False
            if _collides(primary):
                if not _collides(alt):
                    target = alt
                else:
                    continue
            else:
                target = primary
            seen_targets.append(target)
        # 2) Dedupe targets (preserve order of first appearance).
        unique: List[Tuple[float, float]] = []
        seen_keys: set = set()
        for t in seen_targets:
            key = (round(t[0] * 2.0), round(t[1] * 2.0))  # ~0.5 m bucket
            if key in seen_keys:
                continue
            seen_keys.add(key)
            unique.append(t)
        if not unique:
            continue
        # 3) Determine walk direction: vector from BEFORE_RUN vertex
        # to AFTER_RUN vertex.
        before_idx = (run_indices[0] - 1) % n
        after_idx = (run_indices[-1] + 1) % n
        bvx, bvy = coords[before_idx]
        avx, avy = coords[after_idx]
        dirx = avx - bvx
        diry = avy - bvy
        dirlen = math.hypot(dirx, diry)
        if dirlen < 1e-6:
            continue
        # 4) Order targets by progression along walk direction
        # (project onto dir vector starting at BEFORE_RUN vertex).
        def _proj(t):
            return ((t[0] - bvx) * dirx + (t[1] - bvy) * diry) / dirlen
        unique.sort(key=_proj)
        # 5) Keep all unique targets in walk-projection order.  The
        # validity guard at the caller will reject any rewrite that
        # produces an invalid polygon, so we no longer drop
        # backtracking targets here — sharing claimed corners with
        # adjacent junctions is per user 2026-05-02 explicitly OK.
        kept: List[Tuple[float, float]] = list(unique)
        if not kept:
            continue
        # 6) Compose replacement sequence; per-vertex altitudes are
        # inherited from the run's first/last existing entries when
        # available (just to avoid Nones — elevation pipeline will
        # re-interpolate).
        rep_alts: Optional[List[float]] = None
        if node_alts is not None:
            run_alt_avg = (sum(node_alts[i] for i in run_indices)
                           / max(1, len(run_indices)))
            rep_alts = [run_alt_avg] * len(kept)
        replacements.append((set(run_indices), kept, rep_alts))

    if not replacements:
        return None

    # Build new coord list: walk original polygon, dropping in-run
    # vertices and inserting replacement at the END of each run.
    out_pts: List[Tuple[float, float]] = []
    out_alts: Optional[List[float]] = (
        [] if node_alts is not None else None)
    # Index runs by their first index for quick lookup.
    run_by_first: dict = {r[0][0]: r for r in
                          [(run_indices, rep, alts)
                           for run_indices_set, rep, alts in replacements
                           for run_indices in [sorted(run_indices_set)]]}
    # Easier: build a flat per-index drop set + per-first-index insert.
    drop = set()
    insert_at = {}
    for run_set, rep, alts in replacements:
        drop.update(run_set)
        run_sorted = sorted(run_set)
        insert_at[run_sorted[0]] = (rep, alts)
    for i in range(n):
        if i in drop:
            if i in insert_at:
                rep, alts = insert_at[i]
                for k, t in enumerate(rep):
                    out_pts.append(t)
                    if out_alts is not None:
                        out_alts.append(alts[k])
            # Skip this vertex (in-run, dropped).
            continue
        out_pts.append(coords[i])
        if out_alts is not None:
            out_alts.append(node_alts[i])
    return out_pts, out_alts


def _find_circular_runs(flags: Sequence[bool], n: int) -> List[List[int]]:
    """Find contiguous runs of True values in a circular list of
    length n.  Returns each run as a list of indices in walk order.
    Handles wrap-around (a run that crosses the seam between
    index n-1 and index 0)."""
    if n == 0 or not any(flags):
        return []
    if all(flags):
        return [list(range(n))]
    # Find a transition point (False before True).
    start = 0
    for i in range(n):
        if (not flags[i]) and flags[(i + 1) % n]:
            start = (i + 1) % n
            break
    # Walk from start, collecting runs.
    runs: List[List[int]] = []
    cur: List[int] = []
    for offset in range(n):
        idx = (start + offset) % n
        if flags[idx]:
            cur.append(idx)
        else:
            if cur:
                runs.append(cur)
                cur = []
    if cur:
        runs.append(cur)
    return runs


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
    except _GEOM_EXC:
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
        except _GEOM_EXC:
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


# Junction-vertex-outside-pavement thresholds.  Used by the
# regression-tracking test in ``tests/test_junction_rules.py``
# (``test_junction_vertices_outside_pavement``).  No active
# enforcer anymore — the upstream geometry pipeline is the source
# of truth; the test tracks per-airport drift against a baseline.
PAVEMENT_OUTWARD_OFFSET_M = 0.5
PAVEMENT_INSIDE_TOL_M = 0.1


STITCH_PAVEMENT_ROLES = (
    ROLE_JUNCTION,
    ROLE_APRON,
    "primary_parallel",
    "secondary_parallel",
    "stub",
    "cross_connector",
)


def stitch_pavement_to_terminals(
    layout: PavementLayout,
    snap_corner_m: float = 5.0,
    on_edge_tol_m: float = SHARED_VERTEX_TOL_M,
) -> None:
    """Make terminal pads share an identical vertex set with adjacent
    pavement on every shared boundary segment (user 2026-05-04).

    For each pavement vertex (junction / apron / parallel / stub /
    cross_connector) that lies within ``on_edge_tol_m`` of a terminal
    edge interior:

      * If the vertex is within ``snap_corner_m`` of one of that
        edge's endpoints (a terminal corner), rewrite the pavement
        polygon to use the corner instead — the pavement loses a
        vertex and gains exact alignment with the existing terminal
        corner.
      * Otherwise insert the vertex into the terminal polygon at the
        correct position along the edge — the terminal grows a
        vertex so it matches the pavement node.

    Either way both polygons end up with the same vertex sequence on
    the shared segment, so X-Plane renders a seamless meld with no
    sub-metre overlaps.

    ``node_altitudes`` is updated in lockstep with polygon rewrites.
    Run as the last geometry pass before OSM emit.
    """
    terminals = [s for s in layout.shapes if s.role == ROLE_TERMINAL]
    if not terminals:
        return
    pavements = [s for s in layout.shapes
                 if s.role in STITCH_PAVEMENT_ROLES]
    if not pavements:
        return

    # Per-terminal: list of (a_idx, b_idx, ax, ay, bx, by) edges.
    # ``a_idx`` and ``b_idx`` are positions in the terminal's open
    # ring (no closing-duplicate vertex).
    snap_tol2 = snap_corner_m * snap_corner_m
    on_edge_tol2 = on_edge_tol_m * on_edge_tol_m

    # Inserts collected per terminal: edge_idx → list of
    # (frac_along_edge, x, y).  Applied after the pavement-vertex
    # walk so we don't disturb the terminal geometry mid-iteration.
    pending_inserts: dict = {id(t): {} for t in terminals}

    for pav in pavements:
        poly = pav.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        try:
            coords = list(poly.exterior.coords)
        except _GEOM_EXC:
            continue
        ring_closed = (
            len(coords) > 1 and coords[0] == coords[-1])
        coords_open = coords[:-1] if ring_closed else list(coords)
        if len(coords_open) < 3:
            continue

        node_alts = pav.node_altitudes
        # Mirror the closed/open form for altitudes if present.
        alts_open: Optional[List[float]] = None
        if node_alts is not None:
            alts = list(node_alts)
            alts_open = alts[:-1] if (
                ring_closed and len(alts) == len(coords)) else alts
            if len(alts_open) != len(coords_open):
                alts_open = None

        new_coords: List[Tuple[float, float]] = []
        new_alts: Optional[List[float]] = (
            [] if alts_open is not None else None)
        mutated = False
        for vi, (vx, vy) in enumerate(coords_open):
            snapped = False
            for term in terminals:
                tcoords = list(term.polygon.exterior.coords)
                if tcoords and tcoords[0] == tcoords[-1]:
                    tcoords = tcoords[:-1]
                m = len(tcoords)
                if m < 3:
                    continue
                # Skip if vertex matches an existing terminal corner —
                # already shared, no work needed.
                already_corner = False
                for (cx, cy) in tcoords:
                    if (vx - cx) * (vx - cx) \
                            + (vy - cy) * (vy - cy) <= on_edge_tol2:
                        already_corner = True
                        break
                if already_corner:
                    break
                for ei in range(m):
                    ax, ay = tcoords[ei]
                    bx, by = tcoords[(ei + 1) % m]
                    dx = bx - ax
                    dy = by - ay
                    seg2 = dx * dx + dy * dy
                    if seg2 < 1.0:
                        continue
                    t = ((vx - ax) * dx + (vy - ay) * dy) / seg2
                    if t <= 0.0 or t >= 1.0:
                        continue
                    cx = ax + t * dx
                    cy = ay + t * dy
                    d2 = (vx - cx) * (vx - cx) + (vy - cy) * (vy - cy)
                    if d2 > on_edge_tol2:
                        continue
                    # Vertex sits on this terminal edge interior.
                    # Decide: snap to nearest endpoint, or insert.
                    da2 = (vx - ax) * (vx - ax) + (vy - ay) * (vy - ay)
                    db2 = (vx - bx) * (vx - bx) + (vy - by) * (vy - by)
                    if da2 <= snap_tol2 and da2 <= db2:
                        new_coords.append((ax, ay))
                        if new_alts is not None:
                            new_alts.append(alts_open[vi])
                        mutated = True
                        snapped = True
                    elif db2 <= snap_tol2:
                        new_coords.append((bx, by))
                        if new_alts is not None:
                            new_alts.append(alts_open[vi])
                        mutated = True
                        snapped = True
                    else:
                        # Schedule terminal-side insertion at frac t.
                        pending_inserts[id(term)].setdefault(
                            ei, []).append((t, cx, cy))
                        new_coords.append((cx, cy))
                        if new_alts is not None:
                            new_alts.append(alts_open[vi])
                        # Move pavement vertex onto the EXACT edge
                        # geometry so the bucket-intern in to_osm
                        # collapses both to the same nid.
                        if d2 > 1e-9:
                            mutated = True
                        snapped = True
                    break  # done with this pavement vertex
                if snapped:
                    break
            if not snapped:
                new_coords.append((vx, vy))
                if new_alts is not None:
                    new_alts.append(alts_open[vi])

        if not mutated:
            continue
        # Drop consecutive duplicates introduced by snap-to-corner.
        deduped: List[Tuple[float, float]] = []
        deduped_alts: Optional[List[float]] = (
            [] if new_alts is not None else None)
        for k, p in enumerate(new_coords):
            if (deduped
                    and abs(deduped[-1][0] - p[0]) < 1e-6
                    and abs(deduped[-1][1] - p[1]) < 1e-6):
                continue
            deduped.append(p)
            if deduped_alts is not None:
                deduped_alts.append(new_alts[k])
        if len(deduped) < 3:
            continue
        try:
            new_poly = Polygon(deduped + [deduped[0]])
            if not new_poly.is_valid:
                new_poly = new_poly.buffer(0)
            if (new_poly.is_empty
                    or new_poly.geom_type != "Polygon"):
                continue
        except _GEOM_EXC:
            continue
        pav.polygon = new_poly
        if deduped_alts is not None:
            pav.node_altitudes = deduped_alts + [deduped_alts[0]]

    # Apply terminal-side inserts.
    for term in terminals:
        inserts = pending_inserts.get(id(term))
        if not inserts:
            continue
        tcoords = list(term.polygon.exterior.coords)
        ring_closed = (
            len(tcoords) > 1 and tcoords[0] == tcoords[-1])
        tcoords_open = tcoords[:-1] if ring_closed else list(tcoords)
        m = len(tcoords_open)
        out: List[Tuple[float, float]] = []
        for ei in range(m):
            out.append(tcoords_open[ei])
            if ei in inserts:
                # Sort by t ascending; dedup near-duplicates.
                pts = sorted(inserts[ei], key=lambda x: x[0])
                last_t: float = -1.0
                for t, cx, cy in pts:
                    if t - last_t < 1e-4:
                        continue
                    out.append((cx, cy))
                    last_t = t
        if len(out) < 3:
            continue
        try:
            new_poly = Polygon(out + [out[0]])
            if not new_poly.is_valid:
                new_poly = new_poly.buffer(0)
            if (new_poly.is_empty
                    or new_poly.geom_type != "Polygon"):
                continue
        except _GEOM_EXC:
            continue
        term.polygon = new_poly
        # Terminals carry a single ``s.altitude`` (uniform plane); no
        # ``node_altitudes`` to update.


def stitch_pavement_to_flat_runways(
    layout: PavementLayout,
    corner_tol_m: float = SHARED_VERTEX_TOL_M,
    near_edge_snap_m: float = 2.5,
) -> None:
    """Snap pav vertices that sit within ``near_edge_snap_m`` of a
    flat runway edge onto the runway boundary, and schedule a
    matching runway-side insert so the two polygons share a bucket
    vertex at every snap.  Only flat runway shapes (``role=runway``,
    single ``altitude=`` tag, no ``altitude_high``/``altitude_low``)
    are targeted — sloped 4-corner segments must keep their
    canonical [high, low, low, high] vertex order for OSM emit.

    Run BEFORE the per-surface solver so the new shared vertices
    feed HARD anchors.

    Per user 2026-05-12: a prior "Phase B" projected EVERY
    non-coincident pav vertex (including ones far from the runway
    boundary) perpendicular onto the shared edge, creating
    arbitrary corners on flat runway segments at interior-vertex
    projections.  Removed.  Junction interiors are now anchored
    only via the per-surface Jacobi solver's standard propagation
    from real shared corners (segment seams + Phase A snaps); the
    grade test showed no measurable elevation impact when this
    projection pass was disabled.
    """
    flat_runways = [s for s in layout.shapes
                    if s.role == ROLE_RUNWAY
                    and s.polygon is not None
                    and not s.polygon.is_empty
                    and s.polygon.geom_type == "Polygon"
                    and s.altitude is not None
                    and s.altitude_high is None
                    and s.altitude_low is None]
    if not flat_runways:
        return
    pavements = [s for s in layout.shapes
                 if s.role in STITCH_PAVEMENT_ROLES]
    if not pavements:
        return

    corner_tol2 = corner_tol_m * corner_tol_m
    near_snap2 = near_edge_snap_m * near_edge_snap_m

    # rwy id → ei → list of (t_along_edge, x, y) inserts
    rwy_inserts: dict = {id(r): {} for r in flat_runways}

    def _open_ring(coords):
        if coords and coords[0] == coords[-1]:
            return list(coords[:-1])
        return list(coords)

    # ── Phase A: near-edge snap ──────────────────────────────────
    # For each pav vertex within ``near_edge_snap_m`` of a flat
    # runway edge interior, snap it onto the projection AND
    # schedule a runway-side insert at the same point.  Captures
    # cases where a junction's boundary runs nearly-parallel to a
    # runway edge with a small geometric gap (e.g. SPLP -10053
    # vertices at across=24.8 m next to runway edge at across=22.8
    # m).  After snap, the pav and runway share the bucket; the
    # solver HARD-anchors the pav vertex at the runway altitude.
    for pav in pavements:
        poly = pav.polygon
        if poly is None or poly.is_empty or poly.geom_type != "Polygon":
            continue
        try:
            pav_coords = _open_ring(list(poly.exterior.coords))
        except _GEOM_EXC:
            continue
        n_pav = len(pav_coords)
        if n_pav < 3:
            continue
        node_alts = pav.node_altitudes
        ring_closed_alts = (
            node_alts is not None
            and len(node_alts) == n_pav + 1
            and node_alts[0] == node_alts[-1])
        pav_alts: Optional[List[float]] = None
        if node_alts is not None and (
                len(node_alts) == n_pav or ring_closed_alts):
            pav_alts = list(node_alts[:-1] if ring_closed_alts
                             else node_alts)
        snapped = [False] * n_pav
        new_coords: List[Tuple[float, float]] = list(pav_coords)
        for vi, (vx, vy) in enumerate(pav_coords):
            best_d2 = near_snap2
            best: Optional[Tuple[object, int, float, float, float]] = None
            for rwy in flat_runways:
                try:
                    rcoords = _open_ring(list(
                        rwy.polygon.exterior.coords))
                except _GEOM_EXC:
                    continue
                m = len(rcoords)
                if m < 3:
                    continue
                # Skip if vertex is already at a runway corner.
                already_corner = False
                for (cx, cy) in rcoords:
                    if (vx - cx) * (vx - cx) + (vy - cy) * (vy - cy) \
                            <= corner_tol2:
                        already_corner = True
                        break
                if already_corner:
                    break
                for ei in range(m):
                    ax, ay = rcoords[ei]
                    bx, by = rcoords[(ei + 1) % m]
                    dx = bx - ax
                    dy = by - ay
                    seg2 = dx * dx + dy * dy
                    if seg2 < 1.0:
                        continue
                    t = ((vx - ax) * dx + (vy - ay) * dy) / seg2
                    if t <= 0.001 or t >= 0.999:
                        continue
                    cx = ax + t * dx
                    cy = ay + t * dy
                    d2 = (vx - cx) * (vx - cx) + (vy - cy) * (vy - cy)
                    if d2 >= best_d2:
                        continue
                    best_d2 = d2
                    best = (rwy, ei, t, cx, cy)
            if best is not None:
                rwy, ei, t, cx, cy = best
                new_coords[vi] = (cx, cy)
                snapped[vi] = True
                # Per user 2026-05-11: dedup snap targets within 3 m
                # of each other on the same runway edge.  Phase A
                # processes each pav vertex independently; when
                # multiple pav vertices project to nearly-identical
                # points on the runway edge (apron-to-runway throats
                # with several row-110 vertices in a 3 m span), the
                # runway ends up with a cluster of near-duplicate
                # corners that forces the adjacent junction to taper
                # to a needle.  Keep only the first snap target per
                # 3 m cluster; subsequent pav vertices in the cluster
                # snap to the same point so the JUNCTION has merged
                # vertices too.
                existing = rwy_inserts[id(rwy)].setdefault(ei, [])
                duplicate = False
                for et, ex, ey in existing:
                    if (cx - ex) ** 2 + (cy - ey) ** 2 < 9.0:  # 3 m
                        new_coords[vi] = (ex, ey)
                        duplicate = True
                        break
                if not duplicate:
                    existing.append((t, cx, cy))
        if any(snapped):
            try:
                new_poly = Polygon(new_coords + [new_coords[0]])
                if not new_poly.is_valid:
                    new_poly = new_poly.buffer(0)
                if (not new_poly.is_empty
                        and new_poly.geom_type == "Polygon"):
                    pav.polygon = new_poly
                    if pav_alts is not None:
                        pav.node_altitudes = pav_alts + [pav_alts[0]]
            except _GEOM_EXC:
                pass

    # Apply runway-side inserts.
    for rwy in flat_runways:
        inserts = rwy_inserts.get(id(rwy))
        if not inserts:
            continue
        try:
            rcoords_open = _open_ring(list(
                rwy.polygon.exterior.coords))
        except _GEOM_EXC:
            continue
        m = len(rcoords_open)
        if m < 3:
            continue
        # Dedup inserts that land within DEDUP_INSERT_M of each other
        # on the same runway edge.  Phase A snaps each near-runway
        # pav vertex INDIVIDUALLY; when multiple junction vertices
        # project to nearly-identical points on the runway edge
        # (common at apron-to-runway throats where several apt.dat
        # row-110 vertices sit within a 3 m span), the runway ends
        # up with a cluster of 2-3 near-duplicate corners on its
        # apron-side edge — which forces the adjacent junction to
        # taper to a needle at that point.  Per user 2026-05-11:
        # collapse near-duplicate insertions so the runway gets ONE
        # corner per chart-level transition.
        DEDUP_INSERT_M = 3.0
        out: List[Tuple[float, float]] = []
        for ei in range(m):
            out.append(rcoords_open[ei])
            if ei in inserts:
                pts = sorted(inserts[ei], key=lambda x: x[0])
                last_xy: Optional[Tuple[float, float]] = None
                for t, cx, cy in pts:
                    if last_xy is not None:
                        ddx = cx - last_xy[0]
                        ddy = cy - last_xy[1]
                        if ddx * ddx + ddy * ddy < DEDUP_INSERT_M ** 2:
                            continue
                    out.append((cx, cy))
                    last_xy = (cx, cy)
        if len(out) < 3:
            continue
        try:
            new_poly = Polygon(out + [out[0]])
            if not new_poly.is_valid:
                new_poly = new_poly.buffer(0)
            if (new_poly.is_empty
                    or new_poly.geom_type != "Polygon"):
                continue
        except _GEOM_EXC:
            continue
        rwy.polygon = new_poly
        # Flat runway carries a single ``s.altitude`` — uniform
        # plane, no per-vertex node_altitudes to update.


def stitch_pavement_polygons(
    layout: PavementLayout,
    edge_tol_m: float = 0.5,
    snap_corner_m: float = 1.0,
) -> int:
    """Make adjacent pavement polygons share an identical vertex set
    on every shared boundary segment (user 2026-05-08).

    Two pavement polygons (junction / apron) are "adjacent" when one
    polygon's ring vertex lies within ``edge_tol_m`` of another's
    edge interior.  Without this pass, neighbouring junction
    polygons can have parallel edges that don't share OSM nids; the
    grade validator finds 6 m elevation gaps at mid-edge samples
    along those edges (SPLP junction-10053 ↔ junction-10054).

    For each pair where polygon A's vertex V lies on polygon B's
    edge E:
      * If V is within ``snap_corner_m`` of one of E's endpoints, no
        action — the vertex is effectively at a B corner already
        (and ``_enforce_shared_vertex_altitudes`` will average the
        two polygons' z at that bucket).
      * Otherwise insert a new vertex at V's exact (x, y) into B's
        ring with z linearly interpolated from E's endpoints.  After
        insertion, V and the new vertex are at the same XY → same
        OSM-emit bucket → ``_enforce_shared_vertex_altitudes``
        averages their elevations so the rendered surface meets at
        a single z at every shared point.

    Companion to ``stitch_pavement_to_terminals``; runs before the
    final altitude reconciliation chain at the end of
    ``build_airport_pavement``.

    Returns the total number of vertices inserted across all
    polygons.
    """
    PAVEMENT_LIKE = (ROLE_JUNCTION, ROLE_APRON)
    pavements = [
        (i, s) for i, s in enumerate(layout.shapes)
        if s.role in PAVEMENT_LIKE
        and s.polygon is not None
        and not s.polygon.is_empty
        and s.polygon.geom_type == "Polygon"
        and s.node_altitudes
    ]
    if len(pavements) < 2:
        return 0

    edge_tol2 = edge_tol_m * edge_tol_m
    snap_tol2 = snap_corner_m * snap_corner_m

    # Cache per-shape open-ring + altitude views.
    cache: Dict[int, Tuple[
        List[Tuple[float, float]], List[float]]] = {}
    for idx, s in pavements:
        coords = list(s.polygon.exterior.coords)
        ring_closed = (
            len(coords) > 1 and coords[0] == coords[-1])
        coords_open = coords[:-1] if ring_closed else list(coords)
        alts = list(s.node_altitudes)
        alts_open = (alts[:-1] if (
            ring_closed and len(alts) == len(coords)) else alts)
        if len(alts_open) != len(coords_open):
            continue
        cache[idx] = (coords_open, alts_open)

    # Pending inserts per polygon: {b_idx: {edge_idx: [(t, x, y, z), ...]}}
    pending: Dict[int, Dict[int, List[
        Tuple[float, float, float, float]]]] = {}

    for a_idx, A in pavements:
        if a_idx not in cache:
            continue
        a_coords, _ = cache[a_idx]
        for vi, (vx, vy) in enumerate(a_coords):
            for b_idx, B in pavements:
                if b_idx == a_idx or b_idx not in cache:
                    continue
                b_coords, b_alts = cache[b_idx]
                m = len(b_coords)
                if m < 3:
                    continue
                # Skip if V is already at any B corner.
                already_corner = False
                for (cx, cy) in b_coords:
                    if (vx - cx) * (vx - cx) \
                            + (vy - cy) * (vy - cy) <= snap_tol2:
                        already_corner = True
                        break
                if already_corner:
                    continue
                # Find B edge that V projects onto within edge_tol.
                for ei in range(m):
                    ax, ay = b_coords[ei]
                    bx, by = b_coords[(ei + 1) % m]
                    dx = bx - ax
                    dy = by - ay
                    seg2 = dx * dx + dy * dy
                    if seg2 < 1.0:
                        continue
                    t = ((vx - ax) * dx + (vy - ay) * dy) / seg2
                    if t <= 0.001 or t >= 0.999:
                        continue
                    cx = ax + t * dx
                    cy = ay + t * dy
                    d2 = (vx - cx) * (vx - cx) \
                        + (vy - cy) * (vy - cy)
                    if d2 > edge_tol2:
                        continue
                    # Found.  Compute interpolated z along B's edge
                    # at fraction t.
                    z_a = b_alts[ei]
                    z_b = b_alts[(ei + 1) % m]
                    interp_z = z_a * (1.0 - t) + z_b * t
                    # Use V's exact (x, y) so A's vertex and B's new
                    # vertex hash to the same OSM-emit bucket.
                    pending.setdefault(b_idx, {}).setdefault(
                        ei, []).append((t, vx, vy, interp_z))
                    break  # done with this A vertex

    if not pending:
        return 0

    n_inserts = 0
    for b_idx, edge_inserts in pending.items():
        if b_idx not in cache:
            continue
        B = layout.shapes[b_idx]
        b_coords, b_alts = cache[b_idx]
        m = len(b_coords)
        new_coords: List[Tuple[float, float]] = []
        new_alts: List[float] = []
        for ei in range(m):
            new_coords.append(b_coords[ei])
            new_alts.append(b_alts[ei])
            if ei in edge_inserts:
                pts = sorted(edge_inserts[ei], key=lambda x: x[0])
                last_t = -1.0
                for t, ix, iy, iz in pts:
                    if t - last_t < 1e-3:
                        continue
                    new_coords.append((ix, iy))
                    new_alts.append(iz)
                    n_inserts += 1
                    last_t = t
        if len(new_coords) < 3:
            continue
        try:
            new_poly = Polygon(new_coords + [new_coords[0]])
            if not new_poly.is_valid:
                new_poly = new_poly.buffer(0)
            if (new_poly.is_empty
                    or new_poly.geom_type != "Polygon"):
                continue
        except _GEOM_EXC:
            continue
        B.polygon = new_poly
        B.node_altitudes = new_alts + [new_alts[0]]
    return n_inserts


def _vertex_on_any_anchor_edge(
    vx: float, vy: float,
    anchor_edges: Sequence[Tuple[float, float, float, float]],
    tol: float,
) -> bool:
    tol2 = tol * tol
    for ax, ay, bx, by in anchor_edges:
        dx = bx - ax
        dy = by - ay
        seg2 = dx * dx + dy * dy
        if seg2 < 1e-9:
            continue
        t = ((vx - ax) * dx + (vy - ay) * dy) / seg2
        if t < 0.0:
            t = 0.0
        elif t > 1.0:
            t = 1.0
        cx = ax + t * dx
        cy = ay + t * dy
        d2 = (vx - cx) * (vx - cx) + (vy - cy) * (vy - cy)
        if d2 <= tol2:
            return True
    return False
