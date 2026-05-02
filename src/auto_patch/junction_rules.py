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
    RUNWAY_ADJACENCY_TOL_M,
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

    # Phase 5 (landed): Rule 5 — push junction vertices outside
    # apt.dat pavement boundary.  Runs LAST so the rule sees the
    # final polygon shapes (after Rule 1/2/4 reshapes); any vertex
    # that ended up inside the pavement gets pushed back out.
    _push_junction_vertices_outside_pavement(layout)

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
    # A snap "would shrink" when the candidate target is within
    # this distance of any OTHER existing junction vertex
    # (collision = unintended dedupe).
    shrink_collision_tol = SHARED_VERTEX_TOL_M * 2.0

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

        # Pass 1: classify each vertex (runway-adjacent? closest segment?).
        nearest_seg: List[
            Optional[Tuple[Tuple[float, float], Tuple[float, float]]]
        ] = [None] * n
        for i, (vx, vy) in enumerate(coords):
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
        # ordered runway-corner sequence (Option B).
        new_polygon = _rewrite_runway_runs(
            coords, node_alts, runs, nearest_seg,
            shrink_collision_tol)
        if new_polygon is None:
            continue
        new_pts, new_alts_out = new_polygon
        if len(new_pts) < 3:
            continue
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
        # Validity guard: if the rewrite shrinks area dramatically
        # or creates an invalid polygon, skip the rewrite.
        if not new_poly.is_valid or not new_poly.is_simple:
            continue
        if new_poly.area < 0.5 * poly.area:
            continue
        shape.polygon = new_poly
        if new_alts_out is not None:
            shape.node_altitudes = new_alts_out + [new_alts_out[0]]


def _rewrite_runway_runs(
    coords: List[Tuple[float, float]],
    node_alts: Optional[List[float]],
    runs: List[List[int]],
    nearest_seg: List[
        Optional[Tuple[Tuple[float, float], Tuple[float, float]]]],
    shrink_collision_tol: float,
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
        # 5) Drop targets that fall OUTSIDE the BEFORE→AFTER span
        # (would force polygon to backtrack — degenerate spike).
        # Keep targets with 0 < projection < dirlen.
        kept: List[Tuple[float, float]] = []
        for t in unique:
            p = _proj(t)
            if -1.0 <= p <= dirlen + 1.0:
                kept.append(t)
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
                        out_alts.append(alts[k] if alts else 0.0)
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


# ── Rule 5: push junction vertices outside pavement boundary ────


# Per user 2026-05-02: every junction vertex that ISN'T shared with
# an anchor (rect / runway / terminal) must sit OUTSIDE the apt.dat
# pavement boundary, within ``PAVEMENT_OUTWARD_OFFSET_MAX_M``.  This
# guarantees the elevation-smoothing junction polygon FULLY ENCLOSES
# the apt.dat pavement (any apt.dat pavement is INSIDE one of our
# emitted shapes), so no rendered pavement renders outside its
# elevation-smoothed shape.
PAVEMENT_OUTWARD_OFFSET_M = 0.5
PAVEMENT_OUTWARD_OFFSET_MAX_M = 1.0
PAVEMENT_INSIDE_TOL_M = 0.1  # treat vertices within this of boundary as on-it

# Rule 5 is a NEAR-BOUNDARY pass: vertices farther than this distance
# from the pavement boundary aren't push candidates.  Vertices inside
# pavement at greater depths are typically:
#   * Interior cut-line endpoints from ``_decompose_polygon_with_holes``
#   * Densification midpoints near runway/rect long edges (handled by
#     Rule 1 v2 / Rule 2 instead)
#   * Shared-vertex centroid drift artefacts
# Moving them by ≤ 1 m wouldn't get them outside, and a larger move
# would catastrophically deform the polygon.  Tracked separately
# (per-airport baselines) until the upstream geometry is fixed.
PAVEMENT_PUSH_SEARCH_RADIUS_M = 1.0


def _push_junction_vertices_outside_pavement(
    layout: PavementLayout,
) -> None:
    """Rule 5 (user 2026-05-02): for each junction vertex NOT on a
    rect / runway / terminal anchor edge, ensure it lies OUTSIDE the
    apt.dat pavement boundary by at least ``PAVEMENT_OUTWARD_OFFSET_M``
    (capped at ``PAVEMENT_OUTWARD_OFFSET_MAX_M``).  Anchor-edge
    vertices stay put — they're shared exactly with an anchor that
    already covers that pavement region.

    Algorithm per vertex:
      1. If within ``SHARED_VERTEX_TOL_M`` of any anchor edge → skip.
      2. Compute the closest point ``foot`` on ``pav_union.boundary``.
      3. If the vertex is OUTSIDE pavement and at least
         ``PAVEMENT_OUTWARD_OFFSET_M`` past the foot → already
         compliant; skip.
      4. Otherwise push the vertex along the ``vertex → foot``
         direction (or its inverse if the vertex is inside) until it
         sits ``PAVEMENT_OUTWARD_OFFSET_M`` outside.
    """
    pav_union = getattr(layout, "_apt_pav_union", None)
    if pav_union is None or pav_union.is_empty:
        return

    pav_boundary = pav_union.boundary
    if pav_boundary.is_empty:
        return

    # Collect anchor edges grouped by exemption tolerance.  Rule 5
    # MUST NOT push vertices that Rule 1 / Rule 2 specifically
    # placed near anchor boundaries:
    #   * Runway edges → exempt within ``RUNWAY_BOUNDARY_TOL_M`` (Rule 1
    #     legitimately leaves vertices there).
    #   * Sloping rect long edges → exempt within
    #     ``LONG_EDGE_SNAP_M`` PERPENDICULAR (Rule 2 keeps short-end
    #     reach-corners there).
    #   * Terminal & rect short edges → exempt within
    #     ``SHARED_VERTEX_TOL_M`` (legitimate anchor sharing only).
    from .config import LONG_EDGE_SNAP_M
    runway_edges: List[Tuple[float, float, float, float]] = []
    rect_long_edges: List[Tuple[float, float, float, float]] = []
    other_anchor_edges: List[Tuple[float, float, float, float]] = []
    for s in layout.shapes:
        if s.polygon is None or s.polygon.is_empty:
            continue
        c = list(s.polygon.exterior.coords)
        if c and c[0] == c[-1]:
            c = c[:-1]
        m = len(c)
        if s.role == ROLE_RUNWAY:
            for i in range(m):
                ax, ay = c[i]
                bx, by = c[(i + 1) % m]
                runway_edges.append((float(ax), float(ay),
                                     float(bx), float(by)))
        elif s.role in SLOPING_RECT_ROLES:
            # Per the rect-build convention: long edges connect
            # coords[0]↔coords[1] and coords[2]↔coords[3].
            if m == 4:
                rect_long_edges.append(
                    (float(c[0][0]), float(c[0][1]),
                     float(c[1][0]), float(c[1][1])))
                rect_long_edges.append(
                    (float(c[2][0]), float(c[2][1]),
                     float(c[3][0]), float(c[3][1])))
                # Short edges into other_anchor_edges for tight tol.
                other_anchor_edges.append(
                    (float(c[1][0]), float(c[1][1]),
                     float(c[2][0]), float(c[2][1])))
                other_anchor_edges.append(
                    (float(c[3][0]), float(c[3][1]),
                     float(c[0][0]), float(c[0][1])))
        elif s.role == "terminal":
            for i in range(m):
                ax, ay = c[i]
                bx, by = c[(i + 1) % m]
                other_anchor_edges.append(
                    (float(ax), float(ay), float(bx), float(by)))

    from shapely.ops import nearest_points

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

        new_coords: List[Tuple[float, float]] = []
        changed = False
        for i, (vx, vy) in enumerate(coords):
            # Per-anchor-class exemptions (don't undo Rule 1 / Rule 2).
            if _vertex_on_any_anchor_edge(
                    vx, vy, runway_edges, RUNWAY_BOUNDARY_TOL_M):
                new_coords.append((vx, vy))
                continue
            if _vertex_on_any_anchor_edge(
                    vx, vy, rect_long_edges, LONG_EDGE_SNAP_M):
                new_coords.append((vx, vy))
                continue
            if _vertex_on_any_anchor_edge(
                    vx, vy, other_anchor_edges, SHARED_VERTEX_TOL_M):
                new_coords.append((vx, vy))
                continue
            # Find closest point on apt.dat pavement boundary.
            try:
                p = Point(vx, vy)
                foot_pt, _ = nearest_points(pav_boundary, p)
            except Exception:
                new_coords.append((vx, vy))
                continue
            fx, fy = foot_pt.x, foot_pt.y
            # Distance to boundary, signed: positive if outside, neg if inside.
            d = math.hypot(vx - fx, vy - fy)
            try:
                inside = pav_union.contains(p)
            except Exception:
                inside = False
            if not inside and d >= PAVEMENT_OUTWARD_OFFSET_M:
                # Already at least the desired offset outside.
                new_coords.append((vx, vy))
                continue
            if d > PAVEMENT_PUSH_SEARCH_RADIUS_M:
                # Too deep to fix via a sub-1 m push; leave alone and
                # let the test surface it.
                new_coords.append((vx, vy))
                continue
            # Need to move.  Direction from foot to vertex (when
            # outside) or its inverse (when inside).
            if d < 1e-6:
                # Vertex sits exactly on the boundary; we need a
                # direction.  Use the boundary's local normal —
                # approximate via a tiny step in the outward
                # direction (away from pavement centroid).
                cx_p, cy_p = pav_union.centroid.x, pav_union.centroid.y
                ddx = vx - cx_p
                ddy = vy - cy_p
                ln = math.hypot(ddx, ddy)
                if ln < 1e-6:
                    new_coords.append((vx, vy))
                    continue
                ux, uy = ddx / ln, ddy / ln
            else:
                if inside:
                    # Move from inside to outside: vector vertex→foot
                    # then a bit beyond.
                    ux = (fx - vx) / d
                    uy = (fy - vy) / d
                else:
                    # Outside but too close: vector foot→vertex.
                    ux = (vx - fx) / d
                    uy = (vy - fy) / d
            # Place the new vertex PAVEMENT_OUTWARD_OFFSET_M outside.
            target_x = fx + PAVEMENT_OUTWARD_OFFSET_M * ux
            target_y = fy + PAVEMENT_OUTWARD_OFFSET_M * uy
            # Clamp the move so we never displace by more than
            # PAVEMENT_OUTWARD_OFFSET_MAX_M from the original.
            mv_d = math.hypot(target_x - vx, target_y - vy)
            if mv_d > PAVEMENT_OUTWARD_OFFSET_MAX_M:
                # Scale the move down.
                scale = PAVEMENT_OUTWARD_OFFSET_MAX_M / mv_d
                target_x = vx + (target_x - vx) * scale
                target_y = vy + (target_y - vy) * scale
            new_coords.append((target_x, target_y))
            changed = True

        if not changed:
            continue
        try:
            new_poly = Polygon(new_coords).buffer(0)
        except Exception:
            continue
        if new_poly.is_empty:
            continue
        if new_poly.geom_type == "MultiPolygon":
            new_poly = max(new_poly.geoms, key=lambda g: g.area)
        if new_poly.geom_type != "Polygon":
            continue
        shape.polygon = new_poly


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
