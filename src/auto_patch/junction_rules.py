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
from shapely.ops import unary_union

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
    """Per user 2026-05-02: a sloping rect's altitude must vary
    along its source_axis (centerline) only — never perpendicular.
    Per-rect:
      * Project the 4 corners onto source_axis to identify the two
        short-end groups.
      * Compute the axial altitude difference (high-end avg minus
        low-end avg).
      * If |axial diff| < ``SLOPING_RECT_FLAT_THRESHOLD_M``, convert
        to FLAT: set ``altitude = average``, clear ``altitude_high``
        and ``altitude_low``.
      * Otherwise re-set ``altitude_high`` / ``altitude_low`` to
        the axial extremes (drops any perpendicular slope component).

    Flat rects are exempt from sloping-rect connection rules
    (Rule 2 etc.) and may receive junction connections on any
    side.
    """
    rect_roles = ("primary_parallel", "secondary_parallel",
                  "stub", "cross_connector")
    for s in layout.shapes:
        if s.role not in rect_roles:
            continue
        if s.altitude_high is None or s.altitude_low is None:
            continue  # already flat
        if s.source_axis is None or s.source_axis.is_empty:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) != 4:
            continue
        ax_pts = list(s.source_axis.coords)
        if len(ax_pts) < 2:
            continue
        ax_start = ax_pts[0]
        ax_end = ax_pts[-1]
        axdx = ax_end[0] - ax_start[0]
        axdy = ax_end[1] - ax_start[1]
        ax_len2 = axdx * axdx + axdy * axdy
        if ax_len2 < 1e-9:
            continue
        # Per the legacy convention used by triangulation.py:
        # coords[0] and coords[3] = altitude_high end;
        # coords[1] and coords[2] = altitude_low end.
        corner_alt = [s.altitude_high, s.altitude_low,
                      s.altitude_low, s.altitude_high]
        # Project each corner onto source_axis.
        ts = []
        for c in coords:
            t = ((c[0] - ax_start[0]) * axdx
                 + (c[1] - ax_start[1]) * axdy) / ax_len2
            ts.append(t)
        # Group corners into "near-start" and "near-end" by t.
        order = sorted(range(4), key=lambda i: ts[i])
        start_corners = order[:2]
        end_corners = order[2:]
        start_alt = (corner_alt[start_corners[0]]
                     + corner_alt[start_corners[1]]) / 2.0
        end_alt = (corner_alt[end_corners[0]]
                   + corner_alt[end_corners[1]]) / 2.0
        diff = abs(end_alt - start_alt)
        if diff < SLOPING_RECT_FLAT_THRESHOLD_M:
            # Flat in the source_axis direction → use single altitude.
            s.altitude = (start_alt + end_alt) / 2.0
            s.altitude_high = None
            s.altitude_low = None
        # If diff is significant: leave as-is for now (re-aligning
        # the slope direction to source_axis would require also
        # reordering polygon vertex indices to match the convention
        # used by triangulation.py — non-trivial; deferred).


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

    Return format: ``[(p1, p2, corner_a, corner_b), ...]`` (same
    as the legacy ``_rect_long_edges`` API).
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


# Backward-compat alias.
_rect_long_edges = _rect_sloping_edges


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
        # Per user 2026-05-02: flat rects (no altitude_high/low,
        # only ``altitude``) are exempt from sloping-rect connection
        # rules — junctions can connect anywhere on their boundary.
        if (shape.altitude_high is None
                or shape.altitude_low is None):
            continue
        rect = shape.polygon
        if rect is None or rect.is_empty or rect.geom_type != "Polygon":
            continue
        for ax, ay, bx, by in (
            (e[0][0], e[0][1], e[1][0], e[1][1])
            for e in _rect_sloping_edges(rect, shape.source_axis)
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
    except Exception:
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
        if len(coords) != 4:
            continue
        if (s.altitude_high is not None
                and s.altitude_low is not None):
            corner_alts = (
                (coords[0], float(s.altitude_high)),
                (coords[1], float(s.altitude_low)),
                (coords[2], float(s.altitude_low)),
                (coords[3], float(s.altitude_high)),
            )
        elif s.altitude is not None:
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
    except Exception:
        runway_union = None
    return _do_widen(
        layout, chain, corner_index, corner_alt, runway_union)


def _do_widen(
    layout: PavementLayout,
    chain: Sequence[Tuple[float, float]],
    corner_index: dict,
    corner_alt: dict,
    runway_union,
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
        # corner falling within the joining region).
        max_total_shared = 5
        current_shared_count = len(shared_in_poly)
        max_inserts = max(0, max_total_shared - current_shared_count)
        if max_inserts == 0:
            continue

        # Build insertions: list of (insert_at_position, new_vertex,
        # new_altitude).  insert_at_position uses the ORIGINAL coords
        # indices; we'll apply them in reverse order.
        insertions: List[
            Tuple[int, Tuple[float, float], Optional[float]]] = []

        for poly_idx, corner in shared_in_poly:
            if len(insertions) >= max_inserts:
                break
            ci = corner_index[_key(corner)]
            chain_neighbors = (
                chain[(ci - 1) % n_chain],
                chain[(ci + 1) % n_chain],
            )
            prev_v = coords[(poly_idx - 1) % n]
            next_v = coords[(poly_idx + 1) % n]

            for neighbor in chain_neighbors:
                if len(insertions) >= max_inserts:
                    break
                # Skip if already a junction vertex.
                if _key(neighbor) in existing_keys:
                    continue
                # Decide insertion side: BEFORE poly_idx or AFTER.
                # Compare angles from corner to neighbor vs to
                # prev_v / next_v in polygon walk; pick the side
                # whose angle is closer to the neighbor's.
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

                d_to_prev = _ang_diff(a_n, a_p)
                d_to_next = _ang_diff(a_n, a_x)
                if d_to_prev <= d_to_next:
                    insert_at = poly_idx       # before
                    flank_v = prev_v           # vertex before corner
                else:
                    insert_at = poly_idx + 1   # after
                    flank_v = next_v           # vertex after corner

                # Per user 2026-05-02 issue #1/#2: reject the
                # insertion if it creates a JAG (sharp backtrack)
                # in the polygon walk.  At the inserted vertex N,
                # the turn angle = angle between (flank→N) and
                # (N→corner) for "before" insertion, or
                # (corner→N) and (N→flank) for "after".  A turn >
                # JAG_THRESHOLD_DEG means the polygon doubles back
                # — the corner is too far outboard for the
                # polygon body to reach naturally.
                if d_to_prev <= d_to_next:
                    # Walk: flank_v → N → corner
                    e1 = (neighbor[0] - flank_v[0],
                          neighbor[1] - flank_v[1])
                    e2 = (corner[0] - neighbor[0],
                          corner[1] - neighbor[1])
                else:
                    # Walk: corner → N → flank_v
                    e1 = (neighbor[0] - corner[0],
                          neighbor[1] - corner[1])
                    e2 = (flank_v[0] - neighbor[0],
                          flank_v[1] - neighbor[1])
                m1 = math.hypot(*e1)
                m2 = math.hypot(*e2)
                if m1 > 1e-6 and m2 > 1e-6:
                    cos_turn = (e1[0] * e2[0]
                                + e1[1] * e2[1]) / (m1 * m2)
                    # cos_turn near +1 = straight; near -1 = U-turn.
                    # Reject when turn > 135° (cos < -0.7).
                    if cos_turn < -0.7:
                        continue

                alt = corner_alt.get(_key(neighbor))
                insertions.append((insert_at, neighbor, alt))
                # Track this neighbor as "now-in-polygon" so we
                # don't insert duplicates from other shared corners.
                existing_keys.add(_key(neighbor))

        if not insertions:
            continue

        # Apply insertions in REVERSE position order so earlier
        # indices remain valid.
        insertions.sort(key=lambda x: -x[0])
        new_coords = list(coords)
        new_alts: Optional[List[float]] = (
            list(node_alts) if node_alts is not None else None)
        for pos, vert, alt in insertions:
            new_coords.insert(pos, vert)
            if new_alts is not None:
                new_alts.insert(pos, alt if alt is not None else 0.0)

        # Validate polygon.
        if len(new_coords) < 3:
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
        if not new_poly.is_valid or not new_poly.is_simple:
            continue
        # Allow polygon to GROW (Rule 1 v6 widening) but not shrink
        # significantly — a > 50 % shrink suggests degenerate snap.
        if new_poly.area < 0.5 * poly.area:
            continue
        # Overlap rejection: if the widened polygon overlaps the
        # runway OR any other junction polygon by > a tiny noise
        # tolerance, the insertion wrapped the polygon body across
        # a neighbour boundary — revert.
        if runway_union is not None and not runway_union.is_empty:
            try:
                ovl = new_poly.intersection(runway_union).area
            except Exception:
                ovl = 0.0
            if ovl > 1.0:
                continue
        # Check overlap against OTHER junctions (excluding self).
        bad = False
        for other in layout.shapes:
            if other is shape or other.role != ROLE_JUNCTION:
                continue
            if other.polygon is None or other.polygon.is_empty:
                continue
            try:
                if new_poly.intersection(other.polygon).area > 1.0:
                    bad = True
                    break
            except Exception:
                pass
        if bad:
            continue
        shape.polygon = new_poly
        if new_alts is not None:
            shape.node_altitudes = new_alts + [new_alts[0]]


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
        except Exception:
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
            # Sloping vs flat: use source_axis-based detection
            # (per user 2026-05-02 clarification — what matters
            # is direction of slope, not edge length).
            if m == 4:
                sloping = _rect_sloping_edges(s.polygon, s.source_axis)
                sloping_keys = set()
                for sa, sb, _, _ in sloping:
                    rect_long_edges.append(
                        (float(sa[0]), float(sa[1]),
                         float(sb[0]), float(sb[1])))
                    sloping_keys.add((round(sa[0] * 2.0),
                                      round(sa[1] * 2.0),
                                      round(sb[0] * 2.0),
                                      round(sb[1] * 2.0)))
                # Flat edges = the other 2 (short ends in typical
                # case): tight-tol anchor exemption.
                for i in range(4):
                    sa = c[i]
                    sb = c[(i + 1) % 4]
                    key = (round(sa[0] * 2.0), round(sa[1] * 2.0),
                           round(sb[0] * 2.0), round(sb[1] * 2.0))
                    rkey = (round(sb[0] * 2.0), round(sb[1] * 2.0),
                            round(sa[0] * 2.0), round(sa[1] * 2.0))
                    if key in sloping_keys or rkey in sloping_keys:
                        continue
                    other_anchor_edges.append(
                        (float(sa[0]), float(sa[1]),
                         float(sb[0]), float(sb[1])))
                # Skip the legacy index-based block below.
                continue
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
