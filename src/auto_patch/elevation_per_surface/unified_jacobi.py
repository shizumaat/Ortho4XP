"""Per-surface unified Jacobi elevation solver (user 2026-05-03).

Lifted from ``auto_patch.elevation._solve_pavement_elevations_unified``
(commit ``35db401`` baseline) with three targeted changes that
implement the per-axis grade rule:

1. **Rects (taxi roles) get RING EDGES ONLY — no within-shape spatial
   pairs, plus a cross-section flatness constraint.**  Per user
   2026-05-03 terminology: a rect has two "sloping edges" (parallel
   to ``source_axis``, where slope is allowed at ≤ 1.5 % per metre
   of axial travel) and two "axis-end edges" (perpendicular to
   ``source_axis``, which must be EXACTLY FLAT — zero perpendicular
   delta).  Avoid "short" / "long" — a taxi rect can be wider than
   it is long (e.g. SPJC stub F).  The flatness is enforced as an
   equality constraint group on the two corners at each axis-end
   (``rect_flat_groups``), not a 1.5 %-cap edge.

   No diagonals.  No edges from a rect to a perpendicular runway
   150 m away.  Junction vertices touching a rect must coincide
   with the rect's corners only — no intermediate nodes on the
   sloping edge (handled upstream by junction emission rules).

2. **Junctions / aprons / terminals get RING + ALL-PAIR Euclidean
   spatial edges (no radius cap).**  These are multi-directional
   surfaces — a plane can taxi across in any direction, so the
   role's grade cap applies between any two vertices on the same
   polygon, not just ring-adjacent ones.  The legacy 60 m radius
   cap was the source of the F-stub-vs-runway grade violation
   (user 2026-05-03): a junction wider than 60 m had un-constrained
   vertex pairs that ended up at incompatible elevations.

3. **Terminals are SOFT, not HARD-anchored.**  Only runway corners
   (CIFP profile) are immutable.  Terminals enter the solver as
   soft nodes seeded from DEM-median, with a flatness constraint
   (all corners share one value, set to the iteration-average each
   pass).  The pre-solver "max grade-compliant from runway corners
   within 250 m" terminal pin is intentionally bypassed when this
   solver runs.

Cross-shape continuity uses shape-shared vertex buckets (same node
index in the unified graph), which makes elevation continuity at
shared corners automatic.
"""
from __future__ import annotations

import math
import time as _time
from typing import Dict, List, Optional, Tuple

from shapely.errors import GEOSException, TopologicalError

from auto_patch.elevation import APRON_MAX_GRADE, TAXI_MAX_GRADE
from auto_patch.layout import (
    ROLE_APRON, ROLE_BOUNDARY, ROLE_CROSS_CONNECTOR, ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL, ROLE_RUNWAY, ROLE_RUNWAY_CROSSING,
    ROLE_SECONDARY_PARALLEL, ROLE_STUB, ROLE_TERMINAL,
)

# Narrow exception tuple for shapely / numeric-geometry failure
# modes.  Programming errors propagate so they surface immediately.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)


SLOPING_RECT_ROLES = (
    ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
    ROLE_STUB, ROLE_CROSS_CONNECTOR,
)

PAVEMENT_ROLES = {
    ROLE_RUNWAY, *SLOPING_RECT_ROLES,
    ROLE_APRON, ROLE_TERMINAL, ROLE_JUNCTION,
    # Per user 2026-05-18: runway-crossing junctions carry runway-
    # interpolated ``node_altitudes`` from
    # ``_resolve_runway_crossings``.  Treat them as HARD-anchored
    # ring-only edges (same path as ``ROLE_RUNWAY``) so the solver
    # doesn't reshape elevations that the runway-interpolation
    # already established.
    ROLE_RUNWAY_CROSSING,
}

CAP_SWEEPS_PER_ITER = 5


def _role_grade(role: str) -> float:
    """Per-role max grade cap.  All roles now share ``TAXI_MAX_GRADE``
    (1.5 %, user 2026-05-18): the apron-reclassification pipeline
    pass folds true apron-territory pavement into ``ROLE_APRON``,
    and the cap was relaxed from the FAA-1.0 % parking-surface limit
    to match the taxiway cap so reclassified shapes don't trip the
    solver / audit at every other vertex pair.
    """
    if role in (ROLE_RUNWAY, *SLOPING_RECT_ROLES, ROLE_JUNCTION):
        return TAXI_MAX_GRADE
    return APRON_MAX_GRADE


def _open_ring(coords) -> List[Tuple[float, float]]:
    if coords and coords[0] == coords[-1]:
        return list(coords[:-1])
    return list(coords)


def solve(layout, icao: str,
          max_iters: int = 5000, tol_m: float = 0.001,
          dem=None, tile_lat: int = 0, tile_lon: int = 0) -> None:
    """Run the constrained-Laplacian solver and write elevations
    back onto each shape.  Mutates ``layout`` in place.

    ``dem`` is sampled per-vertex during seeding so SOFT nodes
    start at their natural terrain elevation; cap projection then
    pulls them toward HARD anchors only where the per-edge grade
    cap requires it.  Without DEM, soft nodes seed from existing
    layout values (or backfill from nearest HARD), which collapses
    DEM-elevated terrain to runway level.
    """
    t_start = _time.time()
    nodes, bucket_to_idx = _build_node_list(layout)
    if not nodes:
        return
    n = len(nodes)

    elev, is_hard, _have_initial = _seed_elevations(
        layout, nodes, bucket_to_idx,
        dem=dem, tile_lat=tile_lat, tile_lon=tile_lon)
    if not any(is_hard):
        return

    edge_grade, edge_length = _build_edges(
        layout, bucket_to_idx)
    if not edge_grade:
        return

    adj = _build_adjacency(n, edge_grade, edge_length)
    terminal_groups = _build_terminal_groups(
        layout, bucket_to_idx)
    rect_flat_groups = _build_rect_cross_section_groups(
        layout, bucket_to_idx)
    edge_list = list(edge_grade.keys())

    iters_used = _run_jacobi(
        elev, is_hard, adj, edge_list,
        edge_grade, edge_length, terminal_groups,
        rect_flat_groups, max_iters, tol_m)
    n_terms, n_rects, n_juncs = _writeback(
        layout, elev, bucket_to_idx)
    _report(icao, iters_used, max_iters,
             _time.time() - t_start,
             n_terms, n_rects, n_juncs)


# ── Stage 1: build node list ──────────────────────────────────────


def _build_node_list(layout):
    """Assign one node index per unique canonical point across all
    pavement-role shapes.  Returns ``(nodes, bucket_to_idx)`` —
    the dict still names ``bucket_to_idx`` for legacy continuity
    but keys are canonical (x, y) tuples when the layout has a
    registry, else legacy discrete buckets.
    """
    bucket_to_idx: Dict = {}
    nodes: List[Tuple[float, float]] = []
    for s in layout.shapes:
        if s.role not in PAVEMENT_ROLES:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = _open_ring(list(s.polygon.exterior.coords))
        except _GEOM_EXC:
            continue
        for x, y in coords:
            k = layout.canonical_points.get_or_add(float(x), float(y))
            if k not in bucket_to_idx:
                bucket_to_idx[k] = len(nodes)
                nodes.append((float(x), float(y)))
    return nodes, bucket_to_idx


# ── Stage 2: seed initial elevations + HARD anchor flags ─────────


def _seed_elevations(layout, nodes, bucket_to_idx,
                     dem=None, tile_lat: int = 0, tile_lon: int = 0):
    """Returns ``(elev, is_hard, have_initial)``.

    HARD: only CIFP runway corners.  All other nodes are SOFT — even
    terminals and aprons, per user 2026-05-03 ("only the runway ends
    are immutable truth").

    Soft node seeding priority (highest first):
      1. Existing layout altitude_high/low/altitude/node_altitudes
         (warm-start from a previous solver pass).
      2. Per-vertex DEM sample at the node's (x, y).
      3. Nearest-HARD elevation (cheap geometric backfill).

    The DEM step is what lets a soft node settle at its natural
    terrain elevation when the rest of the graph allows it; cap
    projection in subsequent iterations pulls it down toward HARD
    anchors only where the per-edge grade cap is exceeded.
    """
    from auto_patch.elevation import _sample_dem
    n = len(nodes)
    elev: List[float] = [0.0] * n
    is_hard: List[bool] = [False] * n
    have_initial: List[bool] = [False] * n

    # Runway corners — HARD-anchor every runway segment, sloped or
    # flat.  The runway's elevation profile is authoritative truth
    # for adjacent pavement: when a junction shares a vertex with a
    # runway corner, that vertex must adopt the runway's elevation
    # so cap projection can pull the rest of the junction (and its
    # downstream chain of stubs / aprons) up toward it.
    #
    # Sloped segments are 4-corner rects with altitude_high/low.
    # Flat segments use a single ``altitude=`` tag and may carry an
    # arbitrary number of corners — junctions touching the edge
    # interior get inserted as new shared vertices upstream so the
    # solver gets denser HARD anchors along long flat runs (blast
    # pads, runway-interior flats).
    # Two-pass runway HARD seeding: process non-regraded (CIFP only)
    # shapes first, then regraded shapes (those with node_altitudes
    # from the seam pipeline) — the second pass OVERRIDES any shared
    # corner the first pass set.  This ensures that when a runway is
    # segmented into sub-rects and only the seam-crossing sub-rect
    # was regraded, the regraded values propagate to its shared
    # threshold corners with adjacent sub-rects.
    for pass_node_alts in (False, True):
        for s in layout.shapes:
            # ROLE_RUNWAY_CROSSING shares the runway HARD-anchor
            # path: its ``node_altitudes`` come from runway-segment
            # interpolation in ``_resolve_runway_crossings`` and
            # are authoritative; the solver must not reshape them.
            if s.role not in (ROLE_RUNWAY, ROLE_RUNWAY_CROSSING):
                continue
            if s.polygon is None or s.polygon.is_empty:
                continue
            has_node_alts = bool(s.node_altitudes)
            if has_node_alts != pass_node_alts:
                continue
            coords = _open_ring(list(s.polygon.exterior.coords))
            if len(coords) < 3:
                continue
            if s.altitude_high is not None and s.altitude_low is not None:
                if len(coords) != 4:
                    continue
                per = [s.altitude_high, s.altitude_low,
                       s.altitude_low, s.altitude_high]
            elif s.altitude is not None:
                per = [float(s.altitude)] * len(coords)
            elif s.node_altitudes:
                per = [float(a) for a in s.node_altitudes[:len(coords)]]
                if len(per) < len(coords):
                    per += [per[-1]] * (len(coords) - len(per))
            else:
                continue
            for (x, y), a in zip(coords, per):
                k = layout.canonical_points.get_or_add(float(x), float(y))
                idx = bucket_to_idx.get(k)
                if idx is None:
                    continue
                # Pass 1 (CIFP): only set if not already HARD.
                # Pass 2 (regraded): always override.
                if pass_node_alts or not is_hard[idx]:
                    elev[idx] = float(a)
                    is_hard[idx] = True
                    have_initial[idx] = True

    # Per user 2026-05-13: seam vertices are HARD anchors with
    # OVERRIDE priority over runway CIFP corners.  When a runway
    # interior vertex is on a tile-boundary seam, its DEM altitude
    # (already written into node_altitudes by apply_seam_dem_anchors)
    # wins over the CIFP-interpolated value at the same position.
    # Architecturally: seam wins because terrain mesh at the tile
    # boundary is pinned to raw HGT by Ortho4XP's preserve_boundary,
    # and we need pavement to match terrain there to avoid a visible
    # cliff in X-Plane.
    seam_keys = getattr(layout, "_seam_anchor_keys", None) or set()
    if seam_keys:
        for s in layout.shapes:
            if s.polygon is None or s.polygon.is_empty:
                continue
            if not s.node_altitudes:
                continue
            coords = _open_ring(list(s.polygon.exterior.coords))
            if len(coords) < 3:
                continue
            alts = list(s.node_altitudes[:len(coords)])
            for (x, y), a in zip(coords, alts):
                # Match the bucket convention used by seam_anchors.
                from ..layout import SHARED_VERTEX_TOL_M
                bk_s = 1.0 / SHARED_VERTEX_TOL_M
                seam_bk = (int(round(x * bk_s)), int(round(y * bk_s)))
                if seam_bk not in seam_keys:
                    continue
                k = layout.canonical_points.get_or_add(float(x), float(y))
                idx = bucket_to_idx.get(k)
                if idx is None:
                    continue
                # Seam wins: override any existing HARD value too.
                elev[idx] = float(a)
                is_hard[idx] = True
                have_initial[idx] = True

    # Warm-start soft nodes.
    for s in layout.shapes:
        if s.role not in PAVEMENT_ROLES or s.role == ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = _open_ring(list(s.polygon.exterior.coords))
        if (s.altitude_high is not None and s.altitude_low is not None
                and len(coords) == 4):
            per = [s.altitude_high, s.altitude_low,
                   s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            per = [float(s.altitude)] * len(coords)
        elif s.node_altitudes:
            per = [float(a) for a in s.node_altitudes[:len(coords)]]
            if len(per) < len(coords):
                per += [per[-1]] * (len(coords) - len(per))
        else:
            continue
        for (x, y), a in zip(coords, per):
            k = layout.canonical_points.get_or_add(float(x), float(y))
            idx = bucket_to_idx.get(k)
            if idx is None or is_hard[idx] or have_initial[idx]:
                continue
            elev[idx] = float(a)
            have_initial[idx] = True

    # DEM seed for soft nodes that warm-start didn't cover.
    if dem is not None and any(not h for h in have_initial):
        for i in range(n):
            if have_initial[i]:
                continue
            x, y = nodes[i]
            lat, lon = layout.m_to_ll(x, y)
            e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            if e is not None:
                elev[i] = float(e)
                have_initial[i] = True

    # Backfill any node still without an initial value via nearest
    # HARD anchor's elevation (cheap geometric pass).
    if any(not h for h in have_initial):
        hard_pts = [(nodes[i][0], nodes[i][1], elev[i])
                    for i in range(n) if is_hard[i]]
        for i in range(n):
            if have_initial[i]:
                continue
            x, y = nodes[i]
            best_d2 = float("inf")
            best_e = 0.0
            for hx, hy, he in hard_pts:
                d2 = (hx - x) ** 2 + (hy - y) ** 2
                if d2 < best_d2:
                    best_d2 = d2
                    best_e = he
            elev[i] = best_e
            have_initial[i] = True

    return elev, is_hard, have_initial


# ── Stage 3: edge construction (the per-axis rule lives here) ─────


JUNCTION_AXIS_PERP_TOL_M = 15.0  # taxi half-width + small slack


def _collect_junction_axes(layout, polygon):
    """Return every centerline / runway long-axis that passes
    through ``polygon`` — used by ``_build_edges`` to apply
    per-axis grade constraints to a junction.

    Sources:
    * ``layout.apt_taxi_centerlines`` — full apt.dat taxi network.
    * Each runway segment's long-axis (midpoints of its two short
      edges), for runway-crossing junctions.
    """
    from shapely.geometry import LineString
    axes = []
    apt_lines = getattr(layout, "apt_taxi_centerlines", None) or []
    for item in apt_lines:
        ln = item[0] if isinstance(item, tuple) else item
        if ln is None or ln.is_empty:
            continue
        try:
            if polygon.intersects(ln):
                axes.append(ln)
        except _GEOM_EXC:
            continue
    for s2 in layout.shapes:
        if s2.role != ROLE_RUNWAY:
            continue
        if s2.polygon is None or s2.polygon.is_empty:
            continue
        try:
            if not polygon.intersects(s2.polygon):
                continue
            rc = list(s2.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                 0.5 * (rc[0][1] + rc[3][1]))
        b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                 0.5 * (rc[1][1] + rc[2][1]))
        try:
            axes.append(LineString([a_mid, b_mid]))
        except _GEOM_EXC:
            continue
    return axes


def _build_edges(layout, bucket_to_idx
                  ) -> Tuple[Dict[Tuple[int, int], float],
                             Dict[Tuple[int, int], float]]:
    """Build the unified graph's edge list with role-aware geometry.

    For RECT roles (taxi rects, runway segments): ring edges only.
    The within-rect constraint is axial; cross-section flatness
    groups handle the perpendicular dimension.

    For JUNCTION: ring edges + per-axis grade edges.  For each
    apt.dat taxi centerline or runway long-axis that passes
    through the polygon, the vertices within
    ``JUNCTION_AXIS_PERP_TOL_M`` perpendicular of the axis form
    a group; edges between group members use the ALONG-AXIS
    projected distance as the edge length.  Vertices not near any
    axis are bound only by ring continuity.  Per user 2026-05-18:
    a junction may slope in multiple directions along its
    converging centerlines and 1.5 % is enforced ALONG each axis,
    NOT cross-axially.

    For APRON / TERMINAL: ring edges + all-pair Euclidean spatial
    edges.  Aprons must satisfy 1.5 % across the entire interior
    surface (every-direction cap).

    Per-edge cap = role's max grade × edge length.  When two
    shapes contribute to the same vertex pair, the tighter cap
    wins.
    """
    from shapely.geometry import Point
    edge_grade: Dict[Tuple[int, int], float] = {}
    edge_length: Dict[Tuple[int, int], float] = {}

    def _add_edge(ui, uj, length, gr):
        if ui is None or uj is None or ui == uj:
            return
        if length < 0.1:
            return
        key = (ui, uj) if ui < uj else (uj, ui)
        cur_g = edge_grade.get(key, float("inf"))
        if gr < cur_g:
            edge_grade[key] = gr
        cur_l = edge_length.get(key, length)
        edge_length[key] = min(cur_l, length)

    for s in layout.shapes:
        if s.role not in PAVEMENT_ROLES:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = _open_ring(list(s.polygon.exterior.coords))
        if len(coords) < 2:
            continue
        gr = _role_grade(s.role)
        m = len(coords)
        node_idx = [bucket_to_idx.get(layout.canonical_points.get_or_add(float(x), float(y)))
                    for x, y in coords]
        # Ring edges (every shape).
        for i in range(m):
            j = (i + 1) % m
            x1, y1 = coords[i]
            x2, y2 = coords[j]
            length = math.hypot(x2 - x1, y2 - y1)
            _add_edge(node_idx[i], node_idx[j], length, gr)
        # Rects / runways / runway-crossings: ring-only, no spatial
        # pairs.  Runway-crossings are HARD-anchored via the
        # runway-interpolated ``node_altitudes`` seed; spatial
        # edges would constrain them needlessly.
        if (s.role in SLOPING_RECT_ROLES
                or s.role in (ROLE_RUNWAY, ROLE_RUNWAY_CROSSING)):
            continue
        # Junction / apron / terminal: all-pair Euclidean within
        # the polygon.  Per user 2026-05-18: "a junction should not
        # exceed 1.5 % across ANY portion, not just along its
        # edge."  Same all-pair rule as aprons — junction grade
        # holds across the entire interior surface, not only along
        # converging centerlines.
        for i in range(m):
            xi, yi = coords[i]
            for j in range(i + 2, m):
                if i == 0 and j == m - 1:
                    continue  # ring-wrap pair already added
                xj, yj = coords[j]
                length = math.hypot(xj - xi, yj - yi)
                _add_edge(node_idx[i], node_idx[j], length, gr)

    return edge_grade, edge_length


def _build_adjacency(n, edge_grade, edge_length):
    adj: List[List[Tuple[int, float, float]]] = [[] for _ in range(n)]
    for (u, v), gr in edge_grade.items():
        L = edge_length[(u, v)]
        adj[u].append((v, L, gr))
        adj[v].append((u, L, gr))
    return adj


# ── Stage 4: terminal flatness groups ────────────────────────────


def _build_rect_cross_section_groups(layout, bucket_to_idx):
    """Per user 2026-05-03: a taxi rect slopes along its
    ``source_axis`` only — never perpendicular to it.  The two
    corners at each axis-end (the rect's "axis-end edge", or what
    the legacy code called the "short end") share one elevation;
    the cross-section is exactly flat.

    Each rect contributes TWO flatness groups, one per axis-end.
    The solver runs ``_equalize_groups`` on these every iteration,
    same mechanic as terminal flatness.  ``altitude_high`` /
    ``altitude_low`` written by the solver thus correspond to one
    value per axis-end with no perpendicular component.

    Returns ``[[idx_a, idx_b], ...]`` — one group per axis-end (two
    groups per rect).
    """
    from auto_patch.elevation import (
        _corner_elevation_bucket, _short_end_pairs_by_axis,
    )
    groups: List[List[int]] = []
    for s in layout.shapes:
        if s.role not in SLOPING_RECT_ROLES:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = _open_ring(list(s.polygon.exterior.coords))
        if len(coords) != 4:
            continue
        if s.source_axis is None or s.source_axis.is_empty:
            continue
        sp, ep = _short_end_pairs_by_axis(coords, s.source_axis)
        if sp is None:
            continue
        for pair in (sp, ep):
            idxs = []
            for i in pair:
                if 0 <= i < len(coords):
                    k = layout.canonical_points.get_or_add(float(coords[i][0]), float(coords[i][1]))
                    if k in bucket_to_idx:
                        idxs.append(bucket_to_idx[k])
            if len(idxs) >= 2 and idxs[0] != idxs[1]:
                groups.append(idxs)
    return groups


def _build_terminal_groups(layout, bucket_to_idx):
    """Each terminal contributes one group of node indices that
    must share a single elevation (the flatness constraint).
    """
    groups: List[List[int]] = []
    for s in layout.shapes:
        if s.role != ROLE_TERMINAL:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = _open_ring(list(s.polygon.exterior.coords))
        idxs = []
        for x, y in coords:
            k = layout.canonical_points.get_or_add(float(x), float(y))
            if k in bucket_to_idx:
                idxs.append(bucket_to_idx[k])
        if len(idxs) >= 2:
            groups.append(idxs)
    return groups


# ── Stage 5: damped Jacobi + cap projection iteration ────────────


def _equalize_groups(elev, is_hard, groups):
    """Set every member of each group to the group's mean elevation.
    HARD members are immutable; if a HARD member exists, the group
    averages the SOFT members and pulls them toward the HARD value
    (cap projection in subsequent sweeps will pull the rest of the
    graph back into compliance).

    Used for both terminal flatness (one group per terminal, all
    corners) and rect cross-section flatness (one group per rect
    short-end, two corners each).
    """
    for grp in groups:
        if not grp:
            continue
        hard_in_grp = [i for i in grp if is_hard[i]]
        if hard_in_grp:
            target = elev[hard_in_grp[0]]
            for i in grp:
                if not is_hard[i]:
                    elev[i] = target
        else:
            avg = sum(elev[i] for i in grp) / len(grp)
            for i in grp:
                elev[i] = avg


def _run_jacobi(elev, is_hard, adj, edge_list, edge_grade,
                edge_length, terminal_groups,
                rect_flat_groups,
                max_iters, tol_m) -> int:
    """Cap-projection-only relaxation (user 2026-05-03).

    Earlier iterations of this solver included a damped-Jacobi
    neighbour-average step before cap projection.  Jacobi pulls
    every soft node toward the weighted mean of its neighbours,
    which propagates HARD anchor values up through the graph and
    over-flattens DEM-seeded soft nodes (CYXY taxi E ended up at
    700 m next to a 700 m runway, even though DEM said 715 m and
    the grade chain through stubs allowed reaching it).

    Cap-projection-only preserves DEM-seeded values that satisfy
    every per-edge cap.  Soft nodes only move when an edge cap is
    violated, and only by the excess.  Convergence is to
    ``min(DEM_seed, max_reachable_from_HARD)`` per node, which is
    exactly the user's "reach the highs and lows in DEM that are
    possible within grade limits" rule.

    Two equality constraint groups run each iteration: terminals
    (all corners equal) and rect axis-end pairs (per user 2026-
    05-03: rects slope along source_axis only, axis-perpendicular
    is flat).
    """
    n = len(elev)
    for it in range(max_iters):
        prev_elev = list(elev)
        # 1) Multi-sweep edge grade-cap projection — only force
        # acting on soft nodes.  Each sweep visits every edge; an
        # edge is projected (excess split symmetrically for
        # soft-soft, asymmetrically toward the soft side for
        # soft-hard) only when it currently violates its cap.
        for _sweep in range(CAP_SWEEPS_PER_ITER):
            any_proj = False
            for (u, v) in edge_list:
                L = edge_length[(u, v)]
                gr = edge_grade[(u, v)]
                diff = elev[u] - elev[v]
                cap = L * gr
                if abs(diff) <= cap:
                    continue
                excess = abs(diff) - cap
                sign = 1 if diff > 0 else -1
                if is_hard[u] and is_hard[v]:
                    continue
                if is_hard[u]:
                    elev[v] += sign * excess
                    any_proj = True
                elif is_hard[v]:
                    elev[u] -= sign * excess
                    any_proj = True
                else:
                    half = 0.5 * excess * sign
                    elev[u] -= half
                    elev[v] += half
                    any_proj = True
            if not any_proj:
                break
        # 2) Equality constraints — terminal flatness and rect
        # axis-end (cross-section) flatness.
        _equalize_groups(elev, is_hard, terminal_groups)
        _equalize_groups(elev, is_hard, rect_flat_groups)
        # 3) Convergence.
        max_change = 0.0
        for i in range(n):
            if is_hard[i]:
                continue
            d = abs(prev_elev[i] - elev[i])
            if d > max_change:
                max_change = d
        if max_change < tol_m:
            return it + 1
    return max_iters


# ── Stage 6: write elevations back to layout shapes ──────────────


def _writeback(layout, elev, bucket_to_idx):
    """Apply solved elevations to layout shapes.

    For taxi rects: ensure the polygon's vertex order is canonical
    (corners 0, 3 at the higher axis-end; corners 1, 2 at the
    lower).  The OSM emit interpolates altitude_high/low across
    polygon corners via the legacy convention ``[high, low, low,
    high]`` for indices 0..3 — that mapping is wrong for any rect
    whose polygon happens to be ring-rotated relative to canonical,
    leading to a phantom perpendicular slope (the source of the
    user 2026-05-03 SPJC F-stub report).  Rotating the ring at
    writeback aligns the convention with the actual axis-end
    geometry.
    """
    from shapely.geometry import Polygon
    from auto_patch.elevation import (
        _corner_elevation_bucket, _short_end_pairs_by_axis,
    )
    n_terms = n_rects = n_juncs = 0
    for s in layout.shapes:
        if s.role not in PAVEMENT_ROLES:
            continue
        # Runway shapes are normally skipped (their altitudes come
        # from CIFP — HARD-anchored, immutable through the solver).
        # Exceptions where the writeback DOES run:
        #   * Seam-converted runway sub-rects (user 2026-05-13): they
        #     have ``node_altitudes`` set; we write per-vertex
        #     solver-output altitudes so shared corners with adjacent
        #     sub-rects agree on the regraded value.
        #   * Non-4-corner runway shapes (user 2026-05-19): a runway
        #     segment that lost its canonical 4-corner form through
        #     downstream geometry passes (crossing union, snap-to-
        #     corner, etc.) is no longer a sloped rect — its
        #     altitude_high/low tags are stale because X-Plane's
        #     planar 4-corner convention requires exactly 4 corners.
        #     Convert to ``node_altitudes`` so the OSM emit + the
        #     no-vertex-on-sloping-edge invariant treat it as the
        #     non-rect it actually is.
        if s.role == ROLE_RUNWAY and not s.node_altitudes:
            _rc_check = list(s.polygon.exterior.coords) if s.polygon else []
            if _rc_check and _rc_check[0] == _rc_check[-1]:
                _rc_check = _rc_check[:-1]
            if len(_rc_check) == 4:
                continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except _GEOM_EXC:
            continue
        ring_closed = coords and coords[0] == coords[-1]
        coords_open = coords[:-1] if ring_closed else coords
        corner_elevs = _read_corner_elevs(
            coords_open, elev, bucket_to_idx, layout)
        if corner_elevs is None:
            continue
        if s.role == ROLE_TERMINAL:
            # Terminal is FLAT (per user 2026-05-18: a terminal sits
            # on one floor altitude).  The terminal-flatness equality
            # group already enforced this in the solver; average is
            # just a defensive round.
            avg = sum(corner_elevs) / len(corner_elevs)
            s.altitude = round(float(avg), 1)
            s.altitude_high = None
            s.altitude_low = None
            s.node_altitudes = None
            n_terms += 1
        elif s.role == ROLE_APRON:
            # Per user 2026-05-18: aprons are NOT 100 % flat — they
            # satisfy 1.5 % across their surface, NOT zero gradient.
            # Keep the solver's per-corner altitudes (which it
            # already constrained via all-pair Euclidean edges) so
            # adjacent aprons that share corners don't end up at
            # 4-8 m cliff steps (each apron previously averaged to
            # its own single altitude → adjacent aprons diverged).
            alts = [round(float(e), 1) for e in corner_elevs]
            if ring_closed:
                alts.append(alts[0])
            s.node_altitudes = alts
            s.altitude = None
            s.altitude_high = None
            s.altitude_low = None
            n_terms += 1
        elif s.role in SLOPING_RECT_ROLES:
            # Per user 2026-05-13: keep node_altitudes when the shape
            # came in with them — even for 4-corner shapes.  This
            # preserves per-vertex precision for runway sub-rects
            # adjacent to seam-affected sub-rects: their shared
            # corners receive HARD seam altitudes that aren't coplanar
            # with the other 2 CIFP corners, so altitude_high/low
            # (which assumes a planar surface) would average and
            # introduce a > 1 m step at the shared boundary.
            had_node_alts = s.node_altitudes is not None
            if len(coords_open) == 4 and not had_node_alts:
                new_coords, hi, lo = _canonicalise_rect(
                    coords_open, corner_elevs, s.source_axis,
                    _short_end_pairs_by_axis)
                if new_coords is None:
                    continue
                if new_coords != coords_open:
                    s.polygon = Polygon(new_coords + [new_coords[0]])
                s.altitude_high = round(float(hi), 1)
                s.altitude_low = round(float(lo), 1)
                s.altitude = None
                s.node_altitudes = None
                n_rects += 1
            else:
                alts = [round(float(e), 1) for e in corner_elevs]
                if ring_closed:
                    alts.append(alts[0])
                s.node_altitudes = alts
                s.altitude_high = None
                s.altitude_low = None
                s.altitude = None
                n_rects += 1
        elif s.role == ROLE_JUNCTION:
            alts = [round(float(e), 1) for e in corner_elevs]
            if ring_closed:
                alts.append(alts[0])
            s.node_altitudes = alts
            s.altitude = None
            n_juncs += 1
        elif s.role == ROLE_RUNWAY:
            # Seam-converted runway sub-rect — write per-vertex
            # altitudes (the only runway shapes that reach here have
            # node_altitudes pre-set; the skip-guard above filters
            # the CIFP-only altitude_high/low ones).
            alts = [round(float(e), 1) for e in corner_elevs]
            if ring_closed:
                alts.append(alts[0])
            s.node_altitudes = alts
            s.altitude = None
            s.altitude_high = None
            s.altitude_low = None
            n_rects += 1
    return n_terms, n_rects, n_juncs


def _read_corner_elevs(coords_open, elev, bucket_to_idx, layout=None):
    out = []
    for x, y in coords_open:
        idx = bucket_to_idx.get(layout.canonical_points.get_or_add(float(x), float(y)))
        if idx is None:
            return None
        out.append(elev[idx])
    return out


def _canonicalise_rect(coords_open, corner_elevs, source_axis,
                        short_end_pairs_fn):
    """Rotate a rect's 4-vertex ring (and its corner elevations)
    so corners 0, 3 are at the higher axis-end and 1, 2 at the
    lower.  Returns ``(new_coords, hi, lo)`` or ``(None, ...)`` if
    rotation can't be determined.
    """
    sp, ep = short_end_pairs_fn(coords_open, source_axis)
    if sp is None:
        sp, ep = (0, 3), (1, 2)
    a_avg = (corner_elevs[sp[0]] + corner_elevs[sp[1]]) / 2.0
    b_avg = (corner_elevs[ep[0]] + corner_elevs[ep[1]]) / 2.0
    high_pair = sp if a_avg >= b_avg else ep
    hi, lo = max(a_avg, b_avg), min(a_avg, b_avg)
    # Rotation that makes high_pair == (0, 3).
    rotation = _rotation_for_high_pair(high_pair)
    if rotation == 0:
        return list(coords_open), hi, lo
    new_coords = [coords_open[(i - rotation) % 4]
                  for i in range(4)]
    return new_coords, hi, lo


def _rotation_for_high_pair(high_pair) -> int:
    """Return the right-shift k such that rotating the 4-vertex
    ring by k positions makes ``high_pair`` map to ``(0, 3)``.

    Mapping: under right-shift k, old index ``i`` becomes new
    index ``(i + k) % 4``.  We solve for k so that
    ``{(high_pair[0] + k) % 4, (high_pair[1] + k) % 4} == {0, 3}``.
    """
    target = {0, 3}
    a, b = high_pair
    for k in range(4):
        if {(a + k) % 4, (b + k) % 4} == target:
            return k
    return 0


def _report(icao, iters_used, max_iters, elapsed,
             n_terms, n_rects, n_juncs):
    import O4_UI_Utils as UI
    UI.vprint(1,
        f"  [pav-builder] {icao}: per-surface Jacobi solver "
        f"converged in {iters_used}/{max_iters} iters "
        f"({elapsed:.2f} s); applied to {n_terms} terminal/apron(s), "
        f"{n_rects} rect(s), {n_juncs} junction(s).")
