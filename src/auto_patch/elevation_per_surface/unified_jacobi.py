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

from auto_patch.elevation import APRON_MAX_GRADE, TAXI_MAX_GRADE
from auto_patch.layout import (
    ROLE_APRON, ROLE_BOUNDARY, ROLE_CROSS_CONNECTOR, ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL, ROLE_RUNWAY, ROLE_SECONDARY_PARALLEL,
    ROLE_STUB, ROLE_TERMINAL,
)


SLOPING_RECT_ROLES = (
    ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
    ROLE_STUB, ROLE_CROSS_CONNECTOR,
)

PAVEMENT_ROLES = {
    ROLE_RUNWAY, *SLOPING_RECT_ROLES,
    ROLE_APRON, ROLE_TERMINAL, ROLE_JUNCTION,
}

CAP_SWEEPS_PER_ITER = 5


def _role_grade(role: str) -> float:
    """Per-role max grade cap (user 2026-05-02):

    * Runway, taxi rect, **junction**: 1.5 % (``TAXI_MAX_GRADE``).
      Junctions are multi-directional but still capped at the
      taxiway grade — they're transition zones a plane taxis across.
    * Apron, terminal: 1.0 % (``APRON_MAX_GRADE``) — tighter FAA cap
      for parking surfaces.
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
    from auto_patch.elevation import _corner_elevation_bucket

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
    """Assign one node index per unique vertex bucket across all
    pavement-role shapes.  Returns ``(nodes, bucket_to_idx)``.
    """
    from auto_patch.elevation import _corner_elevation_bucket
    bucket_to_idx: Dict[Tuple[int, int], int] = {}
    nodes: List[Tuple[float, float]] = []
    for s in layout.shapes:
        if s.role not in PAVEMENT_ROLES:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = _open_ring(list(s.polygon.exterior.coords))
        except Exception:
            continue
        for x, y in coords:
            b = _corner_elevation_bucket(x, y)
            if b not in bucket_to_idx:
                bucket_to_idx[b] = len(nodes)
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
    from auto_patch.elevation import _corner_elevation_bucket
    n = len(nodes)
    elev: List[float] = [0.0] * n
    is_hard: List[bool] = [False] * n
    have_initial: List[bool] = [False] * n

    # Runway corners — HARD-anchor every runway segment, sloped or
    # flat.  The runway's elevation profile is authoritative truth
    # for adjacent pavement: when a junction shares a vertex with a
    # runway corner, that vertex must adopt the runway's elevation
    # so cap projection can pull the rest of the junction (and its
    # downstream chain of stubs / aprons) up toward it.  The
    # original gate required both altitude_high and altitude_low,
    # which silently dropped flat segments (single ``altitude=``)
    # — leaving long stretches of runway interior with no HARD
    # anchors and adjacent taxiways stuck at terrain.
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = _open_ring(list(s.polygon.exterior.coords))
        if len(coords) != 4:
            continue
        if s.altitude_high is not None and s.altitude_low is not None:
            per = [s.altitude_high, s.altitude_low,
                   s.altitude_low, s.altitude_high]
        elif s.altitude is not None:
            per = [float(s.altitude)] * 4
        else:
            continue
        for (x, y), a in zip(coords, per):
            b = _corner_elevation_bucket(x, y)
            idx = bucket_to_idx.get(b)
            if idx is not None and not is_hard[idx]:
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
            b = _corner_elevation_bucket(x, y)
            idx = bucket_to_idx.get(b)
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


def _build_edges(layout, bucket_to_idx
                  ) -> Tuple[Dict[Tuple[int, int], float],
                             Dict[Tuple[int, int], float]]:
    """Build the unified graph's edge list with role-aware geometry.

    For RECT roles: ring edges only (4 edges per polygon).  No
    spatial pairs — the within-rect constraint is axial, not
    cross-axial.

    For JUNCTION / APRON / TERMINAL roles: ring edges + all-pair
    Euclidean spatial edges within the polygon.  No radius cap —
    every vertex pair on the same multi-directional surface is
    constrained at the role's grade × Euclidean distance.

    Per-edge cap = role's max grade × Euclidean length.  When two
    shapes contribute to the same vertex pair (e.g. a shared edge),
    the tighter cap wins.
    """
    from auto_patch.elevation import _corner_elevation_bucket
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
        node_idx = [bucket_to_idx.get(_corner_elevation_bucket(x, y))
                    for x, y in coords]
        # Ring edges (every shape).
        for i in range(m):
            j = (i + 1) % m
            x1, y1 = coords[i]
            x2, y2 = coords[j]
            length = math.hypot(x2 - x1, y2 - y1)
            _add_edge(node_idx[i], node_idx[j], length, gr)
        # Spatial pairs only for multi-directional roles.
        if s.role in SLOPING_RECT_ROLES or s.role == ROLE_RUNWAY:
            continue
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
                    b = _corner_elevation_bucket(*coords[i])
                    if b in bucket_to_idx:
                        idxs.append(bucket_to_idx[b])
            if len(idxs) >= 2 and idxs[0] != idxs[1]:
                groups.append(idxs)
    return groups


def _build_terminal_groups(layout, bucket_to_idx):
    """Each terminal contributes one group of node indices that
    must share a single elevation (the flatness constraint).
    """
    from auto_patch.elevation import _corner_elevation_bucket
    groups: List[List[int]] = []
    for s in layout.shapes:
        if s.role != ROLE_TERMINAL:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = _open_ring(list(s.polygon.exterior.coords))
        idxs = []
        for x, y in coords:
            b = _corner_elevation_bucket(x, y)
            if b in bucket_to_idx:
                idxs.append(bucket_to_idx[b])
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
        if s.role not in PAVEMENT_ROLES or s.role == ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            coords = list(s.polygon.exterior.coords)
        except Exception:
            continue
        ring_closed = coords and coords[0] == coords[-1]
        coords_open = coords[:-1] if ring_closed else coords
        corner_elevs = _read_corner_elevs(
            coords_open, elev, bucket_to_idx)
        if corner_elevs is None:
            continue
        if s.role == ROLE_TERMINAL or s.role == ROLE_APRON:
            avg = sum(corner_elevs) / len(corner_elevs)
            s.altitude = round(float(avg), 1)
            s.altitude_high = None
            s.altitude_low = None
            s.node_altitudes = None
            n_terms += 1
        elif s.role in SLOPING_RECT_ROLES:
            if len(coords_open) != 4:
                continue
            new_coords, hi, lo = _canonicalise_rect(
                coords_open, corner_elevs, s.source_axis,
                _short_end_pairs_by_axis)
            if new_coords is None:
                continue
            if new_coords != coords_open:
                # Rotate the rect's polygon to canonical order so
                # the OSM-emit convention `[high, low, low, high]`
                # matches the actual axis-end geometry.
                s.polygon = Polygon(new_coords + [new_coords[0]])
            s.altitude_high = round(float(hi), 1)
            s.altitude_low = round(float(lo), 1)
            s.altitude = None
            s.node_altitudes = None
            n_rects += 1
        elif s.role == ROLE_JUNCTION:
            alts = [round(float(e), 1) for e in corner_elevs]
            if ring_closed:
                alts.append(alts[0])
            s.node_altitudes = alts
            s.altitude = None
            n_juncs += 1
    return n_terms, n_rects, n_juncs


def _read_corner_elevs(coords_open, elev, bucket_to_idx):
    from auto_patch.elevation import _corner_elevation_bucket
    out = []
    for x, y in coords_open:
        idx = bucket_to_idx.get(_corner_elevation_bucket(x, y))
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
    try:
        import O4_UI_Utils as UI
        UI.vprint(1,
            f"  [pav-builder] {icao}: per-surface Jacobi solver "
            f"converged in {iters_used}/{max_iters} iters "
            f"({elapsed:.2f} s); applied to {n_terms} terminal/apron(s), "
            f"{n_rects} rect(s), {n_juncs} junction(s).")
    except Exception:
        pass
