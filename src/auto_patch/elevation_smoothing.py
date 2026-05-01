"""2D-Euclidean polygon-grid elevation smoothing.

Implements a per-shape elevation field over an arbitrary polygon
by gridding the polygon, anchoring HARD vertices (rect /
runway / terminal corners on the boundary) at their assigned
altitude, and Jacobi-relaxing free cells to satisfy a 2D
neighbour grade cap.  Used by ``triangulation._triangulate_junctions``
when a junction's anchor set isn't co-planar enough for a planar
fit.

Public API:
    _smooth_polygon_grid(polygon, hard_anchors, dem, tile_lat,
                         tile_lon, layout_anchor, *, grid_step_m,
                         max_iters, convergence_tol_m)
"""
from __future__ import annotations

import math
from typing import List, Tuple

from shapely.geometry import Polygon

from .elevation import (
    ELEVATION_GRID_STEP_M,
    ELEVATION_SMOOTH_CONVERGE_M,
    ELEVATION_SMOOTH_MAX_ITERS,
    TAXI_MAX_GRADE,
    _sample_dem,
)
from .layout import R_EARTH


__all__ = ["_smooth_polygon_grid"]


def _smooth_polygon_grid(
        polygon: Polygon,
        hard_anchors: List[Tuple[float, float, float]],
        dem,
        tile_lat: int,
        tile_lon: int,
        layout_anchor: Tuple[float, float],
        grid_step_m: float = ELEVATION_GRID_STEP_M,
        max_iters: int = ELEVATION_SMOOTH_MAX_ITERS,
        convergence_tol_m: float = ELEVATION_SMOOTH_CONVERGE_M):
    """Smooth a 2D elevation field within ``polygon``.

    ``hard_anchors`` is a list of ``(x, y, elev)`` triples whose
    elevation must be preserved exactly (rect / runway / terminal
    corners; cross-junction shared vertices already locked from a
    previous outer iteration).

    Returns ``(sampler, sample)`` where ``sampler(x, y)`` returns
    the bilinearly-interpolated elevation, or ``None`` if the
    polygon was too small for grid construction.  When no anchors
    are reachable for the polygon, falls back to the global
    elevation graph / DEM at each query point — same path
    ``_vertex_elev_anchored`` step 5 takes today.
    """
    import numpy as np

    # Degenerate polygon — return a graph/DEM sampler.
    if polygon is None or polygon.is_empty:
        return None
    minx, miny, maxx, maxy = polygon.bounds
    span_x = maxx - minx
    span_y = maxy - miny
    if span_x <= 0 or span_y <= 0:
        return None

    # Build grid bbox with a 1-cell margin so boundary vertices land
    # comfortably inside the active region.
    margin = grid_step_m
    minx -= margin
    miny -= margin
    maxx += margin
    maxy += margin
    nx = max(2, int(math.ceil((maxx - minx) / grid_step_m)) + 1)
    ny = max(2, int(math.ceil((maxy - miny) / grid_step_m)) + 1)
    # Cap grid size — should never trigger for realistic airport
    # polygons but prevents pathological all-airport polygons from
    # eating all RAM.
    if nx * ny > 250_000:
        return None

    # Cell-center coords.
    xs = minx + np.arange(nx) * grid_step_m
    ys = miny + np.arange(ny) * grid_step_m

    # ── INSIDE mask ─────────────────────────────────────────────
    # A cell is INSIDE if its center sits inside the polygon OR
    # within ``grid_step_m`` of the polygon boundary (so cells
    # straddling the boundary stay active).  Bilinear sampling at
    # boundary vertices needs the surrounding cells active.
    from shapely.prepared import prep
    prepped = prep(polygon)
    boundary = polygon.boundary
    inside = np.zeros((nx, ny), dtype=bool)
    for i in range(nx):
        for j in range(ny):
            pt = Point(float(xs[i]), float(ys[j]))
            if prepped.contains(pt):
                inside[i, j] = True
            else:
                try:
                    if boundary.distance(pt) < grid_step_m:
                        inside[i, j] = True
                except Exception:
                    pass

    if not inside.any():
        return None

    # ── Pin anchors ─────────────────────────────────────────────
    # Each hard anchor pins its NEAREST INSIDE cell to its elevation.
    # Multiple anchors hitting the same cell are averaged (rare —
    # implies two rect corners within one grid cell, which means
    # they should be the same shared-vertex bucket already).
    pinned = np.zeros((nx, ny), dtype=bool)
    pinned_elev = np.zeros((nx, ny), dtype=float)
    pinned_count = np.zeros((nx, ny), dtype=int)
    eff_anchors: List[Tuple[float, float, float]] = []
    for (ax, ay, az) in hard_anchors:
        i = int(round((ax - minx) / grid_step_m))
        j = int(round((ay - miny) / grid_step_m))
        if 0 <= i < nx and 0 <= j < ny:
            # Snap to the nearest INSIDE cell within a 1-cell radius
            # (anchors at the boundary may map to an INACTIVE cell).
            if not inside[i, j]:
                snapped = False
                for di in range(-1, 2):
                    for dj in range(-1, 2):
                        ii, jj = i + di, j + dj
                        if (0 <= ii < nx and 0 <= jj < ny
                                and inside[ii, jj]):
                            i, j = ii, jj
                            snapped = True
                            break
                    if snapped:
                        break
                if not snapped:
                    continue
            if pinned[i, j]:
                pinned_elev[i, j] = (
                    pinned_elev[i, j] * pinned_count[i, j] + az
                ) / (pinned_count[i, j] + 1)
                pinned_count[i, j] += 1
            else:
                pinned[i, j] = True
                pinned_elev[i, j] = az
                pinned_count[i, j] = 1
            eff_anchors.append((float(xs[i]), float(ys[j]),
                                float(pinned_elev[i, j])))

    # ── Per-cell feasibility band ───────────────────────────────
    INF = float("inf")
    lo = np.full((nx, ny), -INF, dtype=float)
    hi = np.full((nx, ny), INF, dtype=float)
    if eff_anchors:
        XX, YY = np.meshgrid(xs, ys, indexing="ij")
        for (ax, ay, az) in eff_anchors:
            dist = np.hypot(XX - ax, YY - ay)
            band = dist * TAXI_MAX_GRADE
            lo = np.maximum(lo, az - band)
            hi = np.minimum(hi, az + band)

    # ── Initialize cells ────────────────────────────────────────
    # Pinned: their anchor elevation.  Free: prefer the elevation
    # graph's already-grade-compliant network value (smoothed at
    # 1.5 % across the centerline graph), falling back to DEM
    # only when the cell is too far from the network for graph
    # sampling to return.  Using the graph first preserves grade
    # compliance for free cells beyond any anchor cone's reach —
    # raw DEM-init lets real-terrain bumps leak into the polygon's
    # interior, producing wild boundary samples.
    elev = pinned_elev.copy()
    cos0 = math.cos(math.radians(layout_anchor[0]))
    lat0, lon0 = layout_anchor
    for i in range(nx):
        for j in range(ny):
            if not inside[i, j] or pinned[i, j]:
                continue
            d = None
            if dem is not None:
                lat = lat0 + math.degrees(ys[j] / R_EARTH)
                lon = (lon0 + math.degrees(
                    xs[i] / (R_EARTH * cos0)))
                d = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
            l = lo[i, j]
            h = hi[i, j]
            if d is None:
                if l == -INF and h == INF:
                    elev[i, j] = 0.0
                elif l == -INF:
                    elev[i, j] = h
                elif h == INF:
                    elev[i, j] = l
                else:
                    elev[i, j] = 0.5 * (l + h)
            else:
                if l > h:
                    elev[i, j] = 0.5 * (l + h)
                else:
                    elev[i, j] = max(l, min(h, d))

    # ── Iterate Laplacian + grade cap ───────────────────────────
    step_cap = grid_step_m * TAXI_MAX_GRADE
    inside_arr = inside
    pinned_arr = pinned
    free_inside = inside_arr & ~pinned_arr
    damping = 0.5

    for _ in range(max_iters):
        # Laplacian — average of 4 cardinal INSIDE neighbours.
        e_pad = np.pad(elev, 1, mode="edge")
        in_pad = np.pad(inside_arr, 1, mode="constant",
                        constant_values=False)
        n_e = e_pad[2:, 1:-1]
        n_w = e_pad[:-2, 1:-1]
        n_n = e_pad[1:-1, 2:]
        n_s = e_pad[1:-1, :-2]
        m_e = in_pad[2:, 1:-1]
        m_w = in_pad[:-2, 1:-1]
        m_n = in_pad[1:-1, 2:]
        m_s = in_pad[1:-1, :-2]
        cnt = (m_e.astype(np.int8)
               + m_w.astype(np.int8)
               + m_n.astype(np.int8)
               + m_s.astype(np.int8))
        sum_n = (n_e * m_e + n_w * m_w + n_n * m_n + n_s * m_s)
        with np.errstate(divide="ignore", invalid="ignore"):
            mean_n = np.where(cnt > 0,
                              sum_n / np.maximum(cnt, 1),
                              elev)
        target = elev + damping * (mean_n - elev)
        new_elev = np.where(free_inside, target, elev)
        # Clip free cells back into feasibility bands.
        new_elev = np.where(
            free_inside,
            np.maximum(np.minimum(new_elev, hi), lo),
            new_elev)

        # Edge grade cap on each axis-aligned cell pair.
        for axis in (0, 1):
            if axis == 0:
                a = new_elev[:-1, :]
                b = new_elev[1:, :]
                ia = inside_arr[:-1, :]
                ib = inside_arr[1:, :]
                pa = pinned_arr[:-1, :]
                pb = pinned_arr[1:, :]
            else:
                a = new_elev[:, :-1]
                b = new_elev[:, 1:]
                ia = inside_arr[:, :-1]
                ib = inside_arr[:, 1:]
                pa = pinned_arr[:, :-1]
                pb = pinned_arr[:, 1:]
            both_inside = ia & ib
            diff = a - b
            need = both_inside & (np.abs(diff) > step_cap)
            sign = np.sign(diff)
            excess = np.maximum(np.abs(diff) - step_cap, 0.0)
            both_free = (~pa) & (~pb)
            a_pinned = pa & (~pb)
            b_pinned = pb & (~pa)
            half = 0.5 * excess * sign
            corr_a_both = -half
            corr_b_both = half
            corr_b_apf = excess * sign
            corr_a_bpf = -excess * sign
            apply_both = need & both_free
            apply_b_apf = need & a_pinned
            apply_a_bpf = need & b_pinned
            if axis == 0:
                new_elev[:-1, :] = np.where(
                    apply_both, a + corr_a_both, new_elev[:-1, :])
                new_elev[1:, :] = np.where(
                    apply_both, b + corr_b_both, new_elev[1:, :])
                new_elev[1:, :] = np.where(
                    apply_b_apf, b + corr_b_apf, new_elev[1:, :])
                new_elev[:-1, :] = np.where(
                    apply_a_bpf, a + corr_a_bpf, new_elev[:-1, :])
            else:
                new_elev[:, :-1] = np.where(
                    apply_both, a + corr_a_both, new_elev[:, :-1])
                new_elev[:, 1:] = np.where(
                    apply_both, b + corr_b_both, new_elev[:, 1:])
                new_elev[:, 1:] = np.where(
                    apply_b_apf, b + corr_b_apf, new_elev[:, 1:])
                new_elev[:, :-1] = np.where(
                    apply_a_bpf, a + corr_a_bpf, new_elev[:, :-1])
            # Re-pin: pinned cells must never have their elevation
            # mutated by the cap pass.
            new_elev = np.where(pinned_arr, pinned_elev, new_elev)

        max_change = float(np.max(np.abs(new_elev - elev)))
        elev = new_elev
        if max_change < convergence_tol_m:
            break

    # ── Bilinear sampler closure ────────────────────────────────
    def sampler(x: float, y: float) -> Optional[float]:
        fi = (x - minx) / grid_step_m
        fj = (y - miny) / grid_step_m
        i0 = int(math.floor(fi))
        j0 = int(math.floor(fj))
        if i0 < 0:
            i0 = 0
        if j0 < 0:
            j0 = 0
        if i0 > nx - 2:
            i0 = nx - 2
        if j0 > ny - 2:
            j0 = ny - 2
        u = fi - i0
        v = fj - j0
        if u < 0:
            u = 0.0
        elif u > 1:
            u = 1.0
        if v < 0:
            v = 0.0
        elif v > 1:
            v = 1.0
        m00 = inside[i0, j0]
        m10 = inside[i0 + 1, j0]
        m01 = inside[i0, j0 + 1]
        m11 = inside[i0 + 1, j0 + 1]
        if m00 and m10 and m01 and m11:
            c00 = elev[i0, j0]
            c10 = elev[i0 + 1, j0]
            c01 = elev[i0, j0 + 1]
            c11 = elev[i0 + 1, j0 + 1]
            e0 = c00 * (1 - u) + c10 * u
            e1 = c01 * (1 - u) + c11 * u
            return float(e0 * (1 - v) + e1 * v)
        # Fallback: nearest INSIDE cell within a 2-cell radius.
        best_d2 = float("inf")
        best_e: Optional[float] = None
        for di in range(-2, 4):
            for dj in range(-2, 4):
                ii = i0 + di
                jj = j0 + dj
                if (0 <= ii < nx and 0 <= jj < ny
                        and inside[ii, jj]):
                    d2 = (xs[ii] - x) ** 2 + (ys[jj] - y) ** 2
                    if d2 < best_d2:
                        best_d2 = d2
                        best_e = float(elev[ii, jj])
        return best_e

    return sampler


# ── Junction triangulation pass ──────────────────────────────────




