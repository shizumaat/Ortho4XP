"""Phases 3 + 4: junction / apron interior vertex elevations.

For each polygon (junction or apron):

1. Boundary vertices shared with already-solved shapes (rects /
   runways / terminals) inherit those values as anchors.
2. Interior boundary-arc vertices (apt.dat ring points between
   shared corners) are seeded from DEM.
3. Smooth ``node_altitudes`` so that ``|de| ≤ grade × ring_edge_length``
   on consecutive ring edges.  No spatial-pair shortcuts —
   propagation follows the polygon perimeter.

The same routine handles both junctions (1.5 % cap) and aprons
(1.0 % cap); the grade is the only parameter.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from .bfs_propagate import _bucket
from .dem_targets import polygon_vertex_targets

JUNCTION_MAX_GRADE = 0.015
APRON_MAX_GRADE = 0.010


def _ring_edge_caps(coords: List[Tuple[float, float]],
                    max_grade: float) -> List[float]:
    """Per-edge altitude caps for the closed ring of ``coords``."""
    caps = []
    n = len(coords)
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        L = math.hypot(x2 - x1, y2 - y1)
        caps.append(max_grade * L)
    return caps


def _smooth_along_ring(alts: List[float],
                       hard: List[bool],
                       caps: List[float],
                       max_iters: int = 200,
                       tol_m: float = 0.005) -> List[float]:
    """Iterative ring-edge cap projection.  Each iteration sweeps
    every ring edge once; if ``|alts[i] - alts[i+1]|`` exceeds the
    edge's cap, pull the soft endpoint toward the other (or split
    the excess if both are soft).  Hard endpoints never move.
    """
    n = len(alts)
    a = list(alts)
    for _ in range(max_iters):
        max_change = 0.0
        for i in range(n):
            j = (i + 1) % n
            diff = a[i] - a[j]
            cap = caps[i]
            if abs(diff) <= cap:
                continue
            excess = abs(diff) - cap
            sign = 1.0 if diff > 0 else -1.0
            if hard[i] and hard[j]:
                continue
            if hard[i]:
                a[j] += sign * excess
                max_change = max(max_change, excess)
            elif hard[j]:
                a[i] -= sign * excess
                max_change = max(max_change, excess)
            else:
                half = 0.5 * excess * sign
                a[i] -= half
                a[j] += half
                max_change = max(max_change, abs(half))
        if max_change < tol_m:
            break
    return a


def _seed_alts_from_anchors(coords, anchor_map, dem_targets) \
        -> Tuple[List[float], List[bool]]:
    """For each ring vertex, fill from anchor_map if its bucket is
    anchored (returns ``(alt, True)``), else from DEM
    (``(dem, False)``).  Returns ``(alts, hard_flags)``.
    """
    alts: List[float] = []
    hard: List[bool] = []
    for i, (x, y) in enumerate(coords):
        b = _bucket(x, y)
        a = anchor_map.get(b)
        if a is not None:
            alts.append(float(a))
            hard.append(True)
            continue
        d = dem_targets[i]
        if d is None:
            # Fall back to nearest hard neighbour's altitude (filled
            # in below if any hard exists).  Use 0 placeholder.
            alts.append(0.0)
        else:
            alts.append(float(d))
        hard.append(False)
    # Backfill any vertex with no DEM and no anchor by interpolating
    # from the nearest hard neighbour.
    if any(a == 0.0 and not h for a, h in zip(alts, hard)):
        hard_idxs = [i for i, h in enumerate(hard) if h]
        if hard_idxs:
            ref = alts[hard_idxs[0]]
            for i, (a, h) in enumerate(zip(alts, hard)):
                if not h and a == 0.0:
                    alts[i] = ref
    return alts, hard


def solve_polygon_field(layout, shape, dem, tile_lat: int,
                         tile_lon: int,
                         anchor_map: Dict[Tuple[int, int], float],
                         max_grade: float) -> Optional[List[float]]:
    """Compute ``node_altitudes`` for one polygon (junction / apron)
    by ring-edge smoothing from boundary anchors + DEM.

    ``anchor_map`` is the bucket→elev map populated by earlier phases
    (runway corners + rect endpoints + terminal corners).

    Returns the closed-ring altitude list (``alts + [alts[0]]``) on
    success, ``None`` if the polygon is unusable.
    """
    if shape.polygon is None or shape.polygon.is_empty:
        return None
    coords = list(shape.polygon.exterior.coords)
    if not coords:
        return None
    ring_closed = (coords[0] == coords[-1])
    if ring_closed:
        coords = coords[:-1]
    if len(coords) < 3:
        return None
    dem_targets = polygon_vertex_targets(
        layout, dem, tile_lat, tile_lon, coords)
    alts, hard = _seed_alts_from_anchors(coords, anchor_map, dem_targets)
    caps = _ring_edge_caps(coords, max_grade)
    alts = _smooth_along_ring(alts, hard, caps)
    return [round(a, 1) for a in alts] + [round(alts[0], 1)]


def apply_junction_field(layout, dem, tile_lat: int, tile_lon: int,
                          anchor_map) -> int:
    """Phase 3 entry: solve every junction's node_altitudes."""
    n = 0
    for s in layout.shapes:
        if s.role != "junction":
            continue
        result = solve_polygon_field(
            layout, s, dem, tile_lat, tile_lon,
            anchor_map, JUNCTION_MAX_GRADE)
        if result is None:
            continue
        s.node_altitudes = result
        s.altitude = None
        n += 1
    return n


def apply_apron_field(layout, dem, tile_lat: int, tile_lon: int,
                       anchor_map) -> int:
    """Phase 4 entry: solve every apron's node_altitudes (1 % cap)."""
    n = 0
    for s in layout.shapes:
        if s.role != "apron":
            continue
        result = solve_polygon_field(
            layout, s, dem, tile_lat, tile_lon,
            anchor_map, APRON_MAX_GRADE)
        if result is None:
            continue
        s.node_altitudes = result
        s.altitude = None
        n += 1
    return n
