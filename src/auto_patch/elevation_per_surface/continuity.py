"""Phase 5 + shared-vertex continuity reconciliation.

Phase 5 sets terminal altitudes (flat at DEM-median, capped by
the existing terminal-altitude rule); the existing
``elevation.py`` logic already does this and we keep it for now.

Continuity reconciliation: after all phases run, the same physical
vertex shared by two shapes must hold one value.  Resolution
priority (highest wins):

    runway > rect > junction > apron > terminal

Higher-priority shapes' values overwrite lower-priority ones at
shared buckets.  Any conflicts that can't be resolved by ordering
(rare — typically a rect+rect shared vertex) take the average.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

from .bfs_propagate import _bucket

ROLE_PRIORITY = {
    "runway": 100,
    "primary_parallel": 90,
    "secondary_parallel": 90,
    "stub": 90,
    "cross_connector": 90,
    "junction": 70,
    "apron": 60,
    "terminal": 50,
}


def _shape_alt_at_corner(shape, corner_idx: int,
                          coord_xy: Tuple[float, float]
                          ) -> Optional[float]:
    """Return the elevation this shape assigns to a given polygon
    corner.  For sloped rects we use the rect's altitude_high or
    altitude_low depending on the corner's source_axis projection.
    """
    if shape.altitude is not None:
        return float(shape.altitude)
    if shape.node_altitudes:
        if 0 <= corner_idx < len(shape.node_altitudes):
            return float(shape.node_altitudes[corner_idx])
    if shape.altitude_high is not None and shape.altitude_low is not None:
        from auto_patch.elevation import _short_end_pairs_by_axis
        coords = list(shape.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        if len(coords) == 4 and shape.source_axis is not None:
            sp, ep = _short_end_pairs_by_axis(coords, shape.source_axis)
            if sp is not None and corner_idx in sp:
                # source_axis t=0 end — convention says altitude_high
                # (or altitude_low; without direction we use the
                # average of high and low).
                return (float(shape.altitude_high)
                        + float(shape.altitude_low)) / 2.0
            if ep is not None and corner_idx in ep:
                return (float(shape.altitude_high)
                        + float(shape.altitude_low)) / 2.0
    return None


def reconcile_shared_vertices(layout) -> int:
    """Resolve shared-vertex disagreements.

    Pass 1: walk every shape's polygon corners and record (priority,
    elevation) per bucket.  Highest-priority elevation wins.

    Pass 2: walk again and update each shape's stored elevations to
    match the winning bucket value.

    Returns the number of buckets reconciled.
    """
    bucket_winner: Dict[Tuple[int, int],
                        Tuple[int, float]] = {}
    for s in layout.shapes:
        prio = ROLE_PRIORITY.get(s.role)
        if prio is None:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for i, (x, y) in enumerate(coords):
            a = _shape_alt_at_corner(s, i, (x, y))
            if a is None:
                continue
            b = _bucket(x, y)
            cur = bucket_winner.get(b)
            if cur is None or prio > cur[0]:
                bucket_winner[b] = (prio, a)

    n_changed = 0
    for s in layout.shapes:
        if s.role == "runway":
            # Runways are HARD; never overwrite.
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        if s.role != "junction":
            # Only junctions track per-vertex altitudes; rects /
            # aprons / terminals use scalar fields and we leave them
            # alone (the writeback already used reconciled inputs).
            continue
        if not s.node_altitudes:
            continue
        coords = list(s.polygon.exterior.coords)
        ring_closed = (coords and coords[0] == coords[-1])
        coords_open = coords[:-1] if ring_closed else coords
        new_alts = list(s.node_altitudes)
        for i, (x, y) in enumerate(coords_open):
            b = _bucket(x, y)
            w = bucket_winner.get(b)
            if w is None:
                continue
            if i < len(new_alts) and abs(new_alts[i] - w[1]) > 0.05:
                new_alts[i] = w[1]
                n_changed += 1
        if ring_closed and len(new_alts) > 0:
            new_alts[-1] = new_alts[0]
        s.node_altitudes = [round(a, 1) for a in new_alts]
    return n_changed


def collect_anchor_map(layout) -> Dict[Tuple[int, int], float]:
    """Build a bucket→elev anchor map from currently-set values on
    runway segments + rects + terminals.  Used by junction / apron
    phase to seed boundary anchors.
    """
    anchor_map: Dict[Tuple[int, int], float] = {}
    for s in layout.shapes:
        if s.polygon is None or s.polygon.is_empty:
            continue
        coords = list(s.polygon.exterior.coords)
        if coords and coords[0] == coords[-1]:
            coords = coords[:-1]
        for i, (x, y) in enumerate(coords):
            a = _shape_alt_at_corner(s, i, (x, y))
            if a is not None:
                anchor_map.setdefault(_bucket(x, y), a)
    return anchor_map
