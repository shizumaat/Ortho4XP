"""DEM sampling helpers for the per-surface solver.

Thin wrappers around ``auto_patch.elevation._sample_dem`` that return
target elevations at semantically meaningful points: rect short-end
midpoints, rect axial midpoint, polygon vertex.
"""
from __future__ import annotations

import math
from typing import List, Optional, Tuple

from shapely.geometry import LineString


def sample_at_xy(dem, tile_lat: int, tile_lon: int,
                 layout, x: float, y: float) -> Optional[float]:
    """Sample DEM at meter-space (x, y).  Converts to lat/lon via
    the layout's anchor and delegates to ``_sample_dem``.
    """
    from auto_patch.elevation import _sample_dem
    lat, lon = layout.m_to_ll(x, y)
    e = _sample_dem(dem, tile_lat, tile_lon, lat, lon)
    return float(e) if e is not None else None


def rect_axial_targets(layout, dem, tile_lat: int, tile_lon: int,
                       source_axis: LineString
                       ) -> Optional[Tuple[float, float, float]]:
    """DEM samples at axial fractions 0, 0.5, 1 along ``source_axis``.

    Returns ``(start, mid, end)`` triple in meters or ``None`` if the
    axis is unusable or any sample is missing.  Per user 2026-05-02:
    three points (endpoints + midpoint) are sufficient — finer
    sampling would just amplify DEM noise.
    """
    if source_axis is None or source_axis.is_empty:
        return None
    pts = list(source_axis.coords)
    if len(pts) < 2:
        return None
    ax_start = pts[0]
    ax_end = pts[-1]
    mid = ((ax_start[0] + ax_end[0]) / 2.0,
           (ax_start[1] + ax_end[1]) / 2.0)
    e_start = sample_at_xy(dem, tile_lat, tile_lon, layout, *ax_start)
    e_mid = sample_at_xy(dem, tile_lat, tile_lon, layout, *mid)
    e_end = sample_at_xy(dem, tile_lat, tile_lon, layout, *ax_end)
    if e_start is None or e_mid is None or e_end is None:
        return None
    return e_start, e_mid, e_end


def polygon_vertex_targets(layout, dem, tile_lat: int, tile_lon: int,
                           coords) -> List[Optional[float]]:
    """Per-vertex DEM samples for a polygon's exterior ring."""
    return [sample_at_xy(dem, tile_lat, tile_lon, layout, x, y)
            for x, y in coords]
