"""Canonical-point registry for shared-vertex management.

The pavement builder constructs many shapes (rects, junctions, aprons,
boundary pieces) whose perimeters meet at intersection points.  Each
of those meeting points should be a SINGLE geometric vertex shared
by every adjacent shape — at exact floating-point equality, so that
``pav_union.difference(rects)`` and the OSM emitter's vertex
bucketing produce one node ID per real-world point.

Without a shared registry, each shape's corner is snapped
independently to ``pav.boundary`` and ends up at a slightly
different location than its neighbour's "same" corner.  The
sub-millimetre to multi-metre drift cascades downstream:
``buffer(0)`` validity repairs, sliver-corner removal, T-junction
splits, and merge-corner-junctions all exist as workarounds for
this single root cause.

A canonical-point registry replaces independent snapping with a
deterministic ``get_or_add`` lookup.  Every shape that needs a
corner near (x, y) queries the registry; if any prior shape has
already registered a canonical point within
``SHARED_VERTEX_TOL_M``, the same coordinates are returned.
Otherwise (x, y) becomes the new canonical point for that bucket.

Seeded with the input fixed-geometry vertices (apt.dat row-110
pavement polygon vertices + runway corners) so the canonical
set is anchored to real apt.dat data rather than to whatever
rect happened to register a point first.
"""
from __future__ import annotations

import math
from typing import List, Optional, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import Polygon


__all__ = ["CanonicalPointRegistry", "snap_polygon_through_registry"]


_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)


class CanonicalPointRegistry:
    """Spatial-index registry of canonical (x, y) points.

    All shape corners that should be SHARED must go through
    ``get_or_add`` so they pick up the same exact coordinates.

    ``tol_m`` is the bucket radius — points within ``tol_m`` of an
    existing entry resolve to that entry.  Default matches
    ``layout.SHARED_VERTEX_TOL_M`` so registry sharing aligns with
    the OSM emitter's vertex bucketing.
    """

    def __init__(self, tol_m: float = 0.5):
        self.tol_m = tol_m
        # Cell size = tol so neighbours-of-neighbours covers the
        # full lookup radius.
        self._cell = max(tol_m, 0.1)
        self._points: List[Tuple[float, float]] = []
        # cell key (ix, iy) → list of indices into self._points
        self._index: dict = {}

    # ── public API ────────────────────────────────────────────────

    def seed(self, points) -> int:
        """Bulk-add anchor points (apt.dat row-110 + runway corners).

        Duplicates within ``tol_m`` collapse to a single entry.
        Returns the number of NEW canonical points added.
        """
        before = len(self._points)
        for p in points:
            self.get_or_add(float(p[0]), float(p[1]))
        return len(self._points) - before

    def get_or_add(self, x: float, y: float
                    ) -> Tuple[float, float]:
        """Return the canonical (x, y) for the given input point.

        If an existing canonical point sits within ``tol_m`` of
        (x, y), return its coordinates exactly.  Otherwise insert
        (x, y) as a new canonical point and return it.
        """
        nearby = self._find_nearest(x, y, self.tol_m)
        if nearby is not None:
            return nearby
        return self._add(x, y)

    def find_nearest(self, x: float, y: float,
                      max_d: float) -> Optional[Tuple[float, float]]:
        """Find the nearest canonical point within ``max_d`` of
        (x, y).  Does NOT add.  Returns None if no entry qualifies.
        """
        return self._find_nearest(x, y, max_d)

    @property
    def size(self) -> int:
        return len(self._points)

    def points(self) -> List[Tuple[float, float]]:
        """Return a snapshot of every canonical point (insertion
        order).  Caller-owned list."""
        return list(self._points)

    # ── internals ─────────────────────────────────────────────────

    def _cell_key(self, x: float, y: float) -> Tuple[int, int]:
        return (int(math.floor(x / self._cell)),
                int(math.floor(y / self._cell)))

    def _add(self, x: float, y: float) -> Tuple[float, float]:
        idx = len(self._points)
        coords = (float(x), float(y))
        self._points.append(coords)
        self._index.setdefault(
            self._cell_key(x, y), []).append(idx)
        return coords

    def _find_nearest(self, x: float, y: float,
                       max_d: float) -> Optional[Tuple[float, float]]:
        # Number of cells in each direction we have to scan to cover
        # the lookup radius.  +1 to be safe at cell boundaries.
        n_cells = int(math.ceil(max_d / self._cell)) + 1
        cx, cy = self._cell_key(x, y)
        best: Optional[Tuple[float, float]] = None
        best_d = max_d
        for dx in range(-n_cells, n_cells + 1):
            for dy in range(-n_cells, n_cells + 1):
                key = (cx + dx, cy + dy)
                bucket = self._index.get(key)
                if not bucket:
                    continue
                for idx in bucket:
                    px, py = self._points[idx]
                    d = math.hypot(px - x, py - y)
                    if d < best_d:
                        best_d = d
                        best = (px, py)
        return best


def snap_polygon_through_registry(
        poly: Optional[Polygon],
        registry: Optional[CanonicalPointRegistry],
) -> Optional[Polygon]:
    """Route every vertex of ``poly``'s exterior + interior rings
    through ``registry.get_or_add`` so drift introduced by
    ``buffer(0)`` / ``unary_union`` / ``simplify`` resolves to
    canonical (x, y) coordinates shared with adjacent shapes.

    Returns the snapped polygon (a new ``Polygon`` if any vertex
    moved, the input otherwise), or ``None`` if the snap produces
    a degenerate shape.  If either input is ``None`` / empty, the
    input is returned unchanged.
    """
    if registry is None or poly is None or poly.is_empty:
        return poly

    def _snap_ring(coords):
        snapped = []
        for x, y in coords:
            cp = registry.get_or_add(float(x), float(y))
            snapped.append(cp)
        return snapped

    try:
        ext = _snap_ring(poly.exterior.coords)
        if len(set(ext)) < 3:
            return None
        interiors = []
        for ring in poly.interiors:
            ri = _snap_ring(ring.coords)
            if len(set(ri)) < 3:
                continue
            interiors.append(ri)
        snapped_poly = Polygon(ext, interiors)
        if not snapped_poly.is_valid:
            snapped_poly = snapped_poly.buffer(0)
            if (snapped_poly.is_empty
                    or snapped_poly.geom_type not in (
                        "Polygon", "MultiPolygon")):
                return None
            if snapped_poly.geom_type == "MultiPolygon":
                snapped_poly = max(
                    snapped_poly.geoms, key=lambda g: g.area)
        return snapped_poly
    except _GEOM_EXC:
        return None
