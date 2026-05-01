"""Polygon-union helpers shared across the pavement pipeline.

Single-purpose: bridge sub-meter precision gaps between adjacent
apt.dat row-110 polygons that ``shapely.unary_union`` leaves
disjoint.  See ``_merge_near_touching`` for the rule.

Public API:
    _merge_near_touching(geom)
    PAVEMENT_BRIDGE_GAP_M
"""
from __future__ import annotations

from typing import Optional

from shapely.geometry import Polygon


__all__ = [
    "PAVEMENT_BRIDGE_GAP_M",
    "_merge_near_touching",
]


# Tolerance for bridging numerical / sub-meter gaps between near-
# touching apt.dat polygons.  apt.dat at busy airports often stores
# adjacent apron areas as separate row-110 polygons whose shared
# edges have sub-millimeter-to-cm floating-point differences (e.g.
# at SPJC's SE apron a 6.5 m thin "strip" appears between two large
# apron polygons because their shared boundary y-values differ by
# 5 cm).  unary_union doesn't merge these because they don't
# overlap — but for our purposes they ARE one continuous coverage.
# Closing 0.5 m gaps via buffer-shrink merges them while preserving
# real holes (typically meters-wide non-pavement islands).
PAVEMENT_BRIDGE_GAP_M = 0.1


def _merge_near_touching(geom: Optional[Polygon],
                         eps: float = PAVEMENT_BRIDGE_GAP_M
                         ) -> Optional[Polygon]:
    """Force-merge near-touching components of ``geom`` (a possibly
    MultiPolygon) by a buffer-then-shrink.  Returns the same kind
    of geometry (Polygon if single, MultiPolygon if truly disjoint).

    Per user 2026-04-27: we want a single pavement union with holes
    only — apt.dat polygon boundaries shouldn't survive into the
    output.  This helper bridges sub-meter precision gaps between
    apt.dat polygons that ``unary_union`` leaves disjoint.
    """
    if geom is None or geom.is_empty:
        return geom
    try:
        merged = geom.buffer(eps, join_style=2).buffer(
            -eps, join_style=2)
        if merged.is_empty:
            return geom
        if merged.geom_type not in ("Polygon", "MultiPolygon"):
            return geom
        return merged
    except Exception:
        return geom
