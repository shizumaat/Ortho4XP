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
    "_simplify_pavement_polygon",
    # Backwards-compat alias.
    "_drop_close_nonadjacent_pairs",
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


def _simplify_pavement_polygon(geom, tol: float = 1.0):
    """Simplify a pavement polygon: drop sub-``tol`` detail and snip
    sliver-tip corners.

    Per user 2026-05-05: pavement should never have nodes closer than
    1 m, and X-Plane's mesh builder crashes on sub-``SLIVER_ANGLE_
    THRESHOLD_DEG`` (2°) corner spikes.  Apt.dat polygons routinely
    contain both: over-resolved curves stored as 100s of sub-meter
    steps, and 1° needle-tip features that look like real pavement
    tabs but were just floating-point doubled vertices in the
    source data.

    Two passes:

    1. **Douglas-Peucker simplify** (shapely's ``.simplify(tol,
       preserve_topology=True)``).  Drops verts whose perpendicular
       distance to the simplified edge is < ``tol``.  Eliminates
       over-resolution and most close-pair noise.  Preserves
       polygon topology.

    2. **Sliver-tip removal.**  After simplify, any remaining
       vertex with interior angle < ``SLIVER_ANGLE_THRESHOLD_DEG``
       is snipped; the two flanking vertices are joined directly,
       collapsing the needle into a chord.  Iterates so a freshly
       exposed sliver after one drop gets caught on the next pass.

    Returns the simplified polygon (Polygon or MultiPolygon, same
    type as input).  Falls back to the input on any failure.
    """
    if geom is None or geom.is_empty:
        return geom
    if geom.geom_type == "MultiPolygon":
        return type(geom)([_simplify_pavement_polygon(g, tol)
                            for g in geom.geoms])
    if geom.geom_type != "Polygon":
        return geom
    try:
        # Pass 1: DP simplify.
        simp = geom.simplify(tol, preserve_topology=True)
        if (simp.is_empty
                or simp.geom_type not in ("Polygon", "MultiPolygon")):
            return geom
        if simp.geom_type == "MultiPolygon":
            # Topology preservation can split the polygon if a
            # narrow neck collapses; keep the largest piece.
            simp = max(simp.geoms, key=lambda g: g.area)
        # Pass 2: drop sliver-tip corners (re-imported here to avoid
        # circular imports at module load time).
        from .junctions import _drop_sliver_corners
        from shapely.geometry import Polygon as _P

        def _clean_ring(ring_coords):
            coords = list(ring_coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            cleaned = _drop_sliver_corners(coords)
            if len(cleaned) < 3:
                return None
            cleaned.append(cleaned[0])
            return cleaned

        ext = _clean_ring(simp.exterior.coords)
        if ext is None:
            return geom
        ints = []
        for ring in simp.interiors:
            cr = _clean_ring(ring.coords)
            if cr is not None and len(cr) >= 4:
                ints.append(cr)
        out = _P(ext, ints)
        if not out.is_valid:
            out = out.buffer(0)
            if out.geom_type == "MultiPolygon":
                out = max(out.geoms, key=lambda g: g.area)
        if (out.is_valid
                and not out.is_empty
                and out.geom_type == "Polygon"):
            return out
    except Exception:
        pass
    return geom


# Backwards-compatibility alias.  Previous incarnations of this
# helper had narrower behaviour (just non-adjacent close-pair
# removal); the new function does that and more.  Existing call
# sites can use either name.
_drop_close_nonadjacent_pairs = _simplify_pavement_polygon


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
