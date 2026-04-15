"""Decompose a multi-taxiway apt.dat pavement polygon into strips.

X-Plane's ``apt.dat`` format lets a single pavement row describe
arbitrary geometry, and custom scenery packs routinely abuse this
by concatenating an entire network of taxiways into ONE polygon.
At SPJC the "Taxiway V / U / Q / R / L / M" row is an 857 929 m²
blob containing six named taxiways plus their cross-connectors;
"Taxiway Aux B, C, D, E, F, G" is another 160 284 m² multi-blob.
These never fit a single minimum-rotated-rectangle well enough to
be emitted as a clean rect chain, so :mod:`O4_Taxiway_Rects` alone
falls back to apron triangulation — which is the behaviour that
leaves mega-polygons as dense triangle meshes instead of the
handful of rects the user asked for.

The fix is to **decompose** the mega-polygon into individual
strip components, rect-chain each one, and hand the junction
residue (the "hub" where multiple strips meet) back to the apron
triangulator, since a multi-axis junction legitimately needs
triangles.

Algorithm — morphological open-minus decomposition:

1. ``core = polygon.buffer(-half_width).buffer(+half_width)``.  This
   is the "erode-then-dilate" morphological opening: it keeps
   anything wider than ``2 * half_width`` and erases anything
   thinner.  In practice ``core`` is the junction hub — the wide
   regions where strips cross each other or meet a terminal apron.
2. ``branches = polygon - core``.  What's left is the "strip" part
   — the long thin arms of the polygon.  Breaking this into
   connected components via ``unary_union(...).geoms`` gives one
   polygon per strip branch.
3. **Trim each branch** by unioning back the 1 m of ``core`` that
   borders it.  Without this the branch ends abruptly at the
   junction boundary, leaving a visible elevation seam where the
   strip rect meets the hub triangles.  A 1 m expansion gives the
   rect chain a shared vertex with the hub.
4. **Strip-fit check** — each branch is re-tested for MRR fit
   ratio; branches that still fail (e.g. an L-bend inside a single
   connected region) are marked for triangulation.

Returns :class:`Decomposition` with three lists:

* ``strip_polygons`` — disjoint branches suitable for
  :func:`O4_Taxiway_Rects.build_taxiway_rects`.
* ``junction_polygons`` — the ``core`` polygons plus any
  MRR-fit-failing branches, which the caller should triangulate.
* ``used_decomposition`` — ``True`` if morphological opening found
  a non-empty ``core`` (i.e. this is a multi-strip polygon),
  ``False`` if the polygon was already a simple strip and should
  be tried directly with :func:`build_taxiway_rects`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from shapely.geometry import Polygon
from shapely.ops import unary_union

# Default half-width for the morphological opening.  Taxiways on
# SPJC range 20-45 m wide.  A 15 m half-width (30 m diameter) is
# wide enough to survive through a strip narrower than 30 m and
# narrow enough to isolate the junction hubs that are usually
# 60+ m across.
DEFAULT_HALF_WIDTH_M = 15.0

# Branches smaller than this are discarded as noise (a stub pavement
# shard at a concavity corner).
MIN_BRANCH_AREA_M2 = 50.0

# Branches are NOT dilated toward the hub.  An earlier version
# expanded each branch by 1 m so rect emission shared a vertex with
# the junction triangulation, but adjacent branches also got pushed
# into each other at the hub edge, producing cross-rect overlap.
# The MRR rect built from the branch polygon already extends
# slightly into the hub by virtue of covering the branch's full
# bounding rectangle, which is more than enough join slop — see
# O4_Taxiway_Rects.build_taxiway_rects for the rect construction.
BRANCH_JUNCTION_OVERLAP_M = 0.0


@dataclass
class Decomposition:
    """Result of :func:`decompose_multi_taxiway`."""
    strip_polygons: List[Polygon] = field(default_factory=list)
    junction_polygons: List[Polygon] = field(default_factory=list)
    used_decomposition: bool = False


def decompose_multi_taxiway(
    polygon: Polygon,
    half_width_m: float = DEFAULT_HALF_WIDTH_M,
    min_branch_area_m2: float = MIN_BRANCH_AREA_M2,
) -> Decomposition:
    """Split a taxiway polygon into strip branches + junction hubs.

    If the polygon is already a simple strip (morphological opening
    by ``half_width_m`` leaves no wide core), returns a
    :class:`Decomposition` with ``used_decomposition=False`` and the
    original polygon in ``strip_polygons``.  The caller should then
    try :func:`O4_Taxiway_Rects.build_taxiway_rects` directly on
    the whole polygon.

    For a multi-strip polygon this returns the strips and
    junction-hubs as two disjoint lists whose union covers the
    original polygon exactly.
    """
    result = Decomposition()
    if polygon is None or polygon.is_empty:
        return result

    try:
        core = polygon.buffer(-half_width_m).buffer(
            +half_width_m, join_style=2)
    except Exception:
        core = None

    if core is None or core.is_empty or core.area < min_branch_area_m2:
        # Simple strip — no wide hub.  Caller should rect-chain the
        # whole polygon directly.
        result.strip_polygons = [polygon]
        result.used_decomposition = False
        return result

    try:
        branches_geom = polygon.difference(core)
    except Exception:
        branches_geom = None
    if branches_geom is None or branches_geom.is_empty:
        # Whole polygon is the hub — nothing to rectify.  Send it
        # all to the junction triangulator.
        result.junction_polygons = [polygon]
        result.used_decomposition = True
        return result

    # Split the branches into connected components.
    if hasattr(branches_geom, "geoms"):
        branch_list = [g for g in branches_geom.geoms
                       if hasattr(g, "exterior") and not g.is_empty
                       and g.area >= min_branch_area_m2]
    else:
        branch_list = ([branches_geom]
                       if hasattr(branches_geom, "exterior")
                       and branches_geom.area >= min_branch_area_m2
                       else [])

    # Expand each branch by `BRANCH_JUNCTION_OVERLAP_M` toward the
    # hub so the rect-chain emission shares a 1 m overlap vertex
    # with the junction triangulation.  This overlap is handled by
    # the caller's emitted-parts accumulator (the apron path
    # subtracts the rects before triangulating the hub, so the
    # overlap becomes a shared edge).
    trimmed_branches: List[Polygon] = []
    for b in branch_list:
        try:
            expanded = b.buffer(BRANCH_JUNCTION_OVERLAP_M).intersection(
                polygon)
        except Exception:
            expanded = b
        if expanded.is_empty:
            continue
        if hasattr(expanded, "geoms"):
            for g in expanded.geoms:
                if (hasattr(g, "exterior") and not g.is_empty
                        and g.area >= min_branch_area_m2):
                    trimmed_branches.append(g)
        elif (hasattr(expanded, "exterior")
              and expanded.area >= min_branch_area_m2):
            trimmed_branches.append(expanded)

    # Hub = polygon minus the unioned (trimmed) branches.
    hub_polys: List[Polygon] = []
    if trimmed_branches:
        try:
            branch_union = unary_union(trimmed_branches)
            hub = polygon.difference(branch_union)
        except Exception:
            hub = core
    else:
        hub = core

    if not hub.is_empty:
        if hasattr(hub, "geoms"):
            for g in hub.geoms:
                if (hasattr(g, "exterior") and not g.is_empty
                        and g.area >= min_branch_area_m2):
                    hub_polys.append(g)
        elif hasattr(hub, "exterior") and hub.area >= min_branch_area_m2:
            hub_polys.append(hub)

    result.strip_polygons = trimmed_branches
    result.junction_polygons = hub_polys
    result.used_decomposition = True
    return result
