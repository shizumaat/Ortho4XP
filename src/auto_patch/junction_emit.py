"""Junction-polygon emission + pre-elevation geometry finalize.

Run AFTER all rect / terminal / runway shapes are emitted by
``pipeline.build_airport_pavement``.  Computes the residue
(apt.dat row-110 pavement minus rect / terminal / runway coverage),
decomposes it into junction polygons that share boundary vertices
with their rect / terminal / runway neighbours, then enforces the
"no overlap" + "shared vertex exact" invariants across the whole
layout.

Public API:
    emit_junctions_and_finalize(
        layout, *, pav_union, emitted_taxi_rects, terminal_union,
        taxi_rects, icao)

Helpers also exported (used by ``finalize.run_phase2`` for the
post-elevation reclassification pass):

    _aeroway_centerlines_m(layout, taxiway_data=None, to_m=None)
    _reclassify_stranded_junctions(layout, taxiway_data=None, to_m=None)
"""
from __future__ import annotations

import math
from typing import List, Optional, Tuple

import O4_UI_Utils as UI

from shapely.geometry import LineString, Point, Polygon
from shapely.ops import nearest_points, unary_union

from .config import (
    BOUNDARY_SAMPLE_STEP_M,
    EMIT_JUNCTIONS,
    MAX_BOUNDARY_TO_CENTERLINE_M,
)
from .elevation import (
    SHARED_VERTEX_CLUSTER_TOL_M,
    _drop_overlap_against_fixed_shapes,
)
from .layout import (
    BuiltShape,
    ROLE_APRON,
    ROLE_JUNCTION,
    ROLE_RUNWAY,
    ROLE_TERMINAL,
)
from .pavement.centerlines import _insert_points_on_boundary
from .pavement.junctions import (
    _decompose_polygon_with_holes,
    _rect_end_corners,
)
from .pavement.stubs import (
    _add_stub_to_runway_bridges,
    _clip_residue_at_stub_long_edges,
)
from .pavement.union_helpers import _merge_near_touching
from .pavement.vertices import (
    _enforce_shared_vertices,
    _validate_shared_vertex_invariant,
)


from .layout import (
    ROLE_CROSS_CONNECTOR,
    ROLE_PRIMARY_PARALLEL,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
)


__all__ = [
    "_aeroway_centerlines_m",
    "_reclassify_alongside_apron_rects",
    "_reclassify_stranded_junctions",
    "emit_junctions_and_finalize",
]


def _aeroway_centerlines_m(layout, taxiway_data=None, to_m=None):
    """Union of taxi + runway centerlines in layout-local meter
    space.  Used as the "valid junction-boundary anchor" set per
    user 2026-04-30.

    Sources:
      * Taxi rects keep ``source_axis`` (the OSM centerline span
        the rect was built from).
      * Runway shapes' long-axes are derived from the polygon's
        4-corner exterior.
      * Junctions / terminals / aprons / bridges contribute
        nothing — they're regions, not corridors.
      * When ``taxiway_data`` is supplied (Ortho4XP's per-airport
        OSM taxiway centerlines from ``extract_taxiway_info``,
        absolute lat/lon) and ``to_m`` projection is supplied,
        OSM centerlines that did NOT end up as a rect's
        ``source_axis`` are added too.  Catches sub-30 m taxiway
        stubs auto_patch dropped during rect construction; without
        this, a junction whose only nearby centerline came from
        such a stub would look "stranded" at reclassification time.

    Mirrors ``tests/test_junction_invariants.py::_aeroway_centerlines_m``
    so the production reclassifier and the test build the same
    union.
    """
    lines = []
    for s in layout.shapes:
        if s.source_axis is not None and not s.source_axis.is_empty:
            lines.append(s.source_axis)
            continue
        if (s.role == ROLE_RUNWAY
                and s.polygon is not None
                and not s.polygon.is_empty):
            coords = list(s.polygon.exterior.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            if len(coords) == 4:
                a_mid = (0.5 * (coords[0][0] + coords[3][0]),
                         0.5 * (coords[0][1] + coords[3][1]))
                b_mid = (0.5 * (coords[1][0] + coords[2][0]),
                         0.5 * (coords[1][1] + coords[2][1]))
                lines.append(LineString([a_mid, b_mid]))

    # Optional: include OSM taxiway centerlines that didn't survive
    # into a rect.  Each ``taxiway_data`` entry is
    #   {"centerline": [(lon, lat), ...], "wayid": int, "name": str}
    # in absolute coordinates; project to meter space and skip
    # any whose length is largely covered by the existing
    # rect-axis union (avoid double-counting).
    if taxiway_data and to_m is not None:
        existing_union = unary_union(lines) if lines else None
        for entry in taxiway_data:
            cl_coords = entry.get("centerline", [])
            if len(cl_coords) < 2:
                continue
            try:
                m_pts = [to_m(lon, lat) for (lon, lat) in cl_coords]
            except Exception:
                continue
            if len(m_pts) < 2:
                continue
            try:
                ls = LineString(m_pts)
            except Exception:
                continue
            if ls.is_empty or ls.length < 1.0:
                continue
            if existing_union is not None:
                try:
                    cover = ls.intersection(
                        existing_union.buffer(2.0)).length
                    if cover > 0.8 * ls.length:
                        continue
                except Exception:
                    pass
            lines.append(ls)

    if not lines:
        return None
    try:
        return unary_union(lines)
    except Exception:
        return None


def _reclassify_stranded_junctions(layout, taxiway_data=None, to_m=None):
    """Walk every ``ROLE_JUNCTION`` shape; reclassify as
    ``ROLE_APRON`` when its boundary strays farther than
    ``MAX_BOUNDARY_TO_CENTERLINE_M`` from any taxi/runway
    centerline.  Per user 2026-04-30: a valid junction's pavement
    edge is always close to a converging centerline, no matter
    how many taxiways meet there; a junction whose boundary
    strays contains apron-territory pavement that should
    triangulate as apron, not as a multi-directional junction.

    Pure metadata change — geometry untouched.  Runs BEFORE the
    overlap-clip pass so role-based priority decisions in
    ``_drop_overlap_against_fixed_shapes`` see the final
    classification.

    Returns the number of shapes reclassified.
    """
    centers = _aeroway_centerlines_m(layout, taxiway_data, to_m)
    if centers is None or centers.is_empty:
        return 0
    n_reclassified = 0
    for s in layout.shapes:
        if s.role != ROLE_JUNCTION:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        try:
            bnd = s.polygon.boundary
            L = bnd.length
        except Exception:
            continue
        if L <= 0:
            continue
        n_steps = max(2, int(L / BOUNDARY_SAMPLE_STEP_M) + 1)
        max_d = 0.0
        for i in range(n_steps):
            u = min(L, i * BOUNDARY_SAMPLE_STEP_M)
            try:
                p = bnd.interpolate(u)
                d = centers.distance(p)
            except Exception:
                continue
            if d > max_d:
                max_d = d
        if max_d > MAX_BOUNDARY_TO_CENTERLINE_M:
            s.role = ROLE_APRON
            n_reclassified += 1
    return n_reclassified


def _reclassify_alongside_apron_rects(layout):
    """Reclassify sloping rects that are alongside apron-class
    pavement post-emit.

    Per the user 2026-04-30 absorption rule: a sloping rect
    (primary_parallel / secondary_parallel / stub /
    cross_connector) cannot share a long edge with an apron /
    junction polygon, because the rect's straight-line slope
    along the long edge has to match the surrounding pavement's
    natural DEM slope along the seam — which it generally
    won't, producing visible elevation glitches.

    Probe semantics (5 m axial steps, 5 m outside-long-edge
    probe distance, EITHER long edge, ≥ 10 % adjacent run)
    match ``_drop_primary_parallels_embedded_in_pavement``
    exactly — the rule's authoritative semantics, frozen per
    ``feedback_shape_rules.md``.

    The difference vs the legacy Stage-1 absorption: the input
    set is the ACTUAL post-emit ``apron + junction`` shape pool
    rather than the pre-emit ``apt_pav_union − full_runway −
    all_input_taxis`` snapshot.  This eliminates both
    over-absorption (Stage-1 drops based on a too-shrunken
    junction_pav, removing rects that should have survived)
    and under-absorption (Stage-1 misses adjacency the test
    flags because the real apron pool is larger
    post-elevation, post-Phase-B.1 reclass).

    Reclassifies to ROLE_APRON: aligns with the rule's wording
    ("absorbed into the surrounding apron") and keeps the
    polygon's geometry intact (no merge).  Iterated to a fixed
    point: each rect reclassified expands the apron pool, which
    may flag adjacent rects.

    Returns the number of rects reclassified.
    """
    sloping_rect_roles = {
        ROLE_PRIMARY_PARALLEL,
        ROLE_SECONDARY_PARALLEL,
        ROLE_STUB,
        ROLE_CROSS_CONNECTOR,
    }
    SAMPLE_STEP_M = 5.0
    OUTER_PROBE_M = 5.0
    ADJACENCY_FRAC = 0.10
    MIN_AXIS_M = 30.0
    MAX_ITERATIONS = 8

    total_reclassified = 0
    for _it in range(MAX_ITERATIONS):
        # Build apron-class union from the CURRENT shape pool
        # (each iteration sees previous reclassifications).
        other_polys = [s.polygon for s in layout.shapes
                        if s.role in {ROLE_APRON, ROLE_JUNCTION}
                        and s.polygon is not None
                        and not s.polygon.is_empty]
        if not other_polys:
            return total_reclassified
        try:
            other_union = unary_union(other_polys)
        except Exception:
            return total_reclassified

        n_this_iter = 0
        for s in layout.shapes:
            if s.role not in sloping_rect_roles:
                continue
            if s.polygon is None or s.polygon.is_empty:
                continue
            try:
                rc = list(s.polygon.exterior.coords)
            except Exception:
                continue
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            if len(rc) != 4:
                continue
            a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                     0.5 * (rc[0][1] + rc[3][1]))
            b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                     0.5 * (rc[1][1] + rc[2][1]))
            ax = b_mid[0] - a_mid[0]
            ay = b_mid[1] - a_mid[1]
            L = math.hypot(ax, ay)
            if L < MIN_AXIS_M:
                continue
            ux, uy = ax / L, ay / L
            nx, ny = -uy, ux
            half_w = math.hypot(rc[0][0] - a_mid[0],
                                rc[0][1] - a_mid[1])
            if half_w < 1.0:
                continue
            outer = half_w + OUTER_PROBE_M
            # Subtract this rect from the test pool defensively
            # (apron polygons should never include a rect's
            # footprint, but float-noise overlaps could exist).
            try:
                test_pav = other_union.difference(s.polygon)
            except Exception:
                test_pav = other_union
            n_steps = max(2, int(L / SAMPLE_STEP_M) + 1)
            adj = 0
            for i in range(n_steps):
                u = min(L, i * SAMPLE_STEP_M)
                cx = a_mid[0] + u * ux
                cy = a_mid[1] + u * uy
                try:
                    lp = Point(cx + nx * outer,
                                cy + ny * outer)
                    rp = Point(cx - nx * outer,
                                cy - ny * outer)
                    if (test_pav.contains(lp)
                            or test_pav.contains(rp)):
                        adj += 1
                except Exception:
                    continue
            if adj / max(1, n_steps) >= ADJACENCY_FRAC:
                s.role = ROLE_APRON
                n_this_iter += 1

        total_reclassified += n_this_iter
        if n_this_iter == 0:
            break

    return total_reclassified


def emit_junctions_and_finalize(layout, *, pav_union, emitted_taxi_rects,
                                terminal_union, taxi_rects, icao):
    """Emit junction polygons + run final geometry repair passes.

    Mutates layout in place.
    """
    # ── Junction emission (user 2026-04-23): junctions and
    # aprons are treated identically going forward — both will
    # triangulate and slope multi-directionally (unlike rects,
    # which slope only along their axis).  No distinction needed;
    # emit every connected non-rect non-terminal pavement region
    # as a SINGLE junction polygon, with rect + terminal corners
    # injected as exact shared boundary vertices.  Seamless
    # coverage follows by construction.
    MIN_JUNCTION_AREA_M2 = 50.0   # drop only sliver noise from
                                   # rect-snap inexactness
    SIMPLIFY_TOL_M = 2.0          # min-spacing simplification —
                                   # drops apt.dat-curve vertices
                                   # closer than this to their
                                   # neighbours.  Larger ⇒ fewer
                                   # triangles, fewer slivers; rect /
                                   # terminal seam corners survive
                                   # because they're sharp 90°
                                   # turns that DP keeps.
    RECT_CORNER_TOL_M = 2.0       # corner-to-boundary injection range

    taxi_rect_union = (unary_union(emitted_taxi_rects)
                       if emitted_taxi_rects else None)

    # ── Junction polygons: every pavement region not covered by a
    # rect or a terminal.  User 2026-04-23: simplest polygon that
    # covers the remaining pavement and connects to all rects /
    # terminals; each one will triangulate + slope
    # multi-directionally at elevation time.
    if pav_union is not None:
        residue = pav_union
        if taxi_rect_union is not None:
            residue = residue.difference(taxi_rect_union)
        if terminal_union is not None and not terminal_union.is_empty:
            residue = residue.difference(terminal_union)
        # Defensive: subtract runway even though pav_union already
        # had it removed — floating-point boundary artifacts can
        # leave sub-meter residue slivers overlapping runway.  Use
        # the EFFECTIVE runway union (i.e. with apron-merged
        # regions excluded) so the residue covers parts of runways
        # that pass through aprons.
        _eff_rwy = getattr(layout, "_effective_runway_union",
                           layout.runway_union)
        if _eff_rwy is not None and not _eff_rwy.is_empty:
            residue = residue.difference(_eff_rwy)

        # Per user 2026-04-27 invariant: NO polygon along the long
        # edge of a sloping rect.  Even if apt.dat has pavement
        # extending past a stub's long edge (because the boundary
        # bulges between the stub's two short edges), the residue
        # polygon there must NOT become a junction — it would wrap
        # around the stub's long edge.  Subtract a thin strip just
        # OUTSIDE each stub's long edges from the residue so the
        # resulting junctions stop at the stub's short-edge corners.
        try:
            residue = _clip_residue_at_stub_long_edges(
                residue, taxi_rects)
        except Exception:
            pass

        # Per user 2026-04-27 exception: when a stub's runway-facing
        # short edge has a pavement GAP to the runway (apt.dat
        # construction-era discrepancy), project a synthetic
        # quadrilateral from the stub's short edge straight to the
        # runway boundary so the connecting junction has continuous
        # coverage.
        try:
            residue = _add_stub_to_runway_bridges(
                residue, taxi_rects, layout.runway_union)
        except Exception:
            pass

        # Merge near-touching residue parts (post-subtraction
        # MultiPolygon split is often a numerical artifact of
        # apt.dat polygon-boundary precision rather than a real
        # disconnect — see ``_merge_near_touching``).
        residue = _merge_near_touching(residue)
        parts = ([residue] if residue.geom_type == "Polygon"
                 else list(getattr(residue, "geoms", [])))

        # Collect rect corners and terminal corners — both get
        # injected as shared boundary vertices so every junction
        # seams cleanly to its rect/terminal neighbours.
        seam_points: List[Tuple[float, float]] = []
        for rect, axis, role, ref in taxi_rects:
            rc = list(rect.exterior.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            seam_points.extend(rc)
        for shape in layout.shapes:
            if shape.role != ROLE_TERMINAL:
                continue
            tc = list(shape.polygon.exterior.coords)
            if tc and tc[0] == tc[-1]:
                tc = tc[:-1]
            seam_points.extend(tc)
        # Also the runway corners + axis-aligned runway vertices
        # so stub-to-runway edges are shared.
        if layout.runway_union is not None and not layout.runway_union.is_empty:
            ru = layout.runway_union
            geoms = [ru] if ru.geom_type == "Polygon" else list(
                getattr(ru, "geoms", []))
            for g in geoms:
                if g.geom_type != "Polygon":
                    continue
                rc = list(g.exterior.coords)
                if rc and rc[0] == rc[-1]:
                    rc = rc[:-1]
                seam_points.extend(rc)

        # Dedup seam points that coincide within 0.5 m.
        uniq: List[Tuple[float, float]] = []
        for c in seam_points:
            if not any(math.hypot(c[0] - u[0], c[1] - u[1]) < 0.5
                       for u in uniq):
                uniq.append(c)
        seam_points = uniq

        for part in parts:
            if part.geom_type != "Polygon":
                continue
            if part.area < MIN_JUNCTION_AREA_M2:
                continue
            # Inject seam points as exact boundary vertices.
            part = _insert_points_on_boundary(
                part, seam_points, tol=RECT_CORNER_TOL_M)
            # Light simplification of the remaining apt.dat
            # sub-vertex noise.  Tol kept small so injected seam
            # vertices aren't dropped as near-colinear.
            try:
                simp = part.simplify(SIMPLIFY_TOL_M,
                                     preserve_topology=True)
            except Exception:
                simp = part
            if (simp.is_empty
                    or simp.geom_type != "Polygon"
                    or simp.area < MIN_JUNCTION_AREA_M2):
                continue
            # Re-inject seam points in case simplify dropped any
            # (Douglas-Peucker will remove a vertex on an almost-
            # straight segment even when we need it as a seam).
            simp = _insert_points_on_boundary(
                simp, seam_points, tol=RECT_CORNER_TOL_M)
            if EMIT_JUNCTIONS:
                # If the residue polygon wraps a fully-enclosed
                # rect (holes in the polygon), decompose it into
                # multiple simple polygons that join AROUND the
                # rect instead.  Cleaner human-editable output and
                # avoids hole-splicing artefacts in triangulation.
                pieces = _decompose_polygon_with_holes(
                    simp, min_area_m2=MIN_JUNCTION_AREA_M2)
                for piece in pieces:
                    # Re-inject seam points: cut lines added by
                    # the decomposition split may have introduced
                    # new boundary vertices that DON'T align with
                    # any rect/terminal corner; we still need
                    # corners that lie on this piece's boundary.
                    piece = _insert_points_on_boundary(
                        piece, seam_points,
                        tol=RECT_CORNER_TOL_M)
                    if (piece.is_empty
                            or piece.geom_type != "Polygon"
                            or piece.area < MIN_JUNCTION_AREA_M2):
                        continue
                    layout.shapes.append(BuiltShape(
                        polygon=piece, role=ROLE_JUNCTION))

    # ── Runway-taxiway shared-vertex sync ──
    # Stubs that widen into the runway apron (V1-style) need the
    # runway polygon to carry the projection of their outer
    # corners as vertices, so the junction that wraps around the
    # stub has exact vertex coincidence with the runway.
    rwy_vertex_inserts: List[Tuple[float, float]] = []
    if layout.runway_union is not None and not layout.runway_union.is_empty:
        rwy_boundary = layout.runway_union.boundary
        for rect, axis, role, ref in taxi_rects:
            coords_ax = list(axis.coords)
            if len(coords_ax) < 2:
                continue
            for ax_pt in (coords_ax[0], coords_ax[-1]):
                d_rwy = Point(ax_pt).distance(rwy_boundary)
                if d_rwy > 60.0:
                    continue
                pairs = _rect_end_corners(rect, axis)
                if len(pairs) < 2:
                    continue
                end_idx = 0 if ax_pt == coords_ax[0] else 1
                c1, c2 = pairs[end_idx]
                try:
                    rp1 = nearest_points(rwy_boundary, Point(c1))[0]
                    rp2 = nearest_points(rwy_boundary, Point(c2))[0]
                except Exception:
                    continue
                rwy_vertex_inserts.append((rp1.x, rp1.y))
                rwy_vertex_inserts.append((rp2.x, rp2.y))

    if rwy_vertex_inserts:
        for shape in layout.shapes:
            if shape.role != ROLE_RUNWAY:
                continue
            shape.polygon = _insert_points_on_boundary(
                shape.polygon, rwy_vertex_inserts, tol=2.0)

    # ── Global shared-vertex enforcement (user rule 16) ─────────
    # Cluster all emitted-shape vertices within SHARED_VERTEX_TOL_M
    # and replace each with the cluster centroid.  This guarantees
    # adjacent shapes have EXACT vertex coincidence, which target
    # files enforce and compare_target's v_tgt metric measures.
    _enforce_shared_vertices(layout, tol=SHARED_VERTEX_CLUSTER_TOL_M)

    # Pavement layout invariant (user 2026-04-26): NO shape may
    # overlap another, period.  Runs AFTER shared-vertex collapse
    # because the collapse step can shift boundaries by up to
    # SHARED_VERTEX_CLUSTER_TOL_M / 2 and create tiny new overlaps.
    # Re-run shared-vertex collapse afterwards because clipping
    # introduces new vertices at intersection points that may sit
    # within the cluster tol of an existing vertex on the
    # higher-priority shape's edge.
    _drop_overlap_against_fixed_shapes(layout, icao=icao)
    _enforce_shared_vertices(layout, tol=SHARED_VERTEX_CLUSTER_TOL_M)

    # Validate the invariant: every vertex of every shape must either
    # be unique (distance > 2 × tol to any other shape's vertex) OR
    # exactly equal to a vertex on an adjacent shape.  No "close but
    # not equal" drift is permitted.
    _validate_shared_vertex_invariant(layout,
                                      tol=SHARED_VERTEX_CLUSTER_TOL_M)
