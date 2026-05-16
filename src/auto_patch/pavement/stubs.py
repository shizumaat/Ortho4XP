"""Runway-end primary-parallel stub emission.

For each runway end, projects nearby parallel taxi centerlines
into a short connector ("stub") that bridges the parallel to the
runway threshold area.  These appear at most airports as the
A1/F1-style turn-off taxiways near each runway end.

Public API (leading-underscore preserved for backward compatibility
with internal callers in ``O4_Airport_Pavement_Builder``):

    _emit_primary_parallel_runway_stubs
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from shapely.errors import GEOSException, TopologicalError
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from ..layout import ROLE_STUB
from .rects import (
    _extend_rect_corners_perpendicular,
    _natural_half_width,
    _rect_from_axis_extended,
)

# Narrow exception tuple for shapely / numeric-geometry failure
# modes.  Programming errors propagate so they surface immediately.
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)


__all__ = [
    "_add_stub_to_runway_bridges",
    "_clip_residue_at_stub_sloping_edges",
    "_emit_primary_parallel_runway_stubs",
]


def _emit_primary_parallel_runway_stubs(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Tuple[str, List[str], Dict[str, str]]],
    to_m,
    runway_union: Optional[Polygon],
    pav_union: Optional[Polygon],
    apt_vertices: Optional[List[Tuple[float, float]]],
    existing_taxi_rects: List[Tuple[Polygon, LineString, str, str]],
    apt_centerlines: Optional[
        List[Tuple[LineString, str]]] = None,
    rwy_centerlines: Optional[List[LineString]] = None,
) -> List[Tuple[Polygon, LineString, str, str]]:
    """Emit an extra STUB rect at each primary parallel OSM path
    endpoint that terminates INSIDE the runway polygon.

    A / F / L OSM primary taxis at SPJC extend their polyline
    onto the runway pavement itself (the endpoint vertex sits
    inside the runway polygon).  The target hand-drawn OSM has
    a short wide STUB at that transition — the RAMP where the
    taxi meets the runway short-edge.  Detection:

      1. Merge each primary-parallel ref's OSM ways
         (linemerge + gap-bridge, same as the main extraction).
      2. Check each merged polyline's 2 endpoints.
         If an endpoint is INSIDE the runway polygon (or within
         10 m of it), the taxi terminates on runway pavement.
      3. Walk along the polyline from that endpoint toward the
         interior, until the path vertex distance to runway
         boundary exceeds ``STUB_EXIT_D_M`` (~80 m).  The
         "exit" vertex is where the path leaves the runway-
         apron ramp and enters normal taxi corridor.
      4. Emit a STUB rect CENTERED on the exit vertex along the
         local path direction, length ``STUB_LEN_M`` (~80 m).
         The rect width is the perpendicular pav half-width
         at the exit point × 2 (full pav width — the ramp is
         wider than a normal taxi).

    Returns a list of extra (rect, axis, role, ref) tuples to
    append to the main ``taxi_rects`` list.
    """
    if runway_union is None or runway_union.is_empty:
        return []
    if pav_union is None or pav_union.is_empty:
        return []

    # SPJC A/F have loop ramps at runway ends — target stubs sit at
    # the APEX of the loop where the ramp meets the normal taxi
    # corridor (d_rwy ≈ 100 m at a path vertex).  Vertex-based
    # exit at threshold 80 m lands on that apex vertex.  SPLP has
    # no loops — the primary curves smoothly into the runway and
    # target stubs sit MID-CURVE BETWEEN vertices at d_rwy ≈ 75 m.
    # Vertex-based exit over-shoots or under-shoots; interpolate
    # to the exact target d-value.  Separate thresholds:
    STUB_EXIT_D_M = 80.0
    STUB_INTERP_TARGET_D_UNREFED = 75.0
    STUB_LEN_M = 80.0       # target A/F stubs are 79–93 m
    ENDPOINT_INSIDE_TOL_M = 10.0  # allow near-boundary endpoints
    OUTSIDE_NEAR_RWY_M = 135.0  # NE-end endpoint within 135 m of
                                # runway but not inside (SPLP main
                                # taxi ends at 127 m from runway).
                                # SPJC L's internal endpoints sit
                                # at 148-150 m so 135 m excludes
                                # them while catching SPLP's
                                # runway-facing taxi end.
    UNREFED_MIN_LEN_M = 800.0   # unrefed ways that act as long
                                # primary parallels (SPLP main taxi
                                # is 2640 m unrefed)

    # Gather per-ref centerlines.  Prefer apt.dat taxi-network
    # (authoritative, refed at airports like SPLP whose OSM is
    # unrefed) over OSM ways (used at airports without apt.dat
    # coverage).  Per user 2026-05-14: at SPLP the OSM taxiways
    # share an unrefed northeast endpoint, so OSM-based detection
    # emits an unrefed stub at the North runway end while apt.dat
    # (29 'A'-named edges) has the proper ref.  Using apt.dat
    # gives the stub its correct ref while inheriting all the
    # downstream centering / clearance logic.
    by_ref: Dict[str, List[LineString]] = {}
    if apt_centerlines:
        for ls, name in apt_centerlines:
            # Sub-refs (letter+digit) are short connector spurs —
            # never the long primary parallel we're hunting for.
            if name and any(c.isdigit() for c in name):
                continue
            by_ref.setdefault(name, []).append(ls)
    else:
        for wid, nds, tags in ways:
            if tags.get("aeroway") != "taxiway":
                continue
            ref = tags.get("ref", "")
            # Sub-refs (letter+digit, e.g. V1, L3) are short
            # connector spurs at every airport — never the long
            # parallel taxi we're hunting for here.  All other
            # refs (or no ref) are candidates; the length filter
            # at the end (UNREFED_MIN_LEN_M for unrefed; the
            # by-ref endpoint test for everything else) keeps
            # only the long ones that touch the runway.
            if ref and any(c.isdigit() for c in ref):
                continue
            pts = []
            for n in nds:
                if n in nodes:
                    lat, lon = nodes[n]
                    pts.append(to_m(lon, lat))
            if len(pts) >= 2:
                try:
                    by_ref.setdefault(ref, []).append(LineString(pts))
                except _GEOM_EXC:
                    pass

    # Pre-compute existing rect union for overlap detection
    existing_rects_union = None
    if existing_taxi_rects:
        try:
            existing_rects_union = unary_union(
                [r for r, _, _, _ in existing_taxi_rects])
        except _GEOM_EXC:
            existing_rects_union = None

    new_stubs: List[Tuple[Polygon, LineString, str, str]] = []
    rwy_boundary = runway_union.boundary
    emitted_centers: List[Tuple[float, float]] = []
    DEDUP_DIST_M = 50.0  # de-dup stub centers within 50 m
    for ref, lines in by_ref.items():
        # Process INDIVIDUAL OSM ways (not merged).  Shared
        # endpoints between ways (e.g. SPLP -696729 and -696733
        # both end at (226,1182) which is an internal junction
        # to the runway-apron ramp) would become internal
        # vertices after linemerge — and thus hidden from
        # endpoint checking.  Per-way processing exposes each
        # OSM endpoint; the ``emitted_centers`` dedup step
        # coalesces identical runway-facing endpoints.
        if ref == "":
            # Unrefed airports: each individual way must be long
            # enough to represent a primary-parallel taxi.
            # SPLP main taxi (-696729) is 2640 m and -696733 is
            # 1456 m; both exceed the 800 m threshold.
            processed_lines = [ml for ml in lines
                               if ml.length >= UNREFED_MIN_LEN_M]
        else:
            # Refed parallels: each way is already part of a
            # named taxi.  Process all of them.
            processed_lines = list(lines)

        for ml in processed_lines:
            coords = list(ml.coords)
            if len(coords) < 2:
                continue
            # Check both endpoints for runway-terminating condition
            for end_idx in (0, -1):
                ep_pt = Point(coords[end_idx])
                d_ep = ep_pt.distance(rwy_boundary)
                endpoint_inside = (
                    runway_union.contains(ep_pt)
                    or d_ep <= ENDPOINT_INSIDE_TOL_M
                )
                endpoint_outside_near = (
                    not endpoint_inside
                    and d_ep <= OUTSIDE_NEAR_RWY_M
                )
                if not (endpoint_inside or endpoint_outside_near):
                    continue

                if endpoint_inside:
                    # Walk from the endpoint toward the interior
                    # until d_rwy > STUB_EXIT_D_M.  The "exit"
                    # vertex sits just outside the runway-apron
                    # ramp.  For REFED taxis (SPJC A/F) use the
                    # vertex directly as stub center — target
                    # happens to sit at a vertex (apex of loop
                    # ramp).  For UNREFED (SPLP curving primary)
                    # INTERPOLATE back to ``STUB_INTERP_TARGET_D_UNREFED``
                    # because the target sits BETWEEN vertices
                    # in a smooth curve.
                    step = 1 if end_idx == 0 else -1
                    exit_idx = None
                    prev_i = None
                    prev_d = 0.0
                    i = end_idx if end_idx >= 0 else len(coords) - 1
                    while 0 <= i < len(coords):
                        d = Point(coords[i]).distance(rwy_boundary)
                        inside = runway_union.contains(
                            Point(coords[i]))
                        if not inside and d > STUB_EXIT_D_M:
                            exit_idx = i
                            exit_d_val = d
                            break
                        prev_i = i
                        prev_d = d
                        i += step
                    if exit_idx is None:
                        continue
                    # Compute stub center (cx, cy)
                    if (not ref and prev_i is not None
                            and exit_d_val > prev_d
                            and prev_d < STUB_INTERP_TARGET_D_UNREFED
                            < exit_d_val):
                        frac = ((STUB_INTERP_TARGET_D_UNREFED
                                 - prev_d)
                                / (exit_d_val - prev_d))
                        interp_cx = (coords[prev_i][0]
                                     + frac * (coords[exit_idx][0]
                                               - coords[prev_i][0]))
                        interp_cy = (coords[prev_i][1]
                                     + frac * (coords[exit_idx][1]
                                               - coords[prev_i][1]))
                    else:
                        interp_cx = coords[exit_idx][0]
                        interp_cy = coords[exit_idx][1]
                else:
                    # Endpoint is OUTSIDE the runway but within
                    # OUTSIDE_NEAR_RWY_M.  The taxi curves to
                    # runway at this end but doesn't enter runway
                    # pavement (SPLP NE end at 127 m, CYXY F end
                    # at 11 m).  Walk inward along the path until
                    # d_rwy exceeds STUB_EXIT_D_M, same as the
                    # inside-endpoint branch.  Without this walk
                    # the center sits AT the runway-facing
                    # endpoint and the rect ends up touching (or
                    # buried in) the runway.  The Min-clearance
                    # guardrail and gap-midpoint centering below
                    # then handle the geometry uniformly with the
                    # inside case.
                    start_i = (end_idx if end_idx >= 0
                               else len(coords) - 1)
                    step = 1 if end_idx == 0 else -1
                    exit_idx = None
                    i = start_i
                    while 0 <= i < len(coords):
                        d = Point(coords[i]).distance(rwy_boundary)
                        if d > STUB_EXIT_D_M:
                            exit_idx = i
                            break
                        i += step
                    if exit_idx is None:
                        # Path never reaches the widening
                        # threshold — fall back to the original
                        # behaviour (stub at the endpoint) so we
                        # still emit something at the ramp.
                        exit_idx = (0 if end_idx == 0
                                    else len(coords) - 1)
                    interp_cx = coords[exit_idx][0]
                    interp_cy = coords[exit_idx][1]

                # Center stub on interpolated center (interp_cx,
                # interp_cy), length STUB_LEN_M along local path
                # direction (from prev to next of exit vertex).
                prev_idx = max(0, exit_idx - 1)
                next_idx = min(len(coords) - 1, exit_idx + 1)
                dx = coords[next_idx][0] - coords[prev_idx][0]
                dy = coords[next_idx][1] - coords[prev_idx][1]
                mag = math.hypot(dx, dy)
                if mag < 1e-6:
                    continue
                ux, uy = dx / mag, dy / mag
                cx, cy = interp_cx, interp_cy
                # ── Skip perpendicular runway crossings ─────────
                # Per user 2026-05-16: this function emits the
                # "ramp" where a primary parallel MERGES into the
                # runway (i.e. approaches the runway threshold at
                # a shallow angle, parallel-ish to the runway).
                # For a PERPENDICULAR taxiway that simply crosses
                # the runway (CYXY taxiway D crosses 14R/32L at
                # node 135 and 14L/32R at node 121), there is no
                # ramp — the runway pavement takes over.  Without
                # this check the function emits a wrong-shaped
                # stub centered on the runway-adjacent path
                # vertex, planting a junction-area-sized rect
                # right where a clean junction polygon should be.
                # Test: local path direction vs nearest runway
                # bearing.  If perp_diff < 60° (i.e. taxi is
                # within 60° of perpendicular to runway, not
                # within 30° of parallel), skip.
                if rwy_centerlines:
                    local_b = math.degrees(
                        math.atan2(ux, uy)) % 180.0
                    try:
                        nearest_r = min(
                            rwy_centerlines,
                            key=lambda r: Point(cx, cy).distance(r))
                        _rc = list(nearest_r.coords)
                        _rx = _rc[-1][0] - _rc[0][0]
                        _ry = _rc[-1][1] - _rc[0][1]
                        _rmag = math.hypot(_rx, _ry)
                    except _GEOM_EXC:
                        _rmag = 0.0
                    if _rmag > 1e-6:
                        rwy_b = math.degrees(
                            math.atan2(_rx, _ry)) % 180.0
                        db = abs(local_b - rwy_b)
                        db = min(db, 180.0 - db)
                        # db close to 0° = parallel (a ramp)
                        # db close to 90° = perpendicular (a
                        # crossing — skip).  30° is the boundary.
                        if db > 20.0:
                            continue
                # First-pass axis at default length, used only to
                # probe the local pavement width.
                ax_start = (cx - ux * STUB_LEN_M / 2,
                            cy - uy * STUB_LEN_M / 2)
                ax_end = (cx + ux * STUB_LEN_M / 2,
                          cy + uy * STUB_LEN_M / 2)
                try:
                    probe_axis = LineString([ax_start, ax_end])
                except _GEOM_EXC:
                    continue
                _nat, _p90, narrow = _natural_half_width(
                    probe_axis, pav_union)
                if narrow < 3.5 or narrow > 50.0:
                    continue
                width = 2.0 * narrow
                # Width-based stub length (user 2026-04-27 spec):
                # the stub should be roughly square so its long
                # edges sit on the apron-narrowing pavement
                # boundary, corners snap there, and surrounding
                # junctions only connect at the short edges instead
                # of wrapping around the long edges.  Cap range
                # 50..80 m so the stub never collapses to a sliver
                # nor extends beyond the natural runway-apron
                # transition length.  (Diagonal stubs would use a
                # tighter cap; A/F/L at SPJC are perpendicular by
                # construction so the same formula applies.)
                target_len = max(50.0, min(STUB_LEN_M, width + 5.0))
                # ----- Stub area centering -----
                # Distinguish L-style "primary parallel curving
                # into a stub" (narrow connector pavement) from
                # A/F-style "loop ramp" (wide rwy-end pavement).
                # The signal: pavement WIDTH at exit_idx.
                # • Wide  (≥ NARROW_PAV_M): loop ramp — keep stub
                #   centred AT exit_idx (the apex of the loop is
                #   exactly where the user wants the rect).
                # • Narrow (< NARROW_PAV_M): the "stub area" is the
                #   connector pavement between the runway and the
                #   apron-widening point (exit_idx).  Per user
                #   2026-05-14 invariant: the rect sits CENTERED
                #   in the stub area with a junction on either
                #   side (one between rect and runway, one between
                #   rect and apron).  Center the rect at the
                #   midpoint of the gap_curve so both junctions
                #   get equal room.  The MIN_RUNWAY_CLEARANCE_M
                #   guardrail below then trims the axis if either
                #   end is still too close to the runway after
                #   centering.
                NARROW_PAV_M = 60.0
                PULL_BACK_FRAC = 0.5
                # The pull-back applies for narrow connector
                # pavement on BOTH endpoint branches: when the
                # polyline endpoint sits inside the runway
                # (endpoint_inside; e.g. SPJC L/F) the runway
                # anchor is the polyline's runway-boundary
                # crossing; when the endpoint sits outside but
                # within OUTSIDE_NEAR_RWY_M (endpoint_outside_near;
                # e.g. CYXY F at d_rwy ≈ 11 m) the anchor is the
                # endpoint itself.  Wide pavement (loop ramps,
                # SPJC F-style) skips the pull-back to keep the
                # stub centred at the loop apex (exit_idx).
                apply_pullback = (
                    ref and width < NARROW_PAV_M
                    and (endpoint_inside or endpoint_outside_near))
                if apply_pullback:
                    pull_path: List[Tuple[float, float]] = []
                    s = 1 if end_idx == 0 else -1
                    if endpoint_inside:
                        # Find runway-boundary crossing (last
                        # in-rwy vertex → first out-of-rwy
                        # vertex; intersect the connecting
                        # segment with the runway boundary).
                        k = (end_idx if end_idx >= 0
                             else len(coords) - 1)
                        last_in = None
                        while 0 <= k < len(coords):
                            ptk = Point(coords[k])
                            dk = ptk.distance(rwy_boundary)
                            if (runway_union.contains(ptk)
                                    or dk <= ENDPOINT_INSIDE_TOL_M):
                                last_in = k
                                k += s
                            else:
                                break
                        if (last_in is not None
                                and 0 <= k < len(coords)):
                            cross_pt = coords[k]
                            try:
                                seg = LineString(
                                    [coords[last_in], coords[k]])
                                cd = seg.difference(runway_union)
                                if (not cd.is_empty
                                        and cd.geom_type
                                        == "LineString"):
                                    cc = list(cd.coords)
                                    d0 = math.hypot(
                                        cc[0][0] - coords[last_in][0],
                                        cc[0][1] - coords[last_in][1])
                                    d1 = math.hypot(
                                        cc[-1][0] - coords[last_in][0],
                                        cc[-1][1] - coords[last_in][1])
                                    cross_pt = (cc[0] if d0 < d1
                                                else cc[-1])
                            except _GEOM_EXC:
                                pass
                            pull_path.append(
                                (cross_pt[0], cross_pt[1]))
                            m = k
                        else:
                            m = None
                    else:
                        # endpoint_outside_near: anchor at the
                        # OSM endpoint (the closest path point to
                        # the runway) and walk inward to exit_idx.
                        start_i = (end_idx if end_idx >= 0
                                   else len(coords) - 1)
                        pull_path.append(
                            (coords[start_i][0],
                             coords[start_i][1]))
                        m = start_i
                    if m is not None:
                        while True:
                            pull_path.append(
                                (coords[m][0], coords[m][1]))
                            if m == exit_idx:
                                break
                            m += s
                            if not (0 <= m < len(coords)):
                                break
                    if len(pull_path) >= 2:
                        try:
                            gap_curve = LineString(pull_path)
                            gap = gap_curve.length
                        except _GEOM_EXC:
                            gap = 0.0
                        if gap > 30.0:
                            # New centre at PULL_BACK_FRAC along
                            # the gap path measured from the
                            # runway-boundary crossing.  At
                            # PULL_BACK_FRAC = 0.5 this is the
                            # midpoint of the stub area — runway
                            # on one side, apron-widening on the
                            # other, with equal room for a
                            # junction on each side.
                            new_along = PULL_BACK_FRAC * gap
                            cpt = gap_curve.interpolate(new_along)
                            cx, cy = cpt.x, cpt.y
                            # Local tangent at the new centre.
                            eps_ = max(1.0, gap * 0.02)
                            ta = max(0.0, new_along - eps_)
                            tb = min(gap, new_along + eps_)
                            pa = gap_curve.interpolate(ta)
                            pb = gap_curve.interpolate(tb)
                            tdx = pb.x - pa.x
                            tdy = pb.y - pa.y
                            tmag = math.hypot(tdx, tdy)
                            if tmag > 1e-6:
                                ux, uy = tdx / tmag, tdy / tmag
                ax_start = (cx - ux * target_len / 2,
                            cy - uy * target_len / 2)
                ax_end = (cx + ux * target_len / 2,
                          cy + uy * target_len / 2)
                # Per user 2026-05-14: a stub must never touch the
                # runway — the rect sits inside the stub area with
                # a junction between the rect's runway-facing
                # short edge and the runway boundary.  If either
                # axis endpoint is within MIN_RUNWAY_CLEARANCE_M
                # of the runway boundary (after which corner snap
                # would glue the short-edge corners onto the
                # runway boundary), shrink the axis SYMMETRICALLY
                # from both ends to keep the rect centred on
                # (cx, cy) while pulling the runway-facing end
                # back to MIN_RUNWAY_CLEARANCE_M.  If the resulting
                # axis is shorter than MIN_STUB_LEN_M, drop the
                # stub entirely — the connector is too narrow to
                # accommodate both junctions and a meaningful rect.
                # Clearance must exceed RUNWAY_ADJACENCY_TOL_M
                # (config.py = 20 m).  Otherwise the rect's runway-
                # facing corners sit inside the runway-adjacency
                # window and ``_enforce_runway_1to1_sharing`` later
                # snaps them to runway corners — stretching the
                # adjacent junction polygon into the runway-side
                # gap and producing junction/junction overlap.
                # SPJC target B/G runway-end stubs sit at d_rwy
                # ≈ 20-26 m, validating ~22 m as the natural
                # value.
                MIN_RUNWAY_CLEARANCE_M = 22.0
                MIN_STUB_LEN_M = 25.0
                # Use runway_union (a solid) so an axis endpoint
                # inside the runway reads d=0 — boundary distance
                # alone wraps around and reports a small positive
                # value for points just past the boundary into the
                # runway, which would let the shortfall test pass
                # while the rect's runway-facing corners still glue
                # onto the runway boundary in the snap step.
                try:
                    d_start = Point(*ax_start).distance(runway_union)
                    d_end = Point(*ax_end).distance(runway_union)
                except _GEOM_EXC:
                    d_start = d_end = float("inf")
                shortfall = max(
                    MIN_RUNWAY_CLEARANCE_M - d_start,
                    MIN_RUNWAY_CLEARANCE_M - d_end,
                    0.0,
                )
                if shortfall > 0.0:
                    new_len = target_len - 2.0 * shortfall
                    if new_len < MIN_STUB_LEN_M:
                        continue
                    target_len = new_len
                    ax_start = (cx - ux * target_len / 2,
                                cy - uy * target_len / 2)
                    ax_end = (cx + ux * target_len / 2,
                              cy + uy * target_len / 2)
                try:
                    stub_axis = LineString([ax_start, ax_end])
                except _GEOM_EXC:
                    continue
                # Snap against a pav with the runway-clearance
                # buffer carved out so corner snapping and
                # perpendicular extension can't pull the rect's
                # short-edge corners onto the runway boundary.
                # Without this, the runway boundary is the
                # closest pav.boundary point for an axis endpoint
                # placed at d_rwy ≈ MIN_RUNWAY_CLEARANCE_M and
                # the snap glues the corner to the runway.  The
                # carved boundary introduces a new edge at the
                # clearance offset; corners snap there instead,
                # preserving the runway gap that the surrounding
                # junction needs.
                try:
                    pav_for_snap = pav_union.difference(
                        runway_union.buffer(
                            MIN_RUNWAY_CLEARANCE_M))
                    if (pav_for_snap.is_empty
                            or pav_for_snap.geom_type
                            not in ("Polygon", "MultiPolygon")):
                        pav_for_snap = pav_union
                except _GEOM_EXC:
                    pav_for_snap = pav_union
                rect = _rect_from_axis_extended(
                    stub_axis, width, pav_for_snap,
                    apt_vertices=apt_vertices)
                if rect is None or rect.is_empty:
                    continue
                if not rect.is_valid:
                    try:
                        rect = rect.buffer(0)
                    except _GEOM_EXC:
                        continue
                    if (rect.is_empty
                            or rect.geom_type != "Polygon"):
                        continue
                # Per user 2026-04-27: F-style runway-end stubs sit
                # in the wide runway-end ramp; the half-width probe
                # (capped at RAY_CAP_M = 40 m) under-sizes the rect
                # so the ramp extends past the rect's long edges and
                # the surrounding junction wraps around them.  After
                # the rect is built and snapped, ray-cast each corner
                # outward perpendicular to the axis until it hits
                # the apt.dat pavement boundary — turning the rect
                # into a trapezoid that covers the FULL ramp width
                # at each end independently.  Uses the same
                # runway-buffered pav so the extension can't
                # cross into the clearance gap.
                rect = _extend_rect_corners_perpendicular(
                    rect, stub_axis, pav_for_snap)
                if rect is None or rect.is_empty:
                    continue
                # Skip if the stub would overlap an existing rect
                # or a same-ref rect's 30 m buffer (duplicates the
                # L SE stub that the main pipeline already emits).
                skip = False
                if existing_rects_union is not None:
                    try:
                        overlap = rect.intersection(
                            existing_rects_union).area
                        if overlap > rect.area * 0.2:
                            skip = True
                    except _GEOM_EXC:
                        pass
                if not skip and ref:
                    # Same-ref near-duplicate guard: applies only
                    # when ref is non-empty (avoid dropping all
                    # unrefed-airport stubs since every existing
                    # rect has ref="" too).
                    for er, _, _, eref in existing_taxi_rects:
                        if eref != ref:
                            continue
                        try:
                            if er.buffer(30.0).intersects(rect):
                                skip = True
                                break
                        except _GEOM_EXC:
                            pass
                if skip:
                    continue
                # Dedup by stub CENTER position — for unrefed
                # airports, two ways can share the same runway-
                # facing endpoint (e.g. SPLP -696729 and -696733
                # both end at (226,1182)) and produce near-
                # identical stubs.
                c = rect.centroid
                dup = False
                for (ex, ey) in emitted_centers:
                    if math.hypot(c.x - ex, c.y - ey) < DEDUP_DIST_M:
                        dup = True
                        break
                if dup:
                    continue
                emitted_centers.append((c.x, c.y))
                new_stubs.append(
                    (rect, stub_axis, ROLE_STUB, ref))
    return new_stubs


def _clip_residue_at_stub_sloping_edges(
        residue: "Polygon",
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        outer_buffer_m: float = 30.0,
        ) -> "Polygon":
    """Subtract a thin strip just OUTSIDE each STUB rect's sloping
    edges from the residue.  Per user 2026-04-27 invariant: a
    junction polygon must never run along a sloping rect's
    sloping edge ("sloping" = parallel to source_axis = the slope
    direction, regardless of geometric length) — the rect's
    sloping edges are TAXI BOUNDARIES, not junction boundaries.
    When apt.dat pavement bulges past a stub's sloping edge
    between its two cross-edge corners, the bulge becomes a
    "wrap" on the adjacent junction.  Removing the bulge from the
    residue forces the junction to stop at the stub's cross-edge
    corners.

    The strip extends ``outer_buffer_m`` past each sloping edge —
    enough to swallow typical apt.dat curvature noise (the SPJC F
    sloping-edge bulge is 21 m).  Stubs follow
    ``_rect_from_axis_extended``'s corner convention: corners
    [0,1] form one sloping edge "side1", [2,3] form the other.
    """
    if residue is None or residue.is_empty:
        return residue
    for rect, axis, role, ref in taxi_rects:
        if role != ROLE_STUB:
            continue
        try:
            rc = list(rect.exterior.coords)
        except _GEOM_EXC:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Sloping edges per ``_rect_from_axis_extended`` convention.
        sloping_edges = [(rc[0], rc[1]), (rc[2], rc[3])]
        for (e0, e1) in sloping_edges:
            ex = e1[0] - e0[0]
            ey = e1[1] - e0[1]
            mag = math.hypot(ex, ey)
            if mag < 0.5:
                continue
            ux, uy = ex / mag, ey / mag
            # Outward normal (away from the rect centroid).
            nx, ny = -uy, ux
            cx_r, cy_r = rect.centroid.x, rect.centroid.y
            mid_x = 0.5 * (e0[0] + e1[0])
            mid_y = 0.5 * (e0[1] + e1[1])
            if (cx_r - mid_x) * nx + (cy_r - mid_y) * ny > 0:
                nx, ny = -nx, -ny
            # Build a thin rectangle outside the sloping edge.
            o0 = (e0[0] + nx * outer_buffer_m,
                  e0[1] + ny * outer_buffer_m)
            o1 = (e1[0] + nx * outer_buffer_m,
                  e1[1] + ny * outer_buffer_m)
            try:
                strip = Polygon([e0, e1, o1, o0])
                if strip.is_valid and not strip.is_empty:
                    residue = residue.difference(strip)
            except _GEOM_EXC:
                continue
    return residue


def _add_stub_to_runway_bridges(
        residue: "Polygon",
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        runway_union: Optional["Polygon"],
        max_bridge_m: float = 100.0,
        ) -> "Polygon":
    """For each STUB rect whose runway-facing short edge has a
    pavement GAP to the runway, ADD a synthetic quadrilateral from
    the stub's runway-facing edge straight to the runway boundary.

    Per user 2026-04-27: airport apt.dat data is sometimes
    inconsistent with OSM (construction-era discrepancies); when
    we see a stub that should connect to a runway but the apt.dat
    pavement doesn't bridge the gap, project a straight line over
    so the connecting junction has continuous coverage instead of
    a "weird shape with a big gap".

    The runway-facing short edge is identified as whichever of the
    stub's two short edges (corners [0,3] or [1,2]) has its
    midpoint nearest the runway boundary.  If that midpoint sits
    > ``max_bridge_m`` from the runway, we don't bridge (probably
    not a runway-side stub).
    """
    if residue is None or residue.is_empty:
        return residue
    if runway_union is None or runway_union.is_empty:
        return residue
    rwy_boundary = runway_union.boundary
    additions: List[Polygon] = []
    for rect, axis, role, ref in taxi_rects:
        if role != ROLE_STUB:
            continue
        try:
            rc = list(rect.exterior.coords)
        except _GEOM_EXC:
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            continue
        # Two short edges per ``_rect_from_axis_extended``
        # convention: [3,0] and [1,2].
        short_edges = [(rc[3], rc[0]), (rc[1], rc[2])]
        # Pick the one closer to the runway.
        best_edge = None
        best_d = max_bridge_m
        for (e0, e1) in short_edges:
            mid = Point(0.5 * (e0[0] + e1[0]),
                        0.5 * (e0[1] + e1[1]))
            try:
                d = mid.distance(rwy_boundary)
            except _GEOM_EXC:
                continue
            if d < best_d:
                best_d = d
                best_edge = (e0, e1)
        if best_edge is None or best_d <= 1.0:
            # No runway-side short edge in range, OR the stub
            # already touches the runway — nothing to bridge.
            continue
        # Project each short-edge endpoint to the nearest runway-
        # boundary point.
        from shapely.ops import nearest_points
        e0, e1 = best_edge
        try:
            n0, _ = nearest_points(rwy_boundary, Point(e0))
            n1, _ = nearest_points(rwy_boundary, Point(e1))
        except _GEOM_EXC:
            continue
        # Build the bridge quadrilateral: stub edge → runway edge.
        # Order: e0, e1, n1, n0 so the bridge closes properly.
        try:
            bridge = Polygon([e0, e1,
                              (n1.x, n1.y), (n0.x, n0.y)])
            if not bridge.is_valid:
                bridge = bridge.buffer(0)
            if (bridge.is_empty
                    or bridge.geom_type != "Polygon"):
                continue
            # Don't overlap with the rect itself or the runway.
            bridge = bridge.difference(rect)
            if (not bridge.is_empty
                    and bridge.geom_type == "Polygon"):
                bridge = bridge.difference(runway_union)
            if (not bridge.is_empty
                    and bridge.geom_type == "Polygon"
                    and bridge.area >= 1.0):
                additions.append(bridge)
        except _GEOM_EXC:
            continue
    if additions:
        try:
            residue = unary_union([residue] + additions)
        except _GEOM_EXC:
            pass
    return residue
