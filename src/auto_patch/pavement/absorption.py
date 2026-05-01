"""Long-edge absorption rule for taxi rects.

THE RECURRING-REGRESSION HOT SPOT.  Every recent session that
modified this rule introduced a regression at one of the test
airports.  The rule below is reproduced VERBATIM from
``feedback_shape_rules.md`` (memory: user 2026-04-30, AUTHORITATIVE).
Do not edit the rule without explicit user confirmation; do not
add per-ref preservation exceptions; do not switch EITHER -> BOTH;
do not bump the 10 % threshold or the 5 m probe.

----------------------------------------------------------------------
A taxi rect (primary_parallel / secondary_parallel / stub /
cross_connector) has a long edge "covered" by junction-class
pavement (apron / junction / unrolled apt.dat row-110) when that
pavement extends past the long edge for >= 10 % of the rect's axial
length.  When EITHER long edge is covered to >= 10 %, the COVERED
portion is absorbed into the surrounding apron / junction; the
non-covered portion stays as a (possibly shorter) rect.  Partial
absorption -- never absorb the WHOLE rect when only part is
covered.

Implementation:
* Probe at 5 m axial steps; each step is "adjacent" if EITHER
  outside-long-edge sample is inside ``junction_pav`` (apron
  pavement minus runway minus other rects).
* Probe distance = 5 m beyond the rect's long edge (skips sub-1 m
  polygon-imprecision slivers).
* "Adjacent run" must be >= 10 % of axial steps to absorb.
* Kept fragments < 30 m are dropped (avoid emitting micro-stubs).

DO NOT change the absorption behaviour without explicit user
confirmation.  Specifically:
* Don't switch from EITHER to BOTH -- alongside-apron rects must
  absorb.
* Don't add per-ref preservation exceptions -- they override the
  rule.
* Don't bump the 10 % threshold or 5 m probe without confirmation.

The rule is authoritative even when it produces a result that
looks "wrong" at one airport -- propose the change with the
trade-off explained and wait for the user to decide, rather than
edit and report.
----------------------------------------------------------------------

Public API:
    drop_primary_parallels_embedded_in_pavement(
        taxi_rects, apt_pav_union, runway_polys=None, ...)
    split_primary_parallels_at_pavement_boundary(
        taxi_rects, pav_union, ...)

Both retain their leading-underscore names as well for backward
compatibility with internal callers in O4_Airport_Pavement_Builder.
"""
from __future__ import annotations

import math
import sys
from typing import List, Optional, Tuple

import O4_UI_Utils as UI

from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from ..layout import (
    ROLE_CROSS_CONNECTOR,
    ROLE_PRIMARY_PARALLEL,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
)

__all__ = [
    "_drop_primary_parallels_embedded_in_pavement",
    "_split_primary_parallels_at_pavement_boundary",
    "drop_primary_parallels_embedded_in_pavement",
    "split_primary_parallels_at_pavement_boundary",
]


def _drop_primary_parallels_embedded_in_pavement(
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        apt_pav_union: "Optional[Polygon]",
        runway_polys: "Optional[List[Polygon]]" = None,
        adjacency_frac: float = 0.10,
        proximity_m: float = 1.0,
        embed_frac: float = 0.95,         # legacy param, unused
        long_edge_buffer_m: float = 5.0,  # legacy param, unused
        ) -> "List[Tuple[Polygon, LineString, str, str]]":
    """Drop ``primary_parallel`` rects whose long edges sit entirely
    inside apt.dat row-110 pavement.

    Per user 2026-04-27 invariant: a junction polygon must NEVER run
    along a sloping rect's long edge — X-Plane's elevation engine
    over-constrains the junction's spline to F's long-edge altitudes
    and produces visible elevation glitches.

    There are two ways the invariant gets violated for a primary
    parallel:

    1. apt.dat pavement bulges past the long edge between the rect's
       short corners.  The residue then runs along the long edge to
       reach the bulge.  ``_clip_residue_at_stub_long_edges`` already
       handles this for STUBs by carving the bulge away.
    2. The rect sits FULLY INSIDE a paved area (apron, ramp, big
       terminal area).  apt.dat covers the rect's footprint AND
       everything around it; the residue is the apron, and it
       inevitably traces the rect's long edges as it wraps around
       the rect-shaped hole.

    Case 2 doesn't have a "carve away the bulge" fix — there's no
    bulge, the entire surrounding pavement is legitimate apron.  The
    correct behaviour is to NOT emit the rect at all: let the apron
    junction cover the rect's footprint and slope multi-
    directionally.  The primary-parallel grade rule isn't worth
    enforcing for a taxi lane that's just one of many possible
    paths through a paved apron — the apron's own elevation
    propagation produces a perfectly serviceable surface.

    SPJC's F is the canonical example.  F is a 121×60 m primary
    parallel sitting in the middle of the SE apron; row-110 pavement
    covers F + everything around it.  Dropping F lets the SE apron
    (junction -10132) absorb F's footprint instead of wrapping
    around F's long edges.

    Implementation: a primary-parallel rect is considered embedded
    when AT LEAST ONE of its two long edges is ≥ ``embed_frac``
    inside ``apt_pav_union``, AND the strip just outside that edge
    (``long_edge_buffer_m`` wide) is also covered.  Why "at least
    one" instead of "both": a primary parallel can sit on the edge
    of an apron, with one long edge inside the apron (where the
    junction wraps it — the violation we're fixing) and the other
    long edge along the apron's outer pavement boundary (where
    there's no junction at all, so no wrap to worry about).  SPJC's
    F is exactly this: south long edge fully inside the SE apron
    (100 %), north long edge on the apron's outer boundary (27 %).
    Dropping F is correct in both cases — the apron absorbs the
    embedded side, and the boundary side becomes the apron's new
    edge.  Only when NEITHER long edge is embedded does the rect
    actually carry the pavement (a primary parallel through grass)
    — those we must keep.
    """
    if apt_pav_union is None or apt_pav_union.is_empty:
        return taxi_rects
    # Per user 2026-04-28: a sloping rect cannot have a junction or
    # apron polygon running alongside its long edge.  The slope
    # along the rect's long edge is uniform (altitude_high at one
    # short edge, altitude_low at the other); any adjacent junction-
    # class pavement would have to match that slope along the seam,
    # which produces visible elevation glitches in X-Plane when the
    # junction's natural DEM slope differs from the rect's straight-
    # line slope.  Absorb the rect into the junction so the whole
    # region becomes one polygon with per-vertex node_altitudes
    # capturing the natural slope.
    #
    # "Junction-class pavement" = pav_union − runway_polys − OTHER
    # taxi rects.  This excludes:
    #   - Runway adjacency (parallel taxis legitimately run alongside
    #     a runway; the runway has its own slope already).
    #   - Other taxi-rect adjacency (two adjacent taxi rects each
    #     own their slope along their own axes).
    junction_pav = apt_pav_union
    if runway_polys:
        try:
            for r in runway_polys:
                if r is not None and not r.is_empty:
                    junction_pav = junction_pav.difference(r)
        except Exception:
            junction_pav = apt_pav_union
    # Subtract every taxi-rect's footprint (including the rect being
    # tested — that's CORRECT, because we want to know whether the
    # JUNCTION polygon will be adjacent after the rect is emitted
    # and the residue is computed).  After subtraction, the apron
    # polygon containing the rect has a rect-shaped hole; the hole's
    # boundary IS within 1 m of the rect's long edge.
    try:
        all_taxi_polys = [
            r for (r, _ax, _role, _ref) in taxi_rects
            if r is not None and not r.is_empty]
        if all_taxi_polys:
            taxi_union = unary_union(all_taxi_polys)
            if not taxi_union.is_empty:
                junction_pav = junction_pav.difference(taxi_union)
    except Exception:
        pass
    if junction_pav.is_empty:
        return taxi_rects
    # All sloping-rect roles are subject to the absorption rule.
    sloping_rect_roles = {ROLE_PRIMARY_PARALLEL,
                          ROLE_SECONDARY_PARALLEL,
                          ROLE_STUB, ROLE_CROSS_CONNECTOR}
    # Per user 2026-04-30: absorb wherever EITHER long edge has
    # junction-class pavement running alongside it.  Sloping
    # rects cannot share a long edge with an apron/junction
    # polygon — the rect's straight-line slope along the long
    # edge has to match the junction's natural DEM slope along
    # the seam, which it generally won't.  Whether the apron
    # sits on one side or both, the violation is the same.
    #
    # A separate CORRIDOR HEURISTIC above this loop preserves
    # rects that connect to a runway at one of their short
    # edges (e.g. CYXY Taxiway F: ref-tagged corridor running
    # apron→runway), so corridors aren't swallowed even when
    # they're embedded in apron pavement.  Only ``alongside``
    # rects (E primary parallel: both short edges 160+ m from
    # any runway) fall through to absorption.
    #
    # Apply PARTIAL absorption — only the axial range where
    # adjacency holds gets absorbed; the rest of the rect
    # survives as shorter rect(s).
    #
    # Probe semantics: at every 5 m axis step, compute a point
    # ``OUTER_PROBE_M`` METRES OUTSIDE each long edge and ask
    # whether that point lies *directly inside* junction-pav (no
    # buffer).  A real apron extends many meters past the rect's
    # long edge, so a point 2 m beyond the edge will hit it.  The
    # 1 m-and-buffer formulation we replaced was sensitive to
    # sub-1 m slivers caused by apt.dat / DSF polygons that render
    # ~0.5 m wider than the OSM-tagged taxi width — those slivers
    # are polygon-imprecision noise, NOT real apron adjacency, and
    # spuriously absorbed every rect at CYXY.
    #
    # Walk each rect's axis at 5 m steps.  At each step, sample
    # the outside-left and outside-right points.  A step is
    # "adjacent" when EITHER outside point is in junction-pav
    # (rect along an apron edge, on either side).  Find
    # contiguous adjacent runs ≥ 10 % of axial length; keep
    # non-adjacent intervals (≥ 30 m fragments) as new rects.
    #
    # Use cases:
    #   - F primary parallel at CYXY: south 30 % is INSIDE the apron
    #     (apron extends past both long edges) → adjacent →
    #     absorbed.  North 70 % extends out of the apron — only
    #     polygon-imprecision slivers within 0.5 m of the long edge,
    #     no substantive apron at 2 m → not adjacent → kept,
    #     extending from apron edge to the centerline bend.
    #   - E primary parallel at CYXY: apron runs alongside ONE long
    #     side for ~65 % of length, with several meters of apron
    #     past the edge → adjacent → absorbed.  Non-adjacent
    #     fragments (if any ≥ 30 m) survive.
    #   - F embedded in SPJC apron (both sides 100 %): fully
    #     absorbed.
    SAMPLE_STEP_M = 5.0
    MIN_KEPT_M = 30.0
    # Probe distance per user 2026-04-29 (HECA R analysis): a 2 m
    # probe was too sensitive to tile-imprecision slivers up to
    # ~3 m past a rect's long edge — at HECA the DSF pavement
    # along R extends 2–3 m past the rect's natural-half-width
    # edge, the 2 m probe landed inside that sliver, and the
    # absorption fired on a perfectly normal taxiway.  5 m skips
    # those slivers while still catching genuine aprons, which
    # extend many metres past the rect (CYXY's E parallel apron
    # is ~12 m past, SPJC's F apron is ~30 m past — both well
    # within range of the 5 m probe).
    OUTER_PROBE_M = 5.0

    kept: List[Tuple[Polygon, LineString, str, str]] = []
    abs_refs: List[str] = []
    n_full = 0
    n_split = 0
    n_clipped = 0
    for entry in taxi_rects:
        rect, axis, role, ref = entry
        if role not in sloping_rect_roles:
            kept.append(entry)
            continue
        try:
            rc = list(rect.exterior.coords)
        except Exception:
            kept.append(entry)
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            kept.append(entry)
            continue
        # Axis from short-edge midpoints (corners 0,3 → one short
        # edge; 1,2 → other) per ``_rect_from_axis_extended``
        # convention.
        a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                  0.5 * (rc[0][1] + rc[3][1]))
        b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                  0.5 * (rc[1][1] + rc[2][1]))
        ax = b_mid[0] - a_mid[0]
        ay = b_mid[1] - a_mid[1]
        L = math.hypot(ax, ay)
        if L < 30.0:
            kept.append(entry)
            continue
        ux = ax / L
        uy = ay / L
        nx, ny = -uy, ux
        half_w = math.hypot(rc[0][0] - a_mid[0], rc[0][1] - a_mid[1])
        if half_w < 1.0:
            kept.append(entry)
            continue

        n_steps = max(2, int(L / SAMPLE_STEP_M) + 1)
        either_adj = [False] * n_steps
        outer = half_w + OUTER_PROBE_M
        for i in range(n_steps):
            u = min(L, i * SAMPLE_STEP_M)
            cx = a_mid[0] + u * ux
            cy = a_mid[1] + u * uy
            try:
                left_pt = Point(cx + nx * outer,
                                 cy + ny * outer)
                right_pt = Point(cx - nx * outer,
                                  cy - ny * outer)
                either_adj[i] = (
                    bool(junction_pav.contains(left_pt))
                    or bool(junction_pav.contains(right_pt)))
            except Exception:
                continue

        # Find contiguous "either-side adjacent" runs ≥ 10 % of axis.
        min_run_steps = max(1, int(adjacency_frac * n_steps))
        absorbed_intervals: List[Tuple[float, float]] = []
        i = 0
        while i < n_steps:
            if not either_adj[i]:
                i += 1
                continue
            j = i
            while j < n_steps and either_adj[j]:
                j += 1
            if (j - i) >= min_run_steps:
                u_start = i * SAMPLE_STEP_M
                u_end = min(L, j * SAMPLE_STEP_M)
                absorbed_intervals.append((u_start, u_end))
            i = j

        if not absorbed_intervals:
            kept.append(entry)
            continue

        # Compute kept intervals = [0, L] minus absorbed.
        kept_intervals: List[Tuple[float, float]] = []
        prev_end = 0.0
        for s, e in absorbed_intervals:
            if s > prev_end:
                kept_intervals.append((prev_end, s))
            prev_end = e
        if L > prev_end:
            kept_intervals.append((prev_end, L))
        kept_intervals = [(s, e) for s, e in kept_intervals
                          if (e - s) >= MIN_KEPT_M]

        if not kept_intervals:
            n_full += 1
            abs_refs.append(ref or "?")
            continue

        # No-op: kept covers (almost) the full rect.
        if (len(kept_intervals) == 1
                and kept_intervals[0][0] <= SAMPLE_STEP_M
                and kept_intervals[0][1] >= L - SAMPLE_STEP_M):
            kept.append(entry)
            continue

        # Build new rects from kept intervals.  Drop any kept
        # fragment that's apron-interior — i.e. ≥ 2 of its 4 corners
        # are off the pavement boundary.  These are tiny rects
        # floating inside an apron polygon, with no real corridor
        # geometry around them; emitting them just creates a rect-
        # shaped hole in the apron's residue (and a thin diagonal
        # junction polygon connecting the hole to the apron's
        # outer boundary — see CYXY -10005/-10129 regression).
        # Mirrors ``_build_taxi_rects``'s apron-interior check.
        APRON_INTERIOR_TOL_M = 2.0
        try:
            pav_boundary = apt_pav_union.boundary
        except Exception:
            pav_boundary = None
        n_dropped_interior = 0
        new_rects: List[Tuple[Polygon, LineString, str, str]] = []
        for u_lo, u_hi in kept_intervals:
            new_a_mid = (a_mid[0] + u_lo * ux,
                         a_mid[1] + u_lo * uy)
            new_b_mid = (a_mid[0] + u_hi * ux,
                         a_mid[1] + u_hi * uy)
            new_corners = [
                (new_a_mid[0] + nx * half_w,
                 new_a_mid[1] + ny * half_w),
                (new_b_mid[0] + nx * half_w,
                 new_b_mid[1] + ny * half_w),
                (new_b_mid[0] - nx * half_w,
                 new_b_mid[1] - ny * half_w),
                (new_a_mid[0] - nx * half_w,
                 new_a_mid[1] - ny * half_w),
            ]
            try:
                new_rect = Polygon(new_corners)
                if not new_rect.is_valid:
                    new_rect = new_rect.buffer(0)
                if (new_rect.is_empty
                        or new_rect.geom_type != "Polygon"):
                    continue
                # Apron-interior check on the kept fragment.  A
                # kept fragment naturally has 2 corners at the
                # absorbed/kept split (interior to the original
                # rect's footprint, slightly off-boundary by ~1-3 m
                # due to apt.dat boundary imprecision at corridor
                # narrowings).  The discriminator is "all 4 corners
                # off-boundary" — captures truly apron-floating
                # fragments (CYXY -10005: all 4 corners 2.5-23 m
                # off) without dropping normal split-end fragments
                # (CYXY -10004: 3 corners 0-2.9 m off but at most
                # one corner > 3 m off).
                if pav_boundary is not None:
                    n_off = sum(
                        1 for (cx, cy) in new_corners
                        if Point(cx, cy).distance(pav_boundary)
                        > APRON_INTERIOR_TOL_M)
                    max_off = max(
                        (Point(cx, cy).distance(pav_boundary)
                         for (cx, cy) in new_corners),
                        default=0.0)
                    # Drop if all 4 corners off boundary AND at
                    # least one is > 5 m off (indicating real
                    # apron-floating, not boundary imprecision).
                    if n_off == 4 and max_off > 5.0:
                        n_dropped_interior += 1
                        continue
                new_axis = LineString([new_a_mid, new_b_mid])
                new_rects.append((new_rect, new_axis, role, ref))
            except Exception:
                continue
        kept.extend(new_rects)
        if not new_rects:
            n_full += 1
        elif len(new_rects) >= 2:
            n_split += 1
        else:
            n_clipped += 1
        abs_refs.append(
            f"{ref or '?'}"
            f"{'/int=' + str(n_dropped_interior) if n_dropped_interior else ''}")

    if abs_refs:
        try:
            UI.vprint(1,
                f"  [pav-builder] long-edge-adjacent absorption: "
                f"{n_full} dropped, {n_split} split, "
                f"{n_clipped} clipped (refs: "
                f"{', '.join(abs_refs)}).")
        except Exception:
            pass
    return kept


def _split_primary_parallels_at_pavement_boundary(
        taxi_rects: "List[Tuple[Polygon, LineString, str, str]]",
        pav_union: "Optional[Polygon]",
        step_m: float = 10.0,
        min_dropped_m: float = 100.0,
        min_kept_m: float = 50.0,
        ) -> "List[Tuple[Polygon, LineString, str, str]]":
    """Clip primary_parallel rects whose long edge has a contiguous
    embedded sub-range at one end of the rect's axis.

    Per user 2026-04-27: when a taxiway is bounded by apron pavement
    along most of one side, that bounded section should be ABSORBED
    INTO THE APRON (covered by it, removed from the taxi rect).
    The rule fires for FULL embedding via
    ``_drop_primary_parallels_embedded_in_pavement``; this helper
    handles PARTIAL embedding where only a sub-range of the rect's
    axis is bounded.  The unbounded portion stays as a (shorter)
    sloping rect; the absorbed portion's footprint reappears in
    the residue → junction polygon → apron extends to cover it.

    Algorithm:
        1. Walk the rect's axis at ``step_m``-metre steps.  At each
           step ``u``, sample the left and right long-edge midpoints
           at axial offset ``u`` and check whether either falls
           inside ``pav_union``.
        2. Mark each step as "embedded" or not.
        3. Find embedded contiguous prefixes / suffixes ≥
           ``min_dropped_m`` long.
        4. If clipping the embedded prefix and/or suffix leaves a
           kept range ≥ ``min_kept_m``, replace the rect with the
           clipped version.

    The clipped rect uses the same axis direction; its corners are
    at the new axial bounds projected to the original perpendicular
    half-width.

    At CYXY, taxiway E (469 m, NW-SE oriented) has its NW half
    bounded by SW apron pavement on its west side.  This helper
    clips E's NW portion (~290 m) and keeps the SE portion
    (~180 m), leaving the apron polygon to cover E's NW footprint.
    """
    if pav_union is None or pav_union.is_empty:
        return taxi_rects
    out: List[Tuple[Polygon, LineString, str, str]] = []
    n_clipped = 0
    for entry in taxi_rects:
        rect, axis, role, ref = entry
        if role != ROLE_PRIMARY_PARALLEL:
            out.append(entry)
            continue
        try:
            rc = list(rect.exterior.coords)
        except Exception:
            out.append(entry)
            continue
        if rc and rc[0] == rc[-1]:
            rc = rc[:-1]
        if len(rc) != 4:
            out.append(entry)
            continue
        # Axis from short-edge-A midpoint (corners 0,3) to short-
        # edge-B midpoint (corners 1,2).  Per the convention used by
        # ``_rect_from_axis_extended`` / runway-segment emit, corners
        # [0, 1] form one long edge and [2, 3] form the other.
        a_mid = (0.5 * (rc[0][0] + rc[3][0]),
                  0.5 * (rc[0][1] + rc[3][1]))
        b_mid = (0.5 * (rc[1][0] + rc[2][0]),
                  0.5 * (rc[1][1] + rc[2][1]))
        ax = b_mid[0] - a_mid[0]
        ay = b_mid[1] - a_mid[1]
        L = math.hypot(ax, ay)
        if L < 1.0:
            out.append(entry)
            continue
        ux, uy = ax / L, ay / L
        nx, ny = -uy, ux
        half_w = math.hypot(rc[0][0] - a_mid[0], rc[0][1] - a_mid[1])
        if half_w < 1.0:
            out.append(entry)
            continue
        n_steps = max(2, int(math.floor(L / step_m)) + 1)
        embedded: List[bool] = []
        try:
            for i in range(n_steps):
                u = min(L, i * step_m)
                cx = a_mid[0] + u * ux
                cy = a_mid[1] + u * uy
                left = Point(cx + nx * half_w, cy + ny * half_w)
                right = Point(cx - nx * half_w, cy - ny * half_w)
                emb = (pav_union.contains(left)
                       or pav_union.contains(right))
                embedded.append(emb)
        except Exception:
            out.append(entry)
            continue
        # Find embedded prefix length (in steps).
        pfx = 0
        while pfx < n_steps and embedded[pfx]:
            pfx += 1
        # Find embedded suffix length (in steps).
        sfx = 0
        while sfx < n_steps - pfx and embedded[n_steps - 1 - sfx]:
            sfx += 1
        # Convert to metres.  step_m is the axis step; the embedded
        # prefix ends at the LAST embedded step's u-coordinate, so
        # the prefix length in metres is ``pfx_steps × step_m``.
        pfx_m = pfx * step_m
        sfx_m = sfx * step_m
        # Allowable clip: prefix or suffix must be ≥ min_dropped_m,
        # and the kept range must be ≥ min_kept_m.
        clip_pfx = pfx_m >= min_dropped_m
        clip_sfx = sfx_m >= min_dropped_m
        if not (clip_pfx or clip_sfx):
            out.append(entry)
            continue
        u_lo = pfx_m if clip_pfx else 0.0
        u_hi = (L - sfx_m) if clip_sfx else L
        if u_hi - u_lo < min_kept_m:
            out.append(entry)
            continue
        # Build clipped rect: 4 corners at u_lo and u_hi axial
        # offsets, ±half_w perpendicular.
        new_a_mid = (a_mid[0] + u_lo * ux, a_mid[1] + u_lo * uy)
        new_b_mid = (a_mid[0] + u_hi * ux, a_mid[1] + u_hi * uy)
        new_corners = [
            (new_a_mid[0] + nx * half_w, new_a_mid[1] + ny * half_w),
            (new_b_mid[0] + nx * half_w, new_b_mid[1] + ny * half_w),
            (new_b_mid[0] - nx * half_w, new_b_mid[1] - ny * half_w),
            (new_a_mid[0] - nx * half_w, new_a_mid[1] - ny * half_w),
        ]
        try:
            new_rect = Polygon(new_corners)
            if not new_rect.is_valid:
                new_rect = new_rect.buffer(0)
            if (new_rect.is_empty
                    or new_rect.geom_type != "Polygon"):
                out.append(entry)
                continue
            new_axis = LineString([new_a_mid, new_b_mid])
            out.append((new_rect, new_axis, role, ref))
            n_clipped += 1
            try:
                UI.vprint(1,
                    f"  [pav-builder] clipped primary_parallel "
                    f"{ref!r}: {L:.0f}m → {(u_hi - u_lo):.0f}m "
                    f"(dropped "
                    + ("prefix " if clip_pfx else "")
                    + ("suffix " if clip_sfx else "")
                    + "embedded in pavement).")
            except Exception:
                pass
        except Exception:
            out.append(entry)
            continue
    return out




# Public-name aliases (no leading underscore).
drop_primary_parallels_embedded_in_pavement = (
    _drop_primary_parallels_embedded_in_pavement)
split_primary_parallels_at_pavement_boundary = (
    _split_primary_parallels_at_pavement_boundary)
