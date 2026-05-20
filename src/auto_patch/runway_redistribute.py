"""Redistribute runway-segment altitudes after seam DEM anchors enter
the profile (user 2026-05-19).

Background
----------
``runway_segments.generate_patch_osm`` emits an FAA-compliant runway
profile at the time it runs — anchored by CIFP thresholds, cross-runway
projection anchors, and centerline-crossing reconciliation, all passed
through the parabolic envelope clamp + hard-cap + rate-of-change
gates.  But it knows nothing about tile-boundary seams: the seam
pipeline runs later, identifies which runway vertices sit on integer
lat/lon lines, and pins them HARD to ``dem.alt_strict`` (the raw HGT
pixel — required for cross-tile parity / preserve_boundary).

The old ``runway_regrade`` only adjusted the two threshold corners
when seam altitudes were added.  Interior segment-boundary corners
stayed at their emit-time CIFP values, which meant the runway's
combined profile (post-regrade) was no longer FAA-compliant —
adjacent sub-rects could disagree on grade or curvature at the
shared seam, and the per-surface solver had to paper over the
result.

This module finishes the job the regrade started.  For each runway
pair, it:

  1. Pulls the emit-time profile (``layout._runway_profile_state``)
     — the sample grid + anchor flags + phys-end coordinates +
     blast pad lengths.
  2. Folds every seam vertex on that runway into the sample list as
     an additional ANCHORED sample at its DEM altitude.
  3. Re-runs the same FAA gates the emit step ran
     (``runway_segments.faa_joint_solve``) — envelope clamp + hard
     cap + rate-of-grade-change — so the non-anchored interior
     samples shift to honor the new seam anchors while staying
     FAA-compliant.
  4. Evaluates the new profile at every runway sub-rect's vertex
     position (projected onto the runway axis) and writes back per-
     vertex ``node_altitudes``.  Sub-rects that came in as 4-corner
     canonical sloped rects get converted to ``node_altitudes`` if
     any corner moved away from the canonical ``[H, L, L, H]``
     pattern; otherwise they keep their ``altitude_high``/
     ``altitude_low``.

After this pass, every runway vertex carries an FAA-compliant
altitude derived from the same sample grid + same anchor set + same
gates as the emit-time profile.  Adjacent shapes (junctions, aprons,
taxiways) anchor to the runway-emitted altitudes via shared corners
exactly as before.

Cross-tile parity: deterministic from layout geometry + the seam
DEM samples (which preserve_boundary keeps identical between
neighbouring tiles).  Both tile builds compute the same augmented
sample list and run the same iterative passes, so they converge to
the same profile values at every vertex.

Public API: ``redistribute_runway_profile``.
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, List, Tuple

from .layout import ROLE_RUNWAY, SHARED_VERTEX_TOL_M
from .pavement.runway_segments import (
    MAX_RUNWAY_GRADE, MAX_RUNWAY_GRADE_CHANGE_PER_M, faa_joint_solve,
)
from .runway_regrade import regrade_runway, DEFAULT_ARC_K_M


__all__ = ["redistribute_runway_profile"]


def _bucket_key(x: float, y: float) -> Tuple[int, int]:
    s = 1.0 / SHARED_VERTEX_TOL_M
    return (int(round(x * s)), int(round(y * s)))


def _format_ref(desig_a: str, desig_b: str) -> str:
    """Convert a CIFP designator pair to the apt.dat-style ref used as
    ``shape.ref``.  Mirrors the convention in
    ``elevation.py::_ref_from_desig_pair`` — strip the ``RW`` prefix
    so refs match the row-100 format (e.g. ``02/20``).
    """
    def _strip(d):
        if not d:
            return d
        return d[2:] if d.startswith("RW") else d
    a = _strip(desig_a)
    b = _strip(desig_b)
    if a and b:
        return f"{a}/{b}"
    return a or b or ""


def _interp_profile(fractions: List[float], elevs: List[float],
                    t: float) -> float:
    """Linear-interp ``elevs`` at ``t`` over ``fractions`` (assumed
    sorted ascending).  Clamps to endpoints outside the range.
    """
    if not fractions:
        return 0.0
    if t <= fractions[0]:
        return elevs[0]
    if t >= fractions[-1]:
        return elevs[-1]
    # Binary search would be faster but the sample lists are short
    # (≤ 30 entries on typical airports).
    for i in range(len(fractions) - 1):
        f0 = fractions[i]
        f1 = fractions[i + 1]
        if f0 <= t <= f1:
            span = f1 - f0
            if span < 1e-12:
                return elevs[i]
            u = (t - f0) / span
            return elevs[i] + u * (elevs[i + 1] - elevs[i])
    return elevs[-1]


def _shift_thresholds_for_seams(
        fractions: List[float], elevs: List[float],
        anchored: List[bool],
        seam_samples: List[Tuple[float, float]],
        phys_dist: float) -> None:
    """Shift the first/last anchored samples (the runway's two
    threshold ends) via ``regrade_runway`` so they're reachable from
    every interior anchor (existing extras + new seams) within FAA
    grade + K-factor.

    The CIFP threshold altitudes act as a soft preference — the
    optimisation clips them to the feasible band defined by the
    union of interior-anchor envelopes and minimises the shift.
    Mutates ``elevs`` in place at the two threshold indices.
    """
    if not seam_samples:
        return
    # Find first and last anchored samples — the runway's two
    # threshold ends.
    first_i = None
    last_i = None
    for i, a in enumerate(anchored):
        if a:
            if first_i is None:
                first_i = i
            last_i = i
    if first_i is None or last_i is None or first_i == last_i:
        return

    t_a = fractions[first_i]
    t_b = fractions[last_i]
    axis_len = (t_b - t_a) * phys_dist
    if axis_len < 1.0:
        return

    # Build the interior-anchor list for ``regrade_runway``: all
    # currently-anchored samples between thresholds (cross-runway
    # projections + centerline crossings) plus the new seam samples,
    # each as (dist_from_threshold_A, altitude).
    interior: List[Tuple[float, float]] = []
    for i in range(first_i + 1, last_i):
        if anchored[i]:
            d = (fractions[i] - t_a) * phys_dist
            interior.append((d, elevs[i]))
    for t, e in seam_samples:
        d = (t - t_a) * phys_dist
        if 0.5 < d < axis_len - 0.5:
            interior.append((d, e))
    if not interior:
        return

    cifp_a = elevs[first_i]
    cifp_b = elevs[last_i]
    result = regrade_runway(
        cifp_a, cifp_b, axis_len, interior,
        grade_cap=MAX_RUNWAY_GRADE,
        arc_K_m=DEFAULT_ARC_K_M)
    elevs[first_i] = result.threshold_A
    elevs[last_i] = result.threshold_B


def _insert_seam_anchors(fractions: List[float], elevs: List[float],
                          anchored: List[bool],
                          seam_samples: List[Tuple[float, float]]
                          ) -> None:
    """Merge ``seam_samples`` (list of (t, elev)) into the parallel
    arrays, preserving sort order on ``fractions``.  A seam coinciding
    (within 1e-3 in t) with an existing sample takes over that sample
    — sets anchored=True and overrides elev.
    """
    for t, e in sorted(seam_samples):
        if t < 0.0 or t > 1.0:
            continue
        # Find insertion / match position.
        matched = False
        insert_at = len(fractions)
        for i, f in enumerate(fractions):
            if abs(f - t) < 1e-3:
                anchored[i] = True
                elevs[i] = e
                matched = True
                break
            if f > t:
                insert_at = i
                break
        if matched:
            continue
        fractions.insert(insert_at, t)
        elevs.insert(insert_at, e)
        anchored.insert(insert_at, True)


def _find_centerline_boundary_crossings(
        phys_end_a_ll: Tuple[float, float],
        phys_end_b_ll: Tuple[float, float],
        dem,
        tile_lat: int,
        tile_lon: int) -> List[Tuple[float, float]]:
    """Find every integer lat/lon line crossed by the runway's
    centerline and sample DEM at each crossing point.

    Returns ``[(t, altitude), ...]`` where ``t`` is the fraction along
    the centerline from ``phys_end_a`` (matches the convention of the
    sample list stored in ``layout._runway_profile_state``) and
    ``altitude`` is the DEM altitude at the centerline-boundary
    intersection point.

    Why one sample per boundary instead of one per boundary-runway
    vertex (user 2026-05-19): when a runway crosses a tile boundary
    at a shallow angle (e.g. SPLP runway 02/20 is 18° off the
    lon=-77 boundary, slicing diagonally through 148 m of runway
    length), the boundary-runway polygon intersection yields multiple
    vertices fanning across the runway's 45 m width at different
    latitudes.  Each latitude samples a different DEM pixel; on
    rolling terrain, those samples differ by metres even though
    they're all "the same crossing."  Anchoring at the centerline
    crossing gives one FAA-feasible altitude — the value the runway
    actually has where its centerline meets the contour line in the
    real world.

    Cross-tile parity: deterministic from CIFP centerline geometry
    plus the boundary HGT pixel, both of which preserve_boundary
    keeps identical between neighbouring tile builds.
    """
    if dem is None:
        return []
    lat_a, lon_a = phys_end_a_ll
    lat_b, lon_b = phys_end_b_ll
    nodata = getattr(dem, "nodata", -32768)

    crossings: List[Tuple[float, float]] = []

    def _sample(lat_c: float, lon_c: float):
        try:
            v = float(dem.alt_strict(
                (lon_c - tile_lon, lat_c - tile_lat)))
        except _GEOM_EXC:
            return None
        if v != v or v == nodata:
            return None
        return v

    # Integer LATITUDE lines crossed by the centerline (constant-lat).
    if abs(lat_b - lat_a) > 1e-12:
        lat_lo, lat_hi = sorted([lat_a, lat_b])
        n_lo = math.ceil(lat_lo) if lat_lo != math.floor(lat_lo) else int(lat_lo) + 1
        n_hi = math.floor(lat_hi) if lat_hi != math.floor(lat_hi) else int(lat_hi) - 1
        for n in range(int(n_lo), int(n_hi) + 1):
            t = (n - lat_a) / (lat_b - lat_a)
            if not (0.001 < t < 0.999):
                continue
            lat_c = float(n)
            lon_c = lon_a + t * (lon_b - lon_a)
            v = _sample(lat_c, lon_c)
            if v is not None:
                crossings.append((t, v))

    # Integer LONGITUDE lines crossed by the centerline (constant-lon).
    if abs(lon_b - lon_a) > 1e-12:
        lon_lo, lon_hi = sorted([lon_a, lon_b])
        n_lo = math.ceil(lon_lo) if lon_lo != math.floor(lon_lo) else int(lon_lo) + 1
        n_hi = math.floor(lon_hi) if lon_hi != math.floor(lon_hi) else int(lon_hi) - 1
        for n in range(int(n_lo), int(n_hi) + 1):
            t = (n - lon_a) / (lon_b - lon_a)
            if not (0.001 < t < 0.999):
                continue
            lon_c = float(n)
            lat_c = lat_a + t * (lat_b - lat_a)
            v = _sample(lat_c, lon_c)
            if v is not None:
                crossings.append((t, v))

    crossings.sort(key=lambda c: c[0])
    return crossings


_GEOM_EXC = (ValueError, TypeError, IndexError)


def redistribute_runway_profile(
        layout,
        dem=None,
        tile_lat: int = 0,
        tile_lon: int = 0) -> int:
    """Rewrite every runway sub-rect's altitudes by re-running the
    emit-time FAA-compliant profile with tile-boundary seam DEM
    altitudes folded in as additional anchored samples.

    Seam-altitude sampling (user 2026-05-19): one DEM sample per
    centerline-boundary crossing, taken at the geometric centerline
    intersection point with the integer lat/lon line.  This gives a
    single FAA-feasible altitude for the whole boundary cut through
    the runway — required for oblique crossings where the boundary
    diagonally slices across the runway's width and would otherwise
    pick up multiple per-vertex DEM samples at different latitudes
    with inconsistent altitudes.

    Mutates the runway shapes in ``layout`` in place: converts
    4-corner sloped rects to ``node_altitudes`` when any corner
    moves, leaves them canonical otherwise.

    Returns the number of runway shapes touched.
    """
    profile_state = getattr(layout, "_runway_profile_state", None)
    if not profile_state:
        return 0

    # Group runway shapes by ref so each pair's redistribution can be
    # applied to all of its sub-rects.
    shapes_by_ref: Dict[str, list] = defaultdict(list)
    for s in layout.shapes:
        if s.role != ROLE_RUNWAY:
            continue
        if s.polygon is None or s.polygon.is_empty:
            continue
        if not s.ref:
            continue
        shapes_by_ref[s.ref].append(s)

    n_touched = 0
    for pair_key, state in profile_state.items():
        desig_a, desig_b = pair_key
        ref = _format_ref(desig_a, desig_b)
        shapes = shapes_by_ref.get(ref, [])
        if not shapes:
            continue

        # Convert phys-end lat/lon to metre frame.
        phys_a_lat, phys_a_lon = state['phys_end_a_ll']
        phys_b_lat, phys_b_lon = state['phys_end_b_ll']
        ax_a_x, ax_a_y = layout.ll_to_m(phys_a_lat, phys_a_lon)
        ax_b_x, ax_b_y = layout.ll_to_m(phys_b_lat, phys_b_lon)
        ax_dx = ax_b_x - ax_a_x
        ax_dy = ax_b_y - ax_a_y
        ax_len2 = ax_dx * ax_dx + ax_dy * ax_dy
        if ax_len2 < 1.0:
            continue

        # Build the augmented sample list.
        fractions = list(state['fractions'])
        elevs = list(state['elevs'])
        anchored = list(state['anchored'])
        phys_dist = state['phys_dist_m']

        # Find per-boundary centerline crossings.  One DEM sample per
        # boundary line crossed (not per polygon vertex) — the value
        # at the centerline-boundary intersection.
        seam_samples = _find_centerline_boundary_crossings(
            state['phys_end_a_ll'], state['phys_end_b_ll'],
            dem, tile_lat, tile_lon)

        # Step 1: if any new HARD interior anchor entered the profile
        # (centerline-boundary DEMs), the existing CIFP thresholds
        # may no longer be reachable from them within FAA grade +
        # K-factor.  Shift the thresholds via the regrade_runway
        # constrained optimisation — same algorithm the older
        # threshold-only pass used, invoked here as a global step so
        # the WHOLE runway profile can then be smoothed in step 2
        # instead of leaving sub-rect interiors out of sync with the
        # shifted ends.
        #
        # Interior anchors (cross-runway projections, centerline
        # crossings already in the emit-time anchor list) are passed
        # in as immutable constraints — ``regrade_runway`` only moves
        # the two thresholds and keeps every other anchor fixed.
        if seam_samples:
            _shift_thresholds_for_seams(
                fractions, elevs, anchored,
                seam_samples, phys_dist)
            _insert_seam_anchors(fractions, elevs, anchored,
                                  seam_samples)

        # Step 2: re-run the FAA gates on the full sample list.
        # Thresholds are still anchored (at their possibly-shifted
        # values from step 1), seams are anchored at their DEM
        # altitudes, interior samples are free.  ``faa_joint_solve``
        # smooths the interior so every adjacent-edge grade stays
        # within ``MAX_RUNWAY_GRADE`` and every adjacent-triple ΔG
        # stays within the K-factor cap.
        faa_joint_solve(
            fractions, elevs, anchored, phys_dist,
            blast_a=state['blast_a_m'],
            blast_b=state['blast_b_m'],
            grade_cap=MAX_RUNWAY_GRADE,
            max_dg_per_m=MAX_RUNWAY_GRADE_CHANGE_PER_M)

        # Evaluate the new profile at every runway sub-rect's vertex.
        for s in shapes:
            ring = list(s.polygon.exterior.coords)
            ring_closed = bool(ring) and ring[0] == ring[-1]
            ring_open = ring[:-1] if ring_closed else ring
            if len(ring_open) < 3:
                continue
            new_alts: List[float] = []
            for x, y in ring_open:
                vx = x - ax_a_x
                vy = y - ax_a_y
                t = (vx * ax_dx + vy * ax_dy) / ax_len2
                e_new = _interp_profile(fractions, elevs, t)
                new_alts.append(round(e_new, 1))

            # Detect whether the new altitudes still form a canonical
            # ``[H, L, L, H]`` 4-corner sloped rect.  If so AND the
            # shape was originally that form, preserve it (keeps the
            # ``sloping_rect_canonical_form`` invariants happy and
            # avoids unnecessary node_altitudes conversions when
            # nothing actually moved).
            keep_canonical = False
            if (len(ring_open) == 4
                    and s.altitude_high is not None
                    and s.altitude_low is not None
                    and not s.node_altitudes):
                if (abs(new_alts[0] - new_alts[3]) < 0.05
                        and abs(new_alts[1] - new_alts[2]) < 0.05):
                    new_hi = (new_alts[0] + new_alts[3]) / 2.0
                    new_lo = (new_alts[1] + new_alts[2]) / 2.0
                    # Ensure ``hi`` is actually the higher pair (preserve
                    # the canonical convention).
                    if new_hi < new_lo:
                        new_hi, new_lo = new_lo, new_hi
                    if (abs(new_hi - s.altitude_high) < 0.05
                            and abs(new_lo - s.altitude_low) < 0.05):
                        # No-op: emit-time values already match.
                        continue
                    s.altitude_high = round(new_hi, 1)
                    s.altitude_low = round(new_lo, 1)
                    keep_canonical = True
                    n_touched += 1
            if keep_canonical:
                continue

            # Non-canonical (any moved corner or already
            # node_altitudes): write per-vertex.
            closed_alts = new_alts + ([new_alts[0]]
                                       if ring_closed else [])
            s.node_altitudes = closed_alts
            s.altitude = None
            s.altitude_high = None
            s.altitude_low = None
            n_touched += 1

    return n_touched
