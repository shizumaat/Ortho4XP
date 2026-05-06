# Auto-Patch Status — GEOMETRY BASELINE 2026-05-05

## TL;DR

This session reduced the pavement-builder to a clean **geometry-only**
pipeline with a single source of truth for the pavement boundary, and
established the user-approved geometry as a tagged baseline for future
work.  Phase-2 elevations / widening / runway segmentation are
**disabled** on this baseline so the geometry can be reviewed in
isolation; a follow-up agent will re-enable them.

**Branch / tag:** ``dev`` branch, tagged ``geometry-baseline``.

**Build artifact:** [/tmp/SPJC_now.osm](file:///tmp/SPJC_now.osm) —
SPJC built with ``compute_elevations=False``.  All 188 rect corners
land on ``pav_union.boundary`` (worst distance 0.0 m).

## Final SPJC numbers

```
junctions emitted    : 38   (5 sliver triangles dropped at to_osm)
rects (taxi)         : 46   (4 F primary_parallel + 1 F stub absorbed/dropped)
terminals            : 2
runway segments      : 2    (un-segmented; phase-2 disabled)
total OSM ways       : 89
rect corners on pav  : 188 / 188 (worst d = 0.000 m)
```

## Architectural changes this session

### 1. Single ``pav_union`` for everything
Removed the ``pav_union_for_rects`` split.  Both rect building and
junction emit now operate on the same pav_union (pav − effective_runway
− groundside).  Earlier two-variant scheme caused rect corners to
land on a boundary that didn't match the boundary residue subtraction
used at junction emit, producing thin sliver tabs.

* [pipeline.py:712-732](src/auto_patch/pipeline.py:712) — pav_union
  built, simplified, mutated.  ``layout._pav_union_for_rects`` is no
  longer consulted.
* [pavement/rects.py:_build_taxi_rects](src/auto_patch/pavement/rects.py)
  — removed internal ``pav_union − rwy_union`` step; uses pav_union
  directly.

### 2. ``pav_union`` simplified at construction
Apt.dat row-110 polygons routinely contain over-resolved curves
(sub-meter steps) and 0.2 m doubled-vertex needles.  Without
simplification, OSM-emit bucket-collision dedup inflates polygons by
~196K m² and creates sliver corners.

* [pavement/union_helpers.py:_simplify_pavement_polygon](src/auto_patch/pavement/union_helpers.py)
  — DP simplify (1 m tolerance, preserve_topology=True) followed by
  ``_drop_sliver_corners`` to snip needle tips.
* [pipeline.py:725](src/auto_patch/pipeline.py:725) — applied at
  construction so all downstream consumers see the clean polygon.
* Terminals also simplified before subtraction in
  [junction_emit.py](src/auto_patch/junction_emit.py).

### 3. Stub classifier revised
Per user 2026-05-05: a stub can be at a wide range of angles BUT
cannot be directly parallel to the runway.

* [pavement/rects.py:_classify_role](src/auto_patch/pavement/rects.py)
  — db < 20° (parallel) → always primary_parallel or
  secondary_parallel based on distance to runway, never stub.  The
  earlier "short parallel → stub" branch is removed.  Sub-ref
  digit-suffix and diagonal-parent-ref overrides apply only when the
  rect is non-parallel.

Result at SPJC: F was 3 primary_parallel + 2 stubs → now 4 primary_parallel
+ 1 stub.  (The remaining stub was then absorbed by long-edge
adjacency — see #5.)

### 4. Corner-snap node preference (5 m radius)
Per user 2026-05-05: rect corners should prefer a pav.boundary
**vertex** over an arbitrary boundary edge-projection when a vertex
is close, so corners share exact node IDs with pav.

* [pavement/rects.py:_prefer_pav_node](src/auto_patch/pavement/rects.py)
  — helper, 5 m radius.
* Applied in
  [pavement/rects.py:_snap_corners_to_pavement](src/auto_patch/pavement/rects.py)
  and
  [pavement/rects.py:_extend_rect_corners_perpendicular](src/auto_patch/pavement/rects.py).

### 5. Hole-aware sloping-edge snap revised
``_snap_rect_sloping_edges_to_holes`` now validates all 4 corners on
pav.boundary after alignment.  Off-boundary corners snap to nearest
pav node within 5 m; if no node within range, axis is shortened 5 m
and the snap+validate cycle retries (up to 5 iterations).  If all
retries fail, the function falls back to the **original** rect
(pre-alignment) so we never ship a partially-modified rect with
off-boundary corners.

* [pavement/rects.py:_snap_rect_sloping_edges_to_holes](src/auto_patch/pavement/rects.py)
* [pavement/rects.py:_try_align_sloping_to_hole](src/auto_patch/pavement/rects.py)
  — alignment logic factored out for the retry loop.

### 6. ``_extend_rect_corners_perpendicular`` walk-inward fix
When the natural rect corner is OUTSIDE pav (in a void), the function
previously returned ``base_dist`` unchanged, leaving corners 19-22 m
off pav.boundary.  Now walks INWARD until it crosses into pav,
landing on the boundary.

* [pavement/rects.py:_ray_extent](src/auto_patch/pavement/rects.py)
  (inside ``_extend_rect_corners_perpendicular``).

### 7. Long-edge absorption now snaps kept-fragment corners
``_drop_primary_parallels_embedded_in_pavement`` reconstructs corners
along the axis after absorption, which previously left them off
pav.boundary.  Now snapped to pav.boundary; if snap fails (apron-
interior detection: ≥2 corners deep), the fragment is **dropped**
rather than emitted with off-boundary natural corners.

* [pavement/absorption.py:_drop_primary_parallels_embedded_in_pavement](src/auto_patch/pavement/absorption.py)
  — kept-interval reconstruction now calls ``_snap_corners_to_pavement``.

## Critical invariants this baseline maintains

1. **Single boundary source of truth.**  Both rect snapping and
   junction emit consume the same ``pav_union``.
2. **All rect corners on pav.boundary.**  Verified: 0/188 off.
3. **Only terminals are simplified before subtraction in junction
   emit** — rects are aligned by construction (corners snapped to
   simplified pav_union), runway is authoritative apt.dat geometry.
4. **Stub ≠ parallel.**  A rect classified as stub has db ≥ 20° to
   the nearest runway.

## What's NOT in this baseline (intentionally)

The following passes are **disabled** for the geometry baseline so
the geometry can be reviewed in isolation.  A follow-up agent will
re-enable them.

* **Phase-2 elevations**: ``compute_elevations=False`` in build calls.
  No per-vertex altitudes, no Jacobi solver, no within-shape grade
  enforcement.
* **Runway segmentation**: runway emits as 2 polygons (one per
  runway end), not the 50+ segmented chain produced when
  ``compute_elevations=True``.  Re-segmentation depends on
  CIFP-driven splice points and the current widen passes.
* **Junction widening to runway corners**:
  ``widen_junctions_to_runway_corners`` runs only post-elevation in
  pipeline; without phase-2 it doesn't run.
* **Slope alignment / runway 1:1 sharing /
  ``_snap_to_sloping_edge_corners``**: all post-elevation rules in
  ``junction_rules.py`` are skipped.

## Known residual issues at SPJC

1. **5 sliver-corner junction triangles dropped at to_osm.**  Min
   angles 0.5° – 1.9°.  Tiny residue artifacts (3-vert triangles) —
   not user-visible.  Caused by sub-tolerance geometric mismatches at
   rect-rect-pav intersections.  Could be addressed by a residue
   sliver-tip pass; left for follow-up.
2. **F-stub at 71° bearing absorbed away.**  After the long-edge
   absorption pass, the kept fragment of this F-stub had no corners
   snappable to pav.boundary, so it was dropped.  Output has 4 F
   primary_parallel + 0 F stub; user expected 1 F stub.  May or may
   not be a real concern — verify in JOSM.

## Next agent — task spec

**Goal**: re-enable runway segmentation + Phase-2 elevations on the
geometry baseline without regressing the corner-on-boundary invariant.

**Starting point**: tag ``geometry-baseline`` (commit recorded by
``git log`` after this status update).

**What to re-enable**:

1. **Runway segmentation**: in
   [pipeline.py](src/auto_patch/pipeline.py), the runway is currently
   emitted as 2 polygons.  When ``compute_elevations=True`` runs,
   ``run_phase2`` segments the runway based on CIFP breakpoints and
   surrounding pavement boundaries.  Trace where segmentation lives
   (likely in ``finalize.run_phase2`` and the ``_segment_runway*``
   helpers in ``elevation.py``) and ensure it works against the
   simplified pav_union.

2. **Phase-2 elevation pass**: re-run ``run_phase2`` with the
   geometry baseline.  Watch for grade-violation warnings — the
   simplified pav_union may produce slightly different junction
   shapes than the pre-baseline build had.

3. **Junction widening / slope alignment / runway 1:1 sharing**:
   the post-elevation rules in
   [junction_rules.py](src/auto_patch/junction_rules.py) need to
   work against the new geometry.  ``widen_junctions_to_runway_corners``
   in particular has gates the prior baseline tuned (skip-when-aligned,
   multi-step gate, max_total_shared=7, targeted interior prune —
   see ``MEMORY.md`` → feedback_widen_gates) — these may need
   re-tuning.

**Invariants to preserve**:

* All rect corners stay on ``pav_union.boundary``.  Verify with
  the diagnostic in ``/tmp/diag_corners_summary.py`` (regenerate as
  needed).
* No rect's long edge runs alongside a junction polygon (per user
  2026-04-27 invariant).  This is enforced by absorption +
  hole-aware sloping snap.
* Single pav_union — don't reintroduce ``pav_for_rects``.

**Build / verify commands**:

```bash
# Build SPJC with phase-2 elevations enabled
/Users/noah/Ortho4XP-shred86/venv/bin/python -c "
import os, sys
sys.path.insert(0, '/Users/noah/Ortho4XP-shred86/src')
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement('SPJC',
    os.environ.get('XPLANE_ROOT', '/Users/noah/X-Plane 12'),
    compute_elevations=True)
layout.to_osm('/tmp/SPJC_phase2.osm')
print(f'shapes={len(layout.shapes)}')
"

# Run compare-target test
/Users/noah/Ortho4XP-shred86/venv/bin/python -m pytest \
    tests/test_compare_target.py -x

# Verify all rect corners still on pav.boundary
/Users/noah/Ortho4XP-shred86/venv/bin/python /tmp/diag_corners_summary.py
```

**Reference artifacts**:

* [/tmp/SPJC_now.osm](file:///tmp/SPJC_now.osm) — geometry baseline
  (no elevations, no widening)
* [/tmp/SPJC_raw_pav.osm](file:///tmp/SPJC_raw_pav.osm) — raw
  pav_union with only runway + terminals subtracted (visualization
  reference for what pavement looks like)
* ``tests/fixtures/SPJC_target.osm`` — hand-fix target for
  ``test_compare_target_spjc``
