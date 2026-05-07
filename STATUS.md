# Auto-Patch Status — PHASE-2 RE-ENABLED + ELEVATION POLISH 2026-05-07

## TL;DR

Phase-2 is **re-enabled and validated** on top of the
``geometry-baseline`` tag from 2026-05-05.  Runway segmentation,
elevation solver, junction widening, slope alignment, and runway
1:1 sharing all run end-to-end.  Several over-aggressive
simplification passes were disabled / removed.  ``compare-target``
now passes for the first time this branch.

One known visual issue **remains unresolved** on SPJC: a 100k m²
junction polygon over the eastern apron sits over an unfiltered
DEM spike (~9 m above ring elevations, 85k m² footprint), producing
"sharp hills and drops" in X-Plane.  Multiple subdivide approaches
were prototyped during this session (4-way runway-aligned split,
DEM-aware 3-way slice with flat ends + sloping rect middle, raw-DEM
override) — all left in-tree, but **none successfully eliminated
the bumps when tested in X-Plane**.  Spike-slice machinery is
disabled at the call site; the function is preserved in
``junction_repair.py`` for follow-up.  Everything else looks good
in JOSM and in flight.

**Branch:** ``dev``.  **Tag baseline:** ``geometry-baseline``
(commit ``a686e5e``).

**Build artifact:** [/tmp/SPJC_phase2_v17.osm](file:///tmp/SPJC_phase2_v17.osm)
— SPJC with ``compute_elevations=True``, 276 ways, 41 junctions,
phase-2 elevations applied.

## Final SPJC numbers (this session)

```
junctions emitted     : 41   (5 sliver triangles dropped at to_osm)
rects (taxi)          : 50
terminals             : 2
runway segments       : 80   (segmented chain via finalize.run_phase2)
retaining_wall        : 66
tunnel_ramp           : 36
total OSM ways        : 276
rect corners on pav   : 200 / 200  (worst d = 0.000 m)
within-shape grade    : 0 WARN
per-surface solver    : 217 iters / ~7 s
compare-target test   : PASS  (was failing F-stub residual all session)
```

## Changes this session

### 1. Phase-2 re-enabled end-to-end
Verified that ``pipeline.py`` builds correctly with
``compute_elevations=True`` against the geometry baseline.  Order
of operations restored:
``run_phase2`` → ``_align_rect_slope_to_axis`` →
``_snap_to_sloping_edge_corners`` → ``_enforce_runway_1to1_sharing``
→ ``widen_junctions_to_runway_corners`` → ``per_surface_solve`` →
``stitch_pavement_to_terminals`` → ``_report_within_shape_violations``.
All passes converge cleanly.

### 2. Removed apt.dat-arc 4-vertex cap
[pavement/junctions.py:1252](src/auto_patch/pavement/junctions.py:1252)
— ``_build_junction_constructive`` no longer subsamples apt.dat
boundary vertices between rect-corner pairs to 4.  Sharp curves
now preserve every apt.dat boundary vertex.  ``max_arc_vertices``
parameter removed from the function signature.

### 3. Disabled long-edge densification
[triangulation.py:660](src/auto_patch/triangulation.py:660) — the
``_densify_long_boundary_edges`` call inside the per-junction
post-elevation cleanup is **disabled**.  The geometry-baseline
boundary trace is precise enough that synthetic midpoints just
add free vertices the per-surface solver pushes around.  An
area-gated experimental re-enable for >50k m² polygons was tried
mid-session and reverted (didn't measurably reduce the SPJC bumps).

### 4. Disabled colinear-vertex pruning
[triangulation.py:497](src/auto_patch/triangulation.py:497) — the
``_drop_colinear_boundary_vertices`` call is **disabled**.  At its
prior 3 m perpendicular threshold it was flattening real apt.dat
curves; the geometry-baseline pav_union no longer produces the
ear-clip slivers this pass was added to mop up.  Total junction
vertex count fell from 1075 (with both passes on) to 916 (with
both off), with sharp curves preserved exactly.

### 5. Spike subdivide infrastructure (DISABLED)
A full DEM-aware subdivide pass was developed in
[junction_repair.py:_subdivide_dem_spike_junctions](src/auto_patch/junction_repair.py:675).
For each junction whose footprint contains DEM samples >5 m above
its ring max (≥5% of in-polygon samples), the function:
* Computes the elevation gradient direction from ring vertices
  (min-alt → max-alt).
* Cuts perpendicular to the gradient at 25% / 75% positions along
  the axis.
* Builds a 4-corner ``primary_parallel`` rect for the middle 50%
  (clipped to the original polygon for non-convex cases via
  longest-chord-per-cut).
* Pins the LOW and HIGH flat-end altitudes to neighbour-shape
  pin altitudes (rect / runway / terminal corners shared with
  the end region) — falls back to mean of end-region vertex
  elevations when no pins exist.
* Emits LOW and HIGH ends as ``role=junction`` with single
  ``altitude=`` tag, MIDDLE as ``role=primary_parallel`` with
  ``altitude_high`` / ``altitude_low`` matching the flats.

The function works correctly when called.  But the slice CALL is
**disabled** at [pipeline.py:1897](src/auto_patch/pipeline.py:1897)
because the resulting OSM, when rendered in X-Plane, did not
visibly improve the rough-spot bumps the user is trying to fix.
Earlier 4-way runway-aligned variant of the same function was
also tested and replaced — both implementations are visible in
git history (commits before the disable comment).

### 6. Grade-based subdivide re-enabled with lower threshold
[junction_repair.py:349](src/auto_patch/junction_repair.py:349) —
``SUBDIVIDE_VIOLATION_GRADE`` lowered from 10% → **2%** so the
existing ``_subdivide_violating_junctions`` pass can fire on the
1.5–3% violations the per-surface solver leaves on long polygons.
Wired into the per-surface path at
[pipeline.py:1908](src/auto_patch/pipeline.py:1908) (the legacy
unified-solver call site only).  Currently fires zero times on
SPJC because the per-surface solver eliminates ring-pair
violations on its own.

### 7. DEM source — kept as Ortho4XP smoothed
``driver.py`` continues to forward Ortho4XP's
``smooth_raster_over_airports``-processed ``tile.dem`` into
``build_airport_pavement``.  A raw-DEM override was tested on
2026-05-07 and reverted — didn't fix the rough spot and risked
re-introducing terrain artefacts the smoothing was added to remove.

## Investigation notes for the unresolved SPJC apron bumps

The polygon in question is the eastern apron at SPJC, ~100k m²,
ring altitudes 32.4–34.6 m.  Sampled DEM under its footprint:
* p99 = 44.16 m
* max = 44.86 m
* 14.9% of samples >5 m above the ring max.
* Spike footprint ~85k m², 278 m × 308 m, centred at
  ``(-77.1045, -12.0335)``.
* This is almost certainly an unfiltered building (hangar /
  terminal) Ortho4XP's airport smoothing didn't fully flatten
  for this specific airport.

The polygon's perimeter is grade-compliant (max all-pairs ring
grade = 1.57%, 1 marginal pair at ~25 m).  Bumps come from
**Triangle4XP's interior Steiner vertices** which are seeded from
the DEM at .dsf-assembly time and aren't visible to auto_patch's
ring-only checks.

Approaches that DIDN'T work in X-Plane testing:
* **DEM-spike 4-way runway-aligned cut** — produced visually
  similar bumps; sub-polygons still large enough for Steiners.
* **DEM-spike 3-way gradient slice** (flat ends + sloping rect
  middle, with pin-aware altitudes) — produced clean OSM, but
  didn't change the rendered bumps.
* **Raw DEM override** — bumps unchanged or worse; smoothing
  was apparently doing useful work elsewhere.
* **Re-enabling densification only for big polygons** — added
  more ring anchors but didn't tighten the per-axis Steiner
  constraint enough.

Other findings worth carrying forward:
* The polygon is functionally an **apron** (large parking area)
  but the auto_patch classifier tags it as **junction** (residue
  catch-all role).  Apron tagging would apply the tighter 1.0%
  cap and apron-flatness anchors and likely solve the bump
  problem on its own.  Worth investigating
  ``pavement/classifier.py`` for why this polygon falls through
  to junction.
* Manual user fix to apt.dat (split overlapping aprons) on
  2026-05-06 fixed a related rough spot but not this one.
* The user's working theory is that earlier builds rendered this
  area correctly.  No git change directly explains the
  regression.  The smoothed-vs-raw DEM toggle made no visible
  difference here.

## What's NOT in this baseline (intentionally)

* DEM-spike subdivide function (disabled at the call site).
* Long-edge densification for big junctions.
* Colinear-vertex pruning (3 m threshold was too aggressive).
* Raw-DEM source for elevations (kept the smoothed Ortho4XP
  ``tile.dem`` after the test was inconclusive).

## Next agent — task spec

**Goal:** eliminate the SPJC eastern-apron rough spot in X-Plane
without regressing other airports.

**Starting point:** this commit on ``dev``.  Build artifact:
[/tmp/SPJC_phase2_v17.osm](file:///tmp/SPJC_phase2_v17.osm).

**Hypotheses worth chasing (in priority order):**

1. **Reclassify the polygon as apron, not junction.**  Check
   ``pavement/classifier.py`` and the apron extraction in
   ``terminals.py`` to understand why this 100k m² apron-like
   polygon falls through to junction.  The unified Jacobi solver
   applies APRON_MAX_GRADE (1.0%) and an apron-flatness anchor
   for true apron polygons — both would fight the spike harder
   than the current taxi-grade-cap (1.5%) treatment.

2. **Triangle4XP Steiner constraint.**  The auto_patch builder
   doesn't directly control Steiner placement; X-Plane / mesh
   builder inserts them based on its own quality refinement.
   But auto_patch CAN influence outcomes by choosing polygon
   tags carefully.  ``primary_parallel`` rects with altitude_high
   / altitude_low are linearly interpolated between the two
   short ends — Steiners on those polygons get the linear
   interpolation, NOT the DEM.  If we could tag this polygon as
   ``primary_parallel`` (or build several smaller ones) with
   appropriate slope axes, the rendered surface would be slope-
   compliant by construction.  The DEM-spike slice prototyped
   this idea but didn't visibly help — possibly because of how
   X-Plane handles non-rectangular ``primary_parallel`` polygons.

3. **DEM-spike override at .dsf time.**  Investigate whether
   Ortho4XP / Triangle4XP can be told to use auto_patch's
   ring-derived altitudes for Steiners INSIDE patched polygons
   (instead of the raw DEM).  This would be the architecturally
   correct fix — the patch should fully describe the surface,
   not just its boundary.

4. **Per-airport ``smoothing_pix`` tuning.**  SPJC's
   ``dico_airports['SPJC']['smoothing_pix']`` value (or
   ``tile.apt_smoothing_pix`` global) may need a larger kernel
   for this airport so ``smooth_raster_over_airports`` actually
   flattens the building under the apron.

**Re-enable instructions for the disabled spike machinery:**

```python
# In pipeline.py around line 1897, replace the disable stub with:
spike_dem = _load_airport_dem(layout.anchor[0], layout.anchor[1])
n_dem_spike = _subdivide_dem_spike_junctions(
    layout, spike_dem, tile_lat, tile_lon,
    runway_axis_deg=runway_axis_deg)
```

**Build / verify commands:**

```bash
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

/Users/noah/Ortho4XP-shred86/venv/bin/python -m pytest \
    tests/test_compare_target.py -x

/Users/noah/Ortho4XP-shred86/venv/bin/python /tmp/diag_corners_summary.py
```

**Reference artifacts:**

* [/tmp/SPJC_phase2_v17.osm](file:///tmp/SPJC_phase2_v17.osm) —
  current build (slice disabled).
* [/tmp/SPJC_phase2_v15.osm](file:///tmp/SPJC_phase2_v15.osm) —
  build with the 3-way gradient slice enabled (for comparison).
* ``Patches/-20-080/-13-078/SPJC_auto.patch.osm`` — Ortho4XP
  runtime output.

---

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
