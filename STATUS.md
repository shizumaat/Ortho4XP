# Auto-Patch Status — 2026-05-12: TILE-CUT / DIAGONAL STUB / SLOPING-EDGE REFINEMENTS

## TL;DR

**9 commits since the last status.md** (commit `8921aa6`, the
2026-05-11 status entry).  207 of 208 tests pass.  Only the
pre-existing `test_compare_target_spjc` fixture failure remains
(the fixture is OSM-derived; the build switched to apt.dat in the
previous session, so spatial-match positions differ by 1-3 m).

This session was driven by user-directed refinement of the
auto-patch output at SPJC, working through visible issues in
`/tmp/SPJC_v*.osm` outputs.  Five themes:

1. **Cross-tile elevation determinism + regression test** (commit
   `cc3c312`).  Found and fixed a coord/DEM mismatch in
   `per_surface_solve` + `_compute_elevations` for cross-tile
   airports.  Added [tests/test_tile_cut_parity.py](tests/test_tile_cut_parity.py)
   that builds SPLP twice (anchor tile + neighbour tile) and
   asserts elevations agree at the cut (pre-fix worst |dz| =
   4.80 m, post-fix < 1 m).

2. **Phase B asymmetric dedup → removed entirely** (commits
   `05f03f9`, `9bd94b1`).  `stitch_pavement_to_flat_runways`'s
   Phase B perpendicular-projection pass was creating arbitrary-
   looking corners on flat runway segments at perpendicular
   projections of interior junction vertices.  After
   confirming zero grade-test impact when disabled, removed
   entirely (200 lines deleted).  Junction interiors anchored
   to runway altitude via Phase A snaps + per-surface solver.

3. **Sloped end-segment split at blast-pad pav vertex** (commit
   `5bab14c`).  pipeline.py's `pav_runway_intersections`
   capture used row-100 frame for `t` but the rect_boundary it
   tests against uses the blast-extended frame; the mismatch
   filtered out blast-pad-zone pavement vertices.  Extending
   `cl_a/b` by `blast_a/b_m` in pipeline.py captures the SPJC
   V1-throat row-110 vertex (1.79 m off runway boundary, 5.27 m
   into the blast pad).  The runway segmenter then splits the
   sloped end-segment at that fraction, giving the apron
   junction an inline 4-corner sub-rect to snap to.

4. **Diagonal-stub geometry overhaul** (commits `5404f25`,
   `e4f8a9a`, `4f1e0c9`, `ba10a41`):
   * Flipped the 20 % bias from TOWARD-runway to AWAY-from-
     runway, then replaced bias entirely with a centered 30 m
     fixed margin on each end of the gap.
   * Extended the 30 m runway-buffer pull-back from
     perpendicular-only to all centerlines with
     `perp_diff < 70°` (diagonal stubs were skipping it).
   * Moved the diagonal-parent + sub-ref STUB rules in
     `_classify_role` BEFORE the `db<20°` parallel branch so a
     curving sub-segment of a diagonal taxi (SPJC C had a 75 m
     near-parallel slice at `db_local = 19.2°`) doesn't slip
     through as PRIMARY_PARALLEL.  Added a complementary
     "single-rect" dedup for diagonal-parent stub refs.
   * Post-process drop of thin orphan sliver junctions adjacent
     to stubs/parallels.  SPJC 3 → 0 thin slivers.

5. **Sloped-rect split at junction-vertex sloping-edge
   violations** (commit `41c50cd`).  When a junction polygon has
   a vertex on the sloping (long) edge of a sloped 4-corner
   rect, the rect splits at that vertex's axial position into
   two sloped sub-rects (both 4-corner, altitude_high /
   altitude_low linearly interpolated).  The junction vertex
   now coincides with a sub-rect's short-edge corner instead of
   sitting on a sloping edge — topology violation fixed.  SPJC
   1 → 0 violations.

**Branch:** `dev` (56 commits ahead of `origin/dev`).
**Tag baseline:** `spjc-good-baseline` (commit `4d97f89`,
2026-05-07).

## Commits this session

```
41c50cd Split sloped rects at junction-vertex sloping-edge violations
4f1e0c9 Run diagonal-parent + sub-ref STUB rules before db<20° parallel check
ba10a41 Drop thin orphan sliver junctions adjacent to stub/parallel rects
e4f8a9a Diagonal stub: 30 m fixed margin + 30 m runway-perp buffer trim
5404f25 Flip diagonal-stub margin bias: AWAY from runway, not toward
cc3c312 Fix coord/DEM mismatch in per_surface_solve + _compute_elevations
9bd94b1 Remove Phase B perpendicular-projection from stitch pass
05f03f9 Symmetric Phase B insert dedup in stitch_pavement_to_flat_runways
5bab14c Capture blast-pad pav vertices for sloped end-segment split
```

## Build / verify commands

```bash
# Full test suite (207/208 pass; pre-existing SPJC compare-target fail)
venv/bin/python3 -m pytest tests/ --tb=short -q

# Build SPJC + write OSM for visual review
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement(
    "SPJC", "/Users/noah/X-Plane 12", compute_elevations=True)
layout.to_osm("/tmp/SPJC.osm")
print("wrote /tmp/SPJC.osm")
PY

# Cross-tile regression test (needs SPLP scenery + S13W077 + S13W078 hgt)
venv/bin/python3 -m pytest tests/test_tile_cut_parity.py -v
```

## Geometry state at SPJC v23 (current dev tip)

**Diagonal-stub min-distance-to-runway-edge (m):**

| ref | d_rwy_edge | notes |
|-----|-----------|-------|
| A   | 62.5      | F-style loop ramp |
| B   | 25.7      | diagonal stub |
| C   | 19.8      | diagonal stub |
| D   | 55.8      | apron-internal |
| E   | 20.4      | diagonal stub |
| F   | 71.4      | F-style loop ramp |
| G   | 25.6      | diagonal stub |
| L (×3) | 34.8, 54.8, 136.5 | parallel sub-pieces |
| L3, L5 | 40.8, 52.0 | sub-refs |
| V1, V2, V3, V5 | 44.2, 57.5, runway-crossing, 40.3 | sub-refs |

**Junction count:** 38 (target 44; v20 was 38).
**Multi-node flat runways:** 0.
**Thin orphan slivers:** 0.
**Sloping-edge rule violations:** 0.

## Investigated and abandoned this session

### Width-transition splitting (option A from user's option D)

Implemented [`_find_width_transition_breakpoints`](src/auto_patch/pavement/centerlines.py:769)
(probe perpendicular half-width along each centerline, split at
narrow ↔ wide transitions using widen_factor=1.5).  Successfully
brought SPJC V parallel count from 3 → 5 (matching target's 5).
**But** the new V rects sit in regions target treats as junction
territory; their sloping edges then introduced 33 mid-edge step
violations against adjacent junctions (failing the pavement-grade
test, cap 10).

Disabled in [pipeline.py](src/auto_patch/pipeline.py:1271)
behind a `pass` placeholder, comment block documenting the
finding.  Function `_find_width_transition_breakpoints` remains
exported in case future work needs it.

### Cross-taxi geometric intersection splitting (option B)

Already in the pipeline at
[pipeline.py:1243-1261](src/auto_patch/pipeline.py:1243).
Inspected at length — SPJC's missing V splits (Issue 2 below)
are NOT at any apt.dat-detectable cross-taxi intersection.
Target's splits there come from OSM way fragmentation that
apt.dat doesn't capture.  No improvement available without
changing centerline data source.

## OPEN WORK — where to pick up

### Hot topic: TILE-BOUNDARY SLICE WIDTH 5 m → 0.5 m (USER'S NEXT TASK)

**Current:** [`tile_cut.py:45`](src/auto_patch/tile_cut.py:45)
defaults `half_width_m: float = 5.0`, producing a **10 m gap**
along every integer-lat/lon line passing through the airport
footprint.

**Target:** reduce to `half_width_m = 0.5`, producing a **1 m
gap**.  Tighter gap means less pavement lost to the cut.

**Files to touch:**
* [`src/auto_patch/tile_cut.py:45`](src/auto_patch/tile_cut.py:45)
  — change the default.
* Verify with the [`test_tile_cut_parity.py`](tests/test_tile_cut_parity.py)
  regression test — the test currently uses `NEAR_CUT_M = 15`
  to identify near-cut vertices; tighten if needed to reflect
  the new 1 m gap geometry.
* `cut_polys = [line.buffer(half_width_m, cap_style=2)`
  at [`tile_cut.py:111`](src/auto_patch/tile_cut.py:111).
  Tighter buffer = smaller geometric difference = more
  pavement retained, but also less margin for
  floating-point seam alignment between adjacent tiles.
* Check that the cut still produces clean shape splits at
  0.5 m (shapely's polygon difference with a thin buffer can
  produce sliver geometries that get filtered by
  `min_piece_area_m2: float = 1.0`).

**Validation:**
* Run `tests/test_tile_cut_parity.py` for SPLP cross-tile parity.
* Build SPJC + SPLP and visually inspect the cut seams.
* Ensure no new sliver-merge or sliver-drop log lines flag
  unexpected drops.

### Known limitations (not addressed this session)

* **Issue 2 (v23 -10169 junction has 650 m straight edge across
  grass):** SPJC's V parallel is missing a rect at ax 2090..2394
  (target has -10012 there).  Target splits V there based on OSM
  way fragmentation; apt.dat row-1202 has a single long edge
  through that range with no cross-taxi geometric intersection.
  Width-transition splitting (the user's option A) could fix
  this but breaks the pavement-grade test (new V rects in
  target's junction territory produce 33 mid-edge violations).
  Documented in the commit message for `41c50cd`.

* **`test_compare_target_spjc` fixture mismatch:** the SPJC
  target fixture was generated from OSM-derived geometry before
  the apt.dat-primary switch.  Position differences of 1-3 m
  break the spatial-match floor counts.  Regeneration deferred
  until the diagonal-stub + V-parallel geometry is final.

* **Within-shape grade WARN at SPJC:** 1 violation worst-case
  8.0 % on junction -10159 (35.5 → 33.7, d=22.4 m, de=1.8 m).
  Pre-existing and consistent across sessions; not yet
  investigated.

* **Exception-hardening item 3 (project_todos):** broad-except
  pass is done within `src/auto_patch/` but the wider Ortho4XP
  codebase outside `auto_patch/` may still have similar issues.
  Out of scope unless requested.

## Reference artefacts

* `/tmp/SPJC_v23.osm` — most recent SPJC build (commit `41c50cd`).
* `/tmp/SPJC_v22.osm` — width-split experiment (abandoned).
* `/tmp/SPJC_v20.osm` — pre-classifier-fix state.
* `tests/fixtures/SPJC_target.osm` — OSM-derived baseline gate
  (pre-apt.dat-switch).  Regeneration deferred.
* The un-tracked `+60-140/` directory at the repo root: 24 stray
  `.hgt` elevation-tile files that should live under
  `Elevation_data/+60-140/`.  Per user direction: leave
  untracked, never include in `git add -A`.

## Memory notes for next agent

Saved under `~/.claude/projects/-Users-noah-Ortho4XP-shred86/memory/`:
* No new memory entries this session — all decisions captured in
  commit messages.  Existing memories `feedback_root_cause_only`,
  `feedback_general_solutions`, `feedback_shape_rules`,
  `project_target_osm` remain authoritative.
