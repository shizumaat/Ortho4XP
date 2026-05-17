# Auto-Patch Status — 2026-05-17 session: single-pass absorption + Lima end + 20→16 failures

## TL;DR

Fourteen commits this session. Tests went **20 → 16 failures**.
Three user-reported visual issues fixed (CYXY E primary parallel
chain absorbed into south apron; CYXY phantom stubs / junction
wrapping eliminated; SPJC Lima end junctions now proper quadrilaterals
sharing 2 runway corners).  The biggest architectural change: replaced
the multi-pass probe-based absorption with a single-pass
shared-sloping-edge absorption at the END of the pipeline, using
actual junction polygons.

**Tests: 248 / 264 pass at HEAD.  16 failures = the remaining backlog.**

## Commits this session (in order)

| Hash | What |
|---|---|
| `989b865` | Runway-stub emit: skip non-terminal per-ref endpoints.  Source B fix: stubs.py was firing at apt.dat graph forks (count ≥ 2 endpoints) that aren't true chain termini.  SPLP phantom A stubs gone. |
| `b3206fb` | Centerline split: validate split points by distinct refs / runway.  Source A: _split_centerlines_at_points was over-splitting on same-name graph forks.  Defensive fix. |
| `dd9be70` | Junction widening body-prox cap.  Pre-cut runway snap nodes already mark each junction's natural extent — walking one more chain step adds corners past the body.  SPLP runway-shared counts now match target exactly. |
| `295ac06` | Primary-parallel partial-absorption: relax reject + restore split call.  Commit a9b0ef0 (prior session) had over-strict `_rect_long_edges_at_pavement_boundary` AND removed the `_split_primary_parallels_at_pavement_boundary` call.  Restored both — CYXY E now emits as 3 sub-rects through the south apron. |
| `a34aaaf` | Drop unrefed phantom stubs (short rects + apron-edge runway stubs).  CYXY junction `-10082` no longer wraps phantom unrefed stubs. |
| `3a32768` | Junction flat-edge snap: pull near-corner verts onto sloping-rect corners.  Conservative tolerance (≤ 2 m perp, ≤ 10 m corner) — only "almost-at-the-corner" cases.  Wired in twice (post-elevation snap + post-split). |
| `62a2242` | **test_taxi_rects_not_alongside_apron**: rewrite with shared-edge rule.  Old probe-based rule (5 m perpendicular point inside polygon) flagged legitimate adjacent rects and required airport-specific exemption.  New rule: junction perimeter linestring overlap ≥ 2 m of rect's sloping edge within 1 m perp.  Universal, no exemptions. |
| `21985bf` | Post-split absorption: fully-shared sloping edge → drop + extend junction.  Initial incremental absorption fix (full-edge-shared only). |
| `e2efcae` | Single-pass sloping-edge absorption at end of pipeline.  Replaces incremental approach — runs ONCE after all post-elevation junction-refinement, with actual junction polygons.  Detects axial t-ranges of sloping-edge sharing, builds absorbed strip polygons, extends absorbing junctions via shapely union, clips / drops the rect.  Corner-wrap filter clamps t to [0.05, 0.95]. |
| `52a6d8a` | Preserve junction node_altitudes when extending via absorption.  Bucket-index old junction corners → altitudes; assign nearest-old-corner altitude to new strip corners; build new node_altitudes for the union polygon's vertex set.  Fixes test_pavement_grade[CYXY] (was using median-altitude flat fallback). |
| `1aa5a5a` | Re-run sloping / flat-edge cleanup after end-of-pipeline absorption.  Clipped sub-rects from absorption can have new corners that don't yet align with adjacent junction vertices; re-run _split_sloped_rects_at_violations + _snap_junction_vertices_to_rect_flat_edge_corners after absorption. |
| `71ca504` | Absorption strip reuses rect corners to preserve neighbour sharing.  For full-absorption (t_lo=0, t_hi=1) reuse r.polygon directly so unary_union preserves the rect's exact corner positions — neighbour rects that shared short edges via those corners keep their connection through the extended junction's perimeter. |
| `f69e395` | Bump runway-pavement intersection tolerance 6 → 12 m.  SPJC Lima end fix: apt.dat row-110 boundary at lat -12.036366 / lon -77.107520 sits 11.3 m perpendicular from the row-100 runway 16L rect (the runway shoulder).  At 6 m the vertex was missed, no runway seam fired at Lima's natural termination, Lima's end junction was a triangle with 1 runway corner. |
| `629043c` | Frame INTERSECTION_PROX_M as RUNWAY_SHOULDER_M (7.6 m) + CHART_TOL_M (4.4 m).  Same numeric value (12 m) but expressed in terms of physical runway geometry.  SPJC 16L/34R declares shoulder surface codes 27/28 in apt.dat row 100 — shoulders exist even without explicit width. |

## Visual issues fixed (user-reported)

1. **CYXY E primary parallel chain absorbed into south apron**
   (user 2026-05-16).  Way ID `-10082` was absorbing E's southern
   primary parallel section due to commit a9b0ef0's over-strict
   long-edge-at-boundary check + missing partial-absorption call.
   E now emits as multi-segment chain through the south apron with
   the corridor extending south to runway 02.

2. **CYXY phantom stubs + junction wrapping eliminated**
   (user 2026-05-16).  Junction `-10081` was wrapping ways
   `-10467` and `-10468`.  Both were unrefed apt.dat fragments
   (one short centerline → `_build_taxi_rects` stub, one
   `_emit_primary_parallel_runway_stubs` from unrefed long
   line at apron edge).  Filtered both paths.

3. **SPJC Lima end junctions now proper quadrilaterals**
   (user 2026-05-17).  South and north Lima end junctions were
   triangles with only 1 runway corner shared.  Root cause: apt.dat
   row-110 boundary vertex at Lima south end sits 11.3 m past the
   row-100 runway rect (the runway shoulder); the
   `pav_runway_intersections` proximity tolerance was 6 m so the
   vertex was missed, no runway seam fired at Lima's natural
   termination.  Bumped tolerance to RUNWAY_SHOULDER_M (7.6 m) +
   CHART_TOL_M (4.4 m) = 12 m.  Both Lima end junctions now have
   4 verts with 2 runway-shared corners.

## Visual eval output

`/tmp/{SPJC,SPLP,CYXY}_current.osm` written at session end.  Verify
in JOSM:

* **SPJC Lima end junctions** (south at lat -12.0365, north at lat
  -12.008): both should be quadrilaterals with 2 runway-shared
  corners (via OSM way IDs `-10193` south, `-10187` north).
* **CYXY taxiway E**: south-end primary parallel chain should
  extend south through the apron to runway 02, with the apron's
  junction polygon (`-10082` area) absorbing only the embedded
  portion.
* **CYXY junction `-10081`**: should NOT wrap any small unrefed
  stub rects — those phantom shapes are gone.

## The remaining backlog — 16 failing tests

User direction (2026-05-17 session-end):

> Continue getting all tests passing and SPJC, SPLP, CYXY rendering
> correctly.

### Remaining test failures

| Test | Detail |
|---|---|
| `test_compare_target_spjc` | `primary_parallel: matched=25 < floor=26; stub: matched=18 < floor=19` — same shortfall as previous session.  Geometry shift from a9b0ef0 lost 1 primary + 1 stub vs target.  Either fix the geometry regression or re-anchor floor (after visual confirmation that current geometry is OK). |
| `test_compare_target_splp[-13--78-baseline1-196]` | Surfaced this session, likely from one of the SPLP-affecting absorption changes.  Look at A primary parallel match counts. |
| `test_junction_boundary_near_centerline[SPJC,SPLP,CYXY]` | Pre-existing structural pattern — junctions whose perimeter strays > MAX_BOUNDARY_TO_CENTERLINE_M from any centerline (apron-territory pavement misclassified as junction).  Per status-comment: "should be re-classified as ``role=apron`` or split." |
| `test_junction_vertex_count_bounded[SPLP]` | 1 violation: `#31(junction) verts=33` (cap 30).  Edge case. |
| `test_junction_vertex_count_bounded[CYXY]` | 5 violations: `#76(junction) verts=103` (the biggest), `#99 verts=49`, `#73 verts=40`, `#92 verts=37`, `#83 verts=32`.  Mostly the big runway-crossing junction and a few near-cap cases. |
| `test_junction_neighbour_corners_shared[CYXY]` | 1 violation: `#76(junction) ⟂ #466(groundside_pavement/groundside) at (-528.5, 558.5) miss=0.55m`.  Pre-existing 0.55 m groundside-vertex orphan from status.md backlog. |
| `test_taxi_rects_not_alongside_apron[SPLP]` | 1 violation: `#1(primary_parallel/A) sloping edge shared 26.1m with #32(junction)`.  Partial-share case the end-of-pipeline absorption doesn't yet clip (the absorbed sub-rect's neighbours' altitudes would need careful re-derivation).  See "Partial-absorption clip" note below. |
| `test_junction_no_long_edge_proximity[SPJC]` | Pre-existing (5–18 m perp distances).  Different pattern from CYXY/SPLP cases that c1e2ad0 fixed. |
| `test_large_junction_axis_aligned_borders[CYXY]` | Pre-existing (96 violations, mostly bearing misalignment 32–41° at junctions around the runway-crossing area). |
| `test_junction_vertices_outside_pavement[CYXY]` | 34 violations (most are inside-pav by 0.4–25 m, not outside).  Test name is misleading; mostly inside-pav verts that need to snap to pavement boundary. |
| `test_no_vertex_on_sloping_rect_edge[SPLP]` | Phenomenon A — untagged clip residue: SPLP `runway(02/20)` 3-corner sliver + 5-corner pentagon.  Source: `_drop_overlap_against_fixed_shapes` (`elevation.py`) producing non-4-corner sloping rects via `_clip_keep_largest`. |
| `test_rect_short_edges_connect[SPLP]` | Cross (short) edge of a stub not connecting to neighbour properly. |
| `test_pavement_grade[SPJC]` | Cross-shape proximity violations from a9b0ef0 fallout.  Geometry shifted shared corners off alignment. |
| `test_pavement_grade[SPLP]` | 38 within-shape grade/plane violations (cap 30).  Worst 30.76 % at junction.  Apron-junction grade rule violation (1.0 % cap exceeded by big factor at one apron). |

### Suggested attack order

1. **CYXY mega-junction `#76` (103 verts)** — the runway-crossing
   junction is far over the 30-vertex cap.  Splitting it or
   reclassifying parts as apron would resolve:
   - `test_junction_vertex_count_bounded[CYXY]` (5 violations,
     #76 worst at 103)
   - `test_junction_vertices_outside_pavement[CYXY]` (34
     inside-pav verts likely in this junction's perimeter)
   - `test_large_junction_axis_aligned_borders[CYXY]` (96
     bearing-misalignment violations, mostly junction #76 edges)
   - `test_junction_neighbour_corners_shared[CYXY]` (junction #76
     vs groundside)
   - `test_junction_boundary_near_centerline[CYXY]`
   These 5 failures all touch the same root junction — fixing
   it is the highest-leverage move on CYXY.

2. **SPLP `_drop_overlap_against_fixed_shapes` clip artifacts**
   (`test_no_vertex_on_sloping_rect_edge[SPLP]`).  Per status.md
   backlog: prevent the clip from producing non-rect output, or
   absorb / demote the post-clip artifacts.

3. **SPLP `#31(junction) verts=33`** vertex cap — 3 over the
   limit, probably a localized issue.

4. **test_compare_target_*** (SPJC, SPLP) — once visual geometry
   confirmed correct, re-anchor target fixture floors.

5. **SPLP A 26 m partial share**
   (`test_taxi_rects_not_alongside_apron[SPLP]`) — needs partial
   absorption support in `_absorb_rects_at_junction_perimeters`.
   Currently the function detects both full and partial sharing
   in t-range computation but only fully-drops or fully-keeps
   based on whether ALL of [0, 1] − absorbed is < min_kept.
   For SPLP A: 26 m absorbed of 1000 m total → kept range = ~974 m,
   well above min_kept, so SHOULD be clipped.  Investigate why
   it doesn't fire — possibly the partial-sharing detection's
   corner-wrap clamp ([0.05, 0.95]) excludes the 26 m portion
   because it's near a corner.

## Architecture: single-pass absorption (post-emit)

`_absorb_rects_at_junction_perimeters(layout, icao)` in
`src/auto_patch/junction_repair.py:1631`.  Runs at the END of
the pipeline (`src/auto_patch/pipeline.py:2287` area, after
`_split_sloped_rects_at_violations` +
`_snap_junction_vertices_to_rect_flat_edge_corners`).

For each sloping rect:
  1. Get 2 sloping edges (corners 0-1 and 2-3 per
     `_rect_from_axis_extended` convention).
  2. For each junction polygon, intersect `junction.boundary` with
     `sloping_edge.buffer(perp_tol_m=0.5)`; project shared
     LineString endpoints to axial t-range.
  3. Clamp t-range to [0.05, 0.95] to exclude junction-wrap-at-
     short-edge-corner false positives.
  4. Merge all absorbed t-ranges; compute kept = [0, 1] − absorbed.
  5. Filter kept by `min_kept_m = 20`; fold short kept ranges into
     absorbed.
  6. Build absorbed strip polygons; extend each absorbing
     junction via `unary_union([junction.polygon, strip])`.
  7. Build new sub-rect from each kept t-range; drop original.

Strip-building (`_strip_polygon`): for full-absorption uses
`r.polygon` directly so the union preserves the rect's exact
corner positions (so neighbour rects keep their short-edge
connections through the extended junction's perimeter).

Altitude handling on extended junctions: bucket-index old corners
→ altitudes; for new vertices use nearest-known-corner altitude.

## Key memory files (read before continuing)

* `feedback_general_solutions.md` — every fix must work at all
  baseline airports.
* `feedback_root_cause_only.md` — fix root causes; ASK before
  band-aid post-process clean-up.
* `feedback_shape_rules.md` — authoritative rect + junction
  construction rules.
* `feedback_grade_rules.md` — apron / junction grade rule.
* `project_refactor_state.md` — extraction pattern.

## How to verify / reproduce

```bash
# Full test suite (~3 min):
venv/bin/python3 -m pytest tests/ --tb=short -q

# Single-airport build for visual inspection:
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
for icao in ("SPJC", "SPLP", "CYXY"):
    layout = build_airport_pavement(
        icao, "/Users/noah/X-Plane 12", compute_elevations=True)
    layout.to_osm(f"/tmp/{icao}_current.osm")
PY
```

## Files touched this session

Source:
* `src/auto_patch/pavement/stubs.py` — count-1-only filter for
  runway-stub emit (Source B).
* `src/auto_patch/pavement/centerlines.py` — split-point
  intersection validator (Source A).
* `src/auto_patch/junction_rules.py` — body-prox widening cap +
  conservative flat-edge corner snap.
* `src/auto_patch/pavement/rects.py` — drop unrefed phantom
  stubs + relax `_rect_long_edges_at_pavement_boundary` to
  reject only when BOTH long edges majority-embedded.
* `src/auto_patch/pavement/absorption.py` — restore
  `_split_primary_parallels_at_pavement_boundary` (just call
  re-added in pipeline; function unchanged).
* `src/auto_patch/junction_repair.py` — added
  `_absorb_rect_into_junction`,
  `_drop_rects_with_shared_sloping_edge_and_absorb` (initial
  incremental), then
  `_absorb_rects_at_junction_perimeters` (single-pass at end).
* `src/auto_patch/pipeline.py` — wired the absorption call, the
  re-run-cleanup-after-absorption, the flat-edge snap, the
  partial-absorption restore, and bumped `INTERSECTION_PROX_M`
  to `RUNWAY_SHOULDER_M + CHART_TOL_M`.

Tests:
* `tests/test_junction_invariants.py` — rewrote
  `test_taxi_rects_not_alongside_apron` with shared-edge rule;
  removed `TAXI_RECT_ADJACENCY_REGRESSION_BASELINE`.

## Build commands (unchanged)

```bash
venv/bin/python3 -m pytest tests/ --tb=short -q

# Single-airport build for visual inspection:
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
for icao in ("SPJC", "SPLP", "CYXY"):
    layout = build_airport_pavement(
        icao, "/Users/noah/X-Plane 12", compute_elevations=True)
    layout.to_osm(f"/tmp/{icao}.osm")
PY
```
