# Auto-Patch Status — 2026-05-16 (session 2 end): canonical-node bridge rewrite + per-vertex ribbon perp + 23→20 failures

## TL;DR

This session continued the Tier 1 work from the 2026-05-16 morning
handoff (`4ae3230`).  Five commits landed taking the test suite
from **23 → 20 failures**; the bridge feature was rewritten to
build polygons from canonical nodes (boundary + pavement_union
vertices) instead of `parallel_offset()` + layered `difference()`;
the airport-boundary ribbon was fixed to share corners between
adjacent rects in the chain.

**Committed this session (in order):**

1. `7b65cd4` — *Exempt flat shapes from sloping-rect 4-corner
   invariant test.*  Per user 2026-05-09 flat-shape rule: a single
   `altitude=` tag legitimately allows variable node count; the
   test was over-strict on flat sub-rects (CYXY 6-corner runways
   were flat sub-rects from `_split_sloped_rects_at_violations`).
2. `6b61019` — *Accept flat-rect corners as valid snap targets in
   Rule 2 test (`test_junction_no_long_edge_proximity`).*  Same
   rationale — flat-rect corners are equally valid as sloped-rect
   corners for junction-vertex coincidence.  SPLP Rule 2
   violations 3 → 1 (the remaining 1 was the split-but-no-move
   pattern fixed in c1e2ad0).
3. `c1e2ad0` — *Pull junction vertex onto new sub-rect corner
   during sloping-edge split.*  `_split_sloped_rects_at_violations`
   projects the violating junction vertex onto the rect's long
   edge to choose the split parameter and creates a new sub-rect
   corner there — but the junction vertex itself was never moved.
   Result: junction vertex sat ~on_edge_tol_m perpendicular off
   the new sub-rect's long edge AND ~on_edge_tol_m from its
   corner (Tests 2 & 3's same-geometry pattern).  Tracks
   `(j_idx, v_idx)` through clustering and applies the move per
   junction polygon.  SPLP+CYXY Test 2 fully passes; CYXY Test 3
   14 → 10 violations.
4. `3a21c89` — *Rewrite `_emit_boundary_dem_bridge` to use
   canonical nodes.*  Replaces `parallel_offset()` + 5+ layered
   `difference()` ops with a sequential walk: outer side offsets
   the boundary inward by `strip_half_width_m`; inner side walks
   `pavement_union`'s outer ring backward from the bridge run's
   end vertex toward the start, stopping when more than
   `runway_clamp_radius_m` (400 m) from the start.  Every bridge
   vertex is canonical (boundary inner-offset or pavement_union
   outer-ring vertex).  Final small `difference()` cleanup
   against non-runway pavement + ribbon + sloping-rects-buffered
   (1 m, per B12) trims sub-metre residual overlap.  Replaces
   the `_clip_boundary_bridges_against_pavement` re-clip with a
   re-emit after `_split_sloped_rects_at_violations`.  CYXY
   Test 3 10 → 3 violations (remaining 3 are groundside + stub,
   not bridge — different categories).
5. `586b135` — *Boundary ribbon: use per-vertex perpendiculars in
   rect chain.*  Visual report: adjacent boundary rects in the
   chain didn't share their cross-edge nodes (each rect used its
   own per-segment perp; perp at the shared boundary vertex
   differed between adjacent rects).  Switch to per-vertex
   perp (average of incoming + outgoing segment perps) so
   adjacent rects share both inner and outer corners at every
   boundary vertex exactly.

**Tests: 244 / 264 pass at HEAD.  20 failures = the remaining backlog.**

## Visual eval output

`/tmp/{SPJC,SPLP,CYXY}_current.osm` written at session end (after
`586b135`).  Verify visually in JOSM:
* **Boundary ribbon**: every consecutive pair of `airport_boundary`
  shapes should share its cross-edge nodes exactly.  No more
  gaps/overlaps at bends.
* **CYXY bridge polygons**: 2 canonical-node bridges, inner edge
  hugs junction/terminal/rect perimeters; outer edge sits at the
  ribbon's interior side (2.5 m inward of boundary line).  No
  bridge sweeps across pavement.

## The remaining backlog — 20 failing tests

User direction (2026-05-16, session-end):

> Next task is to continue getting all tests passing, and SPJC,
> SPLP, and CYXY airports generating correctly.

So: continue Tier 1 work, drive failures to 0, and verify visually
that each baseline airport renders correctly.

### Remaining Tier 1 invariant failures

`test_junction_neighbour_corners_shared` at all 3 airports still
fails — but the violations are NOT bridge-related anymore (the
rewrite eliminated all bridge orphans at CYXY).  Remaining
patterns:

| Airport | Count | Categories |
|---|---|---|
| **SPJC** | 6 (cap 5) | 5 vs `airport_boundary` (miss 1.4–7 m); 1 vs `runway/16L/34R` (miss 18.76 m, outlier) |
| **SPLP** | 6 (cap 0) | All 6 vs `stub/A` sub-rect corners (miss 4.96–8.51 m) — twin-corner orphans from `_split_sloped_rects_at_violations` (c1e2ad0 fixed the source-vertex side; the OTHER long-edge corner has no matching junction vertex) |
| **CYXY** | 3 (cap 0) | 1 vs `groundside_pavement` (0.55 m); 2 vs `stub` at (-333.5, 194.3) miss 22.48 m — pre-existing structural pattern around the runway crossing |

**Suggested next direction** for Test 3:
* **SPLP twin-corner**: extend `_split_sloped_rects_at_violations`
  to also handle the OPPOSITE-side corner — when the source
  vertex moves to one new sub-rect corner, find adjacent
  junctions whose perimeter passes near the TWIN corner and
  insert/snap there too.  Scoped insertion (only at split-
  introduced corners, not arbitrary neighbours) avoids the
  cascade-and-baseline-bloat of the reverted broad insertion
  pass (`bf37313` reverted in `0c733ba`).
* **CYXY groundside (0.55 m)**: junction#76 vertex sits 0.55 m
  from groundside_pavement#463 vertex at (-528.5, 558.5).  A
  terminal vertex 0.05 m from junction#76 was already
  coincident.  Multi-shape near-coincidence at a single
  junction corner.  Either snap groundside vertex to junction
  (changes groundside shape slightly) or accept via baseline.
* **CYXY 22 m stub orphans**: pre-existing pattern around the
  runway crossing.  Needs targeted diagnosis.
* **SPJC outliers**: 18.76 m miss vs runway/16L/34R is a far
  outlier (junction vertex placed 18 m off any runway corner).
  Pre-existing pre-`a9b0ef0` (per session-1 bisect).

### Other remaining failures

| Test | Airports | Pattern |
|---|---|---|
| `test_compare_target_*` | SPJC, SPLP | Target-fixture rect-match floors not met.  Per `a9b0ef0` commit message: SPJC stub matched=18<floor=19, primary_parallel matched=25<floor=26.  Geometry shift from CYXY-focused changes.  Either re-anchor fixture floors (after confirming shifts are improvements) or locate / fix the regression. |
| `test_junction_no_long_edge_proximity[SPJC]` | SPJC | 6 violations (pre-existing pre-`a9b0ef0`, 5–18 m perp distances).  Different pattern from CYXY/SPLP cases that c1e2ad0 fixed.  Separate diagnosis needed. |
| `test_no_vertex_on_sloping_rect_edge` | SPJC, SPLP | Phenomenon A — untagged clip residue: SPJC `stub(C)` n_corners=5 (overlap-clip pass produced non-rect); SPLP `runway(02/20)` 3-corner sliver + 5-corner pentagon.  Source: `_drop_overlap_against_fixed_shapes` (`elevation.py`) producing non-4-corner sloping rects via `_clip_keep_largest`.  Fix: either prevent clip from producing non-rect output, or absorb/demote-role the post-clip artifacts. |
| `test_pavement_grade[SPJC]` | SPJC | 2 cross-shape proximity violations (`a9b0ef0` fallout).  Geometry shifted shared corners off alignment. |
| `test_pavement_grade[SPLP]` | SPLP | 26 mid-edge steps > 0.5 m (cap 20), worst 1.90 m at junction -10034 / stub A boundary.  `a9b0ef0` fallout. |
| `test_taxi_rects_not_alongside_apron[SPLP/CYXY]` | SPLP, CYXY | Junction residue around runway-crossing area (CYXY) or sub-rect E-D area; stub diagonals flanked by junction pavement. |
| `test_junction_vertex_count_bounded[SPLP/CYXY]` | SPLP, CYXY | Junctions exceed vertex count caps.  May correlate with split-sloped-rect cascades. |
| `test_junction_boundary_near_centerline[SPLP/CYXY]` | SPLP, CYXY | Junction polygon's boundary too close to a taxi centerline (would render as visible overlap). |
| `test_large_junction_axis_aligned_borders[CYXY]` | CYXY | Junction polygon long edges aren't aligned to adjacent rect axes. |
| `test_no_narrow_neck_junctions[CYXY]` | CYXY | Junctions with thin necks (residue around runway crossings has narrow fingers). |
| `test_junction_vertices_outside_pavement[CYXY]` | CYXY | Junction ring has vertices outside `pav_union` (runway-crossing residue subtraction). |
| `test_rect_short_edges_connect[SPLP]` | SPLP | Cross (short) edge of a stub not connecting to neighbour properly. |

### Visual generation issues to verify

User asked for visual eval after tests pass.  Build with the
provided one-liner below and inspect in JOSM:

```bash
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
for icao in ("SPJC", "SPLP", "CYXY"):
    layout = build_airport_pavement(
        icao, "/Users/noah/X-Plane 12", compute_elevations=True)
    layout.to_osm(f"/tmp/{icao}_current.osm")
PY
```

Known visual checkpoints (user-reported across sessions):
* SPJC + SPLP target fixtures — used to match before `a9b0ef0`,
  now off by 1–2 rect matches.
* CYXY D-E junction polygon size, D-west rect, primary E south
  overshoot, rect-rect overlaps — these were the focus of
  `a9b0ef0` and reportedly substantively fixed per the
  morning STATUS handoff.
* CYXY bridge polygons — new shapes from `3a21c89`; verify
  they hug pavement on the inner side without crossing
  pavement.

## Architecture notes

### Bridge architecture (post-`3a21c89`)

The bridge is constructed in two passes:

1. **`finalize.py:262` — initial emit** during the elevation
   phase.  Uses the current pavement_union state.
2. **`pipeline.py` (after `_split_sloped_rects_at_violations`,
   line ~2253)** — drops finalize-time bridges and re-emits with
   final pavement state.  Replaces the old
   `_clip_boundary_bridges_against_pavement` re-clip (which used
   `difference()` and re-introduced intersection vertices).

If a future pass moves pavement vertices AFTER the re-emit, the
bridges become stale.  Currently `_split_sloped_rects_at_violations`
is the last pass that does this; the bridge re-emit runs right
after it.

### Boundary ribbon architecture (post-`586b135`)

The airport boundary ribbon is a chain of 4-corner rectangles
(one per densified boundary segment).  Each rect uses
**per-vertex perpendiculars** (averaged across incoming +
outgoing segment perps at each boundary vertex), so consecutive
rects share their cross-edge nodes at the shared vertex.  My
bridge rewrite uses the same per-vertex perp algorithm for its
outer edge — so bridge outer vertices coincide with ribbon
inner corners.

### Test 3 / canonical-node invariant

`test_junction_neighbour_corners_shared` enforces: every
non-junction polygon vertex within `ORPHAN_NEAR_PERIMETER_M`
(1.0 m) of a junction's perimeter MUST coincide with a junction
vertex (within `ORPHAN_SAME_VERTEX_TOL_M` 0.10 m).  After
`3a21c89`, bridge vertices satisfy this BY CONSTRUCTION (they're
either canonical pavement nodes or 2.5 m off the boundary line
where no pavement is present, beyond the 1 m perimeter check
for any junction).  Remaining Test 3 failures involve other
polygon classes (groundside, airport_boundary, stub sub-rects).

## How to verify / reproduce

```bash
# Full test suite (slow — ~3 min):
venv/bin/python3 -m pytest tests/ --tb=short -q

# Run just the Test 3 family at all 3 airports:
venv/bin/python3 -m pytest \
  "tests/test_junction_invariants.py::test_junction_neighbour_corners_shared" \
  --tb=long -v

# Investigate any specific failure with detailed output:
venv/bin/python3 -m pytest \
  "tests/test_pavement_geometry.py::test_no_vertex_on_sloping_rect_edge[SPJC]" \
  --tb=long -v
```

## Recommended next session

Order of attack:

1. **SPLP `_split_sloped_rects_at_violations` twin-corner**
   (6 Test 3 violations).  Extend `c1e2ad0` to also handle the
   OPPOSITE-side new sub-rect corner: when the source junction
   vertex moves to one new corner, find adjacent junctions
   whose perimeter passes near the TWIN corner and insert/snap.
   Scoped to split-introduced corners only — no broader
   insertion pass (the broad approach was tried in `bf37313`
   and reverted in `0c733ba` due to cascade complexity).
2. **Phenomenon A (untagged clip residue)** at SPJC + SPLP.
   `_drop_overlap_against_fixed_shapes` produces non-rect
   sloping rects.  Either prevent the clip (reject non-rect
   output, absorb instead) or demote the clipped shape's role
   to non-sloping.  See `feedback_root_cause_only.md` —
   prefer source fix to post-process cleanup.
3. **`test_compare_target_*` (SPJC, SPLP)** — investigate the
   1-2 rect matches lost since `a9b0ef0`.  If shifts are real
   improvements, re-anchor floor.  If regressions, locate +
   fix.  Cross-reference with the visual eval for "did the
   geometry get better or worse?"
4. **`test_pavement_grade[SPJC, SPLP]`** — geometry shift
   from `a9b0ef0` moved shared corners off elevation alignment.
   Likely a small refinement to the chart-junction margin
   (introduced in `a9b0ef0`'s change #4) — currently uses
   `CHART_JUNCTION_MARGIN_M = 25` unconditionally.
5. **CYXY-specific (Tier 2)**: `test_taxi_rects_not_alongside_apron`,
   `test_large_junction_axis_aligned_borders`,
   `test_no_narrow_neck_junctions`,
   `test_junction_vertices_outside_pavement`.  All related to
   runway-crossing residue at CYXY.  Diagnose collectively;
   likely shared root cause in junction emit around runway
   crossings.

## Memory notes — read before continuing

* `feedback_general_solutions.md` — every fix must work at all
  baseline airports.  No airport-specific code.  Per user
  2026-05-16: "any changes we make to fix bugs at any
  particular airport always have to take into account whether
  the change will work across all airports."
* `feedback_root_cause_only.md` — fix root causes; ASK before
  band-aid post-process clean-up.
* `feedback_flat_segment_node_count.md` — `altitude=` tag ⇒
  flat shape, any node count; `altitude_high/low` ⇒ 4-corner
  sloped rect.  Tests 1 & 2 fixes in this session respect this.
* `feedback_grade_rules.md` — apron / junction grade rule;
  all-pair Euclidean only between pairs whose straight line
  stays inside the polygon (user 2026-05-14 refinement).
* `feedback_shape_rules.md` — authoritative rect + junction
  construction rules.
* `project_refactor_state.md` — extraction pattern.

## Files touched this session

Source:
* `src/auto_patch/boundary.py` — `_emit_boundary_dem_bridge`
  rewritten (canonical-node sequential walk); per-vertex perp
  in `_emit_airport_boundary_shape`.
* `src/auto_patch/junction_repair.py` —
  `_split_sloped_rects_at_violations` extended to move source
  junction vertex onto new sub-rect corner.
* `src/auto_patch/pipeline.py` — bridge re-emit after
  `_split_sloped_rects_at_violations` (replaces re-clip pass).

Tests:
* `tests/test_pavement_geometry.py` — flat-shape exemption in
  `test_no_vertex_on_sloping_rect_edge`.
* `tests/test_junction_rules.py` — flat-rect corner exemption
  in `test_junction_no_long_edge_proximity`.

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
