# Auto-Patch Status — handoff 2026-05-05

## TL;DR

This session built two passes on the runway-junction widening
problem:

1. (now committed) **boundary-traced widening with multi-step
   walking** — extends each junction to runway corners up to 2
   chain steps from any original runway-shared anchor, with
   pavement-boundary waypoints between the new corner and the
   junction's flank vertex.

2. (planned next) **apt.dat-aware runway segmenter** — instead of
   chasing the gap between fixed-interval runway seams and
   apt.dat-pavement boundaries with boundary-trace waypoints,
   place the seams DIRECTLY at every pavement-boundary
   intersection.  This is a root-cause fix that obsoletes most of
   what (1) does.

Today's tests: 218 passed, 2 failed (`test_pavement_grade[SPJC|
SPLP]`, both out-of-scope).  No new regressions.

---

## Pass 1 — Boundary-traced widening (committed)

### Why

User provided a hand-edited target OSM
(`/private/tmp/SPJC_root_fix_EDITED.osm`) showing two runway
junctions with extra runway corners + apt.dat boundary
waypoints.  The auto build was missing those — the existing
widening rule walked exactly one chain step from each ORIGINAL
runway-shared corner and rejected anything that produced a near-
180° spike at the chain corner.

### Files changed (committed in this session)

* `src/auto_patch/junction_rules.py`:
  * Added `_trace_pav_boundary_waypoints` — picks a single
    apt.dat-pavement-boundary ring vertex that lies between a
    runway chain corner and the polygon's flank vertex, with
    perpendicular offset > 3 m (to avoid 180° spikes) and within
    25 m perpendicular of the chord.
  * Added `_on_sloping_edge_interior` — filters waypoints that
    would land within `SLOPING_EDGE_SNAP_M` (20 m) of a sloping-
    rect / runway-segment edge interior.  Keeps Rule 2 and
    `test_no_vertex_on_sloping_rect_edge` clean.
  * Refactored `_do_widen` to a queue-based walker (multi-step,
    capped at 2 rounds — 1 step from each original anchor + 1
    step from each newly-inserted corner).
  * Added `_validate_trial` with a sliver-corner pre-emption
    check (rejects any trial whose ring would be dropped at OSM
    emit time as a sliver).
  * Added interior-vert pruning at end of widening — drops
    polygon verts > 5 m from any pav_union/runway boundary that
    aren't anchored (chain corners or new waypoints).
  * Replaced the U-turn cos < -0.99 short-circuit with: try
    boundary-trace fallback when cos < -0.95; if no usable
    trace, plain-insert is attempted only when cos >= -0.95.
* `src/auto_patch/junction_rules.py` — also added
  `SLIVER_ANGLE_THRESHOLD_DEG` import from config.

### What's NOT changed today

Nothing else.  No config, no pipeline, no other rule passes.

### Result on SPJC vs EDITED target

* `J-10136`: 9 verts vs EDITED's 7.  Has both runway-end-zone
  corners (B[3], B[5]), one apt.dat-boundary waypoint per arm.
  Polygon shape correct; 2 leftover boundary verts the user
  manually dropped.
* `J-10131`: 13 verts vs EDITED's 12.  Multi-step reaches B[9]
  via chain[A[7] → A[8] → B[9]] plus 1 boundary waypoint.  Has
  one extra chain corner (NEW[6], one more chain step north past
  A[5]) and 1 fewer boundary waypoint than EDITED's two.  User
  said "going 2 steps is fine" so this matches their stated
  preference.

### Test baseline (committed state)

```
O4_TEST_AIRPORTS=SPJC ./venv/bin/pytest tests/
→ 218 passed, 2 failed
   (test_pavement_grade[SPJC] + [SPLP] — out of scope)
```

`tests/test_compare_target.py` baselines unchanged from
2026-05-04 session (`SPJC junction=24, total=63`).

---

## Pass 2 — Apt.dat-aware runway segmenter (PLANNED, NOT
STARTED)

### Why this exists

The current runway segmenter (`_runway_rect_m` plus the
splitting that produces 73 runway shapes for SPJC) places seam
corners at fixed intervals (~100 m).  Apt.dat-pavement
boundaries meet the runway at arbitrary points along its length
that don't generally coincide with these seams.  This forces the
junction-widening pass to "chase" the apt.dat boundary by
inserting boundary-trace waypoints — exactly what Pass 1 does.

### What changes

Walk every apt.dat pavement polygon's boundary; collect points
where it touches the runway boundary (within ε ≈ 0.5 m).  These
points define the runway's seam corners.  Variable-length
segments are fine — elevation/grade interpolate linearly along
the centerline parameter regardless of segment length.

User-stated rule (2026-05-05): if two incoming intersections are
< 2 m apart, collapse to a single seam corner.  The junction can
span a 2 m gap without needing a node there.

### What it obsoletes

* `_trace_pav_boundary_waypoints` and the boundary-trace insert
  logic in `_do_widen` — no longer needed when chain corners
  already align with apt.dat boundary points.
* `_on_sloping_edge_interior` — same.
* The interior-vert pruning at end of widening — same.
* The sliver pre-emption check — likely still useful as a guard
  but no longer the load-bearing piece.
* The multi-step walker queue and `WIDEN_MAX_ROUNDS` — single-
  step walking from each runway-shared anchor will reach all the
  right corners directly.
* The Rule-2 (`SLOPING_EDGE_SNAP_M`) waypoint filter — boundary
  waypoints disappear; chain corners are by definition rect
  corners.

What stays in `widen_junctions_to_runway_corners`: the original
"insert chain[ci-1] / chain[ci+1] adjacent to each runway-shared
junction vertex" logic, with the existing per-insert
validation.

### Implementation sketch

Today's segmenter (in `pavement/runways.py`):
```
for r in apt.runways:
    rect = _runway_rect_m(r, to_m)   # 4-corner quad spanning
                                      # end-to-end + blast pads
    runway_polys.append(rect)
    layout.shapes.append(BuiltShape(polygon=rect, role=ROLE_RUNWAY,
                                    ref=ref))
```
The 4-corner rect later gets densified-and-split into 73 short
segments (need to find where this happens — likely
`finalize.run_phase2` / `_compute_elevations` / segment-emit).

New segmenter:
1. After building the 4-corner runway rect AND the apt.dat
   pav_polys, collect each pav_poly's vertices that lie on the
   runway boundary (within 0.5 m).  Project each onto the
   runway's centerline → centerline parameter `t ∈ [0, 1]`.
2. Sort `t` values, dedup within 2 m centerline distance.
3. Add `t = 0` and `t = 1` (runway ends).
4. Build segments between consecutive `t` values.  Each segment
   is a quad whose short edges sit at the parameter values; the
   long edges follow the runway's left/right side.
5. Set `altitude_high` and `altitude_low` on each segment from
   linear interpolation of the runway's high/low along the
   centerline parameter.

### Tests to revise

`tests/test_compare_target.py`:
  * Per-airport runway segment count baselines (if any) — likely
    will change.

`tests/test_pavement_geometry.py`:
  * `test_no_vertex_on_sloping_rect_edge` — expected to PASS
    more cleanly after refactor (the load-bearing fix).

`tests/test_junction_rules.py`:
  * `test_junction_no_long_edge_proximity` — same.
  * `test_junction_runway_node_sharing` — same.
  * `RULE2_REGRESSION_BASELINE` may shrink (good).

`tests/test_junction_invariants.py`:
  * Vertex-count and boundary-distance baselines may need
    updating.

`tests/test_pavement_grade.py`:
  * Grade computation should be unchanged (linear interpolation
    along centerline).  But if it samples per-segment, variable
    segment lengths may shift the sample distribution.

### Risks

* Adjacent runway segments get out of order if `t` dedup is too
  aggressive.
* Pav-poly vertices NEAR the runway but not ON it (within 0.5 m
  but actually 2-3 m off) get treated as runway-boundary
  intersections — causes spurious seams.  Tighten to actual
  intersection (line-line crossing or shared vertex), not
  proximity.
* Per-airport runway segment counts shift, breaking baselines.

---

## Pipeline order reference (unchanged)

```
build_airport_pavement(icao, xplane_root, compute_elevations=True)
├── load_airport (apt_dat_reader)
├── pav_union = unary_union(pav_polys) - runway - groundside
├── emit rect shapes
├── junction_emit.emit_junctions_and_finalize
│   └── apply_junction_rules
│       ├── _align_rect_slope_to_axis
│       ├── _snap_to_sloping_edge_corners
│       ├── _enforce_runway_1to1_sharing
│       ├── widen_junctions_to_runway_corners  ← Pass 1 lives here
│       └── _push_junction_vertices_outside_pavement
└── if compute_elevations:
    ├── finalize.run_phase2 → _compute_elevations
    │   ├── _push_junction_vertices_off_taxi_rect_edges
    │   └── _triangulate_junctions
    ├── Post-elevation rule passes (rerun in pipeline.py):
    │   _align_rect_slope_to_axis
    │   _snap_to_sloping_edge_corners
    │   _enforce_runway_1to1_sharing
    │   widen_junctions_to_runway_corners
    │   _push_junction_vertices_outside_pavement
    ├── per_surface_solve (Jacobi)
    ├── stitch_pavement_to_terminals
    └── _report_within_shape_violations (WARN audit)
```

---

## Outstanding items (carried over from 2026-05-04)

### §1 — Apron / junction must NOT join sloping edge of taxiways

User reported (2026-05-04 morning): "we need to not allow
junctions or aprons to join the sloping edge of taxiways. This
is resulting in losing some grass areas on east side of terminal
2 along taxiways M and L."

Investigation done: scanned for J/A vertices within 20 m
perpendicular of M/L sloping edges east of terminal 2; found
none.  Either issue is hidden by the post-snap deformation, or
manifests via a different geometric pattern.  Needs fresh eyes
+ user pointing at a specific shape ID.

Pass 2 may resolve this incidentally if the violations come from
runway-adjacent regions.

### §2 — Apron↔sloping-rect smooth join near terminal 1

User mentioned (2026-05-04): "F segment just south of terminal
1 apron, 30.7→33.5, two places where it leaves a cliff."
Specific F segment IDs unconfirmed.  User confirmation needed.

### §3 — Things to NOT do

(carried over)

1. Don't add a "search-and-repair" pass for the wrong-side /
   bump pattern — user prefers root-cause fix.
2. Don't change the absorption rule semantics.
3. Don't reclassify junctions to apron based on centerline-
   distance.
4. Don't whole-rect-flip rects to apron.
5. Don't re-disable bridges/tunnels.
6. Don't re-instate the altitude gate in
   `_snap_to_sloping_edge_corners` or in
   `_densify_long_boundary_edges`'s `sloping_rect_edges` collection.
7. Don't reduce the 20 m runway-snap radius without user
   discussion.

---

## Memory pointers (loaded into every agent session)

* `feedback_shape_rules.md` — authoritative shape rules.
* `feedback_extraction_pattern.md` — refactor extraction recipe.
* `feedback_general_solutions.md` — no airport-specific fixes.
* `project_target_osm.md` — role schema, vertex tolerances.
* `feedback_grade_rules.md` — TAXI 1.5 % / APRON 1.0 % per-axis.
* `project_refactor_state.md` — historical refactor state.

### Recommended new memory entry (Pass 2)

* **Apt.dat-aware runway segmenter** (user 2026-05-05): runway
  segments must have seams at every pavement-boundary
  intersection (within 0.5 m), with seams < 2 m apart collapsed
  to one.  This eliminates the need for boundary-trace
  waypoints when widening junctions to runway corners — chain
  corners (= segment seams) already align with apt.dat
  boundaries.

---

## Stray leftover

`+60-140/` directory at repo root — 24 .hgt elevation tiles,
untracked.  Don't `git add -A`.

`/tmp/SPJC_geometry_dump.osm`, `/tmp/SPJC_compare_J131_J136.osm`,
`/tmp/SPJC_my_build.osm`, `/tmp/SPJC_after_widen_fix.osm` —
visualization / debug files used during today's work.
