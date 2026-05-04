# Auto-Patch Status — handoff 2026-05-05

## TL;DR

Three passes landed today, in commit order:

1. **Pass 1 — boundary-traced widening with multi-step walking**
   (commit 6b3eb4e).  Extended each junction to runway corners up
   to 2 chain steps from any original runway-shared anchor, with
   pavement-boundary waypoints between the new corner and the
   junction's flank vertex.  Worked, but most of it became dead
   weight after Pass 2 landed.

2. **Pass 2 — apt.dat-aware runway segmenter** (commit 15d3a89).
   Root-cause fix: place runway-segment seam corners at every
   point where an apt.dat row-110 pavement polygon's boundary
   touches the runway boundary (within 0.5 m, deduped within 2 m
   centerline distance), instead of only at fixed 100 m
   intervals.  Chain corners now naturally align with apt.dat-
   pavement boundary intersections.

3. **Pass 3 — strip Pass 1's complexity** (commit 9cefc66).  With
   chain corners aligned to apt.dat boundaries, the boundary-
   trace logic + interior-vert pruning + multi-step queue are no
   longer load-bearing.  Removed all of that from
   ``_do_widen``; kept only the original single-step walker plus
   a sliver-corner pre-emption check inside ``_attempt_insert``.

Today's tests (SPJC only):
```
O4_TEST_AIRPORTS=SPJC ./venv/bin/pytest tests/
→ 218 passed, 2 failed (test_pavement_grade SPJC + SPLP — out
  of scope per prior STATUS).  Zero new regressions.
```

---

## Where the new code lives

* `src/auto_patch/pipeline.py:340-411` — collects pav-runway
  intersection points after shoulder absorption.  Walks each
  apt.dat row-110 pavement polygon's exterior; for each vertex
  within 0.5 m of any runway 4-corner rect AND projecting to
  centerline ``t`` in (5 m end-skirt, 1 - 5 m end-skirt), records
  ``(t, lat, lon)``.  Sorts and dedups within 2 m centerline
  distance.  Stashes on ``layout._pav_runway_intersections`` keyed
  by all four designator orderings.

* `src/auto_patch/elevation.py:374-385` — reads the stash, passes
  to ``generate_patch_osm``.

* `src/auto_patch/pavement/runway_segments.py:404-440` — new
  ``pav_intersections`` parameter.  Inside the per-runway loop,
  projects each (lat, lon) onto the phys-end-A → phys-end-B
  centerline, derives ``t``, inserts into the ``fractions`` list
  after threshold injection (with the same 2 m dedup against
  existing fractions).  These t-values are NON-ANCHORED — they're
  segment seams, not elevation constraints.

* `src/auto_patch/junction_rules.py:_attempt_insert` — sliver-
  corner pre-emption check: rejects any trial polygon whose ring
  would have an interior angle below
  ``SLIVER_ANGLE_THRESHOLD_DEG`` (≈ 2°).  Without this guard the
  OSM emitter silently drops the whole junction at emit time
  whenever a chain step lands on a near-180° spike.

## Effects on SPJC

* 73 → 85 runway segments (24 new seams from 12 intersections per
  runway × 2 runways).
* Chain corners now align with apt.dat-pavement boundary
  intersections — chain-walk-based passes (Rule 1 v4 sharing,
  widening) reach the right corners naturally.
* Junction widening rule is back to the simple single-step walker
  it was at 5a50d00.

## Outstanding

* User's hand-edited target ``/private/tmp/SPJC_root_fix_EDITED
  .osm`` had B[3], B[5], B[9] as junction vertices.  These were
  100 m-fixed seam corners from the OLD segmenter; the NEW
  segmenter places different seam corners (at apt.dat-pavement
  intersections).  My build's J-10131 / J-10136 share corners
  with the NEW chain at apt.dat intersections — visually the
  polygons should look right, but they don't match EDITED
  position-for-position.  Verify with JOSM at next session.

* Multi-airport tests (CYXY/HECA/KBNA/SPLP) all show pre-existing
  failures unrelated to today's work (verified by checking out
  5a50d00 — same 23 failures).  Default ``O4_TEST_AIRPORTS=SPJC``
  is what we test against.

---

## Pass 1 — Boundary-traced widening (committed but mostly reverted)

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
