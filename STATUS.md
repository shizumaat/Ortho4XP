# Auto-Patch Status — handoff 2026-05-05 (evening)

## TL;DR

Five passes landed today on the runway-junction widening problem.
SPJC's runway junctions are now all clean per user verification.

| commit | summary |
|---|---|
| 6b3eb4e | Pass 1: boundary-traced widening + multi-step + sliver pre-emption (mostly reverted later) |
| 15d3a89 | **Pass 2: apt.dat-aware runway segmenter** — seam corners at every apt.dat-pavement intersection |
| 9cefc66 | Pass 3: strip Pass 1's complexity (boundary-trace, interior pruning) post-Pass 2 |
| 7d856f7 | Pass 4: smart side-picker + max_total_shared 5→7 + segmenter dedup (apt.dat replaces uniform seams within 12 m) |
| f41a09b | Pass 5: gate multi-step on whether polygon body has a non-anchor far from anchors |
| 8f6d517 | **Pass 6 (final): skip widening when polygon naturally aligned + targeted interior-vert prune** |

Today's tests (SPJC only):
```
O4_TEST_AIRPORTS=SPJC ./venv/bin/pytest tests/
→ 218 passed, 2 failed (test_pavement_grade SPJC + SPLP — out
  of scope per prior STATUS).  Zero new regressions.
```

User confirmed all SPJC runway junctions look correct at HEAD.
Latest build: `/tmp/SPJC_v8.osm`.

---

## What landed (in commit order, current code state)

### Pass 2 — Apt.dat-aware runway segmenter

Place runway-segment seam corners at every point where an apt.dat
row-110 pavement polygon's boundary touches the runway boundary.
Variable-length segments — elevation/grade interpolation works
regardless of segment length.

**Where the code lives:**

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
  centerline, derives ``t``.  When an apt.dat ``t`` lands within
  12 m of a NON-anchored uniform seam, REPLACE the uniform seam
  with the apt.dat position (= dedup; eliminated all sub-12 m
  sliver segments).  Anchor fractions (CIFP thresholds, physical
  ends) are protected — apt.dat points within 2 m of an anchor
  are dropped instead.

### Junction widening rule (current state)

`src/auto_patch/junction_rules.py::_do_widen` — three coupled
gates determine widening behavior per junction:

1. **Skip widening when polygon naturally aligned** — if the
   polygon's runway-shared anchors form a chain-adjacent
   contiguous run AND no non-anchor sits more than 70 m from any
   anchor, skip widening entirely.  These polygons sit tight
   against the runway with their runway-side edge already
   aligned with runway corners; round-0 widening would extend
   them past the body's natural extent.

2. **Multi-step gate** — round 0 (single-step from each original
   anchor) always runs.  Round 1 (multi-step from newly-inserted
   corners) runs only when the polygon has a non-anchor vertex
   more than 50 m from any anchor.  Otherwise round 0 only.
   Cap: 2 rounds total.

3. **max_total_shared = 7** — caps how many runway-shared
   corners a junction can have post-widen (originals + inserts).
   Bumped from 5 because the apt.dat-aware chain has more
   corners adjacent to each junction's runway-side edge, and
   J-10131 / way -10138 needs 6 to reach B[9].

After widening, **targeted interior-vert prune**: drop polygon
vertices that sit > 5 m from any pav_union/runway boundary AND
are adjacent (in walk order) to a newly-inserted chain corner.
The new chain-corner arm supersedes the old interior detour.
Validates the pruned polygon (valid + simple + area ≥ 0.5 ×
original) before committing.

**Smart side-picker** in the inner walk loop: when picking which
side (BEFORE / AFTER the chain corner in polygon walk order) to
insert a chain neighbor, compute cos_turn for both sides and pick
whichever produces the LESS-acute corner (= higher cos).
Previously picked by bearing-closeness which arbitrarily chose
spike-prone sides at near-180° flank angles.

**Sliver-corner pre-emption** inside ``_attempt_insert``:
rejects any trial polygon whose ring would have an interior
angle below ``SLIVER_ANGLE_THRESHOLD_DEG`` (≈ 2°).  Without this
guard the OSM emitter silently drops the whole junction at emit
time.

---

## Effect on SPJC

* 73 → 80 runway segments (apt.dat intersections add ~12
  per-runway, dedup with uniform seams nets the count up by 7).
* Chain corners aligned with apt.dat-pavement boundary
  intersections — chain-walk-based passes (Rule 1 v4 sharing,
  widening) reach the right corners naturally.
* User-verified runway junctions:
  * way -10138 (J-10131): 12 verts, 6 runway-shared incl. B[9],
    A9_INT correctly pruned post-widen.
  * way -10137: 4 verts, 2 runway-shared (back to pre-widen
    state — no over-extension).
  * way -10136: 8 verts, 2 runway-shared (back to pre-widen
    state — no over-extension).

---

## Multi-airport status

Default ``O4_TEST_AIRPORTS=SPJC`` is the only set we run.
Multi-airport tests (CYXY/HECA/KBNA/SPLP) all show 23 pre-
existing failures unrelated to today's work — verified by
checking out 5a50d00 (pre-this-session commit) and seeing the
same failure set.  Don't worry about other airports unless the
user asks.

---

## Outstanding items (next session focus per user 2026-05-05)

### §1 — Holes in the apron

User flagged for next session: holes appearing in the apron.
Specifics TBD — user will point at examples in JOSM at next
session.  Likely root cause candidates:
* Over-aggressive subtraction (taxi rects / runway / groundside
  taking too much from pav_union).
* Sliver gaps between adjacent shapes that don't share vertices.
* DSF overlay polygons missing from pav_union.

### §2 — Junctions / aprons bordering the SLOPING edge of taxi rects

This is §1 from the prior STATUS, restated by user 2026-05-05:
junctions and aprons that share boundary with the SLOPING edge
of a taxiway rect cause grass areas to be lost (e.g. east of
terminal 2 along taxiways M and L at SPJC).

Existing rule passes:
* ``_snap_to_sloping_edge_corners`` snaps junction vertices
  WITHIN ``SLOPING_EDGE_SNAP_M`` (= 20 m) of a sloping rect's
  sloping edge to the nearest cross-edge corner.
* ``_push_junction_vertices_outside_pavement`` (Rule 5) pushes
  junction vertices outside apt.dat pavement.

These run BEFORE / AFTER widening.  But junctions can still
end up sharing boundary with sloping edges (not just vertices).
The fix likely needs a polygon-edge-vs-sloping-edge check, not
just vertex-distance.

Investigation suggestions:
* Walk each junction's polygon edges; for each edge that runs
  parallel to and within ``SLOPING_EDGE_SNAP_M`` of a sloping
  rect's sloping edge, log it as a violation.
* Visual inspection in JOSM: open ``/tmp/SPJC_v8.osm`` at the
  east-of-terminal-2 region (around lat -12.020, lon -77.114?).

### §3 — Things to NOT do (carried over)

1. Don't add "search-and-repair" passes for visual artefacts —
   user prefers root-cause fixes.
2. Don't change the absorption rule semantics.
3. Don't reclassify junctions to apron based on centerline-
   distance.
4. Don't whole-rect-flip rects to apron.
5. Don't re-disable bridges/tunnels.
6. Don't re-instate the altitude gate in
   ``_snap_to_sloping_edge_corners`` or in
   ``_densify_long_boundary_edges``'s ``sloping_rect_edges``
   collection.
7. Don't reduce the 20 m runway-snap radius without user
   discussion.
8. **Don't multi-airport-test as a regression gate** — pre-
   existing failures at CYXY/HECA/KBNA/SPLP are not in scope.
   Use ``O4_TEST_AIRPORTS=SPJC`` only.

---

## Pipeline order reference (unchanged)

```
build_airport_pavement(icao, xplane_root, compute_elevations=True)
├── load_airport (apt_dat_reader)
├── runway_polys built (4-corner rect each)
├── pav_polys built (apt.dat row-110 in meter space)
├── runway-shoulder absorption
├── ★ pav_runway_intersections collected (Pass 2 hook)
├── pav_union = unary_union(pav_polys) - runway - groundside
├── emit rect shapes
├── junction_emit.emit_junctions_and_finalize
│   └── apply_junction_rules
│       ├── _align_rect_slope_to_axis
│       ├── _snap_to_sloping_edge_corners
│       ├── _enforce_runway_1to1_sharing
│       ├── widen_junctions_to_runway_corners  ← Pass 6 lives here
│       └── _push_junction_vertices_outside_pavement
└── if compute_elevations:
    ├── finalize.run_phase2 → _compute_elevations
    │   ├── ★ generate_patch_osm with pav_intersections (Pass 2)
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

## Memory pointers (loaded into every agent session)

* `feedback_shape_rules.md` — authoritative shape rules.
* `feedback_extraction_pattern.md` — refactor extraction recipe.
* `feedback_general_solutions.md` — no airport-specific fixes.
* `project_target_osm.md` — role schema, vertex tolerances.
* `feedback_grade_rules.md` — TAXI 1.5 % / APRON 1.0 % per-axis.
* `project_refactor_state.md` — historical refactor state.

### New memory entry recommended (Pass 2 — landed)

* **Apt.dat-aware runway segmenter** (user 2026-05-05): runway
  segment seam corners are placed at every apt.dat row-110
  pavement-boundary intersection with the runway (within 0.5 m).
  Apt.dat intersections REPLACE non-anchor uniform 100 m seams
  within 12 m (eliminates sliver segments).  CIFP thresholds and
  physical ends are anchored — apt.dat intersections within 2 m
  of an anchor are dropped instead.  Effect: chain corners (=
  segment seams) naturally align with apt.dat-pavement
  boundaries, so chain-walk-based passes (Rule 1 sharing, widen)
  reach the right corners without boundary-trace waypoints.

### New memory entry recommended (Pass 6 — junction widening gates)

* **Junction widening conditional gates** (user 2026-05-05):
  widen_junctions_to_runway_corners now SKIPS widening entirely
  for junctions whose anchors are chain-adjacent AND no non-
  anchor reaches > 70 m off any anchor (the polygon is naturally
  aligned post-Rule-1; round-0 walking would extend it past its
  body extent).  Multi-step (round 1) runs only when a
  non-anchor reaches > 50 m off any anchor (the wrap case at
  J-10131).  Targeted interior-vert prune drops adjacent
  interior verts after a chain-step insert (e.g., A9_INT in
  -10138).

---

## Files of interest (refer at next-session start)

* User's hand-edited target: `/private/tmp/SPJC_root_fix_EDITED.osm`
  (no longer matches HEAD position-for-position because the new
  segmenter places seam corners at apt.dat intersections, not at
  100 m-fixed positions; SHAPES are now correct per user).

* Latest auto-patch build: `/tmp/SPJC_v8.osm`.

* Geometry dumps for visual inspection:
  - `/tmp/SPJC_geometry_dump.osm` — pav_union + runway_union
    outlines + probe points.
  - `/tmp/SPJC_compare_J131_J136.osm` — overlay of mine vs
    EDITED for the two hand-edited junctions.

* Production patch path (rebuild via
  `tools/build_target_osm.py SPJC --out <path>` if you want the
  X-Plane-installed file refreshed):
  `/Users/noah/Ortho4XP-shred86/Patches/-20-080/-13-078/SPJC_auto.patch.osm`.

---

## Stray leftover

`+60-140/` directory at repo root — 24 .hgt elevation tiles,
untracked.  Don't `git add -A`.
