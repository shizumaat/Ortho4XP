# Auto-Patch Status — handoff 2026-05-01

## TL;DR

The auto-patch feature has been **fully refactored, packaged, and
unit-tested**, but **24 pre-existing pavement-geometry tests are
still failing** at our four target airports.  The next agent's job
is to work through Phases A–E in
[`zany-dancing-moonbeam.md`](https://example.com) (the active
plan) to close them.  No structural work remaining.

**Active plan:** `/Users/noah/.claude/plans/zany-dancing-moonbeam.md`

---

## Where the code lives

Everything auto-patch is in `src/auto_patch/`:

```
src/auto_patch/
├── __init__.py            ← public API: generate_auto_patches
├── driver.py              ← tile-level orchestrator (Vector_Map calls this)
├── pipeline.py            ← build_airport_pavement (per-airport)
├── finalize.py            ← Phase-2 elevation + feature emit orchestration
├── junction_emit.py       ← residue → junction polys + geometry-finalize
├── osm_load.py            ← OSM/apt.dat loaders
├── osm_aeroway.py         ← OSM data extraction (Vector_Map calls this)
├── apt_dat_reader.py
├── cifp_reader.py
├── dsf_reader.py
├── layout.py              ← BuiltShape, PavementLayout, role consts (foundation)
├── config.py              ← shared tunables (foundation)
├── elevation.py           ← phase-2 altitude solver
├── elevation_smoothing.py ← _smooth_polygon_grid
├── triangulation.py       ← _triangulate_junctions (per-vertex altitude)
├── junction_repair.py     ← clamp + subdivide + merge-sliver
├── boundary.py            ← airport perimeter ribbon
├── bridges.py             ← taxi/road bridges (gated off via EMIT_BRIDGES_AND_TUNNELS)
├── groundside.py          ← curbside / drop-off pavement
├── terminals.py           ← OSM terminal building pads
└── pavement/              ← airside paved-surface construction (15 files)
    ├── absorption.py      ← THE HOT SPOT — long-edge absorption rule
    ├── centerlines.py     ← OSM aeroway centerline extraction + splitting
    ├── classifier.py      ← apt.dat pavement classifier
    ├── junctions.py       ← junction polys (decomposition + densification)
    ├── rects.py           ← centerline → taxi-rect construction
    ├── runways.py         ← runway rect, crossings, shoulders
    ├── runway_geometry.py ← runway corner / pair primitives
    ├── runway_segments.py ← CIFP-driven runway-profile segments
    ├── strips.py          ← strip decomposition + adjacency
    ├── stubs.py           ← runway-end primary-parallel stubs
    ├── taxiway_decompose.py
    ├── taxiway_rects.py
    ├── taxiway_skeleton.py
    ├── union_helpers.py   ← _merge_near_touching
    └── vertices.py        ← shared-vertex enforcement
```

External entry points: `O4_Vector_Map.py` imports
`generate_auto_patches` and `osm_aeroway` from `auto_patch`.  Tests
in `tests/` import via `from auto_patch.X import Y`.

---

## Test baseline (must not regress)

```
O4_TEST_AIRPORTS=CYXY pytest tests/
→ 5 failed, 207 passed

O4_TEST_AIRPORTS=KPHX,SPJC,SPLP,CYXY pytest \
    tests/test_pavement_geometry.py \
    tests/test_junction_invariants.py \
    tests/test_pavement_grade.py
→ 22 failed, 12 passed
```

The 5 CYXY failures (and 22 multi-airport equivalents) are **the
targets of Phases A–E** in the plan.  Don't touch them as
prerequisites; the plan handles them in dependency order.

The 207 passing tests include:
* 70 new unit tests added 2026-05-01 covering the 4 highest-value
  modules without prior coverage (union_helpers, runway_geometry,
  absorption, layout).
* All 6 pre-existing module-level test files
  (apt_dat_reader, cifp_reader integration, classifier, strips,
  taxiway_decompose, taxiway_rects).
* Surface-mesh tests + auto-patch integration tests at CYXY
  (geometry, junction-invariants except the 5 known failures).

---

## What's pending (in plan order)

| Phase | Targets | Notes |
|---|---|---|
| A | `test_junction_vertex_count_bounded` (4 airports) | Add RDP simplification before junction decomposition.  Pavement_Junctions per-arc densification cap of 4 already landed; this addresses the input boundary trace. |
| B | `test_junction_boundary_near_centerline` (4 airports) | Reclassify junctions whose boundary is > 20 m from any taxi/runway centerline as `role=apron`. |
| C | `test_taxi_rects_not_alongside_apron` (+ `test_rect_short_edges_connect[KPHX]` + `test_coverage_within_source_envelope[KPHX]`) | Diagnose-then-fix.  The absorption rule is correct (DO NOT change semantics — it's the recurring-regression hot spot).  The job is to identify why the rule's input set diverges from the test's view at the failing airports. |
| D | `test_no_self_overlap` (4 airports) | Lower the noise threshold in `_drop_overlap_against_fixed_shapes`; add `ROLE_BOUNDARY` to the priority list. |
| E | `test_junction_neighbour_corners_shared` (4 airports) | New helper to insert neighbour vertices into junction perimeters before final shared-vertex enforcement. |

`test_pavement_grade[SPJC|SPLP]` (2 failures) is **out of scope**
for this plan per user 2026-05-01 — slope is enforced at taxiway
emit time.  Track only if the geometry phases above don't fix it
incidentally.

---

## Useful commands

```bash
# Run full suite at CYXY
O4_TEST_AIRPORTS=CYXY venv/bin/python -m pytest tests/

# Run all 4 target airports against the geometry/junction/grade tests
O4_TEST_AIRPORTS=KPHX,SPJC,SPLP,CYXY venv/bin/python -m pytest \
    tests/test_pavement_geometry.py \
    tests/test_junction_invariants.py \
    tests/test_pavement_grade.py

# Build one airport standalone to inspect the .osm output
PYTHONPATH=src venv/bin/python tools/build_target_osm.py SPJC \
    --xplane "/Users/noah/X-Plane 12" \
    --out /tmp/SPJC_auto.osm

# Compare with the reference target OSM
venv/bin/python tools/compare_target.py \
    /tmp/SPJC_auto.osm tests/fixtures/SPJC_target.osm

# Skip integration tests for shipping
O4_SHIP_MODE=1 venv/bin/python -m pytest tests/
```

---

## Gotchas + invariants for the next agent

1. **Don't touch the absorption rule's semantics.**
   `src/auto_patch/pavement/absorption.py` is the recurring-regression
   hot spot — its docstring (lines 31-43) explicitly forbids:
   - Switching EITHER → BOTH side detection.
   - Adding per-ref preservation exceptions.
   - Changing the 10 % adjacency threshold or 5 m probe distance.
   Phase C is about fixing the rule's INPUT set, not its semantics.

2. **Bridge / tunnel emit is gated off** via
   `EMIT_BRIDGES_AND_TUNNELS = False` in `src/auto_patch/config.py`.
   Re-engagement is a separate post-stabilisation feature pass —
   each emitter must clip its footprint against airside / groundside
   pavement first or `test_no_self_overlap` will fail.

3. **Same-baseline gate**: every commit in Phases A–E must keep the
   currently-passing 207 tests green and target only the failing
   ones.  Run the multi-airport suite before commit.

4. **Stray `+60-140/`** at the repo root is leftover .hgt elevation
   tiles from a misplaced run (24 files, ~24 MB).  Not in `.gitignore`;
   already escaped one accidental `git add -A`.  Don't add it.

5. **Cherry-picked from upstream** in this session: PR #88
   (`shred86/Ortho4XP`) — overpass refactor, GUI memory leak fix,
   version bump to 1.40.13.  When merging future upstream PRs, the
   `auto_patch/` package should not conflict (it's all fork-only files).

---

## Memory pointers

The following memory files are loaded into every session via
`MEMORY.md`:

- **`feedback_extraction_pattern.md`** — orphan-constant audit
  recipe used during the refactor.  Phases A–E are no longer
  extractions but the audit pattern is still useful when adding
  new module-level constants.
- **`feedback_shape_rules.md`** — authoritative rect + junction
  construction rules.  The absorption rule lives here too.
- **`feedback_general_solutions.md`** — every fix must work at any
  airport, no hardcoded ICAOs/refs.
- **`project_target_osm.md`** — role tags + shared-vertex invariants
  for the SPJC/SPLP reference targets.
- **`project_refactor_state.md`** — historical snapshot of the
  monolith refactor (slices 0–3i).  Mostly historical now.

---

## Latent issues found and fixed during this session

1. `auto_patch/pavement/runway_geometry.py` was missing
   `import os` — `parse_aptdat_runway_widths` would NameError but
   the wrapping try/except swallowed the exception silently
   (returns empty dict).  Caught by the new unit test, fixed in
   commit `25fa16e`.

That's the only one.  No other latent issues surfaced.
