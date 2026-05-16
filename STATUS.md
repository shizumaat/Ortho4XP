# Auto-Patch Status — 2026-05-16 (end of session): baseline-test expansion + 23 cross-airport invariants exposed

## TL;DR

This session landed CYXY geometry fixes AND expanded the invariant /
grade test suite to run on `SPJC + SPLP + CYXY` unconditionally
(`baseline_airports()` in `tests/conftest.py`).  The expansion
surfaced **23 pre-existing invariant violations** across all three
airports that env-gating had been hiding.  These are the work
backlog.  Per user direction, every fix from here on must work
across all three baseline airports — no airport-specific patches.

**Committed this session (in order):**

1. `a9b0ef0` — *CYXY: pre-split apt.dat polylines at junction nodes
   + chart-junction margin.*  Six targeted geometry fixes: apt.dat
   pre-split at junction nodes, `at_endpoint` check in
   `split_merged_centerline`, fallback proportional margin in
   `_split_centerlines_at_points`, `CHART_JUNCTION_MARGIN_M = 25`
   replacing bend-shared override at multi-ref endpoints,
   perpendicular-skip in `_emit_primary_parallel_runway_stubs`,
   `_rect_long_edges_at_pavement_boundary` rect-validity check.
2. `c646a1a` — *Run invariant + grade tests on baseline_airports()
   unconditionally.*  Adds `conftest.baseline_airports()`, lifts
   env-gating on every invariant test, expands `test_pavement_grade`
   to baseline airports with CYXY caps.

**Tests: 241 / 264 pass at HEAD.  23 failures = the backlog.**

## CYXY geometry result (after `a9b0ef0`)

* **D-E junction polygon now substantial** — idx 45 at (259, -938),
  10 332 m², covering the E + D + runway intersection area (was
  three separate junctions totalling ~6 000 m² + the bad U-shape).
* **No mis-sized D-west rect.**  The 64 m × 75 m rect that was
  covering the junction area is gone.  Pavement there stays as
  junction residue.
* **No rect-rect overlaps** in either compute_elevations=True or
  False builds.
* **-10003 (primary E south) top corners** moved from 29-32 m →
  80-101 m from node 141.  Junction polygon has room to form.
* **-10007 (stub E) SE corners** at 42-44 m from node 141 (was
  36-38 m at HEAD).  Modest improvement but still triggers
  `test_taxi_rects_not_alongside_apron[CYXY]` — see Tier 2 below.

The user's primary CYXY complaints (E-D junction polygon size,
mis-sized D-west, -10003 overshoot, rect-rect overlaps) are
substantively addressed.  Remaining issues are in the test backlog.

## The work backlog — 23 failing tests across SPJC / SPLP / CYXY

User direction (2026-05-16):

> "Any changes we make to fix bugs at any particular airport always
> have to take into account whether the change will work across all
> airports."

> "I am concerned pursuing fixing CYXY without maintaining SPJC and
> SPLP along the way will result in too many regressions and possible
> architecture changes that are not sustainable."

So the working principle for every fix below:

1. Diagnose the failure across ALL three baseline airports first.
2. Identify the shared root cause when there is one.
3. Propose ONE fix that works at SPJC, SPLP, and CYXY.
4. Verify against the full test suite before committing.
5. Update per-airport regression baselines (e.g.
   `SELF_OVERLAP_BASELINE_M2`, `MID_EDGE_CAP`) where appropriate.

### Tier 1 — multi-airport invariant violations (high leverage)

| Test | SPJC | SPLP | CYXY | Likely shared root cause |
|---|---|---|---|---|
| `test_no_vertex_on_sloping_rect_edge` | ✗ | ✗ | ✗ | junction / runway-crossing emits a vertex on the interior of a rect's sloping edge |
| `test_junction_no_long_edge_proximity` | ✗ | ✗ | ✗ | junction polygon runs alongside a rect's long edge |
| `test_junction_neighbour_corners_shared` | ✗ | ✗ | ✗ | adjacent junctions disagree on shared-corner positions |

Three tests, 9 failures total.  All three test the same family of
invariants ("junction geometry must not violate rect-axis or shared-
corner constraints").  A single root-cause fix in the right
emission-stage layer probably resolves all 9.

**Recommended starting point.**  Spawn an Explore agent to read each
test's body to understand exact assertion logic, then dump the
failing geometry at each airport, then look for the shared pattern
in `pavement/junction_emit.py`, `pavement/junction_repair.py`, and
the post-junction-emit cleanup in `pipeline.py`.

### Tier 2 — CYXY-specific issues exposed/created this session

* **`test_taxi_rects_not_alongside_apron[CYXY]`** — stub E (-10007)
  has 95 m of its 94 m axis adjacent to apron / junction pavement
  (100 %).  Same test fails at SPLP.  Root cause: the diagonal stub
  E between primary E north (-10002) and primary E south (-10003)
  has its long edges flanked by junction residue from the
  surrounding E-D-runway intersection area.  User's original
  complaint about "junction on the sloping edge of a rect."
* **`test_large_junction_axis_aligned_borders[CYXY]`** — junction
  polygons whose long edges aren't aligned to adjacent rect axes.
* **`test_no_narrow_neck_junctions[CYXY]`** — junctions with thin
  necks (the residue around the runway crossings has narrow
  fingers).
* **`test_junction_vertices_outside_pavement[CYXY]`** — junction
  ring has vertices outside `pav_union` (probably from the
  runway-crossing residue subtraction).

These are all flavours of "junction emission produces messy
residue around the E-D-runway intersection."  Tier 1 fixes may
partially resolve them.

### Tier 3 — grade-test caps + structural-fidelity targets

* **`test_pavement_grade[SPJC]`** — 2 cross-shape proximity
  violations (worst 0.60 m).  Cross-shape proximity is a HARD fail
  (cap = 0); shared corners must agree on elevation.  Root cause:
  geometry shift from `a9b0ef0` moved some shared corners off
  alignment.  Need to identify which corners and either fix the
  geometry or re-record the cross-shape baseline.
* **`test_pavement_grade[SPLP]`** — 26 mid-edge steps > 0.5 m
  (cap 20), worst 1.90 m at junction -10034 / stub A boundary.
  Geometry shift effect on the SPLP threshold-pad area.
* **`test_pavement_grade[CYXY]`** — currently passing with
  CYXY-specific caps recorded at session-end state (`WITHIN_SHAPE_CAP
  = 50`, `MID_EDGE_CAP = 5`).  Tighten as geometry improves.
* **`test_compare_target_spjc`** — `stub: matched = 18 < floor 19`,
  `primary_parallel: matched = 25 < floor 26`.  Two rects shifted
  enough to fail target-fixture match.  Investigate which two; if
  the shift represents a real improvement, re-anchor the target
  floors; otherwise locate and fix the regression.
* **`test_compare_target_splp`** — `primary_parallel: matched <
  floor` for the SPLP tile -13/-78 fixture.  Same investigation.

## Test infrastructure (post `c646a1a`)

* `tests/conftest.py` — `baseline_airports()` returns `("SPJC",
  "SPLP", "CYXY")`.  Phase 2 will add KBNA + HECA once Phase 1
  baselines stabilise.
* Each invariant test file (`test_junction_invariants.py`,
  `test_junction_rules.py`, `test_pavement_geometry.py`) has a
  local `_test_airports()` returning `baseline_airports() +
  airports_under_test()` deduplicated.
* `test_pavement_grade.py` parametrises over `baseline_airports()`;
  per-airport caps in `WITHIN_SHAPE_CAP` and `MID_EDGE_CAP`.
* `test_compare_target.py` stays SPJC + SPLP only — those are the
  hand-drawn-target structural-fidelity regression fixtures.  Adding
  CYXY here would require a CYXY target fixture, which we don't
  have.

`SELF_OVERLAP_BASELINE_M2` in `test_pavement_geometry.py` already
tracks per-airport overlap ceilings (SPJC 1700, SPLP 60, CYXY 200).
These should be a model for how to baseline other tests as we work
through the backlog.

## How to verify / reproduce

```bash
# Build CYXY and dump shapes near the E-D junction area:
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement(
    "CYXY", "/Users/noah/X-Plane 12", compute_elevations=True)
layout.to_osm("/tmp/CYXY_current.osm")
PY

# Run the full test suite (slow — ~3 min):
venv/bin/python3 -m pytest tests/ --tb=short -q

# Run just the failing invariant tests at a single airport:
venv/bin/python3 -m pytest \
  "tests/test_junction_invariants.py::test_taxi_rects_not_alongside_apron[CYXY]" \
  "tests/test_junction_rules.py::test_no_narrow_neck_junctions[CYXY]" \
  -v --tb=long

# Investigate overlap details at any airport (passes baseline 200 m²
# cap at CYXY but shows boundary/boundary ribbon overlaps):
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
from shapely.strtree import STRtree
layout = build_airport_pavement(
    "CYXY", "/Users/noah/X-Plane 12", compute_elevations=True)
polys = [(i, s.role, s.ref, s.polygon) for i, s in enumerate(layout.shapes)
         if s.polygon is not None and not s.polygon.is_empty]
tree = STRtree([p for _, _, _, p in polys])
pairs = []
for i, (idx_a, role_a, ref_a, pa) in enumerate(polys):
    for j in tree.query(pa):
        if j <= i: continue
        idx_b, role_b, ref_b, pb = polys[j]
        try: inter = pa.intersection(pb)
        except Exception: continue
        if inter.is_empty: continue
        a = inter.area
        if a > 0.0:
            pairs.append((a, idx_a, role_a, ref_a, idx_b, role_b, ref_b))
pairs.sort(reverse=True)
for p in pairs[:20]:
    print(p)
PY
```

## Recommended next session — Tier 1 root-cause investigation

The three Tier 1 tests fail across all three airports.  My guess
(unverified) is they share a root cause in junction-emission or
post-emission cleanup.  Strategy for the next agent:

1. **Read each failing test's body** to capture exact assertion
   semantics:
   * `test_no_vertex_on_sloping_rect_edge` — for each rect, every
     non-rect polygon vertex within `TOL_M` of the rect must land
     on a corner, not edge-interior.
   * `test_junction_no_long_edge_proximity` — junction polygons
     must not run along (≥ 10 % of axial length) a rect's long
     edge.
   * `test_junction_neighbour_corners_shared` — adjacent junctions
     sharing a boundary must agree on every vertex along it.

2. **Build each baseline airport with elevations enabled** and dump
   the offending geometry (which junctions, which rects, which
   vertices) at each airport.  Look for a shared pattern.

3. **Locate where the offending geometry is emitted.**  Likely
   candidates (in priority order):
   * `pavement/junction_emit.py::emit_junctions_and_finalize` —
     where junction residue is finalised after rect emit.
   * `pavement/junction_repair.py` — post-emission cleanup.
   * The `overlap-clip pass` and `merged sliver junctions` passes
     mentioned in the build log.
   * `_resolve_runway_crossings` — for the runway-crossing case.

4. **Propose ONE fix that resolves all 9 Tier 1 failures across all
   3 airports simultaneously.**  Validate against the full test
   suite.  If the fix can't be made universal, identify what
   makes each airport's case special and consider whether a
   per-airport regression baseline is more appropriate than a
   code change.

5. **Update `status.md` again** with results and remaining backlog.

## Memory notes — read before continuing

* `feedback_general_solutions.md` — every fix must work at all
  baseline airports.  No airport-specific code.  Per user
  2026-05-16: "any changes we make to fix bugs at any particular
  airport always have to take into account whether the change will
  work across all airports."
* `feedback_root_cause_only.md` — fix root causes, ASK before
  band-aid post-process clean-up.
* `feedback_grade_rules.md` — apron / junction grade rule, with
  the user 2026-05-14 refinement that all-pair only applies
  between pairs whose straight line stays inside the polygon.
* `feedback_shape_rules.md` — authoritative rect + junction
  construction rules.
* `project_refactor_state.md` — extraction-pattern recipe for
  pulling code out of the pavement monolith.

## Files touched this session

Source:
* `src/auto_patch/apt_dat_reader.py` — `taxi_centerlines` pre-splits
  at junction nodes via new `_split_polyline_at_junction_vertices`
  helper.
* `src/auto_patch/pavement/centerlines.py` — `at_endpoint` check in
  `split_merged_centerline`; `CHART_JUNCTION_MARGIN_M` + `_is_chart_
  junction` in `_split_centerlines_at_points`; fallback margin for
  short diagonal stubs.
* `src/auto_patch/pavement/rects.py` —
  `_rect_long_edges_at_pavement_boundary` helper called from
  `_build_taxi_rects`.
* `src/auto_patch/pavement/stubs.py` — perpendicular-skip in
  `_emit_primary_parallel_runway_stubs` (new `rwy_centerlines`
  parameter).
* `src/auto_patch/pipeline.py` — passes `rwy_centerlines` to the
  stubs emitter; reverted the chart-junction-share overlap-drop
  exemption that was tried mid-session.

Tests:
* `tests/conftest.py` — `baseline_airports()` helper.
* `tests/test_junction_invariants.py` — `_test_airports()` + flip
  all four `airports_under_test() or [skip]` parametrisations.
* `tests/test_junction_rules.py` — same.
* `tests/test_pavement_geometry.py` — same; `_BASELINE_AIRPORTS`
  removed (lifted to conftest).
* `tests/test_pavement_grade.py` — parametrise over
  `baseline_airports()`; CYXY caps added.

## Build commands (unchanged)

```bash
venv/bin/python3 -m pytest tests/ --tb=short -q

# Single-airport build for visual inspection:
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement(
    "CYXY", "/Users/noah/X-Plane 12", compute_elevations=True)
layout.to_osm("/tmp/CYXY.osm")
PY
```
