# Auto-Patch Status — 2026-05-18 session: canonical-point registry + apron classification + DEM seam infeasibility

## TL;DR

**14 commits this session.  Tests went 16 → 11 failures (248/264 → 253/264).**

Major architectural shift: introduced a **single canonical-point
registry** as the source of truth for vertex identity across every
shape-mutating pass and the OSM emitter.  Replaces a tangle of
discrete-bucket lookups (`_corner_elevation_bucket` /
`int(round(x/0.5))`) that suffered from bucket-boundary aliasing.

Five other related architectural changes:

* **Apron reclassification pass** — junctions whose boundary
  strays > 55 m from any taxi/runway centerline are reclassified
  as `role=apron`.  Threshold widened from 20 m → 55 m to
  accommodate normal fillet-curve geometry at large airports.
* **`ROLE_RUNWAY_CROSSING`** — new role for the polygon
  `_resolve_runway_crossings` builds at runway-runway
  intersections.  Distinct from generic junction: HARD-anchored
  via runway-interpolated `node_altitudes`, never reclassified to
  apron, never subdivided, counts as a source for adjacent
  shapes.
* **Per-corner apron writeback** — aprons keep the solver's
  per-vertex `node_altitudes` (which the solver constrained to
  ≤ 1.5 % via all-pair Euclidean edges), instead of averaging to
  a single flat altitude.  The previous averaging caused
  adjacent aprons to diverge into 4-8 m cliffs at their shared
  corners.
* **`APRON_MAX_GRADE` = 1.5 %** (was 1.0 %).  Aligned with
  taxi grade so reclassified shapes stay feasible.
* **`test_pavement_grade` cap = 0** (was per-airport baselines).
  Hard fail on any within-shape grade violation.

**Test result: 253 / 264 pass at HEAD.  11 failures remain.**
3 of those are airport-specific compare-target fixtures (by
design).  Of the other 8:

* 3 are `test_pavement_grade` failures — the underlying root
  cause is DEM-noisy SEAM altitudes at SPLP (see below).
* 5 are pre-existing backlog items unchanged from prior session.

## Commits this session (in order)

| Hash | What |
|---|---|
| `fc5ed6c` | Apron reclassification pass + shared-vertex altitude invariants.  Junctions whose boundary strays > 20 m (then bumped to 55 m) reclassified as apron.  OSM emitter altitude-aware: same XY bucket but Δalt > 1 m → separate node IDs (cliff preserved). |
| `ff57ee1` | Preserve apt.dat row-110 taxi-network centerlines on `layout.apt_taxi_centerlines`.  Surviving rect `source_axis` covers only ~40 % of original centerlines at SPJC; without the apt.dat set, legitimate junctions sitting on absorbed centerlines got misflagged as apron. |
| `0bc751d` | Per-centerline junction grade + 55 m apron threshold.  *Reverted in `9b3a77e`*. |
| `1670a56` | Recognise apt.dat row-110 vertices as legitimate junction sources.  Captures `layout.apt_pavement_vertices` and `layout.apt_pavement_boundary` at pavement-union build time; `test_junction_vertices_have_source` accepts row-110 vertices + on-row-110-edge points.  SPJC orphans 39 → 2. |
| `1cfea41` | Add canonical-point registry.  `CanonicalPointRegistry` in `canonical_points.py` with spatial-index `get_or_add`/`find_nearest` at `tol = SHARED_VERTEX_TOL_M`.  Rect-corner snapping routed through it. |
| `ffb8ffc` | Route every shape-mutating pass through the registry.  Sites updated: `junction_emit.py:emit_junctions_and_finalize`, `junction_repair.py:_absorb_rects_at_junction_perimeters`, `_subdivide_violating_junctions`, `_split_sloped_rects_at_violations` (both sub-rect-corner and junction-vertex-move sites + modified-junction-polygon re-route). |
| `d44c28e` | Route `_resolve_runway_crossings` output through registry.  Runway-crossing junction's corners now canonical so adjacent rects/junctions snap to them. |
| `c30de9a` | Add `ROLE_RUNWAY_CROSSING`.  `_resolve_runway_crossings` emits the new role.  Solver treats it as HARD-anchored ring-only (runway-tier).  Overlap-clip puts it in the RUNWAY tier.  Test `SOURCE_ROLES` includes it. |
| `1097101` + `df65f41` | Remove accidentally-committed CYXY DEM tiles; add `.gitignore` for `+lat-lon/` cache dirs. |
| `9b3a77e` | Revert junction grade to all-pair (user 2026-05-18 clarification).  "Junctions should not exceed 1.5 % across ANY portion, not just along an edge."  Reverts the per-centerline edge generation in `unified_jacobi._build_edges` and the matching changes in `_subdivide_violating_junctions._worst_grade` and `_report_within_shape_violations`. |
| `e5eb3b3` | Apron writeback per-corner.  `_writeback` now writes `node_altitudes` for `ROLE_APRON` (keeps per-vertex altitudes the solver enforced ≤ 1.5 %) instead of averaging to a single flat altitude.  Terminals stay flat-averaged.  SPJC adjacent-apron 8.6 m cliffs eliminated. |
| `d310dcd` | Solver uses canonical-point registry for node identity.  `_build_node_list` and every other solver site (`_seed_elevations` runway-HARD pass, seam pass, warm-start; `_build_edges`; `_build_rect_cross_section_groups`; `_build_terminal_groups`; `_read_corner_elevs`) routes its key lookup through `layout.canonical_points.get_or_add`. |
| `03a5b77` | Remove `_solver_key` fallback.  Bucket fallback was dead code preserving the buggy discrete-bucket system.  Direct `layout.canonical_points.get_or_add` calls everywhere; no `_solver_key` indirection. |
| `72efbce` | Revert flatten-on-reclassification.  Per-pass altitude tracing revealed `_reclassify_apron_junctions` was AVERAGING each junction's per-corner field to a single flat altitude on role change — adjacent junctions sharing a corner that the solver had at consistent 26.1 m diverged to 22.2 / 30.8 / 26.6 m when each averaged independently.  Cross-shape violations dropped from 42/8/92 at SPJC/SPLP/CYXY to 0/0/1. |
| `fe58280` | Migrate 3 high-priority bucket sites to canonical-point registry.  (1) OSM emitter `_intern` in `layout.py` — every emitted shape vertex now uses registry; emitted lat/lon is the CANONICAL coords (not the input).  (2) `_snap_junction_altitudes_to_rect_corners` in `elevation.py` — `rwy_corner_alt` keyed by canonical.  (3) `_corner_elev_map` in `elevation.py` — consumed by `triangulation.py:_resolve_corner_elev`, now canonical-keyed. |

## Architectural foundation now in place

### Canonical-point registry (`canonical_points.py`)

`CanonicalPointRegistry`:

* Created once at pipeline init (`pipeline.py:1762`), seeded with
  apt.dat row-110 pavement polygon vertices + runway corners.
* Stored on `layout.canonical_points`.
* `get_or_add(x, y)` returns the canonical (x, y) for an input
  point — proximity-based lookup at `tol = SHARED_VERTEX_TOL_M`
  (= 0.5 m).  Adds a new canonical entry only if no existing
  point is within tol.  Two corners 0.002 m apart that would
  have landed in adjacent discrete buckets now resolve to the
  same canonical point.

Single source of truth: every site that needs to identify
"is this vertex the same as that one" goes through the registry.

`snap_polygon_through_registry(poly, registry)` helper used by
junction-mutating passes to route their output polygons' vertices
through the registry (preserving shared corners after `buffer(0)`
/ `unary_union` / interpolation).

### Vertex routing — sites converted

Solver:
* `_build_node_list` — node identity is canonical-point tuple
* `_seed_elevations` — HARD-anchor seed (runway/runway_crossing,
  seam) keys by canonical
* `_build_edges` — node-index lookup canonical
* `_build_rect_cross_section_groups` — canonical
* `_build_terminal_groups` — canonical
* `_read_corner_elevs` — canonical
* OSM emit `_intern` — canonical, emits canonical lat/lon
* `_snap_junction_altitudes_to_rect_corners` — canonical
* `_corner_elev_map` (consumed by `triangulation.py`) — canonical

Junction-mutating passes (route output polygons through
`snap_polygon_through_registry`):
* `junction_emit.emit_junctions_and_finalize` (sliver-corner
  cleanup after `buffer(0)`)
* `_absorb_rects_at_junction_perimeters` (junction extension via
  `unary_union`)
* `_subdivide_violating_junctions` (sub-polygon validation)
* `_split_sloped_rects_at_violations` (sub-rect corners
  registered, junction-vertex move targets registered, modified
  junction polygon re-routed)
* `_resolve_runway_crossings` (runway-crossing-junction's corners
  become canonical sources)

Bucket sites that legitimately stay discrete (audit complete):
* `groundside.py:340` — internal to one polygon
* `junction_repair.py:122-279` (`_clamp_junction_free_vertices`)
  — spatial grid for radius-based clamping
* `triangulation.py:111, 607-627` — spatial indexing for
  cross-junction smoothing
* `junction_rules.py:1377, 1445` — internal lookup tables
* Several `elevation.py` sites — same-shape contexts

## Investigation outcomes

### SPJC cross-shape mismatches resolved

Initial test_pavement_grade[SPJC] showed 81 within-shape and 42
cross-shape violations.  Per-pass altitude tracing identified the
root cause: `_reclassify_apron_junctions` was averaging each
junction's `node_altitudes` to a single altitude on role change.
At a shared corner that the solver had given altitude 26.1 m to
every adjacent junction, each junction's INDEPENDENT average
produced different values (22.2, 30.8, 26.6 m) → 4-8 m cliffs.

Fix: reclassification is now role-only — the solver's per-corner
field is preserved.  SPJC went 42 → 0 cross-shape, 81 → 2
within-shape.

### SPLP within-shape: DEM seam infeasibility

The remaining SPLP `test_pavement_grade` failure is a **DEM
input problem**, not a solver bug.  Instrumentation traced
junction #58 (160 k m², ancestor of the failing apron #33):

* 13 of its 157 corners are HARD-anchored as **seam vertices**
  (tile-boundary DEM samples, per the user's "tile boundary is
  HARD" rule).
* The seam altitudes themselves violate 1.5 % between adjacent
  HARD vertices:
  * y=877 alt 67.7 m
  * y=967 alt 74.0 m  (+6.3 m over 90 m = 7 %)
  * y=982 alt 67.7 m  (-6.3 m over 15 m = **42 %**)
  * y=1026 alt 74.0 m  (+6.3 m over 44 m = 14 %)
* The solver propagates 1.5 % cap from these HARD anchors but
  cannot reconcile them with each other — the inter-anchor
  grade is fixed.
* Subdivision (`_subdivide_violating_junctions`) splits the
  junction but each sub-piece inherits the seam HARD anchors at
  its end and can't satisfy 1.5 % internally.

The DEM at SPLP's tile boundary is genuinely noisy (6.3 m
altitude jumps over 15 m at lat=-13 / lon=-77 boundary).  Three
real fixes for this:

1. **Pre-smooth seam altitudes** along the seam to satisfy 1.5 %
   between adjacent seam vertices before HARD-anchoring.  Must
   be coordinated with Ortho4XP's terrain mesh which also pins
   the seam from DEM — otherwise creates a pavement-vs-terrain
   cliff at the tile boundary.
2. **Don't HARD-anchor seam vertices** — treat them as SOFT,
   let the solver pull them toward feasibility.  Pavement at
   the seam may end up ≠ terrain at the seam (visible step).
3. **Better DEM source** at SPLP.  6 m altitude jumps over 15 m
   on an airport are unrealistic.

Has not been chosen yet — needs user direction.

## Remaining 11 failures — classification

### Category A: DEM seam infeasibility (3 tests)

* `test_pavement_grade[SPJC]` — 2 within-shape violations.
  Worst 5.36 % / 13 m at apron `-10190` near
  lat=-12.0249, lon=-77.1228.  Smaller-scale seam-style issue.
* `test_pavement_grade[SPLP]` — 38 within-shape violations.
  Worst 30.76 % / 14 m at apron `-10032` (former junction #58).
  Direct symptom of the SPLP DEM seam issue described above.
* `test_pavement_grade[CYXY]` — 1 cross-shape violation,
  1.7 m step at the runway / runway-crossing boundary.  Likely
  unrelated to the seam issue.

### Category B: Pre-existing SPLP geometric backlog (3 tests, unchanged from prior session)

* `test_taxi_rects_not_alongside_apron[SPLP]` — primary parallel
  A shares 26.1 m of sloping edge with apron #32 at lat ≈
  -12.156, lon ≈ -76.999.  Partial-absorption corner-wrap clamp
  excludes it.
* `test_no_vertex_on_sloping_rect_edge[SPLP]` — runway 02/20
  has a 3-corner sliver + 5-corner pentagon at lat ≈
  -12.1644.  From `_drop_overlap_against_fixed_shapes._clip_keep_largest`.
* `test_rect_short_edges_connect[SPLP]` — stub A end_B
  unshared at lat=-12.161382 and -12.161122 (both ends
  disconnected from any neighbour).

### Category C: Junction-rule edge case (2 tests)

* `test_junction_no_long_edge_proximity[SPJC]` — 1 junction
  vertex within ~5 m perpendicular of a sloping rect's long
  edge interior.  The `_split_sloped_rects_at_violations` pass
  should catch this; investigate why it doesn't fire.
* `test_junction_vertices_have_source[SPJC]` — 2 orphans, both
  at the same shared corner (527.35, -819.01) between
  junctions #172 and #174.  Sits 0.76 m off the apt.dat
  pavement boundary — leftover from a downstream pass's
  `buffer(0)` rounding.

### Category D: Airport-specific compare-target fixtures (3 tests, by design)

* `test_compare_target_spjc`
* `test_compare_target_splp[-13--77-baseline0-150]`
* `test_compare_target_splp[-13--78-baseline1-196]`

These compare role/count distributions against hand-written
`*_target.osm` fixtures.  Failing because apron reclassification
shifted junction → apron counts.  Need either fixture re-anchor,
or refactor to a parametric `@parametrize("(icao, target_path)")`
form (you said earlier you want NO airport-specific tests).

## Key memory files

* `feedback_general_solutions.md` — every fix must work at all
  baseline airports.
* `feedback_root_cause_only.md` — fix root causes; ASK before
  band-aid post-process clean-up.
* `feedback_shape_rules.md` — authoritative rect + junction
  construction rules.
* `feedback_grade_rules.md` — grade rule (since 2026-05-18: all
  pavement roles share 1.5 %, all-pair within-shape).
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
    layout.to_osm(f"/tmp/{icao}.osm")
PY
```

## Files touched this session

Source:
* `src/auto_patch/canonical_points.py` (NEW) —
  `CanonicalPointRegistry` + `snap_polygon_through_registry`.
* `src/auto_patch/layout.py` — `ROLE_RUNWAY_CROSSING`;
  `apt_pavement_vertices` / `apt_pavement_boundary` /
  `canonical_points` / `apt_taxi_centerlines` on
  `PavementLayout`; OSM emitter `_intern` routed through
  registry; apron writeback per-corner.
* `src/auto_patch/pipeline.py` — registry creation +
  seeding before rect build; apt.dat row-110 vertex + boundary
  capture; apron-reclassification call site; passing registry
  to rect builder.
* `src/auto_patch/junction_emit.py` — sliver-corner cleanup
  routed through registry.
* `src/auto_patch/junction_repair.py` — apron-reclassification
  function; cross-shape clamp restricted to same-shape;
  `_absorb_rects` / `_subdivide_violating_junctions` /
  `_split_sloped_rects_at_violations` routed through registry.
* `src/auto_patch/pavement/rects.py` — registry seeded with
  pav_union + runway corners; corner-snap routes through it.
* `src/auto_patch/pavement/runways.py` —
  `_resolve_runway_crossings` emits `ROLE_RUNWAY_CROSSING`;
  output polygon routed through registry.
* `src/auto_patch/elevation.py` — `APRON_MAX_GRADE` 1.0 % →
  1.5 %; `_snap_junction_altitudes_to_rect_corners` and
  `_corner_elev_map` canonical-keyed; `_clamp_junction_free_vertices`
  cross-shape branch removed; runway-crossing in overlap-clip's
  runway tier.
* `src/auto_patch/elevation_per_surface/unified_jacobi.py` —
  every solver site canonical-keyed (no bucket fallback);
  `PAVEMENT_ROLES` includes `ROLE_RUNWAY_CROSSING`; apron
  writeback per-corner; HARD-anchor seed for runway +
  runway_crossing.
* `src/auto_patch/elevation_per_surface/solver.py` — wraps
  `_jacobi_solve` (no behaviour change this session).
* `src/auto_patch/triangulation.py` — corner-elev lookup
  canonical-keyed.

Tests:
* `tests/test_junction_invariants.py` — `_aeroway_centerlines_m`
  delegates to pipeline helper; `test_junction_vertices_have_source`
  replaces `test_junction_vertex_count_bounded`; row-110
  vertices + boundary accepted as sources; SOURCE_ROLES includes
  `runway_crossing`.
* `tests/test_elevation_terrain_following.py` —
  `APRON_MAX_GRADE` updated to 1.5 %.
* `tests/test_pavement_grade.py` — `WITHIN_SHAPE_CAP = 0` for
  every airport (no soft cap, no per-airport baseline).

## What's next (for handover)

The architectural foundation is in place.  The biggest remaining
decision is the **SPLP DEM seam handling** (Category A above):
pre-smooth seam altitudes, leave seam SOFT, or improve DEM
source.  The other 8 failures are smaller-scope fixes (backlog
items, fixture re-anchor, junction-rule edge cases).

Priority order if continuing:

1. **DEM seam smoothing** (Category A) — pick one of the three
   approaches in the SPLP within-shape section above.  Likely
   unblocks all 3 `test_pavement_grade` failures.
2. **`test_junction_no_long_edge_proximity[SPJC]`** — quick fix
   (1-2 h) once you trace why `_split_sloped_rects_at_violations`
   doesn't fire on this specific junction-vertex.
3. **SPJC orphan corner** at (527.35, -819.01) — 0.76 m drift
   from buffer rounding in a downstream pass.  Snap-to-pavement-
   boundary fix likely.
4. **SPLP backlog** (3 tests) — pre-existing items, separate
   work.
5. **`test_compare_target_*` fixture handling** — your call:
   re-anchor floors, refactor to parametric, or delete.
