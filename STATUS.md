# Auto-Patch Refactor — Status

**Current state:** New emitter `src/O4_Airport_Pavement_Builder.py`
(Phase 1: shapes + roles only) being built to reproduce hand-drawn
target OSMs at `tests/fixtures/{SPJC,SPLP}_target.osm`.  Target
tolerance is **5 m vertex-match for runway** (already passing at
<0.3 m after target resnap) and **1 m for non-runway shapes with
freedom to snap to apt.dat pavement edges**.

This supersedes the apron-deformation model below; the elevation
work is Phase 2.

## Phase 1 — target-driven layout builder

### Infrastructure
- **`tools/compare_target.py`** — parses target + output OSMs,
  reports per-role IoU matching + symmetric vertex-match distance
  distributions.  Success metric is: max vertex distance per matched
  shape ≤ tolerance.
- **`tools/build_target_osm.py`** — driver that runs
  `build_airport_pavement(icao, xplane_root)` and writes OSM.
- **`src/O4_Airport_Pavement_Builder.py`** — the new builder.
  Current implementation:
  1. Emit runway rects (apt.dat row 100 + blast pads).
  2. Load OSM taxi centerlines from tile cache.
  3. Clip each centerline to pavement (pavement union minus runway).
  4. Probe rect width as 2× median distance-to-boundary along axis
     (ray-cast is unusable at SPJC where pavement is one big blob —
     perpendicular rays overshoot into apron and saturate at 120 m).
  5. Extend rect endpoints outward until they leave the pavement,
     so corners land on apt.dat pavement boundaries.
  6. Classify role from angle-to-runway + length + proximity.
  7. Dedup overlapping rects (≥ 70 % inside emitted union).
  8. Aprons = pavement union minus emitted taxis, per connected
     component.

### Target file conventions
- 8 role values (+ 3 unused on SPLP):
  `runway, primary_parallel, secondary_parallel, stub,
   cross_connector, junction, apron, terminal`.
- **Shared vertices pairwise only** (162 of 470 SPJC nodes shared
  by exactly 2 ways, 0 by ≥ 3).  Pairwise shared-vertex invariant
  enforces geometric consistency between adjacent shapes.
- Altitude tags (`altitude`, `altitude_high`, `altitude_low`) are
  Phase 2 — 9 elevation-only ways in SPJC target ignored by Phase 1.
- Runway corners in both targets **were resnapped to apt.dat
  row-100 + blast-pad corners** (user authorized; runways are
  reference only since `generate_patch_osm` emits sloped segmented
  runways in production).

### Current scores (commit `1bdae1c` as baseline, HEAD as current)

```
                        SPJC (72 target shapes)      SPLP (29 target shapes)
                        baseline  current            baseline  current
runway                  0/2       2/2 (IoU 1.00)     0/1       1/1 (IoU 0.99)
primary_parallel        1/9 (0.11) 5/9 (0.40)        0/3       0/3
secondary_parallel      1/4 (0.20) 0/4               0/0       n/a
cross_connector         1/2 (0.12) 2/2 (0.56)        0/2       2/2 (0.24)
stub                    5/15(0.43) 9/15 (0.23)       0/6       2/6 (0.19)
apron                   0/3       3/3 (0.57)         0/3       1/3 (0.33)
junction                0/35      0/35 (not emitted) 0/14      0/14
terminal                0/2       0/2 (not emitted)  —         —
TOTAL matched           8         21                 —         6
```

Runways pass at 0.2–0.3 m.  Everything else has max-vertex-
distance of **60–180 m** on the best matches — well above the 1 m
target.  Main remaining gaps (next-session priorities):

1. **Ref consolidation.** SPJC has sub-refs A1–A6, D1/D2, F1, M1–M3,
   R1/R2 that should merge into their parent refs.  Currently each
   sub-ref emits its own rect, creating 46 stubs vs 15 target.
2. **Rect endpoint trimming.** Rects extend to pavement edge even
   when the target rect ends at a neighbour (stub vs parallel).
   Primary_parallel V: 61 m max vertex distance because my rect
   extends full length of pavement while target stops at stub
   intersection.
3. **Junction emission.** 0/35 SPJC + 0/14 SPLP.  These are the
   filler polygons between rects — target has them explicit.  Need
   to compute the gap polygon and emit with shared corner vertices.
4. **Apron consolidation.** 59 output aprons vs 3 target — each
   apron-residue fragment is emitting separately instead of unioning
   into the big terminal/engine-test regions.
5. **Terminal emission.** From OSM building data; 0 emitted so far.
6. **SPLP primary parallels.** Only 1 emitted from OSM, target has
   3.  SPLP OSM has 15 taxi centerlines, none with refs — need
   smarter selection (length, parallelism, distance to runway).

### Files touched this session
- New: `src/O4_Airport_Pavement_Builder.py`, `tools/compare_target.py`,
  `tools/build_target_osm.py`.
- Modified: `src/O4_Apt_Dat_Reader.py` (X-Plane 12 `Global Scenery/
  Global Airports` path added to `find_airport_apt_dat`).
- Modified: `tests/fixtures/SPJC_target.osm`, `SPLP_target.osm` —
  runway corner nodes resnapped to apt.dat coordinates.  No other
  nodes touched.

### How to resume
```
./venv/bin/python3 tools/build_target_osm.py SPJC
./venv/bin/python3 tools/build_target_osm.py SPLP
./venv/bin/python3 tools/compare_target.py \
    tests/fixtures/SPJC_target.osm /tmp/SPJC_auto.osm --tol 1 -v
./venv/bin/python3 tools/compare_target.py \
    tests/fixtures/SPLP_target.osm /tmp/SPLP_auto.osm --tol 1 -v
```

## Legacy status (Phase 2 territory, not currently active)

(original content below — apron-deformation model, bottom-up merge,
legacy pipeline history — all still relevant for Phase 2 elevation
work.)

## Active plan — apron-deformation reactive model

The user converged on this model over the course of this session;
it supersedes the prior junction-carving model (§5 below) and
every approach in §1-4.  Implementation is staged.

### Rules (all user-authored, non-negotiable)

1. **No strip buffers.**  Each apt.dat pavement polygon IS a
   shape.  Narrow pavement (short side ≤ 30 m) classifies as
   **taxi**; wider pavement classifies as **apron**; runway
   polygons are handled by `generate_patch_osm`, not here.
2. **Taxi shapes have a 1D axial slope** (≤ 1.5 %, relaxable to
   3.0 % when otherwise infeasible).  Axis = MRR midline for
   simple-strip polygons, Voronoi centerline(s) for mega-polys.
3. **Apron shapes are 2D deformable surfaces** at ≤ 1.0 % in any
   direction.  No transition wedges are emitted at apron-apron or
   apron-taxi or apron-runway joins — the apron itself reshapes
   so its boundary elevations match every neighbor exactly.
4. **Shared-boundary vertices have identical elevations.**  All
   vertex altitudes are quantised to **0.5 m increments**; at any
   point where two shapes share a vertex, both emit the same
   quantised value.  "0 m tolerance" means emitted-string
   equality, not float equality.
5. **The last 30 m of a taxi polygon may be a triangle zone** for
   compound-slope blending with an intersecting shape whose axis
   differs.  Triangle zones use ≤ 1.0 % grade per the transition-
   triangle rule in STATUS invariant 5.
6. **Yielding priority:** runway never yields; everything else
   negotiates within its grade cap.  Buildings yield within
   `±BLDG_ADJUST_MAX`; aprons yield via 2D deformation; taxis
   yield via axial profile.  Terminal pads force unrealistic
   flat zones — taxi grade is allowed up to 3.0 % to absorb
   this pressure.
7. **Build order (unchanged from STATUS):** runways → buildings
   → aprons → taxis → reconciliation.

### Data model

- `Shape(polygon, kind, axis, width_m)` — one per (sub-)polygon
  extracted from apt.dat; mega-polys may produce multiple Shapes
  via Voronoi branching.  `kind ∈ {"taxi", "apron"}`.  `axis` is
  `None` for aprons.
- `Adjacency(shape_a, shape_b, shared_boundary)` — the planar
  subdivision's edge graph.  Derived from polygon-boundary
  intersection, not buffered proximity.
- **No Junction type.**  No strip buffers.  No pre-carved discs.
- Shared boundaries are the apt.dat polygon edges themselves
  (possibly after a single unioning pass to merge co-edged
  polygons of the same kind).

### Solver sketch

1. Fix runway segment elevations via CIFP (existing
   `generate_patch_osm`).
2. Fix building pad elevations (DEM centroid, adjustable
   `± BLDG_ADJUST_MAX`).
3. **Solve apron surfaces:** for each connected apron region,
   compute a 2D elevation field satisfying ≤ 1.0 % grade in any
   direction and matching exact elevations at every shared
   boundary with a runway, building pad, or taxi endpoint anchor.
   Interior seed = DEM, smoothed to grade.
4. **Solve taxi axial profiles:** each taxi's 1D profile along
   its axis, endpoints pinned to the apron-edge or runway-segment
   elevation at each axis terminal.  Grade cap 1.5 %; if
   infeasible, relax to 3.0 %.
5. **Reconciliation:** if a taxi still infeasible at 3.0 %,
   relax the adjacent building pad's elevation toward the
   constraint.  Runways never relax.
6. **Quantisation:** round every vertex altitude to 0.5 m.  When
   two shapes share a vertex, they read the same quantised value.

### Emission (role-based pipeline, user-specified staging)

After decomposition, shapes are classified by role
(`classify_shape_roles`) into one of five categories, each with
its own emission treatment:

1. **Primary parallel** (length ≥ 0.25 × runway length, within
   300 m of runway, parallel ± 15°): handed off to the runway
   emitter.  Start / end elevations match the runway's same-side
   elevation at each axis endpoint.  Emitted as a **multi-segment
   sloped rect chain** (same shape type as runway segments), with
   DEM-driven undulations between the CIFP anchors honouring the
   1.5 % grade rule.
2. **Stub** (short taxi adjacent to a primary parallel, axis ≤
   300 m): single sloped rect from the primary parallel edge to
   the runway edge, **leaving an elevation-dependent gap** at both
   ends for the later triangle-join pass.
3. **Secondary parallel** (parallel to runway, but too short /
   too far to be primary): DEM-driven segmented rect chain, same
   pattern as primary parallel but anchored to neighbour shapes
   instead of CIFP.
4. **Cross-connector** (perpendicular to runway, adjacent to two
   parallels): single sloped rect, **stopping short** of both
   parallel intersections.  Gap sized for the triangle-join pass.
5. **Apron:** single flat N-gon when solved elevation range
   < 0.5 m; otherwise triangulated with per-vertex
   `node_altitudes`.  Footprint = polygon minus all primary /
   stub / secondary / cross rects minus their safety gaps.

All shared vertices produce identical `(lat, lon, altitude)`
strings after 0.5 m quantisation.  Gaps between shape categories
will be filled in a later "join" pass by either:

- a flat poly (if both sides agree on elevation to within 0.5 m),
- a single sloped rect (if the slope is 1-axis along the gap),
- or a triangle wedge (compound slope or runway edge meeting taxi
  at different altitude).

### What this supersedes

The previous model had `Junction` objects with pre-carved discs
between meeting strips, plus "transition rect" logic for elevation
differences.  Both are removed.  The apron-deforms-to-fit rule
makes wedges/junctions unnecessary.  This also resolves STATUS's
previously-unsolved problem #3 (transition-rect cascade) because
the rule becomes "apron absorbs Δelev over its whole surface,"
not "insert rect of length Δelev / 0.015 that may cut through
neighbors."

### Intersection handling (derived, not geometric)

Intersections between shapes are **not pre-identified geometric
objects**.  They emerge from the slope-gradient comparison at
each adjacency:

- Apron on either side → apron 2D-deforms to match the neighbour's
  boundary elevation; no extra shape emitted.
- Taxi-to-taxi with similar axial slopes → shapes meet at shared
  vertices directly, no extra shape.
- Taxi-to-anything with incompatible slope axes or nonzero Δelev
  that the apron can't absorb → the taxi's end becomes a
  **triangle zone** (up to 30 m by default) that accepts compound
  slope at ≤ 1.0 %.  The triangle zone is an INTERNAL
  subdivision of the taxi polygon, not a separately-carved shape.

### Validation + iterative fix loop

Because some airport geometries (e.g. forced-flat terminal pads)
cannot be satisfied with default parameters, the emitter runs a
final validation + fix loop:

```
solve_elevations()
for iteration in range(MAX_FIX_ITER):
    violations = validate_joins()
    if not violations:
        break
    for v in violations:
        apply_geometry_fix(v)        # grow triangle zone, shrink axial region
    re_solve_affected_subgraph()
```

Violation types:

1. **Cliff:** two shapes share a vertex with different emitted
   elevations.  Fix: solver bug — should not happen when shared
   vertices use the identical quantised value.
2. **Grade violation in taxi axial region:** edge's slope exceeds
   1.5 %.  Fix: grow the triangle zone at the violating end,
   shortening the axial region.
3. **Grade violation in triangle zone:** compound-slope region's
   internal gradient exceeds 1.0 %.  Fix: grow the triangle zone
   further (up to the full taxi polygon length).
4. **Infeasible at triangle-zone-covers-whole-shape:** fall back
   in this order:
   a. relax taxi grade cap 1.5 % → 3.0 %
   b. relax adjacent building pad ± ``BLDG_ADJUST_MAX``
   c. log an unresolvable violation (runway anchors never move)

**Why this is different from the old greedy-fix loop** (STATUS
§5, which diverged): the old approach patched *elevations*
locally, which introduced new violations at third-party
neighbours.  The new approach modifies *geometry* (triangle-zone
extent) and then **re-runs the solver** on the affected
subgraph.  The solver converges globally even though the loop
runs locally.

### Staged implementation

1. **Decomposition cleanup:** strip the Junction type and all
   strip-buffer / junction-carving code from
   `O4_Pavement_Strips`.  Output is a flat list of Shapes.
   ✅ done
2. **Per-branch local-width decomposition:** mega-polygons are
   split via Voronoi skeleton + median local-width filter;
   branches wider than 30 m (median) remain apron; narrower
   branches become individual taxi shapes with their Voronoi-
   corridor polygon.  Apron residuals union into connected
   components.  ✅ done (60 shapes at SPJC)
3. **Adjacency graph:** compute polygon-boundary adjacencies
   with a 0.5 m tolerance to absorb apt.dat precision gaps.
   Per-boundary-segment, the nearest-neighbour shape claims the
   segment (no double-counting).  🔨 in progress — tolerance
   matching works, nearest-neighbour partition still needed.
4. **Elevation solver:** apron 2D field + taxi 1D profile with
   the hierarchical order in the build-order section above.
5. **Triangulation emitter:** replace the zero-elevation preview
   with actual elevations + triangulated aprons + triangle zones
   at taxi ends.
6. **Validation + fix loop:** the section above.

**Current state:** step 3 (adjacency graph) done with 0.5 m
boundary tolerance.  Steps 1-2 refined extensively; role
classification + multi-segment emission + 10 m taxi-end gap +
apron carving all working.  **Next:** ground-truth-driven
refinement (see below).

## Ground-truth approach — next phase

The user will hand-craft an OSM file representing the desired
output shapes for SPJC: `tests/fixtures/SPJC_target.osm`.  Each
way or multipolygon relation gets:

```
aeroway=taxiway|apron
role=primary_parallel|secondary_parallel|stub|cross_connector|apron
ref=V|L1|A3|...   (optional, human-readable)
width=NN           (optional, for taxis)
```

No elevations — this is a **geometry target** only.  The user
will draw it in JOSM using the SPJC Jeppesen chart as reference.

### What will be built when the target arrives

1. `tools/compare_target.py` — parses both target and current
   output, matches shapes by best IoU, reports:
   - matched (IoU ≥ 0.8)
   - partially matched (0.3 ≤ IoU < 0.8)
   - missed (no target-side match)
   - spurious (no algorithm-side match)
   - per-role F-score
2. `tests/test_target_match_spjc.py` — runs comparison, asserts
   F-score threshold.  Fails the build on regression.
3. **Iterative tuning:** algorithm has ~10-15 tunable knobs
   (NARROW_WIDTH_M, skeleton tolerances, collinear angles, min
   trunk length, accept_median_factor, etc.).  With ground truth
   we can grid-search or manually iterate against the F-score.

### Data sources the algorithm may draw from

- `apt.dat` pavements (primary source currently)
- `apt.dat` runways (for bearing hints, already used)
- **OSM `aeroway=taxiway` ways** — SPJC has 234 of these with 91
  carrying `ref` tags covering A through V.  Previously rejected
  as primary source; worth using as **secondary signal** to seed
  decomposition with known-good centerlines.  apt.dat pavement
  polygons then constrain footprints.
- OSM `aeroway=apron` / airport boundary polygons — bounding
  context.

### Caveats to keep in mind

- **Overfitting to SPJC.** Once SPJC matches, test on 2-3 other
  airports before tuning is locked.
- **Target must be producible from data.** If a target shape has
  no apt.dat + OSM support, flag the gap, don't try to
  hallucinate geometry.
- **IoU tolerance.** Realistic: 0.85+ for taxis, 0.80+ for
  aprons.  Exact vertex match is not achievable.

## Where the code is at end of this session

### Module `src/O4_Pavement_Strips.py`

- `Shape(polygon, kind, axis, width_m)` — one per decomposition
  unit.  kind ∈ {"taxi", "apron"}.
- `Adjacency(shape_a, shape_b, shared, length_m)` — boundary
  graph.
- `decompose_pavement(taxi_polys, apron_polys, …)` → tuple of
  Shapes.  Classification is GEOMETRIC (based on polygon MRR
  short ≤ NARROW_WIDTH_M = **45 m**) + mega-poly decomposition
  via Voronoi skeleton.
  - Taxi extraction: skeleton → longest-path trunk merging with
    runway-bearing preference → local-width filtering
    (accept_median_factor = 1.0, mix_median_factor = 1.5) →
    claim-based polygon carving.
  - Apron: connected components of apron-class polygons + mega-
    polygon residuals.
- `build_adjacency_graph(shapes, tolerance_m=0.5)` — pairwise
  boundary-to-polygon-within-tolerance adjacency.
- `classify_shape_roles(shapes, adjacencies, runway_cls)` — 5
  roles: primary_parallel (≥0.25× runway length), stub, cross_
  connector, secondary_parallel, apron.
- 39 unit tests pass.

### Emitter `_emit_pavement_strip_model` in `src/O4_Auto_Patch.py`

Gated by `O4_STRIP_MODEL=1`.

- Drops apt.dat polygons whose extent is ≥50 % inside runway.
- Subtracts runway union from both taxi and apron pools.
- Passes runway bearings to decompose_pavement for trunk bias.
- For each taxi shape (sorted longest-first to claim pavement):
  - Trims axis by 10 m at each end (GAP_M).
  - Primary/secondary parallels: subdivide into ~100 m segments.
  - Stubs / cross-connectors: single rect from trimmed axis.
  - Builds clean rectangles from (axis, width).
  - Clips each rect to the shape's claimed polygon.
  - Skips rects that would overlap already-claimed taxi area.
- Apron emission: each input apron polygon MINUS the union of
  emitted taxi rects.  No 10 m buffer on long sides (user rule).
  Apron retains the 10 m gap zones at taxi ends.
- Preview elevation is 0.0 m throughout; solver comes later.

### Known issues at handoff

- Apron still over-fragments (~40 pieces at SPJC) where taxis
  carve apart a single big input polygon.  Consolidation pass
  not yet written; candidate: union adjacent apron pieces within
  ~30 m of each other.
- Taxi "A" (east parallel at SPJC) is absorbed into the terminal
  apron because apt.dat classifies it as apron-type pavement.
  Extraction needs either (a) skeleton pass on the apron pool
  too, or (b) OSM centerline seeding.
- L splits at the engine-test-apron jog into two/three primary/
  secondary fragments rather than one continuous L.  The trunk
  extractor's collinear threshold (60°) is enough for small
  bends but not for an explicit jog.  Acceptable per user.
- The audit script `/tmp/check_shapes.py` uses lossy deg→m area
  conversion and reports false-positive overlaps.  In
  meter-space the decomposition has ZERO overlaps.

### Files produced by end-of-day run

- `/tmp/SPJC_shapes.osm` — preview output from
  `tools/dump_shapes_osm.py` (kept in `/tmp/` for now).
  172 ways + multipolygon relations, tagged with role.
- `/tmp/SPJC_legacy.patch.osm` — full pipeline output with
  `O4_STRIP_MODEL=1`.
- `/tmp/spjc_osm_plot.png` — colour-coded role plot.
- `/tmp/run_legacy.py` — launcher for the pipeline.

### How to reproduce at start of next session

```bash
cd /Users/noah/Ortho4XP-shred86
./venv/bin/python3 -m pytest tests/            # 139 pass
O4_STRIP_MODEL=1 O4_DEBUG_TAXI_ONLY=1 \
  ./venv/bin/python3 /tmp/run_legacy.py        # full pipeline
./venv/bin/python3 /tmp/dump_shapes_osm.py     # shapes-only preview
./venv/bin/python3 /tmp/plot_osm.py            # colour-coded PNG
```

## Earlier bottom-up-merge model (retained for reference)

**Commits on dev branch:**
- `f6d3140` — last stable commit from previous session (legacy
  pipeline, 35 taxi rects at SPJC, 9.6 s runtime)
- `7e6199e` — current session: new model function + gating
  infrastructure + merge-based decomposition

## What was tried this session (and what failed)

### 1. Per-polygon taxi/apron classification (legacy model review)

Investigated why named taxiways (TWA, TWE, TWF, TWV) were missing.
Found that apt.dat polygons don't correspond to named taxiways —
e.g. "Taxiway Aux B, C, D, E, F, G, Main Ramp, RWY 33 Head" is
ONE 160k m² polygon.  "Taxiway V / U / Q / R / L / M" is 858k m².
The morphological decomposition (`buffer(-15).buffer(+15)`) can't
extract individual taxi strips from these mega-blobs.  Taxiway A
doesn't even have its own polygon — it's just paint inside "Base
Ramp."  **Conclusion: per-polygon classification is a dead end.**

### 2. OSM centerline-driven emission

OSM has 234 `aeroway=taxiway` ways for SPJC, 91 with `ref` tags
covering all named taxis (A through V).  Prototyped using these
as the emission source with apt.dat pavement as a clip mask.
**Abandoned** because: OSM may not exist at small airports, OSM
and apt.dat don't always align geometrically, and it re-introduces
the buffer-precision problems that caused the original move to
apt.dat polygons.

### 3. Flat-fill + transition-strip model

Decomposed `pavement_region = apt.dat_union − runways − buildings`
into connected components, emitted each as a flat N-gon with
sloped 4-vertex rects at runway edges.  Problems:
- OSM ways can't represent polygons with holes → flat fills with
  building holes required triangulation → 1010 flat triangles.
- Adjacent shapes at different flat elevations → cliffs.
- Phase E1/C4/coverage-fill noise added hundreds more shapes.
Result: 1912 total ways, no improvement over legacy.

### 4. Recursive top-down split

Start with whole pavement, fit one plane, split if fails.
Splitting is elevation-blind → fragments by geometry, not by
elevation planarity.  Produced 210 shapes at first, but with 2.45x
coverage ratio (OBB waste), random elevation cliffs, and no slope
awareness.  Various tuning attempts (OBB ratio, merge tolerance,
simplification) didn't converge to a clean solution.

### 5. Bottom-up merge (current implementation)

**Algorithm:**
1. PREPARE — union apt.dat pavement, subtract runways,
   morphological-close (25 m) + simplify (15 m).
2. SEED — `shapely.ops.triangulate` of boundary vertices + anchor
   points (no densification, no interior grid).  DEM elevation at
   each vertex, grade-clamped to 1.5 %, rounded to 0.5 m.
3. MERGE — greedily merge adjacent coplanar triangles when combined
   region fits one plane (grade ≤ 1.5 %, vertex residual ≤ 0.5 m).
4. EMIT — flat N-gon / sloped 4-vertex rect / 3-vertex triangle.

**SPJC results:** 688 seed → 286 merged → 406 emitted = 479 total
ways (including 73 runway segments).  17 flat, 93 sloped rect,
369 triangles.  **202 join violations remaining (max 2.0 m).**

**Key findings from merge approach:**
- `adaptive_triangulate` (from `O4_Surface_Mesh`) doesn't guarantee
  full coverage at airport scale — it left 27.5 % gaps.  Switched
  to `shapely.ops.triangulate` + centroid-inside-polygon filter
  which achieves 99.99 % coverage.
- Grade-clamping vertex elevations BEFORE merging (the "FAA
  smoothing") dramatically improves merge convergence: 1470 → 48
  with 30 iterations.  But unconstrained propagation flattens the
  whole airport.  A bounded clamp (±2 m from original DEM, 10
  iterations) is the right balance.
- The ">4 vertex → flat at mean" emit fallback is a **cliff
  factory**: it discards the plane's per-vertex elevations and
  emits at a single mean elevation, creating 5 m cliffs with
  neighbors.  Replaced with fan-triangulation from centroid, which
  preserves per-vertex elevations but adds ~200 triangles.
- A greedy post-emit fix loop (insert transition rects at cliffs)
  DIVERGES: each fix creates new violations with third-party
  neighbors.  Violations grew from 127 to 1104 over 10 rounds.
- A pre-emit fix pass (compute all transition rects before emit)
  catches only 8 of 234 violations because most violations come
  from the emit step, not from merge boundaries.
- SPJC is 70 % continuous slope, 30 % flat.  There are no large
  flat zones to discover.  Elevation-quantized flood-fill at 0.5 m
  produces 685 tiny zone fragments, not useful.

### 6. User's structural insight (not yet implemented)

The user described airport pavement as a **graph of flat zones
connected by sloped strips**:

- **Large flat areas** at runway ends and aprons.
- **Sloped rects** connecting adjacent flat areas (length =
  `elev_diff / max_grade` — this is physics, not a tunable).
- **Compound-slope triangles** (~15 m) at runway edges where the
  runway slope axis differs from the pavement, and at building
  pad edges where pad elevation differs from pavement.
- Shared vertices must match **exactly** (not "within tolerance").

The algorithm should build **outward from fixed constraints**
(runways and buildings), not decompose the surface geometrically:

```
Step 1: Emit compound-slope connectors at runway boundaries
        (grade-compliant in all directions, ~15 m but sized
        to meet grade rule)
Step 2: Emit compound-slope connectors at building/terminal
        pad boundaries where needed
Step 3: Between adjacent zones at different elevations, emit
        sloped rects (length = diff / 0.015)
Step 4: Fill remaining area with flat N-gons
Validate: all shared vertices match exactly, all grades ≤ 1.5 %
```

**This approach has NOT been implemented yet.**  It is the
recommended next step.

## Unsolved problems

1. **Join validation:** 202 violations remain at up to 2.0 m.
   Shared vertices between adjacent shapes don't match.  The
   user requires EXACT match (0.0 m tolerance), not "within 0.5 m."

2. **>4 vertex sloped regions:** can't emit as sloped N-gon
   (not human-editable) or as flat at mean (creates cliffs).
   Fan-triangulation preserves elevations but adds many triangles.
   Need a way to decompose these into clean 4-vertex rects.

3. **Transition rect sizing:** when two shapes differ by X meters,
   the connecting sloped rect must be `X / 0.015` meters long.
   This may extend through multiple neighboring shapes, requiring
   those shapes to be cut and their elevations recalculated.

4. **Building pad interaction:** building pads from Phase C3 emit
   at their own elevations after the pavement model.  No compound
   connectors bridge pad-to-pavement elevation differences.

5. **Overlap:** current output has 2.4 % overlap, mostly from fan
   triangulation on concave polygons and transition rect geometry.

## Implementation details

### `_emit_pavement_new_model()` in `O4_Auto_Patch.py`

Module-level function called from inside `generate_airport_surface_
patches` when `O4_NEW_MODEL=1`.  Takes all pipeline closures as
parameters (DEM sampler, emit helpers, runway lookup, etc.).

**Gating:** `O4_NEW_MODEL=1` disables:
- Phase C0 (apt.dat rect-chain), Phase C1 (OSM taxi buffer)
- Phase C2 (apron triangulation), Phase D (junction triangulation)
- Phase C4 (transition strip detection)
- Phase E1 (boundary band), coverage fill, drainage

**Keeps running:** Phase A0.5–A5 (geometry + elevations),
Phase C3 (building pads), `generate_patch_osm` (runways).

### Pavement preparation

```python
pav_region = unary_union(apt_twy_polys_m + apron_polys_m)
pav_region = pav_region.difference(rwy_union_m)
pav_region = pav_region.buffer(25).buffer(-25)   # morphological close
pav_region = pav_region.simplify(15.0)            # RDP simplify
```

The 25 m closing + 15 m simplify was tuned visually against SPJC:
- Absorbs thin arms, jogs, and small interior blobs
- Produces 2 clean connected components with 6 interior holes
- ~170 exterior vertices, 13 % area overflow vs raw apt.dat

### Seed triangulation

Uses `shapely.ops.triangulate(MultiPoint(pts))` with boundary
vertices + runway/building anchor points.  Triangles outside the
polygon or inside holes are discarded via `polygon.contains(
triangle.representative_point())`.  Produces ~688 seed triangles
at 99.99 % coverage.  **No densification, no interior grid** —
boundary vertices + anchors are sufficient.

### Grade clamping

After seed, before merge.  Builds edge adjacency from triangulation.
For each edge where `|z1 - z2| / distance > 1.5 %`, pulls the
non-anchored vertex toward compliance.  Anchored vertices (runway/
building boundary) don't move.  Max adjustment ±2 m from original
DEM.  10 iterations.  All elevations rounded to 0.5 m.

### Merge

Greedy: for each region pair sharing a boundary, try union + plane
refit.  Accept if grade ≤ 1.5 % and max vertex residual ≤ 0.5 m.
Repeat until no merges happen.  Typically reduces 688 → 286.

### Emit

- z_range < 0.5 m → flat N-gon
- 3 vertices → triangle with `node_altitudes`
- 4 vertices, slope in one direction → sloped rect
- \>4 vertices, sloped → fan-triangulate from centroid

### Validation

Post-emit scan of all adjacent shape pairs.  For each pair,
samples the shared boundary midpoint and compares each shape's
interpolated elevation (IDW from vertex elevations).  Reports
violations > 0.5 m.  **Does not fix violations — that's the
next step.**

## How to run

```bash
# Current merge model (behind flag):
O4_NEW_MODEL=1 O4_DEBUG_TAXI_ONLY=1 ./venv/bin/python3 /tmp/run_legacy.py

# Legacy model (default, no flags):
./venv/bin/python3 /tmp/run_legacy.py

# Audit the output:
./venv/bin/python3 /tmp/audit_legacy.py
```

Output goes to `/tmp/SPJC_legacy.patch.osm`.  Copy to
`Patches/-20-080/-13-078/SPJC_auto.patch.osm` for X-Plane.

## Active surface invariants (apply to all emitted geometry)

These are non-negotiable.  They drove the decision to revert.

1. **No overlapping shapes.**  Every emitted way is interior-disjoint
   from every other emitted way.  Shared boundaries are allowed (and
   required by #2); shared interiors are not.
2. **Continuous surface.**  Runways, taxiways, aprons, terminals, and
   buildings together form one continuous painted surface.  Where two
   features meet they must share an actual geometric boundary.
3. **Equal altitudes at joins.**  Where two shapes share a boundary
   (or even a single vertex) their elevations match exactly.
4. **Triangle wedges for compound slope changes.**  When two adjacent
   shapes have incompatible slope axes, the join is bridged by 3-vertex
   triangle patches (`node_altitudes` tag, exactly 3 unique vertices)
   that smoothly transition the elevation.
5. **All FAA/EASA grade rules apply:**  runways ≤ 1.5 % longitudinal
   plus the vertical-curve rate-of-change rule; taxiways ≤ 1.5 %;
   aprons + transition triangles ≤ 1.0 % in any direction; buildings
   are flat pads whose elevation may be adjusted within
   `±BLDG_ADJUST_MAX` from DEM for reconciliation.

These invariants describe a **planar subdivision** of the airport
surface.  Every point belongs to exactly one cell; cells share edges
in pairs; elevations agree across shared edges.

## Authoritative airport elevation build order

Surface geometry can only be emitted after a whole-airport elevation
model has been resolved in this priority order.  Individual features
can each be internally correct and still fail to join smoothly — the
only way to guarantee smooth joins *and* meet FAA/EASA slope limits
everywhere is to solve the elevation field for the whole airport
first.

**Critical principle: CIFP runway elevations are authoritative for the
runway alone.  They are NOT used as anchors for any other feature.**
Propagating CIFP across the airport flattens everything to the runway
slope, which is wrong — taxiways, aprons, and buildings have their own
elevations driven by the local DEM, with the runway being just one of
several boundaries they negotiate against.

1. **Runways:** CIFP threshold elevations are authoritative.  The
   runway emits at those elevations, plus DEM-driven undulations
   between them constrained to ≤ 1.5 % grade.  CIFP elevations stop
   at the runway boundary.
2. **Buildings/terminals:** flat pads at DEM centroid elevation.
   Independent of the runway.  May be adjusted in step 5.
3. **Aprons (1.0 % grade limit):** connect to building edges using
   the building's pad elevation as a hard anchor.  Interior follows
   DEM, smoothed to ≤ 1.0 % in any direction.  Aprons do NOT use
   CIFP runway elevations for any vertex (only at the rare case where
   an apron polygon literally touches a runway boundary).
4. **Taxiways (1.5 % grade limit):** connect smoothly between aprons
   and runways.  Anchored at apron-touch points (using the apron's
   local elevation at that point) and at runway-touch points (using
   the runway segment elevation at that point — only where the
   centerline crosses the runway boundary).  Interior follows DEM,
   smoothed to ≤ 1.5 %.  Taxiways have their own slopes which may
   differ from runway slopes.
5. **Reconciliation:** if step 3 (apron) or step 4 (taxiway) cannot
   meet its grade limit with the current building elevations, adjust
   the contributing building/terminal pad elevation toward the
   constraint.  Buildings yield, the runway never does.

This reflects the user-confirmed priority and supersedes the earlier
"CIFP-everywhere" formulation in commit 1's STATUS.md, which was
misread by the previous session.

**Implementation status:** the legacy Phase A5 still uses CIFP as a
global anchor for building validation and taxiway elevations.
Migrating away from that is in progress.  Apron triangulation (Phase
C2) is the first step — it now uses DEM + building edges only, no
CIFP queries.  Phase D (junctions) will follow.  The full Phase A5
rewrite is a later commit.

## Why we reverted

The previous session scrapped the legacy and started a new module
(`src/O4_Surface_Patch.py`) on the assumption that the legacy was
beyond repair.  After re-reading the legacy carefully, that decision
was wrong.

**The legacy already implements every active invariant** in a six-phase
pipeline:

| Phase | Lines | Purpose |
|---|---|---|
| **A** | `O4_Auto_Patch.py:2022–2237` | Build all input geometries (boundary, runways, taxiways, aprons, buildings) in meter space |
| **A5** | `2238–2591` | Elevation pipeline: CIFP anchors → building platforms (DEM-sampled, validated against CIFP) → taxiway elevations (grade-constrained against both anchor sources) → CIFP "surface model" (IDW from all anchors) → building back-adjustment for taxiway grade compliance |
| **B** | `2592–2641` | Compute junction zones (taxiway∩runway and taxiway∩taxiway overlaps) — these become triangle wedge candidates |
| **C** | `2643–3425` | Clip and emit non-overlapping shapes.  C1 taxiways clipped against runways+junctions+other-taxiways; C2 aprons (flat or triangulated based on elevation range, with apron-internal grade-constrained Delaunay relaxation); C3 buildings clipped against runways+taxiways+junctions; coverage-fill sweep for any remaining paved area |
| **C4** | `3427–3526` | Detect elevation discontinuities between adjacent flat shapes → mark for transition strips |
| **D** | `3527–3851` | Delaunay-triangulate the junction zone, anchor vertices to nearby flat shapes, grade-constrain via Bellman-style relaxation, emit triangles via `node_altitudes` (legacy already uses this tag for triangles only — exactly what the active invariants now require) |
| **E** | `3854+` | Boundary band + tunnel portals (unfinished, see "Phase E/F status") |
| **F** | `4150+` | Drainage zones (unfinished, see "Phase E/F status") |

The accumulator `all_emitted_parts_m` plus `_emitted_union()` at
`O4_Auto_Patch.py:1973` explicitly tracks every emitted shape so
subsequent phases can subtract it.  That **is** the no-overlap
mechanism.

The new pipeline I wrote during the experiment did NOT implement any
of this — it emitted independent rectangles per feature with no
overlap awareness.  When tested against SPJC it produced 81 619 m² of
cross-feature overlap (taxiway ∩ runway, taxiway ∩ apron, apron ∩
building, etc.) as well as 161 069 m² of intra-taxiway overlap from
adjacent centerline rects sharing corners.  The architecture was
fundamentally wrong for the active invariants.

## Current state (after commit `f6d3140`)

The legacy phase-A-through-F surface generator is the active code
path.  Phases C2 (apron) and D (junctions) are wired to
`O4_Surface_Mesh.adaptive_triangulate`.  `O4_Apt_Dat_Reader` is
integrated at Phase A0.5 and supplies authoritative pavement
geometry; it is also the sole source of truth for runway
FOOTPRINT (lat/lon, width, displaced threshold, blast pad) — CIFP
contributes only threshold elevations.  Phase E (boundary band +
tunnel portals) has been rewritten to emit a segmented sloped
band and FAA-compliant tunnel ramps with U-shaped retaining
walls.  Runway grading uses a wide-window-smoothed DEM profile
with envelope pre-clamp, joint hard-cap / rate-of-change solver,
and the correct FAA runway vertical-curve rule
(305 m per 1 % ΔG).

### Taxiway emission pipeline (Phase C0)

Every apt.dat pavement polygon is classified apron / taxiway by
`O4_Pavement_Classifier`.  Taxiway-classified polygons flow
through a recursive Phase A3 walker (`_recursive_rectify`) that
tries four strategies in order:

1. **Fast path `build_taxiway_rects`:** MRR-aligned rect chain
   when the polygon is a clean strip (aspect ≥ 3.5, width 9–45
   m, aspect-scaled fit ratio).  One rect per straight run.
2. **Voronoi skeleton path (`build_rects_along_centerline`):**
   centerline extraction via `shapely.ops.voronoi_diagram`,
   linemerge, drop < 200 m spurs, RDP 20 m simplification.
   The rect builder walks the centerline at its **own
   RDP-simplified vertices** — NOT uniform sampling — so each
   straight skeleton run emits **exactly one rect**.  Rect
   width = MAX of 5 probe widths + 1.5 m outward padding so
   the rect fully covers the pavement strip; the permitted
   overflow onto neighbouring pavement is subtracted by
   downstream apron/junction emission.
3. **Morphological decomposition (`decompose_multi_taxiway`):**
   `buffer(-15).buffer(+15)` isolates wide junction hubs from
   strip branches; each branch is recursed back through the
   pipeline.  Handles polygons with dramatically wider hubs.
4. **Apron triangulation:** last-resort fall-through via
   `adaptive_triangulate`.

A running **global rect-union** (`global_twy_rect_union`)
deduplicates emissions across every source polygon so sibling
branches at a Y-junction can't overlap each other.

Elevation fidelity: the centerline builder RDPs the
grade-clamped elevation profile with a 3 m tolerance at the
Phase A3 call site — a single straight skeleton run with up to
3 m of DEM bump collapses to one rect, matching the user's
"simplify unless there's more than a meter or so of change"
directive.

### Runway-taxi anchoring

`generate_patch_osm` returns its emitted runway segment chain
alongside the OSM string.  Phase A3 builds a meter-space lookup
(`_runway_elev_lookup`) that projects any (x, y) onto the
nearest runway segment for a byte-identical elevation match.
Centerline endpoints within 25 m of a runway polygon (the
distance allows for Voronoi skeleton leaves sitting ≈ half-
width inside the boundary) are pinned to the runway's own
elevation and propagated through `_clamp_profile_with_anchors`.

apt.dat pavement polygons are simplified to 1 m tolerance at
ingest.

### SPJC invariants (from the current patch)

Zero overlap in the flat category, zero in the triangle
category, **0 m²** intra-sloped (was 1 006 m² before the
vertex-based builder).  **561 m² total overlap (0.02 %)** — all
cross-category, from the intentional MAX-width overflow of taxi
rects into neighbouring apron pavement.  User explicitly allows
this per "shapes may go beyond the apt.dat pavement by a small
amount if needed to simplify coverage."

### SPJC performance

Full pipeline runtime dropped from **45.4 s → 9.6 s** in commit
`f6d3140` via a 4-part profiling pass.  See the "Performance"
section below for the breakdown.

### Modules

| File | Status | Tests |
|---|---|---|
| `src/O4_Auto_Patch.py` | active legacy pipeline + Phase A0.5 + Phase C0 | none yet |
| `src/O4_Surface_Mesh.py` | adaptive triangulation, used by C2 + D | 19 |
| `src/O4_Apt_Dat_Reader.py` | apt.dat parser, wired in Phase A0.5 | 22 |
| `src/O4_Pavement_Classifier.py` | name + shape based apron/taxiway split | 21 |
| `src/O4_Taxiway_Rects.py` | MRR + centerline rect-chain builders; grade clamp with anchor support | 18 |
| `src/O4_Taxiway_Decompose.py` | morphological open-minus decomposition | 8 |
| `src/O4_Taxiway_Skeleton.py` | Voronoi medial-axis centerline extraction | none yet |
| `src/O4_Vector_Map.py` | calls `generate_auto_patches` (unchanged) | none |
| `tests/test_surface_mesh.py` | mesh + grade + anchors | 19 |
| `tests/test_apt_dat_reader.py` | runway/pavement/Bezier/search priority | 22 |
| `tests/test_auto_patch_apt_dat_integration.py` | A0.5 adapter | 11 |
| `tests/test_pavement_classifier.py` | name hints, MRR, classification rules | 21 |
| `tests/test_taxiway_rects.py` | long axis, clamp, RDP, rect chain build | 18 |
| `tests/test_taxiway_decompose.py` | simple strip, T/cross/plus junctions, blob | 8 |
| `tests/fixtures/synthetic_apt.dat` | hand-crafted parser fixture | — |

Total tests: **99** (plus 35 tests from `test_boundary_model.py`
contributed by the parallel boundary-model agent, for 134 collected).
Run with `./venv/bin/python3 -m pytest tests/`.

### SPJC numbers (post commit `f6d3140`)

Combined surface + runway-segment patch (`SPJC_auto.patch.osm`):

| Metric | Commit `693c4f0` | **Commit `f6d3140`** |
|---|---|---|
| Total emitted ways | 2 117 | **1 876** |
| Flat polygons | 824 | **784** |
| Sloped rectangles (runway + boundary + **taxiway**) | 364 | **323** |
| Triangles (apron + transition strips) | 929 | **769** |
| apt.dat taxiway classifications | 9 | **9** |
| apt.dat taxiway rectifiable polygons | 52 | **28** |
| apt.dat taxiway rect segments emitted | 76 | **35** (2 flat, 33 sloped) |
| Rects per rectifiable polygon | 1.46 | **1.25** |
| Cross-category overlap | 0 m² ✓ | **561 m²** (0.02 %) ★ |
| Intra-flat overlap | 0 m² ✓ | **0 m²** ✓ |
| Intra-sloped overlap | 1 006 m² | **0 m²** ✓ |
| Intra-triangle overlap | 0 m² ✓ | **0 m²** ✓ |
| Max runway longitudinal grade | 1.50 % | **1.50 %** ✓ |
| Max taxiway longitudinal grade | 1.50 % | **1.50 %** ✓ |
| Max apron longitudinal grade | 1.00 % | **1.00 %** ✓ |
| Full-pipeline runtime | ≈ 45 s | **9.6 s** |

★ The 561 m² cross-category overlap is intentional — it's the
MAX-width + 1.5 m padding overflow of taxi rects onto adjacent
pavement, explicitly allowed by user direction ("shapes may go
beyond the apt.dat pavement by a small amount if needed to
simplify coverage").  There is NO intra-category overlap; taxi
rects tile cleanly against each other, against runway rects,
and against apron triangles.

### User's shape-count targets (not yet hit)

The user's end-state goal at SPJC is:

* **25-30 taxi rects** (excluding junctions) — currently at **35**,
  close but not quite.  5 extra rects come from elevation-RDP
  splits on longer centerlines where DEM bumps exceed the 3 m
  fidelity tolerance.
* **Junctions** should be flat polys (where elevations allow)
  or triangles (compound slopes) — currently triangulated via
  the apron path.
* **Under 1 000 total shapes** excluding runways — currently at
  **1 876** (1 803 excluding the 73 runway segments).  The
  remaining 1 000 to trim is mostly in:
  - 769 apron triangles (task 2: apron-prefer-flat not yet done)
  - 784 flat polys (313 building pads, ~300 coverage fill, ~130
    boundary band segments, ~30 drainage)

Run `./venv/bin/python3 /tmp/run_legacy.py` then
`./venv/bin/python3 /tmp/audit_legacy.py` to reproduce.  All 100
auto-patch tests pass.

### Performance breakdown

Runtime reductions from the commit `f6d3140` profiling pass:

| Step | Before | After | Fix |
|---|---|---|---|
| `_VertexBag` linear scan | 45.35 s | 32.26 s | Grid-based spatial hash replaces linear scans in `add` / `index_at` / `has_near`.  Adaptive_triangulate was O(n²) in vertices. |
| Phase E1 `unary_union` in loop | 32.26 s | 21.23 s | `subtract_base ∪ e1_emitted_union` was being re-computed 419 times per run even though `subtract_base` is constant.  Replaced with two sequential `difference()` calls. |
| Phase C4 strip subtract | 21.23 s | 16.40 s | Transition-strip inner loop scanned all 800+ emitted flat shapes per pair (260 k `difference()` calls).  Now uses the existing STRtree to query only the 0-20 neighbours that touch the strip's bbox. |
| `_delaunay_clipped` prepared | 16.40 s | **9.60 s** | Prepared geometry for `centroid` containment (5-10× faster than `polygon.contains`), plus a hole-free shortcut that skips the expensive `polygon.intersection(t).area` 99 % check when all 3 triangle vertices are strictly inside and the polygon has no interior rings. |

All four changes are pure optimisation — output is bit-identical
to the `84953a6` pre-speedup patch.  The audit reports the same
counts and the same 561 m² overflow.

### `O4_Apt_Dat_Reader` ready to use

Parses runway, pavement (taxiway/apron/ramp) and boundary geometry
from any X-Plane apt.dat file.  Handles row 1 (header), 100
(runway), 110/111/112/113/114 (pavement with Bezier curves), 130
(boundary).  Public API:

```python
import O4_Apt_Dat_Reader as APR

# Find the most-specific apt.dat.  Searches per-airport Custom Scenery
# packs first, then Custom Scenery/Global Airports, then Resources/
# default scenery.  Returns the file path or None.
path = APR.find_airport_apt_dat("/Users/noah/X-Plane 12", "SPJC")

# Parse one airport block.
apt = APR.load_airport(path, "SPJC")
print(apt.icao, apt.name, apt.reference_elev_m)
print(len(apt.runways), "runways")
print(len(apt.pavements), "pavements")  # each is a shapely Polygon
print(apt.boundary)                      # shapely Polygon or None
```

Verified against the user's `SPJC Lima by Los Flipantes 3.0 Nueva
Terminal XP12` Custom Scenery pack: 51 pavements parsed (810 741
m² total), 2 runways with displaced thresholds, boundary present.

## Plan

Working in commits so each step has a roll-back point.  The plan
shifted between commits 3 and 4: rather than fixing all the legacy
bugs in place and then doing a single big refactor, we now
**refactor as we touch** — every commit that touches a section of
the legacy first extracts the affected helper into its own
single-purpose module with unit tests.  Future fixes land in the
small modules, not the monolith.

1. **Commit 1 — revert flag, restore legacy.** ✅ `ee02cc3`
2. **Commit 2 — A4b building merge fix (no convex_hull).** ✅ `811dda1`
3. **Commit 3 — FAA vertical-curve rule + STATUS corrections.** ✅ `27ad1cd`
4. **Commit 4 — extract `O4_Surface_Mesh.py` + adopt pytest.** ✅ `c68680c`
5. **Commit 5 — wire Phase C2 (apron) to adaptive_triangulate.** ✅ `26b63c1`
6. **Commit 6 — wire Phase D (junction) to adaptive_triangulate.** ✅ `7d89f33`
7. **Commit 7a — Phase D hole-vertex + actual-triangle accumulator.** ✅ `e36e6a1`
8. **Commit 7b — add `O4_Apt_Dat_Reader` module + tests.** ✅ `ce3b97c`
9. **Commit 8 — wire `O4_Apt_Dat_Reader` into legacy (Phase A0.5).** ✅ `2e8ab66`
10. **Commit 9 — pavement-as-taxiway grade, morphological terminal clean, flat buffer ring around buildings, drainage from pavement holes.** ✅ `115be6e`
11. **Commit 10 — zero overlap: precise building cutouts + triangle/ditch containment + Phase E buffered subtraction.** ✅ `6d237d6`
12. **Commit 11 — zero overlap with runway segments: include CIFP runway rect in Phase A1, strict 99 % triangle-in-polygon containment in Surface_Mesh.** ✅ `1bb3e4a`
13. **Commit 12 — generate_patch_osm uses apt.dat row-100 for footprint (lat/lon + width), CIFP only for elevation anchors.** ✅ `d5a0ce4`
14. **Commit 13 — apt.dat is sole source of truth for all runway footprint geometry including displaced thresholds and blast pads; `Runway` dataclass gains `blast_a_m` / `blast_b_m`.** ✅ `f46953f`
15. **Commit 14 — Phase E2 tunnel portals: fix missing SW portal (road-direction probe instead of centroid heuristic), dedupe divided-highway carriageways, portal exclusion disc carved from boundary band.** ✅ `b8a096c`
16. **Commit 15 — Phase E2 tunnel portal rewrite: ramp from OSM tunnel node to airport boundary with inverted slope (depth at portal, surface at node), combined carriageway width, U-shaped retaining wall with 0.5 m gap, four to_ll lat/lon unpacking bugs fixed.** ✅ `148d17b`
17. **Commit 16a — Phase E1 segmented sloped boundary band: walk boundary in 40 m steps, sample CIFP surface model at each endpoint, emit sloped rect per segment with flat corner polys where clipping is needed.** ✅ `b2a931b`
18. **Commit 16b — Runway grading: FAA runway vertical-curve rule (305 m per 1 % ΔG), wide-window DEM smoothing, envelope pre-clamp from all anchors, joint hard-cap + rate-of-change solver, blast-pad boundary condition.** ✅ `53a3d53`
19. **Task 1 — pavement classifier (O4_Pavement_Classifier) with name + shape rules; Phase A3 uses it to split apt.dat pavements into apron / taxiway lists.** ✅ `5e5c3dc`
20. **Task 2 — Phase C0 taxiway rect-chain emission (O4_Taxiway_Rects) with MRR long-axis sampling, combined grade + rate-of-change clamp, and RDP profile simplification.** ✅ `b488350`
21. **Task 3 — apron cleanup: union/difference dedup of overlapping apt.dat pavements (Base Ramp + sub-ramps), MAX_APRON_GRADE = 1.0 %, coverage-fill uses emitted triangles rather than piece outlines, terminal_zone_m removed from "covered" set.** ✅ `239efa9`
22. **Task 4 — apron grade-violation filters: sliver triangle drop (< 5 m², min edge < 2 m) and plane-gradient cap at 2 × MAX_APRON_GRADE.** ✅ `603aa1a`
23. **Progress-line fix — include C0 taxi rects in the final "N taxiway (…)" summary.** ✅ `836d492`
24. **Taxiway decomposition + remove per-piece flat-residue fill — recursive morphological splitting of multi-strip mega-polygons, deletion of the piece-mean elevation bandaid that caused visible multi-metre flat/triangle steps.** ✅ `cd128b1`
25. **Voronoi-skeleton rect chains (O4_Taxiway_Skeleton) — medial-axis extraction via shapely.ops.voronoi_diagram, linemerge to paths, per-path rect-chain emission with local-width sizing.  Global rect-union dedup across all taxi source polygons.  Raised SPJC from 17 → 125 taxi rects.** ✅ `df7aaf6`
26. **Task A — rect consolidation: RDP tolerance 0.3 m → 1.0 m, centerline skeleton simplification 3 m → 5 m, seg length 30 m → 50 m.  Long straight taxiways now collapse to one rect.  SPJC 125 → 59 taxi rects, ratio 1.18 rects per rectifiable polygon.** ✅ `c3e98c1`
27. **Task B initial — within-1-m-of-runway elevation anchoring in build_rects_along_centerline via a runway_anchor callback and `_clamp_profile_with_anchors` helper.** ✅ `fa4bb55`
28. **Task B refined — (a) endpoint-based crossing detection replaces the per-sample distance check, and (b) generate_patch_osm returns its emitted runway-segment chain so the taxi-runway elevation lookup uses the SAME per-segment elevation the runway actually emitted (not raw CIFP linear interpolation).** ✅ `693c4f0`
29. **Revert accidentally-included boundary-model wire-in from commit 28** — three scaffolding blocks (imports + Phase E1 env-load + classifier-dispatch stub) got swept up from an uncommitted parallel session; surgically removed without touching any of the taxi work. ✅ `e72ed0a`
30. **Rect-length cap (100 m) + median local width** as a first attempt to fix rotation misalignment on long curved rects. ✅ `a310e2c` (superseded by the vertex-based builder below)
31. **Vertex-based rect emission + MAX-width-with-overflow + skeleton re-tuning.**  `build_rects_along_centerline` now walks the centerline at its own RDP-simplified VERTICES instead of uniform 50 m sampling — one rect per straight skeleton run.  Width uses MAX of 5 probes + 1.5 m padding (was min → then median), explicitly allowed to overflow the apt.dat polygon per user direction.  Skeleton simplification raised to 20 m and min path length to 200 m.  SPJC: 158 → 35 taxi rects, within the user's 25-30 target.  Intra-sloped overlap 1 370 m² → 0 m². ✅ `84953a6`
32. **Profiling pass: 45.4 s → 9.6 s (4.7× speedup).**  Four bottlenecks addressed without changing any output: (a) `O4_Surface_Mesh._VertexBag` gets a grid-based spatial hash replacing the O(n²) linear scan; (b) Phase E1 boundary band stops re-unioning the static `subtract_base` 419 times per run; (c) Phase C4 transition-strip inner loop uses the existing STRtree to query only neighbouring flat shapes instead of scanning all 800+; (d) `O4_Surface_Mesh._delaunay_clipped` uses `shapely.prepared` geometry for centroid-contains plus a hole-free trivial-contain shortcut.  All 100 tests pass, zero-overlap invariants intact. ✅ `f6d3140` (this commit)

### Commits 11–16 (done): runway footprint, tunnels, boundary band, runway grading

**Commit 11 (`1bb3e4a`)** audited the COMBINED output (surface
patches + runway segments from `generate_patch_osm`) for the
first time and found 2 639 m² of cross-feature overlap that the
surface-only audits had missed.  Two root causes: (1) Phase A1's
`rwy_union_m` contained only the apt.dat runway rectangle, but
`generate_patch_osm` emitted CIFP-threshold-based segments with
a 30 m overrun that poked past the apt.dat rect; (2)
`O4_Surface_Mesh._delaunay_clipped` used centroid-in-polygon
which lets Delaunay triangles bridge concave bays (runway cut
out of apron) and poke into subtracted regions.  Fixes: union
a CIFP-derived rect into `rwy_union_m` as well, and require the
triangle area to be ≥ 99 % inside the polygon rather than just
its centroid.

**Commit 12 (`d5a0ce4`)** made apt.dat row-100 the source of
truth for runway footprint geometry in `generate_patch_osm`.
CIFP still supplies threshold elevations, but runway rectangle
lat/lon and width come from apt.dat.  New `apt_runways` dict
parameter passes `{designator: (lat, lon, width_m)}` pre-loaded
at the call site in `generate_auto_patches`.

**Commit 13 (`f46953f`)** extended this to ALL footprint info:
the `Runway` dataclass gained `blast_a_m` / `blast_b_m` fields
parsed from row-100 field 4 of each end block, and the
`apt_runways` tuple grew to
`(lat, lon, width, displaced_m, blast_m)`.
`generate_patch_osm` now uses apt.dat's displaced-threshold
distances and blast-pad lengths instead of CIFP + the legacy
30 m `OVERRUN_EXTENSION`.

**Commit 14 (`b8a096c`)** fixed Phase E2 tunnel portal
detection at SPJC.  Three bugs: (1) the outward-direction
heuristic used dot(road_dir, centroid-vec), which for
boundary-hugging tunnels is nearly perpendicular and flipped
unpredictably; the SW portal's direction ended up pointing INTO
the airport, so all six grade-down segments were filtered out
and the portal was silently dropped.  Fix: probe 20 m along the
road in each sign and pick whichever lands outside the airport
footprint.  (2) Divided-highway tunnels produced two sets of
overlapping portal rects, one per carriageway — now deduplicated
within 40 m.  (3) The Phase E1 band overlapped the Phase E2
rects at the portal corner; Phase E1 now subtracts a 35 m
exclusion disc around each portal before emitting.

**Commit 15 (`148d17b`)** rewrote Phase E2 entirely to match the
user's mental model of tunnel ramps:
- Ramp spans from the OSM tunnel way's outside endpoint to the
  airport boundary (length = `r_ls.difference(airport_footprint)`
  sub-line length), NOT a fixed 120 m.
- Ramp width combines both divided-highway carriageways into
  one wide rect covering the perpendicular span of all clustered
  portals.
- Slope is flipped: surface DEM elevation at the OSM tunnel
  node (HIGH end), `apt_elev − TUNNEL_DEPTH_DEFAULT` at the
  airport boundary (LOW end).  The ramp descends toward the
  airport rather than away from it.
- Retaining walls form a U-shape around the LOW end (two side
  walls along the ramp plus a cap across the portal) with a
  0.5 m `TUNNEL_WALL_GAP` to the ramp on every edge.  The cap
  length is sized to fit BETWEEN the side walls without
  overlapping them at the corners.
- Also fixed four `to_ll(x, y)` lat/lon unpacking bugs (to_ll
  returns `(lon, lat)`, not `(lat, lon)`) that were producing
  degenerate ~1 m wide ramp rectangles.

**Commit 16a (`b2a931b`)** replaced Phase E1's "one big flat
polygon per band piece" with a segmented sloped band.  The old
approach painted huge areas at a single centroid-sampled
elevation, flattening any sloped airport.  New approach walks
the airport boundary in `BAND_SEG_LENGTH = 40 m` steps; each
step emits one rect `BOUNDARY_BAND_WIDTH` deep (inward from the
boundary), with the two short-edge elevations sampled from
`_cifp_surface_elevation` at the two boundary points.  The rect
is emitted sloped when the endpoints differ, flat otherwise.
Each new rect is clipped against a running subtract set
containing `emitted_union`, `portal_exclusion`, AND all E1 rects
emitted earlier in the walk (so adjacent segments at polygon
corners don't overlap).  When the clipped area is < 99.9 % of
the original, the clipped geometry is emitted as a flat polygon
at the averaged elevation — these are the "small flat polygons
for flat areas" connecting adjacent sloped rects.  SPJC goes
from a handful of ~200 k m² flat pieces to 294 sloped rects +
131 corner flats tracking elevation 5.9 – 45.1 m.

**Commit 16b (`53a3d53`)** rebuilt runway grading in
`generate_patch_osm` to match real FAA-compliant vertical
curves.  Six stages:

1. DEM sampling unchanged (every 100 m along centerline).
2. Wide moving-average smoothing — window = runway length / 4
   on each side — collapses local DEM bumps (hangars,
   surface-model noise, vegetation) into the broad trend.
3. Envelope pre-clamp: for each sample, compute the tightest
   `[lower, upper]` band imposed by every CIFP anchor at
   distance `d` using `|elev − anchor| ≤ d × max_grade`; clamp
   smoothed DEM into that envelope.  Makes the profile
   anchor-consistent in one pass.
4. Joint hard-cap / rate-of-change solver: alternate the two
   until both converge.  The hard cap uses a joint range clamp
   (intersection of both neighbors' max-rise bands, midpoint
   when infeasible).  The rate-of-change pass uses the new
   **`MAX_RUNWAY_GRADE_CHANGE_PER_M = 1/30000`** constant — the
   correct FAA runway rule (305 m per 1 % ΔG) instead of the
   taxiway rule 1/3000 that was being mis-applied.
5. Blast-pad boundary condition: at each end anchor the flat
   blast pad implies `g_left = 0` from outside, so the first
   interior segment is clamped explicitly by
   `max_dg × seg_len`.  Without this, the first runway segment
   spiked straight from 0 % to 1.5 % in 100 m, violating the
   VPI rule by 4.5×.
6. Grade separation.  `MAX_RUNWAY_GRADE_CHANGE_PER_M = 1/30000`
   and `MAX_TAXIWAY_GRADE_CHANGE_PER_M = 1/3000` are now
   separate constants.

SPJC RWY B (16R-34L, 650 m displaced threshold) now emits a
smooth vertical curve: blast pad 0 % → 0.78 % → 1.07 % → 1.17 %
→ 0.88 % → 0.49 % → 0.29 % → ... rising from the flat blast pad
to a peak of ~1.2 % and tapering back down over ~500 m.  Max
absolute grade 1.50 % (at the cap).

### Commits 8–10 (done): apt.dat integration and zero overlap

**Commit 8 (`2e8ab66`)** wired `O4_Apt_Dat_Reader` into the legacy
at a new Phase A0.5.  `xplane_root_from_cifp_path()` derives the
X-Plane root from the CIFP directory; if an apt.dat for the ICAO
is found, `_dico_from_apt_dat()` converts its pavements (tile-
relative), runway rectangle, and boundary into the same
`dico_apt_entry` shape Phases A3/C2/C3 already consume.  All
pavement is merged into the apron field; taxiway is left empty.
When `apt_dat_used` is true Phase A3 skips `APRON_BUFFER` and
preserves interior rings, and Phase C2 skips simplification (which
was erasing small holes).

**Commit 9 (`115be6e`)** adjusted four things:
1. Phase C2 grade for apt.dat pavement is `MAX_TAXIWAY_GRADE`
   (1.5 %) rather than `MAX_APRON_GRADE` — all pavement is
   treated as taxiway.
2. Phase A4b applies morphological opening (buffer -5, +5) to
   merged terminals so < 10 m protrusions disappear, keeping
   clean outlines.
3. `_apron_anchors` adds a `TERMINAL_FLAT_BUFFER_M = 10.0` ring
   around each terminal at building elevation, so the ground
   around a terminal is flat before sloping away to taxiways.
4. Phase F paved-contact filter includes `apron_union_m` so
   grass islands are recognised as drainage candidates when
   the surrounding pavement is apt.dat apron.

**Commit 10 (`6d237d6`)** drove total cross-feature overlap from
5 101 m² to exactly 0 m²:
- Phase A4: `BLDG_PAD = 0.0`, `BLDG_SIMPLIFY = 2.0`.  Buildings
  are cut from the terrain precisely, not padded.
- Phase A4 cleanup pass: sort buildings by area descending,
  subtract each prior building from the current one so no two
  pads overlap.
- Phase C2 apron emission: per-triangle intersection-area check
  against the building union (OVERLAP_EPS = 0.5 m²), catching
  triangles whose edges or area cross a building even when the
  centroid is outside.
- Apron triangles are tracked individually in
  `emitted_apron_triangles_m` and added to `all_emitted_parts_m`
  so downstream phases see the real mesh, not the piece polygons
  with holes.
- Phase F Type B drainage: reject candidates whose ditch
  rectangle escapes the infield polygon by more than 0.5 m².
- Phase E boundary band: buffer `emitted_union` by 1.0 m before
  subtracting, consuming sub-meter slivers around building edges.

Historical notes from the original commit 8 plan below.

### ← Commit 8 (done): wire `O4_Apt_Dat_Reader` into the legacy

This is where the real payoff lives.  apt.dat polygons are the
*authoritative* pavement geometry — they're what X-Plane renders
as the texture, so elevation patches built against them align
perfectly with the ground texture (no visible seams).  Most
importantly, **apt.dat polygons are disjoint by construction**:
no buffering, no precision-drift artefacts, no overlap fixups
needed downstream.

#### What commit 8 does

Add a Phase A0.5 step inside `generate_airport_surface_patches`
in `src/O4_Auto_Patch.py`.  Pseudocode:

```python
# At the top of generate_airport_surface_patches, BEFORE Phase A:
import O4_Apt_Dat_Reader as APR

xplane_root = _xplane_root_from_cifp_path(...)  # or pass in
aptdat_path = APR.find_airport_apt_dat(xplane_root, icao)
apt_data = APR.load_airport(aptdat_path, icao) if aptdat_path else None

if apt_data is not None and apt_data.pavements:
    # Replace dico_apt_entry's OSM-derived pavement shapes with
    # apt.dat-derived ones.  Each Pavement.polygon is in absolute
    # lat/lon and goes into dico_apt_entry in the same shape the
    # legacy expects.
    apt_dat_pavements = _convert_to_dico_format(apt_data, tile)
    dico_apt_entry["taxiway"] = apt_dat_pavements["taxiway"]
    dico_apt_entry["apron"]   = apt_dat_pavements["apron"]
    # boundary may or may not replace dico_apt_entry["boundary"] —
    # OSM boundary is sometimes more accurate, decide per case
else:
    pass  # fall through to OSM data — legacy path
```

The downstream phases C1/C2/C3/D all read from `dico_apt_entry` and
work unchanged on the new inputs.  No phase rewrites needed.

#### What integration needs

1. **Find xplane_root.**  The legacy has `find_aptdat(cifp_path)`
   which traces the xplane root from the CIFP directory.  Reuse
   that logic — derive `xplane_root = os.path.dirname(os.path.dirname(cifp_path))`.

2. **Distinguish taxiway vs apron in apt.dat output.**  apt.dat
   row 110 doesn't actually distinguish taxiway from apron — it
   just has a "name" field that may say things like `"TWY A"`,
   `"RAMP 1"`, `"Aeronaval"`, `"AV"`.  Two reasonable approaches:
   (a) Concatenate all pavements into a single `dico_apt_entry["apron"]`
       MultiPolygon and leave `dico_apt_entry["taxiway"]` empty.
       Phase C2 then triangulates everything as "apron-style"
       (1.0 % grade).  Simple but loses the looser 1.5 % grade
       that taxiways get.
   (b) Use the pavement name to classify: `"TWY"` / `"TAXIWAY"` /
       `"taxi"` / etc. → taxiway; everything else → apron.  More
       work, but preserves the per-feature grade rules.

   Recommendation: **start with (a) and refine later**.  At SPJC
   the difference is small because most pavements are aprons.

3. **Boundary handling.**  apt.dat row 130 is the airport perimeter.
   OSM has a `dico_apt_entry["boundary"]` with the same purpose.
   They usually agree but sometimes disagree.  For commit 8, prefer
   apt.dat when available, fall back to OSM.

4. **Extraction → meter space.**  apt.dat polygons are in absolute
   lat/lon; the existing pipeline expects tile-relative shapes for
   `dico_apt_entry`.  Build a small adapter that subtracts
   `(tile.lon, tile.lat)` from each vertex.

5. **Buildings stay OSM.**  apt.dat doesn't have detailed building
   outlines.  `dico_apt_entry["hangar"]` and the OSM building data
   stream are unchanged.

6. **Fallback path.**  If apt.dat is missing or has no row-110
   records (small private airfields), fall back silently to the
   OSM path.  No regression.

#### Expected SPJC outcome

With apt.dat polygons replacing OSM-derived shapes:

* `flat ∩ flat overlap` should drop to near zero (apt.dat polygons
  are disjoint).
* `flat ∩ triangle overlap` should drop dramatically (junction
  zones become explicit edges between adjacent apt.dat polygons,
  not computed buffer overlaps).
* Total ways could drop further (no more taxiway-buffer fan
  overlaps, no synthetic taxiway flats from polylines).
* Pavement boundaries align exactly with the rendered texture.

Estimated effort: **half a day**.  Implementation is mostly the
pavement-name classification and the lat/lon → tile-relative
adapter.  Unit tests for the adapter live in a new
`tests/test_apt_dat_integration.py`.

### Next task queue

User's four-item list from the session that produced
`693c4f0`, in priority order:

1. **Rect consolidation** ✅ *done in task A (commit `c3e98c1`)* —
   RDP tolerance raised to 1 m; long straight taxi stretches
   between junctions now emit as one rect.
2. **Apron prefers flat; slope minimally to meet taxiways.**
   Not yet implemented.  The plan (from user direction):
   - For each apron component, first try a single flat poly
     at the dominant DEM elevation.
   - If taxiway-touch points at the component's boundary
     exceed the grade budget from that flat value, split the
     component into the minimum number of sloped rects /
     triangles that respect the 1.0 % rule AND match the
     taxi-entry anchor elevations.
   - Use the already-emitted taxi rects as fixed-point
     anchors (they're in `emitted_twy_quads_m` by the time
     Phase C2 runs).
3. **Taxi-runway smooth transition** ✅ *mostly done in tasks
   B initial + refined (commits `fa4bb55`, `693c4f0`).*
   Centerline endpoints within 25 m of a runway polygon are
   pinned to the emitted runway chain's interpolated
   elevation.  2 of 3 close taxi-runway pairs at SPJC match
   exactly; **1 of 3 has a 0.7 m residual** on taxi -10025
   which is rect[1] of a multi-rect chain where the clamp at
   rect[0] didn't fully propagate.  To finish:
   - Either (a) every rect-corner near a runway gets its own
     anchor check (not just centerline endpoints), or
     (b) shorten the chain around anchored endpoints so the
     anchor dominates the nearest rect's corners directly.
4. **Building-pad reconciliation for smooth apron joins.**
   Not yet started.  Depends on task 2 (apron path) being
   done so building-pad edges can be reconciled against the
   apron grade budget.  Per the STATUS.md elevation build
   order: "Reconciliation: if step 3 (apron) or step 4
   (taxiway) cannot meet its grade limit with the current
   building elevations, adjust the contributing
   building/terminal pad elevation toward the constraint."

Other pending refinements not in the user's top-4:

5. **Intra-sloped taxi overlap (1 006 m² at SPJC).**  Two
   adjacent rects in a single centerline chain overlap at a
   curved join by a few m² each.  Fix candidates:
   (a) clip each rect's tracking polygon to its source branch
   before emission, (b) emit a single longer-than-2-rect chain
   segment at curved joins and forgo the corner overlap, or
   (c) drop the clamp on adjacent-rect perpendicular widths
   so the transition is always angular, not parallelogrammic.
6. **Drainage (Phase F refinement).**  Unchanged from the
   pre-task-1 queue — see "Phase F status" section below.
7. **Helper extraction: `O4_Apron_Mesh.py`, `O4_Drainage.py`.**
   Taxi work has already extracted `O4_Pavement_Classifier`,
   `O4_Taxiway_Rects`, `O4_Taxiway_Decompose`, and
   `O4_Taxiway_Skeleton`.  Aprons and drainage still live in
   the `O4_Auto_Patch` monolith.

## How to resume next session

```bash
cd /Users/noah/Ortho4XP-shred86
git log --oneline -15            # confirm we're at f6d3140
./venv/bin/python3 -m pytest tests/   # 100 local tests pass

# Full pipeline sanity run + audit (driver at /tmp/run_legacy.py):
./venv/bin/python3 /tmp/run_legacy.py      # writes /tmp/SPJC_legacy.patch.osm (~9.6 s)
./venv/bin/python3 /tmp/audit_legacy.py    # see the SPJC numbers table above
```

Expected output: 35 taxi rects emitted in ~9.6 s, zero intra-
category overlap, 561 m² cross-category overlap (intentional).
If either number has changed, something in the taxi/apron/runway
path regressed.

### What the user said to do next

Exact user quote at end of session: "Update the status.md file
and prepare for context clear.  Then we will proceed with fixing
more issues with taxiways."

So the next session is going to continue **fixing taxiway
issues** — specifics not yet stated.  The current taxi output at
SPJC is:

* 35 rects classified and emitted
* Alignment with source polygon: ~12 % of rects have > 5°
  rotation error (the Voronoi skeleton centerline is straighter
  than the actual polygon curve on some branches)
* 1 taxi-runway join at ~0.7 m residual mismatch (the rest are
  exact via the unified runway chain lookup)
* Junction regions: currently triangulated by the apron path,
  NOT emitted as the user-preferred "flat poly when elevations
  allow, triangle for compound slopes"

Likely candidates for the next taxi fix pass (user will
confirm):

1. **Remaining rotation misalignment** — for polygons that go
   through the Voronoi skeleton path, each straight skeleton
   segment is already one rect, but the skeleton itself may
   not follow the true local tangent at curves.  Possible
   fix: walk the polygon boundary instead of the skeleton.
2. **Junction polys as flat-or-triangle instead of apron tri
   mesh** — per user: "Junctions should be flat polygons where
   possible, or triangles where a compound slope is needed."
3. **Residual 0.7 m mismatch** at one taxi-runway join.
4. **The 5 extra rects over the 25-30 target** — either raise
   elevation fidelity tolerance further or accept the current
   count as good enough.

### Where the relevant code lives

- **`src/O4_Pavement_Classifier.py`** — name + shape rules.
- **`src/O4_Taxiway_Rects.py`** — both rect builders
  (`build_taxiway_rects` MRR fast path, `build_rects_along_
  centerline` skeleton path), grade clamp variants.
  Key tunables at top of file.
- **`src/O4_Taxiway_Skeleton.py`** — Voronoi medial-axis
  extraction.  `DEFAULT_MIN_PATH_LENGTH_M = 200`,
  `DEFAULT_SIMPLIFY_TOL_M = 20`.
- **`src/O4_Taxiway_Decompose.py`** — morphological opening
  for multi-strip hub splitting (used as fallback after
  skeleton fails).
- **`src/O4_Auto_Patch.py`** —
  - Phase A1 runway union: ~line 2654 (`rwy_union_raw_m` saved
    before the safety inflate)
  - Phase A3 taxi classification + recursive rectify:
    ~line 2790-3180 (`_recursive_rectify`, `_try_skeleton`,
    `_try_rectify`, the `global_twy_rect_union` dedup)
  - Phase C0 emit loop: ~line 3200-3310 (uses cached
    `apt_twy_rect_chains`)
  - Phase C2 apron emission: ~line 4010 (where task 2 will
    live)
  - `_runway_elev_lookup` closure: ~line 2880-2960 (reads from
    `runway_segment_chain` returned by `generate_patch_osm`)

### Still-pending tasks from earlier session (not next priority)

The four-item task list from before the profiling session is
still valid — the user set it aside temporarily to chase visible
issues and performance.

1. ✅ Rect consolidation (task A, commit `c3e98c1` + superseded
   by commit `84953a6`).
2. ❌ Apron prefers flat over triangulation (task 2 of the
   earlier list).  Biggest remaining shape-count lever.
3. ✅ mostly — taxi-runway smooth transition (2/3 exact, 1
   residual 0.7 m).
4. ❌ Building-pad reconciliation for smooth apron joins.

### Parallel work on boundary modelling

A parallel agent is developing `O4_Boundary_Model.py` on branch
`smart_airport_boundary` off commit `e72ed0a`.  That module is
NOT imported by any code on the `dev` branch — commit `e72ed0a`
explicitly backed out the accidental import.  If you see
`O4_Boundary_Model.py` or `test_boundary_model.py` appear on
your working tree, leave them alone — they belong to the other
agent.

## Bugs to fix in the legacy

These are specific, targeted fixes — not architectural changes.  Each
should be a small commit.

### 1. A4b building merge uses `convex_hull` (snowballs to giant blobs)
`O4_Auto_Patch.py:2221` runs `terminal.convex_hull` after each
absorption, and `BLDG_MERGE_DIST` is measured from the growing hull.
Snowballs into 394 000 m² polygons at SPJC that swallow the airfield.
**Fix:** measure radius from the *original* terminal footprint
(buffer once at the start), single-pass merge, replace `convex_hull`
with `unary_union().simplify(BLDG_SIMPLIFY)` so the terminal keeps its
real concavities.

### ~~2. `generate_patch_osm()` emits one rect per runway~~ (already fixed)
**Verified incorrect.**  The legacy `generate_patch_osm()` already
samples DEM at `RUNWAY_SEGMENT_LENGTH` (100 m) intervals, anchors
CIFP elevations at displaced thresholds + physical ends, applies
per-segment grade relaxation, and emits per-segment rectangles plus
flat overruns.  At SPJC it produces 73 segmented runway ways (66
sloped, 7 flat overruns) with elevations covering 5.9–32.2 m and
49 unique values.  The previous-session note claiming "one rect per
runway" was looking at an even-older version, or confused this
function with something else.  No fix needed.

### 3. ~~No FAA vertical-curve rule~~ (added in commit 3)
The legacy didn't apply `MAX_GRADE_CHANGE_PER_M = 1/3000` between
adjacent runway segments.  At the legacy's 100 m segment length and
1.5 % grade cap the rule is auto-satisfied, so this pass is mostly
defensive — but it's now in place so any future tightening of segment
length or grade cap stays compliant.  **Done in commit 3.**

### 4. IDW elevation falloff returns 0 m beyond radius
The legacy `_idw_elev` (or its inline equivalent in the CIFP surface
model) zeros out vertices beyond `IDW_RADIUS`, causing runaway-low
drainage shapes.  **Fix:** add a nearest-anchor fallback that returns
the closest anchor's elevation when no anchor falls within radius.

### ~~5. Displaced threshold positions not anchored as samples~~ (already fixed)
**Verified incorrect.**  `generate_patch_osm()` already inserts the
displaced-threshold fraction into the sample chain (lines ~558-566)
and seeds it with the CIFP elevation in the anchor branch (lines
~575-601).  No fix needed.

### 6. Plane fit / centroid math uses absolute meter coordinates
At tropical latitudes `|x| ≈ 8 × 10⁶`, which blows the condition
number of any normal-equation solve.  **Fix:** replace the projector
with an airport-centered version (`(lon-mid_lon, lat-mid_lat)` × scale)
so meter-space coordinates stay in the ±10 km range.  This was found
during the experiment and is the right fix in the legacy too.

## Phase E status (rewritten in commits 14–16a, current)

- **E1 boundary band**: segmented sloped rects along the boundary
  (commit 16a).  Elevation per rect is sampled at both short-edge
  endpoints from the CIFP surface model (IDW over all
  runway / building / taxiway / road anchors).  Flat polys fill
  the places where clipping against emitted shapes leaves a
  non-rectangular residue.  Tracks the airport slope properly
  at SPJC (5.9 – 45.1 m band elevation range).
- **E2 tunnel portals**: full rewrite in commits 14 and 15.
  Per-portal ramp from the OSM tunnel node to the airport
  boundary with inverted slope (depth at boundary, surface at
  node), combined carriageway width, U-shaped retaining walls
  with a 0.5 m gap.  Probe-based outward-direction detection
  handles the SW SPJC portal that the legacy centroid heuristic
  dropped.

Phase E is considered **functionally complete** for SPJC.
Remaining refinements are grade-precision tuning if/when the
visual output needs it.

## Phase F status (drainage zones, refinement pending)

Partially implemented.  Type A (flat infield depression) and
Type B (sloped drainage ditch) both emit shapes, and the
containment check on Type B ditches (`ditch_poly.difference(ip)
< 0.5 m²`) prevents the old class of overlap with surrounding
pavement.  Known issues:

- Type A flat infield elevations are sometimes more than a real
  drainage ditch depth below the surrounding pavement.
- Type A / Type B selection logic is driven by the elevation
  variance of the surrounding pavement samples and mis-classifies
  some flat islands as sloped-ditch candidates.
- The per-pavement-hole classification lumps multi-feature
  drainage basins together instead of producing one ditch per
  natural sub-basin.

Drainage refinement is the **third item** in the next-task
queue below.

## Recently fixed regressions — DO NOT REVERT

These bugs were found earlier in the iteration history; if any
reappear it means an old approach was copied back in.

### Duplicate runway geometry
`generate_auto_patches()` used to call `generate_patch_osm()`
unconditionally *and* call `generate_airport_surface_patches()`.  The
new revert path now does this on purpose — the legacy intentionally
uses `generate_patch_osm()` as the runway baseline and the surface
generator layers other features on top WITHOUT re-emitting runways.
Make sure that distinction stays clear when refactoring.

### Apron subtraction mangling building edges
An earlier attempt subtracted the apron geometry from each building
to avoid apron-on-building overlap.  This produced ugly notched
building polygons.  Fix is to subtract buildings from the apron at
intake (already done in legacy Phase C2), so the apron never contains
building footprints; buildings are emitted unmodified.
**Do not re-introduce `building.difference(apron_union)` anywhere.**

### Legacy public API must be preserved
`generate_airport_surface_patches(icao, taxiway_data, building_data,
runway_pairs, tile, dico_apt_entry, start_node_id=-10000,
road_data=None)` is called from `O4_Vector_Map.py` and from
`generate_auto_patches()` in the same module.  Any refactor must
preserve this signature (or update both call sites in the same
commit).

## Test airport

All debugging is against SPJC (Jorge Chávez, Lima, Peru).  Inputs
documented in `docs/TEST_PLAN_SPJC.md`.

- CIFP: `/Users/noah/X-Plane 12/Custom Data/CIFP/SPJC.dat`
- apt.dat: `/Users/noah/X-Plane 12/Global Scenery/Global Airports/Earth nav data/apt.dat`
- OSM cache: `OSM_data/-20-080/-13-078/-13-078_airports.osm.bz2`
- Building cache: `OSM_data/-20-080/-13-078/-13-078_apt_bldg_local.osm.bz2`
- DEM: `Elevation_data/-20-080/S13W078.hgt`
- Reference output: `Patches/-20-080/-13-078/SPJC_auto.patch.osm`
