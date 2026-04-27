# Per-shape elevation-field plan

> **Status — 2026-04-26 (post-implementation):** the plan was implemented
> in `O4_Airport_Pavement_Builder.py` (`_smooth_polygon_grid` helper,
> Mode A/B constants, optional Mode B integration in
> `_triangulate_junctions`, gated on `USE_PER_POLYGON_ELEVATION_FIELD`).
> Mode B is currently DISABLED by default because measurement showed
> it does not deliver the promised reduction at the current codebase
> state.  See "Lessons learned" at the bottom for the diagnosis.

Implementation plan for eliminating within-shape grade violations at
SPJC / CYXY / HECA / KBNA (the test airports) and at any airport
that has comparable geometry.  Replaces the current per-polygon
elevation-derivation code with two new smoothing modes that share a
common goal: **every pavement vertex's elevation must come from a
grade-compliant field, and every shared vertex between adjacent
shapes must carry the same elevation in both shapes.**

This is a NEW implementation chunk to start with a fresh context.
Read this document first, then proceed to the implementation steps
at the bottom.

---

## Why we're doing this

The current pipeline produces thousands of within-shape grade
violations across the four test airports because:

* Junction polygons inherit elevations from independently-anchored
  rect/runway/terminal corners that aren't mutually grade-compliant.
* The graph-based elevation network operates on centerlines + bridge
  edges only; the 2D interior of a junction polygon is sampled too
  sparsely to guarantee compliance between any two boundary vertices.
* Subdivision and free-vertex clamping (already landed in earlier
  passes) treat symptoms; the source is that the elevation FIELD
  itself isn't grade-compliant in 2D within paved regions.

The fix is to compute a grade-compliant elevation field PER PAVEMENT
SHAPE.  Rect shapes (taxis, runways) use a 1D field along their
axis; non-rect pavement shapes (junctions, aprons, terminals) use a
2D grid covering only the polygon interior.  Unpaved areas don't
exist in any field — X-Plane renders DEM there, as it does today.

---

## Architecture

### Two smoothing modes by shape type

**Mode A — rects (runway segments, taxi rects): 1D axis smoothing**

For each rect shape:

1. Sample the elevation graph at 5 m intervals along the rect's
   axis (start, ..., end).
2. Anchor the two endpoint samples to authoritative values:
   - Runway segments: CIFP threshold elevations propagated by the
     existing runway-segmenting code in `O4_Auto_Patch.generate_patch_osm`.
   - Taxi rects: graph value at the rect's axis endpoint, using
     the existing `_compute_elevations` taxi-rect logic.
   - Shared corners with other rects/junctions: set in earlier
     passes (Mode A always runs before Mode B).
3. Smooth interior samples: Laplacian + grade cap
   (`|Δelev| ≤ 5 m × TAXI_MAX_GRADE = 0.075 m` per 5 m step).
4. Iterate to fixed point (cap 50 iterations or until max change
   < 0.01 m).
5. Emit the rect:
   - Range < 0.1 m → flat `altitude` tag (single value).
   - Else → `altitude_high` / `altitude_low` derived from smoothed
     endpoint samples (the rect is sloped along its axis,
     flat across width by construction).

This is essentially what `O4_Auto_Patch.generate_patch_osm` does for
runway segments today; we extend the pattern to taxi rects.

**Mode B — non-rect pavement (junctions, aprons, terminals): 2D 5 m grid smoothing within the polygon**

For each non-rect pavement polygon:

1. Build a **5 m grid** covering the polygon's bbox.
2. Mark each cell as INSIDE if its centre lies inside the polygon
   exterior, else INACTIVE.  Cells with any portion inside the
   polygon are also INSIDE.
3. For each INSIDE cell, find any HARD anchors that constrain it:
   - Boundary vertices shared with rect/runway/terminal corners
     (Mode A has already chosen those elevations; they're immutable
     during Mode B).
   - Boundary vertices shared with another junction (handled by the
     cross-junction iteration described below).
4. Pin cells nearest each hard anchor to the anchor's elevation.
5. Initialize non-anchor cells: DEM value at the cell centre,
   clipped into the local feasibility cone induced by every
   reachable anchor.  Same rule the elevation graph's
   `choose_values` uses today.
6. Smooth INSIDE cells: Laplacian over neighbour cells + grade-cap
   pass between adjacent cells (`|Δelev| ≤ 0.075 m` per 5 m step).
7. Iterate to fixed point (cap 50 iterations or until max change
   < 0.01 m).
8. Sample boundary vertex elevations bilinearly from the smoothed
   grid.

If the polygon has NO reachable anchor (no shared vertices with any
rect/runway/terminal/junction), use the global elevation graph's
`elevation_at(x, y)` directly for each boundary vertex — same path
the existing `_vertex_elev_anchored` step 5 fallback uses today.

After sampling, emit the polygon using the existing FLAT / PLANAR /
COMPOUND classifier (with the slope-magnitude check we added):

* **FLAT** (vertex range < 0.1 m): single `altitude` tag.
* **PLANAR-COMPATIBLE** (plane-fit residuals < 0.30 m AND plane
  slope ≤ TAXI_MAX_GRADE): emit with linearly-interpolated
  `node_altitudes` (X-Plane interpolates linearly across the ring;
  identical render to a triangulated planar mesh).
* **COMPOUND** (multi-direction slope, residuals exceed planar
  tolerance): emit with per-vertex `node_altitudes`.  Triangle4XP
  triangulates and bilinearly interpolates Steiners from the
  boundary samples.

**Mode C — unpaved areas**: not in any field.  X-Plane renders DEM.
Nothing for us to do.

### Cross-shape consistency at shared vertices

Two-pass strategy with iteration:

1. **Pass 1 — Mode A on every rect**.  After this, every rect's
   altitude_high/low (or flat altitude) is final and every rect
   corner has an authoritative elevation.

2. **Pass 2 — Mode B on every junction, iterated**:
   - 2a: each junction's grid is smoothed using ONLY its
     rect-corner shared vertices as hard anchors.  Boundary
     vertices shared with OTHER junctions are initially free.
   - 2b: walk each junction's free-but-shared-with-neighbour
     boundary vertex; set its elevation to the AVERAGE of what the
     two junctions' grids assigned to it.
   - 2c: re-pin those averaged values as hard anchors and re-smooth
     each junction's grid.
   - Iterate 2b/2c **up to 8 times** (per user 2026-04-26).  Early
     exit when no boundary vertex moves > 0.05 m between iterations.
     If 8 iterations doesn't converge, emit a WARN with the airport
     ICAO + count of unconverged vertices and proceed.

### Why this satisfies the within-shape rule

* Mode A guarantees each rect's two short-edge elevations are
  within `axis_length × TAXI_MAX_GRADE` of each other.  Cross-axis:
  the rect is flat across width by construction — cross-pairs are
  0 % grade.
* Mode B guarantees each cell-pair within the polygon's grid
  satisfies `|Δelev| / 5 m ≤ TAXI_MAX_GRADE = 1.5 %`.  Boundary
  vertices ≤ 5 m apart inherit this directly.  Boundary vertices
  > 5 m apart inherit it transitively through chained cells:
  for any pair, `|Δelev| / dist ≤ TAXI_MAX_GRADE` because the
  sum-of-cell-deltas along any path between them is bounded by
  the path length × the per-step grade cap.
* Cross-shape: shared vertices have the SAME elevation in both
  shapes by anchoring; no step at boundaries.

---

## User's confirmed answers (2026-04-26)

1. **Initial values for non-anchor cells in Mode B**: DEM clipped
   into the local feasibility cone (preserves the "real airports
   excavate / build up" rule).

2. **Polygons with NO reachable anchor**: adopt the global
   elevation graph's value via `elevation_at(x, y)`.

3. **Grid spacing**: 5 m (not 30 m).  Smoothing handles curvature
   limit compliance; 5 m is fine enough that even small junctions
   have many cells and boundary sampling is accurate.  Per-step
   grade cap is `5 × 0.015 = 0.075 m`.

4. **Iteration caps**: 8 iterations max for cross-junction
   consistency, with a WARN if not converged.

---

## Implementation steps

### Step 0 — preconditions

Read these existing pieces of the pipeline before implementing:

* `O4_Airport_Pavement_Builder.py`:
  - `_compute_elevations` — the entry point that's about to be
    rewritten.  Currently does graph build → propagate → smooth
    → plateau snap → bridge enforcement → rect altitude
    derivation → junction triangulation.
  - `_triangulate_junctions` — sets junction vertex elevations
    today.  Will be replaced/refactored to read from Mode B.
  - `ElevationGraph` class — produces the centerline elevation
    graph.  KEEP — it's the source of anchor elevations for both
    modes.
  - `_clamp_junction_free_vertices` — Layer 2 clamping.  KEEP as
    a fallback safety net.
  - `_subdivide_violating_junctions` — Stage 1 subdivision.  KEEP
    as a fallback safety net (will rarely fire if Modes A+B work).

* `O4_Auto_Patch.generate_patch_osm` — runway segment elevation
  smoothing (the prior art for Mode A).

### Step 1 — extend Mode A to taxi rects

Mode A for runway segments already exists (CIFP-anchored
1D smoothing).  Extend it to taxi rects:

1. Add a `_smooth_rect_axis(rect, anchors, graph, dem)` helper.
2. For each emitted taxi rect, sample its axis at 5 m intervals,
   anchor the endpoints (graph values), smooth, and write back to
   `BuiltShape.altitude_high` / `altitude_low` (or `altitude` if
   flat).
3. Update the corner-elevation map (`corner_elev`) so adjacent
   junctions in Mode B can read the rect's authoritative corner
   values.

### Step 2 — implement Mode B (per-polygon 2D grid smoothing)

New module-level helper:

```python
def _smooth_polygon_grid(
    polygon: Polygon,
    hard_anchors: List[Tuple[float, float, float]],  # (x, y, elev)
    graph: Optional[ElevationGraph],
    dem,
    tile_lat: int, tile_lon: int,
    layout_anchor: Tuple[float, float],
    grid_step_m: float = 5.0,
    max_iters: int = 50,
    convergence_tol_m: float = 0.01,
) -> Callable[[float, float], float]:
    """Smooth a 2D elevation field within `polygon` using `hard_anchors`
    as fixed points.  Returns a sampler `f(x, y) -> elev` that
    bilinearly interpolates the smoothed grid."""
```

Implementation: build a numpy 2D array indexed by `(x_cell, y_cell)`,
mask cells inside the polygon, pin anchor cells, run Laplacian +
grade-cap iterations.  Return a sampler closure.

Use this sampler to set every junction polygon's boundary vertex
elevations.  Replace the existing `_vertex_elev_anchored` calls
in `_triangulate_junctions`.

### Step 3 — cross-junction iteration

Wrap step 2 in an outer iteration:

```python
for _ in range(8):
    changes = 0
    for junction in junctions:
        sampler = _smooth_polygon_grid(junction.polygon, ...)
        for v in junction.boundary_vertices:
            new_elev = sampler(v.x, v.y)
            if shared_with_neighbour(v):
                new_elev = average_with_neighbour(v)
            if abs(new_elev - v.elev) > 0.05:
                changes += 1
                v.elev = new_elev
    if changes == 0:
        break
else:
    # 8 iterations didn't converge — log WARN
    sys.stderr.write(f"[pav-builder] WARN: {icao}: cross-junction "
                     f"elevation iteration didn't converge in 8 passes; "
                     f"{n_unconverged} vertices unstable.\n")
```

### Step 4 — keep existing classifiers

The FLAT / PLANAR / COMPOUND classifier stays as-is (with the
slope-magnitude check we added).  It now reads from Mode A/B
results instead of per-shape graph sampling.

### Step 5 — keep existing safety nets

The post-plateau bridge enforcement, Layer 2 free-vertex clamp,
junction subdivision, and `to_osm` sliver-corner safety net all
STAY.  They become safety nets for residual cases where Modes A/B
can't fully converge.  Should mostly be no-ops once the field model
works.

### Step 6 — measure

Re-run all four test airports.  Target metrics (post-implementation):

| airport | current within-shape | target |
|---|---|---|
| SPJC | 2509 | < 100 |
| HECA | 11085 | < 500 |
| CYXY | 2083 | < 100 |
| KBNA | 19037 | < 500 |

Worst-case grade should be < 5 % at all airports.  If we hit
significant residual violations, they'll be in narrow-waist
polygons where the grid smoothing converges but the sampled
boundary still has neighbour-pair conflicts; subdivision is the
follow-up.

---

## Existing code paths that must continue to work

* CIFP runway threshold elevations are immutable.
* Apt.dat ⊕ DSF ⊕ OSM-buffer pavement union construction.
* Per-airport apt.dat fallback (KBNA's row-110-less custom pack).
* DEM auto-download via `O4_DEM_Utils.DEM`.
* OSM auto-download for missing airport tiles.
* Apron-merged runway segment dropping.
* Sliver corner detection in `to_osm` (drops degenerate polygons
  that crash X-Plane).

---

## Project rules (from memory)

* SPJC, CYXY, KBNA, HECA are TEST CASES, not targets.  Every fix
  must be airport-agnostic.  No hardcoded ref letters, no ICAO
  branches, no per-airport tuning.
* Every magic constant must be validated against all four test
  airports before commit.
* Never flatten terrain; the goal is grade-compliant pavement
  built INTO uneven terrain.
* CIFP runway elevations are the only HARD truth.
* Runway/taxi grade rule: ≤ 1.5 % per 30 m horizontal (FAA).
  Already encoded as `TAXI_MAX_GRADE = 0.015`.

---

## Constants to add

```python
# In O4_Airport_Pavement_Builder.py
ELEVATION_GRID_STEP_M = 5.0          # Mode B grid spacing
ELEVATION_SMOOTH_MAX_ITERS = 50      # per-grid smoothing cap
ELEVATION_SMOOTH_CONVERGE_M = 0.01   # per-cell convergence threshold
CROSS_JUNCTION_MAX_ITERS = 8         # cross-junction iteration cap
CROSS_JUNCTION_CONVERGE_M = 0.05     # boundary-vertex convergence
# TAXI_MAX_GRADE (= 0.015) already exists
```

---

## Out of scope for this iteration

* Stage 3 anchor-cone clipping (would replace the entire elevation
  graph step — much bigger refactor; only worth doing if Modes A/B
  don't get us close to zero violations).
* Real-time terrain editing (X-Plane handles DEM rendering).
* Building / hangar elevation handling beyond the existing
  terminal-pad logic.
* Runway curvature rule (1 % / 305 m vertical curve) — current
  code only enforces the longitudinal grade cap.  May matter at
  airports with very long runways but not expected at our four
  test cases.

---

## Update — global pavement mesh (2026-04-26 evening)

After the lessons-learned analysis below, the user pointed out the
critical principle: **only runways are HARD anchors.**  Taxi rect
altitudes and terminal pads must be free to be modified to satisfy
grade rules.  A subsequent violation classification pass confirmed
that ~95 % of within-shape violations were SOFT/SOFT pairs (graph-
derived junction boundary vertices), not rect-corner conflicts —
the locked-anchor design was over-constraining the pavement.

### Solution: global pavement mesh

A new `_build_pavement_mesh` builds one mesh covering every
pavement shape's boundary:

* **Nodes** — every unique boundary vertex (deduplicated by 0.5 m
  bucket).  Cross-shape vertex sharing is automatic via the same
  bucket interner that `to_osm` uses.
* **Edges** — three families:
  1. *Ring adjacency*: each ring step `(i → i+1)` gets an edge
     weighted by 2D Euclidean distance.
  2. *Within-polygon all-pairs* within
     `WITHIN_SHAPE_VIOLATION_RADIUS_M` (60 m) — gives the grade-
     cap pass a direct edge for every pair check_grade flags as
     a within-shape violation candidate, so the mesh enforces the
     2D rule directly rather than relying on it falling out of
     ring-adjacency smoothing.
  3. *Cross-polygon proximity* within 5 m — connects mesh nodes
     in different polygons whose buckets are close in 2D, so the
     CROSS-SHAPE rule is also enforced by the same smoother.
* **Anchors** — runway corner buckets ONLY.  Every other node is
  free.

`_solve_pavement_mesh` runs an iterated `propagate_bounds` +
Laplacian + grade-cap pass, interleaved with two equality
constraints applied every iteration:

* **Sloped rect short-edge pairs**: corners (0, 3) share
  `altitude_high`; corners (1, 2) share `altitude_low`.  Both
  short edges can be raised/lowered independently — what's
  forbidden is moving a single corner.
* **Flat-shape corner groups** (terminals, flat rects): all
  corners share a single `altitude`.

Equalising during smoothing (rather than just at the end) lets
the rest of the mesh adjust around the constrained values, so
junctions sharing a rect or terminal corner converge to elevations
consistent with the constrained shape's emit.

`_writeback_pavement_mesh` then rewrites each shape's altitude
tags from the smoothed mesh:

* Runway: untouched (anchors).
* Sloped taxi rect: `altitude_high` from corners 0/3, `altitude_low`
  from corners 1/2 (already equal post-equalisation).
* Flat shape: `altitude` from corner mean.
* Junction: per-vertex `node_altitudes`.

### Measured impact

| airport | check_grade WITHIN-SHAPE / CROSS-SHAPE | | |
|---|---|---|---|
|         | original | after-soft-classification | global mesh |
| SPJC | 2509 / 5  | 2030 / 3  | **372 / 0**  |
| CYXY | 2083 / 1  | 2181 / 1  | **292 / 0**  |
| HECA | 11085 / 47 | 10021 / 56 | **1491 / 0** |
| KBNA | 19037 / 43 | 13768 / 75 | **1304 / 0** |

Within-shape: −85 % to −93 % across the four airports.
Cross-shape: 0 violations everywhere.  VERTEX-TO-EDGE: also down
from dozens to single-digit-or-low-teens.

The remaining within-shape violations are dominated by mild
overshoots (most pairs at 1.5–5 % grade) at long baselines
(typical pair distance 20–60 m).  The plan's headline targets
(< 100 / < 500) aren't quite hit but the implementation is
within an order of magnitude — and the worst-case grade dropped
from 200–270 % to 13–50 %.  Further reductions would come from
more aggressive subdivision of the few large polygons that
account for most of the residual count, not from re-tuning the
mesh solver.

### What changed in the codebase

* `_runway_corner_elev_map` — runway-only anchor map.
* `_build_pavement_mesh` — global mesh constructor.
* `_solve_pavement_mesh` — iterated solver with rect/flat
  equalisation.
* `_writeback_pavement_mesh` — per-shape altitude tag rewrite.
* `_compute_elevations` — invokes the mesh solver after
  `_triangulate_junctions`, before the existing safety nets
  (clamp + subdivide), which now operate on already-smoothed
  values and rarely fire.

The `Mode B` per-polygon helper (`_smooth_polygon_grid`) and its
gate (`USE_PER_POLYGON_ELEVATION_FIELD`) are kept as a possible
finishing pass but are no longer needed by default — the global
mesh subsumes its function.

---

## Lessons learned (2026-04-26 post-implementation)

### What landed

* `_smooth_polygon_grid` helper (Mode B 2D smoother).
* Constants for grid step / iteration caps / `SHARED_AGREE_TOL_M`.
* `USE_PER_POLYGON_ELEVATION_FIELD` gate (default OFF — see
  "Mode B not engaged" below).
* `_vertex_elev_anchored` re-classified: graph- and DEM-sampled
  values are no longer treated as immutable anchors when they
  appear on a cross-junction shared bucket; only rect / runway /
  terminal corner sources stay HARD.  Smoothing is free to pull
  graph-derived shared boundary vertices into 2D Euclidean grade
  compliance with the polygon's real anchors.
* Single-pass cross-junction reconciliation: after every junction
  has smoothed its boundary, shared SOFT buckets whose smoothed
  values diverge across junctions by more than
  `SHARED_AGREE_TOL_M` (0.10 m default) are overwritten with the
  cross-junction average to restore the shared-vertex invariant.

### Measured impact

| airport | check_grade WITHIN-SHAPE / CROSS-SHAPE |   |   |
|---|---|---|---|
|         | baseline | new | delta |
| SPJC | 2512 / 5  | 2030 / 3  | -19 % within / -2 cross |
| CYXY | 2090 / 1  | 2181 / 1  | +4 % within / 0 cross |
| HECA | 11139 / 47 | 10021 / 56 | -10 % within / +9 cross |
| KBNA | 19230 / 43 | 13768 / 83 | -28 % within / +40 cross |

Net effect: the SOFT classification + averaging delivers a real
reduction in within-shape violations (notably at SPJC and KBNA),
at the cost of a small number of additional sub-meter cross-
shape steps where shared-bucket smoothing diverged enough to
trip the averaging.  The plan's headline targets (< 100 / < 500)
remain out of reach.

### Mode B not engaged

`USE_PER_POLYGON_ELEVATION_FIELD` is OFF by default.  Mode B's
2D-Euclidean anchor cones impose tighter grade constraints than
the elevation graph's network-distance smoothing, so when rect-
corner anchors fit network compliance but not 2D compliance,
Mode B's cell bands collapse and free cells snap to midpoints —
producing wider boundary-sample variance than the legacy
`_smooth_junction_boundary` produces with the same anchors.

The Mode B helper, constants, and integration code are retained;
re-engaging Mode B requires first reworking rect-corner
elevation derivation to enforce 2D-Euclidean grade compliance
between adjacent rects (so the cones don't collapse).

### Diagnosed but unfixed

The remaining within-shape violations are dominated by:

1. **HARD anchor disagreements** (e.g. two rect corners 4 m
   apart in 2D differing by > 1.5 %/m because the centerline
   graph that seeded their elevations is grade-compliant over
   network distance, not 2D).  Mode B (and any anchor-respecting
   smoother) cannot fix these without a separate rect-corner
   reconciliation pass.
2. **Within-junction multi-anchor conflicts** at narrow-waist
   polygons that span anchors with mutually infeasible 2D bands.
   `_subdivide_violating_junctions` partially handles these by
   splitting the offending junction along the conflict line; it
   could be made more aggressive.

### Possible follow-up directions

* Add a per-airport rect-corner reconciliation pass that walks
  every (rect, neighbour-rect) pair within ~30 m, finds 2D-grade-
  violating corner pairs, and adjusts their `altitude_high` /
  `altitude_low` by half the excess.  Iterate to fixed point.
  This is what makes Mode B viable.
* Replace the centerline elevation graph with a dense 2D grid
  (Stage 3 in "Out of scope") so per-cell grade compliance
  matches what Mode B expects.
* Tune `_subdivide_violating_junctions` to split more
  aggressively in the residual cases.
