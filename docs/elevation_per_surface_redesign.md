# Per-surface elevation redesign (2026-05-02)

## Problem

The current unified Laplacian solver in
`src/auto_patch/elevation.py::_solve_pavement_elevations_unified`
treats all pavement as one connected graph and propagates elevation
caps on per-edge `grade × length` budgets. This is geometrically wrong
because **FAA grade is a per-axis constraint along each surface, not a
graph-distance or Euclidean constraint**.

User's example (2026-05-02): a taxiway parallel to a runway 150 m
away, both 500 m long, taxiway climbing 1.5 % and runway descending
1.5 %. After 500 m the elevation delta is 15 m at 150 m perpendicular
— *no FAA violation* because no plane traverses the perpendicular
grass. The current algorithm forbids this and forces the surfaces to
converge.

**Observed at CYXY (2026-05-02):** taxiway E at the south edge of the
SW apron is at solver-output **707.9 m** but DEM truth is **717.2 m**
(runway nearby is 703.9 m). Taxi E should be ~14 m above the runway,
following the natural terrain; current algorithm has cut it down by
~10 m to satisfy a graph-distance grade cap from runway anchors that
shouldn't exist.

Within-shape audit also reports false positives: `714.1 → 707.9, d=7m,
de=6.2m` (88 % grade) on a junction that legitimately bridges a
high-terrain apron edge to a low-terrain runway-side edge.

## Authoritative rules (user 2026-05-02)

* Runway, taxiway (incl. rects: parallel / secondary_parallel / stub /
  cross_connector): 1.5 % grade **along their own centerline axis**.
* Junction: 1.5 % grade, but multi-directional (the polygon may slope
  in more than one direction across its surface).
* Apron: 1.0 % grade, multi-directional.
* Terminal: flat.
* Cross-surface continuity: only at **shared vertices** (exact
  equality). No Euclidean or graph-distance propagation through
  unrelated surfaces.

A plane only ever travels along one surface's axis at a time. Any
constraint that doesn't model a real travel path is wrong.

## Phased per-surface solver

### Phase 0 — Inputs (unchanged)

* CIFP runway profile → `altitude_high` / `altitude_low` on each
  runway segment (HARD anchors).
* DEM tile loaded via existing `_load_airport_dem` /
  `_sample_dem`.
* Layout has all shapes with polygons, roles, refs, source_axis on
  rects.

### Phase 1 — Taxi rect axial profiles

For each rect (`primary_parallel`, `secondary_parallel`, `stub`,
`cross_connector`):

1. Sample DEM at 3 points along source_axis: `t=0`, `t=0.5`, `t=1`
   (the rect's two short-end midpoints + axial midpoint).
2. If a short end is shared with a runway corner, replace that
   sample with the runway corner's HARD elevation.
3. Smooth the 3-point profile to satisfy `|delta_elev| / axial_dist
   ≤ 1.5 %` between consecutive samples. Where DEM exceeds this,
   clamp toward the anchored end (or toward the midpoint when both
   ends are free).
4. Write `altitude_high = max`, `altitude_low = min`. Use the
   short-end-pair geometry helper (`_short_end_pairs_by_axis`) to
   stay convention-free.

### Phase 2 — BFS propagation from runway

Iterate outward from runway-anchored rects:

1. Initialize a frontier of rects with at least one short end
   touching a runway corner. Their runway-anchored end is HARD.
2. For each frontier rect, propagate its non-runway end altitude to
   the adjacent junction (shared vertex).
3. Add rects on the other side of that junction to the frontier;
   those rects now have ONE end anchored (via the junction).
4. Continue BFS until all rects connected to runway anchors are
   processed.
5. Disconnected rects (no path to a runway) fall back to DEM only.

Per-rect altitude is the **min over** (DEM target, anchored-end +
1.5 % × axial_length). This is the user's "max highs and lows in
DEM that are possible within grade limits."

### Phase 3 — Junction interior elevations

For each junction:

1. Boundary vertices shared with rects/runways = anchored to those
   shapes' values (continuity).
2. Interior boundary-arc vertices (apt.dat polygon points between
   rect-shared corners) = sampled from DEM.
3. Smooth all junction `node_altitudes` to satisfy `|de| ≤ 1.5 % ×
   ring_edge_length` on RING edges only — **no within-shape
   spatial-pair edges**.

The junction may slope in multiple directions; this is correct.

### Phase 4 — Apron interior elevations

Same as Phase 3 but cap = 1.0 % on ring edges.

### Phase 5 — Terminals

Flat at DEM-median (existing rule preserved). Where a terminal
shares a vertex with an apron, that vertex feeds back into Phase 4
as a boundary anchor on the apron.

### Continuity (shared-vertex resolution)

Resolution priority when two shapes share a vertex:

1. Runway (HARD)
2. Rect (Phase 1 + Phase 2)
3. Junction (Phase 3)
4. Apron (Phase 4)
5. Terminal (Phase 5)

Higher-priority shape's value wins. Lower-priority shapes treat the
shared vertex as a boundary anchor.

## Module layout

New code lives in `src/auto_patch/elevation_per_surface/`
(package), keeping each file small:

```
elevation_per_surface/
├── __init__.py             ← public entry: solve(layout, dem, ...)
├── dem_targets.py          ← DEM sampling at 3 points along axis
├── axial_profile.py        ← Phase 1: per-rect axial smooth-to-grade
├── bfs_propagate.py        ← Phase 2: BFS frontier from runway
├── junction_field.py       ← Phase 3 + 4: ring-edge smooth (param cap)
└── continuity.py           ← Phase 5 + shared-vertex reconciliation
```

The existing `_compute_elevations` orchestrator stays. We add a new
branch behind a feature flag `USE_PER_SURFACE_SOLVER` that calls
into the new package instead of `_solve_pavement_elevations_unified`.
Once validated against SPJC + CYXY + HECA, the legacy solver and
its supporting helpers can be removed.

## Tests

* **CYXY taxi E regression** (new): assert taxi E at the south edge
  of the SW apron reaches ≥ 714 m (currently 707.9). Fixed
  numeric goal per user 2026-05-02.
* **Existing test_pavement_grade**: SPJC + SPLP must keep passing
  (or come closer to passing — currently failing for out-of-scope
  reasons).
* **Existing test_junction_invariants + test_compare_target**:
  must remain passing.
* **Existing test_junction_rules** (5 rules): unchanged.

## Out of scope

* Re-orienting polygon vertex order to match the legacy
  `[high, low, low, high]` convention — downstream rendering bugs
  on rotated polygons are pre-existing and tracked separately.
* Rule 1 v6 — 2 SPJC junctions still at 1 shared node. Picked back
  up after the elevation redesign lands.
