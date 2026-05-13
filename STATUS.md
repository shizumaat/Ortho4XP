# Auto-Patch Status — 2026-05-13: SEAM-ANCHOR ARCHITECTURE + DIAGONAL-STUB TRAPEZOID

## TL;DR

**Major architecture change this session.**  Replaced the tile-cut
bridge polygon mechanism with **seam DEM HARD anchors** that win
priority over CIFP runway thresholds.  Pavement at the integer
lat/lon tile boundary is now pinned to the raw HGT pixel value
(deterministic across both adjacent tile builds via SRTM 1-pixel
overlap), and the elevation solver propagates FAA-grade-compliant
altitudes outward from the seam.  X-Plane's terrain mesh — pinned to
the same raw HGT at the seam via Ortho4XP's
``preserve_boundary=True`` smoothing — now matches pavement at the
seam by construction.  No more visible cliff.

Also landed: diagonal-stub **trapezoid** emission for digit-ref
sub-rect refs whose pavement flares wide at the apron end
(SPJC's V3 diagonal connector emits correctly now).

**Test status:** 217 / 217 pass.  Both fixture-based compare-target
gates (SPJC + SPLP — the latter newly added) validate clean.

**Branch:** `dev` (3+ commits ahead of last status).
**Tag baselines:**
* `spjc-good-baseline` (commit ``4d97f89``, 2026-05-07) — pre-seam.
* (planned for this session) — post-seam SPJC/SPLP fixtures.

## Commits this session (since 9b7164d STATUS 2026-05-12)

```
(pending commit — large architectural change)
```

The pending commit contains:
* Diagonal-stub trapezoid fallback in ``_rect_from_axis_extended``
  (accept_asymmetric mode for digit refs with parent db ≥ 20°).
* Classifier tightening: digit→STUB now requires ``db_local ≥ 15°``.
* New module ``src/auto_patch/seam_anchors.py``:
  ``split_pavement_at_seams`` inserts seam vertices, converts seam-
  affected sloped rects (and all sub-rects of an affected runway
  chain) to ``node_altitudes``, records seam-vertex bucket keys;
  ``apply_seam_dem_anchors`` writes ``dem.alt_strict`` values into
  ``node_altitudes`` at every seam vertex.
* New module ``src/auto_patch/runway_regrade.py``:
  Stage A regrade solver.  Re-optimises runway threshold altitudes
  against seam HARD anchors with FAA grade cap (1.5 %) and K-factor
  (305 m, ARC Cat C/D default) constraints.  Closed-form 1D-per-
  threshold optimisation with graceful relaxation.
* Solver changes in ``unified_jacobi.py``:
  ``_seed_elevations`` does 2-pass runway HARD seeding (CIFP first,
  regraded shapes override) + seam-HARD override (seam wins over
  CIFP).  ``_writeback`` preserves ``node_altitudes`` representation
  for shapes that came in with it; handles runway shapes with
  ``node_altitudes`` (the seam-converted ones).
* Pipeline integration in ``pipeline.py``: seam pipeline runs
  between DEM load and ``per_surface_solve``.
* Tile-cut bridge polygon mechanism removed entirely
  (``ROLE_TILE_CUT_BRIDGE`` deleted, bridge emission block removed,
  ``box`` import dropped from ``tile_cut.py``, grade-limit entry
  removed from ``config.py``).
* ``_resample_node_altitudes_nn`` upgraded to **edge interpolation**:
  cut-edge vertices now use linear gradient of the underlying old
  edge instead of nearest-neighbour, eliminating jumpy NN artefacts
  at cut buffer boundaries.
* New unit tests: ``tests/test_runway_regrade.py`` (8 cases, all
  pass) covering no-seam / single-seam / two-seam / K-factor /
  grade-cap edge cases.
* Compare-target test now covers SPLP in addition to SPJC.
  Baselines reset to current canonical builds (regenerated 2026-05-13).
* Tile-cut parity test simplified: bridge-vertex test removed
  (bridges gone); seam parity now guaranteed by construction.

## Build / verify commands

```bash
# Full test suite (217/217 pass).
venv/bin/python3 -m pytest tests/ --tb=short -q

# Build SPJC + write OSM for visual review.
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
for icao in ("SPJC", "SPLP"):
    layout = build_airport_pavement(
        icao, "/Users/noah/X-Plane 12", compute_elevations=True)
    layout.to_osm(f"/tmp/{icao}.osm")
    print(f"wrote /tmp/{icao}.osm — {len(layout.shapes)} shapes")
PY

# Cross-tile regression (needs SPLP + S13W077 + S13W078 hgt).
venv/bin/python3 -m pytest tests/test_tile_cut_parity.py -v

# Compare-target gates (both airports).
venv/bin/python3 -m pytest tests/test_compare_target.py -v
```

## Seam-anchor architecture

The seam-anchor pipeline runs inside ``build_airport_pavement`` after
DEM load and before ``per_surface_solve``:

```
1. split_pavement_at_seams(layout)
   For every pavement shape whose polygon boundary intersects any
   integer lat/lon line in the airport footprint:
     a. Insert ring vertices at the intersection points (deterministic
        from polygon geometry alone, identical in both tile builds).
     b. Convert sloped 4-corner rects to node_altitudes representation
        — seam vertices get linearly-interpolated placeholder altitudes.
     c. For each runway whose chain has ANY seam-crossing sub-rect,
        convert ALL sub-rects in the chain to node_altitudes (so the
        shared corners between regraded and non-regraded sub-rects can
        carry independent per-vertex altitudes).
     d. Record every seam vertex's bucket key in layout._seam_anchor_keys.

2. apply_seam_dem_anchors(layout, dem, tile_lat, tile_lon)
   For every vertex in layout._seam_anchor_keys, overwrite the
   placeholder altitude with dem.alt_strict at that vertex's lat/lon.
   alt_strict (vs alt_nostrict) ensures both tile builds sample the
   identical SRTM pixel — SRTM .hgt files include the boundary pixel
   in both adjacent tiles' arrays.

3. regrade_runways_in_layout(layout, dem, tile_lat, tile_lon)
   For each runway sub-rect with seam vertices:
     a. Identify threshold corners (non-seam vertices) and seam
        vertices, projected onto source_axis.
     b. Call regrade_runway with CIFP-seeded threshold altitudes and
        DEM-anchored seam altitudes.
     c. Stage A QP solves for adjusted threshold altitudes that
        minimise CIFP-deviation subject to grade cap + K-factor.
     d. Write adjusted altitudes back to the shape's node_altitudes.
   Per user 2026-05-13 design: SEAM WINS over CIFP at conflict points.

4. per_surface_solve(layout, dem, tile_lat, tile_lon)
   The unified Jacobi:
     - 2-pass HARD seeding: CIFP-only shapes first, regraded shapes
       (those with node_altitudes from the seam pipeline) second
       with OVERRIDE priority.
     - Seam-HARD override: walks layout._seam_anchor_keys and pins
       matching solver-graph nodes to their DEM altitudes (wins
       over runway CIFP if both apply at the same node).
     - Cap projection propagates SOFT nodes from HARDs respecting
       FAA grade limits.
     - Writeback preserves node_altitudes when the shape came in
       with them; converts to altitude_high/low only for CIFP-only
       4-corner rects.
```

## Geometry state at SPJC + SPLP (current dev tip)

**SPJC (single-tile):**
* 896 total shapes (was 894 pre-V3-diagonal-trapezoid)
* boundary 603, junction 40, runway 91, primary_parallel 28,
  stub 20, cross_connector 9, retaining_wall 66, tunnel_ramp 36,
  secondary_parallel 1, terminal 2
* 1 within-shape grade WARN (pre-existing, junction 35.5→33.7
  d=22.4m de=1.8m → 8.0 %)
* Seam pipeline: 4 seam verts, 4 DEM-anchored, 0 runways regraded
  (airport boundary touches a seam but the runway doesn't)

**SPLP (cross-tile, lon=-77 cut):**
* 155 shapes (boundary 123, runway 23, junction 4,
  primary_parallel 2, stub 3)
* Seam pipeline: 22 seam verts, 25 DEM-anchored, 3 runway
  sub-rects regraded
* Runway 02/20 threshold shifts: A −1.47 m, B +2.70 m
* 14 within-shape grade WARN (all on the big 65-vertex junction
  spanning the coastal-terrain gradient; these are all-pair
  Euclidean violations on naturally varying terrain — solver's
  adjacency-graph grade is satisfied)

## Investigated and resolved this session

### Diagonal V3 stub at SPJC

**Pre-fix:** 5 V3 centerlines (50/53/140/57/53 m at bearings
152/141/126/142/148°) entered ``_build_taxi_rects``; the 140 m
diagonal at bearing 126° failed ``_rect_from_axis_extended``'s
symmetry retry loop (width asymmetry 42-47 %, never converged) and
was silently discarded.  Junction ``-10169`` then drew a 100+ m
straight edge across the diagonal pavement.

**Fix:** Added ``accept_asymmetric=True`` fallback in
``_rect_from_axis_extended`` — returns the most-symmetric snapped
quadrilateral encountered across all retries when the digit-ref's
parent way is diagonal (``db_axis ≥ 20°``).  Also tightened the
digit→STUB classifier at ``rects.py:1251`` to require
``db_local ≥ 15°``, so near-parallel V3 sub-segments classify
correctly.  V3 now emits as a 6461 m² trapezoid at bearing 300°
(split into 2 sub-rects post-emit), and junction ``-10169``
shrinks from 63 corners to 5.

### Cross-tile seam (SPLP)

**Pre-fix:** Tile-cut bridge polygons covered the 10 m gap with
``slope_sampler``-interpolated altitudes.  Empirical measurement
showed the visible step lives outside the gap (at the pavement-to-
DEM transition) where bridges have no effect.

**Fix:** Switched to seam-DEM HARD anchors (see architecture
section above).  Bridges removed.  Cross-tile parity:
* Seam vertices: 0.000 m worst |dz| (by construction).
* Cut-edge vertices (5 m off seam): typically ≤ 1 m
  (edge-interpolation resampler).

## Reference artefacts

* ``tests/fixtures/SPJC_target.osm`` — canonical SPJC build,
  regenerated 2026-05-13.  Forward baseline gate for SPJC.
* ``tests/fixtures/SPLP_target.osm`` — canonical SPLP build,
  regenerated 2026-05-13.  Forward baseline gate for SPLP.
* ``/tmp/SPJC.osm`` / ``/tmp/SPLP.osm`` — most recent local builds.
* ``+60-140/`` directory at repo root: 24 stray ``.hgt`` files
  that should live under ``Elevation_data/+60-140/``.  Per user
  direction: leave untracked.

## OPEN WORK — where to pick up

### CYXY missing junctions

User report 2026-05-13: CYXY appears to be missing almost all
junctions.  Not yet investigated.  Pick this up next session.

### In-sim validation (SPLP)

Rebuild SPLP via the full Ortho4XP pipeline (need
``Patches/-13-077/SPLP_auto.patch.osm`` +
``Patches/-13-078/SPLP_auto.patch.osm`` regenerated post-this-
session's changes) and walk the seam in X-Plane.  Architecturally
the cliff should be gone; visual confirmation pending.

### Within-shape grade WARN at SPLP

14 violations on the big 65-vertex junction that spans steep
coastal terrain (54-80 m across the airport).  These are all-pair
Euclidean violations between non-adjacent vertices on a single
polygon — the solver's adjacency-graph grade IS satisfied locally,
but ``check_grade``'s all-pair check finds far-apart vertex pairs
with naturally varying terrain altitudes.  Two paths if needed:
(a) split the big junction along the steep-terrain gradient lines
so each sub-junction's span is smaller, or (b) change
``check_grade`` to use axis-aligned grade (which is what FAA
actually mandates).  Not urgent — these are WARN, not FAIL.

## Memory notes for next agent

Saved under ``~/.claude/projects/-Users-noah-Ortho4XP-shred86/memory/``:
No new memory entries this session — all decisions captured in this
status file + commit messages.  Existing memories
``feedback_root_cause_only``, ``feedback_general_solutions``,
``feedback_shape_rules``, ``project_target_osm``,
``feedback_boundary_clamp_asymmetric`` remain authoritative.

User design rules established this session:
* **Seam wins over everything** — tile-boundary DEM is the
  canonical hard anchor.  CIFP runway thresholds adjust to satisfy
  seam + grade.  Other pavement adjusts to satisfy runway + seam +
  grade.  Visible cliff in X-Plane is unacceptable.
* **FAA grade rules everywhere when possible** — runway: 1.5 %
  longitudinal + K-factor vertical curves (default K = 305 m
  for ARC Cat C/D).  Taxiway/junction: 1.5 % grade cap.
  Apron/terminal: 1.0 % grade cap.  Relax only when seam HARD
  anchors force it.
