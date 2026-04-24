# Auto-Patch Refactor — Status

**Current state:** The new pavement builder
`src/O4_Airport_Pavement_Builder.py` is **WIRED INTO Ortho4XP's
main pipeline** — `O4_Auto_Patch.generate_auto_patches` now
calls `build_airport_pavement(icao, xplane_root)` and writes
the output via `layout.to_osm()` directly to the tile's
`Patches/` directory.  Phase-1 geometry + Phase-2 elevations
both produced in one pass; the legacy surface generator
(buildings, aprons, drainage, boundary band, tunnels) is no
longer called.

**SESSION 7 (2026-04-22 → 2026-04-24):** Sequence of changes:

1. Residue-driven junctions (simplest polygons, no overlaps).
2. Graph-based grade-compliant elevation network.
3. Runway / taxi / terminal elevation tags, X-Plane convention
   (altitude_high/low on corners 0,3 of the ring, cell_size +
   profile spline tags on sloped rects).
4. Multipolygon relations DROPPED — X-Plane's patch parser
   iterates ways only, so junction tags wouldn't reach the
   outer way.  Junctions now emit as simple closed rings.
5. Junction vertices no longer land mid-edge of any taxi rect
   (2 m gap to runway + vertex-level snap/push for other taxis).

**SPJC:** 155 shapes (73 runway + 33 junction + 47 taxi rects
+ 2 terminals); 0 overlaps; 0 mid-edge junction-vertex inserts
on any rect; grade compliance verified on V1 / A / F stubs.

**SPLP:** 57 shapes; 0 overlaps.

## Session 7 iteration 2026-04-24 (late) — Ortho4XP integration + slope-rect fixes

Five commits on 2026-04-24:

- **`c0e1b24`** Graph-based elevation network + runway multi-
  polygon clip fix.
- **`5cca53a`** Corner ordering for altitude_high/low tags —
  rebuild rect ring from 4 raw corner positions in the X-Plane
  patch convention `[hi_left, lo_left, lo_right, hi_right]`.
  Runway segment construction also reorders when B is higher.
  Result at SPJC: 47 / 47 sloped taxi rects + 67 / 67 runway
  segments correctly oriented.
- **`98628db`** Wired `build_airport_pavement` into
  `O4_Auto_Patch.generate_auto_patches`.  Removed multipolygon
  relations from `to_osm()` output (X-Plane patch parser can't
  route tags from relations to their outer ways).  Added
  `cell_size=2` + `profile=spline` tags on sloped rects to
  match legacy patch format.  `Ortho4XP_AutoPatch` naming kept
  so manual-patch priority still works.
- **`039ad60`** Post-runway-clip with a 2 m outward buffer on
  the new segmented runway union, subtracted from every
  junction — fixes 15 junction vertices that were landing
  mid-edge on runway short edges (would have split the
  runway's 4-corner slope rect at render time).
- **`e9dc9aa`** Vertex-level junction-vs-taxi-rect pass.  For
  each junction vertex: if within 2 m of a taxi rect corner,
  snap exactly to that corner; if within 0.5 m of a rect edge
  (not near a corner), push 1 m perpendicular-outward (tested
  via rect containment for the correct side).  Fixes G stub
  and F primary "junction wrapping around both sides" report.
  11 → 0 mid-edge inserts at SPJC.

### Ortho4XP wiring

[O4_Auto_Patch.py:1393-1440](src/O4_Auto_Patch.py) now:

```
xp_root = xplane_root_from_cifp_path(cifp_path)
from O4_Airport_Pavement_Builder import build_airport_pavement
layout = build_airport_pavement(icao, xp_root)
layout.to_osm(auto_patch_file)
```

Replaces the legacy `generate_patch_osm` + `generate_airport_
surface_patches` chain.  Manual patches still take priority via
the `manual_patches` set in `generate_auto_patches`.  Input
fields (`taxiway_data`, `building_data`, `dico_airports`,
`road_data`) are still accepted for API compatibility but
ignored by the new path — the new builder pulls taxiways,
terminals, runway geometry, and CIFP elevations directly from
disk via `_load_osm_airports`, `_extract_osm_terminals`,
`O4_Apt_Dat_Reader`, and `parse_cifp_file`.

### Out-of-scope this session / known-work-ahead

- **Junction / apron elevations.**  Currently un-elevated —
  X-Plane renders them against the DEM.  User plan: junctions
  will triangulate multi-directionally using shared-vertex
  elevations from adjacent rects.  **Next session.**
- **Building / hangar footprints.**  Only terminals emit as
  role=terminal pads.  General buildings + hangars removed in
  the earlier revert.  Later.
- **Legacy drainage / boundary band / tunnels.**  Not
  reproduced in the new builder.  Later.
- **Runway 1 % / 305 m vertical curve rule.**  Currently only
  the 1.5 % longitudinal cap is enforced; the curvature rule
  would require a second-derivative smoothing pass.  Not
  visible at SPJC / SPLP given the segment lengths used.

### Current constants (post-2026-04-24)

- `NETWORK_DENSIFY_M = 30.0`; `NETWORK_BRIDGE_MAX_M = 60.0`.
- `TAXI_MAX_GRADE = 0.015`; `_RWY_HALF_WIDTH_M = 22.5`.
- Runway-vs-junction clip: 2.0 m buffer on new runway union.
- Junction-vertex-vs-taxi-rect: snap within 2.0 m of a corner;
  push 1.0 m perpendicular-outward if within 0.5 m of an edge.
- Sloped rect tags: `altitude_high`, `altitude_low`,
  `cell_size=2`, `profile=spline`.

---

## Session 7 iteration 2026-04-24 — graph-based elevation network

### User-reviewed issues from 2026-04-23 iteration

1. **Missing junction at F stub / runway-34R end.**  A ~7400 m²
   pavement region between the F walking-exit stub and the new
   runway-34R blast-pad segment was uncovered.  Root cause: the
   Phase-2 runway-clip code subtracted the new (segmented, with
   overruns) runway union from each pre-existing junction
   polygon, then kept only the largest resulting sub-polygon.
   For junctions that straddle the new runway's east edge, the
   clip returns a MultiPolygon — the piece east of the runway
   (where the gap is) was dropped as "not the largest".
2. **Stubs' elevations differed 2-5 m from the runway they
   meet.**  V1 / A / F stubs measured +0.8, +4.2, +2.4 m above
   the runway they join, violating the FAA 1.5 % taxi grade cap
   by wide margins.  Root cause: the per-rect elevation logic
   anchored to runway elevation only when the axis endpoint was
   within 30 m of a runway segment.  Perpendicular stubs' axes
   end 40-150 m off runway (the junction fills the gap), so no
   anchor snap applied and the stubs used raw DEM which diverges
   from the actual airport ground plane.

### Changes

1. **Multi-polygon runway clip fix.**  [_compute_elevations
   runway-replacement block]: after clipping each non-runway /
   non-terminal shape against the new runway union, emit EVERY
   sub-polygon above a 50 m² floor (not just the largest).
   Preserves the original shape's metadata on the largest piece
   and uses `dataclasses.replace` to spawn extra shapes for the
   remaining pieces.  Fixes the F/34R gap immediately (uncovered
   pavement 8447 → 1025 m² at SPJC; F/34R-specific uncovered:
   4848 → 3 m²).

2. **Graph-based elevation network.**  Replaced the per-rect DEM
   sampling with a shared elevation graph.  The pipeline:

   - **`_build_elevation_network`** constructs a densified-to-
     30 m graph from OSM taxi centerlines + CIFP-anchored runway
     segment endpoints, plus **bridge edges** connecting any
     non-anchor taxi graph node to the nearest runway graph node
     if within 60 m (tight enough that interior nodes of parallel
     taxis, typically > 100 m off runway at SPJC, don't get
     falsely anchored).  Bridge-edge length is the straight-line
     distance **minus the 22.5 m runway half-width** since grade
     doesn't accrue crossing the flat runway surface width.
   - **`ElevationGraph.propagate_bounds`**: Dijkstra from each
     hard anchor (CIFP threshold + densified runway centerline
     nodes) with edge weight = length × 0.015, producing a
     per-node feasibility interval = ∩ of (anchor ± distance ×
     0.015) across all reachable anchors.
   - **`ElevationGraph.choose_values`**: pick each node's
     elevation as DEM clipped into its interval (per user
     2026-04-24, DEM is a soft preference — real airports
     excavate / build up, so a node free of runway anchors
     follows DEM only when no grade-compliance rule forces
     otherwise).
   - **`ElevationGraph.smooth_rate_of_change`** (up to 30
     iterations): Laplacian-style moving each non-anchor node
     toward neighbours' mean by damping factor 0.4, re-clipping
     into the feasibility interval, and hard-capping each edge
     at 1.5 % afterward.  Addresses the FAA 1 %/30 m
     rate-of-change rule (curvature) as a follow-up to the
     first-derivative grade cap.
   - **`ElevationGraph.elevation_at(x, y)`** samples the network
     surface at arbitrary meter-space points — used for each
     rect's two axis endpoints.

### Results at SPJC

Graph-based path-distance grade compliance for the three stubs
that were previously badly out:

| stub | stub elev | runway anchor | taxi-path dist | Δ | FAA allowed | |
|---|---|---|---|---|---|---|
| V1 | 6.32 m | 5.88 m | 52 m | 0.44 m | 0.79 m | ✓ |
| A  | 14.28 m | 13.41 m | 66 m | 0.87 m | 0.99 m | ✓ |
| F  | 33.37 m | 32.19 m | 87 m | 1.19 m | 1.30 m | ✓ |

All three now within the FAA 1.5 % cap along their actual
graph-network path to the runway anchor.  Overlap count
remains 0 at both SPJC and SPLP.  Coverage 0.04 % uncovered.

### Still-pending / known

- **Terminals**: still DEM-median.  SPJC terminal1 at 31.1 m
  may still be high vs the airport's 13 m reference — the
  terminal pad's apt.dat polygon extends into hilly ground
  east of the apron, and the DEM median picks that up.  Will
  need a building-pad grade constraint later (pulled toward
  adjacent taxi elevations via a similar graph approach).
- **Junctions / aprons / buildings** still un-elevated per user
  instruction.  When later enabled, each junction vertex
  inherits elevation from the adjacent taxi/runway/terminal
  shapes that share that vertex.
- **Bridge-edge threshold at 60 m** works for SPJC/SPLP where
  no parallel taxi passes closer than ~100 m to the runway.
  If future airports have tighter geometry (parallel within
  60 m), the threshold may need per-airport tuning.
- **Runway half-width hard-coded at 22.5 m** (SPJC/SPLP are
  both 45 m runways).  Should use the actual width from
  apt.dat per runway when we touch airports with non-45 m
  runways.

### Constants (elevation graph)

- `NETWORK_DENSIFY_M = 30.0` (max edge length, matches FAA
  rate-of-change rule 1 % / 30 m).
- `NETWORK_BRIDGE_MAX_M = 60.0` (max dist for taxi→runway
  bridge).
- `TAXI_MAX_GRADE = 0.015` (FAA taxi cap).
- `_RWY_HALF_WIDTH_M = 22.5` (offset in bridge-edge length).
- Smoothing: 30 iters, damping 0.4, convergence tol 1 cm.

---

## Session 7 iteration 2026-04-23 (phase 2 — elevations)

User instruction:
- Add elevations to **terminal pads** (flat).
- Re-enable the existing segmented, sloped runway generation.
- Add elevations to **all taxi rects**.
- Leave junctions / aprons un-elevated this iteration.

### Changes

1. **BuiltShape** extended with optional ``altitude``,
   ``altitude_high``, ``altitude_low`` floats.  ``to_osm`` writes
   them as tags (single ``altitude`` when flat; high+low pair
   when sloped).

2. **New `_compute_elevations(layout, icao, xplane_root, apt)`**
   runs as the final pass of ``build_airport_pavement`` when
   ``compute_elevations=True``.  It:
   - Loads the tile DEM from ``Elevation_data/{group}/{hgt}``
     via ``O4_DEM_Utils.DEM``.
   - Parses CIFP via ``O4_Auto_Patch.parse_cifp_file`` and
     calls the legacy ``pair_runways`` + ``generate_patch_osm``
     with a ``Tile`` stub carrying the DEM.
   - Takes the ``runway_segment_chain`` output (list of
     ``(lat_a, lon_a, elev_a, lat_b, lon_b, elev_b, width_m)``
     tuples), **strips the 3 m per-side ``RUNWAY_MARGIN``** so
     the segmented runways match the apt.dat width our Phase-1
     junctions were built against, and builds per-segment
     BuiltShape rects with ``altitude_high``/``altitude_low``
     (or flat ``altitude`` when |Δ| < 0.1 m).  Replaces the
     original single-rect runway shape.
   - Taxi rects (primary / secondary / stub / cross_connector):
     sample DEM at axis endpoints, snap to the nearest runway
     segment elevation when within ``TAXI_ANCHOR_DIST_M = 30 m``,
     cap longitudinal grade at FAA ``1.5 %`` taxi limit, emit
     ``altitude_high/low`` or flat ``altitude``.
   - Terminal pads: sample DEM at centroid + perimeter,
     emit median as flat ``altitude``.
   - Junctions untouched (no altitude tags).
   - Final **runway-vs-other clip** using
     ``new_runway_union.buffer(0.05)`` removes floating-point
     slivers (e.g. the flat overrun segments extending past the
     original runway endpoints into existing junctions).

3. **Defensive geometry cleanup** in the elevation pass: calls
   ``buffer(0)`` on each polygon before ``difference`` to
   avoid shapely's "side location conflict" errors at
   self-kissing boundaries produced by `_enforce_shared_vertices`.

### Results (2026-04-23 phase 2)

| metric | SPJC | SPLP |
|---|---|---|
| runway shapes | **73** (67 sloped + 6 flat) | **24** (21 sloped + 3 flat) |
| taxi rects elevated | 47 / 47 ✓ | 15 / 15 ✓ |
| terminal pads elevated | 2 / 2 ✓ | 0 (none in target) |
| junctions elevated | 0 (by design) | 0 (by design) |
| overlap pairs (> 1 m²) | **0** ✓ | **0** ✓ |
| runway elevation range | 5.9 – 32.2 m | 65.3 – 77.1 m |

### Deferred / known-quirks

- Per-rect DEM sampling is independent; **cross-rect grade-change
  relaxation (FAA 1 %/30 m vertical curve)** is not applied —
  shared vertices between adjacent taxi rects may carry different
  computed altitudes, X-Plane triangulator will smooth.  If
  visible artefacts appear, a post-pass graph relaxation is the
  fix.
- Terminal elevation is DEM-median — SPJC terminal1 lands at
  31.1 m, terminal2 at 22.9 m.  These look high vs the airport
  reference (~13 m at the N threshold); the SPJC site sits on a
  NW-SE slope and the DEM (SRTM-derived) may include
  building-top returns.  Flag for visual review.
- Junctions / aprons / buildings will get elevations in a later
  iteration.
- ``_validate_shared_vertex_invariant`` is still ACTIVE — pre-
  elevation shapes pass it; runway segments are added AFTER so
  the check runs before the new runways exist.

### Current constants (elevation-related)

- ``TAXI_MAX_GRADE = 0.015`` (FAA taxiway longitudinal cap).
- ``TAXI_ANCHOR_DIST_M = 30.0``.
- Legacy ``RUNWAY_MARGIN = 3.0`` stripped from segmented runway
  width to match apt.dat.
- DEM source: ``Elevation_data/{±NN±NNN}/{S|N}NN{E|W}NNN.hgt``.
- CIFP source: ``{xplane_root}/Custom Data/CIFP/{ICAO}.dat``.

---

## Session 7 iteration 2026-04-23 (late) — pivot to simplest-polygons

### Strategic pivot (user 2026-04-23)

Junctions and aprons will be treated identically at elevation
time: both will triangulate and slope multi-directionally
(unlike rects, which slope only along their axis).  So the
distinction no longer exists geometrically.  New rule: emit the
**simplest polygon that covers all remaining pavement and
connects to neighbouring rects + terminals**.  Don't try to
match target's hand-drawn subdivision any more — it was
over-subdivided for a different elevation model.

### Changes

1. **Cluster-based constructive junction builder DELETED from
   the emission path.**  Seed-driven clusters produced
   overlapping, partial junctions when adjacent clusters each
   captured a different subset of rect corners.  Now residue-
   driven only.

2. **Residue = `pav_union - rect_union - terminal_union -
   runway_union`.**  Each connected component is emitted as
   one junction.  No area threshold needed beyond slivers
   (`MIN_JUNCTION_AREA_M2 = 50`).

3. **Rect + terminal + runway corners injected as shared
   boundary vertices** on every junction.  Uses extended
   `_insert_points_on_boundary` that now preserves interior
   rings (critical for junctions that wrap around interior
   rects, creating holes).

4. **OSM writer supports multipolygon relations** for
   polygons with interior rings.  Junction outer + each inner
   ring emit as separate ways; a relation ties them with
   role=outer / role=inner.  [tools/compare_target.py](tools/compare_target.py)
   updated to parse multipolygons correctly.

5. **OSM writer bucket snapping fixed** to use the first
   vertex's actual coordinates (not bucket center).  Bucket
   centers were up to 0.25 m off the real vertex, creating
   sub-metre overlaps after OSM round-trip.  Now zero.

### Results (2026-04-23 late)

| metric | SPJC | SPLP |
|---|---|---|
| shapes total | 82 | 33 |
| junctions (incl. multipolygons) | 31 | 17 |
| overlap pairs (>1 m²) | **0** ✓ | **0** ✓ |
| uncovered pavement | 0.04 % | similar |
| unmatched rect corners | 0 ✓ | 0 ✓ |
| unmatched terminal corners | 0 ✓ | 0 ✓ |

Target match-counts are lower than previous iterations because
target hand-draws finer subdivisions.  This is expected under
the new approach.

### Next steps

1. **Elevations (Phase 2).**  With the geometry now clean and
   seamless, the next phase is adding elevations to each
   shape.  Rects slope along their axis; junctions / terminals
   triangulate and slope multi-directionally.
2. **Re-enable building outlines** (currently only terminal
   pads emit).  The building emission path already exists via
   `_extract_osm_terminals` but is limited to pad polygons.
3. **SPLP refinement**: 5 spurious junction pieces (small
   slivers where rects don't perfectly snap).

### Current constants

- `EMIT_JUNCTIONS = True`; `EMIT_APRONS = False` (the latter
  is dead code — no path emits aprons any more).
- `MIN_JUNCTION_AREA_M2 = 50.0`.
- `SIMPLIFY_TOL_M = 1.0`.
- `RECT_CORNER_TOL_M = 2.0` (corner-to-boundary injection
  range).
- `_snap_corners_to_pavement`: filters apt.dat vertices to
  those within 0.5 m of `pav_union.boundary`.

---

## Session 7 iteration 2026-04-23 (early) — V-chain review + Q/R end junctions

### User-reviewed issues (SPJC, V-chain N→S)

1. **Duplicate runway-stub trapezoid** at every widening-stub →
   runway location.  Two overlapping junctions per runway end —
   one 4-vert trapezoid, one boundary-tracing residue junction.
   Fixed: dropped trapezoid path at `_build_airport_pavement`
   line ~955; kept `rwy_vertex_inserts` collection so runway
   polygons still pick up shared vertices with the residue
   junction.  Overlap pairs 13 → 1 at SPJC.

2. **V1 rect corner 5.88 m off pav boundary.**  Corner snap
   picked an apt.dat vertex that sits INSIDE the pav_union
   interior (where two pavement polygons overlap).  Fixed:
   `_snap_corners_to_pavement` now pre-filters `apt_vertices`
   to those within 0.5 m of `pav_union.boundary`.  All V-chain
   rects (V, V1, V2, V3, V5) now have all 4 corners exactly on
   boundary.

3. **V2 + V-south rect snap** — same fix as #2.

4. **V/Q/R missing NW extent + mid-edge R connection.**
   Constructive cluster at (327,-149) only captured SE rect
   corners (Q, R) — V-primary-north rect end was >120 m from
   cluster centroid, beyond `max_corner_dist_m`.  **See Q/R
   end radius fix below.**

5. **V3 mega-junction (-803794, 35v, 74668 m²)** not emitted.
   Root cause: L primary isn't subdivided at apron vertices,
   so no rect ends seed a cluster in this region.  Same
   root cause as 9 missed primaries (open item #3 from session
   6).  **Deferred** — needs apron-vertex-based primary split.

6. **V-apron junction (-803796)** extends east past U.
   Tried: chord-bulge cap, residue-hull cap — both regressed
   match count by 4 without cleanly fixing the visual issue.
   Reverted.  Still oversized by ~2× (IoU 0.27).  **Deferred**
   — probably needs proper apron emission to bound the east
   side.

### Q/R end junctions fix (user 2026-04-23)

Target junctions at the N and S ends of Q and R each connect
5 rects.  Their rect ends span 300-450 m diagonally.  The
previous `max_corner_dist_m` default of 80 m (sometimes 120 m
with extent slack) only reached the closest 2-3 rects, so the
junction polygons under-covered.

Changes:

- **Corner radius `CORNER_RADIUS_BASE_M = 220 m`** (was 80)
  at `_build_airport_pavement` cluster-build loop.  Allows
  the constructive builder to project rect corners up to 220 m
  from the cluster centroid — enough for the 5-rect junctions.
- **Local disc `local_radius = corner_radius + 60 m`** so the
  boundary walk has a slightly larger window than the corners.
- **`CLUSTER_MERGE_DIST_M: 40 → 80`** so adjacent seed points
  at a single physical junction merge into ONE cluster rather
  than producing 2-3 overlapping junction builds.
- **Dedup overlap threshold `0.8 → 0.3`** at `_polys_overlap_heavily`
  so bigger clusters producing variant polygons around the
  same target junction consolidate to the largest.

### Results (2026-04-23)

| Target | Before | After (IoU) |
|---|---|---|
| -803800 L/M/Q/R (N end, 5 rects) | 0.24 | **0.75** ✓ |
| -803799 V/Q/R (S end, 5 rects) | 0.31 | **0.57** |
| -803801 Q/R N-mid | ~0.40 | **0.88** ✓ |
| -803802 Q/R S-mid | ~0.40 | **0.89** ✓ |

| | junctions matched | missed | spurious | avgIoU |
|---|---|---|---|---|
| SPJC before | 31/43 | 12 | 18 | 0.54 |
| SPJC after | **38/43** | 5 | 18 | 0.56 |
| SPLP before | 12/13 | 1 | 12 | 0.59 |
| SPLP after | 12/13 | 1 | 5 | 0.58 |

### Still-pending at session-7 close

1. **V-apron overflow** (issue 6).  Needs apron emission or a
   better apron-symmetry cap.
2. **V3 mega missing** (issue 5).  Needs L-primary split at
   apt.dat apron-boundary vertices.
3. **18 spurious SPJC junctions.**  Big clusters now produce
   multiple variants that don't dedup at 0.3 threshold.  Next
   step: either tighter dedup on the dominant-area pair, or
   cluster-consolidation pass after build.

### Current session-7 constants

- `EMIT_JUNCTIONS = True`; `EMIT_APRONS = False`.
- `MIN_JUNCTION_AREA_M2 = 1000.0`.
- `CLUSTER_MERGE_DIST_M = 80.0` (was 40).
- `CORNER_RADIUS_BASE_M = 220.0` (was implicit 80).
- `LOCAL_DISC_EXTRA_M = 60.0`.
- Dedup threshold `0.3` (was 0.8).
- Arc-vertex limit per arc: 4 (unchanged).
- `_snap_corners_to_pavement`: pre-filters apt.dat vertices
  to those within 0.5 m of `pav_union.boundary`.

---

## Session 7 original (2026-04-22) — junction re-enable

### Changes landed

1. **`EMIT_JUNCTIONS = True`** (line 72).  Rects judged within
   tolerance at session 6 close; junction emission re-opened for
   iteration.  Aprons stay off.

2. **Area floor `MIN_JUNCTION_AREA_M2: 80 → 1000`** (line 835) +
   residue-fallback floor `500 → MIN_JUNCTION_AREA_M2` (line 1075).
   Observed target minimums: SPJC 2030 m², SPLP 1069 m².  Dropping
   sub-1000 m² fragments removes 4 SPJC + 5 SPLP spurious without
   touching any real target junction.

3. **Apron-symmetry cap in `_build_junction_constructive`**
   (~line 1896).  Per user (2026-04-22) — when an arc between
   consecutive rect corners crosses an apron-exposed side of the
   cluster, the apt.dat pav boundary balloons into apron
   territory because aprons are part of `pav_union`.  Compute
   `rect_bounded_radius = max(dist from cluster centroid to any
   rect corner)`.  For arcs where chord length >
   `rect_bounded_radius × 1.2` (apron-exposed), drop arc vertices
   whose distance to centroid exceeds `radius × 1.3`.  Factor
   tuned via sweep; trigger-only variant avoids regressing the
   normal rect-to-rect arcs.  Result: spurious 21→18 at SPJC
   with matched count preserved at 31.

### Audit findings

- **Rect corners sit on pav boundary correctly**: 65% of SPJC
  rect corners and 48% of SPLP corners are ≤0.5 m from
  `pav_union.boundary`.  Outliers >10 m are all expected —
  walking-exit stubs (A/F) extending onto runway apron (rule-
  sanctioned widening stubs) and a handful of primary-parallel
  corners at irregular pav edges.  The constructive emitter's
  20 m projection filter accepts ≥94%.

- **Target junction inventory**:

  | | count | vertex-count | area (m²) |
  |---|---|---|---|
  | SPJC | 43 | 4-35 (p50=9) | 2030-74668 (p50=6487) |
  | SPLP | 13 | 5-17 (p50=7) | 1069-7671 (p50=2805) |

  Zero target junctions have zero shared nodes; every one shares
  vertices with neighbours (308 incidences at SPJC, 81 at SPLP).

### Remaining gaps at session-7 close

1. **10 of 12 SPJC missed junctions are "no output nearby"** —
   all apron-adjacent, all have primary_parallel=2 in
   neighbours.  Root cause: target subdivides the L primary at
   apron-boundary points that have no OSM bend or junction-node
   signal.  Without the L-subdivisions, no pair of rect ends
   seeds these junction clusters.  Same root cause as the 9
   missed primary_parallels (open item #3 from session 6).

2. **avgIoU 0.54 at SPJC** — emitted junctions have approximately
   right position/area but shape drifts.  Likely from 4-arc-
   vertex limit being too aggressive on large junctions (p90
   target = 17 vertices) and/or from the apron-symmetry cap
   cutting shapes on the tight side.

3. **Runway-stub trapezoid** at line 987 is still 4-vertex-only;
   target runway-stub junctions have 6-7 vertices (boundary-arc
   walk between stub corners + runway projections).  Not wired
   to use the constructive builder.

### Next-session candidates

1. **Split primary parallels at apt.dat pav-boundary vertices
   where an apron meets the rect long edge.**  Without enabling
   `EMIT_APRONS`, scan `apt_pav_vertices` for points within 20 m
   of a primary-parallel rect edge that also sit on the
   apt.dat boundary between pavement polygons; add those as
   split points.  Unblocks the 10 apron-adjacent missed
   junctions + recovers ~9 missed primaries.

2. **Apply the constructive builder to runway-stub trapezoids**
   so they get boundary-arc walks (~2-3 extra vertices each).

3. **Re-enable `EMIT_APRONS`**, then evaluate the 3 SPJC + 3
   SPLP missed aprons as a separate workstream.

### Current session-7 constants

- `EMIT_JUNCTIONS = True`; `EMIT_APRONS = False`.
- `MIN_JUNCTION_AREA_M2 = 1000.0` (was 80).
- `CLUSTER_MERGE_DIST_M = 40.0` (unchanged).
- Constructive builder: `ARC_CAP_FACTOR = 1.3`,
  `ARC_CHORD_TRIGGER = 1.2`.
- Arc-vertex limit: 4 per arc (session-6 user rule).

---

## Session 6 (2026-04-21 / 2026-04-22) — rect correctness iteration

### Changes landed

1. **Aprons also suppressed** (`EMIT_APRONS = False`).  Junction +
   apron are "treated the same" per user — focus is now purely
   taxiway rect correctness.

2. **Unrefed taxi ways as connector junction points**
   (`_find_junction_points`): short unrefed taxi ways (20-200 m)
   that share a node with a refed taxi way add a `_conn` ref at
   that shared node.  Exposed the six V↔U chart-level
   intersections that OSM marks only via unrefed bridge ways.
   V primary south now splits into 3 rects matching target.

3. **Sub-ref BEND_CLUSTER_M=30 m** (was 100 m).  Exposes the
   straight middle run inside 90° sub-ref curves — V5 now
   emits the (838,-1520) L=72 W=49 rect matching target.

4. **`_perpendicular_half_at`** returns `(left+right)/2` instead
   of `min(left,right)`.  Rect width reflects FULL pavement
   strip even when the OSM axis is off-center on the pavement.
   Corner snap pulls the 4 rect corners onto the pav boundary,
   centering the rect.  V1 W: 30 → 46 (target 47).

5. **Pavement-width midpoint cluster check** in
   `_split_centerlines_at_points` replaces fixed
   `CLOSE_INTERSECTION_M=200 m`.  Two consecutive cut-params
   merge into one junction region only if midpoint half-width
   > `narrow_hw × 1.2` (uses same `(left+right)/2` probe as
   narrow_hw).  Restored Q/R 3-rect splits while preserving
   V3's OSM-fragmented 176 m junction on V.

6. **Cross-connector end-segment margin** 30 % each side (40 %
   rect) for first/last emitted segment when ref ∈ {Q, R, X}.
   Middle segments stay 15 %.  R SW L 118→66 (target 67).
   Q SW L 98→56, Q NE L 97→55, R NE L 108→61.

7. **Non-perpendicular diagonal rule** (25 < perp_diff < 65
   vs runway): rect = 35 % of gap, center biased 25 % of gap
   toward the runway-facing endpoint.  V3 L 159→78 (target 70),
   center (438,-720) vs target (441,-723).

8. **Letter-only non-parallel stubs (B, C, D, E, G)** skip
   bend-split when chord/path ≥ 0.95.  Treat the full OSM way
   as ONE centerline and apply diagonal 35 % rule.  Stub
   dedup extended to keep longest per ref (including
   letter-only stubs).  E went from 4 pieces to 1 (L=229,
   target 181).  B, C, D, G similarly.

9. **40 m post-margin floor** on emitted rect centerlines.
   Drops junction-approach tail fragments (e.g. R switchback
   at the V/Q/R triple junction).

10. **Candidates pre-filter**: segments that cannot emit
    (post-margin < 40 m under any plausible margin) are
    excluded from the end-segment indexing so idx==0 / last
    is stable for the cross-connector 30 % rule.

11. **Primary-parallel runway-end stubs**
    (`_emit_primary_parallel_runway_stubs`): A / F / L OSM
    primary taxis extend their polyline ONTO the runway
    polygon at SPJC (the OSM endpoint is INSIDE the runway).
    Target marks the transition as a wider-than-normal STUB
    (A L=79 W=73, F L=93 W=78 — the runway-apron ramp).
    Post-pass: for each primary parallel merged polyline
    whose endpoint is inside the runway polygon (or within
    10 m), walk the path toward the interior until the
    vertex distance to runway boundary exceeds 80 m ("exit
    point"), emit an 80 m STUB rect centered at the exit
    along the local path direction, width = 2 × narrow-hw
    probe at the exit.  Guards: overlap > 20 % with existing
    rects, or same-ref rect within 30 m buffer → skip.
    Recovered 15/15 stubs at SPJC (A, F added; L already
    emitted as "L7" via the main pipeline).

12. **Diagonal range widened to 20-75** (was 25-65) for
    `_rect_margin_frac_for` perp_diff check.  B / C / E / G
    at SPJC measure perp_diff ≈ 69° — just outside the old
    25-65 window.  Per user (2026-04-21): letter-only
    non-parallel taxis that "connect to the runway at an
    angle" are in the same diagonal-stub family as V3 and
    should get the same 35 % retained + 25 % gap bias
    toward the runway-facing endpoint.  Result: B
    (1214,392) L=116 W=43, C (1496,-197) L=114 W=50,
    E (1719,-661) L=118 W=50, G (1821,-873) L=86 W=47 —
    bias shifted each rect toward its segment's
    runway-facing endpoint.  NOTE: each B/C/E/G centerline
    is still bend-split at its A / L crossing partway
    along, so the "gap" the rule applies to is the
    POST-SPLIT segment (~320 m for B, not full 575 m OSM
    path).  Target L=191 for B would require the FULL
    merged path as gap (atomic emission) — attempted, but
    rects then clip against pav_union-minus-runway and
    regress match count 50→46.  Current compromise: rule
    is applied on bend-split segment, giving ~100-120 m
    rects with correct bias direction.

13. **Diagonal bias 25 % → 20 %** per user (2026-04-21):
    adjust "slightly away from the runway" — biasing the
    rect center 25 % of the gap toward the runway-facing
    endpoint was too aggressive; 20 % leaves more room
    for the runway-ramp widening zone.

14. **Diagonal rule scoped to REFED non-parallel taxis**:
    `_rect_margin_frac_for` now early-returns 0.15 for
    unrefed ways (`not ref`) and for PARALLEL_REFS / Q,R,X.
    SPLP's primary taxis are UNREFED at 19° off runway
    (perp_diff = 71° inside the 20-75 diagonal window) but
    the full 35 % + bias rule shrinks them to junk and
    they should stay at 15 % margin like any primary
    parallel.  Only refed letter-only non-parallel stubs
    (B/C/E/G) and digit-bearing sub-refs (V1, V3, V5, L1,
    L3…) now get the diagonal treatment.

15. **Runway-end stubs extended to unrefed long taxi ways**
    (`_emit_primary_parallel_runway_stubs`): accepts
    unrefed OSM ways whose path length ≥ 800 m and
    processes individual ways (not linemerged) so that
    shared junction endpoints like SPLP's (226,1182) stay
    visible as endpoints.  Threshold for "outside but near
    runway" tuned to 135 m (SPLP main taxi NE end is
    127 m; SPJC L's internal endpoints at 148–150 m are
    excluded).  Same-ref 30 m buffer dedup only applies
    when `ref` is non-empty.  Stub centers within 50 m of
    each other coalesce.  Recovered SPLP A-equivalent
    stubs at (-449,-1084), (-507,-1176), and (226,1183)
    — NE runway-facing ramp — without regressing SPJC.

16. **Short unrefed runway-connecting stubs emitted
    atomically** in `_extract_osm_taxi_centerlines`: when
    ref is empty, way length < 300 m, and one endpoint is
    within 30 m of a runway centerline, skip bend-splitting
    and emit the simplified polyline as ONE centerline.
    Recovers SPLP way -696731 (165 m chord 144 m, ratio
    0.87) which was being bend-split into 16 m + 20 m
    fragments that got dropped by the 40 m post-margin
    filter.  Result: target (-293,-641) L=53 W=22 is now
    matched by output (-297,-635) L=90 W=32.

17. **Reject stubs whose rect touches the runway polygon**
    in the stub runway-connection filter (line ~580).  SPLP
    has a short southernmost diagonal stub whose rect has
    one corner exactly on the runway boundary — user
    wanted "only a single stub at the south end" so this
    runway-touching stub is dropped.  Keeps the
    neighbouring "good" stub up and out.

18. **Short parallel-oriented refless rects classified as
    PRIMARY_PARALLEL** (not stub): `_classify_role` refless
    branch lowered `length >= 80` to `length >= 50` for
    parallel rects.  SPLP's 66–80 m parallel slice between
    two diagonal stubs in the south chain now emits as a
    primary_parallel.

19. **Per-segment margin for short parallel slices of a
    curving unrefed taxi**: inside
    `_split_centerlines_at_points` the sub-segment's own
    bearing is compared to the nearest runway when ref is
    empty and gap < 150 m.  If the slice is nearly parallel
    (perp_diff ≥ 75°), apply 22 % margin each side (56 %
    retained) so the resulting primary is roughly half its
    default 70 % length — avoids visual overlap with
    neighbouring diagonal stubs.  Example SPLP:
    (-500,-1054) L=81 W=37 → (-506,-1059) L=54 W=26.

20. **Diagonal rule extended to short unrefed diagonals**
    (< 250 m): SPLP's (-449,-1084) L=109 W=29 diagonal stub
    (perp_diff=43°) now gets 35 % retained + 20 % bias
    toward runway → (-435,-1110) L=54 W=29.  Long unrefed
    centerlines (primary taxis) still return 0.15 margin.

21. **Iterative symmetric axis trim** in
    `_rect_from_axis_extended` — when corner-snap produces an
    asymmetric rect (trapezoid with unequal end-widths or
    parallelogram with unequal long sides), trim BOTH ends
    by 2.5 % of axis length (5 % total) per iteration and
    rebuild.  Uses RATIO-based tolerance (width Δ /
    max_width > 20 % or length Δ / max_length > 10 %) so
    short stubs with modest m-diffs aren't over-trimmed.
    Max 15 iterations (up to 75 % shrink).  Result: SPLP
    trapezoidal stubs now read as proper rectangles with
    near-equal sides.

22. **30 m runway-buffer trim for perpendicular taxis**:
    before cut-param splitting, trim each perpendicular-to-
    runway centerline (perp_diff < 25°) at a 30 m-buffered
    runway polygon.  Captures the runway-apron widening
    zone extending OUTSIDE the runway polygon proper.
    Example SPLP: (58,467) L 95 → 34 (target 39).  Extended
    later (item 24) to all non-parallel centerlines.

23. **Parallel-centerline buffer trim** for non-parallel
    centerlines (perp_diff < 75°): trim at a 15 m-buffered
    polygon around every parallel-to-runway centerline
    (length ≥ 200 m, perp_diff > 75°).  The 15 m buffer
    approximates a primary's half-width, so the trimmed
    perpendicular/diagonal stub starts at the primary's
    physical EDGE not its axis.  Fixes the "off-center
    toward taxiway" complaint: SPLP perpendicular stub
    centers now match target exactly — (50,476) at
    target (52,473) ✓, (-140,-123) at exact target ✓,
    (-301,-638) at target (-293,-641) ✓.

24. **Split STUB_EXIT_D_M between refed and unrefed**:
    SPJC A/F (refed) stay at 80 m threshold — target sits
    at the loop-ramp apex (d_rwy ≈ 100 m) which lands on a
    path vertex.  SPLP unrefed uses vertex exit at 80 m
    PLUS interpolation back to d_rwy = 75 m (target sits
    BETWEEN vertices mid-curve).  Result: SPLP walking
    runway-curve stub (-494,-1205) now at d=1 m from
    target (-493,-1204) — exact match (was d=22 m).

### Results at SPJC (tol=0.5 m, 2026-04-22)

| role | n_t | n_o | matched | spurious | avgIoU |
|---|---|---|---|---|---|
| runway | 2 | 2 | 2 | 0 | 0.95 |
| primary_parallel | 30 | 21 | 21 | 0 | 0.67 |
| secondary_parallel | 4 | 5 | 4 | 1 | 0.73 |
| cross_connector | 6 | 6 | 6 | 0 | **0.76** |
| stub | 15 | 15 | **15** ✓ | 0 | 0.68 |
| terminal | 2 | 2 | 2 | 0 | 0.29 |
| apron | 3 | 0 | 0 | 0 | — (disabled) |
| junction | 43 | 0 | 0 | 0 | — (disabled) |
| **TOTALS** | **105** | **51** | **50** | **1** | |

**V-family final state:**

| | Target (L W) | Output (L W) |
|---|---|---|
| V1 stub | 63 47 | 83 46 ✓ |
| V2 stub | 66 56 | 81 55 ✓ |
| V3 stub | 70 58 | 78 58 ✓ (biased + 35 %) |
| V5 stub | 66 52 | 74 49 ✓ |
| V primary S #1 | 282 51 | 305 42 |
| V primary S-mid | 104 52 | 147 42 ✓ |
| V primary near-V5 | 130 50 | 147 54 ✓ |

### Results at SPLP (tol=0.5 m, 2026-04-22)

| role | n_t | n_o | matched | spurious | avgIoU |
|---|---|---|---|---|---|
| runway | 1 | 1 | 1 | 0 | 0.87 |
| primary_parallel | 5 | 6 | 5 | 1 | 0.60 |
| secondary_parallel | 0 | 1 | 0 | 1 | — |
| cross_connector | 2 | 2 | 2 | 0 | 0.50 |
| stub | 6 | 6 | 5 | 1 | **0.58** |
| apron | 3 | 0 | 0 | 0 | — (disabled) |
| junction | 13 | 0 | 0 | 0 | — (disabled) |
| **TOTALS** | **30** | **16** | **13** | **3** | |

SPLP (up from 10 → 13 matched this session).

**SPLP perpendicular-stub centers now match target exactly:**

| Stub | Output | Target | Δ pos |
|---|---|---|---|
| N cross (58,467) | L=38 @ (50,476) | L=39 @ (52,473) | 3 m |
| Mid (-140,-123) | L=58 @ (-140,-123) | L=31 @ (-140,-123) | 0 m |
| S-mid (-301,-638) | L=50 @ (-301,-638) | L=53 @ (-293,-641) | 9 m |
| Curving primary stub | L=88 @ (-494,-1205) | L=109 @ (-493,-1204) | **1 m** |

### Open items for next session

1. ~~A / F runway-end stubs~~ — **RESOLVED** via change #11
   above.  A stub (687,1563) L=80 W=72 ≈ target (688,1562)
   L=79 W=73.  F stub (2250,-1662) L=88 W=78 ≈ target
   (2243,-1666) L=93 W=78.

2. **Cross-connector middle rect IoU** — Q/R middle rects
   match position/length but avgIoU 0.73 could improve with
   better corner snapping.

3. **Primary-parallel missed = 9** — V primary north / V
   primary mid (big one at (-64,534) L=915) vs my (-26,459)
   L=750 — off by ~165 m.  Likely the NE-most V segment
   splits at a chart-level point with no OSM signal.

4. **Re-enable junctions + aprons** once rect count +
   positions are fully approved.  The pav-width midpoint
   check + end-segment margin logic should NOT conflict
   with junction constructive builds.

5. **Spurious secondary_parallel (1 at SPJC)** — M stub at
   (943,160) L=78 W=18 is a degenerate rect with W=18 m,
   far under normal taxi width.  Likely from a fragmented
   M centerline at the V-Q-R meeting area.  Add a minimum-
   width filter (e.g. reject rect if W < 30 m for
   secondary/primary parallels).

6. **SPLP diagonal stub (-432,-1115) too short** — target
   (-447,-1087) L=86 vs mine L=40.  Target extends ~30 m
   BEYOND the OSM centerline's extent into apt.dat
   pavement.  Would need a "extend-along-pavement" pass
   scoped to unrefed diagonals (risky for SPJC sub-refs
   whose OSM length is definitive).  Left note in code;
   also an attempt at pavement-edge-aware end-trim was
   reverted (dropped a valid stub by interacting with
   diagonal centerline assembly).

### Current constants (session 6)

- `SIGNIFICANT_BEND_DEG = 5.0`
- `BEND_CLUSTER_M = 100.0` (primaries); `30.0` (sub-refs)
- `CLOSE_INTERSECTION_M` (legacy, still referenced) — replaced
  at runtime by pav-width midpoint check with
  `WIDEN_FACTOR=1.2`, `MIN_ALWAYS_MERGE=25 m`,
  `MAX_CLUSTER_SPAN_M=400 m`.
- `GAP_MARGIN_FRAC = 0.15` (default), `0.30` for cross-
  connector end-segments, `0.35` for diagonal stubs
  (20 < perp_diff < 75) with 20 % gap bias toward runway
  (user 2026-04-21: reduced from 25 %).
- `40 m` minimum post-margin centerline length.
- `RWY_JUNCTION_BUFFER_M = 30` (perp-taxi runway buffer).
- `PARALLEL_BUFFER_M = 15`, `PARALLEL_MIN_LEN_M = 200`
  (non-parallel taxi trim at parallel-centerline buffer).
- `STUB_EXIT_D_M = 80` (walking exit threshold);
  `STUB_INTERP_TARGET_D_UNREFED = 75` (unrefed stub
  interpolates back to this d for exact mid-curve fit).
- `ASYM_WIDTH_RATIO_TOL = 0.20`,
  `ASYM_LENGTH_RATIO_TOL = 0.10` (iterative symmetric trim).
- `EMIT_JUNCTIONS = False`, `EMIT_APRONS = False`.

---

## Session 5 (2026-04-20) — legacy

## Session 5 simplified algorithm (2026-04-21)

After multiple width-based iterations (perpendicular ray-cast
probe, narrow-corridor selection, width-threshold clustering),
user approved a **completely simpler approach**: abandon all
pavement-width analysis for rect extraction.  Rects defined
solely by:

  1. **OSM multi-ref intersection points** (`junction_points`).
  2. **Sharp curves** (bend clusters ≥ `SIGNIFICANT_BEND_DEG`).
  3. **Centerline endpoints** (parent-taxi / runway / apron).

Between consecutive breaks, emit **70% centered rect**
(15% margin on each junction-facing end).

**Code simplification**: `_split_centerlines_at_points` went
from ~200 lines to ~70.  Removed `_sub_ref_narrow_corridor`
width-profile logic entirely.  Sub-refs now go through the
SAME pipeline as primaries (both bend-split + 70% trim).

**V-family state (TGT vs OUT at SPJC, after simplification):**

| | Target | Output |
|---|---|---|
| V1 stub | 98 m | 86 m ✓ close |
| V2 stub | 118 m | 78 m (short) |
| V3 stub | 116 m | 159 m (over) |
| V5 stub | 103 m | 48 m (short, structural) |
| V primary count | 5 | 4 (missing south split) |

**Structural limits identified:**

- **V5** target center (839, -1525) is OUTSIDE OSM V5 way's
  extent (endpoints y=-1519 to -1435).  Target covers
  pavement OSM doesn't tag.  Cannot match without chart
  knowledge or apt.dat pavement-polygon analysis.
- **V primary south split** between target -803810 (162m)
  and -803811 (135m): no OSM intersection in that 639m gap.
  Chart-level decision.
- Aprons emit a big "V-Q-R meeting area" apron because with
  junctions disabled, leftover pavement >25000 m² is
  currently classified as apron.

## Session 5 constants (for reference)

- `SIGNIFICANT_BEND_DEG = 5.0`
- `BEND_CLUSTER_M = 100.0` (groups small bends into one curve)
- `CLOSE_INTERSECTION_M = 200.0` (merge close intersections)
- `SAME_INTERSECTION_M = 60.0` (unconditional cluster below)
- `GAP_MARGIN_FRAC = 0.15` (70% rect centered between breaks)

## Session 5 next priorities

1. **Big apron at V-Q-R meeting area**: should be junction
   territory, not apron.  Fix: raise MIN_APRON_AREA or
   exclude apron emission in "junction-like" residue areas
   when `EMIT_JUNCTIONS=False`.
2. **V3 / V5 / V primary south split**: may be unachievable
   without external signals (chart or apt.dat pavement
   polygon analysis).  Accept structural limits for now.
3. **Re-enable junctions** (`EMIT_JUNCTIONS = True`) once V
   rects are visually approved.  Junctions emit via
   `_build_junction_constructive` (corner + apt.dat-arc
   polygons, max 4 arc vertices between corners).
4. **Extend rules to L, F, A, M, U** after V approval.
5. **Later**: apt.dat taxi route 1202 edges as supplemental
   intersection signal (could help V primary south split
   and V5/V3 stub placement).

## Session 5 structural observations

- Sub-refs often have OSM ways that span BOTH the stub's
  narrow corridor AND the junction-widening/curve on either
  side.  Target's sub-ref rects cover a CENTERED portion
  that's narrower than the full OSM way.
- OSM is fragmentary: V2 has 3 disjoint ways, V3 has 2,
  etc.  Target consolidates; my sub-ref dedup picks the
  longest rect per ref.
- Target splits some primaries at positions with NO OSM
  signal (chart convention).  These are structural-limit
  shapes — ~5-10% of total.

---

## Session 4 (2026-04-19) — legacy

| | SPJC | SPLP | Total |
|---|---|---|---|
| matched | 82/105 (78 %) | 22/30 (73 %) | **104/135 (77 %)** |

Match count up +2 vs session-start 102/135 (76 %), **and vertex-
exact (tol=1 m) accuracy jumped dramatically**, which is the
metric that matters for "near perfect match":

| role | v_tgt@1m start | end | gain |
|---|---|---|---|
| SPJC junction | 7.4 % | **30.4 %** | 4× |
| SPJC cross_connector | 8.3 % | **16.7 %** | 2× |
| SPJC primary_parallel | 1.0 % | **4.8 %** | 5× |
| SPJC apron | 14.7 % | 13.9 % | ≈ |
| SPLP junction | 0.0 % | **30.5 %** | new |
| SPLP cross_connector | 0.0 % | **25.0 %** | new |
| SPLP apron (tol=5m) | 4.8 % | **38.1 %** | 8× |

User goal is 95 % match (129/135) AND near-perfect vertex
alignment.  Vertex goal now well in hand; match count needs +25.

**Review artifacts (end of session):**
- `/tmp/SPJC_auto.osm` (146 shapes: 2 runway, 22 stub, 38 primary,
  4 secondary, 10 cross, 2 terminal, 5 apron, 63 junction)
- `/tmp/SPLP_auto.osm` (43 shapes: 1 runway, 11 stub, 6 primary,
  1 secondary, 2 cross, 21 junction, 1 apron)

## Session 4 changes (2026-04-19)

Implemented per user direction (6-step iterate pattern):

1. **Rect corners snap to apt.dat VERTICES** (not edge points):
   `_snap_corners_to_pavement` now does two-stage snap — first
   tries nearest apt.dat vertex within 8 m (the authoritative
   coord set the target also uses), falls back to nearest boundary
   point within 15 m.  Guards against coincident-corner snap that
   would produce degenerate rects.

2. **Global shared-vertex enforcement**
   (`_enforce_shared_vertices`): after all shapes emitted, clusters
   every output vertex across every shape within 1.5 m and replaces
   each with the cluster centroid.  Produces exact shared vertices
   between adjacent shapes (rule 16).

3. **Invariant validator** (`_validate_shared_vertex_invariant`):
   raises RuntimeError if any two shape vertices are 0.01 < d ≤ 1.5
   m apart.  Forces correctness — build FAILS on violation.

4. **Degenerate rect rejection**: `_rect_from_axis_extended`
   returns None when any two post-snap corners are within 1 m of
   each other (rule 7: corners on pav boundary means 4 distinct
   apt.dat vertices).

5. **Refless classifier** (`_classify_role` angle-only branch):
   now uses db + length + dist-to-runway to distinguish primary /
   secondary / cross_connector / stub at airports without OSM
   refs (SPLP).  Long perpendicular pieces far from runway →
   cross_connector; short perpendicular → stub.

6. **Target snapper** (`tools/snap_target_to_apt_dat.py`): added
   vertex-preferred snap.  Target vertices within 8 m of an
   apt.dat vertex snap exactly to it — aligning with the
   builder's choice of the same vertex set.

7. **Final-coincidence sweep** in `_snap_corners_to_pavement`:
   after both vertex-snap and edge-snap, checks all pairs of
   final rect corners; if any two ended within 1 m (vertex-snap
   to one + edge-snap fallback to the same boundary region),
   reverts the farther-from-original to its pre-snap coord.
   **Recovered SPLP 13 → 22 matches.**  Rejection rate 6/20 →
   0/20.

8. **Junction dedup relaxed**: overlap threshold 50 % → 80 %.
   Keeps more near-duplicate junction candidates (they often
   correspond to distinct target junctions in dense cluster
   regions).  Net +1 match at SPJC.

9. **L7 rule confirmed working** via existing `_refine_roles`:
   demotes short (< 150 m) perpendicular (>= 40° off runway)
   segments of parallel refs to stub.  Corner junction emerges
   automatically from rect-end seeding.  Target's L7 itself is
   a user-added ref not in OSM, but the underlying stub geometry
   is generated under ref=L.

10. **Principled collinear merge** (`_merge_collinear_rects_principled`)
    with pav-width-at-joint uniformity check — added as scaffolding
    but disabled.  Regressed SPJC match count by 2 because target
    subdivides at some uniform-width joints the merger combined.
    Kept in source for next iteration with better signals
    (e.g. junction-membership check).

## NEXT SESSION priorities

1. **Target-specific missed-junction analysis**: 15 SPJC + 2 SPLP
   junctions still missed.  Sizes range 2000–72000 m².  For
   each, determine what apt.dat feature it corresponds to and
   why my algorithm doesn't emit there (or emits a different
   shape).  The big-area misses (> 20000 m²) are likely terminal-
   area junctions that my gap-fill considers apron; the small
   (< 3000 m²) ones are rect-end-corner junctions my constructive
   build produces with different geometry.

2. **Over-emission**: 62 spurious outputs (mostly junctions) —
   my per-cluster constructive builds + gap-fill residue together
   emit more junctions than target has.  Need a consolidation
   pass that merges nearby junctions of similar area, OR drops
   gap-fill residue that's already covered by a constructive
   junction.

3. **Apron undercount**: SPLP has 1/3 apron matched (2 missed).
   My algorithm treats most SPLP residue as junctions (area <
   25000 m² threshold).  The 2 missed SPLP aprons may be smaller
   than threshold — consider lowering apron threshold for small
   airports or detecting apron shape (concavity, MRR ratio)
   instead of area alone.

4. **Principled collinear merge v2** (currently scaffolded): the
   width-uniformity check wasn't enough — try combining with
   junction-membership check (if the joint sits inside a
   current junction polygon, do NOT merge; otherwise merge).

5. **SPLP limitation noted**: apt.dat labels all SPLP taxiways as
   "A" — no per-taxi refs available from either OSM or apt.dat
   for this airport.  Classification must rely fully on geometry.

## Legacy: earlier session priorities (2026-04-18 end-of-day)

1. **Too many rects between junctions** — rule violation.  Target
   has V=5, mine has V=11.  L=11 target, mine 10-14 depending on
   settings.  My `_split_centerlines_at_points` splits at every
   junction_point within 25 m of a centerline, including OSM
   multi-ref nodes that the user consolidated in the target.
   **Fix approach:** after splitting, check if adjacent same-ref
   segments are collinear (small angle change) AND have no
   actual widening between them — if so, merge them back.  The
   target's "between junctions" = between WIDENING points, not
   between every OSM crossing.

2. **Junctions should be larger** — my constructive corner+arc
   polygons are bounded by the rect-end corners (which are trimmed
   inward at `widen_factor=1.01 × narrow_hw`).  Target junctions
   cover MORE area than my polygons.  Possible causes:
   - Arc-vertex cap of 4 may truncate the true pav-boundary curve.
   - Rect corners snap within 15 m of boundary but some land
     short, leaving uncovered strips between corner and actual
     pav edge.
   - Gap-fill threshold 500 m² may be dropping legitimate
     junction area as noise.
   **Fix approach:** before the max-4-arc sub-sample, check the
   total arc length; if > some threshold, keep more vertices.
   Also: walk the rect-end corners along the pav boundary past
   the trim point (outward) to capture the widening zone in the
   junction.

## Rules implemented this session

1. **Rect construction**
   - Only straight segments, cut at curves + cross-ref
     intersections (via `_split_centerlines_at_points`).
   - Width = narrowest probe (p10 half-width).
   - 1 % widening → cut (`widen_factor=1.01 × narrow_hw`).
   - Corners snap onto pavement boundary (15 m safety radius).

2. **Junction construction** (`_build_junction_constructive`)
   - Vertices = incoming rect corners.
   - Arcs between consecutive corners trace apt.dat pav
     vertices on the boundary, sub-sampled to max 4.
   - Same-rect corners (rect's 2 short-end corners) connect
     directly (no arc).
   - Clipped to local pav disc + subtracted from rect-union and
     terminal.
   - Seeded from ALL rect-ends + OSM multi-ref clusters;
     clusters merged within 40 m; polygons deduped when overlap
     ≥ 50 %.

3. **Runway-taxiway connections**
   - Runway has vertices inserted at each widening-stub's
     projection (`_insert_points_on_boundary`).
   - Uniform-width stub (no trim at runway end): direct connect,
     NO junction emitted.
   - Widening stub (trim pulled back from runway edge): 4-corner
     trapezoid junction between stub corners and their runway-
     boundary projections.

4. **Rule 15 (no gaps)**
   - After rects + junctions + aprons + terminals emitted, any
     residue ≥ 500 m² but < 25 000 m² emits as additional
     junction (simplified at 3 m).

5. **Shared-vertex invariant**
   - Runway gets vertices inserted at every stub-widening
     projection.
   - Rect corners snapped to pav boundary serve as shared
     anchors with adjacent junctions.

Session 3 changes (2026-04-18, same-day resume):
- **Rect half-width = narrowest probe (p10)** instead of p90.  Per
  user's authoritative rule: rect covers the section of straight
  pavement at its narrowest width; wider portions belong to the
  adjacent junction.  `_natural_half_width` now returns a third
  `narrow_hw` value; `_build_taxi_rects` uses it for width.
- **Rect corner snap = 15 m** onto pavement boundary (was 5 m).
  Per user rule ("rect corners should always be on a pavement
  boundary"), unbounded snap was too aggressive (warped rects at
  degenerate axis placements and produced self-intersecting
  unions); 15 m covers all non-pathological cases.  Corners now
  land on apt.dat boundary in practice.
- **Junctions = ordered rect corners, direct edges, no arc
  vertices.**  Per user's explicit direction: "implement junctions
  that only connect rect corners directly with no intermediate
  vertices, I want to see the result in JOSM."  The entire residue
  + disc-seed + snap pipeline was removed.  For each cluster of
  OSM junction points (merged at 200 m), the builder gathers outer
  corner pairs of every rect whose axis-end is within the cluster,
  orders them angularly around the cluster centroid, and emits a
  simple corner polygon.  No pav-arc tracing, no buffer smoothing.
  Aprons are now the pav-residue left after rects + junctions,
  giving cleaner apron shapes (IoU 0.57 → 0.65).

Session 2 wins retained:
- **Stubs using ORIGINAL un-simplified polyline** instead of RDP-2-coord
  chord: SPJC stubs 7/15 → 13/15.  OSM V5 had 26 curve nodes but
  RDP 1m collapsed to a chord cutting across the curve, clipping
  most of it out of pavement.  Using the raw polyline keeps the
  curve and lets downstream trim + rect-build work.
- **Clip to full pavement union** (not pav - runway) so stubs
  extending to runway edge keep their full physical length.

## NEXT SESSION: user JOSM review + iterate on junction shape

Resume after user has opened `/tmp/SPJC_auto.osm` and
`/tmp/SPLP_auto.osm` in JOSM alongside the target files and
evaluated the direct-connect corner polygons.  Likely feedback
topics:

1. **Which junction clusters are too big?**  The 200 m cluster
   merge may have combined clusters the user expects separate
   (e.g. two adjacent taxi-taxi crossings 150 m apart become one
   polygon with 8 corners, when user wants them as 2 separate
   4-corner junctions).  Adjust CLUSTER_MERGE_DIST_M based on
   feedback.
2. **Which rect corners are still missing from pav boundary?**
   15 m snap covers most cases; the remaining off-boundary
   corners suggest rects with axis-placement issues (trim stopped
   early, or axis too far from centerline).  User may point at
   specific refs in JOSM for diagnosis.
3. **Where does the direct-edge polygon cut across pavement
   incorrectly?**  E.g. a cluster with 3 rects forms a triangle
   that crosses through a rect's interior.  User may want a
   subtraction step or a pav-clip.
4. **Apron fit.**  Junctions eat into what should be apron
   territory, or vice versa.

`_build_junction_constructive` (with the apt.dat pav-arc walk)
is retained as scaffolding for the next iteration once the user
has given concrete shape feedback on direct-edge corners.

Memory: `~/.claude/projects/-Users-noah-Ortho4XP-shred86/memory/feedback_shape_rules.md`
has the user's authoritative spec.

Memory: `~/.claude/projects/-Users-noah-Ortho4XP-shred86/memory/feedback_shape_rules.md`
has the full spec.

Elevation is Phase 2, not started.

## Algorithm overview

1. **Runway:** apt.dat row 100 + blast pads → 4-vertex rect.  IoU
   0.96 at SPJC.
2. **Terminals:** OSM `aeroway=terminal` (ways + multipolygon
   relations); expanded to the containing apt.dat pavement polygon
   (the "pad").  IoU 0.29.
3. **Taxi rects:** OSM taxi centerlines per ref, linemerge,
   same-ref gap-bridge (up to 120 m), RDP at 1 m, split at bends.
   Parallel refs (A/F/L/V/M/U + Q/R/X) split at bends; other refs
   emit ONE rect per merged polyline.  Rect width = 2× 90th-percentile
   probe distance-to-boundary (full pavement span).  Axis trimmed
   at widening points (widen_factor = 1.05 × natural_hw).
4. **Classification:**
   * `A/F/L/V` → primary_parallel; `M/U` → secondary_parallel.
   * `Q/R/X` → cross_connector when Δ-to-runway > 40°.
   * `letter+digit` and plain-letter non-parallels (B/C/D/E/G) → stub.
   * Post-pass demotes short perpendicular (>40°, <150 m) segments
     of parallel refs to stub (captures "short A-runway connector").
5. **Stub filter:** drop stubs whose raw OSM endpoint is > 80 m
   from any runway (user: "stubs must connect to a runway").
6. **Junction seeding:** at every OSM multi-ref node cluster OR
   geometric crossing between different-ref centerlines, seed a
   50 m disc clipped to pavement.
7. **Junction + apron emission:** pavement residue (pav - rects -
   terminals) plus junction discs → merge at 40 m → junction if
   < 25000 m², apron if ≥ 25000 m².
8. **Junction snap-to-rect-corners:** polygon vertices within 3 m
   of a rect corner snap to that corner (shared-vertex invariant).

## Per-role scores at 5 m tolerance (SPJC, 106 target shapes)

| Role | matched / target | avgIoU |
|---|---|---|
| runway | 2/2 | 1.00 |
| apron | 3/3 | 0.64 |
| terminal | 2/2 | 0.29 |
| primary_parallel | 19/30 | 0.56 |
| secondary_parallel | 3/4 | 0.41 |
| cross_connector | 4/6 | 0.40 |
| stub | 6/15 | 0.29 |
| junction | 34/44 | 0.34 |

## Missing ~30 % gap — what it'd take

1. **Primary parallels (11 missed, 12 spurious):** my RDP bends
   don't coincide with user's segment cut points.  Target's 30
   primary segments for A/F/L/V reflect chart-level judgement I
   can't extract from apt.dat + OSM alone.
2. **Stubs (9 missed, 10 spurious):** specifically L1/L7/V5
   (sometimes missing from OSM, sometimes filtered); and the
   short parallel-to-stub demotions don't always pick the right
   segment.
3. **Junctions (10 missed, 24 spurious):** junction polygon shapes
   differ from target (IoU 0.34) — target polygons trace specific
   vertex paths around rect corners + pavement edges; mine are
   residue-hulls.

The last 30 % is hard without external signals (chart / imagery /
detailed user rules per ref).

## Workflow

1. Edit target freehand: `tests/fixtures/SPJC_target.osm` (or
   SPLP_target.osm).
2. Re-snap: `python3 tools/snap_target_to_apt_dat.py SPJC
   tests/fixtures/SPJC_target.osm --force`
3. Re-run: `python3 tools/build_target_osm.py SPJC`
4. Compare: `python3 tools/compare_target.py
   tests/fixtures/SPJC_target_snapped.osm /tmp/SPJC_auto.osm --tol 5`

The snap tool won't overwrite `_snapped.osm` without `--force` (to
avoid losing user's manual edits to snapped files).

---

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
