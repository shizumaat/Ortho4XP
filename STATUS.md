# Auto-Patch Refactor — Status

**Current direction (as of the latest commit):** we are returning to
the legacy `generate_airport_surface_patches()` in
`src/O4_Auto_Patch.py` as the foundation, fixing its known bugs, and
refactoring it into a clean PR-ready module.  The "new pipeline"
experiment in `src/O4_Surface_Patch.py` is being retired — read the
"Why we reverted" section below for the reasoning.

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

1. **CIFP runway threshold elevations are authoritative and immutable.**
   Everything else bends to them.
2. **Sample the DEM to determine centroid elevations for buildings and
   terminals.**  Initial estimates only — may be adjusted in step 5.
3. **Build taxiways that meet slope requirements** (1.5 % longitudinal,
   FAA vertical-curve rate-of-change rule).  Taxiway elevations are
   constrained by CIFP anchors at runway crossings and by DEM-sampled
   building centroids where the taxiway enters a terminal area.
4. **Join buildings and terminals to their nearest taxiways following
   apron slope rules** (1.0 % any direction).  The apron is the
   connective tissue: it must plane-fit (or triangulate) between the
   terminal edge and the taxiway crossing point with grade ≤ 1.0 %.
5. **If an apron still cannot meet 1.0 % with the current
   building/terminal elevation, adjust the building or terminal
   elevation** to produce the smoothest achievable slope.  This is the
   only situation in which building elevations move from their
   DEM-sampled values.

The legacy `generate_airport_surface_patches` already implements a
version of this ordering in its Phase A5 elevation pipeline — see
`src/O4_Auto_Patch.py:2238`.

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

## Current state

- `USE_NEW_SURFACE_PIPELINE` flag has been removed.
- The delegation shim at the top of `generate_airport_surface_patches`
  has been removed.
- `generate_auto_patches()` no longer branches on the flag — it always
  uses `generate_patch_osm()` for the runway baseline and then calls
  `generate_airport_surface_patches()` to layer on the rest.
- The legacy phase-A-through-F implementation runs unchanged.
- `src/O4_Surface_Patch.py` still exists on disk but is no longer
  called.  Will be retired in an upcoming commit (or repurposed —
  see "Plan").
- `scripts/test_spjc_surface_patch.py` and
  `scripts/analyze_spjc_patch.py` need to be updated to call the
  legacy entry point instead of `O4_Surface_Patch.build_surface_patch`.

## Plan

Working in commits so each step has a roll-back point.

1. **Commit 1 — revert flag, restore legacy.**  Done in the same
   commit as this STATUS.md update.
2. **Commit 2 — fix the six known bugs in the legacy** (listed under
   "Bugs to fix in the legacy" below).  Each bug is a small targeted
   change; bugs land as separate commits where it makes sense.
3. **Commit 3 — refactor for PR quality.**  Split the 2 700-line
   monolithic `generate_airport_surface_patches` into named
   module-level helpers, one per phase.  Either (a) keep everything in
   `O4_Auto_Patch.py`, or (b) move the surface mesh builder to its own
   file with a name that actually distinguishes it (current candidate:
   `O4_Surface_Mesh.py`).  This decision is open — see the question
   asked in the conversation.
4. **Commit 4 — apply learnings from the experiment.**  Specifically:
   airport-centered projection for plane fits (numerical stability),
   `patch_feature` tags on every emitted way (for analyzer + JOSM
   debugging), `min_coverage`-based recursive v-subdivision for any
   apron strip emission paths.
5. **Commit 5+ — Phase E (boundary band, tunnel portals) and Phase F
   (drainage) refinement.**  These were producing useful output but
   were never finished; left in place and marked as such.
6. **Commit N — re-run the test plan** in `docs/TEST_PLAN_SPJC.md`,
   verify all 11+ checks plus the four new invariant checks (12-15).

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

### 2. `generate_patch_osm()` emits one rect per runway
The runway-only entry point used as the baseline by Phase 0 emits a
single sloped rectangle from physical end to physical end, losing all
DEM-driven undulation and ignoring displaced thresholds.
**Fix:** port the segmented runway sample chain from the experimental
`O4_Surface_Patch._build_runway_segments` — DEM samples every
`RUNWAY_SEGMENT_LENGTH` (100 m), CIFP anchors at displaced thresholds
+ physical ends, two-pass relaxation (hard 1.5 % cap + FAA
vertical-curve rate-of-change `L ≥ 30 m per 1 % of grade change`,
then re-cap).

### 3. No FAA vertical-curve rule
Same as #2 — the legacy doesn't apply
`MAX_GRADE_CHANGE_PER_M = 1/3000` between adjacent runway segments.
Lands together with the runway segment fix.

### 4. IDW elevation falloff returns 0 m beyond radius
The legacy `_idw_elev` (or its inline equivalent in the CIFP surface
model) zeros out vertices beyond `IDW_RADIUS`, causing runaway-low
drainage shapes.  **Fix:** add a nearest-anchor fallback that returns
the closest anchor's elevation when no anchor falls within radius.

### 5. Displaced threshold positions not anchored as samples
The runway segment builder needs to insert the displaced-threshold
fraction (`displaced_m / phys_dist`) into its sample list and mark
those samples as anchored at the CIFP elevation.  Lands with #2.

### 6. Plane fit / centroid math uses absolute meter coordinates
At tropical latitudes `|x| ≈ 8 × 10⁶`, which blows the condition
number of any normal-equation solve.  **Fix:** replace the projector
with an airport-centered version (`(lon-mid_lon, lat-mid_lat)` × scale)
so meter-space coordinates stay in the ±10 km range.  This was found
during the experiment and is the right fix in the legacy too.

## Phase E/F status (unfinished, kept in place)

Phase E (boundary band + tunnel portals) and Phase F (drainage zones)
are partially implemented and were producing useful output — the
boundary band especially gives the airport a clean perimeter at its
true elevation.  They have known rough edges:

- Phase E1 boundary band: per-piece centroid elevation sampling
  occasionally picks the wrong CIFP anchor when the band crosses a
  major elevation transition.
- Phase E2 tunnel portals: the grade-down rect generator works for
  simple road-under-airport cases but doesn't handle two-stage
  retaining-wall geometry.
- Phase F drainage: the infield triangulator runs but per-vertex
  elevations still occasionally end up below the surrounding pavement
  by more than a real drainage ditch would.

These will be revisited after the core (A-D) phases are clean.  For
now they remain in the source as-is.

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
