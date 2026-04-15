# Auto-Patch Refactor — Status

**Current direction (head = `693c4f0`):** the legacy
`generate_airport_surface_patches()` in `src/O4_Auto_Patch.py`
remains the active code path, now with a dedicated Phase C0 for
apt.dat-sourced taxiways.  Taxiway polygons are classified and
emitted as chains of sloping rectangles via a four-strategy
pipeline (MRR fast path → Voronoi skeleton → morphological
decomposition → apron triangulation fallback), with runway-edge
anchoring tied to the runway chain actually emitted by
`generate_patch_osm`.  The "new pipeline" experiment in
`src/O4_Surface_Patch.py` is retired — read the "Why we reverted"
section below.  **Most recent session (head = `693c4f0`)
delivered Tasks A + B on the user's four-item taxiway refinement
list; Task 2 (apron-prefer-flat) and Task 4 (building-pad
reconciliation) are the next priorities.**

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

## Current state (after commit `693c4f0`)

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

**Taxiway emission is the major change since commit 16.**  Every
apt.dat pavement polygon is classified apron / taxiway by
O4_Pavement_Classifier, and taxiway polygons flow through a new
Phase C0 that emits them as chains of sloping rectangles via
one of four strategies tried in order:

1. **Fast path:** MRR-aligned rect chain if the polygon is a
   clean strip (aspect ≥ 3.5, width 9–45 m, aspect-scaled fit
   ratio).
2. **Voronoi skeleton:** centerline extraction via
   `shapely.ops.voronoi_diagram`, linemerge, drop < 50 m spurs,
   RDP 5 m simplification.  Each centerline is walked in 50 m
   steps and emitted as a local-width rect chain.  Handles
   curved taxiways and multi-taxi mega-polygons.
3. **Morphological decomposition:** buffer(-15).buffer(+15)
   isolates wide junction hubs from strip branches; each branch
   is re-probed with the fast path.  Handles polygons with
   dramatically wider junction hubs.
4. **Apron triangulation:** last-resort fall-through.

A running global rect-union (`global_twy_rect_union`) deduplicates
emissions across all source polygons so sibling branches at a
Y-junction can't overlap.  Chain elevations are RDP-simplified
with a 1 m tolerance so long straight stretches collapse to one
rect (user directive: "unless there's more than 1 m elevation
change over 30 m, simplify and combine").

Runway-taxi joins are anchored to the runway's own emitted
elevation: `generate_patch_osm` now returns its runway segment
chain alongside the OSM string, and Phase A3 builds a meter-space
lookup (`_runway_elev_lookup`) that interpolates the nearest
runway segment.  Taxi-centerline endpoints within 25 m of a
runway polygon (a distance sized for the fact that a Voronoi-
skeleton leaf sits ≈ half-width inside the polygon boundary) are
pinned to the runway's own elevation and propagated through the
grade clamp via `_clamp_profile_with_anchors`.

apt.dat pavement polygons are simplified to 1 m tolerance at
ingest (user directive: "we can simplify taxiway shapes, anything
less than a meter — too many shapes hurts performance").

At SPJC the full pipeline emits **1 006 m² of intra-sloped
overlap** (0.03 % of total area, where adjacent rects in one
chain meet at curve joins) and **zero cross-category overlap**.
Two of the three close taxi-runway joins at SPJC match the
runway elevation exactly; one has a 0.7 m residual to be fixed.

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

### SPJC numbers (post commit `693c4f0`)

Combined surface + runway-segment patch (`SPJC_auto.patch.osm`):

| Metric | Commit 16 | **Commit `693c4f0`** |
|---|---|---|
| Total emitted ways | 3 009 → 5 184 | **2 117** |
| Flat polygons | 545 → 546 | **824** |
| Sloped rectangles (runway + boundary + **taxiway**) | 294 → 293 | **364** |
| Triangles (apron + transition strips) | 2 170 → 4 345 | **929** |
| apt.dat taxiway classifications | — | **9 classified** (7 decomposed) |
| apt.dat taxiway rectifiable polygons | — | **52** |
| apt.dat taxiway rect segments emitted | — | **76** (4 flat, 72 sloped) |
| Cross-category overlap | 0 m² ✓ | **0 m²** ✓ |
| Intra-sloped overlap | 0 m² ✓ | **1 006 m²** (0.03 %) |
| Intra-triangle overlap | 0 m² claimed | **0 m²** ✓ |
| Max runway longitudinal grade | 1.50 % | **1.50 %** ✓ |
| Max taxiway longitudinal grade | 1.50 % | **1.50 %** ✓ |
| Max apron longitudinal grade | (apt.dat path was 1.5 %) | **1.00 %** ✓ |
| Taxi-runway close joins matching exactly | n/a | **2 of 3** |
| Residual join mismatch (worst case) | n/a | **0.70 m** (1 rect) |

**Note on the commit-16 numbers:** the "0 m² overlap" reported at
commit 16 was a cross-feature check that did not catch the
intra-triangle overlap caused by apt.dat's overlapping
"Base Ramp" + named sub-ramps.  At the start of the 693c4f0
session the audit script showed **928 004 m² intra-triangle
overlap** from that layering — ~35 % of the triangle sum.  Task 3
(commit `239efa9`) deduplicated by unioning the apt.dat pool
before emission, which dropped the intra-triangle overlap to
exactly **0 m²**.

Run `./venv/bin/python3 /tmp/run_legacy.py` then
`./venv/bin/python3 /tmp/audit_legacy.py` to reproduce.  All 99
auto-patch tests pass (134 total counting the boundary agent's
test module, if present).

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
28. **Task B refined — (a) endpoint-based crossing detection replaces the per-sample distance check, and (b) generate_patch_osm returns its emitted runway-segment chain so the taxi-runway elevation lookup uses the SAME per-segment elevation the runway actually emitted (not raw CIFP linear interpolation).** ✅ `693c4f0` (this commit)

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
git log --oneline -15            # confirm we're at 693c4f0
./venv/bin/python3 -m pytest tests/   # 99 local, 134 with boundary-model

# Full pipeline sanity run + audit (driver at /tmp/run_legacy.py):
./venv/bin/python3 /tmp/run_legacy.py      # writes /tmp/SPJC_legacy.patch.osm
./venv/bin/python3 /tmp/audit_legacy.py    # see the SPJC numbers table above
```

Next task: **aprons — prefer flat, slope minimally.**  See
item 2 of the Next Task Queue.  Key files and entry points:

- Phase C2 apron emission: `src/O4_Auto_Patch.py` around line
  4010, search for `for ai, p_m in enumerate(apron_polys_m)`.
- Taxi rects that the apron path should treat as anchors:
  `emitted_twy_quads_m` at Phase C0 (line ~2900).
- Runway segment chain available via `runway_segment_chain`
  parameter / `_runway_segs_m` meter-space index.
- Apron anchor collection: `_apron_anchors` helper at line
  ~3790 already merges building edges + taxi crossings +
  terminal flat rings.  Extend it to also read from
  `emitted_twy_quads_m` so apron triangles see the C0 rects
  as immovable boundary conditions.
- Building pad elevations live in `bldg_elevations` (Phase
  A5) and are currently fixed at DEM-sampled centroid; task 4
  will move them within `±BLDG_ADJUST_MAX` to satisfy the
  apron grade budget.

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
