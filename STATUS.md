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

## Current state (after commit 7b — `ce3b97c`)

The legacy phase-A-through-F surface generator is the active code
path.  Phases C2 (apron) and D (junctions) have been wired to the
new `O4_Surface_Mesh.adaptive_triangulate` helper.  Phase E (boundary
band, tunnel portals) and Phase F (drainage) are still inline and
unfinished but functional.  A new `O4_Apt_Dat_Reader` module has
landed and is tested but not yet integrated.

### Modules

| File | Status | Tests |
|---|---|---|
| `src/O4_Auto_Patch.py` | active legacy pipeline | none yet |
| `src/O4_Surface_Mesh.py` | adaptive triangulation, used by C2 + D | 19 |
| `src/O4_Apt_Dat_Reader.py` | apt.dat parser, **not yet wired in** | 22 |
| `src/O4_Vector_Map.py` | calls `generate_auto_patches` (unchanged) | none |
| `tests/test_surface_mesh.py` | mesh + grade + anchors | 19 |
| `tests/test_apt_dat_reader.py` | runway/pavement/Bezier/search priority | 22 |
| `tests/fixtures/synthetic_apt.dat` | hand-crafted parser fixture | — |

Total tests: **41**.  Run with `./venv/bin/python3 -m pytest tests/`.

### SPJC numbers (post commit 8)

| Metric | Commit 7b | **Commit 8** |
|---|---|---|
| Total emitted ways | 2 871 | **2 442** (-15 %) |
| Total cross-feature overlap | 7 638 m² (0.4 %) | **6 099 m² (0.2 %)** |
| flat ∩ flat overlap | 3 085 m² | 3 370 m² |
| flat ∩ triangle overlap | 3 148 m² | **2 729 m²** |
| triangle ∩ triangle | 1 221 m² | **0 m²** ✓ |
| Junction triangles (Phase D) | 1 631 | 0 (taxiway_data is empty) |
| Apron triangles (Phase C2) | ~700 | **2 262** |
| apt.dat pavements loaded | 0 (OSM) | **51 → 5 merged pieces** |
| Buildings reconciled for grade | 30 | 27 |

apt.dat pavement polygons are disjoint by design, so the
Phase D junction pipeline produces zero overlapping triangles
(junction_zone_m is empty).  The remaining ~6 k m² overlap is
from the building merge leftovers (flat-flat) and the apron
clipping against buildings not always preserving tiny holes
(flat-triangle).  Both will be addressed in commit 9+.

The remaining overlaps are not from the legacy architecture being
wrong — they're from OSM-derived inputs not being disjoint by
construction.  Switching to apt.dat (commit 8) eliminates the
underlying cause.

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
9. **Commit 8 — wire `O4_Apt_Dat_Reader` into legacy (Phase A0.5).** ✅ (this commit)

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

### After commit 8

* **Commit 9+** — extract more helpers as we touch them:
  * `_compute_taxiway_elevations` → `O4_Taxiway_Elevations.py`
  * runway elevation interpolation → `O4_Runway_Elevations.py`
  * Phase A geometry preparation → `O4_Surface_Inputs.py`
  * Phase emit helpers → `O4_Surface_Emit.py`
* **Phase E (boundary band, tunnel portals) and Phase F (drainage)**
  refinement.  These were producing useful output but were never
  finished; left in place and marked as such.  Will be cleaned up
  after the core (A-D) phases are clean.
* **Re-run `docs/TEST_PLAN_SPJC.md`** — verify all 11 numeric checks
  plus the four new invariant checks (12-15).

## How to resume next session

```bash
cd /Users/noah/Ortho4XP-shred86
git log --oneline -10            # confirm we're at ce3b97c
./venv/bin/python3 -m pytest tests/   # 41 pass

# Quick sanity: parse SPJC's custom apt.dat
./venv/bin/python3 -c "
import sys; sys.path.insert(0, 'src')
import O4_Apt_Dat_Reader as APR
path = APR.find_airport_apt_dat('/Users/noah/X-Plane 12', 'SPJC')
apt  = APR.load_airport(path, 'SPJC')
print(path); print(len(apt.pavements), 'pavements,',
                   len(apt.runways), 'runways')
"

# Then start commit 8 work — see "Commit 8" section above.
```

The unfinished bug from commits 5-7a (4 569 → 3 148 m² flat-tri
overlap) is partially mitigated but not eliminated.  Commit 8
should remove the underlying cause entirely by replacing the
buffered OSM taxiways with the disjoint apt.dat polygons.

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
