# SPJC Surface-Patch Test Plan

This document describes how to test `src/O4_Surface_Patch.py`
(the analysis+emit pipeline for runways, buildings, taxiways, aprons)
against the reference airport SPJC (Jorge Chávez, Lima, Peru) in tile
`-13/-078`.  It is the persisted reference for future sessions — when
iterating on `build_surface_patch`, run this plan end-to-end and
verify every check before declaring the change good.

## Hard surface invariants (apply to ALL emitted geometry)

These are non-negotiable constraints that every output must satisfy.
They are stronger than the numbered checks below — a build that
violates any of these is broken regardless of the check counts.

1. **No overlapping shapes.**  Every emitted way is disjoint from
   every other emitted way (interiors do not intersect).  Shared
   boundaries are allowed (and required, see #2); shared interiors
   are not.

2. **Continuous surface.**  Runways, taxiways, aprons, terminals,
   and buildings together form a single continuous painted surface.
   Terminals join to aprons, aprons join to taxiways, taxiways join
   to runways.  Where two features meet they must share an actual
   geometric boundary, not merely abut with a gap.

3. **Equal altitudes at joins.**  Where two shapes share a boundary
   (or even a single vertex) their elevations at that boundary are
   identical.  No cliffs, no steps.

4. **Triangle wedges for perpendicular slope changes.**  When the
   slope direction of one shape is perpendicular (or otherwise
   incompatible) with that of an adjacent shape, the join is bridged
   by triangle wedges that smoothly transition the elevation.  Each
   triangle is a 3-vertex sloped patch.  Triangles are bound by the
   same 1.0 % apron / 1.5 % runway-taxiway grade limits.

5. **All slope/grade rules continue to apply.**
   - Runways:        ≤ 1.5 %  longitudinal,  vertical-curve rule applies
   - Taxiways:       ≤ 1.5 %  longitudinal
   - Aprons + transition triangles:  ≤ 1.0 %  any direction
   - Building/terminal pads:  flat, but their elevation is a free
     variable (within ±BLDG_ADJUST_MAX from DEM) for reconciliation.

These invariants change the surface model from the current
"independent shapes per feature" approach to a **planar subdivision**:
the airport surface is a single complex polygon that is partitioned
into runway / taxiway / apron / terminal / triangle cells, every cell
edge being shared between exactly two cells (or one cell + the
unrendered exterior), and all cells agreeing on elevation along
their shared edges.

## Rationale

The legacy pipeline in `O4_Auto_Patch.generate_airport_surface_patches`
is still in place behind the `USE_NEW_SURFACE_PIPELINE` flag.  The new
module is the active area of work.  Because features (runways,
taxiways, aprons, terminals) all have to join smoothly while meeting
FAA/EASA slope limits, the new pipeline uses a two-phase analysis →
emit structure (see STATUS.md "Airport elevation-build order").  This
test plan exists to verify that structure holds up at SPJC, where the
airport slopes ~39 m to ~8 m and an earlier heuristic apron algorithm
produced catastrophic output (509 562 m² blobs, 296 934 m² apron-apron
overlap).

## Inputs and paths

All paths are absolute on the current machine.

| Input                      | Path                                                              |
|----------------------------|-------------------------------------------------------------------|
| Working directory          | `/Users/noah/Ortho4XP-shred86`                                    |
| Ortho4XP config            | `Ortho4XP.cfg` (reads `custom_scenery_dir`)                       |
| X-Plane root               | `/Users/noah/X-Plane 12/`                                         |
| CIFP directory (preferred) | `/Users/noah/X-Plane 12/Custom Data/CIFP/`                        |
| SPJC CIFP                  | `.../Custom Data/CIFP/SPJC.dat`                                   |
| Global airports `apt.dat`  | `/Users/noah/X-Plane 12/Global Scenery/Global Airports/Earth nav data/apt.dat` |
| OSM airport cache          | `OSM_data/-20-080/-13-078/-13-078_airports.osm.bz2`               |
| OSM building cache         | `OSM_data/-20-080/-13-078/-13-078_apt_bldg_local.osm.bz2`         |
| DEM tile                   | `Elevation_data/-20-080/S13W078.hgt`                              |
| Existing reference output  | `Patches/-20-080/-13-078/SPJC_auto.patch.osm`                     |
| venv python                | `venv/bin/python3`                                                |

How `src/O4_Vector_Map.py` locates the CIFP dir when the config doesn't
set `cifp_data_path` explicitly:
```
xplane_root = dirname(normpath(CFG.custom_scenery_dir))
candidate   = os.path.join(xplane_root, "Custom Data", "CIFP")
```
Both `Custom Data/CIFP` and `Resources/default data/CIFP` contain
`SPJC.dat`.  Prefer Custom Data (newer AIRAC cycle).

Tile coordinates: SPJC lat=-12.0219, lon=-77.1143.  Floor → `tile.lat=-13`,
`tile.lon=-78`.  Patch directory is
`Patches/tile_lat_floor_10/tile_lat_lon_coarse` → `Patches/-20-080/-13-078/`.

## Test artifacts (not committed)

Scratch files that should **not** be checked in:

```
/tmp/SPJC_new.patch.osm                  — output of the standalone test
/tmp/spjc_analysis_report.txt            — analyzer output
scripts/test_spjc_surface_patch.py       — ad-hoc test driver (temp)
scripts/analyze_spjc_patch.py            — ad-hoc analyzer (temp)
```

## Test 1 — Standalone driver

Short-circuit the tile pipeline and call `build_surface_patch` directly
on SPJC.  This is the inner loop: seconds per iteration, no DEM tile
loading for other airports, no full OSM pass.

Sketch:

```python
# scripts/test_spjc_surface_patch.py  (ad-hoc; do not commit)
import os, sys
os.chdir("/Users/noah/Ortho4XP-shred86")
sys.path.insert(0, "src")

import O4_File_Names as FNAMES
import O4_Config_Utils as CFG
import O4_OSM_Utils as OSM
import O4_Airport_Utils as APT
import O4_Auto_Patch as AP
import O4_Surface_Patch as SP
import O4_DEM_Utils as DEM
import O4_Tile_Utils as TILE

# 1. Tile object + DEM.  Use whatever constructor matches the main
#    pipeline (see O4_Vector_Map.build_poly_file for the reference
#    invocation).
tile = TILE.Tile(-13, -78, "default")
tile.dem = DEM.DEM(tile.lat, tile.lon,
                   tile.custom_dem, tile.fill_nodata or "to zero",
                   info_only=False)

# 2. Airport layer + dico_airports (same sequence as build_poly_file)
airport_layer = OSM.OSM_layer()
airport_layer.update_dicosm(
    FNAMES.osm_cached(tile.lat, tile.lon, "airports"),
    {"n": [], "w": [("aeroway", "")], "r": []},
    {"n": [], "w": [("aeroway", "")], "r": []},
)
dico_airports = APT.build_airport_array(tile, airport_layer)
APT.smooth_raster_over_airports(tile, dico_airports)

# 3. Taxiway + building extraction (these reuse legacy helpers)
taxiway_data = AP.extract_taxiway_info(airport_layer, dico_airports, tile)
building_layer = OSM.OSM_layer()
building_layer.update_dicosm(
    FNAMES.osm_cached(tile.lat, tile.lon, "apt_bldg_local"),
    {"n": [], "w": [("building", "")], "r": []},
    {"n": [], "w": [("building", "")], "r": []},
)
building_data = AP.extract_building_info(
    airport_layer, dico_airports, tile, building_layer=building_layer)

# 4. CIFP for SPJC
cifp_path = "/Users/noah/X-Plane 12/Custom Data/CIFP"
runways = AP.parse_cifp_file(os.path.join(cifp_path, "SPJC.dat"))
runway_pairs = AP.pair_runways(runways)

# 5. Run the new pipeline
lines, nid = SP.build_surface_patch(
    icao="SPJC",
    taxiway_data=taxiway_data.get("SPJC", []),
    building_data=building_data.get("SPJC", []),
    runway_pairs=runway_pairs,
    tile=tile,
    dico_apt_entry=dico_airports["SPJC"],
    start_node_id=-10000,
    road_data=None,
)

# 6. Write scratch output (does NOT touch Patches/)
with open("/tmp/SPJC_new.patch.osm", "w") as f:
    f.write("<?xml version='1.0' encoding='UTF-8'?>\n")
    f.write("<osm version='0.6' generator='O4_Surface_Patch'>\n")
    f.write("\n".join(lines))
    f.write("\n</osm>\n")
```

Run:

```bash
./venv/bin/python3 scripts/test_spjc_surface_patch.py
```

**Watch the log lines emitted by `build_surface_patch`:**

```
[surface_patch] SPJC: N CIFP anchors
[surface_patch] SPJC: M aprons, worst anchor residual X.XX m, worst strip grade 0.XXX
[surface_patch] SPJC: emitted N runway segs + N buildings + N taxiway segs + N apron strips (N ways total)
```

The `worst strip grade` must be ≤ 0.010.  The `worst anchor residual`
should shrink over reconciliation iterations — if it stays large (>2 m),
reconciliation is blocked by `BLDG_ADJUST_MAX`.

## Test 2 — Analyzer checks

Parse `/tmp/SPJC_new.patch.osm` and run the following numeric checks.
All must pass.

```python
# scripts/analyze_spjc_patch.py  (ad-hoc; do not commit)
import xml.etree.ElementTree as ET
from shapely.geometry import Polygon
from shapely.ops import unary_union

tree = ET.parse("/tmp/SPJC_new.patch.osm")
nodes = {n.get("id"): (float(n.get("lon")), float(n.get("lat")))
         for n in tree.iter("node")}
ways = []
for w in tree.iter("way"):
    nids = [nd.get("ref") for nd in w.findall("nd")]
    tags = {t.get("k"): t.get("v") for t in w.findall("tag")}
    ring = [nodes[n] for n in nids]
    ways.append((ring, tags))
```

### Check 1 — No `node_altitudes` tags
Hard invariant (STATUS.md).  `grep node_altitudes /tmp/SPJC_new.patch.osm`
must return zero hits.

### Check 2 — Apron area within 10 % of OSM source
Sum the area of every emitted way whose tag set looks apron-ish (flat
polygons below the max runway elevation and above sea level, plus all
sloped rects that are not runways or taxiways).  Compare to the sum of
polygon areas in `dico_airports["SPJC"]["apron"][0]`.  Emitted ≥ 0.90 ×
source (slightly less is fine because we subtract buildings and
runways).

### Check 3 — No monster apron polygons
`max(way.area) ≤ 50_000 m²`.  Catches the 509 562 m² pathology of the
old computed-apron heuristic.

### Check 4 — Apron-apron overlap ≈ 0
For all apron ways: `unary_union(all).area ≥ sum(each.area) − 100 m²`.
The new pipeline dissolves overlapping OSM apron polygons at intake, so
this should be essentially zero.

### Check 5 — No apron overlap with building or runway
For each apron way: `poly.intersection(bldg_union).area < 1 m²` and
`poly.intersection(runway_union).area < 1 m²`.  Verifies intake-time
subtraction stuck.  **If this fails for sloped strips only, switch
sloped-strip emission to use the clipped polygon outline instead of
its rotated bbox** — the rotated bbox is a known approximation that
can extend slightly outside the clipped shape.

### Check 6 — Every apron strip ≤ 1.0 % grade
For each sloped apron way, compute actual slope from its 4 corners:

```python
grade = abs(eh - el) / perpendicular_distance_between_long_edges_m
assert grade <= 0.010 + 1e-3
```

### Check 7 — Building ↔ apron continuity (**load-bearing**)
For every building way whose polygon touches an apron:
- Find the apron strip whose polygon contains the midpoint of the
  shared boundary.
- Evaluate the strip's plane at that midpoint.
- Assert `|plane_value − building_elevation| ≤ 0.1 m`.

This is the core smooth-join check.  If it fails for many buildings,
reconciliation isn't converging — first suspect is that the 0.5 ×
nudge step + 3 iterations is too timid, or `BLDG_ADJUST_MAX = 2.0 m`
is binding on too many buildings.

### Check 8 — Taxiway ↔ apron continuity
For every taxiway segment whose centerline crosses an apron polygon:
- Interpolate `(altitude_high + altitude_low) / 2` at the crossing.
- Evaluate the apron strip's plane at the same crossing point.
- Assert match within 0.1 m.

### Check 9 — Taxiway grade compliance
Every taxiway sloped-rect: `|altitude_high − altitude_low| / length_m ≤
0.015 + 1e-3`.

### Check 10 — Building elevations within ±2.1 m of DEM
Re-sample the DEM at each building centroid (replicating
`_building_elevations` step 1).  For every building:
`|final_elev − initial_dem_elev| ≤ BLDG_ADJUST_MAX + 0.1 m = 2.1 m`.
Catches reconciliation runaway.

### Check 12 — Zero overlapping interiors
For every pair of emitted ways `(A, B)`:
`A.poly.intersection(B.poly).area ≤ 1 m²` (1 m² tolerance for
floating-point boundary slop).  This is the "no overlaps" invariant.

### Check 13 — Continuous surface (no gaps at joins)
Build the union of all painted shapes (`unary_union(all_ways)`).
For each runway, taxiway segment, apron strip, terminal, and
building pad, walk its boundary and verify every point is within
0.5 m of another shape's boundary OR is on the airport's true
exterior boundary.  Any boundary segment > 1 m long that is NOT
shared with another shape and NOT on the exterior is a gap.

### Check 14 — Elevation continuity at shared edges
For every pair of emitted ways that share a boundary edge longer
than 1 m, evaluate each way's elevation at the midpoint of the
shared edge.  The two elevations must agree within 0.05 m.

### Check 15 — Triangle wedges where required
Emit a list of every shared edge where two adjacent shapes have
slope axes that are NOT parallel (within 30°).  Each such edge
must either:
  (a) be bridged by one or more triangle wedge(s) so that no two
      sloped rectangles meet directly across an angled boundary, OR
  (b) the boundary is degenerate (length < 5 m).
A failure list pinpoints "fix here" locations.

### Check 11 — Coverage
Non-zero counts in each category: runway segs, buildings, taxiway
segs, apron strips.  SPJC has OSM aprons; if apron count is zero,
intake is broken.

## Test 3 — Visual inspection in JOSM

Load `/tmp/SPJC_new.patch.osm` in JOSM and verify:

- `altitude_high` / `altitude_low` tags are editable on sloped rects.
- Every way tagged sloped has exactly 4 vertices + closing node.
- No `node_altitudes` tag anywhere.
- Apron strip rectangles visually tile their parent OSM apron polygon
  without protruding into buildings or runways.
- Strip orientation is rotated to match the slope gradient (not
  axis-aligned to N/S/E/W) — this is correct, not a bug.

## Test 4 — End-to-end pipeline run

After tests 1–3 pass, run the real pipeline and compare to the
standalone output:

```bash
./venv/bin/python3 Ortho4XP.py   # select tile -13 -78, Step 1
diff /tmp/SPJC_new.patch.osm Patches/-20-080/-13-078/SPJC_auto.patch.osm
```

The files should be byte-identical or nearly so.  If not, identify
whether the wrapper pipeline passes different `taxiway_data` /
`building_data` / `dico_apt_entry` values than the standalone driver
used.

## Expected failure modes (and what to fix)

1. **Check 7 fails broadly** — reconciliation is too timid.  Raise
   step to 0.8 × delta, raise `RECONCILE_MAX_ITERS` to 5.  If buildings
   are hitting `BLDG_ADJUST_MAX` cap, relax to 3 m.
2. **Check 6 fails occasionally near strip edges** — post-clamp
   assertion should re-clamp; verify the re-center step in
   `_split_apron_into_strips`.
3. **Check 5 fails for sloped strips only** — rotated-bbox emission
   extends outside clipped polygon.  Switch to polygon outline for
   sloped strips (accept that X-Plane may not render slope tags on
   non-rect polygons cleanly; needs experimentation).
4. **Check 2 fails — emitted area too low** — intake subtraction is
   over-eating.  Most likely the `apron.difference(bldg_union)` is
   clipping valid apron because a building polygon is slightly
   oversized from `BLDG_PAD`; consider shrinking pad or subtracting
   un-padded footprints.
5. **Check 4 fails — apron-apron overlap** — `unary_union` at intake
   not running or failing silently; check the `except Exception: pass`
   branches.
6. **`worst anchor residual` in the log stays high** — likely cause:
   the reconciliation loop isn't running, or the strip plane fits
   aren't incorporating building edges.  Instrument
   `_collect_apron_anchors` to log how many anchors each apron got.

## Features intentionally NOT tested here

- Standalone taxiways outside any apron (they're emitted by the new
  pipeline but only exercised incidentally at SPJC).
- Tunnel portals, drainage, edge-aware boundary bands (deferred per
  STATUS.md "Next step").
- Other airports — SPJC until clean, then widen.
- Performance on a full-tile run (profile only if needed).
