# Auto-Patch Tier 2 Implementation Plan

## Overview

Tier 2 builds on the Tier 1 CIFP-based runway slope patches with three features:
airport surface modeling, building flattening, and tunnel/overpass handling.

## Feature A: Airport Surface Elevation Modeling (Taxiway Grade Control)

### Problem

Tier 1 gives accurate runway slopes from CIFP data, but taxiways, aprons, and
ramps still rely on the polynomial-fit + Gaussian blur model. This can produce
unrealistic grades between runways and taxiway intersections.

### Current Data Flow

1. OSM provides taxiway centerlines as LineStrings (buffered to 15m polygons)
2. `encode_runways_taxiways_and_aprons()` fits degree-7 polynomials to DEM
   samples along each centerline
3. `weighted_alt()` blends fits with Gaussian weighting (sigma = runway_width)
4. Every taxiway/apron node gets an interpolated elevation from noisy DEM data

### Proposed Approach

Anchor the taxiway elevation model to known CIFP threshold elevations rather
than fitting to noisy DEM data:

1. **Anchor points** - From CIFP thresholds, derive exact elevations at each
   runway-taxiway intersection point via linear interpolation along the runway
   slope.

2. **Grade-constrained propagation** - Walk outward from anchor points along
   the taxiway centerline graph, enforcing max 1% grade per 25m. Where DEM
   suggests a steeper grade, clamp. Where DEM is consistent, use smoothed DEM.

3. **Apron/ramp fill** - For enclosed areas (aprons, ramps), use TIN
   interpolation from surrounding taxiway node elevations.

### Data Source: OSM vs apt.dat

Ortho4XP currently uses OSM exclusively for airport geometry. Both sources
have tradeoffs:

- **OSM**: Already loaded in pipeline, global coverage, often accurate. Can be
  incomplete or outdated for airports under construction.
- **apt.dat** (row codes 110-114): Very detailed geometry including Bezier
  curves, authoritative for X-Plane rendering. Requires new parser; custom
  scenery apt.dat files may override defaults.

**Proposal**: New tile-level config `apt_data_source` with values `"osm"`
(default), `"aptdat"`, or `"both"` (prefer apt.dat, fall back to OSM).

### Implementation Modules

- `O4_Aptdat_Parser.py` - New: parse apt.dat row codes 100, 110-114, 120-121,
  130 into taxiway LineStrings, pavement polygons, boundaries
- `O4_Auto_Patch.py` - Extend with taxiway grade propagation
- `O4_Cfg_Vars.py` - Add `apt_data_source` tile-level config
- `O4_Vector_Map.py` - Modify `include_airports()` to dispatch based on source

---

## Feature B: Building Flattening Within Airport Boundaries

### Problem

Buildings within airport perimeters need flat ground beneath them. Current
`encode_hangars()` flattens only if DEM range is <= 1.5m and uses raw DEM,
not the CIFP-derived surface model. Too conservative for sloped airports.

### Proposed Approach

After runway slopes (Tier 1) and taxiway grades (Feature A) establish the
correct airport surface model:

1. Query OSM for `building=*` ways within buffered airport boundary
2. For each building footprint, sample the auto-patch-derived surface elevation
   at the building centroid
3. Generate flat `altitude` patch polygons matching building footprints
4. Add to `_auto.patch.osm` alongside runway patches

### Implementation

Extend `generate_patch_osm()` in `O4_Auto_Patch.py` to accept building
footprints and their target elevations.

---

## Feature C: Tunnel and Overpass Elevation Handling

### Problem

Where roads pass under taxiways/runways, or where two highways cross at
different levels, the current system either excludes both (bridge/tunnel tags
cause exclusion from leveling) or forces them to the same elevation. This
produces terrain cliffs or roads floating through runways. Only fix today is
hand-crafted patches.

### Current State in Codebase

- `include_roads()` fetches `bridge` and `tunnel` as `tags_of_interest`
- Immediately excludes ways with those tags from leveling via
  `tags_for_exclusion`
- Tags are preserved in `osm_layer.dicosmtags["w"][wayid]` - data exists,
  just discarded
- OSM `layer` tag is NOT currently fetched but should be added

### Three Intersection Types

1. **Road-under-runway/taxiway tunnels**: Most critical. Road passes beneath
   airport surface. Detect by finding road ways with `tunnel=yes` that
   intersect runway/taxiway polygons. Tunnel road needs depression below
   airport surface by ~5.5m (5m clearance + road thickness).

2. **Highway-over-highway overpasses**: Two roads cross at different levels.
   Detect via intersection of road ways where one has `bridge=yes`/`layer=1`
   and the other has `layer=0`/`layer=-1`/`tunnel=yes`.

3. **Road-over-water/rail bridges**: Less critical but same mechanism.

### Patch Generation for Crossings

For each detected crossing, generate:

- **Approach ramp patches** - On either side, sloped `altitude_high`/
  `altitude_low` polygon transitioning the lower road from surface to tunnel
  depth over 50-100m (5-8% ramp grade).

- **Tunnel floor patch** - Under the crossing, flat `altitude` patch at
  depressed elevation (surface minus clearance).

- **Surface patch** - Upper surface at correct elevation from Tier 1/Feature A
  or DEM.

### Key Technical Consideration

X-Plane's terrain mesh is a single-layer heightfield. Tunnels need the mesh
to follow the UPPER surface (runway/overpass), while the tunnel road is
rendered as a 3D object below the mesh. Our patches ensure the mesh surface
is correct above the tunnel; the tunnel road exclusion from leveling (already
implemented) prevents the mesh from following the lower road.

The main value is in **portal approaches** - the ramps where roads descend
to tunnel entrances. The DEM rarely captures these artificial cuts.

### Implementation

- New function: `generate_crossing_patches()` in `O4_Auto_Patch.py`
- Add `layer` to `tags_of_interest` in `include_roads()`
- Categorize bridges/tunnels instead of blanket exclusion
- New config: `auto_patch_crossings` (bool, default True)

---

---

# Tier 3

## Feature D: Custom Scenery Mesh Import

### Problem

Many X-Plane airports ship with custom overlay DSF files containing
hand-crafted terrain meshes that precisely match their 3D terminal models,
jetways, and ramp areas. When Ortho4XP regenerates the terrain mesh for a
tile, this carefully tuned ground shaping is lost, causing buildings to
float or sink into the terrain.

### Proposed Approach

Allow the user to point to an existing scenery package (or DSF file). Parse
the mesh data from the overlay DSF and convert it into auto-patch polygons
using the `node_altitudes` per-vertex elevation format. This effectively
"imports" the scenery developer's ground work into the Ortho4XP pipeline.

1. **DSF mesh reader** - Parse the binary DSF GEOD atom (LOOP/SCAL vertex
   pools) and CMDS atom (triangle opcodes 23/24) to reconstruct a
   triangulated mesh with real-world lat/lon/elevation per vertex.

2. **Triangle-to-patch converter** - Each mesh triangle becomes a 4-node
   closed way (3 vertices + closing node) with `node_altitudes` carrying
   the three elevation values. The triangulator in `include_patches`
   reconstructs essentially the same mesh the scenery developer created.

3. **Priority and clipping** - If a custom mesh exists for an airport, its
   patches take precedence over auto-generated taxiway/building patches for
   the overlapping area. Clip auto-generated patches to avoid the custom
   mesh footprint, or skip auto-patch generation entirely for that airport.

### Existing Codebase Support

- `O4_DSF_Utils.py` already reads DSF atoms (DEMI/DEMD) in
  `extract_elevation_and_bathymetry_data()` and writes full DSF files
  including vertex pools and triangle commands.
- `O4_Overlay_Utils.py` round-trips DSF files through DSFTool for
  text-based editing.
- Legacy C reader in `src/Unused/C/dsf_io.c` provides a complete reference
  implementation for all DSF atom types.
- Vertex format is 5 x uint16 per vertex (X, Y, Z, normal_X, normal_Y)
  with pool scaling parameters mapping to real-world coordinates.

### Implementation

- New function: `parse_dsf_mesh(dsf_path, clip_boundary=None)` in
  `O4_Auto_Patch.py` - reads vertex pools and triangle commands, returns
  list of (lat, lon, elevation) triangle tuples.
- New function: `generate_mesh_import_patches(triangles, ...)` - converts
  triangles to `node_altitudes` patch ways.
- New config: `custom_mesh_path` (per-tile, path to scenery folder or DSF).
- Integration in `generate_auto_patches()` - if custom mesh present, emit
  imported patches and skip auto-generation for that airport.

---

## Proposed Config Parameters

| Parameter | Level | Type | Default | Purpose |
|-----------|-------|------|---------|---------|
| `apt_data_source` | tile | str | `"osm"` | Airport geometry source |
| `aptdat_path` | app | str | `""` | Path to apt.dat |
| `auto_patch_crossings` | tile | bool | `True` | Tunnel/overpass patches |
| `tunnel_clearance` | tile | float | `5.0` | Vertical clearance (meters) |
| `overpass_clearance` | tile | float | `5.0` | Bridge clearance (meters) |
| `custom_mesh_path` | tile | str | `""` | Path to scenery DSF for mesh import |

## Implementation Order

### Tier 2 (airport surface accuracy)

1. **Feature A (taxiway grades)** - DONE. Extends Tier 1 CIFP foundation
   with grade-constrained taxiway patches using per-node elevations.

2. **Feature B (building flattening)** - DONE. Flat patches under hangars,
   terminals, and general buildings within airport boundaries.

3. **Feature C (crossings)** - Most independent, uses existing OSM road data,
   addresses the most painful manual patching. Bridge/tunnel tags already
   fetched, just not used.

### Tier 3 (external data integration)

4. **Feature D (custom mesh import)** - Import hand-crafted terrain meshes
   from custom scenery overlay DSFs. Leverages `node_altitudes` patch format
   to faithfully reproduce scenery developers' ground shaping.
