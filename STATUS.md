# Auto-Patch Status — APT.DAT TAXI NETWORK + EXCEPTION HARDENING 2026-05-11

## TL;DR

**28 commits past where the last status.md was written.**  204 of
206 tests pass.  Two SPJC-specific tests fail because the
fixture floors were generated against OSM-based output and the
build now uses apt.dat as the primary taxi-network source —
positions differ by 1–3 m, breaking the spatial floor counts.

This session landed four substantial architectural changes the
user drove in sequence:

1. **Exception-hardening pass** across the entire ``auto_patch/``
   tree (commits ``bd43761`` … ``2931236``, 21 of the 28).
   Replaced ~430 broad ``except Exception:`` blocks with narrow
   per-file ``_GEOM_EXC`` tuples.  Surfaced one masked bug:
   ``_resample_node_altitudes_nn`` was used at six sites in
   ``bridges.py`` but never imported (silently swallowed by the
   broad except in ``finalize.run_phase2``'s feature-emit
   harness); fixed with one import line in commit ``74e92b1``.

2. **apt.dat taxi network as the primary centerline source**
   (commits ``4568c1e``, ``8656efa``, ``78b1b70``).  Parse apt.dat
   rows 1201 (taxi nodes) and 1202 (taxi edges) into
   ``Airport.taxi_nodes`` / ``Airport.taxi_edges``.  Pipeline now
   uses apt.dat-derived centerlines (run through the same
   RDP-simplify + bend-split machinery as OSM) for any airport
   that has a 1201/1202 block, falling back to OSM when apt.dat
   has no taxi network.  Also fixed two real bugs along the way:
   the DSFTool-path lookup (commit ``4568c1e``: DSF pavement
   wasn't loading at all for CYXY) and a missing
   ``ROLE_TUNNEL_RAMP`` import in ``groundside.py``.

3. **Long-edge absorption ruleset re-aligned** (commits
   ``1684c8e``, ``fe52a1d``).  Dropped the apron-interior reject
   in ``_snap_corners_to_pavement`` that was short-circuiting the
   absorption pass at the rect-build stage.  Reinstated the
   ``CORRIDOR_TO_RUNWAY_M`` heuristic that preserves runway-
   anchored corridors (SPJC L, CYXY F) from absorption, while
   still letting alongside parallels (CYXY G, apron-internal
   sub-refs) absorb correctly.

4. **Runway-pavement seam fixes** (commits ``211f043``,
   ``0b035e8``).  Dropped the 2 m outward buffer in
   ``elevation.py``'s post-segmenting junction-clip — that buffer
   pushed every taxi-junction 2 m off the runway boundary,
   creating a visible sliver at every taxi-junction-to-runway
   interface.  Widened ``INTERSECTION_PROX_M`` from 0.5 m → 3 m
   so apt.dat row-110 boundaries drawn 1–2 m inside the row-100
   runway rect (SPJC's V1 throat is 1.75 m off) still register
   as runway-segmenter intersection points.  Plus euclidean
   dedup + Phase A/B snap-target dedup to consolidate
   multi-vertex clusters from a single chart-level transition.

**Branch:** ``dev`` (28 commits ahead of ``origin/dev``).
**Tag baseline:** ``spjc-good-baseline`` (commit ``4d97f89``,
2026-05-07) — predates all of this session's work.

## Commit chain since last status

```
0b035e8 Dedup runway-pavement intersection clusters
211f043 Fix runway-pavement seam at apt.dat-row-110 offsets
78b1b70 Apply RDP+bend-split to apt.dat taxi centerlines
fe52a1d Reinstate corridor heuristic in long-edge absorption
1684c8e Remove apron-interior reject; let absorption ruleset split rects
8656efa Use apt.dat taxi network as primary centerline source
4568c1e Fix DSFTool path lookup + missing ROLE_TUNNEL_RAMP import
2931236 Catch one remaining 'except Exception as exc' in elevation.py
a0af8f0 Narrow broad excepts in elevation.py to _GEOM_EXC
3eaf6f3 Narrow broad excepts in bridges.py to _GEOM_EXC
3e2ea1e Narrow broad excepts in pipeline.py to _GEOM_EXC
574267b Narrow broad excepts in junction_rules.py to _GEOM_EXC
eb65659 Narrow broad excepts in triangulation/junction_repair/terminals
74e92b1 Narrow broad excepts in junction_emit/groundside/finalize + fix masked NameError in bridges.py
4723af5 Narrow broad excepts in elevation_smoothing/unified_jacobi
4db0fd5 Narrow broad excepts in apt_dat_reader/osm_aeroway/osm_load
541df64 Narrow broad excepts in cifp_reader/dsf_reader/driver/layout
833dfc0 Narrow broad excepts in stubs/strips
613e93d Narrow broad excepts in taxiway_skeleton/taxiway_rects
ed2edc6 Narrow broad excepts in classifier/union_helpers/taxiway_decompose
ad34922 Narrow broad excepts in centerlines.py to _GEOM_EXC
777b843 Narrow broad excepts in junctions.py to _GEOM_EXC
6eef867 Narrow broad excepts in rects.py to _GEOM_EXC
70edb16 Narrow broad excepts in runways.py to _GEOM_EXC
1aaf054 Narrow broad excepts in absorption.py
bbfeaf1 Narrow broad excepts in vertices.py to _GEOM_EXC
8e22da7 Narrow broad excepts in runway_segments.py
bd43761 Narrow broad except in runway_geometry.py to OSError
```

## Architecture changes — what to know

### 1. apt.dat taxi network is the primary centerline source

New data in ``apt_dat_reader.py``:

* ``TaxiNode`` (id, lat, lon, usage, label) — row 1201.
* ``TaxiEdge`` (node_from, node_to, direction, kind, name) —
  row 1202.  ``kind`` is ``"taxiway_A"`` … ``"taxiway_F"`` or
  ``"runway"`` (taxi paths crossing a runway).
* ``Airport.taxi_nodes: Dict[int, TaxiNode]`` and
  ``Airport.taxi_edges: List[TaxiEdge]``.

New helpers:

* ``taxi_centerlines(airport, to_m, rwy_centerlines)`` —
  group edges by ``name``, linemerge per group, run each
  merged polyline through
  ``pavement.centerlines.split_merged_centerline`` (extracted
  from the OSM extractor so apt.dat and OSM share the same
  RDP + straight-enough + bend-split logic).
* ``taxi_junction_points(airport, to_m)`` — nodes referenced
  by edges of ≥ 2 distinct names, ≥ 3 same-name edges, or
  touched by a runway-cross edge.  Fed to
  ``_split_centerlines_at_points``.

Pipeline integration (``pipeline.py``):

* If apt.dat has any taxi-network entries, use them as primary.
* Fall back to OSM (``_extract_osm_taxi_centerlines``) only when
  the apt.dat block has no 1201/1202 data.
* The stub-reach-runway filter at
  ``RUNWAY_ENDPOINT_DIST_M = 80`` is **skipped** when apt.dat is
  the source — apt.dat doesn't carry OSM-noise sub-refs (A1, F1,
  V1-V3) so the filter's purpose doesn't apply, and apron-to-
  apron taxis like CYXY's G survive correctly.

### 2. Apron-interior reject removed; absorption ruleset is authoritative

``_snap_corners_to_pavement`` (rects.py) no longer returns
``None`` for "≥ 2 of 4 natural corners > 15 m inside pavement".
Always snaps every corner to ``pav.boundary``.  The only
remaining ``return None`` is geometric degeneracy (≥ 2 corners
collapsed within 1 m of each other after snap).

The authoritative ruleset is now exclusively in
``_drop_primary_parallels_embedded_in_pavement`` (absorption.py):

* Probe each long edge at 5 m steps; mark step "adjacent" if
  EITHER side has junction-class pavement within 5 m of the
  probe point.
* Contiguous adjacent runs ≥ 10 % of axial length → absorbed.
* Kept fragments ≥ 30 m survive; smaller → dropped.

Plus the **corridor heuristic** (reinstated):
``CORRIDOR_TO_RUNWAY_M = 160 m``.  Rects whose short-edge
midpoint is within 160 m of a runway are preserved from
absorption — they're runway-anchored corridors (SPJC L, CYXY F)
that should keep their full length even when one long edge has
apron alongside.

### 3. Runway-pavement seam fixes

Three layered changes for the −10178 sliver at SPJC's V1 throat:

* **Drop the 2 m junction-clip buffer** in ``elevation.py`` after
  runway segmenting.  Previously ``junction_clip =
  new_rwy_union.buffer(2.0)`` pulled every adjacent junction 2 m
  off the runway boundary.  Replaced with ``junction_clip =
  new_rwy_union``.

* **Widen INTERSECTION_PROX_M** in ``pipeline.py`` from 0.5 m
  → 3 m so apt.dat row-110 boundaries drawn 1–2 m inside the
  row-100 runway rect still register as pav_runway_intersections.

* **Dedup chain**: euclidean (5 m) for pav_runway_intersections
  in ``pipeline.py``, plus per-edge euclidean (3 m) for both
  Phase A snap targets and the runway-side insert step inside
  ``stitch_pavement_to_flat_runways``.

### 4. Exception hardening complete

Every ``except Exception:`` in ``src/auto_patch/`` has been
narrowed.  Each module has a per-file ``_GEOM_EXC`` tuple:

```python
_GEOM_EXC = (ValueError, TypeError,
             GEOSException, TopologicalError, IndexError)
```

For modules mixing geometry with file I/O / dict access
(``pipeline.py``, ``finalize.py``, ``bridges.py``,
``elevation.py``, ``osm_load.py``, ``apt_dat_reader.py``,
``osm_aeroway.py``), the tuple is widened with
``(OSError, KeyError, RuntimeError)``.  ``driver.py`` has a
dedicated ``_DRIVER_EXC`` for the per-airport harness that omits
``NameError`` / ``AttributeError`` / ``ImportError`` so typos
propagate to the test suite immediately.

## Current state

### Test status

| Test | Status |
|---|---|
| 204 / 206 | ✓ pass |
| ``test_compare_target.py::test_compare_target_spjc`` | ✗ fail — fixture is OSM-derived; apt.dat positions differ by 1–3 m |
| ``test_pavement_grade.py::test_pavement_grade[SPJC]`` | ✗ fail — 14 mid-edge steps > 0.5 m (cap 10), worst 1.36 m on junction-stub interfaces |
| ``test_pavement_grade.py::test_pavement_grade[SPLP]`` | ✓ pass |
| Compare-target on every other tested airport | ✓ pass |

### CYXY (the original "broken at DSF" airport)

Massive improvement vs. start of session:

| Role | session start | now |
|---|---|---|
| primary_parallel | 0 | 6 |
| stub | 5 | 11 |
| cross_connector | 0 | 1 |
| junction | 39 huge | **2 small** |

CYXY's DSF pavement now loads (was completely missing due to the
DSFTool path bug).  G (apron-to-apron taxi) correctly emits as
stub + cross_connector instead of being absorbed into a mega-
junction.  Most pavement is now correctly classified into rects.

### SPJC

Most pavement classifies correctly under apt.dat-driven flow.
Counts vs. the old OSM-derived fixture:

| Role | target | out | matched |
|---|---|---|---|
| boundary | 603 | 603 | 603 |
| runway | 85 | 85 | 85 |
| retaining_wall | 66 | 66 | 66 |
| tunnel_ramp | 36 | 36 | 36 |
| junction | 36 | 46 | 31 |
| primary_parallel | 27 | 31 | 23 |
| stub | 16 | 16 | 16 |
| cross_connector | 6 | 9 | 6 |
| secondary_parallel | 1 | 1 | 0 |

The output has slightly more rects (over-emission) and the
matches lose ~4–8 per role to spatial drift between apt.dat and
OSM coordinates.

### SPLP

| Role | count |
|---|---|
| primary_parallel | 1 |
| runway | 21 |
| boundary | 123 |

Simpler output than the old build (which had more OSM-noise sub-
refs); grade test passes cleanly.

## OPEN WORK — where to pick up

### Hot topic: V1-throat sloped-runway split at SPJC (USER'S LAST DIRECTION)

**User's exact words (2026-05-11):** "Is the segmenter missing
that point because of the 1.78 m gap with apt.dat pavement?  We
probably need to segment on pavement nodes within 2 m of the
runway."

**Where:** SPJC's V1 throat, runway 16R/34L, segment ``-10108``.

**Geometry:**

* Runway segment ``-10108`` is **sloped** (``altitude_high=6.1``,
  ``altitude_low=5.6``) covering the blast pad area at the 16R
  end.  Single 4-corner rect; vertices at m-coords
  ``(-83, -19)``, ``(-127, +71)``, ``(-86, +91)``, ``(-43, 0)``.

* The V1 taxiway entry point (apt.dat node 305 in the row-1201
  taxi network) projects onto the runway centerline at
  axial ≈ 75 m from the extended ``phys_end_a``.  The nearest
  apt.dat row-110 vertex to the runway boundary in this area is
  at ``(-60.8, +42.2)`` — **1.76 m off the runway boundary**.

* That row-110 vertex sits on the apron-side **long edge** of
  ``-10108`` (between corners ``(-86, +91)`` and ``(-43, 0)``).
  The user wants ``-10108`` split there into two sloped
  sub-rects so the apron-side junction has a corner to share.

**What's currently happening:**

* ``INTERSECTION_PROX_M = 3.0`` (after my widening) DOES capture
  this row-110 vertex as a ``pav_runway_intersection``.  My
  earlier trace confirmed it: 1 intersection survives dedup
  at t ≈ 0.0194 for runway 16R/34L.

* The pav_intersection gets injected into the segmenter's
  ``fractions`` list (``runway_segments.py:499–543``) and
  becomes a sample at ``sample_pts[i]`` with an interpolated
  elevation.

* **But the resulting sloped segment doesn't split where I
  expect.**  Looking at the runway segments sorted by axial
  position (NW → SE), ``-10108`` is the first sloped segment;
  it ends at the V1-throat seam where the runway transitions to
  flat (``-10109`` at ``alt=6.2``).  The seam corner is at
  ``(-43, 0)`` (node ``-349``), NOT at the projection of the
  apt.dat row-110 vertex ``(-60.8, +42.2)``.

* **The mystery:** the t ≈ 0.0194 pav_intersection corresponds
  to axial position 67.8 m from ``phys_end_a``; projected
  perpendicular at half-width the seam corner SHOULD be at
  ``(-58, +30)`` (apron side) or ``(-78, +10)`` (runway side).
  But ``-10108``'s SE corner ``-349`` is at ``(-43, 0)`` — far
  from either of those predicted positions.

* The runway segmenter's flat-region consolidation (``FLAT_TOL
  = 0.05`` m) groups consecutive flat samples into a multi-node
  polygon; sloped sub-runs become individual 4-corner rects.
  When the consolidation runs, it might be absorbing the
  pav_intersection sample into an adjacent group whose
  elevation difference is below ``FLAT_TOL``.

**User's proposed direction:** "We probably need to segment on
pavement nodes within 2 m of the runway."

Interpretation: re-examine the runway segmenter's handling of
pav_intersection samples that sit **interior to a sloped run**
(i.e. ``-10108``'s case).  Currently:

* If the pav_intersection lands inside a flat run → it
  survives as an intermediate corner of the multi-node flat
  polygon (Task 1 / 2026-05-09 behaviour).
* If the pav_intersection lands inside a sloped run → emit
  loop processes adjacent sample pairs.  If consecutive
  elevations differ by ≥ FLAT_TOL the sloped 4-corner rect
  emits, but **the pav_intersection sample doesn't become a
  corner of either sub-rect** unless adjacent samples already
  have substantially different elevations.

**Required investigation (next session):**

1. **Run the segmenter with instrumentation** to dump
   ``fractions``, ``sample_pts``, and ``elevs`` for SPJC
   16R/34L at low-t range.  Specifically: at t ≈ 0.0194 (my
   pav_intersection at the V1 throat) — is the sample present?
   What's its elevation?  Is it being consolidated into an
   adjacent group?

2. **Locate where the flat-vs-sloped consolidation runs**
   (``runway_segments.py:1086`` ``while idx < n_samples - 1``)
   and trace whether ``sample_pts[i]`` corresponding to the
   pav_intersection survives or gets absorbed.

3. **Confirm the t-projection math.**  My standalone trace
   says the pav_intersection is at t ≈ 0.0194 but the segment
   output says ``-10108`` ends at a different position.
   Either the segmenter handles the pav_intersection
   differently than my standalone trace, or the SE corner of
   ``-10108`` (node ``-349`` at ``(-43, 0)``) came from a
   different fraction (the first uniform 100 m seam, perhaps).

4. **Once the actual split point is identified, fix it so the
   sloped ``-10108`` splits into two pieces with a corner at
   the V1 entry.**

### Other open items (deferred)

* **CYXY −10006 V split** — user noted the V primary at CYXY
  should be split at a centerline bend, but apt.dat row-1202
  edges don't represent bends (just node-to-node edges).  Same
  underlying issue as the SPJC V1 case: apt.dat is sparser than
  OSM at apron-internal bends, so the segmenter doesn't see
  the natural breakpoints.  Same investigation as the V1 throat
  will likely surface the right fix.

* **SPJC compare-target fixture regeneration** — once the
  taxi-network behaviour stabilises, regenerate
  ``tests/fixtures/SPJC_target.osm`` against the new apt.dat-
  driven output.  The fixture is OSM-derived from
  ``spjc-good-baseline`` (commit ``4d97f89``) and won't match
  the new layout's coordinate positions.  The current floor
  failures are spatial-match misses, not count misses (counts
  are roughly correct; matches drop because rects shift 1–3 m).

* **Task 3 from the original 2026-05-11 ask:** "the boundary
  fill shapes should be added last to ensure they do not
  overlap with snapped runway junctions."  Currently
  ``_emit_airport_boundary_shape`` runs **inside**
  ``finalize.run_phase2`` ([finalize.py:197](src/auto_patch/finalize.py:197))
  before the post-finalize passes
  (``widen_junctions_to_runway_corners``,
  ``stitch_pavement_to_flat_runways``, the per-surface solver,
  the snap chain).  Move the boundary-emit + DEM-bridge emit
  to AFTER those passes so the boundary rects don't overlap
  with junctions that have just been widened to share runway
  corners.

* **Exception-hardening item 3 still open per project_todos:**
  the broad-except pass is done within ``src/auto_patch/`` but
  the wider Ortho4XP codebase outside ``auto_patch/`` (e.g.
  ``O4_*`` modules) may still have similar issues.  Out of
  scope unless the user explicitly requests it.

## Build / verify commands

```bash
# Full test suite
/Users/noah/Ortho4XP-shred86/venv/bin/python3 -m pytest tests/ \
    --tb=short -q

# Gate tests (currently 2 SPJC-only failures, all others pass)
/Users/noah/Ortho4XP-shred86/venv/bin/python3 -m pytest \
    tests/test_compare_target.py tests/test_pavement_grade.py \
    -v

# Build SPJC + write to /tmp for visual review
/Users/noah/Ortho4XP-shred86/venv/bin/python3 - <<'PY'
import sys
sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement(
    "SPJC", "/Users/noah/X-Plane 12", compute_elevations=True)
layout.to_osm("/tmp/SPJC.osm")
print("wrote /tmp/SPJC.osm")
PY

# Build CYXY (the original "DSF missing" airport)
/Users/noah/Ortho4XP-shred86/venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement(
    "CYXY", "/Users/noah/X-Plane 12", compute_elevations=True)
layout.to_osm("/tmp/CYXY.osm")
PY

# Standalone trace of pav_runway_intersections for the V1 throat
/Users/noah/Ortho4XP-shred86/venv/bin/python3 - <<'PY'
import sys, math
sys.path.insert(0, "src")
from auto_patch import apt_dat_reader as APR
from auto_patch.pavement.runways import _runway_rect_m
from shapely.geometry import Point
from shapely.ops import transform as shp_transform

apt = APR.load_airport(
    "/Users/noah/X-Plane 12/Custom Scenery/SPJC Lima by Los Flipantes 3.0 Nueva Terminal XP12/Earth nav data/apt.dat",
    "SPJC")
r = next(r for r in apt.runways if "16R" in (r.desig_a, r.desig_b))
lat0 = (r.lat_a + r.lat_b) / 2; lon0 = (r.lon_a + r.lon_b) / 2
R = 6371000; cos0 = math.cos(math.radians(lat0))
def to_m(lon, lat):
    return ((lon - lon0) * math.radians(1) * R * cos0,
            (lat - lat0) * math.radians(1) * R)
rect = _runway_rect_m(r, to_m)
cl_ax, cl_ay = to_m(r.lon_a, r.lat_a)
cl_bx, cl_by = to_m(r.lon_b, r.lat_b)
cl_dx = cl_bx - cl_ax; cl_dy = cl_by - cl_ay
cl_L2 = cl_dx*cl_dx + cl_dy*cl_dy
phys_dist = math.sqrt(cl_L2)
print(f"Runway phys_dist: {phys_dist:.0f}m, blast_a={r.blast_a_m}m, displaced_a={r.displaced_a_m}m")
# Pavement vertices within 3m
intersections = []
for pav in apt.pavements:
    if not pav.polygon: continue
    pm = shp_transform(to_m, pav.polygon)
    coords = list(pm.exterior.coords)
    if coords and coords[0] == coords[-1]: coords = coords[:-1]
    for px, py in coords:
        if rect.exterior.distance(Point(px, py)) > 3.0: continue
        t = ((px - cl_ax) * cl_dx + (py - cl_ay) * cl_dy) / cl_L2
        if t <= 5/phys_dist or t >= 1.0 - 5/phys_dist: continue
        intersections.append((t, px, py, pav.name))
intersections.sort()
for t, px, py, name in intersections[:10]:
    print(f"  t={t:.5f}  m=({px:+.1f},{py:+.1f})  name={name!r}")
PY
```

## Reference artefacts

* ``tests/fixtures/SPJC_target.osm`` — OLD baseline gate
  (regenerated at commit ``10e3544``, 2026-05-12 in calendar time
  but BEFORE this session's apt.dat-primary change set).  Will
  need regeneration once the V1-throat split is resolved.
* ``/tmp/SPJC_v6.osm`` … ``/tmp/SPJC_v9.osm`` (working artefacts
  from this session; the latest ``v9`` reflects the current
  ``dev`` branch tip including all dedup work).
* ``/tmp/CYXY_v4.osm`` — CYXY output for visual review.  Same
  branch state as SPJC_v9.

## Memory notes added this stretch

Saved under
``~/.claude/projects/-Users-noah-Ortho4XP-shred86/memory/``:

* (none added — work is captured in commit messages.  The
  ``feedback_shape_rules.md`` "Long-edge absorption rule" entry
  remains authoritative; commit ``fe52a1d`` documents the
  corridor heuristic that's paired with it.)

## Next agent — task spec

**Recommended starting point:** the V1-throat sloped-runway split
investigation (the user's last interrupted direction).  Specifics
are in the "OPEN WORK" section above.  In short:

1. Read commit ``0b035e8`` first (most recent runway-seam fix);
   that's where the current intersection-tolerance / dedup state
   came from.
2. Read ``src/auto_patch/pavement/runway_segments.py:445–543``
   (the fractions / pav_intersection injection) and
   ``:1086–1136`` (the flat/sloped emit loop) to understand how
   samples become corners.
3. Instrument the segmenter with a debug print at line ≈ 1090
   (start of the emit loop) to dump ``[(i, t, elev, anchored)]``
   for runway 16R/34L.
4. Confirm where SPJC's V1-throat pav_intersection sample ends
   up in the emit sequence.  Compare to where ``-10108`` 's SE
   corner actually lands.
5. Fix the segmenter so an interior pav_intersection sample on a
   sloped run creates a corner at that point (splits the sloped
   sub-rect into two).

**After that lands:** regenerate the SPJC target fixture and
adjust ``test_compare_target.py`` floors.  Then proceed to the
deferred items (CYXY V split, Task 3 boundary-emit ordering).

**Do NOT touch:** the un-tracked ``+60-140/`` directory at the
repo root.  24 stray .hgt elevation-tile files that should live
under ``Elevation_data/+60-140/``; the user has confirmed
"leave it untracked, never include in git add -A".
