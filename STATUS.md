# Auto-Patch Status — RUNWAY / BOUNDARY / TILE-CUT REFACTOR 2026-05-12

## TL;DR

**Twelve commits past the ``spjc-good-baseline`` tag.**  All 206
tests pass.  SPJC compare-target gate intact under regenerated
fixture.  SPLP grade test passes; runway-area step violations are
**zero** (all remaining steps are 0.5 m sub-metre artefacts at
terminal corners, totally unrelated to the runway).

This stretch landed eight architectural changes the user drove in
sequence:

1. Per-surface solver HARD-anchors flat runway segments (not just
   sloped 4-corner with ``altitude_high/low``).
2. New ``stitch_pavement_to_flat_runways`` pass — perpendicular-
   projection + near-edge snap to insert shared vertices on flat
   runway shapes from adjacent junction boundaries.
3. **Task 1** — runway flat-zone consolidation: consecutive flat
   segments collapse into a single multi-node flat polygon.
   *Sample-corner retention rule* (user 2026-05-12): keep EVERY
   uniform 100 m sample as an intermediate corner, not just
   pav_intersections, so junction-snap chains stay intact.
4. **Task 2** — runway-segment chain extends to cover apt.dat
   blast pads as displaced-threshold continuations.
5. Airport boundary emits as a **chain of 4-corner rectangles**
   (one per ~25 m densified segment), each tagged
   ``altitude_high``/``altitude_low`` for sloped pieces or
   ``altitude=`` for flat — JOSM-readable per-segment altitude
   profile.
6. Asymmetric boundary clamp: pull boundary UP toward runway
   when DEM is below the band; **never pull DOWN from DEM**
   (otherwise Ortho4XP's ``smooth_raster_over_airports`` drags
   the rendered terrain down into a canyon around the airport
   perimeter).
7. ``tile_cut.py`` — cut shapes along integer lat/lon tile lines
   with a 10 m buffer, then **drop pieces outside the current
   tile**.  Neighbour-tile auto_patch runs handle their own
   portion of cross-tile airports.
8. Fix to a year-old silently-disabled boundary runway-elevation
   clamp (missing ``_sample_runway_segment_elev`` import — the
   broad ``except Exception`` mask hid the NameError since the
   slice-5 refactor commit ``3f6bb89`` in 2026-04-26).  Then
   narrowed 28 ``except Exception:`` blocks in ``boundary.py`` to
   a focused ``(ValueError, TypeError, GEOSException,
   TopologicalError, IndexError)`` tuple so future bugs of the
   same shape surface immediately.

**Branch:** ``dev`` (12 commits ahead of ``origin/dev``).
**Tag baseline:** ``spjc-good-baseline`` (commit ``4d97f89``,
2026-05-07).

## Commit chain since baseline

```
10e3544 Drop tile-cut shape pieces outside the current tile
63edd2f Keep all intermediate samples as corners in multi-flat runway
f882250 Emit airport boundary as a chain of 4-corner rectangles
58ff723 Revert runway-chain-aware tile-cut; cut entirely at end
6c3323d Asymmetric boundary runway-clamp: pull UP only, never DOWN
8305131 Cut shapes along integer lat/lon tile boundaries
effacdf Narrow broad except blocks in boundary.py
5033ac2 Fix missing import that silently disabled boundary runway clamp
69929d6 Extend runway segment chain to cover blast pads as displaced thresholds
c04ac53 Consolidate flat runway segments into multi-node polygons
e5cc969 Stitch pavement to flat runway shapes (multi-node flat shapes)
1614676 HARD-anchor flat runway segments in per-surface solver
4d97f89 SPJC apron fix: skip airports subsumed by patches_area  ← spjc-good-baseline
```

## Architectural changes — what to know

### 1. Per-surface solver HARD-anchors flat runway segments (``1614676``)

The HARD-seed in ``elevation_per_surface/unified_jacobi.py`` used
to require BOTH ``altitude_high`` and ``altitude_low`` to be set
on a runway shape.  Flat segments (single ``altitude=``) silently
got zero HARD anchors — adjacent junctions saw the runway as
SOFT and didn't get lifted toward runway elevation.  Now both
forms HARD-anchor.

### 2. ``stitch_pavement_to_flat_runways`` (``e5cc969``)

New pass in ``junction_rules.py`` runs before the per-surface
solver.  Two phases:

* **Phase A — near-edge snap.**  For each pavement vertex within
  ``near_edge_snap_m=2.5 m`` of a flat-runway edge interior, snap
  onto the projection and insert a matching vertex into the
  runway shape.
* **Phase B — coincident-edge perpendicular projection.**  For
  each pavement edge whose endpoints both match flat-runway
  corners (a shared boundary edge), project every other pavement
  vertex perpendicular onto the shared edge.  Inserts at every
  projection in BOTH polygons.

### 3. Multi-flat consolidation with full sample retention (``c04ac53``+``63edd2f``)

The runway segmenter's flat-zone consolidator merges consecutive
flat samples into ONE multi-node polygon (single ``altitude=``
tag).  **All** intermediate samples are retained as corners —
both pav_intersection breakpoints AND uniform 100 m seams.

The "sample-corner retention" rule (``63edd2f``) was the user's
2026-05-12 correction to the original Task 1 implementation,
which had dropped uniform 100 m samples and broke the junction-
snap chain (``widen_junctions_to_runway_corners`` relies on
per-segment corner density to find nearby snap points).

### 4. Chain extends to blast pads (``69929d6``)

apt.dat row-100 ``blast_a`` / ``blast_b`` distances are absorbed
into ``displaced_a`` / ``displaced_b`` at the start of segment-
chain setup; ``blast_a/b`` then get zeroed.  CIFP threshold
elevations stay anchored at the displaced-threshold positions
(now interior to the extended chain).  Blast-pad areas get DEM-
sampled + grade-limited like the runway interior.

The legacy separate flat-rect blast-pad emit (at lines 1007-1038
of the pre-refactor ``runway_segments.py``) silently skips
because ``blast_a/b == 0`` after the extension.

### 5. Boundary as a chain of 4-corner rectangles (``f882250``)

``_emit_airport_boundary_shape`` no longer emits a single buffered
strip polygon with per-vertex ``node_altitudes``.  Instead it
walks the airport-boundary ring, densifies to 25 m, and emits one
4-corner rect per consecutive pair of densified points.  Each
rect is tagged ``altitude_high``/``altitude_low`` (sloped) or
``altitude=`` (flat) for direct readability in JOSM.

Half-width 2.5 m on each side of the boundary line.  Pavement
overlap is handled — rects entirely buried in pavement are
skipped; partial-overlap rects are trimmed against pavement.

### 6. Asymmetric runway clamp (``6c3323d``)

In ``_runway_clamped_alt`` (boundary ribbon) and ``_clamped_alt``
(DEM bridge):

```
band = distance * 0.03           # 3 % grade
lo   = runway_elev - band
if dem < lo:  alt = lo           # lift UP toward runway when DEM is low
else:         alt = dem          # follow DEM (no upper clamp)
```

The earlier symmetric clamp pulled boundary DOWN to ``runway +
band`` when DEM was higher.  Ortho4XP's
``smooth_raster_over_airports`` then dragged the rendered terrain
DOWN with it, producing a visible 20 m-deep canyon around the
airport perimeter at the SPLP north end.  Asymmetric clamp
preserves both rules: lift up when terrain dips below runway;
follow DEM when terrain is at or above runway.

### 7. Tile cut with drop-out-of-tile (``8305131`` + ``10e3544``)

New module ``tile_cut.py``.  Runs as a late post-process in
``pipeline.py``.  Three things:

* For each integer lat/lon line passing through the airport's
  pavement footprint, builds a buffered LineString (5 m each
  side → 10 m gap) and subtracts from every shape.
* Clips ``layout.airport_boundary`` itself to the current tile
  so downstream boundary-ribbon / DEM-bridge emit see only the
  in-tile portion.
* **Drops every shape piece whose representative point falls
  outside the current tile.**  The neighbour-tile auto_patch run
  generates the patch covering its portion.

The "drop" rule (``10e3544``) is the user's 2026-05-12
correction: previously cross-tile boundary rects ended up with
``altitude=0`` because the loaded DEM tile didn't cover their
location.  Now those rects are simply dropped — the neighbour
tile handles them with its own DEM context.

**Earlier chain-aware tile-cut in ``runway_segments.py`` was
reverted** (commit ``58ff723``) per user direction: pre-cut
chain breaking interfered with ``widen_junctions_to_runway_
corners`` (the widen pass promoted the gap-edge runway corners
into adjacent junctions, which then misaligned with the post-
process cut).  Cutting at the end uniformly is cleaner.

### 8. Boundary-clamp import fix + exception narrowing (``5033ac2`` + ``effacdf``)

``_sample_runway_segment_elev`` was used in ``boundary.py`` but
never imported — every call raised ``NameError`` and was silently
swallowed by the surrounding ``except Exception:`` block, leaving
the boundary runway-elevation clamp permanently disabled.  Bug
dated back to the slice-5 refactor (commit ``3f6bb89``,
2026-04-26).

Fixed with one import line.  Then narrowed all 28
``except Exception:`` blocks in ``boundary.py`` to
``_GEOM_EXC = (ValueError, TypeError, GEOSException,
TopologicalError, IndexError)`` so future ``NameError`` /
``ImportError`` / ``AttributeError``-on-typo propagate
immediately.

## Final grade state

### SPJC (compare-target baseline)

```
278 total shapes:
  boundary:           603 (chain of 25 m rects, post-cut, in-tile only)
  cross_connector:      6
  junction:            36
  primary_parallel:    27
  retaining_wall:      66
  runway:              85
  secondary_parallel:   1
  stub:                16
  terminal:             2
  tunnel_ramp:         36

Within-shape grade   : 7 violations, worst 2.2 % (junction)
```

Compare-target fixture at ``tests/fixtures/SPJC_target.osm`` was
regenerated at commit ``10e3544``.  Floor table in
``tests/test_compare_target.py:SPJC_BASELINE`` matches the new
target with ~99 % per-role floors.

### SPLP (cross-tile + custom apt.dat)

```
Within-shape:  4 violations (worst 13 % on small isolated
                junction -10043, the original SPLP residue —
                unchanged across this entire refactor)
Plane gradient: 1 violation at 1.52 % (junction -10025, just
                over the cap)
Cross-shape:   0
Vertex-edge step:    0
Mid-edge step: 2 (worst 0.55 m — terminal/primary_parallel
                interface, unrelated to runway)
```

The runway / blast-pad area is now **completely clean** in terms
of step violations.  The original SPLP 20-step "canyon at SW
threshold" report has dissolved entirely through the chain of
HARD anchors + multi-flat snap corners + asymmetric clamp +
tile-cut drop.

## Outstanding TODOs

### 1. Exception-handling hardening across auto_patch (medium effort)

``auto_patch/`` still contains ~412 ``except Exception:`` blocks
in adjacent files (``elevation.py``, ``junction_rules.py``,
``pavement/*.py``, ``bridges.py``, ``groundside.py``, etc.).
Same anti-pattern that silently masked the 1-year boundary import
bug AND a NameError on ``math.ceil`` / ``math.floor`` in the
tile-cut work (``8305131`` — the import was missing, the
NameError swallowed by the broad except wrapping
``generate_patch_osm`` in ``elevation.py:381``, collapsing the
entire SPJC runway chain to 2 shapes until I noticed).

**Approach:** one file per commit.  Define a per-file
``_GEOM_EXC`` tuple of expected geometry exceptions (model on
``boundary.py``):

```python
_GEOM_EXC = (ValueError, TypeError, GEOSException,
             TopologicalError, IndexError)
```

Replace ``except Exception:`` with ``except _GEOM_EXC:``.

**Acceptance:** no ``except Exception:`` in ``auto_patch/**/*.py``
outside of explicit driver-harness layers; ``NameError`` /
``ImportError`` / ``AttributeError``-on-typo propagate; all 206
tests still pass; SPJC compare-target gate intact.

See ``project_todos.md`` memory note for the full rationale.

### 2. Small isolated junction with 10 % within-shape grade (SPLP)

A 4-corner junction at SPLP (way id ``-10043`` in the current
build, but it renumbers run-to-run) has a 10.5 % within-shape
all-pair grade between two of its corners (~2.8 m delta over
~26 m).  No HARD-anchor neighbours, no pav_intersection
constraints.  The per-surface solver's cap-projection between
soft corners isn't moving them toward each other.

Probably a one-vertex solver convergence issue or a topology
quirk in pav_union → junction extraction.  Worth a separate
investigation.  Doesn't affect rendering (no nearby pavement
to step against — within-shape only).

### 3. Stale ``+60-140/`` directory at repo root

Untracked.  24 elevation tile files (~24 MB) that should live
under ``Elevation_data/+60-140/``.  Mentioned in the refactor-
state memory; leave it untracked, do NOT include in any
``git add -A``.

## Build / verify commands

```bash
# Full test suite
/Users/noah/Ortho4XP-shred86/venv/bin/python3 -m pytest tests/ \
    --tb=short -q

# Gate tests
/Users/noah/Ortho4XP-shred86/venv/bin/python3 -m pytest \
    tests/test_compare_target.py tests/test_pavement_grade.py \
    -v

# Build SPLP + grade check
/Users/noah/Ortho4XP-shred86/venv/bin/python3 - <<'PY'
import sys
sys.path.insert(0, "src"); sys.path.insert(0, "tools")
from pathlib import Path
import check_grade as CG
from auto_patch.pipeline import build_airport_pavement
xplane = "/Users/noah/X-Plane 12"
layout = build_airport_pavement("SPLP", xplane, compute_elevations=True)
out = Path("/tmp/SPLP_check.osm"); layout.to_osm(str(out))
CG.run_checks(out, max_grade_pct=1.5, proximity_m=1.0,
              edge_search_m=5.0, edge_step_m=0.5, top_n=5)
PY

# Build SPJC + dump shape counts
/Users/noah/Ortho4XP-shred86/venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from collections import Counter
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement("SPJC", "/Users/noah/X-Plane 12",
                                  compute_elevations=True)
print(Counter(s.role for s in layout.shapes))
PY
```

## Reference artefacts

* ``tests/fixtures/SPJC_target.osm`` — canonical SPJC output gate
  (regenerated 2026-05-12 at commit ``10e3544``).  Includes the
  tile-cut drop + boundary rect chain + multi-flat sample
  retention.
* ``Patches/-20-080/-13-077/SPLP_auto.patch.osm`` — SPLP runtime
  output (in tile -13/-77, the airport-anchor tile).
* ``Patches/-20-080/-13-078/SPLP_auto.patch.osm`` — older SPLP
  patch from a separate generation of tile -13/-78.  Will be
  refreshed when that tile is regenerated; current contents
  pre-date the tile-cut drop fix.
* CIFP data: ``/Users/noah/X-Plane 12/Custom Data/CIFP/SPLP.dat``.

## Memory notes added this stretch

Saved under ``~/.claude/projects/-Users-noah-Ortho4XP-shred86/memory/``:

* ``feedback_runway_end_pavement_classification.md`` — within 10 m
  of runway width → runway (hard anchor); wider → apron.
* ``feedback_flat_segment_node_count.md`` — flat shapes keep
  single ``altitude=`` tag, support arbitrary node count.
* ``feedback_boundary_clamp_asymmetric.md`` — clamp pulls UP only,
  never DOWN from DEM.
* ``project_todos.md`` — deferred work items (exception hardening
  is the standing item).

## Next agent — task spec

**Recommended starting point:** the exception-handling hardening
pass.  It's mechanical, low-risk, and the unmasking has already
caught two real bugs in this stretch — there are likely more
hiding behind the remaining 412 broad excepts.

**Order suggestion (smallest to largest):**
1. ``pavement/runway_geometry.py``, ``pavement/runways.py``,
   ``pavement/runway_segments.py`` (already partially audited
   when adding tile-cut intervals — same patterns as boundary.py).
2. ``pavement/vertices.py``, ``pavement/junctions.py``,
   ``pavement/rects.py``, ``pavement/centerlines.py``,
   ``pavement/absorption.py``.
3. ``junction_rules.py``, ``junction_repair.py``,
   ``junction_emit.py``.
4. ``elevation.py`` (largest — be careful, lots of paths).
5. Remaining: ``bridges.py``, ``boundary.py`` (already done),
   ``groundside.py``, ``pipeline.py``, ``finalize.py``,
   ``triangulation.py``, ``apt_dat_reader.py``, etc.

One file per commit.  Run the full test suite + SPJC compare-
target + SPLP grade test after each.  Any test failure means
the broad except was masking a real bug — investigate, fix,
then continue.

**Other potential follow-ups (in priority order):**
* Visual verification of SPLP in X-Plane after all the changes.
  The user has been driving this iteratively but hasn't yet
  confirmed the final state looks right in-sim.  The 20 m canyon
  at SPLP north and the SW threshold step issues should all be
  gone — worth a check.
* The ``-10043`` small-isolated-junction within-shape grade
  violation (see "Outstanding TODOs" above).
* Regenerate the ``Patches/-20-080/-13-078/`` SPLP patch file by
  running ``driver.generate_auto_patches`` for tile -13/-78.
  The current file is stale and predates the tile-cut drop fix
  (so it has 127 altitude-0 boundary rects); the regenerated
  file would have the correct in-tile portion.
