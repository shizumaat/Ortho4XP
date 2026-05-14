# Auto-Patch Status — 2026-05-14 (end of session): handoff for grass-gap regression

## TL;DR

Session worked through a long chain of fixes culminating in a
known regression that needs careful investigation by the next
agent.  All landed work is committed; one investigation is in
progress.

**Committed this session (in order):**
1. `22bacf2` — Stub centering + apt.dat-driven stub emit +
   canonical runway emit (5 sub-fixes; see commit message).
   Triangle-prevention fix from prior session was rolled back at
   user direction.  SPLP -77 fixture regenerated.
2. `b698166` — Per-runway ref on segmented runway shapes.
   Previously every runway segment at CYXY got the merged
   ``02/20/14L/32R/14R/32L`` ref.  Now segments carry their own
   ``14R/32L`` / ``14L/32R`` / ``02/20`` ref via a per-segment
   ``desig_pair`` extension to ``runway_chain`` tuples.

**Tests: 221 / 221 pass at HEAD.**

## In-progress investigation — DO NOT SKIP THIS BEFORE PROCEEDING

User report 2026-05-14: at CYXY, way -10002 (primary parallel E,
top end where it meets the south apron) should reach DEM ~715 m
but the solver caps it at 710.1 m.  This is a regression — the
same vertex reached 717 m at commit ``0231eaa`` (this session,
the seam-anchor commit, before the GeometryCollection fix).

### What was confirmed

* **Bisected.**  Regression appeared at commit ``bc39174`` ("Fix
  CYXY missing junctions: handle GeometryCollection residue in
  junction_emit").  Pre-bc39174, CYXY emitted only 2 junction
  polygons (49 silently dropped).  Post-bc39174, all 51 emit
  correctly — but at least one of the newly-emitted junctions is
  U-shaped and creates a spurious grade-propagation chain from
  the runway HARD anchor back to primary E.

* **Specific spurious edge identified.**  The SE apron junction
  (current ``shape 94``, area 12 632 m², 12 ring vertices) is
  U-shaped.  Its all-pair edge ``v8 → v2`` (from (240, -904) to
  (305, -874), a 72 m straight segment) is **60 % outside the
  polygon** — the line cuts across the grass interior between
  the two prongs of the U.  The solver's ``_build_edges``
  treats it as a tight 72 m × 1.5 % = 1.08 m grade-cap edge and
  propagates runway altitudes back to primary E via this
  shortcut.

* **Per-axis grade rule is implemented correctly otherwise.**
  Solver graph at CYXY has 11 316 edges and **0 orphan edges**
  (every edge has both endpoints in the same shape).  No
  graph-distance / cross-pav_union propagation.  Per-axis rule
  for taxiway rects + within-shape all-pair for junction /
  apron is what the user prescribed and what the code does.

### User direction (this session)

> "If you're measuring grade across a complex polygon, if any
> point along a straight line between any two nodes crosses
> outside the polygon then that doesn't need to comply with
> grade, it would only be between edges.  We are most concerned
> with along axis of taxiways.  You should evaluate how this was
> being done before, when it was working."

> "We can't have airport specific code, must be generalized to
> work anywhere."

Memory ``feedback_grade_rules.md`` (user 2026-05-07) is
**unchanged direction** — apron/junction = 1.5 % all-pair
within the polygon, but the all-pair must be evaluated relative
to the polygon's interior (i.e. with a cover check), not as a
naïve Euclidean line.

### What was tried and reverted

* **Polygon-covers filter in ``_build_edges``** (and matching
  filter in the within-shape audit in ``elevation.py``).
  Conceptually correct: skip the all-pair edge if the segment
  isn't covered by the polygon.

  - Primary E moved 710.1 → 710.7 m (NOT to 717 — the chain
    length even after filtering is only ~471 m, slack 7.06 m).
  - **SPJC / SPLP** ``test_pavement_grade`` regressions: the
    spurious cross-grass edges were silently synchronising
    adjacent-shape altitudes at near-corners that aren't truly
    shared (4-5 m apart, different buckets).  Removing them
    exposed mid-edge step violations of 1.6 m between adjacent
    junction polygons.
  - Reverted; HEAD is back at 221 / 221.

### Recommended next direction — split the U-shaped junction

The polygon-covers filter at the solver edge level is treating a
symptom.  **The root cause is that the junction emission emits
ONE polygon for what is physically two pavement strips with a
grass gap between them.**

Per user (2026-05-13 already in code) the seam-anchor
architecture is the canonical layer for this kind of geometry.
Per user (this session) the fix must be generalised — no
airport-specific shortcuts.

Next agent: **investigate splitting U-shaped junction polygons
at the grass interior at emission time** (in
``junction_emit.py`` ``emit_junctions_and_finalize`` or its
``_decompose_polygon_with_holes`` helper).  The signal: a
junction polygon whose area is significantly less than its
convex-hull area (concavity ratio) OR whose own straight-line
across two opposite ring vertices exits the polygon by > N m.
Split into two polygons at the narrowest neck of the U.

This keeps the solver graph honest (all-pair within each
resulting polygon stays within actual pavement) without needing
a runtime cover check.  And it composes correctly with the
existing per-axis rule for taxi rects.

### Key data points for the next agent

CYXY primary E v0 — the canonical regression vertex:
* DEM altitude: 716.0 m
* Current solver result: 710.1 m
* Shortest grade-weighted path to a runway HARD anchor:
  `primary E (39, -597) → primary E LOW (115, -797) → stub E
  (126, -808) → stub E (211, -881) → primary E shape 2
  (240, -904) → junction 95 ring hops → runway (305, -874)`.
  6 → 9 nodes depending on whether the spurious 72 m edge in
  junction 94 is active.  Total geometric distance 450-471 m.
* User's mental model of the chain: via stub D + taxiway E at
  ~45° to 32L, then turns parallel to runway.  ~733 m+ chain.
  The current ~471 m chain shortcuts through the SE apron
  U-shape that the user considers "across grass".

CYXY runway profile pathology (separate concern, mentioned by
user):
* CIFP threshold elevations: 14R = 694.0, 32L = 706.2, 02 =
  694.3, 20 = 693.7, 14L = 694.3, 32R = 701.3.
* Runway 14R/32L's CIFP profile is grade-infeasible: 12.2 m
  rise over the runway length with the 02/20 crossing
  constraint, the regrade pass logs "anchor pair grade 1.58 %
  > 1.5 %".
* User direction (this session): 14R/32L should have a HARD
  anchor at the 02/20 crossing pinned to 02/20's altitude there
  (~694), then DEM-grade between the threshold and the crossing.
  ``runway_regrade.regrade_runway`` already accepts a
  ``seam_anchors`` list — runway-runway crossings just need to
  be computed before regrade and fed in.  Currently
  ``seam_anchors`` only gets tile-boundary seams (none at CYXY
  since it's single-tile).

These two issues compound — even with the grass-gap split fix,
primary E only reaches DEM-following altitude if the runway
itself is at DEM-following altitude (~710 m near the SE end of
14R/32L), not the CIFP-linear 704 m.

## How to verify the regression and the fix

```bash
# Baseline (HEAD): primary E v0 at 710.1 m, DEM is 716.0 m
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
layout = build_airport_pavement(
    "CYXY", "/Users/noah/X-Plane 12", compute_elevations=True)
for s in layout.shapes:
    if (s.role == "primary_parallel" and getattr(s, "ref", "") == "E"
            and s.altitude_high is not None
            and abs(s.altitude_high - 710.1) < 1):
        c = s.polygon.centroid
        if -800 < c.y < -500:
            print(f"primary E top: alt_high={s.altitude_high}"
                  f" alt_low={s.altitude_low}"
                  f" centroid=({c.x:.1f},{c.y:.1f})")
PY

# Expected after fix: alt_high ≈ 715-716 m (matching DEM,
# subject to the runway profile also being fixed to follow DEM
# between thresholds + crossings).
```

```bash
# Junction 94 (the offending U-shape) — verify its 60 % grass
# coverage and the v8 ↔ v2 spurious edge:
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
from shapely.geometry import LineString
layout = build_airport_pavement(
    "CYXY", "/Users/noah/X-Plane 12", compute_elevations=True)
apt_pav = layout._apt_pav_union
# Find the SE apron junction (concave U-shape, area > 10 000 m²,
# bounding box around (240..390, -1100..-700)):
for i, s in enumerate(layout.shapes):
    if s.role != "junction":
        continue
    c = s.polygon.centroid
    if c.is_empty:
        continue
    if not (280 < c.x < 340 and -1000 < c.y < -880):
        continue
    if s.polygon.area < 10000:
        continue
    inside = s.polygon.intersection(apt_pav).area
    print(f"junction {i} area={s.polygon.area:.0f} "
          f"inside_pav={100*inside/s.polygon.area:.1f}%")
    # The v8↔v2 spurious edge — its line should be ~60 % inside
    # the polygon (40 % across grass):
    coords = list(s.polygon.exterior.coords)
    if coords and coords[0] == coords[-1]: coords = coords[:-1]
    # Find vertices closest to (240, -904) and (305, -874):
    def near(target):
        return min(range(len(coords)),
                   key=lambda k: (coords[k][0]-target[0])**2
                               + (coords[k][1]-target[1])**2)
    a = near((240, -904))
    b = near((305, -874))
    line = LineString([coords[a], coords[b]])
    inside_line = s.polygon.intersection(line).length
    print(f"  v{a}-v{b} segment: total={line.length:.1f} m, "
          f"inside polygon={inside_line:.1f} m "
          f"({100*inside_line/line.length:.1f}%)")
    break
PY
```

## Tests + fixture state

* All 221 tests pass at HEAD (commit ``b698166``).
* SPLP -77 fixture (``tests/fixtures/SPLP_target_tile-13-77.osm``)
  was regenerated in this session to match the new stub-emit +
  canonical-runway-emit + apt.dat-driven stub-emit code.
* SPLP -77 baseline in ``tests/test_compare_target.py``: junction
  floor 6, total floor 150.
* `+60-140/` is an Ortho4XP tile directory leftover from
  development.  Not tracked by git; ignore or clean.

## Memory notes — read before continuing

* ``feedback_grade_rules.md`` (2026-05-07) — apron / junction
  grade rule.  All-pair within polygon, no radius cap, with the
  user's recent (2026-05-14) refinement that the all-pair only
  applies between pairs whose straight line stays inside the
  polygon.
* ``feedback_root_cause_only.md`` — fix root causes, ASK before
  band-aid post-process clean-up.  The U-shape split is a
  root-cause fix at junction emission time, not a band-aid.
* ``feedback_general_solutions.md`` — no airport-specific code.
  The grass-gap split must be a generic geometric test (e.g.
  concavity ratio threshold, or detect-and-split where a
  polygon's straight diagonal exits the polygon by > N m).
* ``project_target_osm.md`` — target OSM schema and SPJC / SPLP
  target fixture invariants.
* ``project_refactor_state.md`` — extraction-pattern recipe for
  pulling code out of the pavement monolith.

## Build commands (unchanged from prior sessions)

```bash
venv/bin/python3 -m pytest tests/ --tb=short -q

# Full per-airport build (visual review):
venv/bin/python3 - <<'PY'
import sys; sys.path.insert(0, "src")
from auto_patch.pipeline import build_airport_pavement
from O4_DEM_Utils import DEM as O4DEM
for icao in ("SPJC", "CYXY"):
    layout = build_airport_pavement(
        icao, "/Users/noah/X-Plane 12", compute_elevations=True)
    layout.to_osm(f"/tmp/{icao}.osm")
for tile_lat, tile_lon in [(-13, -77), (-13, -78)]:
    dem = O4DEM(tile_lat, tile_lon, fill_nodata="to zero")
    layout = build_airport_pavement(
        "SPLP", "/Users/noah/X-Plane 12", compute_elevations=True,
        tile_dem=dem,
        current_tile_lat=tile_lat,
        current_tile_lon=tile_lon)
    layout.to_osm(
        f"/tmp/SPLP_tile{tile_lat:+d}{tile_lon:+d}.osm")
PY
```
