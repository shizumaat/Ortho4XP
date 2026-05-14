# Auto-Patch Status — 2026-05-14: stub-centering + runway-emit fixes

## TL;DR

Session continuation from 2026-05-13.  Triangle-prevention fix was
rolled back at user direction.  Five new fixes landed (all in this
session's uncommitted tree pending the commit at the end):

1. **Stub centering + min-runway-clearance** (`pavement/stubs.py`).
   Per user invariant: stubs never touch the runway, the rect is
   centered in the connector area with a junction on either side.
   Replaces the asymmetric 0.35 × gap pull-back with a 0.5 × gap
   midpoint center, unifies the inside-runway and outside-near-
   runway endpoint branches, and adds a MIN_RUNWAY_CLEARANCE_M
   (22 m) guardrail that shrinks the axis symmetrically when an
   endpoint falls inside the runway-adjacency snap window.

2. **Snap pav with runway carved out** (`pavement/stubs.py`).
   The corner snap and perpendicular-extension passes use
   `pav_union − runway_union.buffer(MIN_CLEARANCE)` so corners
   can't be pulled onto the runway boundary at snap time —
   the dominant residual leak in the v1 stub-centering attempt.

3. **apt.dat-driven runway-end stub emit** (`pipeline.py` +
   `pavement/stubs.py`).  `_emit_primary_parallel_runway_stubs`
   now accepts an `apt_centerlines` list and prefers it over OSM
   ways when available.  Pipeline passes the raw apt.dat
   linemerged-per-name polylines (NOT the `split_merged_centerline`
   output, which curve-skips at runway transitions and loses the
   endpoint).  At SPLP this replaces the unrefed
   convergence-point stub with a proper `A`-ref stub at the apt.dat
   A-taxi's actual runway terminus.

4. **Pavement-runway intersection tolerance widened**
   (`pipeline.py`).  `INTERSECTION_PROX_M` 3 m → 6 m so apt.dat
   pavement boundaries drawn up to 6 m inside (or outside) the
   runway rect register as seam-corner candidates.  Captures
   SPLP's north-end pavement corner (~5 m inside the runway).
   Added along-axis dedup (5 m) alongside the existing 5 m
   euclidean dedup so opposite-side same-transition pavement
   vertices collapse to one seam instead of producing a 4-m
   micro-segment that fails the grade gate at 23 %.

5. **Canonical runway-emit detection** (`layout.py` `to_osm`).
   When `seam_anchors` converts the entire runway chain to
   `node_altitudes`, the emit boundary now detects 4-corner
   shapes whose per-vertex altitudes follow the canonical
   `[H, L, L, H, H]` pattern (within 5 cm) and re-emits them as
   `altitude_high` / `altitude_low` (sloped) or `altitude=` (flat).
   Result: all CYXY (62/62) and SPJC (91/91) runway segments emit
   canonical tags; SPLP -77 emits 21/24 canonical (the 3 outliers
   are a tile-cut sliver triangle, a threshold-end segment with
   twisted altitudes, and a 5-corner seam-injected segment — none
   are canonical 4-corner rects).

**Reverted from prior session:**
- Triangle-prevention fix in `_extend_rect_corners_perpendicular`
  (user direction 2026-05-14: roll back).
- SPLP -77 fixture regenerated to match rolled-back + new fixes;
  baseline updated (junction floor 4 → 6, total 147 → 150).

**Test status: 221 / 221 pass.**

**Branch:** `dev` (3 commits ahead of `origin/dev`, plus the
in-flight changes summarised above pending this session's commit).

## Reach-down behaviour — pending investigation

User report 2026-05-14: at CYXY, way -10002 (an apron / boundary
junction near the SE runway end) appears to be over-clamped — DEM
calls for it to reach ~715 m and the adjacent runway segment way
-10067 should drop closer to ~700 m, but it's not getting there.
Investigation queued for next pass.

## What's still broken

* CYXY F-stub triangle (way -10016): user direction is to leave it
  alone for now; the triangle-prevention fix was rolled back.
* SPLP -77 has 3 non-canonical runway shapes (idx 5/6/7 — the
  south-threshold transition).  Not breaking the no-overlap or
  grade gates, but they emit `node_altitudes` rather than the
  canonical sloped-rect tags.

## Build / verify commands

```bash
venv/bin/python3 -m pytest tests/ --tb=short -q

# Builds for visual review:
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

## Memory notes

Existing memory entries remain authoritative.  Specifically:
- `feedback_root_cause_only` — fix root causes, ask before
  band-aid post-process clean-up.
- `feedback_shape_rules` — runway has vertex at every taxiway
  intersection (formalised here via the 6 m INTERSECTION_PROX_M
  and the apt.dat-driven runway-end stub emit).
- `feedback_general_solutions` — fixes work at any airport, no
  hardcoded ICAOs.

User direction 2026-05-14 (reinforced):
- Prefer apt.dat for taxi-network data; OSM only when apt.dat is
  completely missing taxiways.
- Stub rect sits centred in the stub area, never touches the
  runway; junction on each side fills the gap.
- Runway segments emit as `altitude=` (flat) or
  `altitude_high` / `altitude_low` (sloped 4-corner) — never as
  `node_altitudes` when the canonical pattern fits.
