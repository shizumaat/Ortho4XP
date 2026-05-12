"""Cross-tile elevation parity for the per-tile slicing pipeline.

A cross-tile airport (e.g. SPLP, which spans tiles -13/-77 and
-13/-78) gets processed once per tile Ortho4XP renders.  The driver
passes ``tile_dem`` (the current build tile's DEM) and
``current_tile_lat/lon``; ``cut_layout_at_tile_boundaries`` keeps
only shapes in the current tile.  For the resulting per-tile patches
to meet smoothly at the cut, the pre-slice elevation calculation
must produce consistent values at the boundary regardless of which
tile is currently being built.

Regression: a coord/DEM mismatch in ``pipeline.py`` and
``elevation.py`` (user 2026-05-12) caused ``_sample_dem`` to be
called with ``tile_lat = floor(layout.anchor)`` while the DEM was
the driver-provided current build tile — for cross-tile airports
where anchor and current build tile differ, ``_sample_dem`` indexed
the wrong row of the DEM array.  Fixed by routing
``current_tile_lat/lon`` through to ``_compute_elevations`` and
``per_surface_solve`` whenever a ``tile_dem`` override is provided.
"""
from __future__ import annotations

import math
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))


def _xplane_root() -> str:
    return os.environ.get("O4_TEST_XPLANE_ROOT", "/Users/noah/X-Plane 12")


def _splp_apt_path() -> str:
    return os.path.join(
        _xplane_root(),
        "Custom Scenery",
        "SPLP Test",
        "Earth nav data",
        "apt.dat",
    )


def _have_splp() -> bool:
    return os.path.isfile(_splp_apt_path())


def _have_dem(tile_lat: int, tile_lon: int) -> bool:
    hem_ns = "S" if tile_lat < 0 else "N"
    hem_ew = "W" if tile_lon < 0 else "E"
    group_lat = (tile_lat // 10) * 10
    group_lon = (tile_lon // 10) * 10
    group_dir = (
        f"{'+' if group_lat >= 0 else '-'}{abs(group_lat):02d}"
        f"{'+' if group_lon >= 0 else '-'}{abs(group_lon):03d}"
    )
    fname = (
        f"{hem_ns}{abs(tile_lat):02d}{hem_ew}{abs(tile_lon):03d}.hgt"
    )
    return os.path.isfile(
        os.path.join(ROOT, "Elevation_data", group_dir, fname)
    )


_REQUIRES_FIXTURES = pytest.mark.skipif(
    not (_have_splp()
         and _have_dem(-13, -77)
         and _have_dem(-13, -78)),
    reason=(
        "Cross-tile parity test needs SPLP scenery + both -13/-77 "
        "and -13/-78 raw HGT DEM tiles on disk."
    ),
)


@_REQUIRES_FIXTURES
def test_cross_tile_build_completes_with_sane_elevations():
    """Build SPLP twice, once with each of its two tiles supplied as
    the current build tile, with the matching tile_dem.  Both builds
    must complete and produce elevations in a sensible range for
    SPLP (sea-level coastal Peru, ~50-200 m).

    Pre-fix, the second build would index the wrong DEM row in
    per_surface_solve, producing garbage altitudes (very large or
    very negative) for shapes whose elevation was solver-derived.
    """
    from auto_patch.pipeline import build_airport_pavement
    from O4_DEM_Utils import DEM

    xp_root = _xplane_root()
    # SPLP elevation range — Las Palmas sits ~50-200 m above sea
    # level.  Allow a generous +-200 m of slack.  The bug
    # would manifest as altitudes outside this range.
    ELEV_MIN, ELEV_MAX = -200.0, 400.0

    for tile_lat, tile_lon in [(-13, -77), (-13, -78)]:
        dem = DEM(tile_lat, tile_lon, fill_nodata="to zero")
        layout = build_airport_pavement(
            "SPLP", xp_root, compute_elevations=True,
            tile_dem=dem,
            current_tile_lat=tile_lat,
            current_tile_lon=tile_lon,
        )
        assert layout.shapes, (
            f"tile ({tile_lat},{tile_lon}): no shapes emitted")
        for s in layout.shapes:
            for attr in ("altitude", "altitude_high", "altitude_low"):
                v = getattr(s, attr, None)
                if v is None:
                    continue
                assert ELEV_MIN <= float(v) <= ELEV_MAX, (
                    f"tile ({tile_lat},{tile_lon}): role={s.role} "
                    f"{attr}={v} out of sane range "
                    f"[{ELEV_MIN}, {ELEV_MAX}]")
            if s.node_altitudes:
                for a in s.node_altitudes:
                    if a is None:
                        continue
                    assert ELEV_MIN <= float(a) <= ELEV_MAX, (
                        f"tile ({tile_lat},{tile_lon}): role={s.role} "
                        f"node_altitude={a} out of sane range "
                        f"[{ELEV_MIN}, {ELEV_MAX}]")


@_REQUIRES_FIXTURES
def test_cross_tile_cut_edge_elevations_consistent():
    """SPLP straddles lon=-77.  Build for each tile and compare
    elevations at vertices adjacent to the cut (within 10 m of the
    boundary line).  The two builds keep different pieces, but the
    near-cut elevations should be within ``TOL_M`` because the raw
    HGT data is identical along the shared edge row and Ortho4XP's
    smoothing fades to raw at the tile boundary.

    Currently a loose tolerance — the goal is to catch obvious
    DEM-indexing bugs, not subtle smoothing differences.  Tighten if
    a full MultiTileDEM refactor lands later.
    """
    from auto_patch.pipeline import build_airport_pavement
    from O4_DEM_Utils import DEM

    xp_root = _xplane_root()
    TOL_M = 1.0  # near-cut altitude agreement tolerance
    NEAR_CUT_M = 15.0  # meters of either side of the boundary

    builds = {}
    for tile_lat, tile_lon in [(-13, -77), (-13, -78)]:
        dem = DEM(tile_lat, tile_lon, fill_nodata="to zero")
        layout = build_airport_pavement(
            "SPLP", xp_root, compute_elevations=True,
            tile_dem=dem,
            current_tile_lat=tile_lat,
            current_tile_lon=tile_lon,
        )
        builds[(tile_lat, tile_lon)] = layout

    # Boundary at lon = -77.  Use anchor lat for cos(lat).
    lat0 = builds[(-13, -77)].anchor[0]
    cos0 = math.cos(math.radians(lat0))
    near_cut_deg = NEAR_CUT_M / (111195.0 * cos0)  # meters → deg lon

    def _near_cut_vertices(layout):
        """Yield (lat, lon, alt_m) tuples for every shape vertex
        within ``near_cut_deg`` of lon=-77."""
        out = []
        for s in layout.shapes:
            if s.polygon is None or s.polygon.is_empty:
                continue
            coords = list(s.polygon.exterior.coords)
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
            for i, (mx, my) in enumerate(coords):
                # Convert meter coords back to lat/lon via the
                # layout's helper.
                lat, lon = layout.m_to_ll(mx, my)
                if abs(lon + 77.0) > near_cut_deg:
                    continue
                # Resolve altitude per the canonical convention.
                alt = None
                if s.node_altitudes and i < len(s.node_altitudes):
                    alt = s.node_altitudes[i]
                elif s.altitude is not None:
                    alt = float(s.altitude)
                elif (s.altitude_high is not None
                        and s.altitude_low is not None):
                    alt = 0.5 * (float(s.altitude_high)
                                  + float(s.altitude_low))
                if alt is None:
                    continue
                out.append((lat, lon, float(alt), s.role))
        return out

    a = _near_cut_vertices(builds[(-13, -77)])
    b = _near_cut_vertices(builds[(-13, -78)])
    assert a, "no near-cut vertices in tile -13/-77 build"
    assert b, "no near-cut vertices in tile -13/-78 build"

    # Cut is the 10 m gap centered on lon=-77.  A-side vertices sit
    # at lon ≈ -76.99995 (5 m east); B-side at lon ≈ -77.00005 (5 m
    # west).  Match across the cut by lat only — the corresponding
    # pair (same source shape) has near-identical lat.  Lon
    # difference equals the cut-gap width (~10 m).
    # Tile-cut may insert new vertices at slightly different lats on
    # the two sides (depends on source-polygon edge geometry).  Match
    # within ~5 m of lat — same-source pair always pairs up; cross-
    # source pairs at similar lats also pair up but their altitudes
    # should still agree by virtue of shared neighbouring corners.
    MATCH_LAT_M = 5.0
    match_lat_deg = MATCH_LAT_M / 111195.0
    pairs = []
    for alat, alon, aalt, _ar in a:
        best_dlat = match_lat_deg
        best_balt = None
        for blat, blon, balt, _br in b:
            dlat = abs(blat - alat)
            if dlat < best_dlat:
                best_dlat = dlat
                best_balt = balt
        if best_balt is not None:
            pairs.append((alat, alon, aalt, best_balt,
                           abs(aalt - best_balt)))

    # Diagnostic dump (visible on test failure via repr).
    n_a, n_b, n_pairs = len(a), len(b), len(pairs)
    worst = max((p[4] for p in pairs), default=0.0)
    mean = (sum(p[4] for p in pairs) / n_pairs) if n_pairs else 0.0
    diag = (
        f"near-cut vertices: A={n_a} B={n_b}  matched (lat within "
        f"{MATCH_LAT_M:.1f}m): {n_pairs}.  worst |dz|={worst:.2f}m  "
        f"mean |dz|={mean:.2f}m"
    )
    assert n_pairs >= 4, f"too few matched pairs.  {diag}"
    assert worst <= TOL_M, (
        f"near-cut elevations diverge by > {TOL_M:.1f}m.  {diag}.  "
        f"Worst 3: "
        f"{sorted(pairs, key=lambda p: -p[4])[:3]}"
    )
