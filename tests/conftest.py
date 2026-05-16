"""Pytest configuration shared by all tests.

* Adds the project's ``src/`` directory to ``sys.path`` so tests can
  import the O4 modules directly without installing the package.
* Provides the airport-discovery + ship-mode toggle helpers used by
  the integration test suites (overlap, junction invariants, grade,
  geometry, …) so airports are not hard-coded in any test file.

Environment variables (per user 2026-04-30):
* ``O4_TEST_TILE=lat,lon`` — discover every airport whose runways
  fall in that 1°×1° tile via the project's CIFP scanner.  This
  matches "all airports in the tile being built".
* ``O4_TEST_AIRPORTS=ICAO1,ICAO2,…`` — explicit ICAO list, takes
  precedence over ``O4_TEST_TILE`` when both are set.
* ``O4_SHIP_MODE=1`` — skip every integration test (used at
  shipping time when tests should not run).  All test modules
  collect normally but each item is marked skip.
* ``XPLANE_ROOT`` — X-Plane install path used by every test that
  needs CIFP / DEM / apt.dat data.  Defaults to
  ``/Users/noah/X-Plane 12``.

When neither ``O4_TEST_TILE`` nor ``O4_TEST_AIRPORTS`` is set, the
discovery returns an empty list and parametrised integration tests
collect zero items.  This is intentional: the project should not
ship a hidden hard-coded list of canonical airports — every run
must be explicit about what it tests against.
"""
import os
import sys
from typing import List, Optional

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.normpath(os.path.join(_HERE, "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


def xplane_root() -> str:
    return os.environ.get("XPLANE_ROOT", "/Users/noah/X-Plane 12")


def xplane_available() -> bool:
    root = xplane_root()
    return (os.path.isdir(root)
            and os.path.isdir(os.path.join(root, "Custom Data", "CIFP")))


# Baseline airports (user 2026-05-16): these run unconditionally on
# every invariant / grade / overlap test, providing CI coverage for
# the canonical geometry classes we need to support:
#
#   * SPJC — apron-heavy, multi-ref diagonal stubs, sloped runway,
#     hand-drawn target available for structural regression.
#   * SPLP — single runway with built-up threshold pad, multi-tile
#     output, sloped runway with seam anchors.  Target available.
#   * CYXY — multi-runway crossings (D crosses 14R/32L and 14L/32R),
#     chart-level E-D junction, parallel taxi E with diagonal stub
#     E and primary-parallel south extension into the apron.
#
# Adding KBNA / HECA to this list once Phase-1 baselines are stable
# would expand coverage to (4) terminal-heavy taxi network with
# multiple intersecting parallels, and (5) wide-runway desert
# airport with extensive aprons.
#
# Hand-drawn target comparison stays SPJC/SPLP-only — those are the
# regression-gate fixtures.  Invariant tests apply universally.
_BASELINE_AIRPORTS: tuple = ("SPJC", "SPLP", "CYXY")


def baseline_airports() -> tuple:
    """Return the canonical baseline airport list every invariant
    test should run against unconditionally.  See module-level
    ``_BASELINE_AIRPORTS`` for the rationale."""
    return _BASELINE_AIRPORTS


_AIRPORTS_CACHE: Optional[List[str]] = None


def airports_under_test() -> List[str]:
    """Resolved ICAO list for parametrised integration tests.

    Cached for the lifetime of the pytest session.  See module
    docstring for the env-var contract.
    """
    global _AIRPORTS_CACHE
    if _AIRPORTS_CACHE is not None:
        return _AIRPORTS_CACHE
    explicit = os.environ.get("O4_TEST_AIRPORTS", "").strip()
    if explicit:
        _AIRPORTS_CACHE = sorted({
            a.strip().upper()
            for a in explicit.split(",")
            if a.strip()})
        return _AIRPORTS_CACHE
    tile_str = os.environ.get("O4_TEST_TILE", "").strip()
    if tile_str:
        try:
            parts = [int(p.strip()) for p in tile_str.split(",")]
            assert len(parts) == 2
            lat, lon = parts
        except (ValueError, AssertionError):
            _AIRPORTS_CACHE = []
            return _AIRPORTS_CACHE
        _AIRPORTS_CACHE = _discover_airports_in_tile(lat, lon)
        return _AIRPORTS_CACHE
    _AIRPORTS_CACHE = []
    return _AIRPORTS_CACHE


def _discover_airports_in_tile(lat: int, lon: int) -> List[str]:
    """Return sorted ICAOs whose runways fall in the 1°×1° tile.

    Uses the same CIFP scanner the build pipeline uses
    (``O4_Cifp_Reader.discover_cifp_airports`` +
    ``parse_cifp_file`` + ``airport_in_tile``) to ensure tests run
    against the exact airport set the build pipeline would touch.
    """
    cifp_path = os.path.join(xplane_root(), "Custom Data", "CIFP")
    if not os.path.isdir(cifp_path):
        return []
    from auto_patch.cifp_reader import (
        discover_cifp_airports, parse_cifp_file, airport_in_tile)
    found: List[str] = []
    for icao, filepath in discover_cifp_airports(cifp_path).items():
        # Only test against true 4-letter ICAOs (mirror the build
        # pipeline's ICAO mode).
        if not (len(icao) == 4 and icao.isalpha()):
            continue
        runways = parse_cifp_file(filepath)
        if runways and airport_in_tile(runways, lat, lon):
            found.append(icao)
    return sorted(found)


def pytest_collection_modifyitems(config, items):
    """When ``O4_SHIP_MODE=1``, skip every collected test."""
    if os.environ.get("O4_SHIP_MODE", "0") == "1":
        skip_marker = pytest.mark.skip(
            reason="O4_SHIP_MODE=1 (tests disabled for shipping)")
        for item in items:
            item.add_marker(skip_marker)
