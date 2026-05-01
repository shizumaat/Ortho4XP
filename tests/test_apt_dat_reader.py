"""Unit tests for O4_Apt_Dat_Reader.

Uses ``tests/fixtures/synthetic_apt.dat`` which contains two
hand-crafted airport blocks:

* ``ZZZZ`` — Test Airport One.  One runway, three pavements (a
  plain square, a square with one Bezier corner, and a square with
  an interior hole) and a boundary polygon.
* ``YYYY`` — minimal stub used to verify the parser stops at the
  next airport header.

The parser is exercised against this fixture so the tests don't
depend on the user's X-Plane install.  A separate test module can
do integration-style checks against the real SPJC data later.
"""
import os

import pytest
from shapely.geometry import Polygon

from auto_patch import apt_dat_reader as APR


_FIXTURE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "fixtures",
    "synthetic_apt.dat",
)


# ──────────────────────────────────────────────────────────────────────
# load_airport — top-level
# ──────────────────────────────────────────────────────────────────────
class TestLoadAirport:
    def test_loads_test_airport_one(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        assert apt is not None
        assert apt.icao == "ZZZZ"
        assert apt.name.startswith("Test Airport One")
        assert apt.reference_elev_ft == 100
        assert apt.reference_elev_m == pytest.approx(30.48, abs=0.01)
        assert apt.source_path == _FIXTURE

    def test_loads_other_airport_independently(self):
        """The parser must stop at the next ``1`` header — YYYY's data
        must not bleed into ZZZZ and vice versa.
        """
        zzzz = APR.load_airport(_FIXTURE, "ZZZZ")
        yyyy = APR.load_airport(_FIXTURE, "YYYY")
        assert zzzz is not None and yyyy is not None
        assert len(zzzz.runways) == 1
        assert len(yyyy.runways) == 1
        assert zzzz.runways[0].desig_a == "09"
        assert yyyy.runways[0].desig_a == "18"
        assert len(yyyy.pavements) == 0    # YYYY has no row 110

    def test_unknown_icao_returns_none(self):
        assert APR.load_airport(_FIXTURE, "AAAA") is None

    def test_missing_file_returns_none(self):
        assert APR.load_airport("/nonexistent/path/apt.dat", "ZZZZ") is None

    def test_case_insensitive_icao(self):
        apt = APR.load_airport(_FIXTURE, "zzzz")
        assert apt is not None
        assert apt.icao == "ZZZZ"


# ──────────────────────────────────────────────────────────────────────
# Runway parsing
# ──────────────────────────────────────────────────────────────────────
class TestRunwayParsing:
    def test_runway_fields(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        rwy = apt.runways[0]
        assert rwy.desig_a == "09"
        assert rwy.desig_b == "27"
        assert rwy.lat_a == pytest.approx(-12.0)
        assert rwy.lon_a == pytest.approx(-77.1)
        assert rwy.lat_b == pytest.approx(-12.0)
        assert rwy.lon_b == pytest.approx(-77.09)
        assert rwy.width_m == 45.0
        assert rwy.surface_code == 1
        assert rwy.displaced_a_m == 0.0
        assert rwy.displaced_b_m == 0.0


# ──────────────────────────────────────────────────────────────────────
# Pavement parsing
# ──────────────────────────────────────────────────────────────────────
class TestPavementParsing:
    def test_three_pavements_emitted(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        assert len(apt.pavements) == 3

    def test_pavement_names_preserved(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        names = [p.name for p in apt.pavements]
        assert "SQUARE" in names
        assert "BEZIER_RECT" in names
        assert "NESTED_HOLE" in names

    def test_simple_square_polygon(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        pav = next(p for p in apt.pavements if p.name == "SQUARE")
        # No Bezier curves → exactly 4 unique vertices in the
        # exterior ring (5 with the closing repeat).
        coords = list(pav.polygon.exterior.coords)
        assert len(coords) == 5
        assert coords[0] == coords[-1]
        # The polygon is small but should have non-zero area.
        assert pav.polygon.area > 0
        assert pav.surface_code == 1

    def test_bezier_pavement_has_extra_vertices(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        pav = next(p for p in apt.pavements if p.name == "BEZIER_RECT")
        coords = list(pav.polygon.exterior.coords)
        # 4 base nodes + Bezier subdivision means more than 5 total.
        assert len(coords) > 5
        assert pav.polygon.is_valid
        assert pav.polygon.area > 0
        assert pav.surface_code == 2

    def test_nested_hole_polygon(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        pav = next(p for p in apt.pavements if p.name == "NESTED_HOLE")
        # Outer ring + 1 interior ring (the hole).
        assert len(pav.polygon.interiors) == 1
        # The hole reduces total area but doesn't make the polygon empty.
        assert pav.polygon.area > 0
        assert pav.polygon.is_valid


# ──────────────────────────────────────────────────────────────────────
# Boundary parsing
# ──────────────────────────────────────────────────────────────────────
class TestBoundaryParsing:
    def test_boundary_polygon_present(self):
        apt = APR.load_airport(_FIXTURE, "ZZZZ")
        assert apt.boundary is not None
        assert isinstance(apt.boundary, Polygon)
        assert apt.boundary.area > 0


# ──────────────────────────────────────────────────────────────────────
# Bezier interpolation primitives
# ──────────────────────────────────────────────────────────────────────
class TestBezierInterpolation:
    def test_quadratic_bezier_endpoints(self):
        pts = APR._quadratic_bezier((0, 0), (1, 1), (2, 0), n_segments=4)
        assert len(pts) == 5
        assert pts[0] == (0, 0)
        assert pts[-1] == (2, 0)

    def test_quadratic_bezier_midpoint(self):
        # Quadratic Bezier midpoint is at (P0 + 2*P1 + P2) / 4 for t=0.5
        pts = APR._quadratic_bezier((0, 0), (10, 10), (20, 0), n_segments=2)
        mid = pts[1]
        assert mid[0] == pytest.approx(10.0)
        assert mid[1] == pytest.approx(5.0)

    def test_cubic_bezier_endpoints(self):
        pts = APR._cubic_bezier((0, 0), (1, 1), (2, 1), (3, 0),
                                n_segments=8)
        assert len(pts) == 9
        assert pts[0] == (0, 0)
        assert pts[-1] == (3, 0)

    def test_mirror_through_anchor(self):
        assert APR._mirror((1, 2), (0, 0)) == (-1, -2)
        assert APR._mirror((5, 5), (3, 3)) == (1, 1)


# ──────────────────────────────────────────────────────────────────────
# find_airport_apt_dat — search priority
# ──────────────────────────────────────────────────────────────────────
class TestFindAirportAptDat:
    def test_returns_none_when_root_missing(self):
        assert APR.find_airport_apt_dat("/nonexistent", "SPJC") is None

    def test_returns_none_when_icao_missing(self):
        assert APR.find_airport_apt_dat("/tmp", "") is None

    def test_finds_per_airport_pack_first(self, tmp_path):
        """A per-airport Custom Scenery pack should win over the
        Global Airports pack.  Build a fake X-Plane root with both.
        """
        xp = tmp_path / "X-Plane 12"
        cs = xp / "Custom Scenery"

        # 1) Per-airport pack containing ZZZZ
        per_apt_pack = cs / "ZZZZ Per Airport Pack" / "Earth nav data"
        per_apt_pack.mkdir(parents=True)
        per_apt_path = per_apt_pack / "apt.dat"
        per_apt_path.write_text(
            "A\n"
            "1    100 0 0 ZZZZ Per-Airport Test\n"
            "100 45.00 1 0 0.25 1 1 0 09 -12 -77 0 60 0 0 0 0 27 -12 -77.01 0 60 0 0 0 0\n",
            encoding="utf-8")

        # 2) Global Airports pack also containing ZZZZ
        global_pack = cs / "Global Airports" / "Earth nav data"
        global_pack.mkdir(parents=True)
        global_path = global_pack / "apt.dat"
        global_path.write_text(
            "A\n"
            "1    100 0 0 ZZZZ Global Test\n"
            "100 45.00 1 0 0.25 1 1 0 09 -12 -77 0 60 0 0 0 0 27 -12 -77.01 0 60 0 0 0 0\n",
            encoding="utf-8")

        found = APR.find_airport_apt_dat(str(xp), "ZZZZ")
        assert found == str(per_apt_path)

    def test_falls_back_to_global(self, tmp_path):
        """If no per-airport pack contains the ICAO, fall back to the
        Global Airports pack.
        """
        xp = tmp_path / "X-Plane 12"
        cs = xp / "Custom Scenery"

        # Per-airport pack contains a DIFFERENT airport
        other_pack = cs / "Other Pack" / "Earth nav data"
        other_pack.mkdir(parents=True)
        (other_pack / "apt.dat").write_text(
            "A\n"
            "1    100 0 0 AAAA Other Airport\n",
            encoding="utf-8")

        # Global has ZZZZ
        global_pack = cs / "Global Airports" / "Earth nav data"
        global_pack.mkdir(parents=True)
        global_path = global_pack / "apt.dat"
        global_path.write_text(
            "A\n"
            "1    100 0 0 ZZZZ Global Test\n",
            encoding="utf-8")

        found = APR.find_airport_apt_dat(str(xp), "ZZZZ")
        assert found == str(global_path)

    def test_falls_back_to_default(self, tmp_path):
        """If neither Custom Scenery contains the ICAO, try the
        default scenery's apt.dat.
        """
        xp = tmp_path / "X-Plane 12"
        default_pack = (xp / "Resources" / "default scenery"
                        / "default apt dat" / "Earth nav data")
        default_pack.mkdir(parents=True)
        default_path = default_pack / "apt.dat"
        default_path.write_text(
            "A\n"
            "1    100 0 0 ZZZZ Default Test\n",
            encoding="utf-8")

        # Custom Scenery dir doesn't even exist
        found = APR.find_airport_apt_dat(str(xp), "ZZZZ")
        assert found == str(default_path)

    def test_returns_none_when_icao_nowhere(self, tmp_path):
        xp = tmp_path / "X-Plane 12"
        (xp / "Custom Scenery").mkdir(parents=True)
        assert APR.find_airport_apt_dat(str(xp), "ZZZZ") is None
