"""Tests for the apt.dat → dico_apt_entry adapter in O4_Auto_Patch.

Exercises ``_dico_from_apt_dat`` and ``xplane_root_from_cifp_path``
without requiring the full ``generate_airport_surface_patches``
pipeline (which has too many external dependencies for a unit
test).  The real SPJC apt.dat file is too large to ship as a
fixture, so these tests use the hand-crafted
``tests/fixtures/synthetic_apt.dat`` that backs the main parser
tests, loaded via the real ``O4_Apt_Dat_Reader`` and fed through
the adapter.

The synthetic fixture has:
  * ZZZZ (Test Airport One): 1 runway, 3 pavements (square,
    bezier-rect, nested-hole), 1 boundary
  * YYYY: 1 runway, 0 pavements (used for negative cases)
"""
import os

import pytest

import O4_Apt_Dat_Reader as APR
import O4_Auto_Patch as AP


_FIXTURE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "fixtures",
    "synthetic_apt.dat",
)


# Minimal fake tile with just the attributes the adapter touches.
class _FakeTile:
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon


# ──────────────────────────────────────────────────────────────────────
# xplane_root_from_cifp_path
# ──────────────────────────────────────────────────────────────────────
class TestXplaneRootFromCifpPath:
    def test_standard_custom_data_layout(self, tmp_path):
        """A real X-Plane-style directory tree: Custom Data/CIFP
        under an X-Plane root.  xplane_root should come back as
        the root.
        """
        xp = tmp_path / "X-Plane 12"
        cifp = xp / "Custom Data" / "CIFP"
        cifp.mkdir(parents=True)
        (xp / "Custom Scenery").mkdir()
        assert AP.xplane_root_from_cifp_path(str(cifp)) == str(xp)

    def test_alternative_resources_layout(self, tmp_path):
        """Laminar's default scenery CIFP lives under Resources/."""
        xp = tmp_path / "X-Plane 12"
        cifp = xp / "Custom Data" / "CIFP"
        cifp.mkdir(parents=True)
        (xp / "Resources").mkdir()
        assert AP.xplane_root_from_cifp_path(str(cifp)) == str(xp)

    def test_returns_none_for_missing_path(self):
        assert AP.xplane_root_from_cifp_path(None) is None
        assert AP.xplane_root_from_cifp_path("") is None

    def test_returns_none_when_root_doesnt_look_right(self, tmp_path):
        """If the derived root has neither Custom Scenery nor
        Resources, reject it — the path doesn't look like an
        X-Plane install.
        """
        weird = tmp_path / "not_xplane" / "Custom Data" / "CIFP"
        weird.mkdir(parents=True)
        assert AP.xplane_root_from_cifp_path(str(weird)) is None


# ──────────────────────────────────────────────────────────────────────
# _dico_from_apt_dat
# ──────────────────────────────────────────────────────────────────────
class TestDicoFromAptDat:
    def setup_method(self):
        self.apt = APR.load_airport(_FIXTURE, "ZZZZ")
        assert self.apt is not None
        # Fake tile positioned so ZZZZ's data (-12°, -77°) is near
        # the tile origin.  With tile.lat = -13 and tile.lon = -78,
        # the pavements at lat ~ -12.003 should have tile-relative
        # y ~ 0.997 (positive, close to 1).
        self.tile = _FakeTile(lat=-13, lon=-78)
        self.fallback_dico = {
            "boundary": None,
            "apron": [],
            "taxiway": [],
            "runway": [],
            "hangar": [],
            "key_type": "icao",
        }

    def test_apron_is_set_as_tuple(self):
        new_dico = AP._dico_from_apt_dat(
            self.apt, self.tile, self.fallback_dico)
        apron = new_dico["apron"]
        assert isinstance(apron, tuple)
        assert len(apron) == 2
        # The first element is a shapely geometry.
        assert hasattr(apron[0], "area")
        assert apron[0].area > 0
        # The second element is the (empty) way-id list.
        assert apron[1] == []

    def test_apron_is_tile_relative(self):
        """The ZZZZ pavements are around lat -12.00 (tile-relative
        y ≈ 1.0), lon -77.10 (tile-relative x ≈ 0.9).  With tile
        origin (lat=-13, lon=-78) the pavement bounds should land
        in a small neighbourhood of those values.  The synthetic
        fixture spans lat -12.003 to -11.999 so tile-relative y can
        slightly exceed 1.0.
        """
        new_dico = AP._dico_from_apt_dat(
            self.apt, self.tile, self.fallback_dico)
        apron_geom = new_dico["apron"][0]
        minx, miny, maxx, maxy = apron_geom.bounds
        assert 0.85 < minx < 0.95, f"minx={minx} not in tile range"
        assert 0.85 < maxx < 0.95, f"maxx={maxx} not in tile range"
        assert 0.99 < miny < 1.01, f"miny={miny} not in tile range"
        assert 0.99 < maxy < 1.01, f"maxy={maxy} not in tile range"

    def test_taxiway_is_empty(self):
        """Commit 8 sends all pavement through the apron path."""
        new_dico = AP._dico_from_apt_dat(
            self.apt, self.tile, self.fallback_dico)
        twy = new_dico["taxiway"]
        assert isinstance(twy, tuple)
        # The geometry should be empty (Polygon()).
        assert twy[0].is_empty

    def test_runway_rectangle_built(self):
        """The runway record in the fixture (09/27) should produce
        a rectangle in the runway field.
        """
        new_dico = AP._dico_from_apt_dat(
            self.apt, self.tile, self.fallback_dico)
        rwy = new_dico["runway"]
        assert isinstance(rwy, tuple)
        assert hasattr(rwy[0], "area")
        assert rwy[0].area > 0
        # Width is 45 m, runway is ~1 km long; rectangle area should
        # be in the ballpark of 45 × 1000 = 45,000 m².  In
        # tile-relative degrees² this is tiny; the assertion is
        # "positive and not absurd".
        assert rwy[0].area < 1e-3

    def test_boundary_is_set_from_aptdat(self):
        new_dico = AP._dico_from_apt_dat(
            self.apt, self.tile, self.fallback_dico)
        bnd = new_dico["boundary"]
        assert bnd is not None
        assert hasattr(bnd, "area")
        assert bnd.area > 0

    def test_fallback_fields_preserved(self):
        """Fields that apt.dat doesn't populate should fall through
        unchanged from fallback_dico.
        """
        self.fallback_dico["key_type"] = "icao"
        self.fallback_dico["hangar"] = ["OSM_HANGAR_DATA"]
        new_dico = AP._dico_from_apt_dat(
            self.apt, self.tile, self.fallback_dico)
        assert new_dico["key_type"] == "icao"
        assert new_dico["hangar"] == ["OSM_HANGAR_DATA"]

    def test_airport_with_no_pavements_leaves_apron_untouched(self):
        """YYYY in the fixture has only a runway, no row-110
        pavements.  The adapter should not set an apron key from
        empty data.
        """
        yyyy = APR.load_airport(_FIXTURE, "YYYY")
        fb = {"boundary": None, "apron": "FALLBACK_APRON",
              "taxiway": "FALLBACK_TAXIWAY",
              "runway": "FALLBACK_RUNWAY"}
        new_dico = AP._dico_from_apt_dat(yyyy, self.tile, fb)
        # apron should still be the fallback.
        assert new_dico["apron"] == "FALLBACK_APRON"
        # runway should be the new rectangle from row 100.
        assert new_dico["runway"] != "FALLBACK_RUNWAY"
        # taxiway is always replaced with an empty tuple.
        assert new_dico["taxiway"][0].is_empty
