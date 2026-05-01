"""Tests for ``O4_Cifp_Reader.xplane_root_from_cifp_path``.

Standalone path-resolution tests; require no X-Plane install or
fixture data.
"""
from auto_patch import cifp_reader as AP


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

