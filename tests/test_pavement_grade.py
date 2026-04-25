"""End-to-end grade validation for the pavement builder.

Skipped automatically unless an X-Plane install is available (the
builder needs CIFP + DEM tiles).  When run, it builds SPJC + SPLP,
writes the output OSM, then invokes ``tools.check_grade.run_checks``
to assert:

* No cross-shape proximity violations (shared corners agree on elev).
* No vertex-to-edge steps > 0.5 m (no visible drops between
  adjacent shapes — the user-reported "1 m drop" regression).
* Within-shape grade violations stay below a soft cap (the long-
  thin-apron-triangle case is a known limitation documented for
  follow-up; this test guards against new regressions).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_TOOLS = _HERE.parent / "tools"
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))


def _xplane_root() -> str:
    return os.environ.get("XPLANE_ROOT", "/Users/noah/X-Plane 12")


def _xplane_available() -> bool:
    root = _xplane_root()
    return (Path(root).is_dir()
            and (Path(root) / "Custom Data" / "CIFP").is_dir())


pytestmark = pytest.mark.skipif(
    not _xplane_available(),
    reason="X-Plane install not found (set XPLANE_ROOT to override)",
)


# Soft cap on within-shape violations.  Most are long-thin-triangle
# artefacts of the apron decomposition; we want to catch large
# regressions without failing on the known background level.
WITHIN_SHAPE_CAP = {"SPJC": 100, "SPLP": 30}


@pytest.mark.parametrize("icao", ["SPJC", "SPLP"])
def test_pavement_grade(tmp_path, icao):
    from O4_Airport_Pavement_Builder import build_airport_pavement
    import check_grade

    layout = build_airport_pavement(
        icao, _xplane_root(), compute_elevations=True)
    out = tmp_path / f"{icao}_test.osm"
    layout.to_osm(str(out))

    within, cross, steps = check_grade.run_checks(
        out,
        max_grade_pct=1.5,
        proximity_m=1.0,
        edge_search_m=5.0,
        edge_step_m=0.5,
        top_n=5,
    )
    # Hard fails — cross-shape continuity must be perfect.
    assert not cross, (
        f"{icao}: {len(cross)} cross-shape proximity violations "
        f"(shared corners disagree on elevation).  Worst: "
        f"{max(v.de_m for v in cross):.2f} m step.")
    assert not steps, (
        f"{icao}: {len(steps)} vertex-to-edge steps > 0.5 m.  "
        f"Worst: {max(s.step_m for s in steps):.2f} m step at "
        f"{steps[0].way_v.role}/{steps[0].way_v.ref or steps[0].way_v.wid}.")
    # Soft cap — log warning if exceeded but still fail to surface
    # regressions.
    cap = WITHIN_SHAPE_CAP[icao]
    assert len(within) <= cap, (
        f"{icao}: {len(within)} within-shape grade violations exceeds "
        f"soft cap {cap}.  Worst: {within[0].grade_pct:.2f}% over "
        f"{within[0].distance_m:.1f} m.")
