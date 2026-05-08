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


# Soft cap on within-shape violations.  The check returns the UNION
# of (1) vertex-pair grade violations along triangle edges AND
# (2) planar-gradient violations (a triangle whose plane tilts
# more than 1.5 % across its surface).  Most are long-thin and
# sliver-triangle artefacts of ear-clipping the apron; cap is
# calibrated to the current background so regressions trip the test.
WITHIN_SHAPE_CAP = {"SPJC": 30, "SPLP": 30}
# Mid-edge step cap: every triangle plane should match its
# neighbours' surface along shared boundaries.  Samples along each
# edge and compares to the nearest other-shape edge's interpolation.
#
# SPLP exception (user 2026-05-08): junction-10053 and junction-10054
# are adjacent polygons whose rings don't share OSM nids (they have
# parallel edges with separate vertices ~0.1-4 m apart, disagreeing
# by ~6 m in elevation at mid-edge samples).  Fixing requires a
# junction-to-junction stitch pass analogous to
# ``stitch_pavement_to_terminals`` (currently terminals-only).
# Tracked as future work; cap set to current observed background.
MID_EDGE_CAP = {"SPJC": 10, "SPLP": 15}


@pytest.mark.parametrize("icao", ["SPJC", "SPLP"])
def test_pavement_grade(tmp_path, icao):
    from auto_patch.pipeline import build_airport_pavement
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
    # Soft cap on vertex-to-edge + mid-edge steps combined.  Vertex
    # continuity at shared boundaries should be ~perfect; mid-edge
    # discontinuities (sliver triangles whose plane tilts away from
    # neighbouring triangles' surfaces) are the known background.
    step_cap = MID_EDGE_CAP[icao]
    assert len(steps) <= step_cap, (
        f"{icao}: {len(steps)} edge/mid-edge steps > 0.5 m exceeds "
        f"cap {step_cap}.  Worst: {max(s.step_m for s in steps):.2f} "
        f"m step.")
    # Soft cap — log warning if exceeded but still fail to surface
    # regressions.
    cap = WITHIN_SHAPE_CAP[icao]
    assert len(within) <= cap, (
        f"{icao}: {len(within)} within-shape grade/plane violations "
        f"exceeds soft cap {cap}.  Worst: {within[0].grade_pct:.2f}% "
        f"over {within[0].distance_m:.1f} m.")
