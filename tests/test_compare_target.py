"""Structural-fidelity gates against the reference fixture outputs.

``tests/fixtures/SPJC_target.osm`` and ``tests/fixtures/SPLP_target.osm``
are the canonical builds (regenerated 2026-05-13 with the seam-anchor
+ diagonal-stub trapezoid pipeline).  Every code change must continue
to reproduce these outputs: the tests compare the produced layout
against each target shape-for-shape via
``tools/compare_target.match_by_role`` and assert that each role's
match count stays at or above an established baseline.

If a future change drops matched shapes below the baseline — even
if every invariant test still passes — these gates fail.  That makes
regressions visible the way the comparison tool does manually.

Baseline reset 2026-05-13 after:
  * diagonal V3-style stub trapezoid emission (Approach A) +
    ``db_local ≥ 15°`` digit→STUB classifier tightening.
  * Seam-anchor architecture: ``split_pavement_at_seams`` +
    ``apply_seam_dem_anchors`` + Stage A runway regrade +
    unified-Jacobi seam-HARD override (seam wins).
  * Tile-cut bridge polygons removed.
  * ``_resample_node_altitudes_nn`` upgraded to edge interpolation
    (cut-edge vertices use linear gradient of the underlying old
    edge instead of nearest-neighbour).

Add new airport baselines as ``tests/fixtures/<ICAO>_target.osm``
files come online.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict

import pytest

from conftest import xplane_available, xplane_root


_HERE = Path(__file__).resolve().parent
_TOOLS = _HERE.parent / "tools"
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))


pytestmark = pytest.mark.skipif(
    not xplane_available(),
    reason="X-Plane install not found (set XPLANE_ROOT to override)",
)


# Per-role floors are set ~5 % below target counts to absorb the
# run-to-run non-determinism in node-ID assignment / sliver-drop
# ordering.  A regression that drops more than ~5 % of any role's
# shapes vs target trips the gate.
SPJC_BASELINE: Dict[str, int] = {
    "boundary":          595,   # of 603 target
    "cross_connector":     8,   # of   9 target
    "junction":           38,   # of  40 target
    "primary_parallel":   26,   # of  28 target
    "retaining_wall":     63,   # of  66 target
    "runway":             86,   # of  91 target
    "secondary_parallel":  1,   # of   1 target
    "stub":               19,   # of  20 target
    "terminal":            2,   # of   2 target
    "tunnel_ramp":        34,   # of  36 target
}
SPJC_BASELINE_TOTAL = 872  # of 896 target

SPLP_BASELINE: Dict[str, int] = {
    "boundary":          117,   # of 123 target
    "junction":            4,   # of   4 target
    "primary_parallel":    2,   # of   2 target
    "runway":             22,   # of  23 target
    "stub":                3,   # of   3 target
}
SPLP_BASELINE_TOTAL = 148  # of 155 target


def _build_layout(icao: str):
    from auto_patch.pipeline import build_airport_pavement
    return build_airport_pavement(icao, xplane_root(),
                                   compute_elevations=True)


def _run_compare(tmp_path: Path, icao: str,
                 baseline: Dict[str, int],
                 baseline_total: int) -> None:
    import compare_target as CT

    target_path = _HERE / "fixtures" / f"{icao}_target.osm"
    assert target_path.is_file(), (
        f"{icao} target fixture missing at {target_path}")

    layout = _build_layout(icao)
    out_path = tmp_path / f"{icao}_out.osm"
    layout.to_osm(str(out_path))

    anchor = CT.pick_anchor(target_path)
    target_shapes = CT.load_shapes(target_path, anchor, "target")
    output_shapes = CT.load_shapes(out_path, anchor, "output")
    pairs = CT.match_by_role(target_shapes, output_shapes)

    matched_by_role: Dict[str, int] = {}
    for p in pairs:
        if p.target is None or p.output is None:
            continue
        if p.iou <= 0.0:
            continue
        role = p.target.role
        matched_by_role[role] = matched_by_role.get(role, 0) + 1

    target_counts: Dict[str, int] = {}
    for s in target_shapes:
        target_counts[s.role] = target_counts.get(s.role, 0) + 1
    output_counts: Dict[str, int] = {}
    for s in output_shapes:
        output_counts[s.role] = output_counts.get(s.role, 0) + 1

    summary_lines = []
    for role in sorted(set(target_counts) | set(output_counts)):
        n_t = target_counts.get(role, 0)
        n_o = output_counts.get(role, 0)
        n_m = matched_by_role.get(role, 0)
        floor = baseline.get(role)
        floor_str = f" (floor {floor})" if floor is not None else ""
        summary_lines.append(
            f"  {role:20s} target={n_t:3d}  out={n_o:3d}  "
            f"matched={n_m:3d}{floor_str}")
    summary = "\n".join(summary_lines)

    failures = []
    for role, floor in baseline.items():
        n_m = matched_by_role.get(role, 0)
        if n_m < floor:
            failures.append(
                f"{role}: matched={n_m} < floor={floor}")
    total_matched = sum(matched_by_role.values())
    if total_matched < baseline_total:
        failures.append(
            f"total: matched={total_matched} < "
            f"floor={baseline_total}")

    assert not failures, (
        f"{icao} structural-fidelity regression vs target:\n"
        f"  failures: {'; '.join(failures)}\n"
        f"  per-role detail:\n{summary}")


def test_compare_target_spjc(tmp_path):
    """SPJC structural fidelity vs ``tests/fixtures/SPJC_target.osm``.

    See ``SPJC_BASELINE`` for current per-role floors.
    """
    _run_compare(tmp_path, "SPJC",
                 SPJC_BASELINE, SPJC_BASELINE_TOTAL)


def test_compare_target_splp(tmp_path):
    """SPLP structural fidelity vs ``tests/fixtures/SPLP_target.osm``.

    Cross-tile airport (spans tiles -13/-77 and -13/-78); the
    standalone build (no ``tile_dem`` override) used to generate
    the fixture exercises the same pipeline as the in-tile builds
    that the test_tile_cut_parity tests run with ``tile_dem`` set.
    """
    _run_compare(tmp_path, "SPLP",
                 SPLP_BASELINE, SPLP_BASELINE_TOTAL)
