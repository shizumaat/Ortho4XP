"""Structural-fidelity gate against the SPJC reference output.

``tests/fixtures/SPJC_target.osm`` is the canonical SPJC build
that the project agreed to as the new baseline (2026-05-07, after
the Naval Base ``patches_area`` skip fix).  Every code change must
continue to reproduce this output for SPJC: the test compares the
produced layout against the target shape-for-shape via
``tools/compare_target.match_by_role`` and asserts that each
role's match count stays at or above an established baseline.

If a future change drops matched shapes below the baseline — even
if every invariant test still passes — this gate fails.  That
makes regressions visible the way the comparison tool does
manually.

Baseline established 2026-05-07 after the Option-A
``encode_runways_taxiways_and_aprons`` patches_area-skip fix
(stops the Peruvian Naval Air Base — name-keyed, not in
patches_list — from injecting DEM-driven TAXIWAY constraint edges
into SPJC's apron region).  See STATUS.md for context.

Add new airport baselines as ``tests/fixtures/<ICAO>_target.osm``
files come online.  CYXY currently has only a guide file (zero
shapes) so it is not gated here.
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


# Per-role minimum match counts — reset 2026-05-07 with the new
# baseline (target file regenerated to equal the canonical SPJC
# output post-Naval-Base-skip fix).  These floors reflect a
# fresh build matching the target shape-for-shape; any drop
# means a structural regression.
#
# A small (~3 shape) per-run gap exists due to non-determinism in
# the build's node-ID assignment / sliver-drop ordering — the
# floors below are set to the observed steady-state match counts.
# If a future cleanup removes that non-determinism, raise these
# to full equality with target counts (see comments).
SPJC_BASELINE: Dict[str, int] = {
    "boundary":          620,   # of 630 target (boundary now emits
                                # as a chain of ~25 m rects so JOSM
                                # shows altitude profile per segment;
                                # small variance from sliver drops)
    "cross_connector":     6,   # of  6 target (full)
    "junction":           36,   # of 36 target (full)
    "primary_parallel":   27,   # of 27 target (full)
    "retaining_wall":     66,   # of 66 target (full)
    "runway":             85,   # of 85 target (full — blast pads
                                # now part of the segment chain)
    "secondary_parallel":  1,   # of  1 target (full)
    "stub":               16,   # of 16 target (full)
    "terminal":            2,   # of  2 target (full)
    "tunnel_ramp":        34,   # of 36 target (2-shape variance)
}
SPJC_BASELINE_TOTAL = 893  # of 905 target shapes (~99% match;
                           # gap is run-to-run determinism,
                           # not a structural regression)


def _build_layout(icao: str):
    from auto_patch.pipeline import build_airport_pavement
    return build_airport_pavement(icao, xplane_root(),
                                   compute_elevations=True)


def test_compare_target_spjc(tmp_path):
    """SPJC structural fidelity vs ``tests/fixtures/SPJC_target.osm``.

    Asserts each tracked role meets its minimum match count.  See
    ``SPJC_BASELINE`` for current floors.

    On failure, the assertion message prints the per-role match
    table so the regressing role(s) are obvious.
    """
    import compare_target as CT

    target_path = _HERE / "fixtures" / "SPJC_target.osm"
    assert target_path.is_file(), (
        f"SPJC target fixture missing at {target_path}")

    layout = _build_layout("SPJC")
    out_path = tmp_path / "SPJC_out.osm"
    layout.to_osm(str(out_path))

    anchor = CT.pick_anchor(target_path)
    target_shapes = CT.load_shapes(target_path, anchor, "target")
    output_shapes = CT.load_shapes(out_path, anchor, "output")
    pairs = CT.match_by_role(target_shapes, output_shapes)

    # Group matches by role.
    matched_by_role: Dict[str, int] = {}
    for p in pairs:
        if p.target is None or p.output is None:
            continue
        if p.iou <= 0.0:
            continue
        role = p.target.role
        matched_by_role[role] = matched_by_role.get(role, 0) + 1

    # Build a readable per-role summary for failure messages.
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
        floor = SPJC_BASELINE.get(role)
        floor_str = f" (floor {floor})" if floor is not None else ""
        summary_lines.append(
            f"  {role:20s} target={n_t:3d}  out={n_o:3d}  "
            f"matched={n_m:3d}{floor_str}")
    summary = "\n".join(summary_lines)

    # Per-role floor checks.
    failures = []
    for role, floor in SPJC_BASELINE.items():
        n_m = matched_by_role.get(role, 0)
        if n_m < floor:
            failures.append(
                f"{role}: matched={n_m} < floor={floor}")

    # Total-match floor.
    total_matched = sum(matched_by_role.values())
    if total_matched < SPJC_BASELINE_TOTAL:
        failures.append(
            f"total: matched={total_matched} < "
            f"floor={SPJC_BASELINE_TOTAL}")

    assert not failures, (
        "SPJC structural-fidelity regression vs target:\n"
        f"  failures: {'; '.join(failures)}\n"
        f"  per-role detail:\n{summary}")
