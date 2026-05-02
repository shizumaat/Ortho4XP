"""Structural-fidelity gate against hand-crafted target OSMs.

A lot of work went into ``tests/fixtures/SPJC_target.osm``; the
production output should match it closely.  Earlier tests in this
suite check geometric invariants (no overlap, vertex count
bounded, junction boundary near centerline) which are useful but
have proven to be poor proxies for "is the output structurally
correct."  This test compares the produced layout against the
target shape-for-shape via ``tools/compare_target.match_by_role``
and asserts that each role's match count stays at or above an
established baseline.

If a future change is going to drop matched shapes below the
baseline — even if every invariant test still passes — this gate
fails.  That makes regressions visible the way the comparison
tool does manually.

Baseline established 2026-05-01 (post-revert of Phases B.1 + C):
SPJC produces 71/104 target shapes matched.  The legacy
absorption rule + Apr 29 reference state.

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


# Per-role minimum match counts — established 2026-05-01 by
# running ``compare_target`` against the post-revert HEAD.  Test
# fails if any role's matched count falls below its baseline.
# To intentionally raise a baseline (output improved), update
# both the count and the date stamp here.
SPJC_BASELINE: Dict[str, int] = {
    # 2026-05-02: junction floor lowered 32 → 30 across two phases.
    # Phase 1 (32 → 31) — Rule 1 (junction-runway 1:1 sharing) snaps
    # runway-near junction vertices to the nearest runway segment
    # endpoint; one junction deformed enough to drop below IoU
    # threshold.
    # Phase 2 (31 → 30) — Rule 3 (axis-aligned cut lines) routes
    # decomposition cuts along the runway axis instead of hole-MRR;
    # one further junction's cut topology shifted off-target.
    # Both intentional per user 2026-05-01 — invariants supersede
    # individual ground-truth matches.
    "junction":           30,   # of 43 target
    "primary_parallel":   20,   # of 29 target
    "stub":               15,   # of 15 target (full match)
    "terminal":            2,   # of  2 target (full match)
    "cross_connector":     2,   # of  6 target
    # secondary_parallel:  0/4  — not yet detected; track but don't gate.
    # apron:               0/3  — not yet emitted; track but don't gate.
    # runway:              0/2  — over-segmented (73 segments vs 2 target);
    #                              CIFP-driven segmentation is correct
    #                              behaviour, the target keeps single rects.
}
SPJC_BASELINE_TOTAL = 69  # of 104 target shapes (sum of above + ...)


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
