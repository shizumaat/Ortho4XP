"""Unit tests for auto_patch.pavement.absorption (the hot spot).

The module's docstring warns: "Every recent session that modified
this rule introduced a regression."  These tests pin the rule's
behavior on synthetic inputs so future edits land tested.

The rule (verbatim from the module docstring):

  A taxi rect (primary_parallel / secondary_parallel / stub /
  cross_connector) has a long edge "covered" by junction-class
  pavement when that pavement extends past the long edge for
  >= 10 % of the rect's axial length.  When EITHER long edge is
  covered to >= 10 %, the COVERED portion is absorbed; the
  non-covered portion stays as a (possibly shorter) rect.
  Partial absorption — never absorb the WHOLE rect when only
  part is covered.

  - Probe at 5 m axial steps; each step is "adjacent" if EITHER
    outside-long-edge sample is inside ``junction_pav``.
  - Probe distance = 5 m beyond the rect's long edge.
  - Adjacent run must be >= 10 % of axial steps to absorb.
  - Kept fragments < 30 m are dropped.

Tests use a horizontal rect convention:
    c0 = (0,    +half_w)   ┐ short edge at x=0
    c3 = (0,    -half_w)   ┘
    c1 = (L,    +half_w)   ┐ short edge at x=L
    c2 = (L,    -half_w)   ┘
Axis is along +x; "left" probe is +y, "right" probe is -y.
"""
from shapely.geometry import LineString, Polygon

from auto_patch.layout import (
    ROLE_APRON,
    ROLE_CROSS_CONNECTOR,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RUNWAY,
    ROLE_SECONDARY_PARALLEL,
    ROLE_STUB,
    ROLE_TERMINAL,
)
from auto_patch.pavement.absorption import (
    _drop_primary_parallels_embedded_in_pavement,
    drop_primary_parallels_embedded_in_pavement,
)


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def _rect(length=100.0, half_w=10.0):
    """Horizontal taxi rect from (0,0) to (length, 0).

    Returns (polygon, axis_linestring) ready to drop into a
    ``taxi_rects`` tuple as ``(rect, axis, role, ref)``.
    """
    poly = Polygon([
        (0.0, +half_w),       # c0
        (length, +half_w),    # c1
        (length, -half_w),    # c2
        (0.0, -half_w),       # c3
    ])
    axis = LineString([(0.0, 0.0), (length, 0.0)])
    return poly, axis


def _apron(x_lo, x_hi, y_lo, y_hi):
    """Axis-aligned apron polygon."""
    return Polygon([
        (x_lo, y_lo),
        (x_hi, y_lo),
        (x_hi, y_hi),
        (x_lo, y_hi),
    ])


def _kept_axes(result):
    """Helper: extract surviving rect axis spans as (xmin, xmax)."""
    out = []
    for rect, axis, role, ref in result:
        coords = list(axis.coords)
        xs = [c[0] for c in coords]
        out.append((min(xs), max(xs)))
    return out


# ──────────────────────────────────────────────────────────────────────
# Trivial pass-through cases
# ──────────────────────────────────────────────────────────────────────
def test_no_apron_returns_input_unchanged():
    """When apt_pav_union is None or empty, every rect passes through."""
    poly, axis = _rect()
    rects = [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")]
    assert _drop_primary_parallels_embedded_in_pavement(rects, None) is rects
    empty = Polygon()
    assert _drop_primary_parallels_embedded_in_pavement(rects, empty) is rects


def test_non_sloping_roles_pass_through():
    """The rule only applies to sloping rects.  APRON, TERMINAL,
    RUNWAY, JUNCTION roles must pass through untouched even when
    fully embedded."""
    big_apron = _apron(-50.0, 200.0, -50.0, 50.0)
    poly, axis = _rect(length=80.0)
    for role in (ROLE_APRON, ROLE_TERMINAL, ROLE_RUNWAY):
        result = _drop_primary_parallels_embedded_in_pavement(
            [(poly, axis, role, "X")], big_apron)
        assert len(result) == 1, f"role {role} must pass through"


def test_short_rect_under_min_axis_passes_through():
    """Rects shorter than 30 m are skipped (too short to meaningfully
    measure adjacency)."""
    big_apron = _apron(-50.0, 100.0, -50.0, 50.0)
    poly, axis = _rect(length=25.0)  # < 30 m
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], big_apron)
    assert len(result) == 1


# ──────────────────────────────────────────────────────────────────────
# THE EITHER-SIDE RULE (the recurring-regression hot spot)
# ──────────────────────────────────────────────────────────────────────
def test_either_side_rule_north_only_absorbs():
    """Apron extending only past the NORTH long edge (y > +half_w + 5
    for the full rect length) absorbs.  This is half the SPJC-F case:
    one side embedded, the other on apron's outer boundary."""
    # Rect: half_w=10 (top edge at y=+10).  Apron starts at y=+5
    # (below top edge — covers part of rect) and extends to y=+50.
    # The +y probe at y = +half_w + 5 = +15 lands inside apron.
    apron = _apron(-50.0, 150.0, +5.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    # Either fully absorbed (zero rects) OR no-op.  EITHER-side rule
    # says north adjacency alone is enough to absorb.
    assert len(result) == 0, (
        "EITHER-side rule violated: north-only adjacency should "
        "absorb")


def test_either_side_rule_south_only_absorbs():
    """Apron extending only past the SOUTH long edge absorbs too.
    Symmetric with the north case."""
    apron = _apron(-50.0, 150.0, -50.0, -5.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    assert len(result) == 0, (
        "EITHER-side rule violated: south-only adjacency should "
        "absorb")


def test_neither_side_keeps_rect():
    """When NEITHER long edge has adjacent apron, the rect must
    survive unchanged (the "primary parallel through grass" case)."""
    # Apron is far away from the rect — outside the +5 m probe range
    # on both sides.
    apron = _apron(200.0, 300.0, -50.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    assert len(result) == 1


# ──────────────────────────────────────────────────────────────────────
# PARTIAL ABSORPTION (the never-absorb-whole-rect-when-only-part-covered rule)
# ──────────────────────────────────────────────────────────────────────
def test_partial_adjacency_splits_not_drops():
    """Apron covering ONLY the first 50 m of the long edge should
    absorb that portion and KEEP the remaining 50 m as a shorter
    rect — not drop the whole rect."""
    # Apron on north side, x=0..50.  Rect length 100.
    apron = _apron(-10.0, 50.0, +5.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    assert len(result) == 1, (
        "Partial-absorption rule violated: half-covered rect must "
        "produce one shorter rect, not drop the whole thing")
    # Surviving rect axis should cover roughly x=50..100.
    spans = _kept_axes(result)
    span_lo, span_hi = spans[0]
    assert span_hi - span_lo >= 30.0  # MIN_KEPT_M
    assert span_hi <= 100.0


def test_two_separate_adjacencies_split_into_three():
    """Apron at x=0..30 and x=70..100 (with a clear gap in the
    middle) should leave the middle 30..70 chunk as a rect."""
    # Two apron strips on the north side.
    a1 = _apron(-10.0, 30.0, +5.0, +50.0)
    a2 = _apron(70.0, 110.0, +5.0, +50.0)
    apron = a1.union(a2)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    # The middle (~30..70 = 40 m) survives as one rect.
    assert len(result) == 1
    spans = _kept_axes(result)
    span_lo, span_hi = spans[0]
    assert 25.0 <= span_lo <= 35.0
    assert 65.0 <= span_hi <= 75.0


def test_full_absorption_drops_all():
    """Apron covering the rect's full length on at least one side
    absorbs the whole rect (no kept fragment ≥ 30 m)."""
    apron = _apron(-50.0, 150.0, +5.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    assert len(result) == 0


def test_short_kept_fragment_under_30m_dropped():
    """Kept fragments < 30 m get dropped (avoid emitting micro-stubs).
    Apron from x=0..80 leaves only 20 m (< 30) at the end → drop."""
    apron = _apron(-10.0, 80.0, +5.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], apron)
    # 20 m kept fragment is below MIN_KEPT_M → dropped → 0 rects.
    assert len(result) == 0


# ──────────────────────────────────────────────────────────────────────
# 10 % THRESHOLD (don't absorb on noise)
# ──────────────────────────────────────────────────────────────────────
def test_below_threshold_adjacency_keeps_rect():
    """Adjacent run < 10 % of axial steps (e.g. one 5 m sliver out of
    100 m) is below the noise threshold and should NOT trigger
    absorption."""
    # 4 m of apron near x=0 — well under 10 % of 100 m axis (10 m).
    sliver = _apron(0.0, 4.0, +5.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], sliver)
    assert len(result) == 1
    spans = _kept_axes(result)
    assert spans[0] == (0.0, 100.0), (
        "below-threshold sliver must not trigger any absorption")


def test_at_threshold_adjacency_triggers():
    """An adjacent run that meets the 10 % threshold triggers
    absorption.  Verified by: original rect does NOT survive
    unchanged.  (The kept fragment may further drop via the apron-
    interior check; that's a separate guarantee — what matters here
    is that 10 % adjacency took action.)"""
    # 12 m apron strip near x=0 covers ~3 sample steps (≥ 10 %
    # of the 21 steps along a 100 m axis).
    strip = _apron(-5.0, 12.0, +5.0, +50.0)
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")], strip)
    # If a fragment survives, it's strictly shorter than the input.
    if result:
        spans = _kept_axes(result)
        for lo, hi in spans:
            assert (hi - lo) < 100.0, (
                "Threshold-meeting adjacency must trim the rect, "
                "not return it unchanged")
    # The all-dropped case is also valid (apron-interior check).


# ──────────────────────────────────────────────────────────────────────
# RUNWAY EXCLUSION
# ──────────────────────────────────────────────────────────────────────
def test_runway_alongside_does_not_absorb():
    """A primary parallel running alongside a RUNWAY should NOT be
    absorbed — runway-adjacency is excluded from junction_pav.

    SPJC's primary parallels along the 16/34 runway are the
    canonical case.
    """
    runway = _apron(-50.0, 150.0, +5.0, +30.0)
    # Same geometry as the apron BUT classified as runway via
    # runway_polys.
    poly, axis = _rect(length=100.0, half_w=10.0)
    result = _drop_primary_parallels_embedded_in_pavement(
        [(poly, axis, ROLE_PRIMARY_PARALLEL, "F")],
        apt_pav_union=runway,
        runway_polys=[runway])
    assert len(result) == 1, (
        "Runway-exclusion rule violated: primary parallel along a "
        "runway must not absorb")


# ──────────────────────────────────────────────────────────────────────
# All sloping roles obey the rule
# ──────────────────────────────────────────────────────────────────────
def test_all_sloping_roles_subject_to_rule():
    """Every sloping rect role (primary/secondary parallel, stub,
    cross_connector) absorbs when its long edge has apron adjacency
    ≥ 10 %."""
    apron = _apron(-50.0, 150.0, +5.0, +50.0)
    for role in (ROLE_PRIMARY_PARALLEL, ROLE_SECONDARY_PARALLEL,
                 ROLE_STUB, ROLE_CROSS_CONNECTOR):
        poly, axis = _rect(length=100.0, half_w=10.0)
        result = _drop_primary_parallels_embedded_in_pavement(
            [(poly, axis, role, "X")], apron)
        assert len(result) == 0, (
            f"Sloping role {role!r} must obey the absorption rule")


# ──────────────────────────────────────────────────────────────────────
# Public-name aliases match the underscored ones
# ──────────────────────────────────────────────────────────────────────
def test_public_alias_is_same_function():
    """The non-underscored ``drop_primary_parallels_embedded_in_pavement``
    is the same callable as the underscored private form."""
    assert (drop_primary_parallels_embedded_in_pavement
            is _drop_primary_parallels_embedded_in_pavement)
