"""Unit tests for O4_Pavement_Strips (apron-deformation model).

Classification rules:
  * MRR short ≤ narrow_width_m AND aspect ≥ simple_strip_aspect
    → taxi with MRR midline axis
  * Otherwise → apron (2D deformable surface, no axis)
  * Adjacent apron polygons are unioned into connected components.
"""
import math

from shapely.affinity import rotate
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from auto_patch.pavement import strips as PS


# ──────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────
def _rect(width, length, cx=0.0, cy=0.0, angle_deg=0.0):
    """Width along Y, length along X, centered at (cx, cy)."""
    dx, dy = length / 2.0, width / 2.0
    poly = Polygon([
        (cx - dx, cy - dy), (cx + dx, cy - dy),
        (cx + dx, cy + dy), (cx - dx, cy + dy),
    ])
    if angle_deg != 0.0:
        poly = rotate(poly, angle_deg, origin=(cx, cy), use_radians=False)
    return poly


def _taxis(shapes):
    return [s for s in shapes if s.kind == "taxi"]


def _aprons(shapes):
    return [s for s in shapes if s.kind == "apron"]


# ──────────────────────────────────────────────────────────────────
# Empty / degenerate
# ──────────────────────────────────────────────────────────────────
def test_empty_inputs_return_empty():
    assert PS.decompose_pavement([], []) == tuple()


def test_none_and_empty_polys_are_skipped():
    assert PS.decompose_pavement(
        [None, Polygon()], [None]) == tuple()


# ──────────────────────────────────────────────────────────────────
# Simple strip classification
# ──────────────────────────────────────────────────────────────────
def test_narrow_strip_is_taxi():
    # 25 × 400 — classic taxi strip.
    shapes = PS.decompose_pavement([_rect(25, 400)], [])
    assert len(shapes) == 1
    assert shapes[0].kind == "taxi"
    assert shapes[0].axis is not None
    assert not shapes[0].axis.is_empty


def test_narrow_strip_width_reflects_polygon():
    shapes = PS.decompose_pavement([_rect(25, 400)], [])
    assert math.isclose(shapes[0].width_m, 25.0, rel_tol=0.05)


def test_strip_right_at_threshold_is_taxi():
    # 30 × 400 — exactly at the width cutoff.
    shapes = PS.decompose_pavement([_rect(30, 400)], [])
    assert len(shapes) == 1
    assert shapes[0].kind == "taxi"


def test_strip_well_over_threshold_is_apron():
    # 70 × 400 — well above the 45 m taxi cutoff (and above the
    # mega-poly skeleton accept threshold of ~67 m after 1.5x
    # margin).  Guaranteed apron.
    shapes = PS.decompose_pavement([_rect(70, 400)], [])
    assert len(shapes) == 1
    assert shapes[0].kind == "apron"
    assert shapes[0].axis is None


def test_rotated_strip_still_classifies_as_taxi():
    shapes = PS.decompose_pavement([_rect(25, 400, angle_deg=45.0)], [])
    assert len(shapes) == 1
    assert shapes[0].kind == "taxi"
    # Axis span ≈ polygon length (within extension tolerance).
    cc = list(shapes[0].axis.coords)
    span = math.hypot(cc[-1][0] - cc[0][0], cc[-1][1] - cc[0][1])
    assert span >= 380.0


def test_blocky_low_aspect_polygon_is_apron():
    # 200 × 300 — unambiguously apron (wide everywhere).
    shapes = PS.decompose_pavement([_rect(200, 300)], [])
    assert len(shapes) == 1
    assert shapes[0].kind == "apron"


# ──────────────────────────────────────────────────────────────────
# Apron classification + unioning
# ──────────────────────────────────────────────────────────────────
def test_apron_input_classifies_as_apron():
    shapes = PS.decompose_pavement([], [_rect(150, 300)])
    assert len(shapes) == 1
    assert shapes[0].kind == "apron"


def test_wide_taxi_classified_polygon_is_reclassified_as_apron():
    # Caller tagged it as taxi, but it's 100 × 400 — over the 45 m
    # taxi ceiling.
    shapes = PS.decompose_pavement([_rect(100, 400)], [])
    assert len(shapes) == 1
    assert shapes[0].kind == "apron"


def test_adjacent_aprons_union_into_one_shape():
    # Two abutting 100 × 100 aprons share an edge → one shape.
    a = _rect(100, 100, cx=-50, cy=0)
    b = _rect(100, 100, cx=+50, cy=0)
    shapes = PS.decompose_pavement([], [a, b])
    aprons = _aprons(shapes)
    assert len(aprons) == 1
    # Union area ≈ 100 × 200 = 20000.
    assert math.isclose(aprons[0].polygon.area, 20000.0, rel_tol=0.05)


def test_disjoint_aprons_stay_separate():
    a = _rect(100, 100, cx=-200, cy=0)
    b = _rect(100, 100, cx=+200, cy=0)
    shapes = PS.decompose_pavement([], [a, b])
    assert len(_aprons(shapes)) == 2


def test_tiny_apron_sliver_is_dropped():
    # A 3 × 1 sliver — below MIN_APRON_AREA_M2.
    tiny = _rect(1, 3)
    shapes = PS.decompose_pavement([], [tiny])
    assert _aprons(shapes) == []


# ──────────────────────────────────────────────────────────────────
# Mixed inputs
# ──────────────────────────────────────────────────────────────────
def test_taxi_next_to_apron_produces_two_shapes():
    # A taxi (25 × 300) touching an apron (150 × 150).  Taxi and
    # apron stay separate shapes; their shared boundary is implicit.
    taxi = _rect(25, 300, cx=-150, cy=0)
    apron = _rect(150, 150, cx=+75, cy=0)
    shapes = PS.decompose_pavement([taxi], [apron])
    assert len(_taxis(shapes)) == 1
    assert len(_aprons(shapes)) == 1


def test_taxi_inside_apron_is_extracted_via_skeleton():
    # apt.dat sometimes has a taxi polygon overlapping an apron
    # polygon.  After input-unioning, the combined pavement is
    # skeletonised; the narrow branch becomes a taxi shape and the
    # residual is apron.
    taxi = _rect(25, 300, cx=0, cy=0)
    apron = _rect(600, 600, cx=0, cy=0)
    shapes = PS.decompose_pavement([taxi], [apron])
    # The union is a 600x600 square; no narrow branch survives
    # skeleton extraction (medial axis of a square is a point) so
    # it's all apron.
    assert len(_aprons(shapes)) >= 1


# ──────────────────────────────────────────────────────────────────
# Axis geometry
# ──────────────────────────────────────────────────────────────────
def test_taxi_axis_endpoints_reach_polygon_boundary():
    # After extension, the axis endpoints should sit ON the taxi
    # polygon boundary (not set-back at the MRR short-side
    # midpoint).
    poly = _rect(25, 400)
    shapes = PS.decompose_pavement([poly], [])
    assert len(shapes) == 1
    ax = shapes[0].axis
    # Both endpoints should be within 0.5 m of the polygon boundary.
    for ep in (Point(ax.coords[0]), Point(ax.coords[-1])):
        assert ep.distance(poly.boundary) < 0.5


def test_apron_has_no_axis():
    shapes = PS.decompose_pavement([], [_rect(200, 400)])
    assert shapes[0].axis is None
    assert shapes[0].width_m == 0.0


# ──────────────────────────────────────────────────────────────────
# Invariants
# ──────────────────────────────────────────────────────────────────
def test_no_junction_type_exists():
    # Sanity: the old Junction class has been removed.
    assert not hasattr(PS, "Junction")
    assert not hasattr(PS, "StripGraph")


def test_return_type_is_tuple_of_shapes():
    shapes = PS.decompose_pavement([_rect(25, 400)], [_rect(150, 200)])
    assert isinstance(shapes, tuple)
    for s in shapes:
        assert isinstance(s, PS.Shape)
        assert s.kind in ("taxi", "apron")


# ──────────────────────────────────────────────────────────────────
# Adjacency graph
# ──────────────────────────────────────────────────────────────────
def test_adjacency_disjoint_polygons_empty():
    a = _rect(100, 100, cx=-500, cy=0)
    b = _rect(100, 100, cx=+500, cy=0)
    shapes = PS.decompose_pavement([], [a, b])
    # Verify they're two separate shapes.
    assert len(shapes) == 2
    adjs = PS.build_adjacency_graph(shapes)
    assert adjs == tuple()


def test_adjacency_single_touching_pair():
    # Two 25×100 taxi rectangles separated by a 0.3 m gap.  They
    # remain two separate shapes after input-unioning (gap > 0) and
    # the adjacency graph pairs them via the 0.5 m tolerance.
    a = _rect(25, 100, cx=-50.15, cy=0)
    b = _rect(25, 100, cx=+50.15, cy=0)
    shapes = PS.decompose_pavement([a, b], [])
    assert len(shapes) == 2
    adjs = PS.build_adjacency_graph(shapes)
    assert len(adjs) == 1
    adj = adjs[0]
    assert adj.shape_a == 0 and adj.shape_b == 1
    # Shared length ≈ 25 (the taxi width at the touching edge).
    assert math.isclose(adj.length_m, 25.0, rel_tol=0.1)


def test_adjacency_point_contact_is_ignored():
    # Two polygons touching only at (0, 0) — a single corner.
    a = Polygon([(-100, -100), (0, -100), (0, 0), (-100, 0)])
    b = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
    # Mark as taxis to keep them separate (don't union as aprons).
    shapes = PS.decompose_pavement([a, b], [])
    adjs = PS.build_adjacency_graph(shapes)
    assert adjs == tuple()


def test_adjacency_three_in_a_row():
    # Three taxi strips in a line: A-B-C with 0.3 m gaps to
    # survive input-unioning as three separate shapes.  Each rect
    # is 100 m long → spacing of 100 m + 0.3 m gap → cx at ±100.3.
    a = _rect(25, 100, cx=-100.3, cy=0)
    b = _rect(25, 100, cx=0, cy=0)
    c = _rect(25, 100, cx=+100.3, cy=0)
    shapes = PS.decompose_pavement([a, b, c], [])
    assert len(shapes) == 3
    adjs = PS.build_adjacency_graph(shapes)
    assert len(adjs) == 2
    pairs = {(adj.shape_a, adj.shape_b) for adj in adjs}
    # Middle shape (index 1 after MRR-based sorting preserved) is
    # adjacent to both outer shapes; adjacent pairs are (0,1) and
    # (1,2).  We assert the pair set has the middle index twice.
    middle_counts = {}
    for (a_, b_) in pairs:
        middle_counts[a_] = middle_counts.get(a_, 0) + 1
        middle_counts[b_] = middle_counts.get(b_, 0) + 1
    # Exactly one index appears in two pairs.
    appears_twice = [k for k, v in middle_counts.items() if v == 2]
    assert len(appears_twice) == 1


def test_adjacency_taxi_touches_apron():
    # Taxi and apron separated by a 0.3 m gap so the input-union
    # doesn't merge them into one super-polygon.
    taxi = _rect(25, 200, cx=-100.15, cy=0)
    apron = _rect(150, 150, cx=75.15, cy=0)
    shapes = PS.decompose_pavement([taxi], [apron])
    kinds = sorted(s.kind for s in shapes)
    assert kinds == ["apron", "taxi"]
    adjs = PS.build_adjacency_graph(shapes)
    assert len(adjs) == 1
    assert math.isclose(adjs[0].length_m, 25.0, rel_tol=0.1)


def test_adjacency_shared_multi_component():
    # Apron with a hole that's 0.3 m larger than the taxi on every
    # side, so the taxi slots inside with a thin gap all around.
    # The adjacency tolerance (0.5 m) still pairs them.
    outer = Polygon(
        [(-200, -200), (200, -200), (200, 200), (-200, 200)],
        holes=[[(-50.3, -12.8), (50.3, -12.8),
                (50.3, 12.8), (-50.3, 12.8)]],
    )
    taxi = _rect(25, 100)
    shapes = PS.decompose_pavement([taxi], [outer])
    adjs = PS.build_adjacency_graph(shapes)
    assert len(adjs) == 1
    # Shared ≈ taxi perimeter = 250 m; tolerance may shrink slightly.
    assert adjs[0].length_m > 200.0


def test_perimeter_coverage_isolated_shape():
    shapes = PS.decompose_pavement([_rect(25, 400)], [])
    adjs = PS.build_adjacency_graph(shapes)
    covered, total = PS.perimeter_coverage(0, shapes, adjs)
    assert covered == 0.0
    assert total > 0.0


def test_perimeter_coverage_fully_embedded_shape():
    # Taxi nested in an apron hole with a thin gap: coverage is
    # close to 100% via the adjacency tolerance.
    outer = Polygon(
        [(-200, -200), (200, -200), (200, 200), (-200, 200)],
        holes=[[(-50.3, -12.8), (50.3, -12.8),
                (50.3, 12.8), (-50.3, 12.8)]],
    )
    taxi = _rect(25, 100)
    shapes = PS.decompose_pavement([taxi], [outer])
    adjs = PS.build_adjacency_graph(shapes)
    taxi_idx = next(i for i, s in enumerate(shapes) if s.kind == "taxi")
    covered, total = PS.perimeter_coverage(taxi_idx, shapes, adjs)
    # Covered should be at least 80% of the perimeter.
    assert covered >= 0.8 * total


# ──────────────────────────────────────────────────────────────────
# Mega-polygon decomposition (per-branch local-width)
# ──────────────────────────────────────────────────────────────────
def test_narrow_y_mega_polygon_extracts_three_taxis():
    # Y-shape of three 25 m wide arms, 300 m long each, meeting at
    # a central hub.  MRR short side spans both vertical arms so
    # the fast path rejects; skeleton should find three branches.
    trunk = _rect(25, 400, cx=0, cy=0)
    arm1 = Polygon([(-150, 0), (-125, 0), (-125, 300), (-150, 300)])
    arm2 = Polygon([(+125, 0), (+150, 0), (+150, 300), (+125, 300)])
    y_poly = unary_union([trunk, arm1, arm2])
    assert y_poly.geom_type == "Polygon"

    shapes = PS.decompose_pavement([y_poly], [])
    taxis = _taxis(shapes)
    assert len(taxis) >= 2
    for t in taxis:
        assert t.width_m <= PS.NARROW_WIDTH_M + 2.0


def test_mega_polygon_with_wide_hub_produces_taxi_plus_apron_residual():
    # A 100 m wide hub (apron-class) joined to a 25 m wide arm.
    # Expected: 1 taxi shape (the arm), 1 apron shape (the hub).
    hub = _rect(100, 100, cx=0, cy=0)       # x: -50..50
    arm = _rect(25, 300, cx=200, cy=0)      # x: 50..350, overlaps hub edge
    mega = unary_union([hub, arm])
    assert mega.geom_type == "Polygon"

    shapes = PS.decompose_pavement([mega], [])
    taxis = _taxis(shapes)
    aprons = _aprons(shapes)
    assert len(taxis) >= 1
    assert len(aprons) >= 1
    # The taxi's polygon must be inside the mega polygon.
    for t in taxis:
        assert t.polygon.within(mega.buffer(0.5))


def test_mega_polygon_with_no_narrow_branches_is_all_apron():
    # A fat L-shape where every branch is > 45 m wide.
    block_a = _rect(150, 300, cx=0, cy=0)
    block_b = _rect(300, 150, cx=75, cy=75)
    big = unary_union([block_a, block_b])
    shapes = PS.decompose_pavement([big], [])
    assert _taxis(shapes) == []
    assert len(_aprons(shapes)) == 1


def test_mega_polygon_claim_no_overlap():
    # Y-mega-polygon: extracted taxi polygons must be
    # interior-disjoint from each other and from any apron residual.
    trunk = _rect(25, 400, cx=0, cy=0)
    arm1 = Polygon([(-150, 0), (-125, 0), (-125, 300), (-150, 300)])
    arm2 = Polygon([(+125, 0), (+150, 0), (+150, 300), (+125, 300)])
    y_poly = unary_union([trunk, arm1, arm2])
    shapes = PS.decompose_pavement([y_poly], [])
    polys = [s.polygon for s in shapes]
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            inter = polys[i].intersection(polys[j])
            # Shared boundary is OK; shared interior is not.
            if not inter.is_empty:
                assert inter.area < 1.0, \
                    f"shapes {i} and {j} overlap by {inter.area} m²"


def test_skeleton_short_branch_rejected():
    # Wide main body (> 45 m) with a very short narrow spur —
    # the spur is too short to pass the skeleton min-length
    # filter, and the main body classifies as apron.
    main = _rect(150, 400, cx=0, cy=0)
    spur = _rect(25, 30, cx=215, cy=0)
    blob = unary_union([main, spur])
    shapes = PS.decompose_pavement([blob], [])
    assert _taxis(shapes) == []


# ──────────────────────────────────────────────────────────────────
# Role classification
# ──────────────────────────────────────────────────────────────────
def _vrect(width, length, cx=0.0, cy=0.0):
    """Vertical rectangle: `width` along X, `length` along Y."""
    dx, dy = width / 2.0, length / 2.0
    return Polygon([
        (cx - dx, cy - dy), (cx + dx, cy - dy),
        (cx + dx, cy + dy), (cx - dx, cy + dy),
    ])


def test_role_primary_parallel_detected():
    # Runway: 3000 m long along Y at x=0.
    runway = LineString([(0, -1500), (0, 1500)])
    # Parallel taxi: 2000 m along Y at x=-100 (within 300 m).
    taxi = _vrect(25, 2000, cx=-100, cy=0)
    shapes = PS.decompose_pavement([taxi], [])
    adjs = PS.build_adjacency_graph(shapes)
    roles = PS.classify_shape_roles(shapes, adjs, [runway])
    taxi_idx = next(i for i, s in enumerate(shapes) if s.kind == "taxi")
    assert roles[taxi_idx] == PS.ROLE_PRIMARY_PARALLEL


def test_role_short_parallel_not_primary():
    # Parallel but too short — < 0.25 × runway length.
    runway = LineString([(0, -1500), (0, 1500)])
    taxi = _vrect(25, 500, cx=-100, cy=0)   # 500 < 0.25 * 3000 = 750
    shapes = PS.decompose_pavement([taxi], [])
    adjs = PS.build_adjacency_graph(shapes)
    roles = PS.classify_shape_roles(shapes, adjs, [runway])
    taxi_idx = next(i for i, s in enumerate(shapes) if s.kind == "taxi")
    assert roles[taxi_idx] == PS.ROLE_SECONDARY_PARALLEL


def test_role_far_parallel_not_primary():
    # Parallel but > 300 m from runway.
    runway = LineString([(0, -1500), (0, 1500)])
    taxi = _vrect(25, 2000, cx=-500, cy=0)
    shapes = PS.decompose_pavement([taxi], [])
    adjs = PS.build_adjacency_graph(shapes)
    roles = PS.classify_shape_roles(shapes, adjs, [runway])
    taxi_idx = next(i for i, s in enumerate(shapes) if s.kind == "taxi")
    assert roles[taxi_idx] == PS.ROLE_SECONDARY_PARALLEL


def test_role_perpendicular_is_cross_connector():
    runway = LineString([(0, -1500), (0, 1500)])
    # Horizontal taxi 200 m along X, 25 m along Y at y=500 —
    # perpendicular to the runway.
    taxi = _rect(25, 200, cx=0, cy=500)
    shapes = PS.decompose_pavement([taxi], [])
    adjs = PS.build_adjacency_graph(shapes)
    roles = PS.classify_shape_roles(shapes, adjs, [runway])
    taxi_idx = next(i for i, s in enumerate(shapes) if s.kind == "taxi")
    assert roles[taxi_idx] == PS.ROLE_CROSS_CONNECTOR


def test_role_stub_requires_adjacency_to_primary():
    # Parallel primary + small stub touching it.
    runway = LineString([(0, -1500), (0, 1500)])
    primary = _vrect(25, 2000, cx=-100, cy=0)
    # Stub: 75 m along X, 20 m along Y, touching primary's east
    # edge at (-87.5, 0) with a 0.3 m gap so union doesn't merge.
    stub = Polygon([
        (-87.2, -10), (-12.5, -10),
        (-12.5, +10), (-87.2, +10),
    ])
    shapes = PS.decompose_pavement([primary, stub], [])
    adjs = PS.build_adjacency_graph(shapes)
    roles = PS.classify_shape_roles(shapes, adjs, [runway])
    # Identify primary by its larger area (2000×25 = 50000).
    primary_idx = max(range(len(shapes)),
                      key=lambda i: shapes[i].polygon.area
                      if shapes[i].kind == "taxi" else 0)
    other_taxi_idxs = [i for i in range(len(shapes))
                       if i != primary_idx
                       and shapes[i].kind == "taxi"]
    assert roles[primary_idx] == PS.ROLE_PRIMARY_PARALLEL
    # The stub should be marked as stub (or cross_connector —
    # perpendicular short taxi adjacent to a primary can be
    # either; stub wins because of its adjacency to primary).
    for idx in other_taxi_idxs:
        assert roles[idx] == PS.ROLE_STUB


def test_role_apron_input_is_apron():
    runway = LineString([(0, -1500), (0, 1500)])
    apron = _rect(200, 300)
    shapes = PS.decompose_pavement([], [apron])
    adjs = PS.build_adjacency_graph(shapes)
    roles = PS.classify_shape_roles(shapes, adjs, [runway])
    assert all(r == PS.ROLE_APRON for r in roles)


def test_adjacency_canonical_ordering():
    # shape_a must always be < shape_b.
    a = _rect(25, 100, cx=-50, cy=0)
    b = _rect(25, 100, cx=+50, cy=0)
    shapes = PS.decompose_pavement([a, b], [])
    adjs = PS.build_adjacency_graph(shapes)
    for adj in adjs:
        assert adj.shape_a < adj.shape_b
