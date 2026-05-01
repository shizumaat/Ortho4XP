"""Unit tests for auto_patch.layout — data model + OSM emission.

Covers:
* PavementLayout.m_to_ll / ll_to_m round-trip.
* _projection / _airport_anchor helpers.
* PavementLayout.to_osm: shape emission, shared-vertex IDs, tag
  handling per role and per elevation form.
"""
import math
import re
import tempfile
from pathlib import Path

import pytest
from shapely.geometry import Polygon

from auto_patch.layout import (
    AEROWAY_FOR_ROLE,
    BuiltShape,
    PavementLayout,
    R_EARTH,
    ROLE_APRON,
    ROLE_JUNCTION,
    ROLE_PRIMARY_PARALLEL,
    ROLE_RUNWAY,
    ROLE_TERMINAL,
    SHARED_VERTEX_TOL_M,
    _airport_anchor,
    _projection,
)


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def _make_layout(icao="KFAKE", lat=40.0, lon=-100.0):
    """Empty PavementLayout anchored at the given lat/lon."""
    return PavementLayout(icao=icao, anchor=(lat, lon))


def _square(cx=0.0, cy=0.0, side=10.0):
    """Axis-aligned square in meter space."""
    h = side / 2.0
    return Polygon([
        (cx - h, cy - h),
        (cx + h, cy - h),
        (cx + h, cy + h),
        (cx - h, cy + h),
    ])


def _emit_and_parse(layout):
    """Call to_osm then parse the resulting XML.

    Returns (nodes, ways) where:
      nodes: {nid: (lat, lon)}
      ways:  [(wid, [nid_refs], {tag: val})]
    """
    with tempfile.NamedTemporaryFile(
            mode="r", suffix=".osm", delete=False) as f:
        path = f.name
    try:
        layout.to_osm(path)
        text = Path(path).read_text()
    finally:
        Path(path).unlink()

    node_re = re.compile(
        r"""<node id='(-?\d+)'[^>]*lat='([^']+)' lon='([^']+)'\s*/>""")
    way_open_re = re.compile(r"""<way id='(-?\d+)'""")
    nd_re = re.compile(r"""<nd ref='(-?\d+)'""")
    tag_re = re.compile(r"""<tag k='([^']+)' v='([^']+)'""")

    nodes = {}
    for m in node_re.finditer(text):
        nodes[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))

    ways = []
    way_blocks = re.findall(
        r"<way id='-?\d+'[^>]*>(.*?)</way>", text, flags=re.DOTALL)
    for i, m_open in enumerate(way_open_re.finditer(text)):
        wid = int(m_open.group(1))
        body = way_blocks[i]
        nds = [int(x) for x in nd_re.findall(body)]
        tags = {k: v for k, v in tag_re.findall(body)}
        ways.append((wid, nds, tags))

    return nodes, ways


# ──────────────────────────────────────────────────────────────────────
# Coordinate round-trip
# ──────────────────────────────────────────────────────────────────────
def test_m_to_ll_then_ll_to_m_is_identity():
    """The two coordinate transforms invert each other to sub-mm
    precision around the anchor."""
    layout = _make_layout(lat=37.5, lon=-122.3)
    for x, y in [(0.0, 0.0), (100.0, 200.0), (-50.0, 1500.0),
                 (12345.6, -7890.1)]:
        lat, lon = layout.m_to_ll(x, y)
        x2, y2 = layout.ll_to_m(lat, lon)
        assert abs(x - x2) < 1e-6
        assert abs(y - y2) < 1e-6


def test_m_to_ll_origin_is_anchor():
    """(0, 0) in meter space is the anchor lat/lon exactly."""
    layout = _make_layout(lat=40.0, lon=-100.0)
    lat, lon = layout.m_to_ll(0.0, 0.0)
    assert abs(lat - 40.0) < 1e-9
    assert abs(lon - (-100.0)) < 1e-9


def test_ll_to_m_anchor_is_origin():
    """The anchor lat/lon maps to (0, 0)."""
    layout = _make_layout(lat=40.0, lon=-100.0)
    x, y = layout.ll_to_m(40.0, -100.0)
    assert abs(x) < 1e-9
    assert abs(y) < 1e-9


# ──────────────────────────────────────────────────────────────────────
# _projection helper
# ──────────────────────────────────────────────────────────────────────
def test_projection_returns_callable_to_m():
    """_projection(anchor) returns a (lon, lat[, z]) → (x, y[, z])
    function whose origin is the anchor."""
    to_m = _projection((40.0, -100.0))
    assert callable(to_m)
    x, y = to_m(-100.0, 40.0)
    assert abs(x) < 1e-9
    assert abs(y) < 1e-9


def test_projection_with_z_passes_through():
    """The optional z argument is preserved in the output tuple."""
    to_m = _projection((40.0, -100.0))
    x, y, z = to_m(-100.0, 40.0, 123.4)
    assert z == 123.4


# ──────────────────────────────────────────────────────────────────────
# _airport_anchor
# ──────────────────────────────────────────────────────────────────────
class _StubRunway:
    def __init__(self, lat_a, lon_a, lat_b, lon_b):
        self.lat_a, self.lon_a = lat_a, lon_a
        self.lat_b, self.lon_b = lat_b, lon_b


class _StubAirport:
    def __init__(self, runways=(), boundary=None):
        self.runways = list(runways)
        self.boundary = boundary


def test_airport_anchor_uses_first_runway_midpoint():
    """When runways are present, anchor = midpoint of first runway."""
    apt = _StubAirport(runways=[
        _StubRunway(40.0, -100.0, 40.02, -100.0),
    ])
    lat, lon = _airport_anchor(apt)
    assert abs(lat - 40.01) < 1e-9
    assert abs(lon - (-100.0)) < 1e-9


def test_airport_anchor_falls_back_to_boundary_centroid():
    """No runways → use boundary centroid (note shapely centroid is
    .y=lat, .x=lon)."""
    boundary = Polygon([
        (-100.0, 40.0),
        (-99.99, 40.0),
        (-99.99, 40.02),
        (-100.0, 40.02),
    ])
    apt = _StubAirport(boundary=boundary)
    lat, lon = _airport_anchor(apt)
    assert abs(lat - 40.01) < 1e-9
    assert abs(lon - (-99.995)) < 1e-9


def test_airport_anchor_default_zero():
    """No runways AND no boundary → (0, 0)."""
    apt = _StubAirport()
    assert _airport_anchor(apt) == (0.0, 0.0)


# ──────────────────────────────────────────────────────────────────────
# to_osm: smoke-level round-trip
# ──────────────────────────────────────────────────────────────────────
def test_to_osm_emits_one_way_per_shape():
    """Each non-empty BuiltShape produces one <way> in the output."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY, ref="RW09/RW27"))
    layout.shapes.append(BuiltShape(
        polygon=_square(50, 0, 10), role=ROLE_PRIMARY_PARALLEL, ref="A"))
    nodes, ways = _emit_and_parse(layout)
    assert len(ways) == 2


def test_to_osm_skips_empty_polygons():
    """Empty / None polygons are skipped silently."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(polygon=Polygon(), role=ROLE_RUNWAY))
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY, ref="RW09/RW27"))
    _, ways = _emit_and_parse(layout)
    assert len(ways) == 1


def test_to_osm_emits_aeroway_and_role_tags():
    """Every way has aeroway= and role= tags set per the AEROWAY_FOR_ROLE
    map."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY))
    _, ways = _emit_and_parse(layout)
    _, _, tags = ways[0]
    assert tags["aeroway"] == AEROWAY_FOR_ROLE[ROLE_RUNWAY]
    assert tags["role"] == ROLE_RUNWAY


def test_to_osm_ref_tag_emitted_when_set():
    """``ref`` gets passed through to the OSM tag."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY, ref="RW09/RW27"))
    _, ways = _emit_and_parse(layout)
    assert ways[0][2]["ref"] == "RW09/RW27"


def test_to_osm_no_ref_tag_when_unset():
    """An empty ref is not emitted."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_JUNCTION))
    _, ways = _emit_and_parse(layout)
    assert "ref" not in ways[0][2]


# ──────────────────────────────────────────────────────────────────────
# to_osm: elevation tag formats
# ──────────────────────────────────────────────────────────────────────
def test_to_osm_sloped_rect_emits_high_low_cell_profile():
    """A rect with altitude_high + altitude_low gets the full
    legacy patch-format quad: altitude_high, altitude_low,
    cell_size, profile."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_PRIMARY_PARALLEL,
        altitude_high=100.5, altitude_low=99.0))
    _, ways = _emit_and_parse(layout)
    tags = ways[0][2]
    assert tags["altitude_high"] == "100.5"
    assert tags["altitude_low"] == "99.0"
    assert tags["cell_size"] == "2"
    assert tags["profile"] == "spline"
    # Flat altitude tag must NOT also be emitted.
    assert "altitude" not in tags


def test_to_osm_flat_altitude_emits_single_tag():
    """A polygon with only ``altitude`` set emits a single
    altitude= tag."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_TERMINAL,
        altitude=42.7))
    _, ways = _emit_and_parse(layout)
    tags = ways[0][2]
    assert tags["altitude"] == "42.7"
    assert "altitude_high" not in tags
    assert "altitude_low" not in tags


def test_to_osm_node_altitudes_emits_csv_list():
    """A polygon with ``node_altitudes`` emits a comma-separated
    list with one value per ring vertex INCLUDING the closing
    repeat."""
    layout = _make_layout()
    poly = _square(0, 0, 10)
    # 4 corners + closing repeat = 5 elevations.
    elevs = [10.0, 11.0, 12.0, 13.0, 10.0]
    layout.shapes.append(BuiltShape(
        polygon=poly, role=ROLE_JUNCTION,
        node_altitudes=elevs))
    _, ways = _emit_and_parse(layout)
    tags = ways[0][2]
    assert "node_altitudes" in tags
    csv_vals = tags["node_altitudes"].split(",")
    # Length must equal the emitted ring length (which includes
    # closing repeat).  Square → 4 unique nids + 1 closing = 5.
    assert len(csv_vals) == len(ways[0][1])
    assert csv_vals[0] == "10.0"
    assert csv_vals[1] == "11.0"


def test_to_osm_no_elevation_tags_when_none_set():
    """A geometry-only polygon (all elevation fields None) emits no
    elevation tags."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_APRON))
    _, ways = _emit_and_parse(layout)
    tags = ways[0][2]
    for k in ("altitude", "altitude_high", "altitude_low",
              "node_altitudes"):
        assert k not in tags


# ──────────────────────────────────────────────────────────────────────
# to_osm: shared vertices (the headline invariant)
# ──────────────────────────────────────────────────────────────────────
def test_to_osm_shared_corner_uses_same_node_id():
    """Two adjacent shapes sharing a corner must reference the SAME
    nid in the OSM output.  This is the v_tgt invariant the target
    files measure."""
    layout = _make_layout()
    # Square A: (0,0)-(10,10).  Square B: (10,0)-(20,10).  Shared
    # edge along x=10, y∈[0,10].
    a = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    b = Polygon([(10, 0), (20, 0), (20, 10), (10, 10)])
    layout.shapes.append(BuiltShape(polygon=a, role=ROLE_RUNWAY, ref="A"))
    layout.shapes.append(BuiltShape(polygon=b, role=ROLE_RUNWAY, ref="B"))
    nodes, ways = _emit_and_parse(layout)

    # Find the two shared corners.  In meter space (10, 0) and
    # (10, 10) are present in both shapes.  After interning they
    # should resolve to the same nid.
    way_a_nids = set(ways[0][1])
    way_b_nids = set(ways[1][1])
    shared = way_a_nids & way_b_nids
    # At least the two shared corners must be in the intersection.
    assert len(shared) >= 2, (
        "shared-vertex invariant violated: adjacent shapes must "
        "reference the SAME node ids at coincident corners")


def test_to_osm_distinct_far_corners_get_different_node_ids():
    """Two shapes far apart should have entirely-disjoint node id
    sets (no false sharing)."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY, ref="A"))
    layout.shapes.append(BuiltShape(
        polygon=_square(1000, 1000, 10),  # 1.4 km away
        role=ROLE_RUNWAY, ref="B"))
    _, ways = _emit_and_parse(layout)
    way_a_nids = set(ways[0][1])
    way_b_nids = set(ways[1][1])
    assert not (way_a_nids & way_b_nids)


def test_to_osm_subtol_distance_clusters_nodes():
    """Vertices within SHARED_VERTEX_TOL_M (0.5 m) of each other
    cluster to the same nid."""
    layout = _make_layout()
    # Two squares whose corners are 0.1 m apart — well within tol.
    a = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    b = Polygon([(10.1, 0), (20, 0), (20, 10), (10.1, 10)])
    layout.shapes.append(BuiltShape(polygon=a, role=ROLE_RUNWAY))
    layout.shapes.append(BuiltShape(polygon=b, role=ROLE_RUNWAY))
    _, ways = _emit_and_parse(layout)
    shared = set(ways[0][1]) & set(ways[1][1])
    assert len(shared) >= 2, (
        f"vertices {SHARED_VERTEX_TOL_M} m apart should cluster")


# ──────────────────────────────────────────────────────────────────────
# to_osm: orphan-node filtering
# ──────────────────────────────────────────────────────────────────────
def test_to_osm_only_referenced_nodes_emitted():
    """Every <node> in the output must be referenced by at least one
    <way>.  No floating nodes."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY))
    layout.shapes.append(BuiltShape(
        polygon=_square(50, 0, 10), role=ROLE_RUNWAY))
    nodes, ways = _emit_and_parse(layout)
    referenced = set()
    for _, nids, _ in ways:
        referenced.update(nids)
    assert set(nodes.keys()) == referenced


# ──────────────────────────────────────────────────────────────────────
# to_osm: ring closure
# ──────────────────────────────────────────────────────────────────────
def test_to_osm_emits_closed_ring():
    """OSM closed ways repeat the first nid as the last nid."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY))
    _, ways = _emit_and_parse(layout)
    nids = ways[0][1]
    assert nids[0] == nids[-1], "ring must be closed"
    # Square → 4 unique nodes + closing repeat = 5 entries.
    assert len(nids) == 5


# ──────────────────────────────────────────────────────────────────────
# to_osm: coordinate precision
# ──────────────────────────────────────────────────────────────────────
def test_to_osm_lat_lon_precision_is_11_decimals():
    """Coordinates are emitted at .11f precision (~1 mm at the
    equator) — matches the test-fixture format."""
    layout = _make_layout()
    layout.shapes.append(BuiltShape(
        polygon=_square(0, 0, 10), role=ROLE_RUNWAY))
    with tempfile.NamedTemporaryFile(
            mode="r", suffix=".osm", delete=False) as f:
        path = f.name
    try:
        layout.to_osm(path)
        text = Path(path).read_text()
    finally:
        Path(path).unlink()
    # Find a <node> line and check the precision.
    m = re.search(r"<node[^>]*lat='([^']+)'[^>]*lon='([^']+)'", text)
    assert m is not None
    lat_str, lon_str = m.group(1), m.group(2)
    # Allow optional trailing zeros — but check that there are at
    # least 11 fractional digits.
    assert len(lat_str.split(".")[1]) == 11
    assert len(lon_str.split(".")[1]) == 11


# ──────────────────────────────────────────────────────────────────────
# Sanity: AEROWAY_FOR_ROLE is complete
# ──────────────────────────────────────────────────────────────────────
def test_aeroway_for_role_covers_all_roles():
    """Every emitted role must have an AEROWAY_FOR_ROLE entry.  If a
    new role is added without updating this map, to_osm falls back
    to "taxiway" silently — this test documents the current
    expected coverage."""
    expected_roles = {
        ROLE_RUNWAY, ROLE_PRIMARY_PARALLEL, ROLE_APRON,
        ROLE_TERMINAL, ROLE_JUNCTION,
    }
    for role in expected_roles:
        assert role in AEROWAY_FOR_ROLE
