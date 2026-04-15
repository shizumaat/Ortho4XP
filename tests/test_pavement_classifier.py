"""Unit tests for O4_Pavement_Classifier."""
from shapely.geometry import Polygon

import O4_Pavement_Classifier as PC


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def _rect(width, length, cx=0.0, cy=0.0):
    """Axis-aligned rectangle of given width (y-dim) and length
    (x-dim), centred at (cx, cy).  Returns a shapely Polygon.
    """
    dx = length / 2.0
    dy = width / 2.0
    return Polygon([
        (cx - dx, cy - dy),
        (cx + dx, cy - dy),
        (cx + dx, cy + dy),
        (cx - dx, cy + dy),
    ])


# ──────────────────────────────────────────────────────────────────────
# name_hint
# ──────────────────────────────────────────────────────────────────────
def test_name_hint_twy():
    assert PC.name_hint("TWY A") == "taxiway"
    assert PC.name_hint("twy A") == "taxiway"
    assert PC.name_hint("Taxiway Alpha") == "taxiway"
    assert PC.name_hint("TAXI B1") == "taxiway"


def test_name_hint_apron():
    assert PC.name_hint("RAMP 1") == "apron"
    assert PC.name_hint("Apron Nord") == "apron"
    assert PC.name_hint("Gate A12") == "apron"
    assert PC.name_hint("Terminal Parking") == "apron"
    assert PC.name_hint("Stand 42") == "apron"


def test_name_hint_none():
    assert PC.name_hint("") is None
    assert PC.name_hint("Aeronaval") is None
    assert PC.name_hint("AV") is None


def test_name_hint_word_boundary():
    # "TAXI" must not match substrings like "TAXIRAMP" is ambiguous,
    # but the earlier-match rule picks the taxi hit first.  More
    # importantly, unrelated words like "TERMINATE" must not match
    # "TERMINAL".
    assert PC.name_hint("TERMINATE") is None


def test_name_hint_earliest_wins():
    # When both a taxiway keyword and an apron keyword appear, the
    # one that starts earliest in the string wins.
    assert PC.name_hint("TWY A near RAMP 1") == "taxiway"
    assert PC.name_hint("RAMP 1 connects TWY A") == "apron"


# ──────────────────────────────────────────────────────────────────────
# min_rotated_bbox_dims
# ──────────────────────────────────────────────────────────────────────
def test_mrr_axis_aligned_rect():
    # 20 m wide, 200 m long strip
    poly = _rect(20, 200)
    long_s, short_s = PC.min_rotated_bbox_dims(poly)
    assert abs(long_s - 200.0) < 1e-6
    assert abs(short_s - 20.0) < 1e-6


def test_mrr_rotated_rect():
    import math
    from shapely.affinity import rotate
    poly = rotate(_rect(20, 200), 37.0, origin=(0, 0))
    long_s, short_s = PC.min_rotated_bbox_dims(poly)
    assert abs(long_s - 200.0) < 1e-3
    assert abs(short_s - 20.0) < 1e-3


def test_mrr_empty_polygon():
    assert PC.min_rotated_bbox_dims(Polygon()) == (0.0, 0.0)
    assert PC.min_rotated_bbox_dims(None) == (0.0, 0.0)


# ──────────────────────────────────────────────────────────────────────
# classify_pavement_m
# ──────────────────────────────────────────────────────────────────────
def test_classify_obvious_taxiway_by_shape():
    # 20 m × 400 m strip with no name → shape says taxiway.
    poly = _rect(20, 400)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "taxiway"
    assert result.reason == "shape:taxiway"
    assert result.aspect == 20.0


def test_classify_obvious_apron_by_shape():
    # 100 m × 150 m chunky ramp with no name → shape says apron
    # (aspect 1.5 < 4).
    poly = _rect(100, 150)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "apron"
    assert result.reason == "shape:apron"


def test_classify_long_but_too_wide_is_apron():
    # 60 m × 400 m: aspect 6.7 but short-side 60 m > 45 m envelope.
    # This is an apron entry lane, not a taxiway.
    poly = _rect(60, 400)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "apron"
    assert result.reason == "shape:apron"


def test_classify_long_but_too_narrow_is_apron():
    # 6 m × 400 m: aspect 66 but short-side 6 m < 9 m envelope.
    # This is probably a service road, not a real taxiway.
    poly = _rect(6, 400)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "apron"
    assert result.reason == "shape:apron"


def test_classify_name_overrides_shape_apron_label_on_strip():
    # Long thin strip labelled as a ramp → apron.
    poly = _rect(20, 400)
    result = PC.classify_pavement_m(poly, "RAMP 3")
    assert result.kind == "apron"
    assert result.reason == "name:apron"


def test_classify_name_overrides_shape_taxiway_label_on_blob():
    # Chunky blob labelled as a taxiway → taxiway.
    poly = _rect(100, 150)
    result = PC.classify_pavement_m(poly, "TWY B")
    assert result.kind == "taxiway"
    assert result.reason == "name:twy"


def test_classify_empty_polygon_falls_to_apron():
    result = PC.classify_pavement_m(Polygon(), "")
    assert result.kind == "apron"
    assert result.reason == "shape:apron"


def test_classify_taxiway_minimum_width_boundary():
    # Right at 9 m width — should still classify as taxiway.
    poly = _rect(9.0, 200)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "taxiway"


def test_classify_taxiway_just_below_min_width():
    # 8.9 m width — below the envelope, falls to apron.
    poly = _rect(8.9, 200)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "apron"


def test_classify_taxiway_maximum_width_boundary():
    poly = _rect(45.0, 300)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "taxiway"


def test_classify_taxiway_just_above_max_width():
    poly = _rect(45.1, 300)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "apron"


def test_classify_taxiway_min_aspect_boundary():
    # Exactly 4:1 aspect, in-envelope width.
    poly = _rect(20, 80)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "taxiway"


def test_classify_taxiway_just_below_min_aspect():
    # 3.95:1 — below the threshold, apron.
    poly = _rect(20, 79)
    result = PC.classify_pavement_m(poly, "")
    assert result.kind == "apron"
