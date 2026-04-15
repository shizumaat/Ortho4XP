"""Classify apt.dat pavement polygons as taxiway or apron.

A :class:`O4_Apt_Dat_Reader.Pavement` polygon does not carry a
semantic type — the row 110 header only stores a free-form label like
``"TWY A"``, ``"RAMP 1"``, ``"Aeronaval"``, or nothing at all.  To
drive per-surface grade rules (1.0 % for aprons, 1.5 % for taxiways)
and to choose between triangulated vs rectangular emission, we need a
reliable taxiway/apron tag.

Two signals combine:

1. **Name hint.**  Substrings ``TWY``, ``TAXIWAY``, ``TAXI`` (case-
   insensitive, word-boundaried) strongly suggest taxiway; ``APRON``,
   ``RAMP``, ``STAND``, ``GATE``, ``PARKING``, ``TERMINAL`` suggest
   apron.
2. **Shape hint.**  Fit a minimum-area rotated bounding rectangle;
   measure its long side *L* and short side *W*.  Taxiways are long
   thin strips (``L/W ≥ 4``) whose short side lies within the
   taxiway-class width envelope (9 – 45 m).  Anything that fails
   either check is an apron.

When name and shape disagree, **the name wins** — a pavement labelled
``"RAMP 1"`` is an apron even if it happens to be long and thin, and a
pavement labelled ``"TWY A"`` is a taxiway even if it has branches
that make it look chunky from a bounding-box perspective.

The module is pure: given a polygon (in meter coordinates) and a
name string, it returns a :class:`PavementKind`.  No file I/O, no
side effects.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Optional, Tuple

from shapely.geometry import Polygon


# ──────────────────────────────────────────────────────────────────────
# Tunables
# ──────────────────────────────────────────────────────────────────────

# Taxiway short-side (width) envelope in meters.  FAA ADG-I is 7.5 m,
# ADG-VI is 30.5 m; apt.dat pavements often include the shoulder as
# part of the polygon width, so the envelope is padded slightly.  A
# pavement whose min-rotated-rect short side falls outside this band
# is NOT a taxiway even if its aspect ratio is high.
TAXIWAY_MIN_WIDTH_M = 9.0
TAXIWAY_MAX_WIDTH_M = 45.0

# Minimum aspect ratio (long/short) for a shape-classified taxiway.
# Four is the visual threshold where "strip" becomes obviously
# correct — anything chunkier than 4:1 is almost always a parking
# apron with a narrow entry throat.
TAXIWAY_MIN_ASPECT = 4.0


# Regex patterns.  Word-boundaried so ``TERMINAL`` does not match
# ``TWY-T-EXTENSION`` and ``TAXI`` does not match ``TAXIRAMP`` etc.
_TWY_NAME_PATS = [
    re.compile(r"(?i)\bTWY\b"),
    re.compile(r"(?i)\bTAXIWAY[S]?\b"),
    re.compile(r"(?i)\bTAXI\b"),
]

_APRON_NAME_PATS = [
    re.compile(r"(?i)\bAPRON[S]?\b"),
    re.compile(r"(?i)\bRAMP[S]?\b"),
    re.compile(r"(?i)\bSTAND[S]?\b"),
    re.compile(r"(?i)\bGATE[S]?\b"),
    re.compile(r"(?i)\bPARKING\b"),
    re.compile(r"(?i)\bTERMINAL\b"),
]


# ──────────────────────────────────────────────────────────────────────
# Data class
# ──────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class PavementKind:
    """Result of classifying one pavement polygon."""
    kind: str          # "taxiway" or "apron"
    long_m: float      # long side of the min-rotated-rect
    short_m: float     # short side of the min-rotated-rect
    aspect: float      # long_m / short_m (0 if short_m == 0)
    reason: str        # "name:twy" | "name:apron" |
                       # "shape:taxiway" | "shape:apron"


# ──────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────
def name_hint(name: str) -> Optional[str]:
    """Return ``"taxiway"``, ``"apron"``, or ``None`` from the free-
    form pavement label.  Case-insensitive, word-boundaried.

    If the name mentions both a taxiway and apron keyword, the earlier
    keyword in the string wins (the label format most pavements use is
    ``"<KIND> <ID>"`` so the first match is almost always the intended
    one).
    """
    if not name:
        return None
    earliest_kind = None
    earliest_pos = len(name) + 1
    for pat in _TWY_NAME_PATS:
        m = pat.search(name)
        if m is not None and m.start() < earliest_pos:
            earliest_pos = m.start()
            earliest_kind = "taxiway"
    for pat in _APRON_NAME_PATS:
        m = pat.search(name)
        if m is not None and m.start() < earliest_pos:
            earliest_pos = m.start()
            earliest_kind = "apron"
    return earliest_kind


def min_rotated_bbox_dims(polygon: Polygon) -> Tuple[float, float]:
    """Return ``(long_side, short_side)`` of the minimum-area rotated
    bounding rectangle of ``polygon``, in the polygon's own coordinate
    units.

    Returns ``(0.0, 0.0)`` for an empty or invalid polygon.
    """
    if polygon is None or polygon.is_empty:
        return (0.0, 0.0)
    try:
        mrr = polygon.minimum_rotated_rectangle
    except Exception:
        return (0.0, 0.0)
    if mrr is None or mrr.is_empty:
        return (0.0, 0.0)
    if not hasattr(mrr, "exterior"):
        return (0.0, 0.0)
    coords = list(mrr.exterior.coords)
    if len(coords) < 5:
        return (0.0, 0.0)
    # MRR has 4 unique vertices then closes; edges 0-1 and 1-2 are
    # perpendicular and are the two distinct side lengths.
    e1 = math.hypot(coords[1][0] - coords[0][0],
                    coords[1][1] - coords[0][1])
    e2 = math.hypot(coords[2][0] - coords[1][0],
                    coords[2][1] - coords[1][1])
    long_s = max(e1, e2)
    short_s = min(e1, e2)
    return (long_s, short_s)


def classify_pavement_m(polygon: Polygon, name: str = "") -> PavementKind:
    """Classify one apt.dat pavement polygon.

    Args:
        polygon: pavement polygon in METER coordinates.  Passing a
            lat/lon polygon will give meaningless short/long sides
            and almost always classify as apron (degrees are tiny).
        name: free-form pavement label from the apt.dat row-110
            header.  May be empty.

    Returns:
        A :class:`PavementKind` describing the classification, the
        min-rotated-rect dimensions, and which rule fired.
    """
    long_s, short_s = min_rotated_bbox_dims(polygon)
    aspect = (long_s / short_s) if short_s > 0 else 0.0

    # 1. Name hint wins when present.
    hint = name_hint(name)
    if hint == "taxiway":
        return PavementKind("taxiway", long_s, short_s, aspect, "name:twy")
    if hint == "apron":
        return PavementKind("apron", long_s, short_s, aspect, "name:apron")

    # 2. Shape fallback.
    if (aspect >= TAXIWAY_MIN_ASPECT
            and TAXIWAY_MIN_WIDTH_M <= short_s <= TAXIWAY_MAX_WIDTH_M):
        return PavementKind("taxiway", long_s, short_s, aspect,
                            "shape:taxiway")
    return PavementKind("apron", long_s, short_s, aspect, "shape:apron")
