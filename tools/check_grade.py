"""Validate elevation continuity + grade across an X-Plane patch OSM.

Usage:
    python3 tools/check_grade.py <output.osm> [--max-grade 1.5]
                                                [--proximity-m 1.0]
                                                [--edge-step-m 0.5]
                                                [--top-n 10]
                                                [--strict]

Three checks are run on the per-vertex elevations encoded in the
patch (``altitude`` / ``altitude_high`` + ``altitude_low`` /
``node_altitudes`` tags):

1. **Within-shape grade.**  Every pair of vertices on the same way
   must obey ``|de| / dist <= max_grade%`` (default 1.5%).

2. **Cross-shape proximity.**  Two vertices on different ways that
   sit within ``proximity-m`` of each other should agree on
   elevation: max permitted step is the same grade rule applied to
   the (sub-metre) distance — effectively zero step for shared
   corners.  Catches "shape A's corner says X, shape B's matching
   corner says Y" desyncs.

3. **Vertex-to-edge step.**  For every vertex, find the closest
   edge of any OTHER way within 5 m, project the vertex onto that
   edge, and compute the elevation X-Plane would render for the
   *edge* at that projected position (linear interpolation between
   the edge's two endpoint elevations).  The vertex's elevation
   should match within ``edge-step-m`` (default 0.5 m).  Catches
   the "junction triangle 1 m below the sloped taxi rect next to
   it" case.

Exit code is 1 in ``--strict`` mode if any check has any violation
beyond its threshold; 0 otherwise.  Without ``--strict`` the tool
always exits 0 and only reports counts — useful as an informational
diagnostic.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

R_EARTH = 6_378_137.0


# Import per-role grade limits from the auto_patch package (the
# single source of truth).  ``ROLE_GRADE_LIMITS`` maps role-tag to
# decimal grade (e.g. 0.015 for 1.5 %); a value of ``None`` means
# "skip the within-shape grade check for this role".  Roles
# missing from the dict fall back to ``max_grade`` (the function
# argument, default 1.5 %).
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_SRC_DIR = os.path.join(os.path.dirname(_THIS_DIR), "src")
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)
try:
    from auto_patch.config import ROLE_GRADE_LIMITS
except Exception:
    ROLE_GRADE_LIMITS: Dict[str, Optional[float]] = {}


# ── OSM parsing ─────────────────────────────────────────────────

_NODE_RE = re.compile(
    r"<node id='(-?\d+)'[^>]*lat='([^']+)'[^>]*lon='([^']+)'"
)
_WAY_RE = re.compile(r"<way id='(-?\d+)'[^>]*>(.*?)</way>", re.S)
_ND_RE = re.compile(r"<nd ref='(-?\d+)'")
_TAG_RE = re.compile(r"<tag k='([^']+)' v='([^']+)'")


@dataclass
class Way:
    wid: str
    role: str
    ref: str
    aeroway: str
    nids: List[str]               # closed ring (first repeats at end)
    elevs: List[Optional[float]]  # one per nid (closed-ring length)
    tags: Dict[str, str]


def _parse_osm(path: Path) -> Tuple[Dict[str, Tuple[float, float]],
                                    List[Way]]:
    txt = path.read_text()
    nodes: Dict[str, Tuple[float, float]] = {}
    for m in _NODE_RE.finditer(txt):
        nodes[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    ways: List[Way] = []
    for m in _WAY_RE.finditer(txt):
        wid = m.group(1)
        body = m.group(2)
        nids = _ND_RE.findall(body)
        if len(nids) < 3:
            continue
        tags = dict(_TAG_RE.findall(body))
        elevs = _derive_per_vertex_elevations(nids, tags)
        ways.append(Way(
            wid=wid,
            role=tags.get("role", ""),
            ref=tags.get("ref", ""),
            aeroway=tags.get("aeroway", ""),
            nids=nids,
            elevs=elevs,
            tags=tags,
        ))
    return nodes, ways


def _derive_per_vertex_elevations(nids: List[str], tags: Dict[str, str]
                                  ) -> List[Optional[float]]:
    """Decode the X-Plane patch elevation tags into a per-nid
    elevation list of the same length as ``nids`` (i.e. closed-ring
    length, last entry == first entry's elevation)."""
    n = len(nids)
    if "node_altitudes" in tags:
        try:
            vals = [float(x) for x in tags["node_altitudes"].split(",")]
        except ValueError:
            return [None] * n
        if len(vals) == n:
            return [float(v) for v in vals]
        return [None] * n
    if "altitude_high" in tags and "altitude_low" in tags:
        try:
            ah = float(tags["altitude_high"])
            al = float(tags["altitude_low"])
        except ValueError:
            return [None] * n
        # X-Plane patch convention for a 4-corner rect way:
        # nids[0]=hi-left, [1]=lo-left, [2]=lo-right, [3]=hi-right,
        # [4]=closing repeat of [0].  Way[-2:] = [n3, n0] is the
        # HIGH short edge; way[1:3] = [n1, n2] is the LOW short
        # edge.  See O4_Vector_Map.include_patches() for the parser.
        if n == 5:
            return [ah, al, al, ah, ah]
        # Rectangles with insertions along the long edges are
        # unsupported by X-Plane's altitude_high/low parser
        # (it expects exactly 5 nodes); flag as unknown.
        return [None] * n
    if "altitude" in tags:
        try:
            a = float(tags["altitude"])
        except ValueError:
            return [None] * n
        return [a] * n
    return [None] * n


# ── Coordinate space ────────────────────────────────────────────

def _ll_to_m_factory(nodes: Dict[str, Tuple[float, float]]):
    if not nodes:
        return lambda lat, lon: (0.0, 0.0)
    lats = [v[0] for v in nodes.values()]
    lons = [v[1] for v in nodes.values()]
    lat0 = sum(lats) / len(lats)
    lon0 = sum(lons) / len(lons)
    cos0 = math.cos(math.radians(lat0))

    def _f(lat: float, lon: float) -> Tuple[float, float]:
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return x, y
    return _f


# ── Vertex / edge tables ────────────────────────────────────────

@dataclass
class Vertex:
    way_idx: int
    nid: str
    x: float
    y: float
    elev: Optional[float]


@dataclass
class Edge:
    way_idx: int
    a: Tuple[float, float]
    b: Tuple[float, float]
    ea: float
    eb: float


def _build_vertex_edge_tables(
    nodes: Dict[str, Tuple[float, float]],
    ways: List[Way],
    ll_to_m,
) -> Tuple[List[Vertex], List[Edge]]:
    vertices: List[Vertex] = []
    edges: List[Edge] = []
    for way_idx, w in enumerate(ways):
        # Vertices (skip the closing repeat to avoid double-counting).
        for k, nid in enumerate(w.nids[:-1] if len(w.nids) > 1
                                and w.nids[0] == w.nids[-1]
                                else w.nids):
            if nid not in nodes:
                continue
            lat, lon = nodes[nid]
            x, y = ll_to_m(lat, lon)
            vertices.append(Vertex(
                way_idx=way_idx, nid=nid, x=x, y=y, elev=w.elevs[k]))
        # Edges (use the closed ring so the last edge wraps).
        ring = w.nids
        if len(ring) >= 2 and ring[0] != ring[-1]:
            ring = ring + [ring[0]]
        for k in range(len(ring) - 1):
            a_nid = ring[k]
            b_nid = ring[k + 1]
            if a_nid not in nodes or b_nid not in nodes:
                continue
            a_xy = ll_to_m(*nodes[a_nid])
            b_xy = ll_to_m(*nodes[b_nid])
            ea = w.elevs[k] if k < len(w.elevs) else None
            eb = (w.elevs[k + 1] if (k + 1) < len(w.elevs)
                  else w.elevs[0])
            if ea is None or eb is None:
                continue
            edges.append(Edge(
                way_idx=way_idx, a=a_xy, b=b_xy, ea=ea, eb=eb))
    return vertices, edges


# ── Spatial bucketing ───────────────────────────────────────────

def _bucket_vertices(vertices: List[Vertex], cell_m: float
                     ) -> Dict[Tuple[int, int], List[int]]:
    out: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, v in enumerate(vertices):
        out[(int(math.floor(v.x / cell_m)),
             int(math.floor(v.y / cell_m)))].append(i)
    return out


def _bucket_edges(edges: List[Edge], cell_m: float
                  ) -> Dict[Tuple[int, int], List[int]]:
    """Bucket each edge into every cell its bounding box touches
    (inflated by 1 cell so a query within any neighbour cell finds
    it)."""
    out: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i, e in enumerate(edges):
        x_lo = min(e.a[0], e.b[0])
        x_hi = max(e.a[0], e.b[0])
        y_lo = min(e.a[1], e.b[1])
        y_hi = max(e.a[1], e.b[1])
        cx_lo = int(math.floor(x_lo / cell_m))
        cx_hi = int(math.floor(x_hi / cell_m))
        cy_lo = int(math.floor(y_lo / cell_m))
        cy_hi = int(math.floor(y_hi / cell_m))
        for cx in range(cx_lo, cx_hi + 1):
            for cy in range(cy_lo, cy_hi + 1):
                out[(cx, cy)].append(i)
    return out


# ── Checks ──────────────────────────────────────────────────────

@dataclass
class Violation:
    grade_pct: float       # %
    excess_pct: float      # %
    distance_m: float
    de_m: float
    way_a: Way
    way_b: Way
    pt_a: Tuple[float, float]
    pt_b: Tuple[float, float]
    elev_a: float
    elev_b: float


@dataclass
class EdgeStep:
    step_m: float
    distance_m: float
    way_v: Way
    way_e: Way
    vert_pt: Tuple[float, float]
    proj_pt: Tuple[float, float]
    elev_v: float
    elev_proj: float


ELEV_ROUNDING_NOISE_M = 0.15  # patch elevations are stored at 1
                               # decimal, contributing up to 0.1 m
                               # of paired-rounding noise; the per-
                               # surface solver converges to within
                               # ~0.05 m of its grade cap before
                               # writeback.  0.15 m envelopes both
                               # without masking real (≥ 1.6 %)
                               # violations.


# X-Plane tile seams run along integer latitude / longitude lines.
# auto_patch handles them in two passes:
#
#   1. ``seam_anchors`` inserts vertices on each integer line that
#      crosses the airport and HARD-anchors them to ``dem.alt_strict``
#      — the terrain mesh in the neighbour tile pins the same points
#      to DEM, so they must agree or the patch tears at the boundary.
#   2. ``tile_cut`` later subtracts a ``half_width_m`` strip
#      (default 5 m each side, 10 m total) at every integer line.
#      The pre-cut seam vertex at the integer line is removed; the
#      resulting polygon gets new boundary vertices on the airport
#      side of the strip, exactly ``half_width_m`` away from the
#      integer line.  Those new vertices inherit altitudes by
#      resampling (nearest-neighbour or slope-projected) from the
#      pre-cut DEM-anchored ring — so they're effectively DEM-pinned
#      too, even though they no longer carry the seam tag.
#
# Both classes of vertex are immovable from the solver's POV: their
# altitudes are dictated by the DEM at the tile boundary, and the
# adjacent-tile patch + terrain mesh must agree exactly.  Within-
# shape grade between any pair touching one of these vertices is a
# function of DEM noise at the tile edge, not of solver feasibility,
# so we skip those pairs (and any triangle that touches one).
#
# Cross-shape proximity and edge/mid-edge step checks naturally
# still pass because both adjacent shapes sample the same DEM at
# the same XY, so the tile-edge vertex altitudes agree across
# shapes.
#
# Detection is geometric: a vertex is on the tile seam iff its lat
# OR its lon is within ``_SEAM_LL_TOL_DEG`` of an integer value.
# The tolerance (1e-4 °, ~11 m) covers the 5-m offset of post-cut
# boundary vertices plus slack for projection round-trip drift.
_SEAM_LL_TOL_DEG = 1e-4


def _seam_nids(nodes: Dict[str, Tuple[float, float]]) -> set:
    """Set of nids on a tile-boundary seam (integer lat or lon)."""
    out: set = set()
    for nid, (lat, lon) in nodes.items():
        if (abs(lat - round(lat)) <= _SEAM_LL_TOL_DEG
                or abs(lon - round(lon)) <= _SEAM_LL_TOL_DEG):
            out.add(nid)
    return out


def _check_plane_gradient(ways: List[Way],
                          nodes: Dict[str, Tuple[float, float]],
                          ll_to_m,
                          max_grade: float,
                          seam_nids: Optional[set] = None,
                          ) -> List[Violation]:
    """For each 3-vertex polygon (a triangle, which X-Plane renders
    as a planar surface), compute the plane's elevation gradient
    and flag if its magnitude exceeds ``max_grade``.

    A triangle can pass every vertex-PAIR grade check yet still
    have a steep perpendicular gradient: for (A=20, B=19, C=19.5)
    placed with A 300 m from B (edge grades ~0.3 %), the plane's
    gradient perpendicular to BC may be several %.  This shows up
    as a visible slope inside the triangle even though no vertex
    pair is "too steep".

    Triangles that touch any seam vertex are skipped — their plane
    is dictated by DEM-pinned corners the solver cannot move.
    """
    seam_nids = seam_nids or set()
    out: List[Violation] = []
    for w in ways:
        grade_cap = _role_grade_limit(w, max_grade)
        if grade_cap is None:
            continue
        # Pre-screen ring nids — skip the whole triangle if any
        # vertex lies on the tile seam.
        ring_nids = (w.nids[:-1] if (len(w.nids) > 1
                     and w.nids[0] == w.nids[-1])
                     else w.nids)
        if any(nid in seam_nids for nid in ring_nids):
            continue
        pts: List[Tuple[float, float, float]] = []
        for k, nid in enumerate(w.nids[:-1] if (len(w.nids) > 1
                                and w.nids[0] == w.nids[-1])
                                else w.nids):
            if nid not in nodes:
                continue
            lat, lon = nodes[nid]
            x, y = ll_to_m(lat, lon)
            e = w.elevs[k]
            if e is None:
                continue
            pts.append((x, y, e))
        if len(pts) != 3:
            continue  # only check triangles
        (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = pts
        # Plane normal via cross product of two in-plane vectors.
        ux, uy, uz = x2 - x1, y2 - y1, z2 - z1
        vx, vy, vz = x3 - x1, y3 - y1, z3 - z1
        nx = uy * vz - uz * vy
        ny = uz * vx - ux * vz
        nz = ux * vy - uy * vx
        if abs(nz) < 1e-6:
            continue  # degenerate triangle in xy plane
        # Plane: nx*X + ny*Y + nz*Z = d; dz/dx = -nx/nz, dz/dy = -ny/nz.
        gx = -nx / nz
        gy = -ny / nz
        grad = math.hypot(gx, gy)
        # Project vertices along the gradient direction to get the
        # plane's altitude swing across the triangle.  The
        # gradient check fires only when the swing exceeds the
        # grade cap allowance for that swing distance — matching
        # the within-shape pair check's rounding-noise envelope.
        gnorm = grad
        if gnorm < 1e-9:
            continue
        ghx, ghy = gx / gnorm, gy / gnorm
        proj = [(p[0] * ghx + p[1] * ghy, p[2], p)
                for p in pts]
        proj.sort()
        lo_p, lo_z, lo_pt = proj[0]
        hi_p, hi_z, hi_pt = proj[-1]
        dist_along_grad = hi_p - lo_p
        de_along_grad = abs(hi_z - lo_z)
        allowance = grade_cap * dist_along_grad + ELEV_ROUNDING_NOISE_M
        if de_along_grad <= allowance:
            continue
        out.append(Violation(
            grade_pct=grad * 100,
            excess_pct=(grad - grade_cap) * 100,
            distance_m=dist_along_grad if dist_along_grad > 0.5
                       else 1.0,
            de_m=de_along_grad,
            way_a=w, way_b=w,
            pt_a=(lo_pt[0], lo_pt[1]),
            pt_b=(hi_pt[0], hi_pt[1]),
            elev_a=lo_z, elev_b=hi_z))
    return out


WITHIN_SHAPE_MAX_PAIR_DIST_M = 60.0   # max distance between two
                                        # polygon vertices for the
                                        # pair to be checked for grade.
                                        # Triangle4XP will not create a
                                        # triangle edge between vertices
                                        # much further apart than this
                                        # (interior Steiner refinement
                                        # subdivides any large triangle),
                                        # so far-pair checks would be
                                        # false positives.  60 m covers
                                        # a typical taxi rect's diagonal
                                        # plus margin.


def _role_grade_limit(way: "Way",
                      default_grade: float) -> Optional[float]:
    """Resolve the within-shape grade limit for a way.

    Looks up the way's ``role`` tag in ``ROLE_GRADE_LIMITS`` (the
    single source of truth in ``auto_patch.config``):

    * Returns the role-specific limit (decimal, e.g. 0.015) if
      the role is present.
    * Returns ``None`` if the role is explicitly mapped to ``None``
      (skip the check — vertical structures, terrain-following
      outlines).
    * Falls back to ``default_grade`` when the role is absent
      from the dict (unknown role; use the function-argument
      cap so behaviour stays compatible with un-tagged input).
    """
    role = way.tags.get("role")
    if role in ROLE_GRADE_LIMITS:
        return ROLE_GRADE_LIMITS[role]
    return default_grade


def _pair_grade_limit(way_a: "Way", way_b: "Way",
                      default_grade: float) -> Optional[float]:
    """Resolve the cross-shape grade limit between two ways.

    Returns ``None`` (= skip the pair) when either way's role is
    on the skip-list (boundary, retaining_wall,
    groundside_pavement) — these are intentionally at terrain
    elevations or stacked at different vertical layers (a
    retaining_wall sitting at apt_elev above a tunnel_ramp at
    apt_elev−8m at the same XY is not an elevation
    "disagreement"; the wall and the ramp are different terrain
    layers by design).

    Otherwise returns the more restrictive of the two role's
    grade caps so close-but-not-shared vertices satisfy both
    surfaces' grade rules.
    """
    a = _role_grade_limit(way_a, default_grade)
    if a is None:
        return None
    b = _role_grade_limit(way_b, default_grade)
    if b is None:
        return None
    return min(a, b)


def _check_within_shape(ways: List[Way],
                        nodes: Dict[str, Tuple[float, float]],
                        ll_to_m,
                        max_grade: float,
                        seam_nids: Optional[set] = None,
                        ) -> List[Violation]:
    """Grade check between vertex pairs on the same way.

    The grade limit per way is resolved from
    ``ROLE_GRADE_LIMITS`` (auto_patch.config) — taxiway-class
    surfaces use 1.5 %, tunnel ramps use 4 %, and roles whose
    limit is ``None`` (boundary, retaining_wall,
    groundside_pavement) are skipped entirely.

    For 3-vertex polygons (triangles), every pair IS a triangle
    edge X-Plane will render — check all 3 pairs.

    For 4+-vertex polygons, check every pair within
    ``WITHIN_SHAPE_MAX_PAIR_DIST_M`` of each other.  Triangle4XP
    triangulates polygon interiors with Steiner refinement; any
    boundary vertex pair this close to each other is a plausible
    triangle edge in the resulting mesh.  Far pairs are skipped —
    Triangle4XP would interpose Steiner points and never connect
    them directly.  This catches the dominant within-shape
    failure mode (free-vertex drift near anchored neighbours)
    that the prior consecutive-only check missed.

    Pairs that include a seam vertex are skipped — those endpoints
    are HARD-anchored to DEM by the seam pipeline and the solver
    cannot move them, so a grade violation here reflects DEM noise
    along the tile boundary, not a solver bug.

    A violation requires ``|de| > grade × dist + ELEV_ROUNDING_NOISE_M``
    so single-decimal rounding doesn't produce spurious flags at
    sub-metre distances (where 0.05 m of true error rounds to 0.10 m
    of stored error).
    """
    seam_nids = seam_nids or set()
    out: List[Violation] = []
    for w in ways:
        grade_cap = _role_grade_limit(w, max_grade)
        if grade_cap is None:
            continue  # skip ROLE_GRADE_LIMITS[role] is None
        pts: List[Tuple[float, float, float, bool]] = []
        for k, nid in enumerate(w.nids[:-1] if (len(w.nids) > 1
                                and w.nids[0] == w.nids[-1])
                                else w.nids):
            if nid not in nodes:
                continue
            lat, lon = nodes[nid]
            x, y = ll_to_m(lat, lon)
            e = w.elevs[k]
            if e is None:
                continue
            pts.append((x, y, e, nid in seam_nids))
        n = len(pts)
        if n < 3:
            continue
        # Build the pairs to check.
        if n == 3:
            pairs = [(0, 1), (1, 2), (2, 0)]  # all 3 triangle edges
        else:
            # All pairs within WITHIN_SHAPE_MAX_PAIR_DIST_M.
            pairs = []
            for i in range(n):
                for j in range(i + 1, n):
                    dx = pts[i][0] - pts[j][0]
                    dy = pts[i][1] - pts[j][1]
                    if (dx * dx + dy * dy
                            <= WITHIN_SHAPE_MAX_PAIR_DIST_M
                            * WITHIN_SHAPE_MAX_PAIR_DIST_M):
                        pairs.append((i, j))
        for i, j in pairs:
            xi, yi, ei, si = pts[i]
            xj, yj, ej, sj = pts[j]
            if si or sj:
                continue  # seam-anchored endpoint — DEM controls.
            d = math.hypot(xi - xj, yi - yj)
            if d < 0.5:
                continue
            de = abs(ei - ej)
            allowance = grade_cap * d + ELEV_ROUNDING_NOISE_M
            if de <= allowance:
                continue
            grade = de / d
            out.append(Violation(
                grade_pct=grade * 100,
                excess_pct=(grade - grade_cap) * 100,
                distance_m=d,
                de_m=de,
                way_a=w, way_b=w,
                pt_a=(xi, yi), pt_b=(xj, yj),
                elev_a=ei, elev_b=ej))
    return out


SHARED_NID_TOLERANCE_M = 0.15  # rounding-precision step at a
                                # shared OSM node (1-decimal elev)


def _check_cross_shape_proximity(
    vertices: List[Vertex],
    ways: List[Way],
    proximity_m: float,
    max_grade: float,
) -> List[Violation]:
    """For every pair of vertices on DIFFERENT ways within
    ``proximity_m`` of each other, verify ``|de| / dist <= grade``.
    For sub-metre distances this is essentially "shared corners
    must agree on elevation".

    When two ways reference the SAME OSM node id, the vertices are
    geometrically identical — any non-zero elevation difference is
    a desync (we tolerate up to SHARED_NID_TOLERANCE_M for one-
    decimal rounding noise, then flag).
    """
    out: List[Violation] = []
    cell = max(proximity_m, 0.5)
    grid = _bucket_vertices(vertices, cell)
    for v_idx, v in enumerate(vertices):
        if v.elev is None:
            continue
        cx = int(math.floor(v.x / cell))
        cy = int(math.floor(v.y / cell))
        for dcx in (-1, 0, 1):
            for dcy in (-1, 0, 1):
                bucket = grid.get((cx + dcx, cy + dcy))
                if not bucket:
                    continue
                for u_idx in bucket:
                    if u_idx <= v_idx:
                        continue
                    u = vertices[u_idx]
                    if u.way_idx == v.way_idx:
                        continue
                    if u.elev is None:
                        continue
                    d = math.hypot(v.x - u.x, v.y - u.y)
                    if d > proximity_m:
                        continue
                    way_v = ways[v.way_idx]
                    way_u = ways[u.way_idx]
                    grade_cap = _pair_grade_limit(
                        way_v, way_u, max_grade)
                    if grade_cap is None:
                        continue
                    de = abs(v.elev - u.elev)
                    # Same OSM node referenced by two ways: the
                    # only valid step is rounding noise.  Don't
                    # apply the grade rule (denominator is zero).
                    if v.nid == u.nid or d < 0.05:
                        if de <= SHARED_NID_TOLERANCE_M:
                            continue
                        out.append(Violation(
                            grade_pct=float("inf"),
                            excess_pct=float("inf"),
                            distance_m=d,
                            de_m=de,
                            way_a=way_v,
                            way_b=way_u,
                            pt_a=(v.x, v.y), pt_b=(u.x, u.y),
                            elev_a=v.elev, elev_b=u.elev))
                        continue
                    allowance = grade_cap * d + ELEV_ROUNDING_NOISE_M
                    if de <= allowance:
                        continue
                    grade = de / d
                    out.append(Violation(
                        grade_pct=grade * 100,
                        excess_pct=(grade - grade_cap) * 100,
                        distance_m=d,
                        de_m=de,
                        way_a=way_v,
                        way_b=way_u,
                        pt_a=(v.x, v.y), pt_b=(u.x, u.y),
                        elev_a=v.elev, elev_b=u.elev))
    return out


def _check_vertex_to_edge_step(
    vertices: List[Vertex],
    edges: List[Edge],
    ways: List[Way],
    edge_search_m: float,
    edge_step_m: float,
) -> List[EdgeStep]:
    """For each vertex, find the closest edge of ANY OTHER way
    within ``edge_search_m``.  Project the vertex onto the edge,
    compute interpolated elevation along the edge at that point,
    and report a violation if the vertex's own elevation differs
    by more than ``edge_step_m``."""
    out: List[EdgeStep] = []
    cell = max(edge_search_m, 1.0)
    edge_grid = _bucket_edges(edges, cell)
    for v in vertices:
        if v.elev is None:
            continue
        way_v = ways[v.way_idx]
        if _role_grade_limit(way_v, 1.0) is None:
            continue  # vertex's role is on the skip-list
        cx = int(math.floor(v.x / cell))
        cy = int(math.floor(v.y / cell))
        best_d2 = edge_search_m * edge_search_m
        best: Optional[Tuple[Edge, float, float, float]] = None
        for dcx in (-1, 0, 1):
            for dcy in (-1, 0, 1):
                bucket = edge_grid.get((cx + dcx, cy + dcy))
                if not bucket:
                    continue
                for e_idx in bucket:
                    e = edges[e_idx]
                    if e.way_idx == v.way_idx:
                        continue
                    way_e = ways[e.way_idx]
                    if _role_grade_limit(way_e, 1.0) is None:
                        continue  # edge's role is on the skip-list
                    ax, ay = e.a
                    bx, by = e.b
                    dx = bx - ax
                    dy = by - ay
                    seg2 = dx * dx + dy * dy
                    if seg2 < 0.04:
                        continue
                    t = ((v.x - ax) * dx + (v.y - ay) * dy) / seg2
                    if t < 0.0:
                        t = 0.0
                    elif t > 1.0:
                        t = 1.0
                    px = ax + t * dx
                    py = ay + t * dy
                    d2 = (v.x - px) * (v.x - px) + (v.y - py) * (v.y - py)
                    if d2 < best_d2:
                        best_d2 = d2
                        best = (e, t, px, py)
        if best is None:
            continue
        e, t, px, py = best
        e_proj = e.ea + t * (e.eb - e.ea)
        step = abs(v.elev - e_proj)
        if step > edge_step_m + 1e-5:
            out.append(EdgeStep(
                step_m=step,
                distance_m=math.sqrt(best_d2),
                way_v=ways[v.way_idx],
                way_e=ways[e.way_idx],
                vert_pt=(v.x, v.y),
                proj_pt=(px, py),
                elev_v=v.elev,
                elev_proj=e_proj))
    return out


def _check_edge_midpoint_step(
    edges: List[Edge],
    ways: List[Way],
    edge_search_m: float,
    edge_step_m: float,
    samples_per_edge: int = 5,
) -> List[EdgeStep]:
    """For every edge, sample at ``samples_per_edge`` points
    (including the midpoint), compute the edge's interpolated
    elevation at that sample, then find the closest edge of any
    OTHER way and compare its interpolated elevation at the
    projected point.

    This catches the "two parallel edges drift apart in elevation"
    case that vertex-only checks miss: endpoints may agree but a
    mid-edge sample can still have a visible step if the two
    edges aren't exactly coincident (e.g. a junction edge running
    0.3 m alongside a sloped rect's long edge, with the junction's
    other endpoint dragging the midpoint elevation off the rect's
    slope at that point).
    """
    out: List[EdgeStep] = []
    cell = max(edge_search_m, 1.0)
    edge_grid = _bucket_edges(edges, cell)
    for e1 in edges:
        way_e1 = ways[e1.way_idx]
        if _role_grade_limit(way_e1, 1.0) is None:
            continue  # this edge's role is on the skip-list
        ax, ay = e1.a
        bx, by = e1.b
        dx = bx - ax
        dy = by - ay
        seg_len = math.hypot(dx, dy)
        if seg_len < 1.0:
            continue
        # Sample at interior t = 1/(N+1), 2/(N+1), ... N/(N+1).
        for k in range(1, samples_per_edge + 1):
            t = k / (samples_per_edge + 1)
            sx = ax + t * dx
            sy = ay + t * dy
            s_elev = e1.ea + t * (e1.eb - e1.ea)
            # Find closest OTHER-way edge to this sample.
            cx = int(math.floor(sx / cell))
            cy = int(math.floor(sy / cell))
            best_d2 = edge_search_m * edge_search_m
            best: Optional[Tuple[Edge, float, float, float]] = None
            for dcx in (-1, 0, 1):
                for dcy in (-1, 0, 1):
                    bucket = edge_grid.get((cx + dcx, cy + dcy))
                    if not bucket:
                        continue
                    for e2_idx in bucket:
                        e2 = edges[e2_idx]
                        if e2.way_idx == e1.way_idx:
                            continue
                        way_e2 = ways[e2.way_idx]
                        if _role_grade_limit(way_e2, 1.0) is None:
                            continue  # other edge's role on skip-list
                        e2ax, e2ay = e2.a
                        e2bx, e2by = e2.b
                        e2dx = e2bx - e2ax
                        e2dy = e2by - e2ay
                        e2seg2 = e2dx * e2dx + e2dy * e2dy
                        if e2seg2 < 0.04:
                            continue
                        tt = ((sx - e2ax) * e2dx
                              + (sy - e2ay) * e2dy) / e2seg2
                        if tt < 0.0:
                            tt = 0.0
                        elif tt > 1.0:
                            tt = 1.0
                        px = e2ax + tt * e2dx
                        py = e2ay + tt * e2dy
                        d2 = (sx - px) * (sx - px) + (sy - py) * (sy - py)
                        if d2 < best_d2:
                            best_d2 = d2
                            best = (e2, tt, px, py)
            if best is None:
                continue
            e2, tt, px, py = best
            e2_elev = e2.ea + tt * (e2.eb - e2.ea)
            step = abs(s_elev - e2_elev)
            if step > edge_step_m + 1e-5:
                out.append(EdgeStep(
                    step_m=step,
                    distance_m=math.sqrt(best_d2),
                    way_v=ways[e1.way_idx],
                    way_e=ways[e2.way_idx],
                    vert_pt=(sx, sy),
                    proj_pt=(px, py),
                    elev_v=s_elev,
                    elev_proj=e2_elev))
    return out


# ── Reporting ───────────────────────────────────────────────────

def _label(w: Way) -> str:
    return f"{w.role or '?'}/{w.ref or w.wid}"


def _print_violations(title: str, vios: List[Violation], top_n: int):
    n = len(vios)
    print(f"\n{title}: {n} violation{'s' if n != 1 else ''}")
    if not vios:
        return
    vios.sort(key=lambda v: -v.grade_pct)
    print(f"  worst {min(top_n, n)}:")
    for v in vios[:top_n]:
        print(f"    {v.grade_pct:6.2f}% (excess {v.excess_pct:+5.2f}%) "
              f"d={v.distance_m:6.2f}m |de|={v.de_m:5.2f}m  "
              f"{_label(v.way_a)} ({v.elev_a:.1f}) -> "
              f"{_label(v.way_b)} ({v.elev_b:.1f})")
    # Bucket distribution by excess.
    band_caps = [0.5, 1.0, 2.0, 5.0]
    band_counts = [0] * (len(band_caps) + 1)
    for v in vios:
        for bi, cap in enumerate(band_caps):
            if v.excess_pct < cap:
                band_counts[bi] += 1
                break
        else:
            band_counts[-1] += 1
    print("  excess distribution:")
    last = 0.0
    for bi, cap in enumerate(band_caps):
        print(f"    {last:.1f} → {cap:.1f}% over: {band_counts[bi]}")
        last = cap
    print(f"    {last:.1f}% → ∞ over: {band_counts[-1]}")


def _print_steps(title: str, steps: List[EdgeStep], top_n: int,
                 step_threshold_m: float):
    n = len(steps)
    print(f"\n{title}: {n} step{'s' if n != 1 else ''} > {step_threshold_m}m")
    if not steps:
        return
    steps.sort(key=lambda s: -s.step_m)
    print(f"  worst {min(top_n, n)}:")
    for s in steps[:top_n]:
        print(f"    step={s.step_m:5.2f}m  d={s.distance_m:5.2f}m  "
              f"vert={_label(s.way_v)} ({s.elev_v:.1f}) -> "
              f"edge={_label(s.way_e)} (proj {s.elev_proj:.1f})")


# ── Main ────────────────────────────────────────────────────────

def run_checks(
    osm_path: Path,
    max_grade_pct: float = 1.5,
    proximity_m: float = 1.0,
    edge_search_m: float = 5.0,
    edge_step_m: float = 0.5,
    top_n: int = 10,
) -> Tuple[List[Violation], List[Violation], List[EdgeStep]]:
    nodes, ways = _parse_osm(osm_path)
    ll_to_m = _ll_to_m_factory(nodes)
    vertices, edges = _build_vertex_edge_tables(nodes, ways, ll_to_m)
    max_grade = max_grade_pct / 100.0
    seam_nids = _seam_nids(nodes)

    print(f"=== Grade validation: {osm_path} ===")
    n_with_elev = sum(1 for v in vertices if v.elev is not None)
    print(f"  ways: {len(ways)} | vertices: {len(vertices)} "
          f"({n_with_elev} with elevation) | edges: {len(edges)} "
          f"| seam vertices: {len(seam_nids)}")

    within = _check_within_shape(
        ways, nodes, ll_to_m, max_grade, seam_nids=seam_nids)
    _print_violations(
        f"WITHIN-SHAPE vertex-pair grade > {max_grade_pct}%",
        within, top_n)

    plane = _check_plane_gradient(
        ways, nodes, ll_to_m, max_grade, seam_nids=seam_nids)
    _print_violations(
        f"PLANE GRADIENT (triangle surface) > {max_grade_pct}%",
        plane, top_n)
    within = within + plane

    cross = _check_cross_shape_proximity(
        vertices, ways, proximity_m, max_grade)
    _print_violations(
        f"CROSS-SHAPE proximity (≤ {proximity_m}m) "
        f"grade > {max_grade_pct}%",
        cross, top_n)

    steps = _check_vertex_to_edge_step(
        vertices, edges, ways, edge_search_m, edge_step_m)
    _print_steps(
        f"VERTEX-TO-EDGE step (within {edge_search_m}m of "
        f"another shape)",
        steps, top_n, edge_step_m)

    mid_steps = _check_edge_midpoint_step(
        edges, ways, edge_search_m, edge_step_m)
    _print_steps(
        f"MID-EDGE step (sample along each edge, compare to "
        f"nearest other-shape edge)",
        mid_steps, top_n, edge_step_m)

    return within, cross, steps + mid_steps


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("osm", type=Path,
                   help="Path to an X-Plane patch.osm file.")
    p.add_argument("--max-grade", type=float, default=1.5,
                   help="Max permitted grade in %% (default 1.5)")
    p.add_argument("--proximity-m", type=float, default=1.0,
                   help="Cross-shape proximity radius (default 1.0 m)")
    p.add_argument("--edge-search-m", type=float, default=5.0,
                   help="Vertex-to-edge search radius (default 5.0 m)")
    p.add_argument("--edge-step-m", type=float, default=0.5,
                   help="Max permitted vertex-to-edge step in m "
                        "(default 0.5)")
    p.add_argument("--top-n", type=int, default=10,
                   help="Show this many worst violations per check")
    p.add_argument("--strict", action="store_true",
                   help="Exit 1 if any check has any violation.")
    args = p.parse_args(argv)
    within, cross, steps = run_checks(
        args.osm,
        max_grade_pct=args.max_grade,
        proximity_m=args.proximity_m,
        edge_search_m=args.edge_search_m,
        edge_step_m=args.edge_step_m,
        top_n=args.top_n,
    )
    if args.strict and (within or cross or steps):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
