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
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

R_EARTH = 6_378_137.0


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


ELEV_ROUNDING_NOISE_M = 0.1  # patch elevations are stored at 1
                              # decimal -> worst-case rounding noise
                              # between any two elevs is 0.1 m


def _check_within_shape(ways: List[Way],
                        nodes: Dict[str, Tuple[float, float]],
                        ll_to_m,
                        max_grade: float) -> List[Violation]:
    """Pairwise grade between every two vertices on the same way.

    A violation requires ``|de| > grade × dist + ELEV_ROUNDING_NOISE_M``
    so single-decimal rounding doesn't produce spurious flags at
    sub-metre distances (where 0.05 m of true error rounds to 0.10 m
    of stored error).
    """
    out: List[Violation] = []
    for w in ways:
        # Build (x, y, elev) list for this way.
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
        n = len(pts)
        for i in range(n):
            xi, yi, ei = pts[i]
            for j in range(i + 1, n):
                xj, yj, ej = pts[j]
                d = math.hypot(xi - xj, yi - yj)
                if d < 0.5:
                    continue
                de = abs(ei - ej)
                allowance = max_grade * d + ELEV_ROUNDING_NOISE_M
                if de <= allowance:
                    continue
                grade = de / d
                out.append(Violation(
                    grade_pct=grade * 100,
                    excess_pct=(grade - max_grade) * 100,
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
                            way_a=ways[v.way_idx],
                            way_b=ways[u.way_idx],
                            pt_a=(v.x, v.y), pt_b=(u.x, u.y),
                            elev_a=v.elev, elev_b=u.elev))
                        continue
                    allowance = max_grade * d + ELEV_ROUNDING_NOISE_M
                    if de <= allowance:
                        continue
                    grade = de / d
                    out.append(Violation(
                        grade_pct=grade * 100,
                        excess_pct=(grade - max_grade) * 100,
                        distance_m=d,
                        de_m=de,
                        way_a=ways[v.way_idx],
                        way_b=ways[u.way_idx],
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

    print(f"=== Grade validation: {osm_path} ===")
    n_with_elev = sum(1 for v in vertices if v.elev is not None)
    print(f"  ways: {len(ways)} | vertices: {len(vertices)} "
          f"({n_with_elev} with elevation) | edges: {len(edges)}")

    within = _check_within_shape(ways, nodes, ll_to_m, max_grade)
    _print_violations(
        f"WITHIN-SHAPE grade > {max_grade_pct}%",
        within, top_n)

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

    return within, cross, steps


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
