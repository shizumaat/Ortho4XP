"""Compare a produced OSM layout against a hand-crafted target OSM.

Usage:
    python3 tools/compare_target.py <target.osm> <output.osm> [--anchor LAT,LON]

Both files are expected to have ways tagged with ``role`` and ``aeroway``.
Ways without a ``role`` tag are skipped.  Multipolygon relations are
handled by taking their outer ring (we score per-ring, holes are
reported separately).

The comparator reports:

1. Per-role shape counts (target vs output).
2. Per-role best-IoU matching (Hungarian-like greedy).  Each target
   shape is paired with its highest-IoU same-role output shape.
3. For each matched pair: symmetric vertex-match report.
   * fraction of target vertices with an output vertex within TOL
   * fraction of output vertices with a target vertex within TOL
4. Global pass/fail verdict at the 0.5m tolerance.
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

from shapely.geometry import Polygon
from shapely.ops import unary_union

R_EARTH = 6_378_137.0
DEFAULT_TOL_M = 0.5


# ----------------------------------------------------------------------
# OSM parsing

@dataclass
class OsmShape:
    way_id: str
    role: str
    aeroway: str
    ref: str
    polygon: Polygon
    source: str  # "target" or "output"

    def label(self) -> str:
        return f"{self.role}:{self.ref or self.way_id}"


_NODE_RE = re.compile(
    r"<node id='(-?\d+)'[^>]*lat='([^']+)'[^>]*lon='([^']+)'"
)
_WAY_RE = re.compile(r"<way id='(-?\d+)'[^>]*>(.*?)</way>", re.S)
_ND_RE = re.compile(r"<nd ref='(-?\d+)'")
_TAG_RE = re.compile(r"<tag k='([^']+)' v='([^']+)'")


def _parse_osm(path: Path) -> Tuple[Dict[str, Tuple[float, float]],
                                    List[Tuple[str, List[str], Dict[str, str]]]]:
    txt = path.read_text()
    nodes: Dict[str, Tuple[float, float]] = {}
    for m in _NODE_RE.finditer(txt):
        nodes[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    ways = []
    for m in _WAY_RE.finditer(txt):
        wid = m.group(1)
        body = m.group(2)
        nds = _ND_RE.findall(body)
        tags = dict(_TAG_RE.findall(body))
        ways.append((wid, nds, tags))
    return nodes, ways


def _ll_to_m(lat: float, lon: float, lat0: float, lon0: float) -> Tuple[float, float]:
    cos0 = math.cos(math.radians(lat0))
    x = math.radians(lon - lon0) * R_EARTH * cos0
    y = math.radians(lat - lat0) * R_EARTH
    return x, y


def load_shapes(path: Path, anchor: Tuple[float, float], source: str) -> List[OsmShape]:
    nodes, ways = _parse_osm(path)
    lat0, lon0 = anchor
    shapes: List[OsmShape] = []
    for wid, nds, tags in ways:
        role = tags.get("role")
        if not role:
            continue
        if len(nds) < 4:
            continue
        coords = []
        for nid in nds:
            if nid not in nodes:
                continue
            lat, lon = nodes[nid]
            coords.append(_ll_to_m(lat, lon, lat0, lon0))
        if len(coords) < 4:
            continue
        # Close if not already closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        try:
            poly = Polygon(coords)
        except Exception:
            continue
        if poly.is_empty or not poly.is_valid:
            poly = poly.buffer(0)
            if poly.is_empty:
                continue
        if poly.geom_type != "Polygon":
            if hasattr(poly, "geoms"):
                poly = max(poly.geoms, key=lambda g: g.area)
            else:
                continue
        shapes.append(OsmShape(
            way_id=wid,
            role=role,
            aeroway=tags.get("aeroway", "?"),
            ref=tags.get("ref", ""),
            polygon=poly,
            source=source,
        ))
    return shapes


def pick_anchor(target_path: Path) -> Tuple[float, float]:
    """Use the centroid of target-file nodes as the meter-space anchor."""
    nodes, _ = _parse_osm(target_path)
    if not nodes:
        raise SystemExit(f"No nodes in {target_path}")
    lats = [lat for lat, _ in nodes.values()]
    lons = [lon for _, lon in nodes.values()]
    return sum(lats) / len(lats), sum(lons) / len(lons)


# ----------------------------------------------------------------------
# Metrics

def _iou(a: Polygon, b: Polygon) -> float:
    if a.is_empty or b.is_empty:
        return 0.0
    inter = a.intersection(b).area
    union = a.union(b).area
    if union <= 0:
        return 0.0
    return inter / union


def _vertices(poly: Polygon) -> List[Tuple[float, float]]:
    out = list(poly.exterior.coords)
    if out[0] == out[-1]:
        out = out[:-1]
    for interior in poly.interiors:
        ic = list(interior.coords)
        if ic[0] == ic[-1]:
            ic = ic[:-1]
        out.extend(ic)
    return out


def vertex_match(a: Polygon, b: Polygon, tol: float = DEFAULT_TOL_M
                 ) -> Tuple[float, float, float, float]:
    """Return (frac_a_matched, frac_b_matched, max_a_to_b, max_b_to_a)."""
    va = _vertices(a)
    vb = _vertices(b)
    if not va or not vb:
        return 0.0, 0.0, float("inf"), float("inf")

    def _dists(src, dst):
        out = []
        for sx, sy in src:
            best = min((sx - dx) ** 2 + (sy - dy) ** 2 for dx, dy in dst)
            out.append(math.sqrt(best))
        return out

    da = _dists(va, vb)
    db = _dists(vb, va)
    frac_a = sum(1 for d in da if d <= tol) / len(da)
    frac_b = sum(1 for d in db if d <= tol) / len(db)
    return frac_a, frac_b, max(da), max(db)


@dataclass
class Pair:
    target: Optional[OsmShape]
    output: Optional[OsmShape]
    iou: float = 0.0
    frac_t: float = 0.0
    frac_o: float = 0.0
    max_t_to_o: float = float("inf")
    max_o_to_t: float = float("inf")


def match_by_role(targets: List[OsmShape], outputs: List[OsmShape],
                  iou_floor: float = 0.10) -> List[Pair]:
    """Greedy best-IoU match within same role.  Leftovers reported separately."""
    pairs: List[Pair] = []
    used_outputs = set()

    # Rank candidate pairs by IoU
    candidates = []
    for ti, t in enumerate(targets):
        for oi, o in enumerate(outputs):
            if o.role != t.role:
                continue
            iou = _iou(t.polygon, o.polygon)
            if iou < iou_floor:
                continue
            candidates.append((iou, ti, oi))
    candidates.sort(reverse=True)

    taken_targets = set()
    for iou, ti, oi in candidates:
        if ti in taken_targets or oi in used_outputs:
            continue
        pair = Pair(target=targets[ti], output=outputs[oi], iou=iou)
        pair.frac_t, pair.frac_o, pair.max_t_to_o, pair.max_o_to_t = \
            vertex_match(targets[ti].polygon, outputs[oi].polygon)
        pairs.append(pair)
        taken_targets.add(ti)
        used_outputs.add(oi)

    for ti, t in enumerate(targets):
        if ti not in taken_targets:
            pairs.append(Pair(target=t, output=None))
    for oi, o in enumerate(outputs):
        if oi not in used_outputs:
            pairs.append(Pair(target=None, output=o))
    return pairs


# ----------------------------------------------------------------------
# Reporting

def _summarize_role(role: str, pairs: List[Pair], tol: float) -> str:
    matched = [p for p in pairs if p.target and p.output]
    missed = [p for p in pairs if p.target and not p.output]
    spurious = [p for p in pairs if p.output and not p.target]

    if matched:
        ious = sorted(p.iou for p in matched)
        avg_iou = sum(ious) / len(ious)
        pass_t = sum(1 for p in matched if p.max_t_to_o <= tol)
        pass_o = sum(1 for p in matched if p.max_o_to_t <= tol)
        frac_t = sum(p.frac_t for p in matched) / len(matched)
        frac_o = sum(p.frac_o for p in matched) / len(matched)
    else:
        avg_iou = 0.0
        pass_t = pass_o = 0
        frac_t = frac_o = 0.0

    n_t = len(matched) + len(missed)
    n_o = len(matched) + len(spurious)

    return (
        f"  {role:<22}  target={n_t:>2}  output={n_o:>2}  "
        f"matched={len(matched):>2}  missed={len(missed):>2}  "
        f"spurious={len(spurious):>2}  "
        f"avgIoU={avg_iou:.2f}  "
        f"v_tgt@{tol}m={frac_t*100:5.1f}%  "
        f"v_out@{tol}m={frac_o*100:5.1f}%  "
        f"pass_t={pass_t}/{len(matched)}  "
        f"pass_o={pass_o}/{len(matched)}"
    )


def _detail_pair(p: Pair, tol: float) -> str:
    if p.target and p.output:
        return (
            f"    [{p.iou:.2f}] {p.target.label():<28} <-> "
            f"{p.output.label():<28}  "
            f"max_t_to_o={p.max_t_to_o:5.2f} m  max_o_to_t={p.max_o_to_t:5.2f} m  "
            f"v_t={p.frac_t*100:5.1f}%  v_o={p.frac_o*100:5.1f}%"
        )
    if p.target and not p.output:
        return f"    MISSED   {p.target.label():<28} (area {p.target.polygon.area:.0f} m²)"
    if p.output and not p.target:
        return f"    SPURIOUS {p.output.label():<28} (area {p.output.polygon.area:.0f} m²)"
    return ""


def report(target_shapes: List[OsmShape], output_shapes: List[OsmShape],
           tol: float = DEFAULT_TOL_M, verbose: bool = False) -> Dict[str, List[Pair]]:
    by_role_t = defaultdict(list)
    for s in target_shapes:
        by_role_t[s.role].append(s)
    by_role_o = defaultdict(list)
    for s in output_shapes:
        by_role_o[s.role].append(s)

    all_roles = sorted(set(by_role_t) | set(by_role_o))
    per_role: Dict[str, List[Pair]] = {}

    print(f"=== Per-role summary (tol={tol} m) ===")
    print(f"  {'role':<22}  {'n_t':>5}  {'n_o':>5}  matched  missed  sp.  avgIoU  v_tgt  v_out  pass_t  pass_o")
    total_t = total_o = total_matched = total_pass_t = total_pass_o = 0
    for role in all_roles:
        pairs = match_by_role(by_role_t.get(role, []), by_role_o.get(role, []))
        per_role[role] = pairs
        print(_summarize_role(role, pairs, tol))
        matched = [p for p in pairs if p.target and p.output]
        total_t += sum(1 for p in pairs if p.target)
        total_o += sum(1 for p in pairs if p.output)
        total_matched += len(matched)
        total_pass_t += sum(1 for p in matched if p.max_t_to_o <= tol)
        total_pass_o += sum(1 for p in matched if p.max_o_to_t <= tol)

    print()
    print(f"TOTALS  targets={total_t}  outputs={total_o}  matched={total_matched}  "
          f"pass_t={total_pass_t}  pass_o={total_pass_o}")

    if verbose:
        for role in all_roles:
            print(f"\n--- {role} ---")
            for p in sorted(per_role[role], key=lambda x: -(x.iou if x.target and x.output else 0)):
                line = _detail_pair(p, tol)
                if line:
                    print(line)

    return per_role


# ----------------------------------------------------------------------
# CLI

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("target", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--anchor", type=str, default=None,
                    help="Override meter-space anchor as LAT,LON. "
                         "Default: centroid of target nodes.")
    ap.add_argument("--tol", type=float, default=DEFAULT_TOL_M,
                    help="Vertex-match tolerance in meters (default 0.5)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    if args.anchor:
        lat0, lon0 = (float(x) for x in args.anchor.split(","))
    else:
        lat0, lon0 = pick_anchor(args.target)

    target_shapes = load_shapes(args.target, (lat0, lon0), "target")
    output_shapes = load_shapes(args.output, (lat0, lon0), "output")
    print(f"Target: {args.target} ({len(target_shapes)} shapes)")
    print(f"Output: {args.output} ({len(output_shapes)} shapes)")
    print(f"Anchor: lat={lat0:.6f} lon={lon0:.6f}")
    print()

    report(target_shapes, output_shapes, tol=args.tol, verbose=args.verbose)


if __name__ == "__main__":
    main()
