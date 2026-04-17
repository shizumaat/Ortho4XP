"""Snap every non-runway vertex in a target OSM to its nearest
apt.dat pavement edge.

Writes a new file ``<basename>_snapped.osm`` next to the input so
the original stays intact.  The user reviews the snapped file and,
if valid, we use it as the new reference target.

Rules:
* Runway vertices are NOT moved (they're already at apt.dat
  row-100 + blast-pad corners).
* Every other vertex is snapped to the nearest point on the
  apt.dat pavement MultiPolygon boundary, PROVIDED the snap
  distance is ≤ ``MAX_SNAP_DIST_M`` (default 25 m) — farther
  than that and we assume the user deliberately drew a vertex in
  the interior (e.g. midpoint of a long edge) and leave it alone.
* Shared nodes move ONCE (the node is snapped, every way that
  references it follows).
* A per-node diagnostic (distance, moved/skipped) is printed so
  you can audit.

Usage:
    python3 tools/snap_target_to_apt_dat.py SPJC tests/fixtures/SPJC_target.osm
    python3 tools/snap_target_to_apt_dat.py SPLP tests/fixtures/SPLP_target.osm
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from shapely.geometry import Point, Polygon
from shapely.ops import nearest_points, transform as shp_transform, unary_union

import O4_Apt_Dat_Reader as APR

R_EARTH = 6_378_137.0
MAX_SNAP_DIST_M = 25.0


def _projection(anchor):
    lat0, lon0 = anchor
    cos0 = math.cos(math.radians(lat0))
    def to_m(lon, lat, z=None):
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return (x, y) if z is None else (x, y, z)
    def to_ll(x, y):
        lon = lon0 + math.degrees(x / (R_EARTH * cos0))
        lat = lat0 + math.degrees(y / R_EARTH)
        return lat, lon
    return to_m, to_ll


def _apt_pavement_boundary_m(apt: APR.Airport, to_m):
    """Return the combined pavement boundary (runways + taxiways +
    aprons) as a single geometry in meter space."""
    pav_polys = []
    # Pavements from apt.dat rows 110+
    for pav in apt.pavements:
        if pav.polygon is None or pav.polygon.is_empty:
            continue
        pm = shp_transform(to_m, pav.polygon)
        if pm.is_empty:
            continue
        if pm.geom_type == "Polygon":
            pav_polys.append(pm)
        else:
            pav_polys.extend(g for g in getattr(pm, "geoms", [])
                             if g.geom_type == "Polygon")
    # Runway footprints (row 100 rect + blast pads)
    for r in apt.runways:
        lat0 = (r.lat_a + r.lat_b) / 2
        lon0 = (r.lon_a + r.lon_b) / 2
        # Use the main to_m, not a local projection
        ax, ay = to_m(r.lon_a, r.lat_a)
        bx, by = to_m(r.lon_b, r.lat_b)
        dx, dy = bx - ax, by - ay
        mag = math.hypot(dx, dy)
        if mag < 1.0:
            continue
        ux, uy = dx / mag, dy / mag
        a_extra = r.blast_a_m or 0.0
        b_extra = r.blast_b_m or 0.0
        ax2 = ax - ux * a_extra; ay2 = ay - uy * a_extra
        bx2 = bx + ux * b_extra; by2 = by + uy * b_extra
        px, py = -uy, ux
        half = r.width_m / 2.0
        pav_polys.append(Polygon([
            (ax2 + px*half, ay2 + py*half),
            (bx2 + px*half, by2 + py*half),
            (bx2 - px*half, by2 - py*half),
            (ax2 - px*half, ay2 - py*half),
        ]))
    if not pav_polys:
        raise SystemExit("No pavement polygons parsed from apt.dat")
    pav = unary_union(pav_polys).buffer(0)
    return pav, pav.boundary


def _parse_osm(txt: str):
    node_re = re.compile(
        r"(<node id='(-?\d+)'[^>]*lat=')([^']+)('[^>]*lon=')([^']+)('[^/]*/>)"
    )
    way_re = re.compile(r"<way id='(-?\d+)'[^>]*>(.*?)</way>", re.S)
    nd_re = re.compile(r"<nd ref='(-?\d+)'")
    tag_re = re.compile(r"<tag k='([^']+)' v='([^']+)'")
    nodes = {}
    for m in node_re.finditer(txt):
        nid = m.group(2)
        lat = float(m.group(3))
        lon = float(m.group(5))
        nodes[nid] = (lat, lon)
    ways = []
    for m in way_re.finditer(txt):
        wid = m.group(1)
        body = m.group(2)
        nds = nd_re.findall(body)
        tags = dict(tag_re.findall(body))
        ways.append((wid, nds, tags))
    return nodes, ways


def _anchor(apt: APR.Airport):
    if apt.runways:
        r = apt.runways[0]
        return ((r.lat_a + r.lat_b)/2, (r.lon_a + r.lon_b)/2)
    return (0.0, 0.0)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("icao")
    ap.add_argument("target_path", type=Path)
    ap.add_argument("--xplane", default="/Users/noah/X-Plane 12")
    ap.add_argument("--max-snap", type=float, default=MAX_SNAP_DIST_M)
    ap.add_argument("--out", type=Path, default=None,
                    help="Output path; default: <target>_snapped.osm")
    args = ap.parse_args(argv)

    out_path = args.out or args.target_path.with_name(
        args.target_path.stem + "_snapped.osm")

    apt_path = APR.find_airport_apt_dat(args.xplane, args.icao)
    apt = APR.load_airport(apt_path, args.icao)
    anchor = _anchor(apt)
    to_m, to_ll = _projection(anchor)
    _, boundary = _apt_pavement_boundary_m(apt, to_m)
    print(f"Loaded {args.icao}: boundary length = {boundary.length:.0f} m")

    txt = args.target_path.read_text()
    nodes, ways = _parse_osm(txt)
    print(f"Target: {len(nodes)} nodes, {len(ways)} ways")

    # Find runway node ids — these we DO NOT move
    runway_nodes: Set[str] = set()
    for wid, nds, tags in ways:
        if tags.get("role") == "runway":
            runway_nodes.update(nds)
    print(f"  runway nodes (frozen): {len(runway_nodes)}")

    # Snap each non-runway node
    moves: Dict[str, Tuple[float, float]] = {}  # nid -> (new_lat, new_lon)
    skipped = 0
    hist = {"lt1": 0, "lt5": 0, "lt10": 0, "lt25": 0}
    for nid, (lat, lon) in nodes.items():
        if nid in runway_nodes:
            continue
        x, y = to_m(lon, lat)
        p = Point(x, y)
        np_pt, _ = nearest_points(boundary, p)
        d = p.distance(np_pt)
        if d > args.max_snap:
            skipped += 1
            continue
        new_lat, new_lon = to_ll(np_pt.x, np_pt.y)
        moves[nid] = (new_lat, new_lon)
        if d < 1.0:
            hist["lt1"] += 1
        elif d < 5.0:
            hist["lt5"] += 1
        elif d < 10.0:
            hist["lt10"] += 1
        else:
            hist["lt25"] += 1

    print(f"  snapped: {len(moves)} nodes")
    print(f"    < 1 m  : {hist['lt1']}")
    print(f"    1-5 m  : {hist['lt5']}")
    print(f"    5-10 m : {hist['lt10']}")
    print(f"    10-25 m: {hist['lt25']}")
    print(f"  skipped (> {args.max_snap} m from any edge): {skipped}")

    # Rewrite nodes in txt
    def _sub_node(m):
        nid = m.group(2)
        if nid not in moves:
            return m.group(0)
        lat, lon = moves[nid]
        return (m.group(1) + f"{lat:.11f}" + m.group(4)
                + f"{lon:.11f}" + m.group(6))

    node_re = re.compile(
        r"(<node id='(-?\d+)'[^>]*lat=')([^']+)('[^>]*lon=')([^']+)('[^/]*/>)"
    )
    new_txt = node_re.sub(_sub_node, txt)
    out_path.write_text(new_txt)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
