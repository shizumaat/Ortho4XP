"""Snap every non-runway vertex in a target OSM to its nearest
PAVEMENT edge.

Pavement = apt.dat row-110 polygons ∪ runway footprints (rows
100 + blast pads) ∪ DSF draped pavement (filtered the same way
the build pipeline does — see O4_DSF_Reader and
``DSF_OVERLAY_FRAC``/``DSF_AIRPORT_RADIUS_M``).  Per user
2026-04-27: the build pipeline now incorporates DSF, so the
target snapping must use the same coverage set or the comparison
target won't match what the pipeline could legitimately produce.

Writes a new file ``<basename>_snapped.osm`` next to the input so
the original stays intact.  The user reviews the snapped file and,
if valid, we use it as the new reference target.

Rules:
* Runway vertices are NOT moved (they're already at apt.dat
  row-100 + blast-pad corners).
* Every other vertex is snapped to the nearest point on the
  pavement union boundary, PROVIDED the snap distance is
  ≤ ``MAX_SNAP_DIST_M`` (default 25 m) — farther than that and
  we assume the user deliberately drew a vertex in the interior
  (e.g. midpoint of a long edge) and leave it alone.
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

from auto_patch import apt_dat_reader as APR

try:
    from auto_patch import dsf_reader as _DSFR
except Exception:
    _DSFR = None

R_EARTH = 6_378_137.0
MAX_SNAP_DIST_M = 25.0
# DSF filter constants: mirror the build pipeline so the snapped
# target reflects the same coverage the pipeline can produce.
DSF_OVERLAY_FRAC = 0.80
DSF_AIRPORT_RADIUS_M = 5_000.0


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


def _apt_pavement_boundary_m(apt: APR.Airport, to_m,
                              xplane_root: str, icao: str,
                              anchor: Tuple[float, float]):
    """Return (pav_union, pav_boundary, pav_vertices) in meter
    space.  ``pav_vertices`` is the full set of pavement vertex
    coords (exterior + interior rings of each apt.dat pavement
    polygon, runway corners, AND surviving DSF pavement
    polygons) — used for vertex-preferred snap.

    DSF polygons are filtered the same way the build pipeline
    does (see O4_Airport_Pavement_Builder):
      • distance gate: drop polygons > DSF_AIRPORT_RADIUS_M from
        the airport bbox;
      • overlay drop: a polygon ≥ DSF_OVERLAY_FRAC inside the
        existing apt.dat pavement is decorative paint, not new
        pavement — drop it.
    """
    apt_pav_polys = []
    # Pavements from apt.dat rows 110+
    for pav in apt.pavements:
        if pav.polygon is None or pav.polygon.is_empty:
            continue
        pm = shp_transform(to_m, pav.polygon)
        if pm.is_empty:
            continue
        if pm.geom_type == "Polygon":
            apt_pav_polys.append(pm)
        else:
            apt_pav_polys.extend(g for g in getattr(pm, "geoms", [])
                                 if g.geom_type == "Polygon")
    # Runway footprints (row 100 rect + blast pads)
    for r in apt.runways:
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
        apt_pav_polys.append(Polygon([
            (ax2 + px*half, ay2 + py*half),
            (bx2 + px*half, by2 + py*half),
            (bx2 - px*half, by2 - py*half),
            (ax2 - px*half, ay2 - py*half),
        ]))
    if not apt_pav_polys:
        raise SystemExit("No pavement polygons parsed from apt.dat")
    pav_polys = list(apt_pav_polys)
    apt_pav_union = unary_union(apt_pav_polys).buffer(0)
    # Compute airport bbox for DSF distance gate.
    apt_bbox_m = None
    try:
        bx_min, by_min, bx_max, by_max = apt_pav_union.bounds
        apt_bbox_m = (bx_min - DSF_AIRPORT_RADIUS_M,
                      by_min - DSF_AIRPORT_RADIUS_M,
                      bx_max + DSF_AIRPORT_RADIUS_M,
                      by_max + DSF_AIRPORT_RADIUS_M)
    except Exception:
        pass
    # DSF pavement (apply same overlay/distance filters as build).
    n_dsf_kept = n_dsf_overlay = n_dsf_far = 0
    if _DSFR is not None:
        try:
            seen_dsf = set()
            for ad in APR.find_all_airport_apt_dats(xplane_root, icao):
                dsf = _DSFR.find_associated_dsf(ad, anchor[0], anchor[1])
                if dsf is None or dsf in seen_dsf:
                    continue
                seen_dsf.add(dsf)
                for ring in _DSFR.read_dsf_pavements(dsf):
                    if len(ring) < 3:
                        continue
                    try:
                        poly_ll = Polygon(
                            [(lon, lat) for (lon, lat) in ring])
                        if not poly_ll.is_valid:
                            poly_ll = poly_ll.buffer(0)
                        if (poly_ll.is_empty
                                or poly_ll.geom_type != "Polygon"):
                            continue
                        pm = shp_transform(to_m, poly_ll)
                        if pm.is_empty or pm.geom_type != "Polygon":
                            continue
                        if apt_bbox_m is not None:
                            px_min, py_min, px_max, py_max = pm.bounds
                            if (px_max < apt_bbox_m[0]
                                    or px_min > apt_bbox_m[2]
                                    or py_max < apt_bbox_m[1]
                                    or py_min > apt_bbox_m[3]):
                                n_dsf_far += 1
                                continue
                        try:
                            inter = pm.intersection(apt_pav_union).area
                            if (pm.area > 0
                                    and inter / pm.area
                                    >= DSF_OVERLAY_FRAC):
                                n_dsf_overlay += 1
                                continue
                        except Exception:
                            pass
                        pav_polys.append(pm)
                        n_dsf_kept += 1
                    except Exception:
                        continue
        except Exception:
            pass
    if n_dsf_kept or n_dsf_overlay or n_dsf_far:
        print(f"  DSF pavement: {n_dsf_kept} kept, "
              f"{n_dsf_overlay} dropped as overlay, "
              f"{n_dsf_far} dropped as off-airport.")
    # Collect every pavement-polygon vertex (apt.dat + DSF + runway).
    vertices: List[Tuple[float, float]] = []
    for pp in pav_polys:
        if pp.is_empty or pp.geom_type != "Polygon":
            continue
        ec = list(pp.exterior.coords)
        if ec and ec[0] == ec[-1]:
            ec = ec[:-1]
        vertices.extend(ec)
        for ring in pp.interiors:
            rc = list(ring.coords)
            if rc and rc[0] == rc[-1]:
                rc = rc[:-1]
            vertices.extend(rc)
    pav = unary_union(pav_polys).buffer(0)
    return pav, pav.boundary, vertices


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
    ap.add_argument("--force", action="store_true",
                    help="Overwrite existing output file (default: refuse).")
    args = ap.parse_args(argv)

    out_path = args.out or args.target_path.with_name(
        args.target_path.stem + "_snapped.osm")
    if out_path.exists() and not args.force:
        raise SystemExit(
            f"Refusing to overwrite existing {out_path}.\n"
            f"Pass --force if you really want to clobber it "
            f"(existing manual edits to the snapped file will be lost)."
        )

    apt_path = APR.find_airport_apt_dat(args.xplane, args.icao)
    apt = APR.load_airport(apt_path, args.icao)
    anchor = _anchor(apt)
    to_m, to_ll = _projection(anchor)
    _, boundary, apt_vertices = _apt_pavement_boundary_m(
        apt, to_m, args.xplane, args.icao, anchor)
    print(f"Loaded {args.icao}: boundary length = {boundary.length:.0f} m,"
          f" vertices = {len(apt_vertices)}")
    VERTEX_SNAP_M = 8.0  # prefer vertex within 8 m of target node

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
        # Stage 1: prefer nearest apt.dat VERTEX within VERTEX_SNAP_M
        # (vertex-level match between target & output is the goal).
        best_vx = best_vy = None
        best_vd = VERTEX_SNAP_M
        for (vx, vy) in apt_vertices:
            d = math.hypot(x - vx, y - vy)
            if d < best_vd:
                best_vd = d
                best_vx, best_vy = vx, vy
        if best_vx is not None:
            new_lat, new_lon = to_ll(best_vx, best_vy)
            moves[nid] = (new_lat, new_lon)
            d = best_vd
        else:
            # Stage 2: nearest point on boundary within max_snap.
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
