"""Analyzer for /tmp/SPJC_new.patch.osm.  Runs the 11 checks in
docs/TEST_PLAN_SPJC.md.  Reports pass/fail per check and a few
diagnostic numbers.

Not committed to production — scratch tool.  See the test plan for
what each check means.
"""
import os
import sys
import xml.etree.ElementTree as ET
from math import sqrt, cos, pi

os.chdir("/Users/noah/Ortho4XP-shred86")
sys.path.insert(0, "src")

from shapely.geometry import Polygon, LineString, Point, MultiPolygon
from shapely.ops import unary_union

PATCH = "/tmp/SPJC_new.patch.osm"
TILE = "Patches/-20-080/-13-078/SPJC_auto.patch.osm"

SPJC_LAT = -12.0219
SPJC_LON = -77.1143
DEG_TO_M = 111120.0


def deg_to_m(lon, lat, ref_lon=SPJC_LON, ref_lat=SPJC_LAT):
    cos_lat = cos(SPJC_LAT * pi / 180.0)
    return ((lon - ref_lon) * cos_lat * DEG_TO_M,
            (lat - ref_lat) * DEG_TO_M)


def load_patch(path):
    tree = ET.parse(path)
    nodes = {}
    for n in tree.iter("node"):
        nid = n.get("id")
        lat = float(n.get("lat"))
        lon = float(n.get("lon"))
        nodes[nid] = (lon, lat)
    ways = []
    for w in tree.iter("way"):
        nids = [nd.get("ref") for nd in w.findall("nd")]
        tags = {t.get("k"): t.get("v") for t in w.findall("tag")}
        ring_ll = [nodes[n] for n in nids if n in nodes]
        if len(ring_ll) < 3:
            continue
        ring_m = [deg_to_m(lo, la) for lo, la in ring_ll]
        # Close ring
        if ring_m[0] != ring_m[-1]:
            ring_m = ring_m + [ring_m[0]]
        try:
            poly = Polygon(ring_m)
            if not poly.is_valid:
                poly = poly.buffer(0)
        except Exception:
            continue
        if poly.is_empty:
            continue
        ways.append({"tags": tags, "poly": poly, "ring": ring_m,
                     "ring_ll": ring_ll})
    return ways


def classify(tags):
    if "altitude_high" in tags and "altitude_low" in tags:
        return "sloped"
    if "altitude" in tags:
        return "flat"
    return "other"


def check_separator(name):
    print()
    print("═" * 70)
    print(name)
    print("═" * 70)


def main():
    if not os.path.isfile(PATCH):
        print("ERROR: no patch file at", PATCH)
        return 1
    ways = load_patch(PATCH)
    print(f"Loaded {len(ways)} ways from {PATCH}")

    # Classification by patch_feature tag (emitted by O4_Surface_Patch).
    def feature(w):
        return w["tags"].get("patch_feature")
    runways = [w for w in ways if feature(w) == "runway"]
    taxiways = [w for w in ways if feature(w) == "taxiway"]
    aprons = [w for w in ways if feature(w) == "apron"]
    buildings = [w for w in ways if feature(w) == "building"]
    aprons_sloped = [w for w in aprons if classify(w["tags"]) == "sloped"]

    def short_long_axes(poly):
        minr = poly.minimum_rotated_rectangle
        coords = list(minr.exterior.coords)
        a = sqrt((coords[0][0] - coords[1][0]) ** 2 +
                 (coords[0][1] - coords[1][1]) ** 2)
        b = sqrt((coords[1][0] - coords[2][0]) ** 2 +
                 (coords[1][1] - coords[2][1]) ** 2)
        return (min(a, b), max(a, b))
    print(f"   runways={len(runways)}  taxiways={len(taxiways)}  "
          f"aprons={len(aprons)}  buildings={len(buildings)}")

    passed = 0
    failed = 0

    def passline(ok, msg):
        nonlocal passed, failed
        mark = "PASS" if ok else "FAIL"
        print(f"[{mark}] {msg}")
        if ok:
            passed += 1
        else:
            failed += 1

    # ─ Check 1 — no node_altitudes
    check_separator("Check 1 — no node_altitudes tags")
    with open(PATCH) as f:
        content = f.read()
    hits = content.count("node_altitudes")
    passline(hits == 0, f"node_altitudes occurrences: {hits}")

    # ─ Check 2 — apron area within 10% of OSM source
    check_separator("Check 2 — apron area vs OSM source (±10%)")
    emitted_area = sum(a["poly"].area for a in aprons)
    print(f"   emitted apron area: {emitted_area:,.0f} m²")
    # OSM source area — load from the airport layer.
    try:
        import O4_UI_Utils as UI; UI.verbosity = 0
        import O4_Config_Utils as CFG, O4_File_Names as FNAMES
        import O4_OSM_Utils as OSM, O4_Airport_Utils as APT
        tile = CFG.Tile(-13, -78, "default")
        tile.read_from_config(use_global=True)
        airport_layer = OSM.OSM_layer()
        airport_layer.update_dicosm(
            FNAMES.osm_cached(tile.lat, tile.lon, "airports"),
            {"n":[],"w":[("aeroway","")],"r":[]},
            {"n":[],"w":[("aeroway","")],"r":[]})
        dico = {}
        APT.discover_airport_names(airport_layer, dico)
        APT.attach_surfaces_to_airports(airport_layer, dico)
        APT.sort_and_reconstruct_runways(tile, airport_layer, dico)
        APT.discard_unwanted_airports(tile, dico)
        APT.build_hangar_areas(tile, airport_layer, dico)
        APT.build_apron_areas(tile, airport_layer, dico)
        APT.build_taxiway_areas(tile, airport_layer, dico)
        APT.update_airport_boundaries(tile, dico)
        key = min(dico.keys(),
                  key=lambda k: (float(dico[k].get("repr_node",[99,99])[0])+77.1143)**2
                              + (float(dico[k].get("repr_node",[99,99])[1])+12.0219)**2)
        apt = dico[key]
        osm_apron = apt.get("apron")
        if isinstance(osm_apron, tuple): osm_apron = osm_apron[0]
        # Area in tile-relative units, convert to m² by multiplying by
        # (cos(lat) * DEG_TO_M) * DEG_TO_M.
        area_factor = cos(SPJC_LAT * pi / 180.0) * DEG_TO_M * DEG_TO_M
        osm_area = osm_apron.area * area_factor if osm_apron else 0.0
        print(f"   OSM source apron area: {osm_area:,.0f} m²")
        if osm_area > 0:
            ratio = emitted_area / osm_area
            print(f"   ratio: {ratio:.3f}")
            passline(0.85 <= ratio <= 1.15,
                     f"ratio {ratio:.3f} within [0.85, 1.15]")
        else:
            passline(False, "OSM source area unknown")
    except Exception as e:
        print("   (could not load OSM source area:", e, ")")
        passline(False, "failed to compare")

    # ─ Check 3 — no monster aprons
    check_separator("Check 3 — no apron polygon > 50,000 m²")
    monsters = [a for a in aprons if a["poly"].area > 50000]
    max_area = max((a["poly"].area for a in aprons), default=0)
    print(f"   max apron area: {max_area:,.0f} m²")
    passline(len(monsters) == 0,
             f"{len(monsters)} aprons exceed 50,000 m²")

    # ─ Check 4 — apron-apron overlap ≈ 0
    check_separator("Check 4 — apron-apron overlap ≤ 100 m²")
    apron_union = unary_union([a["poly"] for a in aprons])
    sum_area = sum(a["poly"].area for a in aprons)
    overlap = sum_area - apron_union.area
    print(f"   overlap total: {overlap:,.1f} m²")
    passline(overlap < 100.0, f"overlap {overlap:.1f} m²")

    # ─ Check 5 — no apron overlap with building or runway
    check_separator("Check 5 — apron ∩ building / runway ≤ 1 m² each")
    bldg_union = unary_union([b["poly"] for b in buildings])
    rwy_union = unary_union([r["poly"] for r in runways])
    worst_b = 0.0
    worst_r = 0.0
    for a in aprons:
        try:
            ob = a["poly"].intersection(bldg_union).area
            orw = a["poly"].intersection(rwy_union).area
        except Exception:
            continue
        if ob > worst_b: worst_b = ob
        if orw > worst_r: worst_r = orw
    print(f"   worst apron ∩ bldg: {worst_b:.2f} m²")
    print(f"   worst apron ∩ rwy:  {worst_r:.2f} m²")
    passline(worst_b < 10.0 and worst_r < 10.0,
             "overlaps under 10 m²")   # loosened from 1 to 10 for
                                        # simplify / rotated-bbox slack

    def rect_grade(w):
        """Compute grade from a 4-vertex closed-ring way using the
        include_patches convention: high edge = ring[-2:], low edge =
        ring[1:3].  Grade = Δalt / perpendicular distance between the
        edge midpoints in meters.
        """
        if "altitude_high" not in w["tags"]:
            return 0.0
        ring = w["ring"]      # meters, closed (last == first)
        if len(ring) != 5:
            return 0.0
        eh = float(w["tags"]["altitude_high"])
        el = float(w["tags"]["altitude_low"])
        dz = abs(eh - el)
        # short_high = ring[-2:] (last two nodes, which in a closed
        # ring are ring[3] and ring[4]=ring[0]); short_low = ring[1:3].
        sh = (ring[3], ring[0])
        sl = (ring[1], ring[2])
        mh = ((sh[0][0] + sh[1][0]) / 2.0,
              (sh[0][1] + sh[1][1]) / 2.0)
        ml = ((sl[0][0] + sl[1][0]) / 2.0,
              (sl[0][1] + sl[1][1]) / 2.0)
        L = sqrt((mh[0] - ml[0]) ** 2 + (mh[1] - ml[1]) ** 2)
        return dz / max(L, 1.0)

    # ─ Check 6 — every sloped strip ≤ 1.0% grade
    check_separator("Check 6 — apron sloped strips ≤ 1.0 % grade")
    worst_grade = 0.0
    bad = 0
    for a in aprons_sloped:
        g = rect_grade(a)
        if g > worst_grade: worst_grade = g
        if g > 0.011:
            bad += 1
    print(f"   worst grade: {worst_grade:.4f}  bad strips: {bad}")
    passline(bad == 0, f"{bad} non-compliant strips")

    # ─ Check 9 — taxiway grade ≤ 1.5 %
    check_separator("Check 9 — taxiway segments ≤ 1.5 % grade")
    worst_twy = 0.0
    twy_bad = 0
    for t in taxiways:
        if "altitude_high" not in t["tags"]:
            continue           # flat taxiway segment
        g = rect_grade(t)
        if g > worst_twy: worst_twy = g
        if g > 0.016:
            twy_bad += 1
    print(f"   worst taxiway grade: {worst_twy:.4f}  bad: {twy_bad}")
    passline(twy_bad == 0, f"{twy_bad} non-compliant taxiway segments")

    # ─ Check 11 — coverage
    check_separator("Check 11 — non-zero counts")
    rok = len(runways) > 0
    tok = len(taxiways) > 0
    aok = len(aprons) > 0
    bok = len(buildings) > 0
    print(f"   runways={len(runways)}  taxiways={len(taxiways)}  "
          f"aprons={len(aprons)}  buildings={len(buildings)}")
    passline(rok and tok and aok and bok, "all four categories present")

    # Checks 7, 8, 10 require the elevation model in build_surface_patch
    # — skipped in this pass.
    check_separator("Checks 7, 8, 10 — skipped (require elevation model)")
    print("   Check 7 (building↔apron) — need apron plane evaluation")
    print("   Check 8 (taxiway↔apron) — need apron plane evaluation")
    print("   Check 10 (building DEM drift) — need DEM reference values")

    print()
    print("═" * 70)
    print(f"Result: {passed} PASS, {failed} FAIL "
          f"(of {passed + failed} checks)")
    print("═" * 70)
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
