"""Run the auto-patch pavement builder on an ICAO and dump to OSM.

Usage:
    python3 tools/build_target_osm.py <ICAO> [--xplane PATH] [--out PATH]

Then compare with ``tools/compare_target.py`` against the reference
target OSM.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from auto_patch.pipeline import build_airport_pavement


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("icao")
    ap.add_argument("--xplane", default="/Users/noah/X-Plane 12")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    out = args.out or f"/tmp/{args.icao}_auto.osm"
    layout = build_airport_pavement(args.icao, args.xplane)
    layout.to_osm(out)

    # Quick summary
    from collections import Counter
    roles = Counter(s.role for s in layout.shapes)
    print(f"Wrote {out}")
    print(f"  anchor={layout.anchor}")
    print(f"  shapes={len(layout.shapes)}")
    for r in sorted(roles):
        print(f"    {r:<22} {roles[r]}")


if __name__ == "__main__":
    main()
