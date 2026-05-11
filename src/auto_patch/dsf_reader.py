"""Read draped pavement polygons from a binary X-Plane DSF file.

DSF (Distribution Scenery File) is X-Plane's binary scenery format.
Beyond apt.dat row-110 polygons, scenery packs commonly add airport
pavement as DRAPED POLYGONS in the adjacent DSF — entries that
reference a ``.pol`` definition (typical paths under
``lib/airport/pavement/`` or ``lib/airport/ground/pavement/``).

We don't bother reimplementing DSF binary parsing; X-Plane ships
``DSFTool`` which converts DSF→text losslessly.  Ortho4XP bundles
``DSFTool`` for all three platforms under ``Utils/{lin,mac,win}/``.
This module:

  1. Locates the platform's ``DSFTool`` binary.
  2. Runs ``DSFTool --dsf2text`` on the requested DSF (caching the
     text output so repeated reads are free).
  3. Parses the text for ``POLYGON_DEF`` blocks and ``BEGIN_POLYGON``
     / ``POLYGON_POINT`` instances; filters to entries whose
     ``POLYGON_DEF`` path looks like pavement.
  4. Returns each pavement polygon as a list of (lon, lat) coords.

Returns are in lat/lon (EPSG:4326).  The caller is responsible for
projecting to its local meter coordinate system.
"""
from __future__ import annotations

import os
import platform
import subprocess
import sys
import tempfile
from typing import List, Optional, Tuple

import O4_File_Names as FNAMES
import O4_UI_Utils as UI


# Pavement-detector patterns: a POLYGON_DEF must START with one of
# these prefixes to be admitted as pavement geometry.  These are
# the X-Plane STOCK pavement library paths — bulk taxiway / apron /
# runway base pavement.  Third-party libraries (zannespol/, custom
# .pol files in scenery packs, etc.) are NOT trusted by this filter
# because their .pol entries are usually visual OVERLAYS painted on
# top of base pavement (e.g. ``LAYER_GROUP taxiways +1``,
# ``LAYER_GROUP runways +3``) — emitting them as pavement footprint
# duplicates apt.dat row-110 coverage and pulls non-pavement
# (grass-tinted, grunge, decorative) areas into the layout.
#
# Confirmed at SPJC where ``zannespol/Asphalt_2_Green_T80.pol`` has
# ``LAYER_GROUP taxiways +1`` and contributes 1.45 M m² of "pavement"
# that's actually a green-tinted overlay on the apron — driving the
# overlap + non-pavement coverage regression the user observed in
# JOSM.  Restricting to stock paths drops 1.86 M m² of decorative
# overlay; CYXY (which legitimately uses DSF as its sole pavement
# source) survives because every CYXY DSF pavement def is under
# ``lib/airport/pavement/`` or ``lib/airport/ground/pavement/``.
_PAVEMENT_PREFIXES = (
    "lib/airport/pavement/",
    "lib/airport/ground/pavement/",
)
# Skip patterns inside the accepted prefixes — line markings,
# direction signs, etc. that share path namespace with bulk
# pavement.
_PAVEMENT_SKIP = (
    "/lines/",
    "/markings/",
    "/lights/",
    "/decals/",
    "DirSigns",
    "shoulder",
)


def _dsftool_path() -> Optional[str]:
    """Return the platform's bundled DSFTool binary, or None."""
    # Mirror the layout O4_Mesh_Utils uses for Triangle4XP.
    base = FNAMES.Utils_dir
    sysname = platform.system().lower()
    if sysname == "darwin":
        cand = os.path.join(base, "mac", "DSFTool")
    elif sysname == "windows":
        cand = os.path.join(base, "win", "DSFTool.exe")
    else:
        cand = os.path.join(base, "lin", "DSFTool")
    if os.path.isfile(cand):
        return cand
    return None


def _is_pavement_def(path: str) -> bool:
    """True if the POLYGON_DEF path is bulk pavement geometry.

    Strict: only X-Plane stock pavement library paths are admitted.
    Third-party libraries (e.g. ``zannespol/``, ``CDB-Library/``,
    ``aericaps_collection/``) commonly use the same path tokens
    (``asphalt``, ``concrete``, ``tarmac``) for layered visual
    OVERLAYS rather than base pavement footprint, so admitting them
    by name pulls non-pavement decorative geometry into the layout.
    """
    p = path.lower()
    if not any(p.startswith(prefix) for prefix in _PAVEMENT_PREFIXES):
        return False
    if any(s in p for s in _PAVEMENT_SKIP):
        return False
    return True


def read_dsf_pavements(
    dsf_path: str,
    cache_dir: Optional[str] = None,
) -> List[List[Tuple[float, float]]]:
    """Extract draped pavement polygons from a DSF file.

    Args:
        dsf_path: path to a binary ``.dsf`` file.
        cache_dir: directory to store the converted text file
            (saves a re-run of DSFTool on subsequent reads).
            Defaults to a per-DSF temp file alongside the source.

    Returns:
        A list of pavement polygon rings.  Each ring is a list of
        ``(lon, lat)`` tuples (NOT closed — the first vertex isn't
        repeated).  Returns ``[]`` on any failure (DSFTool missing,
        DSF unreadable, no pavement defs, etc.).
    """
    if not dsf_path or not os.path.isfile(dsf_path):
        return []
    tool = _dsftool_path()
    if tool is None:
        UI.vprint(1,
            "  [dsf-reader] WARN: DSFTool binary not found at "
            f"{os.path.join(FNAMES.Utils_dir, platform.system().lower())}; "
            "DSF pavement will not be loaded.")
        return []

    # Cache the converted text alongside the DSF (or in cache_dir).
    if cache_dir is None:
        cache_dir = os.path.dirname(dsf_path)
    text_path = os.path.join(
        cache_dir,
        os.path.basename(dsf_path) + ".text",
    )
    # Re-convert if text is missing or older than the DSF.
    needs_convert = (not os.path.isfile(text_path)
                     or (os.path.getmtime(text_path)
                         < os.path.getmtime(dsf_path)))
    if needs_convert:
        try:
            # Some platforms don't allow writing into Custom Scenery;
            # fall back to a temp file in /tmp if the cache write fails.
            try:
                subprocess.run(
                    [tool, "--dsf2text", dsf_path, text_path],
                    check=True, capture_output=True, timeout=120,
                )
            except (PermissionError, subprocess.CalledProcessError):
                fallback = tempfile.NamedTemporaryFile(
                    suffix=".dsf.text", delete=False)
                text_path = fallback.name
                fallback.close()
                subprocess.run(
                    [tool, "--dsf2text", dsf_path, text_path],
                    check=True, capture_output=True, timeout=120,
                )
        except (OSError, subprocess.SubprocessError) as exc:
            UI.vprint(1,
                f"  [dsf-reader] WARN: DSFTool failed on "
                f"{os.path.basename(dsf_path)}: {exc}")
            return []

    try:
        with open(text_path, "r", encoding="utf-8",
                  errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return []

    # Pass 1: collect POLYGON_DEFs in order; track which indices
    # are pavement.
    pav_def_idx: set = set()
    def_idx = 0
    for line in lines:
        if line.startswith("POLYGON_DEF"):
            tok = line.strip().split(maxsplit=1)
            path = tok[1] if len(tok) > 1 else ""
            if _is_pavement_def(path):
                pav_def_idx.add(def_idx)
            def_idx += 1
    if not pav_def_idx:
        return []

    # Pass 2: walk BEGIN_POLYGON / END_POLYGON / BEGIN_WINDING /
    # END_WINDING / POLYGON_POINT to build per-instance rings.
    # A polygon may have multiple windings (outer + holes) — we
    # only emit the OUTER (first) winding here; holes are rare for
    # pavement and the caller's polygon-builder treats each ring
    # as its own outer.
    polys: List[List[Tuple[float, float]]] = []
    in_pavement = False
    in_winding = False
    current_outer: Optional[List[Tuple[float, float]]] = None
    pushed_this_polygon = False
    for line in lines:
        if line.startswith("BEGIN_POLYGON"):
            tok = line.split()
            try:
                idx = int(tok[1])
            except (ValueError, IndexError):
                idx = -1
            in_pavement = idx in pav_def_idx
            current_outer = None
            pushed_this_polygon = False
            continue
        if line.startswith("END_POLYGON"):
            in_pavement = False
            in_winding = False
            current_outer = None
            pushed_this_polygon = False
            continue
        if not in_pavement:
            continue
        if line.startswith("BEGIN_WINDING"):
            if pushed_this_polygon:
                # Subsequent winding = inner ring (hole); we ignore
                # it for pavement extraction.
                in_winding = False
            else:
                in_winding = True
                current_outer = []
            continue
        if line.startswith("END_WINDING"):
            if (in_winding and current_outer
                    and len(current_outer) >= 3):
                polys.append(current_outer)
                pushed_this_polygon = True
            in_winding = False
            current_outer = None
            continue
        if in_winding and line.startswith("POLYGON_POINT"):
            tok = line.split()
            try:
                lon = float(tok[1])
                lat = float(tok[2])
            except (ValueError, IndexError):
                continue
            current_outer.append((lon, lat))
    return polys


def find_associated_dsf(apt_dat_path: str,
                        apt_lat: float,
                        apt_lon: float) -> Optional[str]:
    """Locate the DSF file in the same scenery pack as ``apt_dat_path``
    that covers ``(apt_lat, apt_lon)``.

    DSFs are organized by tile (1° × 1°) inside an ``Earth nav data``
    tree.  The apt.dat sits next to a per-tile subtree like
    ``Earth nav data/+60-140/+60-136.dsf``.

    Returns the absolute DSF path or None.
    """
    if not apt_dat_path or not os.path.isfile(apt_dat_path):
        return None
    # apt.dat lives at <pack>/Earth nav data/apt.dat
    end_dir = os.path.dirname(apt_dat_path)
    if os.path.basename(end_dir) != "Earth nav data":
        return None
    # Tile + group dir naming.
    import math
    tile_lat = int(math.floor(apt_lat))
    tile_lon = int(math.floor(apt_lon))
    grp_lat = (tile_lat // 10) * 10
    grp_lon = (tile_lon // 10) * 10

    def _fmt(v: int, pad: int) -> str:
        sign = "+" if v >= 0 else "-"
        return f"{sign}{abs(v):0{pad}d}"

    grp = f"{_fmt(grp_lat, 2)}{_fmt(grp_lon, 3)}"
    tile = f"{_fmt(tile_lat, 2)}{_fmt(tile_lon, 3)}"
    cand = os.path.join(end_dir, grp, tile + ".dsf")
    if os.path.isfile(cand):
        return cand
    return None
