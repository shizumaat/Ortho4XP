"""PavementLayout + BuiltShape data model and serialisation.

Holds the pavement-builder's core data types (the dataclasses every
phase reads + writes), the role-tag vocabulary, the meter-anchored
projection helpers used to construct a layout, and the .osm
serialisation method (``PavementLayout.to_osm``).

Phase-1 (geometry) and Phase-2 (elevation) both populate
``BuiltShape`` instances inside a ``PavementLayout``; the .osm
emission turns the meter-space layout back into JOSM-readable
WGS-84 OSM with shared node IDs.

Public API:
    BuiltShape, PavementLayout              — data classes
    ROLE_*                                  — role-tag constants
    AEROWAY_FOR_ROLE                        — role -> aeroway tag value
    SHARED_VERTEX_TOL_M, R_EARTH            — geometry constants
    airport_anchor(apt), projection(anchor) — meter-space helpers

Used by every O4_Pavement_* module.  Sits at the bottom of the
pavement-builder dependency hierarchy alongside Pavement_Config.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from shapely.geometry import LineString, Polygon

from . import apt_dat_reader as APR
from .pavement import strips as PS

from .config import SLIVER_ANGLE_THRESHOLD_DEG

__all__ = [
    "BuiltShape",
    "PavementLayout",
    "R_EARTH",
    "SHARED_VERTEX_TOL_M",
    "ROLE_RUNWAY",
    "ROLE_PRIMARY_PARALLEL",
    "ROLE_SECONDARY_PARALLEL",
    "ROLE_STUB",
    "ROLE_CROSS_CONNECTOR",
    "ROLE_APRON",
    "ROLE_TERMINAL",
    "ROLE_JUNCTION",
    "ROLE_BOUNDARY",
    "ROLE_TUNNEL_RAMP",
    "ROLE_RETAINING_WALL",
    "ROLE_GROUNDSIDE_PAVEMENT",
    "AEROWAY_FOR_ROLE",
    "_airport_anchor",
    "_projection",
]


# ──────────────────────────────────────────────────────────────────
# Geometry constants
# ──────────────────────────────────────────────────────────────────
R_EARTH = 6_378_137.0
SHARED_VERTEX_TOL_M = 0.5    # snap vertices closer than this together


# ──────────────────────────────────────────────────────────────────
# Role-tag vocabulary
# ──────────────────────────────────────────────────────────────────
ROLE_RUNWAY = "runway"
ROLE_PRIMARY_PARALLEL = PS.ROLE_PRIMARY_PARALLEL
ROLE_SECONDARY_PARALLEL = PS.ROLE_SECONDARY_PARALLEL
ROLE_STUB = PS.ROLE_STUB
ROLE_CROSS_CONNECTOR = PS.ROLE_CROSS_CONNECTOR
ROLE_APRON = PS.ROLE_APRON
ROLE_TERMINAL = "terminal"
ROLE_JUNCTION = "junction"
ROLE_BOUNDARY = "boundary"
# Tunnel portals: ``tunnel_ramp`` is a sloped 4-corner rect from
# outside-DEM down to apt-elev-6m at the portal; ``retaining_wall``
# is a flat polygon at apt-elev forming the U-shape around the
# portal LOW end.
ROLE_TUNNEL_RAMP = "tunnel_ramp"
ROLE_RETAINING_WALL = "retaining_wall"
# Groundside terminal pavement (curbside / drop-off / parking) —
# emitted with per-vertex DEM altitudes and a 0.1 m gap from the
# terminal building so it follows local terrain instead of being
# flattened to airside-apron elevation.
ROLE_GROUNDSIDE_PAVEMENT = "groundside_pavement"

AEROWAY_FOR_ROLE = {
    ROLE_RUNWAY: "runway",
    ROLE_PRIMARY_PARALLEL: "taxiway",
    ROLE_SECONDARY_PARALLEL: "taxiway",
    ROLE_STUB: "taxiway",
    ROLE_CROSS_CONNECTOR: "taxiway",
    ROLE_JUNCTION: "taxiway",
    ROLE_APRON: "apron",
    ROLE_TERMINAL: "building",
    ROLE_BOUNDARY: "aerodrome",
    ROLE_TUNNEL_RAMP: "taxiway",
    ROLE_RETAINING_WALL: "building",
    ROLE_GROUNDSIDE_PAVEMENT: "apron",
}


# ──────────────────────────────────────────────────────────────────
# Data model
# ──────────────────────────────────────────────────────────────────

@dataclass
class BuiltShape:
    """A single emitted shape: polygon + classification tags.

    Polygons live in meter space anchored at the layout's origin.
    ``ref`` is optional (runway designator, taxi ref from OSM, or
    generated label).  ``source_axis`` is kept on taxi rects for
    elevation sampling along their axis.

    Phase-2 elevation: exactly one of these options is set at any
    time:
      * all four None (no elevation yet)
      * only ``altitude`` set (flat polygon at that elevation, m)
      * ``altitude_high`` + ``altitude_low`` set (linearly sloped
        between the two parallel edges).  Rects use this.
      * ``node_altitudes`` set (per-vertex elevation list, one
        value per ring vertex INCLUDING the closing repeat —
        i.e. len(node_altitudes) == len(closed_ring_nids)).
        Used for triangulated junction polygons that slope in
        more than one direction.
    Runway segments carry altitude_high/low per the legacy patch
    convention.
    """
    polygon: Polygon
    role: str
    ref: str = ""
    source_axis: Optional[LineString] = None
    altitude: Optional[float] = None
    altitude_high: Optional[float] = None
    altitude_low: Optional[float] = None
    node_altitudes: Optional[List[float]] = None
    # OSM ``bridge=yes`` flag.  Set on taxi rects whose source
    # OSM way is tagged as a bridge — see ``_emit_taxi_bridges``.
    is_bridge: bool = False



@dataclass
class PavementLayout:
    icao: str
    anchor: Tuple[float, float]          # (lat0, lon0)
    shapes: List[BuiltShape] = field(default_factory=list)
    # ancillary:
    airport_boundary: Optional[Polygon] = None
    runway_union: Optional[Polygon] = None
    # Path to the apt.dat file the layout was built from.  Used by
    # the bridge-detection step to walk the same scenery pack's
    # DSF and check for taxi-bridge OBJ placements.
    apt_dat_path: Optional[str] = None

    # ---- coordinate helpers ------------------------------------------
    def m_to_ll(self, x: float, y: float) -> Tuple[float, float]:
        lat0, lon0 = self.anchor
        cos0 = math.cos(math.radians(lat0))
        lon = lon0 + math.degrees(x / (R_EARTH * cos0))
        lat = lat0 + math.degrees(y / R_EARTH)
        return lat, lon

    def ll_to_m(self, lat: float, lon: float) -> Tuple[float, float]:
        lat0, lon0 = self.anchor
        cos0 = math.cos(math.radians(lat0))
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return x, y

    # ---- serialization -----------------------------------------------
    def to_osm(self, path: str) -> None:
        """Emit to a JOSM-readable OSM file with shared node IDs.

        Vertices within ``SHARED_VERTEX_TOL_M`` are assigned the same
        node id, matching the target-OSM convention.  Polygons with
        interior rings (holes, typical for junction polygons that
        wrap around rects) are emitted as OSM multipolygon
        relations — each ring becomes a closed way, and a relation
        ties them together with role=outer / role=inner.
        """
        bucket_size = SHARED_VERTEX_TOL_M
        node_key_to_id: Dict[Tuple[int, int], int] = {}
        node_id_to_ll: Dict[int, Tuple[float, float]] = {}
        next_nid = [-1]

        def _intern(x: float, y: float) -> int:
            kx = int(round(x / bucket_size))
            ky = int(round(y / bucket_size))
            key = (kx, ky)
            if key in node_key_to_id:
                return node_key_to_id[key]
            nid = next_nid[0]
            next_nid[0] -= 1
            node_key_to_id[key] = nid
            # Use the ACTUAL coordinates of the first vertex that
            # landed in this bucket (not the bucket center).  Bucket
            # centers can lie up to bucket_size/2 away from the real
            # vertex, which introduces sub-metre overlaps between
            # adjacent shapes after OSM round-trip.
            node_id_to_ll[nid] = self.m_to_ll(x, y)
            return nid

        def _ring_to_nids(ring_coords, ring_elevs=None):
            """Build a closed-ring nid list from coords.

            Returns ``(nids, elevs_or_None)``.  ``elevs_or_None`` is
            an aligned per-vertex elevation list when ``ring_elevs``
            is provided (used for ``node_altitudes`` polygons);
            otherwise None.  Both share the same dedup logic so the
            element-count invariant survives.

            Defensive against upstream polygon-build bugs:

            * Drops consecutive duplicate nids (two ring vertices
              colliding in the SHARED_VERTEX_TOL_M bucket — would
              produce a zero-length edge that crashes downstream
              meshers).
            * Drops non-consecutive duplicate nids (a ring revisits
              the same node — figure-8 / self-touching polygon —
              keeping only the first occurrence).
            """
            coords = list(ring_coords)
            elevs = list(ring_elevs) if ring_elevs is not None else None
            if coords and coords[0] == coords[-1]:
                coords = coords[:-1]
                if elevs is not None and len(elevs) > 1 and elevs[0] == elevs[-1]:
                    elevs = elevs[:-1]
            if len(coords) < 3:
                return None, None
            nids = [_intern(x, y) for (x, y) in coords]
            # Dedup any duplicate nid (consecutive OR not).
            seen: set = set()
            deduped_nids: List[int] = []
            deduped_elevs: List[float] = []
            for k, nid in enumerate(nids):
                if nid in seen:
                    continue
                seen.add(nid)
                deduped_nids.append(nid)
                if elevs is not None and k < len(elevs):
                    deduped_elevs.append(elevs[k])
            if len(deduped_nids) < 3:
                return None, None
            deduped_nids.append(deduped_nids[0])
            if elevs is not None:
                deduped_elevs.append(deduped_elevs[0])
                return deduped_nids, deduped_elevs
            return deduped_nids, None

        # Emit one simple way per shape (exterior ring only, with
        # all tags on that way).  Interior rings — which appear
        # on junction polygons that wrap around rect-shaped holes
        # — are dropped for X-Plane patch compatibility: the
        # Ortho4XP patch parser ([O4_Vector_Map.include_patches])
        # iterates ways only, so tags on an OSM multipolygon
        # relation never reach the outer way.  A junction ring
        # emitted without its holes will slightly overlap the
        # rects that used to punch those holes — X-Plane
        # triangulator handles the overlap by seed-region
        # processing; the rects' altitude_high/low tags prevail
        # where they cover.
        way_blocks: List[Tuple[int, List[int], Dict[str, str]]] = []
        next_wid = [-10001]
        for s in self.shapes:
            # Validate the polygon's geometry before emission.
            # Upstream pipeline stages (decomposition, seam-point
            # injection, shared-vertex enforcement) can occasionally
            # produce a self-touching ring that's geometrically
            # invalid; X-Plane's mesh builder crashes on these.
            poly = s.polygon
            if poly is None or poly.is_empty:
                continue
            if not poly.is_valid:
                try:
                    repaired = poly.buffer(0)
                    if (repaired.is_empty
                            or repaired.geom_type
                            not in ("Polygon", "MultiPolygon")):
                        continue
                    if repaired.geom_type == "MultiPolygon":
                        repaired = max(repaired.geoms,
                                       key=lambda g: g.area)
                    if (repaired.is_empty
                            or repaired.geom_type != "Polygon"):
                        continue
                    # node_altitudes from the original ring no longer
                    # aligns with the repaired ring; degrade to a
                    # flat polygon at the mean of the original
                    # vertex elevations to preserve emission.
                    if s.node_altitudes:
                        valid_elevs = [
                            e for e in s.node_altitudes[:-1]]
                        if valid_elevs:
                            s.altitude = round(
                                sum(valid_elevs) / len(valid_elevs),
                                1)
                        s.node_altitudes = None
                    poly = repaired
                except Exception:
                    continue
            # Pass node_altitudes alongside ring coords so dedup of
            # duplicate nids drops the matching elevations too,
            # keeping the per-vertex count invariant.
            ext_nids, ext_elevs = _ring_to_nids(
                poly.exterior.coords,
                s.node_altitudes)
            if ext_nids is None:
                continue
            # Final validity check: rebuild the polygon from the
            # POST-DEDUP lat/lon coords AT THE PRECISION THE OSM
            # FILE WILL CONTAIN (.11f, ≈ 1 mm at the equator).
            # Polygons that are valid at full float precision can
            # become spike-vertex-on-non-adjacent-edge invalid
            # after this truncation; X-Plane's mesh builder
            # crashes on those.  Drop the whole shape rather than
            # ship a polygon X-Plane can't handle.
            try:
                latlon_ring = []
                for nid in ext_nids[:-1]:
                    lat, lon = node_id_to_ll[nid]
                    latlon_ring.append(
                        (float(f"{lat:.11f}"),
                         float(f"{lon:.11f}")))
                check_poly = Polygon(
                    [(lon, lat) for lat, lon in latlon_ring])
                if not check_poly.is_valid:
                    try:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] WARN: dropping "
                            f"invalid polygon (role={s.role}, "
                            f"nids={len(ext_nids) - 1}): "
                            f"X-Plane mesh builder would crash.\n")
                    except Exception:
                        pass
                    continue
                # Sliver-corner safety net: if any interior angle is
                # below SLIVER_ANGLE_THRESHOLD_DEG, drop the polygon.
                # The source-fix _drop_sliver_corners normally
                # eliminates these for junctions, but this catches
                # any path (rect emission, decomposition fragment,
                # repaired-by-buffer(0) polygon) that could
                # reintroduce a needle tip Triangle4XP can't handle.
                ring_m = [self.ll_to_m(lat, lon)
                          for (lat, lon) in latlon_ring]
                cos_thresh = math.cos(
                    math.radians(SLIVER_ANGLE_THRESHOLD_DEG))
                m = len(ring_m)
                worst_ang = None
                for vi in range(m):
                    ax, ay = ring_m[(vi - 1) % m]
                    bx, by = ring_m[vi]
                    cx, cy = ring_m[(vi + 1) % m]
                    v1x, v1y = ax - bx, ay - by
                    v2x, v2y = cx - bx, cy - by
                    n1 = math.hypot(v1x, v1y)
                    n2 = math.hypot(v2x, v2y)
                    if n1 < 1e-9 or n2 < 1e-9:
                        continue
                    cos = (v1x * v2x + v1y * v2y) / (n1 * n2)
                    if cos > cos_thresh:
                        worst_ang = math.degrees(
                            math.acos(max(-1.0, min(1.0, cos))))
                        break
                if worst_ang is not None:
                    try:
                        import sys as _sys
                        _sys.stderr.write(
                            f"  [pav-builder] WARN: dropping "
                            f"sliver-corner polygon (role={s.role}, "
                            f"nids={len(ext_nids) - 1}, "
                            f"min angle {worst_ang:.2f}°): "
                            f"X-Plane mesh builder would crash.\n")
                    except Exception:
                        pass
                    continue
            except Exception:
                continue
            tags = {
                "aeroway": AEROWAY_FOR_ROLE.get(s.role, "taxiway"),
                "role": s.role,
            }
            if s.ref:
                tags["ref"] = s.ref
            # Phase-2 elevation tags.  Sloped rects also carry
            # cell_size + profile so X-Plane uses spline
            # interpolation between the high and low short
            # edges, matching the legacy auto-patch format.
            if s.altitude_high is not None and s.altitude_low is not None:
                tags["altitude_high"] = f"{s.altitude_high:.1f}"
                tags["altitude_low"] = f"{s.altitude_low:.1f}"
                tags["cell_size"] = "2"
                tags["profile"] = "spline"
            elif (ext_elevs is not None
                  and len(ext_elevs) == len(ext_nids)):
                # Per-vertex elevation: comma-separated list, one
                # value per ring nid (including the closing repeat).
                # X-Plane mesh builder triangulates the polygon and
                # interpolates linearly between vertex elevations.
                tags["node_altitudes"] = ",".join(
                    f"{e:.1f}" for e in ext_elevs)
            elif s.altitude is not None:
                tags["altitude"] = f"{s.altitude:.1f}"
            way_blocks.append((next_wid[0], ext_nids, tags))
            next_wid[0] -= 1
        rel_blocks: List[Tuple[int, List[Tuple[int, str]],
                               Dict[str, str]]] = []

        # Determine which interned nodes are actually referenced by
        # any emitted way (via ``way_blocks`` or ``rel_blocks``
        # member ways).  Per user 2026-04-29: discarded ring builds
        # — short rings that ``_ring_to_nids`` returned None for,
        # or rings whose nodes were dedup'd out of the final ring —
        # leave orphan entries in ``node_id_to_ll``.  Emitting
        # those produces "floating nodes" next to a polygon in
        # JOSM that aren't part of any geometry.  Filter to
        # referenced nids only.
        referenced_nids: set = set()
        for _wid, _nids, _tags in way_blocks:
            referenced_nids.update(_nids)
        # rel_blocks is currently empty in this emitter but be
        # forward-compatible if multipolygons return.
        for _rid, _members, _tags in rel_blocks:
            for _mwid, _role in _members:
                # Member ways' nids — find them in way_blocks.
                for w_id, n_list, _t in way_blocks:
                    if w_id == _mwid:
                        referenced_nids.update(n_list)
                        break
        lines = [
            "<?xml version='1.0' encoding='UTF-8'?>",
            "<osm version='0.6' upload='false' generator='O4_Airport_Pavement_Builder'>",
        ]
        for nid, (lat, lon) in sorted(node_id_to_ll.items(), reverse=True):
            if nid not in referenced_nids:
                continue
            lines.append(
                f"  <node id='{nid}' action='modify' visible='true' "
                f"lat='{lat:.11f}' lon='{lon:.11f}' />"
            )
        for wid, nids, tags in way_blocks:
            lines.append(
                f"  <way id='{wid}' action='modify' visible='true'>"
            )
            for nid in nids:
                lines.append(f"    <nd ref='{nid}' />")
            for k, v in sorted(tags.items()):
                lines.append(f"    <tag k='{k}' v='{v}' />")
            lines.append("  </way>")
        for rid, members, tags in rel_blocks:
            lines.append(
                f"  <relation id='{rid}' action='modify' visible='true'>"
            )
            for mwid, role in members:
                lines.append(
                    f"    <member type='way' ref='{mwid}' role='{role}' />"
                )
            for k, v in sorted(tags.items()):
                lines.append(f"    <tag k='{k}' v='{v}' />")
            lines.append("  </relation>")
        lines.append("</osm>")
        Path(path).write_text("\n".join(lines) + "\n")


# ──────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────
# Projection helpers
# ──────────────────────────────────────────────────────────────────

def _projection(anchor: Tuple[float, float]):
    lat0, lon0 = anchor
    cos0 = math.cos(math.radians(lat0))

    def to_m(lon: float, lat: float, z=None):
        x = math.radians(lon - lon0) * R_EARTH * cos0
        y = math.radians(lat - lat0) * R_EARTH
        return (x, y) if z is None else (x, y, z)

    return to_m


def _airport_anchor(apt: APR.Airport) -> Tuple[float, float]:
    if apt.runways:
        r = apt.runways[0]
        return ((r.lat_a + r.lat_b) / 2.0,
                (r.lon_a + r.lon_b) / 2.0)
    if apt.boundary:
        c = apt.boundary.centroid
        return (c.y, c.x)
    return (0.0, 0.0)


# ──────────────────────────────────────────────────────────────────
# Runway rects
# ──────────────────────────────────────────────────────────────────
