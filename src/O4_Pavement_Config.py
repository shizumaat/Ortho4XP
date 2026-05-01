"""Pavement-builder configuration constants.

Single source of truth for every numeric tunable in the airport
pavement builder.  Module-local constants in other O4_Pavement_*
modules should be reserved for values whose meaning is genuinely
specific to that module; anything tuned across the pipeline lives
here so reviewers can audit the whole tuning surface in one place.

For user-facing knobs (GUI / cfg-file persistence) register the
variable in O4_Cfg_Vars.py instead.
"""

__all__ = [
    "RUNWAY_INSIDE_APRON_FRAC",
    "RUNWAY_APRON_AREA_RATIO",
    "LOAD_DSF_PAVEMENT",
    "EMIT_JUNCTIONS",
    "EMIT_APRONS",
    "EMIT_BRIDGES_AND_TUNNELS",
    "JUNCTION_CLUSTER_DIST_M",
    "MIN_SEGMENT_LEN_M",
    "SLIVER_ANGLE_THRESHOLD_DEG",
]


# Drop emitted line/segment fragments shorter than this length.
# Shared across centerline extraction, taxi-rect splitting, and the
# Phase-A apt.dat-rect chain construction.
MIN_SEGMENT_LEN_M = 15.0


# Cluster of junction-corner candidates: any two within this
# distance get merged when computing the residue's seam points.
JUNCTION_CLUSTER_DIST_M = 40.0

# Interior angles below this threshold count as "needle-tip"
# slivers.  Residue construction can leave thin wedges where
# rect / terminal edges meet the apt.dat boundary at near-collinear
# angles.  The polygon is shapely-valid but a sub-2 deg corner
# forces Triangle4XP to emit a near-degenerate triangle there --
# crashes X-Plane's mesh builder.  Caught at junction-emission
# time by _drop_sliver_corners (drops just the tip vertex), and
# again by a to_osm safety net (drops the whole shape if any
# slipped through).
SLIVER_ANGLE_THRESHOLD_DEG = 2.0


# Apron-merged runway detection: when at least this fraction of a
# runway-segment polygon lies inside an apt.dat / DSF apron polygon
# (and that polygon is much larger than the segment — see
# RUNWAY_APRON_AREA_RATIO), the segment is treated as apron-merged
# and the separate rect is dropped.
RUNWAY_INSIDE_APRON_FRAC = 0.95
# The containing apt.dat / DSF polygon must be >= this ratio times
# the segment area to count as an apron.
RUNWAY_APRON_AREA_RATIO = 3.0

# Per user 2026-04-30: bridge / tunnel emission disabled while we
# stabilise the core pavement-geometry pipeline.  Recent KPHX work
# emitted ``tunnel_ramp`` polygons that overlapped junction /
# groundside_pavement geometry (130 overlap pairs, 75K m^2 total at
# KPHX, caught by ``test_no_self_overlap``).
#
# TODO(bridges): once core geometry is stable + refactored, revisit
# the four feature emit calls gated by this flag in
# build_airport_pavement -- _emit_through_airport_depressed_roads,
# _emit_tunnel_portals, _emit_taxi_bridges,
# _emit_underpass_road_approaches.  Each must carve its footprint
# out of overlapping airside / groundside pavement before emitting,
# or test_no_self_overlap will fail again.
EMIT_BRIDGES_AND_TUNNELS = False

# Combine apt.dat with DSF pavement polygons: when True the
# smart-apt.dat selector still runs to choose the best custom-pack
# vs global candidate by OSM coverage; DSF polygons supplement
# whichever apt.dat is picked.
LOAD_DSF_PAVEMENT = True

# Phase-1 emit-suppression toggles (kept from the pre-refactor
# baseline; iteration aids that remain useful).
EMIT_JUNCTIONS = True
EMIT_APRONS = False
