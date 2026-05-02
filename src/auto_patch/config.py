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
    "AXIS_ALIGN_TOL_DEG",
    "LOAD_DSF_PAVEMENT",
    "LONG_EDGE_SNAP_M",
    "EMIT_JUNCTIONS",
    "EMIT_APRONS",
    "EMIT_BRIDGES_AND_TUNNELS",
    "JUNCTION_CLUSTER_DIST_M",
    "MAX_BOUNDARY_EDGE_M",
    "MIN_SEGMENT_LEN_M",
    "NECK_ABSOLUTE_M",
    "NECK_ABSORB_FRAC",
    "NECK_RELATIVE",
    "RUNWAY_ADJACENCY_TOL_M",
    "RUNWAY_BOUNDARY_TOL_M",
    "RUNWAY_INSIDE_APRON_FRAC",
    "RUNWAY_APRON_AREA_RATIO",
    "SLIVER_ANGLE_THRESHOLD_DEG",
]


# ── Junction-refinement rule constants (user 2026-05-01) ─────────
# Plan: ``/Users/noah/.claude/plans/kind-meandering-sifakis.md``.

# Rule 2: junction vertex within this distance of a sloping rect's
# long edge gets snapped to the rect's nearest short-end CORNER.
# Per user 2026-05-01: 10 m chosen because there should never be any
# junction node along a sloping rect's long side; 10 m is a safety
# band that catches densification midpoints + boundary-trace vertices
# that float off-rect due to apt.dat row-110 curvature noise.
LONG_EDGE_SNAP_M = 10.0

# Rule 4: a junction polygon is split at its narrowest cross-section
# when that thickness is below NECK_ABSOLUTE_M (in metres) OR is
# below NECK_RELATIVE × the polygon's MRR long-side length
# (whichever fires first per user 2026-05-01).
NECK_ABSOLUTE_M = 5.0
NECK_RELATIVE = 0.10
# A piece resulting from a neck split is absorbed into a neighbouring
# rect / junction / apron when the neighbour shares more than this
# fraction of the piece's perimeter.
NECK_ABSORB_FRAC = 0.70

# Rule 1: a junction vertex is "on the runway boundary" when within
# this distance of the runway polygon edge.  Used to identify the
# runway-adjacent vertex run that gets replaced with the runway's
# exact node sequence.
RUNWAY_BOUNDARY_TOL_M = 1.5

# Rule 1 v2 (user 2026-05-02): the WIDER tolerance for detecting
# vertices in a junction's runway-facing region.  Vertices in this
# band but outside RUNWAY_BOUNDARY_TOL_M still count as part of the
# junction's joining edge that needs to widen out to the next
# runway node.  Captures densification midpoints at ~2 m and
# boundary-trace vertices that float a few metres off the runway.
RUNWAY_ADJACENCY_TOL_M = 5.0

# Rule 3 test tolerance: a non-pavement, non-anchor junction edge
# must run parallel or perpendicular to the longest runway axis
# within this many degrees.
AXIS_ALIGN_TOL_DEG = 2.0


# Max length of any junction-polygon ring segment.  Long edges get
# subdivided to anchor Triangle4XP's interior triangulation; without
# this the elevation solver leaves the interior un-anchored on long
# straight runs and produces visible cliffs.  Shared between the
# junction-decomposition pass (densification) and the elevation
# layer (vertex-aware grade clamp).
MAX_BOUNDARY_EDGE_M = 30.0


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
