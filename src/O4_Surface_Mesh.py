"""Per-airport surface patch triangulation helpers.

This module is the *patch-mesh* counterpart to ``O4_Mesh_Utils``.  The
two are easy to confuse but solve different problems:

============================  =========================================
Module                        Builds
============================  =========================================
``O4_Mesh_Utils``             The whole-tile terrain mesh.  Wraps the
                              external ``Triangle4XP`` binary that
                              produces the ``.mesh`` file consumed by
                              X-Plane's DSF builder.
``O4_Surface_Mesh`` *(this)*  The per-airport triangulations that go
                              inside ``{ICAO}_auto.patch.osm`` files.
                              Pure Python (numpy + shapely).
============================  =========================================

The main public function is :func:`adaptive_triangulate`.  Given an
input polygon, a set of "anchor" elevations from neighbouring features
(building edges, taxiway crossings, runway-touch points), and a DEM
sampling callback, it returns the **smallest** set of triangles that:

* exactly match the polygon's outline (using its real exterior
  vertices, not a densified grid),
* preserve interior DEM detail to within ``fidelity_tol`` metres,
* satisfy a per-triangle grade limit (1.0 % aprons, 1.5 % junctions
  and taxiways),
* keep their elevations consistent with neighbouring features at the
  anchored vertices.

This replaces the legacy approach in ``O4_Auto_Patch.generate_airport_surface_patches``,
which densified each polygon into a 30 m grid and Delaunay-triangulated
from there — producing 50+ triangles per apron whether the apron
needed them or not.  An adaptively-refined apron typically yields 2-15
triangles, dramatically reducing the emitted shape count and making
each shape much more meaningful to a human reader (and to JOSM).

The algorithm is:

1. Seed the triangulation with the polygon's exterior vertices, plus
   any anchor points that fall in the interior.  Anchored vertices
   take their elevation from the nearest anchor; non-anchored vertices
   sample DEM.
2. Triangulate (Delaunay) the seed point set.
3. **Initial grade clamp** — pull non-anchored vertices toward
   triangle centroids until every triangle's plane gradient is
   ≤ ``max_grade``.  Anchored vertices never move.
4. **Fidelity refinement loop** — for each triangle, find the DEM
   sample (centroid + edge midpoints) whose deviation from the
   current plane is largest.  If that deviation exceeds
   ``fidelity_tol`` AND the candidate point is grade-compatible with
   its surrounding triangle corners, insert it as a new vertex,
   re-triangulate, and re-run the grade clamp.  The
   grade-compatibility check is what keeps the loop from cascading
   when DEM disagrees with anchors.
5. **Final grade clamp** to converge after the last insertions.
6. Return the triangle list.
"""

from math import sqrt
from typing import Callable, Iterable, List, Optional, Sequence, Tuple

from shapely import geometry as shp_geom
from shapely import ops as shp_ops


# Type aliases for readability.
Vertex3D = Tuple[float, float, float]                # (x, y, z)
Triangle = Tuple[Vertex3D, Vertex3D, Vertex3D]
Anchor = Tuple[float, float, float]                  # (x, y, z)
DemSampler = Callable[[float, float], Optional[float]]


# Tolerances and limits — tuned for SPJC-scale airports.  Override per
# call via the function parameters when needed.
EPS_DUP_M = 0.5             # Vertex deduplication tolerance (m).
SAMPLE_GAP_M = 2.0          # Minimum spacing between vertices when
                            # inserting a new fidelity sample.
GRADE_CLAMP_PASSES = 12     # Iterations of the grade clamp loop.
DEFAULT_FIDELITY_TOL_M = 1.0
DEFAULT_MAX_EXTRA_POINTS = 50


# ──────────────────────────────────────────────────────────────────────
# Plane geometry helper
# ──────────────────────────────────────────────────────────────────────
def fit_plane(p0: Vertex3D, p1: Vertex3D,
              p2: Vertex3D) -> Optional[Tuple[float, float, float]]:
    """Return ``(a, b, c)`` such that ``z = a*x + b*y + c`` passes
    through the three input points exactly, or ``None`` if the
    points are vertically degenerate (the cross product's z component
    is near zero — i.e. the three points are colinear in (x, y) space).
    """
    x0, y0, z0 = p0
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    v1x, v1y, v1z = x1 - x0, y1 - y0, z1 - z0
    v2x, v2y, v2z = x2 - x0, y2 - y0, z2 - z0
    nx = v1y * v2z - v1z * v2y
    ny = v1z * v2x - v1x * v2z
    nz = v1x * v2y - v1y * v2x
    if abs(nz) < 1e-9:
        return None
    a = -nx / nz
    b = -ny / nz
    c = z0 - a * x0 - b * y0
    return (a, b, c)


def plane_grade(plane: Tuple[float, float, float]) -> float:
    """Magnitude of the plane's gradient (i.e. the slope as a
    fraction).  A plane ``z = 0.01*x + 0*y + 5`` has grade 0.01.
    """
    a, b, _c = plane
    return sqrt(a * a + b * b)


# ──────────────────────────────────────────────────────────────────────
# Internal: vertex bag with anchor flags
# ──────────────────────────────────────────────────────────────────────
class _VertexBag:
    """A small mutable list of vertices used during triangulation.

    Each entry is ``[x, y, z, anchored]``.  Anchored vertices are
    pinned by neighbour features (building edges, taxiway/runway
    touch points) and must not be moved during grade enforcement.
    """

    def __init__(self, eps: float = EPS_DUP_M):
        self._eps = eps
        self._eps2 = eps * eps
        self._pts: List[List] = []
        # Spatial-hash grid: cell_size >= eps so two points within
        # eps are guaranteed to be in the same or an adjacent cell.
        # cell_size is also tuned to the dedup + sample gap radii
        # (0.5 m and 2 m respectively): making it larger than 2 m
        # means has_near can't finish in O(1) cells, but making it
        # exactly eps means has_near with gap = SAMPLE_GAP_M must
        # scan (2*gap/cell+1)² ≈ 81 cells — still cheap.  The grid
        # replaces the O(n²) linear scan behaviour with O(n) total
        # builds + O(1) average lookups.
        self._cell_size = max(eps, 1.0)
        self._grid: dict = {}

    def _cell(self, x: float, y: float) -> tuple:
        return (int(x / self._cell_size), int(y / self._cell_size))

    def __iter__(self):
        return iter(self._pts)

    def __len__(self):
        return len(self._pts)

    def __getitem__(self, i):
        return self._pts[i]

    def add(self, x: float, y: float, z: float, anchored: bool) -> int:
        """Add a vertex.  If a vertex within EPS_DUP_M already exists,
        return its index instead.  An incoming anchored vertex
        promotes an existing un-anchored duplicate to anchored, but
        an existing anchor is never overwritten.
        """
        # Check the 9 neighbouring grid cells for a duplicate.
        cx, cy = self._cell(x, y)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                bucket = self._grid.get((cx + dx, cy + dy))
                if bucket is None:
                    continue
                for i in bucket:
                    p = self._pts[i]
                    if (p[0] - x) ** 2 + (p[1] - y) ** 2 < self._eps2:
                        if anchored and not p[3]:
                            self._pts[i] = [p[0], p[1], z, True]
                        return i
        idx = len(self._pts)
        self._pts.append([x, y, float(z), anchored])
        self._grid.setdefault((cx, cy), []).append(idx)
        return idx

    def index_at(self, x: float, y: float) -> Optional[int]:
        """Return the index of the vertex closest to (x, y) within
        EPS_DUP_M, or None if no vertex is that close.
        """
        cx, cy = self._cell(x, y)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                bucket = self._grid.get((cx + dx, cy + dy))
                if bucket is None:
                    continue
                for i in bucket:
                    p = self._pts[i]
                    if (p[0] - x) ** 2 + (p[1] - y) ** 2 < self._eps2:
                        return i
        return None

    def has_near(self, x: float, y: float, gap: float) -> bool:
        """True if any existing vertex is within `gap` of (x, y)."""
        gap2 = gap * gap
        cx, cy = self._cell(x, y)
        # Number of grid cells that could contain a point within
        # `gap` of (x, y).  For gap <= cell_size this is just the
        # 9 surrounding cells; wider `gap` walks more cells.
        r = int(gap / self._cell_size) + 1
        for dx in range(-r, r + 1):
            for dy in range(-r, r + 1):
                bucket = self._grid.get((cx + dx, cy + dy))
                if bucket is None:
                    continue
                for i in bucket:
                    p = self._pts[i]
                    if (p[0] - x) ** 2 + (p[1] - y) ** 2 < gap2:
                        return True
        return False


# ──────────────────────────────────────────────────────────────────────
# Internal: seeding and triangulation
# ──────────────────────────────────────────────────────────────────────
def _find_anchor_at(x: float, y: float, anchors: Sequence[Anchor],
                    radius: float) -> Optional[float]:
    """Return the elevation of the nearest anchor within `radius` of
    (x, y), or None.
    """
    best = None
    best_d2 = (radius * 2.0) ** 2
    for ax, ay, az in anchors:
        d2 = (ax - x) ** 2 + (ay - y) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best = az
    return best


def _seed_vertices(polygon: shp_geom.Polygon,
                   anchors: Sequence[Anchor],
                   sample_dem: DemSampler) -> _VertexBag:
    """Build the initial vertex bag: polygon exterior + interior
    ring (hole) vertices + any interior anchor points.

    Polygon vertices that coincide with an anchor adopt that
    anchor's elevation; otherwise they sample DEM.

    INTERIOR RING vertices are critical: a polygon with holes (e.g.
    a junction zone with building-shaped holes carved by
    ``difference``) needs hole vertices in the triangulation seed,
    or Delaunay produces large triangles that span across the holes.
    The centroid-based clip in :func:`_delaunay_clipped` only rejects
    triangles whose centroid lands inside a hole — large triangles
    whose centroid is elsewhere still cover the hole geometry.
    Adding hole vertices forces Delaunay to wrap the hole tightly so
    the centroid check works.
    """
    bag = _VertexBag()
    for x, y in polygon.exterior.coords[:-1]:
        anc_z = _find_anchor_at(x, y, anchors, EPS_DUP_M)
        if anc_z is not None:
            bag.add(x, y, anc_z, anchored=True)
        else:
            dem_z = sample_dem(x, y)
            bag.add(x, y, dem_z if dem_z is not None else 0.0,
                    anchored=False)
    # Interior ring (hole) vertices.
    for hole in polygon.interiors:
        for x, y in hole.coords[:-1]:
            anc_z = _find_anchor_at(x, y, anchors, EPS_DUP_M)
            if anc_z is not None:
                bag.add(x, y, anc_z, anchored=True)
            else:
                dem_z = sample_dem(x, y)
                bag.add(x, y, dem_z if dem_z is not None else 0.0,
                        anchored=False)
    # Interior anchor points (e.g. a taxiway centerline crossing an
    # apron's middle).
    for ax, ay, az in anchors:
        if bag.has_near(ax, ay, EPS_DUP_M):
            continue
        try:
            if polygon.contains(shp_geom.Point(ax, ay)):
                bag.add(ax, ay, az, anchored=True)
        except Exception:
            pass
    return bag


def _delaunay_clipped(bag: _VertexBag,
                      polygon: shp_geom.Polygon
                      ) -> List[shp_geom.Polygon]:
    """Delaunay-triangulate the vertex bag and keep only triangles
    whose centroid lies strictly inside `polygon` AND whose area
    is at least 99 % inside (to reject triangles that span concave
    bays, see commit 6 Phase D rewrite).

    Uses shapely ``prepared`` geometry for the centroid containment
    test (a 5-10× speedup over the raw ``polygon.contains`` path
    which has to re-index the polygon on every call) and a
    strict-vertex shortcut that skips the expensive
    ``polygon.intersection(t)`` check when the triangle's three
    corners are all strictly inside — in that case the 99 %
    containment is trivially satisfied.
    """
    try:
        mp = shp_geom.MultiPoint([(p[0], p[1]) for p in bag])
        raw = shp_ops.triangulate(mp)
    except Exception:
        return []
    try:
        from shapely.prepared import prep
        prepared = prep(polygon)
    except Exception:
        prepared = None
    # If the polygon has no interior rings, a triangle with all 3
    # corners inside the prepared polygon is trivially 100 %
    # contained — no need for the expensive
    # ``polygon.intersection(t).area`` check.  When the polygon
    # has holes, fall through to the full intersection test.
    has_holes = bool(list(polygon.interiors))
    out = []
    for t in raw:
        if t.is_empty or not t.is_valid:
            continue
        try:
            centroid = t.centroid
            if prepared is not None:
                if not prepared.contains(centroid):
                    continue
            else:
                if not polygon.contains(centroid):
                    continue
            if prepared is not None and not has_holes:
                # Hole-free shortcut: check the 3 raw vertex coords
                # with prepared.contains (no new Point creation
                # needed via the direct coordinate-bounds trick —
                # we must still build Points, but this is much
                # cheaper than polygon.intersection(t) because
                # prepared.contains is O(log n) against a
                # pre-indexed polygon).
                coords = list(t.exterior.coords)[:-1]
                if len(coords) == 3:
                    p0 = shp_geom.Point(coords[0])
                    p1 = shp_geom.Point(coords[1])
                    p2 = shp_geom.Point(coords[2])
                    if (prepared.contains(p0)
                            and prepared.contains(p1)
                            and prepared.contains(p2)):
                        out.append(t)
                        continue
            inter = polygon.intersection(t).area
            if inter >= 0.99 * t.area:
                out.append(t)
        except Exception:
            pass
    return out


def _triangle_corner_indices(t: shp_geom.Polygon,
                             bag: _VertexBag
                             ) -> Optional[List[int]]:
    """Return the bag indices of the three corners of `t`, or None if
    any corner is unmappable.
    """
    coords = list(t.exterior.coords)[:-1]
    if len(coords) != 3:
        return None
    idxs = [bag.index_at(x, y) for x, y in coords]
    if any(i is None for i in idxs):
        return None
    return idxs


# ──────────────────────────────────────────────────────────────────────
# Internal: grade enforcement
# ──────────────────────────────────────────────────────────────────────
def _grade_clamp(triangles: List[shp_geom.Polygon],
                 bag: _VertexBag,
                 max_grade: float,
                 passes: int = GRADE_CLAMP_PASSES) -> bool:
    """Triangle-plane-based grade relaxation with per-vertex
    accumulation.

    The "in any direction" grade rule (≤ ``max_grade`` for any slope
    direction in any emitted patch) bounds the plane *gradient*, not
    just the edge gradients along the polygon's sides.  An
    edge-based relaxation can leave the plane gradient up to √2 times
    the edge limit on right-angled triangles, which is too steep.
    So we enforce the plane gradient directly:

    1. For each pass, walk every triangle.  If its plane gradient
       exceeds ``max_grade``, compute the desired ΔZ that pulling
       each non-anchored corner toward the centroid by an exact
       fraction would produce — the fraction
       ``α = 1 − max_grade / current_grade`` is the closed-form
       solution for "scale the gradient down to exactly max_grade
       while keeping the centroid z fixed".
    2. Accumulate the per-vertex ΔZ contributions across ALL
       triangles incident to each vertex (sum of contributions,
       counted once per triangle).
    3. After processing all triangles, average each vertex's
       accumulated ΔZ over the number of triangles that contributed
       and apply once.  This avoids the double-counting bug where
       a vertex shared between N triangles got pulled N times per
       pass and over-flattened the surface.
    4. Repeat until no triangle exceeds ``max_grade``.

    Anchored vertices NEVER receive a delta.  Triangles where all
    three corners are anchored are skipped (the violation is
    physically unresolvable; the caller's reconciliation step is
    expected to adjust building elevations).

    Returns True if any over-steep triangle remains after the clamp
    finishes (i.e. at least one all-anchored triangle is impossible).
    """
    nv = len(bag)
    for _ in range(passes):
        deltas = [0.0] * nv
        weights = [0] * nv
        any_violation = False
        for t in triangles:
            idxs = _triangle_corner_indices(t, bag)
            if idxs is None:
                continue
            p0 = bag[idxs[0]]
            p1 = bag[idxs[1]]
            p2 = bag[idxs[2]]
            plane = fit_plane(
                (p0[0], p0[1], p0[2]),
                (p1[0], p1[1], p1[2]),
                (p2[0], p2[1], p2[2]))
            if plane is None:
                continue
            grade = plane_grade(plane)
            if grade <= max_grade + 1e-9:
                continue
            any_violation = True
            if p0[3] and p1[3] and p2[3]:
                continue   # all anchored — unresolvable
            cz = (p0[2] + p1[2] + p2[2]) / 3.0
            alpha = 1.0 - max_grade / grade  # exact closed-form pull
            for k, idx in enumerate(idxs):
                p = bag[idx]
                if p[3]:
                    continue
                deltas[idx] += alpha * (cz - p[2])
                weights[idx] += 1
        any_changed = False
        for i in range(nv):
            if weights[i] == 0:
                continue
            avg_delta = deltas[i] / weights[i]
            if abs(avg_delta) > 1e-6:
                bag[i][2] += avg_delta
                any_changed = True
        if not any_changed:
            return any_violation
    # Final scan
    for t in triangles:
        idxs = _triangle_corner_indices(t, bag)
        if idxs is None:
            continue
        plane = fit_plane(
            (bag[idxs[0]][0], bag[idxs[0]][1], bag[idxs[0]][2]),
            (bag[idxs[1]][0], bag[idxs[1]][1], bag[idxs[1]][2]),
            (bag[idxs[2]][0], bag[idxs[2]][1], bag[idxs[2]][2]))
        if plane is None:
            continue
        if plane_grade(plane) > max_grade + 1e-6:
            return True
    return False


# ──────────────────────────────────────────────────────────────────────
# Internal: fidelity refinement
# ──────────────────────────────────────────────────────────────────────
def _find_worst_fidelity_sample(
    triangles: List[shp_geom.Polygon],
    bag: _VertexBag,
    polygon: shp_geom.Polygon,
    sample_dem: DemSampler,
    fidelity_tol: float,
) -> Optional[Vertex3D]:
    """Walk every triangle, sample DEM at its centroid + edge
    midpoints, and return the (x, y, z) of the sample whose deviation
    from the current triangle's plane is largest, if it exceeds
    ``fidelity_tol``.

    This runs against the RAW DEM-driven planes BEFORE grade clamping
    — see the algorithm note in :func:`adaptive_triangulate`.  No
    grade compatibility check is needed because the clamp comes later
    as a single shot, not interleaved with refinement.

    Returns None if no triangle has a sample beyond the tolerance.
    """
    worst_dev = fidelity_tol
    worst_pt: Optional[Vertex3D] = None
    for t in triangles:
        idxs = _triangle_corner_indices(t, bag)
        if idxs is None:
            continue
        corners3: List[Vertex3D] = [
            (bag[i][0], bag[i][1], bag[i][2]) for i in idxs
        ]
        plane = fit_plane(*corners3)
        if plane is None:
            continue
        a, b, c = plane
        (x0, y0, _z0) = corners3[0]
        (x1, y1, _z1) = corners3[1]
        (x2, y2, _z2) = corners3[2]
        samples = (
            ((x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0),
            ((x0 + x1) / 2.0, (y0 + y1) / 2.0),
            ((x1 + x2) / 2.0, (y1 + y2) / 2.0),
            ((x2 + x0) / 2.0, (y2 + y0) / 2.0),
        )
        for sx, sy in samples:
            if bag.has_near(sx, sy, SAMPLE_GAP_M):
                continue
            try:
                # `covers` accepts boundary points so a sample on a
                # shared triangle edge (e.g. the midpoint of a Delaunay
                # diagonal) is valid.  `contains` would reject such
                # samples and miss bumps that land on the diagonal.
                if not polygon.covers(shp_geom.Point(sx, sy)):
                    continue
            except Exception:
                continue
            dem_z = sample_dem(sx, sy)
            if dem_z is None:
                continue
            plane_z = a * sx + b * sy + c
            dev = abs(dem_z - plane_z)
            if dev > worst_dev:
                worst_dev = dev
                worst_pt = (sx, sy, float(dem_z))
    return worst_pt


# ──────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────
def adaptive_triangulate(
    polygon: shp_geom.Polygon,
    anchors: Iterable[Anchor],
    sample_dem: DemSampler,
    max_grade: float,
    fidelity_tol: float = DEFAULT_FIDELITY_TOL_M,
    max_extra_points: int = DEFAULT_MAX_EXTRA_POINTS,
) -> List[Triangle]:
    """Triangulate `polygon` into the smallest set of grade-compliant,
    DEM-faithful triangles.  See module docstring for the algorithm.

    Args:
        polygon: shapely Polygon in meter-space to triangulate.
        anchors: iterable of (x, y, z) anchors from neighbouring
            features.  Anchors that lie on the polygon exterior pin
            the corresponding exterior vertex; anchors in the polygon
            interior become additional triangulation vertices.
        sample_dem: callable (x, y) -> float (elevation) or None.
            Called for every non-anchored vertex and for every
            candidate fidelity refinement sample.
        max_grade: maximum allowed plane gradient per triangle (e.g.
            0.010 for aprons, 0.015 for taxiways and junctions).
        fidelity_tol: refine until every triangle's interior agrees
            with DEM within this many metres.
        max_extra_points: cap on the number of fidelity refinement
            inserts to prevent runaway on noisy DEM.

    Returns:
        A list of Triangle tuples ``((x0,y0,z0),(x1,y1,z1),(x2,y2,z2))``
        in meter space.  All triangle planes satisfy ``max_grade``
        unless every corner of the offending triangle is anchored.
    """
    if polygon.is_empty or not hasattr(polygon, "exterior"):
        return []

    anchors_list: List[Anchor] = list(anchors) if anchors else []

    bag = _seed_vertices(polygon, anchors_list, sample_dem)
    if len(bag) < 3:
        return []

    triangles = _delaunay_clipped(bag, polygon)
    if not triangles:
        return []

    # Phase 1: Fidelity refinement against the raw DEM, BEFORE grade
    # clamping.  This decoupling is essential: if we interleaved
    # refinement and clamping, the clamp would distort each plane
    # away from DEM, the next refinement pass would see those
    # distortions as new "deviations" and insert points to chase
    # them, the next clamp would distort further, and so on — a
    # feedback loop that can't converge on uniformly-steep DEM.
    #
    # Running refinement first means it sees stable DEM-driven planes.
    # On a perfectly linear DEM the initial 2 triangles already
    # match exactly and refinement inserts nothing.  On bumpy terrain
    # refinement converges to capture the bumps with the minimum
    # number of inserts.
    extras = 0
    while extras < max_extra_points:
        candidate = _find_worst_fidelity_sample(
            triangles, bag, polygon, sample_dem, fidelity_tol)
        if candidate is None:
            break
        bag.add(candidate[0], candidate[1], candidate[2], anchored=False)
        extras += 1
        new_tris = _delaunay_clipped(bag, polygon)
        if not new_tris:
            break
        triangles = new_tris

    # Phase 2: Single grade-clamping pass (with internal iterations).
    # No more refinement after this point.  The clamp pulls
    # non-anchored vertices toward triangle centroids until every
    # triangle's plane gradient ≤ max_grade; bumps captured by
    # refinement get smoothed to whatever amplitude the grade rule
    # allows, which is exactly the user's "preserve detail but
    # enforce slope" requirement.
    _grade_clamp(triangles, bag, max_grade)

    out: List[Triangle] = []
    for t in triangles:
        idxs = _triangle_corner_indices(t, bag)
        if idxs is None:
            continue
        v0 = (bag[idxs[0]][0], bag[idxs[0]][1], bag[idxs[0]][2])
        v1 = (bag[idxs[1]][0], bag[idxs[1]][1], bag[idxs[1]][2])
        v2 = (bag[idxs[2]][0], bag[idxs[2]][1], bag[idxs[2]][2])
        out.append((v0, v1, v2))
    return out
