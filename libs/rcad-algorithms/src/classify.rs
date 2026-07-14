//! Point and shape classification algorithms (OCCT BRepClass3d equivalent).
//!
//! OCCT-aligned:
//! - `SolidClassifier` (class): classify point relative to a solid using
//!   `new(brep, solid_ref)` / `perform(point, tol)` / `state()` / `is_done()`.
//! - DS-based functions: classify_point, classify_solid_in_solid for the
//!   boolean pipeline. These are rcad-internal and use the DS data model.

use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods;
use std::collections::HashMap;
use std::sync::Arc;

use crate::bopds::ds::*;
use crate::bvh::{Aabb, Bvh};
use crate::inttools;
use crate::tolerance::{
    AdaptiveTolerance, ToleranceContext, ToleranceLevel, TOLERANCE_CLAMP_MIN,
    TOLERANCE_LEN_MIN, TOLERANCE_LEN_SQ_DIV_SAFE, TOLERANCE_MESH_LEGACY,
};

// =============================================================================
// BRepClass3d_SolidClassifier  ?OCCT-aligned class
// =============================================================================

/// OCCT-aligned: BRepClass3d_SolidClassifier  ?classify point relative to a solid.
///
/// Uses ray casting: cast a ray from the point in the +X direction, count
/// face intersections. Odd  ?In, Even  ?Out. Also checks vertex/edge proximity
/// for On classification.
///
/// Usage:
/// ```rust,ignore
/// let mut cls = SolidClassifier::new(&brep, solid_ref);
/// cls.perform(point, tol);
/// match cls.state() { In => ..., Out => ..., On => ... }
/// ```
pub struct SolidClassifier<'a> {
    brep: &'a topods::BRep,
    solid_ref: topods::ShapeRef,
    state: Classification,
    performed: bool,
}

impl<'a> SolidClassifier<'a> {
    /// OCCT-aligned: constructor with solid.
    pub fn new(brep: &'a topods::BRep, solid_ref: topods::ShapeRef) -> Self {
        Self { brep, solid_ref, state: Classification::Out, performed: false }
    }

    /// OCCT-aligned: Perform  ?classify point against the solid with tolerance.
    pub fn perform(&mut self, point: DVec3, tol: f64) {
        // Collect all face ShapeRefs from this solid
        let face_refs = collect_solid_faces(self.brep, self.solid_ref);
        if face_refs.is_empty() {
            self.state = Classification::Out;
            self.performed = true;
            return;
        }

        // 1. Check vertex/edge proximity  ?On
        for &fr in &face_refs {
            if let topods::TShape::Face(fd) = &*self.brep.tshapes[fr.index] {
                // Check if point is near the face surface
                if let Some(ref surf) = fd.surface {
                    let proj = rcad_kernel::projection::closest_point_on_surface(surf, point, 16);
                    if proj.distance < tol {
                        self.state = Classification::On;
                        self.performed = true;
                        return;
                    }
                }
            }
        }

        // 2. Ray casting: cast ray along +X, count face intersections
        let ray_dir = DVec3::X;
        let mut intersections = 0usize;
        for &fr in &face_refs {
            if let topods::TShape::Face(fd) = &*self.brep.tshapes[fr.index] {
                if let Some(ref surf) = fd.surface {
                    if let Some(t) = ray_face_intersect(point, ray_dir, surf, tol) {
                        if t > tol {
                            // Check if the hit point is within the face boundaries
                            let hit_pt = point + ray_dir * t;
                            let proj = rcad_kernel::projection::closest_point_on_surface(surf, hit_pt, 16);
                            // Approximate UV-boundary check via sampling the face wire
                            if proj.distance < tol * 10.0 {
                                intersections += 1;
                            }
                        }
                    }
                }
            }
        }

        self.state = if intersections % 2 == 1 { Classification::In } else { Classification::Out };
        self.performed = true;
    }

    /// OCCT-aligned: Perform with solid set in constructor.
    pub fn perform_with_point(&mut self, point: DVec3, tol: f64) {
        self.perform(point, tol);
    }

    /// OCCT-aligned: State  ?classification result.
    pub fn state(&self) -> Classification { self.state }

    /// OCCT-aligned: IsDone.
    pub fn is_done(&self) -> bool { self.performed }
}

/// Collect all face ShapeRefs from a solid in a topods::BRep.
fn collect_solid_faces(brep: &topods::BRep, solid_ref: topods::ShapeRef) -> Vec<topods::ShapeRef> {
    let mut faces = Vec::new();
    if let topods::TShape::Solid(sd) = &*brep.tshapes[solid_ref.index] {
        for &sh_ref in &sd.shells {
            if let topods::TShape::Shell(shd) = &*brep.tshapes[sh_ref.index] {
                faces.extend(shd.faces.iter().copied());
            }
        }
    }
    faces
}

/// Compute ray-surface intersection parameter t for a ray P + t*D.
fn ray_face_intersect(origin: DVec3, dir: DVec3, surf: &Surface3, _tol: f64) -> Option<f64> {
    match surf {
        Surface3::Plane(p) => {
            let denom = p.normal.dot(dir);
            if denom.abs() < TOLERANCE_CLAMP_MIN { return None; }
            let t = (p.origin - origin).dot(p.normal) / denom;
            if t >= 0.0 { Some(t) } else { None }
        }
        Surface3::Sphere(s) => {
            let oc = origin - s.center;
            let a = dir.dot(dir);
            let b = 2.0 * oc.dot(dir);
            let c = oc.dot(oc) - s.radius * s.radius;
            let disc = b * b - 4.0 * a * c;
            if disc < 0.0 { return None; }
            let sqrt_disc = disc.sqrt();
            let t1 = (-b - sqrt_disc) / (2.0 * a);
            let t2 = (-b + sqrt_disc) / (2.0 * a);
            let t = if t1 >= 0.0 { t1 } else if t2 >= 0.0 { t2 } else { return None };
            Some(t)
        }
        Surface3::Cylinder(c) => {
            // Cylinder-ray intersection using analytic formula
            let (axis, x_axis, y_axis) = orthonormal_frame(c.axis, c.ref_dir);
            let oc = origin - c.origin;
            let dx = dir.dot(x_axis);
            let dy = dir.dot(y_axis);
            let ox = oc.dot(x_axis);
            let oy = oc.dot(y_axis);
            let a2 = dx * dx + dy * dy;
            let b2 = 2.0 * (ox * dx + oy * dy);
            let c2 = ox * ox + oy * oy - c.radius * c.radius;
            if a2.abs() < TOLERANCE_CLAMP_MIN { return None; }
            let disc = b2 * b2 - 4.0 * a2 * c2;
            if disc < 0.0 { return None; }
            let sqrt_disc = disc.sqrt();
            let t1 = (-b2 - sqrt_disc) / (2.0 * a2);
            let t2 = (-b2 + sqrt_disc) / (2.0 * a2);
            let t = if t1 >= 0.0 { t1 } else if t2 >= 0.0 { t2 } else { return None };
            Some(t)
        }
        _ => {
            // Generic: use polyline sampling + intersection
            let domain = surf.default_domain();
            let n = 20usize;
            let u0 = domain[0]; let u1 = domain[1];
            let v0 = domain[2]; let v1 = domain[3];
            let du = (u1 - u0) / n as f64;
            let dv = (v1 - v0) / n as f64;
            for i in 0..n {
                let u = u0 + (i as f64 + 0.5) * du;
                for j in 0..n {
                    let v = v0 + (j as f64 + 0.5) * dv;
                    let pt = surf.point_at(u, v);
                    let to_pt = pt - origin;
                    let t = to_pt.dot(dir);
                    if t < 0.0 { continue; }
                    let lateral = (to_pt - dir * t).length();
                    if lateral < TOLERANCE_MESH_LEGACY { return Some(t); }
                }
            }
            None
        }
    }
}

use rcad_kernel::geom::any_perpendicular;
fn orthonormal_frame(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    let z = axis.normalize_or_zero();
    let mut x = ref_dir - z * ref_dir.dot(z);
    if x.length_squared() < 1e-24 { x = any_perpendicular(z); }
    x = x.normalize();
    let y = z.cross(x);
    (z, x, y)
}

// =============================================================================
// Classification Types
// =============================================================================

/// Classification of a point relative to a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Classification {
    In,
    Out,
    On,
}

impl Classification {
    /// Returns true if the point is inside or on the boundary.
    pub fn is_inside_or_on(self) -> bool {
        matches!(self, Classification::In | Classification::On)
    }

    /// Returns true if the point is strictly inside.
    pub fn is_inside(self) -> bool {
        self == Classification::In
    }

    /// Returns true if the point is on the boundary.
    pub fn is_on(self) -> bool {
        self == Classification::On
    }

    /// Negate the classification (swap In/Out, keep On).
    pub fn negate(self) -> Self {
        match self {
            Classification::In => Classification::Out,
            Classification::Out => Classification::In,
            Classification::On => Classification::On,
        }
    }
}

/// Classification of one solid relative to another.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidClassification {
    /// The solid is entirely outside the other.
    Outside,
    /// The solid is entirely inside the other.
    Inside,
    /// The solids partially overlap.
    Overlapping,
    /// The solids share a boundary but don't overlap in volume.
    Touching,
    /// The solids are identical (or within tolerance).
    Identical,
}

/// Classification of a point relative to a face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceClassification {
    /// Point is inside the face (projected point lies within face boundary).
    Inside,
    /// Point is outside the face.
    Outside,
    /// Point is on the face boundary (within tolerance).
    OnBoundary,
    /// Point is on the face surface (within tolerance).
    OnSurface,
}

/// Classification of a point relative to an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeClassification {
    /// Point is on the edge (within tolerance).
    OnEdge,
    /// Point is near the edge but not on it.
    Near,
    /// Point is far from the edge.
    Off,
}

// =============================================================================
// Classification Context with Caching
// =============================================================================

/// Cached classification data for a solid.
struct SolidClassifyCache {
    /// Face indices for the solid.
    face_indices: Vec<usize>,
    /// Bounding box of the solid.
    aabb: Aabb,
    /// BVH for fast face queries.
    bvh: Option<Bvh>,
    /// Precomputed face AABBs.
    face_aabbs: Vec<Aabb>,
}

/// Classification context with caching for repeated queries.
///
/// This provides significant performance improvements when classifying
/// multiple points against the same solid.
pub struct ClassifyContext {
    ds: Arc<DS>,
    tolerance: AdaptiveTolerance,
    /// Same role as [`ToleranceContext::workspace_fuzzy`]: lower bound on linear bands for
    /// coarse AABB expansion and relaxed-on-surface checks (boolean fuzzy / user workspace).
    workspace_fuzzy: f64,
    /// Cache keyed by a hash of face indices.
    cache: HashMap<u64, SolidClassifyCache>,
}

impl ClassifyContext {
    /// Create a new classification context.
    pub fn new(ds: Arc<DS>) -> Self {
        let tolerance = AdaptiveTolerance::from_scale(ds.model_scale());
        Self {
            ds,
            tolerance,
            workspace_fuzzy: 0.0,
            cache: HashMap::new(),
        }
    }

    /// Create a context with a custom tolerance.
    pub fn with_tolerance(ds: Arc<DS>, tolerance: AdaptiveTolerance) -> Self {
        Self {
            ds,
            tolerance,
            workspace_fuzzy: 0.0,
            cache: HashMap::new(),
        }
    }

    /// Build from a [`ToleranceContext`] (adaptive scale + optional workspace fuzzy).
    ///
    /// Relaxed / coarse linear bases use `max(adaptive(level), workspace_fuzzy)`, consistent with
    /// [`ToleranceContext::workspace_linear`].
    pub fn with_tolerance_context(ds: Arc<DS>, ctx: ToleranceContext) -> Self {
        let workspace_fuzzy = ctx.workspace_fuzzy.max(0.0);
        let out = Self {
            ds,
            tolerance: ctx.adaptive,
            workspace_fuzzy,
            cache: HashMap::new(),
        };
        tracing::debug!(
            target: "rcad.classify",
            model_scale = out.tolerance.model_scale,
            workspace_fuzzy = out.workspace_fuzzy,
            fine_linear = out.tolerance.tolerance(ToleranceLevel::Strict),
            relaxed_linear = out.workspace_linear(ToleranceLevel::Relaxed),
            "ClassifyContext::with_tolerance_context"
        );
        out
    }

    #[inline]
    fn workspace_linear(&self, level: ToleranceLevel) -> f64 {
        self.tolerance.tolerance(level).max(self.workspace_fuzzy)
    }

    /// Get or create cache for a solid.
    fn get_or_create_cache(&mut self, solid_face_indices: &[usize]) -> &SolidClassifyCache {
        // Canonical order so cache keys and `check_point_on_face` iteration are deterministic
        // (same set of face indices, independent of slice order in the input).
        let mut sorted = solid_face_indices.to_vec();
        sorted.sort_unstable();
        let hash = sorted.iter().fold(0u64, |h, &fi| {
            h.wrapping_mul(31).wrapping_add(fi as u64)
        });

        if !self.cache.contains_key(&hash) {
            let aabb = self.compute_solid_aabb(&sorted);
            let face_aabbs = self.compute_face_aabbs(&sorted);
            let bvh = if sorted.len() > 8 {
                Some(self.build_solid_bvh(&sorted, &face_aabbs))
            } else {
                None
            };

            self.cache.insert(
                hash,
                SolidClassifyCache {
                    face_indices: sorted,
                    aabb,
                    bvh,
                    face_aabbs,
                },
            );
        }

        self.cache.get(&hash).unwrap()
    }

    fn compute_solid_aabb(&self, face_indices: &[usize]) -> Aabb {
        let mut aabb = Aabb::empty();
        for &fi in face_indices {
            let Some(face) = self.ds.faces.get(fi) else { continue; };
            for &vi in &face.boundary_verts {
                aabb.expand_point(self.ds.vertex_point(vi));
            }
        }
        aabb
    }

    fn compute_face_aabbs(&self, face_indices: &[usize]) -> Vec<Aabb> {
        face_indices
            .iter()
            .filter_map(|&fi| {
                let face = self.ds.faces.get(fi)?;
                let mut aabb = Aabb::empty();
                for &vi in &face.boundary_verts {
                    aabb.expand_point(self.ds.vertex_point(vi));
                }
                Some(aabb)
            })
            .collect()
    }

    fn build_solid_bvh(&self, _face_indices: &[usize], _face_aabbs: &[Aabb]) -> Bvh {
        // Create a minimal BRep-like structure for BVH building
        // For now, we'll skip BVH and use linear search for simplicity
        // BVH optimization can be added later
        Bvh::build(&rcad_kernel::BRep::default())
    }

    /// Classify a point relative to a solid.
    pub fn classify_point(&mut self, point: DVec3, solid_face_indices: &[usize]) -> Classification {
        if solid_face_indices.is_empty() {
            return Classification::Out;
        }

        // Filter out invalid face indices to prevent cascading panics downstream.
        let valid: Vec<usize> = solid_face_indices.iter().filter(|&&fi| fi < self.ds.faces.len()).copied().collect();
        if valid.is_empty() { return Classification::Out; }

        // Extract tolerance before borrowing
        let mut sorted_for_tol = valid.to_vec();
        sorted_for_tol.sort_unstable();
        let face_geom_max = sorted_for_tol
            .iter()
            .filter_map(|&fi| self.ds.faces.get(fi).map(|f| f.geom_tol))
            .fold(0.0_f64, f64::max);
        let coarse_tol = crate::tolerance::effective_linear_with_geom_tol(
            self.workspace_linear(ToleranceLevel::Coarse),
            face_geom_max,
        );
        let tolerance = self.tolerance;
        let workspace_fuzzy = self.workspace_fuzzy;

        // Clone sorted face list so the cache borrow does not block `&self.ds` for classification.
        let face_indices_owned = {
            let cache = self.get_or_create_cache(&valid);

            // Quick AABB rejection test
            if !cache.aabb.contains_point(point) {
                let expanded_aabb = Aabb {
                    min: cache.aabb.min - DVec3::splat(coarse_tol),
                    max: cache.aabb.max + DVec3::splat(coarse_tol),
                };
                if !expanded_aabb.contains_point(point) {
                    return Classification::Out;
                }
            }

            cache.face_indices.clone()
        };

        classify_point_internal(
            point,
            &face_indices_owned,
            &self.ds,
            tolerance,
            workspace_fuzzy,
        )
    }

    /// Classify multiple points in parallel.
    pub fn classify_points_parallel(
        &mut self,
        points: &[DVec3],
        solid_face_indices: &[usize],
    ) -> Vec<Classification> {
        use std::thread;

        if points.is_empty() {
            return Vec::new();
        }

        // For small batches, use sequential classification
        if points.len() < 4 {
            return points
                .iter()
                .map(|&p| self.classify_point(p, solid_face_indices))
                .collect();
        }

        // Pre-compute cache
        let _cache = self.get_or_create_cache(solid_face_indices);

        // Split work across threads
        let n_threads = thread::available_parallelism().map(|p| p.get()).unwrap_or(4);
        let n_threads = n_threads.min(points.len()); // Don't create more threads than points
        let chunk_size = points.len().div_ceil(n_threads);

        let ds = Arc::clone(&self.ds);
        let tolerance = self.tolerance;
        let workspace_fuzzy = self.workspace_fuzzy;
        let mut face_indices = solid_face_indices.to_vec();
        face_indices.sort_unstable();

        let handles: Vec<_> = (0..n_threads)
            .map(|i| {
                let start = i * chunk_size;
                let end = ((i + 1) * chunk_size).min(points.len());
                let points_chunk = points[start..end].to_vec();
                let ds = Arc::clone(&ds);
                let face_indices = face_indices.clone();

                thread::spawn(move || {
                    points_chunk
                        .iter()
                        .map(|&p| {
                            classify_point_internal(
                                p,
                                &face_indices,
                                &ds,
                                tolerance,
                                workspace_fuzzy,
                            )
                        })
                        .collect::<Vec<_>>()
                })
            })
            .collect();

        let mut results = vec![Classification::Out; points.len()];
        for (i, handle) in handles.into_iter().enumerate() {
            let chunk_results = handle.join().unwrap();
            let start = i * chunk_size;
            for (j, result) in chunk_results.into_iter().enumerate() {
                results[start + j] = result;
            }
        }

        results
    }

    /// Clear the classification cache.
    pub fn clear_cache(&mut self) {
        self.cache.clear();
    }
}

// =============================================================================
// Relaxed linear tolerance with DS-propagated face tolerances (phase B+)
// =============================================================================

#[inline]
fn relaxed_tol_for_face_geom(
    ds: &DS,
    tol: AdaptiveTolerance,
    face_idx: usize,
    workspace_fuzzy: f64,
) -> f64 {
    let base = tol.tolerance(ToleranceLevel::Relaxed).max(workspace_fuzzy);
    ds.faces
        .get(face_idx)
        .map(|f| crate::tolerance::effective_linear_with_geom_tol(base, f.geom_tol))
        .unwrap_or(base)
}

#[inline]
fn relaxed_tol_for_solid_face_set(
    ds: &DS,
    tol: AdaptiveTolerance,
    solid_face_indices: &[usize],
    workspace_fuzzy: f64,
) -> f64 {
    let base = tol.tolerance(ToleranceLevel::Relaxed).max(workspace_fuzzy);
    let geom_max = solid_face_indices
        .iter()
        .filter_map(|&fi| ds.faces.get(fi).map(|f| f.geom_tol))
        .fold(0.0_f64, f64::max);
    crate::tolerance::effective_linear_with_geom_tol(base, geom_max)
}

// =============================================================================
// Core Classification Functions
// =============================================================================

///  ?OCCT-aligned: BRepClass3d_SolidClassifier::Perform (L171-211).
///   Classify a 3D point relative to a solid defined by face indices.
///
/// OCCT flow (BRepClass3d_SClassifier.cxx L203-253):
///   1. L207: SolidExplorer.Reject(P)  ?bounding box rejection  ?Out
///   2. L218-230: UB-tree select for vertex/edge proximity  ?On
///   3. L236+: Ray intersection with face  ?In/Out from face orientation
///
pub fn classify_point(point: DVec3, solid_face_indices: &[usize], ds: &DS) -> Classification {
    if solid_face_indices.is_empty() {
        return Classification::Out;
    }
    // Filter out invalid face indices (can happen when wire building produces
    // unexpected face splits; skip to avoid cascading panics downstream).
    let valid: Vec<usize> = solid_face_indices.iter().filter(|&&fi| fi < ds.faces.len()).copied().collect();
    if valid.is_empty() { return Classification::Out; }
    let mut sorted = valid;
    sorted.sort_unstable();
    let tol = AdaptiveTolerance::from_scale(ds.model_scale());
    let result = classify_point_internal(point, &sorted, ds, tol, 0.0);
    result
}

///  ?OCCT-aligned: BRepClass3d_SClassifier::Perform (L203-253).
///   Classify point via vertex/edge proximity + ray intersection.
///
/// OCCT flow:
///   1. L207: bounding box reject (SolidExplorer.Reject)  ?handled by caller
///   2. L218-230: vertex/edge proximity  ?On
///   3. L236+: ray intersect faces  ?In/Out from face transition
///
/// rcad: vertex/edge proximity (step 2) + face-on check + multi-ray voting (step 3).

/// OCCT Trans() L728-745: apply transition to state, reversing if parmin < 0.
///   tran: 0=IntCurveSurface_In, 1=IntCurveSurface_Out
///   my_state: 2=ON, 3=IN, 4=OUT
fn trans_helper(parmin: f64, tran: &mut i32, my_state: &mut i32) {
    if parmin < 0.0 {
        *tran = if *tran == 1 { 0i32 } else { 1i32 };
    }
    if *tran == 1 {
        *my_state = 3; // IN (line from inside to outside)
    } else {
        *my_state = 4; // OUT (line from outside to inside)
    }
}

fn classify_point_internal(
    point: DVec3,
    solid_face_indices: &[usize],
    ds: &DS,
    tol: AdaptiveTolerance,
    workspace_fuzzy: f64,
) -> Classification {
    let wf = workspace_fuzzy.max(0.0);
    let boundary_tol = relaxed_tol_for_solid_face_set(ds, tol, solid_face_indices, wf);
    let ray_tol = tol.tolerance(ToleranceLevel::Strict);
    let edge_tol = TOLERANCE_MESH_LEGACY.max(boundary_tol);

    // Build edge BVH (OCCT L214-215: aTree + aMapEV from SolidExplorer)
    let mut edge_aabbs: Vec<(usize, Aabb)> = Vec::new();
    let mut seen_edges = std::collections::HashSet::new();
    for &fi in solid_face_indices {
        if fi >= ds.faces.len() { continue; }
        for &ei in &ds.faces[fi].boundary_edges {
            if !seen_edges.insert(ei) { continue; }
            if let Some(edge) = ds.edges.get(ei) {
                let sv = ds.vertices[edge.start_vertex].point;
                let ev = ds.vertices[edge.end_vertex].point;
                let aabb = Aabb::from_points(&[sv, ev]);
                edge_aabbs.push((ei, aabb));
            }
        }
    }
    let edge_indices: Vec<usize> = edge_aabbs.iter().map(|(ei, _)| *ei).collect();
    let aabbs: Vec<Aabb> = edge_aabbs.iter().map(|(_, a)| *a).collect();
    let edge_tree = crate::bvh::DsBvh::build(edge_indices.clone(), aabbs);

    // OCCT L218-230: Vertex/Edge proximity via UB-tree -> On
    let query_aabb = Aabb {
        min: point - DVec3::splat(edge_tol),
        max: point + DVec3::splat(edge_tol),
    };
    for &ei in &edge_tree.query_aabb(&query_aabb) {
        if let Some(edge) = ds.edges.get(ei) {
            let sv = ds.vertices[edge.start_vertex].point;
            let ev = ds.vertices[edge.end_vertex].point;
            let seg = ev - sv;
            let seg_len_sq = seg.length_squared();
            if seg_len_sq < TOLERANCE_LEN_MIN * TOLERANCE_LEN_MIN {
                if point.distance_squared(sv) < edge_tol * edge_tol {
                    return Classification::On;
                }
                continue;
            }
            let t = ((point - sv).dot(seg) / seg_len_sq).clamp(0.0, 1.0);
            let closest = sv + t * seg;
            if point.distance_squared(closest) < edge_tol * edge_tol {
                return Classification::On;
            }
        }
    }
    for &fi in solid_face_indices {
        if fi >= ds.faces.len() { continue; }
        for &vi in &ds.faces[fi].boundary_verts {
            if let Some(v) = ds.vertices.get(vi) {
                if point.distance_squared(v.point) < edge_tol * edge_tol {
                    return Classification::On;
                }
            }
        }
    }

    // OCCT L232-234: mapEF - edge->face adjacency
    let mut map_ef: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for &fi in solid_face_indices {
        if fi >= ds.faces.len() { continue; }
        for &ei in &ds.faces[fi].boundary_edges {
            map_ef.entry(ei).or_default().push(fi);
        }
    }

    // Build direction list (OCCT L257-278: Segment + OtherSegment)
    struct DirEntry { dir: DVec3, max_par: f64 }
    let mut dirs: Vec<DirEntry> = Vec::new();
    for &fi in solid_face_indices {
        let face_verts = ds.face_boundary_points(fi);
        if face_verts.len() < 3 { continue; }
        let centroid = face_verts.iter().sum::<DVec3>() / face_verts.len() as f64;
        let ray_dir = centroid - point;
        let dir_len = ray_dir.length();
        if dir_len > TOLERANCE_LEN_SQ_DIV_SAFE {
            dirs.push(DirEntry { dir: ray_dir / dir_len, max_par: dir_len * 10.0 });
        }
    }
    for &d in &[DVec3::X, DVec3::Y, DVec3::Z, -DVec3::X, -DVec3::Y, -DVec3::Z] {
        dirs.push(DirEntry { dir: d, max_par: 1e10 });
    }

    // OCCT L203-523: BRepClass3d_SClassifier::Perform
    //   myState: 0=unknown, 1=rejected, 2=ON, 3=IN, 4=OUT
    let mut my_state: i32 = 0;
    let mut is_faulty_line = true;
    let mut an_ind_face = 0usize;

    while is_faulty_line && an_ind_face < dirs.len() {
        let rd = &dirs[an_ind_face];

        // OCCT L259-266: Segment / OtherSegment  ?rcad: direction from list
        let i_flag = 0i32; // rcad: always valid direction (0=OK, 1=OnFace, 2=OUT, 3=bad)

        // OCCT L270-278: anIndFace tracking via GetFaceSegmentIndex
        let a_cur_ind = an_ind_face + 1;
        if a_cur_ind > an_ind_face {
            an_ind_face = a_cur_ind;
        } else {
            my_state = 1;
            break;
        }

        // OCCT L280-297: iFlag handling
        match i_flag {
            1 => { my_state = 2; break; } // OnFace -> ON
            2 => { my_state = 4; break; } // Outside -> OUT
            3 => continue,                 // bad face -> skip
            _ => {}
        }

        // OCCT L299-300: reset parmin
        is_faulty_line = false;
        let mut parmin = f64::MAX;
        let mut near_fault_par = f64::MAX;

        // OCCT L302-360: Line-edge proximity (aTree.Select(aSelectorLine))
        //   L310-326: vertex hits -> NearFaultPar
        //   L328-360: edge hits -> GetTransi -> Trans
        let line_ext = edge_tol * 100.0;
        let line_aabb = Aabb {
            min: point - DVec3::splat(line_ext),
            max: point + DVec3::splat(line_ext),
        };
        // OCCT L314-316: collect vertex hits
        let mut lv_ints: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &fi in solid_face_indices {
            for &vi in &ds.faces[fi].boundary_verts {
                if let Some(v) = ds.vertices.get(vi) {
                    let to_ray = (v.point - point).cross(rd.dir).length();
                    if to_ray < edge_tol * 50.0 {
                        let lp = (v.point - point).dot(rd.dir);
                        lv_ints.insert(vi);
                        if lp.abs() < near_fault_par.abs() {
                            near_fault_par = lp;
                        }
                    }
                }
            }
        }
        // OCCT L328-360: edge hits -> GetTransi
        for &ei in &edge_tree.query_aabb(&line_aabb) {
            if let Some(edge) = ds.edges.get(ei) {
                let sv = ds.vertices[edge.start_vertex].point;
                let ev = ds.vertices[edge.end_vertex].point;
                let to_ray = (sv - point).cross(rd.dir).length();
                if to_ray > edge_tol * 10.0 {
                    continue;
                }
                // OCCT L343-347: skip edges whose vertices are also hit
                if lv_ints.contains(&edge.start_vertex) || lv_ints.contains(&edge.end_vertex) {
                    continue;
                }
                // OCCT L349-360: GetTransi(f1, f2, EE, param, L, tran)
                if let Some(faces) = map_ef.get(&ei) {
                    if faces.len() >= 2 {
                        let f1 = faces[0];
                        let f2 = faces[1];
                        let n1 = ds.face_normal(f1);
                        let n2 = ds.face_normal(f2);
                        let ang_tol = 1e-12;
                        // OCCT L677-683: check line orthogonal to normals
                        if rd.dir.dot(n1).abs() < ang_tol || rd.dir.dot(n2).abs() < ang_tol {
                            near_fault_par = 0.0;
                            continue;
                        }
                        // OCCT L685-701: parallel normals via IsParallel
                        let n_dot = n1.dot(n2).abs();
                        let mut tran = 0i32; // 0=In, 1=Out (matches Trans helper)
                        let tst: i32;
                        if n_dot > 1.0 - ang_tol {
                            // OCCT L687-700: nf1 parallel nf2
                            let ang_d = n1.dot(rd.dir);
                            if ang_d.abs() < ang_tol {
                                near_fault_par = 0.0;
                                continue;
                            } else if ang_d > 0.0 {
                                tran = 1; // Out
                            } else {
                                tran = 0; // In
                            }
                            tst = 1;
                        } else {
                            // OCCT L703-723: non-parallel normals
                            let n_cross = n1.cross(n2);
                            let proj_l = n_cross.cross(rd.dir).cross(n_cross);
                            let proj_len = proj_l.length();
                            if proj_len < ang_tol {
                                near_fault_par = 0.0;
                                continue;
                            }
                            let proj_dir = proj_l / proj_len;
                            let f_ad = n1.dot(proj_dir);
                            let s_ad = n2.dot(proj_dir);
                            if f_ad < -ang_tol && s_ad < -ang_tol {
                                tran = 0; // In
                                tst = 1;
                            } else if f_ad > ang_tol && s_ad > ang_tol {
                                tran = 1; // Out
                                tst = 1;
                            } else {
                                tst = 0; // skip
                            }
                        }
                        // OCCT L351-359: apply Trans if valid
                        let edge_mid = (sv + ev) * 0.5;
                        let lpar = (edge_mid - point).dot(rd.dir);
                        if tst == 1 && lpar.abs() < parmin.abs() {
                            parmin = lpar;
                            trans_helper(parmin, &mut tran, &mut my_state);
                        } else if lpar.abs() < near_fault_par.abs() {
                            near_fault_par = lpar;
                        }
                    }
                }
            }
        }

        // OCCT L363-509: Face intersection loop
        for &fi in solid_face_indices {
            if my_state == 2 {
                break;
            }
            let face = &ds.faces[fi];
            // OCCT L366-370: Shell/Face reject (rcad: flat face list, no rejection)
            // OCCT L375-397: Intersector3d.Perform(L, minW, maxW)
            let min_w = -rd.max_par.max(10.0 * ray_tol + 0.01 * rd.max_par);
            let max_w = rd.max_par.min(1e10);
            let result = ray_cast_classify_point_on_face(
                point, rd.dir, fi, ds, boundary_tol, ray_tol,
            );
            match result {
                RayFaceResult::On => {
                    // OCCT L451-455: |parmin| <= Tol -> ON
                    // If the face contains the point, it's ON the boundary
                    my_state = 2;
                    break;
                }
                RayFaceResult::In(t) | RayFaceResult::Out(t) => {
                    // OCCT L444-488: process intersection point
                    // OCCT L446: |WParameter(i)| < |parmin| - PConfusion
                    if t.abs() < parmin.abs() - ray_tol {
                        parmin = t;
                        // OCCT L448: if |parmin| <= Tol -> ON
                        if parmin.abs() <= ray_tol {
                            my_state = 2;
                            break;
                        }
                        // OCCT L458-473: State == IN -> process transition
                        // RayFaceResult::In(t) = ray exits = transition=Out, state=IN
                        // RayFaceResult::Out(t) = ray enters = transition=In, state=IN
                        let mut tran = match result {
                            RayFaceResult::Out(_) => 0i32, // In (entering)
                            _ => 1i32, // Out (exiting)
                        };
                        // OCCT L463-469: TANGENT -> continue
                        // (rcad doesn't produce tangent transitions from face intersection)
                        // OCCT L472: Trans(parmin, tran, myState)
                        trans_helper(parmin, &mut tran, &mut my_state);
                    }
                }
                RayFaceResult::Faulty => {
                    // OCCT L477-480: State == ON -> isFaultyLine
                    // Check if this face has the point on its boundary
                    // (rcad: Faulty can mean parallel or grazes boundary)
                    // The ray grazes this face; we don't break here,
                    // but we track via isFaultyLine if no valid intersection found
                    continue;
                }
            }
        }
        if my_state == 2 {
            break;
        }

        // OCCT L511-515: NearFaultPar vs parmin -> faulty line
        if near_fault_par.is_finite()
            && parmin.abs() >= near_fault_par.abs() - 1e-12
        {
            is_faulty_line = true;
        }
    }

    // OCCT L525-542: convert myState to Classification
    match my_state {
        2 => Classification::On,
        3 => Classification::In,
        4 => Classification::Out,
        _ => Classification::Out,
    }
}

enum RayFaceResult {
    /// Point is inside the solid (ray exits through this face).
    In(f64),
    /// Point is outside the solid.
    Out(f64),
    /// Point is ON the face boundary.
    On,
    /// Ray grazes face boundary  ?retry needed (faulty line).
    Faulty,
}

/// Fire a single ray from P in ray_dir toward face fi.
/// Returns In/Out/On/Faulty, matching OCCT's face intersection + transition logic.
fn ray_cast_classify_point_on_face(
    point: DVec3,
    ray_dir: DVec3,
    fi: usize,
    ds: &DS,
    boundary_tol: f64,
    ray_tol: f64,
) -> RayFaceResult {
    let face = &ds.faces[fi];

    match &face.surface {
        Surface3::Plane(plane) => {
            let denom = ray_dir.dot(plane.normal);
            if denom.abs() < ray_tol {
                return RayFaceResult::Faulty; // parallel ray
            }
            let t = (plane.origin - point).dot(plane.normal) / denom;
            if t < ray_tol {
                return RayFaceResult::Faulty; // intersection behind P
            }
            let hit = point + ray_dir * t;
            let face_verts = ds.face_boundary_points(fi);
            if is_near_polygon_boundary(&hit, &face_verts, plane, boundary_tol) {
                return RayFaceResult::Faulty; // grazes boundary  ?faulty
            }
            if !inttools::edge_face::point_in_planar_face_with_tol(
                hit, plane, &face_verts, boundary_tol,
            ) {
                return RayFaceResult::Faulty; // hit outside face bounds
            }
            // OCCT transition: ray enters solid when denom < 0 (ray opposite
            // to outward normal), exits when denom > 0 (ray along normal).
            // If ray enters: point was OUTSIDE before intersection.
            // If ray exits:  point was INSIDE before intersection.
            if denom < 0.0 {
                RayFaceResult::Out(t)  // ray enters solid  ?point was Out
            } else {
                RayFaceResult::In(t)   // ray exits solid  ?point was In
            }
        }
        Surface3::Sphere(s) => {
            let oc = point - s.center;
            let a = ray_dir.length_squared();
            if a < ray_tol { return RayFaceResult::Faulty; }
            let b = 2.0 * oc.dot(ray_dir);
            let cc = oc.length_squared() - s.radius * s.radius;
            let disc = b * b - 4.0 * a * cc;
            if disc < 0.0 { return RayFaceResult::Faulty; }
            let sq = disc.sqrt();
            let mut nearest = f64::MAX;
            let mut found = false;
            let face_verts = ds.face_boundary_points(fi);
            for &t in &[(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
                if t > ray_tol && t < nearest {
                    let hit = point + ray_dir * t;
                    // OCCT-aligned: UV-space containment (IntTools_FClass2d) instead
                    // of 3D AABB.  For periodic surfaces (sphere, cylinder), 3D
                    // boundary-vertex AABB may under-represent the face extent.
                    let in_face = if let Some(ref uv_bnd) = ds.faces[fi].uv_boundary {
                        let uv = s.world_to_uv(hit);
                        uv_bnd.len() >= 3 && point_in_uv_polygon(uv, uv_bnd)
                    } else if face_verts.len() < 3 {
                        true
                    } else {
                        point_in_face_aabb(hit, &face_verts, boundary_tol)
                    };
                    if in_face {
                        nearest = t;
                        found = true;
                    }
                }
            }
            if found {
                // OCCT transition: ray_enter check.  For sphere, outward
                // normal = (hit - center).normalize().  ray_dir·normal < 0
                // means entering  ?point was Out; > 0 means exiting  ?In.
                let hit = point + ray_dir * nearest;
                let n = (hit - s.center).normalize_or_zero();
                if ray_dir.dot(n) < 0.0 {
                    RayFaceResult::Out(nearest)
                } else {
                    RayFaceResult::In(nearest)
                }
            } else {
                RayFaceResult::Faulty
            }
        }
        _ => RayFaceResult::Faulty, // non-analytic surface  ?skip
    }
}



// =============================================================================
// Point-on-Face Classification
// =============================================================================

/// Check if a point is on a face surface within the face boundary.
pub fn classify_point_on_face(
    point: DVec3,
    face_idx: usize,
    ds: &DS,
    tolerance: f64,
) -> FaceClassification {
    let face = &ds.faces[face_idx];
    let surface = &face.surface;

    // Check distance to surface
    let dist_to_surface = distance_to_surface(point, surface);

    if dist_to_surface > tolerance {
        return FaceClassification::Outside;
    }

    // If on surface, check if within face boundary
    if dist_to_surface <= tolerance {
        // Project point to surface UV space
        let uv = project_point_to_surface_uv(point, surface);

        // Check if UV point is inside the face boundary
        let inside = if let Some(ref uv_boundary) = face.uv_boundary {
            point_in_uv_polygon(uv, uv_boundary)
        } else {
            // Fallback to 3D boundary check for planar faces
            match surface {
                Surface3::Plane(plane) => {
                    let face_verts = ds.face_boundary_points(face_idx);
                    inttools::edge_face::point_in_planar_face_with_tol(
                        point,
                        plane,
                        &face_verts,
                        tolerance,
                    )
                }
                _ => {
                    // For curved faces without UV boundary, use AABB approximation
                    let face_verts = ds.face_boundary_points(face_idx);
                    point_in_face_aabb(point, &face_verts, tolerance)
                }
            }
        };

        if inside {
            if dist_to_surface <= tolerance * 0.1 {
                FaceClassification::OnSurface
            } else {
                FaceClassification::Inside
            }
        } else {
            FaceClassification::Outside
        }
    } else {
        FaceClassification::Outside
    }
}


/// Check if two solids intersect by testing a point from one against the other's faces.
pub fn classify_solid_in_solid(
    solid_a_faces: &[usize],
    solid_b_faces: &[usize],
    ds: &DS,
    tolerance: f64,
) -> SolidClassification {
    if solid_a_faces.is_empty() || solid_b_faces.is_empty() {
        return SolidClassification::Outside;
    }

    // 1. Check bounding box relationship
    let aabb_a = compute_faces_aabb(solid_a_faces, ds);
    let aabb_b = compute_faces_aabb(solid_b_faces, ds);

    // No overlap in AABBs
    if !aabb_a.intersects(&aabb_b) {
        return SolidClassification::Outside;
    }

    // 2. Check if B is entirely inside A by sampling B's vertices
    let mut b_inside_a = true;
    let mut b_outside_a = false;
    let mut b_on_boundary = false;

    for &fi in solid_b_faces {
        let face = &ds.faces[fi];
        for &vi in &face.boundary_verts {
            let point = ds.vertex_point(vi);
            let class = classify_point(point, solid_a_faces, ds);
            match class {
                Classification::In => {}
                Classification::Out => {
                    b_inside_a = false;
                    b_outside_a = true;
                }
                Classification::On => {
                    b_on_boundary = true;
                }
            }
            if b_outside_a && b_on_boundary {
                break;
            }
        }
        if b_outside_a && b_on_boundary {
            break;
        }
    }

    // 3. Check if A has vertices inside B
    let mut a_inside_b = false;
    let mut a_outside_b = false;

    for &fi in solid_a_faces {
        let face = &ds.faces[fi];
        for &vi in &face.boundary_verts {
            let point = ds.vertex_point(vi);
            let class = classify_point(point, solid_b_faces, ds);
            match class {
                Classification::In => {
                    a_inside_b = true;
                }
                Classification::Out => {
                    a_outside_b = true;
                }
                Classification::On => {}
            }
            if a_inside_b && a_outside_b {
                break;
            }
        }
        if a_inside_b && a_outside_b {
            break;
        }
    }

    // 4. Determine relationship
    if b_inside_a && !b_outside_a {
        if b_on_boundary {
            SolidClassification::Touching
        } else {
            SolidClassification::Inside
        }
    } else if a_inside_b && !a_outside_b {
        SolidClassification::Overlapping // A is inside B
    } else if b_outside_a && !a_inside_b {
        if b_on_boundary {
            SolidClassification::Touching
        } else {
            // Need to check for partial overlap
            if aabb_a.intersects(&aabb_b) {
                SolidClassification::Overlapping
            } else {
                SolidClassification::Outside
            }
        }
    } else if b_outside_a && a_inside_b {
        SolidClassification::Overlapping
    } else {
        // Complex case: check face intersections
        let faces_intersect = check_face_intersections(solid_a_faces, solid_b_faces, ds, tolerance);
        if faces_intersect {
            SolidClassification::Overlapping
        } else if b_on_boundary {
            SolidClassification::Touching
        } else {
            SolidClassification::Outside
        }
    }
}

/// Compute AABB for a set of faces.
fn compute_faces_aabb(face_indices: &[usize], ds: &DS) -> Aabb {
    let mut aabb = Aabb::empty();
    for &fi in face_indices {
        let face = &ds.faces[fi];
        for &vi in &face.boundary_verts {
            aabb.expand_point(ds.vertex_point(vi));
        }
    }
    aabb
}

/// Check if any faces from two sets intersect.
fn check_face_intersections(
    faces_a: &[usize],
    faces_b: &[usize],
    ds: &DS,
    tolerance: f64,
) -> bool {
    // Quick AABB check for face pairs
    for &fi_a in faces_a {
        let face_a = &ds.faces[fi_a];
        let mut aabb_a = Aabb::empty();
        for &vi in &face_a.boundary_verts {
            aabb_a.expand_point(ds.vertex_point(vi));
        }

        for &fi_b in faces_b {
            if fi_a == fi_b {
                continue;
            }

            let face_b = &ds.faces[fi_b];
            let mut aabb_b = Aabb::empty();
            for &vi in &face_b.boundary_verts {
                aabb_b.expand_point(ds.vertex_point(vi));
            }

            if aabb_a.intersects(&aabb_b) {
                // Check if any vertex of B is inside A's face
                for &vi in &face_b.boundary_verts {
                    let point = ds.vertex_point(vi);
                    let class = classify_point_on_face(point, fi_a, ds, tolerance);
                    if matches!(class, FaceClassification::Inside | FaceClassification::OnSurface) {
                        return true;
                    }
                }
            }
        }
    }

    false
}

// =============================================================================
// Point-on-Edge Classification
// =============================================================================

/// Classify a point relative to an edge.
pub fn classify_point_on_edge(
    point: DVec3,
    edge_idx: usize,
    ds: &DS,
    tolerance: f64,
) -> EdgeClassification {
    let edge = &ds.edges[edge_idx];
    let curve = &edge.curve;

    // Project point onto curve
    let (closest_point, param) = project_point_to_curve(point, curve);

    // Check if parameter is within edge range
    let t_range = edge.t_range;
    let t_min = t_range[0].min(t_range[1]);
    let t_max = t_range[0].max(t_range[1]);

    let dist = (point - closest_point).length();

    if dist <= tolerance && param >= t_min - tolerance && param <= t_max + tolerance {
        EdgeClassification::OnEdge
    } else if dist <= tolerance * 10.0 {
        EdgeClassification::Near
    } else {
        EdgeClassification::Off
    }
}

/// Project a point onto a curve and return the closest point and parameter.
fn project_point_to_curve(point: DVec3, curve: &Curve3) -> (DVec3, f64) {
    match curve {
        Curve3::Line(line) => {
            let t = (point - line.origin).dot(line.direction);
            let closest = line.origin + line.direction * t;
            (closest, t)
        }
        Curve3::Circle(circle) => {
            let to_point = point - circle.center;
            let in_plane = to_point - circle.normal * to_point.dot(circle.normal);
            let t = in_plane.normalize_or_zero();
            let angle = t.x.atan2(t.y);
            let closest = circle.center + circle.radius * t;
            (closest, angle)
        }
        Curve3::Ellipse(ellipse) => {
            // Approximate projection for ellipse
            let to_point = point - ellipse.center;
            let angle = (to_point.x / ellipse.major_radius).atan2(to_point.y / ellipse.minor_radius);
            let minor_dir = ellipse.normal.cross(ellipse.major_dir).normalize_or_zero();
            let closest = ellipse.center
                + ellipse.major_dir * angle.cos() * ellipse.major_radius
                + minor_dir * angle.sin() * ellipse.minor_radius;
            (closest, angle)
        }
        _ => {
            // Generic case: sample curve and find closest
            let domain = curve.default_domain();
            let n_samples = 100;
            let mut best_dist = f64::INFINITY;
            let mut best_point = point;
            let mut best_param = domain[0];

            for i in 0..n_samples {
                let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
                let p = curve.point_at(t);
                let d = (p - point).length();
                if d < best_dist {
                    best_dist = d;
                    best_point = p;
                    best_param = t;
                }
            }

            (best_point, best_param)
        }
    }
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Compute the angular range [min, max] of a cylinder face in radians.
fn cylinder_face_angle_range(
    c: &CylindricalSurface,
    face_verts: &[DVec3],
    axis: DVec3,
) -> (f64, f64) {
    if face_verts.len() < 2 {
        return (0.0, std::f64::consts::TAU);
    }
    let angles: Vec<f64> = face_verts
        .iter()
        .map(|&v| {
            let radial = v - c.origin - axis * (v - c.origin).dot(axis);
            cylinder_angle(c, radial)
        })
        .collect();
    let mut min_a = angles[0];
    let mut max_a = angles[0];
    for &a in &angles[1..] {
        if a < min_a { min_a = a; }
        if a > max_a { max_a = a; }
    }
    if max_a - min_a < TOLERANCE_MESH_LEGACY {
        return (0.0, std::f64::consts::TAU);
    }
    if max_a - min_a > std::f64::consts::PI {
        let ref_a = angles[0];
        let wrapped: Vec<f64> = angles
            .iter()
            .map(|&a| {
                let mut d = a - ref_a;
                while d < 0.0 { d += std::f64::consts::TAU; }
                while d > std::f64::consts::TAU { d -= std::f64::consts::TAU; }
                d
            })
            .collect();
        let span = wrapped.iter().cloned().fold(0.0_f64, f64::max);
        if span < std::f64::consts::TAU - 0.01 {
            return (ref_a, ref_a + span);
        } else {
            return (0.0, std::f64::consts::TAU);
        }
    }
    if max_a - min_a < TOLERANCE_MESH_LEGACY {
        return (0.0, std::f64::consts::TAU);
    }
    (min_a, max_a)
}

/// Compute the angle of a radial vector relative to the cylinder's reference direction.
fn cylinder_angle(c: &CylindricalSurface, radial: DVec3) -> f64 {
    let axis = c.axis.normalize();
    let ref_dir = any_perpendicular(axis).normalize();
    let perp_dir = axis.cross(ref_dir).normalize();
    let x = radial.dot(ref_dir);
    let y = radial.dot(perp_dir);
    x.atan2(y)
}

/// Check if angle is within [min, max] range (with angular slack).
fn angle_in_range(angle: f64, min_a: f64, max_a: f64, slack: f64) -> bool {
    angle >= min_a - slack && angle <= max_a + slack
}

/// Conservative face containment check using AABB of the face boundary vertices.
fn point_in_face_aabb(point: DVec3, face_verts: &[DVec3], slack: f64) -> bool {
    if face_verts.is_empty() {
        return false;
    }
    let mut mn = face_verts[0];
    let mut mx = face_verts[0];
    for &v in face_verts.iter().skip(1) {
        mn = mn.min(v);
        mx = mx.max(v);
    }
    point.cmpge(mn - DVec3::splat(slack)).all() && point.cmple(mx + DVec3::splat(slack)).all()
}

/// Check if a point is close to any edge of a polygon (within tolerance).
fn is_near_polygon_boundary(point: &DVec3, verts: &[DVec3], plane: &Plane, boundary_tol: f64) -> bool {
    let (u_axis, v_axis) = inttools::edge_face::plane_local_basis(plane);
    let project = |p: DVec3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(u_axis), d.dot(v_axis))
    };

    let (px, py) = project(*point);
    let n = verts.len();
    let tol_sq = boundary_tol * boundary_tol;

    for i in 0..n {
        let j = (i + 1) % n;
        let (ax, ay) = project(verts[i]);
        let (bx, by) = project(verts[j]);

        let dx = bx - ax;
        let dy = by - ay;
        let len_sq = dx * dx + dy * dy;
        if len_sq < tol_sq {
            continue;
        }

        let t = ((px - ax) * dx + (py - ay) * dy) / len_sq;
        let t = t.clamp(0.0, 1.0);
        let cx = ax + t * dx;
        let cy = ay + t * dy;
        let dist_sq = (px - cx) * (px - cx) + (py - cy) * (py - cy);

        if dist_sq < tol_sq {
            return true;
        }
    }

    false
}

/// Compute distance from point to surface.
fn distance_to_surface(point: DVec3, surface: &Surface3) -> f64 {
    match surface {
        Surface3::Plane(plane) => {
            (point - plane.origin).dot(plane.normal).abs()
        }
        Surface3::Sphere(s) => {
            ((point - s.center).length() - s.radius).abs()
        }
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = (v - c.axis * along).length();
            (perp - c.radius).abs()
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis.normalize_or_zero();
            let apex = cone.apex_point();
            let v = point - apex;
            let along = v.dot(axis);
            let perp = (v - axis * along).length();
            let tan_a = cone.half_angle_rad.tan();
            (perp - along.max(0.0) * tan_a).abs()
        }
        Surface3::Torus(t) => {
            let axis = t.axis.normalize_or_zero();
            let delta = point - t.center;
            let z = delta.dot(axis);
            let radial = delta - axis * z;
            let rho = radial.length();
            let tube_dist = ((rho - t.major_radius).powi(2) + z * z).sqrt();
            (tube_dist - t.minor_radius).abs()
        }
        _ => {
            // Generic case: use projection
            let proj = rcad_kernel::projection::closest_point_on_surface(surface, point, 16);
            (point - proj.point).length()
        }
    }
}

/// Project point to surface UV coordinates.
fn project_point_to_surface_uv(point: DVec3, surface: &Surface3) -> glam::DVec2 {
    match surface {
        Surface3::Plane(plane) => {
            let (u_axis, v_axis) = inttools::edge_face::plane_local_basis(plane);
            let d = point - plane.origin;
            glam::DVec2::new(d.dot(u_axis), d.dot(v_axis))
        }
        Surface3::Sphere(s) => s.world_to_uv(point),
        Surface3::Cylinder(c) => {
            let axis = c.axis.normalize();
            let v = point - c.origin;
            let h = v.dot(axis);
            let radial = v - axis * h;
            let phi = cylinder_angle(c, radial);
            glam::DVec2::new(phi, h)
        }
        _ => {
            let proj = rcad_kernel::projection::closest_point_on_surface(surface, point, 16);
            glam::DVec2::new(proj.params.0, proj.params.1)
        }
    }
}

/// Check if a UV point is inside a UV polygon.
fn point_in_uv_polygon(point: glam::DVec2, polygon: &[glam::DVec2]) -> bool {
    if polygon.len() < 3 {
        return false;
    }

    let eps = 1e-12;
    let mut inside = false;
    let n = polygon.len();

    for i in 0..n {
        let j = (i + 1) % n;
        let xi = polygon[i].x;
        let yi = polygon[i].y;
        let xj = polygon[j].x;
        let yj = polygon[j].y;

        if ((yi > point.y) != (yj > point.y))
            && (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi + eps)
        {
            inside = !inside;
        }
    }

    inside
}

// =============================================================================
// Tests
// =============================================================================


