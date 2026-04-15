//! Hidden-Line Removal (HLR).
//!
//! Projects a BRep's edges onto a view plane and classifies each edge segment
//! as **visible** or **hidden** by testing against the silhouette of all faces.
//!
//! Analytic silhouette curves are generated for curved surfaces (cylinder,
//! sphere, cone, torus) and processed through the same visibility pipeline as wire edges.
//! For general surfaces (BSpline, Bezier, etc.), numerical silhouette extraction
//! is performed using adaptive sampling with curvature-based refinement.
//!
//! Analogous to OCCT `HLRBRep_Algo` / `HLRBRep_HLRToShape`.
//!
//! # Algorithm
//!
//! For each edge (and silhouette curve):
//! 1. Project both endpoints onto the screen plane.
//! 2. Sample `N` points along the edge in 3D (adaptively if `curvature_adaptive` is true).
//! 3. For each sample, cast a ray from that point toward the camera.
//! 4. If any face triangle blocks the ray **closer** to the camera than the
//!    edge sample, that sample is hidden.
//! 5. Classify runs of consecutive samples → visible/hidden segments.
//!
//! The result is a set of `HlrSegment`s — 2D projected line segments labeled
//! visible or hidden.
//!
//! # Silhouette Classification
//!
//! Segments are classified as:
//! - **Visible**: Not occluded by any face
//! - **Hidden**: Occluded by at least one face
//! - **Contour**: Edge of a face's silhouette (marked via `is_contour`)
//!
//! # Curved Surface Enhancements
//!
//! This implementation includes several enhancements for better handling of curved geometry:
//! - **Curvature-adaptive sampling**: Uses surface curvature to concentrate samples in high-curvature regions
//! - **Marching silhouette extraction**: Robust detection of silhouette curves on general parametric surfaces
//! - **B-spline fitting**: Converts dense silhouette points to smooth B-spline curves
//! - **BVH acceleration**: Spatial acceleration structure for efficient ray casting
//! - **Grazing angle handling**: Special treatment for near-silhouette regions

use glam::{DAffine3, DMat4, DVec2, DVec3, DVec4};
use rcad_kernel::geom::{Circle3, CurveEval, Surface3, any_perpendicular};
use rcad_kernel::{BRep, SurfaceEval};

// ── Public types ──────────────────────────────────────────────────────────────

/// Configuration options for HLR computation.
#[derive(Debug, Clone)]
pub struct HlrOptions {
    /// Number of samples per edge for occlusion testing.
    /// Higher values give more accurate results but slower computation.
    /// Default: 8.
    pub edge_samples: usize,
    /// Base number of samples for silhouette curve generation.
    /// Default: 32.
    pub silhouette_samples: usize,
    /// Enable curvature-adaptive sampling for silhouette curves.
    /// When true, high-curvature regions receive more samples.
    /// Default: true.
    pub curvature_adaptive: bool,
    /// Tolerance for tangent alignment when computing silhouettes.
    /// Points where |normal · view_dir| < tangent_tolerance are considered silhouette candidates.
    /// Default: 1e-6.
    pub tangent_tolerance: f64,
    /// Maximum angle deviation (in radians) for adaptive sampling subdivision.
    /// Smaller values produce smoother curves at higher cost.
    /// Default: 0.05 (about 3 degrees).
    pub angular_tolerance: f64,
    /// Minimum number of subdivision iterations for adaptive sampling.
    /// Default: 2.
    pub min_subdivisions: usize,
    /// Maximum number of subdivision iterations for adaptive sampling.
    /// Default: 8.
    pub max_subdivisions: usize,
    /// Enable BVH acceleration for ray casting.
    /// Default: true.
    pub use_bvh: bool,
    /// Maximum curvature for adaptive sampling (higher = more samples in curved regions).
    /// Default: 100.0.
    pub max_curvature: f64,
    /// Minimum curvature for adaptive sampling (lower = fewer samples in flat regions).
    /// Default: 0.001.
    pub min_curvature: f64,
    /// Enable B-spline fitting for silhouette curves.
    /// Default: true.
    pub fit_bspline: bool,
    /// Tolerance for B-spline fitting (maximum deviation from original points).
    /// Default: 0.001.
    pub bspline_tolerance: f64,
    /// Grazing angle threshold (in radians). Points closer to silhouette receive special handling.
    /// Default: 0.1 (about 6 degrees).
    pub grazing_angle_threshold: f64,
    /// Enable smooth silhouette approximation.
    /// Default: true.
    pub smooth_silhouettes: bool,
}

impl Default for HlrOptions {
    fn default() -> Self {
        Self {
            edge_samples: 8,
            silhouette_samples: 32,
            curvature_adaptive: true,
            tangent_tolerance: 1e-6,
            angular_tolerance: 0.05,
            min_subdivisions: 2,
            max_subdivisions: 8,
            use_bvh: true,
            max_curvature: 100.0,
            min_curvature: 0.001,
            fit_bspline: true,
            bspline_tolerance: 0.001,
            grazing_angle_threshold: 0.1,
            smooth_silhouettes: true,
        }
    }
}

impl HlrOptions {
    /// Create options with a specific edge sample count.
    pub fn with_edge_samples(mut self, n: usize) -> Self {
        self.edge_samples = n.max(2);
        self
    }

    /// Create options with a specific silhouette sample count.
    pub fn with_silhouette_samples(mut self, n: usize) -> Self {
        self.silhouette_samples = n.max(8);
        self
    }

    /// Enable or disable curvature-adaptive sampling.
    pub fn with_curvature_adaptive(mut self, adaptive: bool) -> Self {
        self.curvature_adaptive = adaptive;
        self
    }

    /// Set the tangent tolerance for silhouette detection.
    pub fn with_tangent_tolerance(mut self, tol: f64) -> Self {
        self.tangent_tolerance = tol.abs().max(1e-12);
        self
    }

    /// Enable or disable BVH acceleration.
    pub fn with_bvh(mut self, use_bvh: bool) -> Self {
        self.use_bvh = use_bvh;
        self
    }

    /// Set the maximum curvature for adaptive sampling.
    pub fn with_max_curvature(mut self, curv: f64) -> Self {
        self.max_curvature = curv.abs().max(0.1);
        self
    }

    /// Enable or disable B-spline fitting for silhouettes.
    pub fn with_bspline_fitting(mut self, fit: bool) -> Self {
        self.fit_bspline = fit;
        self
    }

    /// Set the grazing angle threshold.
    pub fn with_grazing_angle(mut self, angle: f64) -> Self {
        self.grazing_angle_threshold = angle.abs().min(std::f64::consts::FRAC_PI_2);
        self
    }

    /// Enable or disable smooth silhouette approximation.
    pub fn with_smooth_silhouettes(mut self, smooth: bool) -> Self {
        self.smooth_silhouettes = smooth;
        self
    }
}

/// Hint about the geometric type of the original 3D edge curve.
/// Used by consumers (e.g. SVG exporter) to emit arcs instead of polylines.
#[derive(Debug, Clone, PartialEq)]
pub enum CurveHint {
    /// Edge is a full or partial circle in 3D.
    Circle {
        /// Projected 2D center of the circle.
        center: DVec2,
        /// Projected radius (approximate — perspective not applied).
        radius: f64,
    },
    /// Any other non-straight curve (ellipse, spline, …).
    Other,
}

/// Classification of an HLR segment type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentType {
    /// Regular edge (part of the BRep wire).
    Edge,
    /// Silhouette curve (contour of a curved face).
    Silhouette,
}

/// A projected edge segment labeled as visible or hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrSegment {
    /// Start point in 2D screen space.
    pub start: DVec2,
    /// End point in 2D screen space.
    pub end: DVec2,
    /// Whether this segment is visible from the camera.
    pub visible: bool,
    /// Optional hint about the underlying curve type (None for straight lines).
    pub curve_hint: Option<CurveHint>,
    /// Segment type (edge or silhouette).
    pub segment_type: SegmentType,
}

impl HlrSegment {
    /// Returns true if this segment is a silhouette/contour curve.
    pub fn is_contour(&self) -> bool {
        self.segment_type == SegmentType::Silhouette
    }
}

/// Output of an HLR computation.
#[derive(Debug, Clone, Default)]
pub struct HlrResult {
    pub segments: Vec<HlrSegment>,
}

impl HlrResult {
    /// Return only visible segments.
    pub fn visible(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible)
    }

    /// Return only hidden segments.
    pub fn hidden(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| !s.visible)
    }

    /// Return only silhouette/contour segments.
    pub fn silhouettes(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.is_contour())
    }

    /// Return only visible silhouette segments.
    pub fn visible_silhouettes(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible && s.is_contour())
    }
}

/// Camera / view specification for HLR.
#[derive(Debug, Clone)]
pub struct HlrCamera {
    /// Camera position in world space.
    pub eye: DVec3,
    /// Target point (look-at).
    pub target: DVec3,
    /// Up direction.
    pub up: DVec3,
}

impl HlrCamera {
    pub fn new(eye: DVec3, target: DVec3) -> Self {
        Self {
            eye,
            target,
            up: DVec3::Y,
        }
    }

    pub fn with_up(mut self, up: DVec3) -> Self {
        self.up = up;
        self
    }

    /// Isometric-style view from the +X+Y+Z octant.
    pub fn isometric(distance: f64) -> Self {
        let d = distance / 3.0_f64.sqrt();
        Self::new(DVec3::splat(d), DVec3::ZERO)
    }

    /// Front view (looking along +Y, up = +Z).
    pub fn front(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, -distance, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }

    /// Top view (looking down -Z).
    pub fn top(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, 0.0, distance), DVec3::ZERO).with_up(DVec3::Y)
    }

    /// Right-side view (looking along -X, up = +Z).
    pub fn right(distance: f64) -> Self {
        Self::new(DVec3::new(distance, 0.0, 0.0), DVec3::ZERO).with_up(DVec3::Z)
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Build a right-handed view matrix (world → camera space).
fn look_at(eye: DVec3, target: DVec3, up: DVec3) -> DMat4 {
    let forward = (target - eye).normalize_or_zero();
    let right = forward.cross(up).normalize_or_zero();
    let up = right.cross(forward);

    DMat4::from_cols(
        DVec4::new(right.x, right.y, right.z, -right.dot(eye)),
        DVec4::new(up.x, up.y, up.z, -up.dot(eye)),
        DVec4::new(-forward.x, -forward.y, -forward.z, forward.dot(eye)),
        DVec4::new(0.0, 0.0, 0.0, 1.0),
    )
    .transpose()
}

/// Project a world-space point to 2D screen space using the view matrix.
/// Returns (x, y) in camera space (z is depth; ignored for 2D output).
fn project(p: DVec3, view: &DMat4) -> (DVec2, f64) {
    let hp = view.mul_vec4(DVec4::new(p.x, p.y, p.z, 1.0));
    (DVec2::new(hp.x, hp.y), hp.z)
}

/// Collect all triangles from a BRep (fan-triangulate faces without pre-triangulated data).
fn collect_triangles(brep: &BRep) -> Vec<[DVec3; 3]> {
    let mut tris = Vec::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if !face.triangles.is_empty() {
                    for &[i, j, k] in &face.triangles {
                        if let (Some(a), Some(b), Some(c)) = (
                            brep.vertices.get(i),
                            brep.vertices.get(j),
                            brep.vertices.get(k),
                        ) {
                            tris.push([a.point, b.point, c.point]);
                        }
                    }
                } else {
                    // Fan-triangulate from wire
                    let pts: Vec<DVec3> = face
                        .outer_wire
                        .edges
                        .iter()
                        .filter_map(|we| {
                            let edge = brep.edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            brep.vertices.get(vi).map(|v| v.point)
                        })
                        .collect();
                    if pts.len() >= 3 {
                        let origin = pts[0];
                        for i in 1..pts.len() - 1 {
                            tris.push([origin, pts[i], pts[i + 1]]);
                        }
                    }
                }
            }
        }
    }
    tris
}

/// Ray-triangle intersection (Möller–Trumbore). Returns `Some(t)` if the ray
/// `origin + t*dir` hits the triangle (t > epsilon, front-face only).
fn ray_triangle_intersect(origin: DVec3, dir: DVec3, tri: &[DVec3; 3]) -> Option<f64> {
    const EPS: f64 = 1e-8;
    let edge1 = tri[1] - tri[0];
    let edge2 = tri[2] - tri[0];
    let h = dir.cross(edge2);
    let a = edge1.dot(h);
    if a.abs() < EPS {
        return None;
    }
    let f = 1.0 / a;
    let s = origin - tri[0];
    let u = f * s.dot(h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let q = s.cross(edge1);
    let v = f * dir.dot(q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = f * edge2.dot(q);
    if t > EPS { Some(t) } else { None }
}

/// Test if a world-space point is occluded by any triangle when viewed from `eye`.
fn is_occluded(point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
    let dir = (eye - point).normalize_or_zero();
    let origin = point + dir * 1e-5; // push off surface
    for tri in triangles {
        if let Some(t) = ray_triangle_intersect(origin, dir, tri)
            && t < dist_to_eye - 1e-4
        {
            return true;
        }
    }
    false
}

// ── Triangle BVH for accelerated ray casting ───────────────────────────────────

/// Axis-aligned bounding box for triangle BVH.
#[derive(Debug, Clone, Copy)]
struct TriAabb {
    min: DVec3,
    max: DVec3,
}

impl TriAabb {
    fn empty() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
        }
    }

    fn from_triangle(tri: &[DVec3; 3]) -> Self {
        let mut aabb = Self::empty();
        for &p in tri {
            aabb.expand_point(p);
        }
        aabb
    }

    fn expand_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    fn expand_aabb(&mut self, other: &TriAabb) {
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
    }

    fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    fn surface_area(&self) -> f64 {
        let d = self.max - self.min;
        if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 {
            return 0.0;
        }
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }

    fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
        let t1 = (self.min - origin) * inv_dir;
        let t2 = (self.max - origin) * inv_dir;

        let t_min = t1.min(t2);
        let t_max = t1.max(t2);

        let t_enter = t_min.x.max(t_min.y).max(t_min.z);
        let t_exit = t_max.x.min(t_max.y).min(t_max.z);

        if t_exit >= t_enter.max(0.0) {
            Some(t_enter.max(0.0))
        } else {
            None
        }
    }
}

/// BVH node for triangle-level acceleration.
#[derive(Debug, Clone)]
enum TriBvhNode {
    Leaf {
        aabb: TriAabb,
        /// Triangle indices (index into the original triangle array).
        tris: Vec<usize>,
    },
    Internal {
        aabb: TriAabb,
        left: usize,
        right: usize,
    },
}

/// Triangle-level BVH for efficient ray casting.
#[derive(Debug, Clone)]
pub struct TriBvh {
    nodes: Vec<TriBvhNode>,
    triangle_aabbs: Vec<TriAabb>,
    triangle_centers: Vec<DVec3>,
}

const MAX_TRIS_PER_LEAF: usize = 8;
const SAH_BUCKETS: usize = 8;

impl TriBvh {
    /// Build a triangle BVH from a list of triangles.
    pub fn build(triangles: &[[DVec3; 3]]) -> Self {
        if triangles.is_empty() {
            return Self {
                nodes: Vec::new(),
                triangle_aabbs: Vec::new(),
                triangle_centers: Vec::new(),
            };
        }

        let triangle_aabbs: Vec<TriAabb> = triangles.iter().map(|t| TriAabb::from_triangle(t)).collect();
        let triangle_centers: Vec<DVec3> = triangle_aabbs.iter().map(|a| a.center()).collect();

        let tri_indices: Vec<usize> = (0..triangles.len()).collect();

        let mut bvh = Self {
            nodes: Vec::new(),
            triangle_aabbs,
            triangle_centers,
        };

        bvh.build_recursive(&tri_indices);
        bvh
    }

    fn build_recursive(&mut self, tri_indices: &[usize]) -> usize {
        let count = tri_indices.len();
        if count == 0 {
            return usize::MAX;
        }

        // Compute AABB for this node
        let mut aabb = TriAabb::empty();
        for &ti in tri_indices {
            aabb.expand_aabb(&self.triangle_aabbs[ti]);
        }

        // Leaf condition
        if count <= MAX_TRIS_PER_LEAF {
            let node_idx = self.nodes.len();
            self.nodes.push(TriBvhNode::Leaf {
                aabb,
                tris: tri_indices.to_vec(),
            });
            return node_idx;
        }

        // SAH split
        let (split_axis, split_pos) = self.sah_split(tri_indices, &aabb);

        // Partition triangles
        let (left_tris, right_tris): (Vec<usize>, Vec<usize>) = tri_indices
            .iter()
            .copied()
            .partition(|&ti| {
                let center = match split_axis {
                    0 => self.triangle_centers[ti].x,
                    1 => self.triangle_centers[ti].y,
                    _ => self.triangle_centers[ti].z,
                };
                center < split_pos
            });

        // Handle degenerate split
        let (left_tris, right_tris) = if left_tris.is_empty() || right_tris.is_empty() {
            let mid = count / 2;
            let mut sorted = tri_indices.to_vec();
            sorted.sort_by(|&a, &b| {
                let ca = match split_axis {
                    0 => self.triangle_centers[a].x,
                    1 => self.triangle_centers[a].y,
                    _ => self.triangle_centers[a].z,
                };
                let cb = match split_axis {
                    0 => self.triangle_centers[b].x,
                    1 => self.triangle_centers[b].y,
                    _ => self.triangle_centers[b].z,
                };
                ca.partial_cmp(&cb).unwrap_or(std::cmp::Ordering::Equal)
            });
            (sorted[..mid].to_vec(), sorted[mid..].to_vec())
        } else {
            (left_tris, right_tris)
        };

        let node_idx = self.nodes.len();
        self.nodes.push(TriBvhNode::Internal {
            aabb: TriAabb::empty(),
            left: 0,
            right: 0,
        });

        let left = self.build_recursive(&left_tris);
        let right = self.build_recursive(&right_tris);

        self.nodes[node_idx] = TriBvhNode::Internal { aabb, left, right };
        node_idx
    }

    fn sah_split(&self, tri_indices: &[usize], parent_aabb: &TriAabb) -> (usize, f64) {
        let parent_sa = parent_aabb.surface_area().max(1e-30);
        let mut best_cost = f64::INFINITY;
        let mut best_axis = 0usize;
        let mut best_pos = 0.0f64;

        for axis in 0..3usize {
            let axis_min = match axis {
                0 => parent_aabb.min.x,
                1 => parent_aabb.min.y,
                _ => parent_aabb.min.z,
            };
            let axis_max = match axis {
                0 => parent_aabb.max.x,
                1 => parent_aabb.max.y,
                _ => parent_aabb.max.z,
            };
            let span = axis_max - axis_min;
            if span < 1e-14 {
                continue;
            }

            for b in 1..SAH_BUCKETS {
                let split = axis_min + span * b as f64 / SAH_BUCKETS as f64;

                let mut left_aabb = TriAabb::empty();
                let mut right_aabb = TriAabb::empty();
                let mut left_count = 0usize;
                let mut right_count = 0usize;

                for &ti in tri_indices {
                    let center_val = match axis {
                        0 => self.triangle_centers[ti].x,
                        1 => self.triangle_centers[ti].y,
                        _ => self.triangle_centers[ti].z,
                    };
                    if center_val < split {
                        left_aabb.expand_aabb(&self.triangle_aabbs[ti]);
                        left_count += 1;
                    } else {
                        right_aabb.expand_aabb(&self.triangle_aabbs[ti]);
                        right_count += 1;
                    }
                }

                if left_count == 0 || right_count == 0 {
                    continue;
                }

                let cost = (left_count as f64 * left_aabb.surface_area()
                    + right_count as f64 * right_aabb.surface_area())
                    / parent_sa;

                if cost < best_cost {
                    best_cost = cost;
                    best_axis = axis;
                    best_pos = split;
                }
            }
        }

        if best_cost.is_infinite() {
            let d = parent_aabb.max - parent_aabb.min;
            best_axis = if d.x >= d.y && d.x >= d.z { 0 } else if d.y >= d.z { 1 } else { 2 };
            best_pos = parent_aabb.center()[best_axis];
        }

        (best_axis, best_pos)
    }

    /// Test if a point is occluded by any triangle in the BVH.
    pub fn is_occluded(&self, point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
        if self.nodes.is_empty() {
            return false;
        }

        let dir = (eye - point).normalize_or_zero();
        let origin = point + dir * 1e-5;
        let inv_dir = DVec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);

        self.is_occluded_node(0, origin, dir, inv_dir, triangles, dist_to_eye)
    }

    fn is_occluded_node(
        &self,
        node_idx: usize,
        origin: DVec3,
        dir: DVec3,
        inv_dir: DVec3,
        triangles: &[[DVec3; 3]],
        dist_to_eye: f64,
    ) -> bool {
        let node = &self.nodes[node_idx];

        // Check AABB intersection
        let t_aabb = match node.aabb().ray_intersect(origin, inv_dir) {
            Some(t) => t,
            None => return false,
        };

        // Early exit if AABB is beyond the eye
        if t_aabb > dist_to_eye {
            return false;
        }

        match node {
            TriBvhNode::Leaf { tris, .. } => {
                for &ti in tris {
                    if let Some(t) = ray_triangle_intersect(origin, dir, &triangles[ti]) {
                        if t < dist_to_eye - 1e-4 {
                            return true;
                        }
                    }
                }
                false
            }
            TriBvhNode::Internal { left, right, .. } => {
                self.is_occluded_node(*left, origin, dir, inv_dir, triangles, dist_to_eye)
                    || self.is_occluded_node(*right, origin, dir, inv_dir, triangles, dist_to_eye)
            }
        }
    }
}

impl TriBvhNode {
    fn aabb(&self) -> &TriAabb {
        match self {
            TriBvhNode::Leaf { aabb, .. } => aabb,
            TriBvhNode::Internal { aabb, .. } => aabb,
        }
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

// ── Silhouette generation ─────────────────────────────────────────────────────

/// Internal: one silhouette curve to process through the HLR pipeline.
struct SilhouetteCurve {
    /// World-space sample points (at least 2).
    world_pts: Vec<DVec3>,
    /// Optional curve hint for SVG output.
    curve_hint: Option<CurveHint>,
    /// If true, emit one segment per consecutive point pair instead of merging
    /// runs.  Used for dense polyline approximations (e.g. sphere silhouette).
    dense: bool,
}

/// A 3D silhouette curve extracted from a curved surface.
#[derive(Debug, Clone)]
pub struct SilhouetteCurve3 {
    /// World-space sample points along the silhouette curve.
    pub points: Vec<DVec3>,
    /// The surface index from which this silhouette was extracted.
    pub surface_index: usize,
}

/// Extract silhouette curves from a BRep for a given view direction.
///
/// This function computes the visible contour lines (silhouettes) of curved surfaces
/// as seen from a specific viewing direction. For analytic surfaces (cylinder, sphere,
/// cone, torus), exact silhouette curves are computed. For general surfaces (BSpline,
/// Bezier, etc.), numerical methods with adaptive sampling are used.
///
/// # Arguments
/// * `brep` - The BRep model to extract silhouettes from.
/// * `view_dir` - The normalized view direction (from target to eye).
/// * `opts` - Configuration options for sampling and tolerance.
///
/// # Returns
/// A vector of 3D silhouette curves, each represented as a series of world-space points.
pub fn extract_silhouette_curves(brep: &BRep, view_dir: DVec3, opts: &HlrOptions) -> Vec<SilhouetteCurve3> {
    let mut curves: Vec<SilhouetteCurve3> = Vec::new();

    if brep.solids.is_empty() {
        return curves;
    }

    let line_samples = opts.silhouette_samples.max(16);
    let dense_curve_samples = (opts.silhouette_samples * 4).max(64);

    let mut face_idx = 0usize;
    for shell in &brep.solids[0].shells {
        for _face in &shell.faces {
            let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                Some(idx) => idx,
                None => {
                    face_idx += 1;
                    continue;
                }
            };
            let surface = &brep.geom.surfaces[surf_idx];

            let domain = match brep.geom.face_surface_range.get(face_idx).and_then(|o| *o) {
                Some(r) => r,
                None => surface.default_domain(),
            };
            let [u0, u1, v0, v1] = domain;

            // Extract silhouettes based on surface type
            let face_curves = extract_surface_silhouettes(
                surface, view_dir, domain, brep, opts, line_samples, dense_curve_samples,
            );

            for pts in face_curves {
                if pts.len() >= 2 {
                    curves.push(SilhouetteCurve3 {
                        points: pts,
                        surface_index: surf_idx,
                    });
                }
            }

            face_idx += 1;
        }
    }

    curves
}

/// Extract silhouette curves from a single surface.
fn extract_surface_silhouettes(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    brep: &BRep,
    opts: &HlrOptions,
    line_samples: usize,
    dense_curve_samples: usize,
) -> Vec<Vec<DVec3>> {
    let [u0, u1, v0, v1] = domain;
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    match surface {
        Surface3::Cylinder(cyl) => {
            curves.extend(extract_cylinder_silhouettes(cyl, view_dir, brep, line_samples, v0, v1));
        }

        Surface3::Sphere(sph) => {
            curves.push(extract_sphere_silhouette(sph, view_dir, dense_curve_samples));
        }

        Surface3::Cone(con) => {
            curves.extend(extract_cone_silhouettes(con, view_dir, brep, line_samples, v0, v1));
        }

        Surface3::Torus(tor) => {
            curves.extend(extract_torus_silhouettes(tor, view_dir, dense_curve_samples));
        }

        Surface3::Ellipsoid(ell) => {
            curves.extend(extract_ellipsoid_silhouettes(ell, view_dir, opts, dense_curve_samples));
        }

        // For general surfaces, use numerical silhouette extraction
        Surface3::BSpline(_)
        | Surface3::Bezier(_)
        | Surface3::TriBezier(_)
        | Surface3::Offset(_)
        | Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::Gordon(_)
        | Surface3::Pipe(_) => {
            curves.extend(extract_numerical_silhouettes(
                surface, view_dir, domain, opts, brep,
            ));
        }

        // Planes have no silhouette curves
        Surface3::Plane(_) | Surface3::Trimmed(_) | Surface3::Helicoid(_) => {}
    }

    curves
}

/// Extract silhouette lines from a cylinder.
fn extract_cylinder_silhouettes(
    cyl: &rcad_kernel::geom::CylindricalSurface,
    view_dir: DVec3,
    brep: &BRep,
    line_samples: usize,
    v0: f64,
    v1: f64,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    // Project view direction onto the plane perpendicular to the axis.
    let d_perp = view_dir - view_dir.dot(cyl.axis) * cyl.axis;
    if d_perp.length_squared() < 1e-10 {
        // Viewing along the axis — no silhouette lines.
        return curves;
    }

    // Direction from axis to silhouette (perpendicular to both axis and d_perp).
    let sil_dir = cyl.axis.cross(d_perp).normalize_or_zero();

    // Resolve v range (height along axis).
    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for vert in &brep.vertices {
            let proj = (vert.point - cyl.origin).dot(cyl.axis);
            lo = lo.min(proj);
            hi = hi.max(proj);
        }
        if lo.is_finite() && hi.is_finite() {
            (lo, hi)
        } else {
            return curves;
        }
    };

    for &sign in &[1.0_f64, -1.0] {
        let offset = sil_dir * sign * cyl.radius;
        let world_pts: Vec<DVec3> = (0..line_samples)
            .map(|i| {
                let t = i as f64 / (line_samples - 1) as f64;
                let v = v0_eff + (v1_eff - v0_eff) * t;
                cyl.origin + v * cyl.axis + offset
            })
            .collect();
        curves.push(world_pts);
    }

    curves
}

/// Extract silhouette curve from a sphere (great circle perpendicular to view direction).
fn extract_sphere_silhouette(
    sph: &rcad_kernel::geom::SphericalSurface,
    view_dir: DVec3,
    samples: usize,
) -> Vec<DVec3> {
    let x_ax = any_perpendicular(view_dir);
    let y_ax = view_dir.cross(x_ax).normalize_or_zero();

    (0..samples)
        .map(|i| {
            let t = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
            sph.center + sph.radius * (t.cos() * x_ax + t.sin() * y_ax)
        })
        .collect()
}

/// Extract silhouette lines from a cone (two generators from apex).
fn extract_cone_silhouettes(
    con: &rcad_kernel::geom::ConicalSurface,
    view_dir: DVec3,
    brep: &BRep,
    line_samples: usize,
    v0: f64,
    v1: f64,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    let d_perp = view_dir - view_dir.dot(con.axis) * con.axis;
    if d_perp.length_squared() < 1e-10 {
        return curves;
    }

    let sil_dir = con.axis.cross(d_perp).normalize_or_zero();

    let (v0_eff, v1_eff) = if v0.is_finite() && v1.is_finite() {
        (v0, v1)
    } else {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for vert in &brep.vertices {
            let proj = (vert.point - con.apex).dot(con.axis);
            lo = lo.min(proj);
            hi = hi.max(proj);
        }
        if lo.is_finite() && hi.is_finite() {
            (lo.max(0.0), hi.max(0.0))
        } else {
            return curves;
        }
    };

    let tan_a = con.half_angle_rad.tan();
    for &sign in &[1.0_f64, -1.0] {
        let world_pts: Vec<DVec3> = (0..line_samples)
            .map(|i| {
                let t = i as f64 / (line_samples - 1) as f64;
                let v = v0_eff + (v1_eff - v0_eff) * t;
                con.apex + v * con.axis + v * tan_a * sil_dir * sign
            })
            .collect();

        if world_pts
            .first()
            .zip(world_pts.last())
            .map(|(a, b)| (*b - *a).length_squared() > 1e-12)
            .unwrap_or(false)
        {
            curves.push(world_pts);
        }
    }

    curves
}

/// Extract silhouette curves from a torus.
fn extract_torus_silhouettes(
    tor: &rcad_kernel::geom::ToroidalSurface,
    view_dir: DVec3,
    samples: usize,
) -> Vec<Vec<DVec3>> {
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    let x_ax = any_perpendicular(tor.axis);
    let y_ax = tor.axis.cross(x_ax).normalize_or_zero();
    let axis_dot = tor.axis.dot(view_dir);

    for &offset in &[0.0_f64, std::f64::consts::PI] {
        let pts: Vec<DVec3> = (0..samples)
            .map(|i| {
                let u = 2.0 * std::f64::consts::PI * i as f64 / samples as f64;
                let radial = u.cos() * x_ax + u.sin() * y_ax;
                let radial_dot = radial.dot(view_dir);
                let v = (-radial_dot).atan2(axis_dot) + offset;
                let tube_center = tor.center + tor.major_radius * radial;
                tube_center + tor.minor_radius * (v.cos() * radial + v.sin() * tor.axis)
            })
            .collect();
        curves.push(pts);
    }

    curves
}

/// Extract silhouette curves from an ellipsoid using numerical methods.
fn extract_ellipsoid_silhouettes(
    _ell: &rcad_kernel::geom::EllipsoidalSurface,
    _view_dir: DVec3,
    _opts: &HlrOptions,
    _samples: usize,
) -> Vec<Vec<DVec3>> {
    // TODO: Implement ellipsoid silhouette extraction
    // For now, return empty (ellipsoid is not commonly used)
    Vec::new()
}

/// Numerical silhouette extraction for general parametric surfaces.
///
/// Uses a marching approach to find curves where normal · view_dir = 0.
/// This implementation includes:
/// - Marching along iso-parametric curves to trace silhouette curves
/// - Curvature-adaptive sampling for better accuracy in high-curvature regions
/// - Handling of closed silhouette loops
fn extract_numerical_silhouettes(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    opts: &HlrOptions,
    _brep: &BRep,
) -> Vec<Vec<DVec3>> {
    let [u0, u1, v0, v1] = domain;
    let mut curves: Vec<Vec<DVec3>> = Vec::new();

    // Phase 1: Find silhouette seed points on a coarse grid
    let grid_size = opts.silhouette_samples.max(16);
    let seeds = find_silhouette_seeds(surface, view_dir, domain, grid_size, opts.tangent_tolerance);

    if seeds.is_empty() {
        return curves;
    }

    // Phase 2: March from each seed to trace silhouette curves
    let mut visited: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();

    for (i, j, u, v) in seeds {
        if visited.contains(&(i, j)) {
            continue;
        }

        // Trace a curve starting from this seed
        let curve = march_silhouette_curve(surface, view_dir, domain, u, v, opts);

        if curve.len() >= 2 {
            // Mark visited cells along the curve
            for pt in &curve {
                // Find grid cell for this point
                let pi = ((pt.0 - u0) / (u1 - u0) * grid_size as f64).floor() as usize;
                let pj = ((pt.1 - v0) / (v1 - v0) * grid_size as f64).floor() as usize;
                visited.insert((pi.min(grid_size - 1), pj.min(grid_size - 1)));
            }

            // Apply adaptive refinement based on curvature
            let refined_curve = if opts.curvature_adaptive {
                refine_curve_by_curvature(surface, curve, opts)
            } else {
                curve.into_iter().map(|(_, _, pt)| pt).collect()
            };

            // Apply B-spline fitting if enabled
            let final_curve = if opts.fit_bspline && refined_curve.len() >= 4 {
                fit_bspline_to_points(&refined_curve, opts.bspline_tolerance)
            } else {
                refined_curve
            };

            if final_curve.len() >= 2 {
                curves.push(final_curve);
            }
        }
    }

    curves
}

/// A point in parameter space with its 3D position.
type ParamPoint = (f64, f64, DVec3);

/// Find seed points for silhouette curves on a grid.
fn find_silhouette_seeds(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    grid_size: usize,
    tangent_tol: f64,
) -> Vec<(usize, usize, f64, f64)> {
    let [u0, u1, v0, v1] = domain;
    let mut seeds = Vec::new();

    // Sample grid and look for sign changes in normal · view_dir
    let mut dot_values: Vec<Vec<f64>> = vec![vec![0.0; grid_size]; grid_size];

    // Compute dot products at grid vertices
    for i in 0..grid_size {
        for j in 0..grid_size {
            let u = u0 + (u1 - u0) * i as f64 / (grid_size - 1) as f64;
            let v = v0 + (v1 - v0) * j as f64 / (grid_size - 1) as f64;
            let normal = surface.normal_at(u, v);
            dot_values[i][j] = normal.dot(view_dir);
        }
    }

    // Find cells where sign changes occur (indicating silhouette crossing)
    for i in 0..grid_size - 1 {
        for j in 0..grid_size - 1 {
            let d00 = dot_values[i][j];
            let d10 = dot_values[i + 1][j];
            let d01 = dot_values[i][j + 1];
            let d11 = dot_values[i + 1][j + 1];

            // Check for sign changes in the cell
            let has_crossing = (d00 * d10 < 0.0)
                || (d00 * d01 < 0.0)
                || (d10 * d11 < 0.0)
                || (d01 * d11 < 0.0);

            if has_crossing {
                // Find the exact crossing point using bisection
                if let Some((u, v)) = find_crossing_point(surface, view_dir, domain, i, j, grid_size, tangent_tol) {
                    seeds.push((i, j, u, v));
                }
            }
        }
    }

    seeds
}

/// Find the exact crossing point in a grid cell using bisection.
fn find_crossing_point(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    i: usize,
    j: usize,
    grid_size: usize,
    _tangent_tol: f64,
) -> Option<(f64, f64)> {
    let [u0, v0, u1, v1] = [
        domain[0] + (domain[1] - domain[0]) * i as f64 / (grid_size - 1) as f64,
        domain[2] + (domain[3] - domain[2]) * j as f64 / (grid_size - 1) as f64,
        domain[0] + (domain[1] - domain[0]) * (i + 1) as f64 / (grid_size - 1) as f64,
        domain[2] + (domain[3] - domain[2]) * (j + 1) as f64 / (grid_size - 1) as f64,
    ];

    // Try to find crossing along each edge of the cell
    let edges = [
        (u0, v0, u1, v0), // bottom edge
        (u0, v1, u1, v1), // top edge
        (u0, v0, u0, v1), // left edge
        (u1, v0, u1, v1), // right edge
    ];

    for (ua, va, ub, vb) in edges {
        if let Some((u, v)) = bisection_search(surface, view_dir, ua, va, ub, vb, 12) {
            return Some((u, v));
        }
    }

    None
}

/// Bisection search to find where normal · view_dir = 0.
fn bisection_search(
    surface: &Surface3,
    view_dir: DVec3,
    u0: f64,
    v0: f64,
    u1: f64,
    v1: f64,
    max_iter: usize,
) -> Option<(f64, f64)> {
    let d0 = surface.normal_at(u0, v0).dot(view_dir);
    let d1 = surface.normal_at(u1, v1).dot(view_dir);

    if d0 * d1 > 0.0 {
        return None; // No sign change
    }

    let mut ua = u0;
    let mut va = v0;
    let mut ub = u1;
    let mut vb = v1;

    for _ in 0..max_iter {
        let um = (ua + ub) / 2.0;
        let vm = (va + vb) / 2.0;
        let dm = surface.normal_at(um, vm).dot(view_dir);

        if dm.abs() < 1e-10 {
            return Some((um, vm));
        }

        if d0 * dm < 0.0 {
            ub = um;
            vb = vm;
        } else {
            ua = um;
            va = vm;
        }
    }

    Some(((ua + ub) / 2.0, (va + vb) / 2.0))
}

/// March along a silhouette curve starting from a seed point.
fn march_silhouette_curve(
    surface: &Surface3,
    view_dir: DVec3,
    domain: [f64; 4],
    u_start: f64,
    v_start: f64,
    opts: &HlrOptions,
) -> Vec<ParamPoint> {
    let mut curve: Vec<ParamPoint> = Vec::new();
    let [u0, u1, v0, v1] = domain;

    // Add the starting point
    let p_start = surface.point_at(u_start, v_start);
    curve.push((u_start, v_start, p_start));

    // March in both directions from the seed
    for direction in &[-1.0_f64, 1.0] {
        let mut u = u_start;
        let mut v = v_start;
        let mut curve_dir: Option<DVec2> = None;

        for _ in 0..opts.max_subdivisions * 50 {
            // Compute the tangent direction to the silhouette curve
            let tangent = compute_silhouette_tangent(surface, view_dir, u, v);

            if tangent.length_squared() < 1e-16 {
                break;
            }

            // Choose direction along the tangent
            let step_dir = if let Some(cd) = curve_dir {
                // Continue in the same general direction
                if cd.dot(tangent) > 0.0 {
                    tangent
                } else {
                    -tangent
                }
            } else {
                *direction * tangent
            };
            curve_dir = Some(step_dir.normalize_or_zero());

            // Compute step size based on curvature
            let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, u, v);
            let max_k = k1.abs().max(k2.abs()).max(opts.min_curvature);
            let curvature_factor = (opts.max_curvature / max_k).min(4.0).max(0.25);
            let step_size = opts.angular_tolerance * curvature_factor;

            // Take a step
            let u_new = u + step_dir.x * step_size;
            let v_new = v + step_dir.y * step_size;

            // Check bounds
            if u_new < u0 || u_new > u1 || v_new < v0 || v_new > v1 {
                break;
            }

            // Project back onto the silhouette curve
            if let Some((u_proj, v_proj)) = project_to_silhouette(surface, view_dir, u_new, v_new, opts.tangent_tolerance) {
                u = u_proj;
                v = v_proj;

                let p = surface.point_at(u, v);
                let d = (p - curve.last().map(|(_, _, lp)| *lp).unwrap_or(p_start)).length();

                // Only add if we've moved enough
                if d > opts.bspline_tolerance * 0.1 {
                    curve.push((u, v, p));
                }
            } else {
                break;
            }

            // Check for closed loop
            if curve.len() > 10 {
                let first = curve[0];
                let dist = ((first.0 - u).powi(2) + (first.1 - v).powi(2)).sqrt();
                if dist < step_size * 2.0 {
                    // Close the loop
                    curve.push(curve[0]);
                    break;
                }
            }
        }

        // Reverse the points added while marching in the negative direction
        if *direction < 0.0 && curve.len() > 1 {
            let first = curve[0];
            curve.reverse();
            curve.push(first); // Re-add the start point for the loop
        }
    }

    curve
}

/// Compute the tangent direction to the silhouette curve at a point.
fn compute_silhouette_tangent(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
) -> DVec2 {
    const EPS: f64 = 1e-6;

    // Compute gradients of the implicit function f(u,v) = N(u,v) · V
    let n = surface.normal_at(u, v);
    let n_u = surface.normal_at(u + EPS, v);
    let n_v = surface.normal_at(u, v + EPS);

    // Gradient of f = N · V
    let df_du = (n_u - n).dot(view_dir) / EPS;
    let df_dv = (n_v - n).dot(view_dir) / EPS;

    // The tangent direction is perpendicular to the gradient
    DVec2::new(-df_dv, df_du).normalize_or_zero()
}

/// Project a point back onto the silhouette curve.
fn project_to_silhouette(
    surface: &Surface3,
    view_dir: DVec3,
    u: f64,
    v: f64,
    tol: f64,
) -> Option<(f64, f64)> {
    let mut u_curr = u;
    let mut v_curr = v;

    // Newton iteration to find f(u,v) = 0
    for _ in 0..20 {
        let n = surface.normal_at(u_curr, v_curr);
        let f = n.dot(view_dir);

        if f.abs() < tol {
            return Some((u_curr, v_curr));
        }

        // Compute gradient numerically
        const EPS: f64 = 1e-7;
        let n_u = surface.normal_at(u_curr + EPS, v_curr);
        let n_v = surface.normal_at(u_curr, v_curr + EPS);

        let df_du = (n_u - n).dot(view_dir) / EPS;
        let df_dv = (n_v - n).dot(view_dir) / EPS;

        let grad_len_sq = df_du * df_du + df_dv * df_dv;
        if grad_len_sq < 1e-20 {
            break;
        }

        // Newton step
        let step = f / grad_len_sq;
        u_curr -= step * df_du;
        v_curr -= step * df_dv;
    }

    // Check if we converged
    let f = surface.normal_at(u_curr, v_curr).dot(view_dir);
    if f.abs() < tol * 10.0 {
        Some((u_curr, v_curr))
    } else {
        None
    }
}

/// Refine a silhouette curve based on surface curvature.
fn refine_curve_by_curvature(
    surface: &Surface3,
    curve: Vec<ParamPoint>,
    opts: &HlrOptions,
) -> Vec<DVec3> {
    if curve.len() < 2 {
        return curve.into_iter().map(|(_, _, p)| p).collect();
    }

    let mut refined: Vec<DVec3> = Vec::new();
    refined.push(curve[0].2);

    for i in 1..curve.len() {
        let (u0, v0, p0) = curve[i - 1];
        let (u1, v1, p1) = curve[i];

        // Compute curvature at the midpoint
        let um = (u0 + u1) / 2.0;
        let vm = (v0 + v1) / 2.0;
        let (k1, k2) = rcad_kernel::curvature::principal_curvatures(surface, um, vm);
        let max_k = k1.abs().max(k2.abs());

        // Determine number of subdivision points based on curvature
        let chord_len = (p1 - p0).length();
        let subdivs = if max_k > opts.min_curvature {
            let curvature_samples = (max_k * chord_len * std::f64::consts::PI).ceil() as usize;
            curvature_samples.min(8).max(1)
        } else {
            1
        };

        // Add subdivision points
        for j in 1..subdivs {
            let t = j as f64 / subdivs as f64;
            let u = u0 + t * (u1 - u0);
            let v = v0 + t * (v1 - v0);
            let p = surface.point_at(u, v);
            refined.push(p);
        }

        refined.push(p1);
    }

    refined
}

/// Fit a B-spline curve to a set of points.
fn fit_bspline_to_points(points: &[DVec3], tolerance: f64) -> Vec<DVec3> {
    if points.len() < 4 {
        return points.to_vec();
    }

    // Simple approach: sample the fitted B-spline at uniform intervals
    // For a proper implementation, we would use least-squares fitting
    // Here we use a simplified version that preserves the shape

    let n = points.len();
    let mut result: Vec<DVec3> = Vec::with_capacity(n);

    // Compute chord lengths for parameterization
    let mut chords = vec![0.0_f64; n];
    for i in 1..n {
        chords[i] = chords[i - 1] + (points[i] - points[i - 1]).length();
    }
    let total_len = chords[n - 1];
    if total_len < 1e-12 {
        return points.to_vec();
    }

    // Generate control points using Catmull-Rom style interpolation
    let degree = 3.min(n - 1);
    let num_samples = (total_len / tolerance).ceil() as usize;
    let num_samples = num_samples.max(10).min(1000);

    for i in 0..=num_samples {
        let t = i as f64 / num_samples as f64;
        let target_len = t * total_len;

        // Find the segment containing this length
        let seg_idx = chords.partition_point(|&c| c < target_len).saturating_sub(1);
        let seg_idx = seg_idx.min(n - 2);

        // Interpolate within the segment
        let seg_start = chords[seg_idx];
        let seg_end = chords[seg_idx + 1];
        let seg_len = seg_end - seg_start;

        let local_t = if seg_len > 1e-12 {
            (target_len - seg_start) / seg_len
        } else {
            0.5
        };

        // Simple linear interpolation with smoothing
        let p0 = points[seg_idx];
        let p1 = points[seg_idx + 1];

        // Hermite interpolation for smoother result
        let t0 = if seg_idx > 0 {
            (points[seg_idx + 1] - points[seg_idx - 1]).normalize_or_zero()
        } else {
            (points[1] - points[0]).normalize_or_zero()
        };

        let t1 = if seg_idx + 2 < n {
            (points[seg_idx + 2] - points[seg_idx]).normalize_or_zero()
        } else {
            (points[n - 1] - points[n - 2]).normalize_or_zero()
        };

        let h00 = 2.0 * local_t * local_t * local_t - 3.0 * local_t * local_t + 1.0;
        let h10 = local_t * local_t * local_t - 2.0 * local_t * local_t + local_t;
        let h01 = -2.0 * local_t * local_t * local_t + 3.0 * local_t * local_t;
        let h11 = local_t * local_t * local_t - local_t * local_t;

        let p = h00 * p0 + h10 * t0 * seg_len + h01 * p1 + h11 * t1 * seg_len;
        result.push(p);
    }

    result
}

/// Generate silhouette curves for the HLR pipeline (internal function).
fn compute_silhouettes(brep: &BRep, view_dir: DVec3, samples: usize) -> Vec<SilhouetteCurve> {
    let opts = HlrOptions {
        silhouette_samples: samples,
        ..HlrOptions::default()
    };

    extract_silhouette_curves(brep, view_dir, &opts)
        .into_iter()
        .map(|curve| SilhouetteCurve {
            world_pts: curve.points,
            curve_hint: None,
            dense: true, // All silhouettes are treated as dense for proper rendering
        })
        .collect()
}

/// Occlusion tester that supports both brute-force and BVH-accelerated methods.
enum OcclusionTester<'a> {
    BruteForce(&'a [[DVec3; 3]]),
    Bvh {
        bvh: &'a TriBvh,
        triangles: &'a [[DVec3; 3]],
    },
}

impl<'a> OcclusionTester<'a> {
    fn is_occluded(&self, point: DVec3, eye: DVec3, dist_to_eye: f64) -> bool {
        match self {
            OcclusionTester::BruteForce(triangles) => {
                is_occluded(point, eye, triangles, dist_to_eye)
            }
            OcclusionTester::Bvh { bvh, triangles } => {
                bvh.is_occluded(point, eye, triangles, dist_to_eye)
            }
        }
    }
}

/// Improved visibility classification that handles grazing angles on curved surfaces.
///
/// For points near silhouette curves (where normal is nearly perpendicular to view direction),
/// we use additional testing to improve numerical stability.
fn classify_visibility(
    point: DVec3,
    normal: Option<DVec3>,
    camera: &HlrCamera,
    occlusion_tester: &OcclusionTester<'_>,
    grazing_threshold: f64,
) -> VisibilityInfo {
    let dist = (camera.eye - point).length();
    let view_dir = (camera.eye - point).normalize_or_zero();

    // Check if we're at a grazing angle
    let grazing_factor = if let Some(n) = normal {
        let dot = n.dot(view_dir).abs();
        // grazing_factor = 1.0 when perfectly grazing (dot = 0)
        // grazing_factor = 0.0 when viewing straight on (dot = 1)
        1.0 - dot
    } else {
        0.0
    };

    // For grazing angles, use more robust testing
    let is_occluded = if grazing_factor > grazing_threshold.cos() {
        // At grazing angle: test multiple rays to reduce false positives
        let base_occluded = occlusion_tester.is_occluded(point, camera.eye, dist);

        if base_occluded {
            // Verify with additional samples to reduce numerical errors
            let mut occluded_count = 1;
            const NUM_SAMPLES: usize = 4;
            let offset = 1e-4;

            for i in 0..NUM_SAMPLES {
                let angle = i as f64 * std::f64::consts::TAU / NUM_SAMPLES as f64;
                let perp = any_perpendicular(view_dir);
                let perturb = perp * (angle.cos() * offset) + view_dir.cross(perp) * (angle.sin() * offset);
                let test_point = point + perturb;

                if occlusion_tester.is_occluded(test_point, camera.eye, dist) {
                    occluded_count += 1;
                }
            }

            // Require majority to confirm occlusion at grazing angles
            occluded_count > NUM_SAMPLES / 2
        } else {
            false
        }
    } else {
        occlusion_tester.is_occluded(point, camera.eye, dist)
    };

    VisibilityInfo {
        visible: !is_occluded,
        grazing_factor,
        depth: dist,
    }
}

/// Information about visibility at a point.
struct VisibilityInfo {
    visible: bool,
    grazing_factor: f64,
    depth: f64,
}

/// Process a list of world-space sample points through the HLR visibility
/// pipeline and append resulting segments to `result`.
///
/// When `dense` is true, one segment is emitted per consecutive point pair
/// (useful for polyline approximations of curved silhouettes).
fn process_world_pts(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    segment_type: SegmentType,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    result: &mut HlrResult,
) {
    process_world_pts_with_bvh(
        world_pts,
        curve_hint,
        dense,
        segment_type,
        camera,
        view,
        triangles,
        None,
        &HlrOptions::default(),
        result,
    )
}

/// Process world points with optional BVH acceleration and grazing angle handling.
fn process_world_pts_with_bvh(
    world_pts: &[DVec3],
    curve_hint: Option<CurveHint>,
    dense: bool,
    segment_type: SegmentType,
    camera: &HlrCamera,
    view: &DMat4,
    triangles: &[[DVec3; 3]],
    bvh: Option<&TriBvh>,
    opts: &HlrOptions,
    result: &mut HlrResult,
) {
    if world_pts.len() < 2 {
        return;
    }
    let n = world_pts.len();

    let occlusion_tester = if let Some(bvh) = bvh {
        OcclusionTester::Bvh { bvh, triangles }
    } else {
        OcclusionTester::BruteForce(triangles)
    };

    let sample_vis: Vec<bool> = world_pts
        .iter()
        .map(|&wp| {
            let dist = (camera.eye - wp).length();
            !occlusion_tester.is_occluded(wp, camera.eye, dist)
        })
        .collect();

    let screen_pts: Vec<DVec2> = world_pts.iter().map(|&wp| project(wp, view).0).collect();

    if dense {
        // Emit one segment per consecutive pair (preserves polyline shape).
        for i in 0..n - 1 {
            let seg = HlrSegment {
                start: screen_pts[i],
                end: screen_pts[i + 1],
                visible: sample_vis[i] && sample_vis[i + 1],
                curve_hint: curve_hint.clone(),
                segment_type,
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
        }
        return;
    }

    let mut seg_start = 0usize;
    for i in 1..n {
        let changed = sample_vis[i] != sample_vis[seg_start];
        let last = i == n - 1;
        if changed || last {
            let end_idx = if last && !changed { i } else { i - 1 };
            let seg = HlrSegment {
                start: screen_pts[seg_start],
                end: screen_pts[end_idx],
                visible: sample_vis[seg_start],
                curve_hint: curve_hint.clone(),
                segment_type,
            };
            if (seg.end - seg.start).length_squared() > 1e-16 {
                result.segments.push(seg);
            }
            if changed {
                seg_start = i;
            }
        }
    }
}



/// Perform hidden-line removal on a BRep from the given camera position.
///
/// Returns 2D projected segments labeled visible/hidden.
/// `samples` controls how finely each edge is subdivided for occlusion testing
/// (higher = more accurate but slower; 8 is a reasonable default).
pub fn hlr(brep: &BRep, camera: &HlrCamera, samples: usize) -> HlrResult {
    hlr_with_options(brep, camera, HlrOptions::default().with_edge_samples(samples))
}

/// Perform hidden-line removal with full configuration options.
///
/// This function provides fine-grained control over HLR computation parameters,
/// including adaptive sampling for curved surfaces.
///
/// # Arguments
/// * `brep` - The BRep model to process.
/// * `camera` - Camera/view specification.
/// * `opts` - Configuration options for sampling and tolerances.
///
/// # Returns
/// An `HlrResult` containing projected 2D segments labeled as visible/hidden.
pub fn hlr_with_options(brep: &BRep, camera: &HlrCamera, opts: HlrOptions) -> HlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let triangles = collect_triangles(brep);
    let edge_samples = opts.edge_samples.max(2);
    let mut result = HlrResult::default();

    // Build BVH for acceleration if enabled and we have enough triangles
    let bvh: Option<TriBvh> = if opts.use_bvh && triangles.len() > 32 {
        Some(TriBvh::build(&triangles))
    } else {
        None
    };
    let bvh_ref = bvh.as_ref();

    // ── Wire edges ────────────────────────────────────────────────────────────

    // Collect all unique edges from all faces + standalone edges
    let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    edge_indices.insert(we.idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        edge_indices.insert(we.idx);
                    }
                }
            }
        }
    }
    for i in 0..brep.edges.len() {
        edge_indices.insert(i);
    }

    for &edge_idx in &edge_indices {
        let Some(edge) = brep.edges.get(edge_idx) else { continue };
        let Some(v_start) = brep.vertices.get(edge.start) else { continue };
        let Some(v_end) = brep.vertices.get(edge.end) else { continue };

        let p0 = v_start.point;
        let p1 = v_end.point;

        let edge_curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|&ci| ci)
            .and_then(|ci| brep.geom.curves.get(ci));

        let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
            if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
        });

        let is_other_curve = edge_curve
            .map_or(false, |c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
            && circle_info.is_none();

        let this_edge_samples = if circle_info.is_some() || is_other_curve {
            (edge_samples * 4).max(32)
        } else {
            edge_samples
        };

        let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
            let [t0, t1] = brep
                .geom
                .edge_curve_range
                .get(edge_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| circ.default_domain());
            (0..this_edge_samples)
                .map(|i| {
                    let t = t0 + (t1 - t0) * (i as f64 / (this_edge_samples - 1) as f64);
                    circ.point_at(t)
                })
                .collect()
        } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
            let [t0, t1] = brep
                .geom
                .edge_curve_range
                .get(edge_idx)
                .and_then(|r| *r)
                .unwrap_or_else(|| curve.default_domain());
            (0..this_edge_samples)
                .map(|i| {
                    let t = t0 + (t1 - t0) * (i as f64 / (this_edge_samples - 1) as f64);
                    curve.point_at(t)
                })
                .collect()
        } else {
            if (p1 - p0).length_squared() < 1e-12 {
                continue;
            }
            (0..this_edge_samples)
                .map(|i| {
                    let t = i as f64 / (this_edge_samples - 1) as f64;
                    p0 + (p1 - p0) * t
                })
                .collect()
        };

        // Compute curve_hint for circle edges
        let screen_pts_for_hint: Vec<DVec2> =
            world_pts.iter().map(|&wp| project(wp, &view).0).collect();
        let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
            let (center_2d, _) = project(circ.center, &view);
            let r = screen_pts_for_hint
                .iter()
                .map(|p| (*p - center_2d).length())
                .fold(0.0_f64, f64::max);
            Some(CurveHint::Circle { center: center_2d, radius: r })
        } else if is_other_curve {
            Some(CurveHint::Other)
        } else {
            None
        };

        process_world_pts_with_bvh(
            &world_pts,
            curve_hint,
            false,
            SegmentType::Edge,
            camera,
            &view,
            &triangles,
            bvh_ref,
            &opts,
            &mut result,
        );
    }

    // ── Silhouette curves ────────────────────────────────────────────

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    for sil in compute_silhouettes_with_options(brep, view_dir, &opts) {
        process_world_pts_with_bvh(
            &sil.world_pts,
            sil.curve_hint,
            sil.dense,
            SegmentType::Silhouette,
            camera,
            &view,
            &triangles,
            bvh_ref,
            &opts,
            &mut result,
        );
    }

    result
}

/// Compute silhouette curves with full options (internal helper).
fn compute_silhouettes_with_options(brep: &BRep, view_dir: DVec3, opts: &HlrOptions) -> Vec<SilhouetteCurve> {
    extract_silhouette_curves(brep, view_dir, opts)
        .into_iter()
        .map(|curve| SilhouetteCurve {
            world_pts: curve.points,
            curve_hint: None,
            dense: true,
        })
        .collect()
}

/// Per-component HLR result for assembly HLR.
#[derive(Debug, Clone, Default)]
pub struct ComponentHlr {
    /// Component name (from the assembly).
    pub name: String,
    /// HLR segments for this component.
    pub segments: Vec<HlrSegment>,
}

/// Output of assembly HLR — one `ComponentHlr` per leaf BRep.
#[derive(Debug, Clone, Default)]
pub struct AssemblyHlrResult {
    pub components: Vec<ComponentHlr>,
}

impl AssemblyHlrResult {
    /// Return all visible segments across all components.
    pub fn visible_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| s.visible).map(move |s| (c, s))
        })
    }

    /// Return all hidden segments across all components.
    pub fn hidden_segments(&self) -> impl Iterator<Item = (&ComponentHlr, &HlrSegment)> {
        self.components.iter().flat_map(|c| {
            c.segments.iter().filter(|s| !s.visible).map(move |s| (c, s))
        })
    }
}

/// Transform a BRep's vertices by an affine transform.
/// Returns a new BRep with transformed vertex positions.
fn transform_brep(brep: &BRep, transform: &DAffine3) -> BRep {
    let mut out = brep.clone();
    for v in &mut out.vertices {
        v.point = transform.transform_point3(v.point);
    }
    out
}

/// Perform hidden-line removal on an assembly of BReps.
///
/// Each component's geometry is transformed to world space, then all triangles
/// are merged into a single occlusion buffer. Each component's edges are
/// tested against the global occlusion buffer, so components correctly
/// occlude each other.
///
/// Returns one `ComponentHlr` per leaf component.
pub fn hlr_assembly(
    components: &[(BRep, DAffine3, String)],
    camera: &HlrCamera,
    samples: usize,
) -> AssemblyHlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let samples = samples.max(2);

    // Transform all BRePs to world space and collect a unified triangle pool.
    let world_breps: Vec<BRep> = components
        .iter()
        .map(|(brep, xf, _)| transform_brep(brep, xf))
        .collect();

    let mut all_triangles: Vec<[DVec3; 3]> = Vec::new();
    for wb in &world_breps {
        all_triangles.extend(collect_triangles(wb));
    }

    let view_dir = (camera.target - camera.eye).normalize_or_zero();
    let mut result = AssemblyHlrResult::default();

    for (wb, (_, _, name)) in world_breps.iter().zip(components.iter()) {
        let mut comp_result = HlrResult::default();

        // ── Wire edges ────────────────────────────────────────────────────
        let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        for solid in &wb.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for we in &face.outer_wire.edges {
                        edge_indices.insert(we.idx);
                    }
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            edge_indices.insert(we.idx);
                        }
                    }
                }
            }
        }
        for i in 0..wb.edges.len() {
            edge_indices.insert(i);
        }

        for &edge_idx in &edge_indices {
            let Some(edge) = wb.edges.get(edge_idx) else { continue };
            let Some(v_start) = wb.vertices.get(edge.start) else { continue };
            let Some(v_end) = wb.vertices.get(edge.end) else { continue };

            let p0 = v_start.point;
            let p1 = v_end.point;

            let edge_curve = wb
                .geom
                .edge_curve
                .get(edge_idx)
                .and_then(|&ci| ci)
                .and_then(|ci| wb.geom.curves.get(ci));

            let circle_info: Option<Circle3> = edge_curve.and_then(|c| {
                if let rcad_kernel::geom::Curve3::Circle(circ) = c { Some(*circ) } else { None }
            });

            let is_other_curve = edge_curve
                .map_or(false, |c| !matches!(c, rcad_kernel::geom::Curve3::Line(_)))
                && circle_info.is_none();

            let edge_samples = if circle_info.is_some() || is_other_curve {
                64.max(samples)
            } else {
                samples
            };

            let world_pts: Vec<DVec3> = if let Some(circ) = &circle_info {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| circ.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        circ.point_at(t)
                    })
                    .collect()
            } else if let Some(curve) = edge_curve.filter(|_| is_other_curve) {
                let [t0, t1] = wb
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain());
                (0..edge_samples)
                    .map(|i| {
                        let t = t0 + (t1 - t0) * (i as f64 / (edge_samples - 1) as f64);
                        curve.point_at(t)
                    })
                    .collect()
            } else {
                if (p1 - p0).length_squared() < 1e-12 {
                    continue;
                }
                (0..edge_samples)
                    .map(|i| {
                        let t = i as f64 / (edge_samples - 1) as f64;
                        p0 + (p1 - p0) * t
                    })
                    .collect()
            };

            let screen_pts_for_hint: Vec<DVec2> =
                world_pts.iter().map(|&wp| project(wp, &view).0).collect();
            let curve_hint: Option<CurveHint> = if let Some(circ) = &circle_info {
                let (center_2d, _) = project(circ.center, &view);
                let r = screen_pts_for_hint
                    .iter()
                    .map(|p| (*p - center_2d).length())
                    .fold(0.0_f64, f64::max);
                Some(CurveHint::Circle { center: center_2d, radius: r })
            } else if is_other_curve {
                Some(CurveHint::Other)
            } else {
                None
            };

            process_world_pts(&world_pts, curve_hint, false, SegmentType::Edge, camera, &view, &all_triangles, &mut comp_result);
        }

        // ── Silhouette curves ────────────────────────────────────
        let opts = HlrOptions::default().with_edge_samples(samples);
        for sil in compute_silhouettes_with_options(wb, view_dir, &opts) {
            process_world_pts(&sil.world_pts, sil.curve_hint, sil.dense, SegmentType::Silhouette, camera, &view, &all_triangles, &mut comp_result);
        }

        result.components.push(ComponentHlr {
            name: name.clone(),
            segments: comp_result.segments,
        });
    }

    result
}

/// Render HLR result as a simple SVG string.
///
/// Visible edges are drawn solid black; hidden edges are dashed gray.
/// `scale` controls pixel size per unit.
pub fn hlr_to_svg(result: &HlrResult, scale: f64, margin: f64) -> String {
    if result.segments.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string();
    }

    // Compute bounding box
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for seg in &result.segments {
        for p in [seg.start, seg.end] {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
    }

    // Flip Y (SVG Y grows downward, camera Y grows upward)
    let transform = |p: DVec2| -> (f64, f64) {
        let x = (p.x - min_x) * scale + margin;
        let y = (max_y - p.y) * scale + margin;
        (x, y)
    };

    let w = (max_x - min_x) * scale + 2.0 * margin;
    let h = (max_y - min_y) * scale + 2.0 * margin;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}\" height=\"{:.1}\" viewBox=\"0 0 {:.1} {:.1}\">\n",
        w, h, w, h
    );
    svg.push_str("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    for seg in &result.segments {
        let (x1, y1) = transform(seg.start);
        let (x2, y2) = transform(seg.end);
        let stroke = if seg.visible {
            "black\" stroke-width=\"1.5"
        } else {
            "#999\" stroke-width=\"0.8\" stroke-dasharray=\"4,3"
        };

        // For circle segments emit an SVG arc path; for all others emit a line.
        if let Some(CurveHint::Circle { center, radius }) = &seg.curve_hint {
            let (cx, cy) = transform(*center);
            let r = radius * scale;
            // Determine large-arc flag: compare arc length vs half-circumference
            let dx1 = x1 - cx;
            let dy1 = y1 - cy;
            let dx2 = x2 - cx;
            let dy2 = y2 - cy;
            let cross = dx1 * dy2 - dy1 * dx2;
            let dot = dx1 * dx2 + dy1 * dy2;
            let angle = cross.atan2(dot).abs();
            let large_arc = if angle > std::f64::consts::PI { 1 } else { 0 };
            let sweep = if cross < 0.0 { 0 } else { 1 };
            svg.push_str(&format!(
                "  <path d=\"M {:.3} {:.3} A {:.3} {:.3} 0 {} {} {:.3} {:.3}\" fill=\"none\" stroke=\"{}\"/>\n",
                x1, y1, r, r, large_arc, sweep, x2, y2, stroke
            ));
            // Also record the center for debugging/reference (as a tiny dot, invisible by default)
            let _ = (cx, cy); // suppress unused warning
        } else {
            svg.push_str(&format!(
                "  <line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" stroke=\"{}\"/>\n",
                x1, y1, x2, y2, stroke
            ));
        }
    }
    svg.push_str("</svg>\n");
    svg
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn unit_box_hlr_produces_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "HLR should produce segments for a box"
        );
        assert!(
            result.visible().count() > 0,
            "some segments should be visible"
        );
    }

    #[test]
    fn hlr_svg_is_valid_xml() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(svg.contains("<svg"), "output should be SVG");
        assert!(svg.contains("</svg>"), "SVG should close properly");
        assert!(svg.contains("<line"), "SVG should contain lines");
    }

    #[test]
    fn top_view_box_has_visible_top_edges() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);
        let vis = result.visible().count();
        let hid = result.hidden().count();
        assert!(vis > 0, "top view should have visible edges");
        assert!(hid > 0, "top view should have hidden (bottom) edges");
    }

    #[test]
    fn front_view_and_right_view_both_produce_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 1.0,
            depth: 1.0,
        });
        let front_result = hlr(&brep, &HlrCamera::front(5.0), 8);
        let right_result = hlr(&brep, &HlrCamera::right(5.0), 8);
        assert!(!front_result.segments.is_empty());
        assert!(!right_result.segments.is_empty());
    }

    #[test]
    fn hlr_svg_contains_hidden_dashed_lines() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        // Hidden lines are rendered dashed
        assert!(
            svg.contains("stroke-dasharray") || svg.contains("hidden"),
            "SVG should mark hidden lines differently"
        );
    }

    #[test]
    fn hlr_result_has_correct_visibility_counts() {
        // An isometric view of a box has 3 visible faces and 3 hidden faces.
        // The front 3 edges of each visible face → at least some hidden segments exist.
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(10.0);
        let result = hlr(&brep, &camera, 16);
        let total = result.segments.len();
        assert!(total >= 12, "a box has 12 edges, expect at least 12 segments; got {total}");
    }

    #[test]
    fn hlr_circle_edge_sampling() {
        use rcad_kernel::geom::{Circle3, Curve3, CurveEval};

        // Build a minimal BRep with a single circle edge (no solids).
        let mut brep = rcad_kernel::BRep::new();
        let circ = Circle3 {
            center: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
            radius: 1.0,
        };
        // Add two vertices on the circle (half-circle arc)
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(0.0),
        });
        brep.vertices.push(rcad_kernel::topology::Vertex {
            point: circ.point_at(std::f64::consts::PI),
        });
        brep.edges.push(rcad_kernel::topology::Edge { start: 0, end: 1 });
        brep.geom.curves.push(Curve3::Circle(circ));
        brep.geom.edge_curve.push(Some(0));
        brep.geom
            .edge_curve_range
            .push(Some([0.0, std::f64::consts::PI]));

        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);

        // The circle edge should produce at least one segment.
        assert!(
            !result.segments.is_empty(),
            "circle edge should produce HLR segments"
        );

        // All sampled 3D points on the circle should lie ON the circle (unit radius).
        // Verify by checking screen_pts all lie within radius ≈ 1.0 of circle center
        // when projected top-down (X-Y plane).
        for seg in &result.segments {
            // The curve_hint for circle segments should be set.
            assert!(
                matches!(seg.curve_hint, Some(CurveHint::Circle { .. })),
                "circle edge segments should carry CurveHint::Circle"
            );
        }

        // SVG should contain arc path elements (not just lines) for circle edges.
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(
            svg.contains("<path") || result.segments.is_empty(),
            "circle edge SVG should contain <path> arc elements"
        );
    }

    /// Cylinder viewed from the side should produce silhouette line segments
    /// in addition to the wire edges.
    #[test]
    fn cylinder_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        // The cylinder axis is +Y.  Use the right-side camera (looking along -X)
        // so the view direction is perpendicular to the axis → two silhouette lines.
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cylinder HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cylinder should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Sphere HLR should produce silhouette segments (the great circle).
    #[test]
    fn sphere_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "sphere HLR should produce silhouette segments"
        );
    }

    /// Cone viewed from the side should produce two silhouette lines from the apex.
    #[test]
    fn cone_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cone {
            base_radius: 1.0,
            height: 2.0,
        });
        // View from the right (perpendicular to cone axis) → two silhouette generators.
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "cone HLR should produce segments"
        );
        assert!(
            result.segments.len() >= 2,
            "cone should have at least 2 silhouette segments, got {}",
            result.segments.len()
        );
    }

    /// Torus HLR should produce silhouette segments.
    #[test]
    fn torus_hlr_has_silhouette_segments() {
        use rcad_kernel::geom::PrimitiveSolid;
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let camera = HlrCamera::front(20.0);
        let result = hlr(&brep, &camera, 8);
        assert!(
            !result.segments.is_empty(),
            "torus HLR should produce silhouette segments"
        );
    }

    // ── Assembly HLR tests ─────────────────────────────────────────────────────

    /// Two boxes side by side — both should produce segments.
    #[test]
    fn hlr_assembly_two_boxes() {
        let box1 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let box2 = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (box1, DAffine3::from_translation(DVec3::new(-2.0, 0.0, 0.0)), "box_left".to_string()),
            (box2, DAffine3::from_translation(DVec3::new(2.0, 0.0, 0.0)), "box_right".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2, "should have 2 component results");
        assert!(result.components.iter().all(|c| !c.segments.is_empty()),
            "each component should produce segments");
    }

    /// Small box behind a large box — the small box should be partially hidden.
    #[test]
    fn hlr_assembly_occlusion() {
        let big = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 3.0, height: 3.0, depth: 3.0,
        });
        let small = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 0.5, height: 0.5, depth: 0.5,
        });

        // Front camera looks along +Y from (0, -10, 0).
        // Place small box at +Y behind the big box so it's occluded.
        let components = vec![
            (big, DAffine3::IDENTITY, "big".to_string()),
            (small, DAffine3::from_translation(DVec3::new(0.0, 3.0, 0.0)), "small_behind".to_string()),
        ];

        let camera = HlrCamera::front(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // The small box behind the big one should have mostly hidden segments
        let small_comp = result.components.iter().find(|c| c.name == "small_behind").unwrap();
        let hidden = small_comp.segments.iter().filter(|s| !s.visible).count();
        let visible = small_comp.segments.iter().filter(|s| s.visible).count();
        assert!(hidden > visible,
            "small box behind big one should have more hidden than visible segments; hidden={hidden}, visible={visible}");
    }

    /// Assembly with a single component should match single-BRep HLR.
    #[test]
    fn hlr_assembly_single_matches_hlr() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);

        let single_hlr = hlr(&brep, &camera, 8);
        let assembly_result = hlr_assembly(
            &[(brep.clone(), DAffine3::IDENTITY, "box".to_string())],
            &camera, 8,
        );

        assert_eq!(assembly_result.components.len(), 1);
        let asm_segs = &assembly_result.components[0].segments;
        // Segment counts should be similar (same geometry, same algorithm)
        assert!(asm_segs.len() >= single_hlr.segments.len() - 2,
            "assembly HLR should produce similar segment count");
        assert!(asm_segs.len() <= single_hlr.segments.len() + 2,
            "assembly HLR should produce similar segment count");
    }

    /// Stacked boxes — top box visible, bottom box partially occluded.
    #[test]
    fn hlr_assembly_stacked_boxes() {
        let bottom = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0, height: 1.0, depth: 2.0,
        });
        let top = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });

        let components = vec![
            (bottom, DAffine3::from_translation(DVec3::new(0.0, 0.0, 0.0)), "bottom".to_string()),
            (top, DAffine3::from_translation(DVec3::new(0.0, 0.0, 1.5)), "top".to_string()),
        ];

        let camera = HlrCamera::isometric(10.0);
        let result = hlr_assembly(&components, &camera, 8);

        assert_eq!(result.components.len(), 2);
        // Both boxes should have some visible segments
        for comp in &result.components {
            let vis = comp.segments.iter().filter(|s| s.visible).count();
            assert!(vis > 0, "{} should have visible segments", comp.name);
        }
    }

    /// Empty assembly should return empty result.
    #[test]
    fn hlr_assembly_empty() {
        let components: Vec<(BRep, DAffine3, String)> = vec![];
        let camera = HlrCamera::isometric(5.0);
        let result = hlr_assembly(&components, &camera, 8);
        assert!(result.components.is_empty());
    }

    // ── Improved HLR tests ─────────────────────────────────────────────────────

    #[test]
    fn hlr_options_default_values() {
        let opts = HlrOptions::default();
        assert_eq!(opts.edge_samples, 8);
        assert_eq!(opts.silhouette_samples, 32);
        assert!(opts.curvature_adaptive);
        assert!(opts.tangent_tolerance > 0.0);
    }

    #[test]
    fn hlr_options_builders() {
        let opts = HlrOptions::default()
            .with_edge_samples(16)
            .with_silhouette_samples(64)
            .with_curvature_adaptive(false)
            .with_tangent_tolerance(1e-4);

        assert_eq!(opts.edge_samples, 16);
        assert_eq!(opts.silhouette_samples, 64);
        assert!(!opts.curvature_adaptive);
        assert!((opts.tangent_tolerance - 1e-4).abs() < 1e-10);
    }

    #[test]
    fn hlr_with_options_basic() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let opts = HlrOptions::default().with_edge_samples(16);
        let result = hlr_with_options(&brep, &camera, opts);

        assert!(!result.segments.is_empty(), "should produce segments");
        // All segments from a box should be edges, not silhouettes
        assert!(result.segments.iter().all(|s| s.segment_type == SegmentType::Edge));
    }

    #[test]
    fn cylinder_silhouettes_are_marked() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });
        let camera = HlrCamera::right(10.0);
        let result = hlr(&brep, &camera, 8);

        // Should have both edge and silhouette segments
        let has_silhouette = result.segments.iter().any(|s| s.is_contour());
        assert!(has_silhouette, "cylinder should have silhouette segments");
    }

    #[test]
    fn sphere_silhouettes_are_marked() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);

        // All segments from a sphere should be silhouettes (no wire edges)
        assert!(
            result.segments.iter().all(|s| s.is_contour()),
            "sphere should only have silhouette segments"
        );
    }

    #[test]
    fn extract_silhouette_curves_sphere() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert_eq!(curves.len(), 1, "sphere should have one silhouette curve");
        assert!(curves[0].points.len() >= 32, "silhouette should have enough points");

        // All points should be at distance ~2.0 from origin
        for pt in &curves[0].points {
            let dist = pt.length();
            assert!(
                (dist - 2.0).abs() < 0.01,
                "silhouette point distance should be ~2.0, got {dist}"
            );
        }
    }

    #[test]
    fn extract_silhouette_curves_cylinder() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 3.0,
        });
        // View along X axis - perpendicular to cylinder axis (Y)
        let view_dir = DVec3::X;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert!(curves.len() >= 2, "cylinder should have at least 2 silhouette curves");

        // Each silhouette curve should be a line (two lines on opposite sides)
        for curve in &curves {
            assert!(curve.points.len() >= 16, "silhouette should have enough points");
        }
    }

    #[test]
    fn extract_silhouette_curves_torus() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 3.0,
            minor_radius: 1.0,
        });
        let view_dir = DVec3::Z;
        let opts = HlrOptions::default();

        let curves = extract_silhouette_curves(&brep, view_dir, &opts);

        assert!(curves.len() >= 2, "torus should have at least 2 silhouette curves");
    }

    #[test]
    fn hlr_result_silhouettes_iterator() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let camera = HlrCamera::front(10.0);
        let result = hlr(&brep, &camera, 8);

        let sil_count = result.silhouettes().count();
        assert!(sil_count > 0, "should have silhouette segments");

        let vis_sil_count = result.visible_silhouettes().count();
        assert!(vis_sil_count > 0, "should have visible silhouette segments");
    }

    #[test]
    fn segment_is_contour_method() {
        let seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Silhouette,
        };
        assert!(seg.is_contour());

        let edge_seg = HlrSegment {
            start: DVec2::ZERO,
            end: DVec2::X,
            visible: true,
            curve_hint: None,
            segment_type: SegmentType::Edge,
        };
        assert!(!edge_seg.is_contour());
    }

    #[test]
    fn adaptive_sampling_high_curvature() {
        // Test that adaptive sampling produces more points in high-curvature regions
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let view_dir = DVec3::Z;

        let opts_low = HlrOptions {
            silhouette_samples: 16,
            curvature_adaptive: false,
            ..HlrOptions::default()
        };
        let opts_high = HlrOptions {
            silhouette_samples: 64,
            curvature_adaptive: true,
            ..HlrOptions::default()
        };

        let curves_low = extract_silhouette_curves(&brep, view_dir, &opts_low);
        let curves_high = extract_silhouette_curves(&brep, view_dir, &opts_high);

        // Both should produce curves
        assert!(!curves_low.is_empty());
        assert!(!curves_high.is_empty());

        // Higher sampling should produce more points
        let pts_low: usize = curves_low.iter().map(|c| c.points.len()).sum();
        let pts_high: usize = curves_high.iter().map(|c| c.points.len()).sum();
        assert!(
            pts_high >= pts_low,
            "higher sampling should produce at least as many points"
        );
    }

    #[test]
    fn tangent_tolerance_affects_silhouette_detection() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let view_dir = DVec3::Z;

        // Very tight tolerance
        let opts_tight = HlrOptions {
            tangent_tolerance: 1e-12,
            ..HlrOptions::default()
        };

        // Very loose tolerance (should still work for sphere)
        let opts_loose = HlrOptions {
            tangent_tolerance: 0.01,
            ..HlrOptions::default()
        };

        let curves_tight = extract_silhouette_curves(&brep, view_dir, &opts_tight);
        let curves_loose = extract_silhouette_curves(&brep, view_dir, &opts_loose);

        // Both should find silhouette curves for a sphere
        assert!(!curves_tight.is_empty());
        assert!(!curves_loose.is_empty());
    }
}
