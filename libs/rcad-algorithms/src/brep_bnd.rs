//! BRepBndLib-style bounding box utilities for BRep topology.
//!
//! This module provides utilities analogous to OCCT's `BRepBndLib` class:
//!
//! - **AddShape**: Add topology shapes to bounding boxes
//! - **BoundingBox**: Axis-aligned bounding box type with geometric operations
//! - **ClosestPoint**: Compute closest point in bbox and distances
//! - **SurfaceBounds**: Surface and curve-specific bounding box computation
//!
//! # Example
//!
//! ```
//! # use rcad_algorithms::tolerance::*;
//! use rcad_algorithms::brep_bnd::*;
//! use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
//! use rcad_kernel::BRep;
//!
//! let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
//!     width: 1.0, height: 2.0, depth: 3.0
//! });
//!
//! // Compute bounding box for the entire BRep
//! let mut bbox = BoundingBox::new();
//! add_brep_to_bbox(&rcad_kernel::BRep, &mut bbox);
//!
//! assert!(bbox.is_valid());
//! assert!((bbox.size().x - 1.0).abs() < TOLERANCE_MESH_LEGACY);
//! assert!((bbox.size().y - 2.0).abs() < TOLERANCE_MESH_LEGACY);
//! assert!((bbox.size().z - 3.0).abs() < TOLERANCE_MESH_LEGACY);
//! ```

use glam::DVec3;
use rcad_kernel::{topods, Curve3, Surface3};
use rcad_kernel::geom::{CurveEval, SurfaceEval};

// =============================================================================
// Bounding Box Type
// =============================================================================

/// Axis-aligned bounding box for 3D geometry.
///
/// Provides efficient bounds computation, intersection testing,
/// and closest point queries.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    /// Minimum corner of the bounding box.
    pub min: DVec3,
    /// Maximum corner of the bounding box.
    pub max: DVec3,
}

impl Default for BoundingBox {
    fn default() -> Self {
        Self::new()
    }
}

impl BoundingBox {
    /// Create a new empty bounding box.
    ///
    /// An empty bounding box has min > max and contains no points.
    /// Use `add_point` to expand it.
    pub fn new() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
        }
    }

    /// Create a bounding box from a single point.
    pub fn from_point(p: DVec3) -> Self {
        Self { min: p, max: p }
    }

    /// Create a bounding box from two corner points.
    ///
    /// The points do not need to be ordered; the constructor will
    /// compute the correct min and max.
    pub fn from_corners(p1: DVec3, p2: DVec3) -> Self {
        Self {
            min: p1.min(p2),
            max: p1.max(p2),
        }
    }

    /// Create a bounding box from min and max corners.
    ///
    /// # Safety
    ///
    /// The caller must ensure that min <= max component-wise.
    pub fn from_min_max(min: DVec3, max: DVec3) -> Self {
        Self { min, max }
    }

    /// Check if the bounding box is valid (contains at least one point).
    pub fn is_valid(&self) -> bool {
        self.min.x <= self.max.x
            && self.min.y <= self.max.y
            && self.min.z <= self.max.z
            && self.min.x.is_finite()
            && self.max.x.is_finite()
    }

    /// Check if the bounding box is empty (contains no points).
    pub fn is_empty(&self) -> bool {
        !self.is_valid()
    }

    /// Add a point to the bounding box, expanding it if necessary.
    pub fn add_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    /// Add multiple points to the bounding box.
    pub fn add_points(&mut self, points: &[DVec3]) {
        for &p in points {
            self.add_point(p);
        }
    }

    /// Expand to include another bounding box.
    pub fn add_bbox(&mut self, other: &BoundingBox) {
        if other.is_valid() {
            self.min = self.min.min(other.min);
            self.max = self.max.max(other.max);
        }
    }

    /// Get the corner points as `[min, max]`.
    pub fn corners(&self) -> [DVec3; 2] {
        [self.min, self.max]
    }

    /// Get the size (dimensions) of the bounding box.
    ///
    /// Returns a vector where x, y, z are the extents along each axis.
    pub fn size(&self) -> DVec3 {
        self.max - self.min
    }

    /// Get the center point of the bounding box.
    pub fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    /// Get the diagonal length of the bounding box.
    ///
    /// This is the distance from min to max corner.
    pub fn diagonal(&self) -> f64 {
        (self.max - self.min).length()
    }

    /// Get the squared diagonal length.
    pub fn diagonal_squared(&self) -> f64 {
        (self.max - self.min).length_squared()
    }

    /// Get the volume of the bounding box.
    ///
    /// Returns 0.0 for empty or degenerate boxes.
    pub fn volume(&self) -> f64 {
        let size = self.size();
        if size.x > 0.0 && size.y > 0.0 && size.z > 0.0 {
            size.x * size.y * size.z
        } else {
            0.0
        }
    }

    /// Get the surface area of the bounding box.
    ///
    /// Useful for SAH (Surface Area Heuristic) in BVH construction.
    pub fn surface_area(&self) -> f64 {
        let d = self.size();
        if d.x > 0.0 && d.y > 0.0 && d.z > 0.0 {
            2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
        } else {
            0.0
        }
    }

    /// Check if a point is inside the bounding box (inclusive).
    ///
    /// Points exactly on the boundary are considered inside.
    pub fn contains(&self, p: DVec3, tol: f64) -> bool {
        p.x >= self.min.x - tol
            && p.x <= self.max.x + tol
            && p.y >= self.min.y - tol
            && p.y <= self.max.y + tol
            && p.z >= self.min.z - tol
            && p.z <= self.max.z + tol
    }

    /// Check if a point is strictly inside the bounding box.
    ///
    /// Points on the boundary are not considered inside.
    pub fn contains_strict(&self, p: DVec3) -> bool {
        p.x > self.min.x
            && p.x < self.max.x
            && p.y > self.min.y
            && p.y < self.max.y
            && p.z > self.min.z
            && p.z < self.max.z
    }

    /// Check if this bounding box intersects another.
    pub fn intersects(&self, other: &BoundingBox) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// Compute the intersection of two bounding boxes.
    ///
    /// Returns an empty box if there is no intersection.
    pub fn intersection(&self, other: &BoundingBox) -> BoundingBox {
        BoundingBox::from_min_max(
            self.min.max(other.min),
            self.max.min(other.max),
        )
    }

    /// Compute the union of two bounding boxes.
    pub fn union(&self, other: &BoundingBox) -> BoundingBox {
        let mut result = *self;
        result.add_bbox(other);
        result
    }

    /// Enlarge the bounding box by a tolerance value on all sides.
    pub fn enlarge(&mut self, tol: f64) {
        if self.is_valid() {
            self.min -= DVec3::splat(tol);
            self.max += DVec3::splat(tol);
        }
    }

    /// Get an enlarged copy of the bounding box.
    pub fn enlarged(&self, tol: f64) -> BoundingBox {
        let mut result = *self;
        result.enlarge(tol);
        result
    }

    /// Clamp a point to be inside the bounding box.
    pub fn clamp_point(&self, p: DVec3) -> DVec3 {
        p.clamp(self.min, self.max)
    }

    /// Find the closest point in the bounding box to a given point.
    ///
    /// If the point is inside the box, returns the point itself.
    /// Otherwise, returns the closest point on the box surface.
    pub fn closest_point(&self, p: DVec3) -> DVec3 {
        self.clamp_point(p)
    }

    /// Compute the distance from a point to the bounding box.
    ///
    /// Returns 0.0 if the point is inside the box.
    pub fn distance_to(&self, p: DVec3) -> f64 {
        let closest = self.closest_point(p);
        (p - closest).length()
    }

    /// Compute the squared distance from a point to the bounding box.
    ///
    /// More efficient than `distance_to` when only comparing distances.
    pub fn distance_squared_to(&self, p: DVec3) -> f64 {
        let closest = self.closest_point(p);
        (p - closest).length_squared()
    }

    /// Check if a ray intersects the bounding box.
    ///
    /// Returns the parameter t of the first intersection point along the ray,
    /// or None if there is no intersection.
    ///
    /// # Arguments
    ///
    /// * `origin` - Ray origin
    /// * `direction` - Ray direction (does not need to be normalized)
    pub fn ray_intersect(&self, origin: DVec3, direction: DVec3) -> Option<f64> {
        if direction.x == 0.0 && direction.y == 0.0 && direction.z == 0.0 {
            return None;
        }

        let inv_dir = DVec3::new(
            1.0 / direction.x,
            1.0 / direction.y,
            1.0 / direction.z,
        );

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

    /// Get the 8 corner vertices of the bounding box.
    ///
    /// Corners are returned in the order:
    /// [min, (max.x, min.y, min.z), (max.x, max.y, min.z), (min.x, max.y, min.z),
    ///  (min.x, min.y, max.z), (max.x, min.y, max.z), max, (min.x, max.y, max.z)]
    pub fn all_corners(&self) -> [DVec3; 8] {
        [
            self.min,
            DVec3::new(self.max.x, self.min.y, self.min.z),
            DVec3::new(self.max.x, self.max.y, self.min.z),
            DVec3::new(self.min.x, self.max.y, self.min.z),
            DVec3::new(self.min.x, self.min.y, self.max.z),
            DVec3::new(self.max.x, self.min.y, self.max.z),
            self.max,
            DVec3::new(self.min.x, self.max.y, self.max.z),
        ]
    }
}

// =============================================================================
// AddShape - Add BRep topology to bounding boxes
// =============================================================================

/// Add all geometry from a BRep to a bounding box.
///
/// This includes all vertices, and for faces/edges with geometry,
/// samples points on the surfaces/curves to get accurate bounds.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_bnd::*;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Sphere { radius: 1.0 });
/// let mut bbox = BoundingBox::new();
/// add_brep_to_bbox(&rcad_kernel::BRep, &mut bbox);
///
/// // Sphere should be bounded by [-1, 1] in y (poles)
/// assert!(bbox.is_valid());
/// ```
pub fn add_brep_to_bbox(brep: &rcad_kernel::BRep, bbox: &mut BoundingBox) {
    // Add all vertices
    for vertex in &brep.vertices {
        bbox.add_point(vertex.point);
    }

    // Add all faces (sample surfaces for accurate bounds)
    let n_faces = count_brep_faces(brep);
    for face_idx in 0..n_faces {
        add_face_to_bbox(brep, face_idx, bbox);
    }

    // Add edges that might not be part of faces
    for edge_idx in 0..brep.edges.len() {
        add_edge_to_bbox(brep, edge_idx, bbox);
    }
}

/// Add a face to a bounding box.
///
/// Samples the face's surface at grid points to get accurate bounds
/// for curved surfaces.
pub fn add_face_to_bbox(brep: &rcad_kernel::BRep, face_idx: usize, bbox: &mut BoundingBox) {
    // Get the face and its surface
    let (face, surface_idx) = match get_face_and_surface_idx(brep, face_idx) {
        Some(result) => result,
        None => return,
    };

    // Add points from the surface
    if let Some(surface) = brep.geom.surfaces.get(surface_idx) {
        add_surface_to_bbox(surface, brep, face_idx, bbox);
    }

    // Also add vertices from the face boundary
    for wire_edge in &face.outer_wire.edges {
        if wire_edge.idx < brep.edges.len() {
            let edge = &brep.edges[wire_edge.idx];
            if edge.start < brep.vertices.len() {
                bbox.add_point(brep.vertices[edge.start].point);
            }
            if edge.end < brep.vertices.len() {
                bbox.add_point(brep.vertices[edge.end].point);
            }
        }
    }

    // Add vertices from inner wires
    for inner_wire in &face.inner_wires {
        for wire_edge in &inner_wire.edges {
            if wire_edge.idx < brep.edges.len() {
                let edge = &brep.edges[wire_edge.idx];
                if edge.start < brep.vertices.len() {
                    bbox.add_point(brep.vertices[edge.start].point);
                }
                if edge.end < brep.vertices.len() {
                    bbox.add_point(brep.vertices[edge.end].point);
                }
            }
        }
    }
}

/// Add an edge to a bounding box.
///
/// Samples the edge's 3D curve to get accurate bounds for curved edges.
pub fn add_edge_to_bbox(brep: &rcad_kernel::BRep, edge_idx: usize, bbox: &mut BoundingBox) {
    if edge_idx >= brep.edges.len() {
        return;
    }

    let edge = &brep.edges[edge_idx];

    // Add edge vertices
    if edge.start < brep.vertices.len() {
        bbox.add_point(brep.vertices[edge.start].point);
    }
    if edge.end < brep.vertices.len() {
        bbox.add_point(brep.vertices[edge.end].point);
    }

    // Sample points along the curve if available
    if let Some(curve_idx) = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c)
        && let Some(curve) = brep.geom.curves.get(curve_idx) {
            let range = brep.geom.edge_curve_range.get(edge_idx)
                .copied()
                .flatten()
                .unwrap_or_else(|| curve.default_domain());

            // Sample along the curve
            let samples = 10;
            for i in 0..=samples {
                let t = range[0] + (range[1] - range[0]) * (i as f64) / (samples as f64);
                let p = curve.point_at(t);
                if p.is_finite() {
                    bbox.add_point(p);
                }
            }
        }
}

/// Add a vertex to a bounding box.
pub fn add_vertex_to_bbox(brep: &rcad_kernel::BRep, vertex_idx: usize, bbox: &mut BoundingBox) {
    if vertex_idx < brep.vertices.len() {
        bbox.add_point(brep.vertices[vertex_idx].point);
    }
}

/// Add surface geometry to bounding box with face-specific parameter range.
fn add_surface_to_bbox(surface: &Surface3, brep: &rcad_kernel::BRep, face_idx: usize, bbox: &mut BoundingBox) {
    // Get the parameter range for this face
    let domain = brep.geom.face_surface_range.get(face_idx)
        .copied()
        .flatten()
        .unwrap_or_else(|| surface.default_domain());

    let [u_min, u_max, v_min, v_max] = domain;

    // Sample a grid on the surface
    let n_u = 5;
    let n_v = 5;

    for i in 0..=n_u {
        for j in 0..=n_v {
            let u = u_min + (u_max - u_min) * (i as f64) / (n_u as f64);
            let v = v_min + (v_max - v_min) * (j as f64) / (n_v as f64);
            let p = surface.point_at(u, v);
            if p.is_finite() {
                bbox.add_point(p);
            }
        }
    }
}

// =============================================================================
// SurfaceBounds - Surface and curve-specific bounding boxes
// =============================================================================

/// Compute a bounding box for a surface.
///
/// Samples the surface over its default parameter domain.
pub fn surface_bounds(surface: &Surface3) -> BoundingBox {
    surface_bounds_with_domain(surface, surface.default_domain())
}

/// Compute a bounding box for a surface over a specific parameter domain.
///
/// # Arguments
///
/// * `surface` - The surface to bound
/// * `domain` - Parameter domain as `[u_min, u_max, v_min, v_max]`
pub fn surface_bounds_with_domain(surface: &Surface3, domain: [f64; 4]) -> BoundingBox {
    let [u_min, u_max, v_min, v_max] = domain;
    let mut bbox = BoundingBox::new();

    let n_u = 10;
    let n_v = 10;

    for i in 0..=n_u {
        for j in 0..=n_v {
            let u = u_min + (u_max - u_min) * (i as f64) / (n_u as f64);
            let v = v_min + (v_max - v_min) * (j as f64) / (n_v as f64);
            let p = surface.point_at(u, v);
            if p.is_finite() {
                bbox.add_point(p);
            }
        }
    }

    bbox
}

/// Compute a bounding box for a curve over a parameter range.
///
/// # Arguments
///
/// * `curve` - The curve to bound
/// * `t_min` - Start parameter
/// * `t_max` - End parameter
pub fn curve_bounds(curve: &Curve3, t_min: f64, t_max: f64) -> BoundingBox {
    curve_bounds_with_range(curve, [t_min, t_max])
}

/// Compute a bounding box for a curve over a specific parameter range.
///
/// # Arguments
///
/// * `curve` - The curve to bound
/// * `range` - Parameter range as `[t_min, t_max]`
pub fn curve_bounds_with_range(curve: &Curve3, range: [f64; 2]) -> BoundingBox {
    let [t_min, t_max] = range;
    let mut bbox = BoundingBox::new();

    let samples = 20;
    for i in 0..=samples {
        let t = t_min + (t_max - t_min) * (i as f64) / (samples as f64);
        let p = curve.point_at(t);
        if p.is_finite() {
            bbox.add_point(p);
        }
    }

    bbox
}

/// Compute a bounding box for a curve using its default parameter domain.
pub fn curve_bounds_default(curve: &Curve3) -> BoundingBox {
    let range = curve.default_domain();
    curve_bounds_with_range(curve, range)
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Count the total number of faces in a BRep.
fn count_brep_faces(brep: &rcad_kernel::BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Get a face and its surface index by flat index.
fn get_face_and_surface_idx(brep: &rcad_kernel::BRep, face_idx: usize) -> Option<(&rcad_kernel::topology::Face, usize)> {
    let mut current_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            if face_idx < current_idx + shell.faces.len() {
                let local_idx = face_idx - current_idx;
                let face = &shell.faces[local_idx];

                // Get the surface index
                let surface_idx = brep.geom.face_surface.get(face_idx)
                    .and_then(|s| *s)?;

                return Some((face, surface_idx));
            }
            current_idx += shell.faces.len();
        }
    }

    None
}

// =============================================================================
// Tests
// =============================================================================



