//! BRepExtrema-style distance/extrema calculations.
//!
//!
//! - `DistShapeShape` (class): distance between two shapes, with
//!   `perform()` / `is_done()` / `distance()` / `support_on_shape1()` / `support_on_shape2()`
//! - Utility free functions: distance_point_curve, closest_point_on_curve, etc.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, CurveEval, Line3, Surface3, SurfaceEval};
use rcad_kernel::topods;

use crate::boptools::bvh::{Aabb, Bvh};
use crate::tolerance::*;

/// Finite difference step size for derivative computation.
const H: f64 = TOLERANCE_MESH_LEGACY;

// =============================================================================
// BRepExtrema_DistShapeShape 閳?distance between two shapes (class)
// =============================================================================

/// BRepExtrema_DistShapeShape 閳?compute minimum distance
/// between two shapes using BVH-accelerated surface sampling.
///
/// Usage:
/// ```rust,ignore
/// let mut dss = DistShapeShape::new(s1, s2);
/// dss.perform();
/// if dss.is_done() {
///     let dist = dss.distance();
/// }
/// ```
pub struct DistShapeShape {
    shape1: Option<topods::BRep>,
    shape2: Option<topods::BRep>,
    distance: f64,
    support1: Vec<topods::Shape>,
    support2: Vec<topods::Shape>,
    pt1: DVec3,
    pt2: DVec3,
    performed: bool,
}

impl DistShapeShape {
    /// default constructor.
    pub fn new() -> Self {
        Self {
            shape1: None,
            shape2: None,
            distance: f64::INFINITY,
            support1: Vec::new(),
            support2: Vec::new(),
            pt1: DVec3::ZERO,
            pt2: DVec3::ZERO,
            performed: false,
        }
    }

    /// constructor with both shapes.
    pub fn with_shapes(s1: &topods::BRep, s2: &topods::BRep) -> Self {
        let mut dss = Self::new();
        dss.load_s1(s1);
        dss.load_s2(s2);
        dss
    }

    /// LoadS1.
    pub fn load_s1(&mut self, s1: &topods::BRep) {
        self.shape1 = Some(s1.clone());
    }

    /// LoadS2.
    pub fn load_s2(&mut self, s2: &topods::BRep) {
        self.shape2 = Some(s2.clone());
    }

    /// Perform 閳?compute the minimum distance.
    pub fn perform(&mut self) {
        let Some(ref s1) = self.shape1 else { return };
        let Some(ref s2) = self.shape2 else { return };
        let old1 = s1;
        let old2 = s2;
        let (dist, p1, p2) = distance_brep_brep(&old1, &old2);
        self.distance = dist;
        self.pt1 = p1;
        self.pt2 = p2;
        self.performed = true;
    }

    /// IsDone.
    pub fn is_done(&self) -> bool {
        self.performed
    }

    /// Distance.
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// SupportOnShape1 閳?the support elements on shape 1.
    pub fn support_on_shape1(&self) -> &[topods::Shape] {
        &self.support1
    }

    /// SupportOnShape2 閳?the support elements on shape 2.
    pub fn support_on_shape2(&self) -> &[topods::Shape] {
        &self.support2
    }

    /// Convenience: closest point on shape 1.
    pub fn point_on_shape1(&self) -> DVec3 {
        self.pt1
    }

    /// Convenience: closest point on shape 2.
    pub fn point_on_shape2(&self) -> DVec3 {
        self.pt2
    }
}

// =============================================================================
// DistShapeShape - Distance between shapes
// =============================================================================

/// Compute the Euclidean distance between two points.
pub fn distance_point_point(p1: DVec3, p2: DVec3) -> f64 {
    (p2 - p1).length()
}

/// Compute the distance from a point to a curve.
///
/// Returns the distance and the parameter value on the curve where the closest point lies.
/// Uses Newton iteration for accurate parameter refinement.
pub fn distance_point_curve(point: DVec3, curve: &Curve3) -> (f64, f64) {
    let (param, closest_pt) = closest_point_on_curve(curve, point);
    let distance = distance_point_point(point, closest_pt);
    (distance, param)
}

/// Compute the distance from a point to a surface.
///
/// Returns the distance and the UV parameters where the closest point lies.
/// Uses Newton iteration for accurate parameter refinement.
pub fn distance_point_surface(point: DVec3, surface: &Surface3) -> (f64, f64, f64) {
    let (uv, closest_pt) = closest_point_on_surface(surface, point);
    let distance = distance_point_point(point, closest_pt);
    (distance, uv.x, uv.y)
}

/// Compute the minimum distance between two curves.
///
/// Returns the distance and the two parameter values on each curve.
/// Uses sampling to find initial candidates, then Newton refinement.
pub fn distance_curve_curve(curve1: &Curve3, curve2: &Curve3) -> (f64, f64, f64) {
    let domain1 = curve_domain(curve1);
    let domain2 = curve_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 32;
    let mut best_dist = f64::INFINITY;
    let mut best_t1 = domain1[0];
    let mut best_t2 = domain2[0];

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < best_dist {
                best_dist = dist;
                best_t1 = t1;
                best_t2 = t2;
            }
        }
    }

    // Newton refinement
    let (refined_t1, refined_t2) =
        refine_curve_curve_distance(curve1, curve2, domain1, domain2, best_t1, best_t2);
    let p1 = curve1.point_at(refined_t1);
    let p2 = curve2.point_at(refined_t2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_t1, refined_t2)
}

/// Compute the minimum distance between a curve and a surface.
///
/// Returns the distance and the parameter on the curve plus UV on the surface.
/// Uses sampling to find initial candidates, then Newton refinement.
pub fn distance_curve_surface(curve: &Curve3, surface: &Surface3) -> (f64, f64, f64, f64) {
    let curve_domain = curve_domain(curve);
    let surf_domain = surface_domain(surface);

    // Sample curve and surface to find initial candidates
    let n_curve = 24;
    let n_surf = 12;
    let mut best_dist = f64::INFINITY;
    let mut best_t = curve_domain[0];
    let mut best_u = surf_domain[0];
    let mut best_v = surf_domain[2];

    for i in 0..=n_curve {
        let t = curve_domain[0] + (curve_domain[1] - curve_domain[0]) * i as f64 / n_curve as f64;
        let p_curve = curve.point_at(t);

        for j in 0..=n_surf {
            let u = surf_domain[0] + (surf_domain[1] - surf_domain[0]) * j as f64 / n_surf as f64;
            for k in 0..=n_surf {
                let v =
                    surf_domain[2] + (surf_domain[3] - surf_domain[2]) * k as f64 / n_surf as f64;
                let p_surf = surface.point_at(u, v);
                let dist = (p_surf - p_curve).length();

                if dist < best_dist {
                    best_dist = dist;
                    best_t = t;
                    best_u = u;
                    best_v = v;
                }
            }
        }
    }

    // Newton refinement
    let (refined_t, refined_u, refined_v) = refine_curve_surface_distance(
        curve,
        surface,
        curve_domain,
        surf_domain,
        best_t,
        best_u,
        best_v,
    );
    let p_curve = curve.point_at(refined_t);
    let p_surf = surface.point_at(refined_u, refined_v);
    let final_dist = (p_surf - p_curve).length();

    (final_dist, refined_t, refined_u, refined_v)
}

/// Compute the minimum distance between two surfaces.
///
/// Returns the distance and UV parameters on both surfaces.
/// Uses sampling to find initial candidates, then Newton refinement.
pub fn distance_surface_surface(surf1: &Surface3, surf2: &Surface3) -> (f64, f64, f64, f64, f64) {
    let domain1 = surface_domain(surf1);
    let domain2 = surface_domain(surf2);

    // Sample both surfaces to find initial candidates
    let n_samples = 10;
    let mut best_dist = f64::INFINITY;
    let mut best_u1 = domain1[0];
    let mut best_v1 = domain1[2];
    let mut best_u2 = domain2[0];
    let mut best_v2 = domain2[2];

    for i1 in 0..=n_samples {
        let u1 = domain1[0] + (domain1[1] - domain1[0]) * i1 as f64 / n_samples as f64;
        for j1 in 0..=n_samples {
            let v1 = domain1[2] + (domain1[3] - domain1[2]) * j1 as f64 / n_samples as f64;
            let p1 = surf1.point_at(u1, v1);

            for i2 in 0..=n_samples {
                let u2 = domain2[0] + (domain2[1] - domain2[0]) * i2 as f64 / n_samples as f64;
                for j2 in 0..=n_samples {
                    let v2 = domain2[2] + (domain2[3] - domain2[2]) * j2 as f64 / n_samples as f64;
                    let p2 = surf2.point_at(u2, v2);
                    let dist = (p2 - p1).length();

                    if dist < best_dist {
                        best_dist = dist;
                        best_u1 = u1;
                        best_v1 = v1;
                        best_u2 = u2;
                        best_v2 = v2;
                    }
                }
            }
        }
    }

    // Newton refinement
    let (refined_u1, refined_v1, refined_u2, refined_v2) = refine_surface_surface_distance(
        surf1, surf2, domain1, domain2, best_u1, best_v1, best_u2, best_v2,
    );

    let p1 = surf1.point_at(refined_u1, refined_v1);
    let p2 = surf2.point_at(refined_u2, refined_v2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_u1, refined_v1, refined_u2, refined_v2)
}

/// Compute the minimum distance between two BRep shapes.
///
/// Returns the distance and the two closest points.
/// Uses BVH acceleration for efficiency.
pub fn distance_brep_brep(
    brep1: &rcad_kernel::BRep,
    brep2: &rcad_kernel::BRep,
) -> (f64, DVec3, DVec3) {
    let bvh1 = Bvh::build(brep1);
    let bvh2 = Bvh::build(brep2);

    // Get candidate face pairs
    let candidate_pairs = Bvh::candidate_pairs(&bvh1, &bvh2);

    if candidate_pairs.is_empty() {
        // Fallback: compute bounding box centers distance
        let bb1 = compute_brep_aabb(brep1);
        let bb2 = compute_brep_aabb(brep2);
        let center1 = bb1.center();
        let center2 = bb2.center();
        return ((center2 - center1).length(), center1, center2);
    }

    let mut best_dist = f64::INFINITY;
    let mut best_pt1 = DVec3::ZERO;
    let mut best_pt2 = DVec3::ZERO;

    // Check each candidate pair
    for (fi1, fi2) in candidate_pairs {
        let surf1 = get_brep_surface(brep1, fi1);
        let surf2 = get_brep_surface(brep2, fi2);

        if let (Some(s1), Some(s2)) = (surf1, surf2) {
            let (dist, u1, v1, u2, v2) = distance_surface_surface(&s1, &s2);
            if dist < best_dist {
                best_dist = dist;
                best_pt1 = s1.point_at(u1, v1);
                best_pt2 = s2.point_at(u2, v2);
            }
        }
    }

    // Also check vertex-to-vertex distances
    for ts1 in &brep1.tshapes {
        if let topods::TShape::Vertex(v1) = ts1.as_ref() {
            for ts2 in &brep2.tshapes {
                if let topods::TShape::Vertex(v2) = ts2.as_ref() {
                    let dist = (v2.point - v1.point).length();
                    if dist < best_dist {
                        best_dist = dist;
                        best_pt1 = v1.point;
                        best_pt2 = v2.point;
                    }
                }
            }
        }
    }

    (best_dist, best_pt1, best_pt2)
}

// =============================================================================
// Extrema - Find extremum points
// =============================================================================

/// Find the n closest points on a curve to a given point.
///
/// Returns a vector of (parameter, distance) pairs sorted by distance.
pub fn find_closest_points(curve: &Curve3, point: DVec3, n_points: usize) -> Vec<(f64, f64)> {
    let domain = curve_domain(curve);

    // Sample the curve to find local minima
    let n_samples = 100;
    let mut candidates: Vec<(f64, f64)> = Vec::new();

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        candidates.push((t, dist));
    }

    // Find local minima
    let mut local_minima: Vec<(f64, f64)> = Vec::new();
    for i in 1..candidates.len() - 1 {
        if candidates[i].1 < candidates[i - 1].1 && candidates[i].1 < candidates[i + 1].1 {
            // Refine using Newton
            let refined_t = refine_point_curve_distance(curve, domain, point, candidates[i].0);
            let refined_dist = (curve.point_at(refined_t) - point).length();
            local_minima.push((refined_t, refined_dist));
        }
    }

    // Also include endpoints
    let (t0, d0) = candidates[0];
    let (tn, dn) = candidates[candidates.len() - 1];
    local_minima.push((t0, d0));
    local_minima.push((tn, dn));

    // Sort by distance and take top n
    local_minima.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    local_minima.truncate(n_points);

    local_minima
}

/// Find the furthest points on a BRep in a given direction.
///
/// Returns the two points that are furthest apart when projected onto the direction.
pub fn find_furthest_points(brep: &rcad_kernel::BRep, direction: DVec3) -> (DVec3, DVec3) {
    let dir = direction.normalize();
    let mut min_proj = f64::INFINITY;
    let mut max_proj = f64::NEG_INFINITY;
    let mut min_point = DVec3::ZERO;
    let mut max_point = DVec3::ZERO;

    // Check all vertices
    for ts in &brep.tshapes {
        if let topods::TShape::Vertex(v) = ts.as_ref() {
            let proj = v.point.dot(dir);
            if proj < min_proj {
                min_proj = proj;
                min_point = v.point;
            }
            if proj > max_proj {
                max_proj = proj;
                max_point = v.point;
            }
        }
    }

    // Also sample face interiors for more accuracy
    let face_indices = get_all_face_indices(brep);
    for face_idx in face_indices {
        if let Some(surf) = get_brep_surface(brep, face_idx) {
            let domain = surface_domain(&surf);
            let n_samples = 5;
            for i in 0..=n_samples {
                for j in 0..=n_samples {
                    let u = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
                    let v = domain[2] + (domain[3] - domain[2]) * j as f64 / n_samples as f64;
                    let p = surf.point_at(u, v);
                    let proj = p.dot(dir);
                    if proj < min_proj {
                        min_proj = proj;
                        min_point = p;
                    }
                    if proj > max_proj {
                        max_proj = proj;
                        max_point = p;
                    }
                }
            }
        }
    }

    (min_point, max_point)
}

// =============================================================================
// ClosestPoint - Closest point queries
// =============================================================================

/// Find the closest point on a curve to a given point.
///
/// Returns the parameter value and the closest point on the curve.
/// Uses Newton iteration for accuracy.
pub fn closest_point_on_curve(curve: &Curve3, point: DVec3) -> (f64, DVec3) {
    let domain = curve_domain(curve);

    // Initial guess by sampling
    let n_samples = 50;
    let mut best_t = domain[0];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        if dist < best_dist {
            best_dist = dist;
            best_t = t;
        }
    }

    // Newton refinement
    let refined_t = refine_point_curve_distance(curve, domain, point, best_t);
    let closest = curve.point_at(refined_t);

    (refined_t, closest)
}

/// Find the closest point on a surface to a given point.
///
/// Returns the UV parameters and the closest point on the surface.
/// Uses Newton iteration for accuracy.
pub fn closest_point_on_surface(surface: &Surface3, point: DVec3) -> (DVec2, DVec3) {
    let domain = surface_domain(surface);

    // Initial guess by sampling
    let n_samples = 20;
    let mut best_u = domain[0];
    let mut best_v = domain[2];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let u = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        for j in 0..=n_samples {
            let v = domain[2] + (domain[3] - domain[2]) * j as f64 / n_samples as f64;
            let p = surface.point_at(u, v);
            let dist = (p - point).length();
            if dist < best_dist {
                best_dist = dist;
                best_u = u;
                best_v = v;
            }
        }
    }

    // Newton refinement
    let (refined_u, refined_v) =
        refine_point_surface_distance(surface, domain, point, best_u, best_v);
    let closest = surface.point_at(refined_u, refined_v);

    (DVec2::new(refined_u, refined_v), closest)
}

// =============================================================================
// SupportShapes - Find supporting geometry
// =============================================================================

/// Find the face that supports a point (closest face).
///
/// Returns the face index if found.
pub fn find_supporting_face(brep: &rcad_kernel::BRep, point: DVec3) -> Option<usize> {
    let mut best_face = None;
    let mut best_dist = f64::INFINITY;
    let tolerance = TOLERANCE_ABS * 100.0;

    let face_indices = get_all_face_indices(brep);
    for face_idx in face_indices {
        if let Some(surf) = get_brep_surface(brep, face_idx) {
            let (uv, closest) = closest_point_on_surface(&surf, point);
            let dist = (closest - point).length();

            // Check if point is within the face boundary (UV domain)
            let domain = surface_domain(&surf);
            if uv.x >= domain[0] - tolerance
                && uv.x <= domain[1] + tolerance
                && uv.y >= domain[2] - tolerance
                && uv.y <= domain[3] + tolerance
                && dist < best_dist
            {
                best_dist = dist;
                best_face = Some(face_idx);
            }
        }
    }

    best_face
}

/// Find the edge that supports a point (closest edge).
///
/// Returns the edge index if found.
pub fn find_supporting_edge(brep: &rcad_kernel::BRep, point: DVec3) -> Option<usize> {
    let mut best_edge = None;
    let mut best_dist = f64::INFINITY;
    let tolerance = TOLERANCE_ABS * 100.0;

    for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Edge(ed) = ts.as_ref() {
            // Try to get the curve from geometry data, or create a line from vertices
            let curve = if let Some(c) = ed.curve.clone() {
                c
            } else {
                // Create an implicit line from edge vertex positions
                let start_pt = brep.vertex_point(ed.first.index)?;
                let end_pt = brep.vertex_point(ed.last.index)?;
                let dir = (end_pt - start_pt).normalize();
                Curve3::Line(Line3 {
                    origin: start_pt,
                    direction: dir,
                })
            };

            let (t, closest) = closest_point_on_curve(&curve, point);
            let dist = (closest - point).length();

            // Check if parameter is within edge range
            let edge_range = ed.range;

            if t >= edge_range[0] - tolerance && t <= edge_range[1] + tolerance && dist < best_dist
            {
                best_dist = dist;
                best_edge = Some(edge_idx);
            }
        }
    }

    best_edge
}

// =============================================================================
// Internal helper functions
// =============================================================================

/// Get the domain for a curve, handling infinite lines specially.
fn curve_domain(curve: &Curve3) -> [f64; 2] {
    match curve {
        Curve3::Line(_) => [-1e6, 1e6], // Clamp infinite lines to large range
        other => other.default_domain(),
    }
}

/// Get the domain for a surface, handling infinite domains (planes) specially.
fn surface_domain(surface: &Surface3) -> [f64; 4] {
    let domain = surface.default_domain();
    let u0 = if domain[0].is_infinite() {
        -10.0
    } else {
        domain[0]
    };
    let u1 = if domain[1].is_infinite() {
        10.0
    } else {
        domain[1]
    };
    let v0 = if domain[2].is_infinite() {
        -10.0
    } else {
        domain[2]
    };
    let v1 = if domain[3].is_infinite() {
        10.0
    } else {
        domain[3]
    };
    [u0, u1, v0, v1]
}

/// Compute the AABB of a BRep.
fn compute_brep_aabb(brep: &rcad_kernel::BRep) -> Aabb {
    let mut aabb = Aabb::empty();
    for ts in &brep.tshapes {
        if let topods::TShape::Vertex(v) = ts.as_ref() {
            aabb.expand_point(v.point);
        }
    }
    aabb
}

/// Get a surface from a BRep by face tshape index.
fn get_brep_surface(brep: &rcad_kernel::BRep, face_idx: usize) -> Option<Surface3> {
    brep.tshapes.get(face_idx).and_then(|ts| {
        if let topods::TShape::Face(fd) = ts.as_ref() {
            fd.surface.clone()
        } else {
            None
        }
    })
}

/// Get all face tshape indices from a BRep.
fn get_all_face_indices(brep: &rcad_kernel::BRep) -> Vec<usize> {
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(_, ts)| matches!(ts.as_ref(), topods::TShape::Face(_)))
        .map(|(i, _)| i)
        .collect()
}

/// Get a curve from a BRep by edge tshape index.
fn get_brep_curve(brep: &rcad_kernel::BRep, edge_idx: usize) -> Option<Curve3> {
    brep.tshapes.get(edge_idx).and_then(|ts| {
        if let topods::TShape::Edge(ed) = ts.as_ref() {
            ed.curve.clone()
        } else {
            None
        }
    })
}

/// Compute curve derivative via finite differences.
fn curve_derivative(curve: &Curve3, t: f64) -> DVec3 {
    (curve.point_at(t + H) - curve.point_at(t - H)) / (2.0 * H)
}

/// Compute curve second derivative via finite differences.
fn curve_second_derivative(curve: &Curve3, t: f64) -> DVec3 {
    let d_plus = curve_derivative(curve, t + H);
    let d_minus = curve_derivative(curve, t - H);
    (d_plus - d_minus) / (2.0 * H)
}

/// Compute surface partial derivatives via finite differences.
fn surface_derivatives(surface: &Surface3, u: f64, v: f64) -> (DVec3, DVec3) {
    let du = (surface.point_at(u + H, v) - surface.point_at(u - H, v)) / (2.0 * H);
    let dv = (surface.point_at(u, v + H) - surface.point_at(u, v - H)) / (2.0 * H);
    (du, dv)
}

/// Newton refinement for point-to-curve distance.
fn refine_point_curve_distance(
    curve: &Curve3,
    domain: [f64; 2],
    point: DVec3,
    initial_t: f64,
) -> f64 {
    let mut t = initial_t;

    const MAX_ITER: usize = 20;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p = curve.point_at(t);
        let d = curve_derivative(curve, t);

        let diff = p - point;
        let f = diff.dot(d);

        let d2 = curve_second_derivative(curve, t);
        let df = d.dot(d) + diff.dot(d2);

        if df.abs() < TOL {
            break;
        }

        let delta = -f / df;
        t += delta;

        // Clamp to domain
        t = t.clamp(domain[0], domain[1]);

        if delta.abs() < TOL {
            break;
        }
    }

    t
}

/// Newton refinement for point-to-surface distance.
fn refine_point_surface_distance(
    surface: &Surface3,
    domain: [f64; 4],
    point: DVec3,
    initial_u: f64,
    initial_v: f64,
) -> (f64, f64) {
    let mut u = initial_u;
    let mut v = initial_v;

    const MAX_ITER: usize = 20;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p = surface.point_at(u, v);
        let (du, dv) = surface_derivatives(surface, u, v);

        let diff = p - point;

        // Gradient of distance squared
        let fu = diff.dot(du);
        let fv = diff.dot(dv);

        // Hessian approximation using finite differences for second derivatives
        let (du_du, du_dv) = surface_derivatives(surface, u + H, v);
        let (_dv_du, dv_dv) = surface_derivatives(surface, u, v + H);

        let d2uu = (du_du - du) / H;
        let d2vv = (dv_dv - dv) / H;
        let d2uv = (du_dv - du) / H;

        let fuu = du.dot(du) + diff.dot(d2uu);
        let fvv = dv.dot(dv) + diff.dot(d2vv);
        let fuv = du.dot(dv) + diff.dot(d2uv);

        // Solve 2x2 system
        let det = fuu * fvv - fuv * fuv;
        if det.abs() < TOL {
            break;
        }

        let du_param = (-fu * fvv + fv * fuv) / det;
        let dv_param = (-fv * fuu + fu * fuv) / det;

        u += du_param;
        v += dv_param;

        // Clamp to domain
        u = u.clamp(domain[0], domain[1]);
        v = v.clamp(domain[2], domain[3]);

        if du_param.abs() < TOL && dv_param.abs() < TOL {
            break;
        }
    }

    (u, v)
}

/// Newton refinement for curve-to-curve distance.
fn refine_curve_curve_distance(
    curve1: &Curve3,
    curve2: &Curve3,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    let mut t1 = t1;
    let mut t2 = t2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);

        let d1 = curve_derivative(curve1, t1);
        let d2 = curve_derivative(curve2, t2);

        let diff = p1 - p2;

        // Gradient
        let f1 = diff.dot(d1);
        let f2 = -diff.dot(d2);

        // Hessian
        let d1_2 = curve_second_derivative(curve1, t1);
        let d2_2 = curve_second_derivative(curve2, t2);

        let h11 = d1.dot(d1) + diff.dot(d1_2);
        let h22 = d2.dot(d2) - diff.dot(d2_2);
        let h12 = -d1.dot(d2);

        let det = h11 * h22 - h12 * h12;
        if det.abs() < TOL {
            break;
        }

        let dt1 = (-f1 * h22 + f2 * h12) / det;
        let dt2 = (-f2 * h11 + f1 * h12) / det;

        t1 += dt1;
        t2 += dt2;

        t1 = t1.clamp(domain1[0], domain1[1]);
        t2 = t2.clamp(domain2[0], domain2[1]);

        if dt1.abs() < TOL && dt2.abs() < TOL {
            break;
        }
    }

    (t1, t2)
}

/// Newton refinement for curve-to-surface distance.
fn refine_curve_surface_distance(
    curve: &Curve3,
    surface: &Surface3,
    curve_domain: [f64; 2],
    surf_domain: [f64; 4],
    t: f64,
    u: f64,
    v: f64,
) -> (f64, f64, f64) {
    let mut t = t;
    let mut u = u;
    let mut v = v;

    const MAX_ITER: usize = 30;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let pc = curve.point_at(t);
        let ps = surface.point_at(u, v);

        let dc = curve_derivative(curve, t);
        let (ds_u, ds_v) = surface_derivatives(surface, u, v);

        let diff = pc - ps;

        // Gradient
        let ft = diff.dot(dc);
        let fu = -diff.dot(ds_u);
        let fv = -diff.dot(ds_v);

        // Simple gradient descent step (more robust than full Newton for 3D problems)
        let step = 0.5;
        let htt = dc.dot(dc).max(TOL);
        let huu = ds_u.dot(ds_u).max(TOL);
        let hvv = ds_v.dot(ds_v).max(TOL);

        t -= step * ft / htt;
        u -= step * fu / huu;
        v -= step * fv / hvv;

        t = t.clamp(curve_domain[0], curve_domain[1]);
        u = u.clamp(surf_domain[0], surf_domain[1]);
        v = v.clamp(surf_domain[2], surf_domain[3]);

        if ft.abs() < TOL && fu.abs() < TOL && fv.abs() < TOL {
            break;
        }
    }

    (t, u, v)
}

/// Newton refinement for surface-to-surface distance.
fn refine_surface_surface_distance(
    surf1: &Surface3,
    surf2: &Surface3,
    domain1: [f64; 4],
    domain2: [f64; 4],
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
) -> (f64, f64, f64, f64) {
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;
    const STEP: f64 = 0.3;

    for _ in 0..MAX_ITER {
        let p1 = surf1.point_at(u1, v1);
        let p2 = surf2.point_at(u2, v2);

        let (du1, dv1) = surface_derivatives(surf1, u1, v1);
        let (du2, dv2) = surface_derivatives(surf2, u2, v2);

        let diff = p1 - p2;

        // Gradient
        let fu1 = diff.dot(du1);
        let fv1 = diff.dot(dv1);
        let fu2 = -diff.dot(du2);
        let fv2 = -diff.dot(dv2);

        // Simple gradient descent
        u1 -= STEP * fu1 / (du1.dot(du1) + TOL);
        v1 -= STEP * fv1 / (dv1.dot(dv1) + TOL);
        u2 -= STEP * fu2 / (du2.dot(du2) + TOL);
        v2 -= STEP * fv2 / (dv2.dot(dv2) + TOL);

        u1 = u1.clamp(domain1[0], domain1[1]);
        v1 = v1.clamp(domain1[2], domain1[3]);
        u2 = u2.clamp(domain2[0], domain2[1]);
        v2 = v2.clamp(domain2[2], domain2[3]);

        if fu1.abs() < TOL && fv1.abs() < TOL && fu2.abs() < TOL && fv2.abs() < TOL {
            break;
        }
    }

    (u1, v1, u2, v2)
}

// =============================================================================
// Tests
// =============================================================================
