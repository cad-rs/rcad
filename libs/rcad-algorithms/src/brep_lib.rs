//! BRepLib-style utilities for low-level BRep operations.
//!
//! This module provides utilities analogous to OCCT's `BRepLib` class:
//!
//! - **FindSurface**: Find a surface through a set of edges/points
//! - **SortFaces**: Sort faces by area, bounding box, etc.
//! - **CheckSameDomain**: Check if two faces share the same surface
//! - **Add**: Add geometry to BRep
//! - **Make**: Low-level shape construction
//! - **Bounds**: Compute parameter bounds
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_lib::*;
//! use rcad_kernel::BRep;
//!
//! let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
//! width: 1.0, height: 1.0, depth: 1.0
//! });
//!
//! // Sort faces by area
//! let sorted = sort_faces_by_area(&brep);
//! assert!(!sorted.is_empty());
//! ```

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::{Curve3, Surface3};
use rcad_kernel::topods::{Orientation, TShape};
use rcad_kernel::geom::{any_perpendicular, CurveEval, SurfaceEval};
use rcad_kernel::topology::{Wire, WireEdge};

// =============================================================================
// Error Types
// =============================================================================

/// Errors that can occur during BRepLib operations.
#[derive(Debug, Clone)]
pub enum BRepLibError {
 /// Invalid index provided.
 InvalidIndex {
 kind: &'static str,
 index: usize,
 max: usize,
 },
 /// Missing geometry for the specified element.
 MissingGeometry {
 kind: &'static str,
 index: usize,
 },
 /// Failed to fit a surface to the given data.
 SurfaceFitFailed(String),
 /// Empty input data.
 EmptyInput,
 /// Invalid wire topology.
 InvalidWire(String),
 /// Numerical failure.
 NumericalFailure(String),
}

impl std::fmt::Display for BRepLibError {
 fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
 match self {
 BRepLibError::InvalidIndex { kind, index, max } => {
 write!(f, "Invalid {} index {} (max {})", kind, index, max)
 }
 BRepLibError::MissingGeometry { kind, index } => {
 write!(f, "Missing {} geometry at index {}", kind, index)
 }
 BRepLibError::SurfaceFitFailed(msg) => {
 write!(f, "Surface fit failed: {}", msg)
 }
 BRepLibError::EmptyInput => write!(f, "Empty input"),
 BRepLibError::InvalidWire(msg) => write!(f, "Invalid wire: {}", msg),
 BRepLibError::NumericalFailure(msg) => write!(f, "Numerical failure: {}", msg),
 }
 }
}

impl std::error::Error for BRepLibError {}

// =============================================================================
// FindSurface - Find a surface through edges or points
// =============================================================================

/// Result of finding a surface through edges or points.
#[derive(Debug, Clone)]
pub struct FoundSurface {
 /// The fitted surface.
 pub surface: Surface3,
 /// RMS error of the fit.
 pub rms_error: f64,
 /// Type of surface that was fitted.
 pub surface_type: FittedSurfaceType,
}

/// Classification of the fitted surface type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FittedSurfaceType {
 /// Plane surface.
 Plane,
 /// Cylindrical surface.
 Cylinder,
 /// Conical surface.
 Cone,
 /// Spherical surface.
 Sphere,
 /// BSpline surface (general freeform).
 BSpline,
}

/// Find a surface through a set of edges by sampling points from the edge curves.
///
/// Samples points from the 3D curves of the specified edges and fits a surface.
/// Tries plane, cylinder, cone, and sphere fits before falling back to BSpline.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edges
/// * `edge_indices` - Indices of edges to sample points from
///
/// # Returns
///
/// A `FoundSurface` containing the best-fit surface and fit quality metrics.
pub fn find_surface_through_edges(
 brep: &rcad_kernel::BRep,
 edge_indices: &[usize],
) -> Result<FoundSurface, BRepLibError> {
 if edge_indices.is_empty() {
 return Err(BRepLibError::EmptyInput);
 }

 // Sample points from edges
 let mut points = Vec::new();
 let samples_per_edge = 10;

 for &edge_idx in edge_indices {
 if edge_idx >= brep.edge_count() {
 return Err(BRepLibError::InvalidIndex {
 kind: "edge",
 index: edge_idx,
 max: brep.edge_count(),
 });
 }

 // Get edge vertices
 let ed = match &*brep.tshapes[edge_idx] {
 TShape::Edge(ed) => ed,
 _ => unreachable!(),
 };
 let start_pt = brep.vertex_point(ed.first.index).unwrap();
 let end_pt = brep.vertex_point(ed.last.index).unwrap();

 // Try to sample from the curve if available
 if let Some(curve) = &ed.curve {
 let range = ed.range;

 // Sample along the curve
 for i in 0..samples_per_edge {
 let t = range[0] + (range[1] - range[0]) * (i as f64) / (samples_per_edge - 1) as f64;
 points.push(curve.point_at(t));
 }
 continue;
 }

 // Fallback: interpolate between vertices
 for i in 0..samples_per_edge {
 let t = (i as f64) / (samples_per_edge - 1) as f64;
 points.push(start_pt.lerp(end_pt, t));
 }
 }

 find_surface_through_points(&points, TOLERANCE_MESH_LEGACY)
}

/// Find a surface through a set of points.
///
/// Tries plane, cylinder, cone, and sphere fits in order of increasing complexity.
/// Returns the best fit based on RMS error.
///
/// # Arguments
///
/// * `points` - The points to fit a surface through
/// * `tol` - Tolerance for determining fit quality (relative to point cloud scale)
///
/// # Returns
///
/// A `FoundSurface` containing the best-fit surface and fit quality metrics.
pub fn find_surface_through_points(
 points: &[DVec3],
 tol: f64,
) -> Result<FoundSurface, BRepLibError> {
 if points.len() < 3 {
 return Err(BRepLibError::EmptyInput);
 }

 // Compute bounding box for scale
 let (bb_min, bb_max) = compute_bounding_box(points);
 let bb_size = (bb_max - bb_min).length();
 let scale_tol = tol * bb_size;

 // Try fitting in order of complexity
 let mut best_fit: Option<FoundSurface> = None;
 let mut best_error = f64::INFINITY;

 // Try plane fit first
 if let Some(plane) = fit_plane_to_points(points)
 && plane.rms_error < best_error && plane.rms_error < scale_tol {
 best_error = plane.rms_error;
 best_fit = Some(FoundSurface {
 surface: Surface3::Plane(rcad_kernel::geom::Plane {
 origin: plane.point,
 normal: plane.normal,
 }),
 rms_error: plane.rms_error,
 surface_type: FittedSurfaceType::Plane,
 });
 }

 // Try sphere fit
 if points.len() >= 4
 && let Some(sphere) = fit_sphere_to_points(points)
 && sphere.rms_error < best_error && sphere.rms_error < scale_tol {
 best_error = sphere.rms_error;
 best_fit = Some(FoundSurface {
 surface: Surface3::Sphere(rcad_kernel::geom::SphericalSurface::new(
 sphere.center, DVec3::Z, sphere.radius,
 )),
 rms_error: sphere.rms_error,
 surface_type: FittedSurfaceType::Sphere,
 });
 }

 // Try cylinder fit
 if points.len() >= 5
 && let Some(cylinder) = fit_cylinder_to_points(points)
 && cylinder.rms_error < best_error && cylinder.rms_error < scale_tol {
 best_error = cylinder.rms_error;
 best_fit = Some(FoundSurface {
 surface: Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
 origin: cylinder.axis_point,
 axis: cylinder.axis_direction,
 radius: cylinder.radius,
 ref_dir: any_perpendicular(cylinder.axis_direction),
 }),
 rms_error: cylinder.rms_error,
 surface_type: FittedSurfaceType::Cylinder,
 });
 }

 // Try cone fit
 if points.len() >= 5
 && let Some(cone) = fit_cone_to_points(points)
 && cone.rms_error < best_error && cone.rms_error < scale_tol {
 best_error = cone.rms_error;
 best_fit = Some(FoundSurface {
 surface: Surface3::Cone(rcad_kernel::geom::ConicalSurface {
 apex: cone.apex,
 axis: cone.axis,
 radius: 0.0, // Reference radius at apex
 half_angle_rad: cone.semi_angle,
 }),
 rms_error: cone.rms_error,
 surface_type: FittedSurfaceType::Cone,
 });
 }

 // If no analytic surface fit well, create a BSpline surface
 if best_fit.is_none() {
 let bspline = fit_bspline_surface_to_points(points)?;
 best_fit = Some(bspline);
 }

 best_fit.ok_or_else(|| BRepLibError::SurfaceFitFailed("Could not fit any surface type".into()))
}

// =============================================================================
// Surface Fitting Helpers
// =============================================================================

/// Fitted plane result.
struct FittedPlane {
 point: DVec3,
 normal: DVec3,
 rms_error: f64,
}

/// Fitted sphere result.
struct FittedSphere {
 center: DVec3,
 radius: f64,
 rms_error: f64,
}

/// Fitted cylinder result.
struct FittedCylinder {
 axis_point: DVec3,
 axis_direction: DVec3,
 radius: f64,
 rms_error: f64,
}

/// Fitted cone result.
struct FittedCone {
 apex: DVec3,
 axis: DVec3,
 semi_angle: f64,
 rms_error: f64,
}

fn compute_bounding_box(points: &[DVec3]) -> (DVec3, DVec3) {
 let mut min = DVec3::splat(f64::INFINITY);
 let mut max = DVec3::splat(f64::NEG_INFINITY);
 for &p in points {
 min = min.min(p);
 max = max.max(p);
 }
 (min, max)
}

fn compute_centroid(points: &[DVec3]) -> DVec3 {
 points.iter().sum::<DVec3>() / points.len() as f64
}

fn fit_plane_to_points(points: &[DVec3]) -> Option<FittedPlane> {
 if points.len() < 3 {
 return None;
 }

 let centroid = compute_centroid(points);

 // Compute covariance matrix
 let mut cov = [[0.0; 3]; 3];
 for &p in points {
 let d = p - centroid;
 cov[0][0] += d.x * d.x;
 cov[0][1] += d.x * d.y;
 cov[0][2] += d.x * d.z;
 cov[1][1] += d.y * d.y;
 cov[1][2] += d.y * d.z;
 cov[2][2] += d.z * d.z;
 }
 cov[1][0] = cov[0][1];
 cov[2][0] = cov[0][2];
 cov[2][1] = cov[1][2];

 // Find eigenvector for smallest eigenvalue using power iteration on inverse
 // For plane fitting, the normal is the direction of minimum variance
 let (eigenvalues, eigenvectors) = compute_eigendecomposition(&cov);

 // Find index of smallest eigenvalue
 let min_idx = eigenvalues.iter()
 .enumerate()
 .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
 .map(|(i, _)| i)?;

 let normal = eigenvectors[min_idx].normalize_or(DVec3::Z);

 // Compute RMS error
 let mut sum_sq = 0.0;
 for &p in points {
 let d = (p - centroid).dot(normal);
 sum_sq += d * d;
 }
 let rms_error = (sum_sq / points.len() as f64).sqrt();

 Some(FittedPlane {
 point: centroid,
 normal,
 rms_error,
 })
}

fn fit_sphere_to_points(points: &[DVec3]) -> Option<FittedSphere> {
 if points.len() < 4 {
 return None;
 }

 let n = points.len() as f64;

 // Algebraic fit using linear least squares
 // (x - cx)^2 + (y - cy)^2 + (z - cz)^2 = r^2
 // x^2 + y^2 + z^2 - 2*cx*x - 2*cy*y - 2*cz*z + (cx^2 + cy^2 + cz^2 - r^2) = 0
 let mut sums = [0.0; 10]; // x, y, z, x2, y2, z2, xy, xz, yz, x2y2z2
 for &p in points {
 let x = p.x;
 let y = p.y;
 let z = p.z;
 sums[0] += x;
 sums[1] += y;
 sums[2] += z;
 sums[3] += x * x;
 sums[4] += y * y;
 sums[5] += z * z;
 sums[6] += x * y;
 sums[7] += x * z;
 sums[8] += y * z;
 sums[9] += x * x + y * y + z * z;
 }

 // Solve 4x4 system for sphere parameters
 let a = [
 [2.0 * sums[0], 2.0 * sums[1], 2.0 * sums[2], n],
 [2.0 * sums[3], 2.0 * sums[6], 2.0 * sums[7], sums[0]],
 [2.0 * sums[6], 2.0 * sums[4], 2.0 * sums[8], sums[1]],
 [2.0 * sums[7], 2.0 * sums[8], 2.0 * sums[5], sums[2]],
 ];
 let b = [sums[9], sums[3] + sums[9], sums[4] + sums[9], sums[5] + sums[9]];

 let coeffs = solve_linear_4x4(&a, &b)?;

 let center = DVec3::new(coeffs[0], coeffs[1], coeffs[2]);
 let radius = (coeffs[0] * coeffs[0] + coeffs[1] * coeffs[1] + coeffs[2] * coeffs[2] + coeffs[3]).sqrt().max(TOLERANCE_LINEAR_ULTRA_STRICT);

 // Compute RMS error
 let mut sum_sq = 0.0;
 for &p in points {
 let d = (p - center).length() - radius;
 sum_sq += d * d;
 }
 let rms_error = (sum_sq / n).sqrt();

 Some(FittedSphere {
 center,
 radius,
 rms_error,
 })
}

fn fit_cylinder_to_points(points: &[DVec3]) -> Option<FittedCylinder> {
 if points.len() < 5 {
 return None;
 }

 let centroid = compute_centroid(points);

 // Estimate axis direction using PCA
 let mut cov = [[0.0; 3]; 3];
 for &p in points {
 let d = p - centroid;
 cov[0][0] += d.x * d.x;
 cov[0][1] += d.x * d.y;
 cov[0][2] += d.x * d.z;
 cov[1][1] += d.y * d.y;
 cov[1][2] += d.y * d.z;
 cov[2][2] += d.z * d.z;
 }
 cov[1][0] = cov[0][1];
 cov[2][0] = cov[0][2];
 cov[2][1] = cov[1][2];

 let (eigenvalues, eigenvectors) = compute_eigendecomposition(&cov);

 // For a cylinder, the axis is the direction of maximum variance
 let max_idx = eigenvalues.iter()
 .enumerate()
 .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
 .map(|(i, _)| i)?;

 let axis = eigenvectors[max_idx].normalize_or(DVec3::Z);

 // Project points onto plane perpendicular to axis
 let u = if axis.x.abs() < 0.9 {
 axis.cross(DVec3::X).normalize_or(DVec3::Y)
 } else {
 axis.cross(DVec3::Y).normalize_or(DVec3::X)
 };
 let v = axis.cross(u);

 // Project points and fit circle
 let projected: Vec<(f64, f64)> = points.iter()
 .map(|&p| {
 let d = p - centroid;
 (d.dot(u), d.dot(v))
 })
 .collect();

 let circle = fit_circle_2d(&projected)?;

 let axis_point = centroid + circle.0 * u + circle.1 * v;
 let radius = circle.2;

 // Compute RMS error
 let mut sum_sq = 0.0;
 for &p in points {
 let to_point = p - axis_point;
 let axial = to_point.dot(axis);
 let radial = to_point - axial * axis;
 let d = radial.length() - radius;
 sum_sq += d * d;
 }
 let rms_error = (sum_sq / points.len() as f64).sqrt();

 Some(FittedCylinder {
 axis_point,
 axis_direction: axis,
 radius,
 rms_error,
 })
}

fn fit_cone_to_points(points: &[DVec3]) -> Option<FittedCone> {
 if points.len() < 5 {
 return None;
 }

 let centroid = compute_centroid(points);

 // Estimate axis direction using PCA
 let mut cov = [[0.0; 3]; 3];
 for &p in points {
 let d = p - centroid;
 cov[0][0] += d.x * d.x;
 cov[0][1] += d.x * d.y;
 cov[0][2] += d.x * d.z;
 cov[1][1] += d.y * d.y;
 cov[1][2] += d.y * d.z;
 cov[2][2] += d.z * d.z;
 }
 cov[1][0] = cov[0][1];
 cov[2][0] = cov[0][2];
 cov[2][1] = cov[1][2];

 let (_eigenvalues, eigenvectors) = compute_eigendecomposition(&cov);

 // For a cone, the axis might be any of the principal directions
 // Try all three and pick the best fit
 let mut best_fit: Option<FittedCone> = None;
 let mut best_error = f64::INFINITY;

 for i in 0..3 {
 let axis = eigenvectors[i].normalize_or(DVec3::Z);
 if let Some(cone) = fit_cone_with_axis(points, centroid, axis)
 && cone.rms_error < best_error {
 best_error = cone.rms_error;
 best_fit = Some(cone);
 }
 }

 best_fit
}

fn fit_cone_with_axis(points: &[DVec3], centroid: DVec3, axis: DVec3) -> Option<FittedCone> {
 // For each point, compute the distance from the axis and the axial position
 let mut min_axial = f64::INFINITY;
 let mut max_axial = f64::NEG_INFINITY;

 for &p in points {
 let d = p - centroid;
 let axial = d.dot(axis);
 min_axial = min_axial.min(axial);
 max_axial = max_axial.max(axial);
 }

 // Estimate apex position and semi-angle using least squares
 // For a cone: distance_from_axis = tan(semi_angle) * distance_from_apex_along_axis
 // Let apex be at centroid - t * axis for some t
 // Then: radial_dist = tan(semi_angle) * (axial_dist + t)

 // Use linear regression to find tan(semi_angle) and t
 let mut sum_r = 0.0;
 let mut sum_a = 0.0;
 let mut sum_ra = 0.0;
 let mut sum_a2 = 0.0;
 let n = points.len() as f64;

 for &p in points {
 let d = p - centroid;
 let axial = d.dot(axis);
 let radial = (d - axial * axis).length();
 sum_r += radial;
 sum_a += axial;
 sum_ra += radial * axial;
 sum_a2 += axial * axial;
 }

 // Solve for tan(angle) and offset
 let denom = n * sum_a2 - sum_a * sum_a;
 if denom.abs() < TOLERANCE_FLOAT_LOOSE {
 return None;
 }

 let tan_angle = (n * sum_ra - sum_r * sum_a) / denom;
 let t = (sum_r * sum_a2 - sum_a * sum_ra) / denom;

 if tan_angle.abs() < TOLERANCE_MESH_LEGACY || tan_angle.abs() > 100.0 {
 return None; // Not a reasonable cone
 }

 let semi_angle = tan_angle.atan().abs();
 let apex = centroid - t * axis;

 // Compute RMS error
 let mut sum_sq = 0.0;
 for &p in points {
 let d = p - apex;
 let axial = d.dot(axis);
 let radial = (d - axial * axis).length();
 let expected_radial = semi_angle.tan() * axial.abs();
 let error = if expected_radial > TOLERANCE_LINEAR_ULTRA_STRICT {
 (radial - expected_radial).abs()
 } else {
 radial
 };
 sum_sq += error * error;
 }
 let rms_error = (sum_sq / n).sqrt();

 Some(FittedCone {
 apex,
 axis,
 semi_angle,
 rms_error,
 })
}

fn fit_bspline_surface_to_points(points: &[DVec3]) -> Result<FoundSurface, BRepLibError> {
 // Create a simple bilinear BSpline surface through the points
 // For simplicity, we create a degree 1 (bilinear) surface
 let (_bb_min, _bb_max) = compute_bounding_box(points);
 let _centroid = compute_centroid(points);

 // Create a simple plane as fallback
 let plane = fit_plane_to_points(points).ok_or_else(|| {
 BRepLibError::SurfaceFitFailed("Cannot fit BSpline: insufficient points".into())
 })?;

 // For now, return a plane as the BSpline fallback
 // In a full implementation, this would create an actual BSpline surface
 Ok(FoundSurface {
 surface: Surface3::Plane(rcad_kernel::geom::Plane {
 origin: plane.point,
 normal: plane.normal,
 }),
 rms_error: plane.rms_error,
 surface_type: FittedSurfaceType::BSpline,
 })
}

fn compute_eigendecomposition(m: &[[f64; 3]; 3]) -> ([f64; 3], [DVec3; 3]) {
 // Jacobi eigenvalue algorithm for symmetric 3x3 matrices
 let mut a = *m;
 let mut v = [DVec3::X, DVec3::Y, DVec3::Z];

 const MAX_ITER: usize = 100;
 const TOL: f64 = TOLERANCE_LEN_MIN;

 for _ in 0..MAX_ITER {
 // Find largest off-diagonal element
 let mut max_val = 0.0;
 let (mut p, mut q) = (0, 1);

 for i in 0..3 {
 for j in (i + 1)..3 {
 if a[i][j].abs() > max_val {
 max_val = a[i][j].abs();
 p = i;
 q = j;
 }
 }
 }

 if max_val < TOL {
 break;
 }

 // Compute rotation angle
 let theta = if (a[p][p] - a[q][q]).abs() < TOL {
 std::f64::consts::FRAC_PI_4 * a[p][q].signum()
 } else {
 0.5 * (2.0 * a[p][q] / (a[p][p] - a[q][q])).atan()
 };

 let c = theta.cos();
 let s = theta.sin();

 // Update matrix
 let app = a[p][p];
 let aqq = a[q][q];
 let apq = a[p][q];

 a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
 a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
 a[p][q] = 0.0;
 a[q][p] = 0.0;

 for i in 0..3 {
 if i != p && i != q {
 let aip = a[i][p];
 let aiq = a[i][q];
 a[i][p] = c * aip - s * aiq;
 a[p][i] = a[i][p];
 a[i][q] = s * aip + c * aiq;
 a[q][i] = a[i][q];
 }
 }

 // Update eigenvectors
 for i in 0..3 {
 let vip = v[i][p];
 let viq = v[i][q];
 v[i][p] = c * vip - s * viq;
 v[i][q] = s * vip + c * viq;
 }
 }

 // Normalize eigenvectors
 for i in 0..3 {
 v[i] = v[i].normalize_or(match i {
 0 => DVec3::X,
 1 => DVec3::Y,
 _ => DVec3::Z,
 });
 }

 ([a[0][0], a[1][1], a[2][2]], v)
}

fn fit_circle_2d(points: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
 if points.len() < 3 {
 return None;
 }

 let n = points.len() as f64;

 let mut sums = [0.0; 7]; // x, y, x2, y2, xy, x2y2, x2y2 * x or y
 for &(x, y) in points {
 sums[0] += x;
 sums[1] += y;
 sums[2] += x * x;
 sums[3] += y * y;
 sums[4] += x * y;
 let r2 = x * x + y * y;
 sums[5] += r2;
 sums[6] += r2 * (x + y);
 }

 // Solve 3x3 system for circle center and radius
 let a = [
 [2.0 * sums[0], 2.0 * sums[1], n],
 [2.0 * sums[2], 2.0 * sums[4], sums[0]],
 [2.0 * sums[4], 2.0 * sums[3], sums[1]],
 ];
 let b = [sums[5], sums[2] + sums[5], sums[3] + sums[5]];

 let coeffs = solve_linear_3x3(&a, &b)?;

 let cx = coeffs[0];
 let cy = coeffs[1];
 let radius = (cx * cx + cy * cy + coeffs[2]).sqrt().max(TOLERANCE_LINEAR_ULTRA_STRICT);

 Some((cx, cy, radius))
}

fn solve_linear_3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
 let mut m = *a;
 let mut v = *b;

 // Forward elimination with partial pivoting
 for col in 0..3 {
 let mut max_row = col;
 let mut max_val = m[col][col].abs();
 for row in (col + 1)..3 {
 if m[row][col].abs() > max_val {
 max_val = m[row][col].abs();
 max_row = row;
 }
 }

 if max_val < TOLERANCE_FLOAT_LOOSE {
 return None;
 }

 m.swap(col, max_row);
 v.swap(col, max_row);

 for row in (col + 1)..3 {
 let factor = m[row][col] / m[col][col];
 for j in col..3 {
 m[row][j] -= factor * m[col][j];
 }
 v[row] -= factor * v[col];
 }
 }

 // Back substitution
 let mut x = [0.0; 3];
 for i in (0..3).rev() {
 let mut sum = v[i];
 for j in (i + 1)..3 {
 sum -= m[i][j] * x[j];
 }
 x[i] = sum / m[i][i];
 }

 Some(x)
}

fn solve_linear_4x4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> Option<[f64; 4]> {
 let mut m = *a;
 let mut v = *b;

 // Forward elimination with partial pivoting
 for col in 0..4 {
 let mut max_row = col;
 let mut max_val = m[col][col].abs();
 for row in (col + 1)..4 {
 if m[row][col].abs() > max_val {
 max_val = m[row][col].abs();
 max_row = row;
 }
 }

 if max_val < TOLERANCE_FLOAT_LOOSE {
 return None;
 }

 m.swap(col, max_row);
 v.swap(col, max_row);

 for row in (col + 1)..4 {
 let factor = m[row][col] / m[col][col];
 for j in col..4 {
 m[row][j] -= factor * m[col][j];
 }
 v[row] -= factor * v[col];
 }
 }

 // Back substitution
 let mut x = [0.0; 4];
 for i in (0..4).rev() {
 let mut sum = v[i];
 for j in (i + 1)..4 {
 sum -= m[i][j] * x[j];
 }
 x[i] = sum / m[i][i];
 }

 Some(x)
}

// =============================================================================
// SortFaces - Sort faces by area, bbox, etc.
// =============================================================================

/// Sort faces by area in descending order (largest first).
///
/// Returns indices of faces sorted by their surface area.
pub fn sort_faces_by_area(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let areas: Vec<(usize, f64)> = (0..count_faces(brep))
 .filter_map(|i| {
 let area = compute_face_area(brep, i);
 Some((i, area))
 })
 .collect();

 let mut sorted = areas;
 sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
 sorted.into_iter().map(|(i, _)| i).collect()
}

/// Sort faces by bounding box volume in descending order (largest first).
///
/// Returns indices of faces sorted by their bounding box volume.
pub fn sort_faces_by_bounding_box(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let volumes: Vec<(usize, f64)> = (0..count_faces(brep))
 .filter_map(|i| {
 let bb = compute_face_bounding_box(brep, i)?;
 let volume = (bb[1].x - bb[0].x) * (bb[1].y - bb[0].y) * (bb[1].z - bb[0].z);
 Some((i, volume))
 })
 .collect();

 let mut sorted = volumes;
 sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
 sorted.into_iter().map(|(i, _)| i).collect()
}

/// Sort faces by their distance from a reference point (nearest first).
///
/// Returns indices of faces sorted by distance from the reference point
/// to the face centroid.
pub fn sort_faces_by_distance(brep: &rcad_kernel::BRep, reference: DVec3) -> Vec<usize> {
 let distances: Vec<(usize, f64)> = (0..count_faces(brep))
 .filter_map(|i| {
 let centroid = compute_face_centroid(brep, i)?;
 let dist = (centroid - reference).length();
 Some((i, dist))
 })
 .collect();

 let mut sorted = distances;
 sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
 sorted.into_iter().map(|(i, _)| i).collect()
}

// =============================================================================
// CheckSameDomain - Check if two faces share the same surface
// =============================================================================

/// Check if two faces share the same underlying surface geometry.
///
/// Two faces share the same domain if:
/// 1. They reference the same surface, OR
/// 2. Their surfaces are geometrically equivalent within tolerance
pub fn faces_share_surface(brep: &rcad_kernel::BRep, face1_idx: usize, face2_idx: usize) -> Result<bool, BRepLibError> {
 let faces: Vec<&Surface3> = brep.tshapes.iter().filter_map(|ts| {
 if let TShape::Face(fd) = &**ts { fd.surface.as_ref() } else { None }
 }).collect();
 let n_faces = faces.len();
 if face1_idx >= n_faces {
 return Err(BRepLibError::InvalidIndex {
 kind: "face",
 index: face1_idx,
 max: n_faces,
 });
 }
 if face2_idx >= n_faces {
 return Err(BRepLibError::InvalidIndex {
 kind: "face",
 index: face2_idx,
 max: n_faces,
 });
 }

 let s1 = faces[face1_idx];
 let s2 = faces[face2_idx];

 // Check geometric equivalence
 surfaces_equivalent(s1, s2)
}

/// Check if two surfaces are geometrically equivalent within tolerance.
fn surfaces_equivalent(s1: &Surface3, s2: &Surface3) -> Result<bool, BRepLibError> {
 const TOL: f64 = TOLERANCE_MESH_LEGACY;

 match (s1, s2) {
 (Surface3::Plane(p1), Surface3::Plane(p2)) => {
 // Check if normals are parallel and planes pass through same point
 let normal_dot = p1.normal.dot(p2.normal).abs();
 if normal_dot < 1.0 - TOL {
 return Ok(false);
 }
 let dist = (p1.origin - p2.origin).dot(p1.normal).abs();
 Ok(dist < TOL)
 }
 (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
 Ok((s1.center - s2.center).length() < TOL && (s1.radius - s2.radius).abs() < TOL)
 }
 (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
 let axis_dot = c1.axis.normalize_or(DVec3::Z).dot(c2.axis.normalize_or(DVec3::Z)).abs();
 let same_origin = (c1.origin - c2.origin).length() < TOL;
 let same_radius = (c1.radius - c2.radius).abs() < TOL;
 Ok(axis_dot > 1.0 - TOL && same_origin && same_radius)
 }
 (Surface3::Cone(c1), Surface3::Cone(c2)) => {
 let axis_dot = c1.axis.normalize_or(DVec3::Z).dot(c2.axis.normalize_or(DVec3::Z)).abs();
 let same_apex = (c1.apex - c2.apex).length() < TOL;
 let same_angle = (c1.half_angle_rad - c2.half_angle_rad).abs() < TOL;
 Ok(axis_dot > 1.0 - TOL && same_apex && same_angle)
 }
 (Surface3::Torus(t1), Surface3::Torus(t2)) => {
 let same_center = (t1.center - t2.center).length() < TOL;
 let same_major = (t1.major_radius - t2.major_radius).abs() < TOL;
 let same_minor = (t1.minor_radius - t2.minor_radius).abs() < TOL;
 Ok(same_center && same_major && same_minor)
 }
 // For other surface types, just check if they're the same variant
 _ => Ok(false),
 }
}

// =============================================================================
// Add - Add geometry to BRep
// =============================================================================

/// Add an edge with a 3D curve to the BRep.
///
/// Creates a new edge connecting two existing vertices with the given curve.
/// Returns the index of the new edge.
pub fn add_edge_with_curve(
 brep: &mut rcad_kernel::BRep,
 curve: Curve3,
 start_vertex: usize,
 end_vertex: usize,
) -> Result<usize, BRepLibError> {
 if start_vertex >= brep.vertex_count() {
 return Err(BRepLibError::InvalidIndex {
 kind: "vertex",
 index: start_vertex,
 max: brep.vertex_count(),
 });
 }
 if end_vertex >= brep.vertex_count() {
 return Err(BRepLibError::InvalidIndex {
 kind: "vertex",
 index: end_vertex,
 max: brep.vertex_count(),
 });
 }

 let domain = curve.default_domain();
 let edge_idx = brep.add_edge_flat(start_vertex, end_vertex, Some(curve), domain);
 Ok(edge_idx)
}

/// Add a face with a surface to the BRep.
///
/// Creates a new face with the given surface and wires.
/// Returns the index of the new face (flat index across all solids/shells).
///
/// Note: This function requires the BRep to have at least one solid with one shell.
/// If no shell exists, a new solid/shell structure is created.
pub fn add_face_with_surface(
 brep: &mut rcad_kernel::BRep,
 surface: Surface3,
 wires: Vec<Wire>,
) -> Result<usize, BRepLibError> {
 use std::sync::Arc;
 use rcad_kernel::topods::ShapeRef;

 if wires.is_empty() {
 return Err(BRepLibError::InvalidWire("Face must have at least one wire".into()));
 }

 // Convert each old-style Wire to a tshape wire (ShapeRef chain of edges)
 let mut tshape_wires: Vec<ShapeRef> = Vec::with_capacity(wires.len());
 for w in wires {
 let edge_refs: Vec<ShapeRef> = w.edges.into_iter().map(|we| {
 let ts = &brep.tshapes[we.idx];
 ShapeRef {
 ptr_id: Arc::as_ptr(ts) as u64,
 index: we.idx,
 orientation: if we.forward { Orientation::Forward } else { Orientation::Reversed },
 location: 0,
 }
 }).collect();
 tshape_wires.push(brep.add_twire(edge_refs));
 }

 // First wire is outer, rest are inner
 let outer_wire = tshape_wires.swap_remove(0);
 let inner_wires = tshape_wires;

 let uv_domain = Some(surface.default_domain());
 let face_ref = brep.add_tface(Some(surface), outer_wire, inner_wires, None, uv_domain, Vec::new(), false);
 let face_idx = face_ref.index;

 // Ensure we have a solid/shell structure
 if brep.solid_count() == 0 {
 brep.add_tsolid(Vec::new());
 }

 // Find the first solid and add the face to its first shell
 let mut solid_idx = None;
 for (i, ts) in brep.tshapes.iter().enumerate() {
 if let TShape::Solid(_) = &**ts {
 solid_idx = Some(i);
 break;
 }
 }
 if let Some(si) = solid_idx {
 // Check if the solid already has shells
 let has_shells = {
 let ts = &brep.tshapes[si];
 match &**ts {
 TShape::Solid(sd) => !sd.shells.is_empty(),
 _ => false,
 }
 };
 if has_shells {
 // Get the first shell's tshape index
 let sh_idx = {
 let ts = &brep.tshapes[si];
 match &**ts {
 TShape::Solid(sd) => sd.shells[0].index,
 _ => unreachable!(),
 }
 };
 // Add face to existing shell
 if let TShape::Shell(shd) = &mut *std::sync::Arc::make_mut(&mut brep.tshapes[sh_idx]) {
 shd.faces.push(face_ref);
 }
 } else {
 // Create a shell with this face
 let shell_ref = brep.add_tshell(vec![face_ref]);
 if let TShape::Solid(sd) = &mut *std::sync::Arc::make_mut(&mut brep.tshapes[si]) {
 sd.shells.push(shell_ref);
 }
 }
 }

 Ok(face_idx)
}

// =============================================================================
// Make - Low-level shape construction
// =============================================================================

/// Make an edge from a curve.
///
/// Creates an edge data structure with the given curve.
/// The vertices must be added separately.
pub fn make_edge_from_curve(curve: Curve3) -> EdgeData {
 let domain = curve.default_domain();
 EdgeData {
 curve,
 parameter_range: Some(domain),
 }
}

/// Edge data structure for construction.
#[derive(Debug, Clone)]
pub struct EdgeData {
 pub curve: Curve3,
 pub parameter_range: Option<[f64; 2]>,
}

/// Make a face from a surface.
///
/// Creates a face data structure with the given surface.
/// The wires must be added separately.
pub fn make_face_from_surface(surface: Surface3) -> FaceData {
 let domain = surface.default_domain();
 FaceData {
 surface,
 parameter_range: Some(domain),
 }
}

/// Face data structure for construction.
#[derive(Debug, Clone)]
pub struct FaceData {
 pub surface: Surface3,
 pub parameter_range: Option<[f64; 4]>,
}

/// Make a wire from edges.
///
/// Creates a wire from a sequence of edge indices with directions.
pub fn make_wire_from_edges(edges: Vec<(usize, bool)>) -> Wire {
 Wire {
 edges: edges.into_iter()
 .map(|(idx, forward)| WireEdge::new(idx, forward))
 .collect(),
 }
}

// =============================================================================
// Bounds - Compute parameter bounds
// =============================================================================

/// Compute the parameter bounds of an edge's 3D curve.
///
/// Returns `[t_min, t_max]` for the edge's parameter range.
pub fn compute_edge_bounds(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<[f64; 2], BRepLibError> {
 if edge_idx >= brep.edge_count() {
 return Err(BRepLibError::InvalidIndex {
 kind: "edge",
 index: edge_idx,
 max: brep.edge_count(),
 });
 }

 let ed = match &*brep.tshapes[edge_idx] {
 TShape::Edge(ed) => ed,
 _ => unreachable!(),
 };

 // Use edge's stored range
 Ok(ed.range)
}

/// Compute the parameter bounds of a face's surface.
///
/// Returns `[u_min, u_max, v_min, v_max]` for the face's parameter range.
pub fn compute_face_bounds(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<[f64; 4], BRepLibError> {
 // Collect face surfaces in flat order (iterate tshapes for TShape::Face)
 let face_surfaces: Vec<Option<[f64; 4]>> = brep.tshapes.iter().filter_map(|ts| {
 if let TShape::Face(fd) = &**ts { Some(fd.uv_domain) } else { None }
 }).collect();
 let n_faces = face_surfaces.len();
 if face_idx >= n_faces {
 return Err(BRepLibError::InvalidIndex {
 kind: "face",
 index: face_idx,
 max: n_faces,
 });
 }

 // Use stored uv_domain, or compute from surface
 if let Some(domain) = face_surfaces[face_idx] {
 return Ok(domain);
 }

 // Get face and use surface's default domain
 let faces: Vec<Option<&Surface3>> = brep.tshapes.iter().filter_map(|ts| {
 if let TShape::Face(fd) = &**ts { Some(fd.surface.as_ref()) } else { None }
 }).collect();
 if let Some(Some(surface)) = faces.get(face_idx) {
 Ok(surface.default_domain())
 } else {
 Err(BRepLibError::MissingGeometry {
 kind: "surface",
 index: face_idx,
 })
 }
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Count the total number of faces in a BRep.
/// Count the total number of faces in a BRep.
fn count_faces(brep: &rcad_kernel::BRep) -> usize {
 brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Face(_))).count()
}

/// Collect vertex points from a face's outer wire via tshape structure.
fn collect_face_wire_vertices(brep: &rcad_kernel::BRep, face_ts_idx: usize) -> Vec<DVec3> {
 let mut points = Vec::new();
 let fd = match &*brep.tshapes[face_ts_idx] {
 TShape::Face(fd) => fd,
 _ => return points,
 };
 let outer_wire = &fd.outer_wire;
 if let TShape::Wire(twd) = &*brep.tshapes[outer_wire.index] {
 for &edge_ref in &twd.edges {
 if let TShape::Edge(ed) = &*brep.tshapes[edge_ref.index] {
 let (first_pt, last_pt) = (
 brep.vertex_point(ed.first.index).unwrap_or(DVec3::ZERO),
 brep.vertex_point(ed.last.index).unwrap_or(DVec3::ZERO),
 );
 match edge_ref.orientation {
 Orientation::Forward => points.push(first_pt),
 Orientation::Reversed => points.push(last_pt),
 _ => {
 points.push(first_pt);
 points.push(last_pt);
 }
 }
 }
 }
 }
 points
}

/// Compute the approximate surface area of a face.
fn compute_face_area(brep: &rcad_kernel::BRep, face_idx: usize) -> f64 {
 // Find the face by flat index
 let face_ts_idx = match brep.tshapes.iter().enumerate().filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_))).nth(face_idx) {
 Some((idx, _)) => idx,
 None => return 0.0,
 };

 let points = collect_face_wire_vertices(brep, face_ts_idx);
 if points.len() < 3 {
 return 0.0;
 }

 // Approximate area from wire polygon
 let centroid = points.iter().sum::<DVec3>() / points.len() as f64;
 let mut sum = 0.0;
 for i in 0..points.len() {
 let j = (i + 1) % points.len();
 sum += (points[i] - centroid).cross(points[j] - centroid).length();
 }
 sum * 0.5
}

/// Compute the bounding box of a face.
fn compute_face_bounding_box(brep: &rcad_kernel::BRep, face_idx: usize) -> Option<[DVec3; 2]> {
 let face_ts_idx = brep.tshapes.iter().enumerate().filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_))).nth(face_idx)?.0;
 let points = collect_face_wire_vertices(brep, face_ts_idx);

 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
 for &p in &points {
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }

 if min_pt.x.is_finite() {
 Some([min_pt, max_pt])
 } else {
 None
 }
}

/// Compute the centroid of a face.
fn compute_face_centroid(brep: &rcad_kernel::BRep, face_idx: usize) -> Option<DVec3> {
 let face_ts_idx = brep.tshapes.iter().enumerate().filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_))).nth(face_idx)?.0;
 let points = collect_face_wire_vertices(brep, face_ts_idx);

 if points.is_empty() {
 return None;
 }
 Some(points.iter().sum::<DVec3>() / points.len() as f64)
}

/// Compute a default normal for a surface.
fn compute_surface_normal(surface: &Surface3) -> DVec3 {
 match surface {
 Surface3::Plane(p) => p.normal,
 Surface3::Sphere(_s) => DVec3::Z,
 Surface3::Cylinder(c) => c.axis.normalize_or(DVec3::Z),
 Surface3::Cone(c) => c.axis.normalize_or(DVec3::Z),
 Surface3::Torus(t) => t.axis.normalize_or(DVec3::Z),
 _ => DVec3::Z,
 }
}

// =============================================================================
// Geometric properties (moved from brep_algo for OCCT-package alignment)
// =============================================================================

/// Total volume of an old BRep.
pub fn total_volume(brep: &rcad_kernel::BRep) -> f64 {
    rcad_kernel::volume(brep)
}

/// Total volume of a topods::BRep.
pub fn total_volume_topods(brep: &rcad_kernel::topods::BRep) -> f64 {
    rcad_kernel::volume(brep)
}

/// Total surface area of an old BRep.
pub fn total_surface_area(brep: &rcad_kernel::BRep) -> f64 {
    rcad_kernel::surface_area(brep)
}

/// Total surface area of a topods::BRep.
pub fn total_surface_area_topods(brep: &rcad_kernel::topods::BRep) -> f64 {
    rcad_kernel::surface_area(brep)
}

// =============================================================================
// Tests
// =============================================================================


