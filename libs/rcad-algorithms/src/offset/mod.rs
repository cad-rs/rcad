//! Shell and solid offset operations  ?analogous to OCCT `rcad_kernel::BRepOffsetAPI_MakeOffsetShape`.
//!
//! # Overview
//!
//! This module provides algorithms for offsetting shells and solids:
//!
//! - **`offset_shell`**: Offset all faces of a shell along their normals
//! - **`offset_solid`**: Create a new solid by offsetting (positive = outward, negative = inward)
//! - **`hollow_solid`**: Create a thin-walled solid by removing faces and offsetting remaining faces inward
//!
//! # Supported Surfaces
//!
//! Plane, Sphere, Cylinder, Cone, Torus  ?each has a known parallel-surface construction.
//! B-spline and Bezier surfaces use the `OffsetSurface` wrapper.
//!
//! # Join Types
//!
//! - **Intersection**: Sharp corners at edge intersections (default)
//! - **Arc**: Round corners using fillet arcs at edges
//! - **Tangent**: Smooth transitions between adjacent faces
//!
//! # Algorithm
//!
//! 1. Compute offset surfaces for each face
//! 2. Compute offset curves for each edge (intersection of adjacent offset surfaces)
//! 3. Compute offset vertices (intersection of three or more offset curves)
//! 4. Handle edge extension/intersection for gaps
//! 5. Build result shell from offset faces
//! 6. Apply join type handling (arc, intersection, tangent)
//! 7. Check for self-intersection and repair if enabled
//!
//! # Variable Thickness
//!
//! Per-face thickness can be specified for non-uniform offsets:
//! - Thickness interpolation across transition zones
//! - Smooth blending between different thickness regions
//!
//! # References
//!
//! - OCCT `rcad_kernel::BRepOffsetAPI_MakeOffsetShape`
//! - OCCT `rcad_kernel::BRepOffset_MakeOffset`
//! - OCCT `rcad_kernel::BRepOffset_Mode` (join types)

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::{
 SurfaceEval, CurveEval, any_perpendicular,
 geom::{Curve3, Surface3, Line3, Plane, Circle3, Ellipse3, Parabola3, Hyperbola3, CylindricalSurface, SphericalSurface, ConicalSurface, ToroidalSurface, OffsetSurface},
 topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
};
use crate::tolerance::*;
use crate::inttools::cone_cone::{ConeConeResult, intersect_cone_cone};
use crate::inttools::cylinder_cone::{CylinderConeResult, intersect_cylinder_cone};
use crate::inttools::plane_cone::{PlaneConicalResult, intersect_plane_cone};
use crate::inttools::sphere_cone::{SphereConeResult, intersect_sphere_cone_with_tolerance};
use crate::inttools::torus_torus::{TorusTorusResult, intersect_torus_torus};
use crate::inttools::cylinder_torus::{CylinderTorusResult, intersect_cylinder_torus};
use crate::inttools::torus_cone::{TorusConeResult, intersect_torus_cone};
use crate::inttools::plane_torus::{PlaneTorusResult, intersect_plane_torus};
use crate::inttools::sphere_torus::{SphereTorusResult, intersect_sphere_torus};

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Error Types
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Errors that can occur during offset operations.
#[derive(Debug, Clone)]
pub enum OffsetError {
 /// Offset distance is zero.
 ZeroDistance,
 /// Input shape is empty or invalid.
 InvalidInput(&'static str),
 /// Offset would create a degenerate surface (e.g., sphere radius goes negative).
 DegenerateSurface {
 face_index: usize,
 distance: f64,
 },
 /// Self-intersection detected during offset.
 SelfIntersection {
 description: String,
 },
 /// Failed to compute offset edge intersection.
 EdgeIntersectionFailed {
 edge_index: usize,
 },
 /// Failed to compute offset vertex.
 VertexComputationFailed {
 vertex_index: usize,
 },
 /// Geometry not supported for offset.
 UnsupportedGeometry {
 face_index: usize,
 geometry_type: String,
 },
 /// Numerical failure during computation.
 NumericalFailure(&'static str),
 /// Result has no valid faces.
 EmptyResult,
 /// Wall thickness violation.
 WallThicknessViolation {
 minimum: f64,
 actual: f64,
 location: String,
 },
 /// Failed to create join geometry.
 JoinCreationFailed {
 join_type: JoinType,
 edge_index: usize,
 reason: String,
 },
 /// Invalid variable thickness specification.
 InvalidVariableThickness {
 face_index: usize,
 thickness: f64,
 reason: String,
 },
 /// Offset quality check failed.
 QualityCheckFailed {
 metric: String,
 value: f64,
 threshold: f64,
 },
 /// Recovery failed after self-intersection.
 RecoveryFailed {
 attempts: usize,
 last_error: String,
 },
}

impl std::fmt::Display for OffsetError {
 fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
 match self {
 Self::ZeroDistance => write!(f, "offset distance is zero"),
 Self::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
 Self::DegenerateSurface { face_index, distance } => {
 write!(f, "offset distance {} would degenerate face {}", distance, face_index)
 }
 Self::SelfIntersection { description } => {
 write!(f, "self-intersection detected: {}", description)
 }
 Self::EdgeIntersectionFailed { edge_index } => {
 write!(f, "failed to compute offset edge intersection for edge {}", edge_index)
 }
 Self::VertexComputationFailed { vertex_index } => {
 write!(f, "failed to compute offset vertex {}", vertex_index)
 }
 Self::UnsupportedGeometry { face_index, geometry_type } => {
 write!(f, "unsupported geometry '{}' for face {}", geometry_type, face_index)
 }
 Self::NumericalFailure(msg) => write!(f, "numerical failure: {msg}"),
 Self::EmptyResult => write!(f, "offset produced no valid faces"),
 Self::WallThicknessViolation { minimum, actual, location } => {
 write!(f, "wall thickness {} below minimum {} at {}", actual, minimum, location)
 }
 Self::JoinCreationFailed { join_type, edge_index, reason } => {
 write!(f, "failed to create {:?} join at edge {}: {}", join_type, edge_index, reason)
 }
 Self::InvalidVariableThickness { face_index, thickness, reason } => {
 write!(f, "invalid thickness {} for face {}: {}", thickness, face_index, reason)
 }
 Self::QualityCheckFailed { metric, value, threshold } => {
 write!(f, "quality check '{}' failed: {} exceeds threshold {}", metric, value, threshold)
 }
 Self::RecoveryFailed { attempts, last_error } => {
 write!(f, "recovery failed after {} attempts: {}", attempts, last_error)
 }
 }
 }
}

impl std::error::Error for OffsetError {}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Join Types
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Join type for offset operations at edges.
///
/// Determines how adjacent offset faces are connected at their boundaries.
/// Analogous to OCCT `rcad_kernel::BRepOffset_Mode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum JoinType {
 /// Sharp corners at edge intersections.
 ///
 /// The offset surfaces extend until they intersect, creating sharp corners.
 /// This is the default mode and works well for mechanical parts.
 #[default]
 Intersection,

 /// Round corners using fillet arcs at edges.
 ///
 /// Creates cylindrical transition surfaces at edges with radius equal to
 /// the offset distance. Suitable for smooth, organic shapes.
 Arc,

 /// Smooth transitions between adjacent faces.
 ///
 /// Creates tangent-continuous transitions at edges where the adjacent
 /// faces have similar normals. Falls back to intersection join when
 /// the angle between faces is too large.
 Tangent,
}

impl JoinType {
 /// Returns true if this join type requires additional geometry creation.
 pub fn requires_join_geometry(&self) -> bool {
 matches!(self, JoinType::Arc | JoinType::Tangent)
 }

 /// Returns a string representation of the join type.
 pub fn as_str(&self) -> &'static str {
 match self {
 JoinType::Intersection => "intersection",
 JoinType::Arc => "arc",
 JoinType::Tangent => "tangent",
 }
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Variable Thickness
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Thickness specification for a single face in variable-thickness offset.
#[derive(Debug, Clone)]
pub struct FaceThickness {
 /// Index of the face in the shell.
 pub face_index: usize,
 /// Offset distance for this face.
 pub thickness: f64,
}

/// Specification for variable thickness offset.
///
/// Allows different offset distances for different faces, with smooth
/// transitions between regions of different thickness.
#[derive(Debug, Clone, Default)]
pub struct VariableThickness {
 /// Per-face thickness specifications.
 pub face_thicknesses: Vec<FaceThickness>,
 /// Default thickness for faces not explicitly specified.
 pub default_thickness: f64,
 /// Transition zone width for smoothing between different thicknesses.
 /// Set to 0.0 for sharp transitions.
 pub transition_width: f64,
 /// Enable thickness interpolation across faces.
 pub interpolate: bool,
}

impl VariableThickness {
 /// Create a new variable thickness specification with a default thickness.
 pub fn new(default: f64) -> Self {
 Self {
 face_thicknesses: Vec::new(),
 default_thickness: default,
 transition_width: 0.0,
 interpolate: false,
 }
 }

 /// Add a per-face thickness specification.
 pub fn with_face(mut self, face_index: usize, thickness: f64) -> Self {
 self.face_thicknesses.push(FaceThickness { face_index, thickness });
 self
 }

 /// Set the transition zone width for smoothing.
 pub fn with_transition(mut self, width: f64) -> Self {
 self.transition_width = width;
 self
 }

 /// Enable thickness interpolation.
 pub fn with_interpolation(mut self, interpolate: bool) -> Self {
 self.interpolate = interpolate;
 self
 }

 /// Get the thickness for a specific face.
 pub fn thickness_for_face(&self, face_index: usize) -> f64 {
 self.face_thicknesses
 .iter()
 .find(|ft| ft.face_index == face_index)
 .map(|ft| ft.thickness)
 .unwrap_or(self.default_thickness)
 }

 /// Validate the thickness specification.
 pub fn validate(&self, face_count: usize) -> Result<(), OffsetError> {
 for ft in &self.face_thicknesses {
 if ft.face_index >= face_count {
 return Err(OffsetError::InvalidVariableThickness {
 face_index: ft.face_index,
 thickness: ft.thickness,
 reason: format!("face index {} out of range (0..{})", ft.face_index, face_count),
 });
 }
 if ft.thickness.abs() < TOLERANCE_LEN_MIN {
 return Err(OffsetError::InvalidVariableThickness {
 face_index: ft.face_index,
 thickness: ft.thickness,
 reason: "thickness cannot be zero".to_string(),
 });
 }
 }
 if self.default_thickness.abs() < TOLERANCE_LEN_MIN {
 return Err(OffsetError::InvalidVariableThickness {
 face_index: 0,
 thickness: self.default_thickness,
 reason: "default thickness cannot be zero".to_string(),
 });
 }
 Ok(())
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Self-Intersection Handling
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Configuration for self-intersection handling.
#[derive(Debug, Clone)]
pub struct SelfIntersectionConfig {
 /// Whether to detect self-intersections.
 pub detect: bool,
 /// Whether to attempt automatic repair.
 pub auto_repair: bool,
 /// Maximum number of repair attempts.
 pub max_repair_attempts: usize,
 /// Factor to reduce offset distance during repair (0.0 to 1.0).
 pub reduction_factor: f64,
 /// Minimum offset distance to try during repair.
 pub min_offset_distance: f64,
 /// Whether to allow partial results after repair failure.
 pub allow_partial_results: bool,
}

impl Default for SelfIntersectionConfig {
 fn default() -> Self {
 Self {
 detect: true,
 auto_repair: false,
 max_repair_attempts: 5,
 reduction_factor: 0.8,
 min_offset_distance: TOLERANCE_MESH_LEGACY,
 allow_partial_results: false,
 }
 }
}

/// Result of self-intersection detection.
#[derive(Debug, Clone)]
pub struct SelfIntersectionResult {
 /// Whether self-intersection was detected.
 pub has_intersection: bool,
 /// Pairs of faces that intersect.
 pub intersecting_pairs: Vec<(usize, usize)>,
 /// Estimated minimum distance before self-intersection.
 pub min_safe_distance: Option<f64>,
 /// Description of the intersection.
 pub description: String,
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Offset Quality Analysis
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Quality metrics for an offset result.
#[derive(Debug, Clone, Default)]
pub struct OffsetQuality {
 /// Minimum wall thickness found.
 pub min_wall_thickness: f64,
 /// Maximum deviation from expected offset distance.
 pub max_deviation: f64,
 /// Number of degenerate edges created.
 pub degenerate_edge_count: usize,
 /// Number of self-intersections detected.
 pub self_intersection_count: usize,
 /// Face area ratio (offset/original).
 pub face_area_ratio: f64,
 /// Edge length ratio (offset/original).
 pub edge_length_ratio: f64,
 /// Whether the result is valid.
 pub is_valid: bool,
 /// Warnings generated during analysis.
 pub warnings: Vec<String>,
}

impl OffsetQuality {
 /// Check if the quality meets the given thresholds.
 pub fn check_thresholds(&self, config: &QualityThresholds) -> Result<(), OffsetError> {
 if self.min_wall_thickness < config.min_wall_thickness {
 return Err(OffsetError::QualityCheckFailed {
 metric: "min_wall_thickness".to_string(),
 value: self.min_wall_thickness,
 threshold: config.min_wall_thickness,
 });
 }
 if self.max_deviation > config.max_deviation {
 return Err(OffsetError::QualityCheckFailed {
 metric: "max_deviation".to_string(),
 value: self.max_deviation,
 threshold: config.max_deviation,
 });
 }
 if self.self_intersection_count > 0 && !config.allow_self_intersection {
 return Err(OffsetError::QualityCheckFailed {
 metric: "self_intersection_count".to_string(),
 value: self.self_intersection_count as f64,
 threshold: 0.0,
 });
 }
 Ok(())
 }
}

/// Quality thresholds for offset validation.
#[derive(Debug, Clone)]
pub struct QualityThresholds {
 /// Minimum acceptable wall thickness.
 pub min_wall_thickness: f64,
 /// Maximum acceptable deviation from expected offset.
 pub max_deviation: f64,
 /// Whether to allow self-intersections.
 pub allow_self_intersection: bool,
 /// Maximum acceptable ratio of degenerate edges.
 pub max_degenerate_ratio: f64,
}

impl Default for QualityThresholds {
 fn default() -> Self {
 Self {
 min_wall_thickness: TOLERANCE_MESH_LEGACY,
 max_deviation: TOLERANCE_RETRY_LADDER_COARSE,
 allow_self_intersection: false,
 max_degenerate_ratio: 0.1,
 }
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Options
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Options for offset operations.
#[derive(Debug, Clone)]
pub struct OffsetOptions {
 /// Offset distance. Positive = outward, negative = inward.
 pub distance: f64,
 /// Tolerance for geometric computations.
 pub tolerance: f64,
 /// Whether to check for self-intersection after offset.
 pub check_self_intersection: bool,
 /// Whether to attempt to repair self-intersections by reducing offset distance.
 pub auto_repair: bool,
 /// Minimum feature size to preserve (affects vertex handling).
 pub min_feature_size: f64,
 /// Join type for edge transitions.
 pub join_type: JoinType,
 /// Variable thickness specification (optional).
 pub variable_thickness: Option<VariableThickness>,
 /// Self-intersection handling configuration.
 pub self_intersection_config: SelfIntersectionConfig,
 /// Quality thresholds for validation.
 pub quality_thresholds: QualityThresholds,
 /// Approximation tolerance for curved surfaces.
 pub approximation_tolerance: f64,
 /// Whether to perform wall thickness checking.
 pub check_wall_thickness: bool,
 /// Minimum wall thickness requirement.
 pub min_wall_thickness: f64,
}

impl Default for OffsetOptions {
 fn default() -> Self {
 Self {
 distance: 1.0,
 tolerance: TOLERANCE_ABS,
 check_self_intersection: true,
 auto_repair: false,
 min_feature_size: TOLERANCE_MESH_LEGACY,
 join_type: JoinType::default(),
 variable_thickness: None,
 self_intersection_config: SelfIntersectionConfig::default(),
 quality_thresholds: QualityThresholds::default(),
 approximation_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
 check_wall_thickness: false,
 min_wall_thickness: TOLERANCE_MESH_LEGACY,
 }
 }
}

impl OffsetOptions {
 /// Create options with a given distance.
 pub fn new(distance: f64) -> Self {
 Self {
 distance,
 ..Default::default()
 }
 }

 /// Set tolerance.
 pub fn with_tolerance(mut self, tol: f64) -> Self {
 self.tolerance = tol;
 self
 }

 /// Enable or disable self-intersection checking.
 pub fn with_self_intersection_check(mut self, check: bool) -> Self {
 self.check_self_intersection = check;
 self
 }

 /// Enable or disable auto-repair of self-intersections.
 pub fn with_auto_repair(mut self, repair: bool) -> Self {
 self.auto_repair = repair;
 self
 }

 /// Set the join type for edge transitions.
 pub fn with_join_type(mut self, join_type: JoinType) -> Self {
 self.join_type = join_type;
 self
 }

 /// Set variable thickness specification.
 pub fn with_variable_thickness(mut self, thickness: VariableThickness) -> Self {
 self.variable_thickness = Some(thickness);
 self
 }

 /// Set self-intersection handling configuration.
 pub fn with_self_intersection_config(mut self, config: SelfIntersectionConfig) -> Self {
 self.self_intersection_config = config;
 self
 }

 /// Set quality thresholds.
 pub fn with_quality_thresholds(mut self, thresholds: QualityThresholds) -> Self {
 self.quality_thresholds = thresholds;
 self
 }

 /// Set approximation tolerance for curved surfaces.
 pub fn with_approximation_tolerance(mut self, tol: f64) -> Self {
 self.approximation_tolerance = tol;
 self
 }

 /// Enable wall thickness checking with a minimum value.
 pub fn with_wall_thickness_check(mut self, min_thickness: f64) -> Self {
 self.check_wall_thickness = true;
 self.min_wall_thickness = min_thickness;
 self
 }

 /// Get the effective distance for a specific face (for variable thickness).
 pub fn effective_distance_for_face(&self, face_index: usize) -> f64 {
 if let Some(ref vt) = self.variable_thickness {
 vt.thickness_for_face(face_index)
 } else {
 self.distance
 }
 }
}

/// Result of an offset operation.
#[derive(Debug, Clone)]
pub struct OffsetResult {
 /// The resulting rcad_kernel::BRep.
 pub brep: rcad_kernel::BRep,
 /// Number of offset faces created.
 pub offset_faces: usize,
 /// Number of lateral faces created (for hollow operations).
 pub lateral_faces: usize,
 /// Number of join faces created (for arc/tangent joins).
 pub join_faces: usize,
 /// Whether self-intersection was detected.
 pub self_intersection: bool,
 /// Quality analysis of the result.
 pub quality: OffsetQuality,
 /// Warnings generated during the operation.
 pub warnings: Vec<String>,
 /// The effective offset distance used (may differ from requested if auto-repair occurred).
 pub effective_distance: f64,
 /// Number of repair attempts made (0 if no repair was needed).
 pub repair_attempts: usize,
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Surface Offset
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Compute an offset surface at distance `d` along the normal direction.
///
/// Returns `None` if the offset would create a degenerate surface
/// (e.g., sphere with negative radius).
pub fn offset_surface(surf: &Surface3, d: f64) -> Option<Surface3> {
 match surf {
 Surface3::Plane(p) => {
 // Plane offset: translate along normal
 Some(Surface3::Plane(Plane {
 origin: p.origin + p.normal * d,
 normal: p.normal,
 }))
 }

 Surface3::Sphere(s) => {
 // Sphere offset: adjust radius
 let new_radius = s.radius + d;
 if new_radius <= 0.0 {
 return None;
 }
 Some(Surface3::Sphere(SphericalSurface::new(
 s.center, s.axis, new_radius,
 )))
 }

 Surface3::Cylinder(c) => {
 // Cylinder offset: adjust radius
 let new_radius = c.radius + d;
 if new_radius <= 0.0 {
 return None;
 }
 Some(Surface3::Cylinder(CylindricalSurface {
 origin: c.origin,
 axis: c.axis,
 radius: new_radius,
 ref_dir: c.ref_dir,
 }))
 }

 Surface3::Cone(c) => {
 // Cone offset: the parallel surface to a cone is another cone with the
 // same half-angle but its true apex shifts by d/sin(α) along the axis
 // toward the smaller-radius end (opposite to the cone axis direction).
 let sin_a = c.half_angle_rad.sin();

 // Near-cylinder case
 if sin_a.abs() <= TOLERANCE_LINEAR_ULTRA_STRICT {
 let new_radius = c.radius + d;
 if new_radius <= 0.0 {
 return None;
 }
 return Some(Surface3::Cone(ConicalSurface {
 radius: new_radius,
 ..*c
 }));
 }

 let tan_a = c.half_angle_rad.tan();
 let axis_dir = c.axis.normalize_or(DVec3::Y);

 // The true apex (where radius=0) of the original cone
 let true_apex = c.apex - axis_dir * (c.radius / tan_a);

 // Offset shifts the true apex by d/sin(α) in the outward direction
 // (toward smaller end, opposite to cone axis)
 let new_true_apex = true_apex - axis_dir * (d / sin_a);

 // Represent with reference apex at the new true apex, radius=0
 Some(Surface3::Cone(ConicalSurface {
 apex: new_true_apex,
 axis: c.axis,
 radius: 0.0,
 half_angle_rad: c.half_angle_rad,
 }))
 }

 Surface3::Torus(t) => {
 // Torus offset: adjust minor radius
 let new_minor = t.minor_radius + d;
 if new_minor <= 0.0 {
 return None;
 }
 // Check for self-intersection: minor radius > major radius
 // The offset surface is valid but may be self-intersecting
 Some(Surface3::Torus(ToroidalSurface {
 center: t.center,
 axis: t.axis,
 major_radius: t.major_radius,
 minor_radius: new_minor,
 }))
 }

 // For parametric surfaces, use the generic OffsetSurface wrapper
 Surface3::BSpline(_)
 | Surface3::Bezier(_)
 | Surface3::TriBezier(_)
 | Surface3::LinearExtrusion(_)
 | Surface3::Revolution(_)
 | Surface3::Ruled(_)
 | Surface3::Coons(_)
 | Surface3::Ellipsoid(_)
 | Surface3::Helicoid(_)
 | Surface3::Pipe(_) => {
 Some(Surface3::Offset(OffsetSurface {
 basis: Box::new(surf.clone()),
 offset_distance: d,
 }))
 }

 // Trimmed surface: offset the basis
 Surface3::Trimmed(t) => {
 let offset_basis = offset_surface(&t.basis, d)?;
 Some(Surface3::Trimmed(rcad_kernel::geom::TrimmedSurface {
 basis: Box::new(offset_basis),
 trim: t.trim,
 }))
 }

 // Offset surface: compound the offsets
 Surface3::Offset(o) => {
 Some(Surface3::Offset(OffsetSurface {
 basis: o.basis.clone(),
 offset_distance: o.offset_distance + d,
 }))
 }
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Analytical Intersection Curves for Offset Surfaces
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Result of intersecting two offset surfaces.
#[derive(Debug, Clone)]
pub enum OffsetIntersectionCurve {
 /// No intersection (surfaces are disjoint or nested).
 NoIntersection,
 /// Single point intersection (tangent case).
 TangentPoint(DVec3),
 /// Single line intersection.
 Line(Line3),
 /// Two line intersections (parallel cylinder case).
 TwoLines(Line3, Line3),
 /// Single circle intersection.
 Circle(Circle3),
 /// Tangent circle (surfaces touch along a circle).
 TangentCircle(Circle3),
 /// Two circle intersections.
 TwoCircles(Circle3, Circle3),
 /// Ellipse intersection.
 Ellipse(Ellipse3),
 /// Two ellipse intersections.
 TwoEllipses(Ellipse3, Ellipse3),
 /// Parabola intersection (plane-cone edge case).
 Parabola(Parabola3),
 /// Hyperbola intersection (plane-cone edge case).
 Hyperbola(Hyperbola3),
 /// Intersection requires numerical approximation.
 /// Contains a sampled polyline approximation.
 Numerical(Vec<DVec3>),
 /// General intersection (no simple analytic form).
 General,
 /// Surfaces are coincident (infinite intersection).
 Coincident,
}

/// Intersect two offset planes.
///
/// A plane offset by distance `d1` is simply another plane translated along the normal.
/// The intersection of two planes is always a line (if not parallel).
pub fn intersect_offset_plane_plane(
 plane1: &Plane,
 plane2: &Plane,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 // Compute the offset planes
 let offset_plane1 = Plane {
 origin: plane1.origin + plane1.normal * d1,
 normal: plane1.normal,
 };
 let offset_plane2 = Plane {
 origin: plane2.origin + plane2.normal * d2,
 normal: plane2.normal,
 };

 // Use existing plane-plane intersection
 match crate::inttools::plane_plane::intersect_plane_plane(&offset_plane1, &offset_plane2) {
 crate::inttools::plane_plane::PlanePlaneResult::Line(line) => {
 OffsetIntersectionCurve::Line(line)
 }
 crate::inttools::plane_plane::PlanePlaneResult::Parallel => {
 OffsetIntersectionCurve::NoIntersection
 }
 crate::inttools::plane_plane::PlanePlaneResult::Coincident => {
 OffsetIntersectionCurve::Coincident
 }
 }
}

/// Intersect two offset spheres.
///
/// Two spheres offset by distances `d1` and `d2` are equivalent to spheres
/// with adjusted radii. The intersection is a circle (if any).
pub fn intersect_offset_sphere_sphere(
 sphere1: &SphericalSurface,
 sphere2: &SphericalSurface,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 // Compute effective radii after offset
 let r1 = sphere1.radius + d1;
 let r2 = sphere2.radius + d2;

 // Check for degenerate cases
 if r1 <= 0.0 || r2 <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Distance between centers
 let centers_dist = (sphere2.center - sphere1.center).length();

 // No intersection cases
 if centers_dist > r1 + r2 + TOLERANCE_ABS {
 return OffsetIntersectionCurve::NoIntersection;
 }
 if centers_dist < (r1 - r2).abs() - TOLERANCE_ABS {
 return OffsetIntersectionCurve::NoIntersection;
 }
 if centers_dist < TOLERANCE_ABS {
 // Concentric spheres
 if (r1 - r2).abs() < TOLERANCE_ABS {
 return OffsetIntersectionCurve::Coincident;
 }
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Tangent case
 if (centers_dist - (r1 + r2)).abs() < TOLERANCE_ABS {
 // External tangent - single point
 let t = r1 / (r1 + r2);
 let point = sphere1.center + (sphere2.center - sphere1.center) * t;
 // Return a degenerate circle (radius 0)
 return OffsetIntersectionCurve::Circle(Circle3::new(point, (sphere2.center - sphere1.center).normalize(), 0.0,
 ));
 }
 if (centers_dist - (r1 - r2).abs()).abs() < TOLERANCE_ABS {
 // Internal tangent - single point
 let dir = (sphere2.center - sphere1.center).normalize();
 let point = if r1 > r2 {
 sphere1.center + dir * r1
 } else {
 sphere1.center - dir * r1
 };
 return OffsetIntersectionCurve::Circle(Circle3::new(point, dir, 0.0,
 ));
 }

 // Two intersecting spheres produce a circle
 // The circle lies in a plane perpendicular to the line connecting centers
 // Distance from center1 to the plane
 let h = (r1 * r1 - r2 * r2 + centers_dist * centers_dist) / (2.0 * centers_dist);
 // Radius of the intersection circle
 let circle_radius_sq = r1 * r1 - h * h;
 if circle_radius_sq < 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let circle_radius = circle_radius_sq.sqrt();

 // Center of the intersection circle
 let dir = (sphere2.center - sphere1.center).normalize();
 let circle_center = sphere1.center + dir * h;

 OffsetIntersectionCurve::Circle(Circle3::new(circle_center, dir, circle_radius,
 ))
}

/// Intersect two offset cylinders.
///
/// The intersection depends on the relative orientation of the cylinders:
/// - Parallel axes: 0, 1, or 2 generator lines
/// - Perpendicular axes with intersecting axes: circles or ellipses (Steinmetz)
/// - General case: numerical approximation
pub fn intersect_offset_cylinder_cylinder(
 cyl1: &CylindricalSurface,
 cyl2: &CylindricalSurface,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 // Compute effective radii after offset
 let r1 = cyl1.radius + d1;
 let r2 = cyl2.radius + d2;

 // Check for degenerate cases
 if r1 <= 0.0 || r2 <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Create offset cylinders
 let offset_cyl1 = CylindricalSurface {
 origin: cyl1.origin,
 axis: cyl1.axis,
 radius: r1,
 ref_dir: cyl1.ref_dir,
 };
 let offset_cyl2 = CylindricalSurface {
 origin: cyl2.origin,
 axis: cyl2.axis,
 radius: r2,
 ref_dir: cyl2.ref_dir,
 };

 // Use existing cylinder-cylinder intersection
 match crate::inttools::cylinder_cylinder::intersect_cylinder_cylinder(&offset_cyl1, &offset_cyl2) {
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::NoIntersection => {
 OffsetIntersectionCurve::NoIntersection
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::Coaxial => {
 OffsetIntersectionCurve::Coincident
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::OneGeneratorLine(line) => {
 OffsetIntersectionCurve::Line(line)
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::TwoGeneratorLines(l1, l2) => {
 OffsetIntersectionCurve::TwoLines(l1, l2)
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::TwoCircles(c1, c2) => {
 OffsetIntersectionCurve::TwoCircles(c1, c2)
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::TwoEllipses(e1, e2) => {
 OffsetIntersectionCurve::TwoEllipses(e1, e2)
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::SkewQuartic(branches) => {
 // Use analytical quartic solver output (Ferrari method) instead of marching.
 let mut pts = Vec::new();
 for branch in branches { pts.extend(branch); }
 if pts.is_empty() {
 OffsetIntersectionCurve::NoIntersection
 } else {
 OffsetIntersectionCurve::Numerical(pts)
 }
 }
 crate::inttools::cylinder_cylinder::CylinderCylinderResult::PerpendicularOffsetCurves { .. }
 | crate::inttools::cylinder_cylinder::CylinderCylinderResult::General => {
 // Fall back to numerical approximation
 intersect_cylinders_numerical(&offset_cyl1, &offset_cyl2, d1, d2)
 }
 }
}

/// Intersect an offset plane with an offset cylinder.
///
/// The intersection is:
/// - A circle if the plane is perpendicular to the cylinder axis
/// - One or two lines if the plane is parallel to the cylinder axis
/// - An ellipse for oblique intersections
pub fn intersect_offset_plane_cylinder(
 plane: &Plane,
 cyl: &CylindricalSurface,
 d_plane: f64,
 d_cyl: f64,
) -> OffsetIntersectionCurve {
 // Compute the offset plane
 let offset_plane = Plane {
 origin: plane.origin + plane.normal * d_plane,
 normal: plane.normal,
 };

 // Compute effective cylinder radius after offset
 let r = cyl.radius + d_cyl;
 if r <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Create offset cylinder
 let offset_cyl = CylindricalSurface {
 origin: cyl.origin,
 axis: cyl.axis,
 radius: r,
 ref_dir: cyl.ref_dir,
 };

 // Use existing plane-cylinder intersection
 match crate::inttools::plane_cylinder::intersect_plane_cylinder(&offset_plane, &offset_cyl) {
 crate::inttools::plane_cylinder::PlaneCylinderResult::NoIntersection => {
 OffsetIntersectionCurve::NoIntersection
 }
 crate::inttools::plane_cylinder::PlaneCylinderResult::TangentLine(line) => {
 OffsetIntersectionCurve::Line(line)
 }
 crate::inttools::plane_cylinder::PlaneCylinderResult::TwoLines(l1, l2) => {
 OffsetIntersectionCurve::TwoLines(l1, l2)
 }
 crate::inttools::plane_cylinder::PlaneCylinderResult::Circle(circle) => {
 OffsetIntersectionCurve::Circle(circle)
 }
 crate::inttools::plane_cylinder::PlaneCylinderResult::Ellipse(ellipse) => {
 OffsetIntersectionCurve::Ellipse(ellipse)
 }
 }
}

/// Intersect an offset plane with an offset sphere.
///
/// The intersection is a circle (or tangent point, or none).
pub fn intersect_offset_plane_sphere(
 plane: &Plane,
 sphere: &SphericalSurface,
 d_plane: f64,
 d_sphere: f64,
) -> OffsetIntersectionCurve {
 // Compute the offset plane
 let offset_plane = Plane {
 origin: plane.origin + plane.normal * d_plane,
 normal: plane.normal,
 };

 // Compute effective sphere radius after offset
 let r = sphere.radius + d_sphere;
 if r <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Create offset sphere
 let offset_sphere = SphericalSurface::new(sphere.center, sphere.axis, r);
 // Use existing plane-sphere intersection
 match crate::inttools::plane_sphere::intersect_plane_sphere(&offset_plane, &offset_sphere) {
 crate::inttools::plane_sphere::PlaneSphereResult::NoIntersection => {
 OffsetIntersectionCurve::NoIntersection
 }
 crate::inttools::plane_sphere::PlaneSphereResult::TangentPoint(point) => {
 // Return a degenerate circle
 OffsetIntersectionCurve::Circle(Circle3::new(point, plane.normal, 0.0,
 ))
 }
 crate::inttools::plane_sphere::PlaneSphereResult::Circle(circle) => {
 OffsetIntersectionCurve::Circle(circle)
 }
 }
}

/// Intersect an offset cylinder with an offset sphere.
///
/// For the axis-aligned case (sphere center on cylinder axis), the result is
/// one or two circles. For other cases, the intersection is a quartic curve
/// requiring numerical approximation.
pub fn intersect_offset_cylinder_sphere(
 cyl: &CylindricalSurface,
 sphere: &SphericalSurface,
 d_cyl: f64,
 d_sphere: f64,
) -> OffsetIntersectionCurve {
 // Compute effective radii after offset
 let r_cyl = cyl.radius + d_cyl;
 let r_sphere = sphere.radius + d_sphere;

 if r_cyl <= 0.0 || r_sphere <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Create offset surfaces
 let offset_cyl = CylindricalSurface {
 origin: cyl.origin,
 axis: cyl.axis,
 radius: r_cyl,
 ref_dir: cyl.ref_dir,
 };
 let offset_sphere = SphericalSurface::new(sphere.center, sphere.axis, r_sphere);

 // Use existing sphere-cylinder intersection
 match crate::inttools::sphere_cylinder::intersect_sphere_cylinder(&offset_sphere, &offset_cyl) {
 crate::inttools::sphere_cylinder::SphereCylinderResult::NoIntersection => {
 OffsetIntersectionCurve::NoIntersection
 }
 crate::inttools::sphere_cylinder::SphereCylinderResult::TangentCircle(circle) => {
 OffsetIntersectionCurve::Circle(circle)
 }
 crate::inttools::sphere_cylinder::SphereCylinderResult::TwoCircles(c1, c2) => {
 OffsetIntersectionCurve::TwoCircles(c1, c2)
 }
 crate::inttools::sphere_cylinder::SphereCylinderResult::SkewQuartic(branches) => {
 // Use analytical quartic solver output instead of marching.
 let mut pts = Vec::new();
 for branch in branches { pts.extend(branch); }
 if pts.is_empty() {
 OffsetIntersectionCurve::NoIntersection
 } else {
 OffsetIntersectionCurve::Numerical(pts)
 }
 }
 crate::inttools::sphere_cylinder::SphereCylinderResult::General => {
 // Fall back to numerical approximation
 intersect_cylinder_sphere_numerical(&offset_cyl, &offset_sphere, d_cyl, d_sphere)
 }
 }
}

// ── Cone offset helper ─────────────────────────────────────────────────
fn offset_conical_surface(c: &ConicalSurface, d: f64) -> Option<ConicalSurface> {
 let sin_a = c.half_angle_rad.sin();

 // Near-cylinder: just offset the radius
 if sin_a.abs() <= TOLERANCE_LINEAR_ULTRA_STRICT {
 let new_radius = c.radius + d;
 if new_radius <= 0.0 {
 return None;
 }
 return Some(ConicalSurface {
 radius: new_radius,
 ..*c
 });
 }

 let tan_a = c.half_angle_rad.tan();
 let axis_dir = c.axis.normalize_or(DVec3::Y);

 // The true apex (where radius=0) of the original cone
 let true_apex = c.apex - axis_dir * (c.radius / tan_a);

 // Offsetting the conical surface along its outward normal shifts the true
 // apex by d/sin(α) toward the smaller-radius end (opposite to the cone
 // axis direction, since the outward normal's axial component always points
 // toward the smaller end of the frustum).
 let new_true_apex = true_apex - axis_dir * (d / sin_a);

 // Represent the offset cone with the reference apex at the new true apex.
 // The radius at the reference apex depends only on the distance from the
 // true apex along the axis (which is zero here), so radius=0.
 Some(ConicalSurface {
 apex: new_true_apex,
 axis: c.axis,
 radius: 0.0,
 half_angle_rad: c.half_angle_rad,
 })
}

// ── Phase 2: New offset handler functions ──────────────────────────────

pub fn intersect_offset_cone_cone(
 cone1: &ConicalSurface,
 cone2: &ConicalSurface,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 let off_c1 = match offset_conical_surface(cone1, d1) {
 Some(c) => c,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 let off_c2 = match offset_conical_surface(cone2, d2) {
 Some(c) => c,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 match intersect_cone_cone(&off_c1, &off_c2) {
 ConeConeResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 ConeConeResult::Coaxial => OffsetIntersectionCurve::Coincident,
 ConeConeResult::CoaxialCircle(c) => OffsetIntersectionCurve::Circle(c),
 ConeConeResult::CoaxialPoint(p) => OffsetIntersectionCurve::TangentPoint(p),
 ConeConeResult::SkewQuartic(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 ConeConeResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_cylinder_cone(
 cyl: &CylindricalSurface,
 cone: &ConicalSurface,
 d_cyl: f64,
 d_cone: f64,
) -> OffsetIntersectionCurve {
 let r = cyl.radius + d_cyl;
 if r <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_cyl = CylindricalSurface {
 origin: cyl.origin,
 axis: cyl.axis,
 radius: r,
 ref_dir: cyl.ref_dir,
 };
 let offset_cone = match offset_conical_surface(cone, d_cone) {
 Some(c) => c,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 match intersect_cylinder_cone(&offset_cyl, &offset_cone) {
 CylinderConeResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 CylinderConeResult::CoaxialCircle(c) => OffsetIntersectionCurve::Circle(c),
 CylinderConeResult::CoaxialTwoCircles(c1, _c2) => OffsetIntersectionCurve::Circle(c1),
 CylinderConeResult::ParallelOffsetPolyline(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 CylinderConeResult::SkewQuartic(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 CylinderConeResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_plane_cone(
 plane: &Plane,
 cone: &ConicalSurface,
 d_plane: f64,
 d_cone: f64,
) -> OffsetIntersectionCurve {
 let offset_plane = Plane {
 origin: plane.origin + plane.normal * d_plane,
 normal: plane.normal,
 };
 let offset_cone = match offset_conical_surface(cone, d_cone) {
 Some(c) => c,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 match intersect_plane_cone(&offset_plane, &offset_cone) {
 PlaneConicalResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 PlaneConicalResult::Point(p) => OffsetIntersectionCurve::TangentPoint(p),
 PlaneConicalResult::SingleLine(l) => OffsetIntersectionCurve::Line(l),
 PlaneConicalResult::TwoLines(l1, l2) => OffsetIntersectionCurve::TwoLines(l1, l2),
 PlaneConicalResult::Circle(c) => OffsetIntersectionCurve::Circle(c),
 PlaneConicalResult::Ellipse(e) => OffsetIntersectionCurve::Ellipse(e),
 PlaneConicalResult::Parabola(p) => OffsetIntersectionCurve::Parabola(p),
 PlaneConicalResult::Hyperbola(h) => OffsetIntersectionCurve::Hyperbola(h),
 }
}

pub fn intersect_offset_sphere_cone(
 sphere: &SphericalSurface,
 cone: &ConicalSurface,
 d_sphere: f64,
 d_cone: f64,
) -> OffsetIntersectionCurve {
 let r_sphere = sphere.radius + d_sphere;
 if r_sphere <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_sphere = SphericalSurface::new(sphere.center, sphere.axis, r_sphere);
 let offset_cone = match offset_conical_surface(cone, d_cone) {
 Some(c) => c,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 match intersect_sphere_cone_with_tolerance(&offset_sphere, &offset_cone, TOLERANCE_ABS) {
 SphereConeResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 SphereConeResult::SingleCircle(c) => OffsetIntersectionCurve::Circle(c),
 SphereConeResult::TwoCircles(c1, c2) => OffsetIntersectionCurve::TwoCircles(c1, c2),
 SphereConeResult::TangentPoint(p) => OffsetIntersectionCurve::TangentPoint(p),
 SphereConeResult::Polyline(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 SphereConeResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_torus_torus(
 t1: &ToroidalSurface,
 t2: &ToroidalSurface,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 let r1 = t1.minor_radius + d1;
 if r1 <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let r2 = t2.minor_radius + d2;
 if r2 <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let off_t1 = ToroidalSurface {
 center: t1.center,
 axis: t1.axis,
 major_radius: t1.major_radius,
 minor_radius: r1,
 };
 let off_t2 = ToroidalSurface {
 center: t2.center,
 axis: t2.axis,
 major_radius: t2.major_radius,
 minor_radius: r2,
 };
 match intersect_torus_torus(&off_t1, &off_t2) {
 TorusTorusResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 TorusTorusResult::SingleCircle(c) => OffsetIntersectionCurve::Circle(c),
 TorusTorusResult::TwoCircles(c1, c2) => OffsetIntersectionCurve::TwoCircles(c1, c2),
 TorusTorusResult::TangentCircle(c) => OffsetIntersectionCurve::TangentCircle(c),
 TorusTorusResult::Coaxial => OffsetIntersectionCurve::Coincident,
 TorusTorusResult::SkewQuartic(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 TorusTorusResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_cylinder_torus(
 cyl: &CylindricalSurface,
 torus: &ToroidalSurface,
 d_cyl: f64,
 d_torus: f64,
) -> OffsetIntersectionCurve {
 let r_cyl = cyl.radius + d_cyl;
 if r_cyl <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_cyl = CylindricalSurface {
 origin: cyl.origin,
 axis: cyl.axis,
 radius: r_cyl,
 ref_dir: cyl.ref_dir,
 };
 let r_minor = torus.minor_radius + d_torus;
 if r_minor <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_torus = ToroidalSurface {
 center: torus.center,
 axis: torus.axis,
 major_radius: torus.major_radius,
 minor_radius: r_minor,
 };
 match intersect_cylinder_torus(&offset_cyl, &offset_torus) {
 CylinderTorusResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 CylinderTorusResult::TangentCircle(c) => OffsetIntersectionCurve::TangentCircle(c),
 CylinderTorusResult::TwoCircles(c1, c2) => OffsetIntersectionCurve::TwoCircles(c1, c2),
 CylinderTorusResult::SkewQuartic(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 CylinderTorusResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_torus_cone(
 torus: &ToroidalSurface,
 cone: &ConicalSurface,
 d_torus: f64,
 d_cone: f64,
) -> OffsetIntersectionCurve {
 let r_minor = torus.minor_radius + d_torus;
 if r_minor <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_torus = ToroidalSurface {
 center: torus.center,
 axis: torus.axis,
 major_radius: torus.major_radius,
 minor_radius: r_minor,
 };
 let offset_cone = match offset_conical_surface(cone, d_cone) {
 Some(c) => c,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 match intersect_torus_cone(&offset_torus, &offset_cone) {
 TorusConeResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 TorusConeResult::SingleCircle(c) => OffsetIntersectionCurve::Circle(c),
 TorusConeResult::TwoCircles(c1, c2) => OffsetIntersectionCurve::TwoCircles(c1, c2),
 TorusConeResult::TangentCircle(c) => OffsetIntersectionCurve::TangentCircle(c),
 TorusConeResult::SkewQuartic(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 TorusConeResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_plane_torus(
 plane: &Plane,
 torus: &ToroidalSurface,
 d_plane: f64,
 d_torus: f64,
) -> OffsetIntersectionCurve {
 let offset_plane = Plane {
 origin: plane.origin + plane.normal * d_plane,
 normal: plane.normal,
 };
 let r_minor = torus.minor_radius + d_torus;
 if r_minor <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_torus = ToroidalSurface {
 center: torus.center,
 axis: torus.axis,
 major_radius: torus.major_radius,
 minor_radius: r_minor,
 };
 match intersect_plane_torus(&offset_plane, &offset_torus) {
 PlaneTorusResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 PlaneTorusResult::TangentCircle(c) => OffsetIntersectionCurve::TangentCircle(c),
 PlaneTorusResult::TwoCircles(c1, c2) => OffsetIntersectionCurve::TwoCircles(c1, c2),
 PlaneTorusResult::SkewPolyline(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 PlaneTorusResult::General => OffsetIntersectionCurve::General,
 }
}

pub fn intersect_offset_sphere_torus(
 sphere: &SphericalSurface,
 torus: &ToroidalSurface,
 d_sphere: f64,
 d_torus: f64,
) -> OffsetIntersectionCurve {
 let r_sphere = sphere.radius + d_sphere;
 if r_sphere <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let r_minor = torus.minor_radius + d_torus;
 if r_minor <= 0.0 {
 return OffsetIntersectionCurve::NoIntersection;
 }
 let offset_sphere = SphericalSurface::new(sphere.center, sphere.axis, r_sphere);
 let offset_torus = ToroidalSurface {
 center: torus.center,
 axis: torus.axis,
 major_radius: torus.major_radius,
 minor_radius: r_minor,
 };
 match intersect_sphere_torus(&offset_sphere, &offset_torus) {
 SphereTorusResult::NoIntersection => OffsetIntersectionCurve::NoIntersection,
 SphereTorusResult::OneCircle(c) => OffsetIntersectionCurve::Circle(c),
 SphereTorusResult::TwoCircles(c1, c2) => OffsetIntersectionCurve::TwoCircles(c1, c2),
 SphereTorusResult::SkewPolyline(branches) => branches
 .into_iter()
 .next()
 .filter(|b| b.len() >= 2)
 .map(OffsetIntersectionCurve::Numerical)
 .unwrap_or(OffsetIntersectionCurve::General),
 SphereTorusResult::General => OffsetIntersectionCurve::General,
 }
}

/// Intersect two offset surfaces of arbitrary types.
///
/// Dispatches to the appropriate analytical handler based on surface types,
/// falling back to numerical approximation for unsupported combinations.
pub fn intersect_offset_surfaces(
 surf1: &Surface3,
 surf2: &Surface3,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 match (surf1, surf2) {
 (Surface3::Plane(p1), Surface3::Plane(p2)) => {
 intersect_offset_plane_plane(p1, p2, d1, d2)
 }
 (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
 intersect_offset_sphere_sphere(s1, s2, d1, d2)
 }
 (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
 intersect_offset_cylinder_cylinder(c1, c2, d1, d2)
 }
 (Surface3::Plane(p), Surface3::Cylinder(c))
 | (Surface3::Cylinder(c), Surface3::Plane(p)) => {
 intersect_offset_plane_cylinder(p, c, d1, d2)
 }
 (Surface3::Plane(p), Surface3::Sphere(s))
 | (Surface3::Sphere(s), Surface3::Plane(p)) => {
 intersect_offset_plane_sphere(p, s, d1, d2)
 }
 (Surface3::Cylinder(c), Surface3::Sphere(s))
 | (Surface3::Sphere(s), Surface3::Cylinder(c)) => {
 intersect_offset_cylinder_sphere(c, s, d1, d2)
 }
 (Surface3::Plane(p), Surface3::Cone(c))
 | (Surface3::Cone(c), Surface3::Plane(p)) => {
 intersect_offset_plane_cone(p, c, d1, d2)
 }
 (Surface3::Cylinder(c), Surface3::Cone(k))
 | (Surface3::Cone(k), Surface3::Cylinder(c)) => {
 intersect_offset_cylinder_cone(c, k, d1, d2)
 }
 (Surface3::Sphere(s), Surface3::Cone(c))
 | (Surface3::Cone(c), Surface3::Sphere(s)) => {
 intersect_offset_sphere_cone(s, c, d1, d2)
 }
 (Surface3::Cone(k1), Surface3::Cone(k2)) => {
 intersect_offset_cone_cone(k1, k2, d1, d2)
 }
 (Surface3::Plane(p), Surface3::Torus(t))
 | (Surface3::Torus(t), Surface3::Plane(p)) => {
 intersect_offset_plane_torus(p, t, d1, d2)
 }
 (Surface3::Sphere(s), Surface3::Torus(t))
 | (Surface3::Torus(t), Surface3::Sphere(s)) => {
 intersect_offset_sphere_torus(s, t, d1, d2)
 }
 (Surface3::Cylinder(c), Surface3::Torus(t))
 | (Surface3::Torus(t), Surface3::Cylinder(c)) => {
 intersect_offset_cylinder_torus(c, t, d1, d2)
 }
 (Surface3::Torus(t), Surface3::Cone(c))
 | (Surface3::Cone(c), Surface3::Torus(t)) => {
 intersect_offset_torus_cone(t, c, d1, d2)
 }
 (Surface3::Torus(t1), Surface3::Torus(t2)) => {
 intersect_offset_torus_torus(t1, t2, d1, d2)
 }
 // For other combinations, fall back to numerical approximation
 _ => intersect_surfaces_numerical(surf1, surf2, d1, d2),
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Numerical Intersection Fallbacks
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

#[inline]
fn offset_numeric_geom_floor(da: f64, db: f64) -> f64 {
 if !(da.is_finite() && db.is_finite()) {
 return TOLERANCE_ABS;
 }
 da.abs().max(db.abs()).max(TOLERANCE_ABS)
}

/// Numerical approximation of cylinder-cylinder intersection using marching.
fn intersect_cylinders_numerical(
 cyl1: &CylindricalSurface,
 cyl2: &CylindricalSurface,
 offset_d1: f64,
 offset_d2: f64,
) -> OffsetIntersectionCurve {
 let surf1 = Surface3::Cylinder(*cyl1);
 let surf2 = Surface3::Cylinder(*cyl2);

 // Use the existing marching algorithm
 let tol_floor = offset_numeric_geom_floor(offset_d1, offset_d2);
 let result = crate::inttools::intss::intersect_surfaces_with_density_tol(
 &surf1, &surf2, 64, tol_floor,
 );

 if result.curves.is_empty() {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Convert to polyline
 let mut points = Vec::new();
 for curve in &result.curves {
 if let crate::inttools::intss::SurfaceCurve::Polyline(pts) = &curve.curve_3d {
 points.extend(pts.iter().copied());
 }
 }

 if points.is_empty() {
 OffsetIntersectionCurve::NoIntersection
 } else {
 OffsetIntersectionCurve::Numerical(points)
 }
}

/// Numerical approximation of cylinder-sphere intersection using marching.
fn intersect_cylinder_sphere_numerical(
 cyl: &CylindricalSurface,
 sphere: &SphericalSurface,
 offset_d_cyl: f64,
 offset_d_sphere: f64,
) -> OffsetIntersectionCurve {
 let surf1 = Surface3::Cylinder(*cyl);
 let surf2 = Surface3::Sphere(*sphere);

 let tol_floor = offset_numeric_geom_floor(offset_d_cyl, offset_d_sphere);
 let result = crate::inttools::intss::intersect_surfaces_with_density_tol(
 &surf1, &surf2, 64, tol_floor,
 );

 if result.curves.is_empty() {
 return OffsetIntersectionCurve::NoIntersection;
 }

 let mut points = Vec::new();
 for curve in &result.curves {
 if let crate::inttools::intss::SurfaceCurve::Polyline(pts) = &curve.curve_3d {
 points.extend(pts.iter().copied());
 }
 }

 if points.is_empty() {
 OffsetIntersectionCurve::NoIntersection
 } else {
 OffsetIntersectionCurve::Numerical(points)
 }
}

/// Generic numerical intersection for arbitrary offset surfaces.
fn intersect_surfaces_numerical(
 surf1: &Surface3,
 surf2: &Surface3,
 d1: f64,
 d2: f64,
) -> OffsetIntersectionCurve {
 // Compute offset surfaces
 let offset_surf1 = match offset_surface(surf1, d1) {
 Some(s) => s,
 None => return OffsetIntersectionCurve::NoIntersection,
 };
 let offset_surf2 = match offset_surface(surf2, d2) {
 Some(s) => s,
 None => return OffsetIntersectionCurve::NoIntersection,
 };

 // Use existing intersection
 let tol_floor = offset_numeric_geom_floor(d1, d2);
 let result = crate::inttools::intss::intersect_surfaces_with_density_tol(
 &offset_surf1,
 &offset_surf2,
 64,
 tol_floor,
 );

 if result.curves.is_empty() {
 return OffsetIntersectionCurve::NoIntersection;
 }

 // Convert first curve to appropriate type
 for curve in &result.curves {
 match &curve.curve_3d {
 crate::inttools::intss::SurfaceCurve::Line(l) => {
 return OffsetIntersectionCurve::Line(*l);
 }
 crate::inttools::intss::SurfaceCurve::Circle(c) => {
 return OffsetIntersectionCurve::Circle(*c);
 }
 crate::inttools::intss::SurfaceCurve::Ellipse(e) => {
 return OffsetIntersectionCurve::Ellipse(*e);
 }
 crate::inttools::intss::SurfaceCurve::Polyline(pts) => {
 return OffsetIntersectionCurve::Numerical(pts.clone());
 }
 _ => continue,
 }
 }

 OffsetIntersectionCurve::NoIntersection
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// UV Projection for Offset Curves
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Project a 3D point onto a surface and return the UV parameters.
///
/// This function handles trimmed surface boundaries and maintains accuracy
/// for the parameter space curve (PCurve).
///
/// # Algorithm
///
/// For analytical surfaces (Plane, Sphere, Cylinder, Cone, Torus), uses
/// direct analytical formulas with optional Newton-Raphson refinement.
/// For parametric surfaces (B-spline, Bezier), uses Newton-Raphson iteration.
///
/// # Precision
///
/// Target precision is TOLERANCE_LINEAR_ULTRA_STRICT for planes and TOLERANCE_LINEAR_RELAX_8 for curved surfaces.
pub fn project_point_to_surface_uv(
 point: DVec3,
 surf: &Surface3,
 hint_uv: Option<[f64; 2]>,
) -> Option<[f64; 2]> {
 match surf {
 Surface3::Plane(p) => project_point_to_plane_uv(point, p),
 Surface3::Sphere(s) => project_point_to_sphere_uv(point, s, hint_uv),
 Surface3::Cylinder(c) => project_point_to_cylinder_uv(point, c, hint_uv),
 Surface3::Cone(c) => project_point_to_cone_uv(point, c, hint_uv),
 Surface3::Torus(t) => project_point_to_torus_uv(point, t, hint_uv),
 Surface3::BSpline(_)
 | Surface3::Bezier(_)
 | Surface3::TriBezier(_)
 | Surface3::Offset(_)
 | Surface3::Trimmed(_) => {
 // For parametric surfaces, use Newton iteration
 project_point_to_parametric_surface(point, surf, hint_uv)
 }
 _ => {
 // Fallback: use surface evaluation at hint or center
 hint_uv.or(Some([0.5, 0.5]))
 }
 }
}

/// Compute a deterministic orthonormal frame from a normal vector.
///
/// This creates a consistent UV basis from a normal vector by:
/// 1. Finding the smallest component of the normal
/// 2. Cross-producting with the corresponding axis to get u_dir
/// 3. Cross-producting normal u_dir to get v_dir
///
/// This is more numerically stable than `any_perpendicular` and gives
/// consistent results for the same normal vector.
fn orthonormal_basis_from_normal(normal: DVec3) -> (DVec3, DVec3) {
 let n = normal.normalize_or(DVec3::Z);
 // Find the axis most perpendicular to n for numerical stability
 let abs_n = n.abs();
 let axis = if abs_n.x <= abs_n.y && abs_n.x <= abs_n.z {
 DVec3::X
 } else if abs_n.y <= abs_n.x && abs_n.y <= abs_n.z {
 DVec3::Y
 } else {
 DVec3::Z
 };
 let u_dir = n.cross(axis).normalize_or(DVec3::X);
 let v_dir = n.cross(u_dir).normalize_or(DVec3::Y);
 (u_dir, v_dir)
}
include!("e1.rs");
include!("e2.rs");
include!("e3.rs");

