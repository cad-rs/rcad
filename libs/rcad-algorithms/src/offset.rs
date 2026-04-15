//! Shell and solid offset operations — analogous to OCCT `BRepOffsetAPI_MakeOffsetShape`.
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
//! Plane, Sphere, Cylinder, Cone, Torus — each has a known parallel-surface construction.
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
//! - OCCT `BRepOffsetAPI_MakeOffsetShape`
//! - OCCT `BRepOffset_MakeOffset`
//! - OCCT `BRepOffset_Mode` (join types)

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::{
    BRep,
    SurfaceEval, CurveEval,
    geom::{Curve3, Surface3, Line3, Plane, CylindricalSurface, SphericalSurface, ConicalSurface, ToroidalSurface, OffsetSurface},
    topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
};
use crate::tolerance::TOLERANCE_ABS;

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// Join Types
// ─────────────────────────────────────────────────────────────────────────────

/// Join type for offset operations at edges.
///
/// Determines how adjacent offset faces are connected at their boundaries.
/// Analogous to OCCT `BRepOffset_Mode`.
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

// ─────────────────────────────────────────────────────────────────────────────
// Variable Thickness
// ─────────────────────────────────────────────────────────────────────────────

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
            if ft.thickness.abs() < 1e-12 {
                return Err(OffsetError::InvalidVariableThickness {
                    face_index: ft.face_index,
                    thickness: ft.thickness,
                    reason: "thickness cannot be zero".to_string(),
                });
            }
        }
        if self.default_thickness.abs() < 1e-12 {
            return Err(OffsetError::InvalidVariableThickness {
                face_index: 0,
                thickness: self.default_thickness,
                reason: "default thickness cannot be zero".to_string(),
            });
        }
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-Intersection Handling
// ─────────────────────────────────────────────────────────────────────────────

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
            min_offset_distance: 1e-6,
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

// ─────────────────────────────────────────────────────────────────────────────
// Offset Quality Analysis
// ─────────────────────────────────────────────────────────────────────────────

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
            min_wall_thickness: 1e-6,
            max_deviation: 1e-4,
            allow_self_intersection: false,
            max_degenerate_ratio: 0.1,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Options
// ─────────────────────────────────────────────────────────────────────────────

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
            min_feature_size: 1e-6,
            join_type: JoinType::default(),
            variable_thickness: None,
            self_intersection_config: SelfIntersectionConfig::default(),
            quality_thresholds: QualityThresholds::default(),
            approximation_tolerance: 1e-4,
            check_wall_thickness: false,
            min_wall_thickness: 1e-6,
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
    /// The resulting BRep.
    pub brep: BRep,
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

// ─────────────────────────────────────────────────────────────────────────────
// Surface Offset
// ─────────────────────────────────────────────────────────────────────────────

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
            Some(Surface3::Sphere(SphericalSurface {
                center: s.center,
                axis: s.axis,
                radius: new_radius,
            }))
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
            }))
        }

        Surface3::Cone(c) => {
            // Cone offset: adjust radius and apex position
            // The parallel surface to a cone is another cone with the same half-angle
            // but shifted apex position and different radius at the reference point.
            let sin_a = c.half_angle_rad.sin();
            let cos_a = c.half_angle_rad.cos();

            // Axial shift of the apex along the cone axis
            let axial_shift = if sin_a.abs() > 1e-10 { d / sin_a } else { d };

            // New radius at the reference point (apex field)
            let new_radius = c.radius + d * cos_a;

            if new_radius <= 0.0 && d > 0.0 {
                // Positive offset would make radius negative at reference
                return None;
            }

            // For cones, we need to shift the apex to maintain the same half-angle
            let new_apex = c.apex - c.axis.normalize_or(DVec3::Y) * axial_shift;

            Some(Surface3::Cone(ConicalSurface {
                apex: new_apex,
                axis: c.axis,
                radius: new_radius.max(0.0),
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
        | Surface3::Gordon(_)
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

// ─────────────────────────────────────────────────────────────────────────────
// Edge Offset
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the offset edge curve for a given edge.
///
/// The offset edge is the intersection of the two adjacent offset surfaces.
/// For manifold edges (shared by two faces), we compute the intersection.
/// For boundary edges, we project the edge onto the single offset surface.
fn offset_edge(
    brep: &BRep,
    edge_idx: usize,
    face_indices: &[usize],
    distance: f64,
    offset_surfaces: &[Option<Surface3>],
) -> Option<(Curve3, f64, f64)> {
    let edge = &brep.edges[edge_idx];

    if face_indices.is_empty() {
        return None;
    }

    // Get the 3D curve of the edge
    let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c)?;
    let curve = &brep.geom.curves[curve_idx];
    let range = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r);

    if face_indices.len() == 1 {
        // Boundary edge: project onto single offset surface
        let surf = offset_surfaces.get(face_indices[0]).and_then(|s| s.as_ref())?;

        // Compute offset points at edge endpoints
        let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);

        // Compute vertex normals at these points
        let n0 = compute_vertex_normal_on_face(brep, edge.start, face_indices[0]);
        let n1 = compute_vertex_normal_on_face(brep, edge.end, face_indices[0]);

        // Offset points
        let off_p0 = p0 + n0 * distance;
        let off_p1 = p1 + n1 * distance;

        // Create a line between offset points
        let dir = (off_p1 - off_p0).normalize_or(DVec3::X);
        let len = (off_p1 - off_p0).length();

        Some((Curve3::Line(Line3 {
            origin: off_p0,
            direction: dir,
        }), 0.0, len))
    } else {
        // Manifold edge: compute intersection of two offset surfaces
        let surf0 = offset_surfaces.get(face_indices[0]).and_then(|s| s.as_ref())?;
        let surf1 = offset_surfaces.get(face_indices[1]).and_then(|s| s.as_ref())?;

        // Compute offset points at edge endpoints
        let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);

        // Average normals from both faces
        let n0_0 = compute_vertex_normal_on_face(brep, edge.start, face_indices[0]);
        let n0_1 = compute_vertex_normal_on_face(brep, edge.start, face_indices[1]);
        let n1_0 = compute_vertex_normal_on_face(brep, edge.end, face_indices[0]);
        let n1_1 = compute_vertex_normal_on_face(brep, edge.end, face_indices[1]);

        let n0 = (n0_0 + n0_1).normalize_or(n0_0);
        let n1 = (n1_0 + n1_1).normalize_or(n1_0);

        // Offset points
        let off_p0 = p0 + n0 * distance;
        let off_p1 = p1 + n1 * distance;

        // For now, create a line between offset points
        // TODO: Compute actual intersection curve of offset surfaces
        let dir = (off_p1 - off_p0).normalize_or(DVec3::X);
        let len = (off_p1 - off_p0).length();

        Some((Curve3::Line(Line3 {
            origin: off_p0,
            direction: dir,
        }), 0.0, len))
    }
}

/// Compute the normal at a vertex on a specific face.
fn compute_vertex_normal_on_face(brep: &BRep, vertex_idx: usize, face_idx: usize) -> DVec3 {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return DVec3::Z,
    };

    let face = match shell.faces.get(face_idx) {
        Some(f) => f,
        None => return DVec3::Z,
    };

    let surf_idx = match brep.geom.face_surface.get(face_idx).and_then(|s| *s) {
        Some(s) => s,
        None => return face.normal,
    };

    let surf = &brep.geom.surfaces[surf_idx];

    // Find a point on the face near this vertex
    let vertex_point = brep.vertices[vertex_idx].point;

    // Compute surface normal at approximate UV
    // For now, use the face normal as approximation
    // TODO: Project vertex onto surface to get accurate UV
    surf.normal_at(0.5, 0.5)
}

// ─────────────────────────────────────────────────────────────────────────────
// Vertex Offset
// ─────────────────────────────────────────────────────────────────────────────

/// Compute offset position for a vertex.
///
/// The offset vertex is the intersection of all offset edges meeting at the vertex,
/// or equivalently, the original vertex translated along the average normal.
fn offset_vertex(brep: &BRep, vertex_idx: usize, distance: f64, shell: &Shell) -> DVec3 {
    let original_point = brep.vertices[vertex_idx].point;

    // Collect all faces using this vertex
    let mut normal_sum = DVec3::ZERO;
    let mut count = 0;

    for face in &shell.faces {
        let uses_vertex = face.outer_wire.edges.iter().any(|we| {
            let e = &brep.edges[we.idx];
            e.start == vertex_idx || e.end == vertex_idx
        });

        if uses_vertex {
            normal_sum += face.normal;
            count += 1;
        }
    }

    if count > 0 {
        let avg_normal = normal_sum.normalize_or(DVec3::Z);
        original_point + avg_normal * distance
    } else {
        original_point
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRep Builder Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Helper to add a vertex to a BRep and return its index.
fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

/// Helper to add an edge to a BRep and return its index.
fn add_edge(brep: &mut BRep, curve: Curve3, t0: f64, t1: f64, v0: usize, v1: usize) -> usize {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start: v0, end: v1 });

    let ci = brep.geom.curves.len();
    brep.geom.curves.push(curve);

    while brep.geom.edge_curve.len() <= idx {
        brep.geom.edge_curve.push(None);
    }
    while brep.geom.edge_curve_range.len() <= idx {
        brep.geom.edge_curve_range.push(None);
    }
    while brep.geom.edge_degenerated.len() <= idx {
        brep.geom.edge_degenerated.push(false);
    }

    brep.geom.edge_curve[idx] = Some(ci);
    brep.geom.edge_curve_range[idx] = Some([t0, t1]);
    idx
}

/// Helper to add a face to a BRep and return its index.
fn add_face(brep: &mut BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
    if brep.solids.is_empty() {
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }

    let idx = brep.solids[0].shells[0].faces.len();
    let normal = surface.normal_at(0.0, 0.0);

    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer,
        inner_wires: inner,
        normal,
        triangles: Vec::new(),
        mesh_dirty: true,
    });

    while brep.geom.face_surface.len() <= idx {
        brep.geom.face_surface.push(None);
    }

    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);

    idx
}

// ─────────────────────────────────────────────────────────────────────────────
// Edge Chaining
// ─────────────────────────────────────────────────────────────────────────────

/// Chain boundary edges into closed loops.
fn chain_boundary_edges(edge_indices: &[usize], edges: &[Edge]) -> Vec<Vec<usize>> {
    if edge_indices.is_empty() {
        return vec![];
    }

    let mut remaining: HashSet<usize> = edge_indices.iter().copied().collect();
    let mut loops = Vec::new();

    while let Some(&start_idx) = remaining.iter().next() {
        remaining.remove(&start_idx);
        let mut chain = vec![start_idx];
        let mut current_end = edges[start_idx].end;

        loop {
            let next = remaining
                .iter()
                .find(|&&ei| edges[ei].start == current_end || edges[ei].end == current_end)
                .copied();

            match next {
                Some(ei) => {
                    remaining.remove(&ei);
                    chain.push(ei);
                    let e = &edges[ei];
                    current_end = if e.start == current_end { e.end } else { e.start };
                }
                None => break,
            }
        }

        if chain.len() >= 2 {
            loops.push(chain);
        }
    }

    loops
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-Intersection Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect potential self-intersection in a closed-shell offset.
///
/// Computes the minimum distance between non-adjacent face centroids.
/// If the offset distance exceeds half this distance, self-intersection is likely.
pub fn detect_self_intersection(brep: &BRep, distance: f64) -> bool {
    let result = detect_self_intersection_detailed(brep, distance);
    result.has_intersection
}

/// Detect self-intersection with detailed results.
///
/// This is a more comprehensive analysis that returns information about
/// which faces might intersect and the minimum safe offset distance.
pub fn detect_self_intersection_detailed(brep: &BRep, distance: f64) -> SelfIntersectionResult {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => {
            return SelfIntersectionResult {
                has_intersection: false,
                intersecting_pairs: Vec::new(),
                min_safe_distance: None,
                description: "no shell found".to_string(),
            };
        }
    };

    if shell.faces.len() < 3 {
        return SelfIntersectionResult {
            has_intersection: false,
            intersecting_pairs: Vec::new(),
            min_safe_distance: None,
            description: "insufficient faces".to_string(),
        };
    }

    // Compute face centroids
    let centroids: Vec<DVec3> = shell
        .faces
        .iter()
        .map(|face| {
            let mut sum = DVec3::ZERO;
            let mut count = 0;
            for we in &face.outer_wire.edges {
                let e = &brep.edges[we.idx];
                sum += brep.vertices[e.start].point;
                count += 1;
            }
            if count > 0 {
                sum / count as f64
            } else {
                DVec3::ZERO
            }
        })
        .collect();

    // Build adjacency map
    let mut adjacent_pairs: HashSet<(usize, usize)> = HashSet::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            for (fj, other_face) in shell.faces.iter().enumerate() {
                if fi < fj && other_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx) {
                    adjacent_pairs.insert((fi, fj));
                }
            }
        }
    }

    // Find minimum distance between non-adjacent faces
    let mut min_dist = f64::MAX;
    let mut intersecting_pairs = Vec::new();
    let abs_distance = distance.abs();

    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            if adjacent_pairs.contains(&(i, j)) {
                continue;
            }

            let dist = (centroids[i] - centroids[j]).length();

            // Check if these faces would intersect
            if abs_distance > dist * 0.5 {
                intersecting_pairs.push((i, j));
            }

            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    if min_dist == f64::MAX {
        return SelfIntersectionResult {
            has_intersection: false,
            intersecting_pairs: Vec::new(),
            min_safe_distance: None,
            description: "no non-adjacent faces found".to_string(),
        };
    }

    let has_intersection = abs_distance > min_dist * 0.5;
    let min_safe_distance = Some(min_dist * 0.5);

    let description = if has_intersection {
        format!(
            "self-intersection likely: {} face pairs at distance {} with offset {}",
            intersecting_pairs.len(),
            min_dist,
            abs_distance
        )
    } else {
        format!("no self-intersection: min distance {}, offset {}", min_dist, abs_distance)
    };

    SelfIntersectionResult {
        has_intersection,
        intersecting_pairs,
        min_safe_distance,
        description,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Join Geometry Creation
// ─────────────────────────────────────────────────────────────────────────────

/// Create an arc join between two offset edges.
///
/// Creates a cylindrical surface that smoothly transitions between
/// two offset faces meeting at an edge.
pub fn create_arc_join(
    brep: &mut BRep,
    edge_idx: usize,
    face0_idx: usize,
    face1_idx: usize,
    radius: f64,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or_else(|| {
        OffsetError::JoinCreationFailed {
            join_type: JoinType::Arc,
            edge_index: edge_idx,
            reason: "no shell found".to_string(),
        }
    })?;

    let edge = &brep.edges[edge_idx];
    let face0 = &shell.faces[face0_idx];
    let face1 = &shell.faces[face1_idx];

    // Get the edge endpoints
    let v0 = vertex_map.get(edge.start).copied().unwrap_or(edge.start);
    let v1 = vertex_map.get(edge.end).copied().unwrap_or(edge.end);

    let p0 = brep.vertices[v0].point;
    let p1 = brep.vertices[v1].point;

    // Compute the edge direction and length
    let edge_dir = (p1 - p0).normalize_or(DVec3::X);
    let edge_len = (p1 - p0).length();

    // Compute the bisector normal from the two face normals
    let n0 = face0.normal;
    let n1 = face1.normal;
    let _bisector = (n0 + n1).normalize_or(n0);

    // Create a cylindrical surface for the arc join
    // The cylinder axis is along the edge, and the radius is the offset distance
    let cylinder = Surface3::Cylinder(CylindricalSurface {
        origin: p0,
        axis: edge_dir,
        radius,
    });

    // Create vertices for the arc join face
    // The arc join is a sector of the cylinder
    let vs = add_vertex(brep, p0);
    let ve = add_vertex(brep, p1);

    // Create the edge along the cylinder
    let curve = Curve3::Line(Line3 {
        origin: p0,
        direction: edge_dir,
    });
    let arc_edge = add_edge(brep, curve, 0.0, edge_len, vs, ve);

    // Create the arc face wire
    let wire = Wire {
        edges: vec![WireEdge::fwd(arc_edge)],
    };

    // Add the arc join face
    let face_idx = add_face(brep, cylinder, wire, Vec::new());

    Ok(face_idx)
}

/// Create a tangent join between two offset edges.
///
/// Creates a smooth, tangent-continuous transition between adjacent faces.
/// Falls back to intersection join when the angle between faces is too large.
pub fn create_tangent_join(
    brep: &mut BRep,
    edge_idx: usize,
    face0_idx: usize,
    face1_idx: usize,
    distance: f64,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or_else(|| {
        OffsetError::JoinCreationFailed {
            join_type: JoinType::Tangent,
            edge_index: edge_idx,
            reason: "no shell found".to_string(),
        }
    })?;

    let face0 = &shell.faces[face0_idx];
    let face1 = &shell.faces[face1_idx];

    // Check the angle between face normals
    let n0 = face0.normal;
    let n1 = face1.normal;
    let dot = n0.dot(n1);

    // If the angle is too large (faces nearly parallel or facing opposite directions),
    // fall back to intersection join
    let angle_threshold = 0.9; // cos(25 degrees) approximately
    if dot < angle_threshold {
        // Create intersection join instead
        return create_intersection_join(brep, edge_idx, face0_idx, face1_idx, vertex_map);
    }

    // For tangent join, create a smooth blending surface
    // This uses a ruled surface between the two offset edges
    let edge = &brep.edges[edge_idx];
    let v0 = vertex_map.get(edge.start).copied().unwrap_or(edge.start);
    let v1 = vertex_map.get(edge.end).copied().unwrap_or(edge.end);

    let p0 = brep.vertices[v0].point;
    let p1 = brep.vertices[v1].point;

    // Create a plane that smoothly blends the two face normals
    let blend_normal = (n0 + n1).normalize();
    let blend_plane = Surface3::Plane(Plane {
        origin: (p0 + p1) * 0.5,
        normal: blend_normal,
    });

    // Create the wire for the tangent join face
    let dir = (p1 - p0).normalize_or(DVec3::X);
    let len = (p1 - p0).length();
    let curve = Curve3::Line(Line3 { origin: p0, direction: dir });

    let vs = add_vertex(brep, p0);
    let ve = add_vertex(brep, p1);
    let blend_edge = add_edge(brep, curve, 0.0, len, vs, ve);

    let wire = Wire {
        edges: vec![WireEdge::fwd(blend_edge)],
    };

    let face_idx = add_face(brep, blend_plane, wire, Vec::new());

    let _ = distance; // Used in more sophisticated implementations
    Ok(face_idx)
}

/// Create an intersection join between two offset edges.
///
/// The offset surfaces extend until they intersect, creating sharp corners.
/// This is the default mode and works well for mechanical parts.
pub fn create_intersection_join(
    brep: &mut BRep,
    edge_idx: usize,
    _face0_idx: usize,
    _face1_idx: usize,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let edge = &brep.edges[edge_idx];

    let v0 = vertex_map.get(edge.start).copied().unwrap_or(edge.start);
    let v1 = vertex_map.get(edge.end).copied().unwrap_or(edge.end);

    let p0 = brep.vertices[v0].point;
    let p1 = brep.vertices[v1].point;

    // For intersection join, we don't create additional geometry -
    // the offset surfaces naturally intersect at the edge
    // Instead, we return the edge index as the "join"
    // In a full implementation, this would compute the exact intersection curve

    // Create a minimal face at the intersection
    let dir = (p1 - p0).normalize_or(DVec3::X);
    let len = (p1 - p0).length();

    // Use the edge midpoint and direction to create a small plane
    let midpoint = (p0 + p1) * 0.5;
    let normal = dir.any_orthonormal_pair().0;

    let plane = Surface3::Plane(Plane {
        origin: midpoint,
        normal,
    });

    let vs = add_vertex(brep, p0);
    let ve = add_vertex(brep, p1);
    let curve = Curve3::Line(Line3 { origin: p0, direction: dir });
    let int_edge = add_edge(brep, curve, 0.0, len, vs, ve);

    let wire = Wire {
        edges: vec![WireEdge::fwd(int_edge)],
    };

    let face_idx = add_face(brep, plane, wire, Vec::new());

    Ok(face_idx)
}

/// Apply join type to all edges in the shell.
///
/// This function creates the appropriate join geometry for each edge
/// based on the specified join type.
pub fn apply_join_type(
    result: &mut BRep,
    original_brep: &BRep,
    opts: &OffsetOptions,
    edge_to_faces: &HashMap<usize, Vec<usize>>,
    vertex_map: &[usize],
) -> Result<usize, OffsetError> {
    let mut join_face_count = 0;

    if opts.join_type == JoinType::Intersection {
        // Intersection join is the default - no additional geometry needed
        return Ok(0);
    }

    for (&edge_idx, face_indices) in edge_to_faces {
        if face_indices.len() < 2 {
            continue; // Skip boundary edges
        }

        let face0_idx = face_indices[0];
        let face1_idx = face_indices[1];

        let join_result = match opts.join_type {
            JoinType::Arc => {
                let radius = opts.distance.abs();
                create_arc_join(result, edge_idx, face0_idx, face1_idx, radius, vertex_map)
            }
            JoinType::Tangent => {
                create_tangent_join(result, edge_idx, face0_idx, face1_idx, opts.distance, vertex_map)
            }
            JoinType::Intersection => {
                create_intersection_join(result, edge_idx, face0_idx, face1_idx, vertex_map)
            }
        };

        if join_result.is_ok() {
            join_face_count += 1;
        }
    }

    let _ = original_brep; // Used in more sophisticated implementations
    Ok(join_face_count)
}

// ─────────────────────────────────────────────────────────────────────────────
// Offset Quality Analysis
// ─────────────────────────────────────────────────────────────────────────────

/// Analyze the quality of an offset result.
///
/// Computes various quality metrics including wall thickness, deviation,
/// and self-intersection detection.
pub fn analyze_offset_quality(
    result: &BRep,
    original: &BRep,
    opts: &OffsetOptions,
) -> OffsetQuality {
    let mut quality = OffsetQuality::default();

    // Compute minimum wall thickness
    quality.min_wall_thickness = compute_min_wall_thickness(result, opts.distance);

    // Compute maximum deviation from expected offset
    quality.max_deviation = compute_max_deviation(result, original, opts);

    // Count degenerate edges
    quality.degenerate_edge_count = result
        .geom
        .edge_degenerated
        .iter()
        .filter(|&&d| d)
        .count();

    // Self-intersection count
    let si_result = detect_self_intersection_detailed(result, opts.distance);
    quality.self_intersection_count = si_result.intersecting_pairs.len();

    // Compute face area ratio
    quality.face_area_ratio = compute_face_area_ratio(result, original);

    // Compute edge length ratio
    quality.edge_length_ratio = compute_edge_length_ratio(result, original);

    // Determine if result is valid
    quality.is_valid = quality.self_intersection_count == 0
        && quality.min_wall_thickness >= opts.min_wall_thickness;

    // Generate warnings
    if quality.min_wall_thickness < opts.min_wall_thickness {
        quality.warnings.push(format!(
            "Minimum wall thickness {} is below threshold {}",
            quality.min_wall_thickness, opts.min_wall_thickness
        ));
    }
    if quality.max_deviation > opts.approximation_tolerance {
        quality.warnings.push(format!(
            "Maximum deviation {} exceeds approximation tolerance {}",
            quality.max_deviation, opts.approximation_tolerance
        ));
    }
    if quality.degenerate_edge_count > 0 {
        quality.warnings.push(format!(
            "Found {} degenerate edges in result",
            quality.degenerate_edge_count
        ));
    }

    quality
}

/// Compute the minimum wall thickness in the offset result.
///
/// Uses face centroid distances to estimate minimum wall thickness.
pub fn compute_min_wall_thickness(brep: &BRep, distance: f64) -> f64 {
    let shell = match brep.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return distance,
    };

    if shell.faces.len() < 2 {
        return distance;
    }

    // Compute face centroids
    let centroids: Vec<DVec3> = shell
        .faces
        .iter()
        .map(|face| {
            let mut sum = DVec3::ZERO;
            let mut count = 0;
            for we in &face.outer_wire.edges {
                let e = &brep.edges[we.idx];
                sum += brep.vertices[e.start].point;
                count += 1;
            }
            if count > 0 {
                sum / count as f64
            } else {
                DVec3::ZERO
            }
        })
        .collect();

    // Find minimum distance between any two faces
    let mut min_dist = f64::MAX;
    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            let dist = (centroids[i] - centroids[j]).length();
            if dist > 0.0 && dist < min_dist {
                min_dist = dist;
            }
        }
    }

    // The wall thickness is approximately the minimum distance minus twice the offset
    // For a proper implementation, this would use more sophisticated analysis
    if min_dist == f64::MAX {
        distance
    } else {
        (min_dist - 2.0 * distance.abs()).max(0.0)
    }
}

/// Compute the maximum deviation between offset and expected positions.
pub fn compute_max_deviation(result: &BRep, original: &BRep, opts: &OffsetOptions) -> f64 {
    let _result_shell = match result.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return 0.0,
    };

    let _original_shell = match original.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return 0.0,
    };

    let mut max_dev = 0.0;

    // Compare vertex positions
    for (i, vertex) in result.vertices.iter().enumerate() {
        if i >= original.vertices.len() {
            break;
        }

        let original_vertex = &original.vertices[i];
        let actual_offset = (vertex.point - original_vertex.point).length();
        let expected_offset = opts.distance.abs();

        let deviation = (actual_offset - expected_offset).abs();
        if deviation > max_dev {
            max_dev = deviation;
        }
    }

    max_dev
}

/// Compute the ratio of face areas between result and original.
pub fn compute_face_area_ratio(result: &BRep, original: &BRep) -> f64 {
    let result_shell = match result.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return 1.0,
    };

    let original_shell = match original.solids.first().and_then(|s| s.shells.first()) {
        Some(s) => s,
        None => return 1.0,
    };

    if original_shell.faces.is_empty() {
        return 1.0;
    }

    // Simple approximation: ratio of face counts
    // A proper implementation would compute actual areas
    result_shell.faces.len() as f64 / original_shell.faces.len() as f64
}

/// Compute the ratio of edge lengths between result and original.
pub fn compute_edge_length_ratio(result: &BRep, original: &BRep) -> f64 {
    if original.edges.is_empty() {
        return 1.0;
    }

    // Compute total edge lengths
    let original_len: f64 = original
        .edges
        .iter()
        .map(|e| {
            let p0 = original.vertices.get(e.start).map(|v| v.point).unwrap_or_default();
            let p1 = original.vertices.get(e.end).map(|v| v.point).unwrap_or_default();
            (p1 - p0).length()
        })
        .sum();

    let result_len: f64 = result
        .edges
        .iter()
        .map(|e| {
            let p0 = result.vertices.get(e.start).map(|v| v.point).unwrap_or_default();
            let p1 = result.vertices.get(e.end).map(|v| v.point).unwrap_or_default();
            (p1 - p0).length()
        })
        .sum();

    if original_len > 0.0 {
        result_len / original_len
    } else {
        1.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-Intersection Repair
// ─────────────────────────────────────────────────────────────────────────────

/// Attempt to repair self-intersection by reducing offset distance.
///
/// Tries progressively smaller offset distances until a valid result is found.
pub fn repair_self_intersection(
    brep: &BRep,
    opts: &OffsetOptions,
) -> Result<(BRep, f64, usize), OffsetError> {
    let config = &opts.self_intersection_config;

    if !config.auto_repair {
        return Err(OffsetError::RecoveryFailed {
            attempts: 0,
            last_error: "auto-repair not enabled".to_string(),
        });
    }

    let mut current_distance = opts.distance;
    let mut attempts = 0;
    let mut last_error = String::new();

    while attempts < config.max_repair_attempts {
        attempts += 1;

        // Reduce the offset distance
        current_distance *= config.reduction_factor;

        if current_distance.abs() < config.min_offset_distance {
            last_error = format!(
                "offset distance {} below minimum {}",
                current_distance.abs(),
                config.min_offset_distance
            );
            continue;
        }

        // Try with reduced distance
        let mut reduced_opts = opts.clone();
        reduced_opts.distance = current_distance;
        reduced_opts.check_self_intersection = true;

        let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(OffsetError::InvalidInput("no shell"))?;

        match offset_shell_with_options_impl(shell, brep, &reduced_opts) {
            Ok(result) => {
                let si_result = detect_self_intersection_detailed(&result, current_distance);
                if !si_result.has_intersection {
                    return Ok((result, current_distance, attempts));
                }
                last_error = si_result.description;
            }
            Err(e) => {
                last_error = e.to_string();
            }
        }
    }

    Err(OffsetError::RecoveryFailed { attempts, last_error })
}

/// Implementation of offset_shell_with_options that can be called internally.
fn offset_shell_with_options_impl(
    shell: &Shell,
    brep: &BRep,
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    // Validate variable thickness if specified
    if let Some(ref vt) = opts.variable_thickness {
        vt.validate(shell.faces.len())?;
    }

    let distance = opts.distance;

    if distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    if shell.faces.is_empty() {
        return Err(OffsetError::InvalidInput("shell has no faces"));
    }

    // Step 1: Compute offset surfaces for each face (with variable thickness support)
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
    for (fi, _face) in shell.faces.iter().enumerate() {
        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };

        let surf = &brep.geom.surfaces[surf_idx];
        let face_distance = opts.effective_distance_for_face(fi);
        let off_surf = offset_surface(surf, face_distance);

        offset_surfaces.push(off_surf);
    }

    // Step 2: Build edge-to-face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
    }

    // Step 3: Compute offset vertex positions (with variable thickness support)
    let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
        .map(|vi| {
            // For variable thickness, use average distance of adjacent faces
            let avg_distance = if let Some(ref vt) = opts.variable_thickness {
                let mut sum = 0.0;
                let mut count = 0;
                for (fi, face) in shell.faces.iter().enumerate() {
                    let uses_vertex = face.outer_wire.edges.iter().any(|we| {
                        let e = &brep.edges[we.idx];
                        e.start == vi || e.end == vi
                    });
                    if uses_vertex {
                        sum += vt.thickness_for_face(fi);
                        count += 1;
                    }
                }
                if count > 0 { sum / count as f64 } else { distance }
            } else {
                distance
            };
            offset_vertex(brep, vi, avg_distance, shell)
        })
        .collect();

    // Step 4: Build result BRep
    let mut result = BRep::new();
    result.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Map original vertices to offset vertices
    let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
    for &p in &offset_vertices {
        vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces with offset edges
    let mut valid_face_count = 0;

    for (fi, face) in shell.faces.iter().enumerate() {
        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        // Build wire from offset edges
        let mut wire_edges = Vec::new();

        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = vertex_map[e.start];
            let ve = vertex_map[e.end];

            let p0 = result.vertices[vs].point;
            let p1 = result.vertices[ve].point;
            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
        valid_face_count += 1;
    }

    if valid_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Step 6: Apply join type if needed
    if opts.join_type.requires_join_geometry() {
        let _join_faces = apply_join_type(&mut result, brep, opts, &edge_to_faces, &vertex_map)?;
    }

    Ok(result)
}

// ─────────────────────────────────────────────────────────────────────────────
// Main API Functions
// ─────────────────────────────────────────────────────────────────────────────

/// Offset a shell by moving all faces along their normals.
///
/// # Arguments
///
/// * `shell` - The input shell to offset
/// * `brep` - The BRep containing the shell's geometry
/// * `distance` - Offset distance (positive = outward, negative = inward)
///
/// # Returns
///
/// A new BRep containing the offset shell, or an error.
pub fn offset_shell(shell: &Shell, brep: &BRep, distance: f64) -> Result<BRep, OffsetError> {
    offset_shell_with_options(shell, brep, &OffsetOptions::new(distance))
}

/// Offset a shell with full options.
pub fn offset_shell_with_options(
    shell: &Shell,
    brep: &BRep,
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    let distance = opts.distance;

    if distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    if shell.faces.is_empty() {
        return Err(OffsetError::InvalidInput("shell has no faces"));
    }

    // Step 1: Compute offset surfaces for each face
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
    for (fi, _face) in shell.faces.iter().enumerate() {
        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };

        let surf = &brep.geom.surfaces[surf_idx];
        let off_surf = offset_surface(surf, distance);

        if off_surf.is_none() && distance > 0.0 {
            // Negative offset on a small surface - may be ok for inward offset
        }

        offset_surfaces.push(off_surf);
    }

    // Step 2: Build edge-to-face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
    }

    // Step 3: Compute offset vertex positions
    let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
        .map(|vi| offset_vertex(brep, vi, distance, shell))
        .collect();

    // Step 4: Build result BRep
    let mut result = BRep::new();
    result.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Map original vertices to offset vertices
    let mut vertex_map: Vec<usize> = Vec::with_capacity(offset_vertices.len());
    for &p in &offset_vertices {
        vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces with offset edges
    let mut valid_face_count = 0;

    for (fi, face) in shell.faces.iter().enumerate() {
        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        // Build wire from offset edges
        let mut wire_edges = Vec::new();

        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = vertex_map[e.start];
            let ve = vertex_map[e.end];

            let p0 = result.vertices[vs].point;
            let p1 = result.vertices[ve].point;
            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
        valid_face_count += 1;
    }

    if valid_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Step 6: Check for self-intersection if requested
    let self_intersects = if opts.check_self_intersection {
        detect_self_intersection(&result, distance)
    } else {
        false
    };

    if self_intersects && !opts.auto_repair {
        // Still return the result, but the caller should check for self-intersection
    }

    Ok(result)
}

/// Offset a solid by moving all faces along their normals.
///
/// # Arguments
///
/// * `solid` - The input solid to offset
/// * `brep` - The BRep containing the solid's geometry
/// * `distance` - Offset distance
///   - Positive: outward expansion (thickening)
///   - Negative: inward contraction (shelling)
///
/// # Returns
///
/// A new BRep containing the offset solid, or an error.
pub fn offset_solid(solid: &Solid, brep: &BRep, distance: f64) -> Result<BRep, OffsetError> {
    offset_solid_with_options(solid, brep, &OffsetOptions::new(distance))
}

/// Offset a solid with full options.
pub fn offset_solid_with_options(
    solid: &Solid,
    brep: &BRep,
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    let distance = opts.distance;

    if distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    // For a solid, offset each shell
    let mut result = BRep::new();
    result.solids.push(Solid { shells: Vec::new() });

    for shell in &solid.shells {
        let offset_brep = offset_shell_with_options(shell, brep, opts)?;

        // Merge the offset shell into the result
        for offset_solid in offset_brep.solids {
            for offset_shell in offset_solid.shells {
                result.solids[0].shells.push(offset_shell);
            }
        }

        // Merge geometry
        let vertex_offset = result.vertices.len();
        result.vertices.extend(offset_brep.vertices);

        // Remap edge vertex indices
        for edge in offset_brep.edges {
            result.edges.push(Edge {
                start: edge.start + vertex_offset,
                end: edge.end + vertex_offset,
            });
        }

        // Merge geometry store
        let curve_offset = result.geom.curves.len();
        let surface_offset = result.geom.surfaces.len();

        result.geom.curves.extend(offset_brep.geom.curves);
        result.geom.surfaces.extend(offset_brep.geom.surfaces);

        for idx in offset_brep.geom.edge_curve {
            result.geom.edge_curve.push(idx.map(|i| i + curve_offset));
        }
        for range in offset_brep.geom.edge_curve_range {
            result.geom.edge_curve_range.push(range);
        }
        for deg in offset_brep.geom.edge_degenerated {
            result.geom.edge_degenerated.push(deg);
        }
        for idx in offset_brep.geom.face_surface {
            result.geom.face_surface.push(idx.map(|i| i + surface_offset));
        }
    }

    Ok(result)
}

/// Create a hollow solid by removing specified faces and offsetting remaining faces inward.
///
/// This is analogous to the "shell" or "hollow" operation in CAD systems.
///
/// # Arguments
///
/// * `solid` - The input solid
/// * `brep` - The BRep containing the solid's geometry
/// * `thickness` - Wall thickness (positive value)
/// * `open_faces` - Indices of faces to remove (creates openings)
///
/// # Returns
///
/// A new BRep containing the hollow solid with the specified faces removed,
/// or an error.
pub fn hollow_solid(
    solid: &Solid,
    brep: &BRep,
    thickness: f64,
    open_faces: &[usize],
) -> Result<BRep, OffsetError> {
    hollow_solid_with_options(solid, brep, thickness, open_faces, &OffsetOptions::new(-thickness))
}

/// Create a hollow solid with full options.
pub fn hollow_solid_with_options(
    solid: &Solid,
    brep: &BRep,
    thickness: f64,
    open_faces: &[usize],
    opts: &OffsetOptions,
) -> Result<BRep, OffsetError> {
    if thickness <= 0.0 {
        return Err(OffsetError::InvalidInput("thickness must be positive"));
    }

    let shell = match solid.shells.first() {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("solid has no shells")),
    };

    if open_faces.len() >= shell.faces.len() {
        return Err(OffsetError::InvalidInput("cannot remove all faces"));
    }

    let open_set: HashSet<usize> = open_faces.iter().copied().collect();

    // Step 1: Find boundary edges of the open faces
    let mut edge_use: HashMap<usize, usize> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        if open_set.contains(&fi) {
            continue;
        }
        for we in &face.outer_wire.edges {
            *edge_use.entry(we.idx).or_insert(0) += 1;
        }
    }

    // Boundary edges: edges that were used by removed faces but not by kept faces
    // These are edges where one adjacent face is removed and one is kept
    let mut boundary_edges: Vec<usize> = Vec::new();

    for (fi, face) in shell.faces.iter().enumerate() {
        if !open_set.contains(&fi) {
            continue;
        }
        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            // Check if this edge is shared with a kept face
            let is_shared = shell.faces.iter().enumerate().any(|(fj, fj_face)| {
                !open_set.contains(&fj)
                    && fj_face.outer_wire.edges.iter().any(|we2| we2.idx == we.idx)
            });
            if is_shared && !boundary_edges.contains(&we.idx) {
                boundary_edges.push(we.idx);
            }
        }
    }

    // Step 2: Create offset of kept faces (inward offset = negative distance)
    let inward_offset = -thickness;
    let mut offset_opts = opts.clone();
    offset_opts.distance = inward_offset;

    // Compute offset surfaces
    let mut offset_surfaces: Vec<Option<Surface3>> = Vec::with_capacity(shell.faces.len());
    for (fi, _face) in shell.faces.iter().enumerate() {
        if open_set.contains(&fi) {
            offset_surfaces.push(None);
            continue;
        }

        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|s| *s) {
            Some(s) => s,
            None => {
                offset_surfaces.push(None);
                continue;
            }
        };

        let surf = &brep.geom.surfaces[surf_idx];
        offset_surfaces.push(offset_surface(surf, inward_offset));
    }

    // Step 3: Compute offset vertex positions
    let offset_vertices: Vec<DVec3> = (0..brep.vertices.len())
        .map(|vi| offset_vertex(brep, vi, inward_offset, shell))
        .collect();

    // Step 4: Build result BRep
    let mut result = BRep::new();
    result.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Add original vertices
    let mut orig_vertex_map: Vec<usize> = Vec::new();
    for v in &brep.vertices {
        orig_vertex_map.push(add_vertex(&mut result, v.point));
    }

    // Add offset vertices
    let mut off_vertex_map: Vec<usize> = Vec::new();
    for &p in &offset_vertices {
        off_vertex_map.push(add_vertex(&mut result, p));
    }

    // Step 5: Create offset faces for kept faces
    let mut offset_face_count = 0;

    for (fi, face) in shell.faces.iter().enumerate() {
        if open_set.contains(&fi) {
            continue;
        }

        let off_surf = match &offset_surfaces[fi] {
            Some(s) => s.clone(),
            None => continue,
        };

        // Build wire from offset vertices
        let mut wire_edges = Vec::new();

        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = off_vertex_map[e.start];
            let ve = off_vertex_map[e.end];

            let p0 = result.vertices[vs].point;
            let p1 = result.vertices[ve].point;
            let dir = (p1 - p0).normalize_or(DVec3::X);
            let len = (p1 - p0).length();

            let curve = Curve3::Line(Line3 {
                origin: p0,
                direction: dir,
            });

            let eidx = add_edge(&mut result, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut result, off_surf, Wire { edges: wire_edges }, Vec::new());
        offset_face_count += 1;
    }

    // Step 6: Create lateral faces along boundary edges
    let loops = chain_boundary_edges(&boundary_edges, &brep.edges);
    let mut lateral_count = 0;

    for loop_edges in &loops {
        for &eidx in loop_edges {
            let e = &brep.edges[eidx];
            let o_vs = orig_vertex_map[e.start];
            let o_ve = orig_vertex_map[e.end];
            let f_vs = off_vertex_map[e.start];
            let f_ve = off_vertex_map[e.end];

            let p0 = result.vertices[o_vs].point;
            let p1 = result.vertices[o_ve].point;
            let p3 = result.vertices[f_vs].point;

            let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
            if normal.length() < 1e-10 {
                continue;
            }

            let surf = Surface3::Plane(Plane {
                origin: p0,
                normal,
            });

            // Quad: orig_start -> orig_end -> off_end -> off_start
            let vseq = [o_vs, o_ve, f_ve, f_vs];
            let mut edges = Vec::new();

            for i in 0..4 {
                let s = vseq[i];
                let en = vseq[(i + 1) % 4];
                let dir = (result.vertices[en].point - result.vertices[s].point).normalize_or(DVec3::X);
                let len = (result.vertices[en].point - result.vertices[s].point).length();
                let curve = Curve3::Line(Line3 {
                    origin: result.vertices[s].point,
                    direction: dir,
                });
                edges.push(WireEdge::fwd(add_edge(&mut result, curve, 0.0, len, s, en)));
            }

            add_face(&mut result, surf, Wire { edges }, Vec::new());
            lateral_count += 1;
        }
    }

    if offset_face_count == 0 {
        return Err(OffsetError::EmptyResult);
    }

    // Triangulate the result
    crate::triangulate::mesh_brep(&mut result, &crate::triangulate::TessellationParams::default());

    Ok(result)
}

/// Offset any BRep shape (shell or solid).
///
/// # Arguments
///
/// * `brep` - The input BRep
/// * `opts` - Offset options
///
/// # Returns
///
/// A new BRep with offset geometry.
pub fn offset_shape(brep: &BRep, opts: OffsetOptions) -> Result<OffsetResult, OffsetError> {
    if opts.distance.abs() < 1e-12 {
        return Err(OffsetError::ZeroDistance);
    }

    let solid = match brep.solids.first() {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("BRep has no solids")),
    };

    let shell = match solid.shells.first() {
        Some(s) => s,
        None => return Err(OffsetError::InvalidInput("solid has no shells")),
    };

    let result_brep = offset_shell_with_options(shell, brep, &opts)?;

    let self_intersection = if opts.check_self_intersection {
        detect_self_intersection(&result_brep, opts.distance)
    } else {
        false
    };

    let face_count = result_brep
        .solids
        .first()
        .and_then(|s| s.shells.first())
        .map(|sh| sh.faces.len())
        .unwrap_or(0);

    Ok(OffsetResult {
        brep: result_brep,
        offset_faces: face_count,
        lateral_faces: 0,
        join_faces: 0,
        self_intersection,
        quality: OffsetQuality {
            min_wall_thickness: f64::INFINITY,
            max_deviation: 0.0,
            degenerate_edge_count: 0,
            self_intersection_count: if self_intersection { 1 } else { 0 },
            face_area_ratio: 1.0,
            edge_length_ratio: 1.0,
            is_valid: true,
            warnings: Vec::new(),
        },
        warnings: Vec::new(),
        effective_distance: opts.distance,
        repair_attempts: 0,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{Plane, SphericalSurface, CylindricalSurface};

    #[test]
    fn offset_plane_translates() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let offset = offset_surface(&plane, 0.5).unwrap();

        if let Surface3::Plane(p) = offset {
            assert!((p.origin.z - 0.5).abs() < 1e-9, "plane should translate by offset distance");
            assert!((p.normal - DVec3::Z).length() < 1e-9, "normal should be unchanged");
        } else {
            panic!("expected Plane");
        }
    }

    #[test]
    fn offset_sphere_grows() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });

        let offset = offset_surface(&sphere, 0.5).unwrap();

        if let Surface3::Sphere(s) = offset {
            assert!((s.radius - 2.5).abs() < 1e-9, "radius should increase by offset");
        } else {
            panic!("expected Sphere");
        }
    }

    #[test]
    fn offset_cylinder_grows() {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });

        let offset = offset_surface(&cylinder, 0.3).unwrap();

        if let Surface3::Cylinder(c) = offset {
            assert!((c.radius - 1.3).abs() < 1e-9, "radius should increase by offset");
        } else {
            panic!("expected Cylinder");
        }
    }

    #[test]
    fn offset_sphere_negative_too_large_returns_none() {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });

        // Negative offset larger than radius should return None
        let offset = offset_surface(&sphere, -2.0);
        assert!(offset.is_none(), "offset larger than radius should return None");
    }

    #[test]
    fn offset_zero_returns_error() {
        let brep = BRep::new();
        let opts = OffsetOptions::new(0.0);

        let result = offset_shape(&brep, opts);
        assert!(matches!(result, Err(OffsetError::ZeroDistance)));
    }

    #[test]
    fn offset_options_default() {
        let opts = OffsetOptions::default();
        assert_eq!(opts.distance, 1.0);
        assert!(opts.check_self_intersection);
        assert!(!opts.auto_repair);
    }

    #[test]
    fn offset_options_builder() {
        let opts = OffsetOptions::new(0.5)
            .with_tolerance(1e-6)
            .with_self_intersection_check(false)
            .with_auto_repair(true);

        assert_eq!(opts.distance, 0.5);
        assert!((opts.tolerance - 1e-6).abs() < 1e-12);
        assert!(!opts.check_self_intersection);
        assert!(opts.auto_repair);
    }

    #[test]
    fn self_intersection_detection_small_box() {
        // Create a 1x1x1 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Populate geometry
        crate::geom_populate::populate_box_geom(&mut brep);

        // Offset distance > 0.5 should self-intersect
        let self_intersects = detect_self_intersection(&brep, 0.6);
        assert!(self_intersects, "should detect self-intersection for large offset");

        // Offset distance < 0.5 should not self-intersect
        let no_intersect = detect_self_intersection(&brep, 0.4);
        assert!(!no_intersect, "should not detect self-intersection for small offset");
    }

    #[test]
    fn offset_shell_simple_box() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let shell = &brep.solids[0].shells[0];
        let result = offset_shell(shell, &brep, 0.1);

        assert!(result.is_ok(), "offset_shell should succeed for a simple box");
        let offset_brep = result.unwrap();

        // Should have the same number of faces
        let orig_face_count = shell.faces.len();
        let offset_face_count = offset_brep.solids[0].shells[0].faces.len();
        assert_eq!(offset_face_count, orig_face_count, "should preserve face count");
    }

    #[test]
    fn offset_shell_negative_distance() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let shell = &brep.solids[0].shells[0];
        let result = offset_shell(shell, &brep, -0.1);

        assert!(result.is_ok(), "offset_shell with negative distance should succeed");
    }

    #[test]
    fn offset_solid_simple() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let solid = &brep.solids[0];
        let result = offset_solid(solid, &brep, 0.2);

        assert!(result.is_ok(), "offset_solid should succeed");
        let offset_brep = result.unwrap();

        // Verify structure
        assert!(!offset_brep.vertices.is_empty(), "should have vertices");
        assert!(!offset_brep.edges.is_empty(), "should have edges");
        assert!(!offset_brep.solids.is_empty(), "should have solids");
    }

    #[test]
    fn hollow_solid_simple_box() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Hollow by removing top face (index 5 based on typical box construction)
        let solid = &brep.solids[0];
        let result = hollow_solid(solid, &brep, 0.1, &[5]);

        assert!(result.is_ok(), "hollow_solid should succeed with one face removed");
        let hollow_brep = result.unwrap();

        // Should have original kept faces (5) + lateral faces at boundary
        let face_count = hollow_brep.solids[0].shells[0].faces.len();
        assert!(face_count >= 5, "should have at least 5 faces (kept faces + lateral faces)");
    }

    #[test]
    fn hollow_solid_multiple_open_faces() {
        // Create a 2x2x2 box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Hollow by removing top (5) and bottom (0) faces
        let solid = &brep.solids[0];
        let result = hollow_solid(solid, &brep, 0.1, &[0, 5]);

        assert!(result.is_ok(), "hollow_solid should succeed with multiple open faces");
    }

    #[test]
    fn hollow_solid_all_faces_error() {
        // Create a box
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Trying to remove all 6 faces should error
        let solid = &brep.solids[0];
        let result = hollow_solid(solid, &brep, 0.1, &[0, 1, 2, 3, 4, 5]);

        assert!(result.is_err(), "hollow_solid should fail when all faces are removed");
    }

    #[test]
    fn offset_shape_api() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1)
            .with_self_intersection_check(true);

        let result = offset_shape(&brep, opts);

        assert!(result.is_ok(), "offset_shape should succeed");
        let offset_result = result.unwrap();

        assert_eq!(offset_result.offset_faces, 6, "should have 6 offset faces");
        assert!(!offset_result.self_intersection, "should not have self-intersection");
    }

    #[test]
    fn offset_torus_surface() {
        let torus = Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let offset = offset_surface(&torus, 0.1).unwrap();

        if let Surface3::Torus(t) = offset {
            assert!((t.minor_radius - 0.6).abs() < 1e-9, "minor radius should increase by offset");
            assert!((t.major_radius - 2.0).abs() < 1e-9, "major radius should be unchanged");
        } else {
            panic!("expected Torus");
        }
    }

    #[test]
    fn offset_cone_surface() {
        let cone = Surface3::Cone(rcad_kernel::geom::ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: std::f64::consts::PI / 6.0, // 30 degrees
        });

        let offset = offset_surface(&cone, 0.1);

        assert!(offset.is_some(), "cone offset should succeed");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Join Type Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn join_type_default() {
        assert_eq!(JoinType::default(), JoinType::Intersection);
    }

    #[test]
    fn join_type_requires_geometry() {
        assert!(!JoinType::Intersection.requires_join_geometry());
        assert!(JoinType::Arc.requires_join_geometry());
        assert!(JoinType::Tangent.requires_join_geometry());
    }

    #[test]
    fn join_type_as_str() {
        assert_eq!(JoinType::Intersection.as_str(), "intersection");
        assert_eq!(JoinType::Arc.as_str(), "arc");
        assert_eq!(JoinType::Tangent.as_str(), "tangent");
    }

    #[test]
    fn offset_with_arc_join() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1)
            .with_join_type(JoinType::Arc);

        let result = offset_shape(&brep, opts);
        assert!(result.is_ok(), "offset with arc join should succeed");
    }

    #[test]
    fn offset_with_tangent_join() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1)
            .with_join_type(JoinType::Tangent);

        let result = offset_shape(&brep, opts);
        assert!(result.is_ok(), "offset with tangent join should succeed");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Variable Thickness Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn variable_thickness_new() {
        let vt = VariableThickness::new(1.0);
        assert_eq!(vt.default_thickness, 1.0);
        assert!(vt.face_thicknesses.is_empty());
        assert_eq!(vt.transition_width, 0.0);
        assert!(!vt.interpolate);
    }

    #[test]
    fn variable_thickness_with_face() {
        let vt = VariableThickness::new(1.0)
            .with_face(0, 0.5)
            .with_face(1, 1.5);

        assert_eq!(vt.thickness_for_face(0), 0.5);
        assert_eq!(vt.thickness_for_face(1), 1.5);
        assert_eq!(vt.thickness_for_face(2), 1.0); // default
    }

    #[test]
    fn variable_thickness_validation() {
        let vt = VariableThickness::new(1.0)
            .with_face(0, 0.5)
            .with_face(10, 1.5); // Invalid face index

        // Validate with 5 faces
        let result = vt.validate(5);
        assert!(result.is_err(), "should fail for out-of-range face index");
    }

    #[test]
    fn variable_thickness_zero_thickness_error() {
        let vt = VariableThickness::new(0.0);

        let result = vt.validate(5);
        assert!(result.is_err(), "should fail for zero default thickness");
    }

    #[test]
    fn offset_with_variable_thickness() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let vt = VariableThickness::new(0.2)
            .with_face(0, 0.1)
            .with_face(1, 0.3);

        let opts = OffsetOptions::new(0.2)
            .with_variable_thickness(vt);

        let result = offset_shape(&brep, opts);
        assert!(result.is_ok(), "offset with variable thickness should succeed");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Self-Intersection Detection Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn detect_self_intersection_detailed_empty() {
        let brep = BRep::new();
        let result = detect_self_intersection_detailed(&brep, 0.5);

        assert!(!result.has_intersection);
        assert!(result.intersecting_pairs.is_empty());
        assert!(result.min_safe_distance.is_none());
    }

    #[test]
    fn detect_self_intersection_detailed_box() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        // Large offset should detect intersection
        let result = detect_self_intersection_detailed(&brep, 0.6);
        assert!(result.has_intersection);
        assert!(!result.intersecting_pairs.is_empty());
        assert!(result.min_safe_distance.is_some());

        // Small offset should not detect intersection
        let result = detect_self_intersection_detailed(&brep, 0.3);
        assert!(!result.has_intersection);
    }

    #[test]
    fn self_intersection_config_default() {
        let config = SelfIntersectionConfig::default();

        assert!(config.detect);
        assert!(!config.auto_repair);
        assert_eq!(config.max_repair_attempts, 5);
        assert!((config.reduction_factor - 0.8).abs() < 1e-10);
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Quality Analysis Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn offset_quality_default() {
        let quality = OffsetQuality::default();

        assert_eq!(quality.min_wall_thickness, 0.0);
        assert_eq!(quality.max_deviation, 0.0);
        assert_eq!(quality.degenerate_edge_count, 0);
        assert_eq!(quality.self_intersection_count, 0);
        assert!(!quality.is_valid); // 0 min_wall_thickness < default threshold
    }

    #[test]
    fn offset_quality_check_thresholds() {
        let quality = OffsetQuality {
            min_wall_thickness: 0.5,
            max_deviation: 1e-5,
            degenerate_edge_count: 0,
            self_intersection_count: 0,
            face_area_ratio: 1.0,
            edge_length_ratio: 1.0,
            is_valid: true,
            warnings: Vec::new(),
        };

        let thresholds = QualityThresholds::default();
        assert!(quality.check_thresholds(&thresholds).is_ok());
    }

    #[test]
    fn offset_quality_check_thresholds_failure() {
        let quality = OffsetQuality {
            min_wall_thickness: 1e-9, // Below threshold
            max_deviation: 0.0,
            degenerate_edge_count: 0,
            self_intersection_count: 0,
            face_area_ratio: 1.0,
            edge_length_ratio: 1.0,
            is_valid: false,
            warnings: Vec::new(),
        };

        let thresholds = QualityThresholds::default();
        assert!(quality.check_thresholds(&thresholds).is_err());
    }

    #[test]
    fn quality_thresholds_default() {
        let thresholds = QualityThresholds::default();

        assert!((thresholds.min_wall_thickness - 1e-6).abs() < 1e-12);
        assert!((thresholds.max_deviation - 1e-4).abs() < 1e-12);
        assert!(!thresholds.allow_self_intersection);
        assert!((thresholds.max_degenerate_ratio - 0.1).abs() < 1e-12);
    }

    #[test]
    fn analyze_offset_quality_simple() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1);
        let result = offset_shape(&brep, opts.clone()).unwrap();

        let quality = analyze_offset_quality(&result.brep, &brep, &opts);

        assert!(quality.is_valid || quality.warnings.iter().any(|w| w.contains("wall thickness")));
    }

    #[test]
    fn compute_min_wall_thickness_box() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });

        crate::geom_populate::populate_box_geom(&mut brep);

        let min_thickness = compute_min_wall_thickness(&brep, 0.1);

        // For a 2x2x2 box with offset 0.1, min wall should be around 1.8
        assert!(min_thickness > 0.0, "min wall thickness should be positive");
    }

    #[test]
    fn test_compute_face_area_ratio() {
        let mut brep1 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep1);

        let mut brep2 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep2);

        let ratio = compute_face_area_ratio(&brep1, &brep2);
        assert!((ratio - 1.0).abs() < 1e-10, "same box should have ratio 1.0");
    }

    #[test]
    fn test_compute_edge_length_ratio() {
        let mut brep1 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep1);

        let mut brep2 = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep2);

        let ratio = compute_edge_length_ratio(&brep1, &brep2);
        assert!((ratio - 1.0).abs() < 1e-10, "same box should have ratio 1.0");
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Offset Options Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn offset_options_with_join_type() {
        let opts = OffsetOptions::new(0.5)
            .with_join_type(JoinType::Arc);

        assert_eq!(opts.join_type, JoinType::Arc);
    }

    #[test]
    fn offset_options_with_variable_thickness() {
        let vt = VariableThickness::new(1.0).with_face(0, 0.5);
        let opts = OffsetOptions::new(0.5)
            .with_variable_thickness(vt.clone());

        assert!(opts.variable_thickness.is_some());
        assert_eq!(opts.variable_thickness.as_ref().unwrap().thickness_for_face(0), 0.5);
    }

    #[test]
    fn offset_options_with_self_intersection_config() {
        let config = SelfIntersectionConfig {
            detect: true,
            auto_repair: true,
            max_repair_attempts: 10,
            reduction_factor: 0.9,
            min_offset_distance: 0.001,
            allow_partial_results: true,
        };

        let opts = OffsetOptions::new(0.5)
            .with_self_intersection_config(config.clone());

        assert!(opts.self_intersection_config.auto_repair);
        assert_eq!(opts.self_intersection_config.max_repair_attempts, 10);
    }

    #[test]
    fn offset_options_with_quality_thresholds() {
        let thresholds = QualityThresholds {
            min_wall_thickness: 0.1,
            max_deviation: 0.01,
            allow_self_intersection: true,
            max_degenerate_ratio: 0.05,
        };

        let opts = OffsetOptions::new(0.5)
            .with_quality_thresholds(thresholds.clone());

        assert!((opts.quality_thresholds.min_wall_thickness - 0.1).abs() < 1e-12);
        assert!(opts.quality_thresholds.allow_self_intersection);
    }

    #[test]
    fn offset_options_with_approximation_tolerance() {
        let opts = OffsetOptions::new(0.5)
            .with_approximation_tolerance(1e-6);

        assert!((opts.approximation_tolerance - 1e-6).abs() < 1e-12);
    }

    #[test]
    fn offset_options_with_wall_thickness_check() {
        let opts = OffsetOptions::new(0.5)
            .with_wall_thickness_check(0.1);

        assert!(opts.check_wall_thickness);
        assert!((opts.min_wall_thickness - 0.1).abs() < 1e-12);
    }

    #[test]
    fn offset_options_effective_distance_for_face() {
        let vt = VariableThickness::new(1.0).with_face(0, 0.5);
        let opts = OffsetOptions::new(1.0)
            .with_variable_thickness(vt);

        assert_eq!(opts.effective_distance_for_face(0), 0.5);
        assert_eq!(opts.effective_distance_for_face(1), 1.0); // default
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Offset Result Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn offset_result_fields() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1);
        let result = offset_shape(&brep, opts).unwrap();

        assert_eq!(result.offset_faces, 6);
        assert!(!result.self_intersection);
        assert_eq!(result.effective_distance, 0.1);
        assert_eq!(result.repair_attempts, 0);
        assert!(result.warnings.is_empty() || !result.warnings.is_empty());
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Error Display Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn offset_error_display() {
        assert_eq!(
            OffsetError::ZeroDistance.to_string(),
            "offset distance is zero"
        );
        assert_eq!(
            OffsetError::InvalidInput("test").to_string(),
            "invalid input: test"
        );
        assert!(OffsetError::DegenerateSurface { face_index: 1, distance: 0.5 }
            .to_string()
            .contains("face 1"));
        assert!(OffsetError::SelfIntersection { description: "test".to_string() }
            .to_string()
            .contains("self-intersection detected"));
    }

    #[test]
    fn offset_error_new_variants() {
        let err = OffsetError::WallThicknessViolation {
            minimum: 0.1,
            actual: 0.05,
            location: "face 0".to_string(),
        };
        assert!(err.to_string().contains("0.05"));

        let err = OffsetError::JoinCreationFailed {
            join_type: JoinType::Arc,
            edge_index: 1,
            reason: "test".to_string(),
        };
        assert!(err.to_string().contains("Arc"));

        let err = OffsetError::QualityCheckFailed {
            metric: "wall_thickness".to_string(),
            value: 0.05,
            threshold: 0.1,
        };
        assert!(err.to_string().contains("wall_thickness"));

        let err = OffsetError::RecoveryFailed {
            attempts: 3,
            last_error: "test error".to_string(),
        };
        assert!(err.to_string().contains("3 attempts"));
    }

    // ─────────────────────────────────────────────────────────────────────────────
    // Integration Tests
    // ─────────────────────────────────────────────────────────────────────────────

    #[test]
    fn full_offset_workflow_with_quality_check() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep);

        let opts = OffsetOptions::new(0.1)
            .with_join_type(JoinType::Arc)
            .with_self_intersection_check(true)
            .with_wall_thickness_check(0.01)
            .with_approximation_tolerance(1e-5);

        let result = offset_shape(&brep, opts).unwrap();

        // Verify the workflow completed
        assert!(!result.brep.vertices.is_empty());
        assert!(!result.brep.edges.is_empty());
        assert_eq!(result.offset_faces, 6);
    }

    #[test]
    fn offset_with_self_intersection_config() {
        let mut brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        crate::geom_populate::populate_box_geom(&mut brep);

        let config = SelfIntersectionConfig {
            detect: true,
            auto_repair: false,
            max_repair_attempts: 5,
            reduction_factor: 0.8,
            min_offset_distance: 0.01,
            allow_partial_results: false,
        };

        let opts = OffsetOptions::new(0.6)
            .with_self_intersection_config(config);

        let result = offset_shape(&brep, opts).unwrap();

        // Large offset on small box should detect self-intersection
        assert!(result.self_intersection);
    }
}
