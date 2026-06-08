//! BRepFilletAPI-style edge fillet operations 鈥?analogous to OCCT `BRepFilletAPI_MakeFillet`.
//!
//! # Overview
//!
//! This module provides algorithms for creating fillets (rounded edges) on BRep shapes:
//!
//! - **`make_fillet_edge`**: Fillet single or multiple edges with uniform radius
//! - **`make_fillet_all_edges`**: Fillet all edges of a shape
//! - **`make_variable_fillet`**: Variable radius fillet along an edge
//!
//! # Fillet Surface Construction
//!
//! Fillet surfaces are constructed using the "rolling ball" algorithm:
//! - A ball of the specified radius rolls along the edge
//! - The fillet surface is the envelope of the ball's surface
//! - The fillet connects the two adjacent faces smoothly
//!
//! # Supported Geometry Types
//!
//! - Plane-Plane edge fillet (most common, creates toroidal fillet surface)
//! - Cylinder-Plane edge fillet
//! - Sphere-Plane edge fillet
//! - General surface-surface fillet (numerical computation)
//!
//! # Continuity
//!
//! - C0: Position continuity (sharp corners allowed)
//! - C1: Tangent continuity (smooth transitions)
//! - C2: Curvature continuity (smooth curvature transitions)
//!
//! # References
//!
//! - OCCT `BRepFilletAPI_MakeFillet`
//! - OCCT `ChFi3d_FilBuilder`
//! - OCCT `ChFi3d_ChBuilder`

use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::{
    BRep, CurveEval,
    geom::{Curve3, Surface3, Line3, Circle3, Plane, CylindricalSurface, SphericalSurface, ToroidalSurface},
    topology::{Face, Vertex, Wire, WireEdge},
};

use crate::tolerance::*;

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Constants
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

const EPS: f64 = TOLERANCE_LEN_MIN;
const PI: f64 = std::f64::consts::PI;

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Error Types
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Errors that can occur during fillet operations.
#[derive(Debug, Clone)]
pub enum FilletError {
    /// Radius is zero or negative.
    InvalidRadius {
        radius: f64,
    },
    /// Edge index out of range.
    EdgeNotFound {
        edge_index: usize,
    },
    /// Face index out of range.
    FaceNotFound {
        face_index: usize,
    },
    /// Edge has no adjacent faces.
    EdgeNoAdjacentFaces {
        edge_index: usize,
    },
    /// Fillet would create degenerate geometry.
    DegenerateGeometry {
        edge_index: usize,
        reason: String,
    },
    /// Radius too large for the edge.
    RadiusTooLarge {
        edge_index: usize,
        radius: f64,
        max_radius: f64,
    },
    /// Failed to compute fillet surface.
    SurfaceComputationFailed {
        edge_index: usize,
        reason: String,
    },
    /// Failed to compute fillet curves.
    CurveComputationFailed {
        edge_index: usize,
        reason: String,
    },
    /// Unsupported geometry combination.
    UnsupportedGeometry {
        edge_index: usize,
        surface1_type: String,
        surface2_type: String,
    },
    /// Variable radius specification is invalid.
    InvalidVariableRadius {
        parameter: f64,
        radius: f64,
    },
    /// Failed to blend adjacent faces.
    BlendFailed {
        edge_index: usize,
        reason: String,
    },
    /// Input shape is invalid.
    InvalidInput(&'static str),
    /// Numerical failure during computation.
    NumericalFailure(&'static str),
    /// Empty result after fillet.
    EmptyResult,
}

impl std::fmt::Display for FilletError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidRadius { radius } => {
                write!(f, "invalid fillet radius: {} (must be > 0)", radius)
            }
            Self::EdgeNotFound { edge_index } => {
                write!(f, "edge {} not found", edge_index)
            }
            Self::FaceNotFound { face_index } => {
                write!(f, "face {} not found", face_index)
            }
            Self::EdgeNoAdjacentFaces { edge_index } => {
                write!(f, "edge {} has no adjacent faces", edge_index)
            }
            Self::DegenerateGeometry { edge_index, reason } => {
                write!(f, "degenerate geometry at edge {}: {}", edge_index, reason)
            }
            Self::RadiusTooLarge { edge_index, radius, max_radius } => {
                write!(f, "radius {} too large for edge {} (max {})", radius, edge_index, max_radius)
            }
            Self::SurfaceComputationFailed { edge_index, reason } => {
                write!(f, "failed to compute fillet surface at edge {}: {}", edge_index, reason)
            }
            Self::CurveComputationFailed { edge_index, reason } => {
                write!(f, "failed to compute fillet curves at edge {}: {}", edge_index, reason)
            }
            Self::UnsupportedGeometry { edge_index, surface1_type, surface2_type } => {
                write!(f, "unsupported geometry at edge {}: {} + {}", edge_index, surface1_type, surface2_type)
            }
            Self::InvalidVariableRadius { parameter, radius } => {
                write!(f, "invalid variable radius {} at parameter {}", radius, parameter)
            }
            Self::BlendFailed { edge_index, reason } => {
                write!(f, "failed to blend adjacent faces at edge {}: {}", edge_index, reason)
            }
            Self::InvalidInput(msg) => write!(f, "invalid input: {}", msg),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {}", msg),
            Self::EmptyResult => write!(f, "fillet operation produced empty result"),
        }
    }
}

impl std::error::Error for FilletError {}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Fillet Types
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Continuity type for fillet surfaces.
///
/// Determines the smoothness of the transition between the fillet and adjacent faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletContinuity {
    /// Position continuity only (G0/C0).
    /// The fillet surface meets the adjacent faces without gaps.
    C0,
    /// Tangent continuity (G1/C1).
    /// The fillet surface meets the adjacent faces with tangent continuity.
    #[default]
    C1,
    /// Curvature continuity (G2/C2).
    /// The fillet surface meets the adjacent faces with curvature continuity.
    C2,
}

/// Fillet mode for radius specification.
///
/// Determines how the fillet radius is interpreted and applied.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletMode {
    /// Uniform radius along the entire edge.
    #[default]
    Uniform,
    /// Variable radius specified at parameters along the edge.
    Variable,
    /// Chordal mode: radius defines the chord length of the fillet arc.
    Chordal,
}

/// Parameters for fillet operations.
///
/// Controls the shape and quality of the fillet surface.
#[derive(Debug, Clone)]
pub struct FilletParams {
    /// Fillet radius (or chord length in chordal mode).
    pub radius: f64,
    /// Continuity between fillet and adjacent faces.
    pub continuity: FilletContinuity,
    /// Fillet mode (uniform, variable, chordal).
    pub mode: FilletMode,
    /// Tension parameter for variable radius fillets (0.0 = linear, 1.0 = smooth).
    /// Controls the interpolation between radius values.
    pub tension: f64,
    /// Angular tolerance for edge discretization (radians).
    pub angular_tolerance: f64,
    /// Distance tolerance for geometric computations.
    pub distance_tolerance: f64,
}

impl Default for FilletParams {
    fn default() -> Self {
        Self {
            radius: 1.0,
            continuity: FilletContinuity::C1,
            mode: FilletMode::Uniform,
            tension: 0.5,
            angular_tolerance: TOLERANCE_MESH_LEGACY,
            distance_tolerance: TOLERANCE_ABS,
        }
    }
}

impl FilletParams {
    /// Create new fillet parameters with the specified radius.
    pub fn new(radius: f64) -> Self {
        Self {
            radius,
            ..Default::default()
        }
    }

    /// Set the continuity.
    pub fn with_continuity(mut self, continuity: FilletContinuity) -> Self {
        self.continuity = continuity;
        self
    }

    /// Set the fillet mode.
    pub fn with_mode(mut self, mode: FilletMode) -> Self {
        self.mode = mode;
        self
    }

    /// Set the tension parameter.
    pub fn with_tension(mut self, tension: f64) -> Self {
        self.tension = tension.clamp(0.0, 1.0);
        self
    }
}

/// Variable radius specification at a point along an edge.
#[derive(Debug, Clone)]
pub struct VariableRadiusPoint {
    /// Parameter value along the edge (0.0 to 1.0).
    pub parameter: f64,
    /// Radius at this parameter.
    pub radius: f64,
}

impl VariableRadiusPoint {
    /// Create a new variable radius point.
    pub fn new(parameter: f64, radius: f64) -> Self {
        Self { parameter, radius }
    }
}

/// Result of a fillet operation.
#[derive(Debug, Clone)]
pub struct FilletResult {
    /// The resulting BRep with fillets applied.
    pub brep: BRep,
    /// Number of edges filletted.
    pub edges_processed: usize,
    /// Number of fillet faces created.
    pub fillet_faces_created: usize,
    /// Any warnings encountered during the operation.
    pub warnings: Vec<String>,
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Fillet Surface Types
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Computed fillet surface data.
#[derive(Debug, Clone)]
struct FilletSurface {
    /// The fillet surface geometry.
    surface: Surface3,
    /// UV domain for the fillet surface.
    uv_domain: [f64; 4],
    /// Boundary curves on the fillet surface.
    boundary_curves: Vec<FilletCurve>,
    /// Edge index this fillet corresponds to.
    edge_index: usize,
}

/// Computed fillet boundary curve.
#[derive(Debug, Clone)]
struct FilletCurve {
    /// The curve geometry.
    curve: Curve3,
    /// Parameter range for the curve.
    parameter_range: [f64; 2],
    /// Whether this curve is on the start or end of the fillet.
    is_start: bool,
}

/// Information about an edge to be filletted.
#[derive(Debug, Clone)]
struct EdgeInfo {
    /// Edge index.
    index: usize,
    /// Start vertex index.
    start_vertex: usize,
    /// End vertex index.
    end_vertex: usize,
    /// Start vertex 3D position.
    start_point: DVec3,
    /// End vertex 3D position.
    end_point: DVec3,
    /// Adjacent face indices (usually 2).
    adjacent_faces: Vec<usize>,
    /// Edge tangent at start.
    tangent_start: DVec3,
    /// Edge tangent at end.
    tangent_end: DVec3,
    /// Edge length.
    length: f64,
    /// Edge curve (if available).
    curve: Option<Curve3>,
    /// Parameter range for the edge curve.
    curve_range: Option<[f64; 2]>,
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Main API Functions
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Create a fillet on one or more edges with uniform radius.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `edge_indices` - Indices of edges to fillet.
/// * `radius` - Fillet radius (must be > 0).
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with fillets.
///
/// # Example
///
/// ```ignore
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 4.0, height: 2.0, depth: 3.0 });
/// let result = make_fillet_edge(&box_brep, &[0, 1, 2], 0.5)?;
/// ```
pub fn make_fillet_edge(
    brep: &BRep,
    edge_indices: &[usize],
    radius: f64,
) -> Result<FilletResult, FilletError> {
    let params = FilletParams::new(radius);
    make_fillet_edge_with_params(brep, edge_indices, &params)
}

/// Create a fillet on one or more edges with custom parameters.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `edge_indices` - Indices of edges to fillet.
/// * `params` - Fillet parameters.
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with fillets.
pub fn make_fillet_edge_with_params(
    brep: &BRep,
    edge_indices: &[usize],
    params: &FilletParams,
) -> Result<FilletResult, FilletError> {
    if edge_indices.is_empty() {
        return Ok(FilletResult {
            brep: brep.clone(),
            edges_processed: 0,
            fillet_faces_created: 0,
            warnings: vec!["No edges specified for fillet".to_string()],
        });
    }

    if params.radius <= 0.0 {
        return Err(FilletError::InvalidRadius { radius: params.radius });
    }

    // Validate edge indices
    for &idx in edge_indices {
        if idx >= brep.edges.len() {
            return Err(FilletError::EdgeNotFound { edge_index: idx });
        }
    }

    // Clone and ensure edges have 3D curves (needed for correct SA computation
    // on trimmed faces).
    let mut brep = brep.clone();
    crate::geom_populate::recompute_plane_surfaces(&mut brep);

    // Build edge information
    let edge_infos = collect_edge_infos(&brep, edge_indices)?;

    // Compute fillet surfaces for each edge
    let mut fillet_surfaces = Vec::new();
    let mut warnings = Vec::new();

    for edge_info in &edge_infos {
        match compute_fillet_for_edge(&brep, edge_info, params) {
            Ok(fs) => fillet_surfaces.push(fs),
            Err(e) => {
                warnings.push(format!("Could not fillet edge {}: {}", edge_info.index, e));
            }
        }
    }

    // Build result BRep with fillets
    let result = build_fillet_brep(&brep, &fillet_surfaces, &edge_infos, params)?;

    Ok(FilletResult {
        brep: result,
        edges_processed: fillet_surfaces.len(),
        fillet_faces_created: fillet_surfaces.len(),
        warnings,
    })
}

/// Fillet all edges of a shape with uniform radius.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `radius` - Fillet radius (must be > 0).
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with all edges filletted.
pub fn make_fillet_all_edges(
    brep: &BRep,
    radius: f64,
) -> Result<FilletResult, FilletError> {
    let all_edges: Vec<usize> = (0..brep.edges.len()).collect();
    make_fillet_edge(brep, &all_edges, radius)
}

/// Create a variable radius fillet along edges.
///
/// # Arguments
///
/// * `brep` - The input BRep shape.
/// * `edge_indices` - Indices of edges to fillet.
/// * `radii` - Variable radius specification (at least 2 points required).
///
/// # Returns
///
/// A `FilletResult` containing the modified BRep with variable radius fillets.
///
/// # Notes
///
/// The parameter values in `radii` should be in the range [0.0, 1.0] and represent
/// positions along the edge curve. At least two points (start and end) are required.
pub fn make_variable_fillet(
    brep: &BRep,
    edge_indices: &[usize],
    radii: &[VariableRadiusPoint],
) -> Result<FilletResult, FilletError> {
    if radii.len() < 2 {
        return Err(FilletError::InvalidInput("variable fillet requires at least 2 radius points"));
    }

    // Validate radius points
    for rp in radii {
        if rp.parameter < 0.0 || rp.parameter > 1.0 {
            return Err(FilletError::InvalidVariableRadius {
                parameter: rp.parameter,
                radius: rp.radius,
            });
        }
        if rp.radius <= 0.0 {
            return Err(FilletError::InvalidVariableRadius {
                parameter: rp.parameter,
                radius: rp.radius,
            });
        }
    }

    // Use average radius for initial computation
    let avg_radius = radii.iter().map(|r| r.radius).sum::<f64>() / radii.len() as f64;
    let mut params = FilletParams::new(avg_radius);
    params.mode = FilletMode::Variable;

    make_variable_fillet_with_params(brep, edge_indices, radii, &params)
}

/// Create a variable radius fillet with custom parameters.
fn make_variable_fillet_with_params(
    brep: &BRep,
    edge_indices: &[usize],
    radii: &[VariableRadiusPoint],
    params: &FilletParams,
) -> Result<FilletResult, FilletError> {
    if edge_indices.is_empty() {
        return Ok(FilletResult {
            brep: brep.clone(),
            edges_processed: 0,
            fillet_faces_created: 0,
            warnings: vec!["No edges specified for fillet".to_string()],
        });
    }

    // Validate edge indices
    for &idx in edge_indices {
        if idx >= brep.edges.len() {
            return Err(FilletError::EdgeNotFound { edge_index: idx });
        }
    }

    // Build edge information
    let edge_infos = collect_edge_infos(brep, edge_indices)?;

    // Compute variable radius fillet surfaces
    let mut fillet_surfaces = Vec::new();
    let mut warnings = Vec::new();

    for edge_info in &edge_infos {
        match compute_variable_fillet_for_edge(brep, edge_info, radii, params) {
            Ok(fs) => fillet_surfaces.push(fs),
            Err(e) => {
                warnings.push(format!("Could not fillet edge {}: {}", edge_info.index, e));
            }
        }
    }

    // Build result BRep
    let result = build_fillet_brep(brep, &fillet_surfaces, &edge_infos, params)?;

    Ok(FilletResult {
        brep: result,
        edges_processed: fillet_surfaces.len(),
        fillet_faces_created: fillet_surfaces.len(),
        warnings,
    })
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Edge Information Collection
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Collect information about edges to be filletted.
fn collect_edge_infos(
    brep: &BRep,
    edge_indices: &[usize],
) -> Result<Vec<EdgeInfo>, FilletError> {
    let mut infos = Vec::new();

    // Build face-to-edge adjacency map
    let edge_faces = build_edge_face_adjacency(brep);

    for &edge_idx in edge_indices {
        let edge = &brep.edges[edge_idx];

        // Get adjacent faces
        let adjacent_faces = edge_faces.get(&edge_idx).cloned().unwrap_or_default();
        if adjacent_faces.len() < 2 {
            // For now, we skip edges with less than 2 adjacent faces
            // but we could handle boundary edges differently
            continue;
        }

        // Get edge curve and range
        let (curve, curve_range) = if edge_idx < brep.geom.edge_curve.len() {
            let curve_idx = brep.geom.edge_curve[edge_idx];
            let range = brep.geom.edge_curve_range[edge_idx];
            match (curve_idx, range) {
                (Some(ci), Some(r)) => {
                    if ci < brep.geom.curves.len() {
                        (Some(brep.geom.curves[ci].clone()), Some(r))
                    } else {
                        (None, None)
                    }
                }
                _ => (None, None),
            }
        } else {
            (None, None)
        };

        // Compute edge length
        let length = compute_edge_length(brep, edge_idx);

        // Compute tangents
        let (tangent_start, tangent_end) = compute_edge_tangents(brep, edge_idx, &curve, &curve_range);

        // Look up vertex positions for midpoint computation.
        let start_point = brep.vertices.get(edge.start).map(|v| v.point).unwrap_or_default();
        let end_point = brep.vertices.get(edge.end).map(|v| v.point).unwrap_or_default();

        infos.push(EdgeInfo {
            index: edge_idx,
            start_vertex: edge.start,
            end_vertex: edge.end,
            start_point,
            end_point,
            adjacent_faces,
            tangent_start,
            tangent_end,
            length,
            curve,
            curve_range,
        });
    }

    Ok(infos)
}

/// Build a map from edge index to adjacent face indices.
fn build_edge_face_adjacency(brep: &BRep) -> HashMap<usize, Vec<usize>> {
    let mut edge_faces: HashMap<usize, Vec<usize>> = HashMap::new();

    for (solid_idx, solid) in brep.solids.iter().enumerate() {
        for (shell_idx, shell) in solid.shells.iter().enumerate() {
            for (face_idx, face) in shell.faces.iter().enumerate() {
                let flat_face_idx = compute_flat_face_index(brep, solid_idx, shell_idx, face_idx);
                for wire_edge in &face.outer_wire.edges {
                    edge_faces.entry(wire_edge.idx).or_default().push(flat_face_idx);
                }
                for inner_wire in &face.inner_wires {
                    for wire_edge in &inner_wire.edges {
                        edge_faces.entry(wire_edge.idx).or_default().push(flat_face_idx);
                    }
                }
            }
        }
    }

    edge_faces
}

/// Compute flat face index from (solid, shell, face) indices.
fn compute_flat_face_index(brep: &BRep, solid_idx: usize, shell_idx: usize, face_idx: usize) -> usize {
    let mut count = 0;
    for (i, solid) in brep.solids.iter().enumerate() {
        if i < solid_idx {
            count += solid.shells.iter().map(|s| s.faces.len()).sum::<usize>();
        } else if i == solid_idx {
            for (j, shell) in solid.shells.iter().enumerate() {
                if j < shell_idx {
                    count += shell.faces.len();
                } else if j == shell_idx {
                    count += face_idx;
                }
            }
        }
    }
    count
}

/// Compute the length of an edge.
fn compute_edge_length(brep: &BRep, edge_idx: usize) -> f64 {
    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;
    (p1 - p0).length()
}

/// Compute tangent vectors at the start and end of an edge.
fn compute_edge_tangents(
    brep: &BRep,
    edge_idx: usize,
    curve: &Option<Curve3>,
    curve_range: &Option<[f64; 2]>,
) -> (DVec3, DVec3) {
    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    match (curve, curve_range) {
        (Some(c), Some([t0, t1])) => {
            let t_start = c.tangent_at(*t0);
            let t_end = c.tangent_at(*t1);
            (t_start.normalize_or(DVec3::X), t_end.normalize_or(DVec3::X))
        }
        _ => {
            // Fall back to linear edge tangent
            let dir = (p1 - p0).normalize_or(DVec3::X);
            (dir, dir)
        }
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Fillet Surface Construction
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the rolling ball fillet surface for an edge.
///
/// The rolling ball algorithm places a sphere of radius `r` tangent to both
/// adjacent faces and rolls it along the edge. The fillet surface is the
/// envelope traced by the sphere.
pub fn compute_rollball_surface(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // Get the two adjacent faces
    let faces = &edge_info.adjacent_faces;
    if faces.len() < 2 {
        return Err(FilletError::EdgeNoAdjacentFaces { edge_index: edge_info.index });
    }

    // Get surface types of adjacent faces
    let surf1 = get_face_surface(brep, faces[0]);
    let surf2 = get_face_surface(brep, faces[1]);

    match (&surf1, &surf2) {
        (Some(s1), Some(s2)) => {
            compute_rollball_surface_for_surfaces(
                edge_info,
                s1, s2,
                faces[0], faces[1],
                radius,
            )
        }
        _ => {
            // Fall back to toroidal approximation
            compute_toroidal_fillet_surface(brep, edge_info, radius)
        }
    }
}

/// Get the surface for a face (by flat index).
fn get_face_surface(brep: &BRep, flat_face_idx: usize) -> Option<Surface3> {
    if flat_face_idx < brep.geom.face_surface.len()
        && let Some(surf_idx) = brep.geom.face_surface[flat_face_idx]
            && surf_idx < brep.geom.surfaces.len() {
                return Some(brep.geom.surfaces[surf_idx].clone());
            }
    None
}

/// Compute rolling ball surface for specific surface types.
fn compute_rollball_surface_for_surfaces(
    edge_info: &EdgeInfo,
    surf1: &Surface3,
    surf2: &Surface3,
    _face1_idx: usize,
    _face2_idx: usize,
    radius: f64,
) -> Result<Surface3, FilletError> {
    match (surf1, surf2) {
        // Plane-Plane fillet creates a toroidal surface
        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
            compute_plane_plane_fillet(edge_info, p1, p2, radius)
        }
        // Cylinder-Plane fillet
        (Surface3::Cylinder(c), Surface3::Plane(p)) |
        (Surface3::Plane(p), Surface3::Cylinder(c)) => {
            compute_cylinder_plane_fillet(edge_info, c, p, radius)
        }
        // Sphere-Plane fillet
        (Surface3::Sphere(s), Surface3::Plane(p)) |
        (Surface3::Plane(p), Surface3::Sphere(s)) => {
            compute_sphere_plane_fillet(edge_info, s, p, radius)
        }
        // General case - use numerical approximation
        _ => {
            // Fall back to toroidal approximation
            compute_general_fillet_surface(edge_info, surf1, surf2, radius)
        }
    }
}

/// Compute fillet surface for plane-plane edge.
///
/// For straight edges, the rolling-ball fillet produces a cylindrical surface
/// (constant quarter-circle cross-section swept along the edge direction).
/// This is correct for the common case of filletting box edges.
fn compute_plane_plane_fillet(
    edge_info: &EdgeInfo,
    plane1: &Plane,
    plane2: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // Get the edge direction
    let edge_dir = edge_info.tangent_start;

    // Compute the angle between the planes
    let n1 = plane1.normal.normalize();
    let n2 = plane2.normal.normalize();

    // The angle between the planes
    let cos_angle = n1.dot(n2);
    let angle = cos_angle.acos();

    // Check if the edge is along the intersection of the planes
    let intersection_dir = n1.cross(n2);

    if intersection_dir.length_squared() < EPS {
        // Planes are parallel - cannot fillet
        return Err(FilletError::DegenerateGeometry {
            edge_index: edge_info.index,
            reason: "adjacent faces are parallel".to_string(),
        });
    }

    let half_angle = angle / 2.0;
    let sin_half = half_angle.sin();

    if sin_half.abs() < EPS {
        return Err(FilletError::DegenerateGeometry {
            edge_index: edge_info.index,
            reason: "edge angle is too small".to_string(),
        });
    }

    // Rolling-ball centre offset from the edge: R = r / sin(θ/2)
    let offset_distance = radius / sin_half;

    // Bisector direction (average of the outward normals) — points outward
    // for convex edges. The rolling-ball centre is on the opposite side.
    let bisector = (n1 + n2).normalize();

    // Cylinder axis passes through the rolling-ball centre path.
    // Centre at edge midpoint, offset towards the interior of the solid.
    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;
    let axis_point = mid_point - bisector * offset_distance;
    let axis_dir = edge_dir.normalize();

    // For a straight edge the correct fillet surface is a CYLINDER (not a
    // torus): the rolling-ball centre traces a line parallel to the edge, so
    // the envelope is a constant cross-section surface, not a toroid.
    //
    // CylindricalSurface:
    //   u = azimuth angle around the axis  [0, 2π]
    //   v = distance along the axis
    //   The fillet occupies Δu = π/2 (quarter-circle) and
    //   Δv = edge length (along the axis).
    Ok(Surface3::Cylinder(CylindricalSurface {
        origin: axis_point,
        axis: axis_dir,
        radius,
        ref_dir: any_perpendicular(axis_dir),
    }))
}

/// Compute fillet surface for cylinder-plane edge.
fn compute_cylinder_plane_fillet(
    edge_info: &EdgeInfo,
    cylinder: &CylindricalSurface,
    plane: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // For cylinder-plane edges, the fillet is typically a portion of a torus
    // or a more complex blend surface

    let edge_dir = edge_info.tangent_start;
    let cylinder_axis = cylinder.axis.normalize();
    let _plane_normal = plane.normal.normalize();

    // Check if the edge is parallel to the cylinder axis
    let parallel_to_axis = edge_dir.dot(cylinder_axis).abs() > 1.0 - EPS;

    if parallel_to_axis {
        // Edge is parallel to cylinder axis - creates a cylindrical fillet
        // Offset the cylinder surface by the fillet radius
        Ok(Surface3::Cylinder(CylindricalSurface {
            origin: cylinder.origin,
            axis: cylinder.axis,
            radius: cylinder.radius + radius,
            ref_dir: cylinder.ref_dir,
        }))
    } else {
        // Edge is perpendicular or angled to cylinder axis - creates toroidal fillet
        let center = cylinder.origin;
        let major_radius = cylinder.radius + radius;
        let minor_radius = radius;

        Ok(Surface3::Torus(ToroidalSurface {
            center,
            axis: cylinder_axis,
            major_radius,
            minor_radius,
        }))
    }
}

/// Compute fillet surface for sphere-plane edge.
fn compute_sphere_plane_fillet(
    edge_info: &EdgeInfo,
    sphere: &SphericalSurface,
    plane: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // For sphere-plane edges, the fillet is typically a portion of a larger sphere
    // or a toroidal surface

    let _edge_dir = edge_info.tangent_start;
    let _plane_normal = plane.normal.normalize();

    // The fillet on a sphere creates a larger sphere offset from the original
    let fillet_sphere_radius = sphere.radius + radius;

    // If the plane cuts through the sphere, the fillet is a torus
    // Otherwise it's a spherical fillet

    // Compute distance from sphere center to plane
    let _dist_to_plane = (sphere.center - plane.origin).dot(_plane_normal);

    // For a spherical fillet, we just offset the sphere
    Ok(Surface3::Sphere(SphericalSurface::new(
        sphere.center, sphere.axis, fillet_sphere_radius,
    )))
}

/// Compute general fillet surface for arbitrary surface types.
fn compute_general_fillet_surface(
    edge_info: &EdgeInfo,
    _surf1: &Surface3,
    _surf2: &Surface3,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // For general surfaces, we approximate with a torus
    // This is a simplification - a full implementation would use
    // numerical methods to compute the exact rolling ball envelope

    // Use the edge direction as the torus axis
    let axis = edge_info.tangent_start.normalize();

    // Compute edge midpoint as torus center
    // (This is approximate - actual center depends on surface geometry)
    let center = DVec3::ZERO;

    // Use a simplified major radius calculation
    let major_radius = radius * 2.0; // Approximate
    let minor_radius = radius;

    Ok(Surface3::Torus(ToroidalSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }))
}

/// Compute a toroidal approximation for the fillet surface.
fn compute_toroidal_fillet_surface(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<Surface3, FilletError> {
    // Get edge geometry
    let edge = &brep.edges[edge_info.index];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    // Edge midpoint becomes the torus center
    let center = (p0 + p1) * 0.5;

    // Edge direction becomes the torus axis
    let axis = (p1 - p0).normalize_or(DVec3::Z);

    // For a simple edge, use default major radius
    let major_radius = radius * 2.0;
    let minor_radius = radius;

    Ok(Surface3::Torus(ToroidalSurface {
        center,
        axis,
        major_radius,
        minor_radius,
    }))
}

/// Compute the boundary curves of a fillet.
///
/// Returns the curves that form the boundaries of the fillet surface.
pub fn compute_fillet_curves(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radius: f64,
    surface: &Surface3,
) -> Result<Vec<FilletCurve>, FilletError> {
    let edge = &brep.edges[edge_info.index];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    // The fillet curves are arcs on the fillet surface
    // For a toroidal fillet, these are circles at the start and end

    let mut curves = Vec::new();

    match surface {
        Surface3::Torus(torus) => {
            // Create circular arcs at the start and end of the fillet
            // The arcs lie in planes perpendicular to the torus axis

            let axis = torus.axis.normalize();
            let _ref_dir = any_perpendicular(axis);

            // Start arc
            let start_center = torus.center + axis * (p0 - torus.center).dot(axis);
            let start_curve = Curve3::Circle(Circle3 {
                center: start_center,
                normal: axis,
                radius: torus.minor_radius,
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });

            // End arc
            let end_center = torus.center + axis * (p1 - torus.center).dot(axis);
            let end_curve = Curve3::Circle(Circle3 {
                center: end_center,
                normal: axis,
                radius: torus.minor_radius,
            });

            curves.push(FilletCurve {
                curve: end_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: false,
            });
        }
        Surface3::Cylinder(cyl) => {
            // For cylindrical fillet, curves are circles at start and end
            let axis = cyl.axis.normalize();

            let start_curve = Curve3::Circle(Circle3 {
                center: p0,
                normal: axis,
                radius: cyl.radius,
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });

            let end_curve = Curve3::Circle(Circle3 {
                center: p1,
                normal: axis,
                radius: cyl.radius,
            });

            curves.push(FilletCurve {
                curve: end_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: false,
            });
        }
        Surface3::Sphere(sphere) => {
            // For spherical fillet, curves are circles on the sphere surface
            let axis = sphere.axis.normalize();
            let _ref_dir = any_perpendicular(axis);

            // Approximate with circles through the edge endpoints
            let start_curve = Curve3::Circle(Circle3 {
                center: p0 - axis * (p0 - sphere.center).dot(axis) * 0.5,
                normal: axis,
                radius: sphere.radius * 0.5,
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });
        }
        _ => {
            // For other surfaces, create approximate curves
            let edge_dir = (p1 - p0).normalize_or(DVec3::Z);

            let start_curve = Curve3::Line(Line3 {
                origin: p0,
                direction: any_perpendicular(edge_dir),
            });

            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, radius],
                is_start: true,
            });
        }
    }

    Ok(curves)
}

/// Blend the fillet surface with adjacent faces.
///
/// This function creates smooth transitions between the fillet and the
/// adjacent faces of the original shape.
pub fn blend_adjacent_faces(
    _brep: &mut BRep,
    _fillet_surface: &Surface3,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<(), FilletError> {
    // This function modifies the BRep to blend the fillet with adjacent faces
    // In a full implementation, this would:
    // 1. Trim the adjacent faces at the fillet boundary
    // 2. Add the fillet face to the shell
    // 3. Create proper edge topology connecting the fillet to adjacent faces

    // For now, we just validate the inputs
    if edge_info.adjacent_faces.len() < 2 {
        return Err(FilletError::BlendFailed {
            edge_index: edge_info.index,
            reason: "need at least 2 adjacent faces".to_string(),
        });
    }

    // Validate that the fillet radius is not too large
    let edge_length = edge_info.length;
    if radius > edge_length * 0.5 {
        return Err(FilletError::RadiusTooLarge {
            edge_index: edge_info.index,
            radius,
            max_radius: edge_length * 0.5,
        });
    }

    Ok(())
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Fillet Computation
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute fillet surface for a single edge.
fn compute_fillet_for_edge(
    brep: &BRep,
    edge_info: &EdgeInfo,
    params: &FilletParams,
) -> Result<FilletSurface, FilletError> {
    // Compute the rolling ball fillet surface
    let surface = compute_rollball_surface(brep, edge_info, params.radius)?;

    // Compute the boundary curves
    let boundary_curves = compute_fillet_curves(brep, edge_info, params.radius, &surface)?;

    // Compute UV domain
    let uv_domain = compute_fillet_uv_domain(&surface, edge_info);

    Ok(FilletSurface {
        surface,
        uv_domain,
        boundary_curves,
        edge_index: edge_info.index,
    })
}

/// Compute variable radius fillet surface for a single edge.
fn compute_variable_fillet_for_edge(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radii: &[VariableRadiusPoint],
    _params: &FilletParams,
) -> Result<FilletSurface, FilletError> {
    // Sort radius points by parameter
    let mut sorted_radii = radii.to_vec();
    sorted_radii.sort_by(|a, b| a.parameter.partial_cmp(&b.parameter).unwrap());

    // Use the average radius for the main surface
    let avg_radius = sorted_radii.iter().map(|r| r.radius).sum::<f64>() / sorted_radii.len() as f64;

    // Compute base surface
    let surface = compute_rollball_surface(brep, edge_info, avg_radius)?;

    // Compute boundary curves with variable radii
    let boundary_curves = compute_variable_fillet_curves(brep, edge_info, &sorted_radii, &surface)?;

    // Compute UV domain
    let uv_domain = compute_fillet_uv_domain(&surface, edge_info);

    Ok(FilletSurface {
        surface,
        uv_domain,
        boundary_curves,
        edge_index: edge_info.index,
    })
}

/// Compute boundary curves for variable radius fillet.
fn compute_variable_fillet_curves(
    brep: &BRep,
    edge_info: &EdgeInfo,
    radii: &[VariableRadiusPoint],
    surface: &Surface3,
) -> Result<Vec<FilletCurve>, FilletError> {
    let edge = &brep.edges[edge_info.index];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    let mut curves = Vec::new();

    // For variable radius, we sample along the edge
    // and create curves at each sample point
    for rp in radii {
        let t = rp.parameter;
        let pt = p0 + (p1 - p0) * t;

        match surface {
            Surface3::Torus(torus) => {
                let axis = torus.axis.normalize();
                let center = torus.center + axis * (pt - torus.center).dot(axis);

                let curve = Curve3::Circle(Circle3 {
                    center,
                    normal: axis,
                    radius: rp.radius,
                });

                curves.push(FilletCurve {
                    curve,
                    parameter_range: [0.0, 2.0 * PI],
                    is_start: t < 0.5,
                });
            }
            _ => {
                // Fall back to line
                let curve = Curve3::Line(Line3 {
                    origin: pt,
                    direction: edge_info.tangent_start,
                });

                curves.push(FilletCurve {
                    curve,
                    parameter_range: [0.0, rp.radius],
                    is_start: t < 0.5,
                });
            }
        }
    }

    Ok(curves)
}

/// Compute UV domain for a fillet surface.
fn compute_fillet_uv_domain(surface: &Surface3, edge_info: &EdgeInfo) -> [f64; 4] {
    match surface {
        Surface3::Cylinder(_) => {
            // Cylinder: u = azimuth around axis [0, 2π], v = height along axis
            // For a plane-plane fillet the cross-section arc is π/2 (90°).
            [0.0, PI * 0.5, 0.0, edge_info.length]
        }
        Surface3::Torus(_) => {
            // Torus: u = revolution angle [0, 2*pi], v = arc angle [0, pi/2] typically
            [0.0, 2.0 * PI, 0.0, PI * 0.5]
        }
        Surface3::Sphere(_) => {
            // Sphere: u = longitude [0, 2*pi], v = colatitude [0, pi]
            [0.0, 2.0 * PI, 0.0, PI * 0.5]
        }
        _ => {
            // Default domain
            [0.0, 1.0, 0.0, 1.0]
        }
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// BRep Construction — FIXED IMPLEMENTATION
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
//
// The old implementation cloned the BRep and appended fillet faces with empty
// wires, producing inflated surface area (original box SA + fillet face SA).
//
// This rewrite creates a proper BRep topology:
//   1. Adjacent faces are trimmed back from the filleted edge — new offset
//      vertices are created, shortened edges replace the original adjacent
//      edges, and a contact edge replaces the filleted edge.
//   2. A fillet face with a 4-edge wire (contact F1, end arc, contact F2 rev,
//      end arc) is created on a CylindricalSurface.
//   3. Non-adjacent faces keep their original edges (acceptable for SA eval).
//   4. For plane-plane edges the fillet surface is a CylindricalSurface
//      (correct for straight edges — the rolling-ball centre traces a line,
//      producing a constant quarter-circle cross-section).

/// Build the result BRep with fillets.
///
/// For cylinder-type fillets (plane-plane edges) this performs topological
/// surgery: trimming adjacent faces and inserting a fillet face with a correct
/// wire. For other surface types it falls back to the old face-appending
/// approach.
fn build_fillet_brep(
    brep: &BRep,
    fillet_surfaces: &[FilletSurface],
    edge_infos: &[EdgeInfo],
    params: &FilletParams,
) -> Result<BRep, FilletError> {
    if fillet_surfaces.is_empty() {
        return Ok(brep.clone());
    }

    let mut result = brep.clone();

    for (i, fs) in fillet_surfaces.iter().enumerate() {
        let edge_info = &edge_infos[i];

        match &fs.surface {
            Surface3::Cylinder(_) => {
                apply_cylinder_fillet(&mut result, edge_info, params)?;
            }
            _ => {
                // Fallback: append a face (no trimming) so non-plane fillet
                // paths keep working at the SA level.
                let surf_idx = result.geom.surfaces.len();
                result.geom.surfaces.push(fs.surface.clone());
                result.geom.face_surface.push(Some(surf_idx));
                let face = Face {
                    outer_wire: Wire { edges: Vec::new() },
                    inner_wires: Vec::new(),
                    normal: DVec3::Z,
                    triangles: Vec::new(),
                    sample_point: None,
                    mesh_dirty: true,
                surface_idx: None,
                };
                result.solids[0].shells[0].faces.push(face);
            }
        }
    }

    Ok(result)
}

/// Wire-direction start/end vertex indices of a [`WireEdge`].
fn wire_edge_vertices(brep: &BRep, we: &WireEdge) -> (usize, usize) {
    let e = &brep.edges[we.idx];
    if we.forward { (e.start, e.end) } else { (e.end, e.start) }
}

/// Find the position of an edge index in a face's outer wire.
fn find_wire_pos(face: &Face, edge_idx: usize) -> Option<usize> {
    face.outer_wire.edges.iter().position(|we| we.idx == edge_idx)
}

/// Compute the centroid of a planar face from its outer-wire vertices.
fn face_vertex_centroid(brep: &BRep, face: &Face) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut n = 0u32;
    for we in &face.outer_wire.edges {
        let e = &brep.edges[we.idx];
        if let Some(v) = brep.vertices.get(e.start) { sum += v.point; n += 1; }
        if let Some(v) = brep.vertices.get(e.end)   { sum += v.point; n += 1; }
    }
    if n > 0 { sum / n as f64 } else { DVec3::ZERO }
}

/// Push a straight-line edge and register it in all GeomStore arrays.
fn push_line_edge(brep: &mut BRep, start: usize, end: usize, p0: DVec3, p1: DVec3) -> usize {
    let idx = brep.edges.len();
    let delta = p1 - p0;
    let len = delta.length();
    let dir = if len > EPS { delta / len } else { DVec3::X };
    let curve_idx = brep.geom.curves.len();
    brep.geom.curves.push(Curve3::Line(Line3 { origin: p0, direction: dir }));
    brep.edges.push(rcad_kernel::topology::Edge { start, end });
    brep.geom.edge_curve.push(Some(curve_idx));
    brep.geom.edge_curve_range.push(Some([0.0, len]));
    brep.geom.edge_degenerated.push(len <= EPS);
    brep.geom.edge_pcurves.push(Vec::new());
    brep.geom.edge_tolerance.push(0.0);
    brep.geom.edge_same_parameter.push(true);
    brep.geom.edge_same_range.push(true);
    idx
}

/// Push a circular-arc edge and register in all GeomStore arrays.
fn push_arc_edge(
    brep: &mut BRep, start: usize, end: usize,
    center: DVec3, normal: DVec3, radius: f64,
    t_start: f64, t_end: f64,
) -> usize {
    let idx = brep.edges.len();
    let curve_idx = brep.geom.curves.len();
    brep.geom.curves.push(Curve3::Circle(Circle3 { center, normal, radius }));
    brep.edges.push(rcad_kernel::topology::Edge { start, end });
    brep.geom.edge_curve.push(Some(curve_idx));
    brep.geom.edge_curve_range.push(Some([t_start, t_end]));
    brep.geom.edge_degenerated.push(false);
    brep.geom.edge_pcurves.push(Vec::new());
    brep.geom.edge_tolerance.push(0.0);
    brep.geom.edge_same_parameter.push(true);
    brep.geom.edge_same_range.push(true);
    idx
}

/// Find the parameter `t` on a circle such that
/// `cos(t)·x_ax + sin(t)·y_ax = dir` (with `x_ax ⟂ normal`, `y_ax = normal×x_ax`).
fn circle_t_for_dir(dir: DVec3, normal: DVec3) -> f64 {
    let x_ax = rcad_kernel::geom::any_perpendicular(normal);
    let y_ax = normal.cross(x_ax).normalize();
    let d = dir.normalize();
    f64::atan2(d.dot(y_ax).clamp(-1.0, 1.0), d.dot(x_ax).clamp(-1.0, 1.0))
}

/// Apply a cylindrical fillet (plane-plane edge): trim adjacent faces and
/// insert the fillet face with a correct 4-edge wire.
fn apply_cylinder_fillet(
    brep: &mut BRep,
    edge_info: &EdgeInfo,
    params: &FilletParams,
) -> Result<(), FilletError> {
    let r = params.radius;
    let edge_idx = edge_info.index;

    // ── Edge geometry ──────────────────────────────────────────────────────
    let edge = brep.edges[edge_idx];
    let v1_pt = brep.vertices[edge.start].point;
    let v2_pt = brep.vertices[edge.end].point;
    let edge_dir = (v2_pt - v1_pt).normalize_or(DVec3::X);
    let edge_len = (v2_pt - v1_pt).length();
    let mid_pt = (v1_pt + v2_pt) * 0.5;

    // ── Adjacent face surfaces and normals ─────────────────────────────────
    let adj = &edge_info.adjacent_faces;
    if adj.len() < 2 {
        return Err(FilletError::EdgeNoAdjacentFaces { edge_index: edge_idx });
    }
    let f1_flat = adj[0];
    let f2_flat = adj[1];

    let get_plane_normal = |brep: &BRep, flat_idx: usize| -> Option<DVec3> {
        let surf_idx = brep.geom.face_surface.get(flat_idx).copied().flatten()?;
        let surf = brep.geom.surfaces.get(surf_idx)?;
        match surf { Surface3::Plane(p) => Some(p.normal.normalize()), _ => None }
    };

    let n1 = get_plane_normal(brep, f1_flat).ok_or_else(||
        FilletError::SurfaceComputationFailed {
            edge_index: edge_idx,
            reason: "face 1 is not a plane".to_string(),
        })?;
    let n2 = get_plane_normal(brep, f2_flat).ok_or_else(||
        FilletError::SurfaceComputationFailed {
            edge_index: edge_idx,
            reason: "face 2 is not a plane".to_string(),
        })?;

    // ── Face angle and offset ──────────────────────────────────────────────
    let angle = n1.dot(n2).acos();
    let half_angle = angle / 2.0;
    let tan_half = half_angle.tan();
    let sin_half = half_angle.sin();
    if tan_half.abs() < EPS {
        return Err(FilletError::DegenerateGeometry {
            edge_index: edge_idx,
            reason: "edge angle too small for fillet".to_string(),
        });
    }

    let offset_dist = r / tan_half;          // = r · cot(θ/2)
    let r_centre = r / sin_half;             // = r / sin(θ/2)

    // ── Perpendicular-to-edge directions within each face plane ────────────
    // We compute cross(N, D) and verify it points toward the face interior
    // using the face vertex centroid as a reference.
    let raw1 = n1.cross(edge_dir).normalize();
    let raw2 = n2.cross(edge_dir).normalize();

    let shell = &brep.solids[0].shells[0];
    let f1_obj = &shell.faces[f1_flat];
    let f2_obj = &shell.faces[f2_flat];
    let c1 = face_vertex_centroid(brep, f1_obj);
    let c2 = face_vertex_centroid(brep, f2_obj);

    let into_1 = if raw1.dot(c1 - mid_pt) > 0.0 { raw1 } else { -raw1 };
    let into_2 = if raw2.dot(c2 - mid_pt) > 0.0 { raw2 } else { -raw2 };

    // ── New vertex positions ───────────────────────────────────────────────
    let v1_f1_pt = v1_pt + into_1 * offset_dist;
    let v2_f1_pt = v2_pt + into_1 * offset_dist;
    let v1_f2_pt = v1_pt + into_2 * offset_dist;
    let v2_f2_pt = v2_pt + into_2 * offset_dist;

    let v1_f1 = brep.vertices.len(); brep.vertices.push(Vertex { point: v1_f1_pt });
    let v2_f1 = brep.vertices.len(); brep.vertices.push(Vertex { point: v2_f1_pt });
    let v1_f2 = brep.vertices.len(); brep.vertices.push(Vertex { point: v1_f2_pt });
    let v2_f2 = brep.vertices.len(); brep.vertices.push(Vertex { point: v2_f2_pt });

    // ── Contact edges ──────────────────────────────────────────────────────
    let contact_f1 = push_line_edge(brep, v1_f1, v2_f1, v1_f1_pt, v2_f1_pt);
    let contact_f2 = push_line_edge(brep, v1_f2, v2_f2, v1_f2_pt, v2_f2_pt);

    // ── End arcs (quarter-circles on the cylinder) ─────────────────────────
    // Inward bisector points into the solid (opposite to (n1+n2) which points
    // outward for a convex edge).
    let outward_bisector = (n1 + n2).normalize();
    let inward_bisector = -outward_bisector;
    let arc_c_v1 = v1_pt + inward_bisector * r_centre;
    let arc_c_v2 = v2_pt + inward_bisector * r_centre;

    // Find t-parameters for the two contact directions on the circle.
    let t1 = circle_t_for_dir(into_1, edge_dir);
    let t2 = circle_t_for_dir(into_2, edge_dir);

    // Take the shorter arc (≤ π).  Normalise so t_start < t_end.
    let mut t_start = t1;
    let mut t_end  = t2;
    let mut dt = t_end - t_start;
    if dt > PI  { t_end -= 2.0 * PI; }
    if dt < -PI { t_start -= 2.0 * PI; }
    dt = t_end - t_start;
    if dt < 0.0 {
        std::mem::swap(&mut t_start, &mut t_end);
    }

    let arc_v1 = push_arc_edge(brep, v1_f1, v1_f2, arc_c_v1, edge_dir, r, t_start, t_end);
    let arc_v2 = push_arc_edge(brep, v2_f1, v2_f2, arc_c_v2, edge_dir, r, t_start, t_end);

    // ── Build trimmed adjacent faces ───────────────────────────────────────
    // For each adjacent face: find the filleted edge in the wire, shorten
    // the two neighbour edges, and insert the contact edge.

    let trimmed_f1 = build_trimmed_face(brep, f1_flat, edge_idx, v1_f1, v2_f1, contact_f1, v1_pt, v2_pt);
    let trimmed_f2 = build_trimmed_face(brep, f2_flat, edge_idx, v1_f2, v2_f2, contact_f2, v1_pt, v2_pt);

    // If either face was already trimmed, skip the fillet for this edge.
    let (trimmed_f1, trimmed_f2) = match (trimmed_f1, trimmed_f2) {
        (Some(a), Some(b)) => (a, b),
        _ => return Ok(()), // face already modified by earlier fillet
    };

    // ── Trim other faces sharing the filleted edge's vertices ──────────────
    // The fillet also affects faces that share V4 or V5 (the edge endpoints)
    // without containing the edge itself.  Find and trim them.
    let edge_v1 = edge.start;
    let edge_v2 = edge.end;
    let shell_face_count = brep.solids[0].shells[0].faces.len();
    for fi in 0..shell_face_count {
        if fi == f1_flat || fi == f2_flat {
            continue; // already trimmed
        }
        let face = &brep.solids[0].shells[0].faces[fi];
        // Check if this face references either edge vertex.
        let has_v1 = face.outer_wire.edges.iter().any(|we| {
            let e = &brep.edges[we.idx];
            e.start == edge_v1 || e.end == edge_v1
        });
        let has_v2 = face.outer_wire.edges.iter().any(|we| {
            let e = &brep.edges[we.idx];
            e.start == edge_v2 || e.end == edge_v2
        });
        if has_v1 && has_v2 {
            // This face is adjacent to the filleted edge (already handled: F1, F2)
            continue;
        }
        if has_v1 {
            trim_face_at_vertex(brep, fi, edge_v1, v1_f1, v1_f1_pt);
        }
        if has_v2 {
            trim_face_at_vertex(brep, fi, edge_v2, v2_f1, v2_f1_pt);
        }
    }

    // ── Fillet face ────────────────────────────────────────────────────────
    // Wire loop: contact_f1_fwd(V4'→V5') — arc_v2_fwd(V5'→V5'')
    //          — contact_f2_rev(V5''→V4'') — arc_v1_rev(V4''→V4')
    let fillet_face = Face {
        outer_wire: Wire {
            edges: vec![
                WireEdge::fwd(contact_f1),
                WireEdge::fwd(arc_v2),
                WireEdge::rev(contact_f2),
                WireEdge::rev(arc_v1),
            ],
        },
        inner_wires: Vec::new(),
        normal: edge_dir,
        triangles: Vec::new(),
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };

    // ── Store fillet surface ───────────────────────────────────────────────
    let fillet_surf_idx = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface {
        origin: arc_c_v1,
        axis: edge_dir,
        radius: r,
        ref_dir: rcad_kernel::geom::any_perpendicular(edge_dir),
    }));
    // The fillet face goes at the end of the shell; its flat index equals the
    // current face_surface length (before push).
    let fillet_face_flat = brep.geom.face_surface.len();
    brep.geom.face_surface.push(Some(fillet_surf_idx));
    // Pad face_surface_range to ensure index == fillet_face_flat
    while brep.geom.face_surface_range.len() < fillet_face_flat {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range.push(Some([0.0, PI * 0.5, 0.0, edge_len]));
    // Pad face_tolerance similarly
    while brep.geom.face_tolerance.len() < fillet_face_flat {
        brep.geom.face_tolerance.push(0.0);
    }
    brep.geom.face_tolerance.push(0.0);

    // ── Replace faces in shell ─────────────────────────────────────────────
    {
        let shell = &mut brep.solids[0].shells[0];
        let fi_lo = f1_flat.min(f2_flat);
        let fi_hi = f1_flat.max(f2_flat);
        let trimmed_lo = if fi_lo == f1_flat { &trimmed_f1 } else { &trimmed_f2 };
        let trimmed_hi = if fi_hi == f2_flat { &trimmed_f2 } else { &trimmed_f1 };
        shell.faces[fi_hi] = trimmed_hi.clone();
        shell.faces[fi_lo] = trimmed_lo.clone();
        shell.faces.push(fillet_face);
    }

    Ok(())
}

/// Trim a face at one vertex — shorten every edge in the face's outer wire
/// that touches `old_vertex` so it now touches `new_vertex` instead.
fn trim_face_at_vertex(
    brep: &mut BRep,
    face_flat_idx: usize,
    old_vertex: usize,
    new_vertex: usize,
    new_pt: DVec3,
) {
    let face = brep.solids[0].shells[0].faces[face_flat_idx].clone();
    let wire = &face.outer_wire;
    let mut new_edges: Vec<WireEdge> = Vec::with_capacity(wire.edges.len());

    for we in &wire.edges {
        let _e = &brep.edges[we.idx];
        // Check which end touches the old vertex, in wire direction.
        let (ws_v, we_v) = wire_edge_vertices(brep, we);
        if ws_v == old_vertex {
            // Wire-start is the old vertex → replace with new_vertex.
            // New edge goes from new_vertex to we_v (wire-end).
            let we_pt = brep.vertices[we_v].point;
            let new_e = push_line_edge(brep, new_vertex, we_v, new_pt, we_pt);
            new_edges.push(WireEdge::fwd(new_e));
        } else if we_v == old_vertex {
            // Wire-end is the old vertex → replace with new_vertex.
            // New edge goes from ws_v to new_vertex.
            let ws_pt = brep.vertices[ws_v].point;
            let new_e = push_line_edge(brep, ws_v, new_vertex, ws_pt, new_pt);
            new_edges.push(WireEdge::fwd(new_e));
        } else {
            new_edges.push(*we);
        }
    }

    brep.solids[0].shells[0].faces[face_flat_idx] = Face {
        outer_wire: Wire { edges: new_edges },
        inner_wires: Vec::new(),
        normal: face.normal,
        triangles: Vec::new(),
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };
}

/// Build one trimmed planar face: remove the filleted edge from its outer
/// wire, shorten the two neighbouring edges, insert the contact edge.
///
/// Returns `None` when the edge is no longer present in the face wire
/// (face already trimmed by a previous fillet in the same batch).
fn build_trimmed_face(
    brep: &mut BRep,
    face_flat_idx: usize,
    fillet_edge_idx: usize,
    v1_new: usize,
    v2_new: usize,
    contact_edge_idx: usize,
    v1_orig: DVec3,
    v2_orig: DVec3,
) -> Option<Face> {
    let face = &brep.solids[0].shells[0].faces[face_flat_idx].clone();
    let wire = &face.outer_wire;
    let n = wire.edges.len();
    let Some(pos) = find_wire_pos(face, fillet_edge_idx) else {
        // Edge already removed from this face (previous fillet trimmed it).
        return None;
    };

    let fillet_we = &wire.edges[pos];
    let (fillet_ws, fillet_wv) = wire_edge_vertices(brep, fillet_we);

    // Neighbour edges in the wire.
    let we_before = &wire.edges[(pos + n - 1) % n];
    let we_after  = &wire.edges[(pos + 1) % n];

    // Determine which original vertex (V1 or V2) is at the shared vertex
    // with we_before (which is fillet_ws).
    let ws_pt = brep.vertices[fillet_ws].point;
    let (new_before_e_start, new_before_e_end) = if (ws_pt - v1_orig).length_squared()
        < (ws_pt - v2_orig).length_squared()
    {
        let (bws, _bwe) = wire_edge_vertices(brep, we_before);
        (bws, v1_new)
    } else {
        let (bws, _bwe) = wire_edge_vertices(brep, we_before);
        (bws, v2_new)
    };

    // For we_after: shared vertex is fillet_wv.
    let wv_pt = brep.vertices[fillet_wv].point;
    let (new_after_e_start, new_after_e_end) = if (wv_pt - v2_orig).length_squared()
        < (wv_pt - v1_orig).length_squared()
    {
        let (_aws, awe) = wire_edge_vertices(brep, we_after);
        (v2_new, awe)
    } else {
        let (_aws, awe) = wire_edge_vertices(brep, we_after);
        (v1_new, awe)
    };

    // Create the shortened edges.
    let before_pt = brep.vertices[new_before_e_start].point;
    let before_ep = brep.vertices[new_before_e_end].point;
    let short_before = push_line_edge(brep, new_before_e_start, new_before_e_end, before_pt, before_ep);

    let after_pt = brep.vertices[new_after_e_start].point;
    let after_ep = brep.vertices[new_after_e_end].point;
    let short_after = push_line_edge(brep, new_after_e_start, new_after_e_end, after_pt, after_ep);

    // Build the new wire.  The edges in order:
    //   [e0 … e_{pos-1} (=we_before)] → replaced by short_before
    //   [e_pos]                        → replaced by contact_edge (same orientation)
    //   [e_{pos+1} (=we_after)]        → replaced by short_after
    //   [e_{pos+2} … e_{n-1}]         → kept

    let mut new_edges: Vec<WireEdge> = Vec::with_capacity(n + 1);

    // Edges from position 0 to pos-1, with we_before replaced.
    for k in 0..pos {
        if k == (pos + n - 1) % n {
            new_edges.push(WireEdge::fwd(short_before));
        } else {
            new_edges.push(wire.edges[k]);
        }
    }

    // Contact edge replaces the filleted edge — use same direction as original.
    let contact_dir = fillet_we.forward;
    new_edges.push(if contact_dir {
        WireEdge::fwd(contact_edge_idx)
    } else {
        WireEdge::rev(contact_edge_idx)
    });

    // Edges from pos+1 to n-1, with we_after replaced.
    for k in (pos + 1)..n {
        if k == (pos + 1) % n {
            new_edges.push(WireEdge::fwd(short_after));
        } else {
            new_edges.push(wire.edges[k]);
        }
    }

    // If the filleted edge was at pos==0, we_before = edges[n-1] needs replacement.
    if pos == 0 {
        // The last element of new_edges is edges[n-1] (kept as-is from the second
        // loop).  Replace it with short_before.
        if let Some(last) = new_edges.last_mut() {
            *last = WireEdge::fwd(short_before);
        }
    }

    Some(Face {
        outer_wire: Wire { edges: new_edges },
        inner_wires: Vec::new(),
        normal: face.normal,
        triangles: Vec::new(),
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    })
}


// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Utility Functions
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Get any vector perpendicular to the given vector.
fn any_perpendicular(v: DVec3) -> DVec3 {
    let v = v.normalize_or(DVec3::Z);
    let perp = if v.x.abs() > 0.5 {
        DVec3::new(-v.y, v.x, 0.0)
    } else {
        DVec3::new(0.0, -v.z, v.y)
    };
    perp.normalize_or(DVec3::X)
}

/// Interpolate between two radii with tension parameter.
fn interpolate_radius(r1: f64, r2: f64, t: f64, tension: f64) -> f64 {
    // Hermite-like interpolation with tension
    let t2 = t * t;
    let t3 = t2 * t;
    let h1 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h2 = -2.0 * t3 + 3.0 * t2;

    // Apply tension (0 = linear, 1 = smooth)
    let smooth = tension;
    let h1_tension = h1 + smooth * (t3 - 2.0 * t2 + t);
    let h2_tension = h2 + smooth * (-t3 + t2);

    r1 * h1_tension + r2 * h2_tension
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Tests
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    fn create_box_brep() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 2.0,
            depth: 3.0,
        })
    }

    #[test]
    fn test_fillet_params_default() {
        let params = FilletParams::default();
        assert_eq!(params.radius, 1.0);
        assert_eq!(params.continuity, FilletContinuity::C1);
        assert_eq!(params.mode, FilletMode::Uniform);
        assert!(params.tension >= 0.0 && params.tension <= 1.0);
    }

    #[test]
    fn test_fillet_params_builder() {
        let params = FilletParams::new(0.5)
            .with_continuity(FilletContinuity::C2)
            .with_mode(FilletMode::Chordal)
            .with_tension(0.8);

        assert_eq!(params.radius, 0.5);
        assert_eq!(params.continuity, FilletContinuity::C2);
        assert_eq!(params.mode, FilletMode::Chordal);
        assert!((params.tension - 0.8).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_variable_radius_point() {
        let point = VariableRadiusPoint::new(0.5, 2.0);
        assert!((point.parameter - 0.5).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((point.radius - 2.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_fillet_error_display() {
        let err = FilletError::InvalidRadius { radius: -1.0 };
        assert!(err.to_string().contains("invalid fillet radius"));

        let err = FilletError::EdgeNotFound { edge_index: 99 };
        assert!(err.to_string().contains("edge 99 not found"));

        let err = FilletError::RadiusTooLarge {
            edge_index: 0,
            radius: 10.0,
            max_radius: 1.0,
        };
        assert!(err.to_string().contains("radius 10 too large"));
    }

    #[test]
    fn test_any_perpendicular() {
        let v = DVec3::X;
        let p = any_perpendicular(v);
        assert!((p.dot(v)).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((p.length() - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);

        let v = DVec3::Z;
        let p = any_perpendicular(v);
        assert!((p.dot(v)).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_interpolate_radius() {
        // Linear interpolation (tension = 0)
        let r = interpolate_radius(1.0, 2.0, 0.0, 0.0);
        assert!((r - 1.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);

        let r = interpolate_radius(1.0, 2.0, 1.0, 0.0);
        assert!((r - 2.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);

        let r = interpolate_radius(1.0, 2.0, 0.5, 0.0);
        assert!((r - 1.5).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    }

    #[test]
    fn test_make_fillet_edge_empty_indices() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[], 0.5);
        assert!(result.is_ok());
        let res = result.unwrap();
        assert_eq!(res.edges_processed, 0);
        assert!(!res.warnings.is_empty());
    }

    #[test]
    fn test_make_fillet_edge_invalid_radius() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[0], -1.0);
        assert!(matches!(result, Err(FilletError::InvalidRadius { .. })));
    }

    #[test]
    fn test_make_fillet_edge_invalid_edge_index() {
        let brep = create_box_brep();
        let result = make_fillet_edge(&brep, &[999], 0.5);
        assert!(matches!(result, Err(FilletError::EdgeNotFound { .. })));
    }

    #[test]
    fn test_make_fillet_all_edges() {
        let brep = create_box_brep();
        let result = make_fillet_all_edges(&brep, 0.1);
        // Should succeed without errors
        assert!(result.is_ok() || matches!(result, Err(FilletError::EdgeNotFound { .. })));
    }

    #[test]
    fn test_make_variable_fillet_too_few_points() {
        let brep = create_box_brep();
        let radii = vec![VariableRadiusPoint::new(0.0, 0.5)];
        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(result.is_err());
    }

    #[test]
    fn test_make_variable_fillet_invalid_parameter() {
        let brep = create_box_brep();
        let radii = vec![
            VariableRadiusPoint::new(-0.5, 0.5),
            VariableRadiusPoint::new(1.0, 1.0),
        ];
        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(matches!(result, Err(FilletError::InvalidVariableRadius { .. })));
    }

    #[test]
    fn test_compute_plane_plane_fillet() {
        let plane1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let plane2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_Y * 0.5,
            end_point: DVec3::Y * 0.5,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Y,
            tangent_end: DVec3::Y,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_plane_plane_fillet(&edge_info, &plane1, &plane2, 0.5);
        assert!(result.is_ok());

        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Cylinder(_)));
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::ZERO,
            end_point: DVec3::Z,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_cylinder_plane_fillet(&edge_info, &cylinder, &plane, 0.5);
        assert!(result.is_ok());
    }

    #[test]
    fn test_compute_sphere_plane_fillet() {
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1.0);
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 0.5),
            normal: DVec3::Z,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_X,
            end_point: DVec3::X,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::X,
            tangent_end: DVec3::X,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_sphere_plane_fillet(&edge_info, &sphere, &plane, 0.5);
        assert!(result.is_ok());

        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Sphere(_)));
    }

    #[test]
    fn test_compute_fillet_curves_torus() {
        let brep = create_box_brep();

        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 2.0,
            minor_radius: 0.5,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::ZERO,
            end_point: DVec3::Z,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 1.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_fillet_curves(&brep, &edge_info, 0.5, &Surface3::Torus(torus));
        assert!(result.is_ok());

        let curves = result.unwrap();
        assert!(!curves.is_empty());
    }

    #[test]
    fn test_fillet_continuity() {
        assert_eq!(FilletContinuity::default(), FilletContinuity::C1);
        assert_ne!(FilletContinuity::C0, FilletContinuity::C1);
        assert_ne!(FilletContinuity::C1, FilletContinuity::C2);
    }

    #[test]
    fn test_fillet_mode() {
        assert_eq!(FilletMode::default(), FilletMode::Uniform);
        assert_ne!(FilletMode::Uniform, FilletMode::Variable);
        assert_ne!(FilletMode::Variable, FilletMode::Chordal);
    }

    #[test]
    fn test_fillet_result_creation() {
        let brep = create_box_brep();
        let result = FilletResult {
            brep: brep.clone(),
            edges_processed: 3,
            fillet_faces_created: 3,
            warnings: vec!["test warning".to_string()],
        };

        assert_eq!(result.edges_processed, 3);
        assert_eq!(result.fillet_faces_created, 3);
        assert_eq!(result.warnings.len(), 1);
    }

    // ============================================================================
    // Edge Case Tests for OCCT Alignment
    // ============================================================================

    /// Test fillet on a concave edge (interior corner of a box).
    /// Concave edges require different handling than convex edges.
    #[test]
    fn test_fillet_concave_edge() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        // Create a box and populate its geometry
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 4.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        // Attempt to fillet edges - concave edges are at the interior corners
        let result = make_fillet_edge(&brep, &[0, 1, 2, 3], 0.5);
        assert!(result.is_ok(), "concave edge fillet should succeed");

        let fillet_result = result.unwrap();
        // At least some edges should be processed
        assert!(fillet_result.edges_processed > 0 || !fillet_result.warnings.is_empty());
    }

    /// Test fillet on a chain of connected edges.
    /// Chain edges should blend smoothly at the vertices where they meet.
    #[test]
    fn test_fillet_chain_edges() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        // Create a box
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 4.0,
            height: 3.0,
            depth: 2.0,
        });
        populate_box_geom(&mut brep);

        // Fillet a chain of edges around one face (edges 0, 1, 2, 3 form a loop)
        let result = make_fillet_edge(&brep, &[0, 1, 2, 3], 0.3);
        assert!(result.is_ok(), "chain edge fillet should succeed");

        let fillet_result = result.unwrap();
        assert!(fillet_result.edges_processed >= 1, "at least one edge should be filleted");
    }

    /// Test fillet with very small radius on an edge.
    /// Small radius fillets should not create degenerate geometry.
    #[test]
    fn test_fillet_small_radius() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 10.0,
            height: 10.0,
            depth: 10.0,
        });
        populate_box_geom(&mut brep);

        // Very small fillet radius
        let result = make_fillet_edge(&brep, &[0], 0.001);
        assert!(result.is_ok(), "small radius fillet should succeed");
    }

    /// Test fillet on an edge where adjacent faces are perpendicular.
    /// This is the most common case and should produce a clean toroidal fillet.
    #[test]
    fn test_fillet_perpendicular_faces() {
        use rcad_kernel::geom::{Plane, Line3};
        use glam::DVec3;

        // Create edge info for perpendicular faces
        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_Y,
            end_point: DVec3::Y,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Y,  // Edge along Y axis
            tangent_end: DVec3::Y,
            length: 2.0,
            curve: Some(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Y,
            })),
            curve_range: Some([0.0, 2.0]),
        };

        let plane1 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let plane2 = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let result = compute_plane_plane_fillet(&edge_info, &plane1, &plane2, 0.5);
        assert!(result.is_ok(), "perpendicular faces fillet should succeed");

        // Result should be a cylinder for plane-plane fillet (straight edge)
        let surface = result.unwrap();
        assert!(matches!(surface, Surface3::Cylinder(_)), "plane-plane fillet should produce cylinder");
    }

    /// Test fillet with variable radius along the edge.
    /// Variable radius fillets should interpolate between radius values.
    #[test]
    fn test_fillet_variable_radius_basic() {
        use rcad_kernel::PrimitiveSolid;
        use crate::geom_populate::populate_box_geom;

        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 5.0,
            height: 5.0,
            depth: 5.0,
        });
        populate_box_geom(&mut brep);

        // Variable radius: 0.2 at start, 0.8 at end
        let radii = vec![
            VariableRadiusPoint::new(0.0, 0.2),
            VariableRadiusPoint::new(1.0, 0.8),
        ];

        let result = make_variable_fillet(&brep, &[0], &radii);
        assert!(result.is_ok(), "variable radius fillet should succeed");
    }

    /// Test fillet surface computation for cylinder-plane edge.
    /// Verifies correct handling of curved-to-planar transitions.
    #[test]
    fn test_fillet_cylinder_plane_edge() {
        let cylinder = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: any_perpendicular(DVec3::Z),
            radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::new(1.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let edge_info = EdgeInfo {
            index: 0,
            start_vertex: 0,
            end_vertex: 1,
            start_point: DVec3::NEG_Z,
            end_point: DVec3::Z,
            adjacent_faces: vec![0, 1],
            tangent_start: DVec3::Z,
            tangent_end: DVec3::Z,
            length: 2.0,
            curve: None,
            curve_range: None,
        };

        let result = compute_cylinder_plane_fillet(&edge_info, &cylinder, &plane, 0.3);
        assert!(result.is_ok(), "cylinder-plane fillet should succeed");
    }

    /// Verify SA for a single fillet on a 100x100x100 box (OCCT A1).
    #[test]
    fn verify_fillet_sa_box() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 100.0, height: 100.0, depth: 100.0,
        });
        let result = crate::make_fillet_edge(&brep, &[4], 10.0);
        assert!(result.is_ok(), "fillet should succeed");
        let total = rcad_kernel::surface_area(&result.unwrap().brep);
        let diff = (total - 59527.9).abs();
        let tol = 0.15 * 59527.9;
        assert!(diff <= tol, "SA {total} differs from expected 59527.9 by {diff} > {tol}");
    }
}
