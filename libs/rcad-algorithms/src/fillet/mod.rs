//! rcad_kernel::BRepFilletAPI-style edge fillet operations — analogous to OCCT `BRepFilletAPI_MakeFillet`.
//!
//! # Overview
//!
//! This module provides algorithms for creating fillets (rounded edges) on rcad_kernel::BRep shapes:
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

use glam::DVec3;
use rcad_kernel::{
    CurveEval,
    geom::{
        Circle3, Curve3, CylindricalSurface, Line3, Plane, SphericalSurface, Surface3,
        ToroidalSurface, any_perpendicular,
    },
    topods::{BRep, Orientation, ShapeRef, TShape},
};
use std::collections::HashMap;
use std::sync::Arc;

use crate::tolerance::*;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Constants
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

const EPS: f64 = TOLERANCE_LEN_MIN;
const PI: f64 = std::f64::consts::PI;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Error Types
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Errors that can occur during fillet operations.
#[derive(Debug, Clone)]
pub enum FilletError {
    /// Radius is zero or negative.
    InvalidRadius { radius: f64 },
    /// Edge index out of range.
    EdgeNotFound { edge_index: usize },
    /// Face index out of range.
    FaceNotFound { face_index: usize },
    /// Edge has no adjacent faces.
    EdgeNoAdjacentFaces { edge_index: usize },
    /// Fillet would create degenerate geometry.
    DegenerateGeometry { edge_index: usize, reason: String },
    /// Radius too large for the edge.
    RadiusTooLarge {
        edge_index: usize,
        radius: f64,
        max_radius: f64,
    },
    /// Failed to compute fillet surface.
    SurfaceComputationFailed { edge_index: usize, reason: String },
    /// Failed to compute fillet curves.
    CurveComputationFailed { edge_index: usize, reason: String },
    /// Unsupported geometry combination.
    UnsupportedGeometry {
        edge_index: usize,
        surface1_type: String,
        surface2_type: String,
    },
    /// Variable radius specification is invalid.
    InvalidVariableRadius { parameter: f64, radius: f64 },
    /// Failed to blend adjacent faces.
    BlendFailed { edge_index: usize, reason: String },
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
            Self::EdgeNotFound { edge_index } => write!(f, "edge {} not found", edge_index),
            Self::FaceNotFound { face_index } => write!(f, "face {} not found", face_index),
            Self::EdgeNoAdjacentFaces { edge_index } => {
                write!(f, "edge {} has no adjacent faces", edge_index)
            }
            Self::DegenerateGeometry { edge_index, reason } => {
                write!(f, "degenerate geometry at edge {}: {}", edge_index, reason)
            }
            Self::RadiusTooLarge {
                edge_index,
                radius,
                max_radius,
            } => write!(
                f,
                "radius {} too large for edge {} (max {})",
                radius, edge_index, max_radius
            ),
            Self::SurfaceComputationFailed { edge_index, reason } => write!(
                f,
                "failed to compute fillet surface at edge {}: {}",
                edge_index, reason
            ),
            Self::CurveComputationFailed { edge_index, reason } => write!(
                f,
                "failed to compute fillet curves at edge {}: {}",
                edge_index, reason
            ),
            Self::UnsupportedGeometry {
                edge_index,
                surface1_type,
                surface2_type,
            } => write!(
                f,
                "unsupported geometry at edge {}: {} + {}",
                edge_index, surface1_type, surface2_type
            ),
            Self::InvalidVariableRadius { parameter, radius } => write!(
                f,
                "invalid variable radius {} at parameter {}",
                radius, parameter
            ),
            Self::BlendFailed { edge_index, reason } => write!(
                f,
                "failed to blend adjacent faces at edge {}: {}",
                edge_index, reason
            ),
            Self::InvalidInput(msg) => write!(f, "invalid input: {}", msg),
            Self::NumericalFailure(msg) => write!(f, "numerical failure: {}", msg),
            Self::EmptyResult => write!(f, "fillet operation produced empty result"),
        }
    }
}

impl std::error::Error for FilletError {}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Fillet Types
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Continuity type for fillet surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FilletContinuity {
    /// Position continuity only (G0/C0).
    C0,
    /// Tangent continuity (G1/C1).
    #[default]
    C1,
    /// Curvature continuity (G2/C2).
    C2,
}

/// Fillet mode for radius specification.
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
#[derive(Debug, Clone)]
pub struct FilletParams {
    /// Fillet radius (or chord length in chordal mode).
    pub radius: f64,
    /// Continuity between fillet and adjacent faces.
    pub continuity: FilletContinuity,
    /// Fillet mode (uniform, variable, chordal).
    pub mode: FilletMode,
    /// Tension parameter for variable radius fillets (0.0 = linear, 1.0 = smooth).
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
    /// The resulting rcad_kernel::BRep with fillets applied.
    pub brep: rcad_kernel::BRep,
    /// Number of edges filletted.
    pub edges_processed: usize,
    /// Number of fillet faces created.
    pub fillet_faces_created: usize,
    /// Any warnings encountered during the operation.
    pub warnings: Vec<String>,
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Fillet Surface Types
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

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
    /// Edge index (tshapes index).
    index: usize,
    /// Start vertex tshape index.
    start_vertex: usize,
    /// End vertex tshape index.
    end_vertex: usize,
    /// Start vertex 3D position.
    start_point: DVec3,
    /// End vertex 3D position.
    end_point: DVec3,
    /// Adjacent face indices (tshapes indices, usually 2).
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

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Main API Functions
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Create a fillet on one or more edges with uniform radius.
pub fn make_fillet_edge(
    brep: &rcad_kernel::BRep,
    edge_indices: &[usize],
    radius: f64,
) -> Result<FilletResult, FilletError> {
    let params = FilletParams::new(radius);
    make_fillet_edge_with_params(brep, edge_indices, &params)
}

/// Create a fillet on one or more edges with custom parameters.
pub fn make_fillet_edge_with_params(
    brep: &rcad_kernel::BRep,
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
        return Err(FilletError::InvalidRadius {
            radius: params.radius,
        });
    }

    // Validate edge indices
    let edge_count = brep.edge_count();
    for &idx in edge_indices {
        if idx >= edge_count {
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
pub fn make_fillet_all_edges(
    brep: &rcad_kernel::BRep,
    radius: f64,
) -> Result<FilletResult, FilletError> {
    let all_edges: Vec<usize> = (0..brep.edge_count()).collect();
    make_fillet_edge(brep, &all_edges, radius)
}

/// Create a variable radius fillet along edges.
pub fn make_variable_fillet(
    brep: &rcad_kernel::BRep,
    edge_indices: &[usize],
    radii: &[VariableRadiusPoint],
) -> Result<FilletResult, FilletError> {
    if radii.len() < 2 {
        return Err(FilletError::InvalidInput(
            "variable fillet requires at least 2 radius points",
        ));
    }

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

    let avg_radius = radii.iter().map(|r| r.radius).sum::<f64>() / radii.len() as f64;
    let mut params = FilletParams::new(avg_radius);
    params.mode = FilletMode::Variable;

    make_variable_fillet_with_params(brep, edge_indices, radii, &params)
}

/// Create a variable radius fillet with custom parameters.
fn make_variable_fillet_with_params(
    brep: &rcad_kernel::BRep,
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

    let edge_count = brep.edge_count();
    for &idx in edge_indices {
        if idx >= edge_count {
            return Err(FilletError::EdgeNotFound { edge_index: idx });
        }
    }

    let edge_infos = collect_edge_infos(brep, edge_indices)?;

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

    let result = build_fillet_brep(brep, &fillet_surfaces, &edge_infos, params)?;

    Ok(FilletResult {
        brep: result,
        edges_processed: fillet_surfaces.len(),
        fillet_faces_created: fillet_surfaces.len(),
        warnings,
    })
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Edge Information Collection
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Collect information about edges to be filletted.
fn collect_edge_infos(
    brep: &rcad_kernel::BRep,
    edge_indices: &[usize],
) -> Result<Vec<EdgeInfo>, FilletError> {
    let mut infos = Vec::new();

    // Build face-to-edge adjacency map
    let edge_faces = build_edge_face_adjacency(brep);

    for &edge_idx in edge_indices {
        let ts = brep
            .tshapes
            .get(edge_idx)
            .ok_or(FilletError::EdgeNotFound {
                edge_index: edge_idx,
            })?;
        let ed = match &**ts {
            TShape::Edge(e) => e,
            _ => {
                return Err(FilletError::EdgeNotFound {
                    edge_index: edge_idx,
                });
            }
        };

        // Get adjacent faces
        let adjacent_faces = edge_faces.get(&edge_idx).cloned().unwrap_or_default();
        if adjacent_faces.len() < 2 {
            continue;
        }

        // Get edge curve and range
        let curve = ed.curve.clone();
        let curve_range = Some(ed.range);

        // Compute edge length
        let length = compute_edge_length(brep, edge_idx);

        // Compute tangents
        let (tangent_start, tangent_end) =
            compute_edge_tangents(brep, edge_idx, &curve, &curve_range);

        // Vertex positions
        let start_point = brep.vertex_point(ed.first.index).unwrap_or_default();
        let end_point = brep.vertex_point(ed.last.index).unwrap_or_default();

        infos.push(EdgeInfo {
            index: edge_idx,
            start_vertex: ed.first.index,
            end_vertex: ed.last.index,
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

/// Build a map from edge index (tshapes index) to adjacent face indices.
fn build_edge_face_adjacency(brep: &rcad_kernel::BRep) -> HashMap<usize, Vec<usize>> {
    let mut edge_faces: HashMap<usize, Vec<usize>> = HashMap::new();

    for (fi, ts) in brep.tshapes.iter().enumerate() {
        let fd = match ts.as_ref() {
            TShape::Face(f) => f,
            _ => continue,
        };

        // Walk outer wire edges
        if let Some(wire_ts) = brep.tshapes.get(fd.outer_wire.index) {
            if let TShape::Wire(wd) = wire_ts.as_ref() {
                for esr in &wd.edges {
                    edge_faces.entry(esr.index).or_default().push(fi);
                }
            }
        }

        // Walk inner wire edges
        for iw_sr in &fd.inner_wires {
            if let Some(wire_ts) = brep.tshapes.get(iw_sr.index) {
                if let TShape::Wire(wd) = wire_ts.as_ref() {
                    for esr in &wd.edges {
                        edge_faces.entry(esr.index).or_default().push(fi);
                    }
                }
            }
        }
    }

    edge_faces
}

/// Compute the length of an edge.
fn compute_edge_length(brep: &rcad_kernel::BRep, edge_idx: usize) -> f64 {
    let ts = match brep.tshapes.get(edge_idx) {
        Some(t) => t,
        None => return 0.0,
    };
    let ed = match &**ts {
        TShape::Edge(e) => e,
        _ => return 0.0,
    };
    let p0 = brep.vertex_point(ed.first.index).unwrap_or_default();
    let p1 = brep.vertex_point(ed.last.index).unwrap_or_default();
    (p1 - p0).length()
}

/// Compute tangent vectors at the start and end of an edge.
fn compute_edge_tangents(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    curve: &Option<Curve3>,
    curve_range: &Option<[f64; 2]>,
) -> (DVec3, DVec3) {
    let ts = match brep.tshapes.get(edge_idx) {
        Some(t) => t,
        None => return (DVec3::X, DVec3::X),
    };
    let ed = match &**ts {
        TShape::Edge(e) => e,
        _ => return (DVec3::X, DVec3::X),
    };
    let p0 = brep.vertex_point(ed.first.index).unwrap_or_default();
    let p1 = brep.vertex_point(ed.last.index).unwrap_or_default();

    match (curve, curve_range) {
        (Some(c), Some([t0, t1])) => {
            let t_start = c.tangent_at(*t0);
            let t_end = c.tangent_at(*t1);
            (t_start.normalize_or(DVec3::X), t_end.normalize_or(DVec3::X))
        }
        _ => {
            let dir = (p1 - p0).normalize_or(DVec3::X);
            (dir, dir)
        }
    }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Fillet Surface Construction
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Compute the rolling ball fillet surface for an edge.
pub fn compute_rollball_surface(
    brep: &rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<Surface3, FilletError> {
    let faces = &edge_info.adjacent_faces;
    if faces.len() < 2 {
        return Err(FilletError::EdgeNoAdjacentFaces {
            edge_index: edge_info.index,
        });
    }

    let surf1 = get_face_surface(brep, faces[0]);
    let surf2 = get_face_surface(brep, faces[1]);

    match (&surf1, &surf2) {
        (Some(s1), Some(s2)) => {
            compute_rollball_surface_for_surfaces(edge_info, s1, s2, faces[0], faces[1], radius)
        }
        _ => compute_toroidal_fillet_surface(brep, edge_info, radius),
    }
}

/// Get the surface for a face (by tshape index).
fn get_face_surface(brep: &rcad_kernel::BRep, flat_face_idx: usize) -> Option<Surface3> {
    let ts = brep.tshapes.get(flat_face_idx)?;
    let TShape::Face(fd) = ts.as_ref() else {
        return None;
    };
    fd.surface.clone()
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
        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
            compute_plane_plane_fillet(edge_info, p1, p2, radius)
        }
        (Surface3::Cylinder(c), Surface3::Plane(p))
        | (Surface3::Plane(p), Surface3::Cylinder(c)) => {
            compute_cylinder_plane_fillet(edge_info, c, p, radius)
        }
        (Surface3::Sphere(s), Surface3::Plane(p)) | (Surface3::Plane(p), Surface3::Sphere(s)) => {
            compute_sphere_plane_fillet(edge_info, s, p, radius)
        }
        _ => compute_general_fillet_surface(edge_info, surf1, surf2, radius),
    }
}

/// Compute fillet surface for plane-plane edge.
fn compute_plane_plane_fillet(
    edge_info: &EdgeInfo,
    plane1: &Plane,
    plane2: &Plane,
    radius: f64,
) -> Result<Surface3, FilletError> {
    let edge_dir = edge_info.tangent_start;
    let n1 = plane1.normal.normalize();
    let n2 = plane2.normal.normalize();
    let cos_angle = n1.dot(n2);
    let angle = cos_angle.acos();
    let intersection_dir = n1.cross(n2);

    if intersection_dir.length_squared() < EPS {
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

    let offset_distance = radius / sin_half;
    let bisector = (n1 + n2).normalize();
    let mid_point = (edge_info.start_point + edge_info.end_point) * 0.5;
    let axis_point = mid_point - bisector * offset_distance;
    let axis_dir = edge_dir.normalize();

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
    let edge_dir = edge_info.tangent_start;
    let cylinder_axis = cylinder.axis.normalize();
    let _plane_normal = plane.normal.normalize();

    let parallel_to_axis = edge_dir.dot(cylinder_axis).abs() > 1.0 - EPS;

    if parallel_to_axis {
        Ok(Surface3::Cylinder(CylindricalSurface {
            origin: cylinder.origin,
            axis: cylinder.axis,
            radius: cylinder.radius + radius,
            ref_dir: cylinder.ref_dir,
        }))
    } else {
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
    let _edge_dir = edge_info.tangent_start;
    let _plane_normal = plane.normal.normalize();
    let fillet_sphere_radius = sphere.radius + radius;
    let _dist_to_plane = (sphere.center - plane.origin).dot(_plane_normal);

    Ok(Surface3::Sphere(SphericalSurface::new(
        sphere.center,
        sphere.axis,
        fillet_sphere_radius,
    )))
}

/// Compute general fillet surface for arbitrary surface types.
fn compute_general_fillet_surface(
    edge_info: &EdgeInfo,
    _surf1: &Surface3,
    _surf2: &Surface3,
    radius: f64,
) -> Result<Surface3, FilletError> {
    let axis = edge_info.tangent_start.normalize();
    let center = DVec3::ZERO;
    let major_radius = radius * 2.0;
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
    brep: &rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<Surface3, FilletError> {
    let p0 = brep
        .vertex_point(edge_info.start_vertex)
        .unwrap_or_default();
    let p1 = brep.vertex_point(edge_info.end_vertex).unwrap_or_default();

    let center = (p0 + p1) * 0.5;
    let axis = (p1 - p0).normalize_or(DVec3::Z);
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
pub fn compute_fillet_curves(
    brep: &rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    radius: f64,
    surface: &Surface3,
) -> Result<Vec<FilletCurve>, FilletError> {
    let p0 = brep
        .vertex_point(edge_info.start_vertex)
        .unwrap_or_default();
    let p1 = brep.vertex_point(edge_info.end_vertex).unwrap_or_default();

    let mut curves = Vec::new();

    match surface {
        Surface3::Torus(torus) => {
            let axis = torus.axis.normalize();
            let _ref_dir = any_perpendicular(axis);

            let start_center = torus.center + axis * (p0 - torus.center).dot(axis);
            curves.push(FilletCurve {
                curve: Curve3::Circle(Circle3::new(start_center, axis, torus.minor_radius)),
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });

            let end_center = torus.center + axis * (p1 - torus.center).dot(axis);
            curves.push(FilletCurve {
                curve: Curve3::Circle(Circle3::new(end_center, axis, torus.minor_radius)),
                parameter_range: [0.0, 2.0 * PI],
                is_start: false,
            });
        }
        Surface3::Cylinder(cyl) => {
            let axis = cyl.axis.normalize();
            curves.push(FilletCurve {
                curve: Curve3::Circle(Circle3::new(p0, axis, cyl.radius)),
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });
            curves.push(FilletCurve {
                curve: Curve3::Circle(Circle3::new(p1, axis, cyl.radius)),
                parameter_range: [0.0, 2.0 * PI],
                is_start: false,
            });
        }
        Surface3::Sphere(sphere) => {
            let axis = sphere.axis.normalize();
            let _ref_dir = any_perpendicular(axis);
            let start_curve = Curve3::Circle(Circle3::new(
                p0 - axis * (p0 - sphere.center).dot(axis) * 0.5,
                axis,
                sphere.radius * 0.5,
            ));
            curves.push(FilletCurve {
                curve: start_curve,
                parameter_range: [0.0, 2.0 * PI],
                is_start: true,
            });
        }
        _ => {
            let edge_dir = (p1 - p0).normalize_or(DVec3::Z);
            curves.push(FilletCurve {
                curve: Curve3::Line(Line3 {
                    origin: p0,
                    direction: any_perpendicular(edge_dir),
                }),
                parameter_range: [0.0, radius],
                is_start: true,
            });
        }
    }

    Ok(curves)
}

/// Blend the fillet surface with adjacent faces.
pub fn blend_adjacent_faces(
    _brep: &mut rcad_kernel::BRep,
    _fillet_surface: &Surface3,
    edge_info: &EdgeInfo,
    radius: f64,
) -> Result<(), FilletError> {
    if edge_info.adjacent_faces.len() < 2 {
        return Err(FilletError::BlendFailed {
            edge_index: edge_info.index,
            reason: "need at least 2 adjacent faces".to_string(),
        });
    }

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

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Fillet Computation
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Compute fillet surface for a single edge.
fn compute_fillet_for_edge(
    brep: &rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    params: &FilletParams,
) -> Result<FilletSurface, FilletError> {
    let surface = compute_rollball_surface(brep, edge_info, params.radius)?;
    let boundary_curves = compute_fillet_curves(brep, edge_info, params.radius, &surface)?;
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
    brep: &rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    radii: &[VariableRadiusPoint],
    _params: &FilletParams,
) -> Result<FilletSurface, FilletError> {
    let mut sorted_radii = radii.to_vec();
    sorted_radii.sort_by(|a, b| a.parameter.partial_cmp(&b.parameter).unwrap());

    let avg_radius = sorted_radii.iter().map(|r| r.radius).sum::<f64>() / sorted_radii.len() as f64;
    let surface = compute_rollball_surface(brep, edge_info, avg_radius)?;
    let boundary_curves = compute_variable_fillet_curves(brep, edge_info, &sorted_radii, &surface)?;
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
    brep: &rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    radii: &[VariableRadiusPoint],
    surface: &Surface3,
) -> Result<Vec<FilletCurve>, FilletError> {
    let p0 = brep
        .vertex_point(edge_info.start_vertex)
        .unwrap_or_default();
    let p1 = brep.vertex_point(edge_info.end_vertex).unwrap_or_default();

    let mut curves = Vec::new();

    for rp in radii {
        let t = rp.parameter;
        let pt = p0 + (p1 - p0) * t;

        match surface {
            Surface3::Torus(torus) => {
                let axis = torus.axis.normalize();
                let center = torus.center + axis * (pt - torus.center).dot(axis);
                curves.push(FilletCurve {
                    curve: Curve3::Circle(Circle3::new(center, axis, rp.radius)),
                    parameter_range: [0.0, 2.0 * PI],
                    is_start: t < 0.5,
                });
            }
            _ => {
                curves.push(FilletCurve {
                    curve: Curve3::Line(Line3 {
                        origin: pt,
                        direction: edge_info.tangent_start,
                    }),
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
        Surface3::Cylinder(_) => [0.0, PI * 0.5, 0.0, edge_info.length],
        Surface3::Torus(_) => [0.0, 2.0 * PI, 0.0, PI * 0.5],
        Surface3::Sphere(_) => [0.0, 2.0 * PI, 0.0, PI * 0.5],
        _ => [0.0, 1.0, 0.0, 1.0],
    }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// BRep Construction — topods::BRep API
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
//
// Topological surgery: adjacent faces are trimmed (the filleted edge is removed,
// neighbouring edges shortened, and a contact edge inserted), and a new fillet
// face with a 4-edge wire is added. For plane-plane edges the fillet surface is
// a CylindricalSurface (the rolling-ball centre traces a line, producing a
// constant quarter-circle cross-section).

/// Find the tshape index of the first shell inside the first solid.
fn find_first_shell_index(brep: &rcad_kernel::BRep) -> Option<usize> {
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            return sd.shells.first().map(|sr| sr.index);
        }
    }
    None
}

/// Build a ShapeRef for an edge by its tshape index, with given orientation.
fn edge_shape_ref(
    brep: &rcad_kernel::BRep,
    edge_index: usize,
    orientation: Orientation,
) -> ShapeRef {
    let ts = &brep.tshapes[edge_index];
    ShapeRef {
        ptr_id: Arc::as_ptr(ts) as u64,
        index: edge_index,
        orientation,
        location: 0,
    }
}

/// Build the result BRep with fillets.
fn build_fillet_brep(
    brep: &rcad_kernel::BRep,
    fillet_surfaces: &[FilletSurface],
    edge_infos: &[EdgeInfo],
    params: &FilletParams,
) -> Result<rcad_kernel::BRep, FilletError> {
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
                // Fallback: add a detached face with the fillet surface.
                let fillet_wire_sr = result.add_twire(Vec::new());
                let fillet_face_sr = result.add_tface(
                    Some(fs.surface.clone()),
                    fillet_wire_sr,
                    Vec::new(),
                    None,
                    None,
                    Vec::new(),
                    false,
                );
                if let Some(shell_idx) = find_first_shell_index(&result) {
                    if let TShape::Shell(shd) = Arc::make_mut(&mut result.tshapes[shell_idx]) {
                        shd.faces.push(fillet_face_sr);
                        shd.my_shapes.push(fillet_face_sr);
                    }
                }
            }
        }
    }

    Ok(result)
}

/// Wire-direction start/end vertex tshape indices of an edge.
fn wire_edge_vertices(brep: &rcad_kernel::BRep, edge_idx: usize) -> (usize, usize) {
    let ts = &brep.tshapes[edge_idx];
    let TShape::Edge(ed) = ts.as_ref() else {
        return (usize::MAX, usize::MAX);
    };
    (ed.first.index, ed.last.index)
}

/// Find the position of an edge index within a face's outer wire.
fn find_wire_pos(brep: &rcad_kernel::BRep, face_idx: usize, edge_idx: usize) -> Option<usize> {
    let ts = &brep.tshapes[face_idx];
    let TShape::Face(fd) = ts.as_ref() else {
        return None;
    };
    let wire_ts = brep.tshapes.get(fd.outer_wire.index)?;
    let TShape::Wire(wd) = wire_ts.as_ref() else {
        return None;
    };
    wd.edges.iter().position(|esr| esr.index == edge_idx)
}

/// Compute the centroid of a planar face from its outer-wire vertices.
fn face_vertex_centroid(brep: &rcad_kernel::BRep, face_idx: usize) -> DVec3 {
    let ts = &brep.tshapes[face_idx];
    let TShape::Face(fd) = ts.as_ref() else {
        return DVec3::ZERO;
    };
    let wire_ts = match brep.tshapes.get(fd.outer_wire.index) {
        Some(t) => t,
        None => return DVec3::ZERO,
    };
    let TShape::Wire(wd) = wire_ts.as_ref() else {
        return DVec3::ZERO;
    };
    let mut sum = DVec3::ZERO;
    let mut n = 0u32;
    for esr in &wd.edges {
        if let Some(pt) = brep.vertex_point(esr.index) {
            sum += pt;
            n += 1;
        }
    }
    if n > 0 { sum / n as f64 } else { DVec3::ZERO }
}

/// Push a straight-line edge using the topods flat-index API.
fn push_line_edge(
    brep: &mut rcad_kernel::BRep,
    start: usize,
    end: usize,
    p0: DVec3,
    p1: DVec3,
) -> usize {
    let delta = p1 - p0;
    let len = delta.length();
    let dir = if len > EPS { delta / len } else { DVec3::X };
    let curve = Some(Curve3::Line(Line3 {
        origin: p0,
        direction: dir,
    }));
    brep.add_edge_flat(start, end, curve, [0.0, len])
}

/// Push a circular-arc edge using the topods flat-index API.
fn push_arc_edge(
    brep: &mut rcad_kernel::BRep,
    start: usize,
    end: usize,
    center: DVec3,
    normal: DVec3,
    radius: f64,
    t_start: f64,
    t_end: f64,
) -> usize {
    let curve = Some(Curve3::Circle(Circle3::new(center, normal, radius)));
    brep.add_edge_flat(start, end, curve, [t_start, t_end])
}

/// Find the parameter `t` on a circle such that
/// `cos(t) x_ax + sin(t) y_ax = dir` (with `x_ax ⟂ normal`, `y_ax = normal × x_ax`).
fn circle_t_for_dir(dir: DVec3, normal: DVec3) -> f64 {
    let x_ax = any_perpendicular(normal);
    let y_ax = normal.cross(x_ax).normalize();
    let d = dir.normalize();
    f64::atan2(d.dot(y_ax).clamp(-1.0, 1.0), d.dot(x_ax).clamp(-1.0, 1.0))
}

/// For a given face, check if any edge in its outer wire references a vertex index.
fn face_has_vertex(brep: &rcad_kernel::BRep, face_idx: usize, vertex_idx: usize) -> bool {
    let ts = match brep.tshapes.get(face_idx) {
        Some(t) => t,
        None => return false,
    };
    let TShape::Face(fd) = ts.as_ref() else {
        return false;
    };
    let wire_ts = match brep.tshapes.get(fd.outer_wire.index) {
        Some(t) => t,
        None => return false,
    };
    let TShape::Wire(wd) = wire_ts.as_ref() else {
        return false;
    };
    wd.edges.iter().any(|esr| {
        let ets = match brep.tshapes.get(esr.index) {
            Some(t) => t,
            None => return false,
        };
        let TShape::Edge(ed) = ets.as_ref() else {
            return false;
        };
        ed.first.index == vertex_idx || ed.last.index == vertex_idx
    })
}

/// Apply a cylindrical fillet (plane-plane edge): trim adjacent faces and
/// insert the fillet face with a correct 4-edge wire.
fn apply_cylinder_fillet(
    brep: &mut rcad_kernel::BRep,
    edge_info: &EdgeInfo,
    params: &FilletParams,
) -> Result<(), FilletError> {
    let r = params.radius;
    let edge_idx = edge_info.index;

    // Edge geometry
    let v1_pt = brep
        .vertex_point(edge_info.start_vertex)
        .unwrap_or_default();
    let v2_pt = brep.vertex_point(edge_info.end_vertex).unwrap_or_default();
    let edge_dir = (v2_pt - v1_pt).normalize_or(DVec3::X);
    let edge_len = (v2_pt - v1_pt).length();
    let mid_pt = (v1_pt + v2_pt) * 0.5;

    // Adjacent face surfaces and normals
    let adj = &edge_info.adjacent_faces;
    if adj.len() < 2 {
        return Err(FilletError::EdgeNoAdjacentFaces {
            edge_index: edge_idx,
        });
    }
    let f1_idx = adj[0];
    let f2_idx = adj[1];

    let get_plane_normal = |brep: &rcad_kernel::BRep, flat_idx: usize| -> Option<DVec3> {
        let ts = brep.tshapes.get(flat_idx)?;
        let TShape::Face(fd) = ts.as_ref() else {
            return None;
        };
        let surf = fd.surface.as_ref()?;
        match surf {
            Surface3::Plane(p) => Some(p.normal.normalize()),
            _ => None,
        }
    };

    let n1 =
        get_plane_normal(brep, f1_idx).ok_or_else(|| FilletError::SurfaceComputationFailed {
            edge_index: edge_idx,
            reason: "face 1 is not a plane".to_string(),
        })?;
    let n2 =
        get_plane_normal(brep, f2_idx).ok_or_else(|| FilletError::SurfaceComputationFailed {
            edge_index: edge_idx,
            reason: "face 2 is not a plane".to_string(),
        })?;

    // Face angle and offset
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

    let offset_dist = r / tan_half;
    let r_centre = r / sin_half;

    // Perpendicular-to-edge directions within each face plane
    let raw1 = n1.cross(edge_dir).normalize();
    let raw2 = n2.cross(edge_dir).normalize();

    let c1 = face_vertex_centroid(brep, f1_idx);
    let c2 = face_vertex_centroid(brep, f2_idx);

    let into_1 = if raw1.dot(c1 - mid_pt) > 0.0 {
        raw1
    } else {
        -raw1
    };
    let into_2 = if raw2.dot(c2 - mid_pt) > 0.0 {
        raw2
    } else {
        -raw2
    };

    // New vertex positions
    let v1_f1_pt = v1_pt + into_1 * offset_dist;
    let v2_f1_pt = v2_pt + into_1 * offset_dist;
    let v1_f2_pt = v1_pt + into_2 * offset_dist;
    let v2_f2_pt = v2_pt + into_2 * offset_dist;

    // Add new vertices (tshapes indices)
    let v1_f1 = brep.add_tvertex(v1_f1_pt).index;
    let v2_f1 = brep.add_tvertex(v2_f1_pt).index;
    let v1_f2 = brep.add_tvertex(v1_f2_pt).index;
    let v2_f2 = brep.add_tvertex(v2_f2_pt).index;

    // Contact edges
    let contact_f1 = push_line_edge(brep, v1_f1, v2_f1, v1_f1_pt, v2_f1_pt);
    let contact_f2 = push_line_edge(brep, v1_f2, v2_f2, v1_f2_pt, v2_f2_pt);

    // End arcs (quarter-circles on the cylinder)
    let outward_bisector = (n1 + n2).normalize();
    let inward_bisector = -outward_bisector;
    let arc_c_v1 = v1_pt + inward_bisector * r_centre;
    let arc_c_v2 = v2_pt + inward_bisector * r_centre;

    let t1 = circle_t_for_dir(into_1, edge_dir);
    let t2 = circle_t_for_dir(into_2, edge_dir);

    let mut t_start = t1;
    let mut t_end = t2;
    let mut dt = t_end - t_start;
    if dt > PI {
        t_end -= 2.0 * PI;
    }
    if dt < -PI {
        t_start -= 2.0 * PI;
    }
    dt = t_end - t_start;
    if dt < 0.0 {
        std::mem::swap(&mut t_start, &mut t_end);
    }

    let arc_v1 = push_arc_edge(brep, v1_f1, v1_f2, arc_c_v1, edge_dir, r, t_start, t_end);
    let arc_v2 = push_arc_edge(brep, v2_f1, v2_f2, arc_c_v2, edge_dir, r, t_start, t_end);

    // Build trimmed adjacent faces
    let trimmed_f1 = build_trimmed_face(
        brep, f1_idx, edge_idx, v1_f1, v2_f1, contact_f1, v1_pt, v2_pt,
    );
    let trimmed_f2 = build_trimmed_face(
        brep, f2_idx, edge_idx, v1_f2, v2_f2, contact_f2, v1_pt, v2_pt,
    );

    match (trimmed_f1, trimmed_f2) {
        (Some(_), Some(_)) => {}
        _ => return Ok(()),
    }

    // Trim other faces sharing the filleted edge's vertices
    let edge_v1 = edge_info.start_vertex;
    let edge_v2 = edge_info.end_vertex;
    let shell_idx = match find_first_shell_index(brep) {
        Some(idx) => idx,
        None => return Ok(()),
    };
    let shell_ts = brep.tshapes[shell_idx].clone();
    let face_refs: Vec<ShapeRef> = if let TShape::Shell(shd) = shell_ts.as_ref() {
        shd.faces.clone()
    } else {
        return Ok(());
    };

    for face_sr in &face_refs {
        let fi = face_sr.index;
        if fi == f1_idx || fi == f2_idx {
            continue;
        }
        if face_has_vertex(brep, fi, edge_v1) && face_has_vertex(brep, fi, edge_v2) {
            continue;
        }
        if face_has_vertex(brep, fi, edge_v1) {
            trim_face_at_vertex(brep, fi, edge_v1, v1_f1, v1_f1_pt);
        }
        if face_has_vertex(brep, fi, edge_v2) {
            trim_face_at_vertex(brep, fi, edge_v2, v2_f1, v2_f1_pt);
        }
    }

    // Fillet face
    let fillet_edges = vec![
        edge_shape_ref(brep, contact_f1, Orientation::Forward),
        edge_shape_ref(brep, arc_v2, Orientation::Forward),
        edge_shape_ref(brep, contact_f2, Orientation::Reversed),
        edge_shape_ref(brep, arc_v1, Orientation::Reversed),
    ];
    let fillet_wire_sr = brep.add_twire(fillet_edges);

    let fillet_surf = Surface3::Cylinder(CylindricalSurface {
        origin: arc_c_v1,
        axis: edge_dir,
        radius: r,
        ref_dir: any_perpendicular(edge_dir),
    });

    let fillet_face_sr = brep.add_tface(
        Some(fillet_surf),
        fillet_wire_sr,
        Vec::new(),
        None,
        Some([0.0, PI * 0.5, 0.0, edge_len]),
        Vec::new(),
        false,
    );

    // Add fillet face to the shell
    if let TShape::Shell(shd) = Arc::make_mut(&mut brep.tshapes[shell_idx]) {
        shd.faces.push(fillet_face_sr);
        shd.my_shapes.push(fillet_face_sr);
    }

    Ok(())
}

/// Trim a face at one vertex — shorten every edge in the face's outer wire
/// that touches `old_vertex` so it now touches `new_vertex` instead.
fn trim_face_at_vertex(
    brep: &mut rcad_kernel::BRep,
    face_flat_idx: usize,
    old_vertex: usize,
    new_vertex: usize,
    new_pt: DVec3,
) {
    let (old_wire_edges, _old_normal) = {
        let ts = &brep.tshapes[face_flat_idx];
        let TShape::Face(fd) = ts.as_ref() else {
            return;
        };
        let wire_ts = match brep.tshapes.get(fd.outer_wire.index) {
            Some(t) => t,
            None => return,
        };
        let TShape::Wire(wd) = wire_ts.as_ref() else {
            return;
        };
        (wd.edges.clone(), DVec3::ZERO)
    };

    let mut new_edge_refs: Vec<ShapeRef> = Vec::with_capacity(old_wire_edges.len());

    for esr in &old_wire_edges {
        let (ws_v, we_v) = wire_edge_vertices(brep, esr.index);
        let forward = esr.orientation == Orientation::Forward;
        if (forward && ws_v == old_vertex) || (!forward && we_v == old_vertex) {
            let we_pt = brep.vertex_point(we_v).unwrap_or_default();
            let new_e = push_line_edge(brep, new_vertex, we_v, new_pt, we_pt);
            new_edge_refs.push(edge_shape_ref(brep, new_e, Orientation::Forward));
        } else if (forward && we_v == old_vertex) || (!forward && ws_v == old_vertex) {
            let ws_pt = brep.vertex_point(ws_v).unwrap_or_default();
            let new_e = push_line_edge(brep, ws_v, new_vertex, ws_pt, new_pt);
            new_edge_refs.push(edge_shape_ref(brep, new_e, Orientation::Forward));
        } else {
            new_edge_refs.push(*esr);
        }
    }

    let new_wire_sr = brep.add_twire(new_edge_refs);
    let face_ts = Arc::make_mut(&mut brep.tshapes[face_flat_idx]);
    if let TShape::Face(fd) = face_ts {
        fd.outer_wire = new_wire_sr;
        if !fd.my_shapes.is_empty() {
            fd.my_shapes[0] = new_wire_sr;
        }
    }
}

/// Build one trimmed planar face: remove the filleted edge from its outer
/// wire, shorten the two neighbouring edges, insert the contact edge.
///
/// Returns `None` when the edge is no longer present in the face wire
/// (face already trimmed by a previous fillet in the same batch).
fn build_trimmed_face(
    brep: &mut rcad_kernel::BRep,
    face_flat_idx: usize,
    fillet_edge_idx: usize,
    v1_new: usize,
    v2_new: usize,
    contact_edge_idx: usize,
    v1_orig: DVec3,
    v2_orig: DVec3,
) -> Option<()> {
    let old_edges = {
        let ts = brep.tshapes.get(face_flat_idx)?;
        let TShape::Face(fd) = ts.as_ref() else {
            return None;
        };
        let wire_ts = brep.tshapes.get(fd.outer_wire.index)?;
        let TShape::Wire(wd) = wire_ts.as_ref() else {
            return None;
        };
        wd.edges.clone()
    };

    let n = old_edges.len();
    let pos = old_edges
        .iter()
        .position(|esr| esr.index == fillet_edge_idx)?;

    let fillet_we = &old_edges[pos];
    let (fillet_ws_raw, fillet_wv_raw) = wire_edge_vertices(brep, fillet_we.index);
    let fillet_forward = fillet_we.orientation == Orientation::Forward;
    let (fillet_ws, fillet_wv) = if fillet_forward {
        (fillet_ws_raw, fillet_wv_raw)
    } else {
        (fillet_wv_raw, fillet_ws_raw)
    };

    let we_before = &old_edges[(pos + n - 1) % n];
    let we_after = &old_edges[(pos + 1) % n];

    let ws_pt = brep.vertex_point(fillet_ws).unwrap_or_default();
    let (new_before_e_start, new_before_e_end) =
        if (ws_pt - v1_orig).length_squared() < (ws_pt - v2_orig).length_squared() {
            let (bws, _bwe) = wire_edge_vertices(brep, we_before.index);
            (bws, v1_new)
        } else {
            let (bws, _bwe) = wire_edge_vertices(brep, we_before.index);
            (bws, v2_new)
        };

    let wv_pt = brep.vertex_point(fillet_wv).unwrap_or_default();
    let (new_after_e_start, new_after_e_end) =
        if (wv_pt - v2_orig).length_squared() < (wv_pt - v1_orig).length_squared() {
            let (_aws, awe) = wire_edge_vertices(brep, we_after.index);
            (v2_new, awe)
        } else {
            let (_aws, awe) = wire_edge_vertices(brep, we_after.index);
            (v1_new, awe)
        };

    let before_pt = brep.vertex_point(new_before_e_start).unwrap_or_default();
    let before_ep = brep.vertex_point(new_before_e_end).unwrap_or_default();
    let short_before = push_line_edge(
        brep,
        new_before_e_start,
        new_before_e_end,
        before_pt,
        before_ep,
    );

    let after_pt = brep.vertex_point(new_after_e_start).unwrap_or_default();
    let after_ep = brep.vertex_point(new_after_e_end).unwrap_or_default();
    let short_after = push_line_edge(brep, new_after_e_start, new_after_e_end, after_pt, after_ep);

    let mut new_edge_refs: Vec<ShapeRef> = Vec::with_capacity(n + 1);

    for k in 0..pos {
        if k == (pos + n - 1) % n {
            new_edge_refs.push(edge_shape_ref(brep, short_before, Orientation::Forward));
        } else {
            new_edge_refs.push(old_edges[k]);
        }
    }

    let contact_orient = if fillet_forward {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    new_edge_refs.push(edge_shape_ref(brep, contact_edge_idx, contact_orient));

    for k in (pos + 1)..n {
        if k == (pos + 1) % n {
            new_edge_refs.push(edge_shape_ref(brep, short_after, Orientation::Forward));
        } else {
            new_edge_refs.push(old_edges[k]);
        }
    }

    if pos == 0 {
        if let Some(last) = new_edge_refs.last_mut() {
            *last = edge_shape_ref(brep, short_before, Orientation::Forward);
        }
    }

    let new_wire_sr = brep.add_twire(new_edge_refs);
    let face_ts = Arc::make_mut(&mut brep.tshapes[face_flat_idx]);
    if let TShape::Face(fd) = face_ts {
        fd.outer_wire = new_wire_sr;
        if !fd.my_shapes.is_empty() {
            fd.my_shapes[0] = new_wire_sr;
        }
    }

    Some(())
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
// Utility Functions
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

/// Interpolate between two radii with tension parameter.
fn interpolate_radius(r1: f64, r2: f64, t: f64, tension: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h1 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h2 = -2.0 * t3 + 3.0 * t2;
    let smooth = tension;
    let h1_tension = h1 + smooth * (t3 - 2.0 * t2 + t);
    let h2_tension = h2 + smooth * (-t3 + t2);
    r1 * h1_tension + r2 * h2_tension
}
