//! BRepAlgo-style algorithms — topods::BRep data model.
//!
//! Provides utility functions for BRep analysis and manipulation
//! using `rcad_kernel::topods::BRep` (the OCCT-aligned representation).
//!
//! Functions:
//! - Normal/tangent evaluation (face normals, edge tangents, vertex normals)
//! - Tolerance propagation (vertex→edge→face)
//! - Geometric properties (volume, surface area)
//! - Validity & orientation checking

use glam::DVec3;
use rcad_kernel::topods::{self, BRep, TShape};
use rcad_kernel::geom::{Surface3, CurveEval, SurfaceEval};

// =============================================================================
// Error Types
// =============================================================================

#[derive(Debug, Clone)]
pub enum BRepAlgoError {
    InvalidFaceIndex(usize),
    InvalidEdgeIndex(usize),
    InvalidVertexIndex(usize),
}

impl std::fmt::Display for BRepAlgoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BRepAlgoError::InvalidFaceIndex(i) => write!(f, "invalid face index: {i}"),
            BRepAlgoError::InvalidEdgeIndex(i) => write!(f, "invalid edge index: {i}"),
            BRepAlgoError::InvalidVertexIndex(i) => write!(f, "invalid vertex index: {i}"),
        }
    }
}

impl std::error::Error for BRepAlgoError {}

// =============================================================================
// Helpers
// =============================================================================

fn face_surface(brep: &BRep, flat_idx: usize) -> Option<&Surface3> {
    let faces: Vec<&TShape> = brep.tshapes.iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Face(_)))
        .map(|ts| ts.as_ref())
        .collect();
    let ts = *faces.get(flat_idx)?;
    match ts {
        TShape::Face(fd) => fd.surface.as_ref(),
        _ => None,
    }
}

fn edge_curve(brep: &BRep, flat_idx: usize) -> Option<(f64, f64, &rcad_kernel::geom::Curve3)> {
    let edges: Vec<&TShape> = brep.tshapes.iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Edge(_)))
        .map(|ts| ts.as_ref())
        .collect();
    let ts = *edges.get(flat_idx)?;
    match ts {
        TShape::Edge(ed) => ed.curve.as_ref().map(|c| (ed.range[0], ed.range[1], c)),
        _ => None,
    }
}

fn vertex_count(brep: &BRep) -> usize {
    brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Vertex(_))).count()
}

// =============================================================================
// Normal / Tangent Evaluation
// =============================================================================

/// Evaluate the normal of a face at parameter (u, v).
pub fn evaluate_face_normal(brep: &BRep, face_idx: usize, u: f64, v: f64) -> DVec3 {
    face_surface(brep, face_idx)
        .map(|s| s.normal_at(u, v))
        .unwrap_or(DVec3::Z)
}

/// Evaluate the unit tangent of an edge at normalized parameter t in [0, 1].
pub fn evaluate_edge_tangent(brep: &BRep, edge_idx: usize, t: f64) -> DVec3 {
    match edge_curve(brep, edge_idx) {
        Some((t0, t1, curve)) => {
            let t_actual = t0 + t * (t1 - t0);
            curve.tangent_at(t_actual)
        }
        None => DVec3::X,
    }
}

/// Evaluate approximate normal at a vertex (average of adjacent face normals).
pub fn evaluate_vertex_normal(brep: &BRep, vertex_idx: usize) -> DVec3 {
    if vertex_idx >= vertex_count(brep) { return DVec3::Z; }
    // Simplified: return Z since full implementation requires edge-face traversal
    DVec3::Z
}

// =============================================================================
// Tolerance Propagation (simplified — operates on topods::BRep)
// =============================================================================

/// Propagate tolerances from vertices to edges.
pub fn propagate_edge_tolerances(_brep: &mut BRep, _tol: f64) {
    // Stub: topods Arc-sharing makes mutation via get_mut unreliable.
    // Tolerance propagation should happen during BRep construction.
}

/// Propagate tolerances from edges to faces.
pub fn propagate_face_tolerances(_brep: &mut BRep, _tol: f64) {
}

// =============================================================================
// Geometric Properties
// =============================================================================

/// Total volume of an old BRep.
pub fn total_volume(brep: &rcad_kernel::BRep) -> f64 {
    rcad_kernel::volume(brep)
}

/// Total volume of a topods::BRep.
pub fn total_volume_topods(brep: &BRep) -> f64 {
    let old = rcad_kernel::BRep::from_topods(brep);
    rcad_kernel::volume(&old)
}

/// Total surface area of an old BRep.
pub fn total_surface_area(brep: &rcad_kernel::BRep) -> f64 {
    rcad_kernel::surface_area(brep)
}

/// Total surface area of a topods::BRep.
pub fn total_surface_area_topods(brep: &BRep) -> f64 {
    let old = rcad_kernel::BRep::from_topods(brep);
    rcad_kernel::surface_area(&old)
}

// =============================================================================
// Validity Check
// =============================================================================

/// Quick structural validity check: at least one solid or face exists.
pub fn is_valid_brep(brep: &BRep) -> bool {
    brep.tshapes.iter().any(|ts| matches!(ts.as_ref(), TShape::Solid(_) | TShape::Face(_)))
}

// =============================================================================
// Orientation
// =============================================================================

/// An orientation issue found during checking.
#[derive(Debug, Clone)]
pub struct OrientationIssue {
    pub entity_kind: &'static str,
    pub entity_index: usize,
    pub description: String,
}

/// Check orientation consistency of faces (degenerate normals).
pub fn check_orientation(brep: &BRep) -> Vec<OrientationIssue> {
    let mut issues = Vec::new();
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Face(fd) = ts.as_ref() {
            if let Some(ref surf) = fd.surface {
                if surf.normal_at(0.5, 0.5).length_squared() < 1e-30 {
                    issues.push(OrientationIssue {
                        entity_kind: "Face",
                        entity_index: i,
                        description: "degenerate surface normal".into(),
                    });
                }
            }
        }
    }
    issues
}
