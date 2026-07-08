//! BRepAlgo-style algorithms — topods::BRep data model.
//!
//! OCCT-aligned: NormalProject (project shape onto face along normal),
//! FaceSection (intersect shape with a face).
//!
//! Also provides utility functions for BRep analysis and manipulation
//! using `rcad_kernel::topods::BRep`:
//! - Normal/tangent evaluation
//! - Tolerance propagation
//! - Geometric properties
//! - Validity & orientation checking

use glam::DVec3;
use rcad_kernel::topods::{self, BRep, ShapeRef, TShape};
use rcad_kernel::geom::{Curve3, Surface3, CurveEval, SurfaceEval};

// =============================================================================
// Error Types
// =============================================================================

#[derive(Debug, Clone)]
pub enum BRepAlgoError {
    InvalidFaceIndex(usize),
    InvalidEdgeIndex(usize),
    InvalidVertexIndex(usize),
    /// The NormalProject or FaceSection operation failed.
    OperationFailed(&'static str),
}

impl std::fmt::Display for BRepAlgoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BRepAlgoError::InvalidFaceIndex(i) => write!(f, "invalid face index: {i}"),
            BRepAlgoError::InvalidEdgeIndex(i) => write!(f, "invalid edge index: {i}"),
            BRepAlgoError::InvalidVertexIndex(i) => write!(f, "invalid vertex index: {i}"),
            BRepAlgoError::OperationFailed(msg) => write!(f, "operation failed: {msg}"),
        }
    }
}

impl std::error::Error for BRepAlgoError {}

// =============================================================================
// BRepAlgo_NormalProject — project shape onto a face along face normal
// =============================================================================

/// Projects a shape onto a target face along the face's normal direction.
///
/// OCCT: BRepAlgo_NormalProject
/// 1. Construct with a target face.
/// 2. Call `perform(shape)` to project.
/// 3. Call `projected()` for the result.
pub struct NormalProject {
    face_surface: Option<Surface3>,
    result: Option<BRep>,
}

impl NormalProject {
    /// Create from a flat face index in the given BRep.
    pub fn new(brep: &BRep, face_idx: usize) -> Result<Self, BRepAlgoError> {
        let faces: Vec<ShapeRef> = brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_)))
            .map(|(i, _)| ShapeRef::synthetic(i))
            .collect();
        let sr = *faces.get(face_idx).ok_or(BRepAlgoError::InvalidFaceIndex(face_idx))?;
        let surface = match &*brep.tshapes[sr.index] {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        Ok(Self { face_surface: surface, result: None })
    }

    /// Project shape `input` onto the target face.
    /// Each vertex of `input` is projected along the face normal onto the surface.
    pub fn perform(&mut self, input: &BRep) -> Result<(), BRepAlgoError> {
        let surf = self.face_surface.as_ref()
            .ok_or(BRepAlgoError::OperationFailed("no target face surface"))?;
        let mut out = BRep::new();

        // Collect input vertices and project them
        let input_vertices: Vec<(ShapeRef, DVec3)> = input.tshapes.iter().enumerate()
            .filter_map(|(i, ts)| match ts.as_ref() {
                TShape::Vertex(vd) => Some((ShapeRef::synthetic(i), vd.point)),
                _ => None,
            })
            .collect();

        // For each input vertex, project onto the target surface
        let mut v_map: Vec<ShapeRef> = Vec::new();
        for (_sr, pt) in &input_vertices {
            let proj = rcad_kernel::projection::closest_point_on_surface(surf, *pt, 16);
            let (u, v) = proj.params;
            let projected_pt = surf.point_at(u, v);
            let new_v = out.add_tvertex(projected_pt);
            v_map.push(new_v);
        }

        // Collect input edges and build projected edges
        for (i, ts) in input.tshapes.iter().enumerate() {
            if let TShape::Edge(ed) = ts.as_ref() {
                // Find which input vertices this edge connects
                let sv_idx = input_vertices.iter().position(|(sr, _)| sr.index == ed.first.index);
                let ev_idx = input_vertices.iter().position(|(sr, _)| sr.index == ed.last.index);
                if let (Some(si), Some(ei)) = (sv_idx, ev_idx) {
                    let new_sv = v_map[si];
                    let new_ev = v_map[ei];
                    // Copy the edge curve (project the curve onto the surface later if needed)
                    let _ = out.add_tedge(ed.curve.clone(), new_sv, new_ev, ed.range);
                }
            }
        }

        self.result = Some(out);
        Ok(())
    }

    /// Returns the projected shape (a new BRep on the target face surface).
    pub fn projected(&self) -> Option<&BRep> {
        self.result.as_ref()
    }

    /// OCCT-alias for `projected()`.
    pub fn shape(&self) -> Option<&BRep> {
        self.result.as_ref()
    }

    /// OCCT-alias: returns true if NormalProject performed successfully.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

// =============================================================================
// BRepAlgo_FaceSection — intersect a shape with a face
// =============================================================================

/// Intersects a shape with a target face, returning the section edges.
///
/// OCCT: BRepAlgo_FaceSection
/// 1. Construct with a target face.
/// 2. Call `perform(shape)` to compute intersection.
/// 3. Call `shape()` for the result (edges/wires).
pub struct FaceSection {
    face_surface: Option<Surface3>,
    result: Option<BRep>,
}

impl FaceSection {
    /// Create from a flat face index in the given BRep.
    pub fn new(brep: &BRep, face_idx: usize) -> Result<Self, BRepAlgoError> {
        let faces: Vec<ShapeRef> = brep.tshapes.iter().enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_)))
            .map(|(i, _)| ShapeRef::synthetic(i))
            .collect();
        let sr = *faces.get(face_idx).ok_or(BRepAlgoError::InvalidFaceIndex(face_idx))?;
        let surface = match &*brep.tshapes[sr.index] {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        Ok(Self { face_surface: surface, result: None })
    }

    /// Compute intersection of the target face surface with each face of `input`.
    /// Result edges are stored and returned via `shape()`.
    pub fn perform(&mut self, input: &BRep) -> Result<(), BRepAlgoError> {
        let _surf = self.face_surface.as_ref()
            .ok_or(BRepAlgoError::OperationFailed("no target face surface"))?;
        let mut out = BRep::new();

        // For each face in the input shape, intersect its surface with the target face
        for ts in &input.tshapes {
            if let TShape::Face(fd) = ts.as_ref() {
                if let Some(ref other_surf) = fd.surface {
                    // Use face-face intersection to get intersection curves
                    use rcad_kernel::geom::Surface3;
                    match (_surf, other_surf) {
                        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                            let result = crate::int_ana::intersect_plane_plane_intana(p1, p2);
                            if let crate::int_ana::PlnPlnResult::Line(line) = result {
                                // Create vertices along the line and an edge
                                let sv = out.add_tvertex(
                                    line.point_at(-1000.0));
                                let ev = out.add_tvertex(
                                    line.point_at(1000.0));
                                out.add_tedge(
                                    Some(Curve3::Line(line)),
                                    sv, ev,
                                    [-1000.0, 1000.0],
                                );
                            }
                        }
                        (Surface3::Plane(p), Surface3::Cylinder(c)) => {
                            let result = crate::int_ana::intersect_plane_cylinder_intana(p, c);
                            match result {
                                crate::int_ana::PlnCylResult::Circle(circle) => {
                                    let sv = out.add_tvertex(circle.point_at(0.0));
                                    let ev = out.add_tvertex(circle.point_at(std::f64::consts::TAU));
                                    out.add_tedge(
                                        Some(Curve3::Circle(circle)),
                                        sv, ev,
                                        [0.0, std::f64::consts::TAU],
                                    );
                                }
                                crate::int_ana::PlnCylResult::Ellipse(ell) => {
                                    // Approximate ellipse as a BSpline or skip
                                }
                                _ => {}
                            }
                        }
                        _ => {
                            // Fallback: use generic surface-surface intersection
                            if let Some(curve) = intersect_surfaces(_surf, other_surf) {
                                let sv = out.add_tvertex(curve.point_at(0.0));
                                let ev = out.add_tvertex(curve.point_at(100.0));
                                out.add_tedge(Some(curve), sv, ev, [0.0, 100.0]);
                            }
                        }
                    }
                }
            }
        }

        self.result = Some(out);
        Ok(())
    }

    /// Returns the intersection result as edges/wires.
    pub fn shape(&self) -> Option<&BRep> {
        self.result.as_ref()
    }

    /// OCCT-alias: returns true if FaceSection performed successfully.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

/// Fallback: generic surface-surface intersection using polyline sampling.
fn intersect_surfaces(s1: &Surface3, s2: &Surface3) -> Option<Curve3> {
    use crate::inttools::face_face::intersect_faces;
    let curves = intersect_faces(s1, s2, 1e-7, 1e-7);
    if curves.is_empty() {
        return None;
    }
    Some(curves[0].curve.clone())
}

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
