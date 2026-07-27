//! BRepAlgo-style algorithms — classes.
//!
//! NormalProject (project shape onto face along normal),
//! FaceSection (intersect shape with a face).

use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topods::{self, BRep, ShapeRef, TShape};

#[derive(Debug, Clone)]
pub enum BRepAlgoError {
    InvalidFaceIndex(usize),
    InvalidEdgeIndex(usize),
    InvalidVertexIndex(usize),
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
// BRepAlgo_NormalProject
// =============================================================================

pub struct NormalProject {
    face_surface: Option<Surface3>,
    result: Option<BRep>,
}

impl NormalProject {
    pub fn new(brep: &BRep, face_idx: usize) -> Result<Self, BRepAlgoError> {
        let faces: Vec<ShapeRef> = brep
            .tshapes
            .iter()
            .enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_)))
            .map(|(i, _)| ShapeRef::synthetic(i))
            .collect();
        let sr = *faces
            .get(face_idx)
            .ok_or(BRepAlgoError::InvalidFaceIndex(face_idx))?;
        let surface = match &*brep.tshapes[sr.index] {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        Ok(Self {
            face_surface: surface,
            result: None,
        })
    }

    pub fn perform(&mut self, input: &BRep) -> Result<(), BRepAlgoError> {
        let surf = self
            .face_surface
            .as_ref()
            .ok_or(BRepAlgoError::OperationFailed("no target face surface"))?;
        let mut out = BRep::new();
        let input_vertices: Vec<(ShapeRef, DVec3)> = input
            .tshapes
            .iter()
            .enumerate()
            .filter_map(|(i, ts)| match ts.as_ref() {
                TShape::Vertex(vd) => Some((ShapeRef::synthetic(i), vd.point)),
                _ => None,
            })
            .collect();
        let mut v_map: Vec<ShapeRef> = Vec::new();
        for (_sr, pt) in &input_vertices {
            let proj = rcad_kernel::projection::closest_point_on_surface(surf, *pt, 16);
            let (u, v) = proj.params;
            let projected_pt = surf.point_at(u, v);
            let new_v = out.add_tvertex(projected_pt);
            v_map.push(new_v);
        }
        for (i, ts) in input.tshapes.iter().enumerate() {
            if let TShape::Edge(ed) = ts.as_ref() {
                let sv_idx = input_vertices
                    .iter()
                    .position(|(sr, _)| sr.index == ed.first.index);
                let ev_idx = input_vertices
                    .iter()
                    .position(|(sr, _)| sr.index == ed.last.index);
                if let (Some(si), Some(ei)) = (sv_idx, ev_idx) {
                    let new_sv = v_map[si];
                    let new_ev = v_map[ei];
                    let _ = out.add_tedge(ed.curve.clone(), new_sv, new_ev, ed.range);
                }
            }
        }
        self.result = Some(out);
        Ok(())
    }

    pub fn projected(&self) -> Option<&BRep> {
        self.result.as_ref()
    }
    pub fn shape(&self) -> Option<&BRep> {
        self.result.as_ref()
    }
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

// =============================================================================
// BRepAlgo_FaceSection
// =============================================================================

pub struct FaceSection {
    face_surface: Option<Surface3>,
    result: Option<BRep>,
}

impl FaceSection {
    pub fn new(brep: &BRep, face_idx: usize) -> Result<Self, BRepAlgoError> {
        let faces: Vec<ShapeRef> = brep
            .tshapes
            .iter()
            .enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Face(_)))
            .map(|(i, _)| ShapeRef::synthetic(i))
            .collect();
        let sr = *faces
            .get(face_idx)
            .ok_or(BRepAlgoError::InvalidFaceIndex(face_idx))?;
        let surface = match &*brep.tshapes[sr.index] {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        Ok(Self {
            face_surface: surface,
            result: None,
        })
    }

    pub fn perform(&mut self, input: &BRep) -> Result<(), BRepAlgoError> {
        let _surf = self
            .face_surface
            .as_ref()
            .ok_or(BRepAlgoError::OperationFailed("no target face surface"))?;
        let mut out = BRep::new();
        for ts in &input.tshapes {
            if let TShape::Face(fd) = ts.as_ref() {
                if let Some(ref other_surf) = fd.surface {
                    use rcad_kernel::geom::Surface3;
                    match (_surf, other_surf) {
                        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                            let result = crate::int_ana::intersect_plane_plane_intana(p1, p2);
                            if let crate::int_ana::PlnPlnResult::Line(line) = result {
                                let sv = out.add_tvertex(line.point_at(-1000.0));
                                let ev = out.add_tvertex(line.point_at(1000.0));
                                out.add_tedge(Some(Curve3::Line(line)), sv, ev, [-1000.0, 1000.0]);
                            }
                        }
                        (Surface3::Plane(p), Surface3::Cylinder(c)) => {
                            let result = crate::int_ana::intersect_plane_cylinder_intana(p, c);
                            if let crate::int_ana::PlnCylResult::Circle(circle) = result {
                                let sv = out.add_tvertex(circle.point_at(0.0));
                                let ev = out.add_tvertex(circle.point_at(std::f64::consts::TAU));
                                out.add_tedge(
                                    Some(Curve3::Circle(circle)),
                                    sv,
                                    ev,
                                    [0.0, std::f64::consts::TAU],
                                );
                            }
                        }
                        _ => {
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

    pub fn shape(&self) -> Option<&BRep> {
        self.result.as_ref()
    }
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

fn intersect_surfaces(s1: &Surface3, s2: &Surface3) -> Option<Curve3> {
    use crate::inttools::face_face::intersect_faces;
    let curves = intersect_faces(s1, s2, 1e-7, 1e-7);
    curves.into_iter().next().map(|c| c.curve)
}
