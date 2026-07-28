use crate::bopds::ds::DS;
use crate::tolerance::TOLERANCE_CLAMP_MIN;
use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods::{self, BRepTool, Orientation, Shape, ShapeType};
/// OCCT BRepTool adaptor over the existing DS data source.
///
/// During Phase 1 migration, this allows TopoDS-based wire path code to
/// read from the existing DS by wrapping DS indices as Shape handles.
use std::collections::HashMap;

/// Adaptor: wraps DS + face_idx as a BRepTool, mapping Shape.index  ?DS array index.
///
/// Shape values used with this adaptor:
/// - Vertex Shape.index = DS vertex index
/// - Edge Shape.index = DS edge index
/// - Face Shape.index = DS face index
///
/// pcurves and vertex_params are stored per edge in a lookup map built at construction time.
pub(crate) struct DSAsBRep<'a> {
    pub ds: &'a DS,
    pub face_idx: usize,
    /// edge_index  ?(pc_on_face, t_first, t_last)  ?built once from DSCurveRepOnFace
    pub pcurve_cache: HashMap<usize, (Curve2d, f64, f64)>,
    /// edge_index  ?(vertex_index  ?param)  ?built once from DSEdge.vertex_params
    pub vertex_param_cache: HashMap<usize, HashMap<usize, f64>>,
}

impl<'a> DSAsBRep<'a> {
    pub fn new(ds: &'a DS, face_idx: usize) -> Self {
        // Build pcurve cache: for each edge that has a face_rep for this face,
        // extract the pcurve and range.
        let mut pcurve_cache: HashMap<usize, (Curve2d, f64, f64)> = HashMap::new();
        for (ei, edge) in ds.edges.iter().enumerate() {
            for rep in &edge.face_reps {
                if rep.face_idx == face_idx {
                    let pc = rep.pcurve.clone();
                    let t_first = rep.start_param;
                    let t_last = rep.end_param;
                    pcurve_cache.insert(ei, (pc, t_first, t_last));
                }
            }
        }
        // Build vertex param cache from DSEdge.vertex_params
        let vertex_param_cache: HashMap<usize, HashMap<usize, f64>> = ds
            .edges
            .iter()
            .enumerate()
            .map(|(ei, edge)| (ei, edge.vertex_params.clone()))
            .collect();

        DSAsBRep {
            ds,
            face_idx,
            pcurve_cache,
            vertex_param_cache,
        }
    }
}

impl BRepTool for DSAsBRep<'_> {
    fn vertex_position(&self, v: &Shape) -> DVec3 {
        self.ds
            .vertices
            .get(v.index)
            .map(|v| v.point)
            .unwrap_or(DVec3::ZERO)
    }

    fn vertex_tolerance(&self, v: &Shape) -> f64 {
        self.ds
            .vertices
            .get(v.index)
            .map(|v| v.geom_tol)
            .unwrap_or(0.0)
    }

    fn is_edge_degenerated(&self, e: &Shape) -> bool {
        self.ds.is_edge_degenerated(e.index)
    }

    fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape {
        if let Some(e) = self.ds.edges.get(edge.index) {
            if e.start_vertex == v.index {
                Shape::synthetic(e.end_vertex, Orientation::Forward)
            } else {
                Shape::synthetic(e.start_vertex, Orientation::Forward)
            }
        } else {
            v.clone()
        }
    }

    fn first_vertex(&self, edge: &Shape) -> Shape {
        self.ds
            .edges
            .get(edge.index)
            .map(|e| Shape::synthetic(e.start_vertex, Orientation::Forward))
            .unwrap_or_else(|| edge.clone())
    }

    fn last_vertex(&self, edge: &Shape) -> Shape {
        self.ds
            .edges
            .get(edge.index)
            .map(|e| Shape::synthetic(e.end_vertex, Orientation::Forward))
            .unwrap_or_else(|| edge.clone())
    }

    fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape {
        self.ds
            .edges
            .get(edge.index)
            .map(|e| {
                let vi = if orientation == Orientation::Reversed {
                    e.end_vertex
                } else {
                    e.start_vertex
                };
                Shape::synthetic(vi, orientation)
            })
            .unwrap_or_else(|| edge.clone())
    }

    fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64> {
        self.vertex_param_cache
            .get(&edge.index)
            .and_then(|vpm| vpm.get(&vertex.index).copied())
    }

    fn curve_on_surface(&self, edge: &Shape, _face: &Shape) -> Option<&(Curve2d, f64, f64)> {
        self.pcurve_cache.get(&edge.index)
    }

    fn face_surface(&self, _face: &Shape) -> Option<&Surface3> {
        Some(self.ds.face_surface(self.face_idx).unwrap())
    }

    fn vertex_orientation(&self, v: &Shape) -> Orientation {
        let is_internal = self.ds.vertex_is_internal(v.index);
        let is_new = self
            .ds
            .vertices
            .get(v.index)
            .map_or(false, |dv| dv.origin.is_none());
        if is_internal || is_new {
            Orientation::Internal
        } else {
            Orientation::Forward
        }
    }

    fn face_surface_world(&self, _face: &Shape) -> Option<Surface3> {
        Some(self.ds.face_surface(self.face_idx).cloned().unwrap())
    }

    fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])> {
        self.ds
            .edges
            .get(edge.index)
            .map(|e| (e.curve.clone(), e.t_range))
    }

    fn u_resolution(&self, _face: &Shape, tol3d: f64) -> f64 {
        // Fallback: use the face surface from DS
        let surf = self.ds.face_surface(self.face_idx).unwrap();
        match surf {
            Surface3::Sphere(s) => tol3d / s.radius.max(TOLERANCE_CLAMP_MIN),
            Surface3::Cylinder(c) => tol3d / c.radius.max(TOLERANCE_CLAMP_MIN),
            Surface3::Cone(_) => tol3d * 1e-3,
            Surface3::Torus(t) => tol3d / t.major_radius.max(TOLERANCE_CLAMP_MIN),
            _ => tol3d,
        }
    }

    fn v_resolution(&self, _face: &Shape, tol3d: f64) -> f64 {
        let surf = self.ds.face_surface(self.face_idx).unwrap();
        match surf {
            Surface3::Sphere(s) => tol3d / s.radius.max(TOLERANCE_CLAMP_MIN),
            Surface3::Cylinder(_) => tol3d,
            Surface3::Cone(_) => tol3d,
            Surface3::Torus(t) => tol3d / t.minor_radius.max(TOLERANCE_CLAMP_MIN),
            _ => tol3d,
        }
    }

    fn tolerance(&self, _s: &Shape) -> f64 {
        0.0
    }

    fn shape_type(&self, s: &Shape) -> ShapeType {
        if self.ds.vertices.get(s.index).is_some() {
            ShapeType::Vertex
        } else if self.ds.edges.get(s.index).is_some() {
            ShapeType::Edge
        } else {
            ShapeType::Shape
        }
    }

    fn has_flag(&self, _s: &Shape, _flag: u16) -> bool {
        false
    }

    fn edge_data(&self, _e: &Shape) -> Option<&topods::TEdgeData> {
        None
    }

    fn face_data(&self, _f: &Shape) -> Option<&topods::TFaceData> {
        None
    }
}
