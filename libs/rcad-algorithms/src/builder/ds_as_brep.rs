/// OCCT BRepTool adaptor over the existing DS data source.
///
/// During Phase 1 migration, this allows TopoDS-based wire path code to
/// read from the existing DS by wrapping DS indices as ShapeRef handles.
use std::collections::HashMap;
use glam::DVec3;
use crate::bopds::ds::DS;
use rcad_kernel::geom::*;
use rcad_kernel::topods::{BRepTool, ShapeRef, Orientation};

/// Adaptor: wraps DS + face_idx as a BRepTool, mapping ShapeRef.index → DS array index.
///
/// ShapeRef values used with this adaptor:
/// - Vertex ShapeRef.index = DS vertex index
/// - Edge ShapeRef.index = DS edge index
/// - Face ShapeRef.index = DS face index
///
/// pcurves and vertex_params are stored per edge in a lookup map built at construction time.
pub(crate) struct DSAsBRep<'a> {
    pub ds: &'a DS,
    pub face_idx: usize,
    /// edge_index → (pc_on_face, t_first, t_last) — built once from DSCurveRepOnFace
    pub pcurve_cache: HashMap<usize, (Curve2d, f64, f64)>,
    /// edge_index → (vertex_index → param) — built once from DSEdge.vertex_params
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
        let vertex_param_cache: HashMap<usize, HashMap<usize, f64>> = ds.edges.iter()
            .enumerate()
            .map(|(ei, edge)| (ei, edge.vertex_params.clone()))
            .collect();

        DSAsBRep { ds, face_idx, pcurve_cache, vertex_param_cache }
    }
}

impl BRepTool for DSAsBRep<'_> {
    fn vertex_position(&self, v: ShapeRef) -> DVec3 {
        self.ds.vertices.get(v.index).map(|v| v.point).unwrap_or(DVec3::ZERO)
    }

    fn vertex_tolerance(&self, v: ShapeRef) -> f64 {
        self.ds.vertices.get(v.index).map(|v| v.geom_tol).unwrap_or(0.0)
    }

    fn is_edge_degenerated(&self, e: ShapeRef) -> bool {
        self.ds.is_edge_degenerated(e.index)
    }

    fn edge_other_vertex(&self, edge: ShapeRef, v: ShapeRef) -> ShapeRef {
        if let Some(e) = self.ds.edges.get(edge.index) {
            if e.start_vertex == v.index {
                ShapeRef::new(e.end_vertex)
            } else {
                ShapeRef::new(e.start_vertex)
            }
        } else { v }
    }

    fn first_vertex(&self, edge: ShapeRef) -> ShapeRef {
        self.ds.edges.get(edge.index)
            .map(|e| ShapeRef::new(e.start_vertex))
            .unwrap_or(edge)
    }

    fn last_vertex(&self, edge: ShapeRef) -> ShapeRef {
        self.ds.edges.get(edge.index)
            .map(|e| ShapeRef::new(e.end_vertex))
            .unwrap_or(edge)
    }

    fn oriented_first_vertex(&self, edge: ShapeRef, orientation: Orientation) -> ShapeRef {
        self.ds.edges.get(edge.index).map(|e| {
            let vi = if orientation == Orientation::Reversed { e.end_vertex } else { e.start_vertex };
            ShapeRef::new(vi)
        }).unwrap_or(edge)
    }

    fn parameter_on_edge(&self, vertex: ShapeRef, edge: ShapeRef, _face: ShapeRef) -> Option<f64> {
        self.vertex_param_cache.get(&edge.index)
            .and_then(|vpm| vpm.get(&vertex.index).copied())
    }

    fn curve_on_surface(&self, edge: ShapeRef, _face: ShapeRef) -> Option<&(Curve2d, f64, f64)> {
        self.pcurve_cache.get(&edge.index)
    }

    fn face_surface(&self, _face: ShapeRef) -> Option<&Surface3> {
        Some(&self.ds.faces[self.face_idx].surface)
    }

    fn face_surface_world(&self, _face: ShapeRef) -> Option<Surface3> {
        Some(self.ds.faces[self.face_idx].surface.clone())
    }

    fn edge_curve_world(&self, edge: ShapeRef) -> Option<(Curve3, [f64; 2])> {
        self.ds.edges.get(edge.index).map(|e| (e.curve.clone(), e.t_range))
    }

    fn u_resolution(&self, _face: ShapeRef, tol3d: f64) -> f64 {
        // Fallback: use the face surface from DS
        let surf = &self.ds.faces[self.face_idx].surface;
        match surf {
            Surface3::Sphere(s) => tol3d / s.radius.max(1e-15),
            Surface3::Cylinder(c) => tol3d / c.radius.max(1e-15),
            Surface3::Cone(_) => tol3d * 1e-3,
            Surface3::Torus(t) => tol3d / t.major_radius.max(1e-15),
            _ => tol3d,
        }
    }

    fn v_resolution(&self, _face: ShapeRef, tol3d: f64) -> f64 {
        let surf = &self.ds.faces[self.face_idx].surface;
        match surf {
            Surface3::Sphere(s) => tol3d / s.radius.max(1e-15),
            Surface3::Cylinder(_) => tol3d,
            Surface3::Cone(_) => tol3d,
            Surface3::Torus(t) => tol3d / t.minor_radius.max(1e-15),
            _ => tol3d,
        }
    }
}
