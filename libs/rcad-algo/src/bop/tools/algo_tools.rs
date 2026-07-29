// OCCT BOPTools_AlgoTools — utility functions for boolean operations.
//
// OCCT BOPTools_AlgoTools.cxx / _1.cxx / _2.cxx
// Functions translated 1:1 from OCCT.

use crate::bop::ds::DS;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

// ====================================================================
// DTolerance — OCCT BOPTools_AlgoTools::DTolerance()
// ====================================================================
/// Additional tolerance used in Boolean Operations (1.e-12).
pub fn d_tolerance() -> f64 { 1e-12 }

// ====================================================================
// ComputeVV — OCCT BOPTools_AlgoTools::ComputeVV
// ====================================================================
/// Intersects the vertex with a point.
/// Returns 0 if vertex intersects the point (distance <= tolerance sum).
pub fn compute_vv_vertex_point(v_tol: f64, v_pnt: glam::DVec3,
                                p_pnt: glam::DVec3, p_tol: f64) -> i32 {
    if (v_pnt - p_pnt).length() <= v_tol + p_tol { 0 } else { 1 }
}

/// Intersects two vertices. Returns 0 if they interfere.
pub fn compute_vv(v1_tol: f64, v1_pnt: glam::DVec3,
                  v2_tol: f64, v2_pnt: glam::DVec3, fuzz: f64) -> i32 {
    let tol = v1_tol.max(v2_tol).max(fuzz);
    if (v1_pnt - v2_pnt).length() <= tol { 0 } else { 1 }
}

// ====================================================================
// MakeNewVertex — OCCT BOPTools_AlgoTools::MakeNewVertex (2arg)
// ====================================================================
/// Creates a vertex at the midpoint of two vertices with combined tolerance.
pub fn make_new_vertex(v1_pnt: glam::DVec3, v1_tol: f64,
                       v2_pnt: glam::DVec3, v2_tol: f64) -> (glam::DVec3, f64) {
    let mid = (v1_pnt + v2_pnt) * 0.5;
    let tol = (mid - v1_pnt).length() + v1_tol.max(v2_tol);
    (mid, tol)
}

/// Creates a vertex from a point with given tolerance.
pub fn make_new_vertex_point(p: glam::DVec3, tol: f64) -> (glam::DVec3, f64) {
    (p, tol)
}

// ====================================================================
// IsMicroEdge — OCCT BOPTools_AlgoTools::IsMicroEdge (L2075+)
// ====================================================================
/// Checks if the edge is too short (micro edge) — range < 2 * tolerance.
pub fn is_micro_edge(range_len: f64, edge_tol: f64) -> bool {
    range_len < 2.0 * edge_tol.max(rcad_kernel::CONFUSION)
}

// ====================================================================
// ComputeTolerance — OCCT BOPTools_AlgoTools::ComputeTolerance (L1093-1111)
// ====================================================================
/// Computes max deviation of edge's pcurve on face surface vs 3D curve.
/// Returns (max_distance, max_parameter) or None if edge lacks pcurve.
pub fn compute_tolerance(ds: &DS, edge_idx: usize, face_idx: usize) -> Option<(f64, f64)> {
    let curve = ds.edge_curve(edge_idx)?.clone();
    let range = ds.edge_range(edge_idx);
    let surf = ds.face_surface(face_idx)?;
    let shape = &ds.shapes[edge_idx].shape;
    let pcurve_info = shape.as_edge()?.pcurves.get(&face_idx)?.clone();
    let (pcurve, _f, _l) = pcurve_info;
    let n = 23usize;
    let dt = (range[1] - range[0]) / n as f64;
    let mut max_dist = 0.0f64;
    let mut max_par = range[0];
    for i in 0..=n {
        let t = range[0] + i as f64 * dt;
        let p3d = curve.point_at(t);
        let uv = pcurve.point_at(t);
        let p_surf = surf.point_at(uv.x, uv.y);
        let d = (p3d - p_surf).length();
        if d > max_dist { max_dist = d; max_par = t; }
    }
    Some((max_dist, max_par))
}

// ====================================================================
// TreatCompound — OCCT BOPTools_AlgoTools::TreatCompound
// ====================================================================
/// Flattens a compound shape into a list of non-compound sub-shapes.
pub fn treat_compound(shape: &Shape) -> Vec<Shape> {
    let mut result = Vec::new();
    let mut stack = vec![shape.clone()];
    while let Some(s) = stack.pop() {
        match &*s.data {
            TShape::Compound(children) => {
                for child in children.iter().rev() {
                    stack.push(child.clone());
                }
            }
            _ => result.push(s),
        }
    }
    result
}

// ====================================================================
// Dimensions — OCCT BOPTools_AlgoTools::Dimensions (L546)
// ====================================================================
/// Returns min and max dimension of a shape.
/// VERTEX → (0,0), EDGE → (1,1), FACE → (2,2), SHELL/SOLID → (3,3).
pub fn dimensions(st: ShapeType) -> (i32, i32) {
    match st {
        ShapeType::Vertex => (0, 0),
        ShapeType::Edge => (1, 1),
        ShapeType::Wire => (1, 1),
        ShapeType::Face => (2, 2),
        ShapeType::Shell => (3, 3),
        ShapeType::Solid => (3, 3),
        ShapeType::Compound => (0, 3),
        ShapeType::CompSolid => (3, 3),
        _ => (0, 3),
    }
}

// ====================================================================
// IsHole — OCCT BOPTools_AlgoTools::IsHole (L291)
// ====================================================================
/// Checks if the wire is a hole for the face.
/// rcad: uses signed area of projected wire in 2D.
pub fn is_hole(wire_edges: &[usize], ds: &DS, face_idx: usize) -> bool {
    let surf = match ds.face_surface(face_idx) {
        Some(s) => s,
        None => return false,
    };
    use rcad_kernel::geom::SurfaceEval;
    let mut area_2d = 0.0;
    for &ei in wire_edges {
        let curve = match ds.edge_curve(ei) {
            Some(c) => c.clone(),
            None => continue,
        };
        let range = ds.edge_range(ei);
        let n = 8usize;
        let dt = (range[1] - range[0]) / n as f64;
        let mut prev = {
            let t0 = range[0];
            let p3d = curve.point_at(t0);
            let (_u, v) = closest_uv(&surf, p3d);
            (0.0, v)
        };
        for i in 1..=n {
            let t = range[0] + i as f64 * dt;
            let p3d = curve.point_at(t);
            let (u, v) = closest_uv(&surf, p3d);
            area_2d += (prev.1 + v) * (u - prev.0) * 0.5;
            prev = (u, v);
        }
    }
    area_2d < 0.0
}

fn closest_uv(surf: &rcad_kernel::geom::Surface3, p: glam::DVec3) -> (f64, f64) {
    use crate::bop::closest_point_on_surface;
    let (uv, _) = closest_point_on_surface(surf, p);
    (uv.x, uv.y)
}
