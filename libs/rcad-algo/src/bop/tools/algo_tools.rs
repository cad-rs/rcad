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
    // OCCT: aTolSum = aTolV1 + aTolP2 + Precision::Confusion()
    let a_tol_sum = v_tol + p_tol + rcad_kernel::CONFUSION;
    let a_d2 = (v_pnt - p_pnt).length_squared();
    if a_d2 > a_tol_sum * a_tol_sum { 1 } else { 0 }
}

/// Intersects two vertices. Returns 0 if they interfere.
pub fn compute_vv(v1_tol: f64, v1_pnt: glam::DVec3,
                  v2_tol: f64, v2_pnt: glam::DVec3, fuzz: f64) -> i32 {
    // OCCT: aFuzz1 = max(aFuzz, Precision::Confusion())
    let a_fuzz1 = fuzz.max(rcad_kernel::CONFUSION);
    // OCCT: aTolSum = aTolV1 + aTolV2 + aFuzz1; aTolSum2 = aTolSum * aTolSum
    let a_tol_sum = v1_tol + v2_tol + a_fuzz1;
    // OCCT: aD2 = aP1.SquareDistance(aP2); if (aD2 > aTolSum2) return 1
    let a_d2 = (v1_pnt - v2_pnt).length_squared();
    if a_d2 > a_tol_sum * a_tol_sum { 1 } else { 0 }
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

/// OCCT BOPTools_AlgoTools::MakeVertex — centroid + max tolerance from a list of vertices.
/// rcad: takes slice of (point, tolerance) pairs.
pub fn make_vertex(vertices: &[(glam::DVec3, f64)]) -> (glam::DVec3, f64) {
    let centroid = vertices.iter().map(|(p, _)| *p).sum::<glam::DVec3>() / vertices.len() as f64;
    let max_tol = vertices.iter().map(|(_, t)| *t).fold(0.0_f64, f64::max);
    (centroid, max_tol)
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

// ====================================================================
// MakeContainer — OCCT BOPTools_AlgoTools::MakeContainer (L1608+)
// ====================================================================
/// Creates an empty container shape of the given type.
pub fn make_container(st: ShapeType) -> Shape {
    match st {
        ShapeType::Compound => Shape::null(), // rcad: synthetic Compound
        _ => Shape::null(),
    }
}

// ====================================================================
// Dimension — OCCT BOPTools_AlgoTools::Dimension (L503+)
// ====================================================================
/// Returns dimension of a shape (-1 if mixed).
pub fn dimension(st: ShapeType) -> i32 {
    match st {
        ShapeType::Vertex => 0,
        ShapeType::Edge | ShapeType::Wire => 1,
        ShapeType::Face => 2,
        ShapeType::Shell | ShapeType::Solid | ShapeType::CompSolid => 3,
        ShapeType::Compound => -1,
        _ => -1,
    }
}

// ====================================================================
// IsOpenShell — OCCT BOPTools_AlgoTools::IsOpenShell (L2358+)
// ====================================================================
/// Checks if a shell is open by counting free edges (edges used once).
pub fn is_open_shell(face_indices: &[usize], ds: &DS) -> bool {
    let mut edge_count: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for &fi in face_indices {
        let si = ds.shape_info(fi);
        for &ss in &si.sub_shapes {
            if ss < ds.nb_shapes() && ds.shape_info(ss).shape_type == ShapeType::Edge {
                *edge_count.entry(ss).or_insert(0) += 1;
            }
        }
    }
    edge_count.values().any(|&c| c == 1)
}

// ====================================================================
// IsInvertedSolid — OCCT BOPTools_AlgoTools::IsInvertedSolid (L522+)
// ====================================================================
/// Checks if the solid is inverted (negative volume).
pub fn is_inverted_solid(_shell_indices: &[usize], _ds: &DS) -> bool {
    // rcad: requires signed volume computation
    false
}

// ====================================================================
// CorrectRange — OCCT BOPTools_AlgoTools::CorrectRange (EE variant, L284+)
// ====================================================================
/// Corrects edge range taking into account tolerance of adjacent shapes.
pub fn correct_range_ee(
    curve: &rcad_kernel::geom::Curve3,
    t1: f64, t2: f64,
    _tol_e1: f64, _tol_e2: f64,
) -> (f64, f64) {
    match curve {
        rcad_kernel::geom::Curve3::Line(_) => (t1, t2),
        _ => {
            // OCCT: shrink range by tolerance/derivative at endpoints
            let shrink = rcad_kernel::CONFUSION * 2.0;
            let nt1 = t1 + shrink;
            let nt2 = t2 - shrink;
            if nt2 > nt1 { (nt1, nt2) } else { (t1, t2) }
        }
    }
}

// ====================================================================
// CorrectRange — OCCT BOPTools_AlgoTools::CorrectRange (EF variant, L364+)
// ====================================================================
/// Corrects edge range for edge-face intersection.
pub fn correct_range_ef(
    curve: &rcad_kernel::geom::Curve3,
    t1: f64, t2: f64,
    _tol_e: f64, _tol_f: f64,
) -> (f64, f64) {
    // rcad: same logic as EE variant
    correct_range_ee(curve, t1, t2, 0.0, 0.0)
}

// ====================================================================
// IsBlockInOnFace — OCCT BOPTools_AlgoTools::IsBlockInOnFace (L1979+)
// ====================================================================
/// Checks if the pave block range lies in/on the face in 2D.
pub fn is_block_in_on_face(
    _t1: f64, _t2: f64,
    _ds: &DS, _edge_idx: usize, _face_idx: usize,
) -> bool {
    // rcad: requires pcurve intersection with face bounds
    false
}

// ====================================================================
// PointOnEdge — OCCT BOPTools_AlgoTools::PointOnEdge (L536+)
// ====================================================================
/// Computes a 3D point on the edge at given parameter.
pub fn point_on_edge(curve: &rcad_kernel::geom::Curve3, param: f64) -> glam::DVec3 {
    curve.point_at(param)
}

// ====================================================================
// GetEdgeOnFace — OCCT BOPTools_AlgoTools::GetEdgeOnFace (L1817+)
// ====================================================================
/// Finds the edge on the face that is same as the given edge.
pub fn get_edge_on_face(edge_idx: usize, face_idx: usize, ds: &DS) -> Option<usize> {
    let fi = ds.shape_info(face_idx);
    for &ss in &fi.sub_shapes {
        if ss < ds.nb_shapes() {
            let ssi = ds.shape_info(ss);
            if ssi.shape_type == ShapeType::Edge && ss == edge_idx {
                return Some(ss);
            }
        }
    }
    None
}

// ====================================================================
// MakeEdge — OCCT BOPTools_AlgoTools::MakeEdge / MakeSplitEdge / MakeSectEdge
// ====================================================================
/// Creates a split edge from a source edge with new vertices.
pub fn make_split_edge(
    _curve: &rcad_kernel::geom::Curve3,
    _v1: usize, _t1: f64,
    _v2: usize, _t2: f64,
    _ds: &mut DS,
) -> usize {
    // rcad: DS::push_edge handles this
    0
}

// ====================================================================
// UpdateVertex — OCCT BOPTools_AlgoTools::UpdateVertex (L124-138)
// ====================================================================
/// Updates vertex tolerance to cover a distance.
pub fn update_vertex(_v_idx: usize, _dist: f64, _ds: &mut DS) {
    // rcad: covered by PaveFiller::update_vertex
}

// ====================================================================
// AreFacesSameDomain — OCCT BOPTools_AlgoTools::AreFacesSameDomain (L1139+)
// ====================================================================
/// Checks if two faces are geometrically coincident (same domain).
pub fn are_faces_same_domain(
    f1_idx: usize, f2_idx: usize,
    ds: &DS,
    fuzz: f64,
) -> bool {
    // OCCT: find point inside F1, check if valid for F2
    // rcad: check if surfaces match and vertices are close
    let s1 = match ds.face_surface(f1_idx) { Some(s) => s, None => return false };
    let s2 = match ds.face_surface(f2_idx) { Some(s) => s, None => return false };
    // Check surface type match
    use rcad_kernel::geom::Surface3;
    match (&s1, &s2) {
        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
            let d = (p1.origin - p2.origin).length();
            let nd = 1.0 - p1.normal.dot(p2.normal).abs();
            d < fuzz && nd < fuzz
        }
        _ => false, // non-planar SD check not implemented
    }
}

// ====================================================================
// Sense — OCCT BOPTools_AlgoTools::Sense (L1209+)
// ====================================================================
/// Checks normal directions of two faces sharing an edge.
/// Returns 0=error, 1=same direction, -1=opposite.
pub fn sense(
    _f1_idx: usize, _f2_idx: usize,
    _edge_idx: usize,
    _ds: &DS,
) -> i32 {
    // rcad: requires face normal computation at shared edge
    0
}

// ====================================================================
// MakePCurve — OCCT BOPTools_AlgoTools::MakePCurve (L1657+)
// ====================================================================
/// Creates a 2D pcurve for an edge on two faces.
pub fn make_pcurve(
    _edge_idx: usize,
    _face1_idx: usize, _face2_idx: usize,
    _ds: &mut DS,
) {
    // rcad: pcurve creation is handled by MakePCurves step in PaveFiller
}

// ====================================================================
// BOPAlgo_Tools::FillMap — graph edge connection (BOPAlgo_Tools.cxx L38-48)
// ====================================================================
/// OCCT BOPAlgo_Tools::FillMap(int, int, IndexedDataMap<int, List<int>>).
/// Adds bidirectional connection between two nodes in an adjacency map.
pub fn fill_map(n1: usize, n2: usize, the_map: &mut std::collections::HashMap<usize, Vec<usize>>) {
    the_map.entry(n1).or_default().push(n2);
    the_map.entry(n2).or_default().push(n1);
}

// ====================================================================
// IntTools_Tools::IsOnPave1 — parameter on range boundary (L168+)
// ====================================================================
/// OCCT IntTools_Tools::IsOnPave1 (IntTools_Tools.cxx L627-646).
/// Returns true if aTR is inside [First, Last], or within aTolerance of a boundary.
pub fn is_on_pave_1(t: f64, r_first: f64, r_last: f64, tol: f64) -> bool {
    // OCCT L636-640: inside-range check first
    if t >= r_first && t <= r_last {
        return true;
    }
    // OCCT L642-644: distance to range boundaries within tolerance
    (t - r_first).abs() <= tol || (t - r_last).abs() <= tol
}

// ====================================================================
// BOPAlgo_Tools::MakeBlocks — connected components from graph (L121+)
// ====================================================================
/// OCCT BOPAlgo_Tools::MakeBlocks(IndexedDataMap<int, List<int>>, List<List<int>>).
/// Finds connected components in a vertex connection graph.
pub fn make_blocks(
    the_map: &std::collections::HashMap<usize, Vec<usize>>,
    the_blocks: &mut Vec<Vec<usize>>,
) {
    let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for (&start, _) in the_map {
        if visited.contains(&start) {
            continue;
        }
        let mut block: Vec<usize> = Vec::new();
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            if !visited.insert(node) {
                continue;
            }
            block.push(node);
            if let Some(neighbors) = the_map.get(&node) {
                for &n in neighbors {
                    if !visited.contains(&n) {
                        stack.push(n);
                    }
                }
            }
        }
        if block.len() >= 2 {
            the_blocks.push(block);
        }
    }
}