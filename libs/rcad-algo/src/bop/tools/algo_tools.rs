// OCCT BOPTools_AlgoTools — utility functions for boolean operations.
//
// OCCT BOPTools_AlgoTools.cxx / _1.cxx / _2.cxx
// Functions translated 1:1 from OCCT.

use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::int_tools::face_make_curve::intermediate_point;
use crate::topalgo::brep_top_adaptor::fclass2d::FClass2d;
use crate::topalgo::shape_source::ShapeSource;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use std::sync::Arc;

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
/// OCCT BOPTools_AlgoTools::MakeNewVertex(aP1, aP2, aTol1, aTol2, aVnew):
/// midpoint of the two points; tolerance = max(aTol1, aTol2, 0.5 * aDist).
pub fn make_new_vertex(v1_pnt: glam::DVec3, v1_tol: f64,
                       v2_pnt: glam::DVec3, v2_tol: f64) -> (glam::DVec3, f64) {
    let mid = (v1_pnt + v2_pnt) * 0.5;
    let a_d = v1_pnt.distance(v2_pnt);
    let mut tol = v1_tol.max(v2_tol);
    if tol < a_d * 0.5 {
        tol = a_d * 0.5;
    }
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
// CorrectRange — OCCT BOPTools_AlgoTools::CorrectRange (EF variant, L364-434)
// ====================================================================
/// Corrects the edge range for edge-face intersection.
/// OCCT (BOPTools_AlgoTools_2.cxx L364-434):
///   aTolF = BRep_Tool::Tolerance(aF); aRes = aTolF;
///   spline: aRes /= |D1| at the endpoint (or aBC.Resolution(aRes));
///   other:  aRes = aBC.Resolution(aRes);
///   first += aRes; last -= aRes; reset to aSR if the range < PConfusion.
/// (Note: unlike the EE variant, this does NOT early-return for lines.)
pub fn correct_range_ef(
    curve: &rcad_kernel::geom::Curve3,
    t1: f64, t2: f64,
    _tol_e: f64, tol_f: f64,
) -> (f64, f64) {
    let d_t = rcad_kernel::PCONFUSION;
    let mut nt1 = t1;
    let mut nt2 = t2;
    for i in 0..2 {
        let par = if i == 0 { t1 } else { t2 };
        let mut a_res = tol_f;
        let a_mgn = curve.tangent_at(par).length();
        if a_mgn > 1e-12 {
            a_res = a_res / a_mgn;
        } else {
            a_res = curve_resolution_ef(curve, par, a_res);
        }
        if i == 0 {
            nt1 = t1 + a_res;
        } else {
            nt2 = t2 - a_res;
        }
        if (nt2 - nt1) < d_t {
            nt1 = t1;
            nt2 = t2;
        }
    }
    (nt1, nt2)
}

/// OCCT Adaptor3d_Curve::Resolution(tol) = tol / |dP/dt| (with tol fallback
/// when the tangent speed is degenerate).
fn curve_resolution_ef(curve: &rcad_kernel::geom::Curve3, t: f64, tol: f64) -> f64 {
    let speed = curve.tangent_at(t).length();
    if speed < 1e-15 {
        tol
    } else {
        tol / speed
    }
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
// AreFacesSameDomain — OCCT BOPTools_AlgoTools::AreFacesSameDomain
// (BOPTools_AlgoTools.cxx L1139-1205)
// ====================================================================
/// Checks if two faces are geometrically coincident (same domain).
///
/// OCCT finds a point inside the first face (BOPTools_AlgoTools3D::PointInFace)
/// and checks its validity for the second face (IntTools_Context::
/// IsValidPointForFace). If valid — the faces are same domain.
pub fn are_faces_same_domain(
    f1: usize,
    f2: usize,
    context: &mut IntToolsContext,
    ds: &DS,
    fuzz: f64,
) -> bool {
    // OCCT L1149-1155: the idea is to find a point inside the first face and
    // check its validity for the second face. If valid — the faces are SD.
    let mut b_faces_sd = false;
    // OCCT L1162-1168: find point inside the first face.
    let (i_err, a_p1, _a_p2d1) = point_in_face(f1, ds);
    if i_err != 0 {
        // OCCT L1169-1170: unable to find the point.
        return b_faces_sd;
    }
    // OCCT L1176-1179: compute the tolerance of the faces, taking into account
    // the deviation of the edges from the surfaces.
    let mut a_tol_f1 = ds.face_tolerance(f1);
    let mut a_tol_f2 = ds.face_tolerance(f2);
    // OCCT L1181-1198: find maximal tolerance of edges. The faces should have
    // the same boundaries, so it does not matter which face to explore.
    let mut a_tol_e_max = -1.0;
    for a_e in face_edge_indices(ds, f1) {
        if !ds.is_edge_degenerated(a_e) {
            let a_tol_e = ds.edge_tolerance(a_e);
            if a_tol_e > a_tol_e_max {
                a_tol_e_max = a_tol_e;
            }
        }
    }
    if a_tol_e_max > a_tol_f1 {
        a_tol_f1 = a_tol_e_max;
    }
    if a_tol_e_max > a_tol_f2 {
        a_tol_f2 = a_tol_e_max;
    }
    // OCCT L1201-1203: checking criteria.
    let a_tol = a_tol_f1 + a_tol_f2 + fuzz.max(rcad_kernel::CONFUSION);
    // OCCT L1202: project and classify the point on the second face.
    b_faces_sd = context.is_valid_point_for_face(a_p1, f2, ds, a_tol);
    b_faces_sd
}

// ====================================================================
// PointInFace — OCCT BOPTools_AlgoTools3D::PointInFace
// (BOPTools_AlgoTools3D.cxx L906-938)
// ====================================================================
/// Finds a point inside the face `f` (OCCT L906-938 wrapper).
///
/// Returns `(iErr, theP, theP2D)`: `iErr == 0` on success, with `theP` the 3D
/// point and `theP2D` its UV parameters on the face surface.
fn point_in_face(f: usize, ds: &DS) -> (i32, glam::DVec3, glam::DVec2) {
    // OCCT L918-919: theContext->UVBounds(theF, aUMin, aUMax, aVMin, aVMax).
    // rcad: BRepTools::UVBounds equivalent — the boundary-sampled rect.
    let [a_umin, a_umax, a_vmin, a_vmax] = ds.face_actual_uv_bounds(f);
    let _ = (a_vmin, a_vmax); // OCCT: V bounds are computed but unused.
    // OCCT L921-924: aD2D = gp_Dir2d::D::Y; aUx = IntermediatePoint(aUMin, aUMax).
    let mut a_ux = intermediate_point(a_umin, a_umax);
    // OCCT L926-936: two attempts, the second with a translated (mirrored) line.
    let mut i_err = 1;
    let mut a_p = glam::DVec3::ZERO;
    let mut a_p2d = glam::DVec2::ZERO;
    for _ in 0..2 {
        let (err, p, p2d) = point_in_face_line(f, a_ux, ds);
        i_err = err;
        if i_err == 0 {
            a_p = p;
            a_p2d = p2d;
            break;
        }
        // OCCT L931-934: possible reason — incorrect computation of the 2d box
        // of the face. Try the translated (mirrored) line.
        a_ux = a_umax - (a_ux - a_umin);
    }
    (i_err, a_p, a_p2d)
}

// ====================================================================
// PointInFace (line) — OCCT BOPTools_AlgoTools3D::PointInFace
// (BOPTools_AlgoTools3D.cxx L942-988)
// ====================================================================
/// Finds a point inside the face `f` on the vertical 2D line U = aUx.
///
/// OCCT trims the line with the face boundary via the per-face
/// Geom2dHatch_Hatcher (theContext->Hatcher(theF), IntTools_Context.cxx
/// L343-394) and takes the middle of the first inside domain. rcad has no
/// Geom2dHatch subsystem: the equivalent intersects the line with the face's
/// UV boundary polygon (FClass2d sampling of the boundary pcurves) and applies
/// the same parity rule (the line starts OUTSIDE the bounded face at v = -inf).
fn point_in_face_line(f: usize, a_ux: f64, ds: &DS) -> (i32, glam::DVec3, glam::DVec2) {
    // OCCT L981-1011: trim + compute domains; aNbDomains == 0 → iErr = 2.
    let Some(&(a_v1, a_v2)) = hatch_line_intervals(f, a_ux, ds).first() else {
        return (2, glam::DVec3::ZERO, glam::DVec2::ZERO);
    };
    // OCCT L1023-1025: aVx = IntermediatePoint(aV1, aV2) (theDt2D is 0).
    let a_vx = intermediate_point(a_v1, a_v2);
    // OCCT L1027-1028: theL2D->D0(aVx, theP2D); aS->D0(theP2D, theP).
    let a_p2d = glam::DVec2::new(a_ux, a_vx);
    let Some(a_s) = ds.face_surface(f) else {
        return (1, glam::DVec3::ZERO, a_p2d);
    };
    let a_p = a_s.point_at(a_p2d.x, a_p2d.y);
    (0, a_p, a_p2d)
}

// ====================================================================
// Hatch line intervals — rcad equivalent of Geom2dHatch_Hatcher domains
// ====================================================================
/// rcad equivalent of the Geom2dHatch_Hatcher domain computation
/// (Geom2dHatch_Hatcher.cxx Trim L361-1019 + ComputeDomains L1144-1810) for a
/// vertical 2D line U = aUx: returns the inside intervals [v1, v2] of the face
/// crossed by the line.
///
/// The boundary is the FClass2d UV polygon (the face's sampled boundary
/// pcurves). The line starts OUTSIDE the bounded face at v = -inf, so the
/// inside intervals are the pairs of sorted crossing v-values (even-odd rule).
/// The half-open rule `lo.x <= aUx < hi.x` counts each vertex crossing exactly
/// once; degenerate zero-width intervals (tangency at a U-extremum vertex,
/// which the hatcher merges by its 2D confusion aTolHatch2D = 1e-8) are dropped.
fn hatch_line_intervals(f: usize, a_ux: f64, ds: &DS) -> Vec<(f64, f64)> {
    // OCCT IntTools_Context::FClass2d tolerance (face tolerance floored at
    // CONFUSION) — the same classifier tolerance used by the pipeline.
    let fclass = FClass2d::new(ds, f, ds.face_tolerance(f).max(rcad_kernel::CONFUSION));
    let mut vs: Vec<f64> = Vec::new();
    for poly in fclass.uv_polygons() {
        let n = poly.len();
        for i in 0..n {
            let a = poly[i];
            let b = poly[(i + 1) % n];
            let (lo, hi) = if a.x <= b.x { (a, b) } else { (b, a) };
            if lo.x <= a_ux && a_ux < hi.x {
                // v-coordinate of the crossing (the line is (aUx, t)).
                let t = (a_ux - lo.x) / (hi.x - lo.x);
                vs.push(lo.y + t * (hi.y - lo.y));
            }
        }
    }
    if vs.len() < 2 {
        return Vec::new();
    }
    vs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    const HATCH_CONFUSION: f64 = 1e-8; // OCCT aTolHatch2D (IntTools_Context.cxx L350)
    let mut intervals: Vec<(f64, f64)> = Vec::new();
    let mut k = 0;
    while k + 1 < vs.len() {
        let (v1, v2) = (vs[k], vs[k + 1]);
        if v2 - v1 > HATCH_CONFUSION {
            intervals.push((v1, v2));
        }
        k += 2;
    }
    intervals
}

/// DS edge indices of the face's boundary (OCCT TopExp_Explorer(theF,
/// TopAbs_EDGE) in AreFacesSameDomain L1187).
fn face_edge_indices(ds: &DS, f: usize) -> Vec<usize> {
    let mut edges: Vec<usize> = Vec::new();
    let fshape = ds.shape_at(f);
    let face_data = match &*fshape.data {
        TShape::Face(fd) => fd,
        _ => return edges,
    };
    let wires: Vec<Shape> = std::iter::once(face_data.outer_wire.clone())
        .chain(face_data.inner_wires.iter().cloned())
        .collect();
    for ws in wires {
        let Some(wi) = ds.map_shape_index(ws.ptr_id(), ws.location) else { continue };
        if wi >= ds.nb_shapes() || ds.shape_type(wi) != ShapeType::Wire {
            continue;
        }
        let wshape = ds.shape_at(wi);
        let wire_edges = match &*wshape.data {
            TShape::Wire(w) => w.edges.clone(),
            _ => Vec::new(),
        };
        for eshape in wire_edges {
            if let Some(ei) = ds.map_shape_index(eshape.ptr_id(), eshape.location) {
                edges.push(ei);
            }
        }
    }
    edges
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
// IsSplitToReverse — OCCT BOPTools_AlgoTools::IsSplitToReverse (Face, L1324-1436)
// ====================================================================
/// OCCT BOPTools_AlgoTools::IsSplitToReverse(theFSp, theFSr, theContext, theError).
/// Determines whether the split face `theFSp` must be reversed to be oriented
/// consistently with the original face `theFSr`.
///
/// Returns `(bToReverse, errCode)`; `errCode == 0` means the check succeeded.
///
/// Fast path (OCCT L1336-1341): when the two faces share the same surface
/// handle, only their orientations are compared. Otherwise a point inside the
/// split face is found, the surface normals at that point are computed on both
/// surfaces, and the result is whether the normals point in opposite directions
/// (OCCT L1383-1435). The point is the centroid of the face's outer-wire vertex
/// endpoints projected onto the surface (rcad semantic equivalent of
/// BOPTools_AlgoTools3D::PointInFace).
pub fn is_split_to_reverse_face(f_sp: &Shape, f_sr: &Shape) -> (bool, i32) {
    // OCCT L1336-1341: same surface handle -> compare orientations only.
    // rcad: surfaces are value types, so "same surface" means geometric identity
    // of the analytic surface parameters.
    let surf_sp = f_sp.as_face().and_then(|fd| fd.surface.clone());
    let surf_or = f_sr.as_face().and_then(|fd| fd.surface.clone());
    if let (Some(s1), Some(s2)) = (&surf_sp, &surf_or) {
        if surface_same(s1, s2) {
            return (f_sp.orientation != f_sr.orientation, 0);
        }
    }
    // OCCT L1344-1380: find a point inside the split face.
    let p3d = face_reference_point(f_sp);
    let (uv_sp, _) = match surf_sp.as_ref() {
        Some(s) => crate::bop::closest_point_on_surface(s, p3d),
        None => return (false, 1),
    };
    // OCCT L1383-1392: normal direction of the split face at the point.
    let mut dn_sp = match surf_sp.as_ref() {
        Some(s) => s.normal_at(uv_sp.x, uv_sp.y),
        None => return (false, 2),
    };
    if f_sp.orientation == Orientation::Reversed {
        dn_sp = -dn_sp;
    }
    // OCCT L1401-1414: project the point from the split face on the original face.
    let (uv_or, _) = match surf_or.as_ref() {
        Some(s) => crate::bop::closest_point_on_surface(s, p3d),
        None => return (false, 3),
    };
    // OCCT L1418-1431: normal direction for the original face in this point.
    let mut dn_or = match surf_or.as_ref() {
        Some(s) => s.normal_at(uv_or.x, uv_or.y),
        None => return (false, 4),
    };
    if f_sr.orientation == Orientation::Reversed {
        dn_or = -dn_or;
    }
    // OCCT L1434-1435: compare the normals.
    let a_cos = dn_sp.dot(dn_or);
    (a_cos < 0.0, 0)
}

/// OCCT BOPTools_AlgoTools::IsSplitToReverse(edge) (BOPTools_AlgoTools.cxx
/// L1456-1531) — true when the split edge a_sp has the opposite direction to
/// the original edge a_e. rcad: straight/analytic edge translation.
pub fn is_split_to_reverse_edge(a_sp: &Shape, a_e: &Shape) -> (bool, i32) {
    use rcad_kernel::topods::TShape;
    // OCCT L1461-1468: degenerated edges are not processed.
    let sp_degen = a_sp
        .as_edge()
        .map(|ed| ed.degenerated)
        .unwrap_or(true);
    let e_degen = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(true);
    if sp_degen || e_degen {
        return (false, 1);
    }
    // OCCT L1472-1477: same curve handle -> compare orientations only.
    let c_sp = a_sp.as_edge().and_then(|ed| ed.curve.clone());
    let c_or = a_e.as_edge().and_then(|ed| ed.curve.clone());
    if let (Some(c1), Some(c2)) = (&c_sp, &c_or) {
        if curves_same(c1, c2) {
            return (a_sp.orientation != a_e.orientation, 0);
        }
    }
    // OCCT L1480-1531: compare the tangent vectors at a sample point.
    let t_sp = edge_tangent_3d(a_sp);
    let t_or = edge_tangent_3d(a_e);
    match (t_sp, t_or) {
        (Some(vs), Some(vo)) => (vs.dot(vo) < 0.0, 0),
        _ => (false, 4),
    }
}

/// Edge tangent at the parameterization direction (straight edges).
/// OCCT BOPTools_AlgoTools2D::EdgeTangent — the curve derivative, reversed for
/// a REVERSED edge.
pub fn edge_tangent_3d(e: &Shape) -> Option<glam::DVec3> {
    use rcad_kernel::topods::TShape;
    match &*e.data {
        TShape::Edge(ed) => {
            let p1 = match &*ed.first.data {
                TShape::Vertex(vd) => vd.point,
                _ => return None,
            };
            let p2 = match &*ed.last.data {
                TShape::Vertex(vd) => vd.point,
                _ => return None,
            };
            let d = p2 - p1;
            if d.length_squared() < 1e-24 {
                return None;
            }
            let mut t = d.normalize();
            if e.orientation == rcad_kernel::topods::Orientation::Reversed {
                t = -t;
            }
            Some(t)
        }
        _ => None,
    }
}

/// Geometric identity of two 3D curves (rcad equivalent of OCCT's Geom_Curve
/// handle equality in IsSplitToReverse L1474).
fn curves_same(a: &rcad_kernel::geom::Curve3, b: &rcad_kernel::geom::Curve3) -> bool {
    use rcad_kernel::geom::Curve3;
    const TOL: f64 = 1e-9;
    match (a, b) {
        (Curve3::Line(l1), Curve3::Line(l2)) => {
            (l1.origin - l2.origin).length() < TOL && l1.direction.dot(l2.direction) > 1.0 - TOL
        }
        (Curve3::Circle(c1), Curve3::Circle(c2)) => {
            (c1.center - c2.center).length() < TOL
                && c1.normal.dot(c2.normal) > 1.0 - TOL
                && (c1.radius - c2.radius).abs() < TOL
        }
        (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
            (e1.center - e2.center).length() < TOL
                && e1.normal.dot(e2.normal) > 1.0 - TOL
                && (e1.major_radius - e2.major_radius).abs() < TOL
                && (e1.minor_radius - e2.minor_radius).abs() < TOL
        }
        _ => false,
    }
}

/// Geometric identity of two analytic surfaces (rcad equivalent of OCCT's
/// Geom_Surface handle equality in IsSplitToReverse L1338). Same direction for
/// the axis/normal; the reference directions (u_dir/v_dir) are not compared —
/// the parameterization orientation is irrelevant for the orientation decision.
fn surface_same(a: &rcad_kernel::geom::Surface3, b: &rcad_kernel::geom::Surface3) -> bool {
    use rcad_kernel::geom::Surface3;
    const TOL: f64 = 1e-9;
    match (a, b) {
        (Surface3::Plane(p1), Surface3::Plane(p2)) => {
            (p1.origin - p2.origin).length() < TOL
                && p1.normal.dot(p2.normal) > 1.0 - TOL
        }
        (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
            (c1.origin - c2.origin).length() < TOL
                && c1.axis.dot(c2.axis) > 1.0 - TOL
                && (c1.radius - c2.radius).abs() < TOL
        }
        (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
            (s1.center - s2.center).length() < TOL && (s1.radius - s2.radius).abs() < TOL
        }
        (Surface3::Cone(c1), Surface3::Cone(c2)) => {
            (c1.apex - c2.apex).length() < TOL
                && c1.axis.dot(c2.axis) > 1.0 - TOL
                && (c1.radius - c2.radius).abs() < TOL
                && (c1.half_angle_rad - c2.half_angle_rad).abs() < TOL
        }
        (Surface3::Torus(t1), Surface3::Torus(t2)) => {
            (t1.center - t2.center).length() < TOL
                && t1.axis.dot(t2.axis) > 1.0 - TOL
                && (t1.major_radius - t2.major_radius).abs() < TOL
                && (t1.minor_radius - t2.minor_radius).abs() < TOL
        }
        _ => false,
    }
}

/// Reference point of a face: centroid of the outer-wire edge endpoint vertices.
fn face_reference_point(f: &Shape) -> glam::DVec3 {
    match f.as_face() {
        Some(fd) => {
            let mut pts: Vec<glam::DVec3> = Vec::new();
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                for e in &wd.edges {
                    if let TShape::Edge(ed) = &*e.data {
                        if let TShape::Vertex(vd) = &*ed.first.data {
                            pts.push(vd.point);
                        }
                        if let TShape::Vertex(vd) = &*ed.last.data {
                            pts.push(vd.point);
                        }
                    }
                }
            }
            if pts.is_empty() {
                glam::DVec3::ZERO
            } else {
                pts.iter().sum::<glam::DVec3>() / pts.len() as f64
            }
        }
        None => glam::DVec3::ZERO,
    }
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

/// OCCT IntTools_Tools::IsInRange (IntTools_Tools.cxx L650-666).
/// Returns true if either endpoint of the range (r2_first, r2_last) lies
/// within the reference range (r1_first, r1_last) expanded by aTol:
///   aTRef1 -= tol; aTRef2 += tol;
///   bIsIn = (aT1 >= aTRef1 && aT1 <= aTRef2) || (aT2 >= aTRef1 && aT2 <= aTRef2);
pub fn is_in_range(r1_first: f64, r1_last: f64, r2_first: f64, r2_last: f64, tol: f64) -> bool {
    let t_ref1 = r1_first - tol;
    let t_ref2 = r1_last + tol;
    (r2_first >= t_ref1 && r2_first <= t_ref2) || (r2_last >= t_ref1 && r2_last <= t_ref2)
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