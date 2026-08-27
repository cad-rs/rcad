// OCCT BOPTools_AlgoTools — utility functions for boolean operations.
//
// OCCT BOPTools_AlgoTools.cxx / _1.cxx / _2.cxx
// Functions translated 1:1 from OCCT.

use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::int_tools::face_make_curve::intermediate_point;
use crate::topalgo::brep_top_adaptor::fclass2d::{FClass2d, State};
use crate::topalgo::shape_source::ShapeSource;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape};
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
    // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the pcurve key
    // location is L.Predivided(E.Location()).
    let (fid, floc) = ds.face_key(face_idx)?;
    let fkey = (
        fid,
        crate::bop::algo::compose_face_edge_pcurve_location(floc, shape.location, &ds.locations),
    );
    let pcurve_info = shape.as_edge()?.pcurves.get(&fkey)?.clone();
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
// Shape-based AreFacesSameDomain — OCCT BOPTools_AlgoTools::AreFacesSameDomain
// (BOPTools_AlgoTools.cxx L1139-1205) operates on TopoDS_Face directly, so it
// also classifies split face images that are not registered in the DS.
// rcad's DS-index version above cannot address such pieces; these functions
// reproduce the same checks on the Shape geometry (face tolerance + max edge
// tolerance, PointInFace hatcher point, projection + FClass2d on the piece).
// ====================================================================

/// UV boundary polygons of a face Shape (outer + inner wires), built from the
/// wire edge pcurves keyed by the face (the OCCT Geom2dHatch_Hatcher boundary),
/// falling back to projecting the wire edge endpoint vertices onto the face
/// surface when a pcurve is missing — the piece-equivalent of the DS FClass2d
/// UV polygons. The polygon is returned in wire-edge traversal order.
const PCURVE_SAMPLES: usize = 16;
fn shape_uv_polygons(face: &Shape, locations: &[glam::DAffine3]) -> Option<Vec<Vec<glam::DVec2>>> {
    let fd = face.as_face()?;
    let surf = fd.surface.as_ref()?;
    let mut polys: Vec<Vec<glam::DVec2>> = Vec::new();
    let mut wire_poly = |w: &Shape, out: &mut Vec<glam::DVec2>| {
        if let TShape::Wire(wd) = &*w.data {
            // OCCT Geom2dHatch_Hatcher intersects each edge pcurve once, in the
            // wire's traversal order — the polygon sampled from the edge
            // endpoints must follow the same chain. The WireSplitter loops are
            // valid cycles, but their stored edge order starts at an arbitrary
            // joint (G1 A-B-P loop: [P->A, B->P, A->B]), so the raw
            // first+last sampling walks P->A->B->P->A->B and duplicates every
            // crossing in the parity hatch (no interval found). Reorder the
            // edges into the vertex chain first (traversal endpoints from the
            // composed wire x edge orientation).
            let w_or = w.orientation;
            let edges = wd.edges.clone();
            let trav = |e: &Shape| -> Option<((u64, u32), (u64, u32))> {
                let ed = e.as_edge()?;
                let eff = Orientation::compose(w_or, e.orientation);
                let (vf, vl) = ((ed.first.ptr_id(), ed.first.location), (ed.last.ptr_id(), ed.last.location));
                Some(match eff {
                    Orientation::Reversed => (vl, vf),
                    _ => (vf, vl),
                })
            };
            let mut order: Vec<(usize, bool)> = Vec::with_capacity(edges.len()); // (edge index, reversed)
            {
                let mut used = vec![false; edges.len()];
                let mut cur: Option<(u64, u32)> = None;
                while order.len() < edges.len() {
                    let mut progressed = false;
                    for (i, e) in edges.iter().enumerate() {
                        if used[i] {
                            continue;
                        }
                        let Some((s, t)) = trav(e) else { continue };
                        if cur.is_none() {
                            order.push((i, false));
                            used[i] = true;
                            cur = Some(t);
                            progressed = true;
                            break;
                        }
                        if Some(s) == cur {
                            order.push((i, false));
                            used[i] = true;
                            cur = Some(t);
                            progressed = true;
                            break;
                        }
                        if Some(t) == cur {
                            // The edge's traversal start is its END — the edge
                            // participates in the loop reversed.
                            order.push((i, true));
                            used[i] = true;
                            cur = Some(s);
                            progressed = true;
                            break;
                        }
                    }
                    if !progressed {
                        break;
                    }
                }
                if order.len() < edges.len() {
                    order = edges.iter().enumerate().map(|(i, _)| (i, false)).collect();
                }
            }
            // OCCT Geom2dHatch_Hatcher hatches the edge PCURVES of the face
            // (BOPTools_AlgoTools3D.cxx L992-1066), never 3D vertex projections.
            // When every edge of the wire carries a pcurve keyed by this face,
            // the polygon is sampled from those pcurves in the traversal chain
            // order; otherwise the 3D vertex projections below are used.  A
            // u-periodic band (cylinder/cone side) keeps the band_rect fallback
            // instead: its wire walks [cap circle, seam, cap circle, seam] and
            // the seam junction points (u=0 and u=2*PI of the same vertex) are
            // never deduped, so the pcurve polygon retraces and the parity
            // hatch finds no interval (bcut_simple l8: coaxial cylinders).
            let is_u_band = matches!(
                surf,
                rcad_kernel::geom::Surface3::Cylinder(_) | rcad_kernel::geom::Surface3::Cone(_)
            );
            // OCCT BRep_Tool.cxx L345: pcurve-key location is the composed
            // transform VALUE; rcad stores the stable value hash of the
            // face's own transform here.
            let tr_f = if face.location == 0 {
                glam::DAffine3::IDENTITY
            } else {
                locations
                    .get((face.location - 1) as usize)
                    .copied()
                    .unwrap_or(glam::DAffine3::IDENTITY)
            };
            let mykey = (
                face.ptr_id(),
                rcad_kernel::topo::topods::pcurve_location_id(&tr_f),
            );
            let mut ppts: Vec<glam::DVec2> = Vec::new();
            let mut pcurve_ok = true;
            for &(ei, rev) in &order {
                let e = &edges[ei];
                let ed = match &*e.data {
                    TShape::Edge(ed) => ed,
                    _ => {
                        pcurve_ok = false;
                        break;
                    }
                };
                // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L354-361): a
                // closed surface seam edge (CurveOnClosedSurface) returns the
                // second pcurve for a REVERSED edge and the first otherwise —
                // the two wire instances of the seam map to u=2*PI and u=0.
                let (pc, t1, t2) = if let Some((pc1, pc2, range)) = ed.representations.iter().find_map(|r| match r {
                    rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                        face: f,
                        pcurve1,
                        pcurve2,
                        range,
                    } if *f == mykey => Some((pcurve1.clone(), pcurve2.clone(), *range)),
                    _ => None,
                }) {
                    let pc = if e.orientation == rcad_kernel::topods::Orientation::Reversed {
                        pc2
                    } else {
                        pc1
                    };
                    (pc, range[0], range[1])
                } else if let Some(v) = ed.pcurves.get(&mykey) {
                    v.clone()
                } else {
                    pcurve_ok = false;
                    break;
                };
                // Traversal direction: FORWARD -> [t1, t2], REVERSED -> [t2, t1]
                // (TopoDS_Iterator cumOri), then flipped for a loop-reversed
                // edge (its traversal end joins the current vertex).
                let (ta, tb) = match Orientation::compose(w_or, e.orientation) {
                    Orientation::Reversed => (t2, t1),
                    _ => (t1, t2),
                };
                let (ta, tb) = if rev { (tb, ta) } else { (ta, tb) };
                for i in 0..=PCURVE_SAMPLES {
                    let t = ta + (tb - ta) * (i as f64) / (PCURVE_SAMPLES as f64);
                    let uv = rcad_kernel::geom::Curve2dEval::point_at(&pc, t);
                    if ppts.last().map(|q| (*q - uv).length() < 1e-9).unwrap_or(false) {
                        continue;
                    }
                    ppts.push(uv);
                }
            }
            if pcurve_ok && !ppts.is_empty() {
                *out = ppts;
                return;
            }
            // Fall back: sample the two traversal endpoints (start, end) in
            // order: FORWARD -> [first, last], REVERSED -> [last, first]
            // (TopoDS_Iterator cumOri semantics).
            for &(ei, rev) in &order {
                let e = &edges[ei];
                let (_va, _vb) = match trav(e) {
                    Some(v) => v,
                    None => continue,
                };
                let (p0, p1) = match &*e.data {
                    TShape::Edge(ed) => {
                        let loc0 = if e.location != 0 { e.location } else { ed.first.location };
                        let loc1 = if e.location != 0 { e.location } else { ed.last.location };
                        let pt = |vd: &rcad_kernel::topods::TVertexData, loc: u32| -> glam::DVec3 {
                            if loc == 0 {
                                vd.point
                            } else {
                                locations
                                    .get(loc as usize)
                                    .copied()
                                    .unwrap_or(glam::DAffine3::IDENTITY)
                                    .transform_point3(vd.point)
                            }
                        };
                        let a = match &*ed.first.data {
                            TShape::Vertex(vd) => pt(vd, loc0),
                            _ => glam::DVec3::ZERO,
                        };
                        let b = match &*ed.last.data {
                            TShape::Vertex(vd) => pt(vd, loc1),
                            _ => glam::DVec3::ZERO,
                        };
                        match Orientation::compose(w_or, e.orientation) {
                            Orientation::Reversed => (b, a),
                            _ => (a, b),
                        }
                    }
                    _ => (glam::DVec3::ZERO, glam::DVec3::ZERO),
                };
                // A loop-reversed edge (its traversal end joins the current
                // vertex) participates in the chain in the opposite direction.
                let (p0, p1) = if rev { (p1, p0) } else { (p0, p1) };
                for p in [p0, p1] {
                    let (uv, _) = crate::bop::closest_point_on_surface(surf, p);
                    // The wire's consecutive edges share their joint vertex, so
                    // the raw first+last sampling repeats it (a triangle
                    // becomes [P,A,A,B,B,P]). The duplicated point yields
                    // coincident parity crossings in the hatch
                    // (PointInFaceShape) and must be dropped — the boundary
                    // polygon is the set of distinct joint points (OCCT
                    // Geom2dHatch intersects each pcurve once, so no crossing
                    // is ever duplicated).
                    if out.last().map(|q| (*q - uv).length() < 1e-9).unwrap_or(false) {
                        continue;
                    }
                    out.push(uv);
                }
            }
        }
    };
    let mut outer: Vec<glam::DVec2> = Vec::new();
    wire_poly(&fd.outer_wire, &mut outer);
    // A closed edge (circle) stores first == last, so the vertex sampling
    // yields a degenerate polygon (all joint points collapse onto the seam
    // parameter u=0 for a full-revolution cylindrical face). Sample the edge
    // curves instead (the equivalent of OCCT BRepTools::UVBounds on the
    // pcurves). The check must also catch the zero-area case: a polygon whose
    // points all share the same u (or are otherwise collinear) hatches
    // nothing and would make PointInFaceShape fail on a valid cylindrical
    // piece (bcommon_simple J1: coaxial cylinder common).
    let poly_usable = |poly: &Vec<glam::DVec2>| -> bool {
        if poly.len() < 3 {
            return false;
        }
        // At least 3 distinct corner points. A wire bounded by curved edges
        // (e.g. a disk made of two arcs) sampled at its vertices only
        // collapses onto a line (the first and last corners coincide), and
        // the hatch below finds no interval — such wires must fall back to
        // the curve sampling below. OCCT hatches the edge pcurves, whose
        // curved edges contribute real crossings, so a vertex-only polygon
        // is never used there.
        let mut n_distinct = 0usize;
        for p in poly {
            if !poly[..n_distinct].iter().any(|q| (*q - *p).length() < 1e-9) {
                n_distinct += 1;
            }
        }
        if n_distinct < 3 {
            return false;
        }
        let u0 = poly[0].x;
        !poly.iter().all(|p| (p.x - u0).abs() < 1e-9)
    };
    // Sample the edge curves of a u-periodic band face (cylinder/cone side)
    // into a single-axis rectangle: the wire walks [cap circle, seam, cap
    // circle, seam], whose projected polygon is a self-intersecting figure-8
    // that the periodic ray test cannot classify points between the caps as
    // inside. A band's UV domain is the rectangle [umin,umax] x [vmin,vmax];
    // OCCT's BRepClass_FaceExplorer reaches the same domain through the
    // natural boundaries + pcurves (IntTools_Context::IsValidPointForFace).
    let band_rect = |w: &Shape| -> Vec<glam::DVec2> {
        let pts = sample_wire_uv_curve(w, surf, locations);
        if pts.len() < 2 {
            return pts;
        }
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;
        for p in &pts {
            umin = umin.min(p.x);
            umax = umax.max(p.x);
            vmin = vmin.min(p.y);
            vmax = vmax.max(p.y);
        }
        if !(umin < umax && vmin < vmax) {
            return pts;
        }
        // Bring u into one period: the caps wrap, so umin/umax span a full
        // revolution; normalize to [0, 2*PI).
        let tau = std::f64::consts::TAU;
        if umax - umin > tau * 0.9 {
            umin = 0.0;
            umax = tau;
        }
        vec![
            glam::DVec2::new(umin, vmin),
            glam::DVec2::new(umax, vmin),
            glam::DVec2::new(umax, vmax),
            glam::DVec2::new(umin, vmax),
        ]
    };
    let is_band = matches!(
        surf,
        rcad_kernel::geom::Surface3::Cylinder(_) | rcad_kernel::geom::Surface3::Cone(_)
    );
    if !poly_usable(&outer) {
        if is_band {
            outer = band_rect(&fd.outer_wire);
        } else {
            outer = sample_wire_uv_curve(&fd.outer_wire, surf, locations);
        }
    }
    if poly_usable(&outer) {
        polys.push(outer);
    }
    for iw in &fd.inner_wires {
        let mut inner: Vec<glam::DVec2> = Vec::new();
        wire_poly(iw, &mut inner);
        if !poly_usable(&inner) {
            if is_band {
                inner = band_rect(iw);
            } else {
                inner = sample_wire_uv_curve(iw, surf, locations);
            }
        }
        if poly_usable(&inner) {
            polys.push(inner);
        }
    }
    Some(polys)
}

/// Sample the wire's edge 3D curves (with the composed orientation and the
/// edge Location applied) projected onto the face's surface — the fallback UV
/// polygon for wires whose vertex sampling is degenerate (closed edges with
/// first == last, e.g. a circular boundary).
fn sample_wire_uv_curve(
    w: &Shape,
    surf: &rcad_kernel::geom::Surface3,
    locations: &[glam::DAffine3],
) -> Vec<glam::DVec2> {
    let mut out: Vec<glam::DVec2> = Vec::new();
    const N: usize = 16;
    let (w_or, edges) = match &*w.data {
        TShape::Wire(wd) => (w.orientation, wd.edges.clone()),
        _ => return out,
    };
    for e in &edges {
        let ed = match &*e.data {
            TShape::Edge(ed) => ed,
            _ => continue,
        };
        let Some(curve) = &ed.curve else { continue };
        let (t1, t2) = (ed.range[0], ed.range[1]);
        // OCCT TopoDS_Iterator cumOri: the wire orientation composes into the
        // edge; the traversal runs first->last (FWD) or last->first (REV).
        let eff = Orientation::compose(w_or, e.orientation);
        let (ta, tb) = match eff {
            Orientation::Reversed => (t2, t1),
            _ => (t1, t2),
        };
        let loc = e.location;
        for i in 0..=N {
            let t = ta + (tb - ta) * (i as f64) / (N as f64);
            let mut p = curve.point_at(t);
            if loc != 0 {
                if let Some(tr) = locations.get(loc as usize) {
                    p = tr.transform_point3(p);
                }
            }
            let (uv, _) = crate::bop::closest_point_on_surface(surf, p);
            out.push(uv);
        }
    }
    out
}

/// Periodic-ray classification of `uv` against a UV polygon on a U-periodic
/// surface: every edge is normalized so its u-span stays within one period
/// (seam-crossing edges are unwrapped), and the query u is brought into the
/// edge's period before the standard even-odd ray test. Matches the
/// BRepClass_FaceExplorer seam handling on periodic surfaces.
fn point_in_periodic_polygon(uv: glam::DVec2, poly: &[glam::DVec2], period: f64) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    let half = period * 0.5;
    let mut inside = false;
    for i in 0..n {
        let mut a = poly[i];
        let mut b = poly[(i + 1) % n];
        // Unwrap the edge so its u-span stays within one period anchored at
        // a.x (seam-crossing edges are normalized).
        while b.x - a.x > half {
            b.x -= period;
        }
        while b.x - a.x < -half {
            b.x += period;
        }
        let u = uv.x;
        if (a.y > uv.y) != (b.y > uv.y) {
            let x_int = a.x + (uv.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if u < x_int {
                inside = !inside;
            }
        }
    }
    inside
}

/// OCCT BOPTools_AlgoTools3D::PointInFace (BOPTools_AlgoTools3D.cxx
/// L906-938/L942-990) for a split face image: the hatcher point inside the
/// face, found on the vertical 2D line U = IntermediatePoint(UMin, UMax)
/// crossing the face's UV polygons (parity from V = -inf, same as
/// hatch_line_intervals). Returns (iErr, theP, theP2D) with iErr == 0 on
/// success.
pub(crate) fn point_in_face_shape(
    face: &Shape,
    locations: &[glam::DAffine3],
) -> (i32, glam::DVec3, glam::DVec2) {
    let Some(surf) = face.as_face().and_then(|fd| fd.surface.clone()) else {
        return (1, glam::DVec3::ZERO, glam::DVec2::ZERO);
    };
    let Some(polys) = shape_uv_polygons(face, locations) else {
        return (2, glam::DVec3::ZERO, glam::DVec2::ZERO);
    };
    // UV bounds over all polygons (OCCT UVBounds of the face).
    let mut a_umin = f64::INFINITY;
    let mut a_umax = f64::NEG_INFINITY;
    let mut a_vmin = f64::INFINITY;
    let mut a_vmax = f64::NEG_INFINITY;
    for poly in &polys {
        for p in poly {
            a_umin = a_umin.min(p.x);
            a_umax = a_umax.max(p.x);
            a_vmin = a_vmin.min(p.y);
            a_vmax = a_vmax.max(p.y);
        }
    }
    if !(a_umin <= a_umax && a_vmin <= a_vmax) {
        return (2, glam::DVec3::ZERO, glam::DVec2::ZERO);
    }
    let a_ux0 = intermediate_point(a_umin, a_umax);
    // OCCT L919-935: two attempts; the second uses a translated (mirrored)
    // line aUx = aUMax - (aUx - aUMin) in case the 2d box is wrong.
    let mut a_ux = a_ux0;
    let mut i_err = 2;
    let mut a_vx = 0.0;
    for _ in 0..2 {
        // Vertical line crossings (parity rule: OUT at V = -inf).
        let mut vs: Vec<f64> = Vec::new();
        for poly in &polys {
            let n = poly.len();
            for i in 0..n {
                let a = poly[i];
                let b = poly[(i + 1) % n];
                let (lo, hi) = if a.x <= b.x { (a, b) } else { (b, a) };
                if lo.x <= a_ux && a_ux < hi.x {
                    let t = (a_ux - lo.x) / (hi.x - lo.x);
                    vs.push(lo.y + t * (hi.y - lo.y));
                }
            }
        }
        if vs.len() >= 2 {
            vs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            const HATCH_CONFUSION: f64 = 1e-8; // OCCT aTolHatch2D
            let mut k = 0;
            while k + 1 < vs.len() {
                let (v1, v2) = (vs[k], vs[k + 1]);
                if v2 - v1 > HATCH_CONFUSION {
                    a_vx = intermediate_point(v1, v2);
                    i_err = 0;
                    break;
                }
                k += 2;
            }
            if i_err == 0 {
                break;
            }
        }
        // OCCT L931-934: possible reason — incorrect computation of the 2d
        // box of the face; try again with the translated line.
        a_ux = a_umax - (a_ux - a_umin);
    }
    if i_err != 0 {
        return (i_err, glam::DVec3::ZERO, glam::DVec2::ZERO);
    }
    let a_p2d = glam::DVec2::new(a_ux, a_vx);
    let a_p = surf.point_at(a_p2d.x, a_p2d.y);
    (0, a_p, a_p2d)
}

/// OCCT IntTools_Context::IsValidPointForFace (IntTools_Context.cxx L648-674)
/// for a split face image: project the point on the piece's surface, check
/// the distance against aTol, then classify the UV point with the piece's UV
/// polygons (ON allowed — the point is valid when not strictly OUT). On a
/// U-periodic surface the classification uses the periodic ray test (the
/// polygon may cross the seam, e.g. u near 2*PI equivalent to 0).
fn is_valid_point_for_face_shape(
    p: glam::DVec3,
    face: &Shape,
    a_tol: f64,
    locations: &[glam::DAffine3],
) -> bool {
    let Some(surf) = face.as_face().and_then(|fd| fd.surface.clone()) else {
        return false;
    };
    let (uv, proj) = crate::bop::closest_point_on_surface(&surf, p);
    if (proj - p).length() > a_tol {
        return false;
    }
    let Some(polys) = shape_uv_polygons(face, locations) else {
        return false;
    };
    let is_u_per = matches!(
        &surf,
        rcad_kernel::geom::Surface3::Sphere(_)
            | rcad_kernel::geom::Surface3::Cylinder(_)
            | rcad_kernel::geom::Surface3::Cone(_)
            | rcad_kernel::geom::Surface3::Revolution(_)
            | rcad_kernel::geom::Surface3::Torus(_)
    );
    let period = if is_u_per { std::f64::consts::TAU } else { 0.0 };
    // Strictly inside the outer polygon or on its boundary; not strictly
    // inside any hole (hole boundaries count as ON, which is valid).
    for (pi, poly) in polys.iter().enumerate() {
        let n = poly.len();
        if n < 3 {
            continue;
        }
        let on_boundary = (0..n).any(|i| {
            let a = poly[i];
            let b = poly[(i + 1) % n];
            let ab = b - a;
            let len2 = ab.length_squared();
            if len2 < 1e-30 {
                return false;
            }
            let ap = uv - a;
            let t = (ap.dot(ab) / len2).clamp(0.0, 1.0);
            (ap - t * ab).length() <= a_tol
        });
        let inside = if period > 0.0 {
            point_in_periodic_polygon(uv, poly, period)
        } else {
            rcad_kernel::base::gprop::tri::point_in_polygon_2d(poly, uv)
        };
        if pi == 0 {
            // outer wire: must be inside or on.
            if !inside && !on_boundary {
                return false;
            }
        } else {
            // hole: must not be strictly inside.
            if inside {
                return false;
            }
        }
    }
    true
}

/// OCCT BOPTools_AlgoTools::AreFacesSameDomain (BOPTools_AlgoTools.cxx
/// L1139-1205) for split face images (Shapes not registered in the DS).
pub fn are_faces_same_domain_shapes(
    f1: &Shape,
    f2: &Shape,
    fuzz: f64,
    locations: &[glam::DAffine3],
) -> bool {
    // OCCT L1149-1155: find a point inside the first face.
    let (i_err, a_p1, _a_p2d1) = point_in_face_shape(f1, locations);
    if i_err != 0 {
        return false;
    }
    // OCCT L1162-1168: tolerances of the faces.
    let mut a_tol_f1 = f1.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
    let mut a_tol_f2 = f2.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
    // OCCT L1170-1182: maximal tolerance of the edges of the first face.
    let mut a_tol_e_max = -1.0;
    for a_e in face_edges_shapes(f1) {
        if !a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(true) {
            let a_tol_e = a_e.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
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
    // OCCT L1198-1199: checking criteria.
    let a_tol = a_tol_f1 + a_tol_f2 + fuzz.max(rcad_kernel::CONFUSION);
    // OCCT L1202: project and classify the point on the second face.
    is_valid_point_for_face_shape(a_p1, f2, a_tol, locations)
}

// ====================================================================
// PointInFace — OCCT BOPTools_AlgoTools3D::PointInFace
// (BOPTools_AlgoTools3D.cxx L906-938)
// ====================================================================
/// Finds a point inside the face `f` (OCCT L906-938 wrapper).
///
/// Returns `(iErr, theP, theP2D)`: `iErr == 0` on success, with `theP` the 3D
/// point and `theP2D` its UV parameters on the face surface.
pub(crate) fn point_in_face(f: usize, ds: &DS) -> (i32, glam::DVec3, glam::DVec2) {
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

// ====================================================================
// PointInFace (edge) — OCCT BOPTools_AlgoTools3D::PointInFace
// (BOPTools_AlgoTools3D.cxx L942-990)
// ====================================================================
/// The edge's pcurve on the face, keyed by the DS face index with the
/// surface-identity fallback (same semantics as builder.rs edge_pcurve_on_face).
fn edge_pcurve_of_face<'a>(a_e: &'a Shape, f: usize, ds: &DS) -> Option<&'a rcad_kernel::geom::Curve2d> {
    let ed = a_e.as_edge()?;
    // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the pcurve key
    // location is L.Predivided(E.Location()).
    let (fid, floc) = ds.face_key(f)?;
    let fkey = (
        fid,
        crate::bop::algo::compose_face_edge_pcurve_location(floc, a_e.location, &ds.locations),
    );
    let face = ds.shape(f);
    let surf = face.as_face().and_then(|fd| fd.surface.as_ref())?;
    ed.pcurves
        .get(&fkey)
        .or_else(|| {
            ed.pcurves.iter().find_map(|(k, v)| {
                if let Some(&fi) = ds.map_shape_index.get(k) {
                    if let Some(fs) = ds.face_surface(fi) {
                        if surface_same(surf, &fs) {
                            return Some(v);
                        }
                    }
                }
                None
            })
        })
        .map(|(pc, _, _)| pc)
}

/// rcad equivalent of the Geom2dHatch_Hatcher domain computation for a HALF
/// LINE L(s) = p0 + s*dir, s in [0, +inf) — the edge-normal hatcher line of
/// OCCT BOPTools_AlgoTools3D::PointInFace(theF, theE, theT, theDt2D, ...)
/// (BOPTools_AlgoTools3D.cxx L942-990). The line starts at the edge point p0;
/// the inside intervals are the pairs of sorted crossing parameters, with the
/// s = 0 side decided by the FClass2d state of p0. Returns the inside
/// intervals [s1, s2] with s >= 0.
fn hatch_line_intervals_dir(f: usize, p0: glam::DVec2, dir: glam::DVec2, ds: &DS) -> Vec<(f64, f64)> {
    let fclass = FClass2d::new(ds, f, ds.face_tolerance(f).max(rcad_kernel::CONFUSION));
    let mut ts: Vec<f64> = Vec::new();
    for poly in fclass.uv_polygons() {
        let n = poly.len();
        for i in 0..n {
            let a = poly[i];
            let b = poly[(i + 1) % n];
            let e = b - a;
            // Intersection of the boundary segment a + u*e with p0 + t*dir.
            let denom = e.x * dir.y - e.y * dir.x;
            if denom.abs() < 1e-15 {
                continue; // parallel to the line
            }
            let u = ((p0.x - a.x) * dir.y - (p0.y - a.y) * dir.x) / denom;
            if u < 0.0 || u >= 1.0 {
                continue; // half-open: each shared vertex counted once
            }
            let t = ((p0.x - a.x) * e.y - (p0.y - a.y) * e.x) / denom;
            if t >= 0.0 {
                ts.push(t);
            }
        }
    }
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    const HATCH_CONFUSION: f64 = 1e-8; // OCCT aTolHatch2D (IntTools_Context.cxx L350)
    // The line starts at p0: the s = 0 side state decides the first interval.
    let start_in = fclass.perform(ds, p0, true) != State::Out;
    let mut intervals: Vec<(f64, f64)> = Vec::new();
    let mut k = 0usize;
    if start_in {
        if let Some(&t0) = ts.first() {
            if t0 > HATCH_CONFUSION {
                intervals.push((0.0, t0));
            }
            k = 1;
        }
    }
    while k + 1 < ts.len() {
        let (v1, v2) = (ts[k], ts[k + 1]);
        if v2 - v1 > HATCH_CONFUSION {
            intervals.push((v1, v2));
        }
        k += 2;
    }
    intervals
}

/// OCCT BOPTools_AlgoTools3D::PointInFace (BOPTools_AlgoTools3D.cxx L942-990) —
/// finds a point inside the face on the 2D line through the edge's pcurve
/// point at `a_t`, in the direction perpendicular to the edge tangent
/// (reversed for REVERSED edge and face). The line is trimmed to [0, +inf)
/// and hatched against the face boundary; the first inside domain yields the
/// point at aV1 + theDt2D (when the domain is longer than theDt2D) or at the
/// domain middle.
pub(crate) fn point_in_face_edge(
    f: usize,
    a_e: &Shape,
    a_t: f64,
    d_t2d: f64,
    ds: &DS,
) -> (i32, glam::DVec3, glam::DVec2) {
    // OCCT L946-950: aC2D = CurveOnSurface(theE, theF, f, l); null -> iErr = 5.
    let pc = match edge_pcurve_of_face(a_e, f, ds) {
        Some(pc) => pc,
        None => return (5, glam::DVec3::ZERO, glam::DVec2::ZERO),
    };
    // OCCT L952-956: aC2D->D1(aT, aP2D, aV2D); aD2Dx = Dir(aV2D).
    let a_p2d = pc.point_at(a_t);
    let a_v2d = pc.derivative_at(a_t);
    if a_v2d.length_squared() < 1e-24 {
        // OCCT gp_Dir2d(aV2D) raises Standard_ConstructionError on a zero vector.
        return (5, glam::DVec3::ZERO, glam::DVec2::ZERO);
    }
    // OCCT L958-961: aD2D = (-aD2Dx.Y(), aD2Dx.X()).
    let mut a_d2d = glam::DVec2::new(-a_v2d.y, a_v2d.x).normalize();
    // OCCT L963-969: REVERSED edge / face reverses the direction.
    if a_e.orientation == Orientation::Reversed {
        a_d2d = -a_d2d;
    }
    let a_f = ds.shape(f);
    if a_f.orientation == Orientation::Reversed {
        a_d2d = -a_d2d;
    }
    // OCCT L971-990: hatch the trimmed half-line; take the first inside domain.
    let intervals = hatch_line_intervals_dir(f, a_p2d, a_d2d, ds);
    let Some(&(a_v1, a_v2)) = intervals.first() else {
        return (2, glam::DVec3::ZERO, glam::DVec2::ZERO);
    };
    // OCCT L1035: aVx = (theDt2D > 0 && (aV2 - aV1) > theDt2D) ? aV1 + theDt2D
    //                    : IntTools_Tools::IntermediatePoint(aV1, aV2).
    let a_vx = if d_t2d > 0.0 && (a_v2 - a_v1) > d_t2d {
        a_v1 + d_t2d
    } else {
        intermediate_point(a_v1, a_v2)
    };
    // OCCT L1037-1038: theL2D->D0(aVx, theP2D); aS->D0(theP2D, theP).
    let a_p2d_r = a_p2d + a_d2d * a_vx;
    let Some(a_s) = ds.face_surface(f) else {
        return (1, glam::DVec3::ZERO, a_p2d_r);
    };
    let a_p = a_s.point_at(a_p2d_r.x, a_p2d_r.y);
    (0, a_p, a_p2d_r)
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
/// split face is found by PointInFace (hatcher, L1347-1359), with the
/// PointNearEdge fallback over the non-degenerated, non-closed-on-face edges
/// (L1359-1373); the surface normals at that point are computed on both
/// surfaces, and the result is whether the normals point in opposite
/// directions (OCCT L1383-1435).
pub fn is_split_to_reverse_face(f_sp: &Shape, f_sr: &Shape, ds: &DS) -> (bool, i32) {
    // OCCT L1336-1341: same surface handle -> compare orientations only.
    // rcad: surfaces are value types, so "same surface" means geometric identity
    // of the analytic surface parameters.
    let surf_sp = f_sp.as_face().and_then(|fd| fd.surface.clone());
    let surf_or = f_sr.as_face().and_then(|fd| fd.surface.clone());
    if let (Some(s1), Some(s2)) = (&surf_sp, &surf_or) {
        let same = surface_same(s1, s2);
        if same {
            return (f_sp.orientation != f_sr.orientation, 0);
        }
    }
    // OCCT L1347-1359: PointInFace (hatcher) — the point inside the split face.
    // OCCT keeps the 2D point aP2DFSp and passes it directly to
    // GetNormalToSurface (L1384); rcad keeps the same p2d instead of
    // re-projecting the 3D point.
    let mut p3d: Option<glam::DVec3> = None;
    let mut p2d_sp: Option<glam::DVec2> = None;
    // OCCT L1351: PointInFace(theFSp, ...) works on any TopoDS_Face; a split
    // face image not registered in the DS uses the Shape-based hatcher.
    match ds.map_shape_index.get(&(f_sp.ptr_id(), f_sp.location)).copied() {
        Some(fi) => {
            let (err, p, p2d) = crate::bop::tools::algo_tools::point_in_face(fi, ds);
            if err == 0 {
                p3d = Some(p);
                p2d_sp = Some(p2d);
            }
        }
        None => {
            let (err, p, p2d) = crate::bop::tools::algo_tools::point_in_face_shape(f_sp, &ds.locations);
            if err == 0 {
                p3d = Some(p);
                p2d_sp = Some(p2d);
            }
        }
    }
    if p3d.is_none() {
        // OCCT L1359-1373: try to get the point near some not closed and not
        // degenerated edge on the split face. The 5-arg PointNearEdge
        // (AlgoTools3D.cxx L685-694) computes aT = IntermediatePoint(Range)
        // and forwards to the dT2D overload (L614-652), whose
        // IsPointInOnFace + hatcher-inside-point logic is reproduced here.
        let mut found = false;
        for a_es in face_edges_shapes(f_sp) {
            if a_es.as_edge().map(|ed| ed.degenerated).unwrap_or(true) {
                continue;
            }
            if crate::bop::algo::builder_solid::edge_closed_on_face(&a_es, f_sp) {
                continue; // OCCT L1363: !BRep_Tool::IsClosed(aESp, theFSp)
            }
            let (a_t1, a_t2) = match a_es.as_edge() {
                Some(ed) => (ed.range[0], ed.range[1]),
                None => (0.0, 0.0),
            };
            let a_t = intermediate_point(a_t1, a_t2);
            // OCCT L619-641: dT2D = 10 * MinStepIn2d (1e-4), x10 for
            // cylinder/sphere surfaces, max(2*(tolE+tolF)).
            let mut d_t2d = 10.0 * 1e-5;
            let surf = f_sp.as_face().and_then(|fd| fd.surface.clone());
            if matches!(surf, Some(rcad_kernel::geom::Surface3::Cylinder(_))
                | Some(rcad_kernel::geom::Surface3::Sphere(_)))
            {
                d_t2d = 10.0 * d_t2d;
            }
            let a_tol_e = a_es.as_edge().map(|ed| ed.tolerance).unwrap_or(0.0);
            let a_tol_f = f_sp.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
            let d_tx = 2.0 * (a_tol_e + a_tol_f);
            if d_tx > d_t2d {
                d_t2d = d_tx;
            }
            let (near, err6) = crate::bop::algo::builder::point_near_edge(
                &a_es, f_sp, a_t, d_t2d, ds,
            );
            if err6 == 1 {
                continue;
            }
            let (p2d, p3d_near) = near.unwrap_or((glam::DVec2::ZERO, glam::DVec3::ZERO));
            // OCCT L627-641: the point must be inside (or on) the face;
            // otherwise the hatcher inside-point is taken (or iErr = 2).
            let in_face = match ds.map_shape_index.get(&(f_sp.ptr_id(), f_sp.location)).copied() {
                Some(fi) => {
                    let class2d = FClass2d::new(ds, fi, ds.face_tolerance(fi));
                    class2d.perform(ds, p2d, true) != State::Out
                }
                // OCCT L619-641: IsPointInOnFace(aP, aF, aTol) works on any
                // TopoDS_Face via the face's UV polygons — the split-face
                // image is not registered in the DS, so classify it from its
                // own geometry (rcad: is_valid_point_for_face_shape).
                None => is_valid_point_for_face_shape(
                    p3d_near,
                    f_sp,
                    2.0 * (a_tol_e + a_tol_f),
                    &ds.locations,
                ),
            };
            if in_face {
                p3d = Some(p3d_near);
                p2d_sp = Some(p2d);
                found = true;
            } else if let Some(fi) =
                ds.map_shape_index.get(&(f_sp.ptr_id(), f_sp.location)).copied()
            {
                let (err2, p2, p2d2) = crate::bop::tools::algo_tools::point_in_face(fi, ds);
                if err2 == 0 {
                    p3d = Some(p2);
                    p2d_sp = Some(p2d2);
                    found = true;
                }
            }
            if found {
                break;
            }
        }
        if p3d.is_none() {
            // OCCT L1365-1371: the point has not been found (theError = 1).
            return (false, 1);
        }
    }
    let p3d = p3d.unwrap();
    // OCCT L1383-1392: normal direction of the split face at the point —
    // GetNormalToSurface(aSFSp, aP2DFSp.X(), aP2DFSp.Y()) uses the 2D point
    // from PointInFace/PointNearEdge directly.
    let mut dn_sp = match (&surf_sp, p2d_sp) {
        (Some(s), Some(uv)) => s.normal_at(uv.x, uv.y),
        (Some(s), None) => {
            let (uv_sp, _) = crate::bop::closest_point_on_surface(s, p3d);
            s.normal_at(uv_sp.x, uv_sp.y)
        }
        _ => return (false, 2),
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

/// All edge Shapes of a face (outer + inner wires).
fn face_edges_shapes(face: &Shape) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    if let TShape::Face(fd) = &*face.data {
        for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
            if let TShape::Wire(wd) = &*w.data {
                out.extend(wd.edges.iter().cloned());
            }
        }
    }
    out
}

/// OCCT BOPTools_AlgoTools::IsSplitToReverse(edge) (BOPTools_AlgoTools.cxx
/// L1456-1531) — true when the split edge a_sp has the opposite direction to
/// the original edge a_e. Samples up to 11 points on the split edge (after
/// FindValidRange), projects each onto the original edge, and compares the
/// tangent directions (EdgeTangent = curve derivative, reversed for a
/// REVERSED edge); the first point where both tangents are computable
/// decides. anErr mirrors the OCCT error codes (2/3/4).
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
    // OCCT L1480-1486: FindValidRange(theESp, f, l) — shrink the split edge
    // range so the sample point is inside both edges; fall back to the full
    // range when shrinking fails.
    let mut f = a_sp.as_edge().map(|ed| ed.range[0]).unwrap_or(0.0);
    let mut l = a_sp.as_edge().map(|ed| ed.range[1]).unwrap_or(0.0);
    if !find_valid_range_single_edge(a_sp, &mut f, &mut l) {
        if let Some(ed) = a_sp.as_edge() {
            f = ed.range[0];
            l = ed.range[1];
        }
    }
    // OCCT L1490-1524: try up to 11 sample points.
    let mut an_err = 0;
    let a_nb_p = 11;
    let a_dt = (l - f) / a_nb_p as f64;
    for i in 1..a_nb_p {
        let a_tm = f + i as f64 * a_dt;
        // OCCT L1495-1501: EdgeTangent(theESp, aTm) — curve derivative.
        let Some(v_sp_tgt) = edge_tangent_at_param(a_sp, a_tm) else {
            an_err = 2;
            continue;
        };
        // OCCT L1504-1511: ProjectPointOnEdge(aCSp->Value(aTm), theEOr, aTmOr).
        let Some((a_tm_or, _)) = project_point_on_edge(&c_sp.as_ref().unwrap().point_at(a_tm), a_e) else {
            an_err = 3;
            continue;
        };
        // OCCT L1513-1519: EdgeTangent(theEOr, aTmOr).
        let Some(v_or_tgt) = edge_tangent_at_param(a_e, a_tm_or) else {
            an_err = 4;
            continue;
        };
        // OCCT L1521-1522: aCos = aVSpTgt.Dot(aVOrTgt); return (aCos < 0.).
        let a_cos = v_sp_tgt.dot(v_or_tgt);
        return (a_cos < 0.0, 0);
    }
    (false, an_err)
}

/// OCCT BRepLib::FindValidRange(theEdge, theFirst, theLast) (BRepLib_1.cxx
/// L262-299) — shrink the edge range away from both endpoint tolerance
/// spheres using the multi-parameter FindValidRange (L173-258).
fn find_valid_range_single_edge(a_e: &Shape, out_f: &mut f64, out_l: &mut f64) -> bool {
    use rcad_kernel::topods::TShape;
    let Some(ed) = a_e.as_edge() else { return false };
    let Some(curve) = ed.curve.clone() else { return false };
    let [a_par_v0, a_par_v1] = ed.range;
    if a_par_v1 - a_par_v0 < rcad_kernel::PCONFUSION {
        return false;
    }
    // OCCT L270-271: TopExp::Vertices(theEdge, aV[0], aV[1]).
    let verts = [&ed.first, &ed.last];
    // OCCT L277-289: aTolV = Confusion() + Tolerance(aV); aPntV = Pnt(aV);
    // null vertices take the curve point at the range end with aTolE.
    let a_tol_e = ed.tolerance;
    let mut a_tol_v = [rcad_kernel::CONFUSION; 2];
    let mut a_pnt_v = [glam::DVec3::ZERO; 2];
    for i in 0..2 {
        match &*verts[i].data {
            TShape::Vertex(vd) => {
                a_tol_v[i] += vd.tolerance;
                a_pnt_v[i] = vd.point;
            }
            _ => {
                if !rcad_kernel::is_infinite_value(a_par_v0 + i as f64 * (a_par_v1 - a_par_v0)) {
                    a_tol_v[i] += a_tol_e;
                    a_pnt_v[i] = curve.point_at(if i == 0 { a_par_v0 } else { a_par_v1 });
                }
            }
        }
    }
    crate::bop::algo::pave_filler::find_valid_range_params(
        &curve, a_par_v0, a_par_v1, a_tol_e,
        a_pnt_v[0], a_tol_v[0], a_pnt_v[1], a_tol_v[1],
        out_f, out_l,
    )
}

/// OCCT BOPTools_AlgoTools2D::EdgeTangent (BOPTools_AlgoTools2D.cxx L578-607)
/// — the curve derivative D1 at aT, normalized, reversed for a REVERSED edge.
fn edge_tangent_at_param(e: &Shape, a_t: f64) -> Option<glam::DVec3> {
    let ed = e.as_edge()?;
    if ed.degenerated {
        return None;
    }
    let curve = ed.curve.as_ref()?;
    let mut a_tau = curve.derivative_at(a_t);
    let a_mod = a_tau.length();
    // OCCT BOPTools_AlgoTools2D::EdgeTangent (L88-96): if (mod >
    // gp::Resolution()) aTau /= mod; else return false (gp::Resolution() ==
    // 1e-15).
    if a_mod <= 1e-15 {
        return None;
    }
    a_tau /= a_mod;
    if e.orientation == rcad_kernel::topods::Orientation::Reversed {
        a_tau = -a_tau;
    }
    Some(a_tau)
}

/// OCCT IntTools_Context::ProjectPointOnEdge(aP, aE, aT) — project the 3D
/// point onto the edge's curve within the edge range; returns (param, point).
fn project_point_on_edge(a_p: &glam::DVec3, a_e: &Shape) -> Option<(f64, glam::DVec3)> {
    let ed = a_e.as_edge()?;
    let curve = ed.curve.as_ref()?;
    let [t1, t2] = ed.range;
    let proj =
        rcad_kernel::base::geom_api::project::closest_point_on_curve_range(curve, *a_p, t1, t2, 64);
    Some((proj.param, proj.point))
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
                && c1.x_dir.dot(c2.x_dir) > 1.0 - TOL
                && c1.y_dir.dot(c2.y_dir) > 1.0 - TOL
                && (c1.radius - c2.radius).abs() < TOL
        }
        (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
            (e1.center - e2.center).length() < TOL
                && e1.normal.dot(e2.normal) > 1.0 - TOL
                && (e1.major_radius - e2.major_radius).abs() < TOL
                && (e1.minor_radius - e2.minor_radius).abs() < TOL
                && e1.major_dir.dot(e2.major_dir) > 1.0 - TOL
        }
        _ => false,
    }
}

/// Geometric identity of two analytic surfaces (rcad equivalent of OCCT's
/// Geom_Surface handle equality in IsSplitToReverse L1338). Same direction for
/// the axis/normal; the reference directions (u_dir/v_dir) are not compared —
/// the parameterization orientation is irrelevant for the orientation decision.
/// Same-surface comparison helper for the debug trace.
fn surface_same_match(a: &Option<rcad_kernel::geom::Surface3>, b: &Option<rcad_kernel::geom::Surface3>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => surface_same(x, y),
        _ => false,
    }
}

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

// ====================================================================
// BOPAlgo_Tools::FillMap — graph edge connection (BOPAlgo_Tools.cxx L38-48)
// ====================================================================
/// OCCT BOPAlgo_Tools::FillMap(int, int, IndexedDataMap<int, List<int>>).
/// Adds bidirectional connection between two nodes in an adjacency map.
/// OCCT uses NCollection_IndexedDataMap (insertion-ordered keys); rcad uses
/// IndexMap so MakeBlocks iterates keys in OCCT's insertion order.
pub fn fill_map(n1: usize, n2: usize, the_map: &mut indexmap::IndexMap<usize, Vec<usize>>) {
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
/// OCCT BOPAlgo_Tools::MakeBlocks (BOPAlgo_Tools.hxx L46-86): iterate the
/// IndexedDataMap keys in insertion order; for each unvisited key start a
/// chain and BFS-append all connected unvisited nodes (the chain iterator
/// grows as elements are appended).  rcad mirrors that with an IndexMap and
/// an index-walked Vec chain.
pub fn make_blocks(
    the_map: &indexmap::IndexMap<usize, Vec<usize>>,
    the_blocks: &mut Vec<Vec<usize>>,
) {
    let mut a_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for (&n, _) in the_map {
        // OCCT: if (!aMFence.Add(n)) continue;
        if !a_fence.insert(n) {
            continue;
        }
        // Start the chain (OCCT: theMBlocks.Append; aChain.Append(n)).
        let mut a_chain: Vec<usize> = vec![n];
        // Look for connected elements: iterate the chain, appending new ones.
        let mut i = 0;
        while i < a_chain.len() {
            let n1 = a_chain[i];
            if let Some(a_li) = the_map.get(&n1) {
                for &n2 in a_li {
                    if a_fence.insert(n2) {
                        a_chain.push(n2);
                    }
                }
            }
            i += 1;
        }
        the_blocks.push(a_chain);
    }
}
// ====================================================================
// CorrectTolerances — OCCT BOPTools_AlgoTools::CorrectTolerances
// (BOPTools_AlgoTools_1.cxx L309-313) + CorrectShapeTolerances (L389-436)
// ====================================================================

/// In-place tolerance update of a vertex/edge TShape (single-threaded;
/// preserves the Arc sharing exactly like OCCT mutates TShapes in place).
fn update_shape_tolerance(v: &Shape, new_tol: f64, a_map_to_avoid: &std::collections::HashSet<(u64, u32)>) {
    if a_map_to_avoid.contains(&(v.ptr_id(), v.location)) {
        return;
    }
    // OCCT UpdateShape (BOPTools_AlgoTools_1.cxx L1066-1090): BRep_Builder
    // UpdateVertex/UpdateEdge — set the tolerance on the shared TShape.
    let raw = Arc::as_ptr(&v.data) as *mut TShape;
    unsafe {
        match &mut *raw {
            TShape::Vertex(vd) => vd.tolerance = new_tol,
            TShape::Edge(ed) => ed.tolerance = new_tol,
            _ => {}
        }
    }
}

/// OCCT BOPTools_AlgoTools::CheckEdge (BOPTools_AlgoTools_1.cxx L455-520) —
/// corrects the vertex tolerance from its distance to the edge 3D curve at
/// the vertex parameter (endpoints use the curve range bounds).
fn check_edge(brep: &BRep, e_idx: usize, a_max_tol: f64,
              a_map_to_avoid: &std::collections::HashSet<(u64, u32)>) {
    let (a_tol_e, a_c, a_range, a_v1, a_v2) = {
        let t = &brep.tshapes[e_idx];
        match &**t {
            TShape::Edge(ed) => (ed.tolerance, ed.curve.clone(), ed.range, ed.first.clone(), ed.last.clone()),
            _ => return,
        }
    };
    let Some(a_c) = a_c else { return };
    // OCCT L461: aE.Orientation(FORWARD); the stored first/last endpoints
    // correspond to the curve range bounds (First/Last of the GCurve).
    for (a_v, a_t) in [(&a_v1, a_range[0]), (&a_v2, a_range[1])] {
        let (a_pv, a_tol_v) = match &*a_v.data {
            TShape::Vertex(vd) => (vd.point, vd.tolerance),
            _ => continue,
        };
        // OCCT L470-474: aTol = max(aTolV, aTolE); dd = 0.1*aTol; aTol = aTol^2.
        let mut a_tol = a_tol_v.max(a_tol_e);
        let dd = 0.1 * a_tol;
        a_tol *= a_tol;
        // OCCT L500-506: endpoint check — aPC = aC->Value(First/Last).
        let a_pc = a_c.point_at(a_t);
        let a_d2 = a_pv.distance_squared(a_pc);
        if a_d2 > a_tol {
            let a_new_tolerance = a_d2.sqrt() + dd;
            if a_new_tolerance < a_max_tol {
                update_shape_tolerance(a_v, a_new_tolerance, a_map_to_avoid);
            }
        }
    }
}

/// OCCT BOPTools_AlgoTools::CorrectPointOnCurve (BOPTools_AlgoTools_1.cxx
/// L316-335) — builds a BOPTools_CPC per edge and performs it.
pub fn correct_point_on_curve(brep: &mut BRep,
                              a_map_to_avoid: &std::collections::HashSet<(u64, u32)>,
                              a_max_tol: f64) {
    for i in 0..brep.tshapes.len() {
        if matches!(&*brep.tshapes[i], TShape::Edge(_)) {
            check_edge(brep, i, a_max_tol, a_map_to_avoid);
        }
    }
}

/// OCCT BOPTools_AlgoTools::CorrectCurveOnSurface (BOPTools_AlgoTools_1.cxx
/// L337-386) — pcurve-based correction (BOPTools_CWT/CDT with
/// IntersectCurves2d). Pending translation.
pub fn correct_curve_on_surface(_brep: &mut BRep,
                                _a_map_to_avoid: &std::collections::HashSet<(u64, u32)>,
                                _a_max_tol: f64) {
    // OCCT L351-456: BOPTools_CWT/CDT Perform — 2D pcurve intersection based
    // tolerance correction. Pending (needs IntersectCurves2d translation).
}

/// OCCT BOPTools_AlgoTools::CorrectTolerances (BOPTools_AlgoTools_1.cxx
/// L309-313).
pub fn correct_tolerances(brep: &mut BRep,
                          a_map_to_avoid: &std::collections::HashSet<(u64, u32)>,
                          a_max_tol: f64) {
    correct_point_on_curve(brep, a_map_to_avoid, a_max_tol);
    correct_curve_on_surface(brep, a_map_to_avoid, a_max_tol);
}

/// OCCT BOPTools_AlgoTools::CorrectVertexTolerance (BOPTools_AlgoTools_1.cxx
/// L1005-1017) — vertex tolerance raised to the edge tolerance.
fn correct_vertex_tolerance(brep: &mut BRep,
                            a_map_to_avoid: &std::collections::HashSet<(u64, u32)>) {
    for i in 0..brep.tshapes.len() {
        let (a_tol_e, a_v1, a_v2) = {
            let t = &brep.tshapes[i];
            match &**t {
                TShape::Edge(ed) => (ed.tolerance, ed.first.clone(), ed.last.clone()),
                _ => continue,
            }
        };
        for a_v in [a_v1, a_v2] {
            let a_tol_v = match &*a_v.data {
                TShape::Vertex(vd) => vd.tolerance,
                _ => continue,
            };
            // OCCT L1010-1016: if (aTolV < aTolE) UpdateShape(aV, aTolE).
            if a_tol_v < a_tol_e {
                update_shape_tolerance(&a_v, a_tol_e, a_map_to_avoid);
            }
        }
    }
}

/// OCCT BOPTools_AlgoTools::CorrectEdgeTolerance (BOPTools_AlgoTools_1.cxx
/// L1020-1036) — edge tolerance raised to the face tolerance.
fn correct_edge_tolerance(brep: &mut BRep,
                          a_map_to_avoid: &std::collections::HashSet<(u64, u32)>) {
    for i in 0..brep.tshapes.len() {
        let (a_tol_f, wires) = {
            let t = &brep.tshapes[i];
            match &**t {
                TShape::Face(fd) => {
                    let mut ws = vec![fd.outer_wire.clone()];
                    ws.extend(fd.inner_wires.iter().cloned());
                    (fd.tolerance, ws)
                }
                _ => continue,
            }
        };
        for w in wires {
            let edges: Vec<Shape> = match &*w.data {
                TShape::Wire(wd) => wd.edges.clone(),
                _ => continue,
            };
            for a_e in edges {
                let a_tol_e = match &*a_e.data {
                    TShape::Edge(ed) => ed.tolerance,
                    _ => continue,
                };
                // OCCT L1032-1035: if (aTolE < aTolF) UpdateShape(aE, aTolF).
                if a_tol_e < a_tol_f {
                    update_shape_tolerance(&a_e, a_tol_f, a_map_to_avoid);
                }
            }
        }
    }
}

/// OCCT BOPTools_AlgoTools::CorrectShapeTolerances (BOPTools_AlgoTools_1.cxx
/// L389-436).
pub fn correct_shape_tolerances(brep: &mut BRep,
                                a_map_to_avoid: &std::collections::HashSet<(u64, u32)>) {
    correct_vertex_tolerance(brep, a_map_to_avoid);
    correct_edge_tolerance(brep, a_map_to_avoid);
}
