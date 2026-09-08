//! OCCT BRepFill_CompatibleWires (TKBool/BRepFill) — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill_CompatibleWires.cxx
//! (L125-2635) + BRepFill_CompatibleWires.hxx (L34-101).
//!
//! First consumer: BRepOffsetAPI_ThruSections (Stage 2e, the Build L380-509
//! compatibility pass).  The package helpers consumed below —
//! BRepFill::ComputeACR (BRepFill.cxx L1027-1063) and BRepFill::InsertACR
//! (BRepFill.cxx L1067-1200) — are translated in this file too.
//!
//! Reported gaps (plan §0.6, annotated at the call sites):
//! - BRepLib_FindSurface (TKTopAlgo) — used by PlaneOfWire;
//! - BRepGProp::LinearProperties / GProp principal axes — used by PlaneOfWire;
//! - BRepLProp::Continuity — used by WireContinuity;
//! - BRepExtrema_DistShapeShape — used by EdgeIntersectOnWire;
//! - BRepCheck_Wire — used by SameNumberByPolarMethod.
//!
//! Architecture notes:
//! - OCCT map key (TopTools_ShapeMapHasher == IsSame == same TShape) maps to
//!   `Shape::ptr_id()`.
//! - `BRepTools_WireExplorer` maps to the wire-ordered `TWireData::edges`
//!   list; CurrentVertex (BRepTools_WireExplorer.cxx L381/395/739-742) is the
//!   traversal-first vertex of the current edge.
//! - `BRepAdaptor_Curve` maps to direct `TEdgeData` curve/range access.
//! - `TopExp::Vertices(E, V1, V2)` maps to the stored (first, last) pair of
//!   `TEdgeData`, swapped by the edge wrapper orientation when the caller
//!   passes cumOri = true.

use glam::{DVec3, DVec3 as Vec3Alias};
use std::collections::HashMap;

use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{Curve3, CurveEval, Plane};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

use crate::brep_fill::generator::{
    shape_key, shape_oriented, shape_reversed, top_exp_vertices, top_exp_wire_vertices, ShapeKey,
    BRepFillThruSectionErrorStatus,
};

/// OCCT Precision::Confusion().
pub(super) const TOL_CONFUSION: f64 = CONFUSION;
/// OCCT Precision::Angular().
pub(super) const TOL_ANGULAR: f64 = rcad_kernel::core::precision::ANGULAR;
/// OCCT gp::Resolution().
pub(super) const GP_RESOLUTION: f64 = 2.2250738585072014e-308; // DBL_MIN, as in OCCT gp::Resolution
/// OCCT RealLast().
pub(super) const REAL_LAST: f64 = f64::MAX;

/// OCCT GeomAbs_Shape rank for the `<=` / `>=` comparisons (GeomAbs_ ordering
/// C0 < G1 < C1 < G2 < C2 < C3 < CN).
pub(super) fn continuity_rank(s: GeomAbsShapeAlias) -> i32 {
    use GeomAbsShapeAlias as G;
    match s {
        G::C0 => 0,
        G::G1 => 1,
        G::C1 => 2,
        G::G2 => 3,
        G::C2 => 4,
        G::C3 => 5,
        G::CN => 6,
    }
}
pub(super) type GeomAbsShapeAlias = rcad_kernel::topo::topods::GeomAbsShape;

// =============================================================================
// Kernel-mapping helpers
// =============================================================================

/// OCCT TopExp::Vertices(E, V1, V2) with cumOri = false — the stored
/// (FORWARD, REVERSED) vertex wrappers.
pub(super) fn top_exp_vertices_stored(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    let ed = brep.edge(e.clone());
    (ed.first.clone(), ed.last.clone())
}

/// OCCT TopExp::Vertices(E, V1, V2, cumOri = true) — traversal order.
pub(super) fn top_exp_vertices_cumori(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    top_exp_vertices(brep, e)
}

/// The wire-ordered edge list (BRepTools_WireExplorer mapping).
pub(super) fn wire_edges(brep: &BRep, w: &Shape) -> Vec<Shape> {
    match w.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT BRep_Tool::Degenerated.
pub(super) fn is_degenerated(brep: &BRep, e: &Shape) -> bool {
    brep.edge(e.clone()).degenerated
}

/// OCCT BRep_Tool::Pnt.
pub(super) fn vertex_point(brep: &BRep, v: &Shape) -> DVec3 {
    brep.vertex(v.clone()).point
}

/// OCCT BRep_Tool::Tolerance (vertex or edge).
pub(super) fn shape_tolerance(brep: &BRep, s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Parameter(V, E) — the stored vertex parameter.
pub(super) fn brep_tool_parameter(brep: &BRep, v: &Shape, e: &Shape) -> f64 {
    let ed = brep.edge(e.clone());
    ed.vertex_params
        .get(&v.ptr_id())
        .copied()
        .unwrap_or_else(|| {
            // OCCT guarantees the parameter is stored at edge creation; the
            // fallback picks the nearer range end by the vertex point.
            let p = brep.vertex(v.clone()).point;
            let c = match &ed.curve {
                Some(c) => c,
                None => return ed.range[0],
            };
            let d0 = c.point_at(ed.range[0]).distance(p);
            let d1 = c.point_at(ed.range[1]).distance(p);
            if d0 <= d1 {
                ed.range[0]
            } else {
                ed.range[1]
            }
        })
}

/// OCCT BRep_Tool::Range(E, first, last).
pub(super) fn brep_tool_range(brep: &BRep, e: &Shape) -> (f64, f64) {
    let r = brep.edge(e.clone()).range;
    (r[0], r[1])
}

/// OCCT BRepLib_MakeEdge(C, V0, V1, p0, p1) — null vertices are created at
/// the curve points (BRepLib_MakeEdge.cxx).
fn make_edge_curve_verts(
    brep: &mut BRep,
    c: &Curve3,
    v0: &Shape,
    v1: &Shape,
    p0: f64,
    p1: f64,
) -> Shape {
    let first = if v0.is_null() {
        brep.add_tvertex_unique(c.point_at(p0))
    } else {
        shape_oriented(v0, Orientation::Forward)
    };
    let last = if v1.is_null() {
        brep.add_tvertex_unique(c.point_at(p1))
    } else {
        shape_oriented(v1, Orientation::Reversed)
    };
    brep.add_tedge(Some(c.clone()), first, last, [p0.min(p1), p0.max(p1)])
}

// =============================================================================
// Static helpers (BRepFill_CompatibleWires.cxx)
// =============================================================================

/// OCCT static AddNewEdge (BRepFill_CompatibleWires.cxx L131-150): for
/// <theEdge> find all newest edges in <theEdgeNewEdges> recursively.
pub(super) fn add_new_edge(
    the_edge: &Shape,
    the_edge_new_edges: &HashMap<ShapeKey, Vec<Shape>>,
    list_new_edges: &mut Vec<Shape>,
) {
    if let Some(new_edges) = the_edge_new_edges.get(&shape_key(the_edge)) {
        for an_edge in new_edges {
            add_new_edge(an_edge, the_edge_new_edges, list_new_edges);
        }
    } else {
        list_new_edges.push(the_edge.clone());
    }
}

/// OCCT static SeqOfVertices (L152-172): the distinct vertices of a wire in
/// wire order.
pub(super) fn seq_of_vertices(brep: &BRep, w: &Shape) -> Vec<Shape> {
    let mut s: Vec<Shape> = Vec::new();
    // TopExp_Explorer PE(W, TopAbs_VERTEX) — the wire's edge endpoints.
    for e in wire_edges(brep, w) {
        let (v1, v2) = top_exp_vertices_stored(brep, &e);
        for v in [v1, v2] {
            let mut trouve = false;
            for sv in &s {
                if sv.is_same(&v) {
                    trouve = true;
                    break;
                }
            }
            if !trouve {
                s.push(v);
            }
        }
    }
    s
}

/// OCCT static PlaneOfWire (L174-298).  GAP: BRepLib_FindSurface and
/// BRepGProp::LinearProperties are not ported (plan §0.6, reported).
pub(super) fn plane_of_wire(_brep: &BRep, _w: &Shape, _p: &mut Plane) -> bool {
    panic!(
        "GAP: PlaneOfWire requires BRepLib_FindSurface + BRepGProp::LinearProperties \
         (TKTopAlgo, not translated) — see brep_fill_compatible_wires.rs header"
    )
}

/// OCCT static WireContinuity (L300-365): the minimum continuity of the wire.
pub(super) fn wire_continuity(brep: &BRep, w: &Shape) -> GeomAbsShapeAlias {
    use GeomAbsShapeAlias as G;
    let mut cont_w = G::CN;
    let mut is_degen = false;

    let edges = wire_edges(brep, w);
    let nb_edges = edges.len();
    for e in &edges {
        if is_degenerated(brep, e) {
            is_degen = true;
        }
    }

    if !is_degen {
        let mut testconti = true;

        for j in 0..nb_edges {
            let (edge1, edge2) = if j == nb_edges - 1 {
                (edges[nb_edges - 1].clone(), edges[0].clone())
            } else {
                (edges[j].clone(), edges[j + 1].clone())
            };

            // TopExp::Vertices(Edge1, Vbid, V1, true); (Edge2, V2, Vbid, true)
            let (_, v1) = top_exp_vertices_cumori(brep, &edge1);
            let (v2, _) = top_exp_vertices_cumori(brep, &edge2);
            let u1 = brep_tool_parameter(brep, &v1, &edge1);
            let u2 = brep_tool_parameter(brep, &v2, &edge2);
            let eps = shape_tolerance(brep, &v2) + shape_tolerance(brep, &v1);

            if j == nb_edges - 1 {
                let p1 = brep.edge(edge1.clone()).curve.as_ref().map(|c| c.point_at(u1));
                let p2 = brep.edge(edge2.clone()).curve.as_ref().map(|c| c.point_at(u2));
                testconti = match (p1, p2) {
                    (Some(a), Some(b)) => a.distance(b) <= eps,
                    _ => false,
                };
            }

            if testconti {
                // OCCT L357: cont = BRepLProp::Continuity(Curve1, Curve2, U1,
                // U2, Eps, Precision::Angular()).
                // GAP: BRepLProp::Continuity is not ported (plan §0.6,
                // reported) — the continuity update is skipped.
                let _ = (u1, u2, eps, TOL_ANGULAR);
            }
        }
    }
    let _ = cont_w;
    cont_w
}

/// OCCT static TrimEdge (L367-447): cut the edge at the cut values (wire
/// parameter space mapped to the edge range).
pub(super) fn trim_edge(
    brep: &mut BRep,
    current_edge: &Shape,
    cut_values: &[f64],
    t0: f64,
    t1: f64,
    seq_order: bool,
    s: &mut Vec<Shape>,
) {
    s.clear();
    let ndec = cut_values.len();
    let ed = brep.edge(current_edge.clone());
    let (first, last) = (ed.range[0], ed.range[1]);
    let c = match &ed.curve {
        Some(c) => c.clone(),
        None => return,
    };

    let current_orient = current_edge.orientation;
    // TopExp::Vertices(CurrentEdge, Vf, Vl) — non-oriented.
    let (vf, vl) = top_exp_vertices_stored(brep, current_edge);
    let vbid = Shape::null();

    if seq_order {
        // from first to last
        let mut m0 = first;
        let mut v0 = vf;
        for j in 0..ndec {
            // piece of edge
            let m1 = (cut_values[j] - t0) * (last - first) / (t1 - t0) + first;
            if (m0 - m1).abs() < TOL_CONFUSION {
                return;
            }
            let cut_e = make_edge_curve_verts(brep, &c, &v0, &vbid, m0, m1);
            let cut_e = shape_oriented(&cut_e, current_orient);
            s.push(cut_e.clone());
            m0 = m1;
            v0 = top_exp_first_vertex_stored(brep, &cut_e);
            let _ = &v0;
            v0 = top_exp_last_vertex_stored(brep, &cut_e);
            if j == ndec - 1 {
                // last piece
                if (m0 - last).abs() < TOL_CONFUSION {
                    return;
                }
                let last_e = make_edge_curve_verts(brep, &c, &v0, &vl, m0, last);
                let last_e = shape_oriented(&last_e, current_orient);
                s.push(last_e);
            }
        }
    } else {
        // from last to first
        let mut m1 = last;
        let mut v1 = vl;
        for jj in (0..ndec).rev() {
            let j = jj; // CutValues.Value(j) — 1-based in OCCT, 0-based here
            let _ = j;
            // piece of edge
            let m0 = (cut_values[jj] - t0) * (last - first) / (t1 - t0) + first;
            if (m0 - m1).abs() < TOL_CONFUSION {
                return;
            }
            let cut_e = make_edge_curve_verts(brep, &c, &vbid, &v1, m0, m1);
            let cut_e = shape_oriented(&cut_e, current_orient);
            s.push(cut_e.clone());
            m1 = m0;
            v1 = top_exp_first_vertex_stored(brep, &cut_e);
            if jj == 0 {
                // last piece
                if (first - m1).abs() < TOL_CONFUSION {
                    return;
                }
                let last_e = make_edge_curve_verts(brep, &c, &vf, &v1, first, m1);
                let last_e = shape_oriented(&last_e, current_orient);
                s.push(last_e);
            }
        }
    }
}

/// OCCT TopExp::FirstVertex(E) — the stored FORWARD vertex.
pub(super) fn top_exp_first_vertex_stored(brep: &BRep, e: &Shape) -> Shape {
    brep.edge(e.clone()).first.clone()
}

/// OCCT TopExp::LastVertex(E) — the stored REVERSED vertex.
pub(super) fn top_exp_last_vertex_stored(brep: &BRep, e: &Shape) -> Shape {
    brep.edge(e.clone()).last.clone()
}

/// OCCT static SearchRoot (L449-487): find the map key whose list contains V.
pub(super) fn search_root(
    v: &Shape,
    map: &HashMap<ShapeKey, Vec<Shape>>,
    brep: &BRep,
    v_root: &mut Shape,
) -> bool {
    let mut trouve = false;
    *v_root = Shape::null();
    for (k, list) in map {
        let mut ilyest = false;
        for s in list {
            let vcur = s.clone();
            if vcur.is_same(v) {
                ilyest = true;
            }
            if ilyest {
                break;
            }
        }
        if ilyest {
            trouve = true;
            // VRoot = TopoDS::Vertex(it.Key())
            *v_root = find_vertex_by_key(brep, *k);
        }
        if trouve {
            break;
        }
    }
    trouve
}

/// OCCT static SearchVertex (L489-525): find the wire vertex contained in the
/// list.
pub(super) fn search_vertex(brep: &BRep, list: &[Shape], w: &Shape, von_w: &mut Shape) -> bool {
    let mut trouve = false;
    *von_w = Shape::null();
    let seq_v = seq_of_vertices(brep, w);
    for vi in &seq_v {
        let mut ilyest = false;
        for s in list {
            if s.is_same(vi) {
                ilyest = true;
            }
            if ilyest {
                break;
            }
        }
        if ilyest {
            trouve = true;
            *von_w = vi.clone();
        }
        if trouve {
            break;
        }
    }
    trouve
}

/// OCCT static EdgeIntersectOnWire (L527-680).  GAP: the intersection is
/// computed by BRepExtrema_DistShapeShape (TKTopAlgo), not ported (plan
/// §0.6, reported).
#[allow(clippy::too_many_arguments)]
pub(super) fn edge_intersect_on_wire(
    _brep: &mut BRep,
    _p1: DVec3,
    _p2: DVec3,
    _percent: f64,
    _map: &HashMap<ShapeKey, Vec<Shape>>,
    _w: &Shape,
    _vsol: &mut Shape,
    _new_w: &mut Shape,
    _the_edge_new_edges: &mut HashMap<ShapeKey, Vec<Shape>>,
) -> bool {
    panic!(
        "GAP: EdgeIntersectOnWire requires BRepExtrema_DistShapeShape \
         (TKTopAlgo/BRepExtrema, not translated) — see file header"
    )
}

/// OCCT gp_Vec::AngleWithRef(TheV, VRef) — the signed angle from `a` to `b`
/// about `vref`, in [0, 2*pi).
pub(super) fn angle_with_ref(a: DVec3, b: DVec3, vref: DVec3) -> f64 {
    let la = a.length();
    let lb = b.length();
    if la < 1e-300 || lb < 1e-300 {
        return 0.0;
    }
    let mut ang = (a.dot(b) / (la * lb)).clamp(-1.0, 1.0).acos();
    if a.cross(b).dot(vref) < 0.0 {
        ang = 2.0 * std::f64::consts::PI - ang;
    }
    ang
}

/// OCCT gp_Pnt::Rotated(A1, Ang) — rotate about the axis.
pub(super) fn pnt_rotated(p: DVec3, axis_loc: DVec3, axis_dir: DVec3, ang: f64) -> DVec3 {
    let d = p - axis_loc;
    let c = ang.cos();
    let s = ang.sin();
    let u = axis_dir.normalize_or_zero();
    let rotated = d * c + u.cross(d) * s + u * (u.dot(d)) * (1.0 - c);
    axis_loc + rotated
}

/// OCCT static Transform (L682-734): translate P from (Pos1, Ax1) to
/// (Pos2, Ax2), rotating when the axes are not parallel.
pub(super) fn transform(
    with_rotation: bool,
    p: DVec3,
    pos1: DVec3,
    ax1: DVec3,
    pos2: DVec3,
    ax2: DVec3,
    pnew: &mut DVec3,
) {
    // Pnew = P.Translated(Pos1, Pos2)
    *pnew = p + (pos2 - pos1);
    let mut axe1 = ax1;
    let mut axe2 = ax2;
    // gp_Vec::IsParallel with 1.e-4
    if !dir_is_parallel_1e4(axe1, axe2) {
        let vtrans = pos2 - pos1;
        let mut sign = 1.0f64;
        let mut alpha = vtrans.dot(axe1);
        let mut beta = vtrans.dot(axe2);
        if alpha < -1.0e-7 {
            axe1 = -axe1;
        }
        if beta < 1.0e-7 {
            axe2 = -axe2;
        }
        alpha = vtrans.dot(axe1);
        beta = vtrans.dot(axe2);
        let norm2 = axe1.cross(axe2);
        // Vsign.SetLinearForm(Vtrans.Dot(axe1), axe2, -Vtrans.Dot(axe2), axe1)
        let vsign = axe2 * vtrans.dot(axe1) - axe1 * vtrans.dot(axe2);
        alpha = vsign.dot(axe1);
        beta = vsign.dot(axe2);
        let pasnul = alpha.abs() > 1.0e-4 && beta.abs() > 1.0e-4;
        if alpha * beta > 0.0 && pasnul {
            sign = -1.0;
        }
        let mut ang = angle_with_ref(axe1, axe2, norm2);
        if !with_rotation {
            if ang > std::f64::consts::PI / 2.0 {
                ang -= std::f64::consts::PI;
            }
            if ang < -std::f64::consts::PI / 2.0 {
                ang += std::f64::consts::PI;
            }
        }
        ang *= sign;
        *pnew = pnt_rotated(*pnew, pos2, norm2, ang);
    }
}

/// gp_Vec::IsParallel with the 1.e-4 tolerance used by Transform (L693).
pub(super) fn dir_is_parallel_1e4(a: DVec3, b: DVec3) -> bool {
    let la = a.length();
    let lb = b.length();
    if la < 1e-300 || lb < 1e-300 {
        return false;
    }
    let cross = a.cross(b).length() / (la * lb);
    cross <= 1.0e-4_f64.sin().abs() + 1e-18
}

/// OCCT static BuildConnectedEdges (L736-770).
pub(super) fn build_connected_edges(brep: &BRep, a_wire: &Shape, start_edge: &Shape, start_vertex: &Shape, connected_edges: &mut Vec<Shape>) {
    // TopExp::MapShapesAndAncestors(aWire, TopAbs_VERTEX, TopAbs_EDGE, MapVE)
    let map_ve = map_shapes_and_ancestors_ve(brep, a_wire);
    let mut cur_edge = start_edge.clone();
    let mut cur_vertex = start_vertex.clone();
    let (v1, v2) = top_exp_vertices(brep, start_edge);
    let origin = if v1.is_same(start_vertex) { v2 } else { v1 };

    loop {
        // NCollection_List::Iterator itE(MapVE.FindFromKey(CurVertex))
        let list = map_ve
            .iter()
            .find(|(v, _)| v.is_same(&cur_vertex))
            .map(|(_, l)| l.clone())
            .unwrap_or_default();
        for e in &list {
            if !e.is_same(&cur_edge) {
                connected_edges.push(e.clone());
                let (v1, v2) = top_exp_vertices(brep, e);
                cur_vertex = if v1.is_same(&cur_vertex) { v2 } else { v1 };
                cur_edge = e.clone();
                break;
            }
        }
        if cur_vertex.is_same(&origin) {
            break;
        }
    }
}

/// OCCT TopExp::MapShapesAndAncestors(W, VERTEX, EDGE, Map) — the ordered
/// (vertex -> edges) incidence map of a wire.
pub(super) fn map_shapes_and_ancestors_ve(brep: &BRep, w: &Shape) -> Vec<(Shape, Vec<Shape>)> {
    let mut map: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for e in wire_edges(brep, w) {
        let ed = brep.edge(e.clone());
        for v in [&ed.first, &ed.last] {
            match map.iter_mut().find(|(kv, _)| kv.is_same(v)) {
                Some((_, l)) => l.push(e.clone()),
                None => map.push((v.clone(), vec![e.clone()])),
            }
        }
    }
    map
}

/// Recover a shape handle from its map key (OCCT keeps handles directly).
pub(super) fn find_vertex_by_key(brep: &BRep, k: ShapeKey) -> Shape {
    use std::sync::Arc;
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if Arc::as_ptr(ts) as u64 == k.0 {
            return brep.shape_at(i);
        }
    }
    Shape::null()
}

// =============================================================================
// BRepFill::ComputeACR / BRepFill::InsertACR (BRepFill.cxx L1027-1200)
// =============================================================================

/// OCCT BRepFill::ComputeACR (BRepFill.cxx L1027-1063): the reduced
/// curvilinear abscissas ACR(1..nbEdges) and the total length ACR(0).
pub(super) fn compute_acr(brep: &BRep, wire: &Shape) -> Vec<f64> {
    // calculate the reduced curvilinear abscisses and the length of the wire
    let mut nb_edges = 0usize;
    // cumulated lengths (ACR(0..nbE) — index 0 is the total length slot)
    let mut acr: Vec<f64> = vec![0.0; 1];
    for ecur in wire_edges(brep, wire) {
        nb_edges += 1;
        let prev = acr[nb_edges - 1];
        acr.push(prev);
        if !is_degenerated(brep, &ecur) {
            let ed = brep.edge(ecur.clone());
            if let Some(c) = &ed.curve {
                // BRepAdaptor_Curve + GCPnts_AbscissaPoint::Length
                acr[nb_edges] += rcad_kernel::base::gcpnts::abscissa_point::arc_length(
                    c,
                    ed.range[0],
                    ed.range[1],
                );
            }
        }
    }

    // total length of the wire
    acr[0] = acr[nb_edges];

    // reduced curvilinear abscisses
    let total = acr[0];
    if total > TOL_CONFUSION {
        for a in acr.iter_mut().skip(1) {
            *a /= total;
        }
    } else {
        // punctual wire
        acr[nb_edges] = 1.0;
    }
    acr
}

/// OCCT static TrimEdge (BRepFill.cxx L143-198) — the BRepFill.cxx variant
/// (no degenerate-return guards), used by InsertACR.
pub(super) fn trim_edge_brep_fill(
    brep: &mut BRep,
    current_edge: &Shape,
    cut_values: &[f64],
    t0: f64,
    t1: f64,
    seq_order: bool,
    s: &mut Vec<Shape>,
) {
    s.clear();
    let ndec = cut_values.len();
    let ed = brep.edge(current_edge.clone());
    let (first, last) = (ed.range[0], ed.range[1]);
    let c = match &ed.curve {
        Some(c) => c.clone(),
        None => return,
    };

    let current_orient = current_edge.orientation;
    let (vf, vl) = top_exp_vertices_stored(brep, current_edge);
    let vbid = Shape::null();

    if seq_order {
        // from first to last
        let mut m0 = first;
        let mut v0 = vf;
        for j in 0..ndec {
            // piece of edge
            let m1 = (cut_values[j] - t0) * (last - first) / (t1 - t0) + first;
            let cut_e = make_edge_curve_verts(brep, &c, &v0, &vbid, m0, m1);
            let cut_e = shape_oriented(&cut_e, current_orient);
            s.push(cut_e.clone());
            m0 = m1;
            v0 = top_exp_last_vertex_stored(brep, &cut_e);
            if j == ndec - 1 {
                // last piece
                let last_e = make_edge_curve_verts(brep, &c, &v0, &vl, m0, last);
                let last_e = shape_oriented(&last_e, current_orient);
                s.push(last_e);
            }
        }
    } else {
        // from last to first
        let mut m1 = last;
        let mut v1 = vl;
        for jj in (0..ndec).rev() {
            // piece of edge
            let m0 = (cut_values[jj] - t0) * (last - first) / (t1 - t0) + first;
            let cut_e = make_edge_curve_verts(brep, &c, &vbid, &v1, m0, m1);
            let cut_e = shape_oriented(&cut_e, current_orient);
            s.push(cut_e.clone());
            m1 = m0;
            v1 = top_exp_first_vertex_stored(brep, &cut_e);
            if jj == 0 {
                // last piece
                let last_e = make_edge_curve_verts(brep, &c, &vf, &v1, first, m1);
                let last_e = shape_oriented(&last_e, current_orient);
                s.push(last_e);
            }
        }
    }
}

/// OCCT BRepFill::InsertACR (BRepFill.cxx L1067-1200): insert the ACR cuts
/// into the wire.
pub(super) fn insert_acr(brep: &mut BRep, wire: &Shape, acr_cuts: &[f64], prec: f64) -> Shape {
    // calculate ACR of the wire to be cut
    let mut nb_edges = 0usize;
    for _ in wire_edges(brep, wire) {
        nb_edges += 1;
    }
    let acr_wire = compute_acr(brep, wire);

    let nmax = acr_cuts.len();
    let mut paradec: Vec<f64> = vec![0.0; nmax];
    let mut mw_edges: Vec<Shape> = Vec::new(); // BRepLib_MakeWire MW

    let mut t0 = 0.0f64;
    let mut t1 = 0.0f64;
    nb_edges = 0;

    // processing edge by edge
    for e in wire_edges(brep, wire) {
        nb_edges += 1;
        t0 = t1;
        t1 = acr_wire[nb_edges];

        // parameters of cut on this edge
        let mut ndec = 0usize;
        for i in 0..acr_cuts.len() {
            if t0 + prec < acr_cuts[i] && acr_cuts[i] < t1 - prec {
                paradec[ndec] = acr_cuts[i];
                ndec += 1;
            }
        }

        // const TopoDS_Edge& E = anExp.Current();
        let e = e.clone();
        // const TopoDS_Vertex& V = anExp.CurrentVertex();
        let v = top_exp_vertices_cumori(brep, &e).0;

        if ndec == 0 || is_degenerated(brep, &e) {
            // copy the edge
            mw_edges.push(e);
        } else {
            // it is necessary to cut the edge
            // following the direction of parsing of the wire
            let so = v.is_same(&top_exp_first_vertex_stored(brep, &e));
            let mut se: Vec<Shape> = Vec::new();
            let mut sr: Vec<f64> = Vec::new();
            // the wire is always FORWARD
            // it is necessary to modify the parameter of cut6 if the edge is REVERSED
            if e.orientation == Orientation::Forward {
                for j in 0..ndec {
                    sr.push(paradec[j]);
                }
            } else {
                for j in 0..ndec {
                    sr.push(t0 + t1 - paradec[ndec + 1 - j - 1]);
                }
            }
            trim_edge_brep_fill(brep, &e, &sr, t0, t1, so, &mut se);
            for j in 0..se.len() {
                mw_edges.push(se[j].clone());
            }
        }
    }

    // result
    let orien = wire.orientation;
    let wres = brep.add_twire(mw_edges);
    shape_oriented(&wres, orien)
}

// =============================================================================
// BRepFill_CompatibleWires class (BRepFill_CompatibleWires.hxx L34-101)
// =============================================================================

/// OCCT BRepFill_CompatibleWires — constructs a sequence of wires with good
/// orientation and origin, agreed each other so that the surface passing
/// through these sections is not twisted.
#[derive(Debug)]
pub struct CompatibleWires {
    pub(super) my_init: Vec<Shape>,
    pub(super) my_work: Vec<Shape>,
    pub(super) my_percent: f64,
    pub(super) my_degen1: bool,
    pub(super) my_degen2: bool,
    pub(super) my_status: BRepFillThruSectionErrorStatus,
    pub(super) my_map: HashMap<ShapeKey, Vec<Shape>>,
}

impl Default for CompatibleWires {
    fn default() -> Self {
        Self::new()
    }
}

impl CompatibleWires {
    /// OCCT BRepFill_CompatibleWires::BRepFill_CompatibleWires (L774-777).
    pub fn new() -> Self {
        CompatibleWires {
            my_init: Vec::new(),
            my_work: Vec::new(),
            my_percent: 0.0,
            my_degen1: false,
            my_degen2: false,
            my_status: BRepFillThruSectionErrorStatus::NotDone,
            my_map: HashMap::new(),
        }
    }

    /// OCCT BRepFill_CompatibleWires::BRepFill_CompatibleWires(Sections)
    /// (L781-785).
    pub fn new_with_sections(brep: &BRep, sections: &[Shape]) -> Self {
        let mut me = Self::new();
        me.init(brep, sections);
        me
    }

    /// OCCT BRepFill_CompatibleWires::Init (L789-796).
    pub fn init(&mut self, _brep: &BRep, sections: &[Shape]) {
        self.my_init = sections.to_vec();
        self.my_work = sections.to_vec();
        self.my_percent = 0.1;
        self.my_status = BRepFillThruSectionErrorStatus::NotDone;
        self.my_map.clear();
    }

    /// OCCT BRepFill_CompatibleWires::SetPercent (L800-806).
    pub fn set_percent(&mut self, percent: f64) {
        if 0.0 < percent && percent < 1.0 {
            self.my_percent = percent;
        }
    }

    /// OCCT BRepFill_CompatibleWires::IsDone (L810-813).
    pub fn is_done(&self) -> bool {
        self.my_status == BRepFillThruSectionErrorStatus::Done
    }

    /// OCCT BRepFill_CompatibleWires::GetStatus (hxx L53).
    pub fn get_status(&self) -> BRepFillThruSectionErrorStatus {
        self.my_status
    }

    /// OCCT BRepFill_CompatibleWires::Shape (L817-820).
    pub fn shape(&self) -> &Vec<Shape> {
        &self.my_work
    }

    /// OCCT BRepFill_CompatibleWires::GeneratedShapes (L824-837).
    pub fn generated_shapes(&self, sub_section: &Shape) -> Vec<Shape> {
        self.my_map
            .get(&shape_key(sub_section))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT BRepFill_CompatibleWires::Generated (L1000-1004).
    pub fn generated(&self) -> &HashMap<ShapeKey, Vec<Shape>> {
        &self.my_map
    }

    /// OCCT BRepFill_CompatibleWires::IsDegeneratedFirstSection (L841-844).
    pub fn is_degenerated_first_section(&self) -> bool {
        self.my_degen1
    }

    /// OCCT BRepFill_CompatibleWires::IsDegeneratedLastSection (L848-851).
    pub fn is_degenerated_last_section(&self) -> bool {
        self.my_degen2
    }
}

/// OCCT TopoDS_Shape::Closed() — the wire Closed flag.
pub(super) trait WireFlags {
    fn flags_closed(&self, brep: &BRep) -> bool;
}
impl WireFlags for Shape {
    fn flags_closed(&self, brep: &BRep) -> bool {
        match self.data.as_ref() {
            TShape::Wire(wd) => wd.flags & tshape_flags::CLOSED != 0,
            _ => false,
        }
    }
}

/// Set the wire Closed flag (OCCT TopoDS_Wire::Closed(B)).
pub(super) fn set_wire_closed(brep: &mut BRep, w: &Shape, closed: bool) {
    use std::sync::Arc;
    let idx = w.index;
    if let TShape::Wire(wd) = Arc::make_mut(&mut brep.tshapes[idx]) {
        if closed {
            wd.flags |= tshape_flags::CLOSED;
        } else {
            wd.flags &= !tshape_flags::CLOSED;
        }
    }
}

/// ComputeOrigin sample term (forward pairing, OCCT L1945-1966): the
/// distance sum between the sampled points of the two edges.
pub(super) fn sample_distance(brep: &BRep, prev_edge: &Shape, cur_edge: &Shape, nb_samples: usize, offset: DVec3) -> f64 {
    let mut sum = 0.0f64;
    let prev_curve = brep.edge(prev_edge.clone()).curve.clone();
    let cur_curve = brep.edge(cur_edge.clone()).curve.clone();
    let (prev_range, cur_range) = (
        brep.edge(prev_edge.clone()).range,
        brep.edge(cur_edge.clone()).range,
    );
    if let (Some(pc), Some(cc)) = (prev_curve, cur_curve) {
        let sample_on_prev = (prev_range[1] - prev_range[0]) / nb_samples as f64;
        let sample_on_cur = (cur_range[1] - cur_range[0]) / nb_samples as f64;
        for nbs in 1..nb_samples {
            let par_on_prev = if prev_edge.orientation == Orientation::Forward {
                prev_range[0] + nbs as f64 * sample_on_prev
            } else {
                prev_range[0] + (nb_samples - nbs) as f64 * sample_on_prev
            };
            let par_on_cur = if cur_edge.orientation == Orientation::Forward {
                cur_range[0] + nbs as f64 * sample_on_cur
            } else {
                cur_range[0] + (nb_samples - nbs) as f64 * sample_on_cur
            };
            let pon_prev = pc.point_at(par_on_prev);
            let pon_cur = cc.point_at(par_on_cur) + offset;
            sum += pon_prev.distance(pon_cur);
        }
    }
    sum
}

/// ComputeOrigin sample term (backward pairing, OCCT L2032-2042 / L2060-2073):
/// the sample parameter roles are swapped relative to the forward case.
pub(super) fn sample_distance_backward(brep: &BRep, prev_edge: &Shape, cur_edge: &Shape, nb_samples: usize, offset: DVec3) -> f64 {
    let mut sum = 0.0f64;
    let prev_curve = brep.edge(prev_edge.clone()).curve.clone();
    let cur_curve = brep.edge(cur_edge.clone()).curve.clone();
    let (prev_range, cur_range) = (
        brep.edge(prev_edge.clone()).range,
        brep.edge(cur_edge.clone()).range,
    );
    if let (Some(pc), Some(cc)) = (prev_curve, cur_curve) {
        let sample_on_prev = (prev_range[1] - prev_range[0]) / nb_samples as f64;
        let sample_on_cur = (cur_range[1] - cur_range[0]) / nb_samples as f64;
        for nbs in 1..nb_samples {
            let par_on_prev = if prev_edge.orientation == Orientation::Forward {
                prev_range[0] + nbs as f64 * sample_on_prev
            } else {
                prev_range[0] + (nb_samples - nbs) as f64 * sample_on_prev
            };
            let par_on_cur = if cur_edge.orientation == Orientation::Forward {
                cur_range[0] + (nb_samples - nbs) as f64 * sample_on_cur
            } else {
                cur_range[0] + nbs as f64 * sample_on_cur
            };
            let pon_prev = pc.point_at(par_on_prev);
            let pon_cur = cc.point_at(par_on_cur) + offset;
            sum += pon_prev.distance(pon_cur);
        }
    }
    sum
}
