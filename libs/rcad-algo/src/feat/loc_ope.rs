// OCCT LocOpe.hxx L34-52 + LocOpe.cxx L39-235 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe.cxx
//
// OCCT inheritance chain: LocOpe is a package class — every member is a
// static free function. The commented-out LocOpe::IsInside of the OCCT
// source (cxx L237-263) is dead source text and is not translated.
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> — HashMap keyed
//    by (TShape ptr, Location) (the TopTools_ShapeMapHasher identity,
//    orientation ignored). Where the map is iterated (top_exp_vertices_wire,
//    the "open" branch) an insertion-order Vec stands in for the OCCT bucket
//    iterator — the same reduction as feat::loc_ope_glued_shape arch.
//    diff. #1.
// 2. BRep_Tool::CurveOnSurface(edg, face, f, l) is re-hosted below
//    (brep_tool_curve_on_surface): the pcurve lives on TEdgeData keyed by
//    the face identity (same model as feat::loc_ope_gluer).
// 3. BRepAdaptor_Surface / BRepAdaptor_Curve2d (TgtFaces) — evaluated
//    directly on the TFaceData surface (SurfaceEval::derivatives) and the
//    edge pcurve (Curve2dEval::point_at); standalone feat shapes carry
//    identity locations (feat::loc_ope_find_edges arch. diff. #1).
// 4. BRep_Tool::Curve(edg, Loc, f, l) + C->Transformed(Loc.Transformation())
//    is crate::feat::loc_ope_find_edges::brep_tool_curve (pub(crate)
//    re-host); the location transform applies only when the edge carries an
//    identity location (same arch. difference #1).
// 5. BRepLib_MakeWire in Closed(E, F) runs on a local rcad topods::BRep
//    pool (the same vehicle as feat::loc_ope_build_shape::perform arch.
//    diff. #1 there).
// 6. TopExp::Vertices (TopExp.cxx) has no rcad re-host yet — both overloads
//    (edge/wire) are re-hosted below (top_exp_vertices_edge /
//    top_exp_vertices_wire).
//
// first consumer: LocOpe_Prism/Revol/DPrism Curves/BarycCurve (SampleEdges);
// BRepFeat_Form family (3b) for Closed/TgtFaces.

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_find_edges::brep_tool_curve;
use glam::DVec2;
use rcad_kernel::geom::{ Curve2d, Curve2dEval, Curve3, CurveEval, SurfaceEval };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, TShape };
use std::collections::HashMap;

/// OCCT #define NECHANT 10 (LocOpe.cxx L39).
const NECHANT: i32 = 10;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
fn with_orientation(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> glam::DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Range(e, f, l) (BRep_Tool.cxx).
fn brep_tool_range(edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edg, face, f, l) (architecture
/// difference #2) — the pcurve of the edge on the face with its range.
fn brep_tool_curve_on_surface(edg: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT TopExp::Vertices(const TopoDS_Edge&, Vfirst, Vlast, CumOri=false)
/// (TopExp.cxx; architecture difference #6) — the FORWARD/REVERSED stored
/// vertices of the edge (TopoDS_Iterator(E, CumOri=false): stored
/// orientations, not composed with the edge orientation; TEdgeData.first/
/// last carry them).
pub(crate) fn top_exp_vertices_edge(e: &Shape) -> (Option<Shape>, Option<Shape>) {
    // OCCT: isFirstDefined/isLastDefined boolean flags.
    let mut is_first_defined = false;
    let mut is_last_defined = false;
    let mut v_first: Option<Shape> = None;
    let mut v_last: Option<Shape> = None;

    // OCCT: TopoDS_Iterator ite(E, CumOri) over the edge vertices.
    let sub_vertices: Vec<Shape> = match e.data.as_ref() {
        TShape::Edge(ed) => vec![ed.first.clone(), ed.last.clone()],
        _ => Vec::new(),
    };
    for a_v in &sub_vertices {
        if a_v.orientation == Orientation::Forward {
            v_first = Some(a_v.clone());
            is_first_defined = true;
        } else if a_v.orientation == Orientation::Reversed {
            v_last = Some(a_v.clone());
            is_last_defined = true;
        }
    }
    // OCCT: if (!isFirstDefined) Vfirst.Nullify(); (and the Vlast twin).
    if !is_first_defined {
        v_first = None;
    }
    if !is_last_defined {
        v_last = None;
    }
    (v_first, v_last)
}

/// The OCCT vmap toggle of TopExp::Vertices(Wire) —
/// `if (!vmap.Add(V)) vmap.Remove(V)` over the (TShape, Location) key
/// (arch. difference #1: the order Vec is the bucket-iterator stand-in).
fn vmap_toggle(
    vmap: &mut HashMap<(u64, u32), Shape>,
    vmap_order: &mut Vec<(u64, u32)>,
    v: Shape,
) {
    let key = shape_key(&v);
    if vmap.insert(key, v).is_none() {
        vmap_order.push(key);
    } else {
        vmap.remove(&key);
        if let Some(pos) = vmap_order.iter().position(|&k| k == key) {
            vmap_order.remove(pos);
        }
    }
}

/// OCCT TopExp::Vertices(const TopoDS_Wire&, Vfirst, Vlast) (TopExp.cxx;
/// architecture difference #6).
pub(crate) fn top_exp_vertices_wire(w: &Shape) -> (Option<Shape>, Option<Shape>) {
    // OCCT: Vfirst = Vlast = TopoDS_Vertex().
    let mut v_first: Option<Shape> = None;
    let mut v_last: Option<Shape> = None;

    // OCCT: NCollection_Map vmap (arch. diff. #1).
    let mut vmap: HashMap<(u64, u32), Shape> = HashMap::new();
    let mut vmap_order: Vec<(u64, u32)> = Vec::new();

    let mut v1: Option<Shape> = None;
    let mut v2: Option<Shape> = None;

    // OCCT: TopoDS_Iterator it(W) — cumOri=true (composed orientations;
    // crate::feat::brep_feat_builder::explorer composes).
    for s in explorer(w, ShapeType::Edge, ShapeType::Shape) {
        // OCCT: E.Orientation() == TopAbs_REVERSED -> Vertices(E, V2, V1)
        // (the argument swap).
        let (e1, e2) = if s.orientation == Orientation::Reversed {
            let (vf, vl) = top_exp_vertices_edge(&s);
            (vl, vf)
        } else {
            top_exp_vertices_edge(&s)
        };
        v1 = e1;
        v2 = e2;
        // OCCT: V1.Orientation(TopAbs_FORWARD); + add-or-remove in the map.
        // A null vertex carries no TShape identity to key on; LocOpe wires
        // always carry defined vertices.
        if let Some(v) = v1.take() {
            let v = with_orientation(&v, Orientation::Forward);
            v1 = Some(v);
            let stored = v1.clone().expect("just stored");
            vmap_toggle(&mut vmap, &mut vmap_order, stored);
        }
        // OCCT: V2.Orientation(TopAbs_REVERSED); + add-or-remove.
        if let Some(v) = v2.take() {
            let v = with_orientation(&v, Orientation::Reversed);
            v2 = Some(v);
            let stored = v2.clone().expect("just stored");
            vmap_toggle(&mut vmap, &mut vmap_order, stored);
        }
    }

    // OCCT: if (vmap.IsEmpty()) { closed — Vfirst/Vlast from V2 }.
    if vmap.is_empty() {
        if let Some(v2) = &v2 {
            v_first = Some(with_orientation(v2, Orientation::Forward));
            v_last = Some(with_orientation(v2, Orientation::Reversed));
        }
    } else if vmap.len() == 2 {
        // OCCT: open — first FORWARD key then first REVERSED key.
        for &key in &vmap_order {
            let k = vmap.get(&key).expect("vmap entry");
            if k.orientation == Orientation::Forward {
                v_first = Some(k.clone());
                break;
            }
        }
        for &key in &vmap_order {
            let k = vmap.get(&key).expect("vmap entry");
            if k.orientation == Orientation::Reversed {
                v_last = Some(k.clone());
                break;
            }
        }
    }
    (v_first, v_last)
}

/// OCCT LocOpe::Closed(const TopoDS_Wire& W, const TopoDS_Face& F)
/// (cxx L43-115) — Rust has no overloading: the `_wire` suffix carries the
/// OCCT parameter-type distinction.
pub fn closed_wire(w: &Shape, f: &Shape) -> bool {
    // OCCT cxx L45-50.
    let (v_f_opt, v_l_opt) = top_exp_vertices_wire(w);
    let v_f = v_f_opt.unwrap_or_else(Shape::null);
    let v_l = v_l_opt.unwrap_or_else(Shape::null);
    if !v_f.is_same(&v_l) {
        return false;
    }

    // On recherche l'edge contenant Vf FORWARD
    //
    // OCCT cxx L54-68: exp.Init(W.Oriented(TopAbs_FORWARD), TopAbs_EDGE) and
    // the exp2 vertex walk.
    let w_forward = with_orientation(w, Orientation::Forward);
    let mut ef: Option<Shape> = None;
    for cur_edg in explorer(&w_forward, ShapeType::Edge, ShapeType::Shape) {
        // OCCT cxx L57: exp2.Init(exp.Current(), TopAbs_VERTEX).
        let mut hit = false;
        for v in explorer(&cur_edg, ShapeType::Vertex, ShapeType::Shape) {
            // OCCT cxx L59: IsSame(Vf) && Orientation() == TopAbs_FORWARD.
            if v.is_same(&v_f) && v.orientation == Orientation::Forward {
                hit = true;
                break;
            }
        }
        // OCCT cxx L64-67: if (exp2.More()) break.
        if hit {
            ef = Some(cur_edg);
            break;
        }
    }
    // OCCT cxx L69: TopoDS_Edge Ef = TopoDS::Edge(exp.Current()) — the OCCT
    // source reads an exhausted explorer unprotected; the panic keeps that
    // surface.
    let ef = ef.expect("OCCT exp.Current() on exhausted explorer");

    // On recherche l'edge contenant Vl REVERSED
    //
    // OCCT cxx L73-87.
    let mut el: Option<Shape> = None;
    for cur_edg in explorer(&w_forward, ShapeType::Edge, ShapeType::Shape) {
        let mut hit = false;
        for v in explorer(&cur_edg, ShapeType::Vertex, ShapeType::Shape) {
            // OCCT cxx L77: IsSame(Vl) && Orientation() == TopAbs_REVERSED.
            if v.is_same(&v_l) && v.orientation == Orientation::Reversed {
                hit = true;
                break;
            }
        }
        if hit {
            el = Some(cur_edg);
            break;
        }
    }
    // OCCT cxx L87: TopoDS_Edge El = TopoDS::Edge(exp.Current()).
    let el = el.expect("OCCT exp.Current() on exhausted explorer");

    // OCCT cxx L89-108: pf/pl through BRep_Tool::CurveOnSurface.
    let (c2d, fc, lc) =
        brep_tool_curve_on_surface(&ef, f).expect("OCCT BRep_Tool::CurveOnSurface(Ef, F)");
    let pf: DVec2 = if ef.orientation == Orientation::Forward {
        c2d.point_at(fc)
    } else {
        c2d.point_at(lc)
    };
    let (c2d, fc, lc) =
        brep_tool_curve_on_surface(&el, f).expect("OCCT BRep_Tool::CurveOnSurface(El, F)");
    let pl: DVec2 = if el.orientation == Orientation::Forward {
        c2d.point_at(lc)
    } else {
        c2d.point_at(fc)
    };

    // OCCT cxx L110-114: pf.Distance(pl) <= Precision::PConfusion(Confusion).
    if (pf - pl).length()
        <= rcad_kernel::precision::p_confusion_with_tangent(rcad_kernel::precision::CONFUSION)
    {
        return true;
    }
    false
}

/// OCCT LocOpe::Closed(const TopoDS_Edge& E, const TopoDS_Face& F)
/// (cxx L119-126) — Rust has no overloading: the `_edge` suffix.
pub fn closed_edge(e: &Shape, f: &Shape) -> bool {
    // OCCT cxx L121-125: a one-edge wire through BRep_Builder (architecture
    // difference #5 — local pool).
    let mut pool = rcad_kernel::topo::topods::BRep::new();
    let mut b = rcad_kernel::topo::topods::BRepBuilder::new();
    let w = b.build_wire(&mut pool, vec![with_orientation(e, Orientation::Forward)]);
    closed_wire(&w, f)
}

/// OCCT LocOpe::TgtFaces(const TopoDS_Edge& E, const TopoDS_Face& F1,
/// const TopoDS_Face& F2) (cxx L130-191).
pub fn tgt_faces(e: &Shape, f1: &Shape, f2: &Shape) -> bool {
    // OCCT cxx L132: BRepAdaptor_Surface bs(F1, false) — the dead local of
    // the OCCT source is not translated (loc_ope_build_shape arch.
    // diff. #6 convention).
    // OCCT cxx L133-134.
    let mut u;
    let ta = 0.0001;

    // OCCT cxx L136: TopoDS_Edge e = E.
    let mut e = e.clone();

    // OCCT cxx L138-144: BRepAdaptor_Surface/Curve2d handles re-hosted
    // through the TFaceData surface and the edge pcurve (arch. diff. #3).
    e.orientation = Orientation::Forward;
    let (hc2d, _, _) =
        brep_tool_curve_on_surface(&e, f1).expect("OCCT BRepAdaptor_Curve2d(e, F1)");
    let (hc2d2, _, _) =
        brep_tool_curve_on_surface(&e, f2).expect("OCCT BRepAdaptor_Curve2d(e, F2)");

    // OCCT cxx L146-148.
    let rev1 = f1.orientation == Orientation::Reversed;
    let rev2 = f2.orientation == Orientation::Reversed;
    let (mut f, mut l) = brep_tool_range(&e);
    let mut angmin = std::f64::consts::PI;
    let mut angmax = -std::f64::consts::PI;
    let mut ang;

    // OCCT cxx L151-153.
    let eps = (l - f) / 100.;
    f += eps; // pour eviter de faire des calculs sur les
    l -= eps; // pointes des carreaux pointus.

    // OCCT cxx L154-159.
    let mut p: DVec2;
    let mut du: glam::DVec3;
    let mut dv: glam::DVec3;
    let mut d1: glam::DVec3;
    let mut d2: glam::DVec3;
    let mut uu;
    let mut vv;

    // OCCT cxx L162-189.
    for i in 0..=20i32 {
        u = f + (l - f) * (i as f64) / 20.;
        // OCCT cxx L165: HC2d->D0(u, p).
        p = hc2d.point_at(u);
        // OCCT cxx L166: HS1->D1(p.X(), p.Y(), pp1, du, dv).
        let s1 = match f1.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone().expect("F1 surface"),
            _ => panic!("Standard_NoSuchObject"),
        };
        let (_, du1, dv1) = s1.derivatives(p.x, p.y);
        du = du1;
        dv = dv1;
        // OCCT cxx L167-171.
        d1 = du.cross(dv).normalize_or_zero();
        if rev1 {
            d1 = -d1;
        }
        // OCCT cxx L172-175.
        p = hc2d2.point_at(u);
        uu = p.x;
        vv = p.y;
        let s2 = match f2.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone().expect("F2 surface"),
            _ => panic!("Standard_NoSuchObject"),
        };
        let (_, du2, dv2) = s2.derivatives(uu, vv);
        du = du2;
        dv = dv2;
        // OCCT cxx L176-179.
        d2 = du.cross(dv).normalize_or_zero();
        if rev2 {
            d2 = -d2;
        }
        // OCCT cxx L180-188.
        ang = d1.angle_between(d2);
        if ang <= angmin {
            angmin = ang;
        }
        if ang >= angmax {
            angmax = ang;
        }
    }
    // OCCT cxx L190.
    angmax <= ta
}

/// OCCT LocOpe::SampleEdges(const TopoDS_Shape& theShape,
/// NCollection_Sequence<gp_Pnt>& theSeq) (cxx L195-235).
pub fn sample_edges(the_shape: &Shape, the_seq: &mut Vec<glam::DVec3>) {
    // OCCT cxx L197: theSeq.Clear().
    the_seq.clear();
    // OCCT cxx L198: theMap (arch. diff. #1 — dedup only, never iterated).
    let mut the_map: HashMap<(u64, u32), ()> = HashMap::new();

    // Computes points on edge, but does not take the extremities into account
    //
    // OCCT cxx L207-225.
    for edg in explorer(the_shape, ShapeType::Edge, ShapeType::Shape) {
        // OCCT cxx L210-213: if (!theMap.Add(edg)) continue.
        if the_map.insert(shape_key(&edg), ()).is_none() {
            // OCCT cxx L214.
            if !brep_tool_degenerated(&edg) {
                // OCCT cxx L216-217: BRep_Tool::Curve + Transformed
                // (architecture difference #4).
                let Some((c, f, l)) = brep_tool_curve(&edg) else {
                    continue;
                };
                let c: Curve3 = c;
                // OCCT cxx L218-223.
                let delta = (l - f) / (NECHANT as f64) * 0.123456;
                for i in 1..NECHANT {
                    let prm =
                        delta + ((NECHANT - i) as f64 * f + i as f64 * l) / (NECHANT as f64);
                    the_seq.push(c.point_at(prm));
                }
            }
        }
    }

    // Adds every vertex
    //
    // OCCT cxx L228-234.
    for v in explorer(the_shape, ShapeType::Vertex, ShapeType::Shape) {
        if the_map.insert(shape_key(&v), ()).is_none() {
            the_seq.push(brep_tool_pnt(&v));
        }
    }
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
