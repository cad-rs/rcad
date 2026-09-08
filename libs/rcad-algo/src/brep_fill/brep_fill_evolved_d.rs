//! OCCT BRepFill_Evolved — 1:1 translation (part 4: the file statics).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_Evolved.cxx
//! - IsVertical              L219-238      - IsPlanar            L242-261
//! - Side                    L274-305      - IsInversed          L548-610
//! - ConcaveSide             L627-653      - Bubble              L1527-1544
//! - Compare                 L1825-1840    - AddDegeneratedEdge  L2555-2643
//! - TrimFace                L2647-2711    - PutProfilAt         L2715-2766
//! - TrimEdge                L2770-2864    - ComputeIntervals    L2868-2946
//! - Relative                L2954-3004    - PosOnFace           L3017-3043
//! - DoubleOrNotInFace       L3051-3083    - DistanceToOZ        L3087-3091
//! - Altitud                 L3095-3099    - SimpleExpression    L3103-3119
//! - CutEdgeProf             L3126-3267    - CutEdge             L3277-3381
//! - VertexFromNode          L3391-3443
//!
//! Plus the pure-math gp re-hosts (SetTransformation / SetRotation) and the
//! GAP carriers of the statics (Geom2dAPI_ExtremaCurveCurve,
//! BndLib_Add2dCurve, BRepLProp::Continuity — see the
//! brep_fill_evolved.rs header architecture difference #4).

use std::sync::Arc;

use glam::{DAffine3, DMat3, DVec2, DVec3};

use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, INTERSECTION};
use rcad_kernel::geom::{
    Circle3, Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Line3, Plane, TrimmedCurve2,
    TrimmedCurve3,
};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::gp::Ax3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};

use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_curve_on_surface, builder_add_edge_vertex, builder_range_edge,
    builder_set_degenerated, empty_copied, explorer,
};
use crate::brep_fill::brep_fill_evolved::{
    brep_fill_confusion, brep_tool_degenerated, brep_tool_pnt, brep_tool_tolerance,
    builder_make_vertex, edge_vertices, node_distance, node_infinite, wire_edges,
    DataMapOfShapeItem, LocationTable, MapNodeVertex,
};
use crate::brep_fill::brep_fill_trim_edge_tool::Geom2dAdaptorCurve;
use crate::brep_fill::brep_fill_trim_surface_tool::BRepFillTrimSurfaceTool;
use crate::brep_fill::generator::shape_oriented;
use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Curve2dType, TheIntPCurvePCurveOfGInter};
use crate::geomalgo::int_res2d::Domain as Res2dDomain;
use crate::topalgo::bisector::bisector_bisec::BisectorBisec;
use crate::topalgo::bisector::bisector_bisec_ana::{
    geom2d_curve_of, kind_of, BisectorBisecAna,
};
use crate::topalgo::bisector::bisector_curve::{CurveKind, Geom2dCurveHandle};
use crate::topalgo::mat::HandleMatNode;
use crate::topalgo::mat2d::mat2d_cut_curve::Mat2dCutCurve;

// ---------------------------------------------------------------------------
// GAP carriers + re-hosts (architecture difference #4)
// ---------------------------------------------------------------------------

/// GAP: Geom2dAPI_ExtremaCurveCurve (TKGeomAlgo/Geom2dAPI — the
/// Extrema_ExtCC2d engine has no rcad translation; the mat2d_mini_path.rs
/// ExtCC2d GAP precedent).  OCCT anchor:
/// Geom2dAPI_ExtremaCurveCurve.cxx L29-98.
pub(super) struct Geom2dAPIExtremaCurveCurve;

impl Geom2dAPIExtremaCurveCurve {
    /// OCCT Geom2dAPI_ExtremaCurveCurve(C1, C2, U1min, U1max, U2min, U2max).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        _the_c1: &Line2d,
        _the_c2: &Curve2d,
        _the_u1min: f64,
        _the_u1max: f64,
        _the_u2min: f64,
        _the_u2max: f64,
    ) -> Self {
        panic!(
            "GAP: Geom2dAPI_ExtremaCurveCurve (TKGeomAlgo/Geom2dAPI not translated) — \
             see file header"
        )
    }

    /// OCCT Geom2dAPI_ExtremaCurveCurve::NbExtrema().
    pub(super) fn nb_extrema(&self) -> i32 {
        panic!(
            "GAP: Geom2dAPI_ExtremaCurveCurve (TKGeomAlgo/Geom2dAPI not translated) — \
             see file header"
        )
    }

    /// OCCT Geom2dAPI_ExtremaCurveCurve::Parameters(Index, U1, U2).
    pub(super) fn parameters(&self, _index: i32, _u1: &mut f64, _u2: &mut f64) {
        panic!(
            "GAP: Geom2dAPI_ExtremaCurveCurve (TKGeomAlgo/Geom2dAPI not translated) — \
             see file header"
        )
    }
}

/// GAP: BndLib_Add2dCurve (TKMath/BndLib — not translated).  OCCT anchor:
/// BndLib_Add2dCurve.hxx.
pub(super) struct BndLibAdd2dCurve;

impl BndLibAdd2dCurve {
    /// OCCT BndLib_Add2dCurve::Add(C, Tol, B).
    pub(super) fn add(_the_c: &Geom2dAdaptorCurve, _the_tol: f64, _the_b: &mut BndBox2d) {
        panic!("GAP: BndLib_Add2dCurve (TKMath/BndLib not translated) — see file header")
    }
}

/// GAP: BRepLProp::Continuity(C1, C2, U1, U2) (TKTopAlgo/BRepLProp — not
/// translated).  OCCT anchor: BRepLProp.cxx L54-108.
pub(super) fn brep_lprop_continuity(_c1: &Shape, _c2: &Shape, _u1: f64, _u2: f64) -> i32 {
    panic!("GAP: BRepLProp::Continuity (TKTopAlgo/BRepLProp not translated) — see file header")
}

/// OCCT BRep_Builder::SameRange(E, B) — no-op (the brep_lib.rs stub
/// precedent; the rcad edge encodings carry the same_range flag from
/// construction).
pub(super) fn builder_same_range(_the_e: &Shape, _the_flag: bool) {}

/// OCCT BRep_Builder::SameParameter(E, B) — no-op (the brep_lib.rs stub
/// precedent).
pub(super) fn builder_same_parameter(_the_e: &Shape, _the_flag: bool) {}

/// OCCT myBuilder.MakeEdge(CurrentEdge, CBis, Tol) — the edge over a 3d
/// curve with a tolerance (the brep_algo::tool pool-free form).
pub(super) fn builder_make_edge_curve(the_c: &Curve3, the_tol: f64) -> Shape {
    let mut e = crate::brep_algo::tool::builder_make_edge();
    if let TShape::Edge(ed) = Arc::make_mut(&mut e.data) {
        ed.curve = Some(the_c.clone());
        ed.tolerance = ed.tolerance.max(the_tol);
    }
    e
}

/// OCCT BRepLib_MakeEdge(C, V1, V2) — the edge over a curve with the end
/// vertices (the range is the curve natural domain; the tolerance is the
/// kernel confusion — the BRepLib_MakeEdge defaults).
pub(super) fn builder_make_edge_curve_vertices(the_c: &Curve3, v1: &Shape, v2: &Shape) -> Shape {
    let mut e = builder_make_edge_curve(the_c, CONFUSION);
    let dom = CurveEval::default_domain(the_c);
    builder_add_edge_vertex(&mut e, &shape_oriented(v1, Orientation::Forward));
    builder_add_edge_vertex(&mut e, &shape_oriented(v2, Orientation::Reversed));
    builder_range_edge(&mut e, dom[0], dom[1]);
    e
}

/// OCCT BRepLib_MakeVertex(P) — the vertex at a point.
pub(super) fn brep_lib_make_vertex(the_p: DVec3) -> Shape {
    let mut v = builder_make_vertex();
    crate::brep_algo::tool::builder_update_vertex_point_tol(&mut v, the_p, CONFUSION);
    v
}

/// OCCT BRep_Builder::MakeVertex(V, P, Tol) — the vertex at a point with a
/// tolerance.
fn brep_lib_make_vertex_tol(the_p: DVec3, the_tol: f64) -> Shape {
    let mut v = builder_make_vertex();
    crate::brep_algo::tool::builder_update_vertex_point_tol(&mut v, the_p, the_tol);
    v
}

/// OCCT BRepLib_MakeEdge(C2d, S) — the pcurve-only edge (the TEdgeData
/// pcurve encoding bound on the owning face key).
fn builder_make_edge_pcurve_surface(
    the_c2d: &Curve2d,
    _the_s: &rcad_kernel::geom::Surface3,
    f: f64,
    l: f64,
    owner_face: &Shape,
) -> Shape {
    let mut e = crate::brep_algo::tool::builder_make_edge();
    if let TShape::Edge(ed) = Arc::make_mut(&mut e.data) {
        let key = (owner_face.ptr_id(), owner_face.location);
        ed.pcurves.insert(key, (the_c2d.clone(), f, l));
    }
    e
}

/// OCCT BRepTools::OriEdgeInFace(E, F) — the orientation of the edge as it
/// appears in the face's wires (Forward if absent; the pool-free reduction
/// of the fillet/chfi3d_builder_0.rs translation).
pub(super) fn brep_tools_ori_edge_in_face(e: &Shape, f: &Shape) -> Orientation {
    let fd = match f.data.as_ref() {
        TShape::Face(fd) => fd,
        _ => return Orientation::Forward,
    };
    for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
        if let TShape::Wire(wd) = w.data.as_ref() {
            for we in &wd.edges {
                if we.is_same(e) {
                    return we.orientation;
                }
            }
        }
    }
    Orientation::Forward
}

// ---------------------------------------------------------------------------
// Pure-math gp re-hosts (the draft.rs / loc_ope_pipe.rs precedent)
// ---------------------------------------------------------------------------

/// The gp_Ax3 frame affine (local -> world).
fn ax3_frame(a: &Ax3) -> DAffine3 {
    let m = DMat3::from_cols(a.x_direction, a.y_direction, a.axis.direction);
    let mut f = DAffine3::from_mat3(m);
    f.translation = a.axis.location;
    f
}

/// OCCT gp_Trsf::SetTransformation(A3) — the transformation from the A3
/// coordinate system to the absolute one (p' = Frame(A3)^-1 * p).
pub(super) fn trsf_set_transformation_ax3(a3: &Ax3) -> DAffine3 {
    ax3_frame(a3).inverse()
}

/// OCCT gp_Trsf::SetTransformation(FromT1, ToT2) — the transformation from
/// the T1 coordinate system to the T2 one (p' = Frame(T2)^-1 * Frame(T1) * p).
pub(super) fn trsf_set_transformation_between(from_t1: &Ax3, to_t2: &Ax3) -> DAffine3 {
    ax3_frame(to_t2).inverse() * ax3_frame(from_t1)
}

/// OCCT gp_Trsf::SetRotation(Axis, Angle) — the rotation about an axis
/// through the origin (the gp::OZ() forms consumed here).
pub(super) fn trsf_set_rotation(axis: &DVec3, angle: f64) -> DAffine3 {
    DAffine3::from_mat3(DMat3::from_axis_angle(axis.normalize_or_zero(), angle))
}

// ---------------------------------------------------------------------------
// File statics (BRepFill_Evolved.cxx)
// ---------------------------------------------------------------------------

/// OCCT static IsVertical (L219-238).
pub(super) fn is_vertical(e: &Shape) -> bool {
    let (v1, v2) = edge_vertices(e);
    let p1 = brep_tool_pnt(&v1);
    let p2 = brep_tool_pnt(&v2);

    if (p1.y - p2.y).abs() < brep_fill_confusion() {
        // It is a Line ?
        // OCCT L229-234: GC = BRep_Tool::Curve(E, Loc, f, l);
        // GC->DynamicType() == STANDARD_TYPE(Geom_Line).
        if let Some((gc, _f, _l)) = brep_tool_curve(e) {
            if matches!(gc, Curve3::Line(_)) {
                return true;
            }
        }
    }
    false
}

/// OCCT static IsPlanar (L242-261).
pub(super) fn is_planar(e: &Shape) -> bool {
    let (v1, v2) = edge_vertices(e);
    let p1 = brep_tool_pnt(&v1);
    let p2 = brep_tool_pnt(&v2);

    if (p1.z - p2.z).abs() < brep_fill_confusion() {
        // It is a Line ?
        // OCCT L252-257.
        if let Some((gc, _f, _l)) = brep_tool_curve(e) {
            if matches!(gc, Curve3::Line(_)) {
                return true;
            }
        }
    }
    false
}

/// OCCT static Side (L274-305) — determine the position of the profil
/// correspondingly to plane XOZ.
///           Return 1 : MAT_Left.
///           Return 2 : MAT_Left and Planar.
///           Return 3 : MAT_Left and Vertical.
///           Return 4 : MAT_Right.
///           Return 5 : MAT_Right and Planar.
///           Return 6 : MAT_Right and Vertical.
pub(super) fn side(profil: &Shape, tol: f64) -> i32 {
    // Rem : it is enough to test the first edge of the Wire.
    //       ( Correctly cut in PrepareProfil)
    // OCCT L279-282.
    let first_edge = wire_edges(profil)
        .into_iter()
        .next()
        .expect("TopExp_Explorer::Current (the profile wire is empty)");

    let (v1, v2) = edge_vertices(&first_edge);
    let p1 = brep_tool_pnt(&v1);
    let p2 = brep_tool_pnt(&v2);

    let mut the_side: i32;
    if p1.y < -tol || p2.y < -tol {
        the_side = 4;
    } else {
        the_side = 1;
    }
    if is_vertical(&first_edge) {
        the_side += 2;
    } else if is_planar(&first_edge) {
        the_side += 1;
    }
    the_side
}

/// OCCT static IsInversed (L548-610).
pub(super) fn is_inversed(s: &Shape, e1: &Shape, e2: &Shape, inverse: &mut [bool; 2]) {
    inverse[0] = false;
    inverse[1] = false;
    if s.shape_type() != ShapeType::Edge {
        return;
    }

    // OCCT L562-571: BRepAdaptor_Curve CS(TopoDS::Edge(S)); the D1 at the
    // traversal-first parameter.
    let (cs, cs_f, cs_l) = brep_tool_curve(s).expect("BRepAdaptor_Curve");
    let ds = if s.orientation == Orientation::Forward {
        CurveEval::derivative_at(&cs, cs_f)
    } else {
        -CurveEval::derivative_at(&cs, cs_l)
    };

    if !brep_tool_degenerated(e1) {
        let (c1, c1_f, c1_l) = brep_tool_curve(e1).expect("BRepAdaptor_Curve");
        let dc1 = if e1.orientation == Orientation::Forward {
            CurveEval::derivative_at(&c1, c1_f)
        } else {
            -CurveEval::derivative_at(&c1, c1_l)
        };
        inverse[0] = ds.dot(dc1) < 0.0;
    } else {
        inverse[0] = true;
    }

    if !brep_tool_degenerated(e2) {
        let (c2, c2_f, c2_l) = brep_tool_curve(e2).expect("BRepAdaptor_Curve");
        let dc2 = if e2.orientation == Orientation::Forward {
            CurveEval::derivative_at(&c2, c2_f)
        } else {
            -CurveEval::derivative_at(&c2, c2_l)
        };
        inverse[1] = ds.dot(dc2) < 0.0;
    } else {
        inverse[1] = true;
    }
}

/// OCCT static ConcaveSide (L627-653) — determine if the pipes were at the
/// side of the concavity.  WARNING: Not finished. Done only for circles.
pub(super) fn concave_side(s: &Shape, f: &Shape) -> bool {
    if s.shape_type() == ShapeType::Vertex {
        return false;
    }

    if s.shape_type() == ShapeType::Edge {
        // OCCT L638-641: G2d = BRep_Tool::CurveOnSurface(Edge(S), F, f, l);
        // Geom2dAdaptor_Curve AC(G2d, f, l).
        if let Some((g2d, gf, gl)) = brep_tool_curve_on_surface(s, f) {
            let ac = Geom2dAdaptorCurve::new(g2d, gf, gl);
            if Curve2dAdaptor::get_type(&ac) == Curve2dType::Circle {
                // OCCT L644-649: Direct = AC.Circle().IsDirect() — the
                // kernel Circle2d frame sense (x_dir x y_dir > 0) is the
                // direct encoding (architecture difference).
                let circle = Curve2dAdaptor::circle(&ac);
                let mut direct =
                    circle.x_dir.x * circle.y_dir.y - circle.x_dir.y * circle.y_dir.x > 0.0;
                if s.orientation == Orientation::Reversed {
                    direct = !direct;
                }
                return direct;
            }
        }
    }
    false
}

/// OCCT static Bubble (L1527-1544) — order the sequence of points by
/// growing x.
pub(super) fn bubble(seq: &mut Vec<f64>) {
    let mut invert = true;
    let nb_points = seq.len() as i32;

    while invert {
        invert = false;
        for i in 1..nb_points {
            // OCCT L1537: if (Seq.Value(i + 1) < Seq.Value(i)).
            if seq[i as usize] < seq[(i - 1) as usize] {
                seq.swap((i - 1) as usize, i as usize);
                invert = true;
            }
        }
    }
}

/// OCCT static Compare (L1825-1840).
pub(super) fn compare(e1: &Shape, e2: &Shape) -> Orientation {
    let mut oo = Orientation::Forward;
    let (v1_0, v1_1) = edge_vertices(e1);
    let (v2_0, v2_1) = edge_vertices(e2);
    let p1 = brep_tool_pnt(&v1_0);
    let p2 = brep_tool_pnt(&v2_0);
    let p3 = brep_tool_pnt(&v2_1);
    if p1.distance(p3) < p1.distance(p2) {
        oo = Orientation::Reversed;
    }
    oo
}

/// OCCT static AddDegeneratedEdge (L2555-2643) — degenerated edges can be
/// missing in some face; the missing degenerated edges have vertices
/// corresponding to node of the map.
pub(super) fn add_degenerated_edge(f: &Shape, w: &mut Shape) {
    // OCCT L2557-2572: the trimmed-plane / plane early returns.
    let s = match crate::brep_fill::brep_fill_evolved::brep_tool_surface(f) {
        Some(s) => s,
        None => return,
    };
    if let rcad_kernel::geom::Surface3::Trimmed(ts) = &s {
        if matches!(*ts.basis, rcad_kernel::geom::Surface3::Plane(_)) {
            return;
        }
    }
    if matches!(s, rcad_kernel::geom::Surface3::Plane(_)) {
        return;
    }

    // OCCT L2575.
    let tol_conf = 1.0e-4;

    let mut change = true;

    while change {
        change = false;
        // OCCT L2582-2584: BRepTools_WireExplorer WE(W, F).
        let we_edges = wire_edges(w);
        let mut pf = DVec2::ZERO;
        let mut prev_p = DVec2::ZERO;
        let mut vf = Shape::null();
        let mut v2_last = Shape::null();

        for we_current in &we_edges {
            let ce = we_current.clone();
            let (v1, v2) = edge_vertices(&ce);
            v2_last = v2.clone();
            // OCCT L2590-2597: BRep_Tool::UVPoints(CE, F, P2, P1) /
            // (CE, F, P1, P2).
            let (p1, p2) = crate::brep_algo::tool::brep_tool_uv_points(&ce, f);
            if vf.is_null() {
                vf = v1.clone();
                pf = p1;
            } else if p1.distance(prev_p) >= tol_conf {
                // degenerated edge to be inserted.
                change = true;
                let v = p1 - prev_p;
                let c2d = Curve2d::Line(Line2d::new(prev_p, v));
                let fpar = 0.0;
                let lpar = prev_p.distance(p1);
                // OCCT L2612: CT = new Geom2d_TrimmedCurve(C2d, f, l) — CT
                // is created but not consumed (the MakeEdge below takes
                // C2d) — the source as written.
                let _ct = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c2d.clone()),
                    t_min: fpar,
                    t_max: lpar,
                });
                // OCCT L2613-2618.
                let mut ne = builder_make_edge_pcurve_surface(&c2d, &s, fpar, lpar, f);
                builder_set_degenerated(&mut ne, true);
                builder_add_edge_vertex(&mut ne, &shape_oriented(&v1, Orientation::Forward));
                builder_add_edge_vertex(&mut ne, &shape_oriented(&v1, Orientation::Reversed));
                builder_range_edge(&mut ne, fpar, lpar);
                // OCCT L2618: B.Add(W, NE).
                crate::brep_fill::brep_fill_evolved::builder_add_wire_edge(w, &ne);
                break;
            }
            prev_p = p2;
        }
        // OCCT L2624-2641.
        if !change && vf.is_same(&v2_last) {
            // closed
            if pf.distance(prev_p) >= tol_conf {
                // Degenerated edge to be inserted.
                change = true;
                let v = pf - prev_p;
                let c2d = Curve2d::Line(Line2d::new(prev_p, v));
                let fpar = 0.0;
                let lpar = prev_p.distance(pf);
                let _ct = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c2d.clone()),
                    t_min: fpar,
                    t_max: lpar,
                });
                let mut ne = builder_make_edge_pcurve_surface(&c2d, &s, fpar, lpar, f);
                builder_set_degenerated(&mut ne, true);
                builder_add_edge_vertex(&mut ne, &shape_oriented(&vf, Orientation::Forward));
                builder_add_edge_vertex(&mut ne, &shape_oriented(&vf, Orientation::Reversed));
                builder_range_edge(&mut ne, fpar, lpar);
                // OCCT L2639: B.Add(W, NE).
                crate::brep_fill::brep_fill_evolved::builder_add_wire_edge(w, &ne);
            }
        }
    }
}

/// OCCT static TrimFace (L2647-2711).
pub(super) fn trim_face(face: &Shape, the_edges: &mut Vec<Shape>, s: &mut Vec<Shape>) {
    //--------------------------------------
    // Creation of wires limiting faces.
    //--------------------------------------
    let mut nb_edges;
    let mut new_wire;
    let mut add_edge;
    let mut good_wire = Shape::null();

    while !the_edges.is_empty() {
        // OCCT L2665-2669.
        let mut mwire = MakeWire::from_edge(&the_edges.first().cloned().expect("TheEdges.First()"));
        good_wire = mwire.wire();
        the_edges.remove(0);
        nb_edges = the_edges.len() as i32;
        new_wire = false;

        while !new_wire {
            add_edge = false;

            let mut i: i32 = 1;
            while i <= nb_edges && !add_edge {
                let e = the_edges[(i - 1) as usize].clone();
                if brep_tool_degenerated(&e) {
                    // OCCT L2678-2684: the degenerated edge is removed from
                    // the sequence without being added (the source as
                    // written).
                    the_edges.remove((i - 1) as usize);
                    add_edge = true;
                    nb_edges = the_edges.len() as i32;
                    good_wire = mwire.wire();
                } else {
                    // OCCT L2687-2696: MWire.Add(E); if the connection is
                    // successful the edge is removed from the sequence and
                    // one restarts from the beginning.
                    if mwire.add(&e) {
                        the_edges.remove((i - 1) as usize);
                        add_edge = true;
                        nb_edges = the_edges.len() as i32;
                        good_wire = mwire.wire();
                    }
                }
                i += 1;
            }
            new_wire = !add_edge;
        }
        // OCCT L2701-2709.
        let face_cut_tmp = empty_copied(face);
        let mut face_cut = shape_oriented(&face_cut_tmp, Orientation::Forward);
        // OCCT L2705: BRepTools::Update(FaceCut) — the tolerance
        // propagation walk (no rcad equivalent; the no-op stub precedent).
        add_degenerated_edge(&face_cut, &mut good_wire);
        crate::brep_fill::brep_fill_evolved::builder_add_face_wire(&mut face_cut, &good_wire);
        face_cut = shape_oriented(&face_cut, face.orientation);
        s.push(face_cut);
    }
}

/// OCCT BRepLib_MakeWire — the wire builder with the connection-error
/// semantics consumed by TrimFace (the brep_algo BRepLib_MakeWire GAP
/// reduction).
struct MakeWire {
    wire: Shape,
}

impl MakeWire {
    /// OCCT BRepLib_MakeWire(E).
    fn from_edge(e: &Shape) -> Self {
        let mut w = crate::brep_fill::brep_fill_evolved::builder_make_wire();
        crate::brep_fill::brep_fill_evolved::builder_add_wire_edge(&mut w, e);
        MakeWire { wire: w }
    }

    /// OCCT BRepLib_MakeWire::Wire().
    fn wire(&self) -> Shape {
        self.wire.clone()
    }

    /// OCCT BRepLib_MakeWire::Add(E) — true when the connection succeeded
    /// (Error() == BRepLib_WireDone): the edge shares a vertex with the
    /// wire and is not already in it.
    fn add(&mut self, e: &Shape) -> bool {
        let wire_edges_now = wire_edges(&self.wire);
        // Already in the wire (IsSame).
        for we in &wire_edges_now {
            if we.is_same(e) {
                return false;
            }
        }
        // Connectivity: the edge shares an extremity with the wire.
        let (ef, el) = edge_vertices(e);
        for we in &wire_edges_now {
            let (wf, wl) = edge_vertices(we);
            if ef.is_same(&wf) || ef.is_same(&wl) || el.is_same(&wf) || el.is_same(&wl) {
                crate::brep_fill::brep_fill_evolved::builder_add_wire_edge(&mut self.wire, e);
                return true;
            }
        }
        false
    }
}

/// OCCT static PutProfilAt (L2715-2766).
pub(super) fn put_profil_at(
    prof_ref: &Shape,
    axe_ref: &Ax3,
    e: &Shape,
    f: &Shape,
    at_start: bool,
    table: &mut LocationTable,
) -> Shape {
    // OCCT L2727-2731.
    let (c2d, first, last) = brep_tool_curve_on_surface(e, f).expect("CurveOnSurface");

    // OCCT L2733-2755: the D1 at the requested end (the orientation /
    // AtStart branch; the REVERSED branch reverses D1).
    let (p, d1) = if e.orientation == Orientation::Reversed {
        let (p, d) = if !at_start {
            d1_at(&c2d, first)
        } else {
            d1_at(&c2d, last)
        };
        (p, -d)
    } else if !at_start {
        d1_at(&c2d, last)
    } else {
        d1_at(&c2d, first)
    };

    // OCCT L2756-2761.
    let p3d = DVec3::new(p.x, p.y, 0.0);
    let v3d = DVec3::new(d1.x, d1.y, 0.0);

    let ax = Ax3::from_pnt_n_vx(p3d, DVec3::Z, v3d);
    let trans = trsf_set_transformation_between(&ax, axe_ref);

    // OCCT L2762-2765.
    let loc = table.intern(trans);
    crate::brep_fill::brep_fill_evolved::location_shape_moved(table, prof_ref, loc)
}

/// The Geom2d_Curve::D1(P, D1) reduction (the point and the derivative).
fn d1_at(c2d: &Curve2d, u: f64) -> (DVec2, DVec2) {
    (c2d.point_at(u), c2d.derivative_at(u))
}

/// OCCT static TrimEdge (L2770-2864).
pub(super) fn trim_edge(
    edge: &Shape,
    the_edges_controle: &Vec<Shape>,
    the_ver: &mut Vec<Shape>,
    the_par: &mut Vec<f64>,
    s: &mut Vec<Shape>,
) {
    let mut change = true;
    s.clear();
    //------------------------------------------------------------
    // Parse two sequences depending on the parameter on the edge.
    //------------------------------------------------------------
    // OCCT L2782-2794.
    while change {
        change = false;
        for i in 1..the_par.len() {
            // OCCT L2787: if (ThePar.Value(i) > ThePar.Value(i + 1)).
            if the_par[i - 1] > the_par[i] {
                the_par.swap(i - 1, i);
                the_ver.swap(i - 1, i);
                change = true;
            }
        }
    }

    //----------------------------------------------------------
    // If a vertex is not in the proofing point, it is removed.
    //----------------------------------------------------------
    // OCCT L2799-2810.
    if !brep_tool_degenerated(edge) {
        let mut k: i32 = 1;
        while k <= the_ver.len() as i32 {
            if double_or_not_in_face(the_edges_controle, &the_ver[(k - 1) as usize].clone()) {
                the_ver.remove((k - 1) as usize);
                the_par.remove((k - 1) as usize);
                k -= 1;
            }
            k += 1;
        }
    }

    //-------------------------------------------------------------------
    // Processing of double vertices for non-degenerated edges.
    // If a vertex_double appears twice in the edges of control,
    // the vertex is eliminated .
    // otherwise its only representation is preserved.
    //-------------------------------------------------------------------
    // OCCT L2818-2835.
    if !brep_tool_degenerated(edge) {
        let mut k: i32 = 1;
        while k < the_ver.len() as i32 {
            if the_ver[(k - 1) as usize].is_same(&the_ver[k as usize]) {
                the_ver.remove(k as usize);
                the_par.remove(k as usize);
                if double_or_not_in_face(the_edges_controle, &the_ver[(k - 1) as usize].clone()) {
                    the_ver.remove((k - 1) as usize);
                    the_par.remove((k - 1) as usize);
                    //	  k--;
                }
                k -= 1;
            }
            k += 1;
        }
    }

    //-----------------------------------------------------------
    // Creation of edges.
    // the number of vertices should be even. The edges to be created leave
    // from a vertex with uneven index i to vertex i+1;
    //-----------------------------------------------------------
    // OCCT L2842-2863.
    let mut k: i32 = 1;
    while k < the_ver.len() as i32 {
        let new_edge_tmp = empty_copied(edge);
        let mut new_edge = new_edge_tmp;

        if new_edge.orientation == Orientation::Reversed {
            builder_add_edge_vertex(
                &mut new_edge,
                &shape_oriented(&the_ver[(k - 1) as usize], Orientation::Reversed),
            );
            builder_add_edge_vertex(
                &mut new_edge,
                &shape_oriented(&the_ver[k as usize], Orientation::Forward),
            );
        } else {
            builder_add_edge_vertex(
                &mut new_edge,
                &shape_oriented(&the_ver[(k - 1) as usize], Orientation::Forward),
            );
            builder_add_edge_vertex(
                &mut new_edge,
                &shape_oriented(&the_ver[k as usize], Orientation::Reversed),
            );
        }
        builder_range_edge(
            &mut new_edge,
            the_par[(k - 1) as usize],
            the_par[k as usize],
        );
        //  modified by NIZHNY-EAP Wed Dec 22 12:09:48 1999 ___BEGIN___
        // OCCT L2860.
        crate::brep_fill::brep_fill_evolved::brep_lib_update_tolerances(&new_edge, false);
        //  modified by NIZHNY-EAP Wed Dec 22 13:34:19 1999 ___END___
        s.push(new_edge);
        k += 2;
    }
}

/// OCCT static ComputeIntervals (L2868-2946).
#[allow(clippy::too_many_arguments)]
pub(super) fn compute_intervals(
    v_on_f: &Vec<Shape>,
    v_on_l: &Vec<Shape>,
    par_on_f: &Vec<DVec3>,
    par_on_l: &Vec<DVec3>,
    trim: &BRepFillTrimSurfaceTool,
    bis: &Curve2d,
    vs: &Shape,
    ve: &Shape,
    first_par: &mut Vec<f64>,
    last_par: &mut Vec<f64>,
    first_v: &mut Vec<Shape>,
    last_v: &mut Vec<Shape>,
) {
    let mut i_on_f: i32 = 1;
    let mut i_on_l: i32 = 1;
    let mut u1: f64 = 0.0;
    let mut u2: f64;
    let mut v1 = Shape::null();
    let mut v2;

    if !vs.is_null() {
        u1 = first_parameter_of(bis);
        v1 = vs.clone();
    }
    while i_on_f <= v_on_f.len() as i32 || i_on_l <= v_on_l.len() as i32 {
        //---------------------------------------------------------
        // Return the smallest parameter on the bissectrice
        // corresponding to the current positions IOnF,IOnL.
        //---------------------------------------------------------
        // OCCT L2896-2909.
        if i_on_l > v_on_l.len() as i32
            || (i_on_f <= v_on_f.len() as i32
                && par_on_f[(i_on_f - 1) as usize].x < par_on_l[(i_on_l - 1) as usize].x)
        {
            u2 = par_on_f[(i_on_f - 1) as usize].x;
            v2 = v_on_f[(i_on_f - 1) as usize].clone();
            i_on_f += 1;
        } else {
            u2 = par_on_l[(i_on_l - 1) as usize].x;
            v2 = v_on_l[(i_on_l - 1) as usize].clone();
            i_on_l += 1;
        }
        //---------------------------------------------------------------------
        // When V2 and V1 are different the medium point P of the
        // interval is tested compared to the face. If P is in the face the interval
        // is valid.
        //---------------------------------------------------------------------
        // OCCT L2915-2925.
        if !v1.is_null() && !v2.is_same(&v1) {
            let p = bis.point_at((u2 + u1) * 0.5);
            if trim.is_on_face(p) {
                first_par.push(u1);
                last_par.push(u2);
                first_v.push(v1.clone());
                last_v.push(v2.clone());
            }
        }
        u1 = u2;
        v1 = v2;
    }

    // OCCT L2930-2945.
    if !ve.is_null() {
        u2 = last_parameter_of(bis);
        v2 = ve.clone();
        if !v2.is_same(&v1) {
            let p = bis.point_at((u2 + u1) * 0.5);
            if trim.is_on_face(p) {
                first_par.push(u1);
                last_par.push(u2);
                first_v.push(v1.clone());
                last_v.push(v2);
            }
        }
    }
}

/// OCCT static Relative (L2954-3004) — Commun is true if two wires have V
/// in common; returns FORWARD if the wires near the vertex are at the same
/// side, otherwise REVERSED.  (The rcad reduction carries the result in
/// the Commun flag — the orientation return value is consumed by no
/// caller in BRepFill_Evolved.cxx itself.)
pub(super) fn relative(w1: &Shape, w2: &Shape, v: &Shape, commun: &mut bool) {
    let mut e1 = Shape::null();
    let mut e2 = Shape::null();

    // OCCT L2963-2972.
    for exp_current in explorer(w1, ShapeType::Edge, ShapeType::Shape) {
        let e = exp_current;
        let (v1, v2) = edge_vertices(&e);
        if v1.is_same(v) || v2.is_same(v) {
            e1 = e;
            break;
        }
    }
    // OCCT L2973-2982.
    for exp_current in explorer(w2, ShapeType::Edge, ShapeType::Shape) {
        let e = exp_current;
        let (v1, v2) = edge_vertices(&e);
        if v1.is_same(v) || v2.is_same(v) {
            e2 = e;
            break;
        }
    }

    // OCCT L2984-2989.
    if e1.is_null() || e2.is_null() {
        *commun = false;
        return;
    }
    *commun = true;

    // OCCT L2991-3003.
    let ww1 = crate::brep_fill::brep_fill_evolved::brep_lib_make_wire_from_edge(&e1);
    let ww2 = crate::brep_fill::brep_fill_evolved::brep_lib_make_wire_from_edge(&e2);
    let tol = brep_fill_confusion();
    if side(&ww1, tol) < 4 && side(&ww2, tol) < 4 {
        // two to the left — TopAbs_FORWARD
        return;
    }
    if side(&ww1, tol) > 4 && side(&ww2, tol) > 4 {
        // two to the right — TopAbs_FORWARD
        return;
    }
    // TopAbs_REVERSED — carried by the Commun flag (see the header note).
}

/// OCCT static PosOnFace (L3017-3043).
pub(super) fn pos_on_face(d1: f64, d2: f64, d3: f64) -> i32 {
    if (d1 - d2).abs() <= brep_fill_confusion() {
        return 1;
    }
    if (d1 - d3).abs() <= brep_fill_confusion() {
        return 3;
    }

    if d2 < d3 {
        if d1 > (d2 + brep_fill_confusion()) && d1 < (d3 - brep_fill_confusion()) {
            return 2;
        }
    } else if d1 > (d3 + brep_fill_confusion()) && d1 < (d2 - brep_fill_confusion()) {
        return 2;
    }
    0
}

/// OCCT static DoubleOrNotInFace (L3051-3083) — returns True if V appears
/// zero or two times in the sequence of edges EC.
pub(super) fn double_or_not_in_face(ec: &Vec<Shape>, v: &Shape) -> bool {
    let mut vu = false;

    for i in 1..=ec.len() {
        let (v1, v2) = edge_vertices(&ec[i - 1]);
        if v1.is_same(v) {
            if vu {
                return true;
            } else {
                vu = true;
            }
        }
        if v2.is_same(v) {
            if vu {
                return true;
            } else {
                vu = true;
            }
        }
    }
    !vu
}

/// OCCT static DistanceToOZ (L3087-3091).
pub(super) fn distance_to_oz(v: &Shape) -> f64 {
    let pv3d = brep_tool_pnt(v);
    pv3d.y.abs()
}

/// OCCT static Altitud (L3095-3099).
pub(super) fn altitud(v: &Shape) -> f64 {
    let pv3d = brep_tool_pnt(v);
    pv3d.z
}

/// OCCT static SimpleExpression (L3103-3119) — the simple geometric
/// expression of the bissectrice (unwrap the trimmed BisecAna basis).
pub(super) fn simple_expression(b: &BisectorBisec, bis: &mut Option<Curve2d>) {
    // OCCT L3105: Bis = B.Value().
    let bis_trimmed = b.value();
    let trbis = bis_trimmed.read().expect("Bisector_Bisec::Value()");
    // OCCT L3107-3111: BT == Geom2d_TrimmedCurve; BasBis = TrBis->BasisCurve().
    let bas_bis = trbis.basis_curve.clone();
    // OCCT L3112-3113: BT = BasBis->DynamicType().
    let bt = bas_bis
        .as_any()
        .downcast_ref::<BisectorBisecAna>()
        .map(|_| CurveKind::BisecAna)
        .unwrap_or_else(|| {
            bas_bis
                .as_any()
                .downcast_ref::<Geom2dCurveHandle>()
                .map(|h| kind_of(&h.curve))
                .unwrap_or(CurveKind::Other)
        });
    if bt == CurveKind::BisecAna {
        // OCCT L3115: Bis = down_cast<Bisector_BisecAna>(BasBis)->Geom2dCurve().
        if let Some(c) = geom2d_curve_of(&bas_bis) {
            // OCCT L3116: Bis = new Geom2d_TrimmedCurve(Bis,
            // TrBis->FirstParameter(), TrBis->LastParameter()).
            *bis = Some(Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(c),
                t_min: trbis.first_parameter,
                t_max: trbis.last_parameter,
            }));
        }
    } else {
        // OCCT L3118 (implicit fall-through): Bis keeps the trimmed value —
        // the plain Geom2d curve is extracted from the handle basis; a
        // non-BisecAna bisector basis cannot be represented by the kernel
        // Curve2d enum (architecture difference).
        if let Some(h) = bas_bis.as_any().downcast_ref::<Geom2dCurveHandle>() {
            *bis = Some(h.curve.clone());
        }
    }
}

/// OCCT static CutEdgeProf (L3126-3267) — projection and Cut of an edge at
/// extrema of distance to axis OZ.
pub(super) fn cut_edge_prof(
    e: &Shape,
    plane: &Plane,
    line: &Line2d,
    cuts: &mut Vec<Shape>,
    map_ver_ref_moved: &mut DataMapOfShapeItem<Shape>,
) {
    cuts.clear();

    // OCCT L3142-3144: C = BRep_Tool::Curve(E, L, f, l);
    // CT = new Geom_TrimmedCurve(C, f, l); CT->Transform(L.Transformation())
    // — the rcad curve encodings carry the location (the offset_wire.rs
    // architecture difference).
    let (c_raw, f_raw, l_raw) = brep_tool_curve(e).expect("BRep_Tool::Curve");
    let ct = Curve3::Trimmed(TrimmedCurve3::new(c_raw, f_raw, l_raw));

    // project it in the plane and return the associated PCurve
    // OCCT L3147-3149: Normal = Plane->Pln().Axis().Direction();
    // C = GeomProjLib::ProjectOnPlane(CT, Plane, Normal, false);
    // C2d = GeomProjLib::Curve2d(C, Plane).
    let normal = plane.normal;
    let c = geom_proj_lib_project_on_plane(&ct, plane, normal);
    let c2d = crate::brep_fill::brep_fill_trim_surface_tool::geom_proj_lib_curve2d(&c, plane);

    // Calculate the extrema with the straight line
    // OCCT L3152-3163.
    let mut seq: Vec<f64> = Vec::new();

    let mut u1 = -INFINITE_VALUE;
    let mut u2 = INFINITE_VALUE;
    let f = first_parameter_of(&c2d);
    let l = last_parameter_of(&c2d);

    let mut b = BndBox2d::new();
    let ac2d = Geom2dAdaptorCurve::new(c2d.clone(), f, l);
    BndLibAdd2dCurve::add(&ac2d, brep_fill_confusion(), &mut b);
    // OCCT L3163: B.Get(xmin, U1, xmax, U2).
    if let Some((_, b_u1, _, b_u2)) = b.get() {
        u1 = b_u1;
        u2 = b_u2;
    }

    //  modified by NIZHNY-EAP Wed Feb  2 16:32:37 2000 ___BEGIN___
    // no sense if C2 is normal to Line or really is a point
    // OCCT L3167-3177.
    if u1 != u2 {
        let extrema = Geom2dAPIExtremaCurveCurve::new(line, &c2d, u1 - 1.0, u2 + 1.0, f, l);

        let nb = extrema.nb_extrema();
        for i in 1..=nb {
            let mut p_u1 = 0.0;
            let mut p_u2 = 0.0;
            extrema.parameters(i, &mut p_u1, &mut p_u2);
            seq.push(p_u2);
        }
    }
    //  modified by NIZHNY-EAP Wed Feb  2 16:33:05 2000 ___END___

    // On calcule les intersection avec Oy.
    // OCCT L3181-3191.
    let aline_c2d = Curve2d::Line(line.clone());
    let aline = Geom2dAdaptorCurve::new(
        aline_c2d.clone(),
        first_parameter_of(&aline_c2d),
        last_parameter_of(&aline_c2d),
    );
    let tol = INTERSECTION;
    let tol_c = 0.0;

    let mut intersector = TheIntPCurvePCurveOfGInter::new();
    let dom_line = bounded_domain(&aline_c2d, tol_c);
    let dom_ac2d = bounded_domain(&c2d, tol_c);
    intersector.perform(&aline, &dom_line, &ac2d, &dom_ac2d, tol_c, tol);
    let nb = intersector.base.nb_points();

    for i in 1..=nb {
        seq.push(intersector.base.point(i).param_on_second());
    }

    // Compute the new edges.
    // OCCT L3194-3216.
    let (v_rf, v_rl) = edge_vertices(e);

    let mut vf;
    let mut vl;

    if !map_ver_ref_moved.is_bound(&v_rf) {
        // OCCT L3200: Builder.MakeVertex(Vf, C->Value(f), Tolerance(VRf)) —
        // the Value(f) evaluates the projected 3d C at the 2d curve's first
        // parameter (the parametrizations coincide; OCCT source as
        // written).
        vf = brep_lib_make_vertex_tol(
            CurveEval::point_at(&c, f),
            brep_tool_tolerance(&v_rf),
        );
        map_ver_ref_moved.bind(&v_rf, vf.clone());
    } else {
        vf = map_ver_ref_moved.find(&v_rf).clone();
    }

    if !map_ver_ref_moved.is_bound(&v_rl) {
        vl = brep_lib_make_vertex_tol(CurveEval::point_at(&c, l), brep_tool_tolerance(&v_rl));
        map_ver_ref_moved.bind(&v_rl, vl.clone());
    } else {
        vl = map_ver_ref_moved.find(&v_rl).clone();
    }

    // OCCT L3218-3255.
    if !seq.is_empty() {
        bubble(&mut seq);

        let mut empty = false;

        let mut cur_param = f;

        while !empty {
            let param = seq.first().cloned().expect("Seq.First()");
            seq.remove(0);
            empty = seq.is_empty();
            if (param - cur_param).abs() > brep_fill_confusion()
                && (param - l).abs() > brep_fill_confusion()
            {
                let vv = brep_lib_make_vertex(CurveEval::point_at(&c, param));

                let ee_raw = builder_make_edge_curve_vertices(&c, &vf, &vv);
                let mut ee = shape_oriented(&ee_raw, e.orientation);
                let _ = &mut ee;
                if ee.orientation == Orientation::Forward {
                    cuts.push(ee);
                } else {
                    cuts.insert(0, ee);
                }

                // reinitialize
                cur_param = param;
                vf = vv;
            }
        }
    }

    // OCCT L3257-3266.
    let ee_raw = builder_make_edge_curve_vertices(&c, &vf, &vl);
    let ee = shape_oriented(&ee_raw, e.orientation);
    if ee.orientation == Orientation::Forward {
        cuts.push(ee);
    } else {
        cuts.insert(0, ee);
    }
}

/// OCCT GeomProjLib::ProjectOnPlane(C, Plane, Normal) — the projection of
/// a 3d curve on a plane (the analytic dispatch; a trimmed curve projects
/// its basis and keeps the range).  The non-analytic / oblique-circle
/// branches need the ProjLib approximation machinery — GAP (plan §0.6).
fn geom_proj_lib_project_on_plane(c: &Curve3, plane: &Plane, normal: DVec3) -> Curve3 {
    let n = normal.normalize_or_zero();
    let flatten = |p: DVec3| p - n * (p - plane.origin).dot(n);
    match c {
        Curve3::Line(l) => {
            // The projected line: flattened origin, flattened direction.
            let origin = flatten(l.origin);
            let dir = l.direction - n * l.direction.dot(n);
            Curve3::Line(Line3::new(origin, dir))
        }
        Curve3::Circle(ci) => {
            // OCCT: a circle parallel to the plane projects to a circle; an
            // oblique circle needs the ProjLib approximation (GAP).
            let axis = ci.normal.normalize_or_zero();
            if (axis - n).length() < 1.0e-9 || (axis + n).length() < 1.0e-9 {
                Curve3::Circle(Circle3 {
                    center: flatten(ci.center),
                    ..ci.clone()
                })
            } else {
                panic!(
                    "GAP: GeomProjLib::ProjectOnPlane oblique-circle branch needs ProjLib \
                     (TKGeomBase/ProjLib not translated) — see file header"
                )
            }
        }
        Curve3::Trimmed(t) => {
            let base = geom_proj_lib_project_on_plane(&t.curve, plane, normal);
            Curve3::Trimmed(TrimmedCurve3::new(base, t.first, t.last))
        }
        _ => panic!(
            "GAP: GeomProjLib::ProjectOnPlane non-analytic branch needs ProjLib \
             (TKGeomBase/ProjLib not translated) — see file header"
        ),
    }
}

/// The bounded IntRes2d_Domain of a 2d curve (the trim_edge_tool.rs
/// construction).
fn bounded_domain(c: &Curve2d, tol: f64) -> Res2dDomain {
    let dom = c.default_domain();
    let (f, l) = (dom[0], dom[1]);
    Res2dDomain::bounded(c.point_at(f), f, tol, c.point_at(l), l, tol)
}

/// OCCT static CutEdge (L3277-3381) — cut an edge at the extrema of curves
/// and at points of inflexion; closed circles are also cut in two.
pub(super) fn cut_edge(e: &Shape, f: &Shape, cuts: &mut Vec<Shape>) {
    cuts.clear();
    // OCCT L3280: MAT2d_CutCurve Cuter.
    let mut cuter = Mat2dCutCurve::new();

    // OCCT L3285-3290.
    let (v1, v2) = edge_vertices(e);
    let (c2d, f, l) = brep_tool_curve_on_surface(e, f).expect("BRep_Tool::CurveOnSurface");
    let ct2d = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(c2d.clone()),
        t_min: f,
        t_max: l,
    });

    // OCCT L3292: CT2d->BasisCurve()->IsKind(Geom2d_Circle) &&
    // BRep_Tool::IsClosed(E).
    let basis_is_circle = match &ct2d {
        Curve2d::Trimmed(t) => matches!(t.curve.as_ref(), Curve2d::Circle(_)),
        other => matches!(other, Curve2d::Circle(_)),
    };
    if basis_is_circle && crate::brep_algo::tool::brep_tool_is_closed_edge(e) {
        //---------------------------
        // Cut closed circle.
        //---------------------------
        // OCCT L3297-3337.
        let m1 = (2.0 * f + l) / 3.0;
        let m2 = (f + 2.0 * l) / 3.0;
        let p1 = ct2d.point_at(m1);
        let p2 = ct2d.point_at(m2);

        let vl1 = brep_lib_make_vertex(DVec3::new(p1.x, p1.y, 0.0));
        let vl2 = brep_lib_make_vertex(DVec3::new(p2.x, p2.y, 0.0));
        let fe_tmp = empty_copied(e);
        let me_tmp = empty_copied(e);
        let le_tmp = empty_copied(e);
        let mut fe = shape_oriented(&fe_tmp, Orientation::Forward);
        let mut me = shape_oriented(&me_tmp, Orientation::Forward);
        let mut le = shape_oriented(&le_tmp, Orientation::Forward);

        builder_add_edge_vertex(&mut fe, &v1);
        builder_add_edge_vertex(&mut fe, &shape_oriented(&vl1, Orientation::Reversed));
        builder_range_edge(&mut fe, f, m1);

        builder_add_edge_vertex(&mut me, &shape_oriented(&vl1, Orientation::Forward));
        builder_add_edge_vertex(&mut me, &shape_oriented(&vl2, Orientation::Reversed));
        builder_range_edge(&mut me, m1, m2);

        builder_add_edge_vertex(&mut le, &shape_oriented(&vl2, Orientation::Forward));
        builder_add_edge_vertex(&mut le, &v2);
        builder_range_edge(&mut le, m2, l);

        cuts.push(shape_oriented(&fe, e.orientation));
        cuts.push(shape_oriented(&me, e.orientation));
        cuts.push(shape_oriented(&le, e.orientation));
        //--------
        // Return.
        //--------
        return;
    }

    //-------------------------
    // Cut of the curve.
    //-------------------------
    // OCCT L3342.
    cuter.perform(&ct2d);

    // OCCT L3344-3350.
    if cuter.un_modified() {
        //-----------------------------
        // edge not modified => return.
        //-----------------------------
        return;
    }

    //------------------------
    // Creation of cut edges.
    //------------------------
    // OCCT L3352-3380.
    let mut vf = v1;

    for k in 1..=cuter.nb_curves() {
        let cc = cuter.value(k).clone();
        let vl = if k == cuter.nb_curves() {
            v2.clone()
        } else {
            let p = cc.point_at(last_parameter_of(&cc));
            brep_lib_make_vertex(DVec3::new(p.x, p.y, 0.0))
        };
        let ne_tmp = empty_copied(e);
        let mut ne = shape_oriented(&ne_tmp, Orientation::Forward);
        builder_add_edge_vertex(&mut ne, &shape_oriented(&vf, Orientation::Forward));
        builder_add_edge_vertex(&mut ne, &shape_oriented(&vl, Orientation::Reversed));
        builder_range_edge(&mut ne, first_parameter_of(&cc), last_parameter_of(&cc));
        cuts.push(shape_oriented(&ne, e.orientation));
        vf = vl;
    }
}

/// The Geom2d_TrimmedCurve First/LastParameter reduction (the natural
/// domain of a non-trimmed curve).
fn first_parameter_of(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Trimmed(t) => t.t_min,
        other => Curve2dEval::default_domain(other)[0],
    }
}

fn last_parameter_of(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Trimmed(t) => t.t_max,
        other => Curve2dEval::default_domain(other)[1],
    }
}

/// OCCT static VertexFromNode (L3391-3443) — test if the position of aNode
/// correspondingly to the distance to OZ of vertices VF and VL; returns
/// Status; if Status is different from 0 the vertex corresponding to aNode
/// is created.
pub(super) fn vertex_from_node(
    a_node: &HandleMatNode,
    e: &Shape,
    vf: &Shape,
    vl: &Shape,
    map_node_vertex: &mut MapNodeVertex,
    vn: &mut Shape,
) -> i32 {
    let mut shape_on_node = Shape::null();
    let mut status = 0;

    if !node_infinite(a_node) {
        status = pos_on_face(
            node_distance(a_node),
            distance_to_oz(vf),
            distance_to_oz(vl),
        );
    }
    if status == 2 {
        shape_on_node = e.clone();
    } else if status == 1 {
        shape_on_node = vf.clone();
    } else if status == 3 {
        shape_on_node = vl.clone();
    }

    if !shape_on_node.is_null() {
        //-------------------------------------------------
        // the vertex will correspond to a node of the map
        //-------------------------------------------------
        if map_node_vertex.node_shape_is_bound(a_node, &shape_on_node) {
            *vn = map_node_vertex.node_shape_find(a_node, &shape_on_node);
        } else {
            // OCCT L3434: B.MakeVertex(VN) — the unpointed vertex (the
            // point is set by the caller's UpdateVertex).
            *vn = builder_make_vertex();
            map_node_vertex.node_shape_bind(a_node, &shape_on_node, vn);
        }
    }
    status
}
