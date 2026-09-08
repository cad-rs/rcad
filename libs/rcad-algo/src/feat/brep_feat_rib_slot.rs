// OCCT BRepFeat_RibSlot.hxx L17-222 + BRepFeat_RibSlot.cxx L17-2736 +
// BRepFeat_RibSlot.lxx L17-24 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.lxx
//
// OCCT inheritance chain (BRepFeat_RibSlot.hxx L54):
//   BRepFeat_RibSlot : BRepBuilderAPI_MakeShape
//     BRepBuilderAPI_MakeShape : BRepBuilderAPI_Command
//
// NOTE: RibSlot does NOT derive from BRepFeat_Form in OCCT 8.0 (hxx L54) —
// the Form-like fields (myFShape/myLShape/myPerfSelection/myWire/mySbase/
// mySkface/myPbase/myGShape/mySUntil/myGluedF/myNewEdges/myTgtEdges/
// myFacesForDraft/myStatusError) are declared directly in RibSlot.hxx
// L197-217 and are therefore carried here as plain fields. Consequently the
// BRepFeat_Form pure-virtual slots (Curves/BarycCurve, Form.hxx L130-132)
// have NO counterpart on the OCCT RibSlot class and no coverage mapping is
// needed — the slots are consumed only by the BRepFeat_Form family
// (brep_feat_form.rs), not by RibSlot or its Rib/Slot subclasses
// (MakeLinearForm / MakeRevolutionForm, stage 3c).
//
// Rust has no inheritance -> composition + delegation. The MakeShape base
// sub-object is carried by the my_done/my_shape/my_generated fields (the
// same member-by-member mapping as brep_feat_form.rs). The OCCT protected
// members are pub(crate) so the Rib/Slot subclasses (stage 3c) can fill them
// exactly as the C++ subclasses drive the protected API (Add/Perform...).
//
// Field mapping (hxx L197-217 protected + L201-217 private):
//   BRepBuilderAPI_Command::myDone          -> my_done
//   BRepBuilderAPI_MakeShape::myShape       -> my_shape (None = null)
//   BRepBuilderAPI_MakeShape::myGenerated   -> my_generated
//   BRepFeat_RibSlot::myFirstPnt            -> my_first_pnt
//   BRepFeat_RibSlot::myLastPnt             -> my_last_pnt
//   BRepFeat_RibSlot::myFuse                -> my_fuse
//   BRepFeat_RibSlot::mySliding             -> my_sliding
//   BRepFeat_RibSlot::myMap                 -> my_map (DataMap<Shape,
//                                              List<Shape>>; the key shape is
//                                              carried as the tuple head —
//                                              the UpdateDescendants SkipFace
//                                              test reads the key ShapeType)
//   BRepFeat_RibSlot::myLFMap               -> my_lfmap (same carrier)
//   BRepFeat_RibSlot::myFShape              -> my_f_shape
//   BRepFeat_RibSlot::myLShape              -> my_l_shape
//   BRepFeat_RibSlot::myPerfSelection       -> my_perf_selection
//   BRepFeat_RibSlot::myWire                -> my_wire
//   BRepFeat_RibSlot::mySbase               -> my_sbase
//   BRepFeat_RibSlot::mySkface              -> my_skface
//   BRepFeat_RibSlot::myPbase               -> my_pbase
//   BRepFeat_RibSlot::myGShape              -> my_gshape
//   BRepFeat_RibSlot::mySUntil              -> my_suntil
//   BRepFeat_RibSlot::myGluedF              -> my_glued_f (DataMap<Shape,
//                                              Shape>; key shape = tuple head)
//   BRepFeat_RibSlot::myNewEdges            -> my_new_edges
//   BRepFeat_RibSlot::myTgtEdges            -> my_tgt_edges
//   BRepFeat_RibSlot::myFacesForDraft       -> my_faces_for_draft
//   BRepFeat_RibSlot::myStatusError         -> my_status_error
//
// Architecture differences (referenced from the affected functions):
// 1. The HugeSplit file: ExtremeFaces / PtOnEdgeVertex / SlidingProfile /
//    NoSlidingProfile live in brep_feat_rib_slot_b.rs (the 2000-line rule;
//    the impl block continues there — the OCCT class is one).
// 2. GeomAPI::To2d(C, Pln) (GeomAPI_*.cxx) is re-hosted below as
//    geom_api_to_2d on the rcad GeomProjLib vehicle
//    (geom_proj_lib::curve2d_auto) — the exact analytic pcurve for the
//    line/circle/ellipse curves consumed here, sampled BSpline otherwise.
// 3. Geom2dAPI_InterCurveCurve is re-hosted below as Geom2dAPIInterCurveCurve;
//    the line/line case runs the rcad IntAna2d analytic intersection
//    (AnaIntersection2d); the general conic/BSpline cases need
//    IntCurve_IntConicConic / Geom2dInt_GInter (TKGeomBase — not translated)
//    and stop at the GAP marker with the OCCT structure documented.
// 4. GeomAPI_ProjectPointOnCurve is re-hosted below as
//    GeomAPIProjectPointOnCurve over rcad ExtPC (the Initialize+Perform pair
//    of the OCCT constructor collapses into the rcad constructor — same
//    reduction as ParametricMinMax); Distance(N) =
//    sqrt(SquareDistance(N)) (GeomAPI_ProjectPointOnCurve.cxx L204-210).
// 5. BRepTopAdaptor_FClass2d(fac, tol) + Perform(uv, true) is re-hosted below
//    as brep_top_adaptor_fclass2d_perform on the rcad IntTools_FClass2d
//    vehicle (topalgo::brep_top_adaptor::fclass2d) over the standalone
//    FaceShapeSource adapter (face at index 0, identity location).
// 6. BRepIntCurveSurface_Inter is consumed through the
//    IntCurvesFaceIntersector re-host of loc_ope_cs_intersector.rs (the
//    reduced curve-surface vehicle; same GAP on the full
//    IntCurvesFace_Intersector semantics).
// 7. BRepFeat::IsInside (LFPerform L159) is the brep_feat_form_2 re-host
//    (GAP panic — BRepTopAdaptor_FClass2d +
//    GCPnts_QuasiUniformDeflection); the gluing branch keeps the OCCT
//    structure and stops at the same gap as BRepFeat_Form::GlobalPerform.
// 8. BRepAlgo::IsValid (SlidingProfile/NoSlidingProfile) is the
//    brep_feat_form_2 GAP re-host (BRepCheck_Analyzer on a standalone Shape
//    pending).
// 9. GeomLib::ExtendCurveToPoint (EdgeExtention) is re-hosted below as
//    geom_lib_extend_curve_to_point; the body needs
//    GeomConvert_CompCurveToBSplineCurve (TKGeomBase — not translated) and
//    stops at the GAP marker (OCCT L1269-1413 structure documented).
// 10. BRepLib::SameParameter(myShape, 1.e-7, true) (LFPerform L201) — the
//     shape-level BRepLib::SameParameter(S, Tol, EWithShared) edge-iteration
//     variant is not re-hosted (rcad BRepLib carries the per-edge stub); the
//     call spot carries the OCCT anchor and the deferred note.
// 11. BRepAlgoAPI_Cut (NoSlidingProfile L2641) is the CutVehicle of
//     brep_feat_form_2 (the (PaveFiller, Builder) vehicle with the history
//     knob on) — the same re-host as BRepFeat_Form's trP.
// 12. NCollection_DataMap/NCollection_Map map to HashMap keyed by
//     (TShape ptr, Location) — the TopTools_ShapeMapHasher identity; the
//     OCCT bucket iteration order is not reproduced (same reduction as
//     brep_feat_form.rs myMap).
// 13. BRep_Tool::Curve/Pnt/Tolerance/Degenerated/IsClosed are re-hosted
//     below from the TShape carriers (identity-location reduction of the
//     loc_ope modules); TopExp::FirstVertex/LastVertex are re-hosted below
//     from TopExp.cxx L182-210 (the CumOri composition follows the
//     TopoDS_Iterator rule with the rcad TEdgeData traversal pair).
// 14. BRepLib_MakeEdge/MakeWire/MakeFace/MakeVertex are re-hosted below as
//     small pool constructors (make_edge_*/make_vertex) — the rcad
//     add_tedge/add_tvertex/builders carry the OCCT semantics; the MakeEdge
//     (C, P1, P2) parameter projection uses GeomAPI_ProjectPointOnCurve
//     (BRepLib_MakeEdge.cxx command semantics).
// 15. The OCCT_DEBUG trace blocks are not translated (no debug tracing).
// 16. gp_Lin::Rotated(A1, Ang) (ChoiceOfFaces L672) is re-hosted below as
//     geom_line_rotated — the Rodrigues rotation of the gp_Ax1 frame.
// 17. BRepFeat_Builder::Perform (LFPerform L226) is the rcad
//     BRepFeatBuilder::perform_bop; Shape() is the result-pool compound
//     rebuild (builder_result_shape — the same reduction as
//     brep_feat_form.rs GlobalPerform L1220-1228).

use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder, OcctShapeMap};
use crate::feat::brep_feat_form_2::{
    brep_feat_is_inside, builder_result_shape, CutVehicle,
};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use crate::feat::loc_ope_find_edges::{
    elclib_parameter_circle, elclib_parameter_ellipse, elclib_parameter_lin, LocOpeFindEdges,
};
use crate::feat::loc_ope_gluer::LocOpeGluer;
use crate::feat::loc_ope_operation::LocOpeOperation;
use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;
use crate::topalgo::brep_top_adaptor::fclass2d::FClass2d;
use crate::topalgo::shape_source::FaceShapeSource;
use glam::{DAffine3, DQuat, DVec2, DVec3};
use rcad_kernel::base::extrema::ExtPC;
use rcad_kernel::base::geom_proj_lib::curve2d_auto;
use rcad_kernel::base::int_ana2d::AnaIntersection2d;
use rcad_kernel::geom::{
    Curve2d, Curve3, CurveEval, Ellipse3, Line3, Plane, Surface3, TrimmedCurve3,
};
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::topo::topods::{BRep, State, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Shared shape-identity helpers (the TopTools_ShapeMapHasher reduction).
// ---------------------------------------------------------------------------

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
pub(crate) fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::IsSame(S).
pub(crate) fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT NCollection_DataMap<Shape, List<Shape>>::IsBound.
pub(crate) fn data_map_is_bound(
    m: &HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    s: &Shape,
) -> bool {
    m.contains_key(&shape_key(s))
}

/// OCCT NCollection_DataMap::Bind(k, empty-list).
pub(crate) fn data_map_bind(m: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>, s: &Shape) {
    m.insert(shape_key(s), (s.clone(), Vec::new()));
}

/// OCCT NCollection_DataMap::operator() — the list of the bound entry
/// (Standard_NoSuchObject when unbound).
pub(crate) fn data_map_find<'a>(
    m: &'a HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    s: &Shape,
) -> &'a Vec<Shape> {
    match m.get(&shape_key(s)) {
        Some((_, l)) => l,
        None => panic!("Standard_NoSuchObject"),
    }
}

/// OCCT NCollection_DataMap::ChangeFind — the mutable list of the bound
/// entry (Standard_NoSuchObject when unbound).
pub(crate) fn data_map_change_find<'a>(
    m: &'a mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    key: (u64, u32),
) -> &'a mut Vec<Shape> {
    match m.get_mut(&key) {
        Some((_, l)) => l,
        None => panic!("Standard_NoSuchObject"),
    }
}

/// OCCT NCollection_Map::Add — returns true when newly added.
pub(crate) fn map_add(m: &mut OcctShapeMap, s: &Shape) -> bool {
    m.add(shape_key(s), s.clone())
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp re-hosts (architecture difference #13).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Curve(edg, Loc, f, l) — the 3D curve and parameter range
/// of an edge (the location transformation is the identity-location
/// reduction; the loc_ope_find_edges.rs architecture difference #1).
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// OCCT BRep_Tool::Pnt(vtx).
pub(crate) fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(shape) — the per-type tolerance carrier.
pub(crate) fn brep_tool_tolerance(the_s: &Shape) -> f64 {
    match the_s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
pub(crate) fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Surface(fac) with location applied (identity-location
/// reduction).
pub(crate) fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT TopExp::FirstVertex(E, CumOri) (TopExp.cxx L182-194): the iterator
/// exposes the edge's two stored vertices (composed with the occurrence
/// orientation when CumOri) and returns the composed-FORWARD one. The rcad
/// TEdgeData carries the traversal pair [first, last] (the OCCT
/// FORWARD/REVERSED child pair), so a REVERSED occurrence exposes its stored
/// last vertex as the composed-FORWARD child — the same convention as the
/// kernel BRepTool::oriented_first_vertex.
pub(crate) fn top_exp_first_vertex(e: &Shape, cum_ori: bool) -> Shape {
    let (v_first, v_last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    let mut v = if cum_ori && e.orientation == Orientation::Reversed {
        v_last
    } else {
        v_first
    };
    if cum_ori {
        v.orientation = e.orientation.compose(v.orientation);
    }
    v
}

/// OCCT TopExp::LastVertex(E, CumOri) (TopExp.cxx L198-210): the
/// composed-REVERSED counterpart of FirstVertex.
pub(crate) fn top_exp_last_vertex(e: &Shape, cum_ori: bool) -> Shape {
    let (v_first, v_last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => return Shape::null(),
    };
    let mut v = if cum_ori && e.orientation == Orientation::Reversed {
        v_first
    } else {
        v_last
    };
    if cum_ori {
        v.orientation = e.orientation.compose(v.orientation);
    }
    v
}

/// OCCT BRep_Tool::IsClosed(theShape) — the WIRE branch (BRep_Tool.cxx
/// L1707-1756): the vertices of the forward-oriented wire are parity-counted
/// through a map Add/Remove; closed when at least one vertex exists and
/// every vertex occurs an even number of times.
pub(crate) fn brep_tool_is_closed(the_shape: &Shape) -> bool {
    if the_shape.shape_type() == ShapeType::Wire {
        let mut oriented = the_shape.clone();
        oriented.orientation = Orientation::Forward;
        let mut a_map: HashMap<(u64, u32), Shape> = HashMap::new();
        let mut has_bound = false;
        for v in explorer(&oriented, ShapeType::Vertex, ShapeType::Shape) {
            if v.orientation == Orientation::Internal || v.orientation == Orientation::External {
                continue;
            }
            has_bound = true;
            // OCCT: if (!aMap.Add(V)) { aMap.Remove(V); } — the parity toggle.
            if a_map.remove(&shape_key(&v)).is_none() {
                a_map.insert(shape_key(&v), v.clone());
            }
        }
        has_bound && a_map.is_empty()
    } else {
        // The OCCT shape overload only defines the SHELL/WIRE branches; the
        // edge/vertex cases fall through the shell test and return the wire
        // test (both false for the vertices). RibSlot only consumes wires.
        false
    }
}

// ---------------------------------------------------------------------------
// Geometry re-hosts (architecture differences #2/#3/#4/#9/#16).
// ---------------------------------------------------------------------------

/// OCCT GeomAPI::To2d(C, P) — the pcurve of C on the plane (architecture
/// difference #2): the rcad GeomProjLib vehicle (exact analytic projection
/// for line/circle/ellipse curves, sampled BSpline otherwise).
pub(crate) fn geom_api_to_2d(the_c: &Curve3, the_pln: &Plane) -> Option<Curve2d> {
    curve2d_auto(the_c, &Surface3::Plane(*the_pln))
}

/// OCCT Geom_Line carrier of a Curve3::Line (the ln->Position().Direction() /
/// ln->Lin() readers).
pub(crate) fn geom_line_parts(the_c: &Curve3) -> (DVec3, DVec3) {
    match the_c {
        Curve3::Line(l) => (l.origin, l.direction),
        _ => panic!("Geom_Line carrier expected"),
    }
}

/// OCCT gp_Lin::Rotated(A1, Ang) over the Geom_Line carrier (architecture
/// difference #16): the position and the direction are rotated about the
/// axis (the gp_Trsf::SetRotation Rodrigues algebra).
fn geom_line_rotated(the_l: &Curve3, the_axe: &Ax1, the_ang: f64) -> Curve3 {
    let (origin, direction) = geom_line_parts(the_l);
    let axis = the_axe.direction.normalize_or_zero();
    let q = DQuat::from_axis_angle(axis, the_ang);
    Curve3::Line(Line3 {
        origin: the_axe.location + q * (origin - the_axe.location),
        direction: q * direction,
    })
}

/// OCCT Geom_Curve::Reversed() (Geom_Curve.cxx L31-36: Copy + Reverse) over
/// the rcad Curve3 variants. Per-class Reverse (Geom_Line.cxx L65-68 /
/// Geom_Conic.cxx L23-28 / Geom_TrimmedCurve.cxx L74-84):
/// - Line: the direction flips;
/// - conics (Circle/Ellipse/Hyperbola/Parabola): the position direction
///   flips (the X direction is kept, Y recomputed — u -> -u);
/// - Trimmed: the basis is reversed and the bounds swap
///   (Geom_TrimmedCurve::Reverse).
pub(crate) fn geom_curve_reversed(the_c: &Curve3) -> Curve3 {
    match the_c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin,
            direction: -l.direction,
        }),
        Curve3::Circle(c) => {
            let normal = -c.normal;
            Curve3::Circle(rcad_kernel::geom::Circle3 {
                center: c.center,
                normal,
                x_dir: c.x_dir,
                y_dir: normal.cross(c.x_dir),
                radius: c.radius,
            })
        }
        Curve3::Ellipse(e) => {
            let normal = -e.normal;
            Curve3::Ellipse(Ellipse3 {
                center: e.center,
                normal,
                major_dir: e.major_dir,
                major_radius: e.major_radius,
                minor_radius: e.minor_radius,
            })
        }
        Curve3::Trimmed(t) => {
            let basis = geom_curve_reversed(&t.curve);
            Curve3::Trimmed(TrimmedCurve3::new(basis, t.last, t.first))
        }
        Curve3::Hyperbola(h) => {
            let normal = -h.normal;
            Curve3::Hyperbola(rcad_kernel::geom::Hyperbola3 {
                center: h.center,
                normal,
                major_dir: h.major_dir,
                semi_major: h.semi_major,
                semi_minor: h.semi_minor,
            })
        }
        Curve3::Parabola(p) => {
            let normal = -p.normal;
            Curve3::Parabola(rcad_kernel::geom::Parabola3 {
                vertex: p.vertex,
                normal,
                axis_dir: p.axis_dir,
                focal_param: p.focal_param,
            })
        }
        _ => panic!("GAP(BRepFeat_RibSlot): Geom_Curve::Reversed for the BSpline/Bezier/Others variants needs the per-class Reverse translation"),
    }
}

/// OCCT ElCLib::Parameter(gp_Hyperbola, P) (ElCLib.cxx L1253-1261):
/// U = asinh((P - Loc) . YDirection / MinorRadius).
fn elclib_parameter_hyperbola(
    center: DVec3,
    y_dir: DVec3,
    minor_radius: f64,
    the_p: DVec3,
) -> f64 {
    ((the_p - center).dot(y_dir) / minor_radius).asinh()
}

/// OCCT ElCLib::Parameter(gp_Parabola, P) (ElCLib.cxx L1269-1274):
/// U = (P - Loc) . YDirection.
fn elclib_parameter_parabola(y_dir: DVec3, the_p: DVec3) -> f64 {
    the_p.dot(y_dir)
}

/// OCCT Geom2dAPI_InterCurveCurve re-host (architecture difference #3) —
/// the (C1, C2, Tol) constructor computes the intersection immediately.
pub(crate) struct Geom2dAPIInterCurveCurve {
    my_points: Vec<DVec2>,
}

impl Geom2dAPIInterCurveCurve {
    /// OCCT Geom2dAPI_InterCurveCurve(C1, C2, Tol). The line/line case is
    /// the IntAna2d analytic intersection (rcad AnaIntersection2d); the
    /// general conic/BSpline cases need IntCurve_IntConicConic /
    /// Geom2dInt_GInter (TKGeomBase — not translated): GAP.
    pub(crate) fn new(the_c1: &Curve2d, _the_c2: &Curve2d, _the_tol: f64) -> Self {
        let my_points = match (the_c1, _the_c2) {
            (Curve2d::Line(l1), Curve2d::Line(l2)) => {
                let mut a_int = AnaIntersection2d::new();
                a_int.perform_lin_lin(l1, l2);
                let mut pts = Vec::new();
                for n in 1..=a_int.nb_points() {
                    pts.push(a_int.point(n).value());
                }
                pts
            }
            _ => panic!(
                "GAP(BRepFeat_RibSlot): Geom2dAPI_InterCurveCurve needs IntCurve_IntConicConic / Geom2dInt_GInter (pending translation)"
            ),
        };
        Geom2dAPIInterCurveCurve { my_points }
    }

    /// OCCT NbPoints().
    pub(crate) fn nb_points(&self) -> i32 {
        self.my_points.len() as i32
    }

    /// OCCT Point(Index) — gp_Pnt2d.
    pub(crate) fn point(&self, the_index: i32) -> DVec2 {
        self.my_points[(the_index - 1) as usize]
    }
}

/// OCCT GeomAPI_ProjectPointOnCurve re-host (architecture difference #4) —
/// the Initialize(C, First, Last) + Perform(P) pair collapses into the rcad
/// ExtPC constructor.
pub(crate) struct GeomAPIProjectPointOnCurve {
    my_ext_pc: ExtPC,
}

impl GeomAPIProjectPointOnCurve {
    /// OCCT GeomAPI_ProjectPointOnCurve(P, Curve) (the .cxx Init body).
    pub(crate) fn new(the_p: DVec3, the_curve: &Curve3) -> Self {
        let dom = the_curve.default_domain();
        GeomAPIProjectPointOnCurve {
            my_ext_pc: ExtPC::new(the_p, the_curve, rcad_kernel::precision::CONFUSION, dom[0], dom[1]),
        }
    }

    /// OCCT NbPoints() — the ExtPC extrema count when done (the .cxx
    /// L161-172).
    pub(crate) fn nb_points(&self) -> i32 {
        if self.my_ext_pc.is_done() {
            self.my_ext_pc.nb_ext() as i32
        } else {
            0
        }
    }

    /// OCCT Distance(Index) — sqrt(SquareDistance(Index)) (the .cxx
    /// L204-210).
    pub(crate) fn distance(&self, the_index: i32) -> f64 {
        self.my_ext_pc.square_distance(the_index as usize).sqrt()
    }
}

/// OCCT BRepTopAdaptor_FClass2d(fac, tol) + Perform(uv, true) re-host
/// (architecture difference #5) — the rcad IntTools_FClass2d vehicle over
/// the standalone FaceShapeSource adapter (face at index 0, identity
/// location table).
pub(crate) fn brep_top_adaptor_fclass2d_perform(the_fac: &Shape, the_uv: DVec2) -> State {
    let Some(surf) = brep_tool_surface(the_fac) else {
        // The OCCT classifier is built on the face surface; a surface-less
        // face has no OCCT counterpart (the constructor would raise).
        panic!("BRepTopAdaptor_FClass2d: null face surface");
    };
    let locations = [DAffine3::IDENTITY];
    let src = FaceShapeSource::new(the_fac, surf, &locations);
    let a_cl = FClass2d::new(&src, 0, brep_tool_tolerance(the_fac));
    a_cl.perform(&src, the_uv, true)
}

/// OCCT BRepClass3d_SolidClassifier::State() — the rcad classifier code
/// (0 = IN, 1 = OUT, 2 = ON, see solid_classifier.rs state()) mapped to the
/// TopAbs_State.
fn top_abs_state(the_code: u8) -> State {
    match the_code {
        0 => State::In,
        2 => State::On,
        _ => State::Out,
    }
}

// ---------------------------------------------------------------------------
// BRepLib_MakeEdge / MakeVertex re-hosts (architecture difference #14).
// ---------------------------------------------------------------------------

/// OCCT BRepLib_MakeEdge(C, f, l): the vertices at C(f)/C(l), the range
/// [f, l] (BRepLib_MakeEdge.cxx Perform semantics).
pub(crate) fn make_edge_cl(the_brep: &mut BRep, the_c: &Curve3, the_f: f64, the_l: f64) -> Shape {
    let v1 = the_brep.add_tvertex(the_c.point_at(the_f));
    let v2 = the_brep.add_tvertex(the_c.point_at(the_l));
    the_brep.add_tedge(Some(the_c.clone()), v1, v2, [the_f, the_l])
}

/// OCCT BRepLib_MakeEdge(C, P1, P2): the vertex parameters come from the
/// projection of the points on C (BRepLib_MakeEdge.cxx — the
/// GeomAPI_ProjectPointOnCurve computation); the rcad add_tedge recomputes
/// the vertex parameters by position matching, the range carries the
/// projected parameters.
pub(crate) fn make_edge_c_p_p(
    the_brep: &mut BRep,
    the_c: &Curve3,
    the_p1: DVec3,
    the_p2: DVec3,
) -> Shape {
    let dom = the_c.default_domain();
    let par1 = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
        the_c,
        the_p1,
        dom[0],
        dom[1],
        64,
    )
    .param;
    let par2 = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
        the_c,
        the_p2,
        dom[0],
        dom[1],
        64,
    )
    .param;
    let v1 = the_brep.add_tvertex(the_p1);
    let v2 = the_brep.add_tvertex(the_p2);
    the_brep.add_tedge(Some(the_c.clone()), v1, v2, [par1, par2])
}

/// OCCT BRepLib_MakeEdge(V1, V2): the linear edge between two points
/// (BRepLib_MakeEdge.cxx — the Geom_Line carrier).
pub(crate) fn make_edge_p_p(the_brep: &mut BRep, the_p1: DVec3, the_p2: DVec3) -> Shape {
    let c = Curve3::Line(Line3::new(the_p1, the_p2 - the_p1));
    make_edge_c_p_p(the_brep, &c, the_p1, the_p2)
}

/// OCCT BRepLib_MakeEdge(V1, V2): the linear edge between two vertices
/// (BRepLib_MakeEdge.cxx — the Geom_Line carrier, parameters [0, |P2-P1|]).
pub(crate) fn make_edge_v_v(the_brep: &mut BRep, the_v1: &Shape, the_v2: &Shape) -> Shape {
    let p1 = brep_tool_pnt(the_v1);
    let p2 = brep_tool_pnt(the_v2);
    let c = Curve3::Line(Line3::new(p1, p2 - p1));
    the_brep.add_tedge(Some(c.clone()), the_v1.clone(), the_v2.clone(), [0.0, (p2 - p1).length()])
}

/// OCCT BRepLib_MakeEdge(C, V1, V2): the curve edge between two vertices —
/// the parameters come from the projection of the vertex points on C.
pub(crate) fn make_edge_c_v_v(
    the_brep: &mut BRep,
    the_c: &Curve3,
    the_v1: &Shape,
    the_v2: &Shape,
) -> Shape {
    let dom = the_c.default_domain();
    let p1 = brep_tool_pnt(the_v1);
    let p2 = brep_tool_pnt(the_v2);
    let par1 = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
        the_c, p1, dom[0], dom[1], 64,
    )
    .param;
    let par2 = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
        the_c, p2, dom[0], dom[1], 64,
    )
    .param;
    the_brep.add_tedge(Some(the_c.clone()), the_v1.clone(), the_v2.clone(), [par1, par2])
}

/// OCCT BRepLib_MakeVertex(P) — a fresh vertex (the pool constructor).
pub(crate) fn make_vertex(the_brep: &mut BRep, the_p: DVec3) -> Shape {
    the_brep.add_tvertex(the_p)
}

/// OCCT GeomLib::ExtendCurveToPoint(C, Point, Continuity, After) re-host
/// (architecture difference #9) — GAP: the body needs
/// GeomConvert_CompCurveToBSplineCurve + PLib::HermiteCoefficients
/// (GeomLib.cxx L1269-1413, pending translation); the OCCT structure: the
/// curve is converted (Convert_QuasiAngular), the Hermite constraint segment
/// of order Continuity is concatenated at the First/Last parameter side.
pub(crate) fn geom_lib_extend_curve_to_point(
    the_curve: &mut Curve3,
    the_point: DVec3,
    the_continuity: i32,
    the_after: bool,
) {
    let _ = (the_curve, the_point, the_continuity, the_after);
    panic!(
        "GAP(BRepFeat_RibSlot): GeomLib::ExtendCurveToPoint needs GeomConvert_CompCurveToBSplineCurve (pending translation)"
    );
}

// ---------------------------------------------------------------------------
// BRepFeat_RibSlot (hxx L54-218).
// ---------------------------------------------------------------------------

/// OCCT BRepFeat_RibSlot — provides functions to build mechanical features
/// (ribs and slots).
pub struct BRepFeatRibSlot {
    // --- BRepBuilderAPI_Command / MakeShape base ---
    my_done: bool,            // BRepBuilderAPI_Command::myDone
    my_shape: Option<Shape>,  // BRepBuilderAPI_MakeShape::myShape (None = null)
    my_generated: Vec<Shape>, // BRepBuilderAPI_MakeShape::myGenerated
    // --- BRepFeat_RibSlot members (hxx L197-217) ---
    pub(crate) my_first_pnt: DVec3,
    pub(crate) my_last_pnt: DVec3,
    pub(crate) my_fuse: bool,
    pub(crate) my_sliding: bool,
    // OCCT myMap: DataMap<Shape, List<Shape>>; the key shape is carried as
    // the tuple head (the UpdateDescendants(BOP) SkipFace test reads the key
    // ShapeType).
    pub(crate) my_map: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    pub(crate) my_lfmap: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    pub(crate) my_f_shape: Shape,
    pub(crate) my_l_shape: Shape,
    pub(crate) my_perf_selection: BRepFeatPerfSelection,
    pub(crate) my_wire: Shape,
    pub(crate) my_sbase: Shape,
    pub(crate) my_skface: Shape,
    pub(crate) my_pbase: Shape,
    pub(crate) my_gshape: Shape,
    pub(crate) my_suntil: Shape,
    // OCCT myGluedF: DataMap<Shape, Shape>; the key shape is carried as the
    // tuple head (the LFPerform iteration reads the Key and the Value).
    pub(crate) my_glued_f: HashMap<(u64, u32), (Shape, Shape)>,
    pub(crate) my_new_edges: Vec<Shape>,
    pub(crate) my_tgt_edges: Vec<Shape>,
    pub(crate) my_faces_for_draft: Vec<Shape>,
    pub(crate) my_status_error: BRepFeatStatusError,
}

impl BRepFeatRibSlot {
    /// OCCT BRepFeat_RibSlot::BRepFeat_RibSlot() (RibSlot.lxx L19-24).
    pub fn new() -> Self {
        BRepFeatRibSlot {
            my_done: false,
            my_shape: None,
            my_generated: Vec::new(),
            my_first_pnt: DVec3::ZERO,
            my_last_pnt: DVec3::ZERO,
            my_fuse: false,
            my_sliding: false,
            my_map: HashMap::new(),
            my_lfmap: HashMap::new(),
            my_f_shape: Shape::null(),
            my_l_shape: Shape::null(),
            my_perf_selection: BRepFeatPerfSelection::NoSelection,
            my_wire: Shape::null(),
            my_sbase: Shape::null(),
            my_skface: Shape::null(),
            my_pbase: Shape::null(),
            my_gshape: Shape::null(),
            my_suntil: Shape::null(),
            my_glued_f: HashMap::new(),
            my_new_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
            my_faces_for_draft: Vec::new(),
            my_status_error: BRepFeatStatusError::OK,
        }
    }

    // --- BRepBuilderAPI_Command / MakeShape base ---

    /// OCCT BRepBuilderAPI_Command::Done().
    fn done(&mut self) {
        self.my_done = true;
    }

    /// OCCT BRepBuilderAPI_Command::NotDone().
    fn not_done(&mut self) {
        self.my_done = false;
    }

    /// OCCT BRepBuilderAPI_Command::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepBuilderAPI_MakeShape::Shape().
    pub fn shape(&self) -> Option<&Shape> {
        self.my_shape.as_ref()
    }

    // -----------------------------------------------------------------
    // LFPerform (cxx L82-259) — topological reconstruction of ribs.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::LFPerform (cxx L82-259).
    pub fn lf_perform(&mut self) {
        // OCCT L89-100.
        if self.my_sbase.is_null()
            || self.my_pbase.is_null()
            || self.my_skface.is_null()
            || self.my_gshape.is_null()
            || self.my_lfmap.is_empty()
        {
            self.my_status_error = BRepFeatStatusError::NotInitialized;
            self.not_done();
            return;
        }

        // OCCT L102-103: TopExp_Explorer exp, exp2; int theOpe = 2.
        let mut the_ope = 2i32;

        // OCCT L105-108.
        if !self.my_glued_f.is_empty() {
            the_ope = 1;
        }

        // OCCT L110-128: hope that there is just a solid in the result (the
        // exp/exp2 exit state is not consumed afterwards — the OCCT quirk is
        // kept: the loops only run).
        if !self.my_suntil.is_null() {
            'until: for exp2 in explorer(&self.my_suntil, ShapeType::Face, ShapeType::Shape) {
                let funtil = exp2;
                for exp in explorer(&self.my_sbase, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&exp, &funtil) {
                        continue 'until;
                    }
                }
                break 'until;
            }
        }

        // OCCT L130-132: the list iterators; (int sens = 0 — commented out in
        // OCCT).
        //
        // OCCT L134: LocOpe_Gluer theGlue.
        let mut the_glue = LocOpeGluer::new();

        // --- case of gluing (OCCT L138-186) ---
        if the_ope == 1 {
            let mut collage = true;

            // OCCT L142-144: LocOpe_FindEdges theFE; the locmap DataMap of
            // the OCCT source is a dead local in LFPerform (declared at
            // L143-144, never consumed — the same reduction as the dead
            // locals of brep_feat_form.rs).
            let mut the_fe = LocOpeFindEdges::new();
            // OCCT L145.
            the_glue.init(&self.my_sbase, &self.my_gshape);
            // OCCT L146-175: (Key, Value) pairs of myGluedF — glface = Key,
            // fac = Value (the key shape is the tuple head).
            let glued_items: Vec<(Shape, Shape)> = self
                .my_glued_f
                .values()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect();
            for (glface, fac) in glued_items {
                let glface = glface;
                let fac = fac;
                // OCCT L150-156: find glface among the faces of myGShape.
                let mut found = false;
                for cur in explorer(&self.my_gshape, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, &glface) {
                        found = true;
                        break;
                    }
                }
                // OCCT L157-174.
                if found {
                    collage = brep_feat_is_inside(&glface, &fac);
                    if !collage {
                        the_ope = 2;
                        break;
                    } else {
                        the_glue.bind_face(&glface, &fac);
                        the_fe.set(&glface, &fac);
                        the_fe.init_iterator();
                        while the_fe.more() {
                            let ef = the_fe.edge_from();
                            let et = the_fe.edge_to();
                            the_glue.bind_edge(&ef, &et);
                            the_fe.next();
                        }
                    }
                }
            }

            // OCCT L177-185.
            let ope = the_glue.ope_type();
            if ope == LocOpeOperation::Invalid
                || (self.my_fuse && ope != LocOpeOperation::Fuse)
                || (!self.my_fuse && ope != LocOpeOperation::Cut)
                || !collage
            {
                the_ope = 2;
            }
        }

        // --- gluing is always applicable (OCCT L190-210) ---
        if the_ope == 1 {
            // OCCT L192: theGlue.Perform() — deferred dependency
            // (loc_ope_gluer.rs Perform note): the call is carried at its
            // spot; with the deferred body myGluer stays not-done and the
            // OCCT else branch below is the one taken.
            if the_glue.is_done() {
                // OCCT L195.
                self.update_descendants_gluer(&the_glue);
                // OCCT L196-197.
                self.my_new_edges = the_glue.edges().clone();
                self.my_tgt_edges = the_glue.tgt_edges().clone();
                //
                // OCCT L199-200.
                self.done();
                if let Some(shshs) = the_glue.resulting_shape() {
                    self.my_shape = Some(shshs.clone());
                }
                // OCCT L201: BRepLib::SameParameter(myShape, 1.e-7, true) —
                // the shape-level BRepLib::SameParameter(S, Tol, EWithShared)
                // edge-iteration variant is not re-hosted (architecture
                // difference #10); the rcad edges carry same_parameter on
                // construction.
            } else {
                the_ope = 2;
            }
        }

        // --- case without gluing (OCCT L213-258) ---
        if the_ope == 2 {
            // OCCT L215-219: BRepFeat_Builder theBuilder; partsoftool;
            // BRepClass3d_SolidClassifier oussa; bFlag; aIt.
            let mut the_builder = BRepFeatBuilder::new();
            let partsoftool: Vec<Shape>;
            let mut oussa = SolidClassifier::new();

            // OCCT L221.
            let b_flag = self.my_perf_selection != BRepFeatPerfSelection::NoSelection;
            //
            // OCCT L223-224.
            the_builder.init_with_tool(&self.my_sbase, &self.my_gshape);
            the_builder.set_operation_with_flag(self.my_fuse as i32, b_flag);
            //
            // OCCT L226.
            the_builder.perform_bop();
            // OCCT L227-252.
            if b_flag {
                // OCCT L229-230.
                partsoftool = the_builder.parts_of_tool();
                if !partsoftool.is_empty()
                    && self.my_perf_selection != BRepFeatPerfSelection::NoSelection
                {
                    // OCCT L233.
                    let toler = brep_tool_tolerance(&self.my_pbase) * 2.0;
                    //
                    // OCCT L235-247.
                    for it in partsoftool.clone() {
                        oussa.load(&it);
                        oussa.perform(self.my_first_pnt, toler);
                        let sp1 = top_abs_state(oussa.state());
                        oussa.perform(self.my_last_pnt, toler);
                        let sp2 = top_abs_state(oussa.state());
                        if sp1 != State::Out && sp2 != State::Out {
                            // OCCT L244-245: const TopoDS_Shape& S =
                            // aIt.Value(); theBuilder.KeepPart(S).
                            the_builder.keep_part(&it);
                        }
                    }
                }
                //
                // OCCT L250-251.
                the_builder.perform_result();
                self.my_shape = builder_result_shape(&mut the_builder);
            } else {
                // OCCT L255.
                self.my_shape = builder_result_shape(&mut the_builder);
            }
            // OCCT L257.
            self.done();
        }
    }

    // -----------------------------------------------------------------
    // IsDeleted (cxx L263-266).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::IsDeleted(F) (cxx L263-266). The OCCT
    /// myMap(F) Find raises Standard_NoSuchObject when unbound.
    pub fn is_deleted(&self, the_f: &Shape) -> bool {
        data_map_find(&self.my_map, the_f).is_empty()
    }

    // -----------------------------------------------------------------
    // Modified (cxx L270-293).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::Modified(F) (cxx L270-293) — returns the list
    /// of generated faces (the OCCT function-static list is returned by
    /// value).
    pub fn modified(&self, the_f: &Shape) -> Vec<Shape> {
        // OCCT L277-291.
        if data_map_is_bound(&self.my_map, the_f) {
            let mut list: Vec<Shape> = Vec::new();
            for sh in data_map_find(&self.my_map, the_f).clone() {
                if !shape_is_same(&sh, the_f) {
                    list.push(sh);
                }
            }
            return list;
        }
        self.my_generated.clone() // empty list
    }

    // -----------------------------------------------------------------
    // Generated (cxx L297-357).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::Generated(S) (cxx L297-357) — returns a list
    /// of the faces S created in the shape (the OCCT function-static list is
    /// returned by value).
    pub fn generated(&mut self, the_s: &Shape) -> Vec<Shape> {
        // OCCT L304.
        if the_s.shape_type() != ShapeType::Face {
            // OCCT L306.
            self.my_generated.clear();
            // OCCT L307.
            if self.my_lfmap.is_empty() || !data_map_is_bound(&self.my_lfmap, the_s) {
                // OCCT L309-323: check if filter on face or not.
                if data_map_is_bound(&self.my_map, the_s) {
                    let mut list: Vec<Shape> = Vec::new();
                    for sh in data_map_find(&self.my_map, the_s).clone() {
                        if !shape_is_same(&sh, the_s) {
                            list.push(sh);
                        }
                    }
                    return list;
                } else {
                    // OCCT L326.
                    return self.my_generated.clone();
                }
            } else {
                // OCCT L331.
                self.my_generated.clear();
                // OCCT L332-349.
                let mut list: Vec<Shape> = Vec::new();
                for it in data_map_find(&self.my_lfmap, the_s).clone() {
                    if data_map_is_bound(&self.my_map, &it) {
                        for sh in data_map_find(&self.my_map, &it).clone() {
                            if !shape_is_same(&sh, the_s) {
                                list.push(sh);
                            }
                        }
                    }
                }
                return list;
            }
        } else {
            // OCCT L355.
            return self.my_generated.clone();
        }
    }

    // -----------------------------------------------------------------
    // UpdateDescendants(const LocOpe_Gluer&) (cxx L361-386).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::UpdateDescendants(const LocOpe_Gluer& G)
    /// (cxx L361-386).
    fn update_descendants_gluer(&mut self, the_g: &LocOpeGluer) {
        // OCCT L368.
        let keys: Vec<(u64, u32)> = self.my_map.keys().copied().collect();
        for orig_key in keys {
            // OCCT L370-371: orig = itdm.Key(); newdsc map.
            let orig = match self.my_map.get(&orig_key) {
                Some((sh, _)) => sh.clone(),
                None => continue,
            };
            let mut newdsc = OcctShapeMap::new();
            // OCCT L372-379.
            let items = match self.my_map.get(&orig_key) {
                Some((_, l)) => l.clone(),
                None => Vec::new(),
            };
            for it in items {
                let fdsc = it;
                for it2 in the_g.descendant_faces(&fdsc).clone() {
                    map_add(&mut newdsc, &it2);
                }
            }
            // OCCT L380: myMap.ChangeFind(orig).Clear().
            let entry = data_map_change_find(&mut self.my_map, orig_key);
            entry.clear();
            // OCCT L381-384.
            for (_k, itm) in newdsc.iter() {
                let _ = _k;
                entry.push(itm.clone());
            }
            let _ = orig;
        }
    }

    // -----------------------------------------------------------------
    // Accessors (cxx L390-436).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::FirstShape (cxx L390-397). The absent-map case
    /// of the OCCT myMap(myFShape) Find (Standard_NoSuchObject) panics; the
    /// null myFShape returns myGenerated.
    pub fn first_shape(&self) -> &Vec<Shape> {
        if !self.my_f_shape.is_null() {
            return data_map_find(&self.my_map, &self.my_f_shape);
        }
        &self.my_generated // empty list
    }

    /// OCCT BRepFeat_RibSlot::LastShape (cxx L401-408).
    pub fn last_shape(&self) -> &Vec<Shape> {
        if !self.my_l_shape.is_null() {
            return data_map_find(&self.my_map, &self.my_l_shape);
        }
        &self.my_generated // empty list
    }

    /// OCCT BRepFeat_RibSlot::FacesForDraft (cxx L412-415).
    pub fn faces_for_draft(&self) -> &Vec<Shape> {
        &self.my_faces_for_draft
    }

    /// OCCT BRepFeat_RibSlot::NewEdges (cxx L419-422).
    pub fn new_edges(&self) -> &Vec<Shape> {
        &self.my_new_edges
    }

    /// OCCT BRepFeat_RibSlot::TgtEdges (cxx L426-429).
    pub fn tgt_edges(&self) -> &Vec<Shape> {
        &self.my_tgt_edges
    }

    /// OCCT BRepFeat_RibSlot::CurrentStatusError (cxx L433-436).
    pub fn current_status_error(&self) -> BRepFeatStatusError {
        self.my_status_error
    }

    // -----------------------------------------------------------------
    // CheckPoint (cxx L443-474) — proofing point material side.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::CheckPoint(e, bnd, Pln) (cxx L443-474). The
    /// bnd parameter is unnamed (unused) in OCCT.
    pub(crate) fn check_point(&self, e: &Shape, _bnd: f64, pln: &Plane) -> DVec3 {
        // OCCT L456-457: f, l; cc = BRep_Tool::Curve(e, f, l) — the
        // null-handle continuation is not defined in OCCT (deref).
        let Some((cc, f, l)) = brep_tool_curve(e) else {
            panic!("null 3D curve (OCCT null-handle deref)");
        };

        // OCCT L459-462.
        let par = (f + l) / 2.0;

        // OCCT L464: cc->D1(par, pp, tgt).
        let pp = cc.point_at(par);
        let tgt = cc.derivative_at(par);

        // OCCT L466-469.
        let mut tgt = tgt;
        if e.orientation == Orientation::Reversed {
            tgt = -tgt;
        }

        // OCCT L470-471: gp_Vec D = -tgt.Crossed(Pln->Pln().Position()
        // .Direction()) / 10.; pp.Translate(D).
        let d_vec = -(tgt.cross(pln.normal)) / 10.0;
        let pp = pp + d_vec;

        pp
    }

    // -----------------------------------------------------------------
    // Normal (cxx L481-529) — the normal to a face in a point.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::Normal(F, P) (cxx L481-529).
    pub(crate) fn normal(&self, the_f: &Shape, the_p: DVec3) -> DVec3 {
        // OCCT L489-490.
        let u: f64;
        let v: f64;

        // OCCT L492: BRepAdaptor_Surface AS(F, true) — the surface (the
        // rectangular trimmed wrapper is unwrapped once, the same reduction
        // as brep_feat_form.rs TransformShapeFU).
        let Some(mut s) = brep_tool_surface(the_f) else {
            panic!("BRepAdaptor_Surface: null face surface");
        };
        if let Surface3::Trimmed(t) = &s {
            s = (*t.basis).clone();
        }

        // OCCT L494-516: the per-type ElSLib::Parameters.
        match &s {
            Surface3::Plane(a_pl) => {
                let (pu, pv) = rcad_kernel::math::el::elslib_plane_parameters(
                    the_p, a_pl.origin, a_pl.u_dir, a_pl.v_dir,
                );
                u = pu;
                v = pv;
            }
            Surface3::Cylinder(a_cy) => {
                let x_dir = a_cy.ref_dir;
                let y_dir = a_cy.y_dir.unwrap_or_else(|| a_cy.axis.cross(x_dir));
                let (pu, pv) = rcad_kernel::math::el::elslib_cylinder_parameters(
                    the_p, a_cy.origin, x_dir, y_dir, a_cy.axis, a_cy.radius,
                );
                u = pu;
                v = pv;
            }
            Surface3::Cone(a_co) => {
                let x_dir = a_co.ref_dir;
                let y_dir = a_co.axis.cross(a_co.ref_dir);
                let (pu, pv) = rcad_kernel::math::el::elslib_cone_parameters(
                    the_p,
                    a_co.apex,
                    x_dir,
                    y_dir,
                    a_co.axis,
                    a_co.radius,
                    a_co.half_angle_rad,
                );
                u = pu;
                v = pv;
            }
            Surface3::Torus(a_to) => {
                let x_dir = a_to.ref_dir;
                let y_dir = a_to.axis.cross(a_to.ref_dir);
                let (pu, pv) = rcad_kernel::math::el::elslib_torus_parameters(
                    the_p,
                    a_to.center,
                    x_dir,
                    y_dir,
                    a_to.axis,
                    a_to.major_radius,
                    a_to.minor_radius,
                );
                u = pu;
                v = pv;
            }
            _ => {
                // OCCT L513-515: return gp_Dir(gp_Dir::D::X).
                return DVec3::X;
            }
        }

        // OCCT L518-520: AS.D1(U, V, pt, D1U, D1V).
        use rcad_kernel::SurfaceEval;
        let (_pt, d1u, d1v) = s.derivatives(u, v);

        // OCCT L521-523: CSLib::Normal(D1U, D1V, Precision::Confusion(),
        // St, N).
        let (n, _st) =
            rcad_kernel::math::cs_lib::normal_from_derivatives(d1u, d1v, rcad_kernel::precision::CONFUSION);
        let mut n = match n {
            Some(n) => n,
            None => return DVec3::X,
        };
        // OCCT L524-527.
        if the_f.orientation == Orientation::Forward {
            n = -n;
        }
        n
    }

    /// OCCT BRepFeat_RibSlot::IntPar (static) — the parameter of a point on a
    /// curve (cxx L536-575). The GeomAdaptor_Curve load of a Geom_TrimmedCurve
    /// carries the BASIS curve (GeomAdaptor_Curve.cxx L252-254), so the rcad
    /// Trimmed wrapper dispatches to the basis — the same per-type ElCLib
    /// parameter.
    pub fn int_par(the_c: &Curve3, the_p: DVec3) -> f64 {
        // OCCT L539-542: the null-handle check (the rcad reference carrier is
        // never null).
        let the_c = match the_c {
            Curve3::Trimmed(t) => t.curve.as_ref(),
            c => c,
        };
        match the_c {
            // OCCT L550-552: GeomAbs_Line.
            Curve3::Line(l) => elclib_parameter_lin(l, the_p),
            // OCCT L554-556: GeomAbs_Circle.
            Curve3::Circle(c) => elclib_parameter_circle(c, the_p),
            // OCCT L558-560: GeomAbs_Ellipse.
            Curve3::Ellipse(e) => elclib_parameter_ellipse(e, the_p),
            // OCCT L562-564: GeomAbs_Hyperbola.
            Curve3::Hyperbola(h) => {
                let x_dir = h.major_dir;
                let y_dir = h.normal.cross(x_dir);
                elclib_parameter_hyperbola(h.center, y_dir, h.semi_minor, the_p)
            }
            // OCCT L566-568: GeomAbs_Parabola.
            Curve3::Parabola(p) => {
                let x_dir = p.normal.cross(p.axis_dir);
                let y_dir = p.axis_dir.cross(x_dir);
                elclib_parameter_parabola(y_dir, the_p)
            }
            // OCCT L570-571: default — U = 0.
            _ => 0.0,
        }
    }

    // -----------------------------------------------------------------
    // EdgeExtention (cxx L582-638) — extension of an edge by tangence.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::EdgeExtention(e, bnd, FirstLast)
    /// (cxx L582-638). The OCCT `TopoDS_Edge E;` default-null local is the
    /// never-read initial assignment (kept with an allow).
    #[allow(unused_assignments)]
    pub(crate) fn edge_extention(&self, e: &mut Shape, bnd: f64, first_last: bool) {
        // OCCT L589-591: f, l; cu = BRep_Tool::Curve(e, f, l);
        // C = new Geom_TrimmedCurve(cu, f, l).
        let Some((cu, f, l)) = brep_tool_curve(e) else {
            panic!("null 3D curve (OCCT null-handle deref)");
        };
        let mut the_c = Curve3::Trimmed(TrimmedCurve3::new(cu.clone(), f, l));

        // OCCT L593: TopoDS_Edge E.
        let mut new_e = Shape::null();

        // OCCT L595-599: the analytic (conic) types keep their parameter
        // extension.
        let is_analytic = matches!(
            &cu,
            Curve3::Line(_)
                | Curve3::Circle(_)
                | Curve3::Ellipse(_)
                | Curve3::Hyperbola(_)
                | Curve3::Parabola(_)
        );

        let mut pool = BRep::new();

        if is_analytic {
            // OCCT L601-611.
            if first_last {
                new_e = make_edge_cl(&mut pool, &cu, f - bnd / 10.0, l);
            } else {
                new_e = make_edge_cl(&mut pool, &cu, f, l + bnd / 10.0);
            }
        } else {
            // OCCT L612-636: the non-analytic branch — the G1 extension to
            // the tangent-line point (see geom_lib_extend_curve_to_point).
            if first_last {
                // OCCT L620: C->D1(f, pnt, vct).
                let pnt = the_c.point_at(f);
                let vct = the_c.derivative_at(f);
                // OCCT L621-622: ln = new Geom_Line(pnt, -vct);
                // ln->D0(bnd / 1000., Pt).
                let ln_dir = -vct.normalize_or_zero();
                let pt = pnt + ln_dir * (bnd / 1000.0);
                // OCCT L623.
                geom_lib_extend_curve_to_point(&mut the_c, pt, 1, false);
                // OCCT L624-625.
                let p_last = brep_tool_pnt(&top_exp_last_vertex(e, true));
                new_e = make_edge_c_p_p(&mut pool, &the_c, pt, p_last);
            } else {
                // OCCT L629: C->D1(l, pnt, vct).
                let pnt = the_c.point_at(l);
                let vct = the_c.derivative_at(l);
                // OCCT L630-631: ln = new Geom_Line(pnt, vct);
                // ln->D0(bnd / 1000., Pt).
                let ln_dir = vct.normalize_or_zero();
                let pt = pnt + ln_dir * (bnd / 1000.0);
                // OCCT L632.
                geom_lib_extend_curve_to_point(&mut the_c, pt, 1, true);
                // OCCT L633-634.
                let p_first = brep_tool_pnt(&top_exp_first_vertex(e, true));
                new_e = make_edge_c_p_p(&mut pool, &the_c, p_first, pt);
            }
        }
        // OCCT L637: e = E.
        *e = new_e;
    }

    // -----------------------------------------------------------------
    // ChoiceOfFaces (static, cxx L645-704) — choose face of support in case
    // of support on an edge.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::ChoiceOfFaces(faces, cc, par, bnd, Pln)
    /// (cxx L645-704). The bnd parameter is unnamed (unused) in OCCT.
    pub fn choice_of_faces(
        faces: &mut Vec<Shape>,
        the_cc: &Curve3,
        par: f64,
        _bnd: f64,
        the_pln: &Plane,
    ) -> Shape {
        // OCCT L657: TopoDS_Face FFF.
        let mut fff = Shape::null();

        // OCCT L659-662.
        let pp = the_cc.point_at(par);
        let tgt = the_cc.derivative_at(par);

        // OCCT L664: l1 = new Geom_Line(pp, tgt).
        let l1 = Curve3::Line(Line3::new(pp, tgt));

        // OCCT L666-667: NCollection_Sequence<Geom_Curve> scur;
        // Counter = 0.
        let mut scur: Vec<Option<Curve3>> = Vec::new();
        let mut counter = 0i32;

        // OCCT L669: gp_Ax1 Axe(pp, Pln->Position().Direction()).
        let axe = Ax1 {
            location: pp,
            direction: the_pln.normal,
        };
        // OCCT L670-675.
        for i in 1..=8i32 {
            let rotated = geom_line_rotated(&l1, &axe, i as f64 * std::f64::consts::PI / 9.0);
            scur.push(Some(rotated));
            counter += 1;
        }

        // OCCT L677-701.
        let mut par_best = f64::MAX; // OCCT RealLast()
        for it in faces.clone() {
            let f = it;
            let mut asi = LocOpeCSIntersector::with_shape(&f);
            asi.perform_cur(&scur);
            if !asi.is_done() {
                continue;
            }
            for jj in 1..=counter {
                if asi.nb_points(jj) >= 1 {
                    let app = asi.point(jj, 1).parameter();
                    if app >= 0.0 && app < par_best {
                        par_best = app;
                        fff = f.clone();
                    }
                }
            }
        }

        // OCCT L703.
        fff
    }

    // -----------------------------------------------------------------
    // HeightMax (cxx L711-740).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::HeightMax(theSbase, theSUntil, p1, p2)
    /// (cxx L711-740).
    pub(crate) fn height_max(
        &self,
        the_sbase: &Shape,
        the_suntil: &Shape,
        p1: &mut DVec3,
        p2: &mut DVec3,
    ) -> f64 {
        // OCCT L721-726: Bnd_Box Box; BRepBndLib::Add(theSbase, Box)
        // [+ theSUntil] — the rcad shape_box re-host (the Bnd_Box union of
        // the two adds is the min/max component merge).
        let Some((bmin, bmax)) = BRepFeatBuilder::shape_box(the_sbase, &[]) else {
            panic!("Bnd_Box: empty bounding box (the OCCT Get returns the open box)");
        };
        let mut cmin = bmin;
        let mut cmax = bmax;
        if !the_suntil.is_null() {
            if let Some((umin, umax)) = BRepFeatBuilder::shape_box(the_suntil, &[]) {
                cmin = cmin.min(umin);
                cmax = cmax.max(umax);
            }
        }
        // OCCT L727-729: Box.Get(c[0], c[2], c[4], c[1], c[3], c[5]) —
        // c = [xmin, xmax, ymin, ymax, zmin, zmax]; bnd = c[0].
        let c: [f64; 6] = [
            cmin.x, cmax.x, cmin.y, cmax.y, cmin.z, cmax.z,
        ];
        let mut bnd = c[0];
        // OCCT L730-736.
        for i in 0..6 {
            if c[i] > bnd {
                bnd = c[i];
            }
        }
        // OCCT L737-738.
        *p1 = DVec3::new(c[0] - 2.0 * bnd, c[1] - 2.0 * bnd, c[2] - 2.0 * bnd);
        *p2 = DVec3::new(c[3] + 2.0 * bnd, c[4] + 2.0 * bnd, c[5] + 2.0 * bnd);
        // OCCT L739.
        bnd
    }

    // -----------------------------------------------------------------
    // UpdateDescendants(const BRepAlgoAPI_BooleanOperation&, ...)
    // (cxx L2672-2736).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::UpdateDescendants(aBOP, S, SkipFace)
    /// (cxx L2672-2736). The aBOP argument carries the CutVehicle of the
    /// performed cut (architecture difference #11).
    pub(crate) fn update_descendants_bop(
        &mut self,
        a_bop: &CutVehicle,
        the_s: &Shape,
        skip_face: bool,
    ) {
        // OCCT L2682.
        let keys: Vec<(u64, u32)> = self.my_map.keys().copied().collect();
        for orig_key in keys {
            // OCCT L2684.
            let orig = match self.my_map.get(&orig_key) {
                Some((sh, _)) => sh.clone(),
                None => continue,
            };
            // OCCT L2685-2688.
            if skip_face && orig.shape_type() == ShapeType::Face {
                continue;
            }
            // OCCT L2689.
            let mut newdsc = OcctShapeMap::new();

            // OCCT L2693-2720.
            let items = match self.my_map.get(&orig_key) {
                Some((_, l)) => l.clone(),
                None => Vec::new(),
            };
            for it in items {
                let sh = it;
                if sh.shape_type() != ShapeType::Face {
                    continue;
                }
                let fdsc = sh;
                // OCCT L2701-2708: preserved in S?
                let mut preserved = false;
                for cur in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, &fdsc) {
                        preserved = true;
                        map_add(&mut newdsc, &fdsc);
                        break;
                    }
                }
                // OCCT L2709-2719: if (!exp.More()).
                if !preserved {
                    let a_lm = a_bop.modified(&fdsc);
                    for a_it in a_lm {
                        map_add(&mut newdsc, &a_it);
                    }
                }
            }
            // OCCT L2721: myMap.ChangeFind(orig).Clear().
            let entry = data_map_change_find(&mut self.my_map, orig_key);
            entry.clear();
            // OCCT L2722-2734: check the belonging to the shape.
            for (_k, itm) in newdsc.iter() {
                let _ = _k;
                for cur in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, itm) {
                        entry.push(itm.clone());
                        break;
                    }
                }
            }
        }
    }
}

impl Default for BRepFeatRibSlot {
    fn default() -> Self {
        Self::new()
    }
}

// The remaining impl block (ExtremeFaces / PtOnEdgeVertex / SlidingProfile /
// NoSlidingProfile — the OCCT cxx L747-2668 body) lives in
// brep_feat_rib_slot_b.rs (the 2000-line rule; architecture difference #1).
