//! OCCT BRepOffset_Analyse (TKOffset/BRepOffset) — the 1:1 translation of
//! BRepOffset_Analyse.cxx (L1-1131) + BRepOffset_Analyse.hxx (L40-185).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_Analyse.cxx
//!         $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/BRepOffset_Analyse.hxx
//!
//! Architecture differences (the package numbering continues the
//! brep_offset_tool.rs list):
//!  57. NCollection_DataMap / IndexedDataMap / Map of TopoDS_Shape -> the
//!      shared ShapeDataMap / ShapeIndexedDataMap / OcctShapeSet forms of
//!      brep_offset_tool.rs (arch. diff. #21 cluster); the iteration of the
//!      indexed myAncestors map reproduces the OCCT 1-based operator() walk.
//!  58. The OCCT `mutable myDescendants` member (logical const) -> the
//!      RefCell field; `Descendants(S, theUpdate = false)` keeps the OCCT
//!      default-argument form as the 1-arg `descendants` plus the full
//!      `descendants_with_update` translation.
//!  59. BRep_Tool::Continuity(E, F1, F2) — rcad stores no continuity
//!      records; the order reads C0 (the chfi3d.rs is_tangent_faces L165
//!      precedent), so the `> GeomAbs_C0` branches are structurally
//!      present but not taken.
//!  60. BRep_Tool::IsClosed(E, F) (seam probe) -> the shape-level
//!      brep_tool_is_closed_edge_on_face re-host: rcad TEdgeData stores a
//!      single pcurve per face and no CurveOnClosedSurface record, the
//!      TShape CLOSED bit is the seam marker (BRep_Tool.cxx L795-846;
//!      L819-822: planes are never seam-closed).
//!  61. TopExp::MapShapesAndAncestors / MapShapesAndUniqueAncestors ->
//!      the local TopExp.cxx L80-175 re-hosts (no shared TopExp unit in
//!      rcad; the explorer-based translation of the OCCT bodies).
//!  62. BRepAdaptor_Surface(F, false).GetType() -> the Surface3 variant
//!      probe (the inter2d.rs BRepAdaptorSurface re-host precedent).
//!  63. Message_ProgressScope / Message_ProgressRange -> dropped: the
//!      public Perform carries the OCCT default-empty range, so the scopes
//!      never break; the `if (!aPS.More()) return;` guards are kept as
//!      comments (the inter3d.rs #37 convention).
//!  64. BOPTools_AlgoTools::MakeConnexityBlocks(S, VERTEX, EDGE, LCB, MVEMap)
//!      (BOPTools_AlgoTools.cxx L105-154) -> the
//!      crate::bop::algo::wire_splitter::make_connexity_blocks re-host
//!      (the L187-256 regularity-annotated re-host of the same core BFS,
//!      consumed by brep_offset_make_offset_1_f.rs for the same OCCT
//!      family); the block partition is the core BFS partition, and the
//!      theConnectionMap output (aMVEMap) is not consumed downstream.
//!  65. BOPTools_AlgoTools::IsSplitToReverse(Edge, Edge, theContext) ->
//!      crate::bop::tools::algo_tools::is_split_to_reverse_edge (the
//!      OCCT-context-free rcad form; theIntTools context is bound as a_ctx
//!      to keep the call structure).
//!  66. BRepPrimAPI_MakePrism (TKPrim) / BOPTools_AlgoTools3D::
//!      GetNormalToFaceOnEdge + PointNearEdge (TKBO/BOPTools, the rcad bop
//!      forms are BOPDS-index re-hosts private to the Builder) /
//!      LocalAnalysis_SurfaceContinuity (TKGeomAlgo) — GAP carriers; each
//!      keeps the OCCT call form and takes the OCCT failure path (the
//!      IsDone() == false stubs follow the chfi3d.rs LocalAnalysis
//!      precedent) until the units land.
//!  67. gp_Vec / gp_Dir / gp_Pnt / gp_Pnt2d -> DVec3 / DVec2; gp_Vec::
//!      IsOpposite / Angle -> the gp_Vec.hxx L130-134 + gp_Dir.cxx L27-50
//!      re-hosts; gp::Resolution() -> RealSmall() (DBL_MIN).

use std::cell::RefCell;
use std::collections::HashMap;

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, CurveEval as _, Surface3};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::tools::algo_tools::is_split_to_reverse_edge;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_parameter, brep_tool_range_on_face, brep_tool_tolerance,
};
use crate::geomalgo::gtests_stubs::GeomAbsShape;
use crate::fillet::chfi3d::{define_connect_type, is_tangent_faces};
use crate::fillet::chfi3d_builder_0::GeomAbsSurfaceType;
use crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity;

use super::brep_offset_make_offset_1::empty_compound;
use super::brep_offset_tool::{
    brep_tool_curve, edge_vertices, face_surface_of, oriented, set_add, set_contains,
    shape_data_map, shape_indexed_data_map, top_exp_vertices, OcctShapeSet, ShapeDataMap,
    ShapeIndexedDataMap,
};

// OCCT BRepOffset_Interval (TKOffset/BRepOffset/BRepOffset_Interval.hxx
// L30-70 + BRepOffset_Interval.lxx) — the same-package support class of the
// Analyse; translated here so the Analyse carries no GAP for a same-package
// dependency (the tool_d single-field carrier predates this unit).
/// OCCT BRepOffset_Interval (BRepOffset_Interval.hxx L32-67).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BRepOffsetInterval {
    f: f64,
    l: f64,
    r#type: ChFiDS_TypeOfConcavity,
}

impl BRepOffsetInterval {
    /// OCCT BRepOffset_Interval::BRepOffset_Interval() (lxx L24-29).
    pub fn new() -> Self {
        BRepOffsetInterval {
            f: 0.0,
            l: 0.0,
            r#type: ChFiDS_TypeOfConcavity::Other,
        }
    }

    /// OCCT BRepOffset_Interval::BRepOffset_Interval(U1, U2, Type) (lxx L31-36).
    pub fn with_params(u1: f64, u2: f64, the_type: ChFiDS_TypeOfConcavity) -> Self {
        BRepOffsetInterval {
            f: u1,
            l: u2,
            r#type: the_type,
        }
    }

    /// OCCT BRepOffset_Interval::First(U) (lxx L38-41).
    pub fn set_first(&mut self, u: f64) {
        self.f = u;
    }

    /// OCCT BRepOffset_Interval::Last(U) (lxx L43-46).
    pub fn set_last(&mut self, u: f64) {
        self.l = u;
    }

    /// OCCT BRepOffset_Interval::Type(T) (lxx L48-51).
    pub fn set_type(&mut self, t: ChFiDS_TypeOfConcavity) {
        self.r#type = t;
    }

    /// OCCT BRepOffset_Interval::First() (lxx L53-56).
    pub fn first(&self) -> f64 {
        self.f
    }

    /// OCCT BRepOffset_Interval::Last() (lxx L58-61).
    pub fn last(&self) -> f64 {
        self.l
    }

    /// OCCT BRepOffset_Interval::Type() (lxx L63-66).
    pub fn type_of(&self) -> ChFiDS_TypeOfConcavity {
        self.r#type
    }
}

// ---------------------------------------------------------------------------
// File-local math / gp re-hosts (pure tool functions; the OCCT anchors are
// inlined at each site).
// ---------------------------------------------------------------------------

/// OCCT gp::Resolution() = RealSmall() = DBL_MIN (gp.hxx L60,
/// Standard_Real.hxx RealSmall).
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT gp_Dir::Angle (gp_Dir.cxx L27-50) — reached via gp_Vec::Angle
/// (gp_Vec.hxx L216-221, which normalizes into gp_Dir).  The null-magnitude
/// raise of gp_Vec::Angle is guaranteed away by the OCCT callers (the
/// SquareMagnitude guards of TangentEdges run first).
fn gp_vec_angle(the_v: &DVec3, the_other: &DVec3) -> f64 {
    let coord = the_v.normalize();
    let other = the_other.normalize();
    let cosinus = coord.dot(other);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = coord.cross(other).length();
        if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT gp_Vec::IsOpposite (gp_Vec.hxx L130-134):
/// PI - Angle(theOther) <= theAngularTolerance.
fn gp_vec_is_opposite(the_v: &DVec3, the_other: &DVec3, the_angular_tolerance: f64) -> bool {
    let an_ang = std::f64::consts::PI - gp_vec_angle(the_v, the_other);
    an_ang <= the_angular_tolerance
}

/// OCCT TopoDS_Shape::Reverse() — TopAbs::Reverse flips FORWARD/REVERSED and
/// INTERNAL/EXTERNAL (TopAbs.hxx Reverse).
fn shape_reverse(mut the_s: Shape) -> Shape {
    the_s.orientation = match the_s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    };
    the_s
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp / BRepAdaptor shape-level re-hosts local to the file.
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Continuity(E, F1, F2) — rcad stores no continuity records
/// (arch. diff. #59): the order reads C0.
fn brep_tool_continuity_edge_face(
    _the_edge: &Shape,
    _the_face1: &Shape,
    _the_face2: &Shape,
) -> GeomAbsShape {
    GeomAbsShape::C0
}

/// OCCT GeomAbs_Shape enum rank (GeomAbs_Shape.hxx order: C0, G1, C1, G2,
/// C2, C3, CN) — the `>` comparisons of the cxx; the brep_fill/
/// compatible_wires.rs continuity_rank precedent.
fn geom_abs_shape_rank(the_order: GeomAbsShape) -> i32 {
    match the_order {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::G1 => 1,
        GeomAbsShape::C1 => 2,
        GeomAbsShape::G2 => 3,
        GeomAbsShape::C2 => 4,
        GeomAbsShape::C3 => 5,
        GeomAbsShape::CN => 6,
    }
}

/// OCCT BRep_Tool::IsClosed(E, F) (BRep_Tool.cxx L795-846) — the seam probe
/// (arch. diff. #60): the edge carries a CurveOnClosedSurface representation
/// on a non-plane surface; L819-822: IsPlane(S) -> false.  rcad TEdgeData
/// stores one pcurve per face, so the TShape CLOSED bit is the seam marker.
fn brep_tool_is_closed_edge_on_face(the_edge: &Shape, the_face: &Shape) -> bool {
    let Some(ed) = the_edge.as_edge() else {
        return false;
    };
    if ed.flags & tshape_flags::CLOSED == 0 {
        return false;
    }
    !matches!(face_surface_of(the_face), Some(Surface3::Plane(_)))
}

/// OCCT BRepTools::IsReallyClosed(E, F) (BRepTools.cxx L1204-1220).
fn brep_tools_is_really_closed(e: &Shape, f: &Shape) -> bool {
    if !brep_tool_is_closed_edge_on_face(e, f) {
        return false;
    }
    let mut nbocc = 0;
    for a_s in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        if a_s.is_same(e) {
            nbocc += 1;
        }
    }
    nbocc == 2
}

/// OCCT TopExp::LastVertex(E, CumOri = false) (TopExp.cxx) — the raw stored
/// second vertex (the top_exp_vertices re-host of brep_offset_tool.rs).
fn top_exp_last_vertex(e: &Shape) -> Shape {
    let (_v1, v2) = top_exp_vertices(e);
    v2
}

/// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) (TopExp.cxx L80-120)
/// (arch. diff. #61).
fn top_exp_map_shapes_and_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut ShapeIndexedDataMap<Vec<Shape>>,
) {
    // visit ancestors
    for anc in explorer(s, ta, ShapeType::Shape) {
        for a_s in explorer(&anc, ts, ShapeType::Shape) {
            // OCCT L97-101: FindIndex / Add(exs.Current(), empty) / Append(anc).
            let idx = match m.get_full(&bat::shape_key(&a_s)) {
                Some((i, _, _)) => i,
                None => shape_indexed_data_map::add(m, &a_s, Vec::new()),
            };
            shape_indexed_data_map::value_1_mut(m, idx + 1).push(anc.clone());
        }
    }
    // visit shapes not under ancestors (OCCT L106-116: TopExp_Explorer ex(S, TS, TA))
    for a_s in explorer(s, ts, ta) {
        if !m.contains_key(&bat::shape_key(&a_s)) {
            shape_indexed_data_map::add(m, &a_s, Vec::new());
        }
    }
}

/// OCCT TopExp::MapShapesAndUniqueAncestors(S, TS, TA, M, useOrientation =
/// false) (TopExp.cxx L124-175) (arch. diff. #61).
fn top_exp_map_shapes_and_unique_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut ShapeIndexedDataMap<Vec<Shape>>,
) {
    // visit ancestors
    for anc in explorer(s, ta, ShapeType::Shape) {
        for a_s in explorer(&anc, ts, ShapeType::Shape) {
            let idx = match m.get_full(&bat::shape_key(&a_s)) {
                Some((i, _, _)) => i,
                None => shape_indexed_data_map::add(m, &a_s, Vec::new()),
            };
            let a_list = shape_indexed_data_map::value_1_mut(m, idx + 1);
            // OCCT L146-156: check if anc already exists in the list
            // (useOrientation = false -> IsSame).
            if !a_list.iter().any(|x| x.is_same(&anc)) {
                a_list.push(anc.clone());
            }
        }
    }
    // visit shapes not under ancestors (OCCT L163-173)
    for a_s in explorer(s, ts, ta) {
        if !m.contains_key(&bat::shape_key(&a_s)) {
            shape_indexed_data_map::add(m, &a_s, Vec::new());
        }
    }
}

/// OCCT BRepAdaptor_Surface(F, false).GetType() (cxx L93-97) — the surface
/// kind probe over the face surface (arch. diff. #62).
fn brep_adaptor_surface_type(the_face: &Shape) -> GeomAbsSurfaceType {
    match face_surface_of(the_face) {
        Some(Surface3::Plane(_)) => GeomAbsSurfaceType::Plane,
        Some(Surface3::BSpline(_)) => GeomAbsSurfaceType::BSplineSurface,
        Some(Surface3::Bezier(_)) => GeomAbsSurfaceType::BezierSurface,
        Some(_) => GeomAbsSurfaceType::Other,
        None => GeomAbsSurfaceType::Other,
    }
}

// ---------------------------------------------------------------------------
// GAP carriers (other-package dependencies pending their translation units).
// ---------------------------------------------------------------------------

/// OCCT LocalAnalysis_SurfaceContinuity (TKGeomAlgo/LocalAnalysis) — GAP
/// carrier (arch. diff. #66).  The OCCT IsTangentFaces / CheckMixedContinuity
/// sample loops instantiate it per sample; the pending state reports
/// IsDone() == false for every sample, which per the OCCT flow makes the
/// validity count fall short and the mixed-concavity probe fail.
struct LocalAnalysisSurfaceContinuity;

impl LocalAnalysisSurfaceContinuity {
    #[allow(clippy::too_many_arguments)]
    fn new(
        _the_c2d1: &Curve2d,
        _the_c2d2: &Curve2d,
        _the_par: f64,
        _the_surf1: &Surface3,
        _the_surf2: &Surface3,
        _the_order: GeomAbsShape,
        _the_tol3d: f64,
        _the_tol0d: f64,
        _the_tolang1: f64,
        _the_tolang2: f64,
        _the_tolang3: f64,
    ) -> Self {
        LocalAnalysisSurfaceContinuity
    }

    fn is_done(&self) -> bool {
        false
    }

    /// OCCT LocalAnalysis_SurfaceContinuity::IsG1 — unreachable while the
    /// pending state reports IsDone() == false; the panic keeps the OCCT
    /// failure path.
    fn is_g1(&self) -> bool {
        panic!("GAP: LocalAnalysis_SurfaceContinuity::IsG1 (TKGeomAlgo not translated)");
    }

    /// OCCT LocalAnalysis_SurfaceContinuity::C0Value — see IsG1.
    fn c0_value(&self) -> f64 {
        panic!("GAP: LocalAnalysis_SurfaceContinuity::C0Value (TKGeomAlgo not translated)");
    }
}

/// OCCT BRepPrimAPI_MakePrism(S, V) (TKPrim/BRepPrimAPI) — GAP carrier
/// (arch. diff. #66).  The pending state reports IsDone() == false (the
/// chfi3d.rs LocalAnalysis precedent), so the OCCT
/// `if (!aMP.IsDone()) continue;` branch is taken and the sweep-dependent
/// tail stays inert until the TKPrim unit lands.
struct BRepPrimAPIMakePrism {
    my_shape: Shape,
    my_vec: DVec3,
}

impl BRepPrimAPIMakePrism {
    fn new(the_shape: &Shape, the_vec: DVec3) -> Self {
        BRepPrimAPIMakePrism {
            my_shape: the_shape.clone(),
            my_vec: the_vec,
        }
    }

    fn is_done(&self) -> bool {
        let _ = (&self.my_shape, &self.my_vec);
        false
    }

    /// OCCT BRepPrimAPI_MakePrism::Shape — unreachable while IsDone() is
    /// false; the panic keeps the OCCT failure path.
    fn shape(&self) -> Shape {
        panic!("GAP: BRepPrimAPI_MakePrism::Shape (TKPrim not translated)");
    }

    /// OCCT BRepPrimAPI_MakePrism::Generated(S) (TopTools_ListOfShape) —
    /// see Shape.
    fn generated(&self, _the_s: &Shape) -> Vec<Shape> {
        panic!("GAP: BRepPrimAPI_MakePrism::Generated (TKPrim not translated)");
    }
}

/// OCCT BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(E, F, DN)
/// (BOPTools_AlgoTools3D.cxx L351-376) — GAP carrier (arch. diff. #66);
/// the rcad bop form is the BOPDS-index re-host private to
/// bop/algo/builder.rs, pending the BOPTools_AlgoTools3D unit.
fn get_normal_to_face_on_edge(_the_e: &Shape, _the_f: &Shape, _the_dn: &mut DVec3) {
    panic!("GAP: BOPTools_AlgoTools3D::GetNormalToFaceOnEdge (TKBO/BOPTools not translated)");
}

/// OCCT BOPTools_AlgoTools3D::PointNearEdge(E, F, T, Tol, P2d, P)
/// (BOPTools_AlgoTools3D.cxx L525-612, the 6-parameter overload) — GAP
/// carrier (arch. diff. #66); the rcad bop form is the DS-bound re-host
/// of bop/algo/builder.rs, pending the BOPTools_AlgoTools3D unit.
#[allow(clippy::too_many_arguments)]
fn point_near_edge(
    _the_e: &Shape,
    _the_f: &Shape,
    _the_t: f64,
    _the_tol: f64,
    _the_p2d: &mut DVec2,
    _the_p: &mut DVec3,
) {
    panic!("GAP: BOPTools_AlgoTools3D::PointNearEdge (TKBO/BOPTools not translated)");
}

// ---------------------------------------------------------------------------
// Static functions of BRepOffset_Analyse.cxx.
// ---------------------------------------------------------------------------

// OCCT BRepOffset_Analyse.cxx L46-55
fn correct_orientation_of_tangent(tang_vec: &mut DVec3, a_vertex: &Shape, an_edge: &Shape) {
    let v_last = top_exp_last_vertex(an_edge);
    if a_vertex.is_same(&v_last) {
        // OCCT L53: TangVec.Reverse();
        *tang_vec = -*tang_vec;
    }
}

// OCCT BRepOffset_Analyse.cxx L149-272 (the L57-60 forward declaration)
fn check_mixed_continuity(
    the_edge: &Shape,
    the_face1: &Shape,
    the_face2: &Shape,
    the_ang_tol: f64,
) -> bool {
    let mut a_mixed_cont = false;
    // OCCT L155: GeomAbs_Shape aCurrOrder = BRep_Tool::Continuity(theEdge, theFace1, theFace2);
    let a_curr_order = brep_tool_continuity_edge_face(the_edge, the_face1, the_face2);
    if geom_abs_shape_rank(a_curr_order) > geom_abs_shape_rank(GeomAbsShape::C0) {
        // OCCT L158-159: BRep_Tool::Continuity always returns the minimal
        // continuity between the faces, so aCurrOrder > C0 means the faces
        // are tangent along the whole edge.
        return a_mixed_cont;
    }
    // OCCT L162: the C0 result cannot be trusted — the value is the default.
    let tol_c0 = 0.001f64.max(1.5 * brep_tool_tolerance(the_edge));

    // OCCT L165-168: double aFirst; double aLast; handles aC2d1, aC2d2.
    let a_c2d1: Option<(Curve2d, f64, f64)>;
    let a_c2d2: Option<(Curve2d, f64, f64)>;

    if !the_face1.is_same(the_face2)
        && brep_tool_is_closed_edge_on_face(the_edge, the_face1)
        && brep_tool_is_closed_edge_on_face(the_edge, the_face2)
    {
        // OCCT L173-186: find the edge in the face 1 — that occurrence has
        // the correct orientation.
        let mut an_edge_in_face1: Option<Shape> = None;
        let a_face1 = oriented(the_face1, Orientation::Forward);
        for an_edge in explorer(&a_face1, ShapeType::Edge, ShapeType::Shape) {
            if an_edge.is_same(the_edge) {
                an_edge_in_face1 = Some(an_edge);
                break;
            }
        }
        let an_edge_in_face1 = match an_edge_in_face1 {
            Some(e) => e,
            None => return a_mixed_cont,
        };

        a_c2d1 = brep_tool_curve_on_surface(&an_edge_in_face1, &a_face1);
        let a_face2 = oriented(the_face2, Orientation::Forward);
        let an_edge_in_face1 = shape_reverse(an_edge_in_face1);
        a_c2d2 = brep_tool_curve_on_surface(&an_edge_in_face1, &a_face2);
    } else {
        // OCCT L200-208: obtaining the pcurves of the edge on the two faces.
        a_c2d1 = brep_tool_curve_on_surface(the_edge, the_face1);
        // OCCT L202-207: for the case of a seam edge
        let mut ee = the_edge.clone();
        if the_face1.is_same(the_face2) {
            ee = shape_reverse(ee);
        }
        a_c2d2 = brep_tool_curve_on_surface(&ee, the_face2);
    }

    // OCCT L211-214: the (aFirst, aLast) pair is assigned by the last
    // CurveOnSurface call (the second one) before its null check.
    let (Some((a_c2d1, _, _)), Some((a_c2d2, a_first, a_last))) = (a_c2d1, a_c2d2) else {
        return a_mixed_cont;
    };

    // OCCT L216-218: obtaining the two surfaces of the adjacent faces.
    let a_surf1 = face_surface_of(the_face1);
    let a_surf2 = face_surface_of(the_face2);
    let (Some(a_surf1), Some(a_surf2)) = (a_surf1, a_surf2) else {
        return a_mixed_cont;
    };

    let a_nb_samples: i32 = 23;

    // OCCT L227-231: check for mixed concavity — convex in some regions,
    // concave in others.
    let a_delta = (a_last - a_first) / (a_nb_samples as f64 - 1.0);
    let mut a_has_convex = false;
    let mut a_has_concave = false;
    let mut a_nb_valid = 0;

    for i in 1..=a_nb_samples {
        let a_par = if i == a_nb_samples {
            a_last
        } else {
            a_first + (i as f64 - 1.0) * a_delta
        };

        // OCCT L237-247: LocalAnalysis_SurfaceContinuity aCont(...)
        let a_cont = LocalAnalysisSurfaceContinuity::new(
            &a_c2d1,
            &a_c2d2,
            a_par,
            &a_surf1,
            &a_surf2,
            GeomAbsShape::G1,
            0.001,
            tol_c0,
            the_ang_tol,
            the_ang_tol,
            the_ang_tol,
        );
        if !a_cont.is_done() {
            continue;
        }

        a_nb_valid += 1;

        if !a_cont.is_g1() && (!a_has_convex || !a_has_concave) {
            let an_angle = a_cont.c0_value();
            a_has_convex = a_has_convex || (an_angle > std::f64::consts::FRAC_PI_2 + the_ang_tol);
            a_has_concave = a_has_concave || (an_angle < std::f64::consts::FRAC_PI_2 - the_ang_tol);
        }
    }

    if a_nb_valid < a_nb_samples / 2 {
        return a_mixed_cont;
    }

    // OCCT L268-271: mixed connectivity — both convex and concave regions exist.
    a_mixed_cont = a_has_convex && a_has_concave;

    a_mixed_cont
}

// OCCT BRepOffset_Analyse.cxx L81-145 (static EdgeAnalyse; the brep context
// parameter is the rcad architecture form — the ChFi3d re-hosts take the
// BRep pool handle).
fn edge_analyse(
    brep: &BRep,
    e: &Shape,
    f1: &Shape,
    f2: &Shape,
    sin_tol: f64,
    li: &mut Vec<BRepOffsetInterval>,
) {
    // OCCT L87-88: BRep_Tool::Range(E, F1, f, l) — raises when the pcurve is
    // absent (Standard_NoSuchObject).
    let (f, l) = brep_tool_range_on_face(e, f1)
        .expect("BRep_Tool::Range: no pcurve of the edge on the face (Standard_NoSuchObject)");
    let mut i = BRepOffsetInterval::new();
    i.set_first(f);
    i.set_last(l);

    // OCCT L93-97: BRepAdaptor_Surface(F, false).GetType()
    let a_surf_type1 = brep_adaptor_surface_type(f1);

    let a_surf_type2 = brep_adaptor_surface_type(f2);

    let is_two_planes = a_surf_type1 == GeomAbsSurfaceType::Plane && a_surf_type2 == GeomAbsSurfaceType::Plane;

    let mut connect_type = ChFiDS_TypeOfConcavity::Other;

    if is_two_planes {
        // OCCT L103: then use only the strong condition
        // (BRep_Tool::Continuity(E, F1, F2) > GeomAbs_C0 — arch. diff. #59).
        if geom_abs_shape_rank(brep_tool_continuity_edge_face(e, f1, f2))
            > geom_abs_shape_rank(GeomAbsShape::C0)
        {
            connect_type = ChFiDS_TypeOfConcavity::Tangential;
        } else {
            connect_type = define_connect_type(brep, e, f1, f2, sin_tol, false);
        }
    } else {
        let is_two_splines = (a_surf_type1 == GeomAbsSurfaceType::BSplineSurface
            || a_surf_type1 == GeomAbsSurfaceType::BezierSurface)
            && (a_surf_type2 == GeomAbsSurfaceType::BSplineSurface
                || a_surf_type2 == GeomAbsSurfaceType::BezierSurface);
        let mut is_mixed_concavity = false;
        if is_two_splines {
            let an_ang_tol = 0.1;
            is_mixed_concavity = check_mixed_continuity(e, f1, f2, an_ang_tol);
        }

        if !is_mixed_concavity {
            if is_tangent_faces(brep, e, f1, f2, GeomAbsShape::G1) {
                // OCCT L128: the weak condition
                connect_type = ChFiDS_TypeOfConcavity::Tangential;
            } else {
                connect_type = define_connect_type(brep, e, f1, f2, sin_tol, false);
            }
        } else {
            connect_type = ChFiDS_TypeOfConcavity::Mixed;
        }
    }

    i.set_type(connect_type);
    li.push(i);
}

// OCCT BRepOffset_Analyse.cxx L276-284 (static BuildAncestors)
fn build_ancestors(s: &Shape, ma: &mut ShapeIndexedDataMap<Vec<Shape>>) {
    ma.clear();
    top_exp_map_shapes_and_unique_ancestors(s, ShapeType::Vertex, ShapeType::Edge, ma);
    top_exp_map_shapes_and_unique_ancestors(s, ShapeType::Edge, ShapeType::Face, ma);
}

// ---------------------------------------------------------------------------
// The BRepOffset_Analyse class (BRepOffset_Analyse.hxx L41-185).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Analyse — analyses the shape to find the parts of edges
/// connecting the convex, concave or tangent faces.
///
/// `Clone` models the OCCT `const BRepOffset_Analyse* myAnalyzer` pointer
/// storage (BRepOffset_MakeOffset_1 SetAnalysis): the analyse is immutable
/// after Perform, so a copy is read-equivalent to the OCCT alias (the
/// mutable myDescendants cache recomputes lazily in each copy).
#[derive(Clone)]
pub struct BRepOffsetAnalyse {
    // Inputs
    my_shape: Shape,              //< OCCT: myShape — input shape to analyze
    my_angle: f64,                //< OCCT: myAngle — criteria angle to check tangency

    my_offset: f64,               //< OCCT: myOffset — offset value
    my_face_offset_map: ShapeDataMap<f64>, //< OCCT: myFaceOffsetMap

    // Results
    my_done: bool,                //< OCCT: myDone — status of the algorithm

    my_map_edge_type: ShapeDataMap<Vec<BRepOffsetInterval>>, //< OCCT: myMapEdgeType
    my_ancestors: ShapeIndexedDataMap<Vec<Shape>>,           //< OCCT: myAncestors
    my_replacement: ShapeDataMap<ShapeDataMap<Shape>>,       //< OCCT: myReplacement
    my_descendants: RefCell<ShapeDataMap<Vec<Shape>>>,       //< OCCT: mutable myDescendants

    my_new_faces: Vec<Shape>,     //< OCCT: myNewFaces
    my_generated: ShapeDataMap<Shape>, //< OCCT: myGenerated

    /// rcad architecture field (not in OCCT): the BRep pool handle consumed
    /// by the ChFi3d re-hosts (define_connect_type / is_tangent_faces).
    my_brep: BRep,
}

impl BRepOffsetAnalyse {
    /// OCCT BRepOffset_Analyse::BRepOffset_Analyse() (cxx L64-68).
    pub fn new() -> Self {
        BRepOffsetAnalyse {
            my_shape: Shape::null(),
            my_angle: 0.0,
            my_offset: 0.0,
            my_face_offset_map: ShapeDataMap::new(),
            my_done: false,
            my_map_edge_type: ShapeDataMap::new(),
            my_ancestors: ShapeIndexedDataMap::new(),
            my_replacement: ShapeDataMap::new(),
            my_descendants: RefCell::new(ShapeDataMap::new()),
            my_new_faces: Vec::new(),
            my_generated: ShapeDataMap::new(),
            my_brep: BRep::new(),
        }
    }

    /// OCCT BRepOffset_Analyse::BRepOffset_Analyse(S, Angle) (cxx L72-77).
    pub fn with_shape(s: &Shape, angle: f64) -> Self {
        let mut an_analyse = BRepOffsetAnalyse::new();
        an_analyse.perform(s, angle);
        an_analyse
    }

    /// OCCT BRepOffset_Analyse::Perform(S, Angle, theRange = empty)
    /// (cxx L288-367).
    pub fn perform(&mut self, s: &Shape, angle: f64) {
        self.my_shape = s.clone();
        self.my_new_faces.clear();
        self.my_generated.clear();
        self.my_replacement.clear();
        self.my_descendants.borrow_mut().clear();

        self.my_angle = angle;
        let sin_tol = angle.sin().abs();

        // OCCT L302: Build ancestors.
        build_ancestors(s, &mut self.my_ancestors);

        // OCCT L304-308: aLETang; TopExp_Explorer Exp(S.Oriented(TopAbs_FORWARD),
        // TopAbs_EDGE); the Message_ProgressScope pair is the default-empty
        // progress form (arch. diff. #63).
        let mut a_le_tang: Vec<Shape> = Vec::new();
        for e in explorer(
            &oriented(s, Orientation::Forward),
            ShapeType::Edge,
            ShapeType::Shape,
        ) {
            // OCCT L310-313: if (!aPS.More()) { return; } — never taken with
            // the default-empty range.
            if !shape_data_map::is_bound(&self.my_map_edge_type, &e) {
                let li: Vec<BRepOffsetInterval> = Vec::new();
                shape_data_map::bind(&mut self.my_map_edge_type, &e, li);

                // OCCT L320: const NCollection_List<TopoDS_Shape>& L = Ancestors(E);
                let l = self.ancestors(&e);
                if l.is_empty() {
                    continue;
                }

                if l.len() == 2 {
                    let f1 = l[0].clone();
                    let f2 = l[1].clone();
                    // OCCT L330: EdgeAnalyse(E, F1, F2, SinTol, myMapEdgeType(E));
                    edge_analyse(
                        &self.my_brep,
                        &e,
                        &f1,
                        &f2,
                        sin_tol,
                        shape_data_map::change_find(&mut self.my_map_edge_type, &e),
                    );

                    // OCCT L332-334: for tangent faces add an artificial
                    // perpendicular face to close the gap between them (if
                    // they have different offset values)
                    if shape_data_map::value(&self.my_map_edge_type, &e)
                        .last()
                        .map(|i| i.type_of())
                        == Some(ChFiDS_TypeOfConcavity::Tangential)
                    {
                        a_le_tang.push(e.clone());
                    }
                } else if l.len() == 1 {
                    let f = l[0].clone();
                    // OCCT L343: BRep_Tool::Range(E, F, U1, U2) — raises when
                    // the pcurve is absent (Standard_NoSuchObject).
                    let (u1, u2) = brep_tool_range_on_face(&e, &f).expect(
                        "BRep_Tool::Range: no pcurve of the edge on the face (Standard_NoSuchObject)",
                    );
                    let mut inter = BRepOffsetInterval::with_params(
                        u1,
                        u2,
                        ChFiDS_TypeOfConcavity::Other,
                    );

                    if !brep_tools_is_really_closed(&e, &f) {
                        inter.set_type(ChFiDS_TypeOfConcavity::FreeBound);
                    }
                    shape_data_map::change_find(&mut self.my_map_edge_type, &e)
                        .push(inter);
                }
                // OCCT L352-357: else — "edge shared by more than two faces"
                // (the OCCT_DEBUG branch only).
            }
        }

        // OCCT L361: TreatTangentFaces(aLETang, aPSOuter.Next());
        self.treat_tangent_faces(&a_le_tang);
        // OCCT L362-365: if (!aPSOuter.More()) { return; } — never taken.
        self.my_done = true;
    }

    /// OCCT BRepOffset_Analyse::TreatTangentFaces(theLE, theRange)
    /// (cxx L371-806).
    fn treat_tangent_faces(&mut self, the_le: &[Shape]) {
        // OCCT L372: theRange — the default-empty progress form (arch. diff. #63).
        if the_le.is_empty() || self.my_face_offset_map.is_empty() {
            // OCCT L376-378: nothing to do — either there are no tangent
            // faces in the shape or the face offset map has not been provided
            return;
        }

        // OCCT L382-383: select the edges which connect faces with different
        // offset values
        let mut a_ce_tangent = empty_compound();
        // OCCT L385: bind to each tangent edge a max offset value of its faces
        let mut an_edge_offset_map: ShapeDataMap<f64> = ShapeDataMap::new();
        // OCCT L388: bind the vertices of the tangent edges with the connected
        // edges of the face with the smaller offset value
        let mut a_dvemin: ShapeDataMap<Shape> = ShapeDataMap::new();
        // OCCT L389-392: the progress scope pair — default-empty form.
        for a_e in the_le {
            // OCCT L395-398: if (!aPS1.More()) { return; } — never taken.
            let a_la = self.ancestors(a_e);

            let a_f1 = a_la.first().expect("Ancestors First: empty list");
            let a_f2 = a_la.last().expect("Ancestors Last: empty list");

            let p_offset_val1 = shape_data_map::seek(&self.my_face_offset_map, a_f1);
            let p_offset_val2 = shape_data_map::seek(&self.my_face_offset_map, a_f2);
            let an_offset_val1 = p_offset_val1.map(|v| v.abs()).unwrap_or(self.my_offset);
            let an_offset_val2 = p_offset_val2.map(|v| v.abs()).unwrap_or(self.my_offset);
            if an_offset_val1 != an_offset_val2 {
                // OCCT L410: BRep_Builder().Add(aCETangent, aE);
                bat::builder_add_compound_shape(&mut a_ce_tangent, a_e);
                // OCCT L411: anEdgeOffsetMap.Bind(aE, std::max(...));
                shape_data_map::bind(
                    &mut an_edge_offset_map,
                    a_e,
                    an_offset_val1.max(an_offset_val2),
                );

                // OCCT L413: const TopoDS_Shape& aFMin = anOffsetVal1 < anOffsetVal2 ? aF1 : aF2;
                let a_f_min = if an_offset_val1 < an_offset_val2 { a_f1 } else { a_f2 };
                // OCCT L414: for (TopoDS_Iterator itV(aE); itV.More(); itV.Next())
                for a_v in bat::sub_shapes(a_e) {
                    // OCCT L417: if (Ancestors(aV).Extent() == 3)
                    if self.ancestors(&a_v).len() == 3 {
                        // OCCT L419: for (TopExp_Explorer expE(aFMin, TopAbs_EDGE); ...)
                        for a_e_min in explorer(a_f_min, ShapeType::Edge, ShapeType::Shape) {
                            // OCCT L422: if (aEMin.IsSame(aE)) { continue; }
                            if a_e_min.is_same(a_e) {
                                continue;
                            }
                            // OCCT L426: for (TopoDS_Iterator itV1(aEMin); ...)
                            for a_vx in bat::sub_shapes(&a_e_min) {
                                // OCCT L429: if (aV.IsSame(aVx)) { aDMVEMin.Bind(aV, aEMin); }
                                if a_v.is_same(&a_vx) {
                                    shape_data_map::bind(&mut a_dvemin, &a_v, a_e_min.clone());
                                }
                            }
                        }
                    }
                }
            }
        }

        // OCCT L440-443: if (anEdgeOffsetMap.IsEmpty()) { return; }
        if an_edge_offset_map.is_empty() {
            return;
        }

        // OCCT L445-447: create the map of Face ancestors for the vertices on
        // the tangent edges
        let mut a_dmvf_anc: ShapeDataMap<Vec<Shape>> = ShapeDataMap::new();

        // OCCT L449: the progress scope — default-empty form.
        for a_e in the_le {
            // OCCT L452-455: if (!aPS2.More()) { return; } — never taken.
            // OCCT L457-460: if (!anEdgeOffsetMap.IsBound(aE)) { continue; }
            if !shape_data_map::is_bound(&an_edge_offset_map, a_e) {
                continue;
            }

            // OCCT L462-469: NCollection_Map aMFence; add the edge ancestors.
            let mut a_mfence: OcctShapeSet = HashMap::new();
            {
                let a_lea = self.ancestors(a_e);
                for a_s in &a_lea {
                    set_add(&mut a_mfence, a_s);
                }
            }

            // OCCT L471: for (TopoDS_Iterator itV(aE); itV.More(); itV.Next())
            for a_v in bat::sub_shapes(a_e) {
                // OCCT L474: pLFA = aDMVFAnc.Bound(aV, empty)
                let p_lfa = shape_data_map::bound(&mut a_dmvf_anc, &a_v);
                // OCCT L475: const ... aLVA = Ancestors(aV);
                let a_lva = self.ancestors(&a_v);
                for a_ea in &a_lva {
                    // OCCT L479-483: pIntervals = myMapEdgeType.Seek(aEA);
                    // skip when missing/empty.
                    let p_intervals = shape_data_map::seek(&self.my_map_edge_type, a_ea);
                    let Some(p_intervals) = p_intervals else {
                        continue;
                    };
                    if p_intervals.is_empty() {
                        continue;
                    }
                    // OCCT L484-487: skip the tangential first interval.
                    if p_intervals.first().map(|i| i.type_of())
                        == Some(ChFiDS_TypeOfConcavity::Tangential)
                    {
                        continue;
                    }

                    // OCCT L489-497: the face ancestors of the connected edge.
                    let a_lea = self.ancestors(a_ea);
                    for a_fa in &a_lea {
                        if set_add(&mut a_mfence, a_fa) {
                            p_lfa.push(a_fa.clone());
                        }
                    }
                }
            }
        }

        // OCCT L502: occ::handle<IntTools_Context> aCtx = new IntTools_Context();
        let _a_ctx = IntToolsContext::new();
        // OCCT L504: tangency criteria
        let a_sin_tol = self.my_angle.sin().abs();

        // OCCT L506-509: make blocks of connected edges
        // OCCT L511: BOPTools_AlgoTools::MakeConnexityBlocks(aCETangent,
        // TopAbs_VERTEX, TopAbs_EDGE, aLCB, aMVEMap) — re-hosted via the
        // wire_splitter core BFS (arch. diff. #64); theConnectionMap
        // (aMVEMap) is not exported by the re-host and is not consumed
        // downstream in the Analyse.
        let mut _a_mve_map: ShapeIndexedDataMap<Vec<Shape>> = ShapeIndexedDataMap::new();
        let a_ce_tangent_edges = explorer(&a_ce_tangent, ShapeType::Edge, ShapeType::Shape);
        let a_locations = [glam::DAffine3::IDENTITY];
        let a_lcb =
            crate::bop::algo::wire_splitter::make_connexity_blocks(&a_ce_tangent_edges, &a_locations);

        // OCCT L513-516: analyze each block to find the co-planar edges (the
        // aPS3 progress scope is the default-empty form).
        for a_cb in &a_lcb {
            // OCCT L524: const NCollection_List<TopoDS_Shape>& aCB = itLCB.Value();
            let a_cb = &a_cb.shapes;

            // OCCT L526: NCollection_Map aMFence;
            let mut a_mfence: OcctShapeSet = HashMap::new();
            // OCCT L527: for (itCB1 over aCB) — the list iterator form is the
            // index pair (itCB1, itCB2 = itCB1 advanced).
            for i_cb1 in 0..a_cb.len() {
                let a_e1 = &a_cb[i_cb1];
                // OCCT L530-533: if (!aMFence.Add(aE1)) { continue; }
                if !set_add(&mut a_mfence, a_e1) {
                    continue;
                }

                // OCCT L535-537: TopoDS_Compound aBlock; MakeCompound;
                // Add(aBlock, aE1.Oriented(TopAbs_FORWARD));
                let mut a_block = empty_compound();
                bat::builder_add_compound_shape(&mut a_block, &oriented(a_e1, Orientation::Forward));

                // OCCT L539: double anOffset = anEdgeOffsetMap.Find(aE1);
                let mut an_offset = shape_data_map::find(&an_edge_offset_map, a_e1);
                // OCCT L540: const ... aLF1 = Ancestors(aE1);
                let a_lf1 = self.ancestors(a_e1);

                // OCCT L542-543: GetNormalToFaceOnEdge(aE1, Face(aLF1.First()), aDN1)
                let mut a_dn1 = DVec3::ZERO;
                get_normal_to_face_on_edge(
                    a_e1,
                    a_lf1.first().expect("Ancestors First: empty list"),
                    &mut a_dn1,
                );

                // OCCT L545-546: NCollection_List::Iterator itCB2 = itCB1;
                // for (itCB2.Next(); itCB2.More(); itCB2.Next())
                for i_cb2 in (i_cb1 + 1)..a_cb.len() {
                    let a_e2 = &a_cb[i_cb2];
                    // OCCT L549-552: if (aMFence.Contains(aE2)) { continue; }
                    if set_contains(&a_mfence, a_e2) {
                        continue;
                    }

                    // OCCT L554: const ... aLF2 = Ancestors(aE2);
                    let a_lf2 = self.ancestors(a_e2);

                    // OCCT L556-557: GetNormalToFaceOnEdge(aE2, Face(aLF2.First()), aDN2)
                    let mut a_dn2 = DVec3::ZERO;
                    get_normal_to_face_on_edge(
                        a_e2,
                        a_lf2.first().expect("Ancestors First: empty list"),
                        &mut a_dn2,
                    );

                    // OCCT L559: if (aDN1.XYZ().Crossed(aDN2.XYZ()).Modulus() < aSinTol)
                    if a_dn1.cross(a_dn2).length() < a_sin_tol {
                        // OCCT L561-563
                        bat::builder_add_compound_shape(&mut a_block, &oriented(a_e2, Orientation::Forward));
                        set_add(&mut a_mfence, a_e2);
                        an_offset = an_offset.max(shape_data_map::find(&an_edge_offset_map, a_e2));
                    }
                }

                // OCCT L567-568: make the prism
                let a_mp = BRepPrimAPIMakePrism::new(&a_block, a_dn1 * an_offset);
                // OCCT L569-572: if (!aMP.IsDone()) { continue; }
                if !a_mp.is_done() {
                    continue;
                }

                // OCCT L574-579: the prism ancestors map.
                let mut a_prism_ancestors: ShapeIndexedDataMap<Vec<Shape>> =
                    ShapeIndexedDataMap::new();
                let a_mp_shape = a_mp.shape();
                top_exp_map_shapes_and_ancestors(
                    &a_mp_shape,
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut a_prism_ancestors,
                );
                top_exp_map_shapes_and_ancestors(
                    &a_mp_shape,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                    &mut a_prism_ancestors,
                );

                // OCCT L581: for (TopoDS_Iterator itE(aBlock); itE.More(); itE.Next())
                for a_e in bat::sub_shapes(&a_block) {
                    // OCCT L584-585: aLG = aMP.Generated(aE); aFNew = Face(aLG.First());
                    let a_lg = a_mp.generated(&a_e);
                    let mut a_f_new = a_lg
                        .first()
                        .expect("Generated First: empty list")
                        .clone();

                    // OCCT L587: NCollection_List<TopoDS_Shape>& aLA = myAncestors.ChangeFromKey(aE);
                    let a_la = shape_indexed_data_map::change_find(&mut self.my_ancestors, &a_e);

                    // OCCT L589-590: TopoDS_Shape aF1 = aLA.First(); aF2 = aLA.Last();
                    let a_f1 = a_la.first().expect("Ancestors First: empty list").clone();
                    let a_f2 = a_la.last().expect("Ancestors Last: empty list").clone();

                    // OCCT L592-595
                    let p_offset_val1 = shape_data_map::seek(&self.my_face_offset_map, &a_f1);
                    let p_offset_val2 = shape_data_map::seek(&self.my_face_offset_map, &a_f2);
                    let an_offset_val1 = p_offset_val1.map(|v| v.abs()).unwrap_or(self.my_offset);
                    let an_offset_val2 = p_offset_val2.map(|v| v.abs()).unwrap_or(self.my_offset);

                    // OCCT L597-598
                    let a_f_to_remove = if an_offset_val1 > an_offset_val2 {
                        a_f1.clone()
                    } else {
                        a_f2.clone()
                    };
                    let a_f_opposite = if an_offset_val1 > an_offset_val2 {
                        a_f2
                    } else {
                        a_f1
                    };

                    // OCCT L600-635: orient the face so its normal is directed
                    // to the smaller offset face
                    {
                        // OCCT L602-604: get the normal of the new face
                        let mut a_dn = DVec3::ZERO;
                        get_normal_to_face_on_edge(&a_e, &a_f_new, &mut a_dn);

                        // OCCT L606-615: get the bi-normal for the aFOpposite
                        let mut a_e_in_f: Option<Shape> = None;
                        for a_ee in explorer(&a_f_opposite, ShapeType::Edge, ShapeType::Shape) {
                            if a_e.is_same(&a_ee) {
                                a_e_in_f = Some(a_ee);
                                break;
                            }
                        }
                        // OCCT L620 dereferences aEInF without a null check
                        // (BRep_Tool::Curve raises on the null edge — the
                        // expect keeps the same failure class).
                        let a_e_in_f =
                            a_e_in_f.expect("aEInF: the edge is not in aFOpposite (OCCT leaves it unchecked)");

                        let (a_c3d, f, l) = brep_tool_curve(&a_e_in_f)
                            .expect("BRep_Tool::Curve: no 3D curve of the edge");
                        // OCCT L621: aPOnE = aC3D->Value((f + l) / 2.);
                        let a_p_on_e = a_c3d.point_at((f + l) / 2.0);
                        // OCCT L617-627
                        let mut a_p2d = DVec2::ZERO;
                        let mut a_p_in_f = DVec3::ZERO;
                        point_near_edge(
                            &a_e_in_f,
                            &a_f_opposite,
                            (f + l) / 2.0,
                            1.0e-5,
                            &mut a_p2d,
                            &mut a_p_in_f,
                        );

                        // OCCT L629: gp_Vec aBN(aPOnE, aPInF);
                        let a_bn = a_p_in_f - a_p_on_e;

                        // OCCT L631-634
                        if a_bn.dot(a_dn) < 0.0 {
                            a_f_new = shape_reverse(a_f_new);
                        }
                    }

                    // OCCT L637-646: remove the face with the bigger offset
                    // value from the edge ancestors
                    {
                        let a_la = shape_indexed_data_map::change_find(&mut self.my_ancestors, &a_e);
                        // OCCT L638-645: for (itA over aLA) — Remove(itA) + break.
                        if let Some(pos) = a_la
                            .iter()
                            .position(|x| x.is_same(&a_f_to_remove))
                        {
                            a_la.remove(pos);
                        }
                        // OCCT L646: aLA.Append(aFNew);
                        a_la.push(a_f_new.clone());
                    }

                    // OCCT L648: myMapEdgeType(aE).Clear();
                    shape_data_map::change_find(&mut self.my_map_edge_type, &a_e).clear();
                    // OCCT L650: analyze the edge again
                    edge_analyse(
                        &self.my_brep,
                        &a_e,
                        &a_f_opposite,
                        &a_f_new,
                        a_sin_tol,
                        shape_data_map::change_find(&mut self.my_map_edge_type, &a_e),
                    );

                    // OCCT L652-743: analyze the vertices
                    let mut a_f_new_edge_map: OcctShapeSet = HashMap::new();
                    set_add(&mut a_f_new_edge_map, &a_e);
                    // OCCT L655: for (TopoDS_Iterator itV(aE); itV.More(); itV.Next())
                    for a_v in bat::sub_shapes(&a_e) {
                        // OCCT L659: aEG = Edge(aMP.Generated(aV).First());
                        let mut a_eg = a_mp
                            .generated(&a_v)
                            .first()
                            .expect("Generated First: empty list")
                            .clone();
                        // OCCT L660: myGenerated.Bind(aV, aEG);
                        shape_data_map::bind(&mut self.my_generated, &a_v, a_eg.clone());
                        {
                            // OCCT L662-670: take the aFNew-resident twin of
                            // aEG (the correct orientation).
                            for a_ee in explorer(&a_f_new, ShapeType::Edge, ShapeType::Shape) {
                                if a_ee.is_same(&a_eg) {
                                    a_eg = a_ee.clone();
                                    break;
                                }
                            }
                        }

                        // OCCT L672: if (aDMVEMin.IsBound(aV))
                        if shape_data_map::is_bound(&a_dvemin, &a_v) {
                            // OCCT L674: pSA = aDMVFAnc.Seek(aV);
                            if let Some(p_sa) = shape_data_map::seek(&a_dmvf_anc, &a_v) {
                                // OCCT L675: if (pSA && pSA->Extent() == 1)
                                if p_sa.len() == 1 {
                                    // OCCT L677-686: adjust the orientation of
                                    // the generated edge to its new ancestor.
                                    let mut a_e_min = shape_data_map::find(&a_dvemin, &a_v);
                                    for a_ee in
                                        explorer(&p_sa[0], ShapeType::Edge, ShapeType::Shape)
                                    {
                                        if a_ee.is_same(&a_e_min) {
                                            a_e_min = a_ee.clone();
                                            break;
                                        }
                                    }

                                    let mut an_ori_in_e_min = Orientation::Forward;
                                    for a_s in bat::sub_shapes(&a_e_min) {
                                        if a_s.is_same(&a_v) {
                                            an_ori_in_e_min = a_s.orientation;
                                            break;
                                        }
                                    }

                                    let mut an_ori_in_eg = Orientation::Forward;
                                    for a_s in bat::sub_shapes(&a_eg) {
                                        if a_s.is_same(&a_v) {
                                            an_ori_in_eg = a_s.orientation;
                                            break;
                                        }
                                    }

                                    // OCCT L708-711
                                    if an_ori_in_eg == an_ori_in_e_min {
                                        a_eg = shape_reverse(a_eg);
                                    }
                                }
                            }
                        }

                        // OCCT L715-719
                        {
                            let a_lva = shape_indexed_data_map::change_find(
                                &mut self.my_ancestors,
                                &a_v,
                            );
                            if !a_lva.iter().any(|x| x.is_same(&a_eg)) {
                                a_lva.push(a_eg.clone());
                            }
                        }
                        // OCCT L720: aFNewEdgeMap.Add(aEG);
                        set_add(&mut a_f_new_edge_map, &a_eg);

                        // OCCT L722-723: aLEGA = myAncestors(myAncestors.Add(aEG,
                        // aPrismAncestors.FindFromKey(aEG)));
                        let idx = shape_indexed_data_map::add(
                            &mut self.my_ancestors,
                            &a_eg,
                            shape_indexed_data_map::find(&a_prism_ancestors, &a_eg).clone(),
                        );
                        let mut a_lega_len = 0usize;
                        let mut a_lega_first: Option<Shape> = None;
                        let mut a_lega_last: Option<Shape> = None;
                        {
                            let a_lega =
                                shape_indexed_data_map::value_1_mut(&mut self.my_ancestors, idx + 1);
                            // OCCT L724-732: add the ancestors from the shape
                            if let Some(p_sa) = shape_data_map::seek(&a_dmvf_anc, &a_v) {
                                if !p_sa.is_empty() {
                                    let a_lsa = p_sa.clone();
                                    a_lega.extend(a_lsa);
                                }
                            }
                            a_lega_len = a_lega.len();
                            a_lega_first = a_lega.first().cloned();
                            a_lega_last = a_lega.last().cloned();
                        }

                        // OCCT L734: myMapEdgeType.Bind(aEG, empty);
                        shape_data_map::bind(&mut self.my_map_edge_type, &a_eg, Vec::new());
                        // OCCT L735-742: if (aLEGA.Extent() == 2) — analyze aEG.
                        if a_lega_len == 2 {
                            let a_lf_first = a_lega_first.expect("aLEGA First: empty list");
                            let a_lf_last = a_lega_last.expect("aLEGA Last: empty list");
                            edge_analyse(
                                &self.my_brep,
                                &a_eg,
                                &a_lf_first,
                                &a_lf_last,
                                a_sin_tol,
                                shape_data_map::change_find(&mut self.my_map_edge_type, &a_eg),
                            );
                        }
                    }

                    // OCCT L745-754: find the edge opposite to the tangential
                    // one and add the ancestors for it.
                    let mut a_e_opposite: Option<Shape> = None;
                    for a_ee in explorer(&a_f_new, ShapeType::Edge, ShapeType::Shape) {
                        if !set_contains(&a_f_new_edge_map, &a_ee) {
                            a_e_opposite = Some(a_ee.clone());
                            break;
                        }
                    }
                    // OCCT dereferences aEOpposite at L763 without a null check.
                    let mut a_e_opposite = a_e_opposite
                        .expect("aEOpposite: not found in aFNew (OCCT leaves it unchecked)");

                    {
                        // OCCT L756-770: find it in aFToRemove (the OCCT
                        // comment says aFOpposite, the code iterates aFToRemove).
                        for a_e_in_f_to_rem in
                            explorer(&a_f_to_remove, ShapeType::Edge, ShapeType::Shape)
                        {
                            if a_e.is_same(&a_e_in_f_to_rem) {
                                // OCCT L763: IsSplitToReverse(aEOpposite,
                                // aEInFToRem, aCtx) — the edge dispatch (arch.
                                // diff. #65; the rcad form is context-free).
                                if is_split_to_reverse_edge(&a_e_opposite, &a_e_in_f_to_rem).0 {
                                    // OCCT L765: aEOpposite.Reverse();
                                    a_e_opposite = shape_reverse(a_e_opposite);
                                }
                                break;
                            }
                        }
                    }

                    // OCCT L772-775
                    let mut a_lf_opposite: Vec<Shape> = Vec::new();
                    a_lf_opposite.push(a_f_new.clone());
                    a_lf_opposite.push(a_f_to_remove.clone());
                    shape_indexed_data_map::add(&mut self.my_ancestors, &a_e_opposite, a_lf_opposite);
                    // OCCT L776: myMapEdgeType.Bind(aEOpposite, empty);
                    shape_data_map::bind(&mut self.my_map_edge_type, &a_e_opposite, Vec::new());
                    // OCCT L777-781: EdgeAnalyse(aEOpposite, aFNew, Face(aFToRemove), aSinTol, ...)
                    edge_analyse(
                        &self.my_brep,
                        &a_e_opposite,
                        &a_f_new,
                        &a_f_to_remove,
                        a_sin_tol,
                        shape_data_map::change_find(&mut self.my_map_edge_type, &a_e_opposite),
                    );

                    // OCCT L783-791: pEEMap = myReplacement.ChangeSeek(aFToRemove)
                    // / Bound / Bind(aE, aEOpposite).
                    if !shape_data_map::is_bound(&self.my_replacement, &a_f_to_remove) {
                        shape_data_map::bind(
                            &mut self.my_replacement,
                            &a_f_to_remove,
                            ShapeDataMap::new(),
                        );
                    }
                    let p_ee_map =
                        shape_data_map::change_find(&mut self.my_replacement, &a_f_to_remove);
                    shape_data_map::bind(p_ee_map, &a_e, a_e_opposite.clone());

                    // OCCT L793-799: add the ancestors for the vertices.
                    for a_v in bat::sub_shapes(&a_e_opposite) {
                        let a_lva =
                            shape_indexed_data_map::find(&a_prism_ancestors, &a_v).clone();
                        shape_indexed_data_map::add(&mut self.my_ancestors, &a_v, a_lva);
                    }

                    // OCCT L801-802
                    self.my_new_faces.push(a_f_new.clone());
                    shape_data_map::bind(&mut self.my_generated, &a_e, a_f_new);
                }
            }
        }
    }

    /// OCCT BRepOffset_Analyse::EdgeReplacement(theF, theE) (cxx L810-827).
    pub fn edge_replacement(&self, the_f: &Shape, the_e: &Shape) -> Shape {
        // OCCT L813: pEE = myReplacement.Seek(theF);
        let Some(p_ee) = shape_data_map::seek(&self.my_replacement, the_f) else {
            return the_e.clone();
        };

        // OCCT L820: pE = pEE->Seek(theE);
        let Some(p_e) = shape_data_map::seek(p_ee, the_e) else {
            return the_e.clone();
        };

        p_e.clone()
    }

    /// OCCT BRepOffset_Analyse::Generated(theS) (cxx L831-836).
    pub fn generated(&self, the_s: &Shape) -> Shape {
        shape_data_map::seek(&self.my_generated, the_s)
            .cloned()
            .unwrap_or_else(Shape::null)
    }

    /// OCCT BRepOffset_Analyse::Descendants(theS, theUpdate = false)
    /// (cxx L840-870) — the default-argument form.
    pub fn descendants(&self, the_s: &Shape) -> Option<Vec<Shape>> {
        self.descendants_with_update(the_s, false)
    }

    /// OCCT BRepOffset_Analyse::Descendants(theS, theUpdate) (cxx L840-870) —
    /// the full form (theUpdate materializes the OCCT default parameter).
    pub fn descendants_with_update(&self, the_s: &Shape, the_update: bool) -> Option<Vec<Shape>> {
        if self.my_descendants.borrow().is_empty() || the_update {
            self.my_descendants.borrow_mut().clear();
            // OCCT L846: const int aNbA = myAncestors.Extent();
            let a_nb_a = shape_indexed_data_map::extent(&self.my_ancestors);
            for i in 1..=a_nb_a {
                // OCCT L849-850: aSS = myAncestors.FindKey(i); aLA = myAncestors(i);
                let a_ss = shape_indexed_data_map::find_key_1(&self.my_ancestors, i).clone();
                let a_la = shape_indexed_data_map::value_1(&self.my_ancestors, i).clone();

                let mut a_descendants = self.my_descendants.borrow_mut();
                for a_sa in &a_la {
                    // OCCT L856-860: pLD = myDescendants.ChangeSeek(aSA) / Bound.
                    let p_ld = shape_data_map::bound(&mut a_descendants, a_sa);
                    // OCCT L861-864
                    if !p_ld.iter().any(|x| x.is_same(&a_ss)) {
                        p_ld.push(a_ss.clone());
                    }
                }
            }
        }

        // OCCT L869: return myDescendants.Seek(theS);
        shape_data_map::seek(&self.my_descendants.borrow(), the_s).cloned()
    }

    /// OCCT BRepOffset_Analyse::Clear() (cxx L874-885).
    pub fn clear(&mut self) {
        self.my_done = false;
        // OCCT L877: myShape.Nullify();
        self.my_shape = Shape::null();
        self.my_map_edge_type.clear();
        self.my_ancestors.clear();
        self.my_face_offset_map.clear();
        self.my_replacement.clear();
        self.my_descendants.borrow_mut().clear();
        self.my_new_faces.clear();
        self.my_generated.clear();
    }

    /// OCCT BRepOffset_Analyse::Type(E) (cxx L891-894).
    pub fn type_(&self, e: &Shape) -> Vec<BRepOffsetInterval> {
        shape_data_map::value(&self.my_map_edge_type, e).clone()
    }

    /// OCCT BRepOffset_Analyse::Edges(V, T, LE) (cxx L898-930).
    pub fn edges_on_vertex(
        &self,
        v: &Shape,
        t: ChFiDS_TypeOfConcavity,
        le: &mut Vec<Shape>,
    ) {
        le.clear();
        // OCCT L903: const ... L = Ancestors(V);
        let l = self.ancestors(v);
        for a_s in &l {
            let e = a_s;
            // OCCT L909: pIntervals = myMapEdgeType.Seek(E);
            if let Some(p_intervals) = shape_data_map::seek(&self.my_map_edge_type, e) {
                if !p_intervals.is_empty() {
                    // OCCT L912-913: BRepOffset_Tool::EdgeVertices(E, V1, V2);
                    let mut v1 = Shape::null();
                    let mut v2 = Shape::null();
                    edge_vertices(e, &mut v1, &mut v2);
                    // OCCT L914-920
                    if v1.is_same(v) {
                        if p_intervals
                            .last()
                            .map(|i| i.type_of())
                            == Some(t)
                        {
                            le.push(e.clone());
                        }
                    }
                    // OCCT L921-927
                    if v2.is_same(v) {
                        if p_intervals
                            .first()
                            .map(|i| i.type_of())
                            == Some(t)
                        {
                            le.push(e.clone());
                        }
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_Analyse::Edges(F, T, LE) (cxx L934-955).
    pub fn edges_on_face(&self, f: &Shape, t: ChFiDS_TypeOfConcavity, le: &mut Vec<Shape>) {
        le.clear();
        // OCCT L939: TopExp_Explorer exp(F, TopAbs_EDGE);
        for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
            // OCCT L945: const ... Lint = Type(E);
            let lint = self.type_(&e);
            for it in &lint {
                if it.type_of() == t {
                    le.push(e.clone());
                }
            }
        }
    }

    /// OCCT BRepOffset_Analyse::TangentEdges(Edge, Vertex, Edges)
    /// (cxx L959-1001).
    pub fn tangent_edges(&self, edge: &Shape, vertex: &Shape, edges: &mut Vec<Shape>) {
        // OCCT L968: URef = BRep_Tool::Parameter(Vertex, Edge);
        let u_ref = brep_tool_parameter(vertex, edge);
        // OCCT L969-970: C3dRef = BRepAdaptor_Curve(Edge); VRef = C3dRef.DN(URef, 1);
        let (c3d_ref, _, _) = brep_tool_curve(edge)
            .expect("BRepAdaptor_Curve: no 3D curve of the edge");
        let mut v_ref = c3d_ref.derivative_at(u_ref);
        correct_orientation_of_tangent(&mut v_ref, vertex, edge);
        // OCCT L972-975
        if v_ref.length_squared() < GP_RESOLUTION {
            return;
        }

        // OCCT L977: Edges.Clear();
        edges.clear();

        // OCCT L979: const ... Anc = Ancestors(Vertex);
        let anc = self.ancestors(vertex);
        for cur_e in &anc {
            // OCCT L984-987
            if cur_e.is_same(edge) {
                continue;
            }
            // OCCT L988-990: U = BRep_Tool::Parameter(Vertex, CurE);
            // C3d = BRepAdaptor_Curve(CurE); V = C3d.DN(U, 1);
            let u = brep_tool_parameter(vertex, cur_e);
            let (c3d, _, _) = brep_tool_curve(cur_e)
                .expect("BRepAdaptor_Curve: no 3D curve of the edge");
            let mut v = c3d.derivative_at(u);
            correct_orientation_of_tangent(&mut v, vertex, cur_e);
            // OCCT L992-995
            if v.length_squared() < GP_RESOLUTION {
                continue;
            }
            // OCCT L996-999
            if gp_vec_is_opposite(&v, &v_ref, self.my_angle) {
                edges.push(cur_e.clone());
            }
        }
    }

    /// OCCT BRepOffset_Analyse::Explode(List, T) (cxx L1005-1027).
    pub fn explode(&self, list: &mut Vec<Shape>, t: ChFiDS_TypeOfConcavity) {
        list.clear();
        // OCCT L1010: NCollection_Map Map;
        let mut map: OcctShapeSet = HashMap::new();

        // OCCT L1012-1013: TopExp_Explorer Fexp; Fexp.Init(myShape, TopAbs_FACE);
        for a_f in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if set_add(&mut map, &a_f) {
                // OCCT L1017-1020
                let face = a_f.clone();
                let mut co = empty_compound();
                bat::builder_add_compound_shape(&mut co, &face);
                // OCCT L1021-1023: add to Co all the faces from the cloud of
                // faces G1 created from <Face>
                self.add_faces(&face, &mut co, &mut map, t);
                list.push(co);
            }
        }
    }

    /// OCCT BRepOffset_Analyse::Explode(List, T1, T2) (cxx L1031-1054).
    pub fn explode_2_types(
        &self,
        list: &mut Vec<Shape>,
        t1: ChFiDS_TypeOfConcavity,
        t2: ChFiDS_TypeOfConcavity,
    ) {
        list.clear();
        // OCCT L1037: NCollection_Map Map;
        let mut map: OcctShapeSet = HashMap::new();

        // OCCT L1039-1040
        for a_f in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            if set_add(&mut map, &a_f) {
                // OCCT L1044-1047
                let face = a_f.clone();
                let mut co = empty_compound();
                bat::builder_add_compound_shape(&mut co, &face);
                // OCCT L1048-1050
                self.add_faces_2_types(&face, &mut co, &mut map, t1, t2);
                list.push(co);
            }
        }
    }

    /// OCCT BRepOffset_Analyse::AddFaces(Face, Co, Map, T) (cxx L1058-1092).
    pub fn add_faces(
        &self,
        face: &Shape,
        co: &mut Shape,
        map: &mut OcctShapeSet,
        t: ChFiDS_TypeOfConcavity,
    ) {
        // OCCT L1064: pLE = Descendants(Face);
        let Some(p_le) = self.descendants(face) else {
            return;
        };
        for a_e in &p_le {
            let e = a_e;
            // OCCT L1072: LI = Type(E);
            let li = self.type_(e);
            // OCCT L1073: if (!LI.IsEmpty() && LI.First().Type() == T)
            if !li.is_empty() && li.first().map(|i| i.type_of()) == Some(t) {
                // OCCT L1075: so <NewFace> is attached to G1 by <Face>
                let l = self.ancestors(e);
                if l.len() == 2 {
                    // OCCT L1079-1083
                    let mut f1 = l[0].clone();
                    if f1.is_same(face) {
                        f1 = l[1].clone();
                    }
                    // OCCT L1084-1088
                    if set_add(map, &f1) {
                        bat::builder_add_compound_shape(co, &f1);
                        self.add_faces(&f1, co, map, t);
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_Analyse::AddFaces(Face, Co, Map, T1, T2)
    /// (cxx L1096-1131).
    pub fn add_faces_2_types(
        &self,
        face: &Shape,
        co: &mut Shape,
        map: &mut OcctShapeSet,
        t1: ChFiDS_TypeOfConcavity,
        t2: ChFiDS_TypeOfConcavity,
    ) {
        // OCCT L1103: pLE = Descendants(Face);
        let Some(p_le) = self.descendants(face) else {
            return;
        };
        for a_e in &p_le {
            let e = a_e;
            // OCCT L1111: LI = Type(E);
            let li = self.type_(e);
            // OCCT L1112
            if !li.is_empty()
                && (li.first().map(|i| i.type_of()) == Some(t1)
                    || li.first().map(|i| i.type_of()) == Some(t2))
            {
                // OCCT L1114: so <NewFace> is attached to G1 by <Face>
                let l = self.ancestors(e);
                if l.len() == 2 {
                    // OCCT L1118-1122
                    let mut f1 = l[0].clone();
                    if f1.is_same(face) {
                        f1 = l[1].clone();
                    }
                    // OCCT L1123-1127
                    if set_add(map, &f1) {
                        bat::builder_add_compound_shape(co, &f1);
                        self.add_faces_2_types(&f1, co, map, t1, t2);
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_Analyse::IsDone() (hxx L61).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepOffset_Analyse::HasAncestor(S) (hxx L85).
    pub fn has_ancestor(&self, the_s: &Shape) -> bool {
        self.my_ancestors.contains_key(&bat::shape_key(the_s))
    }

    /// OCCT BRepOffset_Analyse::Ancestors(S) (hxx L88-91) — FindFromKey
    /// (asserts when unbound); the rcad form returns the list by value.
    pub fn ancestors(&self, the_s: &Shape) -> Vec<Shape> {
        shape_indexed_data_map::find(&self.my_ancestors, the_s).clone()
    }

    /// OCCT BRepOffset_Analyse::SetOffsetValue(theOffset) (hxx L119).
    pub fn set_offset_value(&mut self, the_offset: f64) {
        self.my_offset = the_offset;
    }

    /// OCCT BRepOffset_Analyse::SetFaceOffsetMap(theMap) (hxx L122-126).
    pub fn set_face_offset_map(&mut self, the_map: &ShapeDataMap<f64>) {
        self.my_face_offset_map = the_map.clone();
    }

    /// OCCT BRepOffset_Analyse::NewFaces() (hxx L130).
    pub fn new_faces(&self) -> Vec<Shape> {
        self.my_new_faces.clone()
    }

    /// OCCT BRepOffset_Analyse::HasGenerated(S) (hxx L137).
    pub fn has_generated(&self, the_s: &Shape) -> bool {
        shape_data_map::seek(&self.my_generated, the_s).is_some()
    }
}

// OCCT BRepOffset_Analyse.cxx L62/L79/L146/L274/L369/L808/L829/L838/L872/
// L887/L896/L932/L957/L1003/L1029/L1056/L1094: the section separators.
#[cfg(test)]
mod tests {
    use super::*;

    // The Analyse shell tests only — the geometry-dependent paths
    // (Perform / TreatTangentFaces) require the shape fixtures and are
    // exercised by the OCCT offset grids, not here.

    #[test]
    fn interval_accessors() {
        // OCCT BRepOffset_Interval lxx L24-66.
        let mut i = BRepOffsetInterval::new();
        assert_eq!(i.first(), 0.0);
        assert_eq!(i.last(), 0.0);
        assert_eq!(i.type_of(), ChFiDS_TypeOfConcavity::Other);
        i.set_first(1.5);
        i.set_last(2.5);
        i.set_type(ChFiDS_TypeOfConcavity::Tangential);
        assert_eq!(i.first(), 1.5);
        assert_eq!(i.last(), 2.5);
        assert_eq!(i.type_of(), ChFiDS_TypeOfConcavity::Tangential);

        let j = BRepOffsetInterval::with_params(0.25, 0.75, ChFiDS_TypeOfConcavity::FreeBound);
        assert_eq!(j.first(), 0.25);
        assert_eq!(j.last(), 0.75);
        assert_eq!(j.type_of(), ChFiDS_TypeOfConcavity::FreeBound);
    }

    #[test]
    fn analyse_default_state() {
        // OCCT cxx L64-68: the empty ctor leaves myDone = false.
        let mut analyse = BRepOffsetAnalyse::new();
        assert!(!analyse.is_done());
        analyse.clear();
        assert!(!analyse.is_done());
    }

    #[test]
    fn gp_vec_helpers() {
        // OCCT gp_Vec.hxx L130-134: IsOpposite of the exact opposite vector
        // is true for any non-negative tolerance.
        let v = DVec3::new(1.0, 0.0, 0.0);
        assert!(gp_vec_is_opposite(&v, &-v, 1e-12));
        assert!(!gp_vec_is_opposite(&v, &v, 1e-12));
        // OCCT gp_Dir.cxx L27-50: Angle of perpendicular vectors.
        assert!((gp_vec_angle(&v, &DVec3::Y) - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn geom_abs_shape_rank_order() {
        // OCCT GeomAbs_Shape.hxx enum order.
        assert!(geom_abs_shape_rank(GeomAbsShape::C0) < geom_abs_shape_rank(GeomAbsShape::G1));
        assert!(geom_abs_shape_rank(GeomAbsShape::C1) < geom_abs_shape_rank(GeomAbsShape::C3));
        assert!(geom_abs_shape_rank(GeomAbsShape::C3) < geom_abs_shape_rank(GeomAbsShape::CN));
    }
}
