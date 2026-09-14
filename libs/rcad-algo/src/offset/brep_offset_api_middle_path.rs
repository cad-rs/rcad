//! OCCT BRepOffsetAPI_MiddlePath — 1:1 translation (statics + class).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MiddlePath.cxx L35-258 (the statics, the
//!         commented-out forms excluded) + BRepOffsetAPI_MiddlePath.hxx
//!         (L30-110).
//! Build (cxx L265-1010) lives in brep_offset_api_middle_path_b.rs.
//!
//! OCCT inheritance chain (hxx L30): BRepOffsetAPI_MiddlePath ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base members
//! (myShape, the Done flag) are plain fields (the Stage 2e facade
//! precedent).
//!
//! Architecture differences:
//! 1. NCollection_Sequence -> Vec; NCollection_Map<TopoDS_Shape> ->
//!    HashMap<u64, Shape> (the TopTools_ShapeMapHasher key).
//! 2. ShapeUpgrade_UnifySameDomain (TKShHealing) — rewired to the real
//!    body crate::shhealing::shape_upgrade::unify_same_domain::
//!    ShapeUpgradeUnifySameDomain (the GAP carrier is retired, Rule 4);
//!    the real Build is the arena form (cxx L4454-4469), so the
//!    constructor's pool becomes the class my_brep (architecture
//!    difference #4).
//! 3. BRepExtrema_DistShapeShape (TKTopAlgo/BRepExtrema) and
//!    GeomAPI_Interpolate (TKGeomAlgo — the kernel carries the fn-form
//!    interpolate_points only, no tangent-Load class form) — GAP.  The
//!    BRepGProp::SurfaceProperties / LinearProperties calls feed
//!    Properties.CentreOfMass() (cxx L875-890) — GAP (the kernel gprop
//!    re-hosts carry the mass only, no GProp_GProps COM).  GeomLib::
//!    Inertia — rewired to the kernel base::geom_lib::inertia real body.
//!    BRepLib::BuildCurve3d — the arena-form real body
//!    (topalgo::brep_lib::BRepLib::build_curve3d) needs the &mut BRep
//!    pool while the Build call sites carry pool-free Arc<TShape>
//!    handles; the brep_offset_offset.rs (E, Tol) re-host serves the
//!    calls (the brep_offset_inter2d.rs:85 convention).
//! 4. GeomAbs_CurveType -> the canonical rcad_kernel::math enum (re-exported
//!    below); the BRepAdaptor_Curve::GetType
//!    reads map to the Curve3 variant discriminant; the Geom_Line/Bezier/
    //!    BSpline analyses of the statics map to the rcad curve data.
//! 5. TangentOfEdge -> the CurveEval tangent_at carrier; TopExp::
//!    FirstVertex/LastVertex(E, CumOri=true) -> the oriented
//!    top_exp_vertices read of brep_fill/generator.
//! 6. BRepLib_MakeWire -> the BRepLibMakeWire carrier of
//!    brep_algo::normal_projection (Add stores, IsDone=false keeps the OCCT
//!    not-done exit); BRepLib_MakeFace(W, OnlyPlane) -> the BRepLibMakeFace
//!    carrier of brep_offset_make_simple_offset.

use std::collections::HashMap;

use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::topo::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

use glam::DVec3;

use crate::brep_algo::tool::{
    brep_tool_curve, builder_add_edge_vertex, builder_make_edge, builder_make_vertex,
    top_exp_vertices_raw,
};
use rcad_kernel::topo::topods::TShape;

/// OCCT GeomAbs_CurveType (TKMath/GeomAbs/GeomAbs_CurveType.hxx L23-31).
/// Canonical nine-variant enum lives in rcad_kernel::math; re-exported here
/// so the historic import path crate::offset::brep_offset_api_middle_path::
/// GeomAbsCurveType (brep_offset_api_middle_path_b.rs) keeps resolving.  The
/// former local copy (with the shortened Bezier variant) was deleted (Rule 4).
pub use rcad_kernel::math::GeomAbsCurveType;

// ---------------------------------------------------------------------------
// GAP carriers (architecture difference #3).
// ---------------------------------------------------------------------------

use crate::shhealing::shape_upgrade::unify_same_domain::ShapeUpgradeUnifySameDomain;

/// OCCT BRepExtrema_SupportType (TKTopAlgo/BRepExtrema_SupportType.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepExtremaSupportType {
    IsVertex,
    IsOnEdge,
}

/// OCCT BRepExtrema_DistShapeShape (TKTopAlgo/BRepExtrema) — the
/// shape-to-shape distance of IsValidEdge (architecture difference #3;
/// GAP: the general shape distance is not translated — the method surface
/// keeps the OCCT form).
pub struct BRepExtremaDistShapeShape;

impl BRepExtremaDistShapeShape {
    /// OCCT BRepExtrema_DistShapeShape(S1, S2).
    pub fn new(_the_s1: &Shape, _the_s2: &Shape) -> Self {
        BRepExtremaDistShapeShape
    }

    /// OCCT BRepExtrema_DistShapeShape::Value() — GAP.
    pub fn value(&self) -> f64 {
        panic!("GAP: BRepExtrema_DistShapeShape::Value (TKTopAlgo/BRepExtrema not translated)")
    }

    /// OCCT BRepExtrema_DistShapeShape::NbSolution() — GAP.
    pub fn nb_solution(&self) -> i32 {
        panic!("GAP: BRepExtrema_DistShapeShape::NbSolution (TKTopAlgo/BRepExtrema not translated)")
    }

    /// OCCT BRepExtrema_DistShapeShape::SupportTypeShape2(i) — GAP.
    pub fn support_type_shape2(&self, _the_i: i32) -> BRepExtremaSupportType {
        panic!("GAP: BRepExtrema_DistShapeShape::SupportTypeShape2 (TKTopAlgo/BRepExtrema not translated)")
    }

    /// OCCT BRepExtrema_DistShapeShape::SupportOnShape2(i) — GAP.
    pub fn support_on_shape2(&self, _the_i: i32) -> Shape {
        panic!("GAP: BRepExtrema_DistShapeShape::SupportOnShape2 (TKTopAlgo/BRepExtrema not translated)")
    }
}

/// OCCT GProp_GProps (TKMath/GProp_GProps) — the properties carrier of the
/// BRepGProp GAP.
pub struct GPropGProps {
    my_centre: DVec3, // OCCT: myG
}

impl GPropGProps {
    /// OCCT GProp_GProps::CentreOfMass().
    pub fn centre_of_mass(&self) -> DVec3 {
        self.my_centre
    }
}

/// OCCT BRepGProp (TKTopAlgo/BRepGProp) — the shape properties driver
/// (architecture difference #3; GAP): the calls of cxx L875-890 fill
/// GProp_GProps and read Properties.CentreOfMass(), while the kernel
/// base::gprop re-hosts carry only the mass (base::gprop::surface::
/// surface_area f64 / base::gprop::linear::linear_properties f64) — no
/// GProp_GProps centre of mass.
pub struct BRepGProp;

impl BRepGProp {
    /// OCCT BRepGProp::SurfaceProperties(S, Props) — GAP (see above).
    pub fn surface_properties(_the_s: &Shape) -> GPropGProps {
        panic!("GAP: BRepGProp::SurfaceProperties (TKTopAlgo/BRepGProp not translated)")
    }

    /// OCCT BRepGProp::LinearProperties(S, Props) — GAP (see above).
    pub fn linear_properties(_the_s: &Shape) -> GPropGProps {
        panic!("GAP: BRepGProp::LinearProperties (TKTopAlgo/BRepGProp not translated)")
    }
}

/// OCCT GeomAPI_Interpolate (TKGeomAlgo/GeomAPI_Interpolate) — the
/// interpolation of the missed middle edges (architecture difference #3;
/// GAP: no rcad translation yet — IsDone() = false keeps the OCCT
/// not-done branch of Build).
pub struct GeomAPIInterpolate;

impl GeomAPIInterpolate {
    /// OCCT GeomAPI_Interpolate(Points, PeriodicFlag, Tolerance).
    pub fn new(_the_points: &[DVec3], _the_periodic_flag: bool, _the_tolerance: f64) -> Self {
        GeomAPIInterpolate
    }

    /// OCCT GeomAPI_Interpolate::Load(Tangents, Flags) — GAP.
    pub fn load(&mut self, _the_tangents: &[DVec3], _the_flags: &[bool]) {
        panic!("GAP: GeomAPI_Interpolate::Load (TKGeomAlgo not translated)")
    }

    /// OCCT GeomAPI_Interpolate::Perform() — GAP.
    pub fn perform(&mut self) {
        panic!("GAP: GeomAPI_Interpolate::Perform (TKGeomAlgo not translated)")
    }

    /// OCCT GeomAPI_Interpolate::IsDone() — GAP (false).
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT GeomAPI_Interpolate::Curve() — GAP.
    pub fn curve(&self) -> Curve3 {
        panic!("GAP: GeomAPI_Interpolate::Curve (TKGeomAlgo not translated)")
    }
}

/// OCCT GeomLib::Inertia(Pnts, Bary, Xdir, Ydir, Xgap, Ygap, Zgap)
/// (TKGeomAlgo/GeomLib, GeomLib.cxx L1976-2093; the call of cxx L1086) —
/// rewired to the kernel base::geom_lib::inertia real body (the exact OCCT
/// out-param signature).
pub fn geom_lib_inertia(
    the_pnts: &[DVec3],
    the_bary: &mut DVec3,
    the_xdir: &mut DVec3,
    the_ydir: &mut DVec3,
    the_xgap: &mut f64,
    the_ygap: &mut f64,
    the_zgap: &mut f64,
) {
    rcad_kernel::base::geom_lib::inertia(
        the_pnts,
        the_bary,
        the_xdir,
        the_ydir,
        the_xgap,
        the_ygap,
        the_zgap,
    )
}

/// OCCT BRepLib::BuildCurve3d(E) (TKBRep/BRepLib; BRepLib.hxx L90-95 — the
/// defaulted form Tolerance = 1.0e-5, Continuity = GeomAbs_C1, MaxDegree =
/// 14, MaxSegment = 0; the calls of cxx L697/L735/L743/L747) — the real
/// arena-form body (topalgo::brep_lib::BRepLib::build_curve3d) needs the
/// &mut BRep pool while the Build call sites carry pool-free Arc<TShape>
/// handles (architecture differences #19/#22 of brep_offset_offset.rs), so
/// the brep_offset_offset.rs (E, Tol) re-host serves the call (the
/// brep_offset_inter2d.rs:85 convention); the OCCT bool result is
/// discarded at these call sites.
pub fn brep_lib_build_curve3d(the_e: &Shape) {
    super::brep_offset_offset::brep_lib_build_curve3d(the_e, 1.0e-5);
}

// ---------------------------------------------------------------------------
// OCCT statics.
// ---------------------------------------------------------------------------

/// OCCT static IsLinear(anEdge, aLine) (cxx L35-76).
pub(crate) fn is_linear(the_edge: &Shape, a_line: &mut Option<(DVec3, DVec3)>) -> bool {
    // OCCT L36-42: the curve read + the trimmed basis unwrap.
    let Some((a_curve_full, _, _)) = brep_tool_curve(the_edge) else {
        return false;
    };
    let mut a_curve = a_curve_full;
    if let Curve3::Trimmed(tc) = &a_curve {
        // OCCT L43-45: aCurve = BasisCurve().
        a_curve = (*tc.curve).clone();
    }

    // OCCT L47-56: Geom_Line.
    if let Curve3::Line(line) = &a_curve {
        // OCCT L51: aLine = Lin() — the (origin, direction) pair form.
        *a_line = Some((line.origin, line.direction));
        return true;
    }
    // OCCT L57-67: Geom_BezierCurve with 2 poles.
    if let Curve3::Bezier(bezier) = &a_curve {
        if bezier.control_points.len() == 2 {
            let pnt1 = bezier.control_points[0];
            let pnt2 = bezier.control_points[1];
            // OCCT L63: gce_MakeLin(Pnt1, Pnt2).
            *a_line = Some((pnt1, (pnt2 - pnt1).normalize()));
            return true;
        }
    }
    // OCCT L68-80: Geom_BSplineCurve with 2 poles.
    if let Curve3::BSpline(bspline) = &a_curve {
        if bspline.control_points.len() == 2 {
            let pnt1 = bspline.control_points[0];
            let pnt2 = bspline.control_points[1];
            // OCCT L74: gce_MakeLin(Pnt1, Pnt2).
            *a_line = Some((pnt1, (pnt2 - pnt1).normalize()));
            return true;
        }
    }

    // OCCT L83.
    false
}

/// OCCT static TypeOfEdge(anEdge) (cxx L78-89).
pub(crate) fn type_of_edge(the_edge: &Shape) -> GeomAbsCurveType {
    // OCCT L80-84.
    let mut a_lin: Option<(DVec3, DVec3)> = None;
    if is_linear(the_edge, &mut a_lin) {
        return GeomAbsCurveType::Line;
    }

    // OCCT L86-88: BRepAdaptor_Curve::GetType — the variant discriminant.
    if let Some((a_curve, _, _)) = brep_tool_curve(the_edge) {
        return match a_curve {
            Curve3::Circle(_) => GeomAbsCurveType::Circle,
            Curve3::Ellipse(_) => GeomAbsCurveType::Ellipse,
            Curve3::Hyperbola(_) => GeomAbsCurveType::Hyperbola,
            Curve3::Parabola(_) => GeomAbsCurveType::Parabola,
            Curve3::Bezier(_) => GeomAbsCurveType::BezierCurve,
            Curve3::BSpline(_) => GeomAbsCurveType::BSplineCurve,
            Curve3::Offset(_) => GeomAbsCurveType::OffsetCurve,
            Curve3::SineWave(_) => GeomAbsCurveType::OtherCurve,
            Curve3::CircularHelix(_) => GeomAbsCurveType::OtherCurve,
            _ => GeomAbsCurveType::OtherCurve,
        };
    }
    GeomAbsCurveType::OtherCurve
}

/// OCCT static TangentOfEdge(aShape, OnFirst) (cxx L91-119).
pub(crate) fn tangent_of_edge(a_shape: &Shape, on_first: bool) -> DVec3 {
    // OCCT L92-94.
    let an_edge = a_shape;
    let an_or = an_edge.orientation;

    // OCCT L96-105.
    let Some((a_curve, fpar, lpar)) = brep_tool_curve(an_edge) else {
        return DVec3::ZERO;
    };
    let the_par = if on_first {
        if an_or == rcad_kernel::topods::Orientation::Forward {
            fpar
        } else {
            lpar
        }
    } else if an_or == rcad_kernel::topods::Orientation::Forward {
        lpar
    } else {
        fpar
    };

    // OCCT L107-110: aCurve->D1(thePar, thePoint, theTangent).
    let the_tangent = a_curve.tangent_at(the_par);
    // OCCT L111-113.
    if an_or == rcad_kernel::topods::Orientation::Reversed {
        return -the_tangent;
    }

    // OCCT L115.
    the_tangent
}

/// OCCT static IsValidEdge(theEdge, theFace) (cxx L121-155).
pub(crate) fn is_valid_edge(the_edge: &Shape, the_face: &Shape) -> bool {
    // OCCT L123-125.
    let (v1, v2) = top_exp_vertices_raw(the_edge);
    let v1 = v1.unwrap_or_else(Shape::null);
    let v2 = v2.unwrap_or_else(Shape::null);

    // OCCT L127-129.
    let tol = rcad_kernel::precision::CONFUSION;

    // OCCT L131-152.
    for an_edge in explorer_shapes(the_face, ShapeType::Edge) {
        // OCCT L134-135: BRepExtrema_DistShapeShape DistMini(theEdge,
        // anEdge).
        let dist_mini = BRepExtremaDistShapeShape::new(the_edge, &an_edge);
        if dist_mini.value() <= tol {
            // OCCT L137-151.
            for i in 1..=dist_mini.nb_solution() {
                let the_type = dist_mini.support_type_shape2(i);
                if the_type == BRepExtremaSupportType::IsOnEdge {
                    return false;
                }
                // theType is "IsVertex" (OCCT L144).
                let a_vertex = dist_mini.support_on_shape2(i);
                if !(a_vertex.is_same(&v1) || a_vertex.is_same(&v2)) {
                    return false;
                }
            }
        }
    }

    // OCCT L154.
    true
}

/// OCCT BRepTools::OuterWire(F) — the rcad TFaceData outer-wire read.
pub(crate) fn outer_wire(the_face: &Shape) -> Shape {
    if let rcad_kernel::topo::topods::TShape::Face(fd) = the_face.data.as_ref() {
        fd.outer_wire.clone()
    } else {
        Shape::null()
    }
}

/// OCCT TopExp_Explorer walk (the tool.rs explorer carrier).
fn explorer_shapes(the_s: &Shape, to_find: ShapeType) -> Vec<Shape> {
    crate::brep_algo::tool::explorer(the_s, to_find, ShapeType::Shape)
}

/// OCCT static GetUnifiedWire(theWire, theUnifier) (cxx L263-293).
pub(crate) fn get_unified_wire(
    the_wire: &Shape,
    the_unifier: &mut ShapeUpgradeUnifySameDomain,
) -> Shape {
    // OCCT L266: BRepLib_MakeWire aWMaker.
    let mut a_wmaker = crate::brep_algo::normal_projection::BRepLibMakeWire::new();
    // OCCT L267-269: BRepTools_WireExplorer wexp + aGeneratedEdges.
    let mut wexp = crate::offset::brep_offset_inter2d::BRepToolsWireExplorer::new();
    wexp.init(the_wire, &Shape::null());
    let mut a_generated_edges: HashMap<u64, Shape> = HashMap::new();
    // OCCT L269-291.
    while wexp.more() {
        let an_edge = wexp.current();
        // OCCT L272: theUnifier.History()->Modified(anEdge) — the History()
        // handle of the algorithm (never null here: the constructor
        // installs its own history; the rcad None is the null-handle case
        // the OCCT dereference would fault on).
        let a_ls = the_unifier
            .history()
            .expect("ShapeUpgrade_UnifySameDomain::History() is null")
            .modified(&an_edge);
        if !a_ls.is_empty() {
            // OCCT L274-285: wire shouldn't contain duplicated generated
            // edges.
            for a_shape in a_ls {
                if !a_generated_edges.contains_key(&a_shape.ptr_id()) {
                    a_generated_edges.insert(a_shape.ptr_id(), a_shape.clone());
                    // OCCT L280: aWMaker.Add(TopoDS::Edge(aShape)).
                    a_wmaker.add(&[a_shape]);
                }
            }
        } else {
            // OCCT L287-289: no change, put original edge.
            a_wmaker.add(&[an_edge]);
        }
        wexp.next();
    }
    // OCCT L292.
    a_wmaker.shape().clone()
}

/// OCCT BRepOffsetAPI_MiddlePath (hxx L30-110).
pub struct BRepOffsetAPIMiddlePath {
    // OCCT BRepBuilderAPI base members.
    pub(crate) my_done: bool, // OCCT BRepBuilderAPI_Command: myDone
    #[allow(dead_code)]
    pub(crate) my_shape: Shape, // OCCT: myShape
    // The rcad arena stand-in (architecture difference #4 precedent).
    pub(crate) my_brep: rcad_kernel::topo::topods::BRep,
    // OCCT private members (hxx L104-110).
    pub(crate) my_initial_shape: Shape,   // OCCT: myInitialShape
    pub(crate) my_start_wire: Shape,      // OCCT: myStartWire
    pub(crate) my_end_wire: Shape,        // OCCT: myEndWire
    pub(crate) my_closed_section: bool,   // OCCT: myClosedSection
    #[allow(dead_code)]
    pub(crate) my_closed_ring: bool,      // OCCT: myClosedRing
    pub(crate) my_start_wire_edges: HashMap<u64, Shape>, // OCCT: myStartWireEdges
    pub(crate) my_end_wire_edges: HashMap<u64, Shape>,   // OCCT: myEndWireEdges
    pub(crate) my_paths: Vec<Vec<Shape>>, // OCCT: myPaths
    /// The ChooseEdge == 1 restart carrier (the OCCT `j = 1` form of cxx
    /// L712 — the for-counter reset).
    #[allow(dead_code)]
    pub(crate) restart_j: bool,
}

impl BRepOffsetAPIMiddlePath {
    /// OCCT BRepOffsetAPI_MiddlePath::BRepOffsetAPI_MiddlePath(aShape,
    /// StartShape, EndShape) (cxx L297-331).
    pub fn new(a_shape: &Shape, start_shape: &Shape, end_shape: &Shape) -> Self {
        // OCCT L301-303: ShapeUpgrade_UnifySameDomain Unifier(aShape) — the
        // hxx L81-86 defaults UnifyEdges = true, UnifyFaces = true,
        // ConcatBSplines = false; Unifier.Build(); myInitialShape =
        // Unifier.Shape().  The real Build is the arena form (cxx
        // L4454-4469): the pool created here becomes the class my_brep
        // arena (architecture difference #4).
        let mut my_brep = rcad_kernel::topo::topods::BRep::new();
        let mut unifier = ShapeUpgradeUnifySameDomain::with_shape(a_shape, true, true, false);
        unifier.build(&mut my_brep);
        let my_initial_shape = unifier.shape().clone();

        // OCCT L305-318: the start/end wire extraction.
        let a_start_wire = if start_shape.shape_type() == ShapeType::Face {
            // OCCT L308-310: BRepTools::OuterWire(StartFace).
            outer_wire(start_shape)
        } else {
            start_shape.clone()
        };

        let an_end_wire = if end_shape.shape_type() == ShapeType::Face {
            // OCCT L315-317.
            outer_wire(end_shape)
        } else {
            end_shape.clone()
        };

        // OCCT L326-327.
        let my_start_wire = get_unified_wire(&a_start_wire, &mut unifier);
        let my_end_wire = get_unified_wire(&an_end_wire, &mut unifier);

        // OCCT L329-330.
        let my_closed_section = crate::brep_algo::tool::shape_is_closed(&my_start_wire);
        let my_closed_ring = my_start_wire.is_same(&my_end_wire);

        BRepOffsetAPIMiddlePath {
            my_done: false,
            my_shape: Shape::null(),
            my_brep,
            my_initial_shape,
            my_start_wire,
            my_end_wire,
            my_closed_section,
            my_closed_ring,
            my_start_wire_edges: HashMap::new(),
            my_end_wire_edges: HashMap::new(),
            my_paths: Vec::new(),
            restart_j: false,
        }
    }

    /// OCCT BRep_Builder::MakeEdge(V1, V2) form (the BRepLib_MakeEdge
    /// carrier of Build).
    pub(crate) fn brep_lib_make_edge_vertices(the_v1: &Shape, the_v2: &Shape) -> Shape {
        let mut e = builder_make_edge();
        builder_add_edge_vertex(&mut e, the_v1);
        builder_add_edge_vertex(&mut e, the_v2);
        e
    }

    /// OCCT BRepLib_MakeEdge(C, f, l) form.
    pub(crate) fn brep_lib_make_edge_curve_range(
        brep: &mut rcad_kernel::topo::topods::BRep,
        the_curve: &Curve3,
        the_f: f64,
        the_l: f64,
    ) -> Shape {
        let null = Shape::null();
        brep.add_tedge(Some(the_curve.clone()), null.clone(), null, [the_f, the_l])
    }

    /// OCCT BRepLib_MakeEdge(P1, P2) — the segment line edge.
    pub(crate) fn brep_lib_make_edge_of_points(
        brep: &mut rcad_kernel::topo::topods::BRep,
        the_p1: &glam::DVec3,
        the_p2: &glam::DVec3,
    ) -> Shape {
        let mut v1 = builder_make_vertex();
        crate::brep_algo::tool::builder_update_vertex_point_tol(&mut v1, *the_p1, rcad_kernel::precision::CONFUSION);
        let mut v2 = builder_make_vertex();
        crate::brep_algo::tool::builder_update_vertex_point_tol(&mut v2, *the_p2, rcad_kernel::precision::CONFUSION);
        let line = rcad_kernel::geom::Line3::new(*the_p1, (*the_p2 - *the_p1).normalize_or_zero());
        brep.add_tedge(
            Some(Curve3::Line(line)),
            v1,
            v2,
            [0.0, the_p1.distance(*the_p2)],
        )
    }

    /// OCCT BRepLib_MakeEdge(C2d, S, V1, V2, f, l) — the pcurve-on-surface
    /// edge form.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn brep_lib_make_edge_pcurve(
        brep: &mut rcad_kernel::topo::topods::BRep,
        the_pcurve: &rcad_kernel::geom::Curve2d,
        _the_surface: &rcad_kernel::geom::Surface3,
        the_v1: &Shape,
        the_v2: &Shape,
        the_f: f64,
        the_l: f64,
    ) -> Shape {
        let mut e = brep.add_tedge(None, the_v1.clone(), the_v2.clone(), [the_f, the_l]);
        if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut e.data) {
            ed.pcurves.insert((0, 0), (the_pcurve.clone(), the_f, the_l));
        }
        e
    }
}
