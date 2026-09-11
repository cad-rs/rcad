//! OCCT BRepOffsetAPI_ThruSections — CreateSmoothed / TotalSurf /
//! EdgeToBSpline / Generated + the parameter accessors.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_ThruSections.cxx L706-1116 (CreateSmoothed),
//!         L1121-1187 (EdgeToBSpline), L1190-1382 (TotalSurf),
//!         L1385-1740 (Generated) + the accessors (L1620-1747).
//! The class body lives in brep_offset_api_thru_sections.rs.
//!
//! Architecture differences:
//! 1. occ::handle<Geom_BSplineSurface> -> the GeomBSplineSurface carrier
//!    below: the constructor (the anApprox.Surf* outputs) and Copy/IsNull
//!    keep the OCCT form over the landed rcad BSplineSurface data (the full
//!    knot vector is the knots x mults expansion); the used
//!    Geom_BSplineSurface methods (LocateU/UKnot/VKnot/Segment/Bounds/VIso/
//!    UIso and the knot-index accessors) are the TKGeomBase GAP leaves —
//!    they sit behind the TotalSurf null-surface exit that the AppSurf GAP
//!    keeps open (the OCCT failure path).
//! 2. GeomFill_SectionGenerator -> the GeomFillProfiler carrier of
//!    brep_fill/generator.rs (AddCurve/Perform); GeomFill_Line and
//!    GeomFill_AppSurf -> the local carriers below (the approximation
//!    engine is not translated — IsDone() = false is the OCCT no-surface
//!    exit; the Surf* result accessors are the never-reached branch).
//! 3. GeomConvert_ApproxCurve -> the local GAP carrier (HasResult() = false
//!    keeps the EdgeToBSpline fall-through to CurveToBSplineCurve);
//!    GeomConvert::CurveToBSplineCurve -> the base::convert carrier (the
//!    loc_ope_pipe.rs precedent); GeomConvert_CompCurveToBSplineCurve ->
//!    the draft_modification_1_b.rs carrier (Add is the GAP there).
//! 4. BSplCLib::Reparametrize -> the landed rcad_kernel::math::bspl_lib
//!    reparametrize; Geom_BSplineCurve::Reverse -> the local carrier (the
//!    pole/knot reversal).
//! 5. BRepAdaptor_Surface::GetType == GeomAbs_Plane -> the Surface3 variant
//!    discriminant.
//! 6. BRep_Builder SameParameter/SameRange flag setters, MakeShell/Add(
//!    shell, face), MakeEdge(E, C, Tol) and the Range(E, F, first, last)
//!    pcurve range form -> the local carriers (the tool.rs style).
//! 7. TopExp::Vertices(E, V1, V2) over the TShape children (no pool
//!    roundtrip) — the local top_exp_vertices_local (the generator.rs
//!    oriented-read form).

use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::geom::{Curve2d, Curve3, Line2d, Surface3, BSplineCurve3, BSplineSurface};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use glam::DVec2;

use crate::brep_algo::tool::{
    brep_tool_pnt, brep_tool_tolerance, brep_tool_uv_points, builder_add_face_wire,
    builder_add_wire_edge, builder_make_wire, builder_set_closed, empty_copied, explorer,
    reversed, shape_is_closed, builder_update_edge_pcurve,
};
use crate::brep_fill::generator::{BRepFillThruSectionErrorStatus, GeomFillProfiler};
use crate::offset::brep_offset_api_thru_sections::{
    brep_tool_degenerated,
    make_solid, precise_upar, ApproxParametrizationType, BRepOffsetAPIThruSections,
    TOPABS_IN,
};
use crate::offset::brep_offset_inter2d::BRepToolsWireExplorer;
use crate::offset::brep_offset_inter2d::brep_lib_same_parameter;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

// ---------------------------------------------------------------------------
// GAP carriers (the PerformPlan face maker; architecture difference #5).
// ---------------------------------------------------------------------------

/// OCCT BRepBuilderAPI_MakeFace(W) — the wire-form face maker (the
/// BRepLibMakeFace carrier form of brep_offset_make_simple_offset.rs).
/// GAP status: OCCT 8.0 has no BRepLib_MakeFaceWire class; the wire form is
/// BRepBuilderAPI_MakeFace(W) -> BRepLib_MakeFace::Init(W, OnlyPlane)
/// (TKTopAlgo/BRepLib/BRepLib_MakeFace.cxx, the wire-Init overload).  The
/// real body belongs to the topalgo brep_lib MakeFace/MakeWire sibling
/// batch (topalgo/brep_lib/); until that batch lands, IsDone() = false
/// carries the OCCT failure path of this call site
/// (BRepOffsetAPI_ThruSections::CreateSmoothed).
pub(crate) struct BRepLibMakeFaceWire {
    my_face: Shape, // OCCT: myFace
}

impl BRepLibMakeFaceWire {
    /// OCCT BRepBuilderAPI_MakeFace::BRepBuilderAPI_MakeFace(W) — GAP: the
    /// wire face maker (BRepLib_MakeFace::Init(W, OnlyPlane)) is not
    /// translated yet (the topalgo brep_lib sibling batch); IsDone() =
    /// false carries the OCCT failure path.
    pub(crate) fn from_wire(_the_w: &Shape) -> Self {
        BRepLibMakeFaceWire {
            my_face: Shape::null(),
        }
    }

    /// OCCT BRepBuilderAPI_MakeFace::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        false
    }

    /// OCCT BRepBuilderAPI_MakeFace::Face().
    pub(crate) fn face(&self) -> Shape {
        self.my_face.clone()
    }
}

// ---------------------------------------------------------------------------
// Geom_BSplineSurface carrier (architecture difference #1).
// ---------------------------------------------------------------------------

/// OCCT Geom_BSplineSurface (TKGeomBase/Geom_BSplineSurface.hxx) — the
/// TotalSurf result carrier.  The constructor keeps the OCCT
/// (poles, weights, uknots, vknots, umults, vmults, udeg, vdeg) form over
/// the rcad BSplineSurface data; the used surface methods are the
/// Geom_BSplineSurface GAP leaves (see the module header).
pub struct GeomBSplineSurface {
    pub(crate) my_surface: BSplineSurface, // the rcad data carrier
}

impl GeomBSplineSurface {
    /// OCCT Geom_BSplineSurface(Poles, Weights, UKnots, VKnots, UMults,
    /// VMults, UDegree, VDegree).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_poles(
        the_poles: &Vec<Vec<glam::DVec3>>,
        the_weights: &Vec<Vec<f64>>,
        the_u_knots: &Vec<f64>,
        the_v_knots: &Vec<f64>,
        the_u_mults: &Vec<usize>,
        the_v_mults: &Vec<usize>,
        the_u_degree: usize,
        the_v_degree: usize,
    ) -> Self {
        let mut knots_u: Vec<f64> = Vec::new();
        for (k, m) in the_u_knots.iter().zip(the_u_mults.iter()) {
            for _ in 0..*m {
                knots_u.push(*k);
            }
        }
        let mut knots_v: Vec<f64> = Vec::new();
        for (k, m) in the_v_knots.iter().zip(the_v_mults.iter()) {
            for _ in 0..*m {
                knots_v.push(*k);
            }
        }
        GeomBSplineSurface {
            my_surface: BSplineSurface {
                degree_u: the_u_degree,
                degree_v: the_v_degree,
                knots_u,
                knots_v,
                control_points: the_poles.clone(),
                weights: the_weights.clone(),
            },
        }
    }

    /// OCCT handle copy form (TS->Copy()).
    pub(crate) fn copy(&self) -> Self {
        GeomBSplineSurface {
            my_surface: self.my_surface.clone(),
        }
    }

    /// OCCT IsNull over the handle.
    #[allow(dead_code)]
    pub(crate) fn is_null(&self) -> bool {
        false
    }

    /// OCCT Geom_BSplineSurface::LocateU(U, Eps, I1, I2) — GAP.
    pub(crate) fn locate_u(&self, _the_u: f64, _the_eps: f64, _i1: &mut i32, _i2: &mut i32) {
        panic!("GAP: Geom_BSplineSurface::LocateU (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::UKnot(I) — GAP.
    pub(crate) fn u_knot(&self, _the_i: i32) -> f64 {
        panic!("GAP: Geom_BSplineSurface::UKnot (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::VKnot(I) — GAP.
    pub(crate) fn v_knot(&self, _the_i: i32) -> f64 {
        panic!("GAP: Geom_BSplineSurface::VKnot (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::FirstVKnotIndex() — GAP.
    pub(crate) fn first_v_knot_index(&self) -> i32 {
        panic!("GAP: Geom_BSplineSurface::FirstVKnotIndex (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::LastVKnotIndex() — GAP.
    pub(crate) fn last_v_knot_index(&self) -> i32 {
        panic!("GAP: Geom_BSplineSurface::LastVKnotIndex (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::Segment(U1, U2, V1, V2) — GAP.
    pub(crate) fn segment(&mut self, _u1: f64, _u2: f64, _v1: f64, _v2: f64) {
        panic!("GAP: Geom_BSplineSurface::Segment (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::Bounds(U1, U2, V1, V2) — GAP.
    pub(crate) fn bounds(&self, _u1: &mut f64, _u2: &mut f64, _v1: &mut f64, _v2: &mut f64) {
        panic!("GAP: Geom_BSplineSurface::Bounds (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::VIso(V) — GAP.
    pub(crate) fn v_iso(&self, _the_v: f64) -> Curve3 {
        panic!("GAP: Geom_BSplineSurface::VIso (TKGeomBase not translated)")
    }

    /// OCCT Geom_BSplineSurface::UIso(U) — GAP.
    pub(crate) fn u_iso(&self, _the_u: f64) -> Curve3 {
        panic!("GAP: Geom_BSplineSurface::UIso (TKGeomBase not translated)")
    }
}

/// OCCT GeomFill_Line (TKGeomAlgo/GeomFill_Line.hxx) — the section-count
/// carrier of the AppSurf Perform.
pub(crate) struct GeomFillLine {
    #[allow(dead_code)]
    my_nb: i32, // OCCT: myNb
}

impl GeomFillLine {
    /// OCCT GeomFill_Line(Nb).
    pub(crate) fn new(the_nb: i32) -> Self {
        GeomFillLine { my_nb: the_nb }
    }
}

/// OCCT GeomFill_AppSurf (TKGeomAlgo) — the surface approximation engine of
/// TotalSurf (architecture difference #2; GAP: no rcad translation yet —
/// IsDone() = false carries the OCCT no-surface exit of TotalSurf; the
/// parameter storage keeps the OCCT constructor form; the Surf* result
/// accessors are the never-reached branch behind IsDone).
pub(crate) struct GeomFillAppSurf {
    my_continuity: rcad_kernel::topods::GeomAbsShape, // OCCT: myContinuity
    #[allow(dead_code)]
    my_par_type: ApproxParametrizationType,           // OCCT: myParType
}

impl GeomFillAppSurf {
    /// OCCT GeomFill_AppSurf(Degmin, Degmax, Tol3d, Tol2d, NbIt).
    pub(crate) fn new(
        _the_degmin: i32,
        _the_degmax: i32,
        _the_tol3d: f64,
        _the_tol2d: f64,
        _the_nb_it: i32,
    ) -> Self {
        GeomFillAppSurf {
            my_continuity: rcad_kernel::topods::GeomAbsShape::C0,
            my_par_type: ApproxParametrizationType::ChordLength,
        }
    }

    /// OCCT GeomFill_AppSurf::SetContinuity(C).
    pub(crate) fn set_continuity(&mut self, the_c: rcad_kernel::topods::GeomAbsShape) {
        self.my_continuity = the_c;
    }

    /// OCCT GeomFill_AppSurf::SetCriteriumWeight(W1, W2, W3).
    pub(crate) fn set_criterium_weight(&mut self, _w1: f64, _w2: f64, _w3: f64) {
        // The criterion weights storage (the smoothing branch).
    }

    /// OCCT GeomFill_AppSurf::SetParType(Type).
    pub(crate) fn set_par_type(&mut self, the_type: ApproxParametrizationType) {
        self.my_par_type = the_type;
    }

    /// OCCT GeomFill_AppSurf::PerformSmoothing(Line, Section) — GAP.
    pub(crate) fn perform_smoothing(
        &mut self,
        _line: &GeomFillLine,
        _section: &GeomFillProfiler,
    ) {
        // The smoothing engine is not translated; the approximator stays
        // not-done (the OCCT no-surface exit of TotalSurf).
    }

    /// OCCT GeomFill_AppSurf::Perform(Line, Section, SpApprox) — GAP.
    pub(crate) fn perform(&mut self, _line: &GeomFillLine, _section: &GeomFillProfiler, _sp_approx: bool) {
        // The approximation engine is not translated; the approximator stays
        // not-done (the OCCT no-surface exit of TotalSurf).
    }

    /// OCCT GeomFill_AppSurf::IsDone() — GAP (false).
    pub(crate) fn is_done(&self) -> bool {
        false
    }

    /// OCCT GeomFill_AppSurf::SurfPoles() — GAP (never reached).
    pub(crate) fn surf_poles(&self) -> Vec<Vec<glam::DVec3>> {
        panic!("GAP: GeomFill_AppSurf::SurfPoles (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::SurfWeights() — GAP (never reached).
    pub(crate) fn surf_weights(&self) -> Vec<Vec<f64>> {
        panic!("GAP: GeomFill_AppSurf::SurfWeights (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::SurfUKnots() — GAP (never reached).
    pub(crate) fn surf_u_knots(&self) -> Vec<f64> {
        panic!("GAP: GeomFill_AppSurf::SurfUKnots (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::SurfVKnots() — GAP (never reached).
    pub(crate) fn surf_v_knots(&self) -> Vec<f64> {
        panic!("GAP: GeomFill_AppSurf::SurfVKnots (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::SurfUMults() — GAP (never reached).
    pub(crate) fn surf_u_mults(&self) -> Vec<usize> {
        panic!("GAP: GeomFill_AppSurf::SurfUMults (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::SurfVMults() — GAP (never reached).
    pub(crate) fn surf_v_mults(&self) -> Vec<usize> {
        panic!("GAP: GeomFill_AppSurf::SurfVMults (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::UDegree() — GAP (never reached).
    pub(crate) fn u_degree(&self) -> usize {
        panic!("GAP: GeomFill_AppSurf::UDegree (TKGeomAlgo not translated)")
    }

    /// OCCT GeomFill_AppSurf::VDegree() — GAP (never reached).
    pub(crate) fn v_degree(&self) -> usize {
        panic!("GAP: GeomFill_AppSurf::VDegree (TKGeomAlgo not translated)")
    }
}

/// OCCT GeomConvert_ApproxCurve (TKTopAlgo/GeomConvert_ApproxCurve) — the
/// conic approximator of EdgeToBSpline (architecture difference #3; GAP:
/// HasResult() = false keeps the OCCT fall-through to
/// GeomConvert::CurveToBSplineCurve).
pub(crate) struct GeomConvertApproxCurve;

impl GeomConvertApproxCurve {
    /// OCCT GeomConvert_ApproxCurve(Curve, Tol, Order, MaxSegments,
    /// MaxDegree).
    pub(crate) fn new(
        _the_curve: &Curve3,
        _the_tol: f64,
        _the_order: rcad_kernel::topods::GeomAbsShape,
        _the_max_segments: i32,
        _the_max_degree: i32,
    ) -> Self {
        GeomConvertApproxCurve
    }

    /// OCCT GeomConvert_ApproxCurve::HasResult() — GAP (false).
    pub(crate) fn has_result(&self) -> bool {
        false
    }

    /// OCCT GeomConvert_ApproxCurve::Curve() — GAP.
    pub(crate) fn curve(&self) -> BSplineCurve3 {
        panic!("GAP: GeomConvert_ApproxCurve::Curve (TKTopAlgo not translated)")
    }
}

/// OCCT GeomConvert::CurveToBSplineCurve(Trimmed) — the base::convert
/// carrier (the loc_ope_pipe.rs geomconvert_curve_to_bspline precedent).
pub(crate) fn geom_convert_curve_to_bspline(
    the_c1: &Curve3,
    the_p1: f64,
    the_p2: f64,
) -> Option<BSplineCurve3> {
    use rcad_kernel::base::convert;
    match the_c1 {
        // OCCT CaseLine (GeomConvert.cxx): two poles, knots at the range
        // bounds.
        Curve3::Line(l) => Some(convert::line_to_bspline_range(l, the_p1, the_p2)),
        // The staged rcad conversion (exact for bspline identity).
        other => {
            let _ = (the_p1, the_p2);
            Some(convert::curve_to_bspline(other, 33))
        }
    }
}

/// OCCT Geom_BSplineCurve::Reverse() — the pole/knot reversal (after the
/// [0,1] reparametrization the knot complement is the mirror).
pub(crate) fn bspline_curve_reverse(the_c: &mut BSplineCurve3) {
    the_c.control_points.reverse();
    the_c.weights.reverse();
    let (u1, u2) = (
        the_c.knots.first().cloned().unwrap_or(0.0),
        the_c.knots.last().cloned().unwrap_or(1.0),
    );
    let reversed_knots: Vec<f64> = the_c.knots.iter().rev().map(|k| (u1 + u2) - k).collect();
    the_c.knots = reversed_knots;
}

/// OCCT TopExp::Vertices(E, V1, V2) over the TShape children (the
/// generator.rs oriented-read form without the pool roundtrip).
pub(crate) fn top_exp_vertices_local(e: &Shape) -> (Shape, Shape) {
    if let TShape::Edge(ed) = e.data.as_ref() {
        if e.orientation == Orientation::Reversed {
            (ed.last.clone(), ed.first.clone())
        } else {
            (ed.first.clone(), ed.last.clone())
        }
    } else {
        (Shape::null(), Shape::null())
    }
}

// ---------------------------------------------------------------------------
// BRep_Builder carriers (the tool.rs style).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::MakeEdge(E, C, Tol).
fn brep_builder_make_edge_curve(brep: &mut BRep, the_curve: &Curve3, _the_tol: f64) -> Shape {
    let null = Shape::null();
    brep.add_tedge(Some(the_curve.clone()), null.clone(), null, [0.0, 1.0])
}

/// OCCT BRep_Builder::MakeShell().
fn brep_builder_make_shell(brep: &mut BRep) -> Shape {
    brep.add_tshell(vec![])
}

/// OCCT BRep_Builder::Add(Shell, Face).
fn brep_builder_add_shell_face(the_shell: &mut Shape, the_face: &Shape) {
    if let TShape::Shell(sd) = Arc::make_mut(&mut the_shell.data) {
        sd.faces.push(the_face.clone());
        sd.my_shapes.push(the_face.clone());
    }
}

/// OCCT BRep_Builder::MakeFace(F, S, Tol) — the surface-face form (the wire
/// is added afterwards by Add(F, W), as in OCCT).
fn brep_builder_make_face_surface(
    brep: &mut BRep,
    the_surface: &GeomBSplineSurface,
    the_tol: f64,
) -> Shape {
    let _ = the_tol;
    brep.add_tface(
        Some(Surface3::BSpline(the_surface.my_surface.clone())),
        Shape::null(),
        vec![],
        None,
        None,
        vec![],
        true,
    )
}

/// OCCT BRep_Builder::Range(E, F, First, Last) — the pcurve range form.
fn brep_builder_range_on_face(the_e: &mut Shape, _the_f: &Shape, f: f64, l: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.range = [f, l];
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, F, Tol) — the closed-surface
/// two-pcurve form.
fn brep_builder_update_edge_pcurves(
    the_e: &mut Shape,
    the_c1: &rcad_kernel::geom::Curve2d,
    _the_c2: &rcad_kernel::geom::Curve2d,
    the_f: &Shape,
    the_tol: f64,
) {
    let key = (the_f.ptr_id(), the_f.location);
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        let (f0, l0) = (ed.range[0], ed.range[1]);
        // OCCT stores both pcurves of the closed surface; the rcad pcurve
        // index is keyed by face — the second curve lands with the Geom
        // batch.
        ed.pcurves.insert(key, (the_c1.clone(), f0, l0));
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::SameRange(E, B).
fn brep_builder_set_same_range(the_e: &mut Shape, flag: bool) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.same_range = flag;
    }
}

/// OCCT BRep_Builder::SameParameter(E, B).
fn brep_builder_set_same_parameter(the_e: &mut Shape, flag: bool) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.same_parameter = flag;
    }
}

/// OCCT GeomConvert_CompCurveToBSplineCurve re-host (the
/// draft_modification_1_b.rs carrier; Add is the GAP there).
pub(crate) struct GeomConvertCompCurveToBSplineCurveCarrier {
    inner: crate::offset::draft_modification_1_b::GeomConvertCompCurveToBSplineCurve,
}

impl GeomConvertCompCurveToBSplineCurveCarrier {
    pub(crate) fn new(the_bspline: &BSplineCurve3) -> Self {
        GeomConvertCompCurveToBSplineCurveCarrier {
            inner:
                crate::offset::draft_modification_1_b::GeomConvertCompCurveToBSplineCurve::new(
                    the_bspline,
                ),
        }
    }

    /// OCCT GeomConvert_CompCurveToBSplineCurve::Add(BSpline, Tol, After,
    /// WithCancellation, Confusion) — the draft_modification_1_b carrier Add
    /// is the GAP (in-place, no result); the OCCT Standard_Boolean result is
    /// discarded at the TotalSurf call site.
    pub(crate) fn add(
        &mut self,
        the_bspline: &BSplineCurve3,
        the_tol: f64,
        _after: bool,
        _with_cancellation: bool,
        _confusion: i32,
    ) {
        self.inner.add(the_bspline, the_tol, true)
    }

    /// OCCT GeomConvert_CompCurveToBSplineCurve::BSplineCurve().
    pub(crate) fn bspline_curve(&self) -> BSplineCurve3 {
        self.inner.bspline_curve()
    }
}

/// OCCT BRepLib::SameParameter(E, Tol, NewTol, WithCorrection) — the 4-arg
/// form (the inter2d stub leaf).
pub(crate) fn brep_lib_same_parameter_with_result(
    _the_e: &Shape,
    _the_tol: f64,
    _the_new_tol: &mut f64,
    _the_with_correction: bool,
) {
    // The brep_lib_same_parameter stub form (the inter2d GAP leaf).
}

impl BRepOffsetAPIThruSections {
    /// OCCT BRepOffsetAPI_ThruSections::FirstShape() (cxx L1622-1625).
    pub fn first_shape(&self) -> &Shape {
        &self.my_first
    }

    /// OCCT BRepOffsetAPI_ThruSections::LastShape() (cxx L1627-1630).
    pub fn last_shape(&self) -> &Shape {
        &self.my_last
    }

    /// OCCT BRepOffsetAPI_ThruSections::GeneratedFace(Edge) (cxx
    /// L1632-1642).
    pub fn generated_face(&self, edge: &Shape) -> Shape {
        // OCCT L1634-1641.
        if let Some(face) = self.my_edge_face.get(&edge.ptr_id()) {
            face.clone()
        } else {
            Shape::null()
        }
    }

    /// OCCT BRepOffsetAPI_ThruSections::CriteriumWeight(W1, W2, W3) (cxx
    /// L1646-1652).
    pub fn criterium_weight(&self) -> (f64, f64, f64) {
        (
            self.my_crit_weights[0],
            self.my_crit_weights[1],
            self.my_crit_weights[2],
        )
    }

    /// OCCT BRepOffsetAPI_ThruSections::SetCriteriumWeight(W1, W2, W3) (cxx
    /// L1654-1666).
    pub fn set_criterium_weight(&mut self, w1: f64, w2: f64, w3: f64) {
        if w1 < 0.0 || w2 < 0.0 || w3 < 0.0 {
            self.my_status = BRepFillThruSectionErrorStatus::Failed;
            return;
        }
        self.my_crit_weights[0] = w1;
        self.my_crit_weights[1] = w2;
        self.my_crit_weights[2] = w3;
    }

    /// OCCT BRepOffsetAPI_ThruSections::SetContinuity(TheCont) (cxx
    /// L1668-1671).
    pub fn set_continuity(&mut self, the_cont: rcad_kernel::topods::GeomAbsShape) {
        self.my_continuity = the_cont;
    }

    /// OCCT BRepOffsetAPI_ThruSections::Continuity() (cxx L1673-1676).
    pub fn continuity(&self) -> rcad_kernel::topods::GeomAbsShape {
        self.my_continuity
    }

    /// OCCT BRepOffsetAPI_ThruSections::SetParType(ParType) (cxx
    /// L1678-1681).
    pub fn set_par_type(&mut self, par_type: ApproxParametrizationType) {
        self.my_param_type = par_type;
    }

    /// OCCT BRepOffsetAPI_ThruSections::ParType() (cxx L1683-1686).
    pub fn par_type(&self) -> ApproxParametrizationType {
        self.my_param_type
    }

    /// OCCT BRepOffsetAPI_ThruSections::SetMaxDegree(MaxDeg) (cxx
    /// L1688-1691).
    pub fn set_max_degree(&mut self, max_deg: i32) {
        self.my_deg_max = max_deg;
    }

    /// OCCT BRepOffsetAPI_ThruSections::MaxDegree() (cxx L1693-1696).
    pub fn max_degree(&self) -> i32 {
        self.my_deg_max
    }

    /// OCCT BRepOffsetAPI_ThruSections::SetSmoothing(UseVar) (cxx
    /// L1698-1701).
    pub fn set_smoothing(&mut self, use_var: bool) {
        self.my_use_smoothing = use_var;
    }

    /// OCCT BRepOffsetAPI_ThruSections::UseSmoothing() (cxx L1703-1706).
    pub fn use_smoothing(&self) -> bool {
        self.my_use_smoothing
    }

    /// OCCT BRepOffsetAPI_ThruSections::SetMutableInput(theIsMutableInput)
    /// (cxx L1708-1711).
    pub fn set_mutable_input(&mut self, the_is_mutable_input: bool) {
        self.my_mutable_input = the_is_mutable_input;
    }

    /// OCCT BRepOffsetAPI_ThruSections::IsMutableInput() (cxx L1713-1716).
    pub fn is_mutable_input(&self) -> bool {
        self.my_mutable_input
    }

    /// OCCT BRepOffsetAPI_ThruSections::Wires() (hxx inline).
    pub fn wires(&self) -> &Vec<Shape> {
        &self.my_input_wires
    }

    /// OCCT BRepOffsetAPI_ThruSections::GetStatus() (hxx inline).
    pub fn get_status(&self) -> BRepFillThruSectionErrorStatus {
        self.my_status
    }
}

/// OCCT static EdgeToBSpline(theEdge) (cxx L1121-1187) — gets the curve of
/// the edge and converts it to a bspline parameterized from 0 to 1 (the
/// same function as in BRepFill_NSections.cxx).
pub(crate) fn edge_to_bspline(the_edge: &Shape) -> Option<BSplineCurve3> {
    let mut a_bs_curve: Option<BSplineCurve3> = None;
    // OCCT L1126: the degenerated edge — the point-curve construction.
    if brep_tool_degenerated(the_edge) {
        // OCCT L1128-1140: aPoles(1) = Pnt(vf); aPoles(2) = Pnt(vl).
        let (vl, vf) = top_exp_vertices_local(the_edge);
        let vf_pnt = brep_tool_pnt(&vf).unwrap_or(glam::DVec3::ZERO);
        let vl_pnt = brep_tool_pnt(&vl).unwrap_or(glam::DVec3::ZERO);
        // OCCT L1142: Geom_BSplineCurve(aPoles, aKnots, aMults, 1).
        a_bs_curve = Some(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![vf_pnt, vl_pnt],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        });
    } else {
        // OCCT L1145-1154: the curve of the edge (null curve -> nullptr).
        let Some((a_curve, a_first, a_last)) = crate::brep_algo::tool::brep_tool_curve(the_edge)
        else {
            return None;
        };

        // OCCT L1158-1163: the trimmed-curve conversion (copies, segments
        // and removes periodicity when needed).
        let a_trim_curve = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            a_curve.clone(),
            a_first,
            a_last,
        ));

        // OCCT L1165-1175: the special treatment of conic curves.
        let is_conic = matches!(
            a_curve,
            Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_)
        );
        if is_conic {
            let an_appr = GeomConvertApproxCurve::new(
                &a_trim_curve,
                CONFUSION,
                rcad_kernel::topods::GeomAbsShape::C1,
                16,
                14,
            );
            if an_appr.has_result() {
                a_bs_curve = Some(an_appr.curve());
            }
        }

        // OCCT L1177-1181: the general case.
        if a_bs_curve.is_none() {
            a_bs_curve = geom_convert_curve_to_bspline(&a_trim_curve, a_first, a_last);
        }

        let Some(a_bs) = a_bs_curve.as_mut() else {
            return None;
        };

        // OCCT L1183-1185: the location transform — the rcad located-curve
        // storage keeps the transform (the draft_modification.rs #4 note).
        let _ = &a_bs;

        // OCCT L1187-1190: reparameterize to [0,1].
        let mut a_knots = a_bs.knots.clone();
        rcad_kernel::math::bspl_lib::reparametrize(0.0, 1.0, &mut a_knots);
        a_bs.knots = a_knots;

        // OCCT L1192-1195: reverse the curve if the edge is reversed.
        if the_edge.orientation == Orientation::Reversed {
            bspline_curve_reverse(a_bs);
        }
    }

    a_bs_curve
}

impl BRepOffsetAPIThruSections {
    /// OCCT BRepOffsetAPI_ThruSections::TotalSurf(shapes, NbSects, NbEdges,
    /// w1Point, w2Point, vClosed) (cxx L1190-1382).
    pub(crate) fn total_surf(
        &self,
        shapes: &[Shape],
        nb_sects: usize,
        nb_edges: usize,
        w1_point: bool,
        w2_point: bool,
        v_closed: bool,
    ) -> Option<GeomBSplineSurface> {
        // OCCT L1192: jdeb = 1, jfin = NbSects (the OCCT 1-based bounds).
        let mut jdeb = 1usize;
        let mut jfin = nb_sects;

        // OCCT L1196: GeomFill_SectionGenerator section.
        let mut section = GeomFillProfiler::new();
        // OCCT L1197-1200.
        let mut surface: Option<GeomBSplineSurface> = None;
        let mut bs1: Option<BSplineCurve3> = None;

        // OCCT L1202-1218: the w1Point point curve.
        if w1_point {
            jdeb += 1;
            let edge = shapes[0].clone();
            let (vl, vf) = top_exp_vertices_local(&edge);
            let vf_pnt = brep_tool_pnt(&vf).unwrap_or(glam::DVec3::ZERO);
            let vl_pnt = brep_tool_pnt(&vl).unwrap_or(glam::DVec3::ZERO);
            let bs_point = BSplineCurve3 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![vf_pnt, vl_pnt],
                weights: vec![1.0, 1.0],
                is_periodic: false,
            };
            section.add_curve(&Curve3::BSpline(bs_point));
        }

        // OCCT L1220-1223.
        if w2_point {
            jfin -= 1;
        }

        // OCCT L1225-1290: the section curves.
        for j in jdeb..=jfin {
            // OCCT L1228-1231: the looping-section case.
            if j == jfin && v_closed {
                if let Some(bs1_c) = bs1.clone() {
                    section.add_curve(&Curve3::BSpline(bs1_c));
                }
            } else {
                // OCCT L1234-1241: read the first edge to initialize CompBS.
                let a_prev_edge = shapes[(j - 1) * nb_edges].clone();
                let Some(curv_bs) = edge_to_bspline(&a_prev_edge) else {
                    return None;
                };

                // OCCT L1245-1246: initialization.
                let mut comp_bs = GeomConvertCompCurveToBSplineCurveCarrier::new(&curv_bs);

                // OCCT L1248-1275.
                for i in 2..=nb_edges {
                    // OCCT L1250-1259: read the edge.
                    let a_next_edge = shapes[(j - 1) * nb_edges + (i - 1)].clone();
                    let mut a_tol_v = CONFUSION;
                    let (vf, vl) = top_exp_vertices_local(&a_next_edge);
                    a_tol_v = a_tol_v.max(brep_tool_tolerance(&vf));
                    a_tol_v = a_tol_v.max(brep_tool_tolerance(&vl));
                    a_tol_v = a_tol_v.min(1.0e-3);
                    let Some(next_bs) = edge_to_bspline(&a_next_edge) else {
                        // OCCT L1261-1264.
                        return None;
                    };

                    // OCCT L1266-1267: concatenation.
                    comp_bs.add(&next_bs, a_tol_v, true, false, 1);
                }

                // OCCT L1277-1280: return the final section.
                let bs = comp_bs.bspline_curve();
                section.add_curve(&Curve3::BSpline(bs.clone()));

                // OCCT L1282-1286: the looping-section case.
                if j == jdeb && v_closed {
                    bs1 = Some(bs);
                }
            }
        }

        // OCCT L1292-1310: the w2Point point curve.
        if w2_point {
            let edge = shapes[nb_sects * nb_edges - 1].clone();
            let (vl, vf) = top_exp_vertices_local(&edge);
            let vf_pnt = brep_tool_pnt(&vf).unwrap_or(glam::DVec3::ZERO);
            let vl_pnt = brep_tool_pnt(&vl).unwrap_or(glam::DVec3::ZERO);
            let bs_point = BSplineCurve3 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![vf_pnt, vl_pnt],
                weights: vec![1.0, 1.0],
                is_periodic: false,
            };
            section.add_curve(&Curve3::BSpline(bs_point));
        }

        // OCCT L1312-1313.
        section.perform(rcad_kernel::precision::PCONFUSION);
        // OCCT L1313: GeomFill_Line line = new GeomFill_Line(NbSects).
        let line = GeomFillLine::new(nb_sects as i32);

        // OCCT L1315-1319.
        let mut nb_it = 3;
        if self.my_pres3d <= 1.0e-3 {
            nb_it = 0;
        }

        // OCCT L1321-1324.
        let degmin = 2i32;
        let degmax = self.my_deg_max.max(degmin);
        let sp_approx = true;

        // OCCT L1326-1330.
        let mut an_approx =
            GeomFillAppSurf::new(degmin, degmax, self.my_pres3d, self.my_pres3d, nb_it);
        an_approx.set_continuity(self.my_continuity);

        if self.my_use_smoothing {
            // OCCT L1333-1335.
            an_approx.set_criterium_weight(
                self.my_crit_weights[0],
                self.my_crit_weights[1],
                self.my_crit_weights[2],
            );
            an_approx.perform_smoothing(&line, &section);
        } else {
            // OCCT L1337-1340.
            an_approx.set_par_type(self.my_param_type);
            an_approx.perform(&line, &section, sp_approx);
        }

        // OCCT L1342-1355.
        if an_approx.is_done() {
            surface = Some(GeomBSplineSurface::from_poles(
                &an_approx.surf_poles(),
                &an_approx.surf_weights(),
                &an_approx.surf_u_knots(),
                &an_approx.surf_v_knots(),
                &an_approx.surf_u_mults(),
                &an_approx.surf_v_mults(),
                an_approx.u_degree(),
                an_approx.v_degree(),
            ));
        }

        // OCCT L1357.
        surface
    }

    /// OCCT BRepOffsetAPI_ThruSections::CreateSmoothed() (cxx L706-1116).
    pub(crate) fn create_smoothed(&mut self) {
        // OCCT L708-709.
        let nb_sects = self.my_wires.len();

        // OCCT L711-718: check if the first wire is punctual.
        let mut w1_point = true;
        {
            let wire = self.my_wires[0].clone();
            let mut an_exp = BRepToolsWireExplorer::new();
            an_exp.init(&wire, &Shape::null());
            while an_exp.more() {
                w1_point = w1_point && brep_tool_degenerated(&an_exp.current());
                an_exp.next();
            }
        }

        // OCCT L720-727: check if the last wire is punctual.
        let mut w2_point = true;
        {
            let wire = self.my_wires[nb_sects - 1].clone();
            let mut an_exp = BRepToolsWireExplorer::new();
            an_exp.init(&wire, &Shape::null());
            while an_exp.more() {
                w2_point = w2_point && brep_tool_degenerated(&an_exp.current());
                an_exp.next();
            }
        }

        // OCCT L729-734: check if the first wire is the same as the last.
        let mut v_closed = false;
        if self.my_wires[0].is_same(&self.my_wires[self.my_wires.len() - 1]) {
            v_closed = true;
        }

        // OCCT L736-751: find the dimension.
        let mut nb_edges = 0usize;
        if !w1_point {
            let wire = self.my_wires[0].clone();
            let mut an_exp = BRepToolsWireExplorer::new();
            an_exp.init(&wire, &Shape::null());
            while an_exp.more() {
                nb_edges += 1;
                an_exp.next();
            }
        } else {
            let wire = self.my_wires[1].clone();
            let mut an_exp = BRepToolsWireExplorer::new();
            an_exp.init(&wire, &Shape::null());
            while an_exp.more() {
                nb_edges += 1;
                an_exp.next();
            }
        }

        // OCCT L753-787: recover the shapes.
        let mut u_closed = true;
        let mut shapes: Vec<Shape> = vec![Shape::null(); nb_sects * nb_edges];
        let mut nb = 0usize;

        for i in 1..=nb_sects {
            let wire = self.my_wires[i - 1].clone();
            if !shape_is_closed(&wire) {
                // OCCT L759-764: check if the vertices are the same.
                let (v1, v2) = top_exp_vertices_local(&wire);
                if !v1.is_same(&v2) {
                    u_closed = false;
                }
            }
            if (i == 1 && w1_point) || (i == nb_sects && w2_point) {
                // OCCT L766-774: if the wire is punctual.
                let wire = self.my_wires[i - 1].clone();
                let mut an_exp = BRepToolsWireExplorer::new();
                an_exp.init(&wire, &Shape::null());
                for _j in 1..=nb_edges {
                    nb += 1;
                    shapes[nb - 1] = an_exp.current();
                    an_exp.next();
                }
            } else {
                // OCCT L775-783: otherwise.
                let wire = self.my_wires[i - 1].clone();
                let mut an_exp = BRepToolsWireExplorer::new();
                an_exp.init(&wire, &Shape::null());
                while an_exp.more() {
                    nb += 1;
                    shapes[nb - 1] = an_exp.current();
                    an_exp.next();
                }
            }
        }

        // OCCT L789-800: create the new surface.
        let mut shell = brep_builder_make_shell(&mut self.my_brep);
        let mut couture = Shape::null();
        let mut vcouture: Vec<Shape> = vec![Shape::null(); nb_edges];

        // OCCT L803-807: newW1/newW2.
        let mut new_w1 = builder_make_wire();
        let mut new_w2 = builder_make_wire();

        // OCCT L809-811: the points grid declaration (the OCCT leftover
        // declaration — the array is not consumed further).
        let nb_pnts = 21;
        let _ = nb_pnts;

        // OCCT L813-818: concatenate each section to get a total surface.
        let ts = self.total_surf_impl(&shapes, nb_sects, nb_edges, w1_point, w2_point, v_closed);

        // OCCT L820-825.
        let Some(ts) = ts else {
            self.my_status = BRepFillThruSectionErrorStatus::Failed;
            return;
        };

        // OCCT L827-1005: the per-edge face construction.
        for i in 1..=nb_edges {
            // OCCT L830-841: segmentation of TS.
            let mut surface = ts.copy();
            let mut ui1 = (i - 1) as f64;
            let mut ui2 = i as f64;
            ui1 = precise_upar(ui1, &surface);
            ui2 = precise_upar(ui2, &surface);
            let v0 = surface.v_knot(surface.first_v_knot_index());
            let v1 = surface.v_knot(surface.last_v_knot_index());
            surface.segment(ui1, ui2, v0, v1);
            let surface_ref = &surface;

            // OCCT L843-851: return vertices.
            let edge = shapes[i - 1].clone();
            let (mut v1f, mut v1l) = top_exp_vertices_local(&edge);
            if edge.orientation == Orientation::Reversed {
                std::mem::swap(&mut v1f, &mut v1l);
            }
            let first_edge = edge;

            // OCCT L853-860.
            let edge2s = shapes[(nb_sects - 1) * nb_edges + (i - 1)].clone();
            let (mut v2f, mut v2l) = top_exp_vertices_local(&edge2s);
            if edge2s.orientation == Orientation::Reversed {
                std::mem::swap(&mut v2f, &mut v2l);
            }

            // OCCT L862-866: make the face and the wire.
            let mut face =
                brep_builder_make_face_surface(&mut self.my_brep, surface_ref, CONFUSION);
            let mut w = builder_make_wire();

            // OCCT L868-871: make the missing edges.
            let mut f1 = 0.0;
            let mut l1 = 0.0;
            let mut f2 = 0.0;
            let mut l2 = 0.0;
            surface_ref.bounds(&mut f1, &mut l1, &mut f2, &mut l2);

            // OCCT L873-890: --- edge 1.
            let mut edge1: Shape;
            if w1_point {
                // OCCT L875-880: copy the degenerated edge.
                edge1 = empty_copied(&shapes[0]);
                edge1.orientation = Orientation::Forward;
            } else {
                // OCCT L882-883: B.MakeEdge(edge1, surface->VIso(f2), tol).
                let iso = surface_ref.v_iso(f2);
                edge1 = brep_builder_make_edge_curve(&mut self.my_brep, &iso, CONFUSION);
            }
            // OCCT L884-887.
            v1f.orientation = Orientation::Forward;
            crate::brep_algo::tool::builder_add_edge_vertex(&mut edge1, &v1f);
            v1l.orientation = Orientation::Reversed;
            crate::brep_algo::tool::builder_add_edge_vertex(&mut edge1, &v1l);
            brep_builder_range_on_face(&mut edge1, &face, f1, l1);
            // OCCT L888-890: processing of looping sections.
            if v_closed {
                vcouture[i - 1] = edge1.clone();
            }

            // OCCT L892-918: --- edge 2.
            let mut edge2: Shape;
            if v_closed {
                edge2 = vcouture[i - 1].clone();
            } else {
                if w2_point {
                    // OCCT L899-904: copy of the degenerated edge.
                    edge2 = empty_copied(&shapes[nb_sects * nb_edges - 1]);
                    edge2.orientation = Orientation::Forward;
                } else {
                    // OCCT L906-908.
                    let iso = surface_ref.v_iso(l2);
                    edge2 = brep_builder_make_edge_curve(&mut self.my_brep, &iso, CONFUSION);
                }
                v2f.orientation = Orientation::Forward;
                crate::brep_algo::tool::builder_add_edge_vertex(&mut edge2, &v2f);
                v2l.orientation = Orientation::Reversed;
                crate::brep_algo::tool::builder_add_edge_vertex(&mut edge2, &v2l);
                brep_builder_range_on_face(&mut edge2, &face, f1, l1);
            }
            // OCCT L917: edge2.Reverse().
            let mut edge2 = reversed(&edge2);

            // OCCT L919-943: --- edge 3.
            let mut edge3: Shape;
            let mut edge4: Shape = Shape::null();
            if i == 1 {
                let iso = surface_ref.u_iso(f1);
                edge3 = brep_builder_make_edge_curve(&mut self.my_brep, &iso, CONFUSION);
                v1f.orientation = Orientation::Forward;
                crate::brep_algo::tool::builder_add_edge_vertex(&mut edge3, &v1f);
                v2f.orientation = Orientation::Reversed;
                crate::brep_algo::tool::builder_add_edge_vertex(&mut edge3, &v2f);
                brep_builder_range_on_face(&mut edge3, &face, f2, l2);
                if u_closed {
                    couture = edge3.clone();
                }
            } else {
                // OCCT L938-941: edge3 = edge4.
                edge3 = edge4.clone();
            }
            // OCCT L942: edge3.Reverse().
            let mut edge3 = reversed(&edge3);

            // OCCT L944-958: --- edge 4.
            if u_closed && i == nb_edges {
                edge4 = couture.clone();
            } else {
                let iso = surface_ref.u_iso(l1);
                edge4 = brep_builder_make_edge_curve(&mut self.my_brep, &iso, CONFUSION);
                v1l.orientation = Orientation::Forward;
                crate::brep_algo::tool::builder_add_edge_vertex(&mut edge4, &v1l);
                v2l.orientation = Orientation::Reversed;
                crate::brep_algo::tool::builder_add_edge_vertex(&mut edge4, &v2l);
                brep_builder_range_on_face(&mut edge4, &face, f2, l2);
            }

            // OCCT L960-964.
            builder_add_face_wire(&mut w, &edge1);
            builder_add_face_wire(&mut w, &edge4);
            builder_add_face_wire(&mut w, &edge2);
            builder_add_face_wire(&mut w, &edge3);

            // OCCT L966-993: set PCurve.
            if v_closed {
                let p1 = Curve2d::Line(Line2d::new(DVec2::new(0.0, f2), DVec2::new(1.0, 0.0)));
                let p2 = Curve2d::Line(Line2d::new(DVec2::new(0.0, l2), DVec2::new(1.0, 0.0)));
                brep_builder_update_edge_pcurves(&mut edge1, &p1, &p2, &face, CONFUSION);
                brep_builder_range_on_face(&mut edge1, &face, f1, l1);
            } else {
                let p1 = Curve2d::Line(Line2d::new(DVec2::new(0.0, f2), DVec2::new(1.0, 0.0)));
                builder_update_edge_pcurve(&mut edge1, &p1, &face, CONFUSION);
                brep_builder_range_on_face(&mut edge1, &face, f1, l1);
                let p2 = Curve2d::Line(Line2d::new(DVec2::new(0.0, l2), DVec2::new(1.0, 0.0)));
                builder_update_edge_pcurve(&mut edge2, &p2, &face, CONFUSION);
                brep_builder_range_on_face(&mut edge2, &face, f1, l1);
            }

            if u_closed && nb_edges == 1 {
                let p1 = Curve2d::Line(Line2d::new(DVec2::new(l1, 0.0), DVec2::new(0.0, 1.0)));
                let p2 = Curve2d::Line(Line2d::new(DVec2::new(f1, 0.0), DVec2::new(0.0, 1.0)));
                brep_builder_update_edge_pcurves(&mut edge3, &p1, &p2, &face, CONFUSION);
                brep_builder_range_on_face(&mut edge3, &face, f2, l2);
            } else {
                let p1 = Curve2d::Line(Line2d::new(DVec2::new(f1, 0.0), DVec2::new(0.0, 1.0)));
                builder_update_edge_pcurve(&mut edge3, &p1, &face, CONFUSION);
                brep_builder_range_on_face(&mut edge3, &face, f2, l2);
                let p2 = Curve2d::Line(Line2d::new(DVec2::new(l1, 0.0), DVec2::new(0.0, 1.0)));
                builder_update_edge_pcurve(&mut edge4, &p2, &face, CONFUSION);
                brep_builder_range_on_face(&mut edge4, &face, f2, l2);
            }
            // OCCT L994-995: B.Add(face, W); B.Add(shell, face).
            builder_add_face_wire(&mut face, &w);
            brep_builder_add_shell_face(&mut shell, &face);

            // OCCT L997-1003: complete newW1 newW2 + the history.
            let edge12 = reversed(&edge1);
            let edge22 = reversed(&edge2);
            builder_add_wire_edge(&mut new_w1, &edge12);
            builder_add_wire_edge(&mut new_w2, &edge22);

            self.my_edge_face.insert(first_edge.ptr_id(), face);
        }

        // OCCT L1007-1010.
        if u_closed && w1_point && w2_point {
            builder_set_closed(&mut shell, true);
        }

        // OCCT L1012-1037.
        if self.my_is_solid {
            if v_closed {
                // OCCT L1016-1032.
                let mut solid = brep_builder_make_solid_of_shell(&mut self.my_brep, &shell);

                // verify the orientation of the solid.
                let mut clas3d =
                    crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(
                        &solid,
                    );
                clas3d.perform_infinite_point(CONFUSION);
                if clas3d.state() == TOPABS_IN {
                    let local = reversed(&shell);
                    solid = brep_builder_make_solid_of_shell(&mut self.my_brep, &local);
                }
                self.my_shape = solid;
            } else {
                // OCCT L1036: MakeSolid(shell, newW1, newW2, myPres3d, ...).
                let mut shell_mut = shell.clone();
                let mut face1 = self.my_first.clone();
                let mut face2 = self.my_last.clone();
                self.my_shape = make_solid(
                    &mut self.my_brep,
                    &mut shell_mut,
                    &new_w1,
                    &new_w2,
                    self.my_pres3d,
                    &mut face1,
                    &mut face2,
                );
                self.my_first = face1;
                self.my_last = face2;
            }
        } else {
            // OCCT L1042-1044.
            self.my_shape = shell;
        }
        // OCCT L1046: Done().
        self.my_done = true;

        // OCCT L1048-1084: the tolerance pass.
        let mut a_vertex_tolerance_map: HashMap<u64, (Shape, f64)> = HashMap::new();
        for a_cur_edge in explorer(&self.my_shape.clone(), ShapeType::Edge, ShapeType::Shape) {
            let mut a_cur_edge = a_cur_edge;
            // OCCT L1053-1055: B.SameRange(E, false); B.SameParameter(E, false).
            brep_builder_set_same_range(&mut a_cur_edge, false);
            brep_builder_set_same_parameter(&mut a_cur_edge, false);
            // OCCT L1056.
            let a_tolerance = brep_tool_tolerance(&a_cur_edge);
            if self.my_mutable_input {
                // OCCT L1058: BRepLib::SameParameter(aCurEdge, aTolerance).
                brep_lib_same_parameter(&a_cur_edge, a_tolerance);
            } else {
                // OCCT L1060-1081: all edges of myShape can be safely updated;
                // all vertices of myShape are part of the original wires.
                let mut a_new_tolerance = -1.0;
                brep_lib_same_parameter_with_result(&a_cur_edge, a_tolerance, &mut a_new_tolerance, true);
                if a_new_tolerance > 0.0 {
                    let (a_vertex1, a_vertex2) = top_exp_vertices_local(&a_cur_edge);
                    if !a_vertex1.is_null() {
                        let update = match a_vertex_tolerance_map.get(&a_vertex1.ptr_id()) {
                            Some((_, old)) => *old < a_new_tolerance,
                            None => true,
                        };
                        if update {
                            a_vertex_tolerance_map
                                .insert(a_vertex1.ptr_id(), (a_vertex1.clone(), a_new_tolerance));
                        }
                    }
                    if !a_vertex2.is_null() {
                        let update = match a_vertex_tolerance_map.get(&a_vertex2.ptr_id()) {
                            Some((_, old)) => *old < a_new_tolerance,
                            None => true,
                        };
                        if update {
                            a_vertex_tolerance_map
                                .insert(a_vertex2.ptr_id(), (a_vertex2.clone(), a_new_tolerance));
                        }
                    }
                }
            }
        }

        // OCCT L1086-1112.
        if !self.my_mutable_input {
            let mut a_reshaper = ShapeBuildReShape::new();
            for (_key, (a_vertex, a_new_tolerance)) in a_vertex_tolerance_map.clone() {
                if brep_tool_tolerance(&a_vertex) < a_new_tolerance {
                    // OCCT L1093-1096: aNnewVertex = EmptyCopied();
                    // B.UpdateVertex(aNnewVertex, aNewTolerance).
                    let mut a_new_vertex = empty_copied(&a_vertex);
                    crate::brep_algo::tool::builder_update_vertex_tol(
                        &mut a_new_vertex,
                        a_new_tolerance,
                    );
                    // OCCT L1097.
                    a_reshaper.replace(&mut self.my_brep, &a_vertex, &a_new_vertex);
                }
            }
            // OCCT L1100: myShape = aReshaper.Apply(myShape).
            self.my_shape = a_reshaper.apply(
                &mut self.my_brep,
                &self.my_shape.clone(),
                ShapeType::Shape,
            );
        }
    }

    /// The TotalSurf call form of CreateSmoothed (cxx L818) — the private
    /// method dispatch.
    fn total_surf_impl(
        &self,
        shapes: &[Shape],
        nb_sects: usize,
        nb_edges: usize,
        w1_point: bool,
        w2_point: bool,
        v_closed: bool,
    ) -> Option<GeomBSplineSurface> {
        self.total_surf(shapes, nb_sects, nb_edges, w1_point, w2_point, v_closed)
    }
}

/// OCCT BRep_Builder::MakeSolid(S) + Add(Solid, Shell) — the shell-only
/// solid form of the vClosed branch.
pub(crate) fn brep_builder_make_solid_of_shell(brep: &mut BRep, the_shell: &Shape) -> Shape {
    brep.add_tsolid(vec![the_shell.clone()])
}

/// The UV-point pair of the pcurve-period adjustment (the OCCT
/// BRep_Tool::UVPoints(E, F, Pf, Pl) call of the wire reconstruction is
/// carried by the tool.rs brep_tool_uv_points).
pub(crate) fn _uv_points_anchor(
    the_e: &Shape,
    the_f: &Shape,
) -> (DVec2, DVec2) {
    brep_tool_uv_points(the_e, the_f)
}
