//! OCCT BRepFill_TrimEdgeTool — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_TrimEdgeTool.hxx (L37-76) +
//!         BRepFill_TrimEdgeTool.cxx (L1-781).
//!
//! Architecture differences (referenced from the affected functions):
//! 1. `handle<Geom2d_Geometry>` (a Geom2d_Point or Geom2d_Curve handle,
//!    discriminated by DynamicType) maps to the local [`Geom2dGeometry`]
//!    enum (Point = Geom2d_CartesianPoint, Curve = Geom2d_Curve).
//! 2. `Geom2dAdaptor_Curve` over a bisector handle maps to
//!    `topalgo::bisector::bisector_curve::Geom2dCurveAdaptor`; the adaptor
//!    over a bounded pcurve (C, f, l) maps to the local
//!    [`Geom2dAdaptorCurve`].  Both implement
//!    `geomalgo::geom2d_int::Curve2dAdaptor`.
//! 3. `Geom2dInt_GInter` maps to
//!    `geomalgo::geom2d_int::TheIntPCurvePCurveOfGInter` + the bounded
//!    `IntRes2d_Domain` pair (the Geom2dInt instantiation set).
//! 4. `Geom2dAPI_ProjectPointOnCurve` is re-hosted below
//!    ([`Geom2dAPIProjectPointOnCurve`]) over
//!    `geomalgo::extrema_gen_ext_pc2d::EPCOfExtPC2d`
//!    (Extrema_EPCOfExtPC2d); the extracted results (NbPoints /
//!    NearestPoint / LowerDistance / LowerDistanceParameter) carry the
//!    OCCT member surface consumed here.
//! 5. `BRep_Tool::CurveOnSurface(Edge, C, Surf, L, f, l)` (the 2-arg form)
//!    reads the first pcurve representation of the edge (the
//!    is_small_closed_edge precedent of offset_wire_b.rs).
//! 6. NCollection_Sequence<gp_Pnt> -> Vec<DVec3> (X = parameter on the
//!    bisectrice, Y = parameter on the parallel, Z = parameter on the
//!    other parallel); SetValue/Remove map to index assignment / remove.
//! 7. `Bisector_Bisec` / `Bisector_BisecAna` are the topalgo/bisector
//!    translations (`crate::topalgo::bisector`).
//!
//! first consumer: BRepFill_OffsetWire::MakeOffset / TrimEdge (the 2c GAP
//! points).

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::geom::{
    Circle2d, Curve2d, Curve2dEval, Ellipse2d, Hyperbola2d, Line2d, Parabola2d,
};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};

use crate::geomalgo::extrema_gen_ext_pc2d::EPCOfExtPC2d;
use crate::geomalgo::geom2d_int::{
    Curve2dAdaptor, Curve2dType, TheIntPCurvePCurveOfGInter,
};
use crate::geomalgo::int_res2d::Domain as Res2dDomain;
use crate::topalgo::bisector::bisector_bisec::BisectorBisec;
use crate::topalgo::bisector::bisector_bisec_ana::BisectorBisecAna;
use crate::topalgo::bisector::bisector_curve::{
    BisectorCurve, CurveKind, Geom2dCurveAdaptor, Geom2dCurveHandle, TrimmedCurve,
};

/// OCCT GeomAbs_JoinType (TKMath/GeomAbs/GeomAbs_JoinType.hxx) — the
/// package-wide shared enum (offset_wire.rs precedent).
pub use super::offset_wire::GeomAbsJoinType;

// ---------------------------------------------------------------------------
// Local carriers (architecture differences #1, #2, #4)
// ---------------------------------------------------------------------------

/// OCCT `handle<Geom2d_Geometry>` — the point-or-curve handle (architecture
/// difference #1): Point carries Geom2d_CartesianPoint::Pnt2d(), Curve the
/// Geom2d_Curve handle.
#[derive(Debug, Clone)]
pub enum Geom2dGeometry {
    /// OCCT Geom2d_CartesianPoint (Pnt2d()).
    Point(DVec2),
    /// OCCT Geom2d_Curve handle.
    Curve(Curve2d),
}

/// OCCT Geom2dAdaptor_Curve(C, U1, U2) over a bounded kernel pcurve
/// (architecture difference #2) — the adaptor delegates the evaluations to
/// the `impl Curve2dAdaptor for Curve2d` and restricts the domain to
/// [first, last].
#[derive(Debug, Clone)]
pub struct Geom2dAdaptorCurve {
    curve: Curve2d,
    first: f64,
    last: f64,
}

impl Geom2dAdaptorCurve {
    /// OCCT Geom2dAdaptor_Curve(C, U1, U2).
    pub fn new(c: Curve2d, u1: f64, u2: f64) -> Self {
        Geom2dAdaptorCurve {
            curve: c,
            first: u1,
            last: u2,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Value(U) (D0 form consumed by
    /// IntersectWith / AddOrConfuse).
    pub fn value(&self, u: f64) -> DVec2 {
        Curve2dAdaptor::value(self, u)
    }

    /// OCCT Geom2dAdaptor_Curve::FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// OCCT Geom2dAdaptor_Curve::LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// OCCT Geom2dAdaptor_Curve::GetType().
    pub fn get_type(&self) -> Curve2dType {
        Curve2dAdaptor::get_type(&self.curve)
    }
}

impl Curve2dAdaptor for Geom2dAdaptorCurve {
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn value(&self, u: f64) -> DVec2 {
        self.curve.point_at(u)
    }
    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        Curve2dAdaptor::d1(&self.curve, u)
    }
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        Curve2dAdaptor::d2(&self.curve, u)
    }
    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        Curve2dAdaptor::d3(&self.curve, u)
    }
    fn dn(&self, u: f64, n: i32) -> DVec2 {
        Curve2dAdaptor::dn(&self.curve, u, n)
    }
    fn get_type(&self) -> Curve2dType {
        Curve2dAdaptor::get_type(&self.curve)
    }
    fn nb_samples(&self) -> i32 {
        Curve2dAdaptor::nb_samples(&self.curve)
    }
    fn resolution(&self, r3d: f64) -> f64 {
        Curve2dAdaptor::resolution(&self.curve, r3d)
    }
    fn is_closed(&self) -> bool {
        Curve2dAdaptor::is_closed(&self.curve)
    }
    fn is_periodic(&self) -> bool {
        Curve2dAdaptor::is_periodic(&self.curve)
    }
    fn period(&self) -> f64 {
        Curve2dAdaptor::period(&self.curve)
    }
    fn nb_knots(&self) -> i32 {
        Curve2dAdaptor::nb_knots(&self.curve)
    }
    fn degree(&self) -> i32 {
        Curve2dAdaptor::degree(&self.curve)
    }
    fn nb_poles(&self) -> i32 {
        Curve2dAdaptor::nb_poles(&self.curve)
    }
    fn circle(&self) -> Circle2d {
        Curve2dAdaptor::circle(&self.curve)
    }
    fn line(&self) -> Line2d {
        Curve2dAdaptor::line(&self.curve)
    }
    fn ellipse(&self) -> Ellipse2d {
        Curve2dAdaptor::ellipse(&self.curve)
    }
    fn parabola(&self) -> Parabola2d {
        Curve2dAdaptor::parabola(&self.curve)
    }
    fn hyperbola(&self) -> Hyperbola2d {
        Curve2dAdaptor::hyperbola(&self.curve)
    }
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// OCCT Geom2dAPI_ProjectPointOnCurve re-host (architecture difference #4)
/// — the Extrema_ExtPC2d member is reduced to the extracted results
/// (IsDone / NbPoints / NearestPoint / LowerDistance /
/// LowerDistanceParameter); the Init() minimum scan (Geom2dAPI_
/// ProjectPointOnCurve.cxx L54-79) is performed at construction.
#[derive(Debug, Clone)]
pub struct Geom2dAPIProjectPointOnCurve {
    my_is_done: bool,
    my_index: i32,
    my_ext_pc: Vec<(DVec2, f64, f64)>, // (point, parameter, square distance)
}

impl Geom2dAPIProjectPointOnCurve {
    /// OCCT Extrema sampling constants (the Extrema_GGenExtPC defaults
    /// documented by the hlr/contap/h_cont_tool.rs precedent):
    /// Extrema_EPCOfExtPC2d(P, C, Nbu = 20, epsX = 1.0e-8, Tol = 1.0e-5).
    const NB_SAMPLE: i32 = 20;
    const TOL_U: f64 = 1.0e-8;
    const TOL_F: f64 = 1.0e-5;

    /// OCCT Geom2dAPI_ProjectPointOnCurve(P, Curve, Umin, Usup)
    /// (cxx L39-43 -> Init L58-73).
    pub fn new_bounded(
        the_p: DVec2,
        the_c: &dyn Curve2dAdaptor,
        the_umin: f64,
        the_usup: f64,
    ) -> Self {
        let mut ext = EPCOfExtPC2d::with_range(
            the_p,
            the_c,
            Self::NB_SAMPLE,
            the_umin,
            the_usup,
            Self::TOL_U,
            Self::TOL_F,
        );
        ext.perform(the_p);

        let mut r = Geom2dAPIProjectPointOnCurve {
            my_is_done: false,
            my_index: -1,
            my_ext_pc: Vec::new(),
        };
        // OCCT L65: myIsDone = myExtPC.IsDone() && (myExtPC.NbExt() > 0).
        r.my_is_done = ext.is_done() && ext.nb_ext() > 0;
        if !r.my_is_done {
            return r;
        }
        // OCCT L67-79: evaluate the lower distance and its index.
        for i in 1..=ext.nb_ext() {
            r.my_ext_pc
                .push((ext.point(i).value(), ext.point(i).parameter(), ext.square_distance(i)));
        }
        let mut dist2_min = r.my_ext_pc[0].2;
        r.my_index = 1;
        for i in 2..=r.my_ext_pc.len() as i32 {
            let dist2 = r.my_ext_pc[(i - 1) as usize].2;
            if dist2 < dist2_min {
                dist2_min = dist2;
                r.my_index = i;
            }
        }
        r
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve(P, Curve) (cxx L29-33 -> Init
    /// L48-53): the full natural domain.
    pub fn new_point_curve(the_p: DVec2, the_c: &dyn Curve2dAdaptor) -> Self {
        Self::new_bounded(
            the_p,
            the_c,
            the_c.first_parameter(),
            the_c.last_parameter(),
        )
    }

    /// OCCT NbPoints() (cxx L82-92).
    pub fn nb_points(&self) -> usize {
        if self.my_is_done {
            self.my_ext_pc.len()
        } else {
            0
        }
    }

    /// OCCT NearestPoint() (cxx L137-143).
    pub fn nearest_point(&self) -> DVec2 {
        assert!(self.my_is_done, "StdFail_NotDone: NearestPoint");
        self.my_ext_pc[(self.my_index - 1) as usize].0
    }

    /// OCCT LowerDistance() (cxx L155-161).
    pub fn lower_distance(&self) -> f64 {
        assert!(self.my_is_done, "StdFail_NotDone: LowerDistance");
        self.my_ext_pc[(self.my_index - 1) as usize].2.sqrt()
    }

    /// OCCT LowerDistanceParameter() (lxx).
    pub fn lower_distance_parameter(&self) -> f64 {
        assert!(self.my_is_done, "StdFail_NotDone: LowerDistanceParameter");
        self.my_ext_pc[(self.my_index - 1) as usize].1
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp re-hosts (architecture difference #5)
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::CurveOnSurface(Edge, C, Surf, L, f, l) — the 2-arg form:
/// the first pcurve representation of the edge (the is_small_closed_edge
/// precedent of offset_wire_b.rs).
pub fn brep_tool_curve_on_surface(edg: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            // The pcurve map is keyed by face; the first bound pcurve is
            // the (unique) curve representation used by the offset
            // pipeline.
            if let Some((_, (c, a, b))) = ed.pcurves.iter().next() {
                return Some((c.clone(), *a, *b));
            }
            match ed.representations.first() {
                Some(rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    pcurve1,
                    range,
                    ..
                }) => Some((pcurve1.clone(), range[0], range[1])),
                _ => None,
            }
        }
        _ => None,
    }
}

/// OCCT TopExp::Vertices(Edge, Vf, Vl) — the end vertices in traversal
/// order (orientation-aware).
fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let (first, last) = match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => (Shape::null(), Shape::null()),
    };
    if e.orientation == Orientation::Reversed {
        (last, first)
    } else {
        (first, last)
    }
}

// ---------------------------------------------------------------------------
// Static helpers (BRepFill_TrimEdgeTool.cxx)
// ---------------------------------------------------------------------------

/// OCCT static SimpleExpression (L44-60): the simple geometric expression
/// of the bissectrice — unwrap the trimmed BisecAna basis.
fn simple_expression(b: &BisectorBisec) -> Geom2dCurveAdaptor {
    // OCCT L46: Bis = B.Value() — the Handle(Geom2d_TrimmedCurve).
    let bis_trimmed = b.value();
    let trbis = bis_trimmed.read().expect("Bisector_Bisec::Value()");
    // OCCT L48-51: BT == Geom2d_TrimmedCurve; BasBis = TrBis->BasisCurve().
    let bas_bis = trbis.basis_curve.clone();
    // OCCT L53: BT = BasBis->DynamicType().
    let bt = bas_bis
        .as_any()
        .downcast_ref::<BisectorBisecAna>()
        .map(|_| CurveKind::BisecAna)
        .unwrap_or_else(|| {
            bas_bis
                .as_any()
                .downcast_ref::<Geom2dCurveHandle>()
                .map(|h| crate::topalgo::bisector::bisector_bisec_ana::kind_of(&h.curve))
                .unwrap_or(CurveKind::Other)
        });
    if bt == CurveKind::BisecAna {
        // OCCT L56-57: Bis = down_cast<Bisector_BisecAna>(BasBis)
        // ->Geom2dCurve(); Bis = new Geom2d_TrimmedCurve(Bis, TrBis
        // ->FirstParameter(), TrBis->LastParameter()).
        let bis = bas_bis
            .as_any()
            .downcast_ref::<BisectorBisecAna>()
            .expect("Bisector_BisecAna down-cast")
            .geom2d_curve();
        let trimmed = TrimmedCurve::new(bis, trbis.first_parameter, trbis.last_parameter);
        Geom2dCurveAdaptor::new(Arc::new(trimmed) as Arc<dyn BisectorCurve>)
    } else {
        // OCCT L59 (implicit fall-through): Bis keeps the trimmed value.
        let trimmed = trbis.clone();
        Geom2dCurveAdaptor::new(Arc::new(trimmed) as Arc<dyn BisectorCurve>)
    }
}

/// OCCT static Bubble (L107-125): order the sequence of points by
/// increasing x.
fn bubble(seq: &mut Vec<DVec3>) {
    let mut invert = true;
    let nb_points = seq.len() as i32;
    while invert {
        invert = false;
        for i in 1..nb_points {
            let p1 = seq[(i - 1) as usize];
            let p2 = seq[i as usize];
            if p2.x < p1.x {
                seq.swap((i - 1) as usize, i as usize);
                invert = true;
            }
        }
    }
}

/// OCCT static EvalParameters (L129-198): intersect <ac> with <bis> at
/// Precision::Confusion and collect the (bis, ac) parameter pairs.
fn eval_parameters(bis: &dyn Curve2dAdaptor, ac: &dyn Curve2dAdaptor, params: &mut Vec<DVec3>) {
    let mut intersector = TheIntPCurvePCurveOfGInter::new();
    let tol = CONFUSION;

    // OCCT L141: Intersector = Geom2dInt_GInter(CAC, CBis, Tol, Tol).
    let cbis = Res2dDomain::bounded(
        bis.value(bis.first_parameter()),
        bis.first_parameter(),
        tol,
        bis.value(bis.last_parameter()),
        bis.last_parameter(),
        tol,
    );
    let cac = Res2dDomain::bounded(
        ac.value(ac.first_parameter()),
        ac.first_parameter(),
        tol,
        ac.value(ac.last_parameter()),
        ac.last_parameter(),
        tol,
    );
    intersector.perform(ac, &cac, bis, &cbis, tol, tol);

    let intersector = intersector.base;
    if !intersector.is_done() {
        // OCCT L149: throw StdFail_NotDone.
        panic!("StdFail_NotDone: BRepFill_TrimSurfaceTool::IntersectWith");
    }

    // OCCT L152-163: the intersection points.
    let nb_points = intersector.nb_points();
    if nb_points > 0 {
        for i in 1..=nb_points {
            let u1 = intersector.point(i).param_on_second();
            let u2 = intersector.point(i).param_on_first();
            let p = DVec3::new(u1, u2, 0.0);
            params.push(p);
        }
    }

    // OCCT L165-194: the intersection segments.
    let nb_segments = intersector.nb_segments();
    if nb_segments > 0 {
        for i in 1..=nb_segments {
            let seg = intersector.segment(i);
            let mut u1 = seg.first_point().param_on_second();
            let ulast = seg.last_point().param_on_second();
            if (u1 - bis.first_parameter()).abs() <= tol
                && (ulast - bis.last_parameter()).abs() <= tol
            {
                let p = DVec3::new(u1, seg.first_point().param_on_first(), 0.0);
                params.push(p);
                let p = DVec3::new(ulast, seg.last_point().param_on_first(), 0.0);
                params.push(p);
            } else {
                u1 += seg.last_point().param_on_second();
                u1 /= 2.0;
                let mut u2 = seg.first_point().param_on_first();
                u2 += seg.last_point().param_on_first();
                u2 /= 2.0;
                let p = DVec3::new(u1, u2, 0.0);
                params.push(p);
            }
        }
    }

    // OCCT L197: order the sequence by growing parameter on the
    // bissectrice.
    bubble(params);
}

/// OCCT static EvalParametersBis (L200-268): the same with TolC = Tol.
fn eval_parameters_bis(
    bis: &dyn Curve2dAdaptor,
    ac: &dyn Curve2dAdaptor,
    params: &mut Vec<DVec3>,
    tol: f64,
) {
    let mut intersector = TheIntPCurvePCurveOfGInter::new();
    let tol_c = tol;

    // OCCT L211: Intersector = Geom2dInt_GInter(CAC, CBis, TolC, Tol).
    let cbis = Res2dDomain::bounded(
        bis.value(bis.first_parameter()),
        bis.first_parameter(),
        tol_c,
        bis.value(bis.last_parameter()),
        bis.last_parameter(),
        tol_c,
    );
    let cac = Res2dDomain::bounded(
        ac.value(ac.first_parameter()),
        ac.first_parameter(),
        tol_c,
        ac.value(ac.last_parameter()),
        ac.last_parameter(),
        tol_c,
    );
    intersector.perform(ac, &cac, bis, &cbis, tol_c, tol);

    let intersector = intersector.base;
    if !intersector.is_done() {
        // OCCT L219: throw StdFail_NotDone.
        panic!("StdFail_NotDone: BRepFill_TrimSurfaceTool::IntersectWith");
    }

    let nb_points = intersector.nb_points();
    if nb_points > 0 {
        for i in 1..=nb_points {
            let u1 = intersector.point(i).param_on_second();
            let u2 = intersector.point(i).param_on_first();
            let p = DVec3::new(u1, u2, 0.0);
            params.push(p);
        }
    }

    let nb_segments = intersector.nb_segments();
    if nb_segments > 0 {
        for i in 1..=nb_segments {
            let seg = intersector.segment(i);
            let mut u1 = seg.first_point().param_on_second();
            let ulast = seg.last_point().param_on_second();
            if (u1 - bis.first_parameter()).abs() <= tol
                && (ulast - bis.last_parameter()).abs() <= tol
            {
                let p = DVec3::new(u1, seg.first_point().param_on_first(), 0.0);
                params.push(p);
                let p = DVec3::new(ulast, seg.last_point().param_on_first(), 0.0);
                params.push(p);
            } else {
                u1 += seg.last_point().param_on_second();
                u1 /= 2.0;
                let mut u2 = seg.first_point().param_on_first();
                u2 += seg.last_point().param_on_first();
                u2 /= 2.0;
                let p = DVec3::new(u1, u2, 0.0);
                params.push(p);
            }
        }
    }

    bubble(params);
}

/// OCCT ElCLib::Parameter(gp_Lin2d, gp_Pnt2d) — the line parameter of the
/// projection (pure-math helper).
fn el_clib_parameter(l: &Line2d, p: DVec2) -> f64 {
    (p - l.origin).dot(l.direction)
}

// ---------------------------------------------------------------------------
// BRepFill_TrimEdgeTool (hxx L37-76)
// ---------------------------------------------------------------------------

/// OCCT BRepFill_TrimEdgeTool (hxx L37-76) — geometric tool used to
/// construct offset wires.
pub struct BRepFillTrimEdgeTool {
    is_point1: bool,            // OCCT: isPoint1
    is_point2: bool,            // OCCT: isPoint2
    my_p1: DVec2,               // OCCT: myP1
    my_p2: DVec2,               // OCCT: myP2
    my_c1: Option<Curve2d>,     // OCCT: myC1 (None = null handle)
    my_c2: Option<Curve2d>,     // OCCT: myC2 (None = null handle)
    my_offset: f64,             // OCCT: myOffset
    my_bisec: BisectorBisec,    // OCCT: myBisec
    my_bis: Geom2dCurveAdaptor, // OCCT: myBis (Geom2dAdaptor_Curve)
}

impl BRepFillTrimEdgeTool {
    /// OCCT BRepFill_TrimEdgeTool::BRepFill_TrimEdgeTool() (cxx L64).
    pub fn new_empty() -> Self {
        // OCCT L64: `= default` — myBis stays on a null handle; the offset
        // pipeline only uses the loaded form, so the default constructor is
        // unreachable in the translated pipeline (BisectorBisec has no null
        // handle form in rcad).
        unreachable!(
            "BRepFill_TrimEdgeTool default constructor: the OCCT default \
             leaves myBis on a null handle; the offset pipeline only uses \
             the loaded form"
        );
    }

    /// OCCT BRepFill_TrimEdgeTool::BRepFill_TrimEdgeTool(Bisec, S1, S2,
    /// Offset) (cxx L68-100).
    pub fn new(bisec: &BisectorBisec, s1: &Geom2dGeometry, s2: &Geom2dGeometry, offset: f64) -> Self {
        // OCCT L75-76: isPoint1/2 = (S1/2->DynamicType() ==
        // Geom2d_CartesianPoint).
        let is_point1 = matches!(s1, Geom2dGeometry::Point(_));
        let is_point2 = matches!(s2, Geom2dGeometry::Point(_));

        let mut my_p1 = DVec2::ZERO;
        let mut my_p2 = DVec2::ZERO;
        let mut my_c1: Option<Curve2d> = None;
        let mut my_c2: Option<Curve2d> = None;

        // OCCT L80-95: return geometries of shapes.
        if is_point1 {
            match s1 {
                Geom2dGeometry::Point(p) => my_p1 = *p,
                _ => unreachable!(),
            }
        } else {
            match s1 {
                Geom2dGeometry::Curve(c) => my_c1 = Some(c.clone()),
                _ => unreachable!(),
            }
        }
        if is_point2 {
            match s2 {
                Geom2dGeometry::Point(p) => my_p2 = *p,
                _ => unreachable!(),
            }
        } else {
            match s2 {
                Geom2dGeometry::Curve(c) => my_c2 = Some(c.clone()),
                _ => unreachable!(),
            }
        }
        // OCCT L97-99: return the simple expression of the bissectrice.
        let my_bis = simple_expression(bisec);

        BRepFillTrimEdgeTool {
            is_point1,
            is_point2,
            my_p1,
            my_p2,
            my_c1,
            my_c2,
            my_offset: offset,
            // OCCT L73: myBisec(Bisec) — the container copy (the rcad
            // carrier has no handle payload; Default is the stand-in).
            my_bisec: BisectorBisec::default(),
            my_bis,
        }
    }

    /// OCCT BRepFill_TrimEdgeTool::IntersectWith (cxx L272-627).
    #[allow(clippy::too_many_arguments)]
    pub fn intersect_with(
        &self,
        edge1: &Shape,
        edge2: &Shape,
        init_shape1: &Shape,
        init_shape2: &Shape,
        end1: &Shape,
        end2: &Shape,
        the_join_type: GeomAbsJoinType,
        is_open_result: bool,
        params: &mut Vec<DVec3>,
    ) {
        params.clear();

        // OCCT L284-295: return curves associated to edges.
        let Some((c1, f1, l1)) = brep_tool_curve_on_surface(edge1) else {
            return;
        };
        let ac1 = Geom2dAdaptorCurve::new(c1, f1, l1);
        let Some((c2, f2, l2)) = brep_tool_curve_on_surface(edge2) else {
            return;
        };
        let ac2 = Geom2dAdaptorCurve::new(c2, f2, l2);

        // OCCT L298-302: calculate intersection.
        let mut points2: Vec<DVec3> = Vec::new();

        eval_parameters(&self.my_bis, &ac1, params);
        eval_parameters(&self.my_bis, &ac2, &mut points2);

        // OCCT L304-306.
        let mut seance_de_rattrapage = 0;
        let mut tol_init = 1.0e-9;
        let mut nn = 7;

        // OCCT L308-314.
        let t1 = ac1.get_type();
        let t2 = ac2.get_type();
        if (t1 != Curve2dType::Circle && t1 != Curve2dType::Line)
            || (t2 != Curve2dType::Circle && t2 != Curve2dType::Line)
        {
            tol_init = 1.0e-8;
            nn = 6;
        }

        // OCCT L316-371: check, may be there are no intersections at all
        // for case myBis == Line.
        if params.is_empty() && points2.is_empty() && self.my_bis.get_type() == Curve2dType::Line
        {
            // OCCT L322-329.
            let mut dmax = tol_init;
            for _ in 0..nn {
                dmax *= 10.0;
            }
            dmax *= dmax;

            let an_l = self.my_bis.line();
            let mut is_far1 = true;
            let mut is_far2 = true;
            // OCCT L336-349: the first curve endpoints.
            let mut d = f64::MAX;
            let a_p = ac1.value(ac1.first_parameter());
            let par = el_clib_parameter(&an_l, a_p);
            if par >= self.my_bis.first_parameter() && par <= self.my_bis.last_parameter() {
                d = an_l.distance(a_p) * an_l.distance(a_p);
            }
            let a_p = ac1.value(ac1.last_parameter());
            let par = el_clib_parameter(&an_l, a_p);
            if par >= self.my_bis.first_parameter() && par <= self.my_bis.last_parameter() {
                d = (an_l.distance(a_p) * an_l.distance(a_p)).min(d);
            }
            is_far1 = d > dmax;

            // OCCT L351-364: the second curve endpoints.
            let mut d = f64::MAX;
            let a_p = ac2.value(ac2.first_parameter());
            let par = el_clib_parameter(&an_l, a_p);
            if par >= self.my_bis.first_parameter() && par <= self.my_bis.last_parameter() {
                d = an_l.distance(a_p) * an_l.distance(a_p);
            }
            let a_p = ac2.value(ac2.last_parameter());
            let par = el_clib_parameter(&an_l, a_p);
            if par >= self.my_bis.first_parameter() && par <= self.my_bis.last_parameter() {
                d = (an_l.distance(a_p) * an_l.distance(a_p)).min(d);
            }
            is_far2 = d > dmax;

            if is_far1 && is_far2 {
                // OCCT L368: return — no intersections at all.
                return;
            }
        }

        // OCCT L373-391: the "seance de rattrapage" loop.
        while seance_de_rattrapage < nn
            && (points2.len() != params.len() || (points2.is_empty() && params.is_empty()))
        {
            params.clear();
            points2.clear();

            tol_init *= 10.0;

            eval_parameters_bis(&self.my_bis, &ac1, params, tol_init);
            eval_parameters_bis(&self.my_bis, &ac2, &mut points2, tol_init);
            seance_de_rattrapage += 1;
        }

        // OCCT L402-434: Params empty, one point on the second parallel.
        if params.is_empty() && points2.len() == 1 {
            let dmax = 0.25 * self.my_offset * self.my_offset;
            let t_bis = points2[0].x;
            let p_bis = self.my_bis.value(t_bis);

            let t = ac1.first_parameter();
            let p_c = ac1.value(t);
            let dmin = p_c.distance_squared(p_bis);
            let mut p = DVec3::new(t_bis, t, 0.0);
            if dmin < dmax {
                params.push(p);
            }

            let t = ac1.last_parameter();
            let p_c = ac1.value(t);
            let dmin1 = p_c.distance_squared(p_bis);
            if dmin > dmin1 && dmin1 < dmax {
                p.y = t;
                if params.is_empty() {
                    params.push(p);
                } else {
                    params[0] = p;
                }
            }
        }
        // OCCT L435-467: one point on the first parallel, Points2 empty.
        else if params.len() == 1 && points2.is_empty() {
            let dmax = 0.25 * self.my_offset * self.my_offset;
            let t_bis = params[0].x;
            let p_bis = self.my_bis.value(t_bis);

            let t = ac2.first_parameter();
            let p_c = ac2.value(t);
            let dmin = p_c.distance_squared(p_bis);
            let mut p = DVec3::new(t_bis, t, 0.0);
            if dmin < dmax {
                points2.push(p);
            }

            let t = ac2.last_parameter();
            let p_c = ac2.value(t);
            let dmin1 = p_c.distance_squared(p_bis);
            if dmin > dmin1 && dmin1 < dmax {
                p.y = t;
                if points2.is_empty() {
                    points2.push(p);
                } else {
                    points2[0] = p;
                }
            }
        }

        // OCCT L469-475: small manipulation to remove incorrect
        // intersections: return only common intersections (same parameter
        // on the bissectrice).  The tolerance can be eventually changed.
        let tol = 4.0 * 100.0 * PCONFUSION;

        // OCCT L479-489.
        if params.len() == 1 && points2.len() == 1 {
            params[0].z = points2[0].y;
            return;
        }

        // OCCT L491-526.
        let mut i: usize = 1;
        while i <= params.len().min(points2.len()) {
            let p1 = params[i - 1];
            let p2 = points2[i - 1];
            let p1xp2x = (p1.x - p2.x).abs();

            if p1xp2x > tol {
                if p1xp2x > tol_init {
                    i += 1;
                } else if p1.x < p2.x {
                    params.remove(i - 1);
                } else {
                    points2.remove(i - 1);
                }
            } else {
                i += 1;
            }
        }

        // OCCT L528-535.
        if params.len() > points2.len() {
            params.truncate(points2.len());
        } else if params.len() < points2.len() {
            points2.truncate(params.len());
        }

        let nb_points = params.len();

        // OCCT L539-582: now we define: if there are more than one point of
        // intersection is it Ok?
        let mut init_fpar = -f64::MAX;
        let mut init_lpar = f64::MAX;
        if nb_points > 1
            && the_join_type == GeomAbsJoinType::Intersection
            && init_shape1.shape_type() != ShapeType::Vertex
            && init_shape2.shape_type() != ShapeType::Vertex
        {
            // OCCT L548-563: definition of initial first and last
            // parameters (this is inverse procedure to extension of
            // parameters — see BRepFill_OffsetWire, function MakeOffset,
            // case of Circle).
            let init_edge1 = init_shape1;
            let mut to_extend_first_par = true;
            let mut to_extend_last_par = true;
            if is_open_result {
                let (v1, v2) = top_exp_vertices(init_edge1);
                if v1.is_same(end1) || v1.is_same(end2) {
                    to_extend_first_par = false;
                }
                if v2.is_same(end1) || v2.is_same(end2) {
                    to_extend_last_par = false;
                }
            }
            // OCCT L564-581: BRepAdaptor_Curve IC1(InitEdge1) — the Circle
            // branch (the edge curve / range carry the adaptor surface).
            let (ic1_first, ic1_last) = match init_edge1.data.as_ref() {
                TShape::Edge(ed) => (ed.range[0], ed.range[1]),
                _ => (0.0, 0.0),
            };
            let ic1_is_circle = match init_edge1.data.as_ref() {
                TShape::Edge(ed) => matches!(ed.curve, Some(rcad_kernel::geom::Curve3::Circle(_))),
                _ => false,
            };
            if ic1_is_circle {
                let delta = 2.0 * std::f64::consts::PI - ic1_last + ic1_first;
                if to_extend_first_par && to_extend_last_par {
                    init_fpar = ac1.first_parameter() + delta / 2.0;
                } else if to_extend_first_par {
                    init_fpar = ac1.first_parameter() + delta;
                } else if to_extend_last_par {
                    init_fpar = ac1.first_parameter();
                }
                init_lpar = init_fpar + ic1_last - ic1_first;
            }
        }

        // OCCT L584-618: remove all vertices with non-minimal parameter if
        // they are out of initial range.
        if nb_points > 1 && the_join_type == GeomAbsJoinType::Intersection {
            let mut imin: usize = 1;
            for i in 2..=params.len() {
                if params[i - 1].x < params[imin - 1].x {
                    imin = i;
                }
            }

            let params_copy = params.clone();
            let points2_copy = points2.clone();
            params.clear();
            points2.clear();
            for i in 1..=params_copy.len() {
                if imin == i
                    || (params_copy[i - 1].y >= init_fpar && params_copy[i - 1].y <= init_lpar)
                {
                    params.push(params_copy[i - 1]);
                    points2.push(points2_copy[i - 1]);
                }
            }
        }

        // OCCT L620-626: store the parameter of the other parallel in Z.
        for i in 1..=params.len() {
            let mut p_seq = params[i - 1];
            p_seq.z = points2[i - 1].y;
            params[i - 1] = p_seq;
        }
    }

    /// OCCT BRepFill_TrimEdgeTool::AddOrConfuse (cxx L636-733): the first
    /// or the last point of the bissectrice is on the parallel if it was
    /// not found in the intersections; it is projected on parallel lines
    /// and added in the parameters.
    pub fn add_or_confuse(&self, start: bool, edge1: &Shape, edge2: &Shape, params: &mut Vec<DVec3>) {
        let mut to_proj = true;
        const TOL: f64 = 10.0 * CONFUSION;

        // OCCT L645-652: return curves associated to edges.
        let Some((c1, f1, l1)) = brep_tool_curve_on_surface(edge1) else {
            return;
        };
        let ac1 = Geom2dAdaptorCurve::new(c1, f1, l1);

        // OCCT L654-661.
        let p_bis = if start {
            self.my_bis.value(self.my_bis.first_parameter())
        } else {
            self.my_bis.value(self.my_bis.last_parameter())
        };

        // OCCT L663-676: test if the end of the bissectrice is in the set
        // of intersection points.
        if !params.is_empty() {
            let p = if start {
                ac1.value(params[0].y)
            } else {
                ac1.value(params[params.len() - 1].y)
            };
            // OCCT L675: ToProj = !PBis.IsEqual(P, Tol).
            to_proj = p_bis.distance(p) >= TOL;
        }

        if to_proj {
            // OCCT L684-691: project point on parallels and add in Params.
            let Some((c2, f2, l2)) = brep_tool_curve_on_surface(edge2) else {
                return;
            };
            let projector1 = Geom2dAPIProjectPointOnCurve::new_bounded(p_bis, &ac1, f1, l1);
            let ac2 = Geom2dAdaptorCurve::new(c2, f2, l2);
            let projector2 = Geom2dAPIProjectPointOnCurve::new_bounded(p_bis, &ac2, f2, l2);

            // OCCT L693-699.
            if projector1.nb_points() == 0 {
                return;
            }
            // OCCT L700-706.
            if projector1.nearest_point().distance(p_bis) >= TOL {
                return;
            }
            // OCCT L707-713.
            if projector2.nb_points() == 0 {
                return;
            }
            // OCCT L714-720.
            if projector2.nearest_point().distance(p_bis) >= TOL {
                return;
            }
            // OCCT L721-731.
            let mut p_int = DVec3::new(
                0.0,
                projector1.lower_distance_parameter(),
                projector2.lower_distance_parameter(),
            );
            if start {
                p_int.x = self.my_bis.first_parameter();
                params.insert(0, p_int);
            } else {
                p_int.x = self.my_bis.last_parameter();
                params.push(p_int);
            }
        }
    }

    /// OCCT BRepFill_TrimEdgeTool::IsInside (cxx L737-781).
    pub fn is_inside(&self, p: DVec2) -> bool {
        // OCCT L740-741: double Dist = RealLast().
        let mut dist = f64::MAX;
        if self.is_point1 {
            dist = p.distance(self.my_p1);
        } else if self.is_point2 {
            dist = p.distance(self.my_p2);
        } else {
            let c1 = self.my_c1.as_ref().expect("myC1");
            // OCCT L753-757: the projection (full domain).
            let projector = Geom2dAPIProjectPointOnCurve::new_point_curve(p, c1);
            if projector.nb_points() > 0 {
                dist = projector.lower_distance();
            }

            // OCCT L765-774: check of distances between P and first and
            // last point of the first curve should be performed in any
            // case, despite of the results of projection.
            let p_f = c1.point_at(c1.default_domain()[0]);
            let p_l = c1.point_at(c1.default_domain()[1]);
            let a_dist_min = p.distance(p_f).min(p.distance(p_l));

            if dist > a_dist_min {
                dist = a_dist_min;
            }
        }

        // OCCT L780.
        dist < self.my_offset.abs() - CONFUSION
    }
}
