//! OCCT Bisector_Curve (base class of the Bisector package) plus the
//! Geom2d_Curve plumbing the package needs, 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_Curve.hxx (L28-52) / .cxx (L21).
//!
//! Architecture differences (Rust has no inheritance, and the kernel
//! `Curve2d` is a closed enum that cannot host Bisector variants):
//! - OCCT `class Bisector_Curve : public Geom2d_Curve` maps to the single
//!   merged trait [`BisectorCurve`] below: the Geom2d_Curve virtuals
//!   (Value/D0/D1/D2/D3/DN/FirstParameter/.../Copy) and the Bisector_Curve
//!   pure virtuals (Parameter/IsExtendAtStart/IsExtendAtEnd/NbIntervals/
//!   IntervalFirst/IntervalLast) live on one trait.
//! - A plain analytic `Geom2d_Curve` (Line/Circle/BSpline/... from the
//!   kernel) is represented by [`Geom2dCurveHandle`]; the Bisector-only pure
//!   virtuals keep the OCCT pure-virtual failure path (`unimplemented!`).
//! - OCCT `Geom2d_TrimmedCurve` over a Bisector_Curve cannot be expressed by
//!   the kernel enum, so it lives here as [`TrimmedCurve`] (same semantics:
//!   raw parameter bounds + delegation to the basis curve).

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Curve2d, Curve2dEval, Ellipse2d, Line2d, Parabola2d, Hyperbola2d};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::geom2d_int::{geom2d_curve_tool, Curve2dAdaptor, Curve2dType};

/// OCCT gp::Resolution() (gp.hxx L102).
pub const GP_RESOLUTION: f64 = 1e-15;

/// OCCT Geom2d_Curve::ResD1 (Geom2d_Curve.hxx) — point + first derivative.
#[derive(Debug, Clone, Copy)]
pub struct ResD1 {
    pub point: DVec2,
    pub d1: DVec2,
}

/// OCCT Geom2d_Curve::ResD2 — point + first/second derivative.
#[derive(Debug, Clone, Copy)]
pub struct ResD2 {
    pub point: DVec2,
    pub d1: DVec2,
    pub d2: DVec2,
}

/// OCCT Geom2d_Curve::ResD3 — point + first/second/third derivative.
#[derive(Debug, Clone, Copy)]
pub struct ResD3 {
    pub point: DVec2,
    pub d1: DVec2,
    pub d2: DVec2,
    pub d3: DVec2,
}

/// Dynamic kind of a curve handle — the OCCT `DynamicType()` /
/// `STANDARD_TYPE(...)` comparisons used inside the package.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    /// OCCT Bisector_BisecAna.
    BisecAna,
    /// OCCT Bisector_BisecCC.
    BisecCC,
    /// OCCT Bisector_BisecPC.
    BisecPC,
    /// OCCT Geom2d_Line.
    Line,
    /// OCCT Geom2d_Circle.
    Circle,
    /// OCCT Geom2d_Ellipse.
    Ellipse,
    /// OCCT Geom2d_Parabola.
    Parabola,
    /// OCCT Geom2d_Hyperbola.
    Hyperbola,
    /// Any other Geom2d_Curve type.
    Other,
}

/// OCCT `Bisector_Curve : Geom2d_Curve` — the merged interface trait
/// (architecture note above).  Handles are `Arc<dyn BisectorCurve>`.
pub trait BisectorCurve: Send + Sync {
    // ------------------------------------------------------------------
    // OCCT Geom2d_Curve interface (the part Bisector uses).
    // ------------------------------------------------------------------

    /// OCCT Value(U) — the point of parameter U (Geom2d_Curve::Value == D0).
    fn value(&self, u: f64) -> DVec2;

    /// OCCT D0(U, P).
    fn d0(&self, u: f64) -> DVec2 {
        self.value(u)
    }

    /// OCCT D1(U, P, V1).
    fn d1(&self, u: f64) -> (DVec2, DVec2);

    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2);

    /// OCCT D3(U, P, V1, V2, V3).
    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2);

    /// OCCT DN(U, N).
    fn dn(&self, u: f64, n: i32) -> DVec2;

    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64;

    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64;

    /// OCCT IsClosed().
    fn is_closed(&self) -> bool;

    /// OCCT IsPeriodic().
    fn is_periodic(&self) -> bool;

    /// OCCT Period() — only valid on periodic curves (conics: 2*pi).
    fn period(&self) -> f64 {
        // OCCT Geom2d_Conic::Period() = 2*M_PI.
        std::f64::consts::TAU
    }

    /// OCCT Continuity().
    fn continuity(&self) -> GeomAbsShape;

    /// OCCT ReversedParameter(U).
    fn reversed_parameter(&self, u: f64) -> f64;

    /// OCCT Reverse().
    fn reverse(&mut self);

    /// OCCT Reversed() — a reversed copy under a new handle.
    fn reversed(&self) -> Arc<dyn BisectorCurve> {
        let mut c = self.copy_curve();
        BisectorCurve::reverse(Arc::get_mut(&mut c).expect("exclusive copy"));
        c
    }

    /// OCCT IsCN(N).
    fn is_cn(&self, n: i32) -> bool;

    /// OCCT Transform(T) — gp_Trsf2d maps to `glam::DAffine2`.
    fn transform(&mut self, t: &glam::DAffine2);

    /// OCCT Copy() — deep copy as a new handle.
    fn copy_curve(&self) -> Arc<dyn BisectorCurve>;

    // ------------------------------------------------------------------
    // OCCT 8.0 Eval* surface (Geom2d_Curve::EvalD0/EvalD1/EvalD2/EvalD3/
    // EvalDN) — defaults delegate to the classic D-interfaces.
    // ------------------------------------------------------------------

    /// OCCT EvalD0(U).
    fn eval_d0(&self, u: f64) -> DVec2 {
        self.d0(u)
    }

    /// OCCT EvalD1(U).
    fn eval_d1(&self, u: f64) -> ResD1 {
        let (point, d1) = self.d1(u);
        ResD1 { point, d1 }
    }

    /// OCCT EvalD2(U).
    fn eval_d2(&self, u: f64) -> ResD2 {
        let (point, d1, d2) = self.d2(u);
        ResD2 { point, d1, d2 }
    }

    /// OCCT EvalD3(U).
    fn eval_d3(&self, u: f64) -> ResD3 {
        let (point, d1, d2, d3) = self.d3(u);
        ResD3 { point, d1, d2, d3 }
    }

    /// OCCT EvalDN(U, N).
    fn eval_dn(&self, u: f64, n: i32) -> DVec2 {
        self.dn(u, n)
    }

    // ------------------------------------------------------------------
    // OCCT Bisector_Curve pure virtuals (Bisector_Curve.hxx L32-49).
    // ------------------------------------------------------------------

    /// OCCT Parameter(const gp_Pnt2d& P) const — pure virtual.
    fn parameter(&self, _p: DVec2) -> f64 {
        // OCCT pure virtual — only reachable for a plain Geom2d_Curve
        // (e.g. Geom2d_TrimmedCurve); OCCT would fail the pure call.
        unimplemented!("Bisector_Curve::Parameter pure virtual")
    }

    /// OCCT IsExtendAtStart() const — pure virtual.
    fn is_extend_at_start(&self) -> bool {
        unimplemented!("Bisector_Curve::IsExtendAtStart pure virtual")
    }

    /// OCCT IsExtendAtEnd() const — pure virtual.
    fn is_extend_at_end(&self) -> bool {
        unimplemented!("Bisector_Curve::IsExtendAtEnd pure virtual")
    }

    /// OCCT NbIntervals() const — pure virtual.
    fn nb_intervals(&self) -> i32 {
        unimplemented!("Bisector_Curve::NbIntervals pure virtual")
    }

    /// OCCT IntervalFirst(const int Index) const — pure virtual.
    fn interval_first(&self, _index: i32) -> f64 {
        unimplemented!("Bisector_Curve::IntervalFirst pure virtual")
    }

    /// OCCT IntervalLast(const int Index) const — pure virtual.
    fn interval_last(&self, _index: i32) -> f64 {
        unimplemented!("Bisector_Curve::IntervalLast pure virtual")
    }

    // ------------------------------------------------------------------
    // Dynamic type mirror.
    // ------------------------------------------------------------------

    /// OCCT DynamicType() — the kind comparisons used inside the package.
    fn kind(&self) -> CurveKind;

    /// OCCT `occ::down_cast<T>` support — the concrete type as Any
    /// (mirror of Curve2dAdaptor::as_any).
    fn as_any(&self) -> &dyn std::any::Any;

    /// OCCT handle down-cast support (Arc<dyn> -> Arc<Concrete>).
    fn as_any_arc(self: Arc<Self>) -> Arc<dyn std::any::Any + Send + Sync>;
}

// ---------------------------------------------------------------------
// OCCT Geom2d_Line / Geom2d_Circle / ... — a plain kernel `Curve2d` used
// as a Geom2d_Curve handle inside the package.
// ---------------------------------------------------------------------

/// A kernel `Curve2d` as a `Geom2d_Curve` handle (architecture note in the
/// module docs).  The Bisector-only pure virtuals keep the OCCT failure path.
#[derive(Debug, Clone)]
pub struct Geom2dCurveHandle {
    pub curve: Curve2d,
}

impl Geom2dCurveHandle {
    /// OCCT `new Geom2d_Line(Lin2d)` etc. — wrap a kernel curve.
    pub fn new(curve: Curve2d) -> Self {
        Geom2dCurveHandle { curve }
    }

    /// OCCT `new Geom2d_Line(gp_Lin2d)`.
    pub fn line(l: &Line2d) -> Self {
        Geom2dCurveHandle {
            curve: Curve2d::Line(*l),
        }
    }

    /// OCCT `new Geom2d_Circle(gp_Circ2d)`.
    pub fn circle(c: &Circle2d) -> Self {
        Geom2dCurveHandle {
            curve: Curve2d::Circle(*c),
        }
    }

    /// OCCT `new Geom2d_Ellipse(gp_Elips2d)`.
    pub fn ellipse(e: &Ellipse2d) -> Self {
        Geom2dCurveHandle {
            curve: Curve2d::Ellipse(*e),
        }
    }

    /// OCCT `new Geom2d_Parabola(gp_Parab2d)`.
    pub fn parabola(p: &Parabola2d) -> Self {
        Geom2dCurveHandle {
            curve: Curve2d::Parabola(*p),
        }
    }

    /// OCCT `new Geom2d_Hyperbola(gp_Hypr2d)`.
    pub fn hyperbola(h: &Hyperbola2d) -> Self {
        Geom2dCurveHandle {
            curve: Curve2d::Hyperbola(*h),
        }
    }
}

impl BisectorCurve for Geom2dCurveHandle {
    fn value(&self, u: f64) -> DVec2 {
        // OCCT Geom2d_Curve::Value -> the kernel evaluation.
        self.curve.point_at(u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        geom2d_curve_tool::d1(&self.curve, u)
    }

    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        geom2d_curve_tool::d2(&self.curve, u)
    }

    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        geom2d_curve_tool::d3(&self.curve, u)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        geom2d_curve_tool::dn(&self.curve, u, n)
    }

    fn first_parameter(&self) -> f64 {
        geom2d_curve_tool::first_parameter(&self.curve)
    }

    fn last_parameter(&self) -> f64 {
        geom2d_curve_tool::last_parameter(&self.curve)
    }

    fn is_closed(&self) -> bool {
        // Kernel `Curve2d` implements both Curve2dEval and Curve2dAdaptor;
        // the evaluation trait is the OCCT Geom2d_Curve mirror.
        Curve2dEval::is_closed(&self.curve)
    }

    fn is_periodic(&self) -> bool {
        Curve2dEval::is_periodic(&self.curve)
    }

    fn continuity(&self) -> GeomAbsShape {
        // OCCT Geom2d_Conic/Bezier::Continuity() == GeomAbs_CN;
        // Geom2d_BSplineCurve carries its construction continuity (C2 for
        // the kernel default).  GAP: kernel Curve2d stores no continuity
        // field; the OCCT-typical values are mirrored here.
        match &self.curve {
            Curve2d::BSpline(_) => GeomAbsShape::C2,
            _ => GeomAbsShape::CN,
        }
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        Curve2dEval::reversed_parameter(&self.curve, u)
    }

    fn reverse(&mut self) {
        self.curve = rcad_kernel::geom::reverse_curve2d(&self.curve);
    }

    fn is_cn(&self, n: i32) -> bool {
        // OCCT: analytic curves are CN (always true); a BSpline is CN(N)
        // when N <= Degree.  GAP: kernel stores no explicit continuity.
        match &self.curve {
            Curve2d::BSpline(b) => n <= b.degree as i32,
            _ => true,
        }
    }

    fn transform(&mut self, t: &glam::DAffine2) {
        // GAP: kernel Curve2d has no general gp_Trsf2d transform
        // (OCCT Geom2d_Geometry::Transform).  Preserves the loud failure
        // path until the kernel transform lands.
        let _ = t;
        unimplemented!("GAP: kernel Curve2d lacks a gp_Trsf2d transform")
    }

    fn copy_curve(&self) -> Arc<dyn BisectorCurve> {
        Arc::new(Geom2dCurveHandle {
            curve: self.curve.clone(),
        })
    }

    fn kind(&self) -> CurveKind {
        match &self.curve {
            Curve2d::Line(_) => CurveKind::Line,
            Curve2d::Circle(_) => CurveKind::Circle,
            Curve2d::Ellipse(_) => CurveKind::Ellipse,
            Curve2d::Parabola(_) => CurveKind::Parabola,
            Curve2d::Hyperbola(_) => CurveKind::Hyperbola,
            _ => CurveKind::Other,
        }
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_arc(self: Arc<Self>) -> Arc<dyn std::any::Any + Send + Sync> {
        self
    }
}

// ---------------------------------------------------------------------
// OCCT Geom2d_TrimmedCurve (the part Bisector uses) — a trimmed curve over
// a polymorphic Bisector/Geom2d curve.  Architecture note in module docs.
// ---------------------------------------------------------------------

/// OCCT Geom2d_TrimmedCurve over a `Geom2d_Curve` handle.  Parameters are
/// stored raw (the OCCT `Sense` constructor flag affects only periodic
/// re-parameterization, which raw storage already represents).
#[derive(Clone)]
pub struct TrimmedCurve {
    pub basis_curve: Arc<dyn BisectorCurve>,
    pub first_parameter: f64,
    pub last_parameter: f64,
}

impl TrimmedCurve {
    /// OCCT Geom2d_TrimmedCurve(C, U1, U2).
    pub fn new(basis_curve: Arc<dyn BisectorCurve>, u1: f64, u2: f64) -> Self {
        TrimmedCurve {
            basis_curve,
            first_parameter: u1,
            last_parameter: u2,
        }
    }

    /// OCCT Geom2d_TrimmedCurve(C, U1, U2, Sense) — the periodic variant used
    /// by BisecAna for circles/ellipses (`Sense = thesense`).  Raw bounds are
    /// kept (see type doc); the sense flag only selects the OCCT periodic
    /// bound shift, already expressed by the passed bounds.
    pub fn new_with_sense(basis_curve: Arc<dyn BisectorCurve>, u1: f64, u2: f64, _sense: bool) -> Self {
        TrimmedCurve::new(basis_curve, u1, u2)
    }

    /// OCCT BasisCurve().
    pub fn basis(&self) -> &Arc<dyn BisectorCurve> {
        &self.basis_curve
    }

    /// OCCT SetTrim(U1, U2) (Geom2d_TrimmedCurve).
    pub fn set_trim(&mut self, u1: f64, u2: f64) {
        // OCCT: theFirst = U1, theLast = U2 (same-range re-evaluation is a
        // no-op for the plain case used by Bisector_BisecAna::SetTrim).
        self.first_parameter = u1;
        self.last_parameter = u2;
    }
}

impl BisectorCurve for TrimmedCurve {
    fn value(&self, u: f64) -> DVec2 {
        // OCCT Geom2d_TrimmedCurve::Value(U) -> BasisCurve()->Value(U).
        self.basis_curve.value(u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        self.basis_curve.d1(u)
    }

    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        self.basis_curve.d2(u)
    }

    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        self.basis_curve.d3(u)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        self.basis_curve.dn(u, n)
    }

    fn first_parameter(&self) -> f64 {
        self.first_parameter
    }

    fn last_parameter(&self) -> f64 {
        self.last_parameter
    }

    fn is_closed(&self) -> bool {
        // OCCT: IsClosed() -> BasisCurve()->IsClosed() within the range.
        self.basis_curve.is_closed()
    }

    fn is_periodic(&self) -> bool {
        self.basis_curve.is_periodic()
    }

    fn period(&self) -> f64 {
        self.basis_curve.period()
    }

    fn continuity(&self) -> GeomAbsShape {
        self.basis_curve.continuity()
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        // OCCT Geom2d_TrimmedCurve::ReversedParameter(U):
        // LastParameter() + FirstParameter() - U.
        self.last_parameter + self.first_parameter - u
    }

    fn reverse(&mut self) {
        // OCCT Geom2d_TrimmedCurve::Reverse(): reverse the basis, then map
        // the trim bounds through ReversedParameter.  (OCCT mutates the
        // shared basis handle; here the basis is cloned and replaced.)
        let mut basis = self.basis_curve.copy_curve();
        BisectorCurve::reverse(Arc::get_mut(&mut basis).expect("exclusive basis"));
        let new_first = basis.reversed_parameter(self.last_parameter);
        let new_last = basis.reversed_parameter(self.first_parameter);
        self.basis_curve = basis;
        self.first_parameter = new_first;
        self.last_parameter = new_last;
    }

    fn is_cn(&self, n: i32) -> bool {
        self.basis_curve.is_cn(n)
    }

    fn transform(&mut self, t: &glam::DAffine2) {
        // OCCT mutates the shared basis handle; here the basis is cloned and
        // replaced (same result for the exclusive owners used in-package).
        let mut basis = self.basis_curve.copy_curve();
        BisectorCurve::transform(Arc::get_mut(&mut basis).expect("exclusive basis"), t);
        self.basis_curve = basis;
    }

    fn copy_curve(&self) -> Arc<dyn BisectorCurve> {
        Arc::new(TrimmedCurve {
            basis_curve: self.basis_curve.copy_curve(),
            first_parameter: self.first_parameter,
            last_parameter: self.last_parameter,
        })
    }

    fn kind(&self) -> CurveKind {
        self.basis_curve.kind()
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_arc(self: Arc<Self>) -> Arc<dyn std::any::Any + Send + Sync> {
        self
    }
}

// ---------------------------------------------------------------------
// OCCT Geom2dAdaptor_Curve — Adaptor2d_Curve2d over a Geom2d_Curve handle.
// (Glue for the Geom2dInt_GInter chain; architecture note in module docs.)
// ---------------------------------------------------------------------

/// OCCT Geom2dAdaptor_Curve over a [`BisectorCurve`] handle ('static: owns
/// the handle, matching OCCT's owned `occ::handle`).
pub struct Geom2dCurveAdaptor {
    pub curve: Arc<dyn BisectorCurve>,
}

impl Geom2dCurveAdaptor {
    /// OCCT Geom2dAdaptor_Curve(C).
    pub fn new(curve: Arc<dyn BisectorCurve>) -> Self {
        Geom2dCurveAdaptor { curve }
    }
}

impl Curve2dAdaptor for Geom2dCurveAdaptor {
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    fn value(&self, u: f64) -> DVec2 {
        self.curve.value(u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        self.curve.d1(u)
    }

    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        self.curve.d2(u)
    }

    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        self.curve.d3(u)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        self.curve.dn(u, n)
    }

    fn get_type(&self) -> Curve2dType {
        match self.curve.kind() {
            CurveKind::Line => Curve2dType::Line,
            CurveKind::Circle => Curve2dType::Circle,
            _ => {
                // GAP: the trimmed-of-bisector wrapper has no OCCT conic
                // type; OtherCurve matches Geom2dAdaptor's fallback.
                Curve2dType::OtherCurve
            }
        }
    }

    fn nb_samples(&self) -> i32 {
        // OCCT Geom2dInt_Geom2dCurveTool::NbSamples default for non-sampled
        // curves (Geom2dInt_Geom2dCurveTool.cxx L72-77).
        10
    }

    fn resolution(&self, _r3d: f64) -> f64 {
        // OCCT Adaptor2d_Curve2d::Resolution default.
        1e-10
    }

    fn is_closed(&self) -> bool {
        self.curve.is_closed()
    }

    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    fn period(&self) -> f64 {
        self.curve.period()
    }

    fn nb_knots(&self) -> i32 {
        0
    }

    fn degree(&self) -> i32 {
        0
    }

    fn nb_poles(&self) -> i32 {
        0
    }

    fn circle(&self) -> Circle2d {
        // OCCT raises for non-conic; the package only calls it after a
        // Line/Circle kind test (IntCurve_IConicTool input).
        unreachable!("Geom2dCurveAdaptor::circle on non-conic")
    }

    fn line(&self) -> Line2d {
        unreachable!("Geom2dCurveAdaptor::line on non-line")
    }

    fn ellipse(&self) -> Ellipse2d {
        unreachable!("Geom2dCurveAdaptor::ellipse on non-conic")
    }

    fn parabola(&self) -> Parabola2d {
        unreachable!("Geom2dCurveAdaptor::parabola on non-conic")
    }

    fn hyperbola(&self) -> Hyperbola2d {
        unreachable!("Geom2dCurveAdaptor::hyperbola on non-conic")
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}
