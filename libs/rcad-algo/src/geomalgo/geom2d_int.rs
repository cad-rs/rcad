//! OCCT Geom2dInt / IntCurve / IntImpParGen — the imp-par (implicit conic vs
//! parametric curve) 2D intersection chain used by Geom2dInt_GInter and by
//! Geom2dGcc_Circ2d2TanRadGeo.
//!
//! 1:1 translations:
//!   - Adaptor2d_OffsetCurve.hxx/.cxx + Geom2d_OffsetCurveUtils.pxx — the
//!     offset-curve adaptor (the ThePCurve side of the intersection).
//!   - Geom2dInt_Geom2dCurveTool.cxx/.lxx — the parametric-curve tool
//!     (NbSamples/EpsX/FirstParameter/LastParameter/Value/D1/D2/D3/DN).
//!   - IntCurve_IConicTool.cxx — the implicit conic tool (the ImpTool side).
//!   - Extrema_GCurveLocator.hxx + Extrema_GenLocateExtPC.hxx +
//!     Extrema_GFuncExtPC.hxx + Geom2dInt_TheProjPCurOfGInter.cxx — the
//!     point-to-curve projection (FindParameter) used by FindV.
//!
//! The generic engines and the package statics live in their own modules:
//! `int_imp_par_gen` (IntImpParGen_Intersector.gxx + IntImpParGen.cxx +
//! MyImpParTool), `int_conic_curve_gen` (IntCurve_IntConicCurveGen
//! gxx/.lxx) and `user_int_conic_curve_gen` (IntCurve_UserIntConicCurveGen);
//! the GInter instantiations (Geom2dInt_Geom2dCurveTool /
//! Geom2dInt_TheProjPCurOfGInter / TheIntersectorOfTheIntConicCurveOfGInter
//! / TheIntConicCurveOfGInter) are the tool markers and concrete wrappers at
//! the end of this file.

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Curve2d, Curve2dEval, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative};

use super::int_res2d::{Domain as Res2dDomain, IntersectionBase};

/// gp::Resolution() — OCCT gp.hxx.
const GP_RESOLUTION: f64 = 1e-15;
/// Precision::Computational() — OCCT Precision.hxx.
const PRECISION_COMPUTATIONAL: f64 = 1e-14;
/// 2·π.
const PI2: f64 = std::f64::consts::TAU;

/// OCCT GeomAbs_CurveType — the 2D curve kinds needed by the adaptor chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Curve2dType {
    Line,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
    BezierCurve,
    BSplineCurve,
    OffsetCurve,
    OtherCurve,
}

/// OCCT Adaptor2d_Curve2d — the abstract 2D curve adaptor interface.
pub trait Curve2dAdaptor {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64;
    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64;
    /// OCCT Value(U).
    fn value(&self, u: f64) -> DVec2;
    /// OCCT D1(U, P, V).
    fn d1(&self, u: f64) -> (DVec2, DVec2);
    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2);
    /// OCCT D3(U, P, V1, V2, V3).
    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2);
    /// OCCT DN(U, N).
    fn dn(&self, u: f64, n: i32) -> DVec2;
    /// OCCT GetType().
    fn get_type(&self) -> Curve2dType;
    /// OCCT NbSamples().
    fn nb_samples(&self) -> i32;
    /// OCCT Resolution(R3d).
    fn resolution(&self, r3d: f64) -> f64;
    /// OCCT IsClosed().
    fn is_closed(&self) -> bool;
    /// OCCT IsPeriodic().
    fn is_periodic(&self) -> bool;
    /// OCCT Period().
    fn period(&self) -> f64;
    /// OCCT NbKnots().
    fn nb_knots(&self) -> i32;
    /// OCCT Degree().
    fn degree(&self) -> i32;
    /// OCCT NbPoles().
    fn nb_poles(&self) -> i32;
    /// OCCT Circle() — valid when GetType() == Circle.
    fn circle(&self) -> Circle2d;
    /// OCCT Line() — valid when GetType() == Line.
    fn line(&self) -> Line2d;
    /// OCCT Ellipse() — valid when GetType() == Ellipse.
    fn ellipse(&self) -> Ellipse2d;
    /// OCCT Parabola() — valid when GetType() == Parabola.
    fn parabola(&self) -> Parabola2d;
    /// OCCT Hyperbola() — valid when GetType() == Hyperbola.
    fn hyperbola(&self) -> Hyperbola2d;
    /// OCCT Adaptor2d_Curve2d::NbIntervals(GeomAbs_C1) — number of C1
    /// continuity intervals over the curve domain. All rcad curve kinds are
    /// C1-continuous across their domain, so the default is 1 (matching
    /// Geom2dAdaptor_Curve for non-composite curves).
    fn nb_intervals_c1(&self) -> i32 {
        1
    }
    /// OCCT `occ::down_cast<T>(handle)` support — the concrete adaptor as
    /// Any (BRepTopAdaptor_TopolTool::Initialize down-casts the restriction
    /// arc handle to BRepAdaptor_Curve2d).
    fn as_any(&self) -> &dyn std::any::Any;
}

/// OCCT GeomAbs_CurveType of a `Curve2d` (Geom2dAdaptor_Curve::GetType).
pub fn curve2d_type_of(c: &Curve2d) -> Curve2dType {
    match c {
        Curve2d::Line(_) => Curve2dType::Line,
        Curve2d::Circle(_) => Curve2dType::Circle,
        Curve2d::Ellipse(_) => Curve2dType::Ellipse,
        Curve2d::Parabola(_) => Curve2dType::Parabola,
        Curve2d::Hyperbola(_) => Curve2dType::Hyperbola,
        Curve2d::Bezier(_) => Curve2dType::BezierCurve,
        Curve2d::BSpline(_) => Curve2dType::BSplineCurve,
        Curve2d::Offset(_) => Curve2dType::OffsetCurve,
        _ => Curve2dType::OtherCurve,
    }
}

/// OCCT Geom2dAdaptor_Curve — adaptor over a Geom2d_Curve (here rcad's
/// `Curve2d`). Evaluations delegate to `Curve2dEval`.
impl Curve2dAdaptor for Curve2d {
    fn first_parameter(&self) -> f64 {
        self.default_domain()[0]
    }
    fn last_parameter(&self) -> f64 {
        self.default_domain()[1]
    }
    fn value(&self, u: f64) -> DVec2 {
        self.point_at(u)
    }
    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        (self.point_at(u), self.derivative_at(u))
    }
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        (self.point_at(u), self.derivative_at(u), self.derivative2_at(u))
    }
    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        (
            self.point_at(u),
            self.derivative_at(u),
            self.derivative2_at(u),
            self.derivative3_at(u),
        )
    }
    fn dn(&self, u: f64, n: i32) -> DVec2 {
        match n {
            1 => self.derivative_at(u),
            2 => self.derivative2_at(u),
            _ => self.derivative3_at(u),
        }
    }
    fn get_type(&self) -> Curve2dType {
        curve2d_type_of(self)
    }
    fn nb_samples(&self) -> i32 {
        let nbs = match self {
            Curve2d::Bezier(b) => 3 + b.control_points.len() as i32,
            Curve2d::BSpline(b) => b.knots.len() as i32 * b.degree as i32,
            _ => 20,
        };
        nbs.max(20).min(300)
    }
    fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }
    fn is_closed(&self) -> bool {
        Curve2dEval::is_closed(self)
    }
    fn is_periodic(&self) -> bool {
        Curve2dEval::is_periodic(self)
    }
    fn period(&self) -> f64 {
        match self {
            Curve2d::Circle(_) | Curve2d::Ellipse(_) => PI2,
            _ => 0.0,
        }
    }
    fn nb_knots(&self) -> i32 {
        match self {
            Curve2d::BSpline(b) => b.knots.len() as i32,
            _ => 0,
        }
    }
    fn degree(&self) -> i32 {
        match self {
            Curve2d::BSpline(b) => b.degree as i32,
            Curve2d::Bezier(b) => (b.control_points.len() as i32 - 1).max(0),
            _ => 0,
        }
    }
    fn nb_poles(&self) -> i32 {
        match self {
            Curve2d::Bezier(b) => b.control_points.len() as i32,
            Curve2d::BSpline(b) => b.control_points.len() as i32,
            _ => 0,
        }
    }
    fn circle(&self) -> Circle2d {
        match self {
            Curve2d::Circle(c) => *c,
            _ => panic!("Standard_NoSuchObject: Curve2d::Circle"),
        }
    }
    fn line(&self) -> Line2d {
        match self {
            Curve2d::Line(l) => *l,
            _ => panic!("Standard_NoSuchObject: Curve2d::Line"),
        }
    }
    fn ellipse(&self) -> Ellipse2d {
        match self {
            Curve2d::Ellipse(e) => *e,
            _ => panic!("Standard_NoSuchObject: Curve2d::Ellipse"),
        }
    }
    fn parabola(&self) -> Parabola2d {
        match self {
            Curve2d::Parabola(p) => *p,
            _ => panic!("Standard_NoSuchObject: Curve2d::Parabola"),
        }
    }
    fn hyperbola(&self) -> Hyperbola2d {
        match self {
            Curve2d::Hyperbola(h) => *h,
            _ => panic!("Standard_NoSuchObject: Curve2d::Hyperbola"),
        }
    }
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// OCCT Adaptor2d_OffsetCurve (Adaptor2d_OffsetCurve.cxx + the
/// Geom2d_OffsetCurveUtils.pxx formulas) — an algorithmic 2D offset curve.
#[derive(Debug, Clone)]
pub struct AdaptorOffsetCurve {
    /// The basis curve (OCCT handle<Adaptor2d_Curve2d> myCurve).
    pub base: Box<Curve2d>,
    /// The offset distance (OCCT myOffset).
    pub offset: f64,
    /// First parameter (OCCT myFirst).
    pub first: f64,
    /// Last parameter (OCCT myLast).
    pub last: f64,
}

impl AdaptorOffsetCurve {
    /// OCCT Adaptor2d_OffsetCurve(C, Offset).
    pub fn new(base: Curve2d, offset: f64) -> Self {
        let dom = base.default_domain();
        AdaptorOffsetCurve {
            base: Box::new(base),
            offset,
            first: dom[0],
            last: dom[1],
        }
    }

    /// OCCT Adaptor2d_OffsetCurve(C, Offset, WFirst, WLast).
    pub fn new_bounded(base: Curve2d, offset: f64, first: f64, last: f64) -> Self {
        AdaptorOffsetCurve {
            base: Box::new(base),
            offset,
            first,
            last,
        }
    }

    /// OCCT Geom2d_OffsetCurveUtils::CalculateD0 — P(u) = p(u) + Offset·N/|N|
    /// with N = (p'(u).Y, -p'(u).X). Returns false when |p'| is degenerate
    /// (the point is then left unchanged).
    fn calculate_d0(p: &mut DVec2, d1: DVec2, offset: f64) -> bool {
        if d1.length_squared() <= GP_RESOLUTION {
            return false;
        }
        let normal = DVec2::new(d1.y, -d1.x).normalize_or_zero();
        *p += normal * offset;
        true
    }

    /// OCCT Geom2d_OffsetCurveUtils::CalculateD1.
    fn calculate_d1(p: &mut DVec2, d1: &mut DVec2, d2: DVec2, offset: f64) -> bool {
        let ndir = DVec2::new(d1.y, -d1.x);
        let dndir = DVec2::new(d2.y, -d2.x);
        let r2 = ndir.length_squared();
        let r = r2.sqrt();
        let r3 = r * r2;
        let dr = ndir.dot(dndir);
        let (mut dn, mut n) = (dndir, ndir);
        if r3 <= GP_RESOLUTION {
            if r2 <= GP_RESOLUTION {
                return false;
            }
            dn = dn * r - n * (dr / r);
            dn *= offset / r2;
        } else {
            dn = dn * (offset / r) - n * (offset * dr / r3);
        }
        n *= offset / r;
        *p += n;
        *d1 += dn;
        true
    }

    /// OCCT Geom2d_OffsetCurveUtils::CalculateD2.
    fn calculate_d2(
        p: &mut DVec2,
        d1: &mut DVec2,
        d2: &mut DVec2,
        d3: DVec2,
        is_dir_change: bool,
        offset: f64,
    ) -> bool {
        let ndir = DVec2::new(d1.y, -d1.x);
        let mut dndir = DVec2::new(d2.y, -d2.x);
        let mut d2ndir = DVec2::new(d3.y, -d3.x);
        let r2 = ndir.length_squared();
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r4 = r2 * r2;
        let r5 = r3 * r2;
        let dr = ndir.dot(dndir);
        let d2r = ndir.dot(d2ndir) + dndir.dot(dndir);
        if r5 <= GP_RESOLUTION {
            if r4 <= GP_RESOLUTION {
                return false;
            }
            d2ndir = d2ndir - dndir * (2.0 * dr / r2) + ndir * ((3.0 * dr * dr) / r4 - d2r / r2);
            d2ndir *= offset / r;
            dndir = dndir * r - ndir * (dr / r);
            dndir *= offset / r2;
        } else {
            d2ndir = d2ndir * (offset / r) - dndir * (2.0 * offset * dr / r3)
                + ndir * (offset * ((3.0 * dr * dr) / r5 - d2r / r3));
            dndir = dndir * (offset / r) - ndir * (offset * dr / r3);
        }
        let n = ndir * (offset / r);
        *p += n;
        *d1 += dndir;
        if is_dir_change {
            *d2 = -*d2;
        }
        *d2 += d2ndir;
        true
    }

    /// OCCT Geom2d_OffsetCurveUtils::CalculateD3.
    fn calculate_d3(
        p: &mut DVec2,
        d1: &mut DVec2,
        d2: &mut DVec2,
        d3: &mut DVec2,
        d4: DVec2,
        is_dir_change: bool,
        offset: f64,
    ) -> bool {
        let ndir = DVec2::new(d1.y, -d1.x);
        let mut dndir = DVec2::new(d2.y, -d2.x);
        let mut d2ndir = DVec2::new(d3.y, -d3.x);
        let mut d3ndir = DVec2::new(d4.y, -d4.x);
        let r2 = ndir.length_squared();
        let r = r2.sqrt();
        let r3 = r2 * r;
        let r4 = r2 * r2;
        let r5 = r3 * r2;
        let r6 = r3 * r3;
        let r7 = r5 * r2;
        let dr = ndir.dot(dndir);
        let d2r = ndir.dot(d2ndir) + dndir.dot(dndir);
        let d3r = ndir.dot(d3ndir) + 3.0 * dndir.dot(d2ndir);
        if r7 <= GP_RESOLUTION {
            if r6 <= GP_RESOLUTION {
                return false;
            }
            d3ndir = d3ndir - d2ndir * (3.0 * dr / r2)
                - dndir * (3.0 * (d2r / r2 + dr * dr / r4))
                + ndir * (6.0 * dr * dr / r4 + 6.0 * dr * d2r / r4
                    - 15.0 * dr * dr * dr / r6 - d3r);
            d3ndir *= offset / r;
            d2ndir = d2ndir - dndir * (2.0 * dr / r2)
                - ndir * ((3.0 * dr * dr / r4) - d2r / r2);
            d2ndir *= offset / r;
            dndir = dndir * r - ndir * (dr / r);
            dndir *= offset / r2;
        } else {
            d3ndir = d3ndir * (offset / r) - d2ndir * (3.0 * offset * dr / r3)
                - dndir * (3.0 * offset * (d2r / r3 + dr * dr / r5))
                + ndir * (offset * (6.0 * dr * dr / r5 + 6.0 * dr * d2r / r5
                    - 15.0 * dr * dr * dr / r7 - d3r));
            d2ndir = d2ndir * (offset / r) - dndir * (2.0 * offset * dr / r3)
                - ndir * (offset * ((3.0 * dr * dr) / r5 - d2r / r3));
            dndir = dndir * (offset / r) - ndir * (offset * dr / r3);
        }
        let n = ndir * (offset / r);
        *p += n;
        *d1 += dndir;
        *d2 += d2ndir;
        if is_dir_change {
            *d3 = -*d3;
        }
        *d3 += d3ndir;
        true
    }
}

impl Curve2dAdaptor for AdaptorOffsetCurve {
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }

    /// OCCT Adaptor2d_OffsetCurve::Value (L295-309).
    fn value(&self, u: f64) -> DVec2 {
        if self.offset != 0.0 {
            let (mut p, v) = self.base.d1(u);
            let _ = AdaptorOffsetCurve::calculate_d0(&mut p, v, self.offset);
            p
        } else {
            self.base.value(u)
        }
    }

    /// OCCT Adaptor2d_OffsetCurve::D1 (L320-332).
    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        if self.offset != 0.0 {
            let (mut p, mut v, v2) = self.base.d2(u);
            let _ = AdaptorOffsetCurve::calculate_d1(&mut p, &mut v, v2, self.offset);
            (p, v)
        } else {
            self.base.d1(u)
        }
    }

    /// OCCT Adaptor2d_OffsetCurve::D2 (L336-348).
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        if self.offset != 0.0 {
            let (mut p, mut v1, mut v2, v3) = self.base.d3(u);
            let _ = AdaptorOffsetCurve::calculate_d2(&mut p, &mut v1, &mut v2, v3, false, self.offset);
            (p, v1, v2)
        } else {
            self.base.d2(u)
        }
    }

    /// OCCT Adaptor2d_OffsetCurve::D3 (L352-368).
    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        if self.offset != 0.0 {
            let v4 = self.base.dn(u, 4);
            let (mut p, mut v1, mut v2, mut v3) = self.base.d3(u);
            let _ = AdaptorOffsetCurve::calculate_d3(&mut p, &mut v1, &mut v2, &mut v3, v4, false, self.offset);
            (p, v1, v2, v3)
        } else {
            self.base.d3(u)
        }
    }

    fn dn(&self, _u: f64, _n: i32) -> DVec2 {
        panic!("Standard_NotImplemented: Adaptor2d_OffsetCurve::DN");
    }

    /// OCCT Adaptor2d_OffsetCurve::GetType (L386-408).
    fn get_type(&self) -> Curve2dType {
        if self.offset == 0.0 {
            self.base.get_type()
        } else {
            match self.base.get_type() {
                Curve2dType::Line => Curve2dType::Line,
                Curve2dType::Circle => Curve2dType::Circle,
                _ => Curve2dType::OffsetCurve,
            }
        }
    }

    /// OCCT Adaptor2d_OffsetCurve::Line (L412-425).
    fn line(&self) -> Line2d {
        if self.get_type() == Curve2dType::Line {
            let (p, v) = self.d1(0.0);
            Line2d::new(p, v)
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Line");
        }
    }

    /// OCCT Adaptor2d_OffsetCurve::Circle (L429-468).
    fn circle(&self) -> Circle2d {
        if self.get_type() == Curve2dType::Circle {
            if self.offset == 0.0 {
                return self.base.circle();
            }
            let c1 = self.base.circle();
            let mut radius = c1.radius;
            let xd = c1.x_dir;
            let yd = c1.y_dir;
            let crossed = xd.x * yd.y - xd.y * yd.x;
            let signe = if crossed > 0.0 { 1.0 } else { -1.0 };
            radius += signe * self.offset;
            if radius > 0.0 {
                Circle2d {
                    center: c1.center,
                    x_dir: xd,
                    y_dir: yd,
                    radius,
                }
            } else if radius < 0.0 {
                Circle2d {
                    center: c1.center,
                    x_dir: -xd,
                    y_dir: yd,
                    radius: -radius,
                }
            } else {
                panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Circle (null radius)");
            }
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Circle");
        }
    }

    fn ellipse(&self) -> Ellipse2d {
        if self.base.get_type() == Curve2dType::Ellipse && self.offset == 0.0 {
            self.base.ellipse()
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Ellipse");
        }
    }

    fn parabola(&self) -> Parabola2d {
        if self.base.get_type() == Curve2dType::Parabola && self.offset == 0.0 {
            self.base.parabola()
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Parabola");
        }
    }

    fn hyperbola(&self) -> Hyperbola2d {
        if self.base.get_type() == Curve2dType::Hyperbola && self.offset == 0.0 {
            self.base.hyperbola()
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Hyperbola");
        }
    }

    fn degree(&self) -> i32 {
        let t = self.base.get_type();
        if (t == Curve2dType::BezierCurve || t == Curve2dType::BSplineCurve) && self.offset == 0.0 {
            self.base.degree()
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::Degree");
        }
    }

    fn nb_poles(&self) -> i32 {
        let t = self.base.get_type();
        if (t == Curve2dType::BezierCurve || t == Curve2dType::BSplineCurve) && self.offset == 0.0 {
            self.base.nb_poles()
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::NbPoles");
        }
    }

    fn nb_knots(&self) -> i32 {
        if self.offset == 0.0 {
            self.base.nb_knots()
        } else {
            panic!("Standard_NoSuchObject: Adaptor2d_OffsetCurve::NbKnots");
        }
    }

    fn is_closed(&self) -> bool {
        if self.offset == 0.0 {
            Curve2dEval::is_closed(&*self.base)
        } else {
            false
        }
    }

    fn is_periodic(&self) -> bool {
        Curve2dEval::is_periodic(&*self.base)
    }

    fn period(&self) -> f64 {
        self.base.period()
    }

    /// OCCT Adaptor2d_OffsetCurve::Resolution (L379-383):
    /// Precision::PConfusion(R3d) = R3d.
    fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }

    /// OCCT static nbPoints (Adaptor2d_OffsetCurve.cxx L585-604).
    fn nb_samples(&self) -> i32 {
        let mut nbs = 20;
        match self.base.get_type() {
            Curve2dType::BezierCurve => nbs = nbs.max(3 + self.base.nb_poles()),
            Curve2dType::BSplineCurve => nbs = nbs.max(self.base.nb_knots() * self.base.degree()),
            _ => {}
        }
        nbs.min(300)
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// OCCT Geom2dInt_Geom2dCurveTool (Geom2dInt_Geom2dCurveTool.cxx/.lxx) — the
/// parametric-curve tool for the Geom2dInt chain.
pub mod geom2d_curve_tool {
    use super::*;

    /// OCCT Geom2dInt_Geom2dCurveTool::GetType (L32-35).
    pub fn get_type(c: &dyn Curve2dAdaptor) -> Curve2dType {
        c.get_type()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::NbSamples(C, U0, U1) (L23-70).
    pub fn nb_samples_2(c: &dyn Curve2dAdaptor, u0: f64, u1: f64) -> i32 {
        let typ_c = c.get_type();
        let mut nbs = c.nb_samples();
        if typ_c == Curve2dType::BSplineCurve {
            let t = c.last_parameter() - c.first_parameter();
            if t > 1.0e-7 {
                let mut t1 = u1 - u0;
                if t1 < 0.0 {
                    t1 = -t1;
                }
                nbs = c.nb_knots();
                nbs *= c.degree();
                let anb = t1 / t * nbs as f64;
                nbs = anb as i32;
                let a_min_pnt_nb = (c.degree() + 1).max(4);
                if nbs < a_min_pnt_nb {
                    nbs = a_min_pnt_nb;
                }
            }
        } else if typ_c == Curve2dType::Circle {
            // Try to reach deflection = eps*R, eps = 0.01
            let min_r = 1.0;
            let r = c.circle().radius;
            if r > min_r {
                let angl = 0.283079; // 2.*acos(1. - eps)
                let n = ((u1 - u0).abs() / angl) as i32;
                nbs = n.max(nbs);
            }
        }
        if nbs > 300 {
            nbs = 300;
        }
        nbs
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::NbSamples(C) (L73-91).
    pub fn nb_samples(c: &dyn Curve2dAdaptor) -> i32 {
        let mut nbs = c.nb_samples();
        let typ_c = c.get_type();
        if typ_c == Curve2dType::Circle {
            let min_r = 1.0;
            let r = c.circle().radius;
            if r > min_r {
                let angl = 0.283079;
                let n = ((c.last_parameter() - c.first_parameter()) / angl) as i32;
                nbs = n.max(nbs);
            }
        }
        nbs
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::EpsX(C) (L134-137) — tolerance used by
    /// the mathematical algorithms.
    pub fn eps_x(_c: &dyn Curve2dAdaptor) -> f64 {
        1.0e-10
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::EpsX(C, Eps_XYZ) (L140-143).
    pub fn eps_x_2(c: &dyn Curve2dAdaptor, eps_xyz: f64) -> f64 {
        c.resolution(eps_xyz)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::FirstParameter (L120-123).
    pub fn first_parameter(c: &dyn Curve2dAdaptor) -> f64 {
        c.first_parameter()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::LastParameter (L126-129).
    pub fn last_parameter(c: &dyn Curve2dAdaptor) -> f64 {
        c.last_parameter()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::NbIntervals
    /// (Geom2dInt_Geom2dCurveTool.lxx L169-178) — the
    /// Adaptor2d_Curve2d::NbIntervals(GeomAbs_C1) count.
    pub fn nb_intervals(c: &dyn Curve2dAdaptor) -> i32 {
        c.nb_intervals_c1()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Intervals (lxx L146-154) — fills an
    /// NCollection_Array1(1, NbIntervals+1). With the single C1 interval the
    /// bounds are [FirstParameter, LastParameter]
    /// (Geom2dAdaptor_Curve::Intervals). `tab` is 0-based storage for the
    /// 1-based array: logical index i lives at tab[i - 1]; when the C1
    /// interval machinery lands (BSpline C1 knots) this grows with it.
    pub fn intervals(c: &dyn Curve2dAdaptor, tab: &mut [f64]) {
        let n = c.nb_intervals_c1() as usize;
        tab[0] = c.first_parameter();
        tab[n] = c.last_parameter();
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::GetInterval (lxx L158-166) —
    /// a = Tab(i), b = Tab(i+1).
    pub fn get_interval(tab: &[f64], i: usize) -> (f64, f64) {
        (tab[i - 1], tab[i])
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Value (L68-71).
    pub fn value(c: &dyn Curve2dAdaptor, u: f64) -> DVec2 {
        c.value(u)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::D0 (L74-77).
    pub fn d0(c: &dyn Curve2dAdaptor, u: f64) -> DVec2 {
        c.value(u)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::D1 (L80-86).
    pub fn d1(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2) {
        c.d1(u)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::D2 (L89-97).
    pub fn d2(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2, DVec2) {
        c.d2(u)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::D3 (L100-109).
    pub fn d3(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        c.d3(u)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::DN (L112-117).
    pub fn dn(c: &dyn Curve2dAdaptor, u: f64, n: i32) -> DVec2 {
        c.dn(u, n)
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Line (L38-41).
    pub fn line(c: &dyn Curve2dAdaptor) -> Line2d {
        c.line()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Circle (L44-47).
    pub fn circle(c: &dyn Curve2dAdaptor) -> Circle2d {
        c.circle()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Ellipse (lxx) — the adaptor accessor.
    pub fn ellipse(c: &dyn Curve2dAdaptor) -> Ellipse2d {
        c.ellipse()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Parabola (lxx) — the adaptor accessor.
    pub fn parabola(c: &dyn Curve2dAdaptor) -> Parabola2d {
        c.parabola()
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::Hyperbola (lxx) — the adaptor accessor.
    pub fn hyperbola(c: &dyn Curve2dAdaptor) -> Hyperbola2d {
        c.hyperbola()
    }
}

/// ElCLib 2D helpers (ElCLib.cxx) — conic value/D1/D2/parameter on a 2D frame.
pub mod elclib2d {
    use super::*;

    /// OCCT normalizeAngle (ElCLib.cxx L56-72).
    pub fn normalize_angle(a: &mut f64) {
        while *a < -PRECISION_COMPUTATIONAL {
            *a += PI2;
        }
        while *a > PI2 * (1.0 + GP_RESOLUTION) {
            *a -= PI2;
        }
        if *a < 0.0 {
            *a = 0.0;
        }
    }

    /// OCCT ElCLib::LineValue (L521-526).
    pub fn line_value(loc: DVec2, dir: DVec2, u: f64) -> DVec2 {
        DVec2::new(u * dir.x + loc.x, u * dir.y + loc.y)
    }

    /// OCCT ElCLib::LineD1 (L592-598).
    pub fn line_d1(loc: DVec2, dir: DVec2, u: f64) -> (DVec2, DVec2) {
        (DVec2::new(u * dir.x + loc.x, u * dir.y + loc.y), dir)
    }

    /// OCCT ElCLib::LineParameter (L1276-1281).
    pub fn line_parameter(loc: DVec2, dir: DVec2, p: DVec2) -> f64 {
        let coord = p - loc;
        coord.dot(dir)
    }

    /// OCCT ElCLib::CircleValue (L530-539).
    pub fn circle_value(loc: DVec2, xdir: DVec2, ydir: DVec2, radius: f64, u: f64) -> DVec2 {
        let a1 = radius * u.cos();
        let a2 = radius * u.sin();
        DVec2::new(a1 * xdir.x + a2 * ydir.x + loc.x, a1 * xdir.y + a2 * ydir.y + loc.y)
    }

    /// OCCT ElCLib::CircleD1 (L602-619).
    pub fn circle_d1(loc: DVec2, xdir: DVec2, ydir: DVec2, radius: f64, u: f64) -> (DVec2, DVec2) {
        let xc = radius * u.cos();
        let yc = radius * u.sin();
        let p = DVec2::new(xc * xdir.x + yc * ydir.x + loc.x, xc * xdir.y + yc * ydir.y + loc.y);
        let v1 = DVec2::new(-yc * xdir.x + xc * ydir.x, -yc * xdir.y + xc * ydir.y);
        (p, v1)
    }

    /// OCCT ElCLib::CircleD2 (L694-716).
    pub fn circle_d2(loc: DVec2, xdir: DVec2, ydir: DVec2, radius: f64, u: f64) -> (DVec2, DVec2, DVec2) {
        let xc = radius * u.cos();
        let yc = radius * u.sin();
        let p = DVec2::new(xc * xdir.x + yc * ydir.x + loc.x, xc * xdir.y + yc * ydir.y + loc.y);
        let v1 = DVec2::new(-yc * xdir.x + xc * ydir.x, -yc * xdir.y + xc * ydir.y);
        let v2 = DVec2::new(-(xc * xdir.x + yc * ydir.x), -(xc * xdir.y + yc * ydir.y));
        (p, v1, v2)
    }

    /// OCCT ElCLib::CircleParameter (L1285-1291) with gp_Vec2d::Angle /
    /// Crossed conventions.
    pub fn circle_parameter(loc: DVec2, xdir: DVec2, ydir: DVec2, p: DVec2) -> f64 {
        let v = p - loc;
        // gp_Vec2d::Angle(a, b) = atan2(a.Crossed(b), a.Dot(b)).
        let mut teta = (xdir.x * v.y - xdir.y * v.x).atan2(xdir.dot(v));
        let crossed = xdir.x * ydir.y - xdir.y * ydir.x;
        if crossed < 0.0 {
            teta = -teta;
        }
        normalize_angle(&mut teta);
        teta
    }

    /// OCCT ElCLib::EllipseValue (L543-555).
    pub fn ellipse_value(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, u: f64) -> DVec2 {
        let a1 = major * u.cos();
        let a2 = minor * u.sin();
        DVec2::new(a1 * xdir.x + a2 * ydir.x + loc.x, a1 * xdir.y + a2 * ydir.y + loc.y)
    }

    /// OCCT ElCLib::EllipseD1 (L623-642).
    pub fn ellipse_d1(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, u: f64) -> (DVec2, DVec2) {
        let xc = u.cos();
        let yc = u.sin();
        let p = DVec2::new(
            xc * major * xdir.x + yc * minor * ydir.x + loc.x,
            xc * major * xdir.y + yc * minor * ydir.y + loc.y,
        );
        let v1 = DVec2::new(-yc * major * xdir.x + xc * minor * ydir.x, -yc * major * xdir.y + xc * minor * ydir.y);
        (p, v1)
    }

    /// OCCT ElCLib::EllipseD2 (L720-746).
    pub fn ellipse_d2(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, u: f64) -> (DVec2, DVec2, DVec2) {
        let xc = u.cos();
        let yc = u.sin();
        let p = DVec2::new(
            xc * major * xdir.x + yc * minor * ydir.x + loc.x,
            xc * major * xdir.y + yc * minor * ydir.y + loc.y,
        );
        let v1 = DVec2::new(-yc * major * xdir.x + xc * minor * ydir.x, -yc * major * xdir.y + xc * minor * ydir.y);
        let v2 = DVec2::new(
            -(xc * major * xdir.x + yc * minor * ydir.x),
            -(xc * major * xdir.y + yc * minor * ydir.y),
        );
        (p, v1, v2)
    }

    /// OCCT ElCLib::EllipseParameter (L1295-1311).
    pub fn ellipse_parameter(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, p: DVec2) -> f64 {
        let op = p - loc;
        let om = xdir * op.dot(xdir) + ydir * (op.dot(ydir) * (major / minor));
        let mut teta = (xdir.x * om.y - xdir.y * om.x).atan2(xdir.dot(om));
        let crossed = xdir.x * ydir.y - xdir.y * ydir.x;
        if crossed < 0.0 {
            teta = -teta;
        }
        normalize_angle(&mut teta);
        teta
    }

    /// OCCT ElCLib::ParabolaValue (L575-588).
    pub fn parabola_value(loc: DVec2, xdir: DVec2, ydir: DVec2, focal: f64, u: f64) -> DVec2 {
        if focal.abs() <= GP_RESOLUTION {
            return DVec2::new(u * xdir.x + loc.x, u * xdir.y + loc.y);
        }
        let a1 = u * u / (4.0 * focal);
        DVec2::new(a1 * xdir.x + u * ydir.x + loc.x, a1 * xdir.y + u * ydir.y + loc.y)
    }

    /// OCCT ElCLib::ParabolaD1 (L669-690).
    pub fn parabola_d1(loc: DVec2, xdir: DVec2, ydir: DVec2, focal: f64, u: f64) -> (DVec2, DVec2) {
        if focal.abs() <= GP_RESOLUTION {
            let v1 = xdir;
            let p = DVec2::new(u * xdir.x + loc.x, u * xdir.y + loc.y);
            return (p, v1);
        }
        let v1 = DVec2::new(u / (2.0 * focal) * xdir.x + ydir.x, u / (2.0 * focal) * xdir.y + ydir.y);
        let p = DVec2::new(
            (u * u) / (4.0 * focal) * xdir.x + u * ydir.x + loc.x,
            (u * u) / (4.0 * focal) * xdir.y + u * ydir.y + loc.y,
        );
        (p, v1)
    }

    /// OCCT ElCLib::ParabolaD2 (L779-805).
    pub fn parabola_d2(loc: DVec2, xdir: DVec2, ydir: DVec2, focal: f64, u: f64) -> (DVec2, DVec2, DVec2) {
        if focal.abs() <= GP_RESOLUTION {
            let v2 = DVec2::ZERO;
            let v1 = xdir;
            let p = DVec2::new(u * xdir.x + loc.x, u * xdir.y + loc.y);
            return (p, v1, v2);
        }
        let v2 = xdir * (1.0 / (2.0 * focal));
        let v1 = u * v2 + ydir;
        let p = DVec2::new(
            u * u / (4.0 * focal) * xdir.x + u * ydir.x + loc.x,
            u * u / (4.0 * focal) * xdir.y + u * ydir.y + loc.y,
        );
        (p, v1, v2)
    }

    /// OCCT ElCLib::ParabolaParameter (L1331-1335).
    pub fn parabola_parameter(loc: DVec2, ydir: DVec2, p: DVec2) -> f64 {
        (p - loc).dot(ydir)
    }

    /// OCCT ElCLib::HyperbolaValue (L559-571).
    pub fn hyperbola_value(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, u: f64) -> DVec2 {
        let a1 = major * u.cosh();
        let a2 = minor * u.sinh();
        DVec2::new(a1 * xdir.x + a2 * ydir.x + loc.x, a1 * xdir.y + a2 * ydir.y + loc.y)
    }

    /// OCCT ElCLib::HyperbolaD1 (L646-665).
    pub fn hyperbola_d1(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, u: f64) -> (DVec2, DVec2) {
        let xc = u.cosh();
        let yc = u.sinh();
        let p = DVec2::new(
            xc * major * xdir.x + yc * minor * ydir.x + loc.x,
            xc * major * xdir.y + yc * minor * ydir.y + loc.y,
        );
        let v1 = DVec2::new(yc * major * xdir.x + xc * minor * ydir.x, yc * major * xdir.y + xc * minor * ydir.y);
        (p, v1)
    }

    /// OCCT ElCLib::HyperbolaD2 (L750-775).
    pub fn hyperbola_d2(loc: DVec2, xdir: DVec2, ydir: DVec2, major: f64, minor: f64, u: f64) -> (DVec2, DVec2, DVec2) {
        let xc = u.cosh();
        let yc = u.sinh();
        let p = DVec2::new(
            xc * major * xdir.x + yc * minor * ydir.x + loc.x,
            xc * major * xdir.y + yc * minor * ydir.y + loc.y,
        );
        let v1 = DVec2::new(yc * major * xdir.x + xc * minor * ydir.x, yc * major * xdir.y + xc * minor * ydir.y);
        let v2 = DVec2::new(
            xc * major * xdir.x + yc * minor * ydir.x,
            xc * major * xdir.y + yc * minor * ydir.y,
        );
        (p, v1, v2)
    }

    /// OCCT ElCLib::HyperbolaParameter (L1315-1327).
    pub fn hyperbola_parameter(loc: DVec2, ydir: DVec2, minor: f64, p: DVec2) -> f64 {
        let sht = (p - loc).dot(ydir) / minor;
        sht.asinh()
    }
}

/// OCCT IntCurve_IConicTool (IntCurve_IConicTool.hxx/.cxx) — implementation of
/// the ImpTool from IntImpParGen for conics of gp.
#[derive(Debug, Clone)]
pub struct IConicTool {
    /// Line: a/b/c coefficients; Ellipse: a/b/c; Circle: r/x0/y0; Parabola:
    /// f/2p; Hyperbola: a/b.
    prm1: f64,
    prm2: f64,
    prm3: f64,
    /// The conic frame (loc, xdir, ydir).
    axis: (DVec2, DVec2, DVec2),
    /// GeomAbs_CurveType of the conic.
    typ: Curve2dType,
    /// Abs_To_Object — maps absolute coordinates into the object frame whose
    /// X axis is the conic X axis (loc, xdir, ydir).
    abs_to_object: (DVec2, DVec2, DVec2),
}

impl IConicTool {
    /// Object-frame coordinates of an absolute point (P.Transform(Abs_To_Object)).
    fn object_coords(&self, p: DVec2) -> DVec2 {
        let (loc, xdir, ydir) = self.abs_to_object;
        let d = p - loc;
        DVec2::new(d.dot(xdir), d.dot(ydir))
    }

    /// Absolute coordinates of an object-frame point (Object_To_Abs).
    fn absolute_from_object(&self, po: DVec2) -> DVec2 {
        let (loc, xdir, ydir) = self.abs_to_object;
        loc + xdir * po.x + ydir * po.y
    }

    /// OCCT IntCurve_IConicTool(gp_Lin2d) (L75-83).
    pub fn new_line(line: &Line2d) -> Self {
        // gp_Lin2d::Coefficients (gp_Lin2d.hxx L91-96).
        let a = line.direction.y;
        let b = -line.direction.x;
        let c = -(a * line.origin.x + b * line.origin.y);
        // gp_Ax22d(Line.Position(), true): xdir = direction, ydir = +90°.
        let ydir = DVec2::new(-line.direction.y, line.direction.x);
        IConicTool {
            prm1: a,
            prm2: b,
            prm3: c,
            axis: (line.origin, line.direction, ydir),
            typ: Curve2dType::Line,
            abs_to_object: (line.origin, line.direction, ydir),
        }
    }

    /// OCCT IntCurve_IConicTool(gp_Circ2d) (L100-111).
    pub fn new_circle(c: &Circle2d) -> Self {
        IConicTool {
            prm1: c.radius,
            prm2: c.center.x,
            prm3: c.center.y,
            axis: (c.center, c.x_dir, c.y_dir),
            typ: Curve2dType::Circle,
            abs_to_object: (c.center, c.x_dir, c.y_dir),
        }
    }

    /// OCCT IntCurve_IConicTool(gp_Elips2d) (L86-97).
    pub fn new_ellipse(e: &Ellipse2d) -> Self {
        let a = e.major_radius;
        let b = e.minor_radius;
        let minor_dir = DVec2::new(-e.major_dir.y, e.major_dir.x);
        IConicTool {
            prm1: a,
            prm2: b,
            prm3: (a * a - b * b).sqrt(),
            axis: (e.center, e.major_dir, minor_dir),
            typ: Curve2dType::Ellipse,
            abs_to_object: (e.center, e.major_dir, minor_dir),
        }
    }

    /// OCCT IntCurve_IConicTool(gp_Parab2d) (L114-124). rcad's focal_param is
    /// the semi-latus rectum p (x = t²/(2p)); OCCT Focal() = p/2.
    pub fn new_parabola(p: &Parabola2d) -> Self {
        let focal = p.focal_param / 2.0;
        let perp = DVec2::new(-p.axis_dir.y, p.axis_dir.x);
        IConicTool {
            prm1: focal,
            prm2: 4.0 * focal,
            prm3: 0.0,
            axis: (p.origin, p.axis_dir, perp),
            typ: Curve2dType::Parabola,
            abs_to_object: (p.origin, p.axis_dir, perp),
        }
    }

    /// OCCT IntCurve_IConicTool(gp_Hypr2d) (L127-137).
    pub fn new_hyperbola(h: &Hyperbola2d) -> Self {
        let minor_dir = DVec2::new(-h.major_dir.y, h.major_dir.x);
        IConicTool {
            prm1: h.semi_major,
            prm2: h.semi_minor,
            prm3: 0.0,
            axis: (h.center, h.major_dir, minor_dir),
            typ: Curve2dType::Hyperbola,
            abs_to_object: (h.center, h.major_dir, minor_dir),
        }
    }

    /// OCCT IntCurve_IConicTool::Value (L140-159).
    pub fn value(&self, x: f64) -> DVec2 {
        let (loc, xdir, ydir) = self.axis;
        match self.typ {
            Curve2dType::Line => elclib2d::line_value(loc, xdir, x),
            Curve2dType::Ellipse => elclib2d::ellipse_value(loc, xdir, ydir, self.prm1, self.prm2, x),
            Curve2dType::Circle => elclib2d::circle_value(loc, xdir, ydir, self.prm1, x),
            Curve2dType::Parabola => elclib2d::parabola_value(loc, xdir, ydir, self.prm1, x),
            Curve2dType::Hyperbola => elclib2d::hyperbola_value(loc, xdir, ydir, self.prm1, self.prm2, x),
            _ => DVec2::ZERO,
        }
    }

    /// OCCT IntCurve_IConicTool::D1 (L162-186).
    pub fn d1(&self, x: f64) -> (DVec2, DVec2) {
        let (loc, xdir, ydir) = self.axis;
        match self.typ {
            Curve2dType::Line => elclib2d::line_d1(loc, xdir, x),
            Curve2dType::Ellipse => elclib2d::ellipse_d1(loc, xdir, ydir, self.prm1, self.prm2, x),
            Curve2dType::Circle => elclib2d::circle_d1(loc, xdir, ydir, self.prm1, x),
            Curve2dType::Parabola => elclib2d::parabola_d1(loc, xdir, ydir, self.prm1, x),
            Curve2dType::Hyperbola => elclib2d::hyperbola_d1(loc, xdir, ydir, self.prm1, self.prm2, x),
            _ => (DVec2::ZERO, DVec2::ZERO),
        }
    }

    /// OCCT IntCurve_IConicTool::D2 (L189-214).
    pub fn d2(&self, x: f64) -> (DVec2, DVec2, DVec2) {
        let (loc, xdir, ydir) = self.axis;
        match self.typ {
            Curve2dType::Line => {
                let (p, t) = elclib2d::line_d1(loc, xdir, x);
                (p, t, DVec2::ZERO)
            }
            Curve2dType::Ellipse => elclib2d::ellipse_d2(loc, xdir, ydir, self.prm1, self.prm2, x),
            Curve2dType::Circle => elclib2d::circle_d2(loc, xdir, ydir, self.prm1, x),
            Curve2dType::Parabola => elclib2d::parabola_d2(loc, xdir, ydir, self.prm1, x),
            Curve2dType::Hyperbola => elclib2d::hyperbola_d2(loc, xdir, ydir, self.prm1, self.prm2, x),
            _ => (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO),
        }
    }

    /// OCCT IntCurve_IConicTool::Distance (L220-279).
    pub fn distance(&self, p: DVec2) -> f64 {
        match self.typ {
            Curve2dType::Line => self.prm1 * p.x + self.prm2 * p.y + self.prm3,
            Curve2dType::Circle => {
                let dx = self.prm2 - p.x;
                let dy = self.prm3 - p.y;
                (dx * dx + dy * dy).sqrt() - self.prm1
            }
            Curve2dType::Ellipse => {
                let po = self.object_coords(p);
                let x = po.x;
                let y = po.y * (self.prm1 / self.prm2);
                (x * x + y * y).sqrt() - self.prm1
            }
            Curve2dType::Parabola => {
                let po = self.object_coords(p);
                po.y * po.y - self.prm2 * po.x
            }
            Curve2dType::Hyperbola => {
                let po = self.object_coords(p);
                if po.x > 0.0 {
                    (po.x * po.x) / (self.prm1 * self.prm1) - (po.y * po.y) / (self.prm2 * self.prm2) - 1.0
                } else {
                    (-po.x * po.x) / (self.prm1 * self.prm1) - (po.y * po.y) / (self.prm2 * self.prm2) - 1.0
                }
            }
            _ => 0.0,
        }
    }

    /// OCCT IntCurve_IConicTool::GradDistance (L281-370).
    pub fn grad_distance(&self, p: DVec2) -> DVec2 {
        match self.typ {
            Curve2dType::Line => DVec2::new(self.prm1, self.prm2),
            Curve2dType::Circle => {
                let po = self.object_coords(p);
                let temp1 = (po.y * po.y + po.x * po.x).sqrt();
                let mut grad = DVec2::ZERO;
                if temp1 != 0.0 {
                    grad = DVec2::new(po.x / temp1, po.y / temp1);
                }
                self.absolute_from_object(grad)
            }
            Curve2dType::Ellipse => {
                let po = self.object_coords(p);
                let x = po.x;
                let y = po.y * (self.prm1 / self.prm2);
                let temp1 = (y * y + x * x).sqrt();
                let mut grad = DVec2::ZERO;
                if temp1 != 0.0 {
                    grad = DVec2::new(x / temp1, (y * (self.prm1 / self.prm2)) / temp1);
                }
                self.absolute_from_object(grad)
            }
            Curve2dType::Parabola => {
                let po = self.object_coords(p);
                self.absolute_from_object(DVec2::new(-self.prm2, po.y + po.y))
            }
            Curve2dType::Hyperbola => {
                let po = self.object_coords(p);
                self.absolute_from_object(DVec2::new(
                    2.0 * po.x.abs() / (self.prm1 * self.prm1),
                    -2.0 * po.y / (self.prm2 * self.prm2),
                ))
            }
            _ => DVec2::ZERO,
        }
    }

    /// OCCT IntCurve_IConicTool::FindParameter (L372-414).
    pub fn find_parameter(&self, p: DVec2) -> f64 {
        let (loc, xdir, ydir) = self.axis;
        let mut param = 0.0;
        match self.typ {
            Curve2dType::Line => param = elclib2d::line_parameter(loc, xdir, p),
            Curve2dType::Circle => {
                param = elclib2d::circle_parameter(loc, xdir, ydir, p);
                if param < 0.0 {
                    param += PI2;
                }
            }
            Curve2dType::Ellipse => {
                param = elclib2d::ellipse_parameter(loc, xdir, ydir, self.prm1, self.prm2, p);
                if param < 0.0 {
                    param += PI2;
                }
            }
            Curve2dType::Parabola => param = elclib2d::parabola_parameter(loc, ydir, p),
            Curve2dType::Hyperbola => {
                param = elclib2d::hyperbola_parameter(loc, ydir, self.prm2, p)
            }
            _ => {}
        }
        param
    }
}

/// OCCT Extrema_GCurveLocator (Extrema_GCurveLocator.hxx L84-129) — Locate:
/// among a set of samples {C(ui)}, find the point closest to P.
fn locate_on_curve(
    p: DVec2,
    c: &dyn Curve2dAdaptor,
    nb_u: i32,
    u_min: f64,
    u_sup: f64,
) -> (f64, DVec2) {
    assert!(nb_u >= 2, "Standard_OutOfRange");
    let a_uinf = geom2d_curve_tool::first_parameter(c);
    let a_ulast = geom2d_curve_tool::last_parameter(c);
    let a_u1 = a_uinf.min(a_ulast);
    let a_u2 = a_uinf.max(a_ulast);
    let mut a_u11 = u_min.min(u_sup);
    let mut a_u12 = u_min.max(u_sup);
    if a_u11 < a_u1 - f64::EPSILON {
        a_u11 = a_u1;
    }
    if a_u12 > a_u2 + f64::EPSILON {
        a_u12 = a_u2;
    }
    let mut a_u = a_u11;
    let a_pas_u = (a_u12 - a_u) / (nb_u - 1) as f64;
    let mut a_dist2_min = f64::MAX;
    let mut a_u_min = 0.0;
    let mut a_pnt_min = DVec2::ZERO;
    for _ in 1..nb_u {
        let a_pt = geom2d_curve_tool::value(c, a_u);
        let a_dist2 = a_pt.distance_squared(p);
        if a_dist2 < a_dist2_min {
            a_dist2_min = a_dist2;
            a_u_min = a_u;
            a_pnt_min = a_pt;
        }
        a_u += a_pas_u;
    }
    (a_u_min, a_pnt_min)
}

/// OCCT Extrema_GFuncExtPC (Extrema_GFuncExtPC.hxx) — the function
/// F(u) = (C(u)-P)·D1c/|D1c| with its derivative, used by math_FunctionRoot.
struct GFuncExtPC<'a> {
    p: DVec2,
    c: &'a dyn Curve2dAdaptor,
    u: f64,
    pc: DVec2,
    d1f: f64,
    sq_dist: Vec<f64>,
    is_min: Vec<i32>,
    points: Vec<(f64, DVec2)>,
    pinit: bool,
    cinit: bool,
    _d1_init: bool,
    tol: f64,
    max_deriv_order: i32,
    uinfium: f64,
    usupremum: f64,
}

const G_FUNC_MIN_TOL: f64 = 1e-20;
const G_FUNC_MIN_STEP: f64 = 1e-7;
const G_FUNC_MAX_ORDER: i32 = 3;

impl<'a> GFuncExtPC<'a> {
    /// OCCT Extrema_GFuncExtPC(P, C) constructor (L76-101).
    fn new(p: DVec2, c: &'a dyn Curve2dAdaptor) -> Self {
        let mut f = GFuncExtPC {
            p,
            c,
            u: 0.0,
            pc: DVec2::ZERO,
            d1f: 0.0,
            sq_dist: Vec::new(),
            is_min: Vec::new(),
            points: Vec::new(),
            pinit: true,
            cinit: true,
            _d1_init: false,
            tol: G_FUNC_MIN_TOL,
            max_deriv_order: 0,
            uinfium: geom2d_curve_tool::first_parameter(c),
            usupremum: geom2d_curve_tool::last_parameter(c),
        };
        f.sub_interval_initialize();
        match geom2d_curve_tool::get_type(c) {
            Curve2dType::BezierCurve
            | Curve2dType::BSplineCurve
            | Curve2dType::OffsetCurve
            | Curve2dType::OtherCurve => {
                f.max_deriv_order = G_FUNC_MAX_ORDER;
                f.tol = f.search_of_tolerance();
            }
            _ => {
                f.max_deriv_order = 0;
                f.tol = G_FUNC_MIN_TOL;
            }
        }
        f
    }

    /// OCCT SubIntervalInitialize (L420-424).
    fn sub_interval_initialize(&mut self) {
        self.uinfium = geom2d_curve_tool::first_parameter(self.c);
        self.usupremum = geom2d_curve_tool::last_parameter(self.c);
    }

    /// OCCT SearchOfTolerance (L428-457).
    fn search_of_tolerance(&mut self) -> f64 {
        let n_point = 10;
        let a_step = (self.usupremum - self.uinfium) / n_point as f64;
        let mut a_num = 0;
        let mut a_max = f64::NEG_INFINITY;
        loop {
            let mut u = self.uinfium + a_num as f64 * a_step;
            if u > self.usupremum {
                u = self.usupremum;
            }
            let (_, v_der) = geom2d_curve_tool::d1(self.c, u);
            if !(v_der.x.is_infinite() || v_der.y.is_infinite()) {
                let vm = v_der.length();
                if vm > a_max {
                    a_max = vm;
                }
            }
            a_num += 1;
            if a_num >= n_point + 1 {
                break;
            }
        }
        (a_max * 1.0e-12).max(G_FUNC_MIN_TOL)
    }

    /// OCCT SetPoint (L133-140).
    fn set_point(&mut self, p: DVec2) {
        self.p = p;
        self.pinit = true;
        self.sq_dist.clear();
        self.is_min.clear();
        self.points.clear();
    }

    /// OCCT Value(U, F) (L146-253).
    fn func_value(&mut self, u: f64) -> Option<f64> {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch: No init");
        }
        self.u = u;
        let (pc, mut d1c) = geom2d_curve_tool::d1(self.c, u);
        self.pc = pc;
        if d1c.x.is_infinite() || d1c.y.is_infinite() {
            return None;
        }
        let mut ndu = d1c.length();
        if self.max_deriv_order != 0 {
            if ndu <= self.tol {
                // Singular case: the derivative is approximated by a Taylor series.
                let division_factor = 1.0e-3;
                let du = if self.usupremum >= f64::MAX || self.uinfium <= f64::MIN {
                    0.0
                } else {
                    self.usupremum - self.uinfium
                };
                let a_delta = (du * division_factor).max(G_FUNC_MIN_STEP);
                let mut n = 1;
                let mut is_deriv_found = false;
                loop {
                    let v = geom2d_curve_tool::dn(self.c, self.u, n + 1);
                    ndu = v.length();
                    is_deriv_found = ndu > self.tol;
                    if is_deriv_found || n >= self.max_deriv_order {
                        if is_deriv_found {
                            let u2 = if self.u - self.uinfium < a_delta {
                                self.u + a_delta
                            } else {
                                self.u - a_delta
                            };
                            let p1 = geom2d_curve_tool::d0(self.c, self.u.min(u2));
                            let p2 = geom2d_curve_tool::d0(self.c, self.u.max(u2));
                            let v1 = p2 - p1;
                            if v.dot(v1) < 0.0 {
                                d1c = -v;
                            } else {
                                d1c = v;
                            }
                        }
                        break;
                    }
                    n += 1;
                }
                if !is_deriv_found {
                    // Derivative approximated by three points.
                    let (p1, p2, p3) = if self.u - self.uinfium < 2.0 * a_delta {
                        (
                            geom2d_curve_tool::d0(self.c, self.u),
                            geom2d_curve_tool::d0(self.c, self.u + a_delta),
                            geom2d_curve_tool::d0(self.c, self.u + 2.0 * a_delta),
                        )
                    } else {
                        (
                            geom2d_curve_tool::d0(self.c, self.u - 2.0 * a_delta),
                            geom2d_curve_tool::d0(self.c, self.u - a_delta),
                            geom2d_curve_tool::d0(self.c, self.u),
                        )
                    };
                    d1c = if self.u - self.uinfium < 2.0 * a_delta {
                        -3.0 * p1 + 4.0 * p2 - p3
                    } else {
                        p1 - 4.0 * p2 + 3.0 * p3
                    };
                }
                ndu = d1c.length();
            }
        }
        if ndu <= G_FUNC_MIN_TOL {
            return None;
        }
        let ppc = self.p - self.pc;
        Some(ppc.dot(d1c) / ndu)
    }

    /// OCCT Values(U, F, DF) (L274-353).
    fn func_values(&mut self, u: f64) -> Option<(f64, f64)> {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch: No init");
        }
        let pc_old = self.pc;
        let p_old = self.p;
        let the_f = self.func_value(u)?;
        self.u = u;
        self.pc = pc_old;
        self.p = p_old;

        let (pc2, d1c, d2c) = geom2d_curve_tool::d2(self.c, u);
        self.pc = pc2;
        let ndu = d1c.length();
        let the_df;
        if ndu <= self.tol {
            // Singular case: the derivative is approximated by three points.
            let division_factor = 0.01;
            let du = if self.usupremum >= f64::MAX || self.uinfium <= f64::MIN {
                0.0
            } else {
                self.usupremum - self.uinfium
            };
            let a_delta = (du * division_factor).max(G_FUNC_MIN_STEP);
            if self.u - self.uinfium < 2.0 * a_delta {
                let u2 = self.u + a_delta;
                let u3 = self.u + a_delta * 2.0;
                let f2 = self.func_value(u2)?;
                let f3 = self.func_value(u3)?;
                the_df = (-3.0 * the_f + 4.0 * f2 - f3) / (2.0 * a_delta);
            } else {
                let u1 = self.u - a_delta * 2.0;
                let u2 = self.u - a_delta;
                let f2 = self.func_value(u2)?;
                let f1 = self.func_value(u1)?;
                the_df = (f1 - 4.0 * f2 + 3.0 * the_f) / (2.0 * a_delta);
            }
            self.u = u;
            self.pc = pc_old;
            self.p = p_old;
        } else {
            let ppc = self.p - self.pc;
            the_df = ndu + (ppc.dot(d2c) / ndu) - the_f * (d1c.dot(d2c)) / (ndu * ndu);
        }
        self.d1f = the_df;
        self._d1_init = true;
        Some((the_f, the_df))
    }

    /// OCCT GetStateNumber (L357-379).
    fn get_state(&mut self) -> i32 {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        self.sq_dist.push(self.pc.distance_squared(self.p));
        self._d1_init = true;
        let _ = self.func_values(self.u);
        let int_val = if self.d1f > 0.0 { 1 } else { 0 };
        self.is_min.push(int_val);
        self.points.push((self.u, self.pc));
        0
    }

    fn is_min(&self, n: usize) -> bool {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        self.is_min[n - 1] == 1
    }

    fn point(&self, n: usize) -> (f64, DVec2) {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        self.points[n - 1]
    }
}

impl FunctionValue for GFuncExtPC<'_> {
    fn value(&mut self, x: f64) -> Option<f64> {
        self.func_value(x)
    }
}

impl FunctionWithDerivative for GFuncExtPC<'_> {
    fn derivative(&mut self, x: f64) -> Option<f64> {
        self.func_values(x).map(|(_, d)| d)
    }
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        self.func_values(x)
    }
    fn get_state_number(&mut self) -> i32 {
        self.get_state()
    }
}

/// OCCT math_FunctionRoot(F, Guess, Tolerance, A, B, NbIterations) — the
/// bounded Newton root via math_FunctionSetRoot (math_FunctionRoot.cxx
/// L95-119).
fn math_function_root_1d(
    f: &mut dyn FunctionWithDerivative,
    guess: f64,
    tolerance: f64,
    a: f64,
    b: f64,
    nb_iterations: i32,
) -> Option<(f64, f64)> {
    struct Adapter<'a>(&'a mut dyn FunctionWithDerivative);
    impl FunctionSetWithDerivatives for Adapter<'_> {
        fn nb_variables(&self) -> usize {
            1
        }
        fn nb_equations(&self) -> usize {
            1
        }
        fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
            match self.0.value(x[0]) {
                Some(v) => {
                    f[0] = v;
                    true
                }
                None => false,
            }
        }
        fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
            match self.0.derivative(x[0]) {
                Some(d) => {
                    df[0][0] = d;
                    true
                }
                None => false,
            }
        }
        fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
            match self.0.values(x[0]) {
                Some((v, d)) => {
                    f[0] = v;
                    df[0][0] = d;
                    true
                }
                None => false,
            }
        }
    }
    let mut adapter = Adapter(f);
    let mut sol = FunctionSetRoot::new(&adapter, &[tolerance], nb_iterations);
    sol.perform(&mut adapter, &[guess], &[a], &[b], false);
    if !sol.is_done() {
        return None;
    }
    // OCCT math_FunctionRoot (L85-87 / L111-113): F.GetStateNumber() records
    // the extremum point after a successful solve.
    f.get_state_number();
    Some((sol.root()[0], sol.derivative()))
}

/// OCCT Extrema_GenLocateExtPC (Extrema_GenLocateExtPC.hxx L108-126) — local
/// extremum of the distance from a point to a curve near a seed parameter.
/// Returns the parameter of the found extremum.
fn gen_locate_ext_pc(
    p: DVec2,
    c: &dyn Curve2dAdaptor,
    u0: f64,
    u_min: f64,
    u_sup: f64,
    tol_u: f64,
) -> Option<(f64, bool)> {
    let mut f = GFuncExtPC::new(p, c);
    f.set_point(p);
    let root = math_function_root_1d(&mut f, u0, tol_u, u_min, u_sup, 100)?;
    let (uu, _) = (root.0, root.1);
    let (param, _pt) = f.point(1);
    let uu = param;
    match f.func_value(uu) {
        Some(ff) => {
            if ff.abs() >= 1.0e-7 {
                None
            } else {
                Some((param, f.is_min(1)))
            }
        }
        None => None,
    }
}

/// OCCT Geom2dInt_TheProjPCurOfGInter (Geom2dInt_TheProjPCurOfGInter.cxx) —
/// projection of a point onto the parametric curve.
pub mod proj_p_cur_of_g_inter {
    use super::*;

    /// OCCT FindParameter(C, P, LowParameter, HighParameter, Tol) (L27-65).
    pub fn find_parameter_bounded(
        c: &dyn Curve2dAdaptor,
        p: DVec2,
        low_parameter: f64,
        high_parameter: f64,
        _tol: f64,
    ) -> f64 {
        let nb_pts = geom2d_curve_tool::nb_samples(c);
        let the_eps_x = geom2d_curve_tool::eps_x(c);
        let (a_u_min, a_pnt_min) = locate_on_curve(p, c, nb_pts, low_parameter, high_parameter);
        let _ = a_pnt_min;
        let default_param = a_u_min;
        let locate = gen_locate_ext_pc(p, c, default_param, low_parameter, high_parameter, the_eps_x);
        let the_param;
        match locate {
            None => {
                the_param = default_param;
            }
            Some((param, is_min)) => {
                if !is_min {
                    the_param = default_param;
                } else {
                    the_param = param;
                }
            }
        }
        the_param
    }

    /// OCCT FindParameter(C, P, Tol) (L67-79).
    pub fn find_parameter_unbounded(c: &dyn Curve2dAdaptor, p: DVec2, tol: f64) -> f64 {
        find_parameter_bounded(
            c,
            p,
            geom2d_curve_tool::first_parameter(c),
            geom2d_curve_tool::last_parameter(c),
            tol,
        )
    }
}


// ---------------------------------------------------------------------------
// OCCT Geom2dInt instantiations of the generic IntImpParGen_Intersector and
// IntCurve_IntConicCurveGen engines (the _0.cxx template-argument sets).
// ---------------------------------------------------------------------------

/// OCCT Geom2dInt_Geom2dCurveTool (hxx) — the parametric-curve tool of the
/// GInter instantiations, served over `dyn Curve2dAdaptor`. It carries the
/// ParTool role of IntImpParGen_Intersector and — as in OCCT, where
/// HLRBRep_CurveTool plays both roles for the CInter instantiation — also
/// the ThePCurveTool role of IntCurve_UserIntConicCurveGen for
/// Adaptor2d-based callers.
pub struct Geom2dCurveTool;

impl<'a> crate::geomalgo::int_imp_par_gen::ParTool<dyn Curve2dAdaptor + 'a> for Geom2dCurveTool {
    fn nb_samples_uv(c: &(dyn Curve2dAdaptor + 'a), u1: f64, u2: f64) -> i32 {
        geom2d_curve_tool::nb_samples_2(c, u1, u2)
    }
    fn eps_x(c: &(dyn Curve2dAdaptor + 'a)) -> f64 {
        geom2d_curve_tool::eps_x(c)
    }
    fn value(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> DVec2 {
        geom2d_curve_tool::value(c, u)
    }
    fn d1(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> (DVec2, DVec2) {
        geom2d_curve_tool::d1(c, u)
    }
    fn d2(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> (DVec2, DVec2, DVec2) {
        geom2d_curve_tool::d2(c, u)
    }
}

impl<'a> crate::geomalgo::user_int_conic_curve_gen::PCurveTool<dyn Curve2dAdaptor + 'a>
    for Geom2dCurveTool
{
    fn nb_intervals(c: &(dyn Curve2dAdaptor + 'a)) -> i32 {
        geom2d_curve_tool::nb_intervals(c)
    }
    fn intervals(c: &(dyn Curve2dAdaptor + 'a), tab: &mut [f64]) {
        geom2d_curve_tool::intervals(c, tab)
    }
    fn get_interval(_c: &(dyn Curve2dAdaptor + 'a), index: usize, tab: &[f64]) -> (f64, f64) {
        geom2d_curve_tool::get_interval(tab, index)
    }
    fn first_parameter(c: &(dyn Curve2dAdaptor + 'a)) -> f64 {
        geom2d_curve_tool::first_parameter(c)
    }
    fn last_parameter(c: &(dyn Curve2dAdaptor + 'a)) -> f64 {
        geom2d_curve_tool::last_parameter(c)
    }
    fn value(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> DVec2 {
        geom2d_curve_tool::value(c, u)
    }
    fn get_type(c: &(dyn Curve2dAdaptor + 'a)) -> Curve2dType {
        geom2d_curve_tool::get_type(c)
    }
    fn line(c: &(dyn Curve2dAdaptor + 'a)) -> Line2d {
        geom2d_curve_tool::line(c)
    }
    fn circle(c: &(dyn Curve2dAdaptor + 'a)) -> Circle2d {
        geom2d_curve_tool::circle(c)
    }
    fn ellipse(c: &(dyn Curve2dAdaptor + 'a)) -> Ellipse2d {
        geom2d_curve_tool::ellipse(c)
    }
    fn parabola(c: &(dyn Curve2dAdaptor + 'a)) -> Parabola2d {
        geom2d_curve_tool::parabola(c)
    }
    fn hyperbola(c: &(dyn Curve2dAdaptor + 'a)) -> Hyperbola2d {
        geom2d_curve_tool::hyperbola(c)
    }
}

/// OCCT Geom2dInt_TheProjPCurOfGInter (cxx) — the point projection on the
/// parametric curve (the ProjectOnPCurveTool role of the GInter
/// instantiations).
pub struct TheProjPCurOfGInter;

impl<'a> crate::geomalgo::int_imp_par_gen::ProjectOnPCurveTool<dyn Curve2dAdaptor + 'a>
    for TheProjPCurOfGInter
{
    fn find_parameter(c: &(dyn Curve2dAdaptor + 'a), p: DVec2, tol: f64) -> f64 {
        proj_p_cur_of_g_inter::find_parameter_unbounded(c, p, tol)
    }
    fn find_parameter_between(
        c: &(dyn Curve2dAdaptor + 'a),
        p: DVec2,
        low: f64,
        high: f64,
        tol: f64,
    ) -> f64 {
        proj_p_cur_of_g_inter::find_parameter_bounded(c, p, low, high, tol)
    }
}

/// OCCT Geom2dInt_TheIntersectorOfTheIntConicCurveOfGInter — the
/// IntImpParGen_Intersector instantiation over the Geom2dInt tools
/// (Geom2dInt_TheIntersectorOfTheIntConicCurveOfGInter_0.cxx: ImpTool =
/// IntCurve_IConicTool, ParTool = Geom2dInt_Geom2dCurveTool,
/// ProjectOnPCurveTool = Geom2dInt_TheProjPCurOfGInter, ParCurve =
/// Adaptor2d_Curve2d).
#[derive(Debug, Clone)]
pub struct TheIntersectorOfTheIntConicCurveOfGInter {
    pub base: IntersectionBase,
}

impl TheIntersectorOfTheIntConicCurveOfGInter {
    pub fn new() -> Self {
        TheIntersectorOfTheIntConicCurveOfGInter {
            base: IntersectionBase::new(),
        }
    }

    /// OCCT IntImpParGen_Intersector::Perform (IntImpParGen_Intersector.gxx
    /// L245-779) via the generic engine.
    pub fn perform(
        &mut self,
        imp_tool: &IConicTool,
        imp_domain: &Res2dDomain,
        par_curve: &dyn Curve2dAdaptor,
        par_domain: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        fn engine<'a>(
        ) -> super::int_imp_par_gen::Intersector<dyn Curve2dAdaptor + 'a, Geom2dCurveTool, TheProjPCurOfGInter>
        {
            super::int_imp_par_gen::Intersector::new()
        }
        let mut myintersection = engine();
        myintersection
            .base
            .set_reversed_parameters(self.base.reversed_parameters());
        myintersection.perform(imp_tool, imp_domain, par_curve, par_domain, tol_conf, tol);
        self.base = myintersection.base;
    }
}

impl Default for TheIntersectorOfTheIntConicCurveOfGInter {
    fn default() -> Self {
        TheIntersectorOfTheIntConicCurveOfGInter::new()
    }
}

/// OCCT Geom2dInt_TheIntConicCurveOfGInter — the IntCurve_IntConicCurveGen
/// instantiation over the same tools
/// (Geom2dInt_TheIntConicCurveOfGInter_0.cxx).
#[derive(Debug, Clone)]
pub struct TheIntConicCurveOfGInter {
    pub base: IntersectionBase,
}

impl TheIntConicCurveOfGInter {
    /// OCCT IntCurve_IntConicCurveGen() (IntCurve_IntConicCurveGen.lxx L30).
    pub fn bare() -> Self {
        TheIntConicCurveOfGInter {
            base: IntersectionBase::new(),
        }
    }

    /// OCCT IntCurve_IntConicCurveGen::Perform(const gp_Lin2d& L, ...)
    /// (IntCurve_IntConicCurveGen.lxx L44-53).
    pub fn perform_line(
        &mut self,
        line: &Line2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.perform_imp(&IConicTool::new_line(line), d1, pcurve, d2, tol_conf, tol);
    }

    /// OCCT IntCurve_IntConicCurveGen::Perform(const gp_Circ2d& C, ...)
    /// (IntCurve_IntConicCurveGen.lxx L56-74).
    pub fn perform_circle(
        &mut self,
        c: &Circle2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            self.perform_imp(&IConicTool::new_circle(c), &d, pcurve, d2, tol_conf, tol);
        } else {
            self.perform_imp(&IConicTool::new_circle(c), d1, pcurve, d2, tol_conf, tol);
        }
    }

    /// OCCT IntCurve_IntConicCurveGen::Perform(const gp_Elips2d& E, ...)
    /// (IntCurve_IntConicCurveGen.lxx L77-94).
    pub fn perform_ellipse(
        &mut self,
        e: &Ellipse2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            self.perform_imp(&IConicTool::new_ellipse(e), &d, pcurve, d2, tol_conf, tol);
        } else {
            self.perform_imp(&IConicTool::new_ellipse(e), d1, pcurve, d2, tol_conf, tol);
        }
    }

    /// OCCT IntCurve_IntConicCurveGen::Perform(const gp_Parab2d& Prb, ...)
    /// (IntCurve_IntConicCurveGen.lxx L97-105).
    pub fn perform_parabola(
        &mut self,
        p: &Parabola2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.perform_imp(&IConicTool::new_parabola(p), d1, pcurve, d2, tol_conf, tol);
    }

    /// OCCT IntCurve_IntConicCurveGen::Perform(const gp_Hypr2d& H, ...)
    /// (IntCurve_IntConicCurveGen.lxx L108-116).
    pub fn perform_hyperbola(
        &mut self,
        h: &Hyperbola2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.perform_imp(&IConicTool::new_hyperbola(h), d1, pcurve, d2, tol_conf, tol);
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Lin2d&, ...) (lxx L33-42).
    pub fn new_line(
        line: &Line2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = TheIntConicCurveOfGInter::bare();
        r.perform_imp(&IConicTool::new_line(line), d1, pcurve, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Circ2d&, ...) (gxx L24-42).
    pub fn new_circle(
        c: &Circle2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = TheIntConicCurveOfGInter::bare();
        let tool = IConicTool::new_circle(c);
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            r.perform_imp(&tool, &d, pcurve, d2, tol_conf, tol);
        } else {
            r.perform_imp(&tool, d1, pcurve, d2, tol_conf, tol);
        }
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Elips2d&, ...) (gxx L46-64).
    pub fn new_ellipse(
        e: &Ellipse2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = TheIntConicCurveOfGInter::bare();
        let tool = IConicTool::new_ellipse(e);
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            r.perform_imp(&tool, &d, pcurve, d2, tol_conf, tol);
        } else {
            r.perform_imp(&tool, d1, pcurve, d2, tol_conf, tol);
        }
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Parab2d&, ...) (gxx L68-77).
    pub fn new_parabola(
        p: &Parabola2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = TheIntConicCurveOfGInter::bare();
        r.perform_imp(&IConicTool::new_parabola(p), d1, pcurve, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Hypr2d&, ...) (gxx L81-90).
    pub fn new_hyperbola(
        h: &Hyperbola2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = TheIntConicCurveOfGInter::bare();
        r.perform_imp(&IConicTool::new_hyperbola(h), d1, pcurve, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntConicCurveGen::Perform(const IntCurve_IConicTool&,
    /// ...) (IntCurve_IntConicCurveGen.lxx L119-130).
    fn perform_imp(
        &mut self,
        imp_tool: &IConicTool,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut myintersection = TheIntersectorOfTheIntConicCurveOfGInter::new();
        myintersection
            .base
            .set_reversed_parameters(self.base.reversed_parameters());
        myintersection.perform(imp_tool, d1, pcurve, d2, tol_conf, tol);
        self.base.set_values(&myintersection.base);
    }
}

impl Default for TheIntConicCurveOfGInter {
    fn default() -> Self {
        TheIntConicCurveOfGInter::bare()
    }
}

// ---------------------------------------------------------------------------
// OCCT Geom2dInt instantiation of IntCurve_IntPolyPolyGen.
// ---------------------------------------------------------------------------

/// OCCT Geom2dInt tool triple over `dyn Curve2dAdaptor`
/// (Geom2dInt_Geom2dCurveTool + Geom2dInt_TheProjPCurOfGInter, bundled for
/// the IntCurve_IntPolyPolyGen template parameters).  The lifetime
/// parameter keeps the curve references lifetime-elastic like the C++
/// `const Adaptor2d_Curve2d&` signatures.
pub struct GInterPolyTool<'a>(std::marker::PhantomData<&'a ()>);

impl<'a> crate::geomalgo::int_curve_generics::ProjPCurveTool for GInterPolyTool<'a> {
    type Curve = dyn Curve2dAdaptor + 'a;
    fn value(c: &Self::Curve, u: f64) -> DVec2 {
        geom2d_curve_tool::value(c, u)
    }
    fn d1(c: &Self::Curve, u: f64) -> (DVec2, DVec2) {
        geom2d_curve_tool::d1(c, u)
    }
    fn eps_x(c: &Self::Curve) -> f64 {
        geom2d_curve_tool::eps_x(c)
    }
}

impl<'a> crate::geomalgo::int_imp_par_gen::ParTool<dyn Curve2dAdaptor + 'a> for GInterPolyTool<'a> {
    fn nb_samples_uv(c: &(dyn Curve2dAdaptor + 'a), u1: f64, u2: f64) -> i32 {
        geom2d_curve_tool::nb_samples_2(c, u1, u2)
    }
    fn eps_x(c: &(dyn Curve2dAdaptor + 'a)) -> f64 {
        geom2d_curve_tool::eps_x(c)
    }
    fn value(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> DVec2 {
        geom2d_curve_tool::value(c, u)
    }
    fn d1(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> (DVec2, DVec2) {
        geom2d_curve_tool::d1(c, u)
    }
    fn d2(c: &(dyn Curve2dAdaptor + 'a), u: f64) -> (DVec2, DVec2, DVec2) {
        geom2d_curve_tool::d2(c, u)
    }
}

impl<'a> crate::geomalgo::int_imp_par_gen::ProjectOnPCurveTool<dyn Curve2dAdaptor + 'a>
    for GInterPolyTool<'a>
{
    fn find_parameter(c: &(dyn Curve2dAdaptor + 'a), p: DVec2, tol: f64) -> f64 {
        proj_p_cur_of_g_inter::find_parameter_unbounded(c, p, tol)
    }
    fn find_parameter_between(
        c: &(dyn Curve2dAdaptor + 'a),
        p: DVec2,
        low: f64,
        high: f64,
        tol: f64,
    ) -> f64 {
        proj_p_cur_of_g_inter::find_parameter_bounded(c, p, low, high, tol)
    }
}

/// OCCT Geom2dInt_TheIntPCurvePCurveOfGInter — the IntCurve_IntPolyPolyGen
/// instantiation over the Geom2dInt tools
/// (Geom2dInt_TheIntPCurvePCurveOfGInter_0.cxx).
#[derive(Debug, Clone)]
pub struct TheIntPCurvePCurveOfGInter {
    pub base: IntersectionBase,
    pub domain_on_curve1: Res2dDomain,
    pub domain_on_curve2: Res2dDomain,
    min_pnt_nb: i32,
}

impl TheIntPCurvePCurveOfGInter {
    /// OCCT IntCurve_IntPolyPolyGen() (gxx L84-89).
    pub fn new() -> Self {
        TheIntPCurvePCurveOfGInter {
            base: IntersectionBase::new(),
            domain_on_curve1: Res2dDomain::infinite(),
            domain_on_curve2: Res2dDomain::infinite(),
            min_pnt_nb: 20,
        }
    }

    /// OCCT GetMinNbSamples() (gxx L1787-1790).
    pub fn get_min_nb_samples(&self) -> i32 {
        self.min_pnt_nb
    }

    /// OCCT SetMinNbSamples(theMinNbSamples) (gxx L1794-1797).
    pub fn set_min_nb_samples(&mut self, the_min_nb_samples: i32) {
        self.min_pnt_nb = the_min_nb_samples;
    }

    /// OCCT Perform(C1, D1, C2, D2, TheTolConf, TheTol) (gxx L93-292).
    pub fn perform(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        fn engine<'a>(
        ) -> super::int_poly_poly_gen::IntPolyPolyGen<dyn Curve2dAdaptor + 'a, GInterPolyTool<'a>>
        {
            super::int_poly_poly_gen::IntPolyPolyGen::new()
        }
        let mut e = engine();
        e.set_min_nb_samples(self.min_pnt_nb);
        e.perform(c1, d1, c2, d2, the_tol_conf, the_tol);
        self.base = e.base;
        self.domain_on_curve1 = e.domain_on_curve1;
        self.domain_on_curve2 = e.domain_on_curve2;
    }

    /// OCCT Perform(C1, D1, TheTolConf, TheTol) (gxx L297-386) — the
    /// auto-intersection form.
    pub fn perform_cd(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        self.perform_auto(c1, d1, the_tol_conf, the_tol);
    }

    /// OCCT Perform(C1, D1, TheTolConf, TheTol) (gxx L297-386) — the
    /// auto-intersection form.
    pub fn perform_auto(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        fn engine<'a>(
        ) -> super::int_poly_poly_gen::IntPolyPolyGen<dyn Curve2dAdaptor + 'a, GInterPolyTool<'a>>
        {
            super::int_poly_poly_gen::IntPolyPolyGen::new()
        }
        let mut e = engine();
        e.set_min_nb_samples(self.min_pnt_nb);
        e.perform_auto(c1, d1, the_tol_conf, the_tol);
        self.base = e.base;
        self.domain_on_curve1 = e.domain_on_curve1;
        self.domain_on_curve2 = e.domain_on_curve2;
    }
}

impl Default for TheIntPCurvePCurveOfGInter {
    fn default() -> Self {
        TheIntPCurvePCurveOfGInter::new()
    }
}

// ---------------------------------------------------------------------------
// OCCT Geom2dInt_GInter — the IntCurve_IntCurveCurveGen instantiation over
// the Adaptor2d_Curve2d curves (Geom2dInt_GInter_0.cxx: TheCurve =
// Adaptor2d_Curve2d, TheCurveTool = Geom2dInt_Geom2dCurveTool,
// IntCurve_TheIntConicCurve = Geom2dInt_TheIntConicCurveOfGInter,
// IntCurve_TheIntPCurvePCurve = Geom2dInt_TheIntPCurvePCurveOfGInter).
// ---------------------------------------------------------------------------

impl crate::geomalgo::int_curve_curve_gen::IntConicCurveMember<dyn Curve2dAdaptor>
    for TheIntConicCurveOfGInter
{
    fn base(&self) -> &IntersectionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut IntersectionBase {
        &mut self.base
    }
    fn perform_line(
        &mut self,
        l: &Line2d,
        d1: &Res2dDomain,
        c: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntConicCurveOfGInter::perform_line(self, l, d1, c, d2, tol_conf, tol)
    }
    fn perform_circle(
        &mut self,
        c: &Circle2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntConicCurveOfGInter::perform_circle(self, c, d1, pcurve, d2, tol_conf, tol)
    }
    fn perform_ellipse(
        &mut self,
        e: &Ellipse2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntConicCurveOfGInter::perform_ellipse(self, e, d1, pcurve, d2, tol_conf, tol)
    }
    fn perform_parabola(
        &mut self,
        p: &Parabola2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntConicCurveOfGInter::perform_parabola(self, p, d1, pcurve, d2, tol_conf, tol)
    }
    fn perform_hyperbola(
        &mut self,
        h: &Hyperbola2d,
        d1: &Res2dDomain,
        pcurve: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntConicCurveOfGInter::perform_hyperbola(self, h, d1, pcurve, d2, tol_conf, tol)
    }
}

impl crate::geomalgo::int_curve_curve_gen::IntPCurvePCurveMember<dyn Curve2dAdaptor>
    for TheIntPCurvePCurveOfGInter
{
    fn base(&self) -> &IntersectionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut IntersectionBase {
        &mut self.base
    }
    fn perform(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntPCurvePCurveOfGInter::perform(self, c1, d1, c2, d2, tol_conf, tol)
    }
    fn perform_cd(&mut self, c: &dyn Curve2dAdaptor, d: &Res2dDomain, tol_conf: f64, tol: f64) {
        TheIntPCurvePCurveOfGInter::perform_cd(self, c, d, tol_conf, tol)
    }
    fn set_min_nb_samples(&mut self, the_min_nb_samples: i32) {
        TheIntPCurvePCurveOfGInter::set_min_nb_samples(self, the_min_nb_samples)
    }
    fn get_min_nb_samples(&self) -> i32 {
        TheIntPCurvePCurveOfGInter::get_min_nb_samples(self)
    }
}

/// OCCT Geom2dInt_GInter (Geom2dInt_GInter_0.cxx) — the concrete
/// IntCurve_IntCurveCurveGen instantiation for the GInter chain.
pub type GInter = crate::geomalgo::int_curve_curve_gen::IntCurveCurveGen<
    dyn Curve2dAdaptor,
    Geom2dCurveTool,
    TheIntConicCurveOfGInter,
    TheIntPCurvePCurveOfGInter,
>;
