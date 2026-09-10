//! OCCT GCPnts curve-interface projection (TKGeomBase/GCPnts).
//!
//! OCCT's GCPnts algorithms are templates over `Adaptor3d_Curve` /
//! `Adaptor2d_Curve2d`, evaluated through the virtual adaptor interface
//! (D0/D1/D2/GetType/Circle/Bezier/BSpline/NbIntervals/Intervals).
//! The 2D overloads of GCPnts_TangentialDeflection (GCPnts_TangentialDeflection.cxx
//! L46-67) lift the 2D evaluations into gp_Pnt/gp_Vec with a zero Z slot; this
//! module encodes the same interface as one Rust trait returning DVec3.
//!
//! Implementations:
//! - [`GCPntsCurve2d`] — the `Adaptor2d_Curve2d` flavor over the kernel
//!   `Curve2d` (the `Geom2dAdaptor_Curve` semantics),
//! - `dyn Adaptor3dCurve` — the `Adaptor3d_Curve` flavor; the kernel adaptor
//!   trait carries no curve-type accessor, so [`GCPntsCurve::get_type`]
//!   conservatively reports `CurveType::Other` (the OCCT
//!   `Adaptor3d_CurveOnSurface` reports `GeomAbs_OtherCurve` for all
//!   non-KPart forms too; the analytic short-circuits of computeType then
//!   fall through to the Gauss-integration path, which computes the same
//!   lengths).

use glam::DVec3;

use rcad_kernel::base::proj_lib::adaptor::Adaptor3dCurve;
use rcad_kernel::base::proj_lib::adaptor::CurveHandle;
use rcad_kernel::base::proj_lib::adaptor::CurveOnSurface;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::Curve2d;
use rcad_kernel::Curve2dEval;

/// OCCT GCPnts/GCPnts_TangentialDeflection template interface over
/// `Adaptor3d_Curve` / `Adaptor2d_Curve2d` (with the L46-67 2D shims).
pub trait GCPntsCurve {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64;
    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64;
    /// OCCT D0(U, P) — the 2D shims (TangentialDeflection.cxx L46-53) lift
    /// the result into a 3D point with Z = 0.
    fn d0(&self, u: f64) -> DVec3;
    /// OCCT D1(U, P, V) — Z = 0 for the 2D flavor.
    fn d1(&self, u: f64) -> (DVec3, DVec3);
    /// OCCT D2(U, P, V1, V2) — Z = 0 for the 2D flavor.
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3);
    /// OCCT GetType().
    fn get_type(&self) -> CurveType;
    /// OCCT NbIntervals(GeomAbs_CN).
    fn nb_intervals_cn(&self) -> usize;
    /// OCCT Intervals(T, GeomAbs_CN) — the interval bounds, first..last.
    fn intervals_cn(&self) -> Vec<f64>;
    /// OCCT Circle().Radius() — valid when GetType() == Circle.
    fn circle_radius(&self) -> f64;
    /// OCCT Degree() — the BSpline/Bezier degree.
    fn curve_degree(&self) -> i32;
    /// OCCT NbPoles() — the BSpline/Bezier pole count.
    fn nb_poles(&self) -> i32;
    /// OCCT IsRational() — the BSpline/Bezier rationality.
    fn is_rational(&self) -> bool;
    /// OCCT DN(U, 1) — the first derivative as a curve-kind accessor.
    fn dn1(&self, u: f64) -> DVec3;
}

fn to_3d(p: glam::DVec2) -> DVec3 {
    DVec3::new(p.x, p.y, 0.0)
}

/// The `Adaptor2d_Curve2d` flavor over the kernel `Curve2d`
/// (`Geom2dAdaptor_Curve` semantics).
pub struct GCPntsCurve2d<'a> {
    pub curve: &'a Curve2d,
}

impl<'a> GCPntsCurve2d<'a> {
    /// OCCT Geom2dAdaptor_Curve(C).
    pub fn new(curve: &'a Curve2d) -> Self {
        GCPntsCurve2d { curve }
    }
}

impl GCPntsCurve for GCPntsCurve2d<'_> {
    /// OCCT Geom2dAdaptor_Curve::FirstParameter.
    fn first_parameter(&self) -> f64 {
        self.curve.default_domain()[0]
    }

    /// OCCT Geom2dAdaptor_Curve::LastParameter.
    fn last_parameter(&self) -> f64 {
        self.curve.default_domain()[1]
    }

    /// OCCT D0 shim (TangentialDeflection.cxx L46-53).
    fn d0(&self, u: f64) -> DVec3 {
        to_3d(self.curve.point_at(u))
    }

    /// OCCT D1 shim.
    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        (
            to_3d(self.curve.point_at(u)),
            to_3d(self.curve.derivative_at(u)),
        )
    }

    /// OCCT D2 shim (TangentialDeflection.cxx L55-67).
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        (
            to_3d(self.curve.point_at(u)),
            to_3d(self.curve.derivative_at(u)),
            to_3d(self.curve.derivative2_at(u)),
        )
    }

    /// OCCT Geom2dAdaptor_Curve::GetType — the basis-curve kind
    /// (the kernel Geom2dCurveAdaptor mapping, adaptor.rs curve2d_type_of).
    fn get_type(&self) -> CurveType {
        match self.curve.inner() {
            Curve2d::Line(_) => CurveType::Line,
            Curve2d::Circle(_) => CurveType::Circle,
            Curve2d::Ellipse(_) => CurveType::Ellipse,
            Curve2d::Parabola(_) => CurveType::Parabola,
            Curve2d::Hyperbola(_) => CurveType::Hyperbola,
            Curve2d::Bezier(_) => CurveType::Bezier,
            Curve2d::BSpline(_) => CurveType::BSpline,
            _ => CurveType::Other,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::NbIntervals(GeomAbs_CN) — 1 for the
    /// non-composite rcad curve kinds (the codebase-wide encoding).
    fn nb_intervals_cn(&self) -> usize {
        1
    }

    /// OCCT Intervals(T, GeomAbs_CN) — the whole domain as one interval.
    fn intervals_cn(&self) -> Vec<f64> {
        vec![self.first_parameter(), self.last_parameter()]
    }

    /// OCCT Circle().Radius().
    fn circle_radius(&self) -> f64 {
        match self.curve.inner() {
            Curve2d::Circle(c) => c.radius,
            _ => panic!("GCPntsCurve2d::Circle on a non-circle curve"),
        }
    }

    /// OCCT Degree() — the BSpline/Bezier degree.
    fn curve_degree(&self) -> i32 {
        match self.curve.inner() {
            Curve2d::BSpline(b) => b.degree as i32,
            Curve2d::Bezier(b) => (b.control_points.len() - 1) as i32,
            _ => panic!("GCPntsCurve2d::Degree on a non-bspline/bezier curve"),
        }
    }

    /// OCCT NbPoles().
    fn nb_poles(&self) -> i32 {
        match self.curve.inner() {
            Curve2d::BSpline(b) => b.control_points.len() as i32,
            Curve2d::Bezier(b) => b.control_points.len() as i32,
            _ => panic!("GCPntsCurve2d::NbPoles on a non-bspline/bezier curve"),
        }
    }

    /// OCCT IsRational().
    fn is_rational(&self) -> bool {
        match self.curve.inner() {
            Curve2d::BSpline(b) => b.weights.iter().any(|&w| w != 1.0),
            Curve2d::Bezier(b) => b.weights.iter().any(|&w| w != 1.0),
            _ => panic!("GCPntsCurve2d::IsRational on a non-bspline/bezier curve"),
        }
    }

    /// OCCT DN(U, 1).
    fn dn1(&self, u: f64) -> DVec3 {
        to_3d(self.curve.derivative_at(u))
    }
}

/// The `Adaptor3d_CurveOnSurface` flavor (the GeomPlate frontiere type) —
/// the same evaluations as the trait-object form; the curve type stays the
/// OCCT `GeomAbs_OtherCurve` default of the non-KPart forms (EvalKPart
/// staged in the kernel adaptor).
impl GCPntsCurve for CurveOnSurface {
    fn first_parameter(&self) -> f64 {
        Adaptor3dCurve::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        Adaptor3dCurve::last_parameter(self)
    }

    fn d0(&self, u: f64) -> DVec3 {
        Adaptor3dCurve::value(self, u)
    }

    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        Adaptor3dCurve::d1(self, u)
    }

    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        Adaptor3dCurve::d2(self, u)
    }

    fn get_type(&self) -> CurveType {
        CurveType::Other
    }

    fn nb_intervals_cn(&self) -> usize {
        Adaptor3dCurve::nb_intervals(self, rcad_kernel::math::GeomAbsShape::CN)
    }

    fn intervals_cn(&self) -> Vec<f64> {
        Adaptor3dCurve::intervals(self, rcad_kernel::math::GeomAbsShape::CN)
    }

    fn circle_radius(&self) -> f64 {
        panic!("GCPntsCurve: Circle() on a non-circle curve");
    }

    fn curve_degree(&self) -> i32 {
        panic!("GCPntsCurve: Degree() on a non-bspline/bezier curve");
    }

    fn nb_poles(&self) -> i32 {
        panic!("GCPntsCurve: NbPoles() on a non-bspline/bezier curve");
    }

    fn is_rational(&self) -> bool {
        panic!("GCPntsCurve: IsRational() on a non-bspline/bezier curve");
    }

    fn dn1(&self, u: f64) -> DVec3 {
        Adaptor3dCurve::d1(self, u).1
    }
}

/// The `Adaptor3d_Curve` flavor — served over the kernel adaptor trait
/// object (the OCCT virtual dispatch through `const Adaptor3d_Curve&`).
/// The handle wrapper keeps the owning Arc and routes the evaluations
/// through the adaptor trait.
pub struct GCPntsCurve3dHandle(pub CurveHandle);

impl GCPntsCurve for GCPntsCurve3dHandle {
    fn first_parameter(&self) -> f64 {
        self.0.first_parameter()
    }

    fn last_parameter(&self) -> f64 {
        self.0.last_parameter()
    }

    fn d0(&self, u: f64) -> DVec3 {
        self.0.value(u)
    }

    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.0.d1(u)
    }

    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.0.d2(u)
    }

    /// See the `dyn Adaptor3dCurve` impl: the kernel adaptor trait carries
    /// no curve-type accessor; `GeomAbs_OtherCurve` is reported.
    fn get_type(&self) -> CurveType {
        CurveType::Other
    }

    fn nb_intervals_cn(&self) -> usize {
        self.0.nb_intervals(rcad_kernel::math::GeomAbsShape::CN)
    }

    fn intervals_cn(&self) -> Vec<f64> {
        self.0.intervals(rcad_kernel::math::GeomAbsShape::CN)
    }

    fn circle_radius(&self) -> f64 {
        panic!("GCPntsCurve: Circle() on a non-circle curve");
    }

    fn curve_degree(&self) -> i32 {
        panic!("GCPntsCurve: Degree() on a non-bspline/bezier curve");
    }

    fn nb_poles(&self) -> i32 {
        panic!("GCPntsCurve: NbPoles() on a non-bspline/bezier curve");
    }

    fn is_rational(&self) -> bool {
        panic!("GCPntsCurve: IsRational() on a non-bspline/bezier curve");
    }

    fn dn1(&self, u: f64) -> DVec3 {
        self.0.d1(u).1
    }
}

/// The `Adaptor3d_Curve` flavor — served over the kernel adaptor trait
/// object (the OCCT virtual dispatch through `const Adaptor3d_Curve&`).
impl GCPntsCurve for dyn Adaptor3dCurve {
    fn first_parameter(&self) -> f64 {
        Adaptor3dCurve::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        Adaptor3dCurve::last_parameter(self)
    }

    /// OCCT D0(U, P) — the kernel adaptor folds Value/D0 into `value`.
    fn d0(&self, u: f64) -> DVec3 {
        Adaptor3dCurve::value(self, u)
    }

    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        Adaptor3dCurve::d1(self, u)
    }

    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        Adaptor3dCurve::d2(self, u)
    }

    /// OCCT GetType() — the kernel adaptor trait carries no curve-type
    /// accessor; `GeomAbs_OtherCurve` is reported (see the module docs).
    fn get_type(&self) -> CurveType {
        CurveType::Other
    }

    /// OCCT NbIntervals(GeomAbs_CN).
    fn nb_intervals_cn(&self) -> usize {
        Adaptor3dCurve::nb_intervals(self, rcad_kernel::math::GeomAbsShape::CN)
    }

    /// OCCT Intervals(T, GeomAbs_CN).
    fn intervals_cn(&self) -> Vec<f64> {
        Adaptor3dCurve::intervals(self, rcad_kernel::math::GeomAbsShape::CN)
    }

    /// OCCT Circle().Radius() — raises NoSuchObject for the non-circle type
    /// (Adaptor3d_CurveOnSurface.cxx L1390-1395); unreachable under the
    /// `Other` type dispatch.
    fn circle_radius(&self) -> f64 {
        panic!("GCPntsCurve: Circle() on a non-circle curve");
    }

    fn curve_degree(&self) -> i32 {
        panic!("GCPntsCurve: Degree() on a non-bspline/bezier curve");
    }

    fn nb_poles(&self) -> i32 {
        panic!("GCPntsCurve: NbPoles() on a non-bspline/bezier curve");
    }

    fn is_rational(&self) -> bool {
        panic!("GCPntsCurve: IsRational() on a non-bspline/bezier curve");
    }

    fn dn1(&self, u: f64) -> DVec3 {
        Adaptor3dCurve::d1(self, u).1
    }
}
