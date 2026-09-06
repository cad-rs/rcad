// OCCT HLRBRep_CurveTool (TKHLR) — the static curve tool over
// HLRBRep_CurvePtr for the CInter instantiation.
//
// HLRBRep_CurveTool.hxx L17-152 + .cxx L19-83 + .lxx (the one-line
// delegates to HLRBRep_Curve).  The marker struct plays the TheCurveTool /
// ParTool / ThePCurveTool template roles of the CInter sub-engines
// (IntImpParGen_Intersector, IntCurve_IntConicCurveGen,
// IntCurve_UserIntConicCurveGen, IntCurve_IntCurveCurveGen).
//
// Documented conventions (the CurveView cannot express them yet):
// - NbIntervals(C1) == 1 for the C1-continuous rcad edge curve kinds
//   (the Geom2dAdaptor_Curve / Curve2dAdaptor::nb_intervals_c1 default),
//   so Intervals fills [FirstParameter, LastParameter];
// - the NbSamples(C, U1, U2) BSpline branch uses the single C0 interval
//   ((NbIntervals(CN)+1) * Degree = 2 * Degree) pending the C1-interval
//   machinery.

use glam::DVec2;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};

use super::curve::Curve;
use crate::geomalgo::geom2d_int::{Curve2dType, LocatorCurveTool};
use crate::geomalgo::int_imp_par_gen::ParTool;
use crate::geomalgo::user_int_conic_curve_gen::PCurveTool;

/// OCCT HLRBRep_CurveTool (the EpsX constant, lxx L241-245).
pub const CURVE_TOOL_EPS_X: f64 = 1.0e-10;

/// OCCT HLRBRep_CurveTool — the static tool marker.
pub struct CurveTool;

/// OCCT GeomAbs_CurveType of the projected curve (HLRBRep_Curve::GetType).
pub fn get_type_2d(t: CurveType) -> Curve2dType {
    match t {
        CurveType::Line => Curve2dType::Line,
        CurveType::Circle => Curve2dType::Circle,
        CurveType::Ellipse => Curve2dType::Ellipse,
        CurveType::Parabola => Curve2dType::Parabola,
        CurveType::Hyperbola => Curve2dType::Hyperbola,
        CurveType::Bezier => Curve2dType::BezierCurve,
        CurveType::BSpline => Curve2dType::BSplineCurve,
        CurveType::Other => Curve2dType::OtherCurve,
    }
}

impl CurveTool {
    /// OCCT FirstParameter (lxx).
    pub fn first_parameter(c: &Curve<'_>) -> f64 {
        c.first_parameter()
    }

    /// OCCT LastParameter (lxx).
    pub fn last_parameter(c: &Curve<'_>) -> f64 {
        c.last_parameter()
    }

    /// OCCT NbIntervals (lxx) — the single C1 interval convention.
    pub fn nb_intervals(_c: &Curve<'_>) -> i32 {
        1
    }

    /// OCCT Intervals (lxx) — fills Tab(1, NbIntervals+1); `tab` is 0-based
    /// storage for the 1-based array.
    pub fn intervals(c: &Curve<'_>, tab: &mut [f64]) {
        let n = CurveTool::nb_intervals(c) as usize;
        tab[0] = c.first_parameter();
        tab[n] = c.last_parameter();
    }

    /// OCCT GetInterval (lxx) — a = Tab(i), b = Tab(i+1).
    pub fn get_interval(_c: &Curve<'_>, i: usize, tab: &[f64]) -> (f64, f64) {
        (tab[i - 1], tab[i])
    }

    /// OCCT IsClosed (lxx).
    pub fn is_closed(c: &Curve<'_>) -> bool {
        c.is_closed()
    }

    /// OCCT IsPeriodic (lxx).
    pub fn is_periodic(c: &Curve<'_>) -> bool {
        c.is_periodic()
    }

    /// OCCT Period (lxx).
    pub fn period(c: &Curve<'_>) -> f64 {
        c.period()
    }

    /// OCCT Value (lxx).
    pub fn value(c: &Curve<'_>, u: f64) -> DVec2 {
        c.value(u)
    }

    /// OCCT D0 (lxx).
    pub fn d0(c: &Curve<'_>, u: f64) -> DVec2 {
        c.value(u)
    }

    /// OCCT D1 (lxx).
    pub fn d1(c: &Curve<'_>, u: f64) -> (DVec2, DVec2) {
        let (mut p, mut v) = (DVec2::ZERO, DVec2::ZERO);
        c.d1_2d(u, &mut p, &mut v);
        (p, v)
    }

    /// OCCT D2 (lxx).
    pub fn d2(c: &Curve<'_>, u: f64) -> (DVec2, DVec2, DVec2) {
        let (mut p, mut v1, mut v2) = (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO);
        c.d2_2d(u, &mut p, &mut v1, &mut v2);
        (p, v1, v2)
    }

    /// OCCT D3 (lxx) — the HLRBRep_Curve::D3 empty body.
    pub fn d3(c: &Curve<'_>, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        let (mut p, mut v1, mut v2, mut v3) = (
            DVec2::ZERO,
            DVec2::ZERO,
            DVec2::ZERO,
            DVec2::ZERO,
        );
        c.d3_2d(u, &mut p, &mut v1, &mut v2, &mut v3);
        (p, v1, v2, v3)
    }

    /// OCCT DN (lxx) — the HLRBRep_Curve::DN zero vector.
    pub fn dn(c: &Curve<'_>, u: f64, n: i32) -> DVec2 {
        c.dn(u, n)
    }

    /// OCCT Resolution (lxx).
    pub fn resolution(c: &Curve<'_>, r3d: f64) -> f64 {
        c.resolution(r3d)
    }

    /// OCCT GetType (lxx).
    pub fn get_type(c: &Curve<'_>) -> Curve2dType {
        get_type_2d(c.get_type())
    }

    /// OCCT Line (lxx).
    pub fn line(c: &Curve<'_>) -> Line2d {
        c.line()
    }

    /// OCCT Circle (lxx).
    pub fn circle(c: &Curve<'_>) -> Circle2d {
        c.circle()
    }

    /// OCCT Ellipse (lxx).
    pub fn ellipse(c: &Curve<'_>) -> Ellipse2d {
        c.ellipse()
    }

    /// OCCT Hyperbola (lxx) — the default-constructed gp_Hypr2d verbatim
    /// (the projected classification never reports Hyperbola).
    pub fn hyperbola(_c: &Curve<'_>) -> Hyperbola2d {
        Hyperbola2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            semi_major: 0.0,
            semi_minor: 0.0,
        }
    }

    /// OCCT Parabola (lxx) — the default-constructed gp_Parab2d verbatim.
    pub fn parabola(_c: &Curve<'_>) -> Parabola2d {
        Parabola2d {
            origin: DVec2::ZERO,
            axis_dir: DVec2::X,
            focal_param: 0.0,
        }
    }

    /// OCCT EpsX (lxx L241-245).
    pub fn eps_x(_c: &Curve<'_>) -> f64 {
        CURVE_TOOL_EPS_X
    }

    /// OCCT Degree (lxx).
    pub fn degree(c: &Curve<'_>) -> i32 {
        c.degree()
    }

    /// OCCT NbSamples(C) (cxx L20-46): Line 2, Bezier 3+NbPoles,
    /// BSpline NbKnots*Degree (min 2), other 10, clamped to 50.
    pub fn nb_samples(c: &Curve<'_>) -> i32 {
        const NBS_OTHER: f64 = 10.0;
        let mut nbs = NBS_OTHER;
        let typ_c = c.get_type();
        if typ_c == CurveType::Line {
            nbs = 2.0;
        } else if typ_c == CurveType::Bezier {
            nbs = 3.0 + c.nb_poles() as f64;
        } else if typ_c == CurveType::BSpline {
            nbs = c.nb_knots() as f64;
            nbs *= c.degree() as f64;
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        if nbs > 50.0 {
            nbs = 50.0;
        }
        nbs as i32
    }

    /// OCCT NbSamples(C, U1, U2) (cxx L48-74): the BSpline branch counts the
    /// C0 intervals of the 3D curve (GeomAdaptor NbIntervals(CN) + 1) — the
    /// single-interval convention gives 2 * Degree.
    pub fn nb_samples_uv(c: &Curve<'_>, _u1: f64, _u2: f64) -> i32 {
        const NBS_OTHER: f64 = 10.0;
        let mut nbs = NBS_OTHER;
        let typ_c = c.get_type();
        if typ_c == CurveType::Line {
            nbs = 2.0;
        } else if typ_c == CurveType::Bezier {
            nbs = 3.0 + c.nb_poles() as f64;
        } else if typ_c == CurveType::BSpline {
            nbs = (1 + 1) as f64 * c.degree() as f64;
            if nbs < 2.0 {
                nbs = 2.0;
            }
        }
        if nbs > 50.0 {
            nbs = 50.0;
        }
        nbs as i32
    }
}

impl<'a> PCurveTool<Curve<'a>> for CurveTool {
    fn nb_intervals(c: &Curve<'a>) -> i32 {
        CurveTool::nb_intervals(c)
    }
    fn intervals(c: &Curve<'a>, tab: &mut [f64]) {
        CurveTool::intervals(c, tab)
    }
    fn get_interval(c: &Curve<'a>, index: usize, tab: &[f64]) -> (f64, f64) {
        CurveTool::get_interval(c, index, tab)
    }
    fn first_parameter(c: &Curve<'a>) -> f64 {
        CurveTool::first_parameter(c)
    }
    fn last_parameter(c: &Curve<'a>) -> f64 {
        CurveTool::last_parameter(c)
    }
    fn value(c: &Curve<'a>, u: f64) -> DVec2 {
        CurveTool::value(c, u)
    }
    fn get_type(c: &Curve<'a>) -> Curve2dType {
        CurveTool::get_type(c)
    }
    fn line(c: &Curve<'a>) -> Line2d {
        CurveTool::line(c)
    }
    fn circle(c: &Curve<'a>) -> Circle2d {
        CurveTool::circle(c)
    }
    fn ellipse(c: &Curve<'a>) -> Ellipse2d {
        CurveTool::ellipse(c)
    }
    fn parabola(c: &Curve<'a>) -> Parabola2d {
        CurveTool::parabola(c)
    }
    fn hyperbola(c: &Curve<'a>) -> Hyperbola2d {
        CurveTool::hyperbola(c)
    }
}

impl<'a> ParTool<Curve<'a>> for CurveTool {
    fn nb_samples_uv(c: &Curve<'a>, u1: f64, u2: f64) -> i32 {
        CurveTool::nb_samples_uv(c, u1, u2)
    }
    fn eps_x(c: &Curve<'a>) -> f64 {
        CurveTool::eps_x(c)
    }
    fn value(c: &Curve<'a>, u: f64) -> DVec2 {
        CurveTool::value(c, u)
    }
    fn d1(c: &Curve<'a>, u: f64) -> (DVec2, DVec2) {
        CurveTool::d1(c, u)
    }
    fn d2(c: &Curve<'a>, u: f64) -> (DVec2, DVec2, DVec2) {
        CurveTool::d2(c, u)
    }
}

impl<'a> LocatorCurveTool<Curve<'a>> for CurveTool {
    fn first_parameter(c: &Curve<'a>) -> f64 {
        CurveTool::first_parameter(c)
    }
    fn last_parameter(c: &Curve<'a>) -> f64 {
        CurveTool::last_parameter(c)
    }
    fn value(c: &Curve<'a>, u: f64) -> DVec2 {
        CurveTool::value(c, u)
    }
    fn d0(c: &Curve<'a>, u: f64) -> DVec2 {
        CurveTool::d0(c, u)
    }
    fn d1(c: &Curve<'a>, u: f64) -> (DVec2, DVec2) {
        CurveTool::d1(c, u)
    }
    fn d2(c: &Curve<'a>, u: f64) -> (DVec2, DVec2, DVec2) {
        CurveTool::d2(c, u)
    }
    fn dn(c: &Curve<'a>, u: f64, n: i32) -> DVec2 {
        CurveTool::dn(c, u, n)
    }
    fn get_type(c: &Curve<'a>) -> Curve2dType {
        CurveTool::get_type(c)
    }
    fn nb_samples(c: &Curve<'a>) -> i32 {
        CurveTool::nb_samples(c)
    }
    fn eps_x(c: &Curve<'a>) -> f64 {
        CurveTool::eps_x(c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::geom::{Circle3, Ellipse3, Line3, Point3, Vec3};

    use crate::hlr::algo::projector::Projector;
    use crate::hlr::brep::b_curve_tool::CurveView;

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double).
    pub(crate) struct LineEdge(pub f64);

    impl CurveView for LineEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.0
        }
        fn d0(&self, u: f64) -> Point3 {
            Point3::new(u, 0.0, 0.0)
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (Point3::new(u, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0))
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (Point3::new(u, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> Line3 {
            Line3::new(Point3::ZERO, Vec3::new(1.0, 0.0, 0.0))
        }
        fn circle(&self) -> Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn degree(&self) -> i32 {
            0
        }
        fn nb_poles(&self) -> i32 {
            0
        }
        fn nb_knots(&self) -> i32 {
            0
        }
        fn is_closed(&self) -> bool {
            false
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            0.0
        }
        fn resolution(&self, r3d: f64) -> f64 {
            r3d
        }
        fn parameter_3d(&self, p2d: f64) -> f64 {
            p2d
        }
        fn poles(&self) -> Vec<Point3> {
            Vec::new()
        }
    }

    /// The top-view projector over z = 0 (identity image).
    pub(crate) fn top_projector() -> Projector {
        Projector::from_ax2(&Ax2::new(
            glam::DVec3::ZERO,
            glam::DVec3::new(0.0, 0.0, 1.0),
            glam::DVec3::new(1.0, 0.0, 0.0),
        ))
    }

    /// A loaded HLRBRep_Curve over a static edge view and projector.
    pub(crate) fn projected_curve(
        proj: &'static Projector,
        edge: &'static LineEdge,
    ) -> Curve<'static> {
        let mut c = Curve::new();
        c.projector(proj);
        c.load(edge);
        c
    }

    /// OCCT anchor: the CurveTool over a projected line — the lxx delegates
    /// (FirstParameter/Value/D1/Line), NbSamples(Line) = 2 (cxx L25-28) and
    /// EpsX = 1e-10 (lxx L241-245).
    #[test]
    fn curve_tool_over_projected_line_anchor() {
        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(top_projector);
        static EDGE: std::sync::OnceLock<LineEdge> = std::sync::OnceLock::new();
        let edge: &'static LineEdge = EDGE.get_or_init(|| LineEdge(3.0));
        let c = projected_curve(proj, edge);
        // The Update pass runs the projected-curve classification
        // (HLRBRep_Curve::Update sets myType).
        let mut c = c;
        c.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);

        assert_eq!(CurveTool::first_parameter(&c), 0.0);
        assert_eq!(CurveTool::last_parameter(&c), 3.0);
        assert_eq!(CurveTool::get_type(&c), Curve2dType::Line);
        let p = CurveTool::value(&c, 2.0);
        assert!((p.x - 2.0).abs() < 1e-12 && p.y.abs() < 1e-12);
        let (pp, v) = CurveTool::d1(&c, 1.0);
        assert!((pp.x - 1.0).abs() < 1e-12);
        assert!((v.x - 1.0).abs() < 1e-12 && v.y.abs() < 1e-12);
        let l = CurveTool::line(&c);
        assert!((l.direction.x.abs() - 1.0).abs() < 1e-12);
        // NbSamples (cxx L20-46): a Line projects to nbs = 2.
        assert_eq!(CurveTool::nb_samples(&c), 2);
        assert_eq!(CurveTool::nb_samples_uv(&c, 0.0, 3.0), 2);
        assert_eq!(CurveTool::eps_x(&c), 1.0e-10);
        // Intervals: the single C1 interval [0, 3].
        let mut tab = [0.0f64; 3];
        CurveTool::intervals(&c, &mut tab);
        assert_eq!(CurveTool::get_interval(&c, 1, &tab), (0.0, 3.0));
    }
}
