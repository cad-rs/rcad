// OCCT HLRBRep_CInter (TKHLR) — the 2D projected-curve intersection for
// the HLR chain: the IntCurve_IntCurveCurveGen instantiation
// (HLRBRep_CInter.hxx L181-232 alias table + HLRBRep_CInter_0.cxx).
//
// The instantiation family (each OCCT `The*OfCInter` class is the generic
// engine compiled over the HLRBRep tool bindings; rcad expresses them as
// type-level instantiations of the landed generic engines):
//   - HLRBRep_CurveTool                          -> [`curve_tool::CurveTool`]
//   - HLRBRep_TheProjPCurOfCInter                -> [`TheProjPCurOfCInter`]
//     (IntCurve_ProjectOnPCurveGen; the shared FindParameter body over the
//     CurveLocator + LocateExtPC machinery)
//   - HLRBRep_TheIntersectorOfTheIntConicCurveOfCInter
//     -> [`TheIntersectorOfCInter`] (IntImpParGen_Intersector)
//   - HLRBRep_TheIntConicCurveOfCInter           -> [`TheIntConicCurveOfCInter`]
//     (IntCurve_IntConicCurveGen)
//   - HLRBRep_IntConicCurveOfCInter              -> [`IntConicCurveOfCInter`]
//     (IntCurve_UserIntConicCurveGen)
//   - HLRBRep_ThePolygon2dOfTheIntPCurvePCurveOfCInter /
//     HLRBRep_TheDistBetweenPCurvesOfTheIntPCurvePCurveOfCInter /
//     HLRBRep_ExactIntersectionPointOfTheIntPCurvePCurveOfCInter — the
//     IntCurve_IntPolyPolyGen helpers, bound inside the generic engine over
//     [`CPolyTool`]
//   - HLRBRep_TheIntPCurvePCurveOfCInter         -> [`TheIntPCurvePCurveOfCInter`]
//   - HLRBRep_CInter                             -> [`CInter`]
//     (IntCurve_IntCurveCurveGen over the members above)

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};

use super::curve::Curve;
use super::curve_tool::CurveTool;
use crate::geomalgo::geom2d_int::{
    proj_cur_find_parameter_bounded, proj_cur_find_parameter_unbounded, Curve2dAdaptor,
};
use crate::geomalgo::int_conic_curve_gen::IntConicCurveGen;
use crate::geomalgo::int_curve_curve_gen::{
    IntConicCurveMember, IntCurveCurveGen, IntPCurvePCurveMember,
};
use crate::geomalgo::int_curve_generics::ProjPCurveTool;
use crate::geomalgo::int_imp_par_gen::{Intersector as ImpParIntersector, ProjectOnPCurveTool};
use crate::geomalgo::int_poly_poly_gen::IntPolyPolyGen;
use crate::geomalgo::int_res2d::{Domain as Res2dDomain, IntersectionBase};
use crate::geomalgo::user_int_conic_curve_gen::UserIntConicCurveGen;

// ---------------------------------------------------------------------------
// HLRBRep_TheProjPCurOfCInter (TheProjPCurOfCInter_0.cxx L24-80)
// ---------------------------------------------------------------------------

/// OCCT HLRBRep_TheProjPCurOfCInter — the point projection on the projected
/// curve (the ProjectOnPCurveTool role of the CInter instantiation).
#[derive(Debug, Clone, Copy, Default)]
pub struct TheProjPCurOfCInter;

impl<'a> ProjectOnPCurveTool<Curve<'a>> for TheProjPCurOfCInter {
    /// OCCT FindParameter(C, P, LowParameter, HighParameter, Tol)
    /// (TheProjPCurOfCInter_0.cxx L24-56).
    fn find_parameter_between(c: &Curve<'a>, p: DVec2, low: f64, high: f64, tol: f64) -> f64 {
        proj_cur_find_parameter_bounded::<Curve<'a>, CurveTool>(c, p, low, high, tol)
    }

    /// OCCT FindParameter(C, P, Tol) (TheProjPCurOfCInter_0.cxx L58-80).
    fn find_parameter(c: &Curve<'a>, p: DVec2, tol: f64) -> f64 {
        proj_cur_find_parameter_unbounded::<Curve<'a>, CurveTool>(c, p, tol)
    }
}

// ---------------------------------------------------------------------------
// HLRBRep_TheIntersectorOfTheIntConicCurveOfCInter
// (TheIntersectorOfTheIntConicCurveOfCInter_0.cxx: ImpTool =
// IntCurve_IConicTool, ParCurve = HLRBRep_CurvePtr, ParTool =
// HLRBRep_CurveTool, ProjectOnPCurveTool = HLRBRep_TheProjPCurOfCInter,
// IntImpParGen_MyImpParTool = HLRBRep_MyImpParToolOfTheIntersector...)
// ---------------------------------------------------------------------------

/// OCCT HLRBRep_TheIntersectorOfTheIntConicCurveOfCInter.
#[derive(Debug, Clone)]
pub struct TheIntersectorOfCInter {
    pub base: IntersectionBase,
}

impl TheIntersectorOfCInter {
    pub fn new() -> Self {
        TheIntersectorOfCInter {
            base: IntersectionBase::new(),
        }
    }

    /// OCCT IntImpParGen_Intersector::Perform (IntImpParGen_Intersector.gxx
    /// L245-779) via the generic engine.
    pub fn perform(
        &mut self,
        imp_tool: &crate::geomalgo::geom2d_int::IConicTool,
        imp_domain: &Res2dDomain,
        par_curve: &Curve<'_>,
        par_domain: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        fn engine<'a>(
        ) -> ImpParIntersector<Curve<'a>, CurveTool, TheProjPCurOfCInter> {
            ImpParIntersector::new()
        }
        let mut myintersection = engine();
        myintersection
            .base
            .set_reversed_parameters(self.base.reversed_parameters());
        myintersection.perform(imp_tool, imp_domain, par_curve, par_domain, tol_conf, tol);
        self.base = myintersection.base;
    }
}

impl Default for TheIntersectorOfCInter {
    fn default() -> Self {
        TheIntersectorOfCInter::new()
    }
}

// ---------------------------------------------------------------------------
// HLRBRep_TheIntConicCurveOfCInter (TheIntConicCurveOfCInter_0.cxx: the
// IntCurve_IntConicCurveGen instantiation) and HLRBRep_IntConicCurveOfCInter
// (IntConicCurveOfCInter_0.cxx: the IntCurve_UserIntConicCurveGen
// instantiation).
// ---------------------------------------------------------------------------

/// OCCT HLRBRep_TheIntConicCurveOfCInter.
pub type TheIntConicCurveOfCInter<'a> = IntConicCurveGen<Curve<'a>, CurveTool, TheProjPCurOfCInter>;

/// OCCT HLRBRep_IntConicCurveOfCInter.
pub type IntConicCurveOfCInter<'a> =
    UserIntConicCurveGen<Curve<'a>, CurveTool, TheProjPCurOfCInter>;

impl<'a> Default for TheIntConicCurveOfCInter<'a> {
    fn default() -> Self {
        TheIntConicCurveOfCInter::bare()
    }
}

impl<'a> IntConicCurveMember<Curve<'a>> for TheIntConicCurveOfCInter<'a> {
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
        c: &Curve<'a>,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        IntConicCurveGen::perform_line(self, l, d1, c, d2, tol_conf, tol)
    }
    fn perform_circle(
        &mut self,
        c: &Circle2d,
        d1: &Res2dDomain,
        pcurve: &Curve<'a>,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        IntConicCurveGen::perform_circle(self, c, d1, pcurve, d2, tol_conf, tol)
    }
    fn perform_ellipse(
        &mut self,
        e: &Ellipse2d,
        d1: &Res2dDomain,
        pcurve: &Curve<'a>,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        IntConicCurveGen::perform_ellipse(self, e, d1, pcurve, d2, tol_conf, tol)
    }
    fn perform_parabola(
        &mut self,
        p: &Parabola2d,
        d1: &Res2dDomain,
        pcurve: &Curve<'a>,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        IntConicCurveGen::perform_parabola(self, p, d1, pcurve, d2, tol_conf, tol)
    }
    fn perform_hyperbola(
        &mut self,
        h: &Hyperbola2d,
        d1: &Res2dDomain,
        pcurve: &Curve<'a>,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        IntConicCurveGen::perform_hyperbola(self, h, d1, pcurve, d2, tol_conf, tol)
    }
}

// ---------------------------------------------------------------------------
// HLRBRep_TheIntPCurvePCurveOfCInter (TheIntPCurvePCurveOfCInter_0.cxx: the
// IntCurve_IntPolyPolyGen instantiation; the Polygon2dGen /
// DistBetweenPCurvesGen / ExactIntersectionPoint helpers of the gxx are the
// HLRBRep_ThePolygon2dOf... / HLRBRep_TheDistBetweenPCurvesOf... /
// HLRBRep_ExactIntersectionPointOf... instantiations, bound through the
// CPolyTool tool triple).
// ---------------------------------------------------------------------------

/// OCCT tool triple of TheIntPCurvePCurveOfCInter_0.cxx (TheCurveTool =
/// HLRBRep_CurveTool, TheProjPCur = HLRBRep_TheProjPCurOfCInter) bundled
/// for the IntCurve_IntPolyPolyGen template parameters — the
/// [`crate::geomalgo::geom2d_int::GInterPolyTool`] precedent.
#[derive(Debug, Clone, Copy, Default)]
pub struct CPolyTool<'a>(std::marker::PhantomData<&'a ()>);

impl<'a> ProjPCurveTool for CPolyTool<'a> {
    type Curve = Curve<'a>;
    fn value(c: &Curve<'a>, u: f64) -> DVec2 {
        CurveTool::value(c, u)
    }
    fn d1(c: &Curve<'a>, u: f64) -> (DVec2, DVec2) {
        CurveTool::d1(c, u)
    }
    fn eps_x(c: &Curve<'a>) -> f64 {
        CurveTool::eps_x(c)
    }
}

impl<'a> crate::geomalgo::int_imp_par_gen::ParTool<Curve<'a>> for CPolyTool<'a> {
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

impl<'a> ProjectOnPCurveTool<Curve<'a>> for CPolyTool<'a> {
    fn find_parameter(c: &Curve<'a>, p: DVec2, tol: f64) -> f64 {
        proj_cur_find_parameter_unbounded::<Curve<'a>, CurveTool>(c, p, tol)
    }
    fn find_parameter_between(c: &Curve<'a>, p: DVec2, low: f64, high: f64, tol: f64) -> f64 {
        proj_cur_find_parameter_bounded::<Curve<'a>, CurveTool>(c, p, low, high, tol)
    }
}

/// OCCT HLRBRep_TheIntPCurvePCurveOfCInter.
#[derive(Debug, Clone)]
pub struct TheIntPCurvePCurveOfCInter {
    pub base: IntersectionBase,
    pub domain_on_curve1: Res2dDomain,
    pub domain_on_curve2: Res2dDomain,
    min_pnt_nb: i32,
}

impl TheIntPCurvePCurveOfCInter {
    /// OCCT IntCurve_IntPolyPolyGen() (gxx L84-89).
    pub fn new() -> Self {
        TheIntPCurvePCurveOfCInter {
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
    pub fn perform<'a>(
        &mut self,
        c1: &Curve<'a>,
        d1: &Res2dDomain,
        c2: &Curve<'a>,
        d2: &Res2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        fn engine<'a>(
        ) -> IntPolyPolyGen<Curve<'a>, CPolyTool<'a>> {
            IntPolyPolyGen::new()
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
    pub fn perform_cd<'a>(
        &mut self,
        c1: &Curve<'a>,
        d1: &Res2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        fn engine<'a>(
        ) -> IntPolyPolyGen<Curve<'a>, CPolyTool<'a>> {
            IntPolyPolyGen::new()
        }
        let mut e = engine();
        e.set_min_nb_samples(self.min_pnt_nb);
        e.perform_auto(c1, d1, the_tol_conf, the_tol);
        self.base = e.base;
        self.domain_on_curve1 = e.domain_on_curve1;
        self.domain_on_curve2 = e.domain_on_curve2;
    }
}

impl Default for TheIntPCurvePCurveOfCInter {
    fn default() -> Self {
        TheIntPCurvePCurveOfCInter::new()
    }
}

impl<'a> IntPCurvePCurveMember<Curve<'a>> for TheIntPCurvePCurveOfCInter {
    fn base(&self) -> &IntersectionBase {
        &self.base
    }
    fn base_mut(&mut self) -> &mut IntersectionBase {
        &mut self.base
    }
    fn perform(
        &mut self,
        c1: &Curve<'a>,
        d1: &Res2dDomain,
        c2: &Curve<'a>,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        TheIntPCurvePCurveOfCInter::perform(self, c1, d1, c2, d2, tol_conf, tol)
    }
    fn perform_cd(&mut self, c: &Curve<'a>, d: &Res2dDomain, tol_conf: f64, tol: f64) {
        TheIntPCurvePCurveOfCInter::perform_cd(self, c, d, tol_conf, tol)
    }
    fn set_min_nb_samples(&mut self, the_min_nb_samples: i32) {
        TheIntPCurvePCurveOfCInter::set_min_nb_samples(self, the_min_nb_samples)
    }
    fn get_min_nb_samples(&self) -> i32 {
        TheIntPCurvePCurveOfCInter::get_min_nb_samples(self)
    }
}

// ---------------------------------------------------------------------------
// HLRBRep_CInter (CInter.hxx L43-179: the IntCurve_IntCurveCurveGen
// instantiation over the members above).
// ---------------------------------------------------------------------------

/// OCCT HLRBRep_CInter.
pub type CInter<'a> = IntCurveCurveGen<
    Curve<'a>,
    CurveTool,
    TheIntConicCurveOfCInter<'a>,
    TheIntPCurvePCurveOfCInter,
>;

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::base::proj_lib::CurveType;
    use rcad_kernel::geom::{Circle3, Ellipse3, Line3, Point3, Vec3};
    use rcad_kernel::math::gp::Ax2;

    use crate::hlr::algo::projector::Projector;
    use crate::hlr::brep::b_curve_tool::CurveView;

    /// A straight 3D segment adaptor on an arbitrary 3D line.
    pub(crate) struct SegEdge {
        pub origin: Point3,
        pub dir: Vec3,
        pub len: f64,
    }

    impl CurveView for SegEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.len
        }
        fn d0(&self, u: f64) -> Point3 {
            self.origin + self.dir * u
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (self.origin + self.dir * u, self.dir)
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (self.origin + self.dir * u, self.dir, Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> Line3 {
            Line3 {
                origin: self.origin,
                direction: self.dir,
            }
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

    fn loaded_curve() -> CInter<'static> {
        CInter::new()
    }

    /// OCCT anchor: CInter::Perform(C1, C2, TolConf, Tol) over two straight
    /// edges seen in the top view — the projected lines (the X-axis and the
    /// vertical x = 0.5) cross at (0.5, 0) with both parameters at the
    /// middle (IntCurveCurveGen.gxx L235-262 dispatches Line x Line to
    /// IntConicConic::Perform(Lin, Lin)).
    #[test]
    fn c_inter_line_line_cross_anchor() {
        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(|| {
            Projector::from_ax2(&Ax2::new(
                DVec3::ZERO,
                DVec3::new(0.0, 0.0, 1.0),
                DVec3::new(1.0, 0.0, 0.0),
            ))
        });
        static E1: std::sync::OnceLock<SegEdge> = std::sync::OnceLock::new();
        let e1: &'static SegEdge = E1.get_or_init(|| SegEdge {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 1.0,
        });
        static E2: std::sync::OnceLock<SegEdge> = std::sync::OnceLock::new();
        let e2: &'static SegEdge = E2.get_or_init(|| SegEdge {
            origin: Point3::new(0.5, -0.5, 0.0),
            dir: Vec3::new(0.0, 1.0, 0.0),
            len: 1.0,
        });

        let mut c1 = Curve::new();
        c1.projector(proj);
        c1.load(e1);
        c1.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);
        let mut c2 = Curve::new();
        c2.projector(proj);
        c2.load(e2);
        c2.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);

        let mut inter = loaded_curve();
        inter.perform_cc(&c1, &c2, 1.0e-9, 1.0e-9);

        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 1, "nb={}", inter.nb_points());
        let p = inter.point(1);
        // The crossing at (0.5, 0): u = 0.5 on C1, v = 0.5 on C2 (both
        // Middles).
        assert!((p.value().x - 0.5).abs() < 1e-9, "x={}", p.value().x);
        assert!((p.value().y - 0.0).abs() < 1e-9, "y={}", p.value().y);
        assert!((p.param_on_first() - 0.5).abs() < 1e-9);
        assert!((p.param_on_second() - 0.5).abs() < 1e-9);
    }

    /// OCCT anchor: CInter over the projections of two non-conic (Bezier)
    /// edges — the Other x Other arm routes through TheIntPCurvePCurve
    /// (IntPolyPolyGen) and still finds the exact crossing point.
    #[test]
    fn c_inter_other_other_polygon_anchor() {
        // The fixture: two quadratic 3D Bezier-like segments (non-conic view
        // kind) whose top-view images are the straight segments
        // (0,0)-(1,1) and (0,1)-(1,0), crossing at (0.5, 0.5).  The view
        // kind OtherCurve forces the intcurvcurv arm with the polygon
        // refinement (CurveTool::NbSamples -> 10 for non-conic kinds).
        struct OtherEdge {
            pub pts: [Point3; 3],
        }
        // A degree-2 curve evaluator over three control points; the
        // projected image is computed through the projector.
        impl CurveView for OtherEdge {
            fn first_parameter(&self) -> f64 {
                0.0
            }
            fn last_parameter(&self) -> f64 {
                1.0
            }
            fn d0(&self, u: f64) -> Point3 {
                let b0 = (1.0 - u) * (1.0 - u);
                let b1 = 2.0 * (1.0 - u) * u;
                let b2 = u * u;
                self.pts[0] * b0 + self.pts[1] * b1 + self.pts[2] * b2
            }
            fn d1(&self, u: f64) -> (Point3, Vec3) {
                let db0 = -2.0 * (1.0 - u);
                let db1 = 2.0 - 4.0 * u;
                let db2 = 2.0 * u;
                (
                    self.d0(u),
                    self.pts[0] * db0 + self.pts[1] * db1 + self.pts[2] * db2,
                )
            }
            fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
                (
                    self.d0(u),
                    self.d1(u).1,
                    (self.pts[0] * 2.0) - (self.pts[1] * 4.0) + (self.pts[2] * 2.0),
                )
            }
            fn get_type(&self) -> CurveType {
                // A projected non-conic curve reports Other (HLRBRep_Curve::
                // Update keeps myType = OtherCurve for non-conic views).
                CurveType::Other
            }
            fn line(&self) -> Line3 {
                panic!("Standard_NoSuchObject");
            }
            fn circle(&self) -> Circle3 {
                panic!("Standard_NoSuchObject");
            }
            fn ellipse(&self) -> Ellipse3 {
                panic!("Standard_NoSuchObject");
            }
            fn degree(&self) -> i32 {
                2
            }
            fn nb_poles(&self) -> i32 {
                3
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
                self.pts.to_vec()
            }
        }

        static PROJ: std::sync::OnceLock<Projector> = std::sync::OnceLock::new();
        let proj: &'static Projector = PROJ.get_or_init(|| {
            Projector::from_ax2(&Ax2::new(
                DVec3::ZERO,
                DVec3::new(0.0, 0.0, 1.0),
                DVec3::new(1.0, 0.0, 0.0),
            ))
        });
        // The curves lie in z = 0, so the top view is the identity: the
        // images are the straight segments (0,0)-(1,1) and (0,1)-(1,0).
        // The poles are slightly off the chord so the image is genuinely a
        // non-line conic-free curve for the polygon machinery.
        static E1: std::sync::OnceLock<OtherEdge> = std::sync::OnceLock::new();
        let e1: &'static OtherEdge = E1.get_or_init(|| OtherEdge {
            pts: [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.5, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
        });
        static E2: std::sync::OnceLock<OtherEdge> = std::sync::OnceLock::new();
        let e2: &'static OtherEdge = E2.get_or_init(|| OtherEdge {
            pts: [
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.5, 0.5, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
        });

        let mut c1 = Curve::new();
        c1.projector(proj);
        c1.load(e1);
        c1.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);
        let mut c2 = Curve::new();
        c2.projector(proj);
        c2.load(e2);
        c2.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);

        let mut inter = IntCurveCurveGen::<
            Curve<'static>,
            CurveTool,
            TheIntConicCurveOfCInter<'static>,
            TheIntPCurvePCurveOfCInter,
        >::new();
        inter.perform_cc(&c1, &c2, 1.0e-7, 1.0e-7);

        assert!(inter.is_done());
        assert!(inter.nb_points() >= 1, "nb={}", inter.nb_points());
        let mut found = false;
        for i in 1..=inter.nb_points() {
            let p = inter.point(i);
            if (p.value().x - 0.5).abs() < 1e-6 && (p.value().y - 0.5).abs() < 1e-6 {
                found = true;
            }
        }
        assert!(found, "the (0.5, 0.5) crossing was not found");
    }

    // The Curve2dAdaptor import keeps the doc link on the alias table (the
    // OCCT TheCurve = HLRBRep_CurvePtr never converts to Curve2dAdaptor).
    #[allow(unused)]
    fn _alias_parity(_c: &dyn Curve2dAdaptor) {}
}
