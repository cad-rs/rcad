// OCCT HLRBRep_BCurveTool (TKHLR) — the static tool over the edge curve
// adaptor (`BRepAdaptor_Curve`).
//
// HLRBRep_BCurveTool.hxx L36-150 + .cxx L25-129 + .lxx (one-line forwards).
// In rcad the OCCT statics map one-to-one onto the [`CurveView`] trait
// methods (the template-parameter-to-trait pattern); this module carries
// the sample-count dispatch ([`nb_samples`], cxx L25-55) and the free
// one-liners ([`b_curve_value`]).
//
// The OCCT Bezier/BSpline accessors (Bezier(C)/BSpline(C), Poles,
// PolesAndWeights, cxx L59-129) are provided through the trait's pole
// queries — the HLRBRep_Curve Poles/Knots/Multiplicities methods consume
// them.

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Circle3, Ellipse3, Line3, Point3, Vec3};

/// The `BRepAdaptor_Curve` interface consumed by HLRBRep_Curve and
/// HLRBRep_Surface::IsAbove — the full HLRBRep_BCurveTool static set.
pub trait CurveView {
    /// OCCT FirstParameter(C).
    fn first_parameter(&self) -> f64;
    /// OCCT LastParameter(C).
    fn last_parameter(&self) -> f64;
    /// OCCT Value(C, U) / D0(C, U, P).
    fn d0(&self, u: f64) -> Point3;
    /// OCCT D1(C, U, P, V).
    fn d1(&self, u: f64) -> (Point3, Vec3);
    /// OCCT D2(C, U, P, V1, V2) — Standard_NoSuchObject on non-C2 curves.
    fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
        panic!("Standard_NoSuchObject: BCurveTool::D2");
    }
    /// OCCT GetType(C).
    fn get_type(&self) -> CurveType;
    /// OCCT Line(C) — valid when GetType() == Line.
    fn line(&self) -> Line3;
    /// OCCT Circle(C) — valid when GetType() == Circle.
    fn circle(&self) -> Circle3;
    /// OCCT Ellipse(C) — valid when GetType() == Ellipse.
    fn ellipse(&self) -> Ellipse3;
    /// OCCT Degree(C).
    fn degree(&self) -> i32;
    /// OCCT NbPoles(C).
    fn nb_poles(&self) -> i32;
    /// OCCT NbKnots(C).
    fn nb_knots(&self) -> i32;
    /// OCCT IsClosed(C).
    fn is_closed(&self) -> bool;
    /// OCCT IsPeriodic(C).
    fn is_periodic(&self) -> bool;
    /// OCCT Period(C).
    fn period(&self) -> f64;
    /// OCCT Resolution(C, R3d).
    fn resolution(&self, r3d: f64) -> f64;
    /// OCCT Parameter3d(P2d) — the 2D -> 3D parameter bridge of
    /// HLRBRep_Curve (carried here so IsAbove keeps its OCCT call shape).
    fn parameter_3d(&self, p2d: f64) -> f64;
    /// OCCT Poles(C, T) (cxx L59-79) — the pole array; empty for the other
    /// curve types.
    fn poles(&self) -> Vec<Point3>;
}

/// OCCT NbSamples(C, U0, U1) (cxx L25-55).
pub fn nb_samples(c: &dyn CurveView, u0: f64, u1: f64) -> i32 {
    let typ_c = c.get_type();
    // static double nbsOther = 10.0;
    let nbs_other = 10.0f64;
    let mut nbs = nbs_other;

    if typ_c == CurveType::Line {
        nbs = 2.0;
    } else if typ_c == CurveType::Bezier {
        nbs = 3.0 + c.nb_poles() as f64;
    } else if typ_c == CurveType::BSpline {
        nbs = c.nb_knots() as f64;
        nbs *= c.degree() as f64;
        nbs *= c.last_parameter() - c.first_parameter();
        nbs /= u1 - u0;
        if nbs < 2.0 {
            nbs = 2.0;
        }
    }
    if nbs > 50.0 {
        nbs = 50.0;
    }
    nbs as i32
}

/// OCCT Value(C, U) — the D0 evaluation one-liner (lxx).
pub fn b_curve_value(c: &dyn CurveView, u: f64) -> Point3 {
    c.d0(u)
}
