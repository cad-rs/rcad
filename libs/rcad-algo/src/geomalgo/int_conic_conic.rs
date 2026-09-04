// OCCT IntCurve_IntConicConic + IntCurve_PConic + IntCurve_PConicTool —
// 1:1 Rust translation (TKGeomAlgo/IntCurve).
//
// IntCurve_IntConicConic (IntCurve_IntConicConic.hxx/.cxx/.lxx) intersects
// two conics. Its member Inter is an IntCurve_IntImpConicParConic (the
// IntImpParGen_Intersector template instantiated with
// ParCurve=IntCurve_PConic / ParTool=IntCurve_PConicTool); rcad keeps a
// single intersector implementation parameterized by the Curve2dAdaptor
// trait, so PConic implements Curve2dAdaptor with the IntCurve_PConicTool
// semantics (Value/D1/D2/EpsX/NbSamples per PConicTool.cxx L24-131).
//
// Ported so far: the Ellipse-Ellipse Perform overload
// (IntCurve_IntConicConic.cxx L915-958) — the path used by
// Geom2dAPI_InterCurveCurve for two Geom2d_Ellipses (OCC29289 anchor).
//
// Deferred overloads (each needs the IntAna2d offset pre-pass +
// SetBinfBsupFromIntAna2d, or the dedicated closed-form implementation):
//   - Perform(Lin,  Parab)    IntCurve_IntConicConic.cxx L109-226
//   - Perform(Lin,  Hypr)     L230-333
//   - Perform(Circ, Parab)    L337-435
//   - Perform(Circ, Elips)    L439-482
//   - Perform(Circ, Hypr)     L486-581
//   - Perform(Parab, Parab)   L585-688
//   - Perform(Elips, Parab)   L692-806
//   - Perform(Parab, Hypr)    L810-911
//   - Perform(Elips, Hypr)    L962-1063
//   - Perform(Hypr,  Hypr)    L1067-1168
//   - Perform(Circ, Circ)     IntCurve_IntConicConic_1.cxx L807-1238
//   - Perform(Lin,  Lin)      IntCurve_IntConicConic_1.cxx L1381-2233
//   - Perform(Lin,  Circ)     IntCurve_IntConicConic_1.cxx L2236-2652
// plus the (L, Elips)/(L, ...)/... constructor overloads in
// IntCurve_IntConicConic.lxx.

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};

use super::geom2d_int::{
    elclib2d, Curve2dAdaptor, Curve2dType, IConicTool, TheIntersectorOfTheIntConicCurveOfGInter,
};
use super::int_res2d::{Domain as Res2dDomain, IntersectionBase};

/// 2*pi (OCCT: M_PI + M_PI).
const PI2: f64 = std::f64::consts::TAU;

// ---------------------------------------------------------------------------
// IntCurve_PConic
// ---------------------------------------------------------------------------

/// OCCT IntCurve_PConic — a conic from gp represented as a parametric curve.
/// "The Conics are manipulated as objects which only depend on three
/// parameters: Axis and two Reals" (IntCurve_PConic.hxx L64-67).
#[derive(Debug, Clone)]
pub struct PConic {
    /// OCCT axe (gp_Ax22d) — origin and X direction of the local frame.
    axe_origin: DVec2,
    axe_xdir: DVec2,
    prm1: f64,
    prm2: f64,
    the_eps_x: f64,
    the_accuracy: i32,
    typ: Curve2dType,
}

impl PConic {
    /// OCCT IntCurve_PConic(const gp_Elips2d& E) (IntCurve_PConic.cxx L28-36).
    pub fn new_ellipse(e: &Ellipse2d) -> Self {
        PConic {
            axe_origin: e.center,
            axe_xdir: e.major_dir,
            prm1: e.major_radius,
            prm2: e.minor_radius,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            typ: Curve2dType::Ellipse,
        }
    }

    /// OCCT IntCurve_PConic(const gp_Hypr2d& H) (IntCurve_PConic.cxx L38-46).
    pub fn new_hyperbola(h: &Hyperbola2d) -> Self {
        PConic {
            axe_origin: h.center,
            axe_xdir: h.major_dir,
            prm1: h.semi_major,
            prm2: h.semi_minor,
            the_eps_x: 0.00000001,
            the_accuracy: 50,
            typ: Curve2dType::Hyperbola,
        }
    }

    /// OCCT IntCurve_PConic(const gp_Circ2d& C) (IntCurve_PConic.cxx L48-56).
    pub fn new_circle(c: &Circle2d) -> Self {
        PConic {
            axe_origin: c.center,
            axe_xdir: c.x_dir,
            prm1: c.radius,
            prm2: 0.0,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            typ: Curve2dType::Circle,
        }
    }

    /// OCCT IntCurve_PConic(const gp_Parab2d& P) (IntCurve_PConic.cxx L58-66).
    pub fn new_parabola(p: &Parabola2d) -> Self {
        PConic {
            axe_origin: p.origin,
            axe_xdir: p.axis_dir,
            prm1: p.focal_param,
            prm2: 0.0,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            typ: Curve2dType::Parabola,
        }
    }

    /// OCCT IntCurve_PConic(const gp_Lin2d& L) (IntCurve_PConic.cxx L68-76).
    pub fn new_line(l: &Line2d) -> Self {
        PConic {
            axe_origin: l.origin,
            axe_xdir: l.direction,
            prm1: 0.0,
            prm2: 0.0,
            the_eps_x: 0.00000001,
            the_accuracy: 20,
            typ: Curve2dType::Line,
        }
    }

    /// OCCT SetEpsX(EpsDist) (IntCurve_PConic.cxx L78-81).
    pub fn set_eps_x(&mut self, epsx: f64) {
        self.the_eps_x = epsx;
    }

    /// OCCT SetAccuracy(Nb) (IntCurve_PConic.cxx L83-86).
    pub fn set_accuracy(&mut self, n: i32) {
        self.the_accuracy = n;
    }

    /// OCCT Accuracy().
    pub fn accuracy(&self) -> i32 {
        self.the_accuracy
    }

    /// OCCT EpsX().
    pub fn eps_x(&self) -> f64 {
        self.the_eps_x
    }

    /// OCCT TypeCurve().
    pub fn type_curve(&self) -> Curve2dType {
        self.typ
    }

    /// OCCT Axis2() — (origin, X direction) of the local frame.
    pub fn axis2(&self) -> (DVec2, DVec2) {
        (self.axe_origin, self.axe_xdir)
    }

    /// OCCT Param1().
    pub fn param1(&self) -> f64 {
        self.prm1
    }

    /// OCCT Param2().
    pub fn param2(&self) -> f64 {
        self.prm2
    }

    /// Y direction of the local frame (direct Ax22d: +90 deg rotation).
    fn axe_ydir(&self) -> DVec2 {
        DVec2::new(-self.axe_xdir.y, self.axe_xdir.x)
    }
}

// ---------------------------------------------------------------------------
// IntCurve_PConicTool — the Curve2dAdaptor binding of PConic
// ---------------------------------------------------------------------------

/// OCCT IntCurve_PConicTool (IntCurve_PConicTool.hxx/.cxx L24-131) expressed
/// through the Curve2dAdaptor interface used by the single rcad
/// IntImpParGen_Intersector implementation.
impl Curve2dAdaptor for PConic {
    /// PConicTool provides no FirstParameter/LastParameter: a PConic is
    /// unbounded and its evaluation bounds always come from the IntRes2d_Domain
    /// handed to the intersector. The value is never consumed on the
    /// IntImpParGen chain (the only reader is the BSpline branch of
    /// NbSamples, and a PConic is never a BSpline).
    fn first_parameter(&self) -> f64 {
        0.0
    }

    fn last_parameter(&self) -> f64 {
        0.0
    }

    /// OCCT IntCurve_PConicTool::Value (PConicTool.cxx L24-44).
    fn value(&self, u: f64) -> DVec2 {
        let ydir = self.axe_ydir();
        match self.typ {
            Curve2dType::Line => elclib2d::line_value(self.axe_origin, self.axe_xdir, u),
            Curve2dType::Circle => {
                elclib2d::circle_value(self.axe_origin, self.axe_xdir, ydir, self.prm1, u)
            }
            Curve2dType::Ellipse => {
                elclib2d::ellipse_value(self.axe_origin, self.axe_xdir, ydir, self.prm1, self.prm2, u)
            }
            Curve2dType::Parabola => {
                elclib2d::parabola_value(self.axe_origin, self.axe_xdir, ydir, self.prm1, u)
            }
            // -- case GeomAbs_Hyperbola:
            _ => elclib2d::hyperbola_value(self.axe_origin, self.axe_xdir, ydir, self.prm1, self.prm2, u),
        }
    }

    /// OCCT IntCurve_PConicTool::D1 (PConicTool.cxx L47-78).
    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        let ydir = self.axe_ydir();
        match self.typ {
            Curve2dType::Line => elclib2d::line_d1(self.axe_origin, self.axe_xdir, u),
            Curve2dType::Circle => {
                elclib2d::circle_d1(self.axe_origin, self.axe_xdir, ydir, self.prm1, u)
            }
            Curve2dType::Ellipse => {
                elclib2d::ellipse_d1(self.axe_origin, self.axe_xdir, ydir, self.prm1, self.prm2, u)
            }
            Curve2dType::Parabola => {
                elclib2d::parabola_d1(self.axe_origin, self.axe_xdir, ydir, self.prm1, u)
            }
            _ => elclib2d::hyperbola_d1(self.axe_origin, self.axe_xdir, ydir, self.prm1, self.prm2, u),
        }
    }

    /// OCCT IntCurve_PConicTool::D2 (PConicTool.cxx L81-114). For the Line case
    /// OCCT zeroes Tan then calls LineD1 and leaves Norm unset (a line tangent
    /// is never degenerate, so Norm is never consumed by the transition
    /// logic) — rcad returns the zero vector for Norm.
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        let ydir = self.axe_ydir();
        match self.typ {
            Curve2dType::Line => {
                let (p, t) = elclib2d::line_d1(self.axe_origin, self.axe_xdir, u);
                (p, t, DVec2::ZERO)
            }
            Curve2dType::Circle => {
                elclib2d::circle_d2(self.axe_origin, self.axe_xdir, ydir, self.prm1, u)
            }
            Curve2dType::Ellipse => {
                elclib2d::ellipse_d2(self.axe_origin, self.axe_xdir, ydir, self.prm1, self.prm2, u)
            }
            Curve2dType::Parabola => {
                elclib2d::parabola_d2(self.axe_origin, self.axe_xdir, ydir, self.prm1, u)
            }
            _ => elclib2d::hyperbola_d2(self.axe_origin, self.axe_xdir, ydir, self.prm1, self.prm2, u),
        }
    }

    /// PConicTool has no D3 — never called on the IntImpParGen chain. Returns
    /// the D2 data with a zero third derivative.
    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        let (p, t, n) = self.d2(u);
        (p, t, n, DVec2::ZERO)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            _ => DVec2::ZERO,
        }
    }

    /// OCCT GetType() — the stored TypeCurve.
    fn get_type(&self) -> Curve2dType {
        self.typ
    }

    /// OCCT IntCurve_PConicTool::NbSamples(C) (PConicTool.cxx L121-124) —
    /// the Accuracy of the PConic.
    fn nb_samples(&self) -> i32 {
        self.the_accuracy
    }

    fn resolution(&self, r3d: f64) -> f64 {
        r3d
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

    fn nb_knots(&self) -> i32 {
        0
    }

    fn degree(&self) -> i32 {
        0
    }

    fn nb_poles(&self) -> i32 {
        0
    }

    /// Reconstructed from the stored axis + two reals.
    fn circle(&self) -> Circle2d {
        Circle2d {
            center: self.axe_origin,
            x_dir: self.axe_xdir,
            y_dir: self.axe_ydir(),
            radius: self.prm1,
        }
    }

    fn line(&self) -> Line2d {
        Line2d::new(self.axe_origin, self.axe_xdir)
    }

    fn ellipse(&self) -> Ellipse2d {
        Ellipse2d {
            center: self.axe_origin,
            major_dir: self.axe_xdir,
            major_radius: self.prm1,
            minor_radius: self.prm2,
        }
    }

    fn parabola(&self) -> Parabola2d {
        Parabola2d {
            origin: self.axe_origin,
            axis_dir: self.axe_xdir,
            focal_param: self.prm1,
        }
    }

    fn hyperbola(&self) -> Hyperbola2d {
        Hyperbola2d {
            center: self.axe_origin,
            major_dir: self.axe_xdir,
            semi_major: self.prm1,
            semi_minor: self.prm2,
        }
    }
}

// ---------------------------------------------------------------------------
// IntCurve_IntConicConic
// ---------------------------------------------------------------------------

/// OCCT IntCurve_IntConicConic — conic x conic intersection
/// (IntCurve_IntConicConic.hxx). The member Inter is the
/// IntImpParGen_Intersector instantiation; rcad reuses
/// TheIntersectorOfTheIntConicCurveOfGInter with a PConic pcurve.
#[derive(Debug, Clone)]
pub struct IntConicConic {
    pub base: IntersectionBase,
    /// OCCT member Inter.
    inter: TheIntersectorOfTheIntConicCurveOfGInter,
}

impl IntConicConic {
    /// OCCT IntCurve_IntConicConic() (lxx L20).
    pub fn new() -> Self {
        IntConicConic {
            base: IntersectionBase::new(),
            inter: TheIntersectorOfTheIntConicCurveOfGInter::new(),
        }
    }

    /// OCCT constructor (E1, D1, E2, D2, TolConf, Tol) (lxx L128-138).
    pub fn new_ellipse_ellipse(
        e1: &Ellipse2d,
        d1: &Res2dDomain,
        e2: &Ellipse2d,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntConicConic::new();
        r.perform_ellipse_ellipse(e1, d1, e2, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntConicConic::Perform(const gp_Elips2d& E1,
    /// const IntRes2d_Domain& DE1, const gp_Elips2d& E2,
    /// const IntRes2d_Domain& DE2, TolConf, Tol)
    /// (IntCurve_IntConicConic.cxx L915-958).
    pub fn perform_ellipse_ellipse(
        &mut self,
        e1: &Ellipse2d,
        de1: &Res2dDomain,
        e2: &Ellipse2d,
        de2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let itool = IConicTool::new_ellipse(e1);
        let mut pcurve = PConic::new_ellipse(e2);
        pcurve.set_accuracy(20);

        self.inter.base.set_reversed_parameters(self.base.reversed_parameters());

        if !de1.is_closed() {
            let mut d1 = de1.clone();
            d1.set_equivalent_parameters(de1.first_parameter(), de1.first_parameter() + PI2);
            if !de2.is_closed() {
                let mut d2 = de2.clone();
                d2.set_equivalent_parameters(de2.first_parameter(), de2.first_parameter() + PI2);
                self.inter.perform(&itool, &d1, &pcurve, &d2, tol_conf, tol);
            } else {
                self.inter.perform(&itool, &d1, &pcurve, de2, tol_conf, tol);
            }
        } else {
            if !de2.is_closed() {
                let mut d2 = de2.clone();
                d2.set_equivalent_parameters(de2.first_parameter(), de2.first_parameter() + PI2);
                self.inter.perform(&itool, de1, &pcurve, &d2, tol_conf, tol);
            } else {
                self.inter.perform(&itool, de1, &pcurve, de2, tol_conf, tol);
            }
        }
        self.base.set_values(&self.inter.base);
    }
}

impl Default for IntConicConic {
    fn default() -> Self {
        IntConicConic::new()
    }
}
