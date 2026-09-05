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
use super::int_res2d::{
    Domain as Res2dDomain, IntersectionBase, IntersectionPoint, IntersectionSegment, Position,
    Situation, Transition, TypeTrans,
};

/// 2*pi (OCCT: M_PI + M_PI).
const PI2: f64 = std::f64::consts::TAU;

/// OCCT IntCurve_IntConicConic_1.cxx L40-44: the file-local angular
/// tolerance (1.e-15 — "at least to make an accordance between transition
/// and position computation").
const TOLERANCE_ANGULAIRE: f64 = 1.0e-15;
/// OCCT Precision::PConfusion().
const PRECISION_P_CONFUSION: f64 = 1.0e-9;

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

    fn as_any(&self) -> &dyn std::any::Any {
        self
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

    // -- Remaining conic x conic Perform overloads ---------------------------
    //
    // IntCurveCurveGen::InternalPerform dispatches into these by curve kind.
    // Each body is the OCCT overload listed in the file header; they are
    // ported in a follow-up unit (they need the IntAna2d offset pre-pass +
    // SetBinfBsupFromIntAna2d, or the dedicated closed-form implementation).

    /// OCCT Perform(const gp_Lin2d& L1, const gp_Lin2d& L2)
    /// (IntCurve_IntConicConic_1.cxx L1381-2233).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_line_line(
        &mut self,
        l1: &Line2d,
        domain1: &Res2dDomain,
        l2: &Line2d,
        domain2: &Res2dDomain,
        _tol_conf: f64,
        tol_r: f64,
    ) {
        self.base.reset_fields();

        //-- Coordonnees du point d intersection sur chacune des 2 droites
        let mut u1 = 0.0;
        let mut u2 = 0.0;
        //-- Nombre de points solution : 1 : Intersection
        //--                             0 : Non Confondues
        //--                             2 : Confondues a la tolerance pres
        let mut nbsol = 0i32;
        let mut pt_seg1 = IntersectionPoint::empty();
        let mut pt_seg2 = IntersectionPoint::empty();
        let mut a_half_sin_l1_l2 = 0.0;
        let mut tol = tol_r;
        if tol < PRECISION_P_CONFUSION {
            tol = PRECISION_P_CONFUSION;
        }

        line_line_geometric_intersection(l1, l2, tol, &mut u1, &mut u2, &mut a_half_sin_l1_l2, &mut nbsol);

        let tan1 = l1.direction;
        let tan2 = l2.direction;

        let a_cos_t1_t2 = tan1.dot(tan2);
        let is_opposite = a_cos_t1_t2 < 0.0;

        self.base.done = true;

        if nbsol == 1 && check_ll_coincidence(l1, l2, domain1, domain2, tol) {
            nbsol = 2;
        }

        if nbsol == 1 {
            //---------------------------------------------------
            //-- d: distance du point I a partir de laquelle  les
            //--  points de parametre U1+d et U2+-d sont ecartes
            //--  d une distance superieure a Tol.
            //---------------------------------------------------
            let mut pos1a = Position::Middle;
            let mut pos2a = Position::Middle;
            let mut pos1b = Position::Middle;
            let mut pos2b = Position::Middle;
            let d = 0.5 * tol / a_half_sin_l1_l2;
            let mut u1inf = u1 - d;
            let mut u1sup = u1 + d;
            let u1mu2 = u1 - u2;
            let u1pu2 = u1 + u2;
            let mut res1inf = 0.0;
            let mut res1sup = 0.0;
            let mut prod_vect_tan;

            //---------------------------------------------------
            //-- On agrandit la zone U1inf U1sup pour tenir compte
            //-- des tolerances des points en bout
            //--
            if domain1.has_first_point() {
                if l2.distance(domain1.first_point()) < domain1.first_tolerance() {
                    if u1inf > domain1.first_parameter() {
                        u1inf = domain1.first_parameter();
                    }
                    if u1sup < domain1.first_parameter() {
                        u1sup = domain1.first_parameter();
                    }
                }
            }
            if domain1.has_last_point() {
                if l2.distance(domain1.last_point()) < domain1.last_tolerance() {
                    if u1inf > domain1.last_parameter() {
                        u1inf = domain1.last_parameter();
                    }
                    if u1sup < domain1.last_parameter() {
                        u1sup = domain1.last_parameter();
                    }
                }
            }
            if domain2.has_first_point() {
                if l1.distance(domain2.first_point()) < domain2.first_tolerance() {
                    let p = elclib2d::line_parameter(l1.origin, l1.direction, domain2.first_point());
                    if u1inf > p {
                        u1inf = p;
                    }
                    if u1sup < p {
                        u1sup = p;
                    }
                }
            }
            if domain2.has_last_point() {
                if l1.distance(domain2.last_point()) < domain2.last_tolerance() {
                    let p = elclib2d::line_parameter(l1.origin, l1.direction, domain2.last_point());
                    if u1inf > p {
                        u1inf = p;
                    }
                    if u1sup < p {
                        u1sup = p;
                    }
                }
            }
            //-----------------------------------------------------------------

            domain_intersection(
                domain1,
                u1inf,
                u1sup,
                &mut res1inf,
                &mut res1sup,
                &mut pos1a,
                &mut pos1b,
            );

            if (res1sup - res1inf) < 0.0 {
                //-- Si l intersection est vide
                //--
            } else {
                //-- (Domain1  INTER   Zone Intersection)    non vide
                prod_vect_tan = tan1.x * tan2.y - tan1.y * tan2.x;

                // #####################################################################
                // ##  Longueur Minimale d un segment    Sur Courbe 1
                // #####################################################################

                let long_mini_seg = tol;

                if ((res1sup - res1inf) <= long_mini_seg)
                    || ((pos1a == pos1b) && (pos1a != Position::Middle))
                {
                    //-------------------------------  Un seul Point -------------------
                    //--- lorsque la longueur du segment est inferieure a ??
                    //--- ou si deux points designent le meme bout
                    // gka #0022833
                    let a_cur_trans = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                        TypeTrans::Out
                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                        TypeTrans::In
                    } else {
                        TypeTrans::Undecided
                    };

                    let mut new_point1 = IntersectionPoint::empty();
                    if compute_int_point(
                        domain1,
                        domain2,
                        l1,
                        l2,
                        a_cos_t1_t2,
                        u1,
                        u2,
                        &mut res1inf,
                        &mut res1sup,
                        1,
                        a_cur_trans,
                        &mut new_point1,
                    ) {
                        self.base.append_point(&new_point1);
                    }

                    //------------------------------------------------------
                } //---------------   Fin du cas  :   1 seul point --------------------
                else {
                    //-- Intersection AND Domain1  --------> Segment ---------------------
                    let mut u2inf;
                    let mut u2sup;
                    let mut res2inf = 0.0;
                    let mut res2sup = 0.0;

                    if is_opposite {
                        u2inf = u1pu2 - res1sup;
                        u2sup = u1pu2 - res1inf;
                    } else {
                        u2inf = res1inf - u1mu2;
                        u2sup = res1sup - u1mu2;
                    }

                    domain_intersection(
                        domain2,
                        u2inf,
                        u2sup,
                        &mut res2inf,
                        &mut res2sup,
                        &mut pos2a,
                        &mut pos2b,
                    );
                    let _ = (u2inf, u2sup);

                    // ####################################################################
                    // ##  Test sur la longueur minimale d un segment sur Ligne2
                    // ####################################################################
                    let res2sup_m_res2inf = res2sup - res2inf;
                    if res2sup_m_res2inf < 0.0 {
                        //-- Pas de solutions On retourne Vide
                    } else if (res2sup_m_res2inf > long_mini_seg)
                        || ((pos2a == pos2b) && (pos2a != Position::Middle))
                    {
                        //----------- Calcul des attributs du segment --------------
                        //-- Attention, les bornes Res1inf(sup) bougent donc il faut
                        //--  eventuellement recalculer les attributs

                        if is_opposite {
                            res1inf = u1pu2 - res2sup;
                            res1sup = u1pu2 - res2inf;
                            let tampon = res2inf;
                            res2inf = res2sup;
                            res2sup = tampon;
                            let pos = pos2a;
                            pos2a = pos2b;
                            pos2b = pos;
                        } else {
                            res1inf = u1mu2 + res2inf;
                            res1sup = u1mu2 + res2sup;
                        }

                        pos1a = find_position_ll(&mut res1inf, domain1);
                        pos1b = find_position_ll(&mut res1sup, domain1);

                        let mut t1a = Transition::empty();
                        let mut t2a = Transition::empty();
                        let mut t1b = Transition::empty();
                        let mut t2b = Transition::empty();

                        if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                            // &&&&&&&&&&&&&&&
                            t1a.set_value_in_out(false, pos1a, TypeTrans::Out);
                            t2a.set_value_in_out(false, pos2a, TypeTrans::In);
                        } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                            t1a.set_value_in_out(false, pos1a, TypeTrans::In);
                            t2a.set_value_in_out(false, pos2a, TypeTrans::Out);
                        } else {
                            t1a.set_value_touch(false, pos1a, Situation::Unknown, is_opposite);
                            t2a.set_value_touch(false, pos2a, Situation::Unknown, is_opposite);
                        }

                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        //~~~~~~~  C O N V E N T I O N    -    S E G M E N T     ~~~~~~~
                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        //~~ On Renvoie un segment dans les cas suivants :            ~~
                        //~~   (1) Extremite L1 L2   ------>    Extremite L1 L2       ~~
                        //~~   (2) Extremite L1 L2   ------>    Intersection          ~~
                        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                        let mut result_is_a_point = false;

                        if ((res1sup - res1inf) <= long_mini_seg)
                            || ((res2sup - res2inf).abs() <= long_mini_seg)
                        {
                            //-- On force la creation d un point
                            result_is_a_point = true;
                        } else {
                            //------------------------------------------------------------
                            //-- On traite les cas ou l intersection est situee du
                            //-- Mauvais cote du domaine
                            //-- Attention : Res2inf <-> Pos2a        Res2sup <-> Pos2b
                            //--  et         Res1inf <-> Pos1a        Res1sup <-> Pos1b
                            //--             avec Res1inf <= Res1sup
                            //------------------------------------------------------------
                            //-- Le point sera : Res1inf,Res2inf,T1a(Pos1a),T2a(Pos2a)
                            //------------------------------------------------------------

                            if pos1a == Position::Head {
                                if pos1b != Position::End && u1 < res1inf {
                                    result_is_a_point = true;
                                    u1 = res1inf;
                                    u2 = res2inf;
                                }
                            }
                            if pos1b == Position::End {
                                if pos1a != Position::Head && u1 > res1sup {
                                    result_is_a_point = true;
                                    u1 = res1sup;
                                    u2 = res2sup;
                                }
                            }

                            if pos2a == Position::Head {
                                if pos2b != Position::End && u2 < res2inf {
                                    result_is_a_point = true;
                                    u2 = res2inf;
                                    u1 = res1inf;
                                }
                            } else if pos2a == Position::End {
                                if pos2b != Position::Head && u2 > res2inf {
                                    result_is_a_point = true;
                                    u2 = res2inf;
                                    u1 = res1inf;
                                }
                            }
                            if pos2b == Position::Head {
                                if pos2a != Position::End && u2 < res2sup {
                                    result_is_a_point = true;
                                    u2 = res2sup;
                                    u1 = res1sup;
                                }
                            } else if pos2b == Position::End {
                                if pos2a != Position::Head && u2 > res2sup {
                                    result_is_a_point = true;
                                    u2 = res2sup;
                                    u1 = res1sup;
                                }
                            }
                        }

                        if (!result_is_a_point) && (pos1a != Position::Middle || pos2a != Position::Middle)
                        {
                            if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                //&&&&&&&&&&&&&&
                                t1b.set_value_in_out(false, pos1b, TypeTrans::Out);
                                t2b.set_value_in_out(false, pos2b, TypeTrans::In);
                            } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                t1b.set_value_in_out(false, pos1b, TypeTrans::In);
                                t2b.set_value_in_out(false, pos2b, TypeTrans::Out);
                            } else {
                                t1b.set_value_touch(false, pos1b, Situation::Unknown, is_opposite);
                                t2b.set_value_touch(false, pos2b, Situation::Unknown, is_opposite);
                            }
                            let ptdebut;
                            if pos1a == Position::Middle {
                                let t3 = if is_opposite {
                                    if pos2a == Position::Head {
                                        res2sup
                                    } else {
                                        res2inf
                                    }
                                } else if pos2a == Position::Head {
                                    res2inf
                                } else {
                                    res2sup
                                };
                                ptdebut = elclib2d::line_value(l2.origin, l2.direction, t3);
                                res1inf = elclib2d::line_parameter(l1.origin, l1.direction, ptdebut);
                            } else {
                                let t4 = if pos1a == Position::Head {
                                    res1inf
                                } else {
                                    res1sup
                                };
                                ptdebut = elclib2d::line_value(l1.origin, l1.direction, t4);
                                res2inf = elclib2d::line_parameter(l2.origin, l2.direction, ptdebut);
                            }
                            pt_seg1.set_values(ptdebut, res1inf, res2inf, t1a, t2a, false);
                            if pos1b != Position::Middle || pos2b != Position::Middle {
                                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                //~~ Ajustement des parametres et du point renvoye
                                let ptfin;
                                if pos1b == Position::Middle {
                                    ptfin = elclib2d::line_value(l2.origin, l2.direction, res2sup);
                                    res1inf = elclib2d::line_parameter(l1.origin, l1.direction, ptfin);
                                } else {
                                    ptfin = elclib2d::line_value(l1.origin, l1.direction, res1sup);
                                    res2inf = elclib2d::line_parameter(l2.origin, l2.direction, ptfin);
                                }
                                pt_seg2.set_values(ptfin, res1sup, res2sup, t1b, t2b, false);
                                let segment =
                                    IntersectionSegment::with_points(&pt_seg1, &pt_seg2, is_opposite, false);
                                self.base.append_segment(&segment);
                            } else {
                                //-- Extremite(L1 ou L2)  ------>   Point Middle(L1 et L2)

                                pos1b = find_position_ll(&mut u1, domain1);
                                pos2b = find_position_ll(&mut u2, domain2);
                                if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                    t1b.set_value_in_out(false, pos1b, TypeTrans::Out);
                                    t2b.set_value_in_out(false, pos2b, TypeTrans::In);
                                } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                    t1b.set_value_in_out(false, pos1b, TypeTrans::In);
                                    t2b.set_value_in_out(false, pos2b, TypeTrans::Out);
                                } else {
                                    t1b.set_value_touch(false, pos1b, Situation::Unknown, is_opposite);
                                    t2b.set_value_touch(false, pos2b, Situation::Unknown, is_opposite);
                                }

                                pt_seg2.set_values(
                                    elclib2d::line_value(l2.origin, l2.direction, u2),
                                    u1,
                                    u2,
                                    t1b,
                                    t2b,
                                    false,
                                );

                                if ((res1inf - u1).abs() > long_mini_seg)
                                    && ((res2inf - u2).abs() > long_mini_seg)
                                {
                                    let segment = IntersectionSegment::with_points(
                                        &pt_seg1,
                                        &pt_seg2,
                                        is_opposite,
                                        false,
                                    );
                                    self.base.append_segment(&segment);
                                } else {
                                    let p = segment_to_point(&pt_seg1, &t1a, &t2a, &pt_seg2, &t1b, &t2b);
                                    self.base.append_point(&p);
                                }
                            }
                        } //-- (Pos1a!=IntRes2d_Middle || Pos2a!=IntRes2d_Middle) --
                        else {
                            //-- Pos1a == Pos2a == Middle
                            if pos1b == Position::Middle {
                                pos1b = pos1a;
                            }
                            if pos2b == Position::Middle {
                                pos2b = pos2a;
                            }
                            if result_is_a_point {
                                //-- Middle sur le segment A
                                //--
                                if pos1b != Position::Middle || pos2b != Position::Middle {
                                    let ptfin;
                                    if pos1b == Position::Middle {
                                        let t2 = if is_opposite {
                                            if pos2b == Position::Head {
                                                res2sup
                                            } else {
                                                res2inf
                                            }
                                        } else if pos2b == Position::Head {
                                            res2inf
                                        } else {
                                            res2sup
                                        };
                                        ptfin = elclib2d::line_value(l2.origin, l2.direction, t2);
                                        res1sup = elclib2d::line_parameter(l1.origin, l1.direction, ptfin);
                                        // modified by NIZHNY-MKK  Tue Feb 15 10:54:51 2000.BEGIN
                                        pos1b = find_position_ll(&mut res1sup, domain1);
                                        // modified by NIZHNY-MKK  Tue Feb 15 10:54:55 2000.END
                                    } else {
                                        let t1 = if pos1b == Position::Head {
                                            res1inf
                                        } else {
                                            res1sup
                                        };
                                        ptfin = elclib2d::line_value(l1.origin, l1.direction, t1);
                                        res2sup = elclib2d::line_parameter(l2.origin, l2.direction, ptfin);
                                        // modified by NIZHNY-MKK  Tue Feb 15 10:55:08 2000.BEGIN
                                        pos2b = find_position_ll(&mut res2sup, domain2);
                                        // modified by NIZHNY-MKK  Tue Feb 15 10:55:11 2000.END
                                    }
                                    if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                        t1b.set_value_in_out(false, pos1b, TypeTrans::Out);
                                        t2b.set_value_in_out(false, pos2b, TypeTrans::In);
                                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                        t1b.set_value_in_out(false, pos1b, TypeTrans::In);
                                        t2b.set_value_in_out(false, pos2b, TypeTrans::Out);
                                    } else {
                                        t1b.set_value_touch(false, pos1b, Situation::Unknown, is_opposite);
                                        t2b.set_value_touch(false, pos2b, Situation::Unknown, is_opposite);
                                    }
                                    pt_seg2.set_values(ptfin, res1sup, res2sup, t1b, t2b, false);
                                    self.base.append_point(&pt_seg2);
                                } else {
                                    pos1b = find_position_ll(&mut u1, domain1);
                                    pos2b = find_position_ll(&mut u2, domain2);

                                    if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                        t1b.set_value_in_out(false, pos1b, TypeTrans::Out);
                                        t2b.set_value_in_out(false, pos2b, TypeTrans::In);
                                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                        t1b.set_value_in_out(false, pos1b, TypeTrans::In);
                                        t2b.set_value_in_out(false, pos2b, TypeTrans::Out);
                                    } else {
                                        t1b.set_value_touch(false, pos1b, Situation::Unknown, is_opposite);
                                        t2b.set_value_touch(false, pos2b, Situation::Unknown, is_opposite);
                                    }
                                    pt_seg1.set_values(
                                        elclib2d::line_value(l2.origin, l2.direction, u2),
                                        u1,
                                        u2,
                                        t1b,
                                        t2b,
                                        false,
                                    );
                                    self.base.append_point(&pt_seg1);
                                }
                            } else {
                                pt_seg1.set_values(
                                    elclib2d::line_value(l2.origin, l2.direction, u2),
                                    u1,
                                    u2,
                                    t1a,
                                    t2a,
                                    false,
                                );

                                if pos1b != Position::Middle || pos2b != Position::Middle {
                                    if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                                        t1b.set_value_in_out(false, pos1b, TypeTrans::Out);
                                        t2b.set_value_in_out(false, pos2b, TypeTrans::In);
                                    } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                                        t1b.set_value_in_out(false, pos1b, TypeTrans::In);
                                        t2b.set_value_in_out(false, pos2b, TypeTrans::Out);
                                    } else {
                                        t1b.set_value_touch(false, pos1b, Situation::Unknown, is_opposite);
                                        t2b.set_value_touch(false, pos2b, Situation::Unknown, is_opposite);
                                    }
                                    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    //~~ Ajustement des parametres et du point renvoye
                                    let ptfin;
                                    if pos1b == Position::Middle {
                                        ptfin = elclib2d::line_value(l2.origin, l2.direction, res2sup);
                                        res1sup = elclib2d::line_parameter(l1.origin, l1.direction, ptfin);
                                    } else {
                                        ptfin = elclib2d::line_value(l1.origin, l1.direction, res1sup);
                                        res2sup = elclib2d::line_parameter(l2.origin, l2.direction, ptfin);
                                    }

                                    pt_seg2.set_values(ptfin, res1sup, res2sup, t1b, t2b, false);

                                    if ((u1 - res1sup).abs() > long_mini_seg)
                                        || ((u2 - res2sup).abs() > long_mini_seg)
                                    {
                                        //-- Modif du 1er Octobre 92 (Pour Composites)

                                        let segment = IntersectionSegment::with_points(
                                            &pt_seg1,
                                            &pt_seg2,
                                            is_opposite,
                                            false,
                                        );
                                        self.base.append_segment(&segment);
                                    } else {
                                        let p = segment_to_point(&pt_seg1, &t1a, &t2a, &pt_seg2, &t1b, &t2b);
                                        self.base.append_point(&p);
                                    }
                                } else {
                                    self.base.append_point(&pt_seg1);
                                }
                            }
                        }
                    } //----- Fin Creation Segment ----(Res2sup-Res2inf>Tol)-------------
                    else {
                        //------ (Intersection And Domain1)  AND  Domain2  --> Point ------
                        //-- Attention Res1sup peut etre  different de  U2
                        //--   Mais on a Res1sup-Res1inf < Tol

                        // gka #0022833
                        let a_cur_trans = if prod_vect_tan >= TOLERANCE_ANGULAIRE {
                            TypeTrans::In
                        } else if prod_vect_tan <= -TOLERANCE_ANGULAIRE {
                            TypeTrans::Out
                        } else {
                            TypeTrans::Undecided
                        };

                        let mut new_point1 = IntersectionPoint::empty();
                        if compute_int_point(
                            domain2,
                            domain1,
                            l2,
                            l1,
                            a_cos_t1_t2,
                            u2,
                            u1,
                            &mut res2inf,
                            &mut res2sup,
                            2,
                            a_cur_trans,
                            &mut new_point1,
                        ) {
                            self.base.append_point(&new_point1);
                        }
                    }
                }
            }

            // #ifdef OCCT_DEBUG (the printf trace) is an OCCT debug block.
        } else if nbsol == 2 {
            //== Droites confondues a la tolerance pres
            //--On traite ici le cas de segments resultats non neccess. bornes
            //--
            //--On prend la droite D1 comme reference ( pour le sens positif )
            //--
            let mut res_has_first_point = 0i32;
            let mut res_has_last_point = 0i32;
            let mut param_start;
            let mut param_start2;
            let mut param_end;
            let mut param_end2;
            let org2_sur_l1 = elclib2d::line_parameter(l1.origin, l1.direction, l2.origin);
            //== 3 : L1 et L2 bornent
            //== 2 :       L2 borne
            //== 1 : L1 borne
            if domain1.has_first_point() {
                res_has_first_point = 1;
            }
            if domain1.has_last_point() {
                res_has_last_point = 1;
            }
            if is_opposite {
                if domain2.has_last_point() {
                    res_has_first_point += 2;
                }
                if domain2.has_first_point() {
                    res_has_last_point += 2;
                }
            } else {
                if domain2.has_last_point() {
                    res_has_last_point += 2;
                }
                if domain2.has_first_point() {
                    res_has_first_point += 2;
                }
            }
            if res_has_first_point == 0 && res_has_last_point == 0 {
                //~~~~ Creation d un segment infini avec Opposite
                self.base.append_segment(&IntersectionSegment::infinite(is_opposite));
            } else {
                //-- On obtient au pire une demi-droite
                match res_has_first_point {
                    1 => {
                        param_start = domain1.first_parameter();
                        param_start2 = if is_opposite {
                            org2_sur_l1 - param_start
                        } else {
                            param_start - org2_sur_l1
                        };
                    }
                    2 => {
                        if is_opposite {
                            param_start2 = domain2.last_parameter();
                            param_start = org2_sur_l1 - param_start2;
                        } else {
                            param_start2 = domain2.first_parameter();
                            param_start = org2_sur_l1 + param_start2;
                        }
                    }
                    3 => {
                        if is_opposite {
                            param_start2 = domain2.last_parameter();
                            param_start = org2_sur_l1 - param_start2;
                            if param_start < domain1.first_parameter() {
                                param_start = domain1.first_parameter();
                                param_start2 = org2_sur_l1 - param_start;
                            }
                        } else {
                            param_start2 = domain2.first_parameter();
                            param_start = org2_sur_l1 + param_start2;
                            if param_start < domain1.first_parameter() {
                                param_start = domain1.first_parameter();
                                param_start2 = param_start - org2_sur_l1;
                            }
                        }
                    }
                    _ => {
                        //~~~ Segment Infini a gauche
                        param_start = 0.0;
                        param_start2 = 0.0;
                    }
                }

                match res_has_last_point {
                    1 => {
                        param_end = domain1.last_parameter();
                        param_end2 = if is_opposite {
                            org2_sur_l1 - param_end
                        } else {
                            param_end - org2_sur_l1
                        };
                    }
                    2 => {
                        if is_opposite {
                            param_end2 = domain2.first_parameter();
                            param_end = org2_sur_l1 - param_end2;
                        } else {
                            param_end2 = domain2.last_parameter();
                            param_end = org2_sur_l1 + param_end2;
                        }
                    }
                    3 => {
                        if is_opposite {
                            param_end2 = domain2.first_parameter();
                            param_end = org2_sur_l1 - param_end2;
                            if param_end > domain1.last_parameter() {
                                param_end = domain1.last_parameter();
                                param_end2 = org2_sur_l1 - param_end;
                            }
                        } else {
                            param_end2 = domain2.last_parameter();
                            param_end = org2_sur_l1 + param_end2;
                            if param_end > domain1.last_parameter() {
                                param_end = domain1.last_parameter();
                                param_end2 = param_end - org2_sur_l1;
                            }
                        }
                    }
                    _ => {
                        //~~~ Segment Infini a droite
                        param_end = 0.0;
                        param_end2 = 0.0;
                    }
                }

                let mut tinf = Transition::empty();
                let mut tsup = Transition::empty();

                if res_has_first_point != 0 {
                    if res_has_last_point != 0 {
                        //~~~ Creation de la borne superieure
                        //~~~ L1 :     |------------->       ou          |-------------->
                        //~~~ L2 : <------------|            ou  <----|
                        if param_end >= (param_start - tol) {
                            //~~~ Creation d un segment
                            let mut pos1;
                            let mut pos2;
                            pos1 = find_position_ll(&mut param_start, domain1);
                            pos2 = find_position_ll(&mut param_start2, domain2);
                            tinf.set_value_touch(true, pos1, Situation::Unknown, is_opposite);
                            tsup.set_value_touch(true, pos2, Situation::Unknown, is_opposite);
                            let p1 = IntersectionPoint::new(
                                elclib2d::line_value(l1.origin, l1.direction, param_start),
                                param_start,
                                param_start2,
                                tinf,
                                tsup,
                                false,
                            );
                            if param_end > (param_start + tol) {
                                //~~~ Le segment est assez long
                                pos1 = find_position_ll(&mut param_end, domain1);
                                pos2 = find_position_ll(&mut param_end2, domain2);
                                tinf.set_value_touch(true, pos1, Situation::Unknown, is_opposite);
                                tsup.set_value_touch(true, pos2, Situation::Unknown, is_opposite);

                                let p2 = IntersectionPoint::new(
                                    elclib2d::line_value(l1.origin, l1.direction, param_end),
                                    param_end,
                                    param_end2,
                                    tinf,
                                    tsup,
                                    false,
                                );
                                let seg = IntersectionSegment::with_points(&p1, &p2, is_opposite, false);
                                self.base.append_segment(&seg);
                            } else {
                                //~~~~ le segment est de longueur inferieure a Tol
                                self.base.append_point(&p1);
                            }
                        } //-- if( ParamEnd >= ...)
                    } else {
                        //~~~ Creation de la demi droite   |----------->
                        let mut pos1 = find_position_ll(&mut param_start, domain1);
                        let mut pos2 = find_position_ll(&mut param_start2, domain2);
                        tinf.set_value_touch(true, pos1, Situation::Unknown, is_opposite);
                        tsup.set_value_touch(true, pos2, Situation::Unknown, is_opposite);

                        let p = IntersectionPoint::new(
                            elclib2d::line_value(l1.origin, l1.direction, param_start),
                            param_start,
                            param_start2,
                            tinf,
                            tsup,
                            false,
                        );
                        let seg = IntersectionSegment::with_one_point(&p, true, is_opposite, false);
                        self.base.append_segment(&seg);
                        let _ = (&mut pos1, &mut pos2);
                    }
                } else {
                    let mut pos1 = find_position_ll(&mut param_end, domain1);
                    let mut pos2 = find_position_ll(&mut param_end2, domain2);
                    tinf.set_value_touch(true, pos1, Situation::Unknown, is_opposite);
                    tsup.set_value_touch(true, pos2, Situation::Unknown, is_opposite);

                    let p2 = IntersectionPoint::new(
                        elclib2d::line_value(l1.origin, l1.direction, param_end),
                        param_end,
                        param_end2,
                        tinf,
                        tsup,
                        false,
                    );
                    let seg = IntersectionSegment::with_one_point(&p2, false, is_opposite, false);
                    self.base.append_segment(&seg);
                    let _ = (&mut pos1, &mut pos2);
                    //~~~ Creation de la demi droite   <-----------|
                }
            }
        }
    }

    // OCCT Perform(const gp_Lin2d& L, const gp_Circ2d& C) lives in
    // `int_conic_conic_lin_circ` (the _1.cxx Lin-Circ section).

    /// OCCT Perform(const gp_Lin2d& L, const gp_Elips2d& E)
    /// (IntCurve_IntConicConic.hxx L112 decl, body in _1.cxx — dedicated
    /// closed-form implementation).
    pub fn perform_line_ellipse(
        &mut self,
        _l: &Line2d,
        _dl: &Res2dDomain,
        _e: &Ellipse2d,
        _de: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Lin, Elips) — not ported yet");
    }

    /// OCCT Perform(const gp_Lin2d& L, const gp_Parab2d& P)
    /// (IntCurve_IntConicConic.cxx L109-226).
    pub fn perform_line_parabola(
        &mut self,
        _l: &Line2d,
        _dl: &Res2dDomain,
        _p: &Parabola2d,
        _dp: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Lin, Parab) — IntCurve_IntConicConic.cxx L109, not ported yet");
    }

    /// OCCT Perform(const gp_Lin2d& L, const gp_Hypr2d& H)
    /// (IntCurve_IntConicConic.cxx L230-333).
    pub fn perform_line_hyperbola(
        &mut self,
        _l: &Line2d,
        _dl: &Res2dDomain,
        _h: &Hyperbola2d,
        _dh: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Lin, Hypr) — IntCurve_IntConicConic.cxx L230, not ported yet");
    }

    /// OCCT Perform(const gp_Circ2d& C1, const gp_Circ2d& C2)
    /// (IntCurve_IntConicConic_1.cxx L807-1238).
    pub fn perform_circle_circle(
        &mut self,
        _c1: &Circle2d,
        _d1: &Res2dDomain,
        _c2: &Circle2d,
        _d2: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Circ, Circ) — IntCurve_IntConicConic_1.cxx L807, not ported yet");
    }

    /// OCCT Perform(const gp_Circ2d& C, const gp_Elips2d& E)
    /// (IntCurve_IntConicConic.cxx L439-482).
    pub fn perform_circle_ellipse(
        &mut self,
        _c: &Circle2d,
        _dc: &Res2dDomain,
        _e: &Ellipse2d,
        _de: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Circ, Elips) — IntCurve_IntConicConic.cxx L439, not ported yet");
    }

    /// OCCT Perform(const gp_Circ2d& C, const gp_Parab2d& P)
    /// (IntCurve_IntConicConic.cxx L337-435).
    pub fn perform_circle_parabola(
        &mut self,
        _c: &Circle2d,
        _dc: &Res2dDomain,
        _p: &Parabola2d,
        _dp: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Circ, Parab) — IntCurve_IntConicConic.cxx L337, not ported yet");
    }

    /// OCCT Perform(const gp_Circ2d& C, const gp_Hypr2d& H)
    /// (IntCurve_IntConicConic.cxx L486-581).
    pub fn perform_circle_hyperbola(
        &mut self,
        _c: &Circle2d,
        _dc: &Res2dDomain,
        _h: &Hyperbola2d,
        _dh: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Circ, Hypr) — IntCurve_IntConicConic.cxx L486, not ported yet");
    }

    /// OCCT Perform(const gp_Elips2d& E, const gp_Parab2d& P)
    /// (IntCurve_IntConicConic.cxx L692-806).
    pub fn perform_ellipse_parabola(
        &mut self,
        _e: &Ellipse2d,
        _de: &Res2dDomain,
        _p: &Parabola2d,
        _dp: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Elips, Parab) — IntCurve_IntConicConic.cxx L692, not ported yet");
    }

    /// OCCT Perform(const gp_Elips2d& E, const gp_Hypr2d& H)
    /// (IntCurve_IntConicConic.cxx L962-1063).
    pub fn perform_ellipse_hyperbola(
        &mut self,
        _e: &Ellipse2d,
        _de: &Res2dDomain,
        _h: &Hyperbola2d,
        _dh: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Elips, Hypr) — IntCurve_IntConicConic.cxx L962, not ported yet");
    }

    /// OCCT Perform(const gp_Parab2d& P1, const gp_Parab2d& P2)
    /// (IntCurve_IntConicConic.cxx L585-688).
    pub fn perform_parabola_parabola(
        &mut self,
        _p1: &Parabola2d,
        _d1: &Res2dDomain,
        _p2: &Parabola2d,
        _d2: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Parab, Parab) — IntCurve_IntConicConic.cxx L585, not ported yet");
    }

    /// OCCT Perform(const gp_Parab2d& P, const gp_Hypr2d& H)
    /// (IntCurve_IntConicConic.cxx L810-911).
    pub fn perform_parabola_hyperbola(
        &mut self,
        _p: &Parabola2d,
        _dp: &Res2dDomain,
        _h: &Hyperbola2d,
        _dh: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Parab, Hypr) — IntCurve_IntConicConic.cxx L810, not ported yet");
    }

    /// OCCT Perform(const gp_Hypr2d& H1, const gp_Hypr2d& H2)
    /// (IntCurve_IntConicConic.cxx L1067-1168).
    pub fn perform_hyperbola_hyperbola(
        &mut self,
        _h1: &Hyperbola2d,
        _d1: &Res2dDomain,
        _h2: &Hyperbola2d,
        _d2: &Res2dDomain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        unimplemented!("IntConicConic::Perform(Hypr, Hypr) — IntCurve_IntConicConic.cxx L1067, not ported yet");
    }
}

impl Default for IntConicConic {
    fn default() -> Self {
        IntConicConic::new()
    }
}

// ---------------------------------------------------------------------------
// The file-static helpers of IntCurve_IntConicConic_1.cxx
// ---------------------------------------------------------------------------

/// OCCT DomainIntersection (_1.cxx L641-727).
#[allow(clippy::too_many_arguments)]
fn domain_intersection(
    domain: &Res2dDomain,
    u1inf: f64,
    u1sup: f64,
    res1inf: &mut f64,
    res1sup: &mut f64,
    pos_inf: &mut Position,
    pos_sup: &mut Position,
) {
    if domain.has_first_point() {
        if u1sup < (domain.first_parameter() - domain.first_tolerance()) {
            *res1inf = 1.0;
            *res1sup = -1.0;
            return;
        }
        if u1inf > (domain.first_parameter() + domain.first_tolerance()) {
            *res1inf = u1inf;
            *pos_inf = Position::Middle;
        } else {
            *res1inf = domain.first_parameter();
            *pos_inf = Position::Head;
        }
    } else {
        *res1inf = u1inf;
        *pos_inf = Position::Middle;
    }

    if domain.has_last_point() {
        if u1inf > (domain.last_parameter() + domain.last_tolerance()) {
            *res1inf = 1.0;
            *res1sup = -1.0;
            return;
        }
        if u1sup < (domain.last_parameter() - domain.last_tolerance()) {
            *res1sup = u1sup;
            *pos_sup = Position::Middle;
        } else {
            *res1sup = domain.last_parameter();
            *pos_sup = Position::End;
        }
    } else {
        *res1sup = u1sup;
        *pos_sup = Position::Middle;
    }
    //-- Si un des points est en bout ,
    //-- on s assure que les parametres sont corrects
    if *res1inf > *res1sup {
        if *pos_sup == Position::Middle {
            *res1sup = *res1inf;
        } else {
            *res1inf = *res1sup;
        }
    }
    //--- Traitement des cas ou une intersection vraie est dans la tolerance
    //--  d un des bouts (the OCCT commented-out block).
}

/// OCCT LineLineGeometricIntersection (_1.cxx L730-775).
fn line_line_geometric_intersection(
    l1: &Line2d,
    l2: &Line2d,
    tol: f64,
    u1: &mut f64,
    u2: &mut f64,
    sin_demi_angle: &mut f64,
    nbsol: &mut i32,
) {
    let u1x = l1.direction.x;
    let u1y = l1.direction.y;
    let u2x = l2.direction.x;
    let u2y = l2.direction.y;
    let uo21x = l2.origin.x - l1.origin.x;
    let uo21y = l2.origin.y - l1.origin.y;

    let mut d = u1y * u2x - u1x * u2y;

    // modified by NIZHNY-MKK  Tue Feb 15 10:54:04 2000.BEGIN
    //    if(std::abs(D)<1e-15) { //-- Droites //
    if d.abs() < TOLERANCE_ANGULAIRE {
        // modified by NIZHNY-MKK  Tue Feb 15 10:54:11 2000.END
        d = u1y * uo21x - u1x * uo21y;
        *nbsol = if d.abs() <= tol { 2 } else { 0 };
    } else {
        *u1 = (uo21y * u2x - uo21x * u2y) / d;
        *u2 = (uo21y * u1x - uo21x * u1y) / d;
        //------------------- Calcul du Sin du demi angle  entre L1 et L2
        //----
        if d < 0.0 {
            d = -d;
        }
        if d > 1.0 {
            d = 1.0; //-- Deja vu !
        }
        *sin_demi_angle = (0.5 * d.asin()).sin();
        *nbsol = 1;
    }
}

/// OCCT FindPositionLL (_1.cxx L1209-1237).
fn find_position_ll(param: &mut f64, domain: &Res2dDomain) -> Position {
    let mut a_dpar = rcad_kernel::precision::INFINITE_VALUE;
    let mut a_pos = Position::Middle;
    let mut a_res_par = *param;
    if domain.has_first_point() {
        a_dpar = (*param - domain.first_parameter()).abs();
        if a_dpar <= domain.first_tolerance() {
            a_res_par = domain.first_parameter();
            a_pos = Position::Head;
        }
    }
    if domain.has_last_point() {
        let a_d2 = (*param - domain.last_parameter()).abs();
        if a_d2 <= domain.last_tolerance() && (a_pos == Position::Middle || a_d2 < a_dpar) {
            a_res_par = domain.last_parameter();
            a_pos = Position::End;
        }
    }
    *param = a_res_par;
    a_pos
}

/// OCCT getDomainParametrs (_1.cxx L1240-1252) — gka 0022833.
fn get_domain_parameters(
    the_domain: &Res2dDomain,
    the_first: &mut f64,
    the_last: &mut f64,
    the_tol1: &mut f64,
    the_tol2: &mut f64,
) {
    *the_first = if the_domain.has_first_point() {
        the_domain.first_parameter()
    } else {
        -rcad_kernel::precision::INFINITE_VALUE
    };
    *the_last = if the_domain.has_last_point() {
        the_domain.last_parameter()
    } else {
        rcad_kernel::precision::INFINITE_VALUE
    };
    *the_tol1 = if the_domain.has_first_point() {
        the_domain.first_tolerance()
    } else {
        0.0
    };
    *the_tol2 = if the_domain.has_last_point() {
        the_domain.last_tolerance()
    } else {
        0.0
    };
}

/// OCCT computeIntPoint (_1.cxx L1254-1356) — the intersection point for the
/// case when the specified domain is less than the intersection tolerance
/// (gka 0022833).
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_arguments)]
fn compute_int_point(
    the_cur_domain: &Res2dDomain,
    the_domain_other: &Res2dDomain,
    the_cur_lin: &Line2d,
    the_other_lin: &Line2d,
    the_cos_t1_t2: f64,
    the_par_cur: f64,
    the_par_other: f64,
    the_res_inf: &mut f64,
    the_res_sup: &mut f64,
    the_num: i32,
    the_cur_trans: TypeTrans,
    the_new_point: &mut IntersectionPoint,
) -> bool {
    if (*the_res_sup - the_par_cur).abs() > (*the_res_inf - the_par_cur).abs() {
        *the_res_sup = *the_res_inf;
    }

    let a_res2 = the_par_other + (*the_res_sup - the_par_cur) * the_cos_t1_t2;

    let mut a_first2 = 0.0;
    let mut a_last2 = 0.0;
    let mut a_tol21 = 0.0;
    let mut a_tol22 = 0.0;
    let mut a_tol11 = 0.0;
    let mut a_tol12 = 0.0;

    get_domain_parameters(
        the_domain_other,
        &mut a_first2,
        &mut a_last2,
        &mut a_tol21,
        &mut a_tol22,
    );

    if a_res2 < a_first2 - a_tol21 || a_res2 > a_last2 + a_tol22 {
        return false;
    }

    //------ compute parameters of intersection point --
    let mut a_t1 = Transition::empty();
    let mut a_t2 = Transition::empty();
    let mut res_sup = *the_res_sup;
    let a_pos1a = find_position_ll(&mut res_sup, the_cur_domain);
    *the_res_sup = res_sup;
    let mut a_res2m = a_res2;
    let a_pos2a = find_position_ll(&mut a_res2m, the_domain_other);
    let an_other_trans = if the_cur_trans == TypeTrans::Out {
        TypeTrans::In
    } else if the_cur_trans == TypeTrans::In {
        TypeTrans::Out
    } else {
        TypeTrans::Undecided
    };

    if the_cur_trans != TypeTrans::Undecided {
        a_t1.set_value_in_out(false, a_pos1a, the_cur_trans);
        a_t2.set_value_in_out(false, a_pos2a, an_other_trans);
    } else {
        let an_opposite = the_cos_t1_t2 < 0.0;
        a_t1.set_value_touch(false, a_pos1a, Situation::Unknown, an_opposite);
        a_t2.set_value_touch(false, a_pos2a, Situation::Unknown, an_opposite);
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //--------------------------------------------------
    // gka bug 0022833
    let mut a_res_u1 = the_par_cur;
    let mut a_res_u2 = the_par_other;

    let mut a_first1 = 0.0;
    let mut a_last1 = 0.0;
    get_domain_parameters(the_cur_domain, &mut a_first1, &mut a_last1, &mut a_tol11, &mut a_tol12);

    let is_inside1 = the_par_cur >= a_first1 && the_par_cur <= a_last1;
    let is_inside2 = the_par_other >= a_first2 && the_par_other <= a_last2;

    if !is_inside1 || !is_inside2 {
        if is_inside1 {
            let pt1 = elclib2d::line_value(the_other_lin.origin, the_other_lin.direction, a_res2);
            a_res_u2 = a_res2;
            let a_par1 = elclib2d::line_parameter(the_cur_lin.origin, the_cur_lin.direction, pt1);
            a_res_u1 = if a_par1 >= a_first1 && a_par1 <= a_last1 {
                a_par1
            } else {
                *the_res_sup
            };
        } else if is_inside2 {
            let a_pt1 =
                elclib2d::line_value(the_cur_lin.origin, the_cur_lin.direction, *the_res_sup);
            a_res_u1 = *the_res_sup;
            let a_par2 = elclib2d::line_parameter(the_other_lin.origin, the_other_lin.direction, a_pt1);
            a_res_u2 = if a_par2 >= a_first2 && a_par2 <= a_last2 {
                a_par2
            } else {
                a_res2
            };
        } else {
            // PKVf
            //  check that parameters are within range on both curves
            if the_par_cur < a_first1 - a_tol11
                || the_par_cur > a_last1 + a_tol12
                || the_par_other < a_first2 - a_tol21
                || the_par_other > a_last2 + a_tol22
            {
                return false;
            }
            // PKVt
            a_res_u1 = *the_res_sup;
            a_res_u2 = a_res2;
        }
    }
    let a_pres = (elclib2d::line_value(the_cur_lin.origin, the_cur_lin.direction, a_res_u1)
        + elclib2d::line_value(the_other_lin.origin, the_other_lin.direction, a_res_u2))
        * 0.5;
    if the_num == 1 {
        the_new_point.set_values(a_pres, a_res_u1, a_res_u2, a_t1, a_t2, false);
    } else {
        the_new_point.set_values(a_pres, a_res_u2, a_res_u1, a_t2, a_t1, false);
    }
    true
}

/// OCCT CheckLLCoincidence (_1.cxx L1363-1377) — returns true if the input
/// are trimmed curves and they coincide within tolerance.
fn check_ll_coincidence(
    l1: &Line2d,
    l2: &Line2d,
    domain1: &Res2dDomain,
    domain2: &Res2dDomain,
    the_tol: f64,
) -> bool {
    let is_first1 = domain1.has_first_point() && l2.distance(domain1.first_point()) < the_tol;
    let is_last1 = domain1.has_last_point() && l2.distance(domain1.last_point()) < the_tol;
    if is_first1 && is_last1 {
        return true;
    }
    let is_first2 = domain2.has_first_point() && l1.distance(domain2.first_point()) < the_tol;
    let is_last2 = domain2.has_last_point() && l1.distance(domain2.last_point()) < the_tol;
    is_first2 && is_last2
}

/// OCCT SegmentToPoint (_1.cxx L2620-2649).
#[allow(clippy::too_many_arguments)]
fn segment_to_point(
    pa: &IntersectionPoint,
    t1a: &Transition,
    t2a: &Transition,
    pb: &IntersectionPoint,
    t1b: &Transition,
    t2b: &Transition,
) -> IntersectionPoint {
    if (t1b.position_on_curve() == Position::Middle) && (t2b.position_on_curve() == Position::Middle)
    {
        return pa.clone();
    }
    if (t1a.position_on_curve() == Position::Middle) && (t2a.position_on_curve() == Position::Middle)
    {
        return pb.clone();
    }

    let mut t1 = *t1a;
    let mut t2 = *t2a;
    let mut u1 = pa.param_on_first();
    let mut u2 = pa.param_on_second();

    if t1.position_on_curve() == Position::Middle {
        t1.set_position(t1b.position_on_curve());
        u1 = pb.param_on_first();
    }
    if t2.position_on_curve() == Position::Middle {
        t2.set_position(t2b.position_on_curve());
        u2 = pb.param_on_second();
    }
    IntersectionPoint::new(pa.value(), u1, u2, t1, t2, false)
}
