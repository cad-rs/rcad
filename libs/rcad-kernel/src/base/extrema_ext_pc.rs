//! OCCT Extrema_GGExtPC (TKGeomBase/Extrema/Extrema_GGExtPC.hxx L47-654) —
//! all the extremum distances between a point and a curve, dispatching on the
//! curve type (elementary curves through `Extrema_ExtPElC`, Bezier and BSpline
//! through interval decomposition over `Extrema_EPCOfExtPC`, general curves
//! through the C2/deflection interval scan).
//!
//! OCCT template instantiation (Extrema_ExtPC.hxx L31-38):
//! `Extrema_ExtPC = Extrema_GGExtPC<Adaptor3d_Curve, Extrema_CurveTool,
//! Extrema_ExtPElC, gp_Pnt, gp_Vec, Extrema_POnCurv,
//! NCollection_Sequence<Extrema_POnCurv>, Extrema_EPCOfExtPC>`, so this single
//! translation unit is both the `Extrema_GGExtPC` body and the
//! `Extrema_ExtPC` alias.
//!
//! Template-parameter mapping: `TheCurve` + `TheCurveTool` map to the
//! [`ExtPCurveTool`] trait below (the rcad `ExtremaCurveTool` facade over the
//! `Adaptor3d_Curve` evaluation traits, extended with the
//! Extrema_CurveTool.hxx L128-141 Bezier/BSpline statics that `Perform`
//! consumes); `TheExtPElC` maps to [`ExtremaExtPElC`]; `ThePoint`/`TheVector`
//! map to `glam::DVec3`; `ThePOnC`/`TheSequenceOfPOnC` map to [`POnCurve`] /
//! `Vec<POnCurve>`; `TheEPC` maps to
//! [`GenExtPC`](crate::base::extrema_gen_ext_pc::GenExtPC).

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema_curve_tool::{CurveToolHandle, ExtremaCurveTool};
use crate::base::extrema_ext_elc::{elclib_in_period, PRECISION_INFINITE, REAL_LAST};
use crate::base::extrema_ext_p_elc::ExtremaExtPElC;
use crate::base::extrema_gen_ext_pc::GenExtPC;
use crate::base::proj_lib::CurveType;
use crate::core::precision::{
    is_infinite_value, is_positive_infinite_value, CONFUSION, PCONFUSION, SQUARE_CONFUSION,
};
use crate::math::GeomAbsShape;

/// OCCT `Extrema_CurveTool` (hxx L38-142) is a SINGLE tool struct: every
/// static the `Extrema_GGExtPC` / `Extrema_GGenExtCC` bodies call lives on it,
/// including the `Bezier` / `BSpline` downcasts of hxx L136-141.  rcad carries
/// them all on [`ExtremaCurveTool`]; this trait marks the `Extrema_ExtPC`
/// instantiation sites — it is the trait-object type `ExtremaExtPC` stores.
pub trait ExtPCurveTool: ExtremaCurveTool {}

impl<T: ExtremaCurveTool + ?Sized> ExtPCurveTool for T {}

/// OCCT `Geom_BSplineCurve` knot/degree view (Extrema_GGExtPC.hxx L191-194,
/// L231): `FirstUKnotIndex()` / `LastUKnotIndex()` / `Knots()` / `Degree()`.
///
/// Architecture glue: OCCT Geom_BSplineCurve stores the reduced knot array
/// plus multiplicities; the rcad `BSplineCurve3` stores the flat knot vector
/// with multiplicities expanded (geom/mod.rs L159-173).  The reduced view is
/// reconstructed on demand; `FirstUKnotIndex()` is 1 and `LastUKnotIndex()`
/// is `knots.Length()` exactly as in Geom_BSplineCurve.  The OCCT array is
/// 1-based: knot i reads `knots[i - 1]`.
pub struct BSplineView {
    /// OCCT Geom_BSplineCurve::Knots() — the reduced (no multiplicity) knots.
    pub knots: Vec<f64>,
    /// OCCT Geom_BSplineCurve::FirstUKnotIndex() (1-based).
    pub first_u_knot_index: usize,
    /// OCCT Geom_BSplineCurve::LastUKnotIndex() (1-based).
    pub last_u_knot_index: usize,
    /// OCCT Geom_BSplineCurve::Degree().
    pub degree: usize,
}

/// OCCT Extrema_CurveTool::DeflCurvIntervals (Extrema_CurveTool.cxx L41-91)
/// — the parameters bounding the intervals of subdivision of the curve
/// according to the curvature deflection.
///
/// GAP carrier: the GCPnts_TangentialDeflection engine the OCCT body builds
/// at cxx L83 is translated 1:1 in
/// `rcad-algo/src/geomalgo/gcpnts_tangential_deflection.rs`, which
/// rcad-kernel cannot reach (crate dependency direction; the geomalgo home
/// note records the pending move into the kernel).  The prologue (cxx L44-81)
/// is translated exactly, including both early-return paths; at cxx L83 the
/// carrier raises Standard_NotImplemented — the OCCT flow is preserved up to
/// the untranslated dependency.
fn extrema_curve_tool_defl_curv_intervals(the_c: &dyn ExtPCurveTool) -> Vec<f64> {
    // cxx L44-48.
    let epsd = 1.0e-3;
    let maxdefl = 1.0e3;
    let mindefl = 1.0e-3;
    let mut nbpnts = 23i32;
    let mut l = 0.0f64;
    let tf = the_c.first_parameter();
    let tl = the_c.last_parameter();
    let a_p = the_c.value(tf);
    // cxx L52-57.
    for i in 2..=nbpnts {
        let t = (tf * ((nbpnts - i) as f64) + ((i - 1) as f64) * tl) / ((nbpnts - 1) as f64);
        let a_p1 = the_c.value(t);
        l += a_p.distance(a_p1);
    }

    // cxx L59-68.
    let d_ldt = l / (tl - tf);
    if l <= CONFUSION || d_ldt < epsd || (tl - tf) > 10000.0 {
        nbpnts = 2;
        let mut intervals = vec![0.0f64; nbpnts as usize];
        intervals[0] = tf;
        intervals[nbpnts as usize - 1] = tl;
        return intervals;
    }

    // cxx L70-78.
    let a_defl = (0.01 * l / (2.0 * std::f64::consts::PI)).max(mindefl);
    if a_defl > maxdefl {
        nbpnts = 2;
        let mut intervals = vec![0.0f64; nbpnts as usize];
        intervals[0] = tf;
        intervals[nbpnts as usize - 1] = tl;
        return intervals;
    }

    // cxx L80-83: aMinLen / aTol feed the GCPnts_TangentialDeflection
    // construction — the untranslated kernel dependency (see the GAP note).
    let a_min_len = (0.00001 * l).max(CONFUSION);
    let a_tol = (0.00001 * (tl - tf)).max(PCONFUSION);
    let _ = (a_min_len, a_tol);
    panic!(
        "Standard_NotImplemented: GCPnts_TangentialDeflection is not available in \
         rcad-kernel (Extrema_CurveTool::DeflCurvIntervals L83); the 1:1 engine \
         lives in rcad-algo/geomalgo/gcpnts_tangential_deflection.rs pending the \
         kernel move"
    );
}

/// OCCT Extrema_GGExtPC (hxx L47-654); the `Extrema_ExtPC` alias.
pub struct ExtremaExtPC<'a> {
    /// hxx L635: TheCurve* myC.
    my_c: Option<&'a dyn ExtPCurveTool>,
    /// hxx L636: ThePoint Pf.
    p_f: DVec3,
    /// hxx L637: ThePoint Pl.
    p_l: DVec3,
    /// hxx L638: TheExtPElC myExtPElC.
    my_ext_p_elc: ExtremaExtPElC,
    /// hxx L639: TheSequenceOfPOnC mypoint.
    my_point: Vec<POnCurve>,
    /// hxx L640: bool mydone.
    my_done: bool,
    /// hxx L641: double mydist1.
    my_dist1: f64,
    /// hxx L642: double mydist2.
    my_dist2: f64,
    /// hxx L643: TheEPC myExtPC.
    my_ext_pc: GenExtPC<'a>,
    /// hxx L644: double mytolu.
    my_tol_u: f64,
    /// hxx L645: double mytolf.
    my_tol_f: f64,
    /// hxx L646: int mysample.
    my_sample: i32,
    /// hxx L647: double myintuinf.
    my_int_u_inf: f64,
    /// hxx L648: double myintusup.
    my_int_u_sup: f64,
    /// hxx L649: double myuinf.
    my_u_inf: f64,
    /// hxx L650: double myusup.
    my_u_sup: f64,
    /// hxx L651: GeomAbs_CurveType type.
    type_: CurveType,
    /// hxx L652: NCollection_Sequence<bool> myismin.
    my_is_min: Vec<bool>,
    /// hxx L653: NCollection_Sequence<double> mySqDist.
    my_sq_dist: Vec<f64>,
}

impl<'a> ExtremaExtPC<'a> {
    /// OCCT Extrema_GGExtPC() default constructor (hxx L61-75).
    pub fn new() -> Self {
        ExtremaExtPC {
            my_c: None,
            p_f: DVec3::ZERO,
            p_l: DVec3::ZERO,
            my_ext_p_elc: ExtremaExtPElC::new(),
            my_point: Vec::new(),
            my_done: false,
            my_dist1: REAL_LAST,
            my_dist2: REAL_LAST,
            my_ext_pc: GenExtPC::new(),
            my_tol_u: 0.0,
            my_tol_f: 0.0,
            my_sample: 17,
            my_int_u_inf: PRECISION_INFINITE,
            my_int_u_sup: PRECISION_INFINITE,
            my_u_inf: PRECISION_INFINITE,
            my_u_sup: PRECISION_INFINITE,
            type_: CurveType::Other,
            my_is_min: Vec::new(),
            my_sq_dist: Vec::new(),
        }
    }

    /// OCCT Extrema_GGExtPC(theP, theC, theUinf, theUsup, theTolF = 1.0e-10)
    /// (hxx L84-92).
    pub fn new_point_curve_ranged(
        the_p: DVec3,
        the_c: &'a dyn ExtPCurveTool,
        the_u_inf: f64,
        the_u_sup: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut this = ExtremaExtPC::new();
        this.initialize(the_c, the_u_inf, the_u_sup, the_tol_f);
        this.perform(the_p);
        this
    }

    /// OCCT Extrema_GGExtPC(theP, theC, theTolF = 1.0e-10) (hxx L99-106) —
    /// the full curve parameter range.
    pub fn new_point_curve(
        the_p: DVec3,
        the_c: &'a dyn ExtPCurveTool,
        the_tol_f: f64,
    ) -> Self {
        let mut this = ExtremaExtPC::new();
        this.initialize(
            the_c,
            the_c.first_parameter(),
            the_c.last_parameter(),
            the_tol_f,
        );
        this.perform(the_p);
        this
    }

    /// OCCT Initialize(theC, theUinf, theUsup, theTolF = 1.0e-10)
    /// (hxx L113-128).
    pub fn initialize(
        &mut self,
        the_c: &'a dyn ExtPCurveTool,
        the_u_inf: f64,
        the_u_sup: f64,
        the_tol_f: f64,
    ) {
        self.my_c = Some(the_c);
        self.my_int_u_inf = the_u_inf;
        self.my_u_inf = the_u_inf;
        self.my_int_u_sup = the_u_sup;
        self.my_u_sup = the_u_sup;
        self.my_tol_f = the_tol_f;
        self.my_tol_u = the_c.resolution(CONFUSION);
        self.type_ = the_c.get_type();
        self.my_done = false;
        self.my_dist1 = REAL_LAST;
        self.my_dist2 = REAL_LAST;
        self.my_sample = 17;
    }

    /// OCCT Perform(theP) (hxx L132-528).
    pub fn perform(&mut self, the_p: DVec3) {
        self.my_sq_dist.clear();
        self.my_point.clear();
        self.my_is_min.clear();
        self.my_sample = 17;
        // hxx L140: constexpr double t3d = Precision::Confusion().
        let t3d = CONFUSION;

        // hxx L158: TheCurve& aCurve = *myC; (the Copy of the &'a reference
        // keeps it usable across the &mut self calls below).
        let a_curve = self.my_c;

        // hxx L142-148.
        if is_infinite_value(self.my_u_inf) {
            self.my_dist1 = REAL_LAST;
        } else {
            let c = a_curve.expect("Extrema_ExtPC: no curve");
            self.p_f = c.value(self.my_u_inf);
            self.my_dist1 = the_p.distance_squared(self.p_f);
        }

        // hxx L150-156.
        if is_infinite_value(self.my_u_sup) {
            self.my_dist2 = REAL_LAST;
        } else {
            let c = a_curve.expect("Extrema_ExtPC: no curve");
            self.p_l = c.value(self.my_u_sup);
            self.my_dist2 = the_p.distance_squared(self.p_l);
        }

        // hxx L160-471: the switch on the curve type.
        match self.type_ {
            CurveType::Circle => {
                // hxx L162-165.
                let circ = a_curve.expect("Extrema_ExtPC: no curve").circle();
                self.my_ext_p_elc
                    .perform_point_circle(the_p, &circ, t3d, self.my_u_inf, self.my_u_sup);
            }
            CurveType::Ellipse => {
                // hxx L166-169.
                let ell = a_curve.expect("Extrema_ExtPC: no curve").ellipse();
                self.my_ext_p_elc
                    .perform_point_ellipse(the_p, &ell, t3d, self.my_u_inf, self.my_u_sup);
            }
            CurveType::Parabola => {
                // hxx L170-173.
                let par = a_curve.expect("Extrema_ExtPC: no curve").parabola();
                self.my_ext_p_elc
                    .perform_point_parabola(the_p, &par, self.my_u_inf, self.my_u_sup);
            }
            CurveType::Hyperbola => {
                // hxx L174-177.
                let hyp = a_curve.expect("Extrema_ExtPC: no curve").hyperbola();
                self.my_ext_p_elc
                    .perform_point_hyperbola(the_p, &hyp, t3d, self.my_u_inf, self.my_u_sup);
            }
            CurveType::Line => {
                // hxx L178-181.
                let lin = a_curve.expect("Extrema_ExtPC: no curve").line();
                self.my_ext_p_elc
                    .perform_point_line(the_p, &lin, t3d, self.my_u_inf, self.my_u_sup);
            }
            CurveType::Bezier => {
                // hxx L182-189.
                let c = a_curve.expect("Extrema_ExtPC: no curve");
                self.my_int_u_inf = self.my_u_inf;
                self.my_int_u_sup = self.my_u_sup;
                self.my_sample = (c.bezier_nb_poles() * 2) as i32;
                self.my_ext_pc.initialize_curve(c);
                self.interval_perform(the_p);
                return;
            }
            CurveType::BSpline => {
                self.perform_bspline(the_p, a_curve);
            }
            _ => {
                self.perform_default(the_p, a_curve);
            }
        }

        // Postprocessing (hxx L473-527).
        // OCCT tests GeomAbs_BSplineCurve || GeomAbs_OffsetCurve ||
        // GeomAbs_OtherCurve; the rcad CurveType enum (proj_lib/mod.rs
        // L44-53) carries no offset value — an OCCT offset curve reports
        // Other at the rcad adaptor layer.
        if self.type_ == CurveType::BSpline || self.type_ == CurveType::Other {
            if self.my_dist1 < SQUARE_CONFUSION || self.my_dist2 < SQUARE_CONFUSION {
                // hxx L478-488.
                let mut is_first_added = false;
                let mut is_last_added = false;
                let a_nb_points = self.my_point.len();
                for i in 1..=a_nb_points {
                    let u = self.my_point[i - 1].param;
                    if (u - self.my_u_inf).abs() < self.my_tol_u {
                        is_first_added = true;
                    } else if (self.my_u_sup - u).abs() < self.my_tol_u {
                        is_last_added = true;
                    }
                }
                // hxx L489-494.
                if !is_first_added && self.my_dist1 < SQUARE_CONFUSION {
                    self.my_sq_dist.insert(0, self.my_dist1);
                    self.my_is_min.insert(0, true);
                    self.my_point
                        .insert(0, POnCurve { param: self.my_u_inf, point: self.p_f });
                }
                // hxx L495-500.
                if !is_last_added && self.my_dist2 < SQUARE_CONFUSION {
                    self.my_sq_dist.push(self.my_dist2);
                    self.my_is_min.push(true);
                    self.my_point
                        .push(POnCurve { param: self.my_u_sup, point: self.p_l });
                }
                // hxx L501.
                self.my_done = true;
            }
        } else {
            // hxx L504-527.
            self.my_done = self.my_ext_p_elc.is_done();
            if self.my_done {
                let nb_ext = self.my_ext_p_elc.nb_ext();
                for i in 1..=nb_ext {
                    // hxx L512: ThePOnC PC = myExtPElC.Point(i); (copy).
                    let mut pc = self.my_ext_p_elc.point(i).clone();
                    let mut u = pc.param;
                    let c = a_curve.expect("Extrema_ExtPC: no curve");
                    // hxx L514-517.
                    if c.is_periodic() {
                        u = elclib_in_period(u, self.my_u_inf, self.my_u_inf + c.period());
                    }
                    // hxx L518-524.
                    if (u >= self.my_u_inf - self.my_tol_u)
                        && (u <= self.my_u_sup + self.my_tol_u)
                    {
                        // hxx L520: PC.SetValues(U, myExtPElC.Point(i).Value()).
                        pc.param = u;
                        pc.point = self.my_ext_p_elc.point(i).point;
                        self.my_sq_dist.push(self.my_ext_p_elc.square_distance(i));
                        self.my_is_min.push(self.my_ext_p_elc.is_min(i));
                        self.my_point.push(pc);
                    }
                }
            }
        }
    }

    /// OCCT IsDone() (hxx L531).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT SquareDistance(theN) (hxx L536-541) — 1-based.
    pub fn square_distance(&self, the_n: usize) -> f64 {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[the_n - 1]
    }

    /// OCCT NbExt() (hxx L545-550).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_sq_dist.len()
    }

    /// OCCT IsMin(theN) (hxx L555-560) — 1-based.
    pub fn is_min(&self, the_n: usize) -> bool {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_is_min[the_n - 1]
    }

    /// OCCT Point(theN) (hxx L565-570) — 1-based.
    pub fn point(&self, the_n: usize) -> &POnCurve {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        &self.my_point[the_n - 1]
    }

    /// OCCT TrimmedSquareDistances(theDist1, theDist2, theP1, theP2)
    /// (hxx L577-586) — the out parameters map to the tuple, in the OCCT
    /// declaration order.
    pub fn trimmed_square_distances(&self) -> (f64, f64, DVec3, DVec3) {
        (self.my_dist1, self.my_dist2, self.p_f, self.p_l)
    }

    /// The GeomAbs_BSplineCurve switch arm (hxx L190-389).
    fn perform_bspline(&mut self, the_p: DVec3, a_curve: Option<&'a dyn ExtPCurveTool>) {
        let c = a_curve.expect("Extrema_ExtPC: no curve");
        // hxx L191-194: TheCurveTool::BSpline(aCurve).
        let a_bspline = c.bspline();
        let a_first_idx = a_bspline.first_u_knot_index as i32;
        let a_last_idx = a_bspline.last_u_knot_index as i32;
        let a_knots = &a_bspline.knots; // 1-based: aKnots(i) -> a_knots[i-1]

        // hxx L196-204.
        let mut a_period_jump = 0.0f64;
        let a_tol_coeff = (self.my_u_sup - self.my_u_inf) * PCONFUSION;
        if c.is_periodic() {
            let mut a_period_shift =
                ((self.my_u_inf - a_knots[(a_first_idx - 1) as usize]) / c.period()) as i32;
            if self.my_u_inf < a_knots[(a_first_idx - 1) as usize] - a_tol_coeff {
                a_period_shift -= 1;
            }
            a_period_jump = c.period() * a_period_shift as f64;
        }

        // hxx L206-215.
        let mut an_idx: i32;
        let mut a_first_used_knot = a_first_idx;
        let mut a_last_used_knot = a_last_idx;
        an_idx = a_first_idx;
        while an_idx <= a_last_idx {
            let a_knot = a_knots[(an_idx - 1) as usize] + a_period_jump;
            if self.my_u_inf >= a_knot - a_tol_coeff {
                a_first_used_knot = an_idx;
            } else {
                break;
            }
            an_idx += 1;
        }
        // hxx L216-223.
        an_idx = a_last_idx;
        while an_idx >= a_first_idx {
            let a_knot = a_knots[(an_idx - 1) as usize] + a_period_jump;
            if self.my_u_sup <= a_knot + a_tol_coeff {
                a_last_used_knot = an_idx;
            } else {
                break;
            }
            an_idx -= 1;
        }

        // hxx L225-229.
        if a_first_used_knot == a_last_used_knot {
            a_first_used_knot = a_first_idx;
            a_last_used_knot = a_first_idx + 1;
        }

        // hxx L231.
        self.my_sample = a_bspline.degree as i32 + 1;

        if self.my_sample == 2 {
            // hxx L233-298.
            let mut a_pmin = DVec3::ZERO;
            let mut tmin = 0.0;
            let mut distmin = REAL_LAST;
            let mut a_min1;
            let mut a_min2 = 0.0;
            self.my_ext_pc.initialize_curve(c);
            an_idx = a_first_used_knot;
            while an_idx < a_last_used_knot {
                // hxx L241.
                let mut a_f = a_knots[(an_idx - 1) as usize] + a_period_jump;
                let mut a_l = a_knots[an_idx as usize] + a_period_jump;

                // hxx L243-246.
                if an_idx == a_first_used_knot {
                    a_f = self.my_u_inf;
                } else if an_idx == a_last_used_knot - 1 {
                    a_l = self.my_u_sup;
                }

                // hxx L248-254.
                let a_p1 = c.d0(a_f);
                let a_p2 = c.d0(a_l);
                let a_base1 = a_p1 - the_p; // TheVector aBase1(theP, aP1)
                let a_base2 = a_p2 - the_p; // TheVector aBase2(theP, aP2)
                let a_v = a_p1 - a_p2; // TheVector aV(aP2, aP1)
                let a_val1 = a_v.dot(a_base1);
                let a_val2 = a_v.dot(a_base2);
                // hxx L255-266.
                if an_idx == a_first_used_knot {
                    a_min1 = the_p.distance_squared(a_p1);
                } else {
                    a_min1 = a_min2;
                    if distmin > a_min1 {
                        distmin = a_min1;
                        tmin = a_f;
                        a_pmin = a_p1;
                    }
                }
                // hxx L267-269.
                a_min2 = the_p.distance_squared(a_p2);
                let a_min_sq_dist = a_min1.min(a_min2);
                let a_min_der = a_val1.abs().min(a_val2.abs());
                // hxx L270-279.
                if !(is_infinite_value(a_val1) || is_infinite_value(a_val2)) {
                    if a_val1 * a_val2 <= 0.0
                        || a_min_sq_dist < 100.0 * SQUARE_CONFUSION
                        || 2.0 * a_min_der < CONFUSION
                    {
                        self.my_int_u_inf = a_f;
                        self.my_int_u_sup = a_l;
                        self.interval_perform(the_p);
                    }
                }
                an_idx += 1;
            }
            // hxx L281-297.
            if !is_infinite_value(distmin) {
                let mut is_to_add = true;
                let nb_ext = self.my_point.len();
                let mut i = 1usize;
                while i <= nb_ext && is_to_add {
                    let t = self.my_point[i - 1].param;
                    is_to_add = (distmin < self.my_sq_dist[i - 1])
                        && ((t - tmin).abs() > self.my_tol_u);
                    i += 1;
                }
                if is_to_add {
                    let pc = POnCurve { param: tmin, point: a_pmin };
                    self.my_sq_dist.push(distmin);
                    self.my_is_min.push(true);
                    self.my_point.push(pc);
                }
            }
        } else {
            // hxx L299-386.
            let mut a_val_idx = 1usize;
            // hxx L302-303: NCollection_Array1 (1, mysample*(L-F) + 1).
            let nb_slots =
                (self.my_sample as usize) * ((a_last_used_knot - a_first_used_knot) as usize) + 1;
            let mut a_val = vec![0.0f64; nb_slots + 1];
            let mut a_param = vec![0.0f64; nb_slots + 1];
            an_idx = a_first_used_knot;
            while an_idx < a_last_used_knot {
                // hxx L306.
                let mut a_f = a_knots[(an_idx - 1) as usize] + a_period_jump;
                let mut a_l = a_knots[an_idx as usize] + a_period_jump;

                // hxx L308-311.
                if an_idx == a_first_used_knot {
                    a_f = self.my_u_inf;
                }
                if an_idx == a_last_used_knot - 1 {
                    a_l = self.my_u_sup;
                }

                // hxx L313-320.
                let a_step = (a_l - a_f) / self.my_sample as f64;
                for a_pnt_idx in 0..self.my_sample {
                    let a_current_param = a_f + a_step * a_pnt_idx as f64;
                    a_val[a_val_idx] = c.value(a_current_param).distance_squared(the_p);
                    a_param[a_val_idx] = a_current_param;
                    a_val_idx += 1;
                }
                an_idx += 1;
            }
            // hxx L322-323.
            a_val[a_val_idx] = c.value(self.my_u_sup).distance_squared(the_p);
            a_param[a_val_idx] = self.my_u_sup;

            // hxx L325.
            self.my_ext_pc.initialize_curve(c);

            // hxx L327-342: for (anIdx = aVal.Lower() + 1; anIdx < aVal.Upper(); anIdx++).
            let mut an_idx2 = 2usize;
            while an_idx2 < nb_slots {
                if a_val[an_idx2] <= SQUARE_CONFUSION {
                    self.my_sq_dist.push(a_val[an_idx2]);
                    self.my_is_min.push(true);
                    self.my_point.push(POnCurve {
                        param: a_param[an_idx2],
                        point: c.value(a_param[an_idx2]),
                    });
                }
                if (a_val[an_idx2] >= a_val[an_idx2 + 1] && a_val[an_idx2] >= a_val[an_idx2 - 1])
                    || (a_val[an_idx2] <= a_val[an_idx2 + 1]
                        && a_val[an_idx2] <= a_val[an_idx2 - 1])
                {
                    self.my_int_u_inf = a_param[an_idx2 - 1];
                    self.my_int_u_sup = a_param[an_idx2 + 1];
                    self.interval_perform(the_p);
                }
                an_idx2 += 1;
            }

            // hxx L344-363.
            if self.my_dist1 > SQUARE_CONFUSION && !is_positive_infinite_value(self.my_dist1) {
                let (a_p1, a_v1) = c.d1(a_param[1]);
                let (a_p2, a_v2) = c.d1(a_param[2]);
                let a_base1 = a_p1 - the_p;
                let a_base2 = a_p2 - the_p;
                let a_val1 = a_v1.dot(a_base1);
                let a_val2 = a_v2.dot(a_base2);
                if !(is_infinite_value(a_val1) || is_infinite_value(a_val2)) {
                    if a_val1 * a_val2 <= 0.0
                        || a_base1.dot(a_base2) <= 0.0
                        || 2.0 * a_val1.abs() < CONFUSION
                    {
                        self.my_int_u_inf = a_param[1];
                        self.my_int_u_sup = a_param[2];
                        self.interval_perform(the_p);
                    }
                }
            }

            // hxx L365-385.
            if self.my_dist2 > SQUARE_CONFUSION && !is_positive_infinite_value(self.my_dist2) {
                let (a_p1, a_v1) = c.d1(a_param[nb_slots - 1]);
                let (a_p2, a_v2) = c.d1(a_param[nb_slots]);
                let a_base1 = a_p1 - the_p;
                let a_base2 = a_p2 - the_p;
                let a_val1 = a_v1.dot(a_base1);
                let a_val2 = a_v2.dot(a_base2);

                if !(is_infinite_value(a_val1) || is_infinite_value(a_val2)) {
                    if a_val1 * a_val2 <= 0.0
                        || a_base1.dot(a_base2) <= 0.0
                        || 2.0 * a_val2.abs() < CONFUSION
                    {
                        self.my_int_u_inf = a_param[nb_slots - 1];
                        self.my_int_u_sup = a_param[nb_slots];
                        self.interval_perform(the_p);
                    }
                }
            }
        }
        // hxx L387.
        self.my_done = true;
    }

    /// The default switch arm (hxx L390-470) — the C2 / curvature-deflection
    /// interval scan over `Extrema_EPCOfExtPC`.
    fn perform_default(&mut self, the_p: DVec3, a_curve: Option<&'a dyn ExtPCurveTool>) {
        let c = a_curve.expect("Extrema_ExtPC: no curve");
        // hxx L391-394.
        let a_max_sample = 17i32;
        let mut int_ext_is_done = false;

        // hxx L395-405.
        let n = c.nb_intervals(GeomAbsShape::C2);
        let the_h_inter: Vec<f64> = if n > 1 {
            // hxx L398-399: Intervals(aCurve, theHInter->ChangeArray1(), GeomAbs_C2).
            c.intervals(GeomAbsShape::C2)
        } else {
            // hxx L403: TheCurveTool::DeflCurvIntervals(aCurve).
            extrema_curve_tool_defl_curv_intervals(c)
        };
        let n = the_h_inter.len() as i32 - 1; // hxx L404 (holds for both paths)

        // hxx L406.
        self.my_sample = (self.my_sample / n).max(a_max_sample);
        // hxx L407-413.
        let mut maxint = 0.0f64;
        for i in 1..=n {
            let dt = the_h_inter[i as usize] - the_h_inter[(i - 1) as usize];
            if maxint < dt {
                maxint = dt;
            }
        }
        // hxx L414-418.
        let is_periodic = c.is_periodic();
        let mut s1 = 0.0f64;
        let mut s2 = 0.0f64;
        // hxx L419.
        self.my_ext_pc.initialize_curve(c);
        // hxx L420-466.
        for i in 1..=n {
            self.my_int_u_inf = the_h_inter[(i - 1) as usize];
            self.my_int_u_sup = the_h_inter[i as usize];
            self.my_sample = ((a_max_sample as f64 * (self.my_int_u_sup - self.my_int_u_inf)
                / maxint) as i32)
                .max(3);

            // hxx L426-427.
            let mut an_inf_to_check = self.my_int_u_inf;
            let mut a_sup_to_check = self.my_int_u_sup;

            // hxx L429-434.
            if is_periodic {
                let a_period = c.period();
                an_inf_to_check =
                    elclib_in_period(self.my_int_u_inf, self.my_u_inf, self.my_u_inf + a_period);
                a_sup_to_check = self.my_int_u_sup + (an_inf_to_check - self.my_int_u_inf);
            }
            // hxx L435-438.
            let int_is_not_valid =
                (self.my_u_inf > a_sup_to_check) || (self.my_u_sup < an_inf_to_check);

            if int_is_not_valid {
                continue;
            }

            // hxx L440-445.
            if self.my_u_inf >= an_inf_to_check {
                an_inf_to_check = self.my_u_inf;
            }
            if self.my_u_sup <= a_sup_to_check {
                a_sup_to_check = self.my_u_sup;
            }
            if (a_sup_to_check - an_inf_to_check) <= self.my_tol_u {
                continue;
            }

            // hxx L447-457.
            if i != 1 {
                let (pp, v1) = c.d1(self.my_int_u_inf);
                s1 = (pp - the_p).dot(v1);
                if s1 * s2 < 0.0 {
                    self.my_sq_dist.push(pp.distance_squared(the_p));
                    self.my_is_min.push(s1 < 0.0);
                    self.my_point.push(POnCurve {
                        param: self.my_int_u_inf,
                        point: pp,
                    });
                }
            }
            // hxx L458-462.
            if i != n {
                let (pp, v1) = c.d1(self.my_int_u_sup);
                s2 = (pp - the_p).dot(v1);
            }

            // hxx L464-465.
            self.interval_perform(the_p);
            int_ext_is_done = int_ext_is_done || self.my_done;
        }

        // hxx L468.
        self.my_done = int_ext_is_done;
    }

    /// OCCT IntervalPerform(theP) (hxx L590-614).
    fn interval_perform(&mut self, the_p: DVec3) {
        self.my_ext_pc.initialize_interval(
            self.my_sample,
            self.my_int_u_inf,
            self.my_int_u_sup,
            self.my_tol_u,
            self.my_tol_f,
        );
        self.my_ext_pc.perform(the_p);
        self.my_done = self.my_ext_pc.is_done();
        if self.my_done {
            let nb_ext = self.my_ext_pc.nb_ext();
            for i in 1..=nb_ext {
                // hxx L602-603.
                let pc = self.my_ext_pc.point(i).clone();
                let mut u = pc.param;
                let c = self.my_c.expect("Extrema_ExtPC: no curve");
                // hxx L604-607.
                if c.is_periodic() {
                    u = elclib_in_period(u, self.my_u_inf, self.my_u_inf + c.period());
                }
                // hxx L608-611.
                if (u >= self.my_u_inf - self.my_tol_u) && (u <= self.my_u_sup + self.my_tol_u)
                {
                    self.add_sol(u, pc.point, self.my_ext_pc.square_distance(i), self.my_ext_pc.is_min(i));
                }
            }
        }
    }

    /// OCCT AddSol(theU, theP, theSqDist, isMin) (hxx L617-632).
    fn add_sol(&mut self, the_u: f64, the_p: DVec3, the_sq_dist: f64, is_min: bool) {
        let nb_ext = self.my_point.len();
        for i in 1..=nb_ext {
            let t = self.my_point[i - 1].param;
            if (t - the_u).abs() <= self.my_tol_u {
                return;
            }
        }
        let pc = POnCurve { param: the_u, point: the_p };
        self.my_sq_dist.push(the_sq_dist);
        self.my_is_min.push(is_min);
        self.my_point.push(pc);
    }
}
