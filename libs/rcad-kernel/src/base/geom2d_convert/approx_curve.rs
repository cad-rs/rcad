//! OCCT Geom2dConvert_ApproxCurve (TKGeomBase/Geom2dConvert).
//!
//! Sources (1:1 statement-mapped):
//! - `Geom2dConvert_ApproxCurve.hxx` L29-93 — the class frame (ctor
//!   defaults, member fields);
//! - `Geom2dConvert_ApproxCurve.cxx` L31-108 — the AdvApprox_EvaluatorFunction
//!   subclass driving Value/D1/D2 through the curve adaptor;
//! - `Geom2dConvert_ApproxCurve.cxx` L110-209 — the two constructors, the
//!   Approximate driver over AdvApprox_ApproxAFunction with the
//!   AdvApprox_PrefAndRec cutting tool, and the accessors.
//!
//! The engines live kernel-side (`math::adv_approx::ApproxAFunction`,
//! `math::adv_approx::PrefAndRec`) — no duplication across crates.
//!
//! The adaptor is the kernel `Geom2dCurveAdaptor`
//! (`base::proj_lib::adaptor`), the Geom2dAdaptor_Curve encoding: the
//! unrestricted constructor form wraps a `Curve2d` and evaluates through
//! `Curve2dEval`, so the `Curve2d::Offset` payload
//! (`OffsetCurve2d { basis, offset_distance }`) carries everything the
//! adaptor's Value/D1/D2 need.

use std::sync::Arc;

use crate::base::proj_lib::adaptor::{Adaptor2dCurve2d, Geom2dCurveAdaptor};
use crate::geom::{BSplineCurve2, Curve2d};
use crate::math::adv_approx::EvaluatorFunction;
use crate::math::adv_approx::{ApproxAFunction, PrefAndRec};
use crate::math::GeomAbsShape;

// ---------------------------------------------------------------------------
// OCCT Geom2dConvert_ApproxCurve.cxx L31-108 — the evaluator.
// ---------------------------------------------------------------------------

/// OCCT Geom2dConvert_ApproxCurve_Eval (cxx L32-54) — holds the curve
/// adaptor and the saved evaluation window.
pub struct ApproxCurveEval {
    /// OCCT: occ::handle<Adaptor2d_Curve2d> fonct.
    fonct: Arc<dyn Adaptor2dCurve2d>,
    /// OCCT: double StartEndSav[2].
    start_end_sav: [f64; 2],
}

impl ApproxCurveEval {
    /// OCCT ctor (cxx L35-42).
    fn new(the_func: Arc<dyn Adaptor2dCurve2d>, first: f64, last: f64) -> Self {
        ApproxCurveEval {
            fonct: the_func,
            start_end_sav: [first, last],
        }
    }
}

impl EvaluatorFunction for ApproxCurveEval {
    /// OCCT Geom2dConvert_ApproxCurve_Eval::Evaluate (cxx L56-108).  The
    /// dimension travels in the `result` slice length (the OCCT
    /// `*Dimension`); the returned value is `*ErrorCode`.
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // cxx L63: *ErrorCode = 0.
        let mut error_code = 0;
        // cxx L64: double par = *Param.
        let par = parameter;

        // cxx L67-70: Dimension is incorrect.
        if result.len() != 2 {
            error_code = 1;
        }
        // cxx L72-75: Parameter is incorrect.
        if par < start_end[0] || par > start_end[1] {
            error_code = 2;
        }
        // cxx L76-81.
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct = self
                .fonct
                .trim(start_end[0], start_end[1], crate::core::precision::PCONFUSION);
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

        // cxx L83-84: gp_Pnt2d pnt; gp_Vec2d v1, v2.
        // cxx L86-107: the derivative-request switch.
        match derivative_request {
            0 => {
                let pnt = self.fonct.value(par);
                result[0] = pnt.x;
                result[1] = pnt.y;
            }
            1 => {
                let (_pnt, v1) = self.fonct.d1(par);
                result[0] = v1.x;
                result[1] = v1.y;
            }
            2 => {
                let (_pnt, _v1, v2) = self.fonct.d2(par);
                result[0] = v2.x;
                result[1] = v2.y;
            }
            _ => {
                result[0] = 0.0;
                result[1] = 0.0;
                error_code = 3;
            }
        }
        error_code
    }
}

// ---------------------------------------------------------------------------
// OCCT Geom2dConvert_ApproxCurve.hxx L29-93 — the class frame.
// ---------------------------------------------------------------------------

/// OCCT Geom2dConvert_ApproxCurve — converts a 2D curve to a BSpline by
/// approximation within a given tolerance.
pub struct Geom2dConvertApproxCurve {
    /// OCCT hxx L89: bool myIsDone.
    my_is_done: bool,
    /// OCCT hxx L90: bool myHasResult.
    my_has_result: bool,
    /// OCCT hxx L91: occ::handle<Geom2d_BSplineCurve> myBSplCurve.
    my_bspl_curve: Option<BSplineCurve2>,
    /// OCCT hxx L92: double myMaxError.
    my_max_error: f64,
}

impl Geom2dConvertApproxCurve {
    /// OCCT ctor from a Geom2d_Curve (cxx L110-118): the curve is wrapped
    /// into a Geom2dAdaptor_Curve and handed to Approximate.
    pub fn new(
        curve: &Curve2d,
        tol2d: f64,
        order: GeomAbsShape,
        max_segments: i32,
        max_degree: i32,
    ) -> Self {
        // OCCT L116: handle<Geom2dAdaptor_Curve> HCurve =
        // new Geom2dAdaptor_Curve(Curve).
        let hcurve: Arc<dyn Adaptor2dCurve2d> = Arc::new(Geom2dCurveAdaptor::new(curve.clone()));
        let mut this = Geom2dConvertApproxCurve {
            my_is_done: false,
            my_has_result: false,
            my_bspl_curve: None,
            my_max_error: 0.0,
        };
        this.approximate(hcurve, tol2d, order, max_segments, max_degree);
        this
    }

    /// OCCT ctor from an Adaptor2d_Curve2d (cxx L120-127).
    pub fn new_adaptor(
        curve: Arc<dyn Adaptor2dCurve2d>,
        tol2d: f64,
        order: GeomAbsShape,
        max_segments: i32,
        max_degree: i32,
    ) -> Self {
        let mut this = Geom2dConvertApproxCurve {
            my_is_done: false,
            my_has_result: false,
            my_bspl_curve: None,
            my_max_error: 0.0,
        };
        this.approximate(curve, tol2d, order, max_segments, max_degree);
        this
    }

    /// OCCT Geom2dConvert_ApproxCurve::Approximate (cxx L129-183) — the
    /// AdvApprox_ApproxAFunction driver with the AdvApprox_PrefAndRec
    /// cutting tool.
    pub fn approximate(
        &mut self,
        the_curve: Arc<dyn Adaptor2dCurve2d>,
        the_tol2d: f64,
        the_order: GeomAbsShape,
        the_max_segments: i32,
        the_max_degree: i32,
    ) {
        // Initialisation of input parameters of AdvApprox (cxx L135).
        // OCCT L137: int Num1DSS = 0, Num2DSS = 1, Num3DSS = 0.
        let num1dss = 0i32;
        let num2dss = 1i32;
        let num3dss = 0i32;
        // OCCT L138-140: OneDTolNul / ThreeDTolNul are null handles;
        // TwoDTol = new HArray1(1, Num2DSS) initialized with theTol2d.
        let one_d_tol: Option<&[f64]> = None;
        let two_d_tol = [the_tol2d];
        let three_d_tol: Option<&[f64]> = None;

        // OCCT L142-143.
        let first = the_curve.first_parameter();
        let last = the_curve.last_parameter();

        // OCCT L145-150: the C2 / C3 cutting point arrays.
        let nb_interv_c2 = the_curve.nb_intervals(GeomAbsShape::C2);
        let cut_pnts_c2 = the_curve.intervals(GeomAbsShape::C2);
        debug_assert_eq!(cut_pnts_c2.len(), nb_interv_c2 + 1);
        let nb_interv_c3 = the_curve.nb_intervals(GeomAbsShape::C3);
        let cut_pnts_c3 = the_curve.intervals(GeomAbsShape::C3);
        debug_assert_eq!(cut_pnts_c3.len(), nb_interv_c3 + 1);
        // OCCT L151: AdvApprox_PrefAndRec CutTool(CutPnts_C2, CutPnts_C3)
        // — the 2-argument form with the default Weight = 5.
        let cut_tool = PrefAndRec::with_default_weight(&cut_pnts_c2, &cut_pnts_c3);

        // OCCT L153: myMaxError = 0.
        self.my_max_error = 0.0;

        // OCCT L155: Geom2dConvert_ApproxCurve_Eval ev(theCurve, First, Last).
        let mut ev = ApproxCurveEval::new(the_curve.clone(), first, last);

        // OCCT L156-168: the AdvApprox_ApproxAFunction construction (the
        // 13-argument ctor running Perform with the user cutting tool).
        let a_approx = ApproxAFunction::with_cut_tool(
            num1dss,
            num2dss,
            num3dss,
            one_d_tol,
            Some(&two_d_tol),
            three_d_tol,
            first,
            last,
            the_order,
            the_max_degree,
            the_max_segments,
            &mut ev,
            &cut_tool,
        );

        // OCCT L170-171.
        self.my_is_done = a_approx.is_done();
        self.my_has_result = a_approx.has_result();

        if self.my_has_result {
            // OCCT L175-179: the result poles / knots / mults / degree.
            let poles_flat = a_approx.poles2d_flat(1);
            let knots = a_approx.knots_vec();
            let mults = a_approx.multiplicities_vec();
            let degree = a_approx.degree();
            // OCCT L180: new Geom2d_BSplineCurve(Poles, Knots, Mults,
            // Degree) — the non-periodic (knots, mults) constructor; the
            // rcad carrier stores the flat (expanded) knot vector (the
            // same bijection as `geom::bspline_ops`), unit weights (the
            // AdvApprox 2D result is non-rational).
            let control_points: Vec<glam::DVec2> = poles_flat
                .chunks_exact(2)
                .map(|ch| glam::DVec2::new(ch[0], ch[1]))
                .collect();
            let weights = vec![1.0; control_points.len()];
            let mut knots_flat = Vec::with_capacity(knots.len());
            for (k, m) in knots.iter().zip(mults.iter()) {
                for _ in 0..*m {
                    knots_flat.push(*k);
                }
            }
            self.my_bspl_curve = Some(BSplineCurve2 {
                degree: degree as usize,
                knots: knots_flat,
                control_points,
                weights,
                is_periodic: false,
            });
            // OCCT L181.
            self.my_max_error = a_approx.max_error_at(2, 1);
        }
    }

    /// OCCT Geom2dConvert_ApproxCurve::Curve() (cxx L185-188).
    pub fn curve(&self) -> Option<BSplineCurve2> {
        self.my_bspl_curve.clone()
    }

    /// OCCT Geom2dConvert_ApproxCurve::IsDone() (cxx L190-193).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT Geom2dConvert_ApproxCurve::HasResult() (cxx L195-198).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT Geom2dConvert_ApproxCurve::MaxError() (cxx L200-203).
    pub fn max_error(&self) -> f64 {
        self.my_max_error
    }

    /// OCCT Geom2dConvert_ApproxCurve::Dump(Standard_OStream& o)
    /// (cxx L205-209) — the stream is the caller's String buffer.
    pub fn dump(&self, o: &mut String) {
        o.push_str("******* Dump of ApproxCurve *******\n");
        o.push_str(format!("******* Error   {}\n", self.max_error()).as_str());
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle2d, OffsetCurve2d};
    use glam::DVec2;

    /// The offset of a circle by d is the concentric circle of radius
    /// R + d; the approximation must land on it within the AdvApprox
    /// tolerance (the offset evaluation rides the kernel Curve2d exact
    /// BSpline derivatives through the adaptor).
    fn offset_circle() -> Curve2d {
        Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(Curve2d::Circle(Circle2d {
                center: DVec2::new(1.0, 2.0),
                x_dir: DVec2::new(1.0, 0.0),
                y_dir: DVec2::new(0.0, 1.0),
                radius: 3.0,
            })),
            offset_distance: 0.5,
        })
    }

    #[test]
    fn approx_offset_circle_lands_on_target() {
        let c = offset_circle();
        let appr = Geom2dConvertApproxCurve::new(
            &c,
            1.0e-4,
            GeomAbsShape::C2,
            16,
            14,
        );
        assert!(appr.has_result(), "HasResult");
        assert!(appr.is_done(), "IsDone");
        let bs = appr.curve().expect("result curve");
        // Sample against the analytic offset circle of radius 3.5.
        for k in 0..=8 {
            let t = std::f64::consts::FRAC_PI_2 * (k as f64) / 8.0;
            let p = crate::geom::Curve2dEval::point_at(&bs, t);
            let r = (p - DVec2::new(1.0, 2.0)).length();
            assert!(
                (r - 3.5).abs() < 1.0e-3,
                "u={t}: radius {r}, want 3.5 (max error {})",
                appr.max_error()
            );
        }
        let mut dump = String::new();
        appr.dump(&mut dump);
        assert!(dump.contains("Dump of ApproxCurve"));
    }

    /// The Geom2dConvert::CurveToBSplineCurve offset branches
    /// (Geom2dConvert.cxx L353-368 trimmed, L425-440 full): the
    /// approximation replaces the offset with a non-rational BSpline.
    #[test]
    fn curve_to_bspline_offset_branches() {
        use crate::base::convert::ConvertParameterisation;
        use crate::geom::{Curve2dEval, TrimmedCurve2};

        // The full offset branch.
        let bs = super::super::curve_to_bspline_curve_2d(
            &offset_circle(),
            ConvertParameterisation::TgtThetaOver2,
        );
        for k in 0..=4 {
            let t = std::f64::consts::FRAC_PI_2 * (k as f64) / 4.0;
            let p = Curve2dEval::point_at(&bs, t);
            let r = (p - DVec2::new(1.0, 2.0)).length();
            assert!((r - 3.5).abs() < 1.0e-3, "full branch u={t}: r={r}");
        }

        // The trimmed offset branch — the trim window restricts the
        // approximation domain.
        let trc = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(offset_circle()),
            t_min: 0.2,
            t_max: 1.1,
        });
        let bs_trim = super::super::curve_to_bspline_curve_2d(
            &trc,
            ConvertParameterisation::TgtThetaOver2,
        );
        for k in 0..=4 {
            let t = 0.2 + 0.9 * (k as f64) / 4.0;
            let p = Curve2dEval::point_at(&bs_trim, t);
            let r = (p - DVec2::new(1.0, 2.0)).length();
            assert!((r - 3.5).abs() < 1.0e-3, "trimmed branch u={t}: r={r}");
        }
    }
}
