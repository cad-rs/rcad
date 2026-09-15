//! OCCT GeomConvert_ApproxCurve (TKGeomBase/GeomConvert).
//!
//! Sources (1:1 statement-mapped), OCCT tree root
//! `src/ModelingData/TKGeomBase/GeomConvert/`:
//! - `GeomConvert_ApproxCurve.hxx` L29-94 — the class frame (no default
//!   arguments in this OCCT version: every caller passes Tol3d / Order /
//!   MaxSegments / MaxDegree explicitly);
//! - `GeomConvert_ApproxCurve.cxx` L31-106 — the AdvApprox_EvaluatorFunction
//!   subclass driving Value/D1/D2 through the curve adaptor;
//! - `GeomConvert_ApproxCurve.cxx` L108-208 — the two constructors, the
//!   Approximate driver over AdvApprox_ApproxAFunction with the
//!   AdvApprox_PrefAndRec cutting tool, and the accessors / Dump.
//!
//! The engines live kernel-side (`math::adv_approx::ApproxAFunction`,
//! `math::adv_approx::PrefAndRec`) — no duplication across crates.  The
//! structural precedent is the 2D sibling
//! `base::geom2d_convert::approx_curve::Geom2dConvertApproxCurve`.
//!
//! The adaptor is the kernel `GeomCurveAdaptor` (`base::proj_lib::geom_adaptor_curve`),
//! the GeomAdaptor_Curve encoding: the plain constructor form wraps a `Curve3`
//! and evaluates through `CurveEval`, so every `Curve3` payload carries what
//! the adaptor's Value/D1/D2 need.

use std::sync::Arc;

use glam::DVec3;

use crate::base::proj_lib::adaptor::{Adaptor3dCurve, CurveHandle};
use crate::base::proj_lib::geom_adaptor_curve::GeomCurveAdaptor;
use crate::geom::{BSplineCurve3, Curve3};
use crate::math::adv_approx::{ApproxAFunction, EvaluatorFunction, PrefAndRec};
use crate::math::GeomAbsShape;

// ---------------------------------------------------------------------------
// OCCT GeomConvert_ApproxCurve.cxx L31-106 — the evaluator.
// ---------------------------------------------------------------------------

/// OCCT GeomConvert_ApproxCurve_Eval (cxx L31-53) — holds the curve adaptor
/// and the saved evaluation window.
pub struct ApproxCurveEval {
    /// OCCT: occ::handle<Adaptor3d_Curve> fonct.
    fonct: Arc<dyn Adaptor3dCurve>,
    /// OCCT: double StartEndSav[2].
    start_end_sav: [f64; 2],
}

impl ApproxCurveEval {
    /// OCCT ctor (cxx L34-41).
    fn new(the_func: Arc<dyn Adaptor3dCurve>, first: f64, last: f64) -> Self {
        ApproxCurveEval {
            fonct: the_func,
            start_end_sav: [first, last],
        }
    }
}

impl EvaluatorFunction for ApproxCurveEval {
    /// OCCT GeomConvert_ApproxCurve_Eval::Evaluate (cxx L55-106).  The
    /// dimension travels in the `result` slice length (the OCCT
    /// `*Dimension`); the returned value is `*ErrorCode`.
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // cxx L62: *ErrorCode = 0.
        let mut error_code = 0;
        // cxx L63: double par = *Param.
        let par = parameter;

        // cxx L66-69: Dimension is incorrect.
        if result.len() != 3 {
            error_code = 1;
        }

        // cxx L71-76.
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct = self
                .fonct
                .trim(start_end[0], start_end[1], crate::core::precision::PCONFUSION);
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

        // cxx L78-79: gp_Pnt pnt; gp_Vec v1, v2.
        // cxx L81-105: the derivative-request switch.
        match derivative_request {
            0 => {
                let pnt = self.fonct.value(par);
                result[0] = pnt.x;
                result[1] = pnt.y;
                result[2] = pnt.z;
            }
            1 => {
                let (_pnt, v1) = self.fonct.d1(par);
                result[0] = v1.x;
                result[1] = v1.y;
                result[2] = v1.z;
            }
            2 => {
                let (_pnt, _v1, v2) = self.fonct.d2(par);
                result[0] = v2.x;
                result[1] = v2.y;
                result[2] = v2.z;
            }
            _ => {
                result[0] = 0.0;
                result[1] = 0.0;
                result[2] = 0.0;
                error_code = 3;
            }
        }
        error_code
    }
}

// ---------------------------------------------------------------------------
// OCCT GeomConvert_ApproxCurve.hxx L29-94 — the class frame.
// ---------------------------------------------------------------------------

/// OCCT GeomConvert_ApproxCurve — converts a 3D curve to a 3D BSpline by
/// approximation within a given tolerance.
pub struct GeomConvertApproxCurve {
    /// OCCT hxx L90: bool myIsDone.
    my_is_done: bool,
    /// OCCT hxx L91: bool myHasResult.
    my_has_result: bool,
    /// OCCT hxx L92: occ::handle<Geom_BSplineCurve> myBSplCurve.
    my_bspl_curve: Option<BSplineCurve3>,
    /// OCCT hxx L93: double myMaxError.
    my_max_error: f64,
}

impl GeomConvertApproxCurve {
    /// OCCT ctor from a Geom_Curve (cxx L108-116): the curve is wrapped into
    /// a GeomAdaptor_Curve and handed to Approximate.
    pub fn new(
        curve: &Curve3,
        tol3d: f64,
        order: GeomAbsShape,
        max_segments: i32,
        max_degree: i32,
    ) -> Self {
        // OCCT L114: handle<GeomAdaptor_Curve> HCurve = new
        // GeomAdaptor_Curve(Curve).
        let hcurve: Arc<dyn Adaptor3dCurve> = Arc::new(GeomCurveAdaptor::new(curve.clone()));
        let mut this = GeomConvertApproxCurve {
            my_is_done: false,
            my_has_result: false,
            my_bspl_curve: None,
            my_max_error: 0.0,
        };
        this.approximate(hcurve, tol3d, order, max_segments, max_degree);
        this
    }

    /// OCCT ctor from an Adaptor3d_Curve (cxx L118-125).
    pub fn new_adaptor(
        curve: Arc<dyn Adaptor3dCurve>,
        tol3d: f64,
        order: GeomAbsShape,
        max_segments: i32,
        max_degree: i32,
    ) -> Self {
        let mut this = GeomConvertApproxCurve {
            my_is_done: false,
            my_has_result: false,
            my_bspl_curve: None,
            my_max_error: 0.0,
        };
        this.approximate(curve, tol3d, order, max_segments, max_degree);
        this
    }

    /// OCCT GeomConvert_ApproxCurve::Approximate (cxx L127-182) — the
    /// AdvApprox_ApproxAFunction driver with the AdvApprox_PrefAndRec
    /// cutting tool.
    pub fn approximate(
        &mut self,
        the_curve: Arc<dyn Adaptor3dCurve>,
        the_tol3d: f64,
        the_order: GeomAbsShape,
        the_max_segments: i32,
        the_max_degree: i32,
    ) {
        // Initialisation of input parameters of AdvApprox (cxx L133).
        // OCCT L135: int Num1DSS = 0, Num2DSS = 0, Num3DSS = 1.
        let num1dss = 0i32;
        let num2dss = 0i32;
        let num3dss = 1i32;
        // OCCT L136-138: OneDTolNul / TwoDTolNul are null handles;
        // ThreeDTol = new HArray1(1, Num3DSS) initialized with theTol3d.
        let one_d_tol: Option<&[f64]> = None;
        let two_d_tol: Option<&[f64]> = None;
        let three_d_tol = [the_tol3d];

        // OCCT L140-141.
        let first = the_curve.first_parameter();
        let last = the_curve.last_parameter();

        // OCCT L143-148: the C2 / C3 cutting point arrays.
        let nb_interv_c2 = the_curve.nb_intervals(GeomAbsShape::C2);
        let cut_pnts_c2 = the_curve.intervals(GeomAbsShape::C2);
        debug_assert_eq!(cut_pnts_c2.len(), nb_interv_c2 + 1);
        let nb_interv_c3 = the_curve.nb_intervals(GeomAbsShape::C3);
        let cut_pnts_c3 = the_curve.intervals(GeomAbsShape::C3);
        debug_assert_eq!(cut_pnts_c3.len(), nb_interv_c3 + 1);
        // OCCT L150: AdvApprox_PrefAndRec CutTool(CutPnts_C2, CutPnts_C3)
        // — the 2-argument form with the default Weight = 5.
        let cut_tool = PrefAndRec::with_default_weight(&cut_pnts_c2, &cut_pnts_c3);

        // OCCT L152: myMaxError = 0.
        self.my_max_error = 0.0;

        // OCCT L154: GeomConvert_ApproxCurve_Eval ev(theCurve, First, Last).
        let mut ev = ApproxCurveEval::new(the_curve.clone(), first, last);

        // OCCT L155-167: the AdvApprox_ApproxAFunction construction (the
        // 13-argument ctor running Perform with the user cutting tool).
        let a_approx = ApproxAFunction::with_cut_tool(
            num1dss,
            num2dss,
            num3dss,
            one_d_tol,
            two_d_tol,
            Some(&three_d_tol),
            first,
            last,
            the_order,
            the_max_degree,
            the_max_segments,
            &mut ev,
            &cut_tool,
        );

        // OCCT L169-170.
        self.my_is_done = a_approx.is_done();
        self.my_has_result = a_approx.has_result();

        if self.my_has_result {
            // OCCT L174-178: the result poles / knots / mults / degree.
            let poles_flat = a_approx.poles_flat(1);
            let knots = a_approx.knots_vec();
            let mults = a_approx.multiplicities_vec();
            let degree = a_approx.degree();
            // OCCT L179: new Geom_BSplineCurve(Poles, Knots, Mults, Degree)
            // — the non-periodic (knots, mults) constructor; the rcad
            // carrier stores the flat (expanded) knot vector (the same
            // bijection as `geom::bspline_ops`), unit weights (the AdvApprox
            // 3D result is non-rational).
            let control_points: Vec<DVec3> = poles_flat
                .chunks_exact(3)
                .map(|ch| DVec3::new(ch[0], ch[1], ch[2]))
                .collect();
            let weights = vec![1.0; control_points.len()];
            let mut knots_flat = Vec::with_capacity(
                knots.len() + mults.iter().map(|m| *m as usize).sum::<usize>(),
            );
            for (k, m) in knots.iter().zip(mults.iter()) {
                for _ in 0..*m {
                    knots_flat.push(*k);
                }
            }
            self.my_bspl_curve = Some(BSplineCurve3 {
                degree: degree as usize,
                knots: knots_flat,
                control_points,
                weights,
                is_periodic: false,
            });
            // OCCT L180.
            self.my_max_error = a_approx.max_error_at(3, 1);
        }
    }

    /// OCCT GeomConvert_ApproxCurve::Curve() (cxx L184-187).
    pub fn curve(&self) -> Option<BSplineCurve3> {
        self.my_bspl_curve.clone()
    }

    /// OCCT GeomConvert_ApproxCurve::IsDone() (cxx L189-192).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT GeomConvert_ApproxCurve::HasResult() (cxx L194-197).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT GeomConvert_ApproxCurve::MaxError() (cxx L199-202).
    pub fn max_error(&self) -> f64 {
        self.my_max_error
    }

    /// OCCT GeomConvert_ApproxCurve::Dump(Standard_OStream& o)
    /// (cxx L204-208) — the stream is the caller's String buffer.
    pub fn dump(&self, o: &mut String) {
        o.push_str("******* Dump of ApproxCurve *******\n");
        o.push_str(format!("*******Error   {}\n", self.max_error()).as_str());
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, CurveEval};
    use crate::math::bspl_lib;

    /// The unit circle about the origin in the XY plane (the parameter is
    /// the angle over [0, 2*pi]).
    fn full_circle() -> Curve3 {
        Curve3::Circle(Circle3 {
            center: DVec3::new(0.0, 0.0, 0.0),
            normal: DVec3::new(0.0, 0.0, 1.0),
            x_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: DVec3::new(0.0, 1.0, 0.0),
            radius: 1.0,
        })
    }

    /// The OCCT-caller parameter set (no ctor defaults exist in this OCCT
    /// version): Tol3d = 1e-4, Continuity = GeomAbs_C1, MaxSegments = 50,
    /// MaxDegree = 14.
    fn approx_circle_defaults() -> GeomConvertApproxCurve {
        GeomConvertApproxCurve::new(
            &full_circle(),
            1.0e-4,
            GeomAbsShape::C1,
            50,
            14,
        )
    }

    /// A full circle approximates to a non-rational BSpline of degree <=
    /// MaxDegree = 14 that is at least C1 everywhere and lands on the
    /// circle within the approximation tolerance.
    #[test]
    fn approx_full_circle_lands_on_circle() {
        let appr = approx_circle_defaults();
        assert!(appr.has_result(), "HasResult");
        assert!(appr.is_done(), "IsDone");
        let bs = appr.curve().expect("result curve");
        assert!(!bs.is_periodic, "non-periodic result");
        assert!(bs.degree >= 1 && bs.degree <= 14, "degree {}", bs.degree);
        assert!(
            appr.max_error() <= 1.0e-4,
            "MaxError {} within Tol3d 1e-4",
            appr.max_error()
        );
        println!(
            "circle: degree = {}, poles = {}, knots(regrouped) = {}, max_error = {:e}",
            bs.degree,
            bs.control_points.len(),
            bspline_knots_mults(&bs).0.len(),
            appr.max_error()
        );

        // C1 continuity: every interior knot multiplicity m keeps
        // degree - m >= 1.
        let (knots, mults) = bspline_knots_mults(&bs);
        for &m in mults.iter().skip(1).take(mults.len() - 2) {
            assert!(
                (m as i32) <= bs.degree as i32 - 1,
                "interior mult {m} breaks C1 at degree {}",
                bs.degree
            );
        }

        // Sampled accuracy against the analytic circle.
        for k in 0..=16 {
            let t = std::f64::consts::TAU * (k as f64) / 16.0;
            let p = CurveEval::point_at(&bs, t);
            let r = p.length();
            assert!(
                (r - 1.0).abs() < 1.0e-3,
                "u={t}: radius {r}, want 1.0 (max error {})",
                appr.max_error()
            );
        }

        let mut dump = String::new();
        appr.dump(&mut dump);
        assert!(dump.contains("Dump of ApproxCurve"));
    }

    /// The Adaptor3d_Curve constructor (cxx L118-125) yields the same
    /// approximation as the Geom_Curve constructor.
    #[test]
    fn ctor_from_adaptor_matches_ctor_from_curve() {
        let direct = approx_circle_defaults();
        let via_adaptor = GeomConvertApproxCurve::new_adaptor(
            Arc::new(GeomCurveAdaptor::new(full_circle())),
            1.0e-4,
            GeomAbsShape::C1,
            50,
            14,
        );
        assert_eq!(direct.has_result(), via_adaptor.has_result());
        assert_eq!(direct.is_done(), via_adaptor.is_done());
        assert_eq!(direct.max_error(), via_adaptor.max_error());
        let a = direct.curve().expect("direct curve");
        let b = via_adaptor.curve().expect("adaptor curve");
        assert_eq!(a.degree, b.degree);
        assert_eq!(a.knots, b.knots);
        for (pa, pb) in a.control_points.iter().zip(b.control_points.iter()) {
            assert_eq!(pa, pb);
        }
    }

    /// A C0-kinked input (the `glue` caller case: the concatenated BSpline
    /// with continuity < C1) is segmented at the discontinuity and the
    /// result follows the polyline within tolerance.
    #[test]
    fn approx_c0_kinked_bspline_is_segmented() {
        // Degree-1 clamped BSpline through (0,0,0)-(1,0,0)-(1,1,0): C0 at
        // u = 0.5 only.
        let kinked = Curve3::BSpline(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(1.0, 1.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
            is_periodic: false,
        });
        let appr = GeomConvertApproxCurve::new(&kinked, 1.0e-7, GeomAbsShape::C1, 16, 14);
        assert!(appr.has_result(), "HasResult");
        let bs = appr.curve().expect("result curve");
        assert!(bs.degree <= 14);
        // The C0 join must survive as a full-multiplicity knot (a C1
        // approximation across a corner would cut it off).
        let dist = |p: DVec3, q: DVec3| (p - q).length();
        for k in 0..=8 {
            let t = (k as f64) / 8.0;
            let p = CurveEval::point_at(&bs, t);
            let w = if t <= 0.5 {
                DVec3::new(t * 2.0, 0.0, 0.0)
            } else {
                DVec3::new(1.0, (t - 0.5) * 2.0, 0.0)
            };
            assert!(
                dist(p, w) < 1.0e-3,
                "u={t}: {p:?} want {w:?} (max error {})",
                appr.max_error()
            );
        }
    }

    /// The (knots, mults) regrouping of the flat carrier knot vector (the
    /// AdvApprox knots are bitwise-stable after expansion, so the equality
    /// window is the BSplCLib resolution constant).
    fn bspline_knots_mults(bs: &BSplineCurve3) -> (Vec<f64>, Vec<i32>) {
        let mut knots: Vec<f64> = Vec::new();
        let mut mults: Vec<i32> = Vec::new();
        for &k in &bs.knots {
            if let Some(i) = knots
                .iter()
                .position(|k2| (*k2 - k).abs() <= bspl_lib::GP_RESOLUTION)
            {
                mults[i] += 1;
            } else {
                knots.push(k);
                mults.push(1);
            }
        }
        (knots, mults)
    }
}
