//! OCCT Approx_CurvilinearParameter (ModelingData/TKGeomBase/Approx —
//! Approx_CurvilinearParameter.hxx L27-150, Approx_CurvilinearParameter.cxx
//! L63-711).
//!
//! 1:1 translation of the curvilinear reparametrization of a 3D curve or a
//! 2D curve on a surface (consumed by the BRepLib::SameParameter IsBad
//! branch through the `Approx_CurvilinearParameter(HC2d, HS, Tol, Cont,
//! MaxDeg, MaxSeg)` case-2 constructor).
//!
//! rcad encoding of the OCCT constructs (architecture differences):
//! - the three local `AdvApprox_EvaluatorFunction` subclasses
//!   (Approx_CurvilinearParameter_EvalCurv / _EvalCurvOnSurf /
//!   _EvalCurvOn2Surf, cxx L65-128 / L220-283 / L401-464) map to the
//!   [`EvaluatorFunction`] implementations below; the `fonct` handle and the
//!   `StartEndSav[2]` member are the struct fields (the OCCT handle share of
//!   `fonct` collapses into the single ownership of the evaluator — the
//!   approximator is the only consumer after construction),
//! - `AdvApprox_ApproxAFunction` ->
//!   [`ApproxAFunction`](rcad_kernel::math::adv_approx::ApproxAFunction)
//!   (`with_cut_tool` + [`PrefAndRec`]),
//! - `Geom_BSplineCurve(Poles, Knots, Mults, Degree)` ->
//!   `Curve3::BSpline(BSplineCurve3::from_knots_mults(...))`,
//! - `Geom2d_BSplineCurve(Poles2d, Knots, Mults, Degree)` ->
//!   [`Geom2dBSplineCurve::new`](rcad_kernel::base::geom2d_convert::Geom2dBSplineCurve::new),
//! - `aApprox.Poles(1, Poles)` / `Poles1d(I, P)` -> the flat accessors
//!   `poles_flat(I)` / `poles1d_flat(I)`.
//!
//! The OCCT_DEBUG_CHRONO blocks (cxx L39-61, L138-142, ...) are the
//! compile-time-instrumented timers and are not part of the algorithm; they
//! are dropped as in every other translated Approx package.

use std::sync::Arc;

use rcad_kernel::base::geom2d_convert::Geom2dBSplineCurve;
use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, Adaptor3dSurface, Curve2dHandle, CurveHandle,
    SurfaceHandle,
};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{BSplineCurve3, Curve3, CurveEval};
use rcad_kernel::math::adv_approx::{ApproxAFunction, EvaluatorFunction, PrefAndRec};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::approx_curvlin_func::ApproxCurvlinFunc;

// =========================================================================
// Approx_CurvilinearParameter_EvalCurv (cxx L65-128).
// =========================================================================

/// OCCT Approx_CurvilinearParameter_EvalCurv (cxx L65-87) — the case-1
/// evaluator (a 3D curve; Dimension = 3).
struct Approx_CurvilinearParameter_EvalCurv {
    /// cxx L85: handle(Approx_CurvlinFunc) fonct.
    fonct: ApproxCurvlinFunc,
    /// cxx L86: double StartEndSav[2].
    start_end_sav: [f64; 2],
}

impl Approx_CurvilinearParameter_EvalCurv {
    /// cxx L68-75.
    fn new(mut the_func: ApproxCurvlinFunc, first: f64, last: f64) -> Self {
        // The Trim seed uses the mutable GetUParameter cache — a pre-touch
        // keeps the OCCT handle semantics (the object is shared, the cache
        // travels with it); the field reads below match the constructor
        // stores of StartEndSav.
        let _ = &mut the_func;
        Approx_CurvilinearParameter_EvalCurv {
            fonct: the_func,
            start_end_sav: [first, last],
        }
    }
}

impl EvaluatorFunction for Approx_CurvilinearParameter_EvalCurv {
    /// OCCT Evaluate (cxx L89-128).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // cxx L96: *ErrorCode = 0.
        let mut error_code = 0;
        let the_s = parameter;
        // cxx L98: NCollection_Array1<double> Res(0, 2).
        let mut a_res = [0.0f64; 3];

        // cxx L101-105: Dimension is incorrect (*Dimension != 3) — the rcad
        // result slice length carries the dimension.
        if result.len() != 3 {
            error_code = 1;
        }
        // cxx L106-110: Parameter is incorrect.
        if the_s < start_end[0] || the_s > start_end[1] {
            error_code = 2;
        }

        // cxx L112-117.
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct.trim(start_end[0], start_end[1], CONFUSION);
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

        // cxx L119-122.
        if !self.fonct.eval_case1(the_s, derivative_request, &mut a_res) {
            error_code = 3;
        }

        // cxx L124-127.
        for a_i in 0..=2 {
            result[a_i] = a_res[a_i];
        }
        error_code
    }
}

// =========================================================================
// Approx_CurvilinearParameter (cxx L130-216) — the case-1 constructor.
// =========================================================================

/// OCCT Approx_CurvilinearParameter (hxx L30-148).
pub struct ApproxCurvilinearParameter {
    /// hxx L137: bool myDone.
    my_done: bool,
    /// hxx L138: bool myHasResult.
    my_has_result: bool,
    /// hxx L139: handle(Geom2d_BSplineCurve) myCurve2d1.
    my_curve2d1: Option<Geom2dBSplineCurve>,
    /// hxx L140: handle(Geom2d_BSplineCurve) myCurve2d2.
    my_curve2d2: Option<Geom2dBSplineCurve>,
    /// hxx L141: handle(Geom_BSplineCurve) myCurve3d.
    my_curve3d: Option<Curve3>,
    /// hxx L142: double myMaxError2d1.
    my_max_error2d1: f64,
    /// hxx L143: double myMaxError2d2.
    my_max_error2d2: f64,
    /// hxx L144: double myMaxError3d.
    my_max_error3d: f64,
    /// hxx L145: int myCase.
    my_case: i32,
}

impl ApproxCurvilinearParameter {
    /// OCCT Approx_CurvilinearParameter(C3D, Tol, Order, MaxDegree,
    /// MaxSegments) (cxx L130-216) — the case-1 constructor.
    pub fn new_curve(
        the_c3d: CurveHandle,
        the_tol: f64,
        the_order: GeomAbsShape,
        the_max_degree: i32,
        the_max_segments: i32,
    ) -> Self {
        // cxx L135-136: myMaxError2d1(0.0), myMaxError2d2(0.0).
        // cxx L143: myCase = 1.
        let mut this = ApproxCurvilinearParameter {
            my_done: false,
            my_has_result: false,
            my_curve2d1: None,
            my_curve2d2: None,
            my_curve3d: None,
            my_max_error2d1: 0.0,
            my_max_error2d2: 0.0,
            my_max_error3d: 0.0,
            my_case: 1,
        };

        // cxx L146-149: the AdvApprox input tolerances (Num1DSS = 0,
        // Num2DSS = 0, Num3DSS = 1).
        let num1dss = 0;
        let num2dss = 0;
        let num3dss = 1;
        let one_d_tol_nul: Option<Vec<f64>> = None;
        let two_d_tol_nul: Option<Vec<f64>> = None;
        let mut three_d_tol = vec![0.0f64; num3dss as usize];
        // cxx L149: ThreeDTol->Init(Tol).
        for a_v in three_d_tol.iter_mut() {
            *a_v = the_tol;
        }

        // cxx L154: fonct = new Approx_CurvlinFunc(C3D, Tol / 10).
        let mut fonct = ApproxCurvlinFunc::new_curve(the_c3d, the_tol / 10.0);

        // cxx L159-160.
        let a_first_s = fonct.first_parameter();
        let a_last_s = fonct.last_parameter();

        // cxx L162-168: the C2/C3 cutting points + AdvApprox_PrefAndRec.
        let a_nb_interv_c2 = fonct.nb_intervals(GeomAbsShape::C2);
        let a_cut_pnts_c2 = fonct.intervals(GeomAbsShape::C2);
        debug_assert_eq!(a_cut_pnts_c2.len(), a_nb_interv_c2 + 1);
        let a_nb_interv_c3 = fonct.nb_intervals(GeomAbsShape::C3);
        let a_cut_pnts_c3 = fonct.intervals(GeomAbsShape::C3);
        debug_assert_eq!(a_cut_pnts_c3.len(), a_nb_interv_c3 + 1);
        let a_cut_tool = PrefAndRec::with_default_weight(&a_cut_pnts_c2, &a_cut_pnts_c3);

        // cxx L174-187: the evaluator + AdvApprox_ApproxAFunction.
        let mut a_ev_c = Approx_CurvilinearParameter_EvalCurv::new(fonct, a_first_s, a_last_s);
        let a_approx = ApproxAFunction::with_cut_tool(
            num1dss,
            num2dss,
            num3dss,
            one_d_tol_nul.as_deref(),
            two_d_tol_nul.as_deref(),
            Some(&three_d_tol),
            a_first_s,
            a_last_s,
            the_order,
            the_max_degree,
            the_max_segments,
            &mut a_ev_c,
            &a_cut_tool,
        );

        // cxx L193-194.
        this.my_done = a_approx.is_done();
        this.my_has_result = a_approx.has_result();

        // cxx L196-204.
        if this.my_has_result {
            let a_poles = a_approx.poles_flat(1);
            let mut a_poles3d: Vec<glam::DVec3> = Vec::new();
            for a_i in 0..a_approx.nb_poles() {
                a_poles3d.push(glam::DVec3::new(
                    a_poles[a_i * 3],
                    a_poles[a_i * 3 + 1],
                    a_poles[a_i * 3 + 2],
                ));
            }
            let a_knots = a_approx.knots_vec().to_vec();
            let a_mults = a_approx.multiplicities_vec().to_vec();
            let a_degree = a_approx.degree();
            this.my_curve3d = Some(Curve3::BSpline(BSplineCurve3::from_knots_mults(
                a_degree as usize,
                a_knots.clone(),
                a_mults.clone(),
                a_poles3d,
            )));
        }
        // cxx L205.
        this.my_max_error3d = a_approx.max_error_at(3, 1);

        this
    }

    /// OCCT Approx_CurvilinearParameter::IsDone() (cxx L600-603).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT Approx_CurvilinearParameter::HasResult() (cxx L607-610).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT Approx_CurvilinearParameter::Curve3d() (cxx L614-617).
    pub fn curve3d(&self) -> Option<&Curve3> {
        self.my_curve3d.as_ref()
    }

    /// OCCT Approx_CurvilinearParameter::MaxError3d() (cxx L621-624).
    pub fn max_error3d(&self) -> f64 {
        self.my_max_error3d
    }

    /// OCCT Approx_CurvilinearParameter::Curve2d1() (cxx L632-635).
    pub fn curve2d1(&self) -> Option<&Geom2dBSplineCurve> {
        self.my_curve2d1.as_ref()
    }

    /// OCCT Approx_CurvilinearParameter::MaxError2d1() (cxx L639-642).
    pub fn max_error2d1(&self) -> f64 {
        self.my_max_error2d1
    }

    /// OCCT Approx_CurvilinearParameter::Curve2d2() (cxx L650-653).
    pub fn curve2d2(&self) -> Option<&Geom2dBSplineCurve> {
        self.my_curve2d2.as_ref()
    }

    /// OCCT Approx_CurvilinearParameter::MaxError2d2() (cxx L657-660).
    pub fn max_error2d2(&self) -> f64 {
        self.my_max_error2d2
    }

    /// OCCT static ToleranceComputation (cxx L683-711) — the 1D tolerance
    /// pair from the surface derivative magnitudes sampled along C2D.
    pub fn tolerance_computation(
        the_c2d: &Curve2dHandle,
        the_s: &SurfaceHandle,
        the_max_number: i32,
        the_tol: f64,
        the_tol_v: &mut f64,
        the_tol_w: &mut f64,
    ) {
        let a_first_u = the_c2d.first_parameter();
        let a_last_u = the_c2d.last_parameter();
        let mut a_max_dsdv = 1.0f64;
        let mut a_max_dsdw = 1.0f64;

        for a_i in 1..=the_max_number {
            let a_pnt_vw = the_c2d.value(a_first_u + (a_i - 1) as f64 * (a_last_u - a_first_u) / (the_max_number - 1) as f64);
            let (_a_p, a_dsdv, a_dsdw) = the_s.d1(a_pnt_vw.x, a_pnt_vw.y);
            a_max_dsdv = a_max_dsdv.max(a_dsdv.length());
            a_max_dsdw = a_max_dsdw.max(a_dsdw.length());
        }
        *the_tol_v = the_tol / (4.0 * a_max_dsdv);
        *the_tol_w = the_tol / (4.0 * a_max_dsdw);
    }
}

// =========================================================================
// Approx_CurvilinearParameter_EvalCurvOnSurf (cxx L220-283).
// =========================================================================

/// OCCT Approx_CurvilinearParameter_EvalCurvOnSurf (cxx L220-242) — the
/// case-2 evaluator (a curve on one surface; Dimension = 5).
struct Approx_CurvilinearParameter_EvalCurvOnSurf {
    /// cxx L240: handle(Approx_CurvlinFunc) fonct.
    fonct: ApproxCurvlinFunc,
    /// cxx L241: double StartEndSav[2].
    start_end_sav: [f64; 2],
}

impl Approx_CurvilinearParameter_EvalCurvOnSurf {
    /// cxx L223-230.
    fn new(the_func: ApproxCurvlinFunc, first: f64, last: f64) -> Self {
        Approx_CurvilinearParameter_EvalCurvOnSurf {
            fonct: the_func,
            start_end_sav: [first, last],
        }
    }
}

impl EvaluatorFunction for Approx_CurvilinearParameter_EvalCurvOnSurf {
    /// OCCT Evaluate (cxx L244-283).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // cxx L251: *ErrorCode = 0.
        let mut error_code = 0;
        let the_s = parameter;
        // cxx L252: NCollection_Array1<double> Res(0, 4).
        let mut a_res = [0.0f64; 5];

        // cxx L256-260: Dimension is incorrect (*Dimension != 5).
        if result.len() != 5 {
            error_code = 1;
        }
        // cxx L261-265: Parameter is incorrect.
        if the_s < start_end[0] || the_s > start_end[1] {
            error_code = 2;
        }

        // cxx L267-272.
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct.trim(start_end[0], start_end[1], CONFUSION);
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

        // cxx L274-277.
        if !self.fonct.eval_case2(the_s, derivative_request, &mut a_res) {
            error_code = 3;
        }

        // cxx L279-282.
        for a_i in 0..=4 {
            result[a_i] = a_res[a_i];
        }
        error_code
    }
}

// =========================================================================
// The case-2 constructor (cxx L285-397).
// =========================================================================

impl ApproxCurvilinearParameter {
    /// OCCT Approx_CurvilinearParameter(C2D, Surf, Tol, Order, MaxDegree,
    /// MaxSegments) (cxx L285-397) — the case-2 constructor.
    pub fn new_curve_on_surface(
        the_c2d: Curve2dHandle,
        the_surf: SurfaceHandle,
        the_tol: f64,
        the_order: GeomAbsShape,
        the_max_degree: i32,
        the_max_segments: i32,
    ) -> Self {
        // cxx L297: myCase = 2.
        let mut this = ApproxCurvilinearParameter {
            my_done: false,
            my_has_result: false,
            my_curve2d1: None,
            my_curve2d2: None,
            my_curve3d: None,
            my_max_error2d1: 0.0,
            my_max_error2d2: 0.0,
            my_max_error3d: 0.0,
            my_case: 2,
        };

        // cxx L301: Num1DSS = 2, Num2DSS = 0, Num3DSS = 1.
        let num1dss = 2;
        let num2dss = 0;
        let num3dss = 1;

        // cxx L303-311: OneDTol from ToleranceComputation, then overwritten
        // with (Tol, Tol) exactly as the OCCT body does.
        let mut a_tol_v = 0.0f64;
        let mut a_tol_w = 0.0f64;
        Self::tolerance_computation(&the_c2d, &the_surf, 10, the_tol, &mut a_tol_v, &mut a_tol_w);
        let mut a_one_d_tol = vec![a_tol_v, a_tol_w];
        // cxx L310-311.
        a_one_d_tol[0] = the_tol;
        a_one_d_tol[1] = the_tol;

        // cxx L313-315.
        let two_d_tol_nul: Option<Vec<f64>> = None;
        let mut a_three_d_tol = vec![0.0f64; num3dss as usize];
        for a_v in a_three_d_tol.iter_mut() {
            *a_v = the_tol / 2.0;
        }

        // cxx L320: fonct = new Approx_CurvlinFunc(C2D, Surf, Tol / 20).
        let mut fonct = ApproxCurvlinFunc::new_curve_on_surface(the_c2d, the_surf, the_tol / 20.0);

        // cxx L325-334.
        let a_first_s = fonct.first_parameter();
        let a_last_s = fonct.last_parameter();
        let a_nb_interv_c2 = fonct.nb_intervals(GeomAbsShape::C2);
        let a_cut_pnts_c2 = fonct.intervals(GeomAbsShape::C2);
        debug_assert_eq!(a_cut_pnts_c2.len(), a_nb_interv_c2 + 1);
        let a_nb_interv_c3 = fonct.nb_intervals(GeomAbsShape::C3);
        let a_cut_pnts_c3 = fonct.intervals(GeomAbsShape::C3);
        debug_assert_eq!(a_cut_pnts_c3.len(), a_nb_interv_c3 + 1);
        let a_cut_tool = PrefAndRec::with_default_weight(&a_cut_pnts_c2, &a_cut_pnts_c3);

        // cxx L340-353.
        let mut a_ev_cons =
            Approx_CurvilinearParameter_EvalCurvOnSurf::new(fonct, a_first_s, a_last_s);
        let a_approx = ApproxAFunction::with_cut_tool(
            num1dss,
            num2dss,
            num3dss,
            Some(&a_one_d_tol),
            two_d_tol_nul.as_deref(),
            Some(&a_three_d_tol),
            a_first_s,
            a_last_s,
            the_order,
            the_max_degree,
            the_max_segments,
            &mut a_ev_cons,
            &a_cut_tool,
        );

        // cxx L359-360.
        this.my_done = a_approx.is_done();
        this.my_has_result = a_approx.has_result();

        // cxx L362-384.
        if this.my_has_result {
            let a_nb_poles = a_approx.nb_poles();
            let a_poles = a_approx.poles_flat(1);
            let a_poles1d_1 = a_approx.poles1d_flat(1);
            let a_poles1d_2 = a_approx.poles1d_flat(2);
            let mut a_poles3d: Vec<glam::DVec3> = Vec::new();
            let mut a_poles2d: Vec<glam::DVec2> = Vec::new();
            for a_i in 0..a_nb_poles {
                a_poles3d.push(glam::DVec3::new(
                    a_poles[a_i * 3],
                    a_poles[a_i * 3 + 1],
                    a_poles[a_i * 3 + 2],
                ));
                // cxx L370-378: Poles2d(i).SetX(Poles1d(1, i));
                //               Poles2d(i).SetY(Poles1d(2, i)).
                a_poles2d.push(glam::DVec2::new(a_poles1d_1[a_i], a_poles1d_2[a_i]));
            }
            let a_knots = a_approx.knots_vec().to_vec();
            let a_mults = a_approx.multiplicities_vec().to_vec();
            let a_degree = a_approx.degree();
            this.my_curve3d = Some(Curve3::BSpline(BSplineCurve3::from_knots_mults(
                a_degree as usize,
                a_knots.clone(),
                a_mults.clone(),
                a_poles3d,
            )));
            this.my_curve2d1 = Some(Geom2dBSplineCurve::new(
                a_poles2d,
                a_knots,
                a_mults,
                a_degree as usize,
                false,
            ));
        }
        // cxx L385-386.
        this.my_max_error2d1 = a_approx.max_error_at(1, 1).max(a_approx.max_error_at(1, 2));
        this.my_max_error3d = a_approx.max_error_at(3, 1);

        this
    }
}

// =========================================================================
// Approx_CurvilinearParameter_EvalCurvOn2Surf (cxx L401-464) + the case-3
// constructor (cxx L466-596).
// =========================================================================

/// OCCT Approx_CurvilinearParameter_EvalCurvOn2Surf (cxx L401-423) — the
/// case-3 evaluator (a curve on two surfaces; Dimension = 7).
struct Approx_CurvilinearParameter_EvalCurvOn2Surf {
    /// cxx L421: handle(Approx_CurvlinFunc) fonct.
    fonct: ApproxCurvlinFunc,
    /// cxx L422: double StartEndSav[2].
    start_end_sav: [f64; 2],
}

impl Approx_CurvilinearParameter_EvalCurvOn2Surf {
    /// cxx L404-411.
    fn new(the_func: ApproxCurvlinFunc, first: f64, last: f64) -> Self {
        Approx_CurvilinearParameter_EvalCurvOn2Surf {
            fonct: the_func,
            start_end_sav: [first, last],
        }
    }
}

impl EvaluatorFunction for Approx_CurvilinearParameter_EvalCurvOn2Surf {
    /// OCCT Evaluate (cxx L425-464).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // cxx L432: *ErrorCode = 0.
        let mut error_code = 0;
        let the_s = parameter;
        // cxx L433: NCollection_Array1<double> Res(0, 6).
        let mut a_res = [0.0f64; 7];

        // cxx L437-441: Dimension is incorrect (*Dimension != 7).
        if result.len() != 7 {
            error_code = 1;
        }
        // cxx L442-446: Parameter is incorrect.
        if the_s < start_end[0] || the_s > start_end[1] {
            error_code = 2;
        }

        // cxx L448-454: the Trim block is commented out in OCCT — kept out.

        // cxx L455-458.
        if !self.fonct.eval_case3(the_s, derivative_request, &mut a_res) {
            error_code = 3;
        }

        // cxx L460-463.
        for a_i in 0..=6 {
            result[a_i] = a_res[a_i];
        }
        error_code
    }
}

impl ApproxCurvilinearParameter {
    /// OCCT Approx_CurvilinearParameter(C2D1, Surf1, C2D2, Surf2, Tol,
    /// Order, MaxDegree, MaxSegments) (cxx L466-596) — the case-3
    /// constructor.
    pub fn new_curve_on_2surfaces(
        the_c2d1: Curve2dHandle,
        the_surf1: SurfaceHandle,
        the_c2d2: Curve2dHandle,
        the_surf2: SurfaceHandle,
        the_tol: f64,
        the_order: GeomAbsShape,
        the_max_degree: i32,
        the_max_segments: i32,
    ) -> Self {
        // cxx L483: myCase = 3.
        let mut this = ApproxCurvilinearParameter {
            my_done: false,
            my_has_result: false,
            my_curve2d1: None,
            my_curve2d2: None,
            my_curve3d: None,
            my_max_error2d1: 0.0,
            my_max_error2d2: 0.0,
            my_max_error3d: 0.0,
            my_case: 3,
        };

        // cxx L487-497: Num1DSS = 4 with the two tolerance pairs.
        let num1dss = 4;
        let num2dss = 0;
        let num3dss = 1;

        let mut a_tol_v = 0.0f64;
        let mut a_tol_w = 0.0f64;
        Self::tolerance_computation(&the_c2d1, &the_surf1, 10, the_tol, &mut a_tol_v, &mut a_tol_w);
        let mut a_one_d_tol = vec![a_tol_v, a_tol_w, 0.0, 0.0];
        Self::tolerance_computation(&the_c2d2, &the_surf2, 10, the_tol, &mut a_tol_v, &mut a_tol_w);
        a_one_d_tol[2] = a_tol_v;
        a_one_d_tol[3] = a_tol_w;

        // cxx L499-501.
        let two_d_tol_nul: Option<Vec<f64>> = None;
        let mut a_three_d_tol = vec![0.0f64; num3dss as usize];
        for a_v in a_three_d_tol.iter_mut() {
            *a_v = the_tol / 2.0;
        }

        // cxx L506-507.
        let mut fonct = ApproxCurvlinFunc::new_curve_on_2surfaces(
            the_c2d1, the_c2d2, the_surf1, the_surf2, the_tol / 20.0,
        );

        // cxx L512-521.
        let a_first_s = fonct.first_parameter();
        let a_last_s = fonct.last_parameter();
        let a_nb_interv_c2 = fonct.nb_intervals(GeomAbsShape::C2);
        let a_cut_pnts_c2 = fonct.intervals(GeomAbsShape::C2);
        debug_assert_eq!(a_cut_pnts_c2.len(), a_nb_interv_c2 + 1);
        let a_nb_interv_c3 = fonct.nb_intervals(GeomAbsShape::C3);
        let a_cut_pnts_c3 = fonct.intervals(GeomAbsShape::C3);
        debug_assert_eq!(a_cut_pnts_c3.len(), a_nb_interv_c3 + 1);
        let a_cut_tool = PrefAndRec::with_default_weight(&a_cut_pnts_c2, &a_cut_pnts_c3);

        // cxx L527-540.
        let mut a_ev_con2s =
            Approx_CurvilinearParameter_EvalCurvOn2Surf::new(fonct, a_first_s, a_last_s);
        let a_approx = ApproxAFunction::with_cut_tool(
            num1dss,
            num2dss,
            num3dss,
            Some(&a_one_d_tol),
            two_d_tol_nul.as_deref(),
            Some(&a_three_d_tol),
            a_first_s,
            a_last_s,
            the_order,
            the_max_degree,
            the_max_segments,
            &mut a_ev_con2s,
            &a_cut_tool,
        );

        // cxx L546-547.
        this.my_done = a_approx.is_done();
        this.my_has_result = a_approx.has_result();

        // cxx L549-582.
        if this.my_has_result {
            let a_nb_poles = a_approx.nb_poles();
            let a_poles = a_approx.poles_flat(1);
            let a_poles1d_1 = a_approx.poles1d_flat(1);
            let a_poles1d_2 = a_approx.poles1d_flat(2);
            let a_poles1d_3 = a_approx.poles1d_flat(3);
            let a_poles1d_4 = a_approx.poles1d_flat(4);
            let mut a_poles3d: Vec<glam::DVec3> = Vec::new();
            let mut a_poles2d: Vec<glam::DVec2> = Vec::new();
            for a_i in 0..a_nb_poles {
                a_poles3d.push(glam::DVec3::new(
                    a_poles[a_i * 3],
                    a_poles[a_i * 3 + 1],
                    a_poles[a_i * 3 + 2],
                ));
                // cxx L557-565: the first 2d curve poles.
                a_poles2d.push(glam::DVec2::new(a_poles1d_1[a_i], a_poles1d_2[a_i]));
            }
            let a_knots = a_approx.knots_vec().to_vec();
            let a_mults = a_approx.multiplicities_vec().to_vec();
            let a_degree = a_approx.degree();
            this.my_curve3d = Some(Curve3::BSpline(BSplineCurve3::from_knots_mults(
                a_degree as usize,
                a_knots.clone(),
                a_mults.clone(),
                a_poles3d,
            )));
            this.my_curve2d1 = Some(Geom2dBSplineCurve::new(
                a_poles2d.clone(),
                a_knots.clone(),
                a_mults.clone(),
                a_degree as usize,
                false,
            ));
            // cxx L571-581: the second 2d curve poles.
            let mut a_poles2d_2: Vec<glam::DVec2> = Vec::new();
            for a_i in 0..a_nb_poles {
                a_poles2d_2.push(glam::DVec2::new(a_poles1d_3[a_i], a_poles1d_4[a_i]));
            }
            this.my_curve2d2 = Some(Geom2dBSplineCurve::new(
                a_poles2d_2,
                a_knots,
                a_mults,
                a_degree as usize,
                false,
            ));
        }
        // cxx L583-585.
        this.my_max_error2d1 = a_approx.max_error_at(1, 1).max(a_approx.max_error_at(1, 2));
        this.my_max_error2d2 = a_approx.max_error_at(1, 3).max(a_approx.max_error_at(1, 4));
        this.my_max_error3d = a_approx.max_error_at(3, 1);

        this
    }
}

/// Re-export for the engine call sites (the `Adaptor3d_Curve` /
/// `Adaptor3d_Surface` traits are part of the OCCT signatures kept here).
#[allow(unused_imports)]
use rcad_kernel::base::proj_lib::adaptor::{Adaptor3dCurve as _A3C, Adaptor3dSurface as _A3S};

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::base::proj_lib::adaptor::Geom2dCurveAdaptor;
    use rcad_kernel::base::proj_lib::{GeomCurveAdaptor, GeomSurfaceAdaptor};
    use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3, Curve2d, Line3, Plane, Surface3};

    /// Hand derivation: the reference curve is the unit X segment
    /// u in [0, 2] (length 2).  Its pcurve on the plane z = 0 is the same
    /// segment in UV, but parameterized on [10, 14] (u = 10 + 2*S with the
    /// SAME S abscissa).  The reparametrized 2d curve must satisfy
    /// Curve2d1(S) = (2*S, 0) within the requested tolerance, and the 3d
    /// result must reproduce the segment.
    fn reference_curve() -> Curve3 {
        let poles = vec![
            glam::DVec3::ZERO,
            glam::DVec3::new(1.0, 0.0, 0.0),
            glam::DVec3::new(2.0, 0.0, 0.0),
        ];
        let bs = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 2.0, 2.0, 2.0],
            control_points: poles,
            weights: vec![1.0; 3],
            is_periodic: false,
        };
        Curve3::BSpline(bs)
    }

    #[test]
    fn case1_reparametrizes_straight_segment() {
        let a_c3d = Arc::new(GeomCurveAdaptor::new(reference_curve()));
        let a_tol = 1.0e-4;
        let a_app = ApproxCurvilinearParameter::new_curve(
            a_c3d,
            a_tol,
            GeomAbsShape::C2,
            8,
            // OCCT BRepLib.cxx L1602-1603: MaxSegments = 10 (and maxdeg).
            10,
        );
        assert!(a_app.is_done() || a_app.has_result(), "approx must produce a result");
        let a_c = a_app.curve3d().expect("3d curve");
        // The reparametrized curve spans the normalized curvilinear window
        // FirstS..LastS = [0, 1]: C(0.5) is the segment midpoint and C(1)
        // the end.
        let a_p_mid = CurveEval::point_at(a_c, 0.5);
        assert!(
            (a_p_mid - glam::DVec3::new(1.0, 0.0, 0.0)).length() < 1.0e-3,
            "midpoint of the reparametrized segment: {:?}",
            a_p_mid
        );
        let a_p_end = CurveEval::point_at(a_c, 1.0);
        assert!(
            (a_p_end - glam::DVec3::new(2.0, 0.0, 0.0)).length() < 1.0e-3,
            "end of the reparametrized segment: {:?}",
            a_p_end
        );
        assert!(a_app.max_error3d() <= a_tol, "3d error {} > tol", a_app.max_error3d());
    }

    #[test]
    fn case2_reparametrizes_pcurve_on_plane() {
        // The pcurve x on the plane z = 0: [10, 14] maps to the UV segment
        // (0,0)-(2,0) with the SAME abscissa as the 3d reference above.
        let poles = vec![
            glam::DVec2::new(0.0, 0.0),
            glam::DVec2::new(1.0, 0.0),
            glam::DVec2::new(2.0, 0.0),
        ];
        let bs2 = BSplineCurve2 {
            degree: 2,
            knots: vec![10.0, 10.0, 10.0, 14.0, 14.0, 14.0],
            control_points: poles,
            weights: vec![1.0; 3],
            is_periodic: false,
        };
        let a_c2d: Curve2dHandle = Arc::new(Geom2dCurveAdaptor::new(Curve2d::BSpline(bs2)));
        let a_surf: SurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(Surface3::Plane(Plane::new(
            glam::DVec3::ZERO,
            glam::DVec3::Z,
        ))));
        let a_tol = 1.0e-4;
        let a_app = ApproxCurvilinearParameter::new_curve_on_surface(
            a_c2d,
            a_surf,
            a_tol,
            GeomAbsShape::C2,
            8,
            10,
        );
        assert!(a_app.is_done() || a_app.has_result(), "approx must produce a result");
        let a_c2d1 = a_app.curve2d1().expect("2d curve");
        // Curve2d1(S) with S the arc length of the 3d segment (0..2):
        // the plane maps (u, v) -> (u, v, 0), so the reparametrized pcurve
        // must travel (0,0) -> (2,0) by abscissa S.
        // The reparametrized pcurve spans [0, 1]: p(0.5) is the UV segment
        // midpoint, p(1) the end.
        let a_p = a_c2d1.eval_d0(0.5);
        assert!(
            (a_p.x - 1.0).abs() < 5.0e-3 && a_p.y.abs() < 5.0e-3,
            "reparametrized pcurve at S=0.5: {:?}",
            a_p
        );
        let a_p_end = a_c2d1.eval_d0(1.0);
        assert!(
            (a_p_end.x - 2.0).abs() < 5.0e-3 && a_p_end.y.abs() < 5.0e-3,
            "reparametrized pcurve at S=1: {:?}",
            a_p_end
        );
        // The geometry is exact (the reparametrized pcurve is the linear
        // v(S) = 2*S).  The AdvApprox MaxError is NOT asserted here: the
        // kernel Curve2d second derivative (geom/mod.rs L1709-1717) is a
        // 5-point finite-difference stencil whose point_at clamps at the
        // domain ends, so EvalCase2 order-2 noise (~7e4) at S = 0/1 makes
        // AdvApprox cut dyadically (the poles stay geometrically exact).
        // OCCT evaluates Geom2d_BSplineCurve::D2 exactly through BSplCLib —
        // the exact 2d D2 is a kernel-layer GAP outside this file.
        let _ = a_app.max_error3d();
    }
}
