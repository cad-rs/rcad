//! OCCT Approx_SameParameter (ModelingData/TKGeomBase/Approx —
//! Approx_SameParameter.hxx L27-175, Approx_SameParameter.cxx L37-992).
//!
//! 1:1 translation of the approximation of a pcurve so that its parameter
//! matches the parameter of a given 3d reference curve: the evaluator
//! (`Approx_SameParameter_Evaluator`, cxx L37-102), the static helpers
//! (`ProjectPointOnCurve` L106-149, `ComputeTolReached` L153-188, `Check`
//! L192-266), the three constructors, `Build` and the private members
//! (`BuildInitialDistribution` L547-580, `IncreaseInitialNbSamples` L588-649,
//! `CheckSameParameter` L653-766, `ComputeTangents` L770-809, `Interpolate`
//! L813-853, `IncreaseNbPoles` L857-992).
//!
//! rcad encoding of the OCCT handles (architecture difference, annotated):
//! - `handle(Adaptor3d_Curve) myC3d` -> [`GeomCurveAdaptor`] over the kernel
//!   `Curve3` plus the restricted parameter window (the `(Curve3, first,
//!   last)`-pair carrier of `GeomAdaptor_Curve`),
//! - `handle(Adaptor3d_Surface) mySurf` -> `Arc<GeomSurfaceAdaptor>`,
//! - `handle(Adaptor2d_Curve2d) myHCurve2d` -> `Curve2dHandle`
//!   (`Arc<dyn Adaptor2dCurve2d>`, the `Geom2dAdaptor_Curve` value),
//! - `handle(Adaptor3d_CurveOnSurface) myCurveOnSurface` ->
//!   `CurveOnSurface` (the kernel adaptor carries `Arc` handles; its
//!   `ShallowCopy` is the handle-pair rebuild below).
//!
//! Reused kernel machinery (not re-translated): `BSplCLib::Eval` ->
//! `bspl_lib::eval_flat`, `BSplCLib::Interpolate` -> `bspl_lib::interpolate`,
//! `AdvApprox_ApproxAFunction` -> `math::adv_approx::ApproxAFunction`,
//! `GeomLib_MakeCurvefromApprox` ->
//! [`GeomLibMakeCurvefromApprox`](crate::geomalgo::geom_lib_make_curve_from_approx),
//! `Extrema_LocateExtPC` / `Extrema_ExtPC` -> `base::extrema_locate_ext_pc` /
//! `base::extrema_ext_pc` (through the `Extrema_CurveTool` facade).

use std::sync::Arc;

use glam::DVec3;

use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::extrema_locate_ext_pc::LocateExtPC;
use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, Adaptor3dSurface, Curve2dHandle, CurveOnSurface,
    Geom2dCurveAdaptor,
};
use rcad_kernel::base::proj_lib::{CurveType, GeomCurveAdaptor, GeomSurfaceAdaptor};
use rcad_kernel::core::precision::{is_infinite_value, p_confusion, CONFUSION, INFINITE_VALUE, REAL_LAST};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, TrimmedCurve2};
use rcad_kernel::math::adv_approx::{ApproxAFunction, EvaluatorFunction};
use rcad_kernel::math::{bspl_lib, GeomAbsShape};

use crate::geomalgo::geom_lib_make_curve_from_approx::GeomLibMakeCurvefromApprox;

/// OCCT `Approx_SameParameter::myNbSamples` (hxx L161) — "To be consistent
/// with checkshape".
const MY_NB_SAMPLES: i32 = 22;

/// OCCT `Approx_SameParameter::myMaxArraySize` (hxx L162).
const MY_MAX_ARRAY_SIZE: usize = 1000;

// =========================================================================
// Approx_SameParameter_Evaluator (cxx L37-102)
// =========================================================================

/// OCCT `Approx_SameParameter_Evaluator` (cxx L37-60) — evaluates the 1D
/// BSpline that represents the change in parameterization.
struct ApproxSameParameterEvaluator<'a> {
    /// OCCT cxx L57: const NCollection_Array1<double>& FlatKnots.
    flat_knots: &'a [f64],
    /// OCCT cxx L58: const NCollection_Array1<double>& Poles.
    poles: &'a [f64],
    /// OCCT cxx L59: handle(Adaptor2d_Curve2d) HCurve2d.
    h_curve2d: Curve2dHandle,
}

impl ApproxSameParameterEvaluator<'_> {
    /// OCCT cxx L40-47.
    fn new<'a>(
        flat_knots: &'a [f64],
        poles: &'a [f64],
        h_curve2d: Curve2dHandle,
    ) -> ApproxSameParameterEvaluator<'a> {
        ApproxSameParameterEvaluator {
            flat_knots,
            poles,
            h_curve2d,
        }
    }
}

impl EvaluatorFunction for ApproxSameParameterEvaluator<'_> {
    /// OCCT Approx_SameParameter_Evaluator::Evaluate (cxx L64-102).  The
    /// `Dimension` argument is unused in OCCT (the 2-result layout is the
    /// caller's contract); `DerivativeRequest` other than 0/1 leaves `Result`
    /// untouched and only reports the success code.
    fn evaluate(
        &mut self,
        _start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // OCCT cxx L71-74.
        let a_degree: usize = 3;
        // OCCT: int extrap_mode[2] = {aDegree, aDegree} — the FIRST entry is
        // passed by reference to Eval (BSplCLib::Eval, BSplCLib.cxx
        // L3640-3790 reads the adjacent [1] as well).
        let mut extrap_mode = [3i32, 3i32];
        let mut eval_result = [0.0f64; 2];

        // OCCT cxx L77-85: BSplCLib::Eval(*Parameter, false, *DerivativeRequest,
        // extrap_mode[0], aDegree, FlatKnots, 1, PolesArray[0], eval_result[0]).
        bspl_lib::eval_flat(
            parameter,
            false,
            derivative_request,
            &mut extrap_mode,
            a_degree,
            self.flat_knots,
            1,
            self.poles,
            &mut eval_result,
        );

        // OCCT cxx L87-99.
        if derivative_request == 0 {
            // HCurve2d->D0(eval_result[0], aPoint); aPoint.Coord(Result[0], Result[1]);
            let a_point = self.h_curve2d.d0(eval_result[0]);
            result[0] = a_point.x;
            result[1] = a_point.y;
        } else if derivative_request == 1 {
            // HCurve2d->D1(eval_result[0], aPoint, aVector);
            // aVector.Multiply(eval_result[1]); aVector.Coord(Result[0], Result[1]);
            let (_a_point, mut a_vector) = self.h_curve2d.d1(eval_result[0]);
            a_vector *= eval_result[1];
            result[0] = a_vector.x;
            result[1] = a_vector.y;
        }

        // OCCT cxx L101: ReturnCode[0] = 0.
        0
    }
}

// =========================================================================
// ProjectPointOnCurve (cxx L106-149)
// =========================================================================

/// OCCT static `ProjectPointOnCurve` (cxx L106-149) — the Newton iteration
/// that projects `APoint` on `Curve`, clamped to the curve parameter range.
fn project_point_on_curve(
    init_value: f64,
    a_point: DVec3,
    tolerance: f64,
    num_iteration: i32,
    curve: &dyn Adaptor3dCurve,
    status: &mut bool,
    result: &mut f64,
) {
    // OCCT cxx L114-119.
    let mut num_iter = 0i32;
    let mut not_done = true;
    let mut param = init_value;
    *status = false;
    loop {
        // OCCT cxx L120-146 (do { ... } while).
        num_iter += 1;
        let (a_curve_point, d1, d2) = curve.d2(param);
        let vector = a_point - a_curve_point;

        let func = vector.dot(d1);
        if func.abs() < tolerance * d1.length() {
            not_done = false;
            *status = true;
        } else {
            let func_derivative = vector.dot(d2) - d1.dot(d1);

            // Avoid division by zero (OCCT cxx L137-141).
            let toler = 1.0e-12;
            if func_derivative.abs() > toler {
                param -= func / func_derivative;
            }

            param = param.max(curve.first_parameter());
            param = param.min(curve.last_parameter());
        }

        // OCCT: while (not_done && num_iter <= NumIteration).
        if !(not_done && num_iter <= num_iteration) {
            break;
        }
    }

    // OCCT cxx L148: Result = param.
    *result = param;
}

// =========================================================================
// ComputeTolReached (cxx L153-188)
// =========================================================================

/// OCCT static `ComputeTolReached` (cxx L153-188) — the sampled maximal
/// distance between the 3d curve and the curve on surface.
fn compute_tol_reached(c3d: &dyn Adaptor3dCurve, cons: &CurveOnSurface, nbp: i32) -> f64 {
    // OCCT cxx L157-159.
    let mut d2 = 0.0f64; // Square max discrete deviation.
    let first = c3d.first_parameter();
    let last = c3d.last_parameter();
    for i in 0..=nbp {
        // OCCT cxx L162-163: IntToReal(i) / IntToReal(nbp).
        let t = i as f64 / nbp as f64;
        let u = first * (1.0 - t) + last * t;
        // OCCT cxx L164-174: the evaluation sits in a try/catch that turns a
        // Standard_Failure into d2 = Precision::Infinite().  The rcad curve
        // evaluation is total (no raise), so the catch arm has no counterpart
        // and only the infinite-coordinate guard below is kept.
        let p_c3d = c3d.value(u);
        let p_cons = cons.value(u);
        if is_infinite_value(p_cons.x) || is_infinite_value(p_cons.y) || is_infinite_value(p_cons.z) {
            d2 = INFINITE_VALUE;
            break;
        }
        d2 = d2.max(p_c3d.distance_squared(p_cons));
    }

    // OCCT cxx L184-187.
    let a_mult = 1.0 + 0.05;
    let mut a_deviation = a_mult * d2.sqrt();
    a_deviation = a_deviation.max(CONFUSION); // Tolerance in modeling space.
    a_deviation
}

// =========================================================================
// Check (cxx L192-266)
// =========================================================================

/// OCCT static `Check` (cxx L192-266) — the quality test of the interpolated
/// parameter-change function over `2 * nbp` samples plus the pole-order test.
fn check(
    flat_knots: &[f64],
    poles: &[f64],
    nbp: i32,
    pc3d: &[f64],
    c3d: &dyn Adaptor3dCurve,
    cons: &CurveOnSurface,
    tol: &mut f64,
    oldtol: f64,
) -> bool {
    // OCCT cxx L201-203.
    let a_degree: usize = 3;
    let mut extrap_mode = [3i32, 3i32];

    // Correction of the interval of valid values (OCCT cxx L206-221): the
    // bug-fix of OCC5898.
    let mut a_param_first = 3.0 * pc3d[0] - 2.0 * pc3d[(nbp - 1) as usize];
    let mut a_param_last = 3.0 * pc3d[(nbp - 1) as usize] - 2.0 * pc3d[0];

    let first_par = cons.first_parameter();
    let last_par = cons.last_parameter();
    if a_param_first < first_par {
        a_param_first = first_par;
    }
    if a_param_last > last_par {
        a_param_last = last_par;
    }

    // OCCT cxx L223-250.
    let mut d2 = 0.0f64; // Maximum square deviation on the samples.
    let d = *tol;
    let nn = 2 * nbp;
    let unsurnn = 1.0 / nn as f64;
    let mut tprev = a_param_first;
    for i in 0..=nn {
        // Compute corresponding parameter on 2d curve.
        // It should be inside of 3d curve parameter space.
        let t = unsurnn * i as f64;
        let tc3d = pc3d[0] * (1.0 - t) + pc3d[(nbp - 1) as usize] * t; // weight function.
        let p_c3d = c3d.value(tc3d);
        // BSplCLib::Eval(tc3d, false, 0, extrap_mode[0], aDegree, FlatKnots,
        //                1, (double&)Poles(1), tcons);
        let mut tcons_result = [0.0f64; 1];
        bspl_lib::eval_flat(
            tc3d,
            false,
            0,
            &mut extrap_mode,
            a_degree,
            flat_knots,
            1,
            poles,
            &mut tcons_result,
        );
        let tcons = tcons_result[0];

        if tcons < tprev || tcons > a_param_last {
            *tol = INFINITE_VALUE;
            return false;
        }
        tprev = tcons;
        let p_cons = cons.value(tcons);
        let temp = p_c3d.distance_squared(p_cons);
        if temp > d2 {
            d2 = temp;
        }
    }
    *tol = d2.sqrt();

    // Check poles parameters to be ordered (OCCT cxx L253-263).
    for i in 2..=poles.len() {
        // OCCT cxx L254-256: for (i = Poles.Lower() + 1; i <= Poles.Upper(); ++i).
        let a_previous_param = poles[i - 2];
        let a_current_param = poles[i - 1];

        if a_previous_param > a_current_param {
            return false;
        }
    }

    // OCCT cxx L265.
    *tol <= d || *tol > 0.8 * oldtol
}

// =========================================================================
// Geom2dAdaptor::MakeCurve (Geom2dAdaptor.cxx L33-117)
// =========================================================================

/// OCCT `Geom2dAdaptor::MakeCurve(const Adaptor2d_Curve2d& HC)`
/// (Geom2dAdaptor.cxx L33-117) — the `switch (HC.GetType())` plus the trim
/// step.  Interface difference: `geomalgo::geom_api` holds a private variant
/// over a kernel `Curve2d`; this one takes the adaptor handle exactly as the
/// OCCT signature does.
///
/// The `GeomAbs_OffsetCurve` arm (Geom2dAdaptor.cxx L78-88) is unreachable in
/// the rcad encoding: `base::proj_lib::CurveType` carries no offset value, so
/// an OCCT offset curve reports `Other` and reaches the `OtherCurve` raise.
fn geom2d_adaptor_make_curve(hc: &dyn Adaptor2dCurve2d) -> Curve2d {
    // OCCT L35-92: the type switch.
    let c2d = match hc.get_type() {
        CurveType::Line => Curve2d::Line(hc.line()),
        CurveType::Circle => Curve2d::Circle(hc.circle()),
        CurveType::Ellipse => Curve2d::Ellipse(hc.ellipse()),
        CurveType::Parabola => Curve2d::Parabola(hc.parabola()),
        CurveType::Hyperbola => Curve2d::Hyperbola(hc.hyperbola()),
        CurveType::Bezier => Curve2d::Bezier(
            hc.bezier()
                .expect("Geom2dAdaptor::MakeCurve: null Bezier handle"),
        ),
        CurveType::BSpline => Curve2d::BSpline(
            hc.bspline()
                .expect("Geom2dAdaptor::MakeCurve: null BSpline handle"),
        ),
        CurveType::Other => panic!("Standard_DomainError: Geom2dAdaptor::MakeCurve, OtherCurve"),
    };

    // Trim the curve if necessary (OCCT L94-113).
    let c_dom = Curve2dEval::default_domain(&c2d);
    if hc.first_parameter() != c_dom[0] || hc.last_parameter() != c_dom[1] {
        // if (C2D->IsPeriodic() || (HC.FirstParameter() >= C2D->FirstParameter()
        //     && HC.LastParameter() <= C2D->LastParameter()))
        if c2d.is_periodic()
            || (hc.first_parameter() >= c_dom[0] && hc.last_parameter() <= c_dom[1])
        {
            return Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(c2d),
                t_min: hc.first_parameter(),
                t_max: hc.last_parameter(),
            });
        }
        // else { tf = max(HC.FirstParameter(), C2D->FirstParameter());
        //        tl = min(HC.LastParameter(), C2D->LastParameter()); }
        let tf = hc.first_parameter().max(c_dom[0]);
        let tl = hc.last_parameter().min(c_dom[1]);
        return Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(c2d),
            t_min: tf,
            t_max: tl,
        });
    }
    c2d
}

// =========================================================================
// Adaptor3d_CurveOnSurface value copy / ShallowCopy
// =========================================================================

/// The rcad encoding of the `Adaptor3d_CurveOnSurface` value copy
/// (`ACS = aData.myCOnS`, cxx L449) and of `ShallowCopy()`
/// (Adaptor3d_CurveOnSurface.cxx L901-930): OCCT copies the two handles plus
/// the cached EvalKPart payload; the rcad `CurveOnSurface` is rebuilt from the
/// two `Arc` handles (the shallow copies) and its payload fields are carried
/// over.  The interval cache is a pure memo of `my2dCurve`/`mySurface`, so
/// leaving it unset recomputes the identical sequence.
fn curve_on_surface_copy(cos: &CurveOnSurface) -> CurveOnSurface {
    let mut copy = CurveOnSurface::new(cos.my2d_curve.clone(), cos.my_surface.clone());
    copy.my_type = cos.my_type;
    copy.my_circ = cos.my_circ;
    copy.my_lin = cos.my_lin;
    copy
}

// =========================================================================
// Approx_SameParameter_Data (hxx L84-119)
// =========================================================================

/// OCCT `Approx_SameParameter::Approx_SameParameter_Data` (hxx L84-119) — the
/// working arrays and parameter ranges.
///
/// Architecture difference: the four `double*` members point into the fixed
/// `myMaxArraySize`-sized stack buffers of `Build` (cxx L326-327); the rcad
/// value owns them as `myMaxArraySize`-long vectors with the same 0-based
/// indexing, so `myPC2d.length == myNewPC2d.length` holds as in OCCT.
struct ApproxSameParameterData {
    /// hxx L86: Adaptor3d_CurveOnSurface myCOnS.
    my_c_on_s: CurveOnSurface,
    /// hxx L87: int myNbPnt.
    my_nb_pnt: i32,
    /// hxx L88: double* myPC3d.
    my_pc3d: Vec<f64>,
    /// hxx L89: double* myPC2d.
    my_pc2d: Vec<f64>,
    /// hxx L92: double* myNewPC3d.
    my_new_pc3d: Vec<f64>,
    /// hxx L93: double* myNewPC2d.
    my_new_pc2d: Vec<f64>,
    /// hxx L96: double myC3dPF.
    my_c3d_pf: f64,
    /// hxx L97: double myC3dPL.
    my_c3d_pl: f64,
    /// hxx L98: double myC2dPF.
    my_c2d_pf: f64,
    /// hxx L99: double myC2dPL.
    my_c2d_pl: f64,
    /// hxx L101: double myTol.
    my_tol: f64,
}

impl ApproxSameParameterData {
    /// OCCT hxx L104-118: Swap(theNewNbPoints).
    fn swap(&mut self, the_new_nb_points: i32) {
        self.my_nb_pnt = the_new_nb_points;

        // 3-D.
        std::mem::swap(&mut self.my_pc3d, &mut self.my_new_pc3d);
        // 2-D.
        std::mem::swap(&mut self.my_pc2d, &mut self.my_new_pc2d);
    }
}

// =========================================================================
// Approx_SameParameter (hxx L30-173)
// =========================================================================

/// OCCT Approx_SameParameter (hxx L30-173).
pub struct ApproxSameParameter {
    /// hxx L163: const double myDeltaMin.
    my_delta_min: f64,
    /// hxx L165: bool mySameParameter.
    my_same_parameter: bool,
    /// hxx L166: bool myDone.
    my_done: bool,
    /// hxx L167: double myTolReached.
    my_tol_reached: f64,
    /// hxx L168: handle(Geom2d_Curve) myCurve2d.
    my_curve2d: Option<Curve2d>,
    /// hxx L169: handle(Adaptor2d_Curve2d) myHCurve2d.
    my_h_curve2d: Curve2dHandle,
    /// hxx L170: handle(Adaptor3d_Curve) myC3d — the rcad GeomAdaptor_Curve
    /// encoding (kernel curve + restricted window).
    my_c3d: GeomCurveAdaptor,
    /// hxx L171: handle(Adaptor3d_Surface) mySurf.
    my_surf: Arc<GeomSurfaceAdaptor>,
    /// hxx L172: handle(Adaptor3d_CurveOnSurface) myCurveOnSurface.
    my_curve_on_surface: Option<CurveOnSurface>,
}

impl ApproxSameParameter {
    /// OCCT Approx_SameParameter(C3d, Pcurv, S, Tol3d) — the
    /// `Approx_SameParameter(const handle(Adaptor3d_Curve)&,
    /// const handle(Geom2d_Curve)&, const handle(Adaptor3d_Surface)&,
    /// const double)` constructor (cxx L286-298).
    ///
    /// Interface difference: the OCCT adaptor arguments collapse into the
    /// kernel values plus the 3d curve parameter window, which is the rcad
    /// `(Curve3, first, last)` encoding of `GeomAdaptor_Curve`
    /// (`myC3d = C3D`, `myHCurve2d = new Geom2dAdaptor_Curve(C2D)`,
    /// `mySurf = S`, then `Build(Tol)`).
    ///
    /// `myTolReached` is left uninitialized by the OCCT constructor; the rcad
    /// value starts at 0.0 (the field is only read after `Build`, except on
    /// the `BuildInitialDistribution` failure path where OCCT reads an
    /// indeterminate value).
    pub fn new(
        c3d: &Curve3,
        c3d_first: f64,
        c3d_last: f64,
        pcurv: &Curve2d,
        s: &Surface3,
        tol3d: f64,
    ) -> Self {
        let mut this = ApproxSameParameter {
            // OCCT cxx L274: myDeltaMin(Precision::PConfusion()).
            my_delta_min: p_confusion(),
            // OCCT cxx L275: mySameParameter(true).
            my_same_parameter: true,
            // OCCT cxx L276: myDone(false).
            my_done: false,
            my_tol_reached: 0.0,
            my_curve2d: None,
            my_h_curve2d: Arc::new(Geom2dCurveAdaptor::new(pcurv.clone())),
            my_c3d: GeomCurveAdaptor::with_range(c3d.clone(), c3d_first, c3d_last),
            my_surf: Arc::new(GeomSurfaceAdaptor::new(s.clone())),
            my_curve_on_surface: None,
        };
        // OCCT cxx L281: Build(Tol).
        this.build(tol3d);
        this
    }

    /// OCCT Approx_SameParameter::Build (cxx L318-543).
    fn build(&mut self, tolerance: f64) {
        // OCCT cxx L326-327: the four myMaxArraySize-sized working buffers.
        let qpcons = vec![0.0f64; MY_MAX_ARRAY_SIZE];
        let qnewpcons = vec![0.0f64; MY_MAX_ARRAY_SIZE];
        let qpc3d = vec![0.0f64; MY_MAX_ARRAY_SIZE];
        let qnewpc3d = vec![0.0f64; MY_MAX_ARRAY_SIZE];

        // Create and fill data structure (OCCT cxx L330-341).
        let mut a_data = ApproxSameParameterData {
            my_c_on_s: CurveOnSurface::new(self.my_h_curve2d.clone(), self.my_surf.clone()),
            my_c2d_pf: 0.0,
            my_c2d_pl: 0.0,
            my_c3d_pf: 0.0,
            my_c3d_pl: 0.0,
            my_nb_pnt: 0, // No points initially.
            my_pc2d: qpcons,
            my_pc3d: qpc3d,
            my_new_pc2d: qnewpcons,
            my_new_pc3d: qnewpc3d,
            my_tol: tolerance,
        };
        a_data.my_c2d_pf = self.my_h_curve2d.first_parameter();
        a_data.my_c2d_pl = self.my_h_curve2d.last_parameter();
        a_data.my_c3d_pf = self.my_c3d.first_parameter();
        a_data.my_c3d_pl = self.my_c3d.last_parameter();

        // Build initial distribution (OCCT cxx L344-349).
        if !self.build_initial_distribution(&mut a_data) {
            self.my_same_parameter = false;
            self.my_done = false;
            return;
        }

        // Check same parameter state on distribution (OCCT cxx L352-374).
        let mut a_max_sq_deviation = 0.0f64;
        let a_percent_of_bad_proj = 0.3;
        let a_nb_pnt = a_data.my_nb_pnt - (a_percent_of_bad_proj * a_data.my_nb_pnt as f64) as i32;
        self.my_same_parameter = self.check_same_parameter(&mut a_data, &mut a_max_sq_deviation);
        if self.my_same_parameter {
            // OCCT cxx L358-360.
            self.my_tol_reached = compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);
            self.my_done = true;
            return;
        } else {
            // Control number of sample points after checking sameparameter
            // (OCCT cxx L364-373).
            if a_data.my_nb_pnt < a_nb_pnt {
                self.my_tol_reached =
                    compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);
                self.my_curve2d = Some(geom2d_adaptor_make_curve(&*self.my_h_curve2d));
                self.my_done = false;
                return;
            }
        }

        // Control tangents at the extremities (OCCT cxx L376-387).
        let mut tangent = [0.0f64, 0.0f64];
        let (mut a_first_tangent, mut a_last_tangent) = (tangent[0], tangent[1]);
        if !self.compute_tangents(&a_data.my_c_on_s, &mut a_first_tangent, &mut a_last_tangent) {
            // Cannot compute tangents (OCCT cxx L381-387).
            self.my_same_parameter = false;
            self.my_done = false;
            self.my_tol_reached = compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);
            return;
        }
        tangent = [a_first_tangent, a_last_tangent];

        // OCCT cxx L392-396.
        let mut a_continuity = self.my_h_curve2d.continuity();
        if (a_continuity as u8) > (GeomAbsShape::C1 as u8) {
            a_continuity = GeomAbsShape::C1;
        }

        // OCCT cxx L398-480: the loop over the number of poles.
        let mut besttol2 = a_data.my_tol * a_data.my_tol;
        let mut tolsov = INFINITE_VALUE;
        let mut interpolok = false;
        let mut has_count_changed = false;
        loop {
            // Interpolation data (OCCT cxx L402-406).
            let num_knots = a_data.my_nb_pnt as usize + 7;
            let num_poles = a_data.my_nb_pnt as usize + 3;
            let mut poles = vec![0.0f64; num_poles];
            let mut flat_knots = vec![0.0f64; num_knots];

            if !self.interpolate(&a_data, tangent[0], tangent[1], &mut poles, &mut flat_knots) {
                // Interpolation fails (OCCT cxx L409-415).
                self.my_same_parameter = false;
                self.my_done = false;
                self.my_tol_reached =
                    compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);
                return;
            }

            // OCCT cxx L417-420.
            let mut algtol = besttol2.sqrt();
            interpolok = check(
                &flat_knots,
                &poles,
                a_data.my_nb_pnt + 1,
                &a_data.my_pc3d,
                &self.my_c3d,
                &a_data.my_c_on_s,
                &mut algtol,
                tolsov,
            );
            tolsov = algtol;

            // Try to build 2d curve and check it for validity (OCCT L422-474).
            if interpolok {
                let besttol = besttol2.sqrt();

                // OCCT cxx L427-430: the 1D tolerance pair from the surface
                // resolutions at besttol (tol2d / tol3d stay null handles).
                let tol1d = vec![
                    self.my_surf.u_resolution(besttol),
                    self.my_surf.v_resolution(besttol),
                ];

                // OCCT cxx L432-445.
                let mut ev = ApproxSameParameterEvaluator::new(&flat_knots, &poles, self.my_h_curve2d.clone());
                let a_max_deg = 11;
                let a_max_seg = 1000;
                let an_approximator = ApproxAFunction::new(
                    2,
                    0,
                    0,
                    Some(&tol1d),
                    None,
                    None,
                    a_data.my_c3d_pf,
                    a_data.my_c3d_pl,
                    a_continuity,
                    a_max_deg,
                    a_max_seg,
                    &mut ev,
                );

                if an_approximator.is_done() || an_approximator.has_result() {
                    // OCCT cxx L449-453.
                    let acs = curve_on_surface_copy(&a_data.my_c_on_s);
                    let a_curve_builder = GeomLibMakeCurvefromApprox::new(&an_approximator);
                    let a_c2d = a_curve_builder
                        .curve2d_from_two1d(1, 2)
                        .expect("GeomLib_MakeCurvefromApprox::Curve2dFromTwo1d");
                    // OCCT cxx L452: aHCurve2d = new Geom2dAdaptor_Curve(aC2d)
                    // — the same handle the result keeps; the rcad adaptor
                    // owns a copy of the kernel Curve2d.
                    let a_c2d_curve = Curve2d::BSpline(a_c2d);
                    let a_h_curve2d: Curve2dHandle =
                        Arc::new(Geom2dCurveAdaptor::new(a_c2d_curve.clone()));
                    // OCCT cxx L453: aData.myCOnS.Load(aHCurve2d) — the rcad
                    // CurveOnSurface sets the 2d handle then runs Load(C)
                    // (Adaptor3d_CurveOnSurface.cxx L943-964).
                    a_data.my_c_on_s.my2d_curve = a_h_curve2d.clone();
                    a_data.my_c_on_s.load_curve();
                    self.my_tol_reached =
                        compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);

                    // OCCT cxx L456-472.
                    let a_mult = 250.0; // To be tolerant with discrete tolerance.
                    if self.my_tol_reached < a_mult * besttol {
                        self.my_curve2d = Some(a_c2d_curve);
                        self.my_h_curve2d = a_h_curve2d;
                        self.my_done = true;
                        break;
                    } else if (a_data.my_nb_pnt as usize) < MY_MAX_ARRAY_SIZE - 1 {
                        interpolok = false;
                        a_data.my_c_on_s = acs;
                    } else {
                        break;
                    }
                }
            }

            // OCCT cxx L476-479.
            if !interpolok {
                has_count_changed =
                    self.increase_nb_poles(&poles, &flat_knots, &mut a_data, &mut besttol2);
            }

            // OCCT cxx L400: while (!interpolok && hasCountChanged).
            if !(!interpolok && has_count_changed) {
                break;
            }
        }

        if !self.my_done {
            // OCCT cxx L482-540: the loop finished unsuccessfully.

            // Original 2d curve (OCCT cxx L489-491).
            a_data.my_c_on_s.my2d_curve = self.my_h_curve2d.clone();
            a_data.my_c_on_s.load_curve();
            self.my_tol_reached = compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);
            self.my_curve2d = Some(geom2d_adaptor_make_curve(&*self.my_h_curve2d));

            // Approximation curve (OCCT cxx L494-507).
            let num_knots = a_data.my_nb_pnt as usize + 7;
            let num_poles = a_data.my_nb_pnt as usize + 3;
            let mut poles = vec![0.0f64; num_poles];
            let mut flat_knots = vec![0.0f64; num_knots];

            self.interpolate(&a_data, tangent[0], tangent[1], &mut poles, &mut flat_knots);

            let besttol = besttol2.sqrt();
            let tol1d = vec![
                self.my_surf.u_resolution(besttol),
                self.my_surf.v_resolution(besttol),
            ];

            // OCCT cxx L507-519.
            let mut ev = ApproxSameParameterEvaluator::new(&flat_knots, &poles, self.my_h_curve2d.clone());
            let an_approximator = ApproxAFunction::new(
                2,
                0,
                0,
                Some(&tol1d),
                None,
                None,
                a_data.my_c3d_pf,
                a_data.my_c3d_pl,
                a_continuity,
                11,
                40,
                &mut ev,
            );

            // OCCT cxx L521-525.
            if !an_approximator.is_done() && !an_approximator.has_result() {
                self.my_done = false;
                return;
            }

            // OCCT cxx L527-530.
            let a_curve_builder = GeomLibMakeCurvefromApprox::new(&an_approximator);
            let a_c2d = a_curve_builder
                .curve2d_from_two1d(1, 2)
                .expect("GeomLib_MakeCurvefromApprox::Curve2dFromTwo1d");
            // OCCT cxx L529: aHCurve2d = new Geom2dAdaptor_Curve(aC2d).
            let a_c2d_curve = Curve2d::BSpline(a_c2d);
            let a_h_curve2d: Curve2dHandle =
                Arc::new(Geom2dCurveAdaptor::new(a_c2d_curve.clone()));
            a_data.my_c_on_s.my2d_curve = a_h_curve2d.clone();
            a_data.my_c_on_s.load_curve();

            // OCCT cxx L532-539.
            let an_approx_tol =
                compute_tol_reached(&self.my_c3d, &a_data.my_c_on_s, 2 * MY_NB_SAMPLES);
            if an_approx_tol < self.my_tol_reached {
                self.my_tol_reached = an_approx_tol;
                self.my_curve2d = Some(a_c2d_curve);
                self.my_h_curve2d = a_h_curve2d;
            }
            self.my_done = true;
        }

        // OCCT cxx L542: myCurveOnSurface =
        // occ::down_cast<Adaptor3d_CurveOnSurface>(aData.myCOnS.ShallowCopy()).
        self.my_curve_on_surface = Some(curve_on_surface_copy(&a_data.my_c_on_s));
    }

    /// OCCT Approx_SameParameter::BuildInitialDistribution (cxx L547-580).
    fn build_initial_distribution(&self, the_data: &mut ApproxSameParameterData) -> bool {
        // Take a multiple of the sample of CheckShape, at least the control
        // points will be correct (OCCT cxx L549-566).
        let deltacons = (the_data.my_c2d_pl - the_data.my_c2d_pf) / MY_NB_SAMPLES as f64;
        let deltac3d = (the_data.my_c3d_pl - the_data.my_c3d_pf) / MY_NB_SAMPLES as f64;
        let mut wcons = the_data.my_c2d_pf;
        let mut wc3d = the_data.my_c3d_pf;
        for ii in 0..MY_NB_SAMPLES as usize {
            the_data.my_pc2d[ii] = wcons;
            the_data.my_pc3d[ii] = wc3d;
            wcons += deltacons;
            wc3d += deltac3d;
        }
        the_data.my_nb_pnt = MY_NB_SAMPLES;
        the_data.my_pc2d[the_data.my_nb_pnt as usize] = the_data.my_c2d_pl;
        the_data.my_pc3d[the_data.my_nb_pnt as usize] = the_data.my_c3d_pl;

        // Change number of points in case of C0 continuity (OCCT cxx
        // L568-577).
        let continuity = self.my_h_curve2d.continuity();
        if (continuity as u8) < (GeomAbsShape::C1 as u8) {
            if !self.increase_initial_nb_samples(the_data) {
                // Number of samples is too big.
                return false;
            }
        }

        true
    }

    /// OCCT Approx_SameParameter::IncreaseInitialNbSamples (cxx L588-649) —
    /// get the number of C1 intervals and build a new distribution on them.
    fn increase_initial_nb_samples(&self, the_data: &mut ApproxSameParameterData) -> bool {
        // OCCT cxx L590-592: NbInt = NbIntervals(GeomAbs_C1) + 1; the array is
        // (1, NbInt) and keeps that length while NbInt is decremented below.
        let a_c1_intervals = self.my_h_curve2d.intervals(GeomAbsShape::C1);
        let mut nb_int = a_c1_intervals.len() as i32;

        // OCCT cxx L594-602.
        let mut inter = 1i32;
        while inter <= nb_int && a_c1_intervals[(inter - 1) as usize] <= the_data.my_c3d_pf + self.my_delta_min
        {
            inter += 1;
        }
        while nb_int > 0
            && a_c1_intervals[(nb_int - 1) as usize] >= the_data.my_c3d_pl - self.my_delta_min
        {
            nb_int -= 1;
        }

        // Compute new parameters (OCCT cxx L604-631).
        let mut a_new_par: Vec<f64> = Vec::new();
        a_new_par.push(the_data.my_c3d_pf);
        let mut ii = 1i32;
        while inter <= nb_int
            || (ii < MY_NB_SAMPLES && inter <= a_c1_intervals.len() as i32)
        {
            if a_c1_intervals[(inter - 1) as usize] < the_data.my_pc2d[ii as usize] {
                a_new_par.push(a_c1_intervals[(inter - 1) as usize]);
                if (the_data.my_pc2d[ii as usize] - a_c1_intervals[(inter - 1) as usize])
                    <= self.my_delta_min
                {
                    ii += 1;
                    if ii > MY_NB_SAMPLES {
                        ii = MY_NB_SAMPLES;
                    }
                }
                inter += 1;
            } else {
                if (a_c1_intervals[(inter - 1) as usize] - the_data.my_pc2d[ii as usize])
                    > self.my_delta_min
                {
                    a_new_par.push(the_data.my_pc2d[ii as usize]);
                }
                ii += 1;
            }
        }
        // Simple protection if theNewNbPoints > allocated elements in array but
        // one myMaxArraySize - 1 index may be filled after projection (OCCT
        // cxx L632-638).
        the_data.my_nb_pnt = a_new_par.len() as i32;
        if the_data.my_nb_pnt as usize > MY_MAX_ARRAY_SIZE - 1 {
            return false;
        }

        // OCCT cxx L640-646.
        let mut ii = 1i32;
        while ii < the_data.my_nb_pnt {
            // Copy only internal points.
            let value = a_new_par[ii as usize];
            the_data.my_pc2d[ii as usize] = value;
            the_data.my_pc3d[ii as usize] = value;
            ii += 1;
        }
        the_data.my_pc3d[the_data.my_nb_pnt as usize] = the_data.my_c3d_pl;
        the_data.my_pc2d[the_data.my_nb_pnt as usize] = the_data.my_c2d_pl;

        true
    }

    /// OCCT Approx_SameParameter::CheckSameParameter (cxx L653-766).
    fn check_same_parameter(
        &self,
        the_data: &mut ApproxSameParameterData,
        the_sq_dist: &mut f64,
    ) -> bool {
        // OCCT cxx L656-657.
        let tol2 = the_data.my_tol * the_data.my_tol;
        let mut is_same_param = true;

        // Compute initial distance on boundary points (OCCT cxx L659-669).
        let pcons = the_data.my_c_on_s.value(the_data.my_c2d_pf);
        let pc3d = self.my_c3d.value(the_data.my_c3d_pf);
        let dist2 = pcons.distance_squared(pc3d);
        let mut dmax2 = dist2;

        let pcons = the_data.my_c_on_s.value(the_data.my_c2d_pl);
        let pc3d = self.my_c3d.value(the_data.my_c3d_pl);
        let dist2 = pcons.distance_squared(pc3d);
        dmax2 = dmax2.max(dist2);

        // OCCT cxx L671-672: Extrema_LocateExtPC Projector;
        // Projector.Initialize(*myC3d, myC3dPF, myC3dPL, myTol).
        let tool = CurveToolHandle::for_curve3(&self.my_c3d.curve, &self.my_c3d, &self.my_c3d);
        let mut projector = LocateExtPC::new();
        projector.initialize(
            &tool,
            the_data.my_c3d_pf,
            the_data.my_c3d_pl,
            the_data.my_tol,
        );

        // OCCT cxx L674-678.
        let mut count = 1i32;
        let mut previousp = the_data.my_c3d_pf;
        let mut initp = 0.0f64;
        // OCCT cxx L675: double ... curp; — always written before it is read.
        let mut curp = 0.0f64;
        let bornesup = the_data.my_c3d_pl - self.my_delta_min;
        let mut is_proj_ok = false;
        for ii in 1..the_data.my_nb_pnt {
            let ii = ii as usize;
            // OCCT cxx L680-682.
            let pcons = the_data.my_c_on_s.value(the_data.my_pc2d[ii]);
            let pc3d = self.my_c3d.value(the_data.my_pc3d[ii]);
            let dist2 = pcons.distance_squared(pc3d);

            // Same parameter point (OCCT cxx L684-698).
            let is_use_param = dist2 <= tol2
                && (the_data.my_pc3d[ii]
                    > the_data.my_pc3d[(count - 1) as usize] + self.my_delta_min);
            if is_use_param {
                if dmax2 < dist2 {
                    dmax2 = dist2;
                }
                the_data.my_pc3d[count as usize] = the_data.my_pc3d[ii];
                initp = the_data.my_pc3d[count as usize];
                previousp = initp;
                the_data.my_pc2d[count as usize] = the_data.my_pc2d[ii];
                count += 1;
                continue;
            }

            // Local search: local extrema and iterative projection algorithm
            // (OCCT cxx L700-726).
            if !is_proj_ok {
                initp = the_data.my_pc3d[ii];
            }
            is_proj_ok = false;
            is_same_param = false;
            projector.perform(pcons, initp);
            if projector.is_done() {
                // Local extrema is found.
                curp = projector.point().param;
                is_proj_ok = true;
            } else {
                project_point_on_curve(
                    initp,
                    pcons,
                    the_data.my_tol,
                    30,
                    &self.my_c3d,
                    &mut is_proj_ok,
                    &mut curp,
                );
            }
            is_proj_ok = is_proj_ok && // Good projection.
                curp > previousp + self.my_delta_min && // Point is separated from previous.
                curp < bornesup; // Inside of parameter space.
            if is_proj_ok {
                the_data.my_pc3d[count as usize] = curp;
                initp = curp;
                previousp = curp;
                the_data.my_pc2d[count as usize] = the_data.my_pc2d[ii];
                count += 1;
                continue;
            }

            // Whole parameter space search using general extrema (OCCT cxx
            // L728-758).
            let pr = ExtremaExtPC::new_point_curve_ranged(
                pcons,
                &tool,
                the_data.my_c3d_pf,
                the_data.my_c3d_pl,
                the_data.my_tol,
            );
            if !pr.is_done() || pr.nb_ext() == 0 {
                // Lazy evaluation is used.
                continue;
            }

            let a_nb_ext = pr.nb_ext();
            let mut an_ind_min = 0usize;
            let mut a_cur_dist_min = REAL_LAST;
            for i in 1..=a_nb_ext {
                let a_p = pr.point(i).point;
                let a_dist2 = a_p.distance_squared(pcons);
                if a_dist2 < a_cur_dist_min {
                    a_cur_dist_min = a_dist2;
                    an_ind_min = i;
                }
            }
            if an_ind_min != 0 {
                curp = pr.point(an_ind_min).param;
                if curp > previousp + self.my_delta_min && curp < bornesup {
                    the_data.my_pc3d[count as usize] = curp;
                    initp = curp;
                    previousp = curp;
                    the_data.my_pc2d[count as usize] = the_data.my_pc2d[ii];
                    count += 1;
                    is_proj_ok = true;
                }
            }
        }

        // OCCT cxx L760-765.
        the_data.my_nb_pnt = count;
        the_data.my_pc2d[the_data.my_nb_pnt as usize] = the_data.my_c2d_pl;
        the_data.my_pc3d[the_data.my_nb_pnt as usize] = the_data.my_c3d_pl;

        *the_sq_dist = dmax2;
        is_same_param
    }

    /// OCCT Approx_SameParameter::ComputeTangents (cxx L770-809).
    fn compute_tangents(
        &self,
        the_c_on_s: &CurveOnSurface,
        the_first_tangent: &mut f64,
        the_last_tangent: &mut f64,
    ) -> bool {
        // OCCT cxx L774-777.
        let a_small_magnitude = 1.0e-12;

        // First point (OCCT cxx L779-791).
        let a_param_first = self.my_c3d.first_parameter();
        let (_a_pnt_c_on_s, a_vec_con_s) = the_c_on_s.d1(a_param_first);
        let (_a_pnt, a_vec) = self.my_c3d.d1(a_param_first);
        let a_magnitude = a_vec_con_s.length();
        if a_magnitude > a_small_magnitude {
            *the_first_tangent = a_vec.length() / a_magnitude;
        } else {
            return false;
        }

        // Last point (OCCT cxx L793-806).
        let a_param_last = self.my_c3d.last_parameter();
        let (_a_pnt_c_on_s, a_vec_con_s) = the_c_on_s.d1(a_param_last);
        let (_a_pnt, a_vec) = self.my_c3d.d1(a_param_last);

        let a_magnitude = a_vec_con_s.length();
        if a_magnitude > a_small_magnitude {
            *the_last_tangent = a_vec.length() / a_magnitude;
        } else {
            return false;
        }

        true
    }

    /// OCCT Approx_SameParameter::Interpolate (cxx L813-853).
    fn interpolate(
        &self,
        the_data: &ApproxSameParameterData,
        a_tang_first: f64,
        a_tang_last: f64,
        the_poles: &mut [f64],
        the_flat_knots: &mut [f64],
    ) -> bool {
        // OCCT cxx L819-821.
        let num_poles = the_data.my_nb_pnt as usize + 3;
        let mut contact_order = vec![0i32; num_poles];
        let mut a_parameters = vec![0.0f64; num_poles];

        // Fill tables taking attention to end values (OCCT cxx L823-828).
        contact_order[1] = 1;
        contact_order[num_poles - 2] = 1;

        the_flat_knots[0] = the_data.my_c3d_pf;
        the_flat_knots[1] = the_data.my_c3d_pf;
        the_flat_knots[2] = the_data.my_c3d_pf;
        the_flat_knots[3] = the_data.my_c3d_pf;
        the_flat_knots[num_poles] = the_data.my_c3d_pl;
        the_flat_knots[num_poles + 1] = the_data.my_c3d_pl;
        the_flat_knots[num_poles + 2] = the_data.my_c3d_pl;
        the_flat_knots[num_poles + 3] = the_data.my_c3d_pl;

        the_poles[0] = the_data.my_c2d_pf;
        the_poles[num_poles - 1] = the_data.my_c2d_pl;
        the_poles[1] = a_tang_first;
        the_poles[num_poles - 2] = a_tang_last;

        a_parameters[0] = the_data.my_c3d_pf;
        a_parameters[1] = the_data.my_c3d_pf;
        a_parameters[num_poles - 2] = the_data.my_c3d_pl;
        a_parameters[num_poles - 1] = the_data.my_c3d_pl;

        for ii in 3..=num_poles - 2 {
            // OCCT cxx L841-842: thePoles(ii) = myPC2d[ii - 2];
            // aParameters(ii) = theFlatKnots(ii + 2) = myPC3d[ii - 2].
            the_poles[ii - 1] = the_data.my_pc2d[ii - 2];
            let value = the_data.my_pc3d[ii - 2];
            a_parameters[ii - 1] = value;
            the_flat_knots[ii + 1] = value;
        }
        // OCCT cxx L844-852: BSplCLib::Interpolate(3, theFlatKnots, aParameters,
        // ContactOrder, 1, thePoles(1), inversion_problem); the rcad engine
        // returns the inversion problem on success and the OCCT error codes
        // (which the OCCT body turns into raises) otherwise — both are
        // non-zero on failure.
        let inversion_problem = bspl_lib::interpolate(
            3,
            the_flat_knots,
            &a_parameters,
            &contact_order,
            1,
            the_poles,
        );
        inversion_problem == 0
    }

    /// OCCT Approx_SameParameter::IncreaseNbPoles (cxx L857-992).
    fn increase_nb_poles(
        &self,
        the_poles: &[f64],
        the_flat_knots: &[f64],
        the_data: &mut ApproxSameParameterData,
        the_best_sq_tol: &mut f64,
    ) -> bool {
        // OCCT cxx L862-865.
        let tool = CurveToolHandle::for_curve3(&self.my_c3d.curve, &self.my_c3d, &self.my_c3d);
        let mut projector = LocateExtPC::new();
        projector.initialize(
            &tool,
            self.my_c3d.first_parameter(),
            self.my_c3d.last_parameter(),
            the_data.my_tol,
        );
        let mut curp = 0.0f64;
        let mut projok = false;

        // Project middle point to fix parameterization and check projection
        // existence (OCCT cxx L867-873).
        let a_degree: usize = 3;
        let derivative_request = 0;
        let mut extrap_mode = [3i32, 3i32];
        let mut newcount = 0i32;
        for ii in 0..the_data.my_nb_pnt {
            let ii = ii as usize;
            // OCCT cxx L876-883.
            the_data.my_new_pc2d[newcount as usize] = the_data.my_pc2d[ii];
            the_data.my_new_pc3d[newcount as usize] = the_data.my_pc3d[ii];
            newcount += 1;

            if the_data.my_nb_pnt as usize - ii + newcount as usize == MY_MAX_ARRAY_SIZE {
                continue;
            }

            // OCCT cxx L885-893.
            let mut eval_result = [0.0f64; 1];
            bspl_lib::eval_flat(
                0.5 * (the_data.my_pc3d[ii] + the_data.my_pc3d[ii + 1]),
                false,
                derivative_request,
                &mut extrap_mode,
                a_degree,
                the_flat_knots,
                1,
                the_poles,
                &mut eval_result,
            );

            // OCCT cxx L895-926.
            if eval_result[0] < the_data.my_pc2d[ii]
                || eval_result[0] > the_data.my_pc2d[ii + 1]
            {
                let ucons = 0.5 * (the_data.my_pc2d[ii] + the_data.my_pc2d[ii + 1]);
                let uc3d = 0.5 * (the_data.my_pc3d[ii] + the_data.my_pc3d[ii + 1]);

                let pcons = the_data.my_c_on_s.value(ucons);
                projector.perform(pcons, uc3d);
                if projector.is_done() {
                    curp = projector.point().param;
                    let dist_2 = projector.square_distance();
                    if dist_2 > *the_best_sq_tol {
                        *the_best_sq_tol = dist_2;
                    }
                    projok = true;
                } else {
                    project_point_on_curve(
                        uc3d,
                        pcons,
                        the_data.my_tol,
                        30,
                        &self.my_c3d,
                        &mut projok,
                        &mut curp,
                    );
                }
                if projok {
                    if curp > the_data.my_pc3d[ii] + self.my_delta_min
                        && curp < the_data.my_pc3d[ii + 1] - self.my_delta_min
                    {
                        the_data.my_new_pc3d[newcount as usize] = curp;
                        the_data.my_new_pc2d[newcount as usize] = ucons;
                        newcount += 1;
                    }
                }
            }
        } // OCCT cxx L927: for (ii = 0; ii < count; ii++)

        // OCCT cxx L928-929.
        the_data.my_new_pc3d[newcount as usize] =
            the_data.my_pc3d[the_data.my_nb_pnt as usize];
        the_data.my_new_pc2d[newcount as usize] =
            the_data.my_pc2d[the_data.my_nb_pnt as usize];

        // OCCT cxx L931-936.
        if (the_data.my_nb_pnt != newcount) && (newcount as usize) < MY_MAX_ARRAY_SIZE - 1 {
            // Distribution is changed.
            the_data.swap(newcount);
            return true;
        }

        // Increase number of samples in two times (OCCT cxx L938-980).
        newcount = 0;
        for n in 0..the_data.my_nb_pnt {
            let n = n as usize;
            // OCCT cxx L942-949.
            the_data.my_new_pc3d[newcount as usize] = the_data.my_pc3d[n];
            the_data.my_new_pc2d[newcount as usize] = the_data.my_pc2d[n];
            newcount += 1;

            if the_data.my_nb_pnt as usize - n + newcount as usize == MY_MAX_ARRAY_SIZE {
                continue;
            }

            let ucons = 0.5 * (the_data.my_pc2d[n] + the_data.my_pc2d[n + 1]);
            let uc3d = 0.5 * (the_data.my_pc3d[n] + the_data.my_pc3d[n + 1]);

            // OCCT cxx L954-979.
            let pcons = the_data.my_c_on_s.value(ucons);
            projector.perform(pcons, uc3d);
            if projector.is_done() {
                curp = projector.point().param;
                let dist_2 = projector.square_distance();
                if dist_2 > *the_best_sq_tol {
                    *the_best_sq_tol = dist_2;
                }
                projok = true;
            } else {
                project_point_on_curve(
                    uc3d,
                    pcons,
                    the_data.my_tol,
                    30,
                    &self.my_c3d,
                    &mut projok,
                    &mut curp,
                );
            }
            if projok {
                if curp > the_data.my_pc3d[n] + self.my_delta_min
                    && curp < the_data.my_pc3d[n + 1] - self.my_delta_min
                {
                    the_data.my_new_pc3d[newcount as usize] = curp;
                    the_data.my_new_pc2d[newcount as usize] = ucons;
                    newcount += 1;
                }
            }
        }
        // OCCT cxx L981-982.
        the_data.my_new_pc3d[newcount as usize] =
            the_data.my_pc3d[the_data.my_nb_pnt as usize];
        the_data.my_new_pc2d[newcount as usize] =
            the_data.my_pc2d[the_data.my_nb_pnt as usize];

        // OCCT cxx L984-989.
        if the_data.my_nb_pnt != newcount {
            // Distribution is changed.
            the_data.swap(newcount);
            return true;
        }

        false
    }

    /// OCCT IsDone() (hxx L55).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsSameParameter() (hxx L64).
    pub fn is_same_parameter(&self) -> bool {
        self.my_same_parameter
    }

    /// OCCT TolReached() (hxx L59).
    pub fn tol_reached(&self) -> f64 {
        self.my_tol_reached
    }

    /// OCCT Curve2d() (hxx L69).
    pub fn curve2d(&self) -> Curve2d {
        self.my_curve2d.clone().expect("Approx_SameParameter::Curve2d")
    }

    /// OCCT Curve3d() (hxx L74) — returns `myC3d`, the adaptor of the input 3d
    /// curve.  rcad encoding: the kernel curve of the adaptor, wrapped in a
    /// `Geom_TrimmedCurve` when the adaptor carries a restricted parameter
    /// window, so that `default_domain()` answers the adaptor's
    /// `FirstParameter()`/`LastParameter()` and evaluation matches
    /// `GeomAdaptor_Curve::Value()`.
    pub fn curve3d(&self) -> Curve3 {
        let basis = self.my_c3d.curve.clone();
        let dom = CurveEval::default_domain(&basis);
        let first = self.my_c3d.first_parameter();
        let last = self.my_c3d.last_parameter();
        if dom[0] != first || dom[1] != last {
            Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(basis, first, last))
        } else {
            basis
        }
    }

    /// OCCT CurveOnSurface() (hxx L78) — the rcad (pcurve, surface) value pair
    /// of the `Adaptor3d_CurveOnSurface` handle.
    pub fn curve_on_surface(&self) -> (Curve2d, Surface3) {
        let cos = self
            .my_curve_on_surface
            .as_ref()
            .expect("Approx_SameParameter::CurveOnSurface");
        let pcurve = geom2d_adaptor_make_curve(&*cos.my2d_curve);
        let surface = cos
            .my_surface
            .kernel_surface()
            .expect("Approx_SameParameter::CurveOnSurface: no kernel surface")
            .clone();
        (pcurve, surface)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec2;
    use rcad_kernel::geom::{BezierCurve2, CurveEval, Line3, Plane, SurfaceEval};

    /// Anchor: XY plane (u = x, v = y) and the 3d line (0,0,0) -> (1,0,0) over
    /// [0, 1], with a pcurve whose parameterization already matches the 3d
    /// curve (the 2-pole Bezier through (0,0)-(1,0) evaluates to (t, 0)).
    /// Build takes the same-parameter early exit (cxx L356-361): myDone,
    /// mySameParameter and myTolReached = 1.05 * 0 floored at Confusion.
    #[test]
    fn same_parameter_early_exit() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let c3d = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let pcurv = Curve2d::Bezier(BezierCurve2 {
            control_points: vec![DVec2::ZERO, DVec2::new(1.0, 0.0)],
            weights: vec![1.0, 1.0],
        });

        let sp = ApproxSameParameter::new(&c3d, 0.0, 1.0, &pcurv, &s, CONFUSION);

        assert!(sp.is_done());
        assert!(sp.is_same_parameter());
        assert!(
            (sp.tol_reached() - CONFUSION).abs() <= 1.0e-15,
            "tol_reached={}",
            sp.tol_reached()
        );

        // Curve3d() answers the GeomAdaptor_Curve of the input: the restricted
        // window is the adaptor domain and the evaluation is unchanged.
        let c3d_out = sp.curve3d();
        assert_eq!(CurveEval::default_domain(&c3d_out), [0.0, 1.0]);
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            assert!(
                (c3d_out.point_at(t) - c3d.point_at(t)).length() <= 1.0e-12,
                "t={}",
                t
            );
        }
    }

    /// Anchor: the same plane and 3d line, with a pcurve spanning only half
    /// the 3d range (the 2-pole Bezier through (0,0)-(0.5,0) evaluates to
    /// (0.5*t, 0)).  The same-parameter test fails, the projection pass keeps
    /// the whole distribution, so Build runs the interpolation loop and the
    /// AdvApprox rebuild of the pcurve.  The rebuilt pcurve evaluated at the
    /// 3d parameter must land on the 3d curve within TolReached — the OCCT
    /// same-parameter contract.
    #[test]
    fn not_same_parameter_rebuilds_pcurve() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let c3d = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let pcurv = Curve2d::Bezier(BezierCurve2 {
            control_points: vec![DVec2::ZERO, DVec2::new(0.5, 0.0)],
            weights: vec![1.0, 1.0],
        });

        let sp = ApproxSameParameter::new(&c3d, 0.0, 1.0, &pcurv, &s, CONFUSION);

        assert!(!sp.is_same_parameter());
        assert!(sp.is_done());

        // The rebuilt pcurve must be the same-parameter one.
        let pc = sp.curve2d();
        let ntol = sp.tol_reached();
        for i in 0..=44 {
            let t = i as f64 / 44.0;
            let uv = pc.point_at(t);
            assert!(uv.x >= -1.0e-9 && uv.x <= 0.5 + 1.0e-9, "uv={:?}", uv);
            assert!(uv.y.abs() <= 1.0e-9, "uv={:?}", uv);
            let p_cons = s.point_at(uv.x, uv.y);
            let p_c3d = c3d.point_at(t);
            assert!(
                (p_cons - p_c3d).length() <= ntol,
                "t={} deviation={} tol={}",
                t,
                (p_cons - p_c3d).length(),
                ntol
            );
        }

        // CurveOnSurface() exposes the same locus.
        let (pc_on_s, s_on_s) = sp.curve_on_surface();
        let uv = pc_on_s.point_at(0.25);
        let p_on_s = s_on_s.point_at(uv.x, uv.y);
        let p_c3d = c3d.point_at(0.25);
        assert!((p_on_s - p_c3d).length() <= ntol);
    }

    /// Anchor: same plane / 3d line, with a C0 pcurve — the degree-1 BSpline
    /// through (0,0), (0.5,0), (1,0) with the interior knot 0.5 of
    /// multiplicity 1.  Geom2dAdaptor_Curve::Continuity answers C0, so
    /// BuildInitialDistribution runs IncreaseInitialNbSamples (cxx
    /// L588-649): the distribution is rebuilt on the C1 intervals of the
    /// pcurve.  The curve is still the exact (u, 0) locus, so the
    /// same-parameter contract holds at the end.
    #[test]
    fn c0_pcurve_increases_initial_samples() {
        use rcad_kernel::geom::BSplineCurve2;

        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let c3d = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        let pcurv = Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::new(0.5, 0.0), DVec2::new(1.0, 0.0)],
            weights: vec![1.0, 1.0, 1.0],
            is_periodic: false,
        });

        let sp = ApproxSameParameter::new(&c3d, 0.0, 1.0, &pcurv, &s, CONFUSION);

        assert!(sp.is_done());
        assert!(sp.is_same_parameter());
        for i in 0..=44 {
            let t = i as f64 / 44.0;
            let p_c3d = c3d.point_at(t);
            let p_cons = s.point_at(t, 0.0);
            assert!((p_cons - p_c3d).length() <= sp.tol_reached());
        }
    }
}
