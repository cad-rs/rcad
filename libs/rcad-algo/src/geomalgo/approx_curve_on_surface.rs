// OCCT Approx_CurveOnSurface.hxx L17-141 + Approx_CurveOnSurface.cxx L44-813
// (ModelingData/TKGeomBase/Approx).
//
// 1:1 translation of the approximation of a curve on surface: the two
// constructors, Perform, the isoline fast path (isIsoLine / buildC3dOnIsoLine)
// and the result accessors, plus the three AdvApprox_EvaluatorFunction
// bindings (Approx_CurveOnSurface_Eval / _Eval3d / _Eval2d, cxx L46-306).
//
// Dependency encodings:
// - handle(Adaptor2d_Curve2d) / handle(Adaptor3d_Surface) ->
//   rcad_kernel::base::proj_lib::adaptor::{Curve2dHandle, SurfaceHandle},
// - Adaptor3d_CurveOnSurface -> the CurveOnSurface encoding in the same
//   kernel module,
// - AdvApprox_ApproxAFunction / AdvApprox_DichoCutting /
//   AdvApprox_EvaluatorFunction -> rcad_kernel::math::adv_approx.
//
// GAP (staged, kernel adv_approx): AdvApprox_PrefAndRec and the
// cut-tool-parameterized ApproxAFunction constructor are not translated —
// the OCCT CutTool selection below keeps the OCCT form and the PrefAndRec
// branches panic until the kernel closes the gap; ApproxAFunction::new runs
// Perform with its internal dichotomy cutting.
// GAP (staged, kernel adv_approx): the 1D subspace results (Poles1d,
// MaxError(1, N)) are not stored by the kernel ApproxAFunction, so the 2D
// result extraction of Perform panics (the !theOnly3d arm).
// GAP (staged, kernel geom): Geom_Surface::UIso/VIso and
// Geom_RectangularTrimmedSurface are not translated, so the
// buildC3dOnIsoLine isoline extraction panics at its first GAP leaf.

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::proj_lib::adaptor::{
    is_surf_g1, Adaptor2dCurve2d, Adaptor3dCurve, Curve2dHandle, CurveOnSurface, SurfaceHandle,
};
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::core::precision::{p_confusion, ANGULAR, CONFUSION};
use rcad_kernel::geom::{Curve2d, Curve3, CurveEval, Surface3, SurfaceEval, BSplineCurve3, TrimmedCurve3};
use rcad_kernel::base::proj_lib::adaptor::GeomAbsSurfaceType;
use rcad_kernel::math::adv_approx::{ApproxAFunction, Cutting, DichoCutting, EvaluatorFunction};
use rcad_kernel::math::bspl_lib;
use rcad_kernel::math::GeomAbsShape;

// ---------------------------------------------------------------------------
// Approx_CurveOnSurface_Eval (cxx L46-142)
// ---------------------------------------------------------------------------

/// OCCT Approx_CurveOnSurface_Eval (cxx L46-71) — evaluates the 5D function
/// (u, v, x, y, z) of the curve-on-surface.
struct CurveOnSurfaceEval {
    /// OCCT: handle(Adaptor3d_Curve) fonct.
    fonct: Arc<dyn Adaptor3dCurve>,
    /// OCCT: handle(Adaptor2d_Curve2d) fonct2d.
    fonct2d: Curve2dHandle,
    /// OCCT: StartEndSav[2].
    start_end_sav: [f64; 2],
}

impl CurveOnSurfaceEval {
    /// OCCT ctor (cxx L49-58).
    fn new(
        the_func: Arc<dyn Adaptor3dCurve>,
        the_func2d: Curve2dHandle,
        first: f64,
        last: f64,
    ) -> Self {
        CurveOnSurfaceEval {
            fonct: the_func,
            fonct2d: the_func2d,
            start_end_sav: [first, last],
        }
    }
}

impl EvaluatorFunction for CurveOnSurfaceEval {
    /// OCCT Approx_CurveOnSurface_Eval::Evaluate (cxx L73-142).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        let mut error_code = 0;
        let par = parameter;

        // Dimension is incorrect (cxx L84-87).
        if result.len() != 5 {
            error_code = 1;
        }

        // Parameter is incorrect (cxx L90-96): re-trim both functions.
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct = self.fonct.trim(start_end[0], start_end[1], p_confusion());
            self.fonct2d = self.fonct2d.trim(start_end[0], start_end[1], p_confusion());
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

        match derivative_request {
            0 => {
                let pnt2d = self.fonct2d.d0(par);
                let pnt = self.fonct.value(par);
                result[0] = pnt2d.x;
                result[1] = pnt2d.y;
                result[2] = pnt.x;
                result[3] = pnt.y;
                result[4] = pnt.z;
            }
            1 => {
                let (_pnt2d, v21) = self.fonct2d.d1(par);
                let (_pnt, v1) = self.fonct.d1(par);
                result[0] = v21.x;
                result[1] = v21.y;
                result[2] = v1.x;
                result[3] = v1.y;
                result[4] = v1.z;
            }
            2 => {
                let (_pnt2d, _v21, v22) = self.fonct2d.d2(par);
                let (_pnt, _v1, v2) = self.fonct.d2(par);
                result[0] = v22.x;
                result[1] = v22.y;
                result[2] = v2.x;
                result[3] = v2.y;
                result[4] = v2.z;
            }
            _ => {
                result[0] = 0.0;
                result[1] = 0.0;
                result[2] = 0.0;
                result[3] = 0.0;
                result[4] = 0.0;
                error_code = 3;
            }
        }
        error_code
    }
}

// ---------------------------------------------------------------------------
// Approx_CurveOnSurface_Eval3d (cxx L146-225)
// ---------------------------------------------------------------------------

/// OCCT Approx_CurveOnSurface_Eval3d (cxx L146-168) — evaluates the 3D
/// function of the curve-on-surface.
struct CurveOnSurfaceEval3d {
    /// OCCT: handle(Adaptor3d_Curve) fonct.
    fonct: Arc<dyn Adaptor3dCurve>,
    start_end_sav: [f64; 2],
}

impl EvaluatorFunction for CurveOnSurfaceEval3d {
    /// OCCT Approx_CurveOnSurface_Eval3d::Evaluate (cxx L170-225).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        let mut error_code = 0;
        let par = parameter;

        // Dimension is incorrect (cxx L181-184).
        if result.len() != 3 {
            error_code = 1;
        }

        // Parameter is incorrect (cxx L187-192).
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct = self.fonct.trim(start_end[0], start_end[1], p_confusion());
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

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
// Approx_CurveOnSurface_Eval2d (cxx L229-306)
// ---------------------------------------------------------------------------

/// OCCT Approx_CurveOnSurface_Eval2d (cxx L229-251) — evaluates the 2D
/// function of the pcurve.
struct CurveOnSurfaceEval2d {
    /// OCCT: handle(Adaptor2d_Curve2d) fonct2d.
    fonct2d: Curve2dHandle,
    start_end_sav: [f64; 2],
}

impl EvaluatorFunction for CurveOnSurfaceEval2d {
    /// OCCT Approx_CurveOnSurface_Eval2d::Evaluate (cxx L253-306).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        let mut error_code = 0;
        let par = parameter;

        // Dimension is incorrect (cxx L264-267).
        if result.len() != 2 {
            error_code = 1;
        }

        // Parameter is incorrect (cxx L270-275).
        if start_end[0] != self.start_end_sav[0] || start_end[1] != self.start_end_sav[1] {
            self.fonct2d = self.fonct2d.trim(start_end[0], start_end[1], p_confusion());
            self.start_end_sav[0] = start_end[0];
            self.start_end_sav[1] = start_end[1];
        }

        match derivative_request {
            0 => {
                let pnt = self.fonct2d.d0(par);
                result[0] = pnt.x;
                result[1] = pnt.y;
            }
            1 => {
                let (_pnt, v1) = self.fonct2d.d1(par);
                result[0] = v1.x;
                result[1] = v1.y;
            }
            2 => {
                let (_pnt, _v1, v2) = self.fonct2d.d2(par);
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
// Approx_CurveOnSurface
// ---------------------------------------------------------------------------

/// OCCT Approx_CurveOnSurface (hxx L28-139) — approximation of a curve on
/// surface.
pub struct ApproxCurveOnSurface {
    /// Input curve (hxx L118: const handle(Adaptor2d_Curve2d) myC2D).
    my_c2d: Curve2dHandle,
    /// Input surface (hxx L121: const handle(Adaptor3d_Surface) mySurf).
    my_surf: SurfaceHandle,
    /// First parameter of the result (hxx L124).
    my_first: f64,
    /// Last parameter of the result (hxx L127).
    my_last: f64,
    /// Tolerance (hxx L130).
    my_tol: f64,
    /// hxx L132: handle(Geom2d_BSplineCurve) myCurve2d.
    my_curve2d: Option<Curve2d>,
    /// hxx L133: handle(Geom_BSplineCurve) myCurve3d.
    my_curve3d: Option<BSplineCurve3>,
    /// hxx L134: bool myIsDone.
    my_is_done: bool,
    /// hxx L135: bool myHasResult.
    my_has_result: bool,
    /// hxx L136: double myError3d.
    my_error3d: f64,
    /// hxx L137: double myError2dU.
    my_error2d_u: f64,
    /// hxx L138: double myError2dV.
    my_error2d_v: f64,
}

impl ApproxCurveOnSurface {
    /// OCCT deprecated constructor calling Perform (cxx L310-332).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_perform(
        c2d: Curve2dHandle,
        surf: SurfaceHandle,
        first: f64,
        last: f64,
        tol: f64,
        continuity: GeomAbsShape,
        max_degree: i32,
        max_segments: i32,
        only3d: bool,
        only2d: bool,
    ) -> Self {
        let mut this = ApproxCurveOnSurface {
            my_c2d: c2d,
            my_surf: surf,
            my_first: first,
            my_last: last,
            my_tol: tol,
            my_curve2d: None,
            my_curve3d: None,
            my_is_done: false,
            my_has_result: false,
            my_error3d: 0.0,
            my_error2d_u: 0.0,
            my_error2d_v: 0.0,
        };
        this.perform(max_segments, max_degree, continuity, only3d, only2d);
        this
    }

    /// OCCT constructor that does not call Perform (cxx L336-352).
    pub fn new(
        the_c2d: Curve2dHandle,
        the_surf: SurfaceHandle,
        the_first: f64,
        the_last: f64,
        the_tol: f64,
    ) -> Self {
        ApproxCurveOnSurface {
            my_c2d: the_c2d,
            my_surf: the_surf,
            my_first: the_first,
            my_last: the_last,
            my_tol: the_tol,
            my_curve2d: None,
            my_curve3d: None,
            my_is_done: false,
            my_has_result: false,
            my_error3d: 0.0,
            my_error2d_u: 0.0,
            my_error2d_v: 0.0,
        }
    }

    /// OCCT Approx_CurveOnSurface::Perform (cxx L356-552).
    pub fn perform(
        &mut self,
        the_max_segments: i32,
        the_max_degree: i32,
        the_continuity: GeomAbsShape,
        the_only3d: bool,
        the_only2d: bool,
    ) {
        self.my_is_done = false;
        self.my_has_result = false;
        self.my_error2d_u = 0.0;
        self.my_error2d_v = 0.0;
        self.my_error3d = 0.0;

        if the_only3d && the_only2d {
            // OCCT: throw Standard_ConstructionError();
            panic!("Standard_ConstructionError");
        }

        // OCCT L373-385: normalize the continuity.  (GeomAbs_G1 / GeomAbs_G2
        // are not representable in rcad's GeomAbsShape; the >C2 restriction
        // maps C3/CN to C2.)
        let a_continuity = match the_continuity {
            GeomAbsShape::C0 => GeomAbsShape::C0,
            GeomAbsShape::C1 => GeomAbsShape::C1,
            GeomAbsShape::C2 => GeomAbsShape::C2,
            // The GeomAbs_G1 -> GeomAbs_C1 and GeomAbs_G2 -> GeomAbs_C2 arms
            // sit here in OCCT.
            GeomAbsShape::C3 | GeomAbsShape::CN => {
                GeomAbsShape::C2 // Restriction of AdvApprox_ApproxAFunction
            }
        };

        // OCCT L387: TrimmedC2D = myC2D->Trim(myFirst, myLast,
        // Precision::PConfusion()).
        let trimmed_c2d: Curve2dHandle = self
            .my_c2d
            .trim(self.my_first, self.my_last, p_confusion());

        // OCCT L389-399: the isoline fast path (only3d only).
        let mut is_u = false;
        let mut a_param = 0.0f64;
        let mut is_forward = false;
        if the_only3d
            && Self::is_iso_line(
                trimmed_c2d.as_ref(),
                &mut is_u,
                &mut a_param,
                &mut is_forward,
            )
        {
            if self.build_c3d_on_iso_line(trimmed_c2d.as_ref(), is_u, a_param, is_forward) {
                self.my_is_done = true;
                self.my_has_result = true;
                return;
            }
        }

        // OCCT L401: HCOnS = new Adaptor3d_CurveOnSurface(TrimmedC2D, mySurf).
        let hcons: Arc<dyn Adaptor3dCurve> = Arc::new(CurveOnSurface::new(
            trimmed_c2d.clone(),
            self.my_surf.clone(),
        ));

        // OCCT L403-406: Num1DSS / Num2DSS / Num3DSS and the tolerance arrays.
        let mut num1dss = 0i32;
        let num2dss = 0i32;
        let mut num3dss = 0i32;
        let mut one_d_tol: Vec<f64> = Vec::new();
        let two_d_tol_nul: Vec<f64> = Vec::new();
        let mut three_d_tol: Vec<f64> = Vec::new();

        // create evaluators and choose appropriate one (cxx L408-424).
        let mut eval3d = CurveOnSurfaceEval3d {
            fonct: hcons.clone(),
            start_end_sav: [self.my_first, self.my_last],
        };
        let mut eval2d = CurveOnSurfaceEval2d {
            fonct2d: trimmed_c2d.clone(),
            start_end_sav: [self.my_first, self.my_last],
        };
        let mut eval = CurveOnSurfaceEval::new(
            hcons.clone(),
            trimmed_c2d.clone(),
            self.my_first,
            self.my_last,
        );
        let eval_ptr: &mut dyn EvaluatorFunction = if the_only3d {
            &mut eval3d
        } else if the_only2d {
            &mut eval2d
        } else {
            &mut eval
        };

        // Initialization for 2d approximation (cxx L427-463).
        if !the_only3d {
            num1dss = 2;

            let mut tol_u = self.my_surf.u_resolution(self.my_tol) / 2.0;
            let mut tol_v = self.my_surf.v_resolution(self.my_tol) / 2.0;

            if self.my_surf.u_continuity() == GeomAbsShape::C0 {
                if !is_surf_g1(self.my_surf.as_ref(), true, ANGULAR) {
                    tol_u = 1.0e-3f64.min(1.0e3 * tol_u);
                }
                if !is_surf_g1(self.my_surf.as_ref(), true, CONFUSION) {
                    tol_u = 1.0e-3f64.min(1.0e2 * tol_u);
                }
            }

            if self.my_surf.v_continuity() == GeomAbsShape::C0 {
                if !is_surf_g1(self.my_surf.as_ref(), false, ANGULAR) {
                    tol_v = 1.0e-3f64.min(1.0e3 * tol_v);
                }
                if !is_surf_g1(self.my_surf.as_ref(), false, CONFUSION) {
                    tol_v = 1.0e-3f64.min(1.0e2 * tol_v);
                }
            }

            one_d_tol = vec![tol_u, tol_v];
        }

        if !the_only2d {
            num3dss = 1;
            three_d_tol = vec![self.my_tol / 2.0];
        }

        // OCCT L472-500: the cutting tool selection.
        if a_continuity <= self.my_c2d.continuity()
            && a_continuity <= self.my_surf.u_continuity()
            && a_continuity <= self.my_surf.v_continuity()
        {
            // OCCT: CutTool = new AdvApprox_DichoCutting();
            let cut_tool = DichoCutting;
            self.run_approxa_function(
                &cut_tool,
                num1dss,
                num2dss,
                num3dss,
                &one_d_tol,
                &two_d_tol_nul,
                &three_d_tol,
                a_continuity,
                the_max_degree,
                the_max_segments,
                the_only3d,
                the_only2d,
                eval_ptr,
            );
        } else if a_continuity == GeomAbsShape::C1 {
            // OCCT L479-489: NbIntervals/Intervals C1 + C2, then
            // CutTool = new AdvApprox_PrefAndRec(CutPnts_C1, CutPnts_C2).
            let nb_interv_c1 = hcons.nb_intervals(GeomAbsShape::C1);
            let mut cut_pnts_c1: Vec<f64> = hcons.intervals(GeomAbsShape::C1);
            cut_pnts_c1.resize(nb_interv_c1 + 1, 0.0);
            let nb_interv_c2 = hcons.nb_intervals(GeomAbsShape::C2);
            let mut cut_pnts_c2: Vec<f64> = hcons.intervals(GeomAbsShape::C2);
            cut_pnts_c2.resize(nb_interv_c2 + 1, 0.0);
            let _ = (cut_pnts_c1, cut_pnts_c2);
            panic!("GAP: AdvApprox_PrefAndRec not translated");
        } else {
            // OCCT L491-499: NbIntervals/Intervals C2 + C3, then
            // CutTool = new AdvApprox_PrefAndRec(CutPnts_C2, CutPnts_C3).
            let nb_interv_c2 = hcons.nb_intervals(GeomAbsShape::C2);
            let mut cut_pnts_c2: Vec<f64> = hcons.intervals(GeomAbsShape::C2);
            cut_pnts_c2.resize(nb_interv_c2 + 1, 0.0);
            let nb_interv_c3 = hcons.nb_intervals(GeomAbsShape::C3);
            let mut cut_pnts_c3: Vec<f64> = hcons.intervals(GeomAbsShape::C3);
            cut_pnts_c3.resize(nb_interv_c3 + 1, 0.0);
            let _ = (cut_pnts_c2, cut_pnts_c3);
            panic!("GAP: AdvApprox_PrefAndRec not translated");
        }
    }

    /// The OCCT L502-551 tail of Perform: construct the
    /// AdvApprox_ApproxAFunction and unpack its result.
    ///
    /// GAP (staged): the OCCT constructor receives `*CutTool`; the kernel
    /// ApproxAFunction runs its internal dichotomy cutting, so the selected
    /// tool is carried here for form until the kernel constructor lands.
    #[allow(clippy::too_many_arguments)]
    fn run_approxa_function(
        &mut self,
        _cut_tool: &dyn Cutting,
        num1dss: i32,
        num2dss: i32,
        num3dss: i32,
        one_d_tol: &[f64],
        two_d_tol_nul: &[f64],
        three_d_tol: &[f64],
        a_continuity: GeomAbsShape,
        the_max_degree: i32,
        the_max_segments: i32,
        the_only3d: bool,
        the_only2d: bool,
        eval_ptr: &mut dyn EvaluatorFunction,
    ) {
        // OCCT L502-514.
        let a_approx = ApproxAFunction::new(
            num1dss,
            num2dss,
            num3dss,
            Some(one_d_tol),
            Some(two_d_tol_nul),
            Some(three_d_tol),
            self.my_first,
            self.my_last,
            a_continuity,
            the_max_degree,
            the_max_segments,
            eval_ptr,
        );

        // OCCT L516: delete CutTool (RAII here).

        self.my_is_done = a_approx.is_done();
        self.my_has_result = a_approx.has_result();

        if self.my_has_result {
            // OCCT L523-525: Knots / Multiplicities / Degree.
            let knots = a_approx.knots_vec().to_vec();
            let mults = a_approx.multiplicities_vec().to_vec();
            let degree = a_approx.degree();

            let a_nb_poles = a_approx.nb_poles();
            if !the_only2d {
                // OCCT L530-533: Poles(1, Poles); myCurve3d = new
                // Geom_BSplineCurve(Poles, Knots->Array1(), Mults->Array1(),
                // Degree).
                let poles_flat = a_approx.poles_flat(1);
                let poles: Vec<DVec3> = (0..a_nb_poles)
                    .map(|i| {
                        DVec3::new(
                            poles_flat[i * 3],
                            poles_flat[i * 3 + 1],
                            poles_flat[i * 3 + 2],
                        )
                    })
                    .collect();
                self.my_curve3d = Some(BSplineCurve3 {
                    degree: degree as usize,
                    // Architecture difference: rcad BSplineCurve3 stores the
                    // full (multiplicity-expanded) knot vector; OCCT passes
                    // the distinct knots + multiplicities to the constructor.
                    knots: full_knots(&knots, &mults),
                    control_points: poles,
                    weights: vec![1.0; a_nb_poles],
                    is_periodic: false,
                });
                self.my_error3d = a_approx.max_error_at(3, 1);
            }
            if !the_only3d {
                // OCCT L537-550: Poles1dU/Poles1dV extraction, the
                // Geom2d_BSplineCurve construction and
                // myError2dU = MaxError(1, 1); myError2dV = MaxError(1, 2).
                // GAP: the kernel ApproxAFunction stores no 1D subspace
                // poles/errors (staged in rcad-kernel math/adv_approx).
                panic!("GAP: AdvApprox_ApproxAFunction 1D subspace storage (Poles1d / MaxError(1,N)) not translated");
            }
        }
    }

    /// OCCT Approx_CurveOnSurface::IsDone (cxx L554-557).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT Approx_CurveOnSurface::HasResult (cxx L559-562).
    pub fn has_result(&self) -> bool {
        self.my_has_result
    }

    /// OCCT Approx_CurveOnSurface::Curve3d (cxx L564-567).
    pub fn curve3d(&self) -> Option<BSplineCurve3> {
        self.my_curve3d.clone()
    }

    /// OCCT Approx_CurveOnSurface::Curve2d (cxx L569-572).
    pub fn curve2d(&self) -> Option<Curve2d> {
        self.my_curve2d.clone()
    }

    /// OCCT Approx_CurveOnSurface::MaxError3d (cxx L574-577).
    pub fn max_error3d(&self) -> f64 {
        self.my_error3d
    }

    /// OCCT Approx_CurveOnSurface::MaxError2dU (cxx L579-582).
    pub fn max_error2d_u(&self) -> f64 {
        self.my_error2d_u
    }

    /// OCCT Approx_CurveOnSurface::MaxError2dV (cxx L584-587).
    pub fn max_error2d_v(&self) -> f64 {
        self.my_error2d_v
    }

    /// OCCT Approx_CurveOnSurface::isIsoLine (cxx L591-675) — checks whether
    /// the 2d curve is a horizontal or vertical isoline.
    fn is_iso_line(
        the_c2d: &dyn Adaptor2dCurve2d,
        the_is_u: &mut bool,
        the_param: &mut f64,
        the_is_forward: &mut bool,
    ) -> bool {
        // These variables are used to check line state (vertical or
        // horizontal).
        let mut is_appropriate_type = false;
        let mut a_loc2d = DVec2::ZERO;
        let mut a_dir2d = DVec2::ZERO;

        // Test type (cxx L602-649).
        let a_type = the_c2d.get_type();
        if a_type == CurveType::Line {
            let a_lin2d = the_c2d.line();
            a_loc2d = a_lin2d.origin;
            a_dir2d = a_lin2d.direction;
            is_appropriate_type = true;
        } else if a_type == CurveType::BSpline {
            // OCCT: aBSpline2d = theC2D->BSpline() — the downcast; the rcad
            // trait models the null handle as None.
            let Some(a_bspline2d) = the_c2d.bspline() else {
                return false;
            };
            if a_bspline2d.degree != 1 || a_bspline2d.control_points.len() != 2 {
                return false; // Not a line or uneven parameterization.
            }

            a_loc2d = a_bspline2d.control_points[0];

            // Vector should be non-degenerated.
            let a_vec2d = a_bspline2d.control_points[1] - a_bspline2d.control_points[0];
            if a_vec2d.length_squared() < CONFUSION {
                return false; // Degenerated spline.
            }
            a_dir2d = a_vec2d;

            is_appropriate_type = true;
        } else if a_type == CurveType::Bezier {
            let Some(a_bezier2d) = the_c2d.bezier() else {
                return false;
            };
            // OCCT: aBezier2d->Degree() != 1 — a 2-pole Bezier is degree 1.
            if a_bezier2d.control_points.len() != 2 {
                return false; // Not a line or uneven parameterization.
            }

            a_loc2d = a_bezier2d.control_points[0];

            // Vector should be non-degenerated.
            let a_vec2d = a_bezier2d.control_points[1] - a_bezier2d.control_points[0];
            if a_vec2d.length_squared() < CONFUSION {
                return false; // Degenerated spline.
            }
            a_dir2d = a_vec2d;

            is_appropriate_type = true;
        }

        if !is_appropriate_type {
            return false;
        }

        // Check line to be vertical or horizontal (cxx L657-674).  OCCT:
        // aDir2d.IsParallel(gp::DX2d(), Precision::Angular()).
        if dir2d_is_parallel(a_dir2d, DVec2::X, ANGULAR) {
            // Horizontal line. V = const.
            *the_is_u = false;
            *the_param = a_loc2d.y;
            *the_is_forward = a_dir2d.dot(DVec2::X) > 0.0;
            return true;
        } else if dir2d_is_parallel(a_dir2d, DVec2::Y, ANGULAR) {
            // Vertical line. U = const.
            *the_is_u = true;
            *the_param = a_loc2d.x;
            *the_is_forward = a_dir2d.dot(DVec2::Y) > 0.0;
            return true;
        }

        false
    }

    /// OCCT Approx_CurveOnSurface::buildC3dOnIsoLine (cxx L681-813).
    fn build_c3d_on_iso_line(
        &mut self,
        the_c2d: &dyn Adaptor2dCurve2d,
        the_is_u: bool,
        the_param: f64,
        the_is_forward: bool,
    ) -> bool {
        // Convert adapter to the appropriate type (cxx L687-691): the OCCT
        // down_cast to GeomAdaptor_Surface maps to checking the adaptor wraps
        // a kernel surface.
        let Some(kernel_surf) = self.my_surf.kernel_surface() else {
            return false;
        };
        let mut a_surf: Surface3 = kernel_surf.clone();

        if self.my_surf.get_type() == GeomAbsSurfaceType::Sphere {
            return false;
        }

        // Extract isoline (cxx L699-703).
        let a_f2d = the_c2d.value(the_c2d.first_parameter());
        let a_l2d = the_c2d.value(the_c2d.last_parameter());

        let mut is_to_trim = true;
        // OCCT L706-707: aSurf->Bounds(U1, U2, V1, V2).
        let domain = a_surf.default_domain();
        let (u1, u2, v1, v2) = (domain[0], domain[1], domain[2], domain[3]);

        let a_c3d: Curve3;
        if the_is_u {
            let mut a_v1_param = a_f2d.y.min(a_l2d.y);
            let mut a_v2_param = a_f2d.y.max(a_l2d.y);
            if a_v2_param < v1 - self.my_tol || a_v1_param > v2 + self.my_tol {
                return false;
            // OCCT Approx_CurveOnSurface.cxx L717:
            // Precision::IsInfinite(V1) || Precision::IsInfinite(V2).
            } else if rcad_kernel::precision::is_infinite_value(v1)
                || rcad_kernel::precision::is_infinite_value(v2)
            {
                if (a_v2_param - a_v1_param).abs() < p_confusion() {
                    return false;
                }
                // OCCT L723: aSurf = new Geom_RectangularTrimmedSurface(aSurf,
                // U1, U2, aV1Param, aV2Param); isToTrim = false.
                a_surf = surface_rectangular_trimmed(&a_surf, u1, u2, a_v1_param, a_v2_param);
                is_to_trim = false;
            } else {
                a_v1_param = a_v1_param.max(v1);
                a_v2_param = a_v2_param.min(v2);
                if (a_v2_param - a_v1_param).abs() < p_confusion() {
                    return false;
                }
            }
            // OCCT L735-739: aC3d = aSurf->UIso(theParam);
            // if (isToTrim) aC3d = new Geom_TrimmedCurve(aC3d, aV1Param,
            // aV2Param).
            let iso = surface_u_iso(&a_surf, the_param);
            a_c3d = if is_to_trim {
                Curve3::Trimmed(TrimmedCurve3::new(iso, a_v1_param, a_v2_param))
            } else {
                iso
            };
        } else {
            let mut a_u1_param = a_f2d.x.min(a_l2d.x);
            let mut a_u2_param = a_f2d.x.max(a_l2d.x);
            if a_u2_param < u1 - self.my_tol || a_u1_param > u2 + self.my_tol {
                return false;
            // OCCT Approx_CurveOnSurface.cxx L749:
            // Precision::IsInfinite(U1) || Precision::IsInfinite(U2).
            } else if rcad_kernel::precision::is_infinite_value(u1)
                || rcad_kernel::precision::is_infinite_value(u2)
            {
                if (a_u2_param - a_u1_param).abs() < p_confusion() {
                    return false;
                }
                // OCCT L755: aSurf = new Geom_RectangularTrimmedSurface(aSurf,
                // aU1Param, aU2Param, V1, V2); isToTrim = false.
                a_surf = surface_rectangular_trimmed(&a_surf, a_u1_param, a_u2_param, v1, v2);
                is_to_trim = false;
            } else {
                a_u1_param = a_u1_param.max(u1);
                a_u2_param = a_u2_param.min(u2);
                if (a_u2_param - a_u1_param).abs() < p_confusion() {
                    return false;
                }
            }
            // OCCT L767-771: aC3d = aSurf->VIso(theParam);
            // if (isToTrim) aC3d = new Geom_TrimmedCurve(aC3d, aU1Param,
            // aU2Param).
            let iso = surface_v_iso(&a_surf, the_param);
            a_c3d = if is_to_trim {
                Curve3::Trimmed(TrimmedCurve3::new(iso, a_u1_param, a_u2_param))
            } else {
                iso
            };
        }

        // OCCT L775: myCurve3d = GeomConvert::CurveToBSplineCurve(aC3d,
        // Convert_QuasiAngular).  Architecture difference: the kernel
        // converter is sampling-based (the exact GeomConvert conversion is
        // staged); the path is behind the UIso/VIso GAP above.
        self.my_curve3d = Some(geom_convert_curve_to_bspline(&a_c3d));
        if !the_is_forward {
            // OCCT L778: myCurve3d->Reverse().
            reverse_bspline3(self.my_curve3d.as_mut().unwrap());
        }

        // OCCT L783-785: rebuild the parameterization to match the 2d curve:
        // BSplCLib::Reparametrize(First, Last, Knots); SetKnots.
        // Architecture difference: the kernel curve stores the full
        // (multiplicity-expanded) knot vector; the linear map applies
        // entry-wise to it.
        let mut a_knots = self.my_curve3d.as_ref().unwrap().knots.clone();
        bspl_lib::reparametrize(the_c2d.first_parameter(), the_c2d.last_parameter(), &mut a_knots);
        self.my_curve3d.as_mut().unwrap().knots = a_knots;

        // Evaluate error (cxx L788-806).
        self.my_error3d = 0.0;

        let a_par_f = self.my_first;
        let a_par_l = self.my_last;
        let a_nb_pnt = 23usize;
        for an_idx in 0..=a_nb_pnt {
            let a_par = a_par_f + ((a_par_l - a_par_f) * an_idx as f64) / a_nb_pnt as f64;

            let a_pnt2d = the_c2d.value(a_par);

            let a_pnt_c3d = self.my_curve3d.as_ref().unwrap().point_at(a_par);
            let a_pnt_c2d = a_surf.point_at(a_pnt2d.x, a_pnt2d.y);

            let a_sq_deviation = a_pnt_c3d.distance_squared(a_pnt_c2d);
            self.my_error3d = a_sq_deviation.max(self.my_error3d);
        }

        self.my_error3d = self.my_error3d.sqrt();

        // Target tolerance is not obtained. This situation happens for
        // isolines on the sphere (cxx L808-812).
        self.my_error3d <= self.my_tol
    }
}

// ---------------------------------------------------------------------------
// GAP leaves of the isoline path (staged in the kernel geom package)
// ---------------------------------------------------------------------------

/// OCCT Geom_RectangularTrimmedSurface construction (cxx L723 / L755).
fn surface_rectangular_trimmed(
    _surf: &Surface3,
    _u1: f64,
    _u2: f64,
    _v1: f64,
    _v2: f64,
) -> Surface3 {
    panic!("GAP: Geom_RectangularTrimmedSurface not translated")
}

/// OCCT Geom_Surface::UIso — the virtual dispatch of the concrete surface
/// classes: Geom_Plane / Geom_CylindricalSurface / Geom_ConicalSurface /
/// Geom_SphericalSurface / Geom_ToroidalSurface (all via the ElSLib
/// constructors) and Geom_SurfaceOfRevolution (cxx L372-379: a rotated copy of
/// the basis curve).
fn surface_u_iso(surf: &Surface3, param: f64) -> Curve3 {
    use rcad_kernel::base::proj_lib::elslib_iso as el;
    match surf {
        Surface3::Plane(p) => Curve3::Line(el::elslib_plane_u_iso(
            &el::Ax3View::from_axes(p.origin, p.normal, p.u_dir),
            param,
        )),
        Surface3::Cylinder(c) => Curve3::Line(el::elslib_cylinder_u_iso(
            &el::Ax3View::from_axes(c.origin, c.axis, c.ref_dir),
            c.radius,
            param,
        )),
        Surface3::Cone(c) => Curve3::Line(el::elslib_cone_u_iso(
            &el::Ax3View::from_axes(c.apex, c.axis, c.ref_dir),
            c.radius,
            c.half_angle_rad,
            param,
        )),
        Surface3::Sphere(s) => Curve3::Circle(el::elslib_sphere_u_iso(
            &el::Ax3View::from_axes(s.center, s.axis, s.ref_dir),
            s.radius,
            param,
        )),
        Surface3::Torus(t) => Curve3::Circle(el::elslib_torus_u_iso(
            &el::Ax3View::from_axes(t.center, t.axis, t.ref_dir),
            t.major_radius,
            t.minor_radius,
            param,
        )),
        // OCCT Geom_SurfaceOfRevolution::UIso (cxx L372-379):
        //   C = basisCurve->Copy(); C->Rotate(Ax1(loc, direction), U); return C.
        Surface3::Revolution(r) => rotate_curve_about_axis(
            &r.profile,
            r.axis_origin,
            r.axis_dir,
            param,
        ),
        // OCCT Geom_BSplineSurface::UIso (Geom_BSplineSurface_1.cxx L598-635):
        // BSplSLib::Iso on the U direction; the result curve's knots and
        // degree are the V direction's.
        Surface3::BSpline(b) => {
            let (cpoles, cweights) = bspl_lib::bspl_slib_iso(
                param,
                true,
                b.degree_u,
                &b.knots_u,
                &b.control_points,
                &b.weights,
            );
            Curve3::BSpline(BSplineCurve3 {
                degree: b.degree_v,
                knots: b.knots_v.clone(),
                control_points: cpoles,
                weights: cweights,
                is_periodic: false,
            })
        }
        _ => panic!("GAP: Geom_Surface::UIso not translated for this surface type"),
    }
}

/// OCCT Geom_Surface::VIso — the same virtual dispatch, the V-isoparametric
/// counterpart (Geom_SurfaceOfRevolution::VIso cxx L383-410: the parallel
/// circle of the basis point through the axis).
fn surface_v_iso(surf: &Surface3, param: f64) -> Curve3 {
    use rcad_kernel::base::proj_lib::elslib_iso as el;
    match surf {
        Surface3::Plane(p) => Curve3::Line(el::elslib_plane_v_iso(
            &el::Ax3View::from_axes(p.origin, p.normal, p.u_dir),
            param,
        )),
        Surface3::Cylinder(c) => Curve3::Circle(el::elslib_cylinder_v_iso(
            &el::Ax3View::from_axes(c.origin, c.axis, c.ref_dir),
            c.radius,
            param,
        )),
        Surface3::Cone(c) => Curve3::Circle(el::elslib_cone_v_iso(
            &el::Ax3View::from_axes(c.apex, c.axis, c.ref_dir),
            c.radius,
            c.half_angle_rad,
            param,
        )),
        Surface3::Sphere(s) => Curve3::Circle(el::elslib_sphere_v_iso(
            &el::Ax3View::from_axes(s.center, s.axis, s.ref_dir),
            s.radius,
            param,
        )),
        Surface3::Torus(t) => Curve3::Circle(el::elslib_torus_v_iso(
            &el::Ax3View::from_axes(t.center, t.axis, t.ref_dir),
            t.major_radius,
            t.minor_radius,
            param,
        )),
        // OCCT Geom_BSplineSurface::VIso (Geom_BSplineSurface_1.cxx L775-812):
        // BSplSLib::Iso on the V direction; the result curve's knots and
        // degree are the U direction's.
        Surface3::BSpline(b) => {
            let (cpoles, cweights) = bspl_lib::bspl_slib_iso(
                param,
                false,
                b.degree_v,
                &b.knots_v,
                &b.control_points,
                &b.weights,
            );
            Curve3::BSpline(BSplineCurve3 {
                degree: b.degree_u,
                knots: b.knots_u.clone(),
                control_points: cpoles,
                weights: cweights,
                is_periodic: false,
            })
        }
        // OCCT Geom_SurfaceOfRevolution::VIso (cxx L383-410): the circle of the
        // basis point at V about the axis.  Rad = distance from the axis; the
        // circle frame is gp_Ax2(C, direction, D) where C is the projection of
        // the basis point onto the axis and D the unit vector from C to it.
        Surface3::Revolution(r) => {
            let pc = r.profile.point_at(param);
            let d = pc - r.axis_origin;
            let rad = (d - r.axis_dir * d.dot(r.axis_dir)).length();
            let c = r.axis_origin + r.axis_dir * d.dot(r.axis_dir);
            let radial = pc - c;
            let (normal, x_dir) = if rad > rcad_kernel::precision::CONFUSION {
                (r.axis_dir, radial.normalize_or_zero())
            } else {
                // OCCT: Rep = gp_Ax2(C, direction) — the zero-radius case uses
                // the frame's default X direction.
                let x = rcad_kernel::geom::any_perpendicular(r.axis_dir);
                (r.axis_dir, x)
            };
            Curve3::Circle(rcad_kernel::geom::Circle3 {
                center: c,
                normal,
                x_dir,
                y_dir: normal.cross(x_dir).normalize_or_zero(),
                radius: rad,
            })
        }
        _ => panic!("GAP: Geom_Surface::VIso not translated for this surface type"),
    }
}

/// OCCT Geom_Geometry::Rotate (Geom_Geometry.cxx L47-55) applied to a curve:
/// the per-type `Rotate(Ax1, Ang)` overrides.  The rotation is a rigid motion
/// so a TrimmedCurve rotates its basis and keeps its parameter range
/// (Geom_TrimmedCurve::Transform).
fn rotate_curve_about_axis(
    curve: &Curve3,
    axis_loc: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Curve3 {
    let axis_dir = axis_dir.normalize_or_zero();
    let trsf = glam::DAffine3::from_translation(axis_loc)
        * glam::DAffine3::from_axis_angle(axis_dir, angle)
        * glam::DAffine3::from_translation(-axis_loc);
    match curve {
        Curve3::Trimmed(t) => Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            rotate_curve_about_axis(&t.curve, axis_loc, axis_dir, angle),
            t.first,
            t.last,
        )),
        other => rcad_kernel::geom::transform_curve(other, &trsf),
    }
}

/// OCCT GeomConvert::CurveToBSplineCurve(C, Convert_QuasiAngular)
/// (GeomConvert.cxx L157-209).
///
/// The kernel's exact conversion has landed the Trimmed(Line) and
/// Trimmed(Circle) arms (GeomConvert.cxx L200-282).  The isoline path only
/// ever converts `Geom_TrimmedCurve`s here (BuildC3dOnIsoLine L735-739 wraps
/// `aSurf->UIso/VIso` in a trimmed curve), and the rcad TrimmedCurve3
/// `map_param` is the identity, so the outermost trim supplies the domain and
/// an inner trim is inert.  Inputs outside the landed arms keep the kernel's
/// sampling conversion (the pre-existing stand-in for the staged arms).
fn geom_convert_curve_to_bspline(c: &Curve3) -> BSplineCurve3 {
    if let Curve3::Trimmed(tc) = c {
        let mut basis = tc.basis_curve();
        let mut depth = 0;
        while let Curve3::Trimmed(inner) = basis {
            basis = inner.basis_curve();
            depth += 1;
            if depth > 8 {
                break;
            }
        }
        if matches!(basis, Curve3::Line(_) | Curve3::Circle(_)) {
            return rcad_kernel::base::convert::geom_convert_curve_to_bspline_curve(
                &Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
                    basis.clone(),
                    tc.first,
                    tc.last,
                )),
                rcad_kernel::base::convert::ConvertParameterisation::QuasiAngular,
            );
        }
    }
    rcad_kernel::base::convert::curve_to_bspline(c, 23)
}

/// OCCT Geom_BSplineCurve::Reverse (Geom_BSplineCurve_1.cxx): reverse the
/// parameterization — poles reversed, knots negated and mirrored.
///
/// Architecture difference: carried here because the kernel BSplineCurve3 has
/// no Reverse member yet (staged to move into the kernel geom package); the
/// full (multiplicity-expanded) knot vector is negated and reversed as a
/// whole, which matches the OCCT distinct-knot operation.
fn reverse_bspline3(c: &mut BSplineCurve3) {
    c.control_points.reverse();
    c.knots = c.knots.iter().rev().map(|k| -k).collect();
}

// ---------------------------------------------------------------------------
// Local helpers (architecture differences)
// ---------------------------------------------------------------------------

/// Expand the distinct knots + multiplicities into the full knot vector the
/// rcad BSplineCurve representation stores (the OCCT constructor consumes the
/// pair directly).
fn full_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut out = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            out.push(*k);
        }
    }
    out
}

/// OCCT gp_Dir2d::IsParallel(Other, AngularTolerance) — true when the angular
/// distance between the two directions is within the tolerance (mod pi).
fn dir2d_is_parallel(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    let ang_a = a.y.atan2(a.x);
    let ang_b = b.y.atan2(b.x);
    let diff = (ang_a - ang_b).abs();
    let diff = diff.min(std::f64::consts::PI * 2.0 - diff);
    diff <= angular_tolerance || (diff - std::f64::consts::PI).abs() <= angular_tolerance
}
