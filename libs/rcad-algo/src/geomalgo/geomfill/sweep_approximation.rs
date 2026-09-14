//! OCCT Approx_SweepApproximation (TKGeomBase/Approx) — 1:1 translation of
//! Approx_SweepApproximation.hxx (L50-227) + Approx_SweepApproximation.cxx
//! (whole file L31-872) plus the .lxx accessors, over
//! AdvApprox_ApproxAFunction.
//!
//! File split: this engine used to live inside `sweep.rs` (GeomFill_Sweep);
//! it is wired back as a `#[path]` child module of `sweep.rs` and
//! re-exported from there, so the consumers (`GeomFill_Pipe::Perform` in
//! pipe.rs, `GeomFill_Sweep::BuildAll` in sweep.rs) keep their original
//! paths.  Pure file split — nothing changed semantically.
//!
//! Architecture differences:
//! - `gp_GTrsf2d` maps to the local [`GpGTrsf2d`] re-host over the members
//!   consumed by the approximation (identity constructor, SetValue, Invert,
//!   Inverted, Transforms).
//! - `Approx_SweepApproximation_Eval` (cxx L31-58) holds a reference back
//!   to the approximation; the Rust form holds the mutable borrow,
//!   constructed only inside [`ApproxSweepApproximation::approximation`]
//!   around the AdvApprox_ApproxAFunction construction.
//! - `handle(Approx_SweepFunction) myFunc` (hxx L188) maps to an owned box:
//!   the function object is created by the consumer solely for this
//!   approximation (the only consumers, GeomFill_Pipe::Perform
//!   L813-820 and GeomFill_Sweep::BuildAll, never touch their local
//!   function object after handing it over, so the shared-handle mutation
//!   visibility is preserved).

use glam::{DVec2, DVec3};

use rcad_kernel::math::adv_approx::{
    ApproxAFunction, Cutting, DichoCutting, EvaluatorFunction, PrefAndRec,
};
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::math::GeomAbsShape;

use super::super::approx_sweep_function::ApproxSweepFunction;
// ---------------------------------------------------------------------------
// OCCT Approx_SweepApproximation (TKGeomBase/Approx)
// ---------------------------------------------------------------------------

/// OCCT gp_GTrsf2d — local re-host of the members consumed by
/// Approx_SweepApproximation: the identity constructor (gp_GTrsf2d.hxx
/// L51-58), SetValue (hxx L255-270), Invert (gp_GTrsf2d.cxx L56-70),
/// Inverted (hxx L144-150), Transforms (hxx L287-297), over
/// gp_Mat2d::Invert (gp_Mat2d.cxx L125-140) and gp_XY::Multiply(gp_Mat2d)
/// (gp_XY.hxx L401-406).
#[derive(Clone, Copy)]
struct GpGTrsf2d {
    /// OCCT gp_Mat2d matrix (the vectorial part).
    matrix: [[f64; 2]; 2],
    /// OCCT gp_XY loc (the translation part).
    loc: DVec2,
    /// OCCT shape == gp_Other (SetValue switches the general form on).
    is_other: bool,
    /// OCCT double scale.
    scale: f64,
}

impl GpGTrsf2d {
    /// OCCT gp_GTrsf2d() — identity (gp_GTrsf2d.hxx L51-58).
    fn new() -> Self {
        GpGTrsf2d {
            matrix: [[1.0, 0.0], [0.0, 1.0]],
            loc: DVec2::ZERO,
            is_other: false,
            scale: 1.0,
        }
    }

    /// OCCT gp_GTrsf2d::SetValue (hxx L255-270) — columns 1-2 go to the
    /// matrix, column 3 to the translation part; the form becomes gp_Other.
    fn set_value(&mut self, the_row: usize, the_col: usize, the_value: f64) {
        if the_col == 3 {
            self.loc[the_row - 1] = the_value;
        } else {
            self.matrix[the_row - 1][the_col - 1] = the_value;
        }
        self.is_other = true;
    }

    /// OCCT gp_XY::Multiply(gp_Mat2d) (gp_XY.hxx L401-406) — row-vector
    /// form: (x y) * M.
    fn xy_multiplied_matrix(&self, coord: DVec2) -> DVec2 {
        DVec2::new(
            self.matrix[0][0] * coord.x + self.matrix[0][1] * coord.y,
            self.matrix[1][0] * coord.x + self.matrix[1][1] * coord.y,
        )
    }

    /// OCCT gp_GTrsf2d::Invert (gp_GTrsf2d.cxx L56-70) — the gp_Other
    /// branch; the Trsf2d fallback (cxx L64-68) is unreachable here: the
    /// AAffin entries always carry gp_Other after Perform (SetValue was
    /// applied in step 2.2).
    fn invert(&mut self) {
        if self.is_other {
            // OCCT gp_Mat2d::Invert (gp_Mat2d.cxx L125-140).
            let a_new_mat = [
                [self.matrix[1][1], -self.matrix[0][1]],
                [-self.matrix[1][0], self.matrix[0][0]],
            ];
            let mut a_det = a_new_mat[0][0] * a_new_mat[1][1] - a_new_mat[0][1] * a_new_mat[1][0];
            if a_det.abs() <= GP_RESOLUTION {
                panic!("gp_Mat2d::Invert() - matrix has zero determinant");
            }
            a_det = 1.0 / a_det;
            self.matrix[0][0] = a_new_mat[0][0] * a_det;
            self.matrix[1][0] = a_new_mat[1][0] * a_det;
            self.matrix[0][1] = a_new_mat[0][1] * a_det;
            self.matrix[1][1] = a_new_mat[1][1] * a_det;
            // OCCT: loc.Multiply(matrix); loc.Reverse().
            self.loc = self.xy_multiplied_matrix(self.loc);
            self.loc = -self.loc;
        } else {
            panic!("unreachable: gp_GTrsf2d::Invert Trsf2d branch (gp_GTrsf2d.cxx L64-68)");
        }
    }

    /// OCCT gp_GTrsf2d::Inverted (hxx L144-150).
    fn inverted(&self) -> GpGTrsf2d {
        let mut a_t = *self;
        a_t.invert();
        a_t
    }

    /// OCCT gp_GTrsf2d::Transforms(gp_XY&) (hxx L287-297).
    fn transforms(&self, coord: DVec2) -> DVec2 {
        let mut coord = self.xy_multiplied_matrix(coord);
        if !self.is_other && self.scale != 1.0 {
            coord *= self.scale;
        }
        coord + self.loc
    }
}

/// OCCT Approx_SweepApproximation_Eval (cxx L31-58) — the AdvApprox
/// evaluator forwarding Evaluate to Approx_SweepApproximation::Eval.  The
/// OCCT class holds a reference back to the approximation; the Rust form
/// holds the mutable borrow (constructed only inside Approximation, around
/// the AdvApprox_ApproxAFunction construction).
struct ApproxSweepApproximationEval<'a> {
    tool: &'a mut ApproxSweepApproximation,
}

impl EvaluatorFunction for ApproxSweepApproximationEval<'_> {
    /// OCCT Approx_SweepApproximation_Eval::Evaluate (cxx L50-58) — the
    /// `result` slice plays the role of Result (Result[0] is the double&
    /// passed down to Eval); the returned code is Eval's ier.
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        self.tool
            .eval(parameter, derivative_request, start_end[0], start_end[1], result)
    }
}

/// OCCT Approx_SweepApproximation (Approx_SweepApproximation.hxx L50-227) —
/// approximation of a swept surface S(u,v) (and associated 2d curves)
/// defined by a section law; drives AdvApprox_ApproxAFunction through the
/// Approx_SweepApproximation_Eval evaluator.
pub struct ApproxSweepApproximation {
    /// OCCT handle(Approx_SweepFunction) myFunc (hxx L188) — an owned box:
    /// the function object is created by the consumer solely for this
    /// approximation (architecture note for the shared handle).
    my_func: Box<dyn ApproxSweepFunction>,
    /// OCCT bool done (hxx L189).
    done: bool,
    /// OCCT int Num1DSS, Num2DSS, Num3DSS (hxx L190-192).
    num1dss: usize,
    num2dss: usize,
    num3dss: usize,
    /// OCCT int udeg, vdeg, deg2d (hxx L193-195).
    udeg: usize,
    vdeg: usize,
    deg2d: usize,
    /// OCCT handle(HArray2<gp_Pnt>) tabPoles / HArray2<double> tabWeights
    /// (hxx L196-197) — row-major [section][pole].
    tab_poles: Vec<Vec<DVec3>>,
    tab_weights: Vec<Vec<f64>>,
    /// OCCT tabUKnots / tabVKnots / tab2dKnots (hxx L198-200).
    tab_u_knots: Vec<f64>,
    tab_v_knots: Vec<f64>,
    tab_2d_knots: Vec<f64>,
    /// OCCT tabUMults / tabVMults / tab2dMults (hxx L201-203).
    tab_u_mults: Vec<i32>,
    tab_v_mults: Vec<i32>,
    tab_2d_mults: Vec<i32>,
    /// OCCT NCollection_Sequence<HArray1<gp_Pnt2d>> seqPoles2d (hxx L204).
    seq_poles2d: Vec<Vec<DVec2>>,
    /// OCCT MError1d / tab2dError / MError3d / AError1d / Ave2dError /
    /// AError3d (hxx L205-210).
    m_error1d: Vec<f64>,
    tab2d_error: Vec<f64>,
    m_error3d: Vec<f64>,
    a_error1d: Vec<f64>,
    ave2d_error: Vec<f64>,
    a_error3d: Vec<f64>,
    /// OCCT handle(HArray1<gp_GTrsf2d>) AAffin + HArray1<double> COnSurfErr
    /// (hxx L211-212).
    a_affin: Vec<GpGTrsf2d>,
    c_on_surf_err: Vec<f64>,
    /// OCCT gp_Vec Translation (hxx L213).
    translation: DVec3,
    /// OCCT myPoles / myPoles2d / myWeigths (hxx L214-216) — the D0 scratch.
    my_poles: Vec<DVec3>,
    my_poles2d: Vec<DVec2>,
    my_weigths: Vec<f64>,
    /// OCCT myDPoles / myD2Poles / myDPoles2d / myD2Poles2d (hxx L217-220).
    my_dpoles: Vec<DVec3>,
    my_d2poles: Vec<DVec3>,
    my_dpoles2d: Vec<DVec2>,
    my_d2poles2d: Vec<DVec2>,
    /// OCCT myDWeigths / myD2Weigths (hxx L221-222).
    my_dweigths: Vec<f64>,
    my_d2weigths: Vec<f64>,
    /// OCCT int myOrder / double myParam / first / last (hxx L223-226) — the
    /// control variables.
    my_order: i32,
    my_param: f64,
    first: f64,
    last: f64,
}

impl ApproxSweepApproximation {
    /// OCCT Approx_SweepApproximation::Approx_SweepApproximation
    /// (cxx L60-69).
    ///
    /// Architecture note: the OCCT handle shares Func with the caller; the
    /// Rust form clones the concrete function into the engine — the only
    /// consumers (GeomFill_Pipe::Perform L813-820, GeomFill_Sweep::BuildAll)
    /// never touch their local function object after handing it over, so
    /// the shared-handle mutation visibility is preserved.
    pub fn new<F: ApproxSweepFunction + Clone + 'static>(func: &F) -> Self {
        Self::with_owned(Box::new(func.clone()))
    }

    /// The handle-of-abstract constructor form (OCCT
    /// `Approx_SweepApproximation(const handle(Approx_SweepFunction)&)`)
    /// for an already-boxed function — used by GeomFill_Sweep::BuildAll.
    pub fn with_owned(func: Box<dyn ApproxSweepFunction>) -> Self {
        //  Init of variables of control
        ApproxSweepApproximation {
            my_func: func,
            done: false,
            num1dss: 0,
            num2dss: 0,
            num3dss: 0,
            udeg: 0,
            vdeg: 0,
            deg2d: 0,
            tab_poles: Vec::new(),
            tab_weights: Vec::new(),
            tab_u_knots: Vec::new(),
            tab_v_knots: Vec::new(),
            tab_2d_knots: Vec::new(),
            tab_u_mults: Vec::new(),
            tab_v_mults: Vec::new(),
            tab_2d_mults: Vec::new(),
            seq_poles2d: Vec::new(),
            m_error1d: Vec::new(),
            tab2d_error: Vec::new(),
            m_error3d: Vec::new(),
            a_error1d: Vec::new(),
            ave2d_error: Vec::new(),
            a_error3d: Vec::new(),
            a_affin: Vec::new(),
            c_on_surf_err: Vec::new(),
            translation: DVec3::ZERO,
            my_poles: Vec::new(),
            my_poles2d: Vec::new(),
            my_weigths: Vec::new(),
            my_dpoles: Vec::new(),
            my_d2poles: Vec::new(),
            my_dpoles2d: Vec::new(),
            my_d2poles2d: Vec::new(),
            my_dweigths: Vec::new(),
            my_d2weigths: Vec::new(),
            my_order: -1,
            my_param: 0.0,
            first: 1.0e100,
            last: -1.0e100,
        }
    }

    /// OCCT Approx_SweepApproximation::Perform (cxx L71-274).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        first: f64,
        last: f64,
        tol3d: f64,
        bound_tol: f64,
        tol2d: f64,
        tol_angular: f64,
        continuity: GeomAbsShape,
        degmax: i32,
        segmax: i32,
    ) {
        let mut continuity = continuity;
        // OCCT: int NbPolSect, NbKnotSect; double Tol, Tol3dMin = Tol3d,
        // The3D2DTol = 0.
        let mut nb_pol_sect = 0usize;
        let mut nb_knot_sect = 0usize;
        let mut tol3d_min = tol3d;
        let mut the3d2dtol = 0.0f64;

        // (1) Characteristics of a section
        self.my_func.section_shape(&mut nb_pol_sect, &mut nb_knot_sect, &mut self.udeg);
        self.num2dss = self.my_func.nb_2d_curves();
        self.tab_u_knots = vec![0.0; nb_knot_sect];
        self.tab_u_mults = vec![0; nb_knot_sect];
        self.my_func.knots(&mut self.tab_u_knots);
        self.my_func.mults(&mut self.tab_u_mults);

        // (2) Decompositition into sub-spaces
        let two_d_tol: Option<Vec<f64>>;
        self.num3dss = nb_pol_sect;

        // (2.1) Tolerance 3d and 1d
        let mut one_d_tol = vec![0.0f64; self.num3dss];
        let mut three_d_tol = vec![0.0f64; self.num3dss];

        self.my_func
            .get_tolerance(bound_tol, tol3d, tol_angular, &mut three_d_tol);

        for ii in 1..=self.num3dss {
            if three_d_tol[ii - 1] < tol3d_min {
                tol3d_min = three_d_tol[ii - 1];
            }
        }

        if self.my_func.is_rational() {
            self.num1dss = nb_pol_sect;
            let mut wmin = vec![0.0f64; self.num1dss];
            self.my_func.get_minimal_weight(&mut wmin);
            let size = self.my_func.maximal_section();
            self.translation = self.my_func.barycentre_of_surf();
            for ii in 1..=self.num3dss {
                let mut tol = three_d_tol[ii - 1] / 2.0; // To take account of the error on the final result.
                one_d_tol[ii - 1] = tol * wmin[ii - 1] / size;
                tol *= wmin[ii - 1]; // Factor of projection
                three_d_tol[ii - 1] = tol.max(1.0e-20);
            }
        } else {
            self.num1dss = 0;
        }

        // (2.2) Tolerance and Transformation 2d.
        if self.num2dss == 0 {
            two_d_tol = None; // OCCT: TwoDTol.Nullify().
        } else {
            // for 2d define affinity using resolutions, to
            // avoid homogeneous tolerance of approximation (u/v and 2d/3d)
            let mut two_d_tol_v = vec![0.0f64; self.num2dss];
            self.a_affin = vec![GpGTrsf2d::new(); self.num2dss];
            the3d2dtol = 0.9 * bound_tol; // 10% of security
            for ii in 1..=self.num2dss {
                let mut tolu = 0.0f64;
                let mut tolv = 0.0f64;
                self.my_func.resolution(ii, the3d2dtol, &mut tolu, &mut tolv);
                let res;
                if tolu > tolv {
                    res = tolv;
                    self.a_affin[ii - 1].set_value(1, 1, tolv / tolu);
                } else {
                    res = tolu;
                    self.a_affin[ii - 1].set_value(2, 2, tolu / tolv);
                }
                two_d_tol_v[ii - 1] = tol2d.min(res);
            }
            two_d_tol = Some(two_d_tol_v);
        }

        // (3) Approximation
        // Init
        self.my_poles = vec![DVec3::ZERO; self.num3dss];
        self.my_dpoles = vec![DVec3::ZERO; self.num3dss];
        self.my_d2poles = vec![DVec3::ZERO; self.num3dss];

        self.my_weigths = vec![0.0; self.num3dss];
        self.my_dweigths = vec![0.0; self.num3dss];
        self.my_d2weigths = vec![0.0; self.num3dss];

        if self.num2dss > 0 {
            self.my_poles2d = vec![DVec2::ZERO; self.num2dss];
            self.my_dpoles2d = vec![DVec2::ZERO; self.num2dss];
            self.my_d2poles2d = vec![DVec2::ZERO; self.num2dss];
            self.c_on_surf_err = vec![0.0; self.num2dss];
        } else {
            self.my_poles2d = Vec::new();
            self.my_dpoles2d = Vec::new();
            self.my_d2poles2d = Vec::new();
            self.c_on_surf_err = Vec::new();
        }

        // Checks if myFunc->D2 is implemented
        if continuity >= GeomAbsShape::C2 {
            let b = self.my_func.d2(
                first,
                first,
                last,
                &mut self.my_poles,
                &mut self.my_dpoles,
                &mut self.my_d2poles,
                &mut self.my_poles2d,
                &mut self.my_dpoles2d,
                &mut self.my_d2poles2d,
                &mut self.my_weigths,
                &mut self.my_dweigths,
                &mut self.my_d2weigths,
            );
            if !b {
                continuity = GeomAbsShape::C1;
            }
        }
        // Checks if myFunc->D1 is implemented
        if continuity == GeomAbsShape::C1 {
            let b = self.my_func.d1(
                first,
                first,
                last,
                &mut self.my_poles,
                &mut self.my_dpoles,
                &mut self.my_poles2d,
                &mut self.my_dpoles2d,
                &mut self.my_weigths,
                &mut self.my_dweigths,
            );
            if !b {
                continuity = GeomAbsShape::C0;
            }
        }

        // So that F was at least 20 times more exact than its approx
        self.my_func.set_tolerance(tol3d_min / 20.0, tol2d / 20.0);

        let nb_interval_c2 = self.my_func.nb_intervals(GeomAbsShape::C2);
        let nb_interval_c3 = self.my_func.nb_intervals(GeomAbsShape::C3);

        if nb_interval_c3 > 1 {
            // (3.1) Approximation with preferential cut
            let mut param_de_decoupe_c2 = vec![0.0f64; nb_interval_c2 + 1];
            self.my_func.intervals(&mut param_de_decoupe_c2, GeomAbsShape::C2);
            let mut param_de_decoupe_c3 = vec![0.0f64; nb_interval_c3 + 1];
            self.my_func.intervals(&mut param_de_decoupe_c3, GeomAbsShape::C3);

            let preferentiel =
                PrefAndRec::with_default_weight(&param_de_decoupe_c2, &param_de_decoupe_c3);
            // OCCT: Approx_SweepApproximation_Eval ev(*this) is passed as
            // TheApproxFunction; the Rust form constructs the evaluator
            // inside Approximation (see the note there).
            self.approximation(
                Some(&one_d_tol),
                two_d_tol.as_deref(),
                &three_d_tol,
                the3d2dtol,
                first,
                last,
                continuity,
                degmax,
                segmax,
                &preferentiel,
            );
        } else {
            // (3.2) Approximation without preferential cut
            let dichotomie = DichoCutting;
            self.approximation(
                Some(&one_d_tol),
                two_d_tol.as_deref(),
                &three_d_tol,
                the3d2dtol,
                first,
                last,
                continuity,
                degmax,
                segmax,
                &dichotomie,
            );
        }
    }

    /// OCCT Approx_SweepApproximation::Approximation (cxx L281-409) — call
    /// F(t) and store the results.
    ///
    /// Architecture note: the OCCT signature carries TheApproxFunction (an
    /// Approx_SweepApproximation_Eval aliasing *this); Rust cannot take a
    /// second &mut to self, so the evaluator is constructed here around the
    /// ApproxAFunction construction (the Perform call sites build `ev` the
    /// same way in OCCT).
    #[allow(clippy::too_many_arguments)]
    fn approximation(
        &mut self,
        one_d_tol: Option<&[f64]>,
        two_d_tol: Option<&[f64]>,
        three_d_tol: &[f64],
        bound_tol: f64,
        first: f64,
        last: f64,
        continuity: GeomAbsShape,
        degmax: i32,
        segmax: i32,
        the_cutting_tool: &dyn Cutting,
    ) {
        let (num1dss, num2dss, num3dss) = (
            self.num1dss as i32,
            self.num2dss as i32,
            self.num3dss as i32,
        );
        // OCCT: AdvApprox_ApproxAFunction Approx(Num1DSS, Num2DSS, Num3DSS,
        // ..., TheApproxFunction, TheCuttingTool).
        let approx = {
            let mut ev = ApproxSweepApproximationEval { tool: self };
            ApproxAFunction::with_cut_tool(
                num1dss,
                num2dss,
                num3dss,
                one_d_tol,
                two_d_tol,
                Some(three_d_tol),
                first,
                last,
                continuity,
                degmax,
                segmax,
                &mut ev,
                the_cutting_tool,
            )
        };
        self.done = approx.has_result();

        if self.done {
            // --> Fill Champs of the surface ----
            self.vdeg = approx.degree() as usize;
            let a_nb_poles = approx.nb_poles();
            // Unfortunately Adv_Approx stores the transposition of the required
            // so, writing tabPoles = Approx.Poles() will give an erroneous result
            // It is only possible to allocate and recopy term by term...
            self.tab_poles = vec![vec![DVec3::ZERO; a_nb_poles]; self.num3dss];
            self.tab_weights = vec![vec![0.0; a_nb_poles]; self.num3dss];

            if self.num1dss == self.num3dss {
                for ii in 1..=self.num3dss {
                    // OCCT: P = Approx.Poles()->Value(jj, ii);
                    // wpoid = Approx.Poles1d()->Value(jj, ii) — the column
                    // accessors give the transposed (jj, ii) reads.
                    let poles_col = approx.poles_flat(ii);
                    let poles1d_col = approx.poles1d_flat(ii);
                    for jj in 1..=a_nb_poles {
                        let mut p = DVec3::new(
                            poles_col[(jj - 1) * 3],
                            poles_col[(jj - 1) * 3 + 1],
                            poles_col[(jj - 1) * 3 + 2],
                        );
                        let wpoid = poles1d_col[jj - 1];
                        p /= wpoid; // It is necessary to divide poles by weight
                        p += self.translation; // OCCT: P.Translate(Translation)
                        self.tab_poles[ii - 1][jj - 1] = p;
                        self.tab_weights[ii - 1][jj - 1] = wpoid;
                    }
                }
            } else {
                for row in self.tab_weights.iter_mut() {
                    for w in row.iter_mut() {
                        *w = 1.0; // OCCT: tabWeights->Init(1)
                    }
                }
                for ii in 1..=self.num3dss {
                    let poles_col = approx.poles_flat(ii);
                    for jj in 1..=a_nb_poles {
                        self.tab_poles[ii - 1][jj - 1] = DVec3::new(
                            poles_col[(jj - 1) * 3],
                            poles_col[(jj - 1) * 3 + 1],
                            poles_col[(jj - 1) * 3 + 2],
                        );
                    }
                }
            }

            // this is better
            self.tab_v_knots = approx.knots_vec().to_vec();
            self.tab_v_mults = approx.multiplicities_vec().to_vec();

            // --> Filling of curves 2D  ----------
            if self.num2dss > 0 {
                self.deg2d = self.vdeg;
                self.tab_2d_knots = approx.knots_vec().to_vec();
                self.tab_2d_mults = approx.multiplicities_vec().to_vec();

                for ii in 1..=self.num2dss {
                    let trsf_inv = self.a_affin[ii - 1].inverted();
                    let mut p2d_col = approx.poles2d_flat(ii);
                    // do not forget to apply inverted homothety.
                    // OCCT: TrsfInv.Transforms(P2d->ChangeValue(jj).ChangeCoord()).
                    for jj in 0..a_nb_poles {
                        let coord = trsf_inv.transforms(DVec2::new(p2d_col[jj * 2], p2d_col[jj * 2 + 1]));
                        p2d_col[jj * 2] = coord.x;
                        p2d_col[jj * 2 + 1] = coord.y;
                    }
                    self.seq_poles2d.push(
                        (0..a_nb_poles)
                            .map(|jj| DVec2::new(p2d_col[jj * 2], p2d_col[jj * 2 + 1]))
                            .collect(),
                    );
                }
            }
            // ---> Filling of errors
            self.m_error3d = vec![0.0; self.num3dss];
            self.a_error3d = vec![0.0; self.num3dss];
            for ii in 1..=self.num3dss {
                self.m_error3d[ii - 1] = approx.max_error_at(3, ii);
                self.a_error3d[ii - 1] = approx.average_error_at(3, ii);
            }

            if self.my_func.is_rational() {
                self.m_error1d = vec![0.0; self.num3dss];
                self.a_error1d = vec![0.0; self.num3dss];
                for ii in 1..=self.num1dss {
                    self.m_error1d[ii - 1] = approx.max_error_at(1, ii);
                    self.a_error1d[ii - 1] = approx.average_error_at(1, ii);
                }
            }

            if self.num2dss > 0 {
                let two_d_tol_v = two_d_tol.expect("TwoDTol");
                self.tab2d_error = vec![0.0; self.num2dss];
                self.ave2d_error = vec![0.0; self.num2dss];
                self.c_on_surf_err = vec![0.0; self.num2dss];
                for ii in 1..=self.num2dss {
                    self.tab2d_error[ii - 1] = approx.max_error_at(2, ii);
                    self.ave2d_error[ii - 1] = approx.average_error_at(2, ii);
                    self.c_on_surf_err[ii - 1] =
                        (self.tab2d_error[ii - 1] / two_d_tol_v[ii - 1]) * bound_tol;
                }
            }
        }
    }

    /// OCCT Approx_SweepApproximation::Eval (cxx L411-433) — the
    /// EvaluatorFunction from AdvApprox; `result` plays the role of the
    /// double& Result (the LocalResult of the OCCT evaluator).
    pub fn eval(
        &mut self,
        parameter: f64,
        derivative_request: i32,
        first: f64,
        last: f64,
        result: &mut [f64],
    ) -> i32 {
        match derivative_request {
            0 => (!self.d0(parameter, first, last, result)) as i32,
            1 => (!self.d1(parameter, first, last, result)) as i32,
            2 => (!self.d2(parameter, first, last, result)) as i32,
            _ => 2,
        }
    }

    /// OCCT Approx_SweepApproximation::D0 (cxx L435-501).
    fn d0(&mut self, param: f64, first: f64, last: f64, result: &mut [f64]) -> bool {
        let mut ok = true;

        // Management of limits
        if (self.first != first) || (last != self.last) {
            self.my_func.set_interval(first, last);
        }

        if (param != self.my_param) || (self.my_order < 0) || (self.first != first) || (last != self.last) {
            // Positioning in case when the last operation is not repeated.
            ok = self.my_func.d0(
                param,
                first,
                last,
                &mut self.my_poles,
                &mut self.my_poles2d,
                &mut self.my_weigths,
            );

            //  poles3d are multiplied by weight after translation.
            for ii in 1..=self.num1dss {
                self.my_poles[ii - 1] -= self.translation;
                self.my_poles[ii - 1] *= self.my_weigths[ii - 1];
            }

            //  The transformation is applied to poles 2d.
            for ii in 1..=self.num2dss {
                self.my_poles2d[ii - 1] = self.a_affin[ii - 1].transforms(self.my_poles2d[ii - 1]);
            }

            // Update variables of control and return
            self.first = first;
            self.last = last;
            self.my_order = 0;
            self.my_param = param;
        }

        //  Extraction of results
        let mut index = 0usize;
        for ii in 1..=self.num1dss {
            result[index] = self.my_weigths[ii - 1];
            index += 1;
        }
        for ii in 1..=self.num2dss {
            result[index] = self.my_poles2d[ii - 1].x;
            result[index + 1] = self.my_poles2d[ii - 1].y;
            index += 2;
        }
        for ii in 1..=self.num3dss {
            result[index] = self.my_poles[ii - 1].x;
            result[index + 1] = self.my_poles[ii - 1].y;
            result[index + 2] = self.my_poles[ii - 1].z;
            index += 3;
        }

        ok
    }

    /// OCCT Approx_SweepApproximation::D1 (cxx L503-583).
    #[allow(clippy::too_many_arguments)]
    fn d1(&mut self, param: f64, first: f64, last: f64, result: &mut [f64]) -> bool {
        let mut ok = true;

        if (self.first != first) || (last != self.last) {
            self.my_func.set_interval(first, last);
        }

        if (param != self.my_param) || (self.my_order < 1) || (self.first != first) || (last != self.last) {
            // Positioning
            ok = self.my_func.d1(
                param,
                first,
                last,
                &mut self.my_poles,
                &mut self.my_dpoles,
                &mut self.my_poles2d,
                &mut self.my_dpoles2d,
                &mut self.my_weigths,
                &mut self.my_dweigths,
            );

            //  Take into account the multiplication of poles3d by weights.
            //  and the translation.
            for ii in 1..=self.num1dss {
                // Translation on the section
                self.my_poles[ii - 1] -= self.translation;
                // Homothety on all.
                let a_weight = self.my_weigths[ii - 1];
                self.my_dpoles[ii - 1] *= a_weight;
                let vaux = self.my_poles[ii - 1];
                self.my_dpoles[ii - 1] += self.my_dweigths[ii - 1] * vaux;
                self.my_poles[ii - 1] *= a_weight; // for the cash
            }

            //  Apply transformation 2d to suitable vectors
            for ii in 1..=self.num2dss {
                let vcoord = self.a_affin[ii - 1].transforms(self.my_dpoles2d[ii - 1]);
                self.my_dpoles2d[ii - 1] = vcoord;
                self.my_poles2d[ii - 1] = self.a_affin[ii - 1].transforms(self.my_poles2d[ii - 1]);
            }

            // Update control variables and return
            self.first = first;
            self.last = last;
            self.my_order = 1;
            self.my_param = param;
        }

        //  Extraction of results
        let mut index = 0usize;
        for ii in 1..=self.num1dss {
            result[index] = self.my_dweigths[ii - 1];
            index += 1;
        }
        for ii in 1..=self.num2dss {
            result[index] = self.my_dpoles2d[ii - 1].x;
            result[index + 1] = self.my_dpoles2d[ii - 1].y;
            index += 2;
        }
        for ii in 1..=self.num3dss {
            result[index] = self.my_dpoles[ii - 1].x;
            result[index + 1] = self.my_dpoles[ii - 1].y;
            result[index + 2] = self.my_dpoles[ii - 1].z;
            index += 3;
        }
        ok
    }

    /// OCCT Approx_SweepApproximation::D2 (cxx L585-678).
    #[allow(clippy::too_many_arguments)]
    fn d2(&mut self, param: f64, first: f64, last: f64, result: &mut [f64]) -> bool {
        let mut ok = true;

        // management of limits
        if (self.first != first) || (last != self.last) {
            self.my_func.set_interval(first, last);
        }

        if (param != self.my_param) || (self.my_order < 2) || (self.first != first) || (last != self.last) {
            // Positioning in case when the last operation is not repeated
            ok = self.my_func.d2(
                param,
                first,
                last,
                &mut self.my_poles,
                &mut self.my_dpoles,
                &mut self.my_d2poles,
                &mut self.my_poles2d,
                &mut self.my_dpoles2d,
                &mut self.my_d2poles2d,
                &mut self.my_weigths,
                &mut self.my_dweigths,
                &mut self.my_d2weigths,
            );

            //  Multiply poles3d by the weight after translations.
            for ii in 1..=self.num1dss {
                // First translate
                self.my_poles[ii - 1] -= self.translation;

                // Calculate the second derivative
                self.my_d2poles[ii - 1] *= self.my_weigths[ii - 1];
                let vaux = self.my_dpoles[ii - 1];
                self.my_d2poles[ii - 1] += (2.0 * self.my_dweigths[ii - 1]) * vaux;
                let vaux = self.my_poles[ii - 1];
                self.my_d2poles[ii - 1] += self.my_d2weigths[ii - 1] * vaux;

                // Then the remainder for the cash
                self.my_dpoles[ii - 1] *= self.my_weigths[ii - 1];
                let vaux = self.my_poles[ii - 1];
                self.my_dpoles[ii - 1] += self.my_dweigths[ii - 1] * vaux;
                self.my_poles[ii - 1] *= self.my_weigths[ii - 1];
            }

            //  Apply transformation to poles 2d.
            for ii in 1..=self.num2dss {
                let vcoord = self.a_affin[ii - 1].transforms(self.my_d2poles2d[ii - 1]);
                self.my_d2poles2d[ii - 1] = vcoord;
                let vcoord = self.a_affin[ii - 1].transforms(self.my_dpoles2d[ii - 1]);
                self.my_dpoles2d[ii - 1] = vcoord;
                self.my_poles2d[ii - 1] = self.a_affin[ii - 1].transforms(self.my_poles2d[ii - 1]);
            }

            // Update variables of control and return
            self.first = first;
            self.last = last;
            self.my_order = 2;
            self.my_param = param;
        }

        //  Extraction of results
        let mut index = 0usize;
        for ii in 1..=self.num1dss {
            result[index] = self.my_d2weigths[ii - 1];
            index += 1;
        }
        for ii in 1..=self.num2dss {
            result[index] = self.my_d2poles2d[ii - 1].x;
            result[index + 1] = self.my_d2poles2d[ii - 1].y;
            index += 2;
        }
        for ii in 1..=self.num3dss {
            result[index] = self.my_d2poles[ii - 1].x;
            result[index + 1] = self.my_d2poles[ii - 1].y;
            result[index + 2] = self.my_d2poles[ii - 1].z;
            index += 3;
        }

        ok
    }

    /// OCCT Approx_SweepApproximation::SurfShape (cxx L680-697).
    #[allow(clippy::too_many_arguments)]
    pub fn surf_shape(
        &self,
        u_degree: &mut usize,
        v_degree: &mut usize,
        nb_u_poles: &mut usize,
        nb_v_poles: &mut usize,
        nb_u_knots: &mut usize,
        nb_v_knots: &mut usize,
    ) {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        *u_degree = self.udeg;
        *v_degree = self.vdeg;
        *nb_u_poles = self.tab_poles.len(); // OCCT: tabPoles->ColLength().
        *nb_v_poles = self.tab_poles.first().map(|r| r.len()).unwrap_or(0); // RowLength().
        *nb_u_knots = self.tab_u_knots.len();
        *nb_v_knots = self.tab_v_knots.len();
    }

    /// OCCT Approx_SweepApproximation::Surface (cxx L699-716).
    #[allow(clippy::too_many_arguments)]
    pub fn surface(
        &self,
        t_poles: &mut [Vec<DVec3>],
        t_weights: &mut [Vec<f64>],
        t_u_knots: &mut [f64],
        t_v_knots: &mut [f64],
        t_u_mults: &mut [i32],
        t_v_mults: &mut [i32],
    ) {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        t_poles.clone_from_slice(&self.tab_poles);
        t_weights.clone_from_slice(&self.tab_weights);
        t_u_knots.copy_from_slice(&self.tab_u_knots);
        t_u_mults.copy_from_slice(&self.tab_u_mults);
        t_v_knots.copy_from_slice(&self.tab_v_knots);
        t_v_mults.copy_from_slice(&self.tab_v_mults);
    }

    /// OCCT Approx_SweepApproximation::MaxErrorOnSurf (cxx L718-753).
    pub fn max_error_on_surf(&self) -> f64 {
        let mut max_error = 0.0f64;
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }

        if self.my_func.is_rational() {
            let mut wmin = vec![0.0f64; self.num1dss];
            self.my_func.get_minimal_weight(&mut wmin);
            let size = self.my_func.maximal_section();
            for ii in 1..=self.num3dss {
                let err = (size * self.m_error1d[ii - 1] + self.m_error3d[ii - 1]) / wmin[ii - 1];
                if err > max_error {
                    max_error = err;
                }
            }
        } else {
            for ii in 1..=self.num3dss {
                let err = self.m_error3d[ii - 1];
                if err > max_error {
                    max_error = err;
                }
            }
        }
        max_error
    }

    /// OCCT Approx_SweepApproximation::AverageErrorOnSurf (cxx L755-784).
    pub fn average_error_on_surf(&self) -> f64 {
        let mut moy_error = 0.0f64;
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }

        if self.my_func.is_rational() {
            let mut wmin = vec![0.0f64; self.num1dss];
            self.my_func.get_minimal_weight(&mut wmin);
            let size = self.my_func.maximal_section();
            for ii in 1..=self.num3dss {
                let err = (size * self.a_error1d[ii - 1] + self.a_error3d[ii - 1]) / wmin[ii - 1];
                moy_error += err;
            }
        } else {
            for ii in 1..=self.num3dss {
                let err = self.a_error3d[ii - 1];
                moy_error += err;
            }
        }
        moy_error / self.num3dss as f64
    }

    /// OCCT Approx_SweepApproximation::Curves2dShape (cxx L786-799).
    pub fn curves2d_shape(&self, degree: &mut usize, nb_poles: &mut usize, nb_knots: &mut usize) {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: Approx_SweepApproximation");
        }
        *degree = self.deg2d;
        *nb_poles = self.seq_poles2d[0].len();
        *nb_knots = self.tab_2d_knots.len();
    }

    /// OCCT Approx_SweepApproximation::Curve2d (cxx L801-817).
    #[allow(clippy::too_many_arguments)]
    pub fn curve2d(
        &self,
        index: usize,
        t_poles: &mut [DVec2],
        t_knots: &mut [f64],
        t_mults: &mut [i32],
    ) {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: Approx_SweepApproximation");
        }
        t_poles.copy_from_slice(&self.seq_poles2d[index - 1]);
        t_knots.copy_from_slice(&self.tab_2d_knots);
        t_mults.copy_from_slice(&self.tab_2d_mults);
    }

    /// OCCT Approx_SweepApproximation::Max2dError (cxx L819-826).
    pub fn max2d_error(&self, index: usize) -> f64 {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.tab2d_error[index - 1]
    }

    /// OCCT Approx_SweepApproximation::Average2dError (cxx L828-835).
    pub fn average2d_error(&self, index: usize) -> f64 {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.ave2d_error[index - 1]
    }

    /// OCCT Approx_SweepApproximation::TolCurveOnSurf (cxx L837-844).
    pub fn tol_curve_on_surf(&self, index: usize) -> f64 {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.c_on_surf_err[index - 1]
    }

    /// OCCT Approx_SweepApproximation::Dump (cxx L846-872).
    pub fn dump(&self, o: &mut dyn std::fmt::Write) {
        let _ = writeln!(o, "Dump of SweepApproximation");
        if self.done {
            let _ = writeln!(o, "Error 3d = {}", self.max_error_on_surf());

            if self.num2dss > 0 {
                let _ = write!(o, "Error 2d = ");
                for ii in 1..=self.num2dss {
                    let _ = write!(o, "{}", self.max2d_error(ii));
                    if ii < self.num2dss {
                        let _ = write!(o, " , \n");
                    }
                }
                println!(); // OCCT: this newline goes to std::cout (cxx L864).
            }
            let _ = writeln!(
                o,
                "{} Segment(s) of degree {}",
                self.tab_v_knots.len() - 1,
                self.vdeg
            );
        } else {
            println!(" Not Done "); // OCCT: std::cout (cxx L870).
        }
    }

    // -----------------------------------------------------------------
    // OCCT Approx_SweepApproximation.lxx (inline accessors)
    // -----------------------------------------------------------------

    /// OCCT IsDone (lxx L26-29).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT UDegree (lxx L31-38).
    pub fn u_degree(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.udeg
    }

    /// OCCT VDegree (lxx L40-47).
    pub fn v_degree(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.vdeg
    }

    /// OCCT SurfPoles (lxx L49-56) — the OCCT const-ref return maps to a
    /// clone (the frozen GeomFill_Pipe consumer builds the BSplineSurface
    /// by value).
    pub fn surf_poles(&self) -> Vec<Vec<DVec3>> {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.tab_poles.clone()
    }

    /// OCCT SurfWeights (lxx L58-65) — cloned (see surf_poles).
    pub fn surf_weights(&self) -> Vec<Vec<f64>> {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.tab_weights.clone()
    }

    /// OCCT SurfUKnots (lxx L67-74).
    pub fn surf_u_knots(&self) -> &[f64] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        &self.tab_u_knots
    }

    /// OCCT SurfVKnots (lxx L76-83).
    pub fn surf_v_knots(&self) -> &[f64] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        &self.tab_v_knots
    }

    /// OCCT SurfUMults (lxx L85-92).
    pub fn surf_u_mults(&self) -> &[i32] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        &self.tab_u_mults
    }

    /// OCCT SurfVMults (lxx L94-101).
    pub fn surf_v_mults(&self) -> &[i32] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        &self.tab_v_mults
    }

    /// OCCT NbCurves2d (lxx L103-110).
    pub fn nb_curves2d(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        self.num2dss
    }

    /// OCCT Curves2dDegree (lxx L112-123).
    pub fn curves2d_degree(&self) -> i32 {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: Approx_SweepApproximation");
        }
        self.deg2d as i32
    }

    /// OCCT Curve2dPoles (lxx L125-137).
    pub fn curve2d_poles(&self, index: usize) -> &[DVec2] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: Approx_SweepApproximation");
        }
        &self.seq_poles2d[index - 1]
    }

    /// OCCT Curves2dKnots (lxx L139-150).
    pub fn curves2d_knots(&self) -> &[f64] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: Approx_SweepApproximation");
        }
        &self.tab_2d_knots
    }

    /// OCCT Curves2dMults (lxx L152-163).
    pub fn curves2d_mults(&self) -> &[i32] {
        if !self.done {
            panic!("StdFail_NotDone: Approx_SweepApproximation");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: Approx_SweepApproximation");
        }
        &self.tab_2d_mults
    }

    // OCCT TolReached(Tol3d, Tol2d) is commented out in
    // Approx_SweepApproximation.lxx (L165-169) — not translated.
}
