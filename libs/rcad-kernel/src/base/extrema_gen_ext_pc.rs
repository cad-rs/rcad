//! OCCT Extrema_GGenExtPC (TKGeomBase/Extrema/Extrema_GGenExtPC.hxx L37-234)
//! — generic search of all parameter values u where the distance function
//! F(u) = distance(P, C(u)) has an extremum (dF/du = 0), driven by
//! math_FunctionRoots.
//!
//! OCCT template instantiation (Extrema_EPCOfExtPC.hxx L28-32):
//! `Extrema_EPCOfExtPC = Extrema_GGenExtPC<Adaptor3d_Curve, Extrema_CurveTool,
//! Extrema_POnCurv, gp_Pnt, Extrema_PCFOfEPCOfExtPC>`, so this single
//! translation unit is both the `Extrema_GGenExtPC` body and the
//! `Extrema_EPCOfExtPC` alias.
//!
//! Template-parameter mapping: `TheCurve` + `TheTool` map to the
//! [`ExtPCurveTool`] trait (see `extrema_ext_pc`); `ThePOnC` maps to
//! [`POnCurve`]; `ThePoint` maps to `glam::DVec3`; `ThePCF` maps to
//! [`FuncExtPC`](crate::base::extrema_func_ext_pc::FuncExtPC).

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema_ext_pc::ExtPCurveTool;
use crate::base::extrema_func_ext_pc::FuncExtPC;
use crate::math::root::FunctionRoots;

/// OCCT Extrema_GGenExtPC (hxx L38-234); the `Extrema_EPCOfExtPC` alias.
pub struct GenExtPC<'a> {
    /// hxx L226: bool myDone.
    my_done: bool,
    /// hxx L227: bool myInit — written by the OCCT body, never read there.
    #[allow(dead_code)]
    my_init: bool,
    /// hxx L228: int mynbsample.
    my_nb_sample: i32,
    /// hxx L229: double myumin.
    my_umin: f64,
    /// hxx L230: double myusup.
    my_usup: f64,
    /// hxx L231: double mytolu.
    my_tol_u: f64,
    /// hxx L232: double mytolF.
    my_tol_f: f64,
    /// hxx L233: ThePCF myF.
    my_f: FuncExtPC<'a>,
}

impl<'a> GenExtPC<'a> {
    /// OCCT Extrema_GGenExtPC() default constructor (hxx L44-53).
    pub fn new() -> Self {
        GenExtPC {
            my_done: false,
            my_init: false,
            my_nb_sample: 0,
            my_umin: 0.0,
            my_usup: 0.0,
            my_tol_u: 0.0,
            my_tol_f: 0.0,
            my_f: FuncExtPC::new(),
        }
    }

    /// OCCT Extrema_GGenExtPC(theP, theC, theNbSample, theTolU, theTolF)
    /// (hxx L61-70).
    pub fn new_point_curve(
        the_p: DVec3,
        the_c: &'a dyn ExtPCurveTool,
        the_nb_sample: i32,
        the_tol_u: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut this = GenExtPC {
            my_done: false,
            my_init: false,
            my_nb_sample: 0,
            my_umin: 0.0,
            my_usup: 0.0,
            my_tol_u: 0.0,
            my_tol_f: 0.0,
            my_f: FuncExtPC::new_point_curve(the_p, the_c),
        };
        this.initialize_full(the_c, the_nb_sample, the_tol_u, the_tol_f);
        this.perform(the_p);
        this
    }

    /// OCCT Extrema_GGenExtPC(theP, theC, theNbSample, theUmin, theUsup,
    /// theTolU, theTolF) (hxx L80-91).
    pub fn new_point_curve_ranged(
        the_p: DVec3,
        the_c: &'a dyn ExtPCurveTool,
        the_nb_sample: i32,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut this = GenExtPC {
            my_done: false,
            my_init: false,
            my_nb_sample: 0,
            my_umin: 0.0,
            my_usup: 0.0,
            my_tol_u: 0.0,
            my_tol_f: 0.0,
            my_f: FuncExtPC::new_point_curve(the_p, the_c),
        };
        this.initialize_ranged(the_c, the_nb_sample, the_umin, the_usup, the_tol_u, the_tol_f);
        this.perform(the_p);
        this
    }

    /// OCCT Initialize(theC, theNbU, theTolU, theTolF) (hxx L98-110).
    pub fn initialize_full(
        &mut self,
        the_c: &'a dyn ExtPCurveTool,
        the_nb_u: i32,
        the_tol_u: f64,
        the_tol_f: f64,
    ) {
        self.my_init = true;
        self.my_nb_sample = the_nb_u;
        self.my_tol_u = the_tol_u;
        self.my_tol_f = the_tol_f;
        self.my_f.initialize(the_c);
        self.my_umin = the_c.first_parameter();
        self.my_usup = the_c.last_parameter();
    }

    /// OCCT Initialize(theC, theNbU, theUmin, theUsup, theTolU, theTolF)
    /// (hxx L119-133).
    pub fn initialize_ranged(
        &mut self,
        the_c: &'a dyn ExtPCurveTool,
        the_nb_u: i32,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
        the_tol_f: f64,
    ) {
        self.my_init = true;
        self.my_nb_sample = the_nb_u;
        self.my_tol_u = the_tol_u;
        self.my_tol_f = the_tol_f;
        self.my_f.initialize(the_c);
        self.my_umin = the_umin;
        self.my_usup = the_usup;
    }

    /// OCCT Initialize(theNbU, theUmin, theUsup, theTolU, theTolF)
    /// (hxx L141-152).
    pub fn initialize_interval(
        &mut self,
        the_nb_u: i32,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
        the_tol_f: f64,
    ) {
        self.my_nb_sample = the_nb_u;
        self.my_tol_u = the_tol_u;
        self.my_tol_f = the_tol_f;
        self.my_umin = the_umin;
        self.my_usup = the_usup;
    }

    /// OCCT Initialize(theC) (hxx L156).
    pub fn initialize_curve(&mut self, the_c: &'a dyn ExtPCurveTool) {
        self.my_f.initialize(the_c);
    }

    /// OCCT Perform(theP) (hxx L160-173).
    pub fn perform(&mut self, the_p: DVec3) {
        self.my_f.set_point(the_p);
        self.my_f.sub_interval_initialize(self.my_umin, self.my_usup);
        self.my_done = false;

        // hxx L166: math_FunctionRoots S(myF, myumin, myusup, mynbsample,
        // mytolu, mytolF, mytolF) — the trailing K keeps the OCCT default 0.
        let s = FunctionRoots::new(
            &mut self.my_f,
            self.my_umin,
            self.my_usup,
            self.my_nb_sample,
            self.my_tol_u,
            self.my_tol_f,
            self.my_tol_f,
            0.0,
        );
        if !s.is_done() || s.is_all_null() {
            return;
        }

        self.my_done = true;
    }

    /// OCCT IsDone() (hxx L176).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT NbExt() (hxx L180-187).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_f.nb_ext()
    }

    /// OCCT SquareDistance(theN) (hxx L192-199) — 1-based.
    pub fn square_distance(&self, the_n: usize) -> f64 {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_f.square_distance(the_n)
    }

    /// OCCT IsMin(theN) (hxx L204-211) — 1-based.
    pub fn is_min(&self, the_n: usize) -> bool {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_f.is_min(the_n)
    }

    /// OCCT Point(theN) (hxx L216-223) — 1-based.
    pub fn point(&self, the_n: usize) -> &POnCurve {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_f.point(the_n)
    }
}
