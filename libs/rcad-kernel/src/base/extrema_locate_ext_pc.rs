//! OCCT Extrema_LocateExtPC (TKGeomBase/Extrema) — the template-alias chain
//! `Extrema_LocateExtPC = Extrema_GLocateExtPC<Adaptor3d_Curve,
//! Extrema_CurveTool, gp_Pnt, gp_Vec, Extrema_POnCurv,
//! Extrema_ELPCOfLocateExtPC, Extrema_LocEPCOfLocateExtPC>`
//! (Extrema_LocateExtPC.hxx L28-36), consumed by ChFi3d_Builder_2.cxx
//! L43/L338/L3825 and IntTools_Context.cxx L829/L890.  Template-parameter
//! mapping: TheELPC = Extrema_ELPCOfLocateExtPC = Extrema_GGExtPC (the kernel
//! [`ExtremaExtPC`](crate::base::extrema_ext_pc::ExtremaExtPC)); TheLocEPC =
//! Extrema_LocEPCOfLocateExtPC = Extrema_GenLocateExtPC over ThePCLocF =
//! Extrema_PCLocFOfLocEPCOfLocateExtPC = Extrema_GFuncExtPC (the kernel
//! [`FuncExtPC`](crate::base::extrema_func_ext_pc::FuncExtPC)); the bounded
//! math_FunctionRoot ctor maps to the kernel
//! [`NewtonFunctionRoot::new_bounded`](crate::math::newton_function_root).

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema_ext_pc::{ExtPCurveTool, ExtremaExtPC};
use crate::base::extrema_func_ext_pc::FuncExtPC;
use crate::base::proj_lib::CurveType;
use crate::core::precision::{is_infinite_value, CONFUSION};
use crate::math::newton_function_root::NewtonFunctionRoot;
use crate::math::root::FunctionValue;
use crate::math::GeomAbsShape;

/// OCCT Extrema_GenLocateExtPC (Extrema_GenLocateExtPC.hxx L44-168); the
/// `Extrema_LocEPCOfLocateExtPC` alias.
pub struct GenLocateExtPC<'a> {
    /// hxx L163: bool myDone.
    my_done: bool,
    /// hxx L164: double mytolU.
    my_tol_u: f64,
    /// hxx L165: double myumin.
    my_umin: f64,
    /// hxx L166: double myusup.
    my_usup: f64,
    /// hxx L167: ThePCLocF myF.
    my_f: FuncExtPC<'a>,
}

impl<'a> GenLocateExtPC<'a> {
    /// OCCT Extrema_GenLocateExtPC() default constructor (hxx L49-55).
    pub fn new() -> Self {
        GenLocateExtPC {
            my_done: false,
            my_tol_u: 0.0,
            my_umin: 0.0,
            my_usup: 0.0,
            my_f: FuncExtPC::new(),
        }
    }

    /// OCCT Initialize(theC, theUmin, theUsup, theTolU) (hxx L118-126).
    pub fn initialize(
        &mut self,
        the_c: &'a dyn ExtPCurveTool,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
    ) {
        // hxx L120.
        self.my_done = false;
        self.my_f.initialize(the_c);
        self.my_umin = the_umin;
        self.my_usup = the_usup;
        self.my_tol_u = the_tol_u;
    }

    /// OCCT Perform(theP, theU0) (hxx L131-155).
    pub fn perform(&mut self, the_p: DVec3, the_u0: f64) {
        // hxx L135.
        self.my_f.set_point(the_p);
        // hxx L138: math_FunctionRoot S(myF, theU0, mytolU, myumin, myusup)
        // — the OCCT default NbIterations is 100.
        let s = NewtonFunctionRoot::new_bounded(
            &mut self.my_f,
            the_u0,
            self.my_tol_u,
            self.my_umin,
            self.my_usup,
            100,
        );
        // hxx L139.
        self.my_done = s.is_done();
        if self.my_done {
            // hxx L140-147.
            let uu = self.point().param;
            match self.my_f.value(uu) {
                Some(ff) => {
                    if ff.abs() >= 1.0e-07 {
                        self.my_done = false;
                    }
                }
                None => self.my_done = false,
            }
        }
    }

    /// OCCT IsDone() (hxx L158).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT SquareDistance() (hxx L161-168) — myF.SquareDistance(1).
    pub fn square_distance(&self) -> f64 {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_f.square_distance(1)
    }

    /// OCCT IsMin() (hxx L171-178) — myF.IsMin(1).
    pub fn is_min(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_f.is_min(1)
    }

    /// OCCT Point() (hxx L181-188) — myF.Point(1).
    pub fn point(&self) -> &POnCurve {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_f.point(1)
    }
}

/// OCCT Extrema_GLocateExtPC (Extrema_GLocateExtPC.hxx L54-381); the
/// `Extrema_LocateExtPC` alias.
pub struct LocateExtPC<'a> {
    /// hxx L374: ThePOnC mypp.
    my_pp: POnCurve,
    /// hxx L375: TheCurve* myC.
    my_c: Option<&'a dyn ExtPCurveTool>,
    /// hxx L376: double mydist2.
    my_dist2: f64,
    /// hxx L377: bool myismin.
    my_ismin: bool,
    /// hxx L378: bool myDone.
    my_done: bool,
    /// hxx L379: double myumin.
    my_umin: f64,
    /// hxx L380: double myusup.
    my_usup: f64,
    /// hxx L381: double mytol.
    my_tol: f64,
    /// hxx L382: TheLocEPC myLocExtPC.
    my_loc_ext_pc: GenLocateExtPC<'a>,
    /// hxx L383: TheELPC myExtremPC.
    my_extrem_pc: ExtremaExtPC<'a>,
    /// hxx L384: GeomAbs_CurveType type.
    type_: CurveType,
    /// hxx L385: int numberext.
    numberext: usize,
}

impl<'a> LocateExtPC<'a> {
    /// OCCT Extrema_GLocateExtPC() default constructor (hxx L59-73).
    pub fn new() -> Self {
        LocateExtPC {
            my_pp: POnCurve {
                param: 0.0,
                point: DVec3::ZERO,
            },
            my_c: None,
            my_dist2: 0.0,
            my_ismin: false,
            my_done: false,
            my_umin: 0.0,
            my_usup: 0.0,
            my_tol: 0.0,
            my_loc_ext_pc: GenLocateExtPC::new(),
            my_extrem_pc: ExtremaExtPC::new(),
            type_: CurveType::Other,
            numberext: 0,
        }
    }

    /// OCCT Extrema_GLocateExtPC(theP, theC, theU0, theTolF) (hxx L82-92):
    /// Initialize over the full curve range, then Perform.
    pub fn new_point_curve_seed(
        the_p: DVec3,
        the_c: &'a dyn ExtPCurveTool,
        the_u0: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut this = LocateExtPC::new();
        this.initialize(
            the_c,
            the_c.first_parameter(),
            the_c.last_parameter(),
            the_tol_f,
        );
        this.perform(the_p, the_u0);
        this
    }

    /// OCCT Extrema_GLocateExtPC(theP, theC, theU0, theUmin, theUsup,
    /// theTolF) (hxx L94-101): Initialize over [theUmin, theUsup], then
    /// Perform.
    pub fn new_point_curve_seed_ranged(
        the_p: DVec3,
        the_c: &'a dyn ExtPCurveTool,
        the_u0: f64,
        the_umin: f64,
        the_usup: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut this = LocateExtPC::new();
        this.initialize(the_c, the_umin, the_usup, the_tol_f);
        this.perform(the_p, the_u0);
        this
    }

    /// OCCT Initialize(theC, theUmin, theUsup, theTolF) (hxx L99-119).
    pub fn initialize(
        &mut self,
        the_c: &'a dyn ExtPCurveTool,
        the_umin: f64,
        the_usup: f64,
        the_tol_f: f64,
    ) {
        // hxx L101-107.
        self.my_c = Some(the_c);
        self.my_tol = the_tol_f;
        self.my_umin = the_umin;
        self.my_usup = the_usup;
        self.type_ = the_c.get_type();
        // hxx L109: double tolu = TheCurveTool::Resolution(theC,
        // Precision::Confusion()).
        let tolu = the_c.resolution(CONFUSION);
        // hxx L110-118: the BSpline/Bezier/Offset/Other types drive the local
        // Newton locator; the analytic types drive the full Extrema_ExtPC.
        // (the rcad CurveType enum carries no offset value — an OCCT offset
        // curve reports Other at the rcad adaptor layer.)
        match self.type_ {
            CurveType::BSpline | CurveType::Bezier | CurveType::Other => {
                self.my_loc_ext_pc.initialize(the_c, the_umin, the_usup, tolu);
            }
            _ => {
                self.my_extrem_pc.initialize(the_c, the_umin, the_usup, tolu);
            }
        }
    }

    /// OCCT Perform(theP, theU0) (hxx L124-345).
    pub fn perform(&mut self, the_p: DVec3, the_u0: f64) {
        // hxx L125-130.
        let mut local_u0 = the_u0;
        match self.type_ {
            // hxx L132-321: the Other/Offset/BSpline arm — the search runs
            // interval by continuous-C2 interval.
            CurveType::Other | CurveType::BSpline => {
                let c = self.my_c.expect("Extrema_GLocateExtPC: no curve");
                // hxx L134-137: the C2 interval array (1, n + 1).
                let n = c.nb_intervals(GeomAbsShape::C2) as i32;
                let the_inter = c.intervals(GeomAbsShape::C2);
                // hxx L139-145: be gentle with the caller.
                if local_u0 < self.my_umin {
                    local_u0 = self.my_umin;
                } else if local_u0 > self.my_usup {
                    local_u0 = self.my_usup;
                }
                // hxx L129: double myintuinf = 0, myintusup = 0 — declared
                // outside the search loop; the values are loop-carried and
                // reused after the loop (the last tried interval when U0 is
                // not found inside any interval).
                let mut my_int_u_inf = 0.0;
                let mut my_int_u_sup = 0.0;
                // hxx L147-158: search for the interval containing U0.
                let mut found = false;
                let mut inter: i32 = 1;
                while !found && inter <= n {
                    my_int_u_inf = the_inter[(inter - 1) as usize].max(self.my_umin);
                    my_int_u_sup = the_inter[inter as usize].min(self.my_usup);
                    if (local_u0 >= my_int_u_inf) && (local_u0 < my_int_u_sup) {
                        found = true;
                    }
                    inter += 1;
                }
                // hxx L160-161: IFV 16.06.00 — inter is increased after found.
                if found {
                    inter -= 1;
                }
                // hxx L164-171: try on the found interval — the OCCT code
                // reuses the loop-carried myintuinf/myintusup here.
                self.my_loc_ext_pc
                    .initialize(c, my_int_u_inf, my_int_u_sup, self.my_tol);
                self.my_loc_ext_pc.perform(the_p, local_u0);
                self.my_done = self.my_loc_ext_pc.is_done();
                if self.my_done {
                    // hxx L168-170.
                    self.my_pp = self.my_loc_ext_pc.point().clone();
                    self.my_ismin = self.my_loc_ext_pc.is_min();
                    self.my_dist2 = self.my_loc_ext_pc.square_distance();
                } else {
                    // hxx L173-312: try on the neighboring intervals.
                    let mut k: i32 = 1;
                    let mut i1: i32 = inter;
                    let mut i2: i32 = inter;
                    let mut s1_inf;
                    let mut s2_inf;
                    let mut s1_sup;
                    let mut s2_sup;
                    // hxx L181-186.
                    let (p1, v1) = c.d1(my_int_u_inf);
                    s2_inf = (p1 - the_p).dot(v1);
                    let (p1, v1) = c.d1(my_int_u_sup);
                    s1_sup = (p1 - the_p).dot(v1);

                    while !self.my_done && i2 > 0 && i1 <= n {
                        i1 = inter + k;
                        i2 = inter - k;
                        // hxx L196-224: the interval above.
                        if i1 <= n {
                            let my_int_u_inf = the_inter[(i1 - 1) as usize].max(self.my_umin);
                            let my_int_u_sup = the_inter[i1 as usize].min(self.my_usup);
                            if my_int_u_inf < my_int_u_sup {
                                // hxx L201-203.
                                let (p1, v1) = c.d1(my_int_u_inf);
                                s2_sup = (p1 - the_p).dot(v1);
                                if is_infinite_value(s2_sup) || is_infinite_value(s1_sup) {
                                    break;
                                }
                                // hxx L206-214: s1sup * s2sup <= RealEpsilon —
                                // the extremum on the boundary.
                                if s1_sup * s2_sup <= f64::EPSILON {
                                    self.my_done = true;
                                    self.my_pp =
                                        POnCurve { param: my_int_u_inf, point: p1 };
                                    self.my_ismin = s1_sup <= 0.0;
                                    self.my_dist2 = the_p.distance_squared(p1);
                                    break;
                                }
                                // hxx L216-217.
                                let (p1, v1) = c.d1(my_int_u_sup);
                                s1_sup = (p1 - the_p).dot(v1);
                                // hxx L218-224.
                                self.my_loc_ext_pc
                                    .initialize(c, my_int_u_inf, my_int_u_sup, self.my_tol);
                                self.my_loc_ext_pc
                                    .perform(the_p, (my_int_u_inf + my_int_u_sup) * 0.5);
                                self.my_done = self.my_loc_ext_pc.is_done();
                                if self.my_done {
                                    self.my_pp = self.my_loc_ext_pc.point().clone();
                                    self.my_ismin = self.my_loc_ext_pc.is_min();
                                    self.my_dist2 = self.my_loc_ext_pc.square_distance();
                                    break;
                                }
                            }
                        }
                        // hxx L226-306: the interval below.
                        if i2 > 0 {
                            let my_int_u_inf = the_inter[(i2 - 1) as usize].max(self.my_umin);
                            let my_int_u_sup = the_inter[i2 as usize].min(self.my_usup);
                            if my_int_u_inf < my_int_u_sup {
                                // hxx L232-234.
                                let (p1, v1) = c.d1(my_int_u_sup);
                                s1_inf = (p1 - the_p).dot(v1);
                                if is_infinite_value(s2_inf) || is_infinite_value(s1_inf) {
                                    break;
                                }
                                // hxx L236-246.
                                if s1_inf * s2_inf <= f64::EPSILON {
                                    self.my_done = true;
                                    self.my_pp =
                                        POnCurve { param: my_int_u_sup, point: p1 };
                                    self.my_ismin = s1_inf <= 0.0;
                                    self.my_dist2 = the_p.distance_squared(p1);
                                    break;
                                }
                                // hxx L248-249.
                                let (p1, v1) = c.d1(my_int_u_inf);
                                s2_inf = (p1 - the_p).dot(v1);
                                // hxx L250-256.
                                self.my_loc_ext_pc
                                    .initialize(c, my_int_u_inf, my_int_u_sup, self.my_tol);
                                self.my_loc_ext_pc
                                    .perform(the_p, (my_int_u_inf + my_int_u_sup) * 0.5);
                                self.my_done = self.my_loc_ext_pc.is_done();
                                if self.my_done {
                                    self.my_pp = self.my_loc_ext_pc.point().clone();
                                    self.my_ismin = self.my_loc_ext_pc.is_min();
                                    self.my_dist2 = self.my_loc_ext_pc.square_distance();
                                    break;
                                }
                            }
                        }
                        // hxx L311.
                        k += 1;
                    }
                }
            }

            // hxx L323-329: the Bezier arm.
            CurveType::Bezier => {
                self.my_loc_ext_pc.perform(the_p, the_u0);
                self.my_done = self.my_loc_ext_pc.is_done();
            }

            // hxx L331-344: the default arm — the full Extrema_ExtPC, the
            // extremum nearest U0 wins.
            _ => {
                self.my_extrem_pc.perform(the_p);
                self.numberext = 0;
                if self.my_extrem_pc.is_done() {
                    let mut val_u2 = f64::MAX; // RealLast()
                    for i in 1..=self.my_extrem_pc.nb_ext() {
                        let par = self.my_extrem_pc.point(i).param;
                        let val_u = (par - the_u0).abs();
                        if val_u <= val_u2 {
                            val_u2 = val_u;
                            self.numberext = i;
                            self.my_done = true;
                        }
                    }
                }
                if self.numberext == 0 {
                    self.my_done = false;
                }
            }
        }
    }

    /// OCCT IsDone() (hxx L347).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT SquareDistance() (hxx L350-372).
    pub fn square_distance(&self) -> f64 {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        // hxx L355-370: the type dispatch.
        if self.type_ == CurveType::Bezier {
            self.my_loc_ext_pc.square_distance()
        } else if self.type_ == CurveType::BSpline || self.type_ == CurveType::Other {
            self.my_dist2
        } else if self.numberext != 0 {
            self.my_extrem_pc.square_distance(self.numberext)
        } else {
            0.0
        }
    }

    /// OCCT IsMin() (hxx L375-397 in the source order of the getters).
    pub fn is_min(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        if self.type_ == CurveType::Bezier {
            self.my_loc_ext_pc.is_min()
        } else if self.type_ == CurveType::BSpline || self.type_ == CurveType::Other {
            self.my_ismin
        } else if self.numberext != 0 {
            self.my_extrem_pc.is_min(self.numberext)
        } else {
            false
        }
    }

    /// OCCT Point() (hxx L400-410).
    pub fn point(&self) -> &POnCurve {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        if self.type_ == CurveType::Bezier {
            self.my_loc_ext_pc.point()
        } else if self.type_ == CurveType::BSpline || self.type_ == CurveType::Other {
            &self.my_pp
        } else {
            self.my_extrem_pc.point(self.numberext)
        }
    }
}
