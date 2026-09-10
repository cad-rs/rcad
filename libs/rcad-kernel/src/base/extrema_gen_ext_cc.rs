//! OCCT Extrema_GGenExtCC (TKGeomBase/Extrema/Extrema_GGenExtCC.hxx L55-920)
//! — the template class that computes the extremal distances between two
//! curves with Evtushenko's global optimization solver.  `Extrema_ECC`
//! (Extrema_ECC.hxx L28-34) is the alias
//! `Extrema_GGenExtCC<Adaptor3d_Curve, Extrema_CurveTool, Adaptor3d_Curve,
//! Extrema_CurveTool, Extrema_POnCurv, gp_Pnt, Extrema_ExtPC>`, so this single
//! translation unit is also the `Extrema_ECC` body.
//!
//! Also carries the file-local helpers:
//! - `Extrema_GGenExtCC_comp` (hxx L135-150) — the gp_XY comparator,
//! - `Extrema_GGenExtCC_ChangeIntervals` (hxx L152-202) — interval refinement,
//! - `Extrema_GGenExtCC_PointsInspector` (hxx L204-246),
//! - `Extrema_GGenExtCC_ProjPOnC` (hxx L248-263).
//!
//! Template-parameter mapping: `TheCurve1` / `TheCurve2` (Adaptor3d_Curve) and
//! `TheCurveTool1` / `TheCurveTool2` (Extrema_CurveTool, the static facade over
//! that curve) both map to the `ExtremaCurveTool` trait; `ThePOnC` maps to
//! `POnCurve`; `TheExtPC` maps to `crate::base::extrema::ExtPC`.
//!
//! rcad note: OCCT lets `aFunc` (the objective) and `aFinder` (the solver)
//! coexist because `math_GlobOptMin` holds a raw pointer; Rust forbids the
//! aliasing, so the objective is re-created after the solver is dropped.  The
//! objective is stateless, so the call graph is unchanged.

use glam::DVec2;
use glam::DVec3;

use crate::base::extrema::ExtPC;
use crate::base::extrema::POnCurve;
use crate::base::extrema_curve_tool::ExtremaCurveTool;
use crate::base::extrema_glob_opt_func_cc::GlobOptFuncCCC2;
use crate::base::extrema_math_glob_opt_min::{MathGlobOptMin, MathMultiVarFunc};
use crate::base::extrema_math_opt::{CellFilter, CellFilterAction};
use crate::base::proj_lib::CurveType;
use crate::core::precision::{is_infinite_value, PCONFUSION, CONFUSION};
use crate::math::math_bfgs::{MultipleVarFunction, MultipleVarFunctionWithGradient};
use crate::math::math_matrix::Vector;
use crate::math::GeomAbsShape;

/// OCCT M_SQRT2 (Extrema_GGenExtCC.hxx L40-42).
const M_SQRT2: f64 = 1.41421356237309504880168872420969808;

/// OCCT RealLast().
const REAL_LAST: f64 = f64::MAX;

/// OCCT Extrema_GGenExtCC_comp (hxx L135-150) — comparator used in std::sort.
fn ggen_ext_cc_comp(the_a: &DVec2, the_b: &DVec2) -> std::cmp::Ordering {
    if the_a.x < the_b.x {
        std::cmp::Ordering::Less
    } else if the_a.x == the_b.x {
        if the_a.y < the_b.y {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    } else {
        std::cmp::Ordering::Greater
    }
}

/// OCCT Extrema_GGenExtCC_ChangeIntervals(theInts, theNbInts) (hxx L152-202).
///
/// OCCT indexes the array 1..Length(); rcad keeps it 0-based, so every OCCT
/// `Value(i)` maps to `the_ints[i - 1]`.
fn ggen_ext_cc_change_intervals(the_ints: &mut Vec<f64>, the_nb_ints: usize) {
    let a_nb_ints = the_ints.len() - 1;
    let mut a_nb_add = the_nb_ints as i64 - a_nb_ints as i64;
    let mut a_new_ints = vec![0.0f64; the_nb_ints + 1];
    let a_nb_last = the_ints.len();
    if a_nb_ints == 1 {
        // hxx L161-173.
        a_new_ints[0] = the_ints[0];
        a_new_ints[the_nb_ints] = the_ints[the_ints.len() - 1];
        let dt = (the_ints[the_ints.len() - 1] - the_ints[0]) / the_nb_ints as f64;
        let mut t = the_ints[0] + dt;
        let mut i = 2usize;
        while i <= the_nb_ints {
            a_new_ints[i - 1] = t;
            t += dt;
            i += 1;
        }
        *the_ints = a_new_ints;
        return;
    }
    // hxx L174-200.
    for i in 1..=a_nb_last {
        a_new_ints[i - 1] = the_ints[i - 1];
    }
    let mut a_nb_last = a_nb_last;
    while a_nb_add > 0 {
        let mut an_l_int_max = -1.0;
        let mut a_max_ind: i64 = -1;
        for i in 1..a_nb_last {
            let an_l = a_new_ints[i] - a_new_ints[i - 1];
            if an_l > an_l_int_max {
                an_l_int_max = an_l;
                a_max_ind = i as i64;
            }
        }

        let t = (a_new_ints[a_max_ind as usize] + a_new_ints[a_max_ind as usize - 1]) / 2.0;
        let mut i = a_nb_last;
        while i > a_max_ind as usize {
            a_new_ints[i] = a_new_ints[i - 1];
            i -= 1;
        }
        a_nb_last += 1;
        a_nb_add -= 1;
        a_new_ints[a_max_ind as usize] = t;
    }
    *the_ints = a_new_ints;
}

/// OCCT Extrema_GGenExtCC_PointsInspector (hxx L204-246).
pub struct GGenCcPointsInspector {
    /// hxx L207: static constexpr int Dimension = 2.
    dimension: usize,
    /// hxx L243: double myTol.
    my_tol: f64,
    /// hxx L244: gp_XY myCurrent.
    my_current: DVec2,
    /// hxx L245: bool myIsFind.
    my_is_find: bool,
}

impl GGenCcPointsInspector {
    /// OCCT Extrema_GGenExtCC_PointsInspector(theTol) (hxx L219-223).
    pub fn new(the_tol: f64) -> Self {
        GGenCcPointsInspector {
            dimension: 2,
            my_tol: the_tol * the_tol,
            my_current: DVec2::ZERO,
            my_is_find: false,
        }
    }

    /// OCCT Coord(i, thePnt) (hxx L212).
    pub fn coord(i: usize, the_pnt: &DVec2) -> f64 {
        match i {
            0 => the_pnt.x,
            _ => the_pnt.y,
        }
    }

    /// OCCT Shift(thePnt, theTol) (hxx L214-217).
    pub fn shift(&self, the_pnt: &DVec2, the_tol: f64) -> DVec2 {
        DVec2::new(the_pnt.x + the_tol, the_pnt.y + the_tol)
    }

    /// OCCT ClearFind() (hxx L225).
    pub fn clear_find(&mut self) {
        self.my_is_find = false;
    }

    /// OCCT isFind() (hxx L227).
    pub fn is_find(&self) -> bool {
        self.my_is_find
    }

    /// OCCT SetCurrent(theCurPnt) (hxx L229).
    pub fn set_current(&mut self, the_cur_pnt: &DVec2) {
        self.my_current = *the_cur_pnt;
    }

    /// OCCT Inspect(theObject) (hxx L231-240).
    pub fn inspect(&mut self, the_object: &DVec2) -> CellFilterAction {
        let a_pt = self.my_current - *the_object;
        let a_sq_dist = a_pt.length_squared();
        if a_sq_dist < self.my_tol {
            self.my_is_find = true;
        }
        CellFilterAction::Keep
    }

    /// OCCT `Dimension`.
    pub fn dimension(&self) -> usize {
        self.dimension
    }
}

/// OCCT Extrema_GGenExtCC_ProjPOnC(theP, theProjTool) (hxx L248-263).
///
/// The rcad `ExtPC` takes the projection range per `Perform`, so the range the
/// OCCT `Extrema_ExtPC::Initialize` fixed is passed through.
fn ggen_ext_cc_proj_p_on_c<F>(
    the_p: DVec3,
    the_proj_tool: &mut ExtPC,
    eval: &F,
    uinf: f64,
    usup: f64,
) -> f64
where
    F: Fn(f64) -> DVec3,
{
    let mut a_dist = REAL_LAST;
    the_proj_tool.perform_fn(the_p, eval, uinf, usup);
    if the_proj_tool.is_done() && the_proj_tool.nb_ext() > 0 {
        for i in 1..=the_proj_tool.nb_ext() {
            let a_d = the_proj_tool.square_distance(i);
            if a_d < a_dist {
                a_dist = a_d;
            }
        }
    }
    a_dist
}

/// OCCT `Extrema_ECC` = `Extrema_GGenExtCC<...>` (Extrema_GGenExtCC.hxx
/// L55-126).
pub struct GenExtCC<'a> {
    /// hxx L117: bool myIsFindSingleSolution.
    my_is_find_single_solution: bool,
    /// hxx L118: bool myParallel.
    my_parallel: bool,
    /// hxx L119: double myCurveMinTol.
    my_curve_min_tol: f64,
    /// hxx L120: math_Vector myLowBorder.
    my_low_border: Vector,
    /// hxx L121: math_Vector myUppBorder.
    my_upp_border: Vector,
    /// hxx L122: NCollection_Sequence<double> myPoints1.
    my_points1: Vec<f64>,
    /// hxx L123: NCollection_Sequence<double> myPoints2.
    my_points2: Vec<f64>,
    /// hxx L124: void* myC[2] — the two curves.
    my_c1: &'a dyn ExtremaCurveTool,
    my_c2: &'a dyn ExtremaCurveTool,
    /// hxx L125: bool myDone.
    my_done: bool,
}

impl<'a> GenExtCC<'a> {
    /// OCCT Extrema_GGenExtCC(C1, C2) (hxx L300-315).
    pub fn new(c1: &'a dyn ExtremaCurveTool, c2: &'a dyn ExtremaCurveTool) -> Self {
        let mut this = GenExtCC {
            my_is_find_single_solution: false,
            my_parallel: false,
            my_curve_min_tol: PCONFUSION,
            my_low_border: Vector::new(1, 2),
            my_upp_border: Vector::new(1, 2),
            my_points1: Vec::new(),
            my_points2: Vec::new(),
            my_c1: c1,
            my_c2: c2,
            my_done: false,
        };
        this.my_low_border.set(1, c1.first_parameter());
        this.my_low_border.set(2, c2.first_parameter());
        this.my_upp_border.set(1, c1.last_parameter());
        this.my_upp_border.set(2, c2.last_parameter());
        this
    }

    /// OCCT Extrema_GGenExtCC(C1, C2, Uinf, Usup, Vinf, Vsup) (hxx L326-346).
    pub fn with_params(
        c1: &'a dyn ExtremaCurveTool,
        c2: &'a dyn ExtremaCurveTool,
        the_uinf: f64,
        the_usup: f64,
        the_vinf: f64,
        the_vsup: f64,
    ) -> Self {
        let mut this = GenExtCC::new(c1, c2);
        this.set_params(c1, c2, the_uinf, the_usup, the_vinf, the_vsup);
        this
    }

    /// OCCT SetParams (hxx L357-376).
    pub fn set_params(
        &mut self,
        c1: &'a dyn ExtremaCurveTool,
        c2: &'a dyn ExtremaCurveTool,
        the_uinf: f64,
        the_usup: f64,
        the_vinf: f64,
        the_vsup: f64,
    ) {
        self.my_c1 = c1;
        self.my_c2 = c2;
        self.my_low_border.set(1, the_uinf);
        self.my_low_border.set(2, the_vinf);
        self.my_upp_border.set(1, the_usup);
        self.my_upp_border.set(2, the_vsup);
    }

    /// OCCT SetTolerance (hxx L387-396).
    pub fn set_tolerance(&mut self, the_tol: f64) {
        self.my_curve_min_tol = the_tol;
    }

    /// OCCT SetSingleSolutionFlag (hxx L407-416).
    pub fn set_single_solution_flag(&mut self, the_flag: bool) {
        self.my_is_find_single_solution = the_flag;
    }

    /// OCCT GetSingleSolutionFlag (hxx L427-436).
    pub fn get_single_solution_flag(&self) -> bool {
        self.my_is_find_single_solution
    }

    /// OCCT Perform() (hxx L447-802).
    pub fn perform(&mut self) {
        self.my_done = false;
        self.my_parallel = false;

        let c1 = self.my_c1;
        let c2 = self.my_c2;

        // hxx L461-471.
        let mut a_nb_inter = [
            c1.nb_intervals(GeomAbsShape::C2),
            c2.nb_intervals(GeomAbsShape::C2),
        ];
        let mut a_continuity = GeomAbsShape::C2;

        if a_nb_inter[0] as i64 * a_nb_inter[1] as i64 > 100 {
            a_continuity = GeomAbsShape::C1;
            a_nb_inter[0] = c1.nb_intervals(a_continuity);
            a_nb_inter[1] = c2.nb_intervals(a_continuity);
        }

        // hxx L473-511.
        let mut an_l = [0.0f64; 2];
        let mut indmax: i64 = -1;
        let mut indmin: i64 = -1;
        let mult = 20.0;
        let params_finite = !(is_infinite_value(c1.first_parameter())
            || is_infinite_value(c1.last_parameter())
            || is_infinite_value(c2.first_parameter())
            || is_infinite_value(c2.last_parameter()));
        if params_finite {
            // OCCT GCPnts_AbscissaPoint::Length(C1) / (C2) through the
            // Extrema_CurveTool facade.
            an_l[0] = c1.abscissa_length();
            an_l[1] = c2.abscissa_length();
            if an_l[0] >= 0.0 && an_l[1] >= 0.0 {
                if an_l[0] / a_nb_inter[0] as f64 > mult * an_l[1] / a_nb_inter[1] as f64 {
                    indmax = 0;
                    indmin = 1;
                } else if an_l[1] / a_nb_inter[1] as f64 > mult * an_l[0] / a_nb_inter[0] as f64 {
                    indmax = 1;
                    indmin = 0;
                }
            }
        }
        let mut a_nb_int_opt: i64 = 0;
        if indmax >= 0 {
            a_nb_int_opt = (an_l[indmax as usize] * a_nb_inter[indmin as usize] as f64
                / an_l[indmin as usize]
                / (mult / 4.0)) as i64
                + 1;
            if a_nb_int_opt > 100 || a_nb_int_opt < a_nb_inter[indmax as usize] as i64 {
                indmax = -1;
            } else if a_nb_int_opt * a_nb_inter[indmin as usize] as i64 > 100 {
                a_nb_int_opt = 100 / a_nb_inter[indmin as usize] as i64;
                if a_nb_int_opt < a_nb_inter[indmax as usize] as i64 {
                    indmax = -1;
                }
            }
        }

        // hxx L513-541.
        let mut an_intervals1: Vec<f64> = c1.intervals(a_continuity);
        let mut an_intervals2: Vec<f64> = c2.intervals(a_continuity);
        if indmax >= 0 {
            if indmax == 0 {
                ggen_ext_cc_change_intervals(&mut an_intervals1, a_nb_int_opt as usize);
                a_nb_inter[0] = an_intervals1.len() as i32 - 1;
            } else {
                ggen_ext_cc_change_intervals(&mut an_intervals2, a_nb_int_opt as usize);
                a_nb_inter[1] = an_intervals2.len() as i32 - 1;
            }
        }
        if c1.is_closed() && a_nb_inter[0] == 1 {
            ggen_ext_cc_change_intervals(&mut an_intervals1, 3);
            a_nb_inter[0] = an_intervals1.len() as i32 - 1;
        }
        if c2.is_closed() && a_nb_inter[1] == 1 {
            ggen_ext_cc_change_intervals(&mut an_intervals2, 3);
            a_nb_inter[1] = an_intervals2.len() as i32 - 1;
        }

        // hxx L543-579: the Lipschitz constant.
        let a_max_lc = 10000.0;
        let mut a_lc = 100.0f64;
        let a_max_der1 = 1.0 / c1.resolution(1.0);
        let a_max_der2 = 1.0 / c2.resolution(1.0);
        let mut a_max_der = a_max_der1.max(a_max_der2) * M_SQRT2;
        if a_lc > a_max_der {
            a_lc = a_max_der;
        }

        let mut is_const_locked_flag = false;
        let a_cr = 0.001;
        if a_max_der1 / a_max_der < a_cr || a_max_der2 / a_max_der < a_cr {
            is_const_locked_flag = true;
        }
        if a_max_der > a_max_lc {
            a_lc = a_max_lc;
            is_const_locked_flag = true;
        }
        if c1.get_type() == CurveType::Line {
            a_max_der = 1.0 / c2.resolution(1.0);
            if a_lc > a_max_der {
                is_const_locked_flag = true;
                a_lc = a_max_der;
            }
        }
        if c2.get_type() == CurveType::Line {
            a_max_der = 1.0 / c1.resolution(1.0);
            if a_lc > a_max_der {
                is_const_locked_flag = true;
                a_lc = a_max_der;
            }
        }

        // hxx L581-616: probe the largest gradient modulus on a 21x21 grid.
        {
            let mut a_func = GlobOptFuncCCC2::new_3d(c1, c2);
            if a_lc < a_max_lc || a_max_der > a_max_lc {
                let mut a_t = Vector::new(1, 2);
                let mut a_g = Vector::new(1, 2);
                let mut a_max_g = 0.0f64;
                let n1 = 21;
                let n2 = 21;
                let dt1 = (c1.last_parameter() - c1.first_parameter()) / (n1 - 1) as f64;
                let dt2 = (c2.last_parameter() - c2.first_parameter()) / (n2 - 1) as f64;
                let mut t1 = c1.first_parameter();
                for _i1 in 1..=n1 {
                    a_t.set(1, t1);
                    let mut t2 = c2.first_parameter();
                    for _i2 in 1..=n2 {
                        a_t.set(2, t2);
                        let mut a_f = 0.0;
                        a_func.values(&a_t, &mut a_f, &mut a_g);
                        let a_mod = a_g.get(1) * a_g.get(1) + a_g.get(2) * a_g.get(2);
                        a_max_g = a_max_g.max(a_mod);
                        t2 += dt2;
                    }
                    t1 += dt1;
                }
                a_max_g = a_max_g.sqrt();
                if a_max_g > a_max_der {
                    a_lc = a_max_g.min(a_max_lc);
                    is_const_locked_flag = true;
                }
                if a_max_g > 100.0 * a_max_lc {
                    a_lc = 100.0 * a_max_lc;
                    is_const_locked_flag = true;
                } else if a_max_g < 0.1 * a_max_der {
                    is_const_locked_flag = true;
                }
            }
        }

        // hxx L626-631: the duplicate-detection cell size and the filter.
        let a_cell_size = ((an_intervals1[an_intervals1.len() - 1] - an_intervals1[0])
            .max(an_intervals2[an_intervals2.len() - 1] - an_intervals2[0])
            * PCONFUSION
            / (2.0 * M_SQRT2))
            .max(PCONFUSION);

        let mut a_pnts: Vec<DVec2> = Vec::new();
        let mut a_f = REAL_LAST;

        // hxx L634-688: the global search over the interval boxes.
        {
            let mut a_func = GlobOptFuncCCC2::new_3d(c1, c2);
            let mut a_finder = MathGlobOptMin::new(
                MathMultiVarFunc::WithHessian(&mut a_func),
                &self.my_low_border.clone(),
                &self.my_upp_border.clone(),
                a_lc,
                1.0e-2,
                1.0e-7,
            );
            a_finder.set_lip_const_state(is_const_locked_flag);
            a_finder.set_continuity(if a_continuity == GeomAbsShape::C2 {
                2
            } else {
                1
            });
            let a_disc_tol = 1.0e-2;
            let a_value_tol = 1.0e-2;
            let a_same_tol = self.my_curve_min_tol / a_disc_tol;
            a_finder.set_tol(a_disc_tol, a_same_tol);
            a_finder.set_functional_minimal_value(0.0);

            let mut an_inspector = GGenCcPointsInspector::new(a_cell_size);
            let mut a_filter: CellFilter<DVec2, GGenCcPointsInspector> =
                CellFilter::new(2, a_cell_size);

            let mut a_first_border_interval = Vector::new(1, 2);
            let mut a_second_border_interval = Vector::new(1, 2);
            for i in 1..=a_nb_inter[0] {
                for j in 1..=a_nb_inter[1] {
                    a_first_border_interval.set(1, an_intervals1[(i - 1) as usize]);
                    a_first_border_interval.set(2, an_intervals2[(j - 1) as usize]);
                    a_second_border_interval.set(1, an_intervals1[i as usize]);
                    a_second_border_interval.set(2, an_intervals2[j as usize]);

                    a_finder.set_local_params(&a_first_border_interval, &a_second_border_interval);
                    a_finder.perform(self.get_single_solution_flag());

                    let a_curr_f = a_finder.get_f();
                    if a_curr_f >= a_f + a_same_tol * a_value_tol {
                        continue;
                    }

                    if a_curr_f > a_f - a_same_tol * a_value_tol {
                        if a_curr_f < a_f {
                            a_f = a_curr_f;
                        }
                    } else {
                        a_f = a_curr_f;
                        a_filter.reset_scalar(a_cell_size);
                        a_pnts.clear();
                    }

                    let mut sol = Vector::new(1, 2);
                    for k in 1..=a_finder.nb_extrema() {
                        a_finder.points(k, &mut sol);
                        let a_pnt2d = DVec2::new(sol.get(1), sol.get(2));

                        let a_xy_min = an_inspector.shift(&a_pnt2d, -a_cell_size);
                        let a_xy_max = an_inspector.shift(&a_pnt2d, a_cell_size);

                        an_inspector.clear_find();
                        an_inspector.set_current(&a_pnt2d);
                        let mut found = false;
                        {
                            let mut inspector_ref = &mut an_inspector;
                            a_filter.inspect_range(
                                &a_xy_min,
                                &a_xy_max,
                                &GGenCcPointsInspector::coord,
                                &mut |obj: &DVec2| {
                                    let r = inspector_ref.inspect(obj);
                                    if inspector_ref.is_find() {
                                        found = true;
                                    }
                                    r
                                },
                            );
                        }
                        if !found {
                            a_filter.add(&a_pnt2d, &a_pnt2d, &GGenCcPointsInspector::coord);
                            a_pnts.push(a_pnt2d);
                        }
                    }
                }
            }
        }

        // hxx L690-705.
        let a_nb_sol = a_pnts.len();
        if a_nb_sol == 0 {
            self.my_done = false;
            return;
        }

        self.my_done = true;

        if a_nb_sol == 1 {
            let a_sol = a_pnts[0];
            self.my_points1.push(a_sol.x);
            self.my_points2.push(a_sol.y);
            return;
        }

        // hxx L707: std::sort with the file-local comparator.
        a_pnts.sort_by(ggen_ext_cc_comp);

        // hxx L709-759.
        let mut a_solutions: Vec<usize> = Vec::new();
        let mut b_save_solution = true;
        let mut b_dirs_coincide = true;
        let mut b_different_solutions = false;

        let mut is_parallel = true;
        {
            let mut a_func = GlobOptFuncCCC2::new_3d(c1, c2);
            let mut a_val = 0.0f64;
            let mut a_vec = Vector::new_init(1, 2, 0.0);

            for an_idx in 0..(a_nb_sol - 1) {
                let a_current = a_pnts[an_idx];
                let a_next = a_pnts[an_idx + 1];

                a_vec.set(1, (a_current.x + a_next.x) * 0.5);
                a_vec.set(2, (a_current.y + a_next.y) * 0.5);

                a_func.value(&a_vec, &mut a_val);
                if (a_val - a_f).abs() < CONFUSION {
                    if b_save_solution {
                        a_solutions.push(an_idx);
                        b_save_solution = false;
                    }
                } else {
                    is_parallel = false;
                    a_solutions.push(an_idx);
                    b_save_solution = true;
                }

                if !b_different_solutions && a_next.x > a_current.x {
                    if a_next.y > a_current.y {
                        b_different_solutions = true;
                        b_dirs_coincide = true;
                    } else if a_next.y < a_current.y {
                        b_different_solutions = true;
                        b_dirs_coincide = false;
                    }
                }
            }
        }
        a_solutions.push(a_nb_sol - 1);

        if !b_different_solutions {
            is_parallel = false;
        }

        // hxx L765-783: the parallel check with the point-curve extrema.
        if is_parallel {
            let a_t1 = [self.my_low_border.get(1), self.my_upp_border.get(1)];
            let a_t2 = [
                if b_dirs_coincide {
                    self.my_low_border.get(2)
                } else {
                    self.my_upp_border.get(2)
                },
                if b_dirs_coincide {
                    self.my_upp_border.get(2)
                } else {
                    self.my_low_border.get(2)
                },
            ];

            // OCCT: TheExtPC anExtPC1, anExtPC2 with Initialize(C, Low, Upp).
            let mut an_ext_pc1 = ExtPC::new_fn(
                DVec3::ZERO,
                self.my_curve_min_tol,
                self.my_low_border.get(1),
                self.my_upp_border.get(1),
                &|u| c1.value(u),
            );
            let mut an_ext_pc2 = ExtPC::new_fn(
                DVec3::ZERO,
                self.my_curve_min_tol,
                self.my_low_border.get(2),
                self.my_upp_border.get(2),
                &|u| c2.value(u),
            );

            for i_t in 0..2 {
                if !is_parallel {
                    break;
                }
                let a_dist1 = ggen_ext_cc_proj_p_on_c(
                    c1.value(a_t1[i_t]),
                    &mut an_ext_pc2,
                    &|u| c2.value(u),
                    self.my_low_border.get(2),
                    self.my_upp_border.get(2),
                );
                let a_dist2 = ggen_ext_cc_proj_p_on_c(
                    c2.value(a_t2[i_t]),
                    &mut an_ext_pc1,
                    &|u| c1.value(u),
                    self.my_low_border.get(1),
                    self.my_upp_border.get(1),
                );
                is_parallel = (a_dist1.min(a_dist2) - a_f * a_f).abs() < CONFUSION;
            }
        }

        // hxx L785-801.
        if is_parallel {
            let a_sol = a_pnts[0];
            self.my_points1.push(a_sol.x);
            self.my_points2.push(a_sol.y);
            self.my_parallel = true;
        } else {
            for idx in &a_solutions {
                let a_sol = a_pnts[*idx];
                self.my_points1.push(a_sol.x);
                self.my_points2.push(a_sol.y);
            }
        }
    }

    /// OCCT IsDone() (hxx L813-822).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsParallel() (hxx L833-844).
    pub fn is_parallel(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_parallel
    }

    /// OCCT NbExt() (hxx L855-866).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_points1.len()
    }

    /// OCCT SquareDistance(N) (hxx L877-892) — 1-based.
    pub fn square_distance(&self, the_n: usize) -> f64 {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        let p1 = self.my_c1.value(self.my_points1[the_n - 1]);
        let p2 = self.my_c2.value(self.my_points2[the_n - 1]);
        (p1 - p2).length_squared()
    }

    /// OCCT Points(N, P1, P2) (hxx L903-918) — 1-based.
    pub fn points(&self, the_n: usize, the_p1: &mut POnCurve, the_p2: &mut POnCurve) {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        let t1 = self.my_points1[the_n - 1];
        let t2 = self.my_points2[the_n - 1];
        the_p1.param = t1;
        the_p1.point = self.my_c1.value(t1);
        the_p2.param = t2;
        the_p2.point = self.my_c2.value(t2);
    }
}
