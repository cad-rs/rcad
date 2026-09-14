// OCCT Approx_ComputeLine (TKGeomBase Approx) — 1:1 Rust translation of
// `Approx_ComputeLine.gxx` (L134-1750) as instantiated by
// `AppDef_Compute_0.cxx` (L29-53):
//   Approx_ComputeLine                = AppDef_Compute (rcad `Compute`)
//   MultiLine                         = AppDef_MultiLine (rcad app_def::MultiLine)
//   LineTool                          = AppDef_MyLineTool (rcad app_def::my_line_tool)
//   Approx_MyGradient                 = AppDef_MyGradientOfCompute
//     (= AppParCurves_Gradient.gxx, rcad app_par_curves::Gradient)
//   Approx_ParLeastSquareOfMyGradient = AppDef_ParLeastSquareOfMyGradientOfCompute
//     (= AppParCurves_LeastSquare.gxx, rcad app_par_curves::LeastSquare)
//
// The struct keeps the AppDef instantiation name (Compute in module
// app_def_compute); the gxx is shared with the Approx_ComputeLine and the
// GeomInt/BRepApprox ComputeLineBezier bindings, which will extract the
// generic engine when their consumers land (the bspl_compute_line
// precedent: concrete first, genericize at the second instantiation).
//
// Scope note: Approx_ComputeLine.gxx contains no AppDef_Variational /
// smoothing branch anywhere (checked L134-1750) — the only smoothing user
// of this machinery family is Approx_BSplComputeLine.gxx, translated
// separately in bspl_compute_line.rs. There is no GAP branch.
//
// OCCT quirk kept as a comment: `ComputeCurve` is declared in
// AppDef_Compute.hxx (L168-170) AND defined here (gxx L1443-1690) — unlike
// the stale BSplineCompute declaration, this one has a body.

use glam::{DVec2, DVec3};
use rcad_kernel::core::precision::{ANGULAR, CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use rcad_kernel::math::math_matrix::Vector as RVector;
use rcad_kernel::math::VecD;

use super::app_def::{my_line_tool, MultiLine};
use super::app_par_curves::{Gradient, LeastSquare};
use super::approx_int::{
    ApproxParamType, ApproxStatus, AppParConstraint, ConstraintCouple, MCurvesToBSpCurve,
    MultiBSpCurve, MultiCurve, MultiPoint,
};

/// OCCT Standard_Real RealLast() (Standard_Real.hxx).
const REAL_LAST: f64 = f64::MAX;

/// OCCT gp::Resolution() == RealSmall() (gp.hxx L59-60) — the smallest
/// positive normalized double.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT Standard::Epsilon(1.0) (Standard_Real.hxx L242-246) —
/// nextafter(1.0, RealLast()) - 1.0 == DBL_EPSILON.
const EPSILON_1: f64 = f64::EPSILON;

/// OCCT gp_Dir::Angle (gp_Dir.cxx L27-50) — the angle in [0, PI] between
/// the (unit) directions, computed with the acos/asin switch at 45 degrees.
fn gp_dir_angle(coord: DVec3, other: DVec3) -> f64 {
    let cosinus = coord.dot(other);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = coord.cross(other).length();
        if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488-493) — raises
/// VectorWithNullMagnitude when a magnitude <= gp::Resolution(); the
/// rcad mapping of the raise is a panic.
fn gp_vec_angle(coord: DVec3, other: DVec3) -> f64 {
    assert!(
        coord.length() > GP_RESOLUTION && other.length() > GP_RESOLUTION,
        "gp_VectorWithNullMagnitude: gp_Vec::Angle"
    );
    // gp_Dir(coord) normalizes the coordinates.
    gp_dir_angle(coord.normalize(), other.normalize())
}

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx L136-141).
fn gp_vec_is_parallel(coord: DVec3, other: DVec3, angular_tolerance: f64) -> bool {
    let an_ang = gp_vec_angle(coord, other);
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT gp_Vec2d::Angle (gp_Vec2d.cxx L47-87) — the signed angle in
/// (-PI, PI]; raises VectorWithNullMagnitude when a magnitude <=
/// gp::Resolution().
fn gp_vec2d_angle(coord: DVec2, other: DVec2) -> f64 {
    let a_norm = coord.length();
    let an_other_norm = other.length();
    assert!(
        a_norm > GP_RESOLUTION && an_other_norm > GP_RESOLUTION,
        "gp_VectorWithNullMagnitude: gp_Vec2d::Angle"
    );

    let a_d = a_norm * an_other_norm;
    let a_cosinus = coord.dot(other) / a_d;
    // OCCT coord.Crossed(theOther.coord) — the 2d cross product scalar.
    let a_sinus = coord.perp_dot(other) / a_d;

    // Use M_SQRT1_2 (1/sqrt(2)) for better readability and precision.
    const A_COS_45_DEG: f64 = std::f64::consts::FRAC_1_SQRT_2;

    if a_cosinus > -A_COS_45_DEG && a_cosinus < A_COS_45_DEG {
        // For angles near +/-90 degrees, use acos for better precision.
        if a_sinus > 0.0 {
            a_cosinus.acos()
        } else {
            -a_cosinus.acos()
        }
    } else {
        // For angles near 0 degrees or +/-180 degrees, use asin for better
        // precision.
        if a_cosinus > 0.0 {
            a_sinus.asin()
        } else if a_sinus > 0.0 {
            std::f64::consts::PI - a_sinus.asin()
        } else {
            -std::f64::consts::PI - a_sinus.asin()
        }
    }
}

/// OCCT gp_Vec2d::IsParallel (gp_Vec2d.hxx L372-376).
fn gp_vec2d_is_parallel(coord: DVec2, other: DVec2, angular_tolerance: f64) -> bool {
    let an_ang = gp_vec2d_angle(coord, other).abs();
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT CheckMultiCurve static (Approx_ComputeLine.gxx L134-428) — detects
/// "loop" results (the approximated Bezier folding back onto itself) and
/// reports a bad point index. Returns true when the MultiCurve is
/// acceptable.
fn check_multi_curve(
    the_multi_curve: &MultiCurve,
    the_line: &MultiLine,
    the_indfirst: i32,
    the_indlast: i32,
    the_indbad: &mut i32,
) -> bool {
    let nbp3d = my_line_tool::nb_p3d(the_line) as i32;
    let nbp2d = my_line_tool::nb_p2d(the_line) as i32;

    let coeff = 4.0f64; // 2*2

    if nbp3d > 1 {
        // only simple cases are analysed
        return true;
    }

    let min_scal_prod = -0.9f64;
    let sq_tol3d = SQUARE_CONFUSION;

    *the_indbad = 0;
    // OCCT `int indbads[4]` scans indexes 1..3; the 2d branch writes
    // indbads[indcur] up to NbCur, which is an out-of-bounds write in OCCT
    // for NbCur > 3 (UB) — rcad sizes the array to stay in bounds.
    let nb_cur = the_multi_curve.nb_curves() as i32;
    let mut indbads = vec![0i32; (nb_cur as usize + 1).max(4)];
    let mut loop_found = false;

    let a_nb_p3d = nbp3d.max(1);
    let a_nb_p2d = nbp2d.max(1);

    let mut tab_p = vec![DVec3::ZERO; a_nb_p3d as usize];
    let mut tab_p2d = vec![DVec2::ZERO; a_nb_p2d as usize];

    if the_multi_curve.dimension(1) == 3 {
        // 3d case (gxx L169-300)
        let mut a_poles: Vec<DVec3> = Vec::new();
        the_multi_curve.curve(1, &mut a_poles);
        let mut first_vec = DVec3::ZERO;
        let mut second_vec;
        let mut indp = 2usize;
        while indp <= a_poles.len() {
            // OCCT: FirstVec = gp_Vec(aPoles(1), aPoles(indp++)) — the
            // subscript uses indp before the post-increment.
            first_vec = a_poles[indp - 1] - a_poles[0];
            indp += 1;
            let a_length = first_vec.length();
            if a_length > GP_RESOLUTION {
                first_vec /= a_length;
                break;
            }
        }
        // OCCT: MidPnt = aPoles(indp - 1) — 1-based index indp-1.
        let mut mid_pnt = a_poles[indp - 2];
        while indp <= a_poles.len() {
            second_vec = a_poles[indp - 1] - mid_pnt;
            let a_length = second_vec.length();
            if a_length <= GP_RESOLUTION {
                indp += 1;
                continue;
            }
            second_vec /= a_length;
            let scal_prod = first_vec.dot(second_vec);
            if scal_prod < min_scal_prod {
                loop_found = true;
                break;
            }
            first_vec = second_vec;
            mid_pnt = a_poles[indp - 1];
            indp += 1;
        }
        // Check: may be it is a real loop (gxx L217-257)
        if loop_found {
            for first_ind in the_indfirst..=(the_indlast - 2) {
                my_line_tool::value_3d(the_line, first_ind as usize, &mut tab_p);
                let first_pnt = tab_p[0];
                for k in (first_ind + 1)..the_indlast {
                    my_line_tool::value_3d(the_line, k as usize, &mut tab_p);
                    let pnt1 = tab_p[0];
                    my_line_tool::value_3d(the_line, (k + 1) as usize, &mut tab_p);
                    let pnt2 = tab_p[0];
                    if first_pnt.distance_squared(pnt1) <= sq_tol3d
                        || first_pnt.distance_squared(pnt2) <= sq_tol3d
                    {
                        loop_found = false;
                        break;
                    }
                    // OCCT Vec1.Normalize() — raises on a null vector; the
                    // guard above rules the null case out.
                    let vec1 = (pnt1 - first_pnt).normalize();
                    let vec2 = (pnt2 - first_pnt).normalize();
                    let scal_prod = vec1.dot(vec2);
                    if scal_prod < min_scal_prod {
                        loop_found = false;
                        break;
                    }
                }
                if !loop_found {
                    break;
                }
            }
        }
        if loop_found {
            // search <indbad> (gxx L258-299)
            let mut max_sq_dist = 0.0f64;
            let mut min_sq_dist = REAL_LAST;
            for k in (the_indfirst + 1)..=the_indlast {
                my_line_tool::value_3d(the_line, (k - 1) as usize, &mut tab_p);
                let prev_pnt = tab_p[0];
                my_line_tool::value_3d(the_line, k as usize, &mut tab_p);
                let cur_pnt = tab_p[0];
                let a_sq_dist = prev_pnt.distance_squared(cur_pnt);
                if a_sq_dist > max_sq_dist {
                    max_sq_dist = a_sq_dist;
                    indbads[1] = k;
                }
                if a_sq_dist > GP_RESOLUTION && a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                }
            }
            let relation = max_sq_dist / min_sq_dist;
            if relation < coeff {
                loop_found = false;
            } else {
                for indcur in 2..=nb_cur {
                    max_sq_dist = 0.0;
                    for k in (the_indfirst + 1)..=the_indlast {
                        // OCCT quirk (gxx L287-290): the 2d Value overload is
                        // called here even in the 3d case.
                        my_line_tool::value_2d(the_line, (k - 1) as usize, &mut tab_p2d);
                        let prev_pnt = tab_p2d[(indcur - 2) as usize];
                        my_line_tool::value_2d(the_line, k as usize, &mut tab_p2d);
                        let cur_pnt = tab_p2d[(indcur - 2) as usize];
                        let a_sq_dist = prev_pnt.distance_squared(cur_pnt);
                        if a_sq_dist > max_sq_dist {
                            max_sq_dist = a_sq_dist;
                            indbads[indcur as usize] = k;
                        }
                    }
                }
            }
        }
    } else {
        // 2d case (gxx L301-414)
        let mut a_poles2d: Vec<DVec2> = Vec::new();
        the_multi_curve.curve2d(1, &mut a_poles2d);
        let a_sq_norm_toler = EPSILON_1 * EPSILON_1;
        let mut first_vec = a_poles2d[1] - a_poles2d[0];
        let mut second_vec;
        let mut a_vec_sq_norm = first_vec.length_squared();
        if a_vec_sq_norm < a_sq_norm_toler {
            *the_indbad = the_indfirst + 1;
            return false;
        }

        // OCCT divides by std::sqrt(aSqNormToler), not by the actual norm.
        first_vec /= a_sq_norm_toler.sqrt();
        let mut mid_pnt = a_poles2d[1];
        for k in 3..=a_poles2d.len() {
            second_vec = a_poles2d[k - 1] - mid_pnt;
            a_vec_sq_norm = second_vec.length_squared();
            if a_vec_sq_norm < a_sq_norm_toler {
                *the_indbad = the_indfirst + k as i32 - 1;
                return false;
            }

            second_vec /= a_vec_sq_norm.sqrt();
            let scal_prod = first_vec.dot(second_vec);
            if scal_prod < min_scal_prod {
                loop_found = true;
                break;
            }
            first_vec = second_vec;
            mid_pnt = a_poles2d[k - 1];
        }
        // Check: may be it is a real loop (gxx L346-386)
        if loop_found {
            for first_ind in the_indfirst..=(the_indlast - 2) {
                my_line_tool::value_2d(the_line, first_ind as usize, &mut tab_p2d);
                let first_pnt = tab_p2d[0];
                for k in (first_ind + 1)..the_indlast {
                    my_line_tool::value_2d(the_line, k as usize, &mut tab_p2d);
                    let pnt1 = tab_p2d[0];
                    my_line_tool::value_2d(the_line, (k + 1) as usize, &mut tab_p2d);
                    let pnt2 = tab_p2d[0];
                    if first_pnt.distance_squared(pnt1) <= sq_tol3d
                        || first_pnt.distance_squared(pnt2) <= sq_tol3d
                    {
                        loop_found = false;
                        break;
                    }
                    let vec1 = (pnt1 - first_pnt).normalize();
                    let vec2 = (pnt2 - first_pnt).normalize();
                    let scal_prod = vec1.dot(vec2);
                    if scal_prod < min_scal_prod {
                        loop_found = false;
                        break;
                    }
                }
                if !loop_found {
                    break;
                }
            }
        }
        if loop_found {
            // search <indbad> (gxx L387-413)
            for indcur in 1..=nb_cur {
                let mut max_sq_dist = 0.0f64;
                let mut min_sq_dist = REAL_LAST;
                for k in (the_indfirst + 1)..=the_indlast {
                    my_line_tool::value_2d(the_line, (k - 1) as usize, &mut tab_p2d);
                    let prev_pnt = tab_p2d[(indcur - 1) as usize];
                    my_line_tool::value_2d(the_line, k as usize, &mut tab_p2d);
                    let cur_pnt = tab_p2d[(indcur - 1) as usize];
                    let a_sq_dist = prev_pnt.distance_squared(cur_pnt);
                    if a_sq_dist > max_sq_dist {
                        max_sq_dist = a_sq_dist;
                        indbads[indcur as usize] = k;
                    }
                    if a_sq_dist > GP_RESOLUTION && a_sq_dist < min_sq_dist {
                        min_sq_dist = a_sq_dist;
                    }
                }
                let relation = max_sq_dist / min_sq_dist;
                if relation < coeff {
                    loop_found = false;
                }
            }
        }
    }

    // Define <indbad> (gxx L416-422)
    for i in 1..=3 {
        if indbads[i] != 0 {
            *the_indbad = indbads[i];
            break;
        }
    }

    if !loop_found {
        *the_indbad = 0;
    }

    !loop_found
}

/// OCCT AppDef_Compute (AppDef_Compute.hxx L44-221; body
/// Approx_ComputeLine.gxx L430-1750).
pub struct Compute {
    /// OCCT myMultiCurves.
    my_multi_curves: Vec<MultiCurve>,
    /// OCCT TheMultiCurve — the default AppParCurves_MultiCurve() has a null
    /// pole array (NbCurves() == 0); rcad uses one empty MultiPoint, whose
    /// nb_curves() is 0 as well.
    the_multi_curve: MultiCurve,
    /// OCCT myspline.
    myspline: MultiBSpCurve,
    alldone: bool,
    /// OCCT leaves tolreached uninitialized in the constructors (gxx
    /// L717-830); rcad initializes it to false.
    tolreached: bool,
    /// OCCT Par.
    par: ApproxParamType,
    /// OCCT myParameters (HArray1, null until Compute stores the best fit;
    /// the reads behind the NbCurves() guards dereference it unconditionally
    /// in OCCT).
    my_parameters: Option<RVector>,
    /// OCCT myfirstParam (HArray1, null unless a Parameters-taking
    /// constructor was used; Nullify()-ed after its first use in Perform).
    myfirst_param: Option<RVector>,
    /// OCCT myPar.
    my_par: Vec<RVector>,
    /// OCCT Tolers3d.
    tolers3d: Vec<f64>,
    /// OCCT Tolers2d.
    tolers2d: Vec<f64>,
    /// OCCT myConstraints (HArray1(1, 2)).
    my_constraints: [ConstraintCouple; 2],
    mydegremin: i32,
    mydegremax: i32,
    mytol3d: f64,
    mytol2d: f64,
    /// OCCT leaves currenttol3d/currenttol2d uninitialized in the
    /// constructors; Compute assigns them (gxx L1353) before any read.
    currenttol3d: f64,
    currenttol2d: f64,
    mycut: bool,
    mysquares: bool,
    myitermax: i32,
    /// OCCT myfirstC.
    myfirstc: AppParConstraint,
    /// OCCT mylastC.
    mylastc: AppParConstraint,
    /// OCCT myMultiLineNb.
    my_multi_line_nb: i32,
    /// OCCT myIsClear.
    my_is_clear: bool,
}

impl Compute {
    /// The shared member initialization of the four OCCT constructors
    /// (gxx L726-746 / L757-777 / L787-802 / L813-829). Every constructor
    /// sets myMultiLineNb(0), myIsClear(false), myConstraints(1, 2),
    /// myfirstC/mylastC = TangencyPoint, alldone = false and the algorithm
    /// parameters; they differ only in Par (IsoParametric for the
    /// Parameters-taking constructors, parametrization for the two others)
    /// and in myfirstParam. OCCT leaves tolreached, currenttol3d,
    /// currenttol2d and myspline uninitialized — rcad gives them the
    /// documented defaults.
    fn make_members(
        par: ApproxParamType,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        Compute {
            my_multi_curves: Vec::new(),
            the_multi_curve: MultiCurve::new(1, 0, 0),
            myspline: MultiBSpCurve::new_nbpol(0),
            alldone: false,
            tolreached: false,
            par,
            my_parameters: None,
            myfirst_param: None,
            my_par: Vec::new(),
            tolers3d: Vec::new(),
            tolers2d: Vec::new(),
            my_constraints: [
                ConstraintCouple {
                    index: 0,
                    constraint: AppParConstraint::NoConstraint,
                },
                ConstraintCouple {
                    index: 0,
                    constraint: AppParConstraint::NoConstraint,
                },
            ],
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            currenttol3d: REAL_LAST,
            currenttol2d: REAL_LAST,
            mycut: cutting,
            mysquares: squares,
            myitermax: nb_iterations,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            my_multi_line_nb: 0,
            my_is_clear: false,
        }
    }

    /// OCCT Approx_ComputeLine(Line, Parameters, degreemin..Squares)
    /// (gxx L717-747) — Par is forced to Approx_IsoParametric; the MultiLine
    /// is approximated right away.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_line_and_parameters(
        line: &MultiLine,
        parameters: &RVector,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        let mut r = Compute::make_members(
            ApproxParamType::IsoParametric,
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            squares,
        );
        // myfirstParam = new HArray1(Parameters.Lower(), Parameters.Upper());
        let mut fp = RVector::new(parameters.lower(), parameters.upper());
        for i in parameters.lower()..=parameters.upper() {
            fp.set(i, parameters.get(i));
        }
        r.myfirst_param = Some(fp);
        r.perform(line);
        r
    }

    /// OCCT Approx_ComputeLine(Parameters, degreemin..Squares)
    /// (gxx L749-777) — Par is forced to Approx_IsoParametric.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_parameters(
        parameters: &RVector,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        squares: bool,
    ) -> Self {
        let mut r = Compute::make_members(
            ApproxParamType::IsoParametric,
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            squares,
        );
        let mut fp = RVector::new(parameters.lower(), parameters.upper());
        for i in parameters.lower()..=parameters.upper() {
            fp.set(i, parameters.get(i));
        }
        r.myfirst_param = Some(fp);
        r
    }

    /// OCCT Approx_ComputeLine(degreemin, degreemax, Tolerance3d, Tolerance2d,
    /// NbIterations, cutting, parametrization, Squares) (gxx L779-802).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) -> Self {
        Compute::make_members(
            parametrization,
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            squares,
        )
    }

    /// OCCT Approx_ComputeLine(Line, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, NbIterations, cutting, parametrization, Squares)
    /// (gxx L804-830) — the MultiLine is approximated right away.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_line(
        line: &MultiLine,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) -> Self {
        let mut r = Compute::make_members(
            parametrization,
            degreemin,
            degreemax,
            tolerance3d,
            tolerance2d,
            nb_iterations,
            cutting,
            squares,
        );
        r.perform(line);
        r
    }

    /// OCCT FirstTangencyVector(Line, index, V) (gxx L430-512).
    fn first_tangency_vector(&self, line: &MultiLine, index: i32, v: &mut RVector) {
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_v = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_v2d = vec![DVec2::ZERO; mynb_p2d as usize];

        let ok = if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::tangency_3d_2d(line, index as usize, &mut tab_v, &mut tab_v2d)
        } else if nb_p2d != 0 {
            my_line_tool::tangency_2d(line, index as usize, &mut tab_v2d)
        } else {
            my_line_tool::tangency_3d(line, index as usize, &mut tab_v)
        };

        if ok {
            if nb_p3d != 0 {
                let mut j = 1;
                for t in tab_v.iter().take(nb_p3d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    v.set(j + 2, t.z);
                    j += 3;
                }
            }
            if nb_p2d != 0 {
                let mut j = nb_p3d * 3 + 1;
                for t in tab_v2d.iter().take(nb_p2d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    j += 2;
                }
            }
        } else {
            // Search for a tangent vector by construction of a parabola:
            let first_c = AppParConstraint::PassPoint;
            let last_c = AppParConstraint::PassPoint;
            let nbpoles = 3;
            let mut mypar = RVector::new(index, index + 2);
            self.parameters_compute(line, index, index + 2, &mut mypar);
            let mut lsq = LeastSquare::new(
                line,
                index,
                index + 2,
                first_c,
                last_c,
                &rvector_to_vecd(&mypar),
                nbpoles,
            );
            let c = lsq.bezier_value();

            let mut j = 1;
            for i in 1..=nb_p3d {
                let (_p, v1) = c.d1(i as usize, 0.0);
                v.set(j, v1.x);
                v.set(j + 1, v1.y);
                v.set(j + 2, v1.z);
                j += 3;
            }
            let mut j = nb_p3d * 3 + 1;
            for i in (nb_p3d + 1)..=(nb_p3d + nb_p2d) {
                let (_p2d, v2d) = c.d1(i as usize, 0.0);
                v.set(j, v2d.x);
                v.set(j + 1, v2d.y);
                j += 2;
            }
        }
    }

    /// OCCT LastTangencyVector(Line, index, V) (gxx L514-595) — same shape
    /// with the D1 evaluated at parameter 1.0 on the last parabola segment.
    fn last_tangency_vector(&self, line: &MultiLine, index: i32, v: &mut RVector) {
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_v = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_v2d = vec![DVec2::ZERO; mynb_p2d as usize];

        let ok = if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::tangency_3d_2d(line, index as usize, &mut tab_v, &mut tab_v2d)
        } else if nb_p2d != 0 {
            my_line_tool::tangency_2d(line, index as usize, &mut tab_v2d)
        } else {
            my_line_tool::tangency_3d(line, index as usize, &mut tab_v)
        };

        if ok {
            if nb_p3d != 0 {
                let mut j = 1;
                for t in tab_v.iter().take(nb_p3d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    v.set(j + 2, t.z);
                    j += 3;
                }
            }
            if nb_p2d != 0 {
                let mut j = nb_p3d * 3 + 1;
                for t in tab_v2d.iter().take(nb_p2d as usize) {
                    v.set(j, t.x);
                    v.set(j + 1, t.y);
                    j += 2;
                }
            }
        } else {
            // Search for a tangent vector by construction of a parabola:
            let first_c = AppParConstraint::PassPoint;
            let last_c = AppParConstraint::PassPoint;
            let nbpoles = 3;
            let mut mypar = RVector::new(index - 2, index);
            self.parameters_compute(line, index - 2, index, &mut mypar);
            let mut lsq = LeastSquare::new(
                line,
                index - 2,
                index,
                first_c,
                last_c,
                &rvector_to_vecd(&mypar),
                nbpoles,
            );
            let c = lsq.bezier_value();

            let mut j = 1;
            for i in 1..=nb_p3d {
                let (_p, v1) = c.d1(i as usize, 1.0);
                v.set(j, v1.x);
                v.set(j + 1, v1.y);
                v.set(j + 2, v1.z);
                j += 3;
            }
            let mut j = nb_p3d * 3 + 1;
            for i in (nb_p3d + 1)..=(nb_p3d + nb_p2d) {
                let (_p2d, v2d) = c.d1(i as usize, 1.0);
                v.set(j, v2d.x);
                v.set(j + 1, v2d.y);
                j += 2;
            }
        }
    }

    /// OCCT SearchFirstLambda(Line, TheParam, V, index) (gxx L597-655).
    fn search_first_lambda(
        &self,
        line: &MultiLine,
        the_param: &RVector,
        v: &RVector,
        index: i32,
    ) -> f64 {
        // dq/dw = lambda* V = (p2-p1)/(u2-u1)
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_p1 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p2 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p12d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_p22d = vec![DVec2::ZERO; mynb_p2d as usize];

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, index as usize, &mut tab_p1, &mut tab_p12d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, index as usize, &mut tab_p12d);
        } else if nb_p3d != 0 {
            my_line_tool::value_3d(line, index as usize, &mut tab_p1);
        }

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, (index + 1) as usize, &mut tab_p2, &mut tab_p22d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, (index + 1) as usize, &mut tab_p22d);
        } else if nb_p3d != 0 {
            my_line_tool::value_3d(line, (index + 1) as usize, &mut tab_p2);
        }

        let u1 = the_param.get(index);
        let u2 = the_param.get(index + 1);
        let low = v.lower(); // OCCT V.Lower().

        let (lambda, s);
        if nb_p3d != 0 {
            let p1 = tab_p1[0];
            let p2 = tab_p2[0];
            let p1p2 = p2 - p1;
            let my_v = DVec3::new(v.get(low), v.get(low + 1), v.get(low + 2));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        } else {
            let p12d = tab_p12d[0];
            let p22d = tab_p22d[0];
            let p1p2 = p22d - p12d;
            let my_v = DVec2::new(v.get(low), v.get(low + 1));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        }
        s * lambda
    }

    /// OCCT SearchLastLambda(Line, TheParam, V, index) (gxx L657-715).
    fn search_last_lambda(
        &self,
        line: &MultiLine,
        the_param: &RVector,
        v: &RVector,
        index: i32,
    ) -> f64 {
        // dq/dw = lambda* V = (p2-p1)/(u2-u1)
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }
        let mut tab_p = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p2 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p2d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_p22d = vec![DVec2::ZERO; mynb_p2d as usize];

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, (index - 1) as usize, &mut tab_p, &mut tab_p2d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, (index - 1) as usize, &mut tab_p2d);
        } else if nb_p3d != 0 {
            my_line_tool::value_3d(line, (index - 1) as usize, &mut tab_p);
        }

        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, index as usize, &mut tab_p2, &mut tab_p22d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, index as usize, &mut tab_p22d);
        } else if nb_p3d != 0 {
            my_line_tool::value_3d(line, index as usize, &mut tab_p2);
        }

        let u1 = the_param.get(index - 1);
        let u2 = the_param.get(index);
        let low = v.lower();

        let (lambda, s);
        if nb_p3d != 0 {
            let p1 = tab_p[0];
            let p2 = tab_p2[0];
            let p1p2 = p2 - p1;
            let my_v = DVec3::new(v.get(low), v.get(low + 1), v.get(low + 2));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        } else {
            let p12d = tab_p2d[0];
            let p22d = tab_p22d[0];
            let p1p2 = p22d - p12d;
            let my_v = DVec2::new(v.get(low), v.get(low + 1));
            lambda = p1p2.length() / (my_v.length() * (u2 - u1));
            s = if p1p2.dot(my_v) > 0.0 { 1.0 } else { -1.0 };
        }
        s * lambda
    }

    /// OCCT Perform(Line) (gxx L832-1219).
    pub fn perform(&mut self, line: &MultiLine) {
        if !self.my_is_clear {
            self.my_multi_curves.clear();
            self.my_par.clear();
            self.tolers3d.clear();
            self.tolers2d.clear();
            self.my_multi_line_nb = 0;
        } else {
            self.my_is_clear = false;
        }

        let mut nbp;
        let mut oldlastpt;
        let mut finish = false;
        let mut begin = true;
        let mut ok = false;
        let mut go_up = false;
        let mut my_status;
        // thetol3d/thetol2d are uninitialized locals in OCCT, only written
        // through Compute's output references; Perform never reads them.
        let mut thetol3d = 0.0f64;
        let mut thetol2d = 0.0f64;

        let thefirstpt = my_line_tool::first_point(line) as i32;
        let thelastpt = my_line_tool::last_point(line) as i32;
        let mut myfirstpt = thefirstpt;
        let mut mylastpt = thelastpt;

        let mut my_couple1 = ConstraintCouple {
            index: myfirstpt,
            constraint: self.myfirstc,
        };
        let mut my_couple2 = ConstraintCouple {
            index: mylastpt,
            constraint: self.mylastc,
        };
        self.my_constraints[0] = my_couple1;
        self.my_constraints[1] = my_couple2;

        let mut the_param = RVector::new(thefirstpt, thelastpt);

        if !self.mycut {
            // Case where no additional points are desired (gxx L867-915).
            match &self.myfirst_param {
                None => {
                    self.parameters_compute(line, thefirstpt, thelastpt, &mut the_param);
                }
                Some(fp) => {
                    for i in fp.lower()..=fp.upper() {
                        the_param.set(i + thefirstpt - 1, fp.get(i));
                    }
                }
            }
            self.the_multi_curve = MultiCurve::new(1, 0, 0);
            let mut an_other_line0 = MultiLine::new();
            let mut is_other_line0_made = false;
            let mut indbad = 0i32;
            self.alldone = self.compute(
                line,
                myfirstpt,
                mylastpt,
                &mut the_param,
                &mut thetol3d,
                &mut thetol2d,
                &mut indbad,
            );
            if indbad != 0 {
                is_other_line0_made = my_line_tool::make_ml_one_more_point(
                    line,
                    myfirstpt as usize,
                    mylastpt as usize,
                    indbad as usize,
                    &mut an_other_line0,
                );
            }
            if is_other_line0_made {
                self.my_is_clear = true;
                //++myMultiLineNb;
                self.perform(&an_other_line0);
                self.alldone = true;
            }
            if !self.alldone && self.the_multi_curve.nb_curves() > 0 {
                self.my_multi_curves.push(self.the_multi_curve.clone());
                self.tolers3d.push(self.currenttol3d);
                self.tolers2d.push(self.currenttol2d);
                let mylen = mylastpt - myfirstpt + 1;
                let my_par_len = self.my_parameters.as_ref().expect("myParameters").length();
                let a_len = if my_par_len > mylen { my_par_len } else { mylen };
                let mut the_par = RVector::new(myfirstpt, myfirstpt + a_len - 1);
                if let Some(mp) = self.my_parameters.as_ref() {
                    for i in 0..a_len {
                        the_par.set(myfirstpt + i, mp.get(mp.lower() + i));
                    }
                }
                self.my_par.push(the_par);
            }
        } else {
            // Cutting into troncons (gxx L916-1218).
            while !finish {
                oldlastpt = mylastpt;
                // Management of the multiline splitting for approximation:
                if !begin {
                    if !go_up {
                        if ok {
                            // Calcul de la partie a approximer.
                            myfirstpt = mylastpt;
                            mylastpt = thelastpt;
                            if myfirstpt == thelastpt {
                                finish = true;
                                self.alldone = true;
                                return;
                            }
                        } else {
                            nbp = mylastpt - myfirstpt + 1;
                            my_status = my_line_tool::what_status(
                                line,
                                myfirstpt as usize,
                                mylastpt as usize,
                            );
                            if my_status == ApproxStatus::NoPointsAdded
                                && nbp <= self.mydegremax + 1
                            {
                                let interpol = self.compute_curve(line, myfirstpt, mylastpt);
                                if interpol {
                                    if mylastpt == thelastpt {
                                        finish = true;
                                        self.alldone = true;
                                        return;
                                    }
                                }
                            }
                            mylastpt = (myfirstpt + mylastpt) / 2;
                        }
                    }
                    go_up = false;
                }

                // Verification du nombre de points restants par rapport au
                // degre demande.
                // ========================================================
                nbp = mylastpt - myfirstpt + 1;
                my_status =
                    my_line_tool::what_status(line, myfirstpt as usize, mylastpt as usize);
                if nbp <= self.mydegremax + 5 {
                    // Rajout necessaire de points si possible.
                    // ========================================
                    go_up = false;
                    ok = true;
                    if my_status == ApproxStatus::PointsAdded {
                        // Appel recursif du decoupage:
                        go_up = true;

                        let an_other_line1 = my_line_tool::make_ml_between(
                            line,
                            myfirstpt as usize,
                            mylastpt as usize,
                            (nbp - 1) as usize,
                        );

                        let nbpdsotherligne = my_line_tool::first_point(&an_other_line1) as i32
                            - my_line_tool::last_point(&an_other_line1) as i32;

                        // -- If MakeML failed, return an empty line
                        if nbpdsotherligne == 0 || self.my_multi_line_nb >= 3 {
                            if myfirstpt == mylastpt {
                                break; // Pour etre sur de ne pas
                                       // planter la station !!
                            }
                            my_couple1.index = myfirstpt;
                            my_couple2.index = mylastpt;
                            self.my_constraints[0] = my_couple1;
                            self.my_constraints[1] = my_couple2;

                            let mut param = RVector::new(myfirstpt, mylastpt);
                            let save_par = self.par;
                            self.par = ApproxParamType::IsoParametric;
                            self.parameters_compute(line, myfirstpt, mylastpt, &mut param);
                            self.the_multi_curve = MultiCurve::new(1, 0, 0);
                            let mut an_other_line2 = MultiLine::new();
                            let mut is_other_line2_made = false;
                            let mut indbad = 0i32;
                            ok = self.compute(
                                line,
                                myfirstpt,
                                mylastpt,
                                &mut param,
                                &mut thetol3d,
                                &mut thetol2d,
                                &mut indbad,
                            );
                            if indbad != 0 {
                                is_other_line2_made = my_line_tool::make_ml_one_more_point(
                                    line,
                                    myfirstpt as usize,
                                    mylastpt as usize,
                                    indbad as usize,
                                    &mut an_other_line2,
                                );
                            }
                            if is_other_line2_made {
                                self.my_is_clear = true;
                                //++myMultiLineNb;
                                self.par = save_par;
                                self.perform(&an_other_line2);
                                ok = true;
                            }

                            if !ok {
                                let tt3d = self.currenttol3d;
                                let tt2d = self.currenttol2d;
                                let save_parameters = self.my_parameters.clone();
                                let save_multi_curve = self.the_multi_curve.clone();

                                if save_par != ApproxParamType::IsoParametric {
                                    self.par = save_par;
                                } else {
                                    self.par = ApproxParamType::ChordLength;
                                }

                                self.parameters_compute(line, myfirstpt, mylastpt, &mut param);
                                is_other_line2_made = false;
                                indbad = 0;
                                ok = self.compute(
                                    line,
                                    myfirstpt,
                                    mylastpt,
                                    &mut param,
                                    &mut thetol3d,
                                    &mut thetol2d,
                                    &mut indbad,
                                );
                                if indbad != 0 {
                                    is_other_line2_made = my_line_tool::make_ml_one_more_point(
                                        line,
                                        myfirstpt as usize,
                                        mylastpt as usize,
                                        indbad as usize,
                                        &mut an_other_line2,
                                    );
                                }
                                if is_other_line2_made {
                                    self.my_is_clear = true;
                                    //++myMultiLineNb;
                                    self.perform(&an_other_line2);
                                    ok = true;
                                }

                                if !ok && tt3d <= self.currenttol3d && tt2d <= self.currenttol2d
                                {
                                    self.currenttol3d = tt3d;
                                    self.currenttol2d = tt2d;
                                    self.my_parameters = save_parameters;
                                    self.the_multi_curve = save_multi_curve;
                                }
                            }
                            self.par = save_par;
                            if myfirstpt == thelastpt {
                                finish = true;
                                self.alldone = true;
                                return;
                            }

                            oldlastpt = mylastpt;
                            if !ok {
                                self.tolreached = false;
                                if self.the_multi_curve.nb_curves() == 0 {
                                    self.my_multi_curves.clear();
                                    return;
                                }
                                let mut an_other_line3 = MultiLine::new();
                                let mut is_other_line3_made = false;
                                let mut indbad2 = 0i32;
                                if !check_multi_curve(
                                    &self.the_multi_curve,
                                    line,
                                    myfirstpt,
                                    mylastpt,
                                    &mut indbad2,
                                ) {
                                    is_other_line3_made = my_line_tool::make_ml_one_more_point(
                                        line,
                                        myfirstpt as usize,
                                        mylastpt as usize,
                                        indbad2 as usize,
                                        &mut an_other_line3,
                                    );
                                }
                                if is_other_line3_made {
                                    self.my_is_clear = true;
                                    //++myMultiLineNb;
                                    self.perform(&an_other_line3);
                                    myfirstpt = mylastpt;
                                    mylastpt = thelastpt;
                                } else {
                                    self.my_multi_curves.push(self.the_multi_curve.clone());
                                    self.tolers3d.push(self.currenttol3d);
                                    self.tolers2d.push(self.currenttol2d);
                                    let mylen = oldlastpt - myfirstpt + 1;
                                    let my_par_len = self
                                        .my_parameters
                                        .as_ref()
                                        .expect("myParameters")
                                        .length();
                                    let a_len = if my_par_len > mylen { my_par_len } else { mylen };
                                    let mut the_par = RVector::new(myfirstpt, myfirstpt + a_len - 1);
                                    if let Some(mp) = self.my_parameters.as_ref() {
                                        for i in 0..a_len {
                                            the_par.set(myfirstpt + i, mp.get(mp.lower() + i));
                                        }
                                    }
                                    self.my_par.push(the_par);
                                }
                            }
                            myfirstpt = oldlastpt;
                            mylastpt = thelastpt;
                        } else {
                            self.my_is_clear = true;
                            self.my_multi_line_nb += 1;
                            self.perform(&an_other_line1);
                            myfirstpt = mylastpt;
                            mylastpt = thelastpt;
                        }
                    }

                    if my_status == ApproxStatus::NoPointsAdded && !begin {
                        // On rend la meilleure approximation obtenue
                        // precedemment.
                        // ================================================
                        go_up = true;
                        self.tolreached = false;
                        if self.the_multi_curve.nb_curves() == 0 {
                            self.my_multi_curves.clear();
                            return;
                        }
                        self.my_multi_curves.push(self.the_multi_curve.clone());
                        self.tolers3d.push(self.currenttol3d);
                        self.tolers2d.push(self.currenttol2d);
                        let mylen = oldlastpt - myfirstpt + 1;
                        let my_par_len =
                            self.my_parameters.as_ref().expect("myParameters").length();
                        let a_len = if my_par_len > mylen { my_par_len } else { mylen };
                        let mut the_par = RVector::new(myfirstpt, myfirstpt + a_len - 1);
                        if let Some(mp) = self.my_parameters.as_ref() {
                            for i in 0..a_len {
                                the_par.set(myfirstpt + i, mp.get(mp.lower() + i));
                            }
                        }
                        self.my_par.push(the_par);

                        myfirstpt = oldlastpt;
                        mylastpt = thelastpt;
                    } else if my_status == ApproxStatus::NoApproximation {
                        // No approximation is done between myfirstpt and
                        // mylastpt. Store to be able to inform the user.
                        // ==================================================
                        go_up = true;
                        myfirstpt = mylastpt;
                        mylastpt = thelastpt;
                    }
                }

                if myfirstpt == thelastpt {
                    finish = true;
                    self.alldone = true;
                    return;
                }
                if !go_up {
                    if myfirstpt == mylastpt {
                        break; // Pour etre sur de ne pas
                               // planter la station !!
                    }
                    my_couple1.index = myfirstpt;
                    my_couple2.index = mylastpt;
                    self.my_constraints[0] = my_couple1;
                    self.my_constraints[1] = my_couple2;

                    // Calculation of parameters on this new interval.
                    // Retrieve the initial parameters during splitting.
                    let mut param = RVector::new(myfirstpt, mylastpt);
                    if begin {
                        match &self.myfirst_param {
                            None => {
                                self.parameters_compute(line, myfirstpt, mylastpt, &mut param);
                            }
                            Some(fp) => {
                                for i in fp.lower()..=fp.upper() {
                                    param.set(i, fp.get(i));
                                }
                            }
                        }
                        self.myfirst_param = None; // myfirstParam.Nullify();
                        the_param = param.clone(); // TheParam = Param;
                        begin = false;
                    } else {
                        let pfirst = the_param.get(myfirstpt);
                        let plast = the_param.get(mylastpt);
                        for i in myfirstpt..=mylastpt {
                            param.set(i, (the_param.get(i) - pfirst) / (plast - pfirst));
                        }
                    }

                    self.the_multi_curve = MultiCurve::new(1, 0, 0);
                    let mut indbad = 0i32;
                    ok = self.compute(
                        line,
                        myfirstpt,
                        mylastpt,
                        &mut param,
                        &mut thetol3d,
                        &mut thetol2d,
                        &mut indbad,
                    );
                    if myfirstpt == thelastpt {
                        finish = true;
                        self.alldone = true;
                        return;
                    }
                }
            }
        }
    }

    /// OCCT Parameters(Index) (gxx L1221-1224) — the parameters of the
    /// approximation corresponding to the multicurve Index (1-based).
    pub fn parameters(&self, index: usize) -> &RVector {
        &self.my_par[index - 1] // OCCT myPar.Value(Index)->Array1()
    }

    /// OCCT NbMultiCurves() (gxx L1226-1229).
    pub fn nb_multi_curves(&self) -> usize {
        self.my_multi_curves.len()
    }

    /// OCCT ChangeValue(Index) (gxx L1231-1234).
    pub fn change_value(&mut self, index: usize) -> &mut MultiCurve {
        &mut self.my_multi_curves[index - 1]
    }

    /// OCCT Value(Index) (gxx L1236-1239).
    pub fn value(&self, index: usize) -> &MultiCurve {
        &self.my_multi_curves[index - 1]
    }

    /// OCCT SplineValue() (gxx L1241-1247) — concatenates the Bezier
    /// multicurves into one MultiBSpCurve.
    pub fn spline_value(&mut self) -> &MultiBSpCurve {
        // OCCT: Approx_MCurvesToBSpCurve Trans; Trans.Perform(myMultiCurves);
        //       myspline = Trans.Value();
        // The rcad MCurvesToBSpCurve receives the curves via append() before
        // perform().
        let mut trans = MCurvesToBSpCurve::new();
        for cu in &self.my_multi_curves {
            trans.append(cu.clone());
        }
        trans.perform();
        self.myspline = trans.value();
        &self.myspline
    }

    /// OCCT Parameters(Line, firstP, LastP, TheParameters) (gxx L1249-1322).
    fn parameters_compute(
        &self,
        line: &MultiLine,
        first_p: i32,
        last_p: i32,
        the_parameters: &mut RVector,
    ) {
        if self.par == ApproxParamType::ChordLength || self.par == ApproxParamType::Centripetal {
            let nb_p3d = my_line_tool::nb_p3d(line) as i32;
            let nb_p2d = my_line_tool::nb_p2d(line) as i32;
            let mut mynb_p3d = nb_p3d;
            let mut mynb_p2d = nb_p2d;
            if nb_p3d == 0 {
                mynb_p3d = 1;
            }
            if nb_p2d == 0 {
                mynb_p2d = 1;
            }

            the_parameters.set(first_p, 0.0);
            let mut dist = 0.0f64;
            let mut tab_p = vec![DVec3::ZERO; mynb_p3d as usize];
            let mut tab_pp = vec![DVec3::ZERO; mynb_p3d as usize];
            let mut tab_p2d = vec![DVec2::ZERO; mynb_p2d as usize];
            let mut tab_pp2d = vec![DVec2::ZERO; mynb_p2d as usize];

            for i in (first_p + 1)..=last_p {
                if nb_p3d != 0 && nb_p2d != 0 {
                    my_line_tool::value_3d_2d(line, (i - 1) as usize, &mut tab_p, &mut tab_p2d);
                } else if nb_p2d != 0 {
                    my_line_tool::value_2d(line, (i - 1) as usize, &mut tab_p2d);
                } else if nb_p3d != 0 {
                    my_line_tool::value_3d(line, (i - 1) as usize, &mut tab_p);
                }

                if nb_p3d != 0 && nb_p2d != 0 {
                    my_line_tool::value_3d_2d(line, i as usize, &mut tab_pp, &mut tab_pp2d);
                } else if nb_p2d != 0 {
                    my_line_tool::value_2d(line, i as usize, &mut tab_pp2d);
                } else if nb_p3d != 0 {
                    my_line_tool::value_3d(line, i as usize, &mut tab_pp);
                }
                dist = 0.0;
                for j in 1..=nb_p3d {
                    let a_p1 = tab_p[(j - 1) as usize];
                    let a_p2 = tab_pp[(j - 1) as usize];
                    dist += a_p2.distance_squared(a_p1);
                }
                for j in 1..=nb_p2d {
                    let a_p12d = tab_p2d[(j - 1) as usize];
                    let a_p22d = tab_pp2d[(j - 1) as usize];
                    dist += a_p22d.distance_squared(a_p12d);
                }

                dist = dist.sqrt();
                if self.par == ApproxParamType::ChordLength {
                    the_parameters.set(i, the_parameters.get(i - 1) + dist);
                } else {
                    // Par == Approx_Centripetal
                    the_parameters.set(i, the_parameters.get(i - 1) + dist.sqrt());
                }
            }
            for i in first_p..=last_p {
                let v = the_parameters.get(i) / the_parameters.get(last_p);
                the_parameters.set(i, v);
            }
        } else {
            for i in first_p..=last_p {
                the_parameters
                    .set(i, (i as f64 - first_p as f64) / (last_p as f64 - first_p as f64));
            }
        }
    }

    /// OCCT Compute(Line, fpt, lpt, Para, TheTol3d, TheTol2d, indbad)
    /// (gxx L1324-1441).
    #[allow(clippy::too_many_arguments)]
    fn compute(
        &mut self,
        line: &MultiLine,
        fpt: i32,
        lpt: i32,
        para: &mut RVector,
        the_tol3d: &mut f64,
        the_tol2d: &mut f64,
        indbad: &mut i32,
    ) -> bool {
        *indbad = 0;
        let mut mydone;
        // OCCT Fv is an uninitialized local written by SQ.Error.
        let mut fv = 0.0f64;
        let nbp = lpt - fpt + 1;

        let mut par_sav = RVector::new(para.lower(), para.upper());
        for i in para.lower()..=para.upper() {
            par_sav.set(i, para.get(i));
        }
        let mut mdegmax = self.mydegremax;
        if nbp < mdegmax + 5 && self.mycut {
            mdegmax = nbp - 5;
        }
        if mdegmax < self.mydegremin {
            mdegmax = self.mydegremin;
        }

        self.currenttol3d = REAL_LAST;
        self.currenttol2d = REAL_LAST;
        let mut deg = std::cmp::min(nbp - 1, self.mydegremin);
        while deg <= mdegmax {
            let mut my_scu = MultiCurve::new((deg + 1) as usize, 0, 0);
            if self.mysquares {
                // OCCT passes Para as const math_Vector& — the rcad
                // LeastSquare takes the flat 1-based copy.
                let a_params = rvector_to_vecd(para);
                let mut sq = LeastSquare::new(
                    line,
                    fpt,
                    lpt,
                    self.myfirstc,
                    self.mylastc,
                    &a_params,
                    deg + 1,
                );
                mydone = sq.is_done();
                my_scu = sq.bezier_value();
                sq.error(&mut fv, the_tol3d, the_tol2d);
            } else {
                // OCCT passes Para as math_Vector& — GRAD updates it in
                // place; the rcad Gradient takes the flat 1-based copy and
                // the optimized parameters are written back into Para.
                let mut a_params = rvector_to_vecd(para);
                let mut grad = Gradient::new(
                    line,
                    fpt,
                    lpt,
                    &self.my_constraints,
                    &mut a_params,
                    deg,
                    self.mytol3d,
                    self.mytol2d,
                    self.myitermax,
                );
                for i in fpt..=lpt {
                    para.set(i, a_params.get((i - fpt + 1) as usize));
                }
                mydone = grad.is_done();
                my_scu = grad.value();
                if my_scu.nb_curves() == 0 {
                    deg += 1;
                    continue;
                }
                *the_tol3d = grad.max_error_3d();
                *the_tol2d = grad.max_error_2d();
            }
            let mut uu1 = para.get(para.lower());
            let mut restau = false;
            for i in (para.lower() + 1)..=para.upper() {
                let uu2 = para.get(i);
                if uu2 <= uu1 {
                    restau = true;
                    break;
                }
                uu1 = uu2;
            }
            if restau {
                for i in para.lower()..=para.upper() {
                    para.set(i, par_sav.get(i));
                }
            }
            if mydone {
                if *the_tol3d <= self.mytol3d && *the_tol2d <= self.mytol2d {
                    // Stockage de la multicurve approximee.
                    self.tolreached = true;
                    if !check_multi_curve(&my_scu, line, fpt, lpt, indbad) {
                        return false;
                    } else {
                        self.my_multi_curves.push(my_scu.clone());
                        // Storage of parameters of the approximated
                        // MultiLine part: A ameliorer !! (bq trop de
                        // recopies)
                        let mut the_par = RVector::new(para.lower(), para.upper());
                        for i in para.lower()..=para.upper() {
                            the_par.set(i, para.get(i));
                        }
                        self.my_par.push(the_par);
                        self.tolers3d.push(*the_tol3d);
                        self.tolers2d.push(*the_tol2d);
                        return true;
                    }
                }
            }

            if *the_tol3d <= self.currenttol3d && *the_tol2d <= self.currenttol2d {
                self.the_multi_curve = my_scu;
                self.currenttol3d = *the_tol3d;
                self.currenttol2d = *the_tol2d;
                let mut mp = RVector::new(para.lower(), para.upper());
                for i in para.lower()..=para.upper() {
                    mp.set(i, para.get(i));
                }
                self.my_parameters = Some(mp);
            }
            deg += 1;
        }

        false
    }

    /// OCCT ComputeCurve(Line, firstpt, lastpt) (gxx L1443-1690) — the
    /// interpolation fallback for a troncon with too few points to cut
    /// further.
    fn compute_curve(&mut self, line: &MultiLine, firstpt: i32, lastpt: i32) -> bool {
        let mut mydone = false;
        let myfirstpt = firstpt;
        let mylastpt = lastpt;
        let nbp = lastpt - firstpt + 1;
        let mut para = RVector::new(firstpt, lastpt);

        self.parameters_compute(line, firstpt, lastpt, &mut para);

        let nb_p3d = my_line_tool::nb_p3d(line) as i32;
        let nb_p2d = my_line_tool::nb_p2d(line) as i32;
        let mut mynb_p3d = nb_p3d;
        let mut mynb_p2d = nb_p2d;
        if nb_p3d == 0 {
            mynb_p3d = 1;
        }
        if nb_p2d == 0 {
            mynb_p2d = 1;
        }

        let mut tab_v1 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_v2 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p1 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p2 = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_p = vec![DVec3::ZERO; mynb_p3d as usize];
        let mut tab_v12d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_v22d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_p12d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_p22d = vec![DVec2::ZERO; mynb_p2d as usize];
        let mut tab_p2d = vec![DVec2::ZERO; mynb_p2d as usize];

        let tangent1;
        let tangent2;
        if nb_p3d != 0 && nb_p2d != 0 {
            my_line_tool::value_3d_2d(line, myfirstpt as usize, &mut tab_p1, &mut tab_p12d);
            my_line_tool::value_3d_2d(line, mylastpt as usize, &mut tab_p2, &mut tab_p22d);
            tangent1 =
                my_line_tool::tangency_3d_2d(line, myfirstpt as usize, &mut tab_v1, &mut tab_v12d);
            tangent2 =
                my_line_tool::tangency_3d_2d(line, mylastpt as usize, &mut tab_v2, &mut tab_v22d);
        } else if nb_p2d != 0 {
            my_line_tool::value_2d(line, myfirstpt as usize, &mut tab_p12d);
            my_line_tool::value_2d(line, mylastpt as usize, &mut tab_p22d);
            tangent1 = my_line_tool::tangency_2d(line, myfirstpt as usize, &mut tab_v12d);
            tangent2 = my_line_tool::tangency_2d(line, mylastpt as usize, &mut tab_v22d);
        } else {
            my_line_tool::value_3d(line, myfirstpt as usize, &mut tab_p1);
            my_line_tool::value_3d(line, mylastpt as usize, &mut tab_p2);
            tangent1 = my_line_tool::tangency_3d(line, myfirstpt as usize, &mut tab_v1);
            tangent2 = my_line_tool::tangency_3d(line, mylastpt as usize, &mut tab_v2);
        }
        if nbp == 2 {
            // If there are only 2 points, we still verify that the tangents
            // are aligned (gxx L1494-1563).
            if tangent1 {
                for i in 1..=nb_p3d {
                    let p_vec = tab_p2[(i - 1) as usize] - tab_p1[(i - 1) as usize];
                    let v13d = tab_v1[(i - 1) as usize];
                    if !gp_vec_is_parallel(p_vec, v13d, ANGULAR) {
                        break;
                    }
                }
                for i in 1..=nb_p2d {
                    let p_vec2d = tab_p22d[(i - 1) as usize] - tab_p12d[(i - 1) as usize];
                    let v12d = tab_v12d[(i - 1) as usize];
                    if !gp_vec2d_is_parallel(p_vec2d, v12d, ANGULAR) {
                        break;
                    }
                }
            }

            if tangent2 {
                for i in 1..=nb_p3d {
                    let p_vec = tab_p2[(i - 1) as usize] - tab_p1[(i - 1) as usize];
                    let v23d = tab_v2[(i - 1) as usize];
                    if !gp_vec_is_parallel(p_vec, v23d, ANGULAR) {
                        break;
                    }
                }
                for i in 1..=nb_p2d {
                    let p_vec2d = tab_p22d[(i - 1) as usize] - tab_p12d[(i - 1) as usize];
                    let v22d = tab_v22d[(i - 1) as usize];
                    if !gp_vec2d_is_parallel(p_vec2d, v22d, ANGULAR) {
                        break;
                    }
                }
            }

            let mut my_scu = MultiCurve::new((self.mydegremin + 1) as usize, 0, 0);
            if nb_p3d != 0 && nb_p2d != 0 {
                let m_pole1 = MultiPoint::new_tab_p3d_p2d(&tab_p1, &tab_p12d);
                let m_pole2 = MultiPoint::new_tab_p3d_p2d(&tab_p2, &tab_p22d);
                my_scu.set_value(1, m_pole1);
                my_scu.set_value((self.mydegremin + 1) as usize, m_pole2);
                for i in 2..=self.mydegremin {
                    for j in 1..=nb_p3d {
                        let p1 = tab_p1[(j - 1) as usize];
                        let p2 = tab_p2[(j - 1) as usize];
                        tab_p[(j - 1) as usize] =
                            p1 + (i - 1) as f64 * (p2 - p1) / self.mydegremin as f64;
                    }
                    for j in 1..=nb_p2d {
                        let p12d = tab_p12d[(j - 1) as usize];
                        let p22d = tab_p22d[(j - 1) as usize];
                        tab_p2d[(j - 1) as usize] =
                            p12d + (i - 1) as f64 * (p22d - p12d) / self.mydegremin as f64;
                    }
                    let m_pole = MultiPoint::new_tab_p3d_p2d(&tab_p, &tab_p2d);
                    my_scu.set_value(i as usize, m_pole);
                }
            } else if nb_p3d != 0 {
                let m_pole1 = MultiPoint::new_tab_p3d(&tab_p1);
                let m_pole2 = MultiPoint::new_tab_p3d(&tab_p2);
                my_scu.set_value(1, m_pole1);
                my_scu.set_value((self.mydegremin + 1) as usize, m_pole2);
                for i in 2..=self.mydegremin {
                    for j in 1..=nb_p3d {
                        let p1 = tab_p1[(j - 1) as usize];
                        let p2 = tab_p2[(j - 1) as usize];
                        tab_p[(j - 1) as usize] =
                            p1 + (i - 1) as f64 * (p2 - p1) / self.mydegremin as f64;
                    }
                    let m_pole = MultiPoint::new_tab_p3d(&tab_p);
                    my_scu.set_value(i as usize, m_pole);
                }
            } else if nb_p2d != 0 {
                let m_pole1 = MultiPoint::new_tab_p2d(&tab_p12d);
                let m_pole2 = MultiPoint::new_tab_p2d(&tab_p22d);
                my_scu.set_value(1, m_pole1);
                my_scu.set_value((self.mydegremin + 1) as usize, m_pole2);
                for i in 2..=self.mydegremin {
                    for j in 1..=nb_p2d {
                        let p12d = tab_p12d[(j - 1) as usize];
                        let p22d = tab_p22d[(j - 1) as usize];
                        tab_p2d[(j - 1) as usize] =
                            p12d + (i - 1) as f64 * (p22d - p12d) / self.mydegremin as f64;
                    }
                    let m_pole = MultiPoint::new_tab_p2d(&tab_p2d);
                    my_scu.set_value(i as usize, m_pole);
                }
            }
            mydone = true;
            // Stockage de la multicurve approximee.
            self.tolreached = true;
            self.my_multi_curves.push(my_scu);
            let mut the_par = RVector::new(para.lower(), para.upper());
            for i in para.lower()..=para.upper() {
                the_par.set(i, para.get(i));
            }
            self.my_par.push(the_par);
            self.tolers3d.push(CONFUSION);
            self.tolers2d.push(PCONFUSION);
            return mydone;
        }

        // with the tangents (gxx L1645-1690).
        let deg = nbp + 1;
        let mut my_scu = MultiCurve::new((deg + 1) as usize, 0, 0);
        let cons = AppParConstraint::TangencyPoint;
        let mut v1 = RVector::new(1, nb_p3d * 3 + nb_p2d * 2);
        let mut v2 = RVector::new(1, nb_p3d * 3 + nb_p2d * 2);
        self.first_tangency_vector(line, myfirstpt, &mut v1);
        let mut lambda1 = self.search_first_lambda(line, &para, &v1, myfirstpt);

        self.last_tangency_vector(line, mylastpt, &mut v2);
        let mut lambda2 = self.search_last_lambda(line, &para, &v2, mylastpt);

        let mut lsq = LeastSquare::new(
            line,
            myfirstpt,
            mylastpt,
            cons,
            cons,
            &rvector_to_vecd(&para),
            deg + 1,
        );

        lambda1 /= deg as f64;
        lambda2 /= deg as f64;
        lsq.perform_v1tv2t(
            &rvector_to_vecd(&para),
            &rvector_to_vecd(&v1),
            &rvector_to_vecd(&v2),
            lambda1,
            lambda2,
        );
        mydone = lsq.is_done();
        my_scu = lsq.bezier_value();

        if mydone {
            let mut fv = 0.0f64;
            let mut the_tol3d = 0.0f64;
            let mut the_tol2d = 0.0f64;
            lsq.error(&mut fv, &mut the_tol3d, &mut the_tol2d);

            // Stockage de la multicurve approximee.
            self.tolreached = true;
            self.my_multi_curves.push(my_scu);
            let mut the_par = RVector::new(para.lower(), para.upper());
            for i in para.lower()..=para.upper() {
                the_par.set(i, para.get(i));
            }
            self.my_par.push(the_par);
            self.tolers3d.push(the_tol3d);
            self.tolers2d.push(the_tol2d);
            return true;
        }
        mydone
    }

    /// OCCT Init(degreemin..Squares) (gxx L1692-1709).
    #[allow(clippy::too_many_arguments)]
    pub fn init(
        &mut self,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) {
        self.mydegremin = degreemin;
        self.mydegremax = degreemax;
        self.mytol3d = tolerance3d;
        self.mytol2d = tolerance2d;
        self.par = parametrization;
        self.mysquares = squares;
        self.mycut = cutting;
        self.myitermax = nb_iterations;
    }

    /// OCCT SetDegrees(degreemin, degreemax) (gxx L1711-1715).
    pub fn set_degrees(&mut self, degreemin: i32, degreemax: i32) {
        self.mydegremin = degreemin;
        self.mydegremax = degreemax;
    }

    /// OCCT SetTolerances(Tolerance3d, Tolerance2d) (gxx L1717-1721).
    pub fn set_tolerances(&mut self, tolerance3d: f64, tolerance2d: f64) {
        self.mytol3d = tolerance3d;
        self.mytol2d = tolerance2d;
    }

    /// OCCT SetConstraints(FirstC, LastC) (gxx L1723-1728).
    pub fn set_constraints(&mut self, first_c: AppParConstraint, last_c: AppParConstraint) {
        self.myfirstc = first_c;
        self.mylastc = last_c;
    }

    /// OCCT IsAllApproximated() (gxx L1730-1733).
    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }

    /// OCCT IsToleranceReached() (gxx L1735-1738).
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }

    /// OCCT Error(Index, tol3d, tol2d) (gxx L1740-1744) — the tolerances 2d
    /// and 3d of the Index MultiCurve, returned as a pair.
    pub fn error(&self, index: usize) -> (f64, f64) {
        (self.tolers3d[index - 1], self.tolers2d[index - 1])
    }

    /// OCCT Parametrization() (gxx L1746-1749).
    pub fn parametrization(&self) -> ApproxParamType {
        self.par
    }
}

/// OCCT math_Vector -> VecD conversion for the LeastSquare/Gradient calls
/// (the rcad solvers take the flat 1-based copy: VecD(k) = Para(lower+k-1)).
fn rvector_to_vecd(v: &RVector) -> VecD {
    let mut r = VecD::new(v.length() as usize);
    for i in 1..=v.length() {
        r.set(i as usize, v.get(i));
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT anchor: a quarter unit circle sampled at 17 points approximated
    /// with cutting (Perform, gxx L832-1219) — the splitting loop terminates
    /// with alldone and one stored MultiCurve covering the samples.
    #[test]
    fn compute_quarter_circle() {
        let n = 17;
        let pts: Vec<DVec3> = (0..n)
            .map(|i| {
                let t = 0.5 * std::f64::consts::FRAC_PI_2 * (i as f64) / ((n - 1) as f64);
                DVec3::new(t.cos(), t.sin(), 0.0)
            })
            .collect();
        let line = MultiLine::new_tab_p3d(&pts);

        let mut c = Compute::new(
            3,
            8,
            1.0e-3,
            1.0e-6,
            5,
            true,
            ApproxParamType::ChordLength,
            false,
        );
        c.perform(&line);

        assert!(c.is_all_approximated(), "alldone");
        assert_eq!(c.nb_multi_curves(), 1);
        let cu_nb_poles = c.value(1).nb_poles();
        // The per-pole multipoint coordinate count (1 = 3D-only line).
        assert_eq!(c.value(1).nb_curves(), 1);
        // The stored parameters of the approximated part exist (myPar).
        let pars = c.parameters(1);
        assert_eq!(pars.length(), n as i32);
        // The spline concatenation of the single Bezier piece.
        let spl = c.spline_value();
        assert_eq!(spl.poles.len(), cu_nb_poles);
    }

    /// OCCT anchor: a two-point MultiLine with NoConstraint ends is fitted
    /// exactly by Compute's degree loop (gxx L1354-1424) — at deg =
    /// min(nbp-1, mydegremin) = 1 the GRAD fit (nbpoles = 2, Householder
    /// branch, gxx L536-545) interpolates the two points with zero error,
    /// Compute stores the curve and returns true (gxx L1394-1424), and
    /// Perform terminates with alldone at the myfirstpt == Thelastpt check
    /// (gxx L1160-1164). ComputeCurve is not reached for this input.
    ///
    /// Note on the constrained (default TangencyPoint) variant of this
    /// input: OCCT's deg=2 GRAD would run AppParCurves_LeastSquare with
    /// FirstP = TheFirstPoint(TangencyPoint, 1) = 2 > LastP =
    /// TheLastPoint(TangencyPoint, 2) = 1 (gxx L1465-1483), hence
    /// Nlignes = NA * Neq = 0 (gxx L500-503) and MakeTAA constructing the
    /// EMPTY math_Vector myB(myfirst, mylast) with mylast = myfirst - 1
    /// (gxx L1532-1542) — legal in OCCT (math_VectorBase.lxx L36-50 builds
    /// a length-0 array) and every consuming loop (FirstP..=LastP) is empty.
    /// The rcad kernel math_Vector asserts i2 >= i1 instead
    /// (math_matrix.rs:208), so that variant panics in rcad — a kernel-level
    /// divergence outside this file's scope; the test uses NoConstraint,
    /// whose FirstP/LastP (gxx L1465-1471 NoConstraint branch) keep the
    /// ranges well-formed in both implementations.
    #[test]
    fn compute_two_points_exact_fit() {
        let pts = vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 2.0, 3.0)];
        let line = MultiLine::new_tab_p3d(&pts);

        let mut c = Compute::new(
            3,
            8,
            1.0e-3,
            1.0e-6,
            5,
            true,
            ApproxParamType::ChordLength,
            false,
        );
        c.set_constraints(AppParConstraint::NoConstraint, AppParConstraint::NoConstraint);
        c.perform(&line);

        assert!(c.is_all_approximated(), "alldone");
        assert!(c.is_tolerance_reached(), "tolreached");
        let (tol3d, tol2d) = c.error(1);
        assert!(tol3d <= 1.0e-12, "tol3d={}", tol3d);
        assert!(tol2d <= 1.0e-12, "tol2d={}", tol2d);
        let cu = c.value(1);
        // nbpoles = deg + 1 = 2 at the first loop degree (gxx L1354-1356).
        assert_eq!(cu.nb_poles(), 2, "nbpoles");
        let p0 = cu.poles[0].point(1);
        let p1 = cu.poles[1].point(1);
        assert!(p0.distance(pts[0]) < 1.0e-12, "p0={:?}", p0);
        assert!(p1.distance(pts[1]) < 1.0e-12, "p1={:?}", p1);
    }
}
