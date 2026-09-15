//! OCCT AppDef_Variational (ModelingData/TKGeomBase/AppDef).
//!
//! 1:1 translation of `AppDef_Variational.hxx` (L17-330) and
//! `AppDef_Variational.cxx` (L1-3403).  This file carries the class
//! declaration mirror, the constructor (cxx L71-219), `Init` (L224-484),
//! `Approximate` (L489-660), the accessors / `Dump` (L665-947) and the
//! setters (L951-1156).  The private engine methods are split into
//! `#[path]` child modules:
//! - [`variational_b`]: TheMotor / Optimization / Project / ACR /
//!   NearIndex / GettingKnots / SplitCurve (cxx L1160-2126),
//! - [`variational_c`]: InitSmoothCriterion / InitParameters /
//!   InitCriterionEstimations / EstTangent / EstSecnd / InitCutting /
//!   Adjusting / AssemblingConstraints / InitTthetaF (cxx L2129-3403).
//!
//! Handle / architecture mapping:
//! - `occ::handle<AppDef_SmoothCriterion>` ->
//!   [`CriterionHandle`] (`Rc<RefCell<dyn SmoothCriterion>>`; the concrete
//!   `AppDef_LinearCriteria` lives behind it, mirroring the CurveHandle
//!   convention of app_def_smooth_criterion.rs);
//! - `occ::handle<FEmTool_Curve>` -> [`CurveHandle`]
//!   (`Rc<RefCell<Curve>>`, shared mutable aliasing preserved);
//! - `occ::handle<NCollection_HArray1<double>>` -> [`RealArray1`] by value
//!   (the OCCT Array1 assignment `ChangeArray1() = Array1()` copies bounds
//!   and data - NCollection_Array1.hxx L273/L325 - so the Rust clone is the
//!   exact counterpart);
//! - `occ::handle<NCollection_HArray1<int>>` ->
//!   `rcad_kernel::math::fem_tool::IntArray1`;
//! - `occ::handle<NCollection_HArray1<AppParCurves_ConstraintCouple>>` ->
//!   the local [`CoupleArray1`] mirror;
//! - `gp_Pnt / gp_Pnt2d / gp_Vec / gp_Vec2d` -> `glam::DVec3 / DVec2` (the
//!   gp free behaviors used here - Vec2d::Angle, Vec::CrossMagnitude,
//!   Vec::IsNormal, Dir::Angle - are mirrored by the module-private helpers
//!   [`gp_vec2d_angle`], [`gp_cross_magnitude`], [`gp_dir_angle`],
//!   [`gp_vec_is_normal`], [`not_parallel`]);
//! - the OCCT `math_Vector` algebraic operators are mirrored by the
//!   module-private helpers [`v_sub`], [`v_add`], [`v_mul`], [`v_norm`],
//!   [`v_norm2`], [`v_div_assign`] (math_VectorBase.lxx semantics).

use std::cell::RefCell;
use std::rc::Rc;

use glam::{DVec2, DVec3};

use rcad_kernel::math::convert_comp_polynomial_to_poles::ConvertCompPolynomialToPoles;
use rcad_kernel::math::fem_tool::IntArray1;
use rcad_kernel::math::math_matrix::{Matrix, Vector};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::app_def::my_line_tool;
use crate::geomalgo::app_def::MultiLine;
use crate::geomalgo::app_def_linear_criteria::LinearCriteria;
use crate::geomalgo::app_def_smooth_criterion::{CurveHandle, RealArray1, SmoothCriterion};
use crate::geomalgo::approx_int::{AppParConstraint, ConstraintCouple, MultiBSpCurve, MultiPoint};

#[path = "app_def_variational_b.rs"]
mod variational_b;

#[path = "app_def_variational_c.rs"]
mod variational_c;

/// OCCT `occ::handle<AppDef_SmoothCriterion>` - shared, mutable criterion
/// (the concrete `AppDef_LinearCriteria` lives behind it).
pub type CriterionHandle = Rc<RefCell<dyn SmoothCriterion>>;

/// Mirror of `NCollection_HArray1<AppParCurves_ConstraintCouple>`: the
/// constraint couples with an arbitrary (1-based in practice) lower bound.
#[derive(Debug, Clone)]
pub struct CoupleArray1 {
    lower: i32,
    data: Vec<ConstraintCouple>,
}

impl CoupleArray1 {
    /// OCCT new NCollection_HArray1<AppParCurves_ConstraintCouple>(Lower,
    /// Upper).
    pub fn new(lower: i32, upper: i32) -> Self {
        assert!(upper >= lower - 1, "NCollection_HArray1: bad range");
        CoupleArray1 {
            lower,
            data: vec![
                ConstraintCouple {
                    index: 0,
                    constraint: AppParConstraint::NoConstraint,
                };
                (upper - lower + 1).max(0) as usize
            ],
        }
    }

    /// OCCT Lower().
    #[inline]
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// OCCT Upper().
    #[inline]
    pub fn upper(&self) -> i32 {
        self.lower + self.data.len() as i32 - 1
    }

    /// OCCT Length().
    #[inline]
    pub fn length(&self) -> i32 {
        self.data.len() as i32
    }

    /// OCCT Value(Index).
    #[inline]
    pub fn value(&self, index: i32) -> ConstraintCouple {
        self.data[(index - self.lower) as usize]
    }

    /// OCCT SetValue(Index, Value).
    #[inline]
    pub fn set_value(&mut self, index: i32, value: ConstraintCouple) {
        self.data[(index - self.lower) as usize] = value;
    }

    /// OCCT ChangeValue(Index).
    #[inline]
    pub fn change_value(&mut self, index: i32) -> &mut ConstraintCouple {
        &mut self.data[(index - self.lower) as usize]
    }
}

// ---------------------------------------------------------------------------
// gp / math_Vector helper mirrors (pure utilities over glam / rcad Vector)
// ---------------------------------------------------------------------------

/// OCCT gp_Vec2d::Angle (gp_Vec2d.cxx L47-88) - angle between two 2d vectors
/// in [-PI, PI] (acos near +/-90 deg, asin near 0/180 deg).
fn gp_vec2d_angle(v: DVec2, other: DVec2) -> f64 {
    use std::f64::consts::FRAC_1_SQRT_2;
    let a_norm = v.length();
    let other_norm = other.length();
    assert!(
        a_norm > 0.0 && other_norm > 0.0,
        "gp_VectorWithNullMagnitude: gp_Vec2d::Angle"
    );
    let d = a_norm * other_norm;
    let cosinus = v.dot(other) / d;
    let sinus = (v.x * other.y - v.y * other.x) / d;

    if cosinus > -FRAC_1_SQRT_2 && cosinus < FRAC_1_SQRT_2 {
        // For angles near +/-90 degrees, use acos for better precision.
        if sinus > 0.0 {
            cosinus.acos()
        } else {
            -cosinus.acos()
        }
    } else {
        // For angles near 0 degrees or +/-180 degrees, use asin for better
        // precision.
        if cosinus > 0.0 {
            sinus.asin()
        } else if sinus > 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            -std::f64::consts::PI - sinus.asin()
        }
    }
}

/// OCCT gp_Vec::CrossMagnitude (gp_Vec.hxx L270-273, gp_XYZ::CrossMagnitude)
/// - the magnitude of the cross product || V ^ R ||.
fn gp_cross_magnitude(v: DVec3, right: DVec3) -> f64 {
    v.cross(right).length()
}

/// OCCT gp_Dir::Angle (gp_Dir.hxx) - ACos of the dot product of the two
/// (normalized) directions.
fn gp_dir_angle(a: DVec3, b: DVec3) -> f64 {
    a.normalize().dot(b.normalize()).acos()
}

/// OCCT gp_Vec::IsNormal (gp_Vec.hxx L480-484) - Abs(PI/2 - Angle) <=
/// AngularTolerance.
fn gp_vec_is_normal(v: DVec3, other: DVec3, angular_tolerance: f64) -> bool {
    let an_ang = (std::f64::consts::FRAC_PI_2 - gp_dir_angle(v, other)).abs();
    an_ang <= angular_tolerance
}

/// OCCT static NotParallel (AppDef_Variational.cxx L2976-2991) - returns
/// true as soon as V (built from T by unit coordinate shifts) is not
/// parallel to T.
fn not_parallel(t: DVec3, v: &mut DVec3) -> bool {
    *v = t;
    v.x += 1.0;
    if gp_cross_magnitude(*v, t) > 1.0e-12 {
        return true;
    }
    *v = t;
    v.y += 1.0;
    if gp_cross_magnitude(*v, t) > 1.0e-12 {
        return true;
    }
    *v = t;
    v.z += 1.0;
    gp_cross_magnitude(*v, t) > 1.0e-12
}

/// OCCT math_VectorBase operator- (math_VectorBase.lxx) - a - b.
fn v_sub(a: &Vector, b: &Vector) -> Vector {
    let mut r = Vector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set(i, a.get(i) - b.get(i));
    }
    r
}

/// OCCT math_VectorBase operator+ (math_VectorBase.lxx) - a + b.
fn v_add(a: &Vector, b: &Vector) -> Vector {
    let mut r = Vector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set(i, a.get(i) + b.get(i));
    }
    r
}

/// OCCT math_VectorBase operator*(scalar) (math_VectorBase.lxx) - a * s.
fn v_mul(a: &Vector, s: f64) -> Vector {
    let mut r = Vector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set(i, a.get(i) * s);
    }
    r
}

/// OCCT math_Vector::Norm (math_VectorBase.lxx) - sqrt(Norm2()).
fn v_norm(a: &Vector) -> f64 {
    v_norm2(a).sqrt()
}

/// OCCT math_Vector::Norm2 (math_VectorBase.lxx) - sum of squares.
fn v_norm2(a: &Vector) -> f64 {
    let mut s = 0.0;
    for i in a.lower()..=a.upper() {
        s += a.get(i) * a.get(i);
    }
    s
}

/// OCCT math_VectorBase operator/=(scalar) (math_VectorBase.lxx).
fn v_div_assign(v: &mut Vector, s: f64) {
    for i in v.lower()..=v.upper() {
        let x = v.get(i) / s;
        v.set(i, x);
    }
}

/// Build a RealArray1 with lower bound 1 from a data slice (the mirror of
/// `new NCollection_HArray1<double>(Array1)` - OCCT copies the array with
/// its bounds; every construction site here uses the 1-based storage).
fn real_array1_from(data: &[f64]) -> RealArray1 {
    let mut a = RealArray1::new(1, data.len() as i32);
    for (k, v) in data.iter().enumerate() {
        a.set_value(1 + k as i32, *v);
    }
    a
}

/// std::stable_sort over the whole RealArray1 storage (the mirror of
/// `std::stable_sort(CurrentTi->begin(), CurrentTi->end())`; RealArray1
/// owns its storage by value, so the sort round-trips through a Vec).
fn stable_sort_real_array(a: &mut RealArray1) {
    let mut v: Vec<f64> = (a.lower()..=a.upper()).map(|i| a.value(i)).collect();
    v.sort_by(|x, y| x.partial_cmp(y).unwrap());
    for (k, val) in v.iter().enumerate() {
        a.set_value(a.lower() + k as i32, *val);
    }
}

// ---------------------------------------------------------------------------
// AppDef_Variational
// ---------------------------------------------------------------------------

/// OCCT AppDef_Variational - smooths N points with constraints by
/// minimization of quadratic criterium but also variational criterium in
/// order to obtain "fair Curve" (AppDef_Variational.hxx L45-328).
pub struct AppDefVariational {
    /// OCCT mySSP.
    my_ssp: MultiLine,
    /// OCCT myNbP3d.
    my_nb_p3d: i32,
    /// OCCT myNbP2d.
    my_nb_p2d: i32,
    /// OCCT myDimension.
    my_dimension: i32,
    /// OCCT myFirstPoint.
    my_first_point: i32,
    /// OCCT myLastPoint.
    my_last_point: i32,
    /// OCCT myNbPoints.
    my_nb_points: i32,
    /// OCCT myTabPoints.
    my_tab_points: RealArray1,
    /// OCCT myConstraints.
    my_constraints: CoupleArray1,
    /// OCCT myNbConstraints.
    my_nb_constraints: i32,
    /// OCCT myTabConstraints.
    my_tab_constraints: RealArray1,
    /// OCCT myNbPassPoints.
    my_nb_pass_points: i32,
    /// OCCT myNbTangPoints.
    my_nb_tang_points: i32,
    /// OCCT myNbCurvPoints.
    my_nb_curv_points: i32,
    /// OCCT myTypConstraints.
    my_typ_constraints: IntArray1,
    /// OCCT myTtheta.
    my_ttheta: RealArray1,
    /// OCCT myTfthet.
    my_tfthet: RealArray1,
    /// OCCT myMaxDegree.
    my_max_degree: i32,
    /// OCCT myMaxSegment.
    my_max_segment: i32,
    /// OCCT myNbIterations.
    my_nb_iterations: i32,
    /// OCCT myTolerance.
    my_tolerance: f64,
    /// OCCT myContinuity.
    my_continuity: GeomAbsShape,
    /// OCCT myNivCont.
    my_niv_cont: i32,
    /// OCCT myWithMinMax.
    my_with_min_max: bool,
    /// OCCT myWithCutting.
    my_with_cutting: bool,
    /// OCCT myPercent[3].
    my_percent: [f64; 3],
    /// OCCT myCriterium[4].
    my_criterium: [f64; 4],
    /// OCCT mySmoothCriterion.
    my_smooth_criterion: CriterionHandle,
    /// OCCT myParameters.
    my_parameters: RealArray1,
    /// OCCT myKnots.
    my_knots: RealArray1,
    /// OCCT myMBSpCurve.
    my_mbsp_curve: MultiBSpCurve,
    /// OCCT myMaxError.
    my_max_error: f64,
    /// OCCT myMaxErrorIndex.
    my_max_error_index: i32,
    /// OCCT myAverageError.
    my_average_error: f64,
    /// OCCT myIsCreated.
    my_is_created: bool,
    /// OCCT myIsDone.
    my_is_done: bool,
    /// OCCT myIsOverConstr.
    my_is_over_constr: bool,
}

impl AppDefVariational {
    /// OCCT AppDef_Variational(SSP, FirstPoint, LastPoint, TheConstraints)
    /// with the declared default arguments (hxx L62-73): MaxDegree = 14,
    /// MaxSegment = 100, Continuity = GeomAbs_C2, WithMinMax = false,
    /// WithCutting = true, Tolerance = 1.0, NbIterations = 2.
    pub fn new(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: CoupleArray1,
    ) -> Self {
        Self::new_full(
            ssp,
            first_point,
            last_point,
            the_constraints,
            14,
            100,
            GeomAbsShape::C2,
            false,
            true,
            1.0,
            2,
        )
    }

    /// OCCT AppDef_Variational::AppDef_Variational (cxx L71-219).
    #[allow(clippy::too_many_arguments)]
    pub fn new_full(
        ssp: &MultiLine,
        first_point: i32,
        last_point: i32,
        the_constraints: CoupleArray1,
        max_degree: i32,
        max_segment: i32,
        continuity: GeomAbsShape,
        with_min_max: bool,
        with_cutting: bool,
        tolerance: f64,
        nb_iterations: i32,
    ) -> Self {
        // OCCT member init list (cxx L83-93).
        let my_ssp = ssp.clone();
        let my_first_point = first_point;
        let my_last_point = last_point;
        let my_constraints = the_constraints;
        let mut my_max_degree = max_degree;
        let my_max_segment = max_segment;
        let my_nb_iterations = nb_iterations;
        let my_tolerance = tolerance;
        let my_continuity = continuity;
        let my_with_min_max = with_min_max;
        let my_with_cutting = with_cutting;

        // Verifications:
        if my_max_degree < 1 {
            panic!("Standard_DomainError: AppDef_Variational");
        }
        my_max_degree = my_max_degree.min(30); // std::min(30, myMaxDegree)
        //
        if my_max_segment < 1 {
            panic!("Standard_DomainError: AppDef_Variational");
        }
        //
        // (myWithMinMax != 0 && myWithMinMax != 1 - a bool is always 0 or 1.)
        //
        let my_is_over_constr = false;
        let my_is_created = false;
        let my_is_done = false;
        let my_niv_cont;
        match my_continuity {
            GeomAbsShape::C0 => my_niv_cont = 0,
            GeomAbsShape::C1 => my_niv_cont = 1,
            GeomAbsShape::C2 => my_niv_cont = 2,
            _ => panic!("Standard_ConstructionError: AppDef_Variational"),
        }
        //
        let my_nb_p2d = my_line_tool::nb_p2d(ssp) as i32;
        let my_nb_p3d = my_line_tool::nb_p3d(ssp) as i32;
        let my_dimension = 2 * my_nb_p2d + 3 * my_nb_p3d;
        //
        let my_percent = [0.4, 0.2, 0.4];
        let mut my_knots = RealArray1::new(1, 2);
        my_knots.set_value(1, 0.0);
        my_knots.set_value(2, 1.0);

        //  Declaration
        //
        let my_smooth_criterion: CriterionHandle =
            Rc::new(RefCell::new(LinearCriteria::new(ssp, my_first_point, my_last_point)));
        let my_parameters = RealArray1::new(my_first_point, my_last_point);
        let my_nb_points = my_last_point - my_first_point + 1;
        if my_nb_points <= 0 {
            panic!("Standard_ConstructionError: AppDef_Variational");
        }
        //
        let mut my_tab_points = RealArray1::new(1, my_dimension * my_nb_points);
        //
        //  Table of Points initialization
        //
        let mut tab_p3d = vec![DVec3::ZERO; 1.max(my_nb_p3d) as usize];
        let mut tab_p2d = vec![DVec2::ZERO; 1.max(my_nb_p2d) as usize];
        let mut index = 1;

        for ipoint in my_first_point..=my_last_point {
            if my_nb_p2d != 0 && my_nb_p3d == 0 {
                my_line_tool::value_2d(ssp, ipoint as usize, &mut tab_p2d);

                for jp2d in 1..=my_nb_p2d {
                    let p2d = tab_p2d[(jp2d - 1) as usize];

                    my_tab_points.set_value(index, p2d.x);
                    index += 1;
                    my_tab_points.set_value(index, p2d.y);
                    index += 1;
                }
            }
            if my_nb_p3d != 0 && my_nb_p2d == 0 {
                my_line_tool::value_3d(ssp, ipoint as usize, &mut tab_p3d);

                for jp3d in 1..=my_nb_p3d {
                    let p3d = tab_p3d[(jp3d - 1) as usize];

                    my_tab_points.set_value(index, p3d.x);
                    index += 1;
                    my_tab_points.set_value(index, p3d.y);
                    index += 1;
                    my_tab_points.set_value(index, p3d.z);
                    index += 1;
                }
            }
            if my_nb_p3d != 0 && my_nb_p2d != 0 {
                my_line_tool::value_3d_2d(ssp, ipoint as usize, &mut tab_p3d, &mut tab_p2d);

                for jp3d in 1..=my_nb_p3d {
                    let p3d = tab_p3d[(jp3d - 1) as usize];

                    my_tab_points.set_value(index, p3d.x);
                    index += 1;
                    my_tab_points.set_value(index, p3d.y);
                    index += 1;
                    my_tab_points.set_value(index, p3d.z);
                    index += 1;
                }
                for jp2d in 1..=my_nb_p2d {
                    let p2d = tab_p2d[(jp2d - 1) as usize];

                    my_tab_points.set_value(index, p2d.x);
                    index += 1;
                    my_tab_points.set_value(index, p2d.y);
                    index += 1;
                }
            }
        }

        let mut this = AppDefVariational {
            my_ssp,
            my_nb_p3d,
            my_nb_p2d,
            my_dimension,
            my_first_point,
            my_last_point,
            my_nb_points,
            my_tab_points,
            my_constraints,
            my_nb_constraints: 0,
            my_tab_constraints: RealArray1::new(1, 1),
            my_nb_pass_points: 0,
            my_nb_tang_points: 0,
            my_nb_curv_points: 0,
            my_typ_constraints: IntArray1::new(1, 1),
            my_ttheta: RealArray1::new(1, 1),
            my_tfthet: RealArray1::new(1, 1),
            my_max_degree,
            my_max_segment,
            my_nb_iterations,
            my_tolerance,
            my_continuity,
            my_niv_cont,
            my_with_min_max,
            my_with_cutting,
            my_percent,
            my_criterium: [0.0; 4],
            my_smooth_criterion,
            my_parameters,
            my_knots,
            my_mbsp_curve: MultiBSpCurve::new_nbpol(0),
            my_max_error: 0.0,
            my_max_error_index: 0,
            my_average_error: 0.0,
            my_is_created,
            my_is_done,
            my_is_over_constr,
        };
        // OCCT L218: Init();
        this.init();
        this
    }

    /// OCCT AppDef_Variational::Init (cxx L224-484).
    fn init(&mut self) {
        let mut cur_multy_point;
        let mut tab_v3d = vec![DVec3::ZERO; 1.max(self.my_nb_p3d) as usize];
        let mut tab_v2d = vec![DVec2::ZERO; 1.max(self.my_nb_p2d) as usize];
        let mut tab_v3dcurv = vec![DVec3::ZERO; 1.max(self.my_nb_p3d) as usize];
        let mut tab_v2dcurv = vec![DVec2::ZERO; 1.max(self.my_nb_p2d) as usize];

        self.my_nb_constraints = self.my_constraints.length();
        if self.my_nb_constraints < 0 {
            panic!("Standard_ConstructionError: AppDef_Variational::Init");
        }

        self.my_typ_constraints = IntArray1::new(1, 1.max(2 * self.my_nb_constraints));
        self.my_tab_constraints =
            RealArray1::new(1, 1.max(2 * self.my_dimension * self.my_nb_constraints));
        self.my_ttheta = RealArray1::new(
            1,
            1.max((2 * self.my_nb_p2d + 6 * self.my_nb_p3d) * self.my_nb_constraints),
        );
        self.my_tfthet = RealArray1::new(
            1,
            1.max((2 * self.my_nb_p2d + 6 * self.my_nb_p3d) * self.my_nb_constraints),
        );

        //
        // Table of types initialization
        let mut index = 1;
        let mut jndex = 1;
        cur_multy_point = 1;
        self.my_nb_pass_points = 0;
        self.my_nb_tang_points = 0;
        self.my_nb_curv_points = 0;
        let mut valcontr;

        let (cons_lower, cons_upper) = (self.my_constraints.lower(), self.my_constraints.upper());
        for iconstr in cons_lower..=cons_upper {
            let couple = self.my_constraints.value(iconstr);
            let ipoint = couple.index;

            valcontr = couple.constraint;
            match valcontr {
                AppParConstraint::NoConstraint => {
                    cur_multy_point -= self.my_nb_p3d * 6 + self.my_nb_p2d * 2;
                }
                AppParConstraint::PassPoint => {
                    self.my_typ_constraints.set_value(index, ipoint);
                    index += 1;
                    self.my_typ_constraints.set_value(index, 0);
                    index += 1;
                    self.my_nb_pass_points += 1;
                    if self.my_nb_p2d != 0 {
                        jndex += 4 * self.my_nb_p2d;
                    }
                    if self.my_nb_p3d != 0 {
                        jndex += 6 * self.my_nb_p3d;
                    }
                }
                AppParConstraint::TangencyPoint => {
                    self.my_typ_constraints.set_value(index, ipoint);
                    index += 1;
                    self.my_typ_constraints.set_value(index, 1);
                    index += 1;
                    self.my_nb_tang_points += 1;
                    if self.my_nb_p2d != 0 && self.my_nb_p3d == 0 {
                        if !my_line_tool::tangency_2d(&self.my_ssp, ipoint as usize, &mut tab_v2d)
                        {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        for jp2d in 1..=self.my_nb_p2d {
                            let mut vt2d = tab_v2d[(jp2d - 1) as usize];
                            vt2d = vt2d.normalize();
                            self.my_tab_constraints.set_value(jndex, vt2d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt2d.y);
                            jndex += 1;
                            jndex += 2;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + (jp2d - 1) * 2,
                                jndex - 4,
                            );
                        }
                    }
                    if self.my_nb_p3d != 0 && self.my_nb_p2d == 0 {
                        if !my_line_tool::tangency_3d(&self.my_ssp, ipoint as usize, &mut tab_v3d)
                        {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        for jp3d in 1..=self.my_nb_p3d {
                            let mut vt3d = tab_v3d[(jp3d - 1) as usize];
                            vt3d = vt3d.normalize();
                            self.my_tab_constraints.set_value(jndex, vt3d.x);
                            jndex += 1;

                            self.my_tab_constraints.set_value(jndex, vt3d.y);
                            jndex += 1;

                            self.my_tab_constraints.set_value(jndex, vt3d.z);
                            jndex += 1;
                            jndex += 3;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + (jp3d - 1) * 6,
                                jndex - 6,
                            );
                        }
                    }
                    if self.my_nb_p3d != 0 && self.my_nb_p2d != 0 {
                        if !my_line_tool::tangency_3d_2d(
                            &self.my_ssp,
                            ipoint as usize,
                            &mut tab_v3d,
                            &mut tab_v2d,
                        ) {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        for jp3d in 1..=self.my_nb_p3d {
                            let mut vt3d = tab_v3d[(jp3d - 1) as usize];
                            vt3d = vt3d.normalize();
                            self.my_tab_constraints.set_value(jndex, vt3d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt3d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt3d.z);
                            jndex += 1;
                            jndex += 3;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + (jp3d - 1) * 6,
                                jndex - 6,
                            );
                        }

                        for jp2d in 1..=self.my_nb_p2d {
                            let mut vt2d = tab_v2d[(jp2d - 1) as usize];
                            vt2d = vt2d.normalize();
                            self.my_tab_constraints.set_value(jndex, vt2d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt2d.y);
                            jndex += 1;
                            jndex += 2;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + (jp2d - 1) * 2 + self.my_nb_p3d * 6,
                                jndex - 4,
                            );
                        }
                    }
                }
                AppParConstraint::CurvaturePoint => {
                    self.my_typ_constraints.set_value(index, ipoint);
                    index += 1;
                    self.my_typ_constraints.set_value(index, 2);
                    index += 1;
                    self.my_nb_curv_points += 1;
                    if self.my_nb_p2d != 0 && self.my_nb_p3d == 0 {
                        if !my_line_tool::tangency_2d(&self.my_ssp, ipoint as usize, &mut tab_v2d)
                        {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        if !my_line_tool::curvature_2d(
                            &self.my_ssp,
                            ipoint as usize,
                            &mut tab_v2dcurv,
                        ) {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        for jp2d in 1..=self.my_nb_p2d {
                            let mut vt2d = tab_v2d[(jp2d - 1) as usize];
                            vt2d = vt2d.normalize();
                            let vc2d = tab_v2dcurv[(jp2d - 1) as usize];
                            if (gp_vec2d_angle(vc2d, vt2d).abs() - std::f64::consts::FRAC_PI_2)
                                .abs()
                                > rcad_kernel::core::precision::ANGULAR
                            {
                                panic!("Standard_ConstructionError: AppDef_Variational::Init");
                            }
                            self.my_tab_constraints.set_value(jndex, vt2d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt2d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc2d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc2d.y);
                            jndex += 1;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + (jp2d - 1) * 2,
                                jndex - 4,
                            );
                        }
                    }

                    if self.my_nb_p3d != 0 && self.my_nb_p2d == 0 {
                        if !my_line_tool::tangency_3d(&self.my_ssp, ipoint as usize, &mut tab_v3d)
                        {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        if !my_line_tool::curvature_3d(
                            &self.my_ssp,
                            ipoint as usize,
                            &mut tab_v3dcurv,
                        ) {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        for jp3d in 1..=self.my_nb_p3d {
                            let mut vt3d = tab_v3d[(jp3d - 1) as usize];
                            vt3d = vt3d.normalize();
                            let vc3d = tab_v3dcurv[(jp3d - 1) as usize];
                            if !gp_vec_is_normal(
                                vc3d.normalize(),
                                vt3d,
                                rcad_kernel::core::precision::ANGULAR,
                            ) {
                                panic!("Standard_ConstructionError: AppDef_Variational::Init");
                            }
                            self.my_tab_constraints.set_value(jndex, vt3d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt3d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt3d.z);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc3d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc3d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc3d.z);
                            jndex += 1;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + (jp3d - 1) * 6,
                                jndex - 6,
                            );
                        }
                    }
                    if self.my_nb_p3d != 0 && self.my_nb_p2d != 0 {
                        if !my_line_tool::tangency_3d_2d(
                            &self.my_ssp,
                            ipoint as usize,
                            &mut tab_v3d,
                            &mut tab_v2d,
                        ) {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        if !my_line_tool::curvature_3d_2d(
                            &self.my_ssp,
                            ipoint as usize,
                            &mut tab_v3dcurv,
                            &mut tab_v2dcurv,
                        ) {
                            panic!("Standard_ConstructionError: AppDef_Variational::Init");
                        }
                        for jp3d in 1..=self.my_nb_p3d {
                            let mut vt3d = tab_v3d[(jp3d - 1) as usize];
                            vt3d = vt3d.normalize();
                            let vc3d = tab_v3dcurv[(jp3d - 1) as usize];
                            if !gp_vec_is_normal(
                                vc3d.normalize(),
                                vt3d,
                                rcad_kernel::core::precision::ANGULAR,
                            ) {
                                panic!("Standard_ConstructionError: AppDef_Variational::Init");
                            }
                            self.my_tab_constraints.set_value(jndex, vt3d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt3d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt3d.z);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc3d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc3d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc3d.z);
                            jndex += 1;
                            self.init_ttheta_f(
                                3,
                                valcontr,
                                cur_multy_point + (jp3d - 1) * 6,
                                jndex - 6,
                            );
                        }
                        for jp2d in 1..=self.my_nb_p2d {
                            let mut vt2d = tab_v2d[(jp2d - 1) as usize];
                            vt2d = vt2d.normalize();
                            let vc2d = tab_v2dcurv[(jp2d - 1) as usize];
                            if (gp_vec2d_angle(vc2d, vt2d).abs() - std::f64::consts::FRAC_PI_2)
                                .abs()
                                > rcad_kernel::core::precision::ANGULAR
                            {
                                panic!("Standard_ConstructionError: AppDef_Variational::Init");
                            }
                            self.my_tab_constraints.set_value(jndex, vt2d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vt2d.y);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc2d.x);
                            jndex += 1;
                            self.my_tab_constraints.set_value(jndex, vc2d.y);
                            jndex += 1;
                            self.init_ttheta_f(
                                2,
                                valcontr,
                                cur_multy_point + (jp2d - 1) * 2 + self.my_nb_p3d * 6,
                                jndex - 4,
                            );
                        }
                    }
                }
            }
            cur_multy_point += self.my_nb_p3d * 6 + self.my_nb_p2d * 2;
        }
        // OverConstraint Detection
        let max_seg;
        if self.my_with_cutting {
            max_seg = self.my_max_segment;
        } else {
            max_seg = 1;
        }
        if ((self.my_max_degree - self.my_niv_cont) * max_seg
            - self.my_nb_pass_points
            - 2 * self.my_nb_tang_points
            - 3 * self.my_nb_curv_points)
            < 0
        {
            self.my_is_over_constr = true;
            self.my_is_created = false;
        } else {
            self.init_smooth_criterion();
            self.my_is_created = true;
        }
    }

    /// OCCT AppDef_Variational::Approximate (cxx L489-660).
    pub fn approximate(&mut self) {
        if !self.my_is_created {
            panic!("StdFail_NotDone: AppDef_Variational::Approximate");
        }

        let mut w_quadratic = 0.0;
        let mut w_quality = 0.0;

        let mut ecarts = RealArray1::new(self.my_first_point, self.my_last_point);

        self.my_smooth_criterion
            .borrow()
            .get_weight(&mut w_quadratic, &mut w_quality);

        let mut the_curve: Option<CurveHandle> = None;

        self.my_smooth_criterion.borrow().get_curve(&mut the_curve);
        let mut the_curve = the_curve.expect("AppDef_Variational::Approximate: null TheCurve");

        //---------------------------------------------------------------------

        let mut j = self.my_smooth_criterion.clone();
        self.the_motor(&j, w_quadratic, w_quality, &mut the_curve, &mut ecarts);

        if self.my_with_min_max && self.my_tolerance < self.my_max_error {
            self.adjusting(&mut j, &mut w_quadratic, &mut w_quality, &mut the_curve, &mut ecarts);
            // OCCT Adjusting rebinds the caller's handle (J = JNew,
            // cxx L2962); the field write-back mirrors that reference
            // semantics (Rust has no aliased handle& parameter).
            self.my_smooth_criterion = j.clone();
        }

        //---------------------------------------------------------------------

        let nb_elem = the_curve.borrow().nb_elements();

        let mut tab_p3d = vec![DVec3::ZERO; 1.max(self.my_nb_p3d) as usize];
        let mut tab_p2d = vec![DVec2::ZERO; 1.max(self.my_nb_p2d) as usize];
        // double debfin[2] = {-1., 1};
        let debfin = [-1.0f64, 1.0f64];

        {
            // OCCT: PolynomialIntervalsPtr = new HArray2<double>(1, NbElem, 1, 2)
            // - mirrored by a row-major flat Vec (NbElem rows of 2 columns).
            let mut polynomial_intervals = vec![0.0f64; (nb_elem * 2) as usize];

            // OCCT: NbCoeffPtr = new HArray1<int>(1, myMaxSegment).
            let mut nb_coeff: Vec<i32> = vec![0; self.my_max_segment as usize];

            let size = self.my_max_segment * (self.my_max_degree + 1) * self.my_dimension;
            // OCCT: CoeffPtr = new HArray1<double>(1, size); CoeffPtr->Init(0.);
            let mut coeff = vec![0.0f64; size as usize];

            // OCCT: IntervallesPtr = new HArray1<double>(1, NbElem + 1);
            //       IntervallesPtr->ChangeArray1() = TheCurve->Knots();
            let intervalles: Vec<f64> = the_curve.borrow_mut().knots().to_vec();

            the_curve.borrow_mut().get_polynom(&mut coeff);

            for ii in 1..=nb_elem {
                nb_coeff[(ii - 1) as usize] = the_curve.borrow().degree(ii) + 1;
            }

            for ii in 0..nb_elem {
                // PolynomialIntervalsPtr->SetValue(ii, 1, debfin[0]); ... (ii, 2, debfin[1])
                polynomial_intervals[(ii * 2) as usize] = debfin[0];
                polynomial_intervals[(ii * 2 + 1) as usize] = debfin[1];
            }

            // OCCT: Convert_CompPolynomialToPoles AConverter(NbElem, myNivCont,
            //       myDimension, myMaxDegree, NbCoeffPtr, CoeffPtr,
            //       PolynomialIntervalsPtr, IntervallesPtr);
            // (the OCCT Continuity array is only read at entries 2..NbCurves;
            // it is filled with myNivCont.)
            let continuity = vec![self.my_niv_cont; (nb_elem + 1) as usize];
            let aconverter = ConvertCompPolynomialToPoles::from_arrays(
                nb_elem as usize,
                self.my_dimension as usize,
                self.my_max_degree as usize,
                &continuity,
                &nb_coeff,
                &coeff,
                &polynomial_intervals,
                &intervalles,
            );
            if aconverter.is_done() {
                let nb_poles = aconverter.nb_poles() as i32;
                // const NCollection_Array2<double>& aPoles = AConverter.Poles()
                // - the rcad mirror is a row-major flat (NbPoles x Dimension)
                // slice.
                let a_poles = aconverter.poles();
                let poles_dim = self.my_dimension as usize;
                let knots_data = aconverter.knots_vec();
                self.my_knots = real_array1_from(&knots_data);
                let mults_data = aconverter.multiplicities_vec();
                let mut mults = IntArray1::new(1, mults_data.len() as i32);
                for (k, m) in mults_data.iter().enumerate() {
                    mults.set_value(1 + k as i32, *m);
                }

                let mut tab_mu: Vec<MultiPoint> = Vec::with_capacity(nb_poles as usize);
                for ipole in 1..=nb_poles {
                    let mut index = 0usize; // aPoles.LowerCol() - 1 (row-major flat)
                    if self.my_nb_p3d != 0 {
                        for jp3d in 1..=self.my_nb_p3d {
                            let mut p3d = DVec3::ZERO;
                            p3d.x = a_poles[(ipole - 1) as usize * poles_dim + index];
                            index += 1;
                            p3d.y = a_poles[(ipole - 1) as usize * poles_dim + index];
                            index += 1;
                            p3d.z = a_poles[(ipole - 1) as usize * poles_dim + index];
                            index += 1;
                            tab_p3d[(jp3d - 1) as usize] = p3d;
                        }
                    }
                    if self.my_nb_p2d != 0 {
                        for jp2d in 1..=self.my_nb_p2d {
                            let mut p2d = DVec2::ZERO;
                            p2d.x = a_poles[(ipole - 1) as usize * poles_dim + index];
                            index += 1;
                            p2d.y = a_poles[(ipole - 1) as usize * poles_dim + index];
                            index += 1;
                            tab_p2d[(jp2d - 1) as usize] = p2d;
                        }
                    }
                    if self.my_nb_p2d != 0 && self.my_nb_p3d != 0 {
                        tab_mu.push(MultiPoint::new_tab_p3d_p2d(&tab_p3d, &tab_p2d));
                    } else if self.my_nb_p2d != 0 {
                        tab_mu.push(MultiPoint::new_tab_p2d(&tab_p2d));
                    } else {
                        tab_mu.push(MultiPoint::new_tab_p3d(&tab_p3d));
                    }
                }
                // OCCT: AppParCurves_MultiBSpCurve aCurve(TabMU, myKnots->Array1(),
                //       Mults->Array1()); - the ctor computes
                //       myDegree = ComputeDegree(Mults, NbPoles)
                //                = sum(Mults) - NbPoles - 1
                //       (AppParCurves_MultiBSpCurve.cxx L85-98).
                let mults_vec: Vec<usize> =
                    (mults.lower()..=mults.upper()).map(|i| mults.value(i) as usize).collect();
                let knots_vec: Vec<f64> =
                    (self.my_knots.lower()..=self.my_knots.upper())
                        .map(|i| self.my_knots.value(i))
                        .collect();
                let degree = mults_vec.iter().sum::<usize>() - nb_poles as usize - 1;
                let a_curve = MultiBSpCurve::from_multipoints(tab_mu, knots_vec, mults_vec, degree);
                self.my_mbsp_curve = a_curve;
                self.my_is_done = true;
            }
        }
    }

    /// OCCT AppDef_Variational::IsCreated (cxx L665-668).
    pub fn is_created(&self) -> bool {
        self.my_is_created
    }

    /// OCCT AppDef_Variational::IsDone (cxx L673-676).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT AppDef_Variational::IsOverConstrained (cxx L681-684).
    pub fn is_over_constrained(&self) -> bool {
        self.my_is_over_constr
    }

    /// OCCT AppDef_Variational::Value (cxx L687-694).
    pub fn value(&self) -> MultiBSpCurve {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::Value");
        }
        self.my_mbsp_curve.clone()
    }

    /// OCCT AppDef_Variational::MaxError (cxx L699-706).
    pub fn max_error(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::MaxError");
        }
        self.my_max_error
    }

    /// OCCT AppDef_Variational::MaxErrorIndex (cxx L711-718).
    pub fn max_error_index(&self) -> i32 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::MaxErrorIndex");
        }
        self.my_max_error_index
    }

    /// OCCT AppDef_Variational::QuadraticError (cxx L723-730).
    pub fn quadratic_error(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::QuadraticError");
        }
        self.my_criterium[0]
    }

    /// OCCT AppDef_Variational::Distance (cxx L735-787).
    pub fn distance(&self, mat: &mut Matrix) {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::Distance");
        }
        let mut tab_p3d = vec![DVec3::ZERO; 1.max(self.my_nb_p3d) as usize];
        let mut tab_p2d = vec![DVec2::ZERO; 1.max(self.my_nb_p2d) as usize];
        let j0 = mat.lower_col() - self.my_first_point;

        let mut pt3d;
        let mut pt2d;

        for ipoint in self.my_first_point..=self.my_last_point {
            let mut index = 1;
            if self.my_nb_p3d != 0 {
                my_line_tool::value_3d(&self.my_ssp, ipoint as usize, &mut tab_p3d);

                for jp3d in 1..=self.my_nb_p3d {
                    let p3d = tab_p3d[(jp3d - 1) as usize];
                    pt3d = self.mbsp_value_p3d(index, self.my_parameters.value(ipoint));
                    let v = p3d.distance(pt3d);
                    mat.set(index, j0 + ipoint, v);
                    index += 1;
                }
            }
            if self.my_nb_p2d != 0 {
                if self.my_nb_p3d == 0 {
                    my_line_tool::value_2d(&self.my_ssp, ipoint as usize, &mut tab_p2d);
                } else {
                    my_line_tool::value_3d_2d(
                        &self.my_ssp,
                        ipoint as usize,
                        &mut tab_p3d,
                        &mut tab_p2d,
                    );
                }
                for jp2d in 1..=self.my_nb_p2d {
                    let p2d = tab_p2d[(jp2d - 1) as usize];
                    pt2d = self.mbsp_value_p2d(index, self.my_parameters.value(ipoint));
                    let v = p2d.distance(pt2d);
                    mat.set(index, j0 + ipoint, v);
                    index += 1;
                }
            }
        }
    }

    /// OCCT AppParCurves_MultiBSpCurve::Value(CuIndex, U, Pt)
    /// (AppParCurves_MultiBSpCurve.cxx L123-143) - D0 evaluation of the
    /// CuIndex-th curve of myMBSpCurve at U through BSplCLib::D0
    /// (mirrored by `bspl_lib::eval_flat`, the flat-knots Eval).
    fn mbsp_value_p3d(&self, cu_index: i32, u: f64) -> DVec3 {
        if self.my_mbsp_curve.poles[0].dimension(cu_index as usize) != 3 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::Value");
        }
        let mut poles: Vec<DVec3> = Vec::new();
        self.my_mbsp_curve.curve(cu_index as usize, &mut poles);
        let mut flat_poles: Vec<f64> = Vec::with_capacity(poles.len() * 3);
        for p in &poles {
            flat_poles.extend_from_slice(&[p.x, p.y, p.z]);
        }
        let mults_i32: Vec<i32> = self.my_mbsp_curve.mults.iter().map(|&m| m as i32).collect();
        let degree = self.my_mbsp_curve.degree;
        let nflat = rcad_kernel::math::bspl_lib::knot_sequence_length(&mults_i32, degree, false);
        let mut flat_knots = vec![0.0f64; nflat];
        rcad_kernel::math::bspl_lib::knot_sequence(
            &self.my_mbsp_curve.knots,
            &mults_i32,
            degree,
            false,
            &mut flat_knots,
        );
        let mut results = vec![0.0f64; 3];
        let mut extrap_mode = [0i32; 2]; // BSplCLib::D0: no extrapolation
        rcad_kernel::math::bspl_lib::eval_flat(
            u,
            false,
            0,
            &mut extrap_mode,
            degree,
            &flat_knots,
            3,
            &flat_poles,
            &mut results,
        );
        DVec3::new(results[0], results[1], results[2])
    }

    /// OCCT AppParCurves_MultiBSpCurve::Value(CuIndex, U, Pt2d)
    /// (AppParCurves_MultiBSpCurve.cxx L147-167) - D0 evaluation of the
    /// CuIndex-th (2d) curve of myMBSpCurve at U through BSplCLib::D0
    /// (mirrored by `bspl_lib::eval_flat`, the flat-knots Eval).
    fn mbsp_value_p2d(&self, cu_index: i32, u: f64) -> DVec2 {
        if self.my_mbsp_curve.poles[0].dimension(cu_index as usize) != 2 {
            panic!("Standard_OutOfRange: AppParCurves_MultiBSpCurve::Value");
        }
        let mut poles: Vec<DVec2> = Vec::new();
        self.my_mbsp_curve.curve2d(cu_index as usize, &mut poles);
        let mut flat_poles: Vec<f64> = Vec::with_capacity(poles.len() * 2);
        for p in &poles {
            flat_poles.extend_from_slice(&[p.x, p.y]);
        }
        let mults_i32: Vec<i32> = self.my_mbsp_curve.mults.iter().map(|&m| m as i32).collect();
        let degree = self.my_mbsp_curve.degree;
        let nflat = rcad_kernel::math::bspl_lib::knot_sequence_length(&mults_i32, degree, false);
        let mut flat_knots = vec![0.0f64; nflat];
        rcad_kernel::math::bspl_lib::knot_sequence(
            &self.my_mbsp_curve.knots,
            &mults_i32,
            degree,
            false,
            &mut flat_knots,
        );
        let mut results = vec![0.0f64; 2];
        let mut extrap_mode = [0i32; 2]; // BSplCLib::D0: no extrapolation
        rcad_kernel::math::bspl_lib::eval_flat(
            u,
            false,
            0,
            &mut extrap_mode,
            degree,
            &flat_knots,
            2,
            &flat_poles,
            &mut results,
        );
        DVec2::new(results[0], results[1])
    }

    /// OCCT AppDef_Variational::AverageError (cxx L792-799).
    pub fn average_error(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::AverageError");
        }
        self.my_average_error
    }

    /// OCCT AppDef_Variational::Parameters (cxx L804-811).
    pub fn parameters(&self) -> &RealArray1 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::Parameters");
        }
        &self.my_parameters
    }

    /// OCCT AppDef_Variational::Knots (cxx L816-823).
    pub fn knots(&self) -> &RealArray1 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::Knots");
        }
        &self.my_knots
    }

    /// OCCT AppDef_Variational::Criterium (cxx L828-839).
    pub fn criterium(&self) -> (f64, f64, f64) {
        if !self.my_is_done {
            panic!("StdFail_NotDone: AppDef_Variational::Criterium");
        }
        (
            self.my_criterium[1],
            self.my_criterium[2],
            self.my_criterium[3],
        )
    }

    /// OCCT AppDef_Variational::CriteriumWeight (cxx L844-849).
    pub fn criterium_weight(&self) -> (f64, f64, f64) {
        (self.my_percent[0], self.my_percent[1], self.my_percent[2])
    }

    /// OCCT AppDef_Variational::MaxDegree (cxx L854-857).
    pub fn max_degree(&self) -> i32 {
        self.my_max_degree
    }

    /// OCCT AppDef_Variational::MaxSegment (cxx L862-865).
    pub fn max_segment(&self) -> i32 {
        self.my_max_segment
    }

    /// OCCT AppDef_Variational::Continuity (cxx L870-873).
    pub fn continuity(&self) -> GeomAbsShape {
        self.my_continuity
    }

    /// OCCT AppDef_Variational::WithMinMax (cxx L878-881).
    pub fn with_min_max(&self) -> bool {
        self.my_with_min_max
    }

    /// OCCT AppDef_Variational::WithCutting (cxx L886-889).
    pub fn with_cutting(&self) -> bool {
        self.my_with_cutting
    }

    /// OCCT AppDef_Variational::Tolerance (cxx L894-897).
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT AppDef_Variational::NbIterations (cxx L902-905).
    pub fn nb_iterations(&self) -> i32 {
        self.my_nb_iterations
    }

    /// OCCT AppDef_Variational::Dump (cxx L910-947) - prints on stdout the
    /// information on the current state of the object (the OCCT
    /// scientific/setprecision(3)/setw(9) formatting is mirrored with the
    /// Rust `{:>9.3e}` form).
    pub fn dump(&self) {
        println!(" \nVariational Smoothing ");
        println!(" Number of multipoints                   {}", self.my_nb_points);
        println!(" Number of 2d par multipoint {}", self.my_nb_p2d);
        println!(" Number of 3d per multipoint {}", self.my_nb_p3d);
        println!(" Number of PassagePoint      {}", self.my_nb_pass_points);
        println!(" Number of TangencyPoints    {}", self.my_nb_tang_points);
        println!(" Number of CurvaturePoints   {}", self.my_nb_curv_points);
        println!(" \nTolerance {:>9.3e}", self.my_tolerance);
        if self.with_min_max() {
            println!("  as Max Error.");
        } else {
            println!("  as size Error.");
        }
        println!(
            "CriteriumWeights : {} , {} , {}",
            self.my_percent[0], self.my_percent[1], self.my_percent[2]
        );

        if self.my_is_done {
            println!(" MaxError             {:>9.3e}", self.my_max_error);
            println!(" Index of  MaxError   {}", self.my_max_error_index);
            println!(" Average Error        {:>9.3e}", self.my_average_error);
            println!(" Quadratic Error      {:>9.3e}", self.my_criterium[0]);
            println!(" Tension Criterium    {:>9.3e}", self.my_criterium[1]);
            println!(" Flexion  Criterium   {:>9.3e}", self.my_criterium[2]);
            println!(" Jerk  Criterium      {:>9.3e}", self.my_criterium[3]);
            println!(" NbSegments           {}", self.my_knots.length() - 1);
        } else if self.my_is_over_constr {
            println!(" The problem is overconstraint");
        } else {
            println!(" Error in approximation");
        }
    }

    /// OCCT AppDef_Variational::SetConstraints (cxx L951-958).
    pub fn set_constraints(&mut self, a_constraint: CoupleArray1) -> bool {
        self.my_constraints = a_constraint;
        self.init();
        !self.my_is_over_constr
    }

    /// OCCT AppDef_Variational::SetParameters (cxx L963-966) - the OCCT
    /// `ChangeArray1() = Array1()` copies bounds and data -> the Rust clone.
    pub fn set_parameters(&mut self, param: &RealArray1) {
        self.my_parameters = param.clone();
    }

    /// OCCT AppDef_Variational::SetKnots (cxx L971-975).
    pub fn set_knots(&mut self, knots: &RealArray1) -> bool {
        self.my_knots = knots.clone();
        true
    }

    /// OCCT AppDef_Variational::SetMaxDegree (cxx L980-996).
    pub fn set_max_degree(&mut self, degree: i32) -> bool {
        if ((degree - self.my_niv_cont) * self.my_max_segment
            - self.my_nb_pass_points
            - 2 * self.my_nb_tang_points
            - 3 * self.my_nb_curv_points)
            < 0
        {
            false
        } else {
            self.my_max_degree = degree;

            self.init_smooth_criterion();

            true
        }
    }

    /// OCCT AppDef_Variational::SetMaxSegment (cxx L1001-1015).
    pub fn set_max_segment(&mut self, nb_segment: i32) -> bool {
        if self.my_with_cutting
            && ((self.my_max_degree - self.my_niv_cont) * nb_segment
                - self.my_nb_pass_points
                - 2 * self.my_nb_tang_points
                - 3 * self.my_nb_curv_points)
                < 0
        {
            false
        } else {
            self.my_max_segment = nb_segment;
            true
        }
    }

    /// OCCT AppDef_Variational::SetContinuity (cxx L1020-1051).
    pub fn set_continuity(&mut self, c: GeomAbsShape) -> bool {
        let niv_cont;
        match c {
            GeomAbsShape::C0 => niv_cont = 0,
            GeomAbsShape::C1 => niv_cont = 1,
            GeomAbsShape::C2 => niv_cont = 2,
            _ => panic!("Standard_ConstructionError: AppDef_Variational::SetContinuity"),
        }
        if ((self.my_max_degree - niv_cont) * self.my_max_segment
            - self.my_nb_pass_points
            - 2 * self.my_nb_tang_points
            - 3 * self.my_nb_curv_points)
            < 0
        {
            false
        } else {
            self.my_continuity = c;
            self.my_niv_cont = niv_cont;

            self.init_smooth_criterion();
            true
        }
    }

    /// OCCT AppDef_Variational::SetWithMinMax (cxx L1056-1061).
    pub fn set_with_min_max(&mut self, min_max: bool) {
        self.my_with_min_max = min_max;

        self.init_smooth_criterion();
    }

    /// OCCT AppDef_Variational::SetWithCutting (cxx L1066-1098).
    pub fn set_with_cutting(&mut self, cutting: bool) -> bool {
        if !cutting {
            if ((self.my_max_degree - self.my_niv_cont) * self.my_knots.length()
                - self.my_nb_pass_points
                - 2 * self.my_nb_tang_points
                - 3 * self.my_nb_curv_points)
                < 0
            {
                false
            } else {
                self.my_with_cutting = cutting;
                self.init_smooth_criterion();
                true
            }
        } else {
            if ((self.my_max_degree - self.my_niv_cont) * self.my_max_segment
                - self.my_nb_pass_points
                - 2 * self.my_nb_tang_points
                - 3 * self.my_nb_curv_points)
                < 0
            {
                false
            } else {
                self.my_with_cutting = cutting;
                self.init_smooth_criterion();
                true
            }
        }
    }

    /// OCCT AppDef_Variational::SetCriteriumWeight(Percent1, Percent2,
    /// Percent3) (cxx L1103-1117).
    pub fn set_criterium_weight(&mut self, percent1: f64, percent2: f64, percent3: f64) {
        if percent1 < 0.0 || percent2 < 0.0 || percent3 < 0.0 {
            panic!("Standard_DomainError: AppDef_Variational::SetCriteriumWeight");
        }
        let total = percent1 + percent2 + percent3;
        self.my_percent[0] = percent1 / total;
        self.my_percent[1] = percent2 / total;
        self.my_percent[2] = percent3 / total;

        self.init_smooth_criterion();
    }

    /// OCCT AppDef_Variational::SetCriteriumWeight(Order, Percent)
    /// (cxx L1122-1139).
    pub fn set_criterium_weight_order(&mut self, order: i32, percent: f64) {
        if percent < 0.0 {
            panic!("Standard_DomainError: AppDef_Variational::SetCriteriumWeight");
        }
        if !(1..=3).contains(&order) {
            panic!("Standard_ConstructionError: AppDef_Variational::SetCriteriumWeight");
        }
        self.my_percent[(order - 1) as usize] = percent;
        let total = self.my_percent[0] + self.my_percent[1] + self.my_percent[2];
        self.my_percent[0] /= total;
        self.my_percent[1] /= total;
        self.my_percent[2] /= total;

        self.init_smooth_criterion();
    }

    /// OCCT AppDef_Variational::SetTolerance (cxx L1144-1148).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.my_tolerance = tol;
        self.init_smooth_criterion();
    }

    /// OCCT AppDef_Variational::SetNbIterations (cxx L1153-1156).
    pub fn set_nb_iterations(&mut self, iter: i32) {
        self.my_nb_iterations = iter;
    }
}
