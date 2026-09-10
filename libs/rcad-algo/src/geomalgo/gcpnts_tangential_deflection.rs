//! OCCT GCPnts_TangentialDeflection + GCPnts_DistFunction
//! (TKGeomBase/GCPnts) — 1:1 port.
//!
//! - GCPnts_DistFunction.hxx/.cxx: the curve/line square-distance function
//!   used by the deflection estimator (both the math_Function and the
//!   math_MultipleVarFunction views).
//! - GCPnts_TangentialDeflection.cxx: the tangential (angular +
//!   curvature) discretisation of a curve over an adaptor.
//!
//! The C++ template parameter `TheCurve` (Adaptor3d_Curve /
//! Adaptor2d_Curve2d with the L36-67 2D shims) is the
//! [`GCPntsCurve`](super::gcpnts_curve::GCPntsCurve) trait; a caller picks
//! the dimension by constructing the adaptor flavor it holds
//! ([`GCPntsCurve2d`](super::gcpnts_curve::GCPntsCurve2d) for 2D curves,
//! the `&dyn Adaptor3dCurve` blanket impl for 3D curves) and calling
//! [`TangentialDeflection::new`].
//!
//! Indexing convention: the OCCT lists are 1-based; the loop variables of
//! the translated bodies keep the OCCT numeric values and index the rcad
//! Vecs through the `-1` offset, exactly as printed in the OCCT source.
//!
//! Home note: OCCT GCPnts lives in ModelingData/TKGeomBase and the rcad
//! kernel hosts the GCPnts package (rcad-kernel/src/base/gcpnts,
//! rcad-kernel/src/math/gcpnts); the kernel tree is outside this batch's
//! exclusive file domain, so the adaptor-based GCPnts translations land in
//! geomalgo pending the move (see the interface-change report).

use glam::DVec3;

use rcad_kernel::core::precision::{ANGULAR, CONFUSION};
use rcad_kernel::math::opt::{pso_minimize, BrentMinimum};
use rcad_kernel::math::root::FunctionValue;

use super::gcpnts_curve::GCPntsCurve;

/// OCCT GCPnts_TangentialDeflection.cxx L34.
const US_3: f64 = 0.3333333333333333333333333333;

/// OCCT gp::Resolution() — the smallest positive double.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

// ---------------------------------------------------------------------------
// GCPnts_DistFunction (GCPnts_DistFunction.hxx L26-72, .cxx L27-86)
// ---------------------------------------------------------------------------

/// OCCT GCPnts_DistFunction — the square distance between the point on the
/// curve C(u), U1 <= u <= U2 and the line through C(U1) and C(U2), fed to
/// the minimization algorithms as a one-variable function (negated).
pub struct DistFunction<'a> {
    my_curve: &'a dyn GCPntsCurve,
    my_lin_p: DVec3,
    my_lin_dir: DVec3,
    my_u1: f64,
    my_u2: f64,
}

impl<'a> DistFunction<'a> {
    /// OCCT GCPnts_DistFunction(theCurve, U1, U2) (GCPnts_DistFunction.cxx
    /// L27-45).
    pub fn new(the_curve: &'a dyn GCPntsCurve, u1: f64, u2: f64) -> Self {
        let p1 = the_curve.d0(u1);
        let mut p2 = the_curve.d0(u2);
        let dir = if (p1 - p2).length_squared() > GP_RESOLUTION {
            p2 - p1
        } else {
            // For #28812.
            p2 = the_curve.d0(u1 + 0.01 * (u2 - u1));
            p2 - p1
        };
        DistFunction {
            my_curve: the_curve,
            my_lin_p: p1,
            my_lin_dir: dir,
            my_u1: u1,
            my_u2: u2,
        }
    }

    /// OCCT GCPnts_DistFunction::Value (GCPnts_DistFunction.cxx L49-57).
    pub fn value(&mut self, x: f64) -> Option<f64> {
        if x < self.my_u1 || x > self.my_u2 {
            return None;
        }
        // F = -myLin.SquareDistance(myCurve.Value(X)).
        let p = self.my_curve.d0(x);
        let d = p - self.my_lin_p;
        let t = d.dot(self.my_lin_dir) / self.my_lin_dir.length_squared();
        let q = d - self.my_lin_dir * t;
        Some(-q.length_squared())
    }
}

/// OCCT GCPnts_DistFunctionMV (GCPnts_DistFunction.hxx L51-71) — the
/// multi-variable view used by math_PSO (the BrentMinimum step goes through
/// the same math_Function interface).
struct DistFunctionMV<'a> {
    my_max_curv_lin_dist: &'a std::cell::RefCell<DistFunction<'a>>,
}

impl FunctionValue for DistFunctionMV<'_> {
    /// OCCT GCPnts_DistFunctionMV::Value (GCPnts_DistFunction.cxx L70-76) —
    /// forwards to the one-variable function at X(1).
    fn value(&mut self, x: f64) -> Option<f64> {
        self.my_max_curv_lin_dist.borrow_mut().value(x)
    }
}

// ---------------------------------------------------------------------------
// free helpers
// ---------------------------------------------------------------------------

/// OCCT EstimAngl (GCPnts_TangentialDeflection.cxx L69-81).
fn estim_angl(p1: DVec3, pm: DVec3, p2: DVec3) -> f64 {
    let v1 = pm - p1;
    let v2 = p2 - pm;
    let l = v1.length() * v2.length();
    if l > GP_RESOLUTION {
        v1.cross(v2).length() / l
    } else {
        0.0
    }
}

/// OCCT getIntervalIdx (GCPnts_TangentialDeflection.cxx L85-99) — the
/// continuity interval holding `the_param` (search from `the_previous_idx`).
/// The OCCT indices are 1-based; the returned index follows the same
/// convention (0 stands for the first interval).
fn get_interval_idx(the_param: f64, the_intervs: &[f64], the_previous_idx: usize) -> usize {
    // for (anIdx = thePreviousIdx; anIdx < theIntervs.Upper(); anIdx++)
    //   if (theParam >= theIntervs(anIdx) && theParam <= theIntervs(anIdx+1)) break;
    let mut an_idx = the_previous_idx;
    while an_idx + 1 < the_intervs.len() {
        if the_param >= the_intervs[an_idx] && the_param <= the_intervs[an_idx + 1] {
            break;
        }
        an_idx += 1;
    }
    an_idx
}

// ---------------------------------------------------------------------------
// GCPnts_TangentialDeflection
// ---------------------------------------------------------------------------

/// OCCT GCPnts_TangentialDeflection — discretisation of a curve with points
/// restricted by the maximal curvature deflection and by the angle between
/// tangents (GCPnts_TangentialDeflection.hxx L32-118).
pub struct TangentialDeflection {
    /// hxx L120: Standard_Real myAngularDeflection.
    my_angular_deflection: f64,
    /// hxx L121: Standard_Real myCurvatureDeflection.
    my_curvature_deflection: f64,
    /// hxx L122: Standard_Real myUTol.
    my_utol: f64,
    /// hxx L123: Standard_Integer myMinNbPnts.
    my_min_nb_pnts: i32,
    /// hxx L124: Standard_Real myMinLen.
    my_min_len: f64,
    /// hxx L126: NCollection_List<Standard_Real> myParameters.
    my_parameters: Vec<f64>,
    /// hxx L129: NCollection_List<gp_Pnt> myPoints.
    my_points: Vec<DVec3>,
    /// hxx L131: Standard_Real myLastU.
    my_last_u: f64,
    /// hxx L133: Standard_Real myFirstu.
    my_firstu: f64,
}

impl TangentialDeflection {
    /// OCCT GCPnts_TangentialDeflection() (cxx L104-113) — the empty
    /// constructor.
    pub fn empty() -> Self {
        TangentialDeflection {
            my_angular_deflection: 0.0,
            my_curvature_deflection: 0.0,
            my_utol: 0.0,
            my_min_nb_pnts: 0,
            my_min_len: 0.0,
            my_parameters: Vec::new(),
            my_points: Vec::new(),
            my_last_u: 0.0,
            my_firstu: 0.0,
        }
    }

    /// OCCT GCPnts_TangentialDeflection(theC, theAngularDeflection,
    /// theCurvatureDeflection, theMinimumOfPoints, theUTol, theMinLen)
    /// (cxx L117-137 / L169-189) — the full-domain constructors (the 2d and
    /// 3d flavors share the body; the dimension is carried by the adaptor).
    pub fn new(
        the_c: &dyn GCPntsCurve,
        the_angular_deflection: f64,
        the_curvature_deflection: f64,
        the_minimum_of_points: i32,
        the_utol: f64,
        the_min_len: f64,
    ) -> Self {
        let mut r = TangentialDeflection::empty();
        r.initialize(
            the_c,
            the_angular_deflection,
            the_curvature_deflection,
            the_minimum_of_points,
            the_utol,
            the_min_len,
        );
        r
    }

    /// OCCT GCPnts_TangentialDeflection(theC, theFirstParameter,
    /// theLastParameter, ...) (cxx L141-165 / L193-217) — the
    /// restricted-range constructors.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_range(
        the_c: &dyn GCPntsCurve,
        the_first_parameter: f64,
        the_last_parameter: f64,
        the_angular_deflection: f64,
        the_curvature_deflection: f64,
        the_minimum_of_points: i32,
        the_utol: f64,
        the_min_len: f64,
    ) -> Self {
        let mut r = TangentialDeflection::empty();
        r.initialize_with_range(
            the_c,
            the_first_parameter,
            the_last_parameter,
            the_angular_deflection,
            the_curvature_deflection,
            the_minimum_of_points,
            the_utol,
            the_min_len,
        );
        r
    }

    /// OCCT Initialize(theC, ...) (cxx L221-236 / L240-255) — the
    /// full-domain delegates.
    pub fn initialize(
        &mut self,
        the_c: &dyn GCPntsCurve,
        the_angular_deflection: f64,
        the_curvature_deflection: f64,
        the_minimum_of_points: i32,
        the_utol: f64,
        the_min_len: f64,
    ) {
        self.initialize_with_range(
            the_c,
            the_c.first_parameter(),
            the_c.last_parameter(),
            the_angular_deflection,
            the_curvature_deflection,
            the_minimum_of_points,
            the_utol,
            the_min_len,
        );
    }

    /// OCCT Initialize(theC, theFirstParameter, theLastParameter, ...)
    /// (cxx L259-297) — forwards to the private template `initialize`.
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_with_range(
        &mut self,
        the_c: &dyn GCPntsCurve,
        the_first_parameter: f64,
        the_last_parameter: f64,
        the_angular_deflection: f64,
        the_curvature_deflection: f64,
        the_minimum_of_points: i32,
        the_utol: f64,
        the_min_len: f64,
    ) {
        self.initialize_template(
            the_c,
            the_first_parameter,
            the_last_parameter,
            the_angular_deflection,
            the_curvature_deflection,
            the_minimum_of_points,
            the_utol,
            the_min_len,
        );
    }

    /// OCCT EvaluateDu (cxx L301-322).
    fn evaluate_du(
        &self,
        the_c: &dyn GCPntsCurve,
        the_u: f64,
        the_p: &mut DVec3,
        the_du: &mut f64,
        the_not_done: &mut bool,
    ) {
        let (p, t, n) = the_c.d2(the_u);
        *the_p = p;
        let lt = t.length();
        const L_TOL: f64 = CONFUSION;
        if lt > L_TOL && n.length() > L_TOL {
            let lc = n.cross(t).length();
            let ln = lc / lt;
            if ln > L_TOL {
                *the_du = (8.0 * self.my_curvature_deflection.max(self.my_min_len) / ln).sqrt();
                *the_not_done = false;
            }
        }
    }

    /// OCCT PerformLinear (cxx L326-348).
    fn perform_linear(&mut self, the_c: &dyn GCPntsCurve) {
        let p = the_c.d0(self.my_firstu);
        self.my_parameters.push(self.my_firstu);
        self.my_points.push(p);
        if self.my_min_nb_pnts > 2 {
            let du = (self.my_last_u - self.my_firstu) / self.my_min_nb_pnts as f64;
            let mut u = self.my_firstu + du;
            for _i in 2..self.my_min_nb_pnts {
                let p = the_c.d0(u);
                self.my_parameters.push(u);
                self.my_points.push(p);
                u += du;
            }
        }
        let p = the_c.d0(self.my_last_u);
        self.my_parameters.push(self.my_last_u);
        self.my_points.push(p);
    }

    /// OCCT PerformCircular (cxx L352-381).
    fn perform_circular(&mut self, the_c: &dyn GCPntsCurve) {
        // akm 8/01/02 : check the radius before divide by it.
        let df_r = the_c.circle_radius();
        let mut du = TangentialDeflection::arc_angular_step(
            df_r,
            self.my_curvature_deflection,
            self.my_angular_deflection,
            self.my_min_len,
        );

        let a_diff = self.my_last_u - self.my_firstu;
        // Round up number of points to satisfy curvatureDeflection more precisely.
        let mut nb_points = (a_diff / du).ceil().min(1.0e+6) as i32;
        nb_points = nb_points.max(self.my_min_nb_pnts - 1);
        du = a_diff / nb_points as f64;

        let mut p;
        let mut u = self.my_firstu;
        for _i in 1..=nb_points {
            p = the_c.d0(u);
            self.my_parameters.push(u);
            self.my_points.push(p);
            u += du;
        }

        p = the_c.d0(self.my_last_u);
        self.my_parameters.push(self.my_last_u);
        self.my_points.push(p);
    }

    /// OCCT initialize (the private template, cxx L385-454).
    #[allow(clippy::too_many_arguments)]
    fn initialize_template(
        &mut self,
        the_c: &dyn GCPntsCurve,
        the_first_parameter: f64,
        the_last_parameter: f64,
        the_angular_deflection: f64,
        the_curvature_deflection: f64,
        the_minimum_of_points: i32,
        the_utol: f64,
        the_min_len: f64,
    ) {
        use rcad_kernel::base::proj_lib::CurveType;
        if the_curvature_deflection < CONFUSION || the_angular_deflection < ANGULAR {
            panic!("GCPnts_TangentialDeflection::Initialize - Zero Deflection");
        }
        self.my_parameters.clear();
        self.my_points.clear();
        if the_first_parameter < the_last_parameter {
            self.my_firstu = the_first_parameter;
            self.my_last_u = the_last_parameter;
        } else {
            self.my_last_u = the_first_parameter;
            self.my_firstu = the_last_parameter;
        }
        self.my_utol = the_utol;
        self.my_angular_deflection = the_angular_deflection;
        self.my_curvature_deflection = the_curvature_deflection;
        self.my_min_nb_pnts = the_minimum_of_points.max(2);
        self.my_min_len = the_min_len.max(CONFUSION);

        match the_c.get_type() {
            CurveType::Line => {
                self.perform_linear(the_c);
            }
            CurveType::Circle => {
                self.perform_circular(the_c);
            }
            CurveType::BSpline => {
                // Handle(BSplineCurve) aBS = theC.BSpline();
                // if (aBS->NbPoles() == 2).
                if the_c.nb_poles() == 2 {
                    self.perform_linear(the_c);
                } else {
                    self.perform_curve(the_c);
                }
            }
            CurveType::Bezier => {
                // Handle(BezierCurve) aBZ = theC.Bezier();
                // if (aBZ->NbPoles() == 2).
                if the_c.nb_poles() == 2 {
                    self.perform_linear(the_c);
                } else {
                    self.perform_curve(the_c);
                }
            }
            _ => {
                self.perform_curve(the_c);
            }
        }
    }

    /// OCCT AddPoint (cxx L458-491) — inserts (or replaces) a point keeping
    /// the parameter order; returns its 1-based index.
    pub fn add_point(&mut self, the_pnt: DVec3, the_param: f64, the_is_replace: bool) -> i32 {
        const TOL: f64 = 1.0e-12; // Precision::PConfusion()
        let mut index: i32 = -1;
        let nb = self.my_parameters.len();
        for i in 0..nb {
            if index != -1 {
                break;
            }
            let dist = self.my_parameters[i] - the_param;
            if dist.abs() <= TOL {
                index = i as i32 + 1;
                if the_is_replace {
                    self.my_points[i] = the_pnt;
                    self.my_parameters[i] = the_param;
                }
            } else if dist > TOL {
                self.my_points.insert(i, the_pnt);
                self.my_parameters.insert(i, the_param);
                index = i as i32 + 1;
            }
        }
        if index == -1 {
            self.my_points.push(the_pnt);
            self.my_parameters.push(the_param);
            index = self.my_parameters.len() as i32;
        }
        index
    }

    /// OCCT ArcAngularStep (static, cxx L495-518).
    pub fn arc_angular_step(
        the_radius: f64,
        the_linear_deflection: f64,
        the_angular_deflection: f64,
        the_min_length: f64,
    ) -> f64 {
        if the_radius < 0.0 {
            panic!("Negative radius");
        }

        const A_PRECISION: f64 = CONFUSION;

        let mut du = 0.0;
        let mut a_min_size_ang = 0.0;
        if the_radius > A_PRECISION {
            du = (1.0 - (the_linear_deflection / the_radius)).max(0.0);

            // It is not suitable to consider min size greater than 1/4 arc len.
            if the_min_length > A_PRECISION {
                a_min_size_ang = (the_min_length / the_radius).min(std::f64::consts::FRAC_PI_2);
            }
        }
        du = 2.0 * du.acos();
        du = du.min(the_angular_deflection).max(a_min_size_ang);
        du
    }

    /// OCCT PerformCurve (cxx L522-956).
    fn perform_curve(&mut self, the_c: &dyn GCPntsCurve) {
        let mut v1: DVec3;
        let mut v2: DVec3;
        // The C++ locals are assigned before every read; the Rust
        // flow checker needs an initializer.
        let mut middle_point = DVec3::ZERO;
        let mut current_point = DVec3::ZERO;
        let last_point;

        let mut u1 = self.my_firstu;
        // LTol — protection against zero length.
        const L_TOL: f64 = CONFUSION;
        let mut atol = 1.0e-2 * self.my_angular_deflection;
        if atol > 1.0e-2 {
            atol = 1.0e-2;
        } else if atol < 1.0e-7 {
            atol = 1.0e-7;
        }

        last_point = the_c.d0(self.my_last_u);

        // Initialization of the computation.
        let mut not_done = true;
        let mut dusave = (self.my_last_u - self.my_firstu) * US_3;
        let mut du = dusave;
        current_point = DVec3::ZERO;
        self.evaluate_du(the_c, u1, &mut current_point, &mut du, &mut not_done);
        self.my_parameters.push(u1);
        self.my_points.push(current_point);

        // Used to detect "isLine" current bspline and in Du computation in
        // general handling.
        let nb_interv = the_c.nb_intervals_cn();
        let mut intervs = the_c.intervals_cn();

        if not_done || du > 5.0 * dusave {
            // It is either a line or a singularity:
            v1 = last_point - current_point;
            let mut l1 = v1.length();
            if l1 > L_TOL {
                // If it is a line, we verify by computing minNbPoints:
                let mut is_line = true;
                let mut nb_points = if self.my_min_nb_pnts > 3 {
                    self.my_min_nb_pnts
                } else {
                    3
                };
                use rcad_kernel::base::proj_lib::CurveType;
                match the_c.get_type() {
                    CurveType::BSpline => {
                        // NbPoints = std::max(BS->Degree() + 1, NbPoints).
                        nb_points = nb_points.max(the_c.curve_degree() + 1);
                    }
                    CurveType::Bezier => {
                        nb_points = nb_points.max(the_c.curve_degree() + 1);
                    }
                    _ => {}
                }
                let mut param = 0.0f64;
                let mut i = 1usize; // OCCT 1-based interval index
                while i <= nb_interv && is_line {
                    // Avoid usage intervals out of [myFirstu, myLastU].
                    if intervs[i] < self.my_firstu || intervs[i - 1] > self.my_last_u {
                        i += 1;
                        continue;
                    }

                    // Fix border points in applicable intervals, to avoid
                    // being out of the target interval.
                    if intervs[i - 1] < self.my_firstu && intervs[i] > self.my_firstu {
                        intervs[i - 1] = self.my_firstu;
                    }
                    if intervs[i - 1] < self.my_last_u && intervs[i] > self.my_last_u {
                        intervs[i] = self.my_last_u;
                    }

                    let delta = (intervs[i] - intervs[i - 1]) / nb_points as f64;
                    let mut j = 1i32;
                    while j <= nb_points && is_line {
                        param = intervs[i - 1] + j as f64 * delta;
                        middle_point = the_c.d0(param);
                        v2 = middle_point - current_point;
                        let l2 = v2.length();
                        if l2 > L_TOL {
                            let a_angle = v2.cross(v1).length() / (l1 * l2);
                            is_line = a_angle < atol;
                        }
                        j += 1;
                    }
                    i += 1;
                }

                if is_line {
                    self.my_parameters.clear();
                    self.my_points.clear();

                    self.perform_linear(the_c);
                    return;
                } else {
                    // It was a singularity, continue:
                    // Du = Dusave;
                    let mut du_loc = du;
                    let mut not_done_loc = true;
                    self.evaluate_du(the_c, param, &mut middle_point, &mut du_loc, &mut not_done_loc);
                    du = du_loc;
                }
            } else {
                du = (self.my_last_u - self.my_firstu) / 2.1;
                let middle_u = self.my_firstu + du;
                middle_point = the_c.d0(middle_u);
                v1 = middle_point - current_point;
                l1 = v1.length();
                if l1 < L_TOL {
                    // L1 < LTol: this is a zero-length curve, computation is
                    // finished: return a segment of 2 points (protection).
                    self.my_parameters.push(self.my_last_u);
                    self.my_points.push(last_point);
                    return;
                }
            }
        }

        if du > dusave {
            du = dusave;
        } else {
            dusave = du;
        }

        if du < self.my_utol {
            du = self.my_last_u - self.my_firstu;
            if du < self.my_utol {
                self.my_parameters.push(self.my_last_u);
                self.my_points.push(last_point);
                return;
            }
        }

        // Normal processing for a curve.
        let mut more_points = true;
        let mut u2 = self.my_firstu;
        let angle_max = self.my_angular_deflection * 0.5; // because we take the midpoint
        // Indexes of intervals of U1 and U2, used to handle non-uniform case.
        let mut a_idx: [usize; 2] = [0, 0];
        let mut is_need_to_check = false;
        let mut a_prev_point = *self.my_points.last().unwrap();

        while more_points {
            let mut coef = 0.0f64;
            let mut acoef = 0.0f64;
            let mut fcoef = 0.0f64;

            a_idx[0] = get_interval_idx(u1, &intervs, a_idx[0]);
            u2 += du;

            if u2 >= self.my_last_u {
                // End of curve
                u2 = self.my_last_u;
                current_point = last_point;
                du = u2 - u1;
                dusave = du;
            } else {
                current_point = the_c.d0(u2); // Next point
            }

            let mut too_large = false;
            let mut correction = true;
            let mut too_small = false;

            while correction {
                // Adjustment of Du
                if is_need_to_check {
                    a_idx[1] = get_interval_idx(u2, &intervs, a_idx[0]);
                    if a_idx[1] > a_idx[0] {
                        // Jump to another polynom.
                        // Set Du to the smallest value and check deflection
                        // on it.
                        if du > (intervs[a_idx[0] + 1] - intervs[a_idx[0]]) * US_3 {
                            du = (intervs[a_idx[0] + 1] - intervs[a_idx[0]]) * US_3;
                            u2 = u1 + du;
                            if u2 > self.my_last_u {
                                u2 = self.my_last_u;
                            }
                            current_point = the_c.d0(u2);
                        }
                    }
                }
                let middle_u = (u1 + u2) * 0.5; // Verification at the midpoint
                middle_point = the_c.d0(middle_u);

                v1 = current_point - a_prev_point; // Deflection criterion
                v2 = middle_point - a_prev_point;
                let l1 = v1.length();

                fcoef = if l1 > self.my_min_len {
                    v1.cross(v2).length() / (l1 * self.my_curvature_deflection)
                } else {
                    // FCoef = 0.0
                    0.0
                };

                v1 = current_point - middle_point; // Angular criterion
                let l1 = v1.length();
                let l2 = v2.length();
                if l1 > self.my_min_len && l2 > self.my_min_len {
                    let angg = v1.cross(v2).length() / (l1 * l2);
                    acoef = angg / angle_max;
                } else {
                    acoef = 0.0;
                }

                // Keep the most penalizing coefficient
                coef = acoef.max(fcoef);

                if is_need_to_check && coef < 0.55 {
                    is_need_to_check = false;
                    du = dusave;
                    u2 = u1 + du;
                    if u2 > self.my_last_u {
                        u2 = self.my_last_u;
                    }
                    current_point = the_c.d0(u2);
                    continue;
                }

                if coef <= 1.0 {
                    if (self.my_last_u - u2).abs() < self.my_utol {
                        self.my_parameters.push(self.my_last_u);
                        self.my_points.push(last_point);
                        more_points = false;
                        correction = false;
                    } else {
                        if coef >= 0.55 || too_large {
                            self.my_parameters.push(u2);
                            self.my_points.push(current_point);
                            a_prev_point = current_point;
                            correction = false;
                            is_need_to_check = true;
                        } else if too_small {
                            correction = false;
                            a_prev_point = current_point;
                        } else {
                            too_small = true;
                            // double UUU2 = U2;
                            du += ((u2 - u1) * (1.0 - coef)).min(du * US_3);

                            u2 = u1 + du;
                            if u2 > self.my_last_u {
                                u2 = self.my_last_u;
                            }
                            current_point = the_c.d0(u2);
                        }
                    }
                } else {
                    if coef >= 1.5 {
                        if a_prev_point.distance(*self.my_points.last().unwrap()) > CONFUSION {
                            self.my_parameters.push(u1);
                            self.my_points.push(a_prev_point);
                        }
                        u2 = middle_u;
                        du = u2 - u1;
                        current_point = middle_point;
                    } else {
                        du *= 0.9;
                        u2 = u1 + du;
                        current_point = the_c.d0(u2);
                        too_large = true;
                    }
                }
            }

            du = u2 - u1;

            if more_points {
                if u1 > self.my_firstu {
                    if fcoef > acoef {
                        // Deflection is the splitting criterion
                        let mut not_done_loc = true;
                        self.evaluate_du(the_c, u2, &mut current_point, &mut du, &mut not_done_loc);
                        if not_done_loc {
                            du += (du - dusave) * (du / dusave);
                            if du > 1.5 * dusave {
                                du = 1.5 * dusave;
                            }
                            if du < 0.75 * dusave {
                                du = 0.75 * dusave;
                            }
                        }
                    } else {
                        // Angle is the splitting criterion
                        du += (du - dusave) * (du / dusave);
                        if du > 1.5 * dusave {
                            du = 1.5 * dusave;
                        }
                        if du < 0.75 * dusave {
                            du = 0.75 * dusave;
                        }
                    }
                }

                if du < self.my_utol {
                    du = self.my_last_u - u2;
                    if du < self.my_utol {
                        self.my_parameters.push(self.my_last_u);
                        self.my_points.push(last_point);
                        more_points = false;
                    } else if du * US_3 > self.my_utol {
                        du *= US_3;
                    }
                }
                u1 = u2;
                dusave = du;
            }
        }
        // Readjustment of the second to last point:
        let i = self.my_points.len() as i32 - 1;
        if i >= 2 {
            // MiddleU = myParameters(i - 1)  (1-based).
            let mut middle_u = self.my_parameters[(i - 1 - 1) as usize];
            middle_u = (self.my_last_u + middle_u) * 0.5;
            middle_point = the_c.d0(middle_u);
            self.my_parameters[(i - 1) as usize] = middle_u;
            self.my_points[(i - 1) as usize] = middle_point;
        }

        // Add points at midpoints of segments if the minimum
        // number of points is not yet reached.
        let mut nbp = self.my_points.len() as i32;

        while nbp < self.my_min_nb_pnts {
            let mut i = 2usize; // OCCT 1-based
            while i <= nbp as usize {
                let middle_u =
                    (self.my_parameters[i - 1 - 1] + self.my_parameters[i - 1]) * 0.5;
                middle_point = the_c.d0(middle_u);
                self.my_parameters.insert(i - 1, middle_u);
                self.my_points.insert(i - 1, middle_point);
                nbp += 1;
                i += 2;
            }
        }
        // Additional check for intervals
        let min_len2 = self.my_min_len * self.my_min_len;
        let max_nbp = 10 * nbp;
        let mut i = 1usize; // OCCT 1-based
        while i < nbp as usize {
            let u1 = self.my_parameters[i - 1];
            let u2 = self.my_parameters[i];

            if u2 - u1 <= self.my_utol {
                i += 1;
                continue;
            }

            // Check maximal deflection on interval;
            let mut dmax = 0.0f64;
            let mut umax = 0.0f64;
            self.estim_defl_entry(the_c, u1, u2, &mut dmax, &mut umax);
            let p1 = self.my_points[i - 1];
            let p2 = self.my_points[i];
            middle_point = the_c.d0(umax);
            let amax = estim_angl(p1, middle_point, p2);
            if dmax > self.my_curvature_deflection || amax > angle_max {
                if umax - u1 > self.my_utol && u2 - umax > self.my_utol {
                    if p1.distance_squared(middle_point) > min_len2
                        && p2.distance_squared(middle_point) > min_len2
                    {
                        // myParameters.InsertAfter(i, umax) — 1-based i maps
                        // to the 0-based insert position i.
                        self.my_parameters.insert(i, umax);
                        self.my_points.insert(i, middle_point);
                        nbp += 1;
                        // --i: to compensate ++i in the loop header, i must
                        // point to the first part of the split interval.
                        i -= 1;
                        if nbp > max_nbp {
                            break;
                        }
                    }
                }
            }
            i += 1;
        }
    }

    /// OCCT EstimDefl (cxx L960-1010) — the maximal deviation of the curve
    /// from the [theU1, theU2] chord.
    fn estim_defl_entry(
        &self,
        the_c: &dyn GCPntsCurve,
        the_u1: f64,
        the_u2: f64,
        the_max_defl: &mut f64,
        the_umax: &mut f64,
    ) {
        estim_defl(
            &self.my_utol,
            &self.my_last_u,
            &self.my_firstu,
            the_c,
            the_u1,
            the_u2,
            the_max_defl,
            the_umax,
        );
    }

    /// OCCT Parameters() (hxx accessor).
    pub fn parameters(&self) -> &Vec<f64> {
        &self.my_parameters
    }

    /// OCCT Points() (hxx accessor).
    pub fn points(&self) -> &Vec<DVec3> {
        &self.my_points
    }
}

/// OCCT EstimDefl (GCPnts_TangentialDeflection.cxx L960-1010) — the free
/// form carrying the member context explicitly.
#[allow(clippy::too_many_arguments)]
fn estim_defl(
    my_utol: &f64,
    my_last_u: &f64,
    my_firstu: &f64,
    the_c: &dyn GCPntsCurve,
    the_u1: f64,
    the_u2: f64,
    the_max_defl: &mut f64,
    the_umax: &mut f64,
) {
    let du = my_last_u - my_firstu;
    //
    // typename GCPnts_TCurveTypes<TheCurve>::DistFunction aFunc(theC, theU1, theU2);
    let a_func = std::cell::RefCell::new(DistFunction::new(the_c, the_u1, the_u2));
    //
    let a_nb_iter = 100;
    let a_rel_tol = (1.0e-3f64).max(2.0 * my_utol / (the_u1.abs() + the_u2.abs()));
    //
    let mut an_opt_loc = BrentMinimum::new(a_rel_tol, a_nb_iter, *my_utol);
    an_opt_loc.perform(
        &mut DistFunctionMV {
            my_max_curv_lin_dist: &a_func,
        },
        the_u1,
        (the_u1 + the_u2) / 2.0,
        the_u2,
    );
    if an_opt_loc.is_done() {
        *the_max_defl = (-an_opt_loc.minimum()).sqrt();
        *the_umax = an_opt_loc.location();
        return;
    }
    //
    let a_steps = (0.1 * du).max(100.0 * my_utol);
    let a_nb_particles = (32.0 * (the_u2 - the_u1) / du) as i32;
    let a_low_border = [the_u1];
    let a_upp_border = [the_u2];
    //
    //
    // math_PSO aFinder(&aFuncMV, aLowBorder, aUppBorder, aSteps, aNbParticles);
    // aFinder.Perform(aSteps, aValue, aT);
    let a_t = pso_minimize(
        |x: &[f64]| a_func.borrow_mut().value(x[0]).unwrap_or(0.0),
        &a_low_border,
        &a_upp_border,
        a_nb_particles.max(8) as usize,
        100,
        *my_utol,
    );
    let a_value = a_func.borrow_mut().value(a_t[0]).unwrap_or(0.0);
    //
    an_opt_loc.perform(
        &mut DistFunctionMV {
            my_max_curv_lin_dist: &a_func,
        },
        (a_t[0] - a_steps).max(the_u1),
        a_t[0],
        (a_t[0] + a_steps).min(the_u2),
    );
    if an_opt_loc.is_done() {
        *the_max_defl = (-an_opt_loc.minimum()).sqrt();
        *the_umax = an_opt_loc.location();
        return;
    }

    *the_max_defl = (-a_value).sqrt();
    *the_umax = a_t[0];
}

