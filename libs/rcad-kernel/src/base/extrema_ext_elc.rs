//! OCCT Extrema_ExtElC (TKGeomBase/Extrema/Extrema_ExtElC.hxx L35-93 and
//! Extrema_ExtElC.cxx L256-1197) — all the extremum distances between two
//! elementary curves (line-line, line-circle, line-ellipse, line-hyperbola,
//! line-parabola, circle-circle).
//!
//! Also carries the two file-local statics of the same translation unit:
//! - `ExtremaExtElC_TrigonometricRoots` (cxx L47-252) — the sorted roots of
//!   the trigonometric polynomial, code duplicated from IntAna_IntQuadQuad,
//! - `RefineDir` (cxx L1171-1197) — the direction snap to a coordinate axis.
//!
//! Result encoding: `mySqDist[6]` / `myPoint[6][2]` are fixed arrays exactly
//! as in the OCCT header (L91-92).

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::core::precision::{is_infinite_value, ANGULAR, CONFUSION, INFINITE_VALUE};
use crate::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};
use crate::math::direct_polynomial_roots::DirectPolynomialRoots;
use crate::math::root::trig_function_roots;

/// OCCT gp::Resolution() (gp.cxx) — 2.2250738585072014e-308.
pub(crate) const GP_RESOLUTION: f64 = 2.2250738585072014e-308;

/// OCCT RealEpsilon() — DBL_EPSILON.
pub(crate) const REAL_EPSILON: f64 = f64::EPSILON;

/// OCCT RealLast().
pub(crate) const REAL_LAST: f64 = f64::MAX;

/// OCCT RealFirst().
pub(crate) const REAL_FIRST: f64 = -f64::MAX;

/// OCCT `Precision::Infinite()`.
pub(crate) const PRECISION_INFINITE: f64 = INFINITE_VALUE;

// =============================================================================
// ExtremaExtElC_TrigonometricRoots (cxx L47-252)
// =============================================================================

/// OCCT ExtremaExtElC_TrigonometricRoots (cxx L47-117) — internal class
/// providing the sorted roots of a trigonometric polynomial.
pub(crate) struct ExtremaExtElCTrigonometricRoots {
    /// cxx L50: double Roots[4].
    roots: [f64; 4],
    /// cxx L51: bool done.
    done: bool,
    /// cxx L52: int NbRoots.
    nb_roots: usize,
    /// cxx L53: bool infinite_roots.
    infinite_roots: bool,
}

impl ExtremaExtElCTrigonometricRoots {
    /// OCCT ExtremaExtElC_TrigonometricRoots(CC, SC, C, S, Cte, Binf, Bsup)
    /// (cxx L121-252).
    pub(crate) fn new(
        cc_in: f64,
        sc_in: f64,
        c_in: f64,
        s_in: f64,
        cte_in: f64,
        binf: f64,
        bsup: f64,
    ) -> Self {
        // cxx L131-140.
        let mut nbessai = 1;
        let mut cc = cc_in;
        let mut sc = sc_in;
        let mut c = c_in;
        let mut s = s_in;
        let mut cte = cte_in;
        let mut this = ExtremaExtElCTrigonometricRoots {
            roots: [0.0; 4],
            done: false,
            nb_roots: 0,
            infinite_roots: false,
        };

        // cxx L141: while (nbessai <= 2 && !done)
        while nbessai <= 2 && !this.done {
            // OCCT: math_TrigonometricFunctionRoots MTFR(cc, sc, c, s, cte,
            // Binf, Bsup).
            let mtfr = trig_function_roots(cc, sc, c, s, cte, binf, bsup);
            if mtfr.done {
                this.done = true;
                if mtfr.infinite {
                    this.infinite_roots = true;
                } else {
                    // cxx L155-223 (else #1).
                    let a_two_pi = std::f64::consts::PI + std::f64::consts::PI;
                    this.nb_roots = mtfr.roots.len();
                    for i in 0..this.nb_roots {
                        this.roots[i] = mtfr.roots[i];
                        if this.roots[i] < 0.0 {
                            this.roots[i] += a_two_pi;
                        }
                        if this.roots[i] > a_two_pi {
                            this.roots[i] -= a_two_pi;
                        }
                    }

                    // cxx L174-179: reliability check on the direct search.
                    let mut a_max_coef = cc.max(sc);
                    a_max_coef = a_max_coef.max(c);
                    a_max_coef = a_max_coef.max(s);
                    a_max_coef = a_max_coef.max(cte);
                    let a_precision = 1.0e-8f64.max(1.0e-12 * a_max_coef);

                    let sv_nb_roots = this.nb_roots;
                    for i in 0..sv_nb_roots {
                        let co = this.roots[i].cos();
                        let si = this.roots[i].sin();
                        let y = co * (cc * co + (sc + sc) * si + c) + s * si + cte;
                        // cxx L189-193.
                        if y.abs() > a_precision {
                            this.nb_roots -= 1;
                            this.roots[i] = 1000.0;
                        }
                    }

                    // cxx L196-211: bubble sort.
                    loop {
                        let mut triee = true;
                        let mut i = 1usize;
                        let mut j = 0usize;
                        while i < sv_nb_roots {
                            if this.roots[i] < this.roots[j] {
                                triee = false;
                                this.roots.swap(i, j);
                            }
                            i += 1;
                            j += 1;
                        }
                        if triee {
                            break;
                        }
                    }

                    // cxx L213-223.
                    this.infinite_roots = false;
                    if this.nb_roots == 0 {
                        if cc.abs() + sc.abs() + c.abs() + s.abs() < 1e-10 {
                            if cte.abs() < 1e-10 {
                                this.infinite_roots = true;
                            }
                        }
                    }
                }
            } else {
                // cxx L227-250: try to set very small coefficients to ZERO.
                if cc_in.abs() < 1e-10 {
                    cc = 0.0;
                }
                if sc_in.abs() < 1e-10 {
                    sc = 0.0;
                }
                if c_in.abs() < 1e-10 {
                    c = 0.0;
                }
                if s_in.abs() < 1e-10 {
                    s = 0.0;
                }
                if cte_in.abs() < 1e-10 {
                    cte = 0.0;
                }
                nbessai += 1;
            }
        }
        this
    }

    /// OCCT IsDone() (cxx L65).
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsARoot(U) (cxx L68-86).
    pub(crate) fn is_a_root(&self, u: f64) -> bool {
        let a_eps = REAL_EPSILON;
        let pi_p_pi = std::f64::consts::PI + std::f64::consts::PI;
        for i in 0..self.nb_roots {
            if (u - self.roots[i]).abs() <= a_eps {
                return true;
            }
            if (u - self.roots[i] - pi_p_pi).abs() <= a_eps {
                return true;
            }
        }
        false
    }

    /// OCCT NbSolutions() (cxx L89-95).
    pub(crate) fn nb_solutions(&self) -> usize {
        self.nb_roots
    }

    /// OCCT InfiniteRoots() (cxx L99-106).
    pub(crate) fn infinite_roots(&self) -> bool {
        self.infinite_roots
    }

    /// OCCT Value(N) (cxx L109-116) — 1-based.
    pub(crate) fn value(&self, n: usize) -> f64 {
        self.roots[n - 1]
    }
}

// =============================================================================
// RefineDir (cxx L1171-1197)
// =============================================================================

/// OCCT static RefineDir(gp_Dir&) (cxx L1171-1197) — snap a direction whose
/// component is within RealEpsilon() of ±1 onto the exact axis.
pub(crate) fn refine_dir(a_dir: DVec3) -> DVec3 {
    let i_k = 3usize;
    let a_eps = REAL_EPSILON;
    let mut a_cx = [a_dir.x, a_dir.y, a_dir.z];

    for i in 0..i_k {
        let a_one = if a_cx[i] > 0.0 { 1.0 } else { -1.0 };
        let a_x1 = a_one - a_eps;
        let a_x2 = a_one + a_eps;
        if a_cx[i] > a_x1 && a_cx[i] < a_x2 {
            let j = (i + 1) % i_k;
            let k = (i + 2) % i_k;
            a_cx[i] = a_one;
            a_cx[j] = 0.0;
            a_cx[k] = 0.0;
            return DVec3::new(a_cx[0], a_cx[1], a_cx[2]);
        }
    }
    a_dir
}

// =============================================================================
// ElCLib evaluation on the elementary curves (ElCLib.cxx)
// =============================================================================

/// OCCT ElCLib::LineValue(U, Pos) / ElCLib::Value(U, gp_Lin).
pub(crate) fn elclib_line_value(u: f64, lin: &Line3) -> DVec3 {
    lin.origin + u * lin.direction
}

/// OCCT ElCLib::CircleValue(U, Pos, R) / ElCLib::Value(U, gp_Circ).
pub(crate) fn elclib_circle_value(u: f64, circ: &Circle3) -> DVec3 {
    circ.center + circ.radius * (u.cos() * circ.x_dir + u.sin() * circ.y_dir)
}

/// OCCT ElCLib::EllipseValue(U, Pos, A, B) / ElCLib::Value(U, gp_Elips).
pub(crate) fn elclib_ellipse_value(u: f64, elips: &Ellipse3) -> DVec3 {
    let y_dir = elips.normal.cross(elips.major_dir);
    elips.center
        + elips.major_radius * u.cos() * elips.major_dir
        + elips.minor_radius * u.sin() * y_dir
}

/// OCCT ElCLib::HyperbolaValue(U, Pos, A, B) / ElCLib::Value(U, gp_Hypr).
pub(crate) fn elclib_hyperbola_value(u: f64, hypr: &Hyperbola3) -> DVec3 {
    let y_dir = hypr.normal.cross(hypr.major_dir);
    hypr.center
        + hypr.semi_major * u.cosh() * hypr.major_dir
        + hypr.semi_minor * u.sinh() * y_dir
}

/// OCCT ElCLib::ParabolaValue(U, Pos, P) / ElCLib::Value(U, gp_Parab).
pub(crate) fn elclib_parabola_value(u: f64, parab: &Parabola3) -> DVec3 {
    parab.vertex
        + (u * u / (2.0 * parab.focal_param)) * parab.axis_dir
        + u * parab.normal.cross(parab.axis_dir).normalize_or_zero()
}

/// OCCT ElCLib::Parameter(gp_Lin, gp_Pnt) — the parameter of the projection of
/// P on the line.
pub(crate) fn elclib_line_parameter(lin: &Line3, p: DVec3) -> f64 {
    (p - lin.origin).dot(lin.direction)
}

/// OCCT ElCLib::Parameter(gp_Circ, gp_Pnt).
pub(crate) fn elclib_circle_parameter(circ: &Circle3, p: DVec3) -> f64 {
    let v = p - circ.center;
    let x = v.dot(circ.x_dir);
    let y = v.dot(circ.y_dir);
    y.atan2(x)
}

/// OCCT gp_Lin::SquareDistance(gp_Pnt) — squared distance from the point to
/// the line.
pub(crate) fn line_square_distance(lin: &Line3, p: DVec3) -> f64 {
    let v = p - lin.origin;
    let t = v.dot(lin.direction);
    (v - t * lin.direction).length_squared()
}

/// OCCT ElCLib::InPeriod(U, UFirst, ULast) — reuse of the kernel encoding.
pub(crate) fn elclib_in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    crate::math::el::in_period(u, u_first, u_last)
}

/// OCCT ElCLib::AdjustPeriodic (ElCLib.cxx L115-149).
pub(crate) fn elclib_adjust_periodic(
    u_first: f64,
    u_last: f64,
    preci: f64,
    u1: &mut f64,
    u2: &mut f64,
) {
    if is_infinite_value(u_first) || is_infinite_value(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    if a_period < epsilon_of(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}

/// OCCT Epsilon(theValue) (Precision.hxx) —
/// `Max(Abs(theValue), RealSmall()) * RealEpsilon()`.
pub(crate) fn epsilon_of(the_value: f64) -> f64 {
    the_value.abs().max(2.2250738585072014e-308) * f64::EPSILON
}

// =============================================================================
// Extrema_ExtElC (hxx L35-93, cxx L256-1197)
// =============================================================================

/// OCCT Extrema_ExtElC (hxx L35-93) — all the distances between two
/// elementary curves.
#[derive(Debug, Clone)]
pub struct ExtremaExtElC {
    /// hxx L88: bool myDone.
    my_done: bool,
    /// hxx L89: bool myIsPar.
    my_is_par: bool,
    /// hxx L90: int myNbExt.
    my_nb_ext: usize,
    /// hxx L91: double mySqDist[6].
    my_sq_dist: [f64; 6],
    /// hxx L92: Extrema_POnCurv myPoint[6][2].
    my_point: [[POnCurve; 2]; 6],
}

/// The OCCT default-constructed Extrema_POnCurv (parameter 0, point (0,0,0)).
pub(crate) fn default_p_on_curve() -> POnCurve {
    POnCurve {
        param: 0.0,
        point: DVec3::ZERO,
    }
}

impl ExtremaExtElC {
    /// OCCT Extrema_ExtElC() (cxx L256-265).
    pub fn new() -> Self {
        let mut this = ExtremaExtElC {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 6],
            my_point: [
                [default_p_on_curve(), default_p_on_curve()],
                [default_p_on_curve(), default_p_on_curve()],
                [default_p_on_curve(), default_p_on_curve()],
                [default_p_on_curve(), default_p_on_curve()],
                [default_p_on_curve(), default_p_on_curve()],
                [default_p_on_curve(), default_p_on_curve()],
            ],
        };
        // cxx L261-264.
        for an_idx in 0..6 {
            this.my_sq_dist[an_idx] = REAL_LAST;
        }
        this
    }

    /// The OCCT member initializer repeated at the head of every computing
    /// constructor (cxx L315-321, L477-483, L660-666, L793-799, L895-901,
    /// L962-968).
    fn reset(&mut self) {
        self.my_is_par = false;
        self.my_done = false;
        self.my_nb_ext = 0;
        for an_idx in 0..6 {
            self.my_sq_dist[an_idx] = REAL_LAST;
        }
    }

    /// OCCT Extrema_ExtElC(const gp_Lin& C1, const gp_Lin& C2, AngTol)
    /// (cxx L268-357).
    pub fn line_line(the_c1: &Line3, the_c2: &Line3, _ang_tol: f64) -> Self {
        let mut this = ExtremaExtElC::new();
        this.reset();

        // cxx L323-339.
        let a_d1 = the_c1.direction;
        let a_d2 = the_c2.direction;
        let a_cos_a = a_d1.dot(a_d2);
        let a_sq_sin_a = 1.0 - a_cos_a * a_cos_a;
        let mut a_u1 = 0.0;
        let mut a_u2 = 0.0;
        if a_sq_sin_a < GP_RESOLUTION || a_cos_a.abs() > ANGULAR.cos() {
            this.my_is_par = true;
        } else {
            let a_l1l2 = the_c2.origin - the_c1.origin;
            let a_d1l = a_d1.dot(a_l1l2);
            let a_d2l = a_d2.dot(a_l1l2);
            a_u1 = (a_d1l - a_cos_a * a_d2l) / a_sq_sin_a;
            a_u2 = (a_cos_a * a_d1l - a_d2l) / a_sq_sin_a;
            this.my_is_par = is_infinite_value(a_u1) || is_infinite_value(a_u2);
        }

        // cxx L341-347.
        if this.my_is_par {
            this.my_sq_dist[0] = line_square_distance(the_c2, the_c1.origin);
            this.my_nb_ext = 1;
            this.my_done = true;
            return this;
        }

        // cxx L351-356.
        let a_p1 = elclib_line_value(a_u1, the_c1);
        let a_p2 = elclib_line_value(a_u2, the_c2);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve {
            param: a_u2,
            point: a_p2,
        };
        this.my_nb_ext = 1;
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC(const gp_Lin& C1, const gp_Circ& C2, Tol)
    /// (cxx L471-623) with the protected `PlanarLineCircleExtrema`
    /// (cxx L361-439).
    pub fn line_circle(c1: &Line3, c2: &Circle3, _tol: f64) -> Self {
        let mut this = ExtremaExtElC::new();
        this.reset();

        // cxx L485-488.
        if this.planar_line_circle_extrema(c1, c2) {
            return this;
        }

        // cxx L490-503: T1 in the reference of the circle.
        let d = c1.direction;
        let d1 = d;
        let x2 = c2.x_dir;
        let y2 = c2.y_dir;
        let z2 = c2.normal;
        let mut dx = d.dot(x2);
        let mut dy = d.dot(y2);
        let mut dz = d.dot(z2);
        // cxx L500-502: D.SetCoord(Dx, Dy, Dz); RefineDir(D); D.Coord(Dx, Dy, Dz).
        let a_d = refine_dir(DVec3::new(dx, dy, dz));
        dx = a_d.x;
        dy = a_d.y;
        dz = a_d.z;

        // cxx L504-523: V in the reference frame of the circle.
        let o1 = c1.origin;
        let o2 = c2.center;
        let o2o1_raw = o1 - o2;
        let a_tol_ro2o1 = GP_RESOLUTION;
        let a_ro2o1 = o2o1_raw.length();
        let o2o1 = if a_ro2o1 > a_tol_ro2o1 {
            // cxx L513-518.
            let a_do2o1 = refine_dir(DVec3::new(
                o2o1_raw.dot(x2),
                o2o1_raw.dot(y2),
                o2o1_raw.dot(z2),
            ));
            a_ro2o1 * a_do2o1
        } else {
            // cxx L522.
            DVec3::new(o2o1_raw.dot(x2), o2o1_raw.dot(y2), o2o1_raw.dot(z2))
        };

        // cxx L525.
        let a_d_xyz = DVec3::new(dx, dy, dz);
        let v_xyz = a_d_xyz * o2o1.dot(a_d_xyz) - o2o1;

        // cxx L556-585.
        let a_tol = 1.0e-12;
        let r = c2.radius;
        let mut a5 = r * dx * dy;
        let mut a1 = -2.0 * a5;
        let mut a2 = 0.5 * r * (dx * dx - dy * dy);
        let mut a3 = v_xyz.y;
        let mut a4 = -v_xyz.x;
        if a1 >= -a_tol && a1 <= a_tol {
            a1 = 0.0;
        }
        if a2 >= -a_tol && a2 <= a_tol {
            a2 = 0.0;
        }
        if a3 >= -a_tol && a3 <= a_tol {
            a3 = 0.0;
        }
        if a4 >= -a_tol && a4 <= a_tol {
            a4 = 0.0;
        }
        if a5 >= -a_tol && a5 <= a_tol {
            a5 = 0.0;
        }

        // cxx L588-600.
        let sol = ExtremaExtElCTrigonometricRoots::new(
            a1,
            a2,
            a3,
            a4,
            a5,
            0.0,
            std::f64::consts::PI + std::f64::consts::PI,
        );
        if !sol.is_done() {
            return this;
        }
        if sol.infinite_roots() {
            this.my_is_par = true;
            this.my_sq_dist[0] = r * r;
            this.my_nb_ext = 1;
            this.my_done = true;
            return this;
        }

        // cxx L601-622: storage of solutions.
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let u2 = sol.value(no_sol);
            let p2 = elclib_circle_value(u2, c2);
            let u1 = (p2 - o1).dot(d1);
            let p1 = elclib_line_value(u1, c1);
            this.my_sq_dist[this.my_nb_ext] = p1.distance_squared(p2);
            this.my_point[this.my_nb_ext][0] = POnCurve {
                param: u1,
                point: p1,
            };
            this.my_point[this.my_nb_ext][1] = POnCurve {
                param: u2,
                point: p2,
            };
            this.my_nb_ext += 1;
        }
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC::PlanarLineCircleExtrema (cxx L361-439) — GAP
    /// carrier.
    ///
    /// GAP dependency: the body needs `Extrema_ExtElC2d(gp_Lin2d, gp_Circ2d)`
    /// (Extrema_ExtElC2d.cxx) and `IntAna2d_AnaIntersection(gp_Lin2d, gp_Circ2d)`
    /// — neither is translated into the Extrema-language 2D form yet (only the
    /// line/ellipse form lives in the geomalgo IntAna consumers).
    ///
    /// OCCT failure path preserved (cxx L365-368): when the line is not
    /// parallel to the circle plane the method returns `false` and the caller
    /// (cxx L485) falls through into the general 3D trigonometric solve; the
    /// GAP returns `false` for every input, i.e. exactly the OCCT branch taken
    /// whenever the planarity test fails, and the general solve then runs
    /// (mathematically the same stationary-distance equation, so only the
    /// planar short-circuit's numerically-safer root set is lost).
    fn planar_line_circle_extrema(&mut self, _the_lin: &Line3, _the_circ: &Circle3) -> bool {
        // cxx L361-368.
        false
    }

    /// OCCT Extrema_ExtElC(const gp_Lin& C1, const gp_Elips& C2)
    /// (cxx L627-753).
    pub fn line_ellipse(c1: &Line3, c2: &Ellipse3) -> Self {
        let mut this = ExtremaExtElC::new();
        this.reset();

        // cxx L669-678.
        let d = c1.direction;
        let d1 = d;
        let x2 = c2.major_dir;
        let y2 = c2.normal.cross(c2.major_dir);
        let z2 = c2.normal;
        let dx = d.dot(x2);
        let dy = d.dot(y2);
        let _dz = d.dot(z2);

        // cxx L681-685.
        let o1 = c1.origin;
        let o2 = c2.center;
        let o2o1_raw = o1 - o2;
        let o2o1 = DVec3::new(o2o1_raw.dot(x2), o2o1_raw.dot(y2), o2o1_raw.dot(z2));
        let d_ref = DVec3::new(dx, dy, d.dot(z2));
        let v_xyz = d_ref * o2o1.dot(d_ref) - o2o1;

        // cxx L688-719.
        let maj_r = c2.major_radius;
        let min_r = c2.minor_radius;
        let mut a5 = maj_r * min_r * dx * dy;
        let mut a1 = -2.0 * a5;
        let r2 = maj_r * maj_r;
        let r2_minor = min_r * min_r;
        let mut a2 = (r2 * dx * dx - r2_minor * dy * dy - r2 + r2_minor) / 2.0;
        let mut a3 = min_r * v_xyz.y;
        let mut a4 = -maj_r * v_xyz.x;
        let a_eps = 1.0e-12;
        if a5.abs() <= a_eps {
            a5 = 0.0;
        }
        if a1.abs() <= a_eps {
            a1 = 0.0;
        }
        if a2.abs() <= a_eps {
            a2 = 0.0;
        }
        if a3.abs() <= a_eps {
            a3 = 0.0;
        }
        if a4.abs() <= a_eps {
            a4 = 0.0;
        }

        // cxx L721-735.
        let sol = ExtremaExtElCTrigonometricRoots::new(
            a1,
            a2,
            a3,
            a4,
            a5,
            0.0,
            std::f64::consts::PI + std::f64::consts::PI,
        );
        if !sol.is_done() {
            return this;
        }
        if sol.infinite_roots() {
            this.my_is_par = true;
            let a_p = elclib_ellipse_value(0.0, c2);
            this.my_sq_dist[0] = (a_p - c1.origin).length_squared();
            this.my_nb_ext = 1;
            this.my_done = true;
            return this;
        }

        // cxx L737-752.
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let u2 = sol.value(no_sol);
            let p2 = elclib_ellipse_value(u2, c2);
            let u1 = (p2 - o1).dot(d1);
            let p1 = elclib_line_value(u1, c1);
            this.my_sq_dist[this.my_nb_ext] = p1.distance_squared(p2);
            this.my_point[this.my_nb_ext][0] = POnCurve {
                param: u1,
                point: p1,
            };
            this.my_point[this.my_nb_ext][1] = POnCurve {
                param: u2,
                point: p2,
            };
            this.my_nb_ext += 1;
        }
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC(const gp_Lin& C1, const gp_Hypr& C2)
    /// (cxx L757-858).
    pub fn line_hyperbola(c1: &Line3, c2: &Hyperbola3) -> Self {
        let mut this = ExtremaExtElC::new();
        this.reset();

        // cxx L801-811.
        let d = c1.direction;
        let d1 = d;
        let x2 = c2.major_dir;
        let y2 = c2.normal.cross(c2.major_dir);
        let z2 = c2.normal;
        let dx = d.dot(x2);
        let dy = d.dot(y2);
        let dz = d.dot(z2);

        // cxx L813-820.
        let o1 = c1.origin;
        let o2 = c2.center;
        let o2o1_raw = o1 - o2;
        let o2o1 = DVec3::new(o2o1_raw.dot(x2), o2o1_raw.dot(y2), o2o1_raw.dot(z2));
        let d_ref = DVec3::new(dx, dy, dz);
        let v_xyz = d_ref * o2o1.dot(d_ref) - o2o1;
        let vx = v_xyz.x;
        let vy = v_xyz.y;

        // cxx L822-830.
        let r = c2.semi_major;
        let r_minor = c2.semi_minor;
        let a = -2.0 * r * r_minor * dx * dy;
        let b = -r * r * dx * dx - r_minor * r_minor * dy * dy + r * r + r_minor * r_minor;
        let a1 = a + b;
        let a2 = 2.0 * r * vx + 2.0 * r_minor * vy;
        let a4 = -2.0 * r * vx + 2.0 * r_minor * vy;
        let a5 = a - b;

        // OCCT: math_DirectPolynomialRoots Sol(A1, A2, 0.0, A4, A5).
        let sol = DirectPolynomialRoots::new_quartic(a1, a2, 0.0, a4, a5);
        if !sol.is_done() {
            return this;
        }

        // cxx L838-857.
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let v = sol.value(no_sol);
            if v > 0.0 {
                let u2 = v.ln();
                let p2 = elclib_hyperbola_value(u2, c2);
                let u1 = (p2 - o1).dot(d1);
                let p1 = elclib_line_value(u1, c1);
                this.my_sq_dist[this.my_nb_ext] = p1.distance_squared(p2);
                this.my_point[this.my_nb_ext][0] = POnCurve {
                    param: u1,
                    point: p1,
                };
                this.my_point[this.my_nb_ext][1] = POnCurve {
                    param: u2,
                    point: p2,
                };
                this.my_nb_ext += 1;
            }
        }
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC(const gp_Lin& C1, const gp_Parab& C2)
    /// (cxx L862-951).
    pub fn line_parabola(c1: &Line3, c2: &Parabola3) -> Self {
        let mut this = ExtremaExtElC::new();
        this.reset();

        // cxx L903-913.
        let d = c1.direction;
        let d1 = d;
        let x2 = c2.axis_dir;
        let y2 = c2.normal.cross(c2.axis_dir);
        let z2 = c2.normal;
        let dx = d.dot(x2);
        let dy = d.dot(y2);
        let dz = d.dot(z2);

        // cxx L915-920.
        let o1 = c1.origin;
        let o2 = c2.vertex;
        let o2o1_raw = o1 - o2;
        let o2o1 = DVec3::new(o2o1_raw.dot(x2), o2o1_raw.dot(y2), o2o1_raw.dot(z2));
        let d_ref = DVec3::new(dx, dy, dz);
        let v_xyz = d_ref * o2o1.dot(d_ref) - o2o1;

        // cxx L922-927.
        let p_par = c2.focal_param;
        let a1 = (1.0 - dx * dx) / (2.0 * p_par * p_par);
        let a2 = -3.0 * dx * dy / (2.0 * p_par);
        let a3 = 1.0 - dy * dy + v_xyz.x / p_par;
        let a4 = v_xyz.y;

        // OCCT: math_DirectPolynomialRoots Sol(A1, A2, A3, A4).
        let sol = DirectPolynomialRoots::new_cubic(a1, a2, a3, a4);
        if !sol.is_done() {
            return this;
        }

        // cxx L935-950.
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let u2 = sol.value(no_sol);
            let p2 = elclib_parabola_value(u2, c2);
            let u1 = (p2 - o1).dot(d1);
            let p1 = elclib_line_value(u1, c1);
            this.my_sq_dist[this.my_nb_ext] = p1.distance_squared(p2);
            this.my_point[this.my_nb_ext][0] = POnCurve {
                param: u1,
                point: p1,
            };
            this.my_point[this.my_nb_ext][1] = POnCurve {
                param: u2,
                point: p2,
            };
            this.my_nb_ext += 1;
        }
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC(const gp_Circ& C1, const gp_Circ& C2)
    /// (cxx L955-1113).
    pub fn circle_circle(c1: &Circle3, c2: &Circle3) -> Self {
        let mut this = ExtremaExtElC::new();
        this.reset();

        // cxx L970-972.
        let a_tol_a = ANGULAR;
        let a_tol_d = CONFUSION;
        let a_tol_d2 = a_tol_d * a_tol_d;

        // cxx L974-985.
        let a_pc1 = c1.center;
        let a_dc1 = c1.normal;
        let a_pc2 = c2.center;
        let a_dc2 = c2.normal;
        // OCCT: gp_Pln aPlc1(aPc1, aDc1); aD2 = aPlc1.SquareDistance(aPc2).
        let a_d2 = ((a_pc2 - a_pc1).dot(a_dc1)).powi(2);
        let b_is_same_plane = a_dc1.dot(a_dc2).abs() >= (1.0 - a_tol_a) && a_d2 < a_tol_d2;
        if !b_is_same_plane {
            return this;
        }

        // cxx L991-1002.
        let a_dc2 = (a_pc2 - a_pc1).length_squared();
        let b_is_same_axe = a_dc2 < a_tol_d2;
        if b_is_same_axe {
            this.my_is_par = true;
            this.my_nb_ext = 1;
            this.my_done = true;
            let a_dr = c1.radius - c2.radius;
            this.my_sq_dist[0] = a_dr * a_dr;
            return this;
        }

        // cxx L1004-1028.
        this.my_done = true;
        let a_r1_0 = c1.radius;
        let a_r2_0 = c2.radius;
        let mut j1 = 0usize;
        let mut j2 = 1usize;
        let mut a_c1 = c1.clone();
        let mut a_c2 = c2.clone();
        if a_r2_0 > a_r1_0 {
            j1 = 1;
            j2 = 0;
            a_c1 = c2.clone();
            a_c2 = c1.clone();
        }
        let a_r1 = a_c1.radius;
        let a_r2 = a_c2.radius;

        // cxx L1030-1035.
        let a_pc1 = a_c1.center;
        let a_pc2 = a_c2.center;
        let a_d12 = (a_pc1 - a_pc2).length();
        let a_dir12 = (a_pc2 - a_pc1).normalize_or_zero();

        // cxx L1037-1067: 1. Four common solutions.
        this.my_nb_ext = 4;
        let a_p11 = a_pc1 - a_r1 * a_dir12;
        let a_p12 = a_pc1 + a_r1 * a_dir12;
        let a_p21 = a_pc2 - a_r2 * a_dir12;
        let a_p22 = a_pc2 + a_r2 * a_dir12;
        let a_t11 = elclib_circle_parameter(&a_c1, a_p11);
        let a_t12 = elclib_circle_parameter(&a_c1, a_p12);
        let a_t21 = elclib_circle_parameter(&a_c2, a_p21);
        let a_t22 = elclib_circle_parameter(&a_c2, a_p22);

        this.my_point[0][j1] = POnCurve {
            param: a_t11,
            point: a_p11,
        };
        this.my_point[0][j2] = POnCurve {
            param: a_t21,
            point: a_p21,
        };
        this.my_sq_dist[0] = a_p11.distance_squared(a_p21);
        this.my_point[1][j1] = POnCurve {
            param: a_t11,
            point: a_p11,
        };
        this.my_point[1][j2] = POnCurve {
            param: a_t22,
            point: a_p22,
        };
        this.my_sq_dist[1] = a_p11.distance_squared(a_p22);
        this.my_point[2][j1] = POnCurve {
            param: a_t12,
            point: a_p12,
        };
        this.my_point[2][j2] = POnCurve {
            param: a_t21,
            point: a_p21,
        };
        this.my_sq_dist[2] = a_p12.distance_squared(a_p21);
        this.my_point[3][j1] = POnCurve {
            param: a_t12,
            point: a_p12,
        };
        this.my_point[3][j2] = POnCurve {
            param: a_t22,
            point: a_p22,
        };
        this.my_sq_dist[3] = a_p12.distance_squared(a_p22);

        // cxx L1069-1112: 2. Check for intersections.
        let b_out = a_d12 > (a_r1 + a_r2 + a_tol_d);
        let b_in = a_d12 < (a_r1 - a_r2 - a_tol_d);
        if !b_out && !b_in {
            let a_alpha = 0.5 * (a_r1 * a_r1 - a_r2 * a_r2 + a_d12 * a_d12) / a_d12;
            let mut a_val = a_r1 * a_r1 - a_alpha * a_alpha;
            if a_val < 0.0 {
                a_val = -a_val;
            }
            let a_beta = a_val.sqrt();
            let a_pt = a_pc1 + a_alpha * a_dir12;

            let a_d_lt = a_dc1.cross(a_dir12);
            let a_pl1 = a_pt + a_beta * a_d_lt;
            let a_pl2 = a_pt - a_beta * a_d_lt;

            let a_dist2 = a_pl1.distance_squared(a_pl2);
            let b_nb_ext6 = a_dist2 > a_tol_d2;

            this.my_nb_ext = 5;
            let mut a_t = [0.0f64; 2];
            a_t[j1] = elclib_circle_parameter(&a_c1, a_pl1);
            a_t[j2] = elclib_circle_parameter(&a_c2, a_pl1);
            this.my_point[4][j1] = POnCurve {
                param: a_t[j1],
                point: a_pl1,
            };
            this.my_point[4][j2] = POnCurve {
                param: a_t[j2],
                point: a_pl1,
            };
            this.my_sq_dist[4] = 0.0;

            if b_nb_ext6 {
                this.my_nb_ext = 6;
                a_t[j1] = elclib_circle_parameter(&a_c1, a_pl2);
                a_t[j2] = elclib_circle_parameter(&a_c2, a_pl2);
                this.my_point[5][j1] = POnCurve {
                    param: a_t[j1],
                    point: a_pl2,
                };
                this.my_point[5][j2] = POnCurve {
                    param: a_t[j2],
                    point: a_pl2,
                };
                this.my_sq_dist[5] = 0.0;
            }
        }
        this
    }

    /// OCCT IsDone() (cxx L1117-1120).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsParallel() (cxx L1124-1131).
    pub fn is_parallel(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_is_par
    }

    /// OCCT NbExt() (cxx L1135-1142).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_nb_ext
    }

    /// OCCT SquareDistance(N) (cxx L1146-1154) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT Points(N, P1, P2) (cxx L1158-1167) — 1-based.
    pub fn points(&self, n: usize, p1: &mut POnCurve, p2: &mut POnCurve) {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        *p1 = self.my_point[n - 1][0].clone();
        *p2 = self.my_point[n - 1][1].clone();
    }
}

impl Default for ExtremaExtElC {
    fn default() -> Self {
        ExtremaExtElC::new()
    }
}
