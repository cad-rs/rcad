//! OCCT Extrema_ExtPElC2d (TKGeomBase/Extrema/Extrema_ExtPElC2d.hxx L33-161
//! and Extrema_ExtPElC2d.cxx L34-405) — all the extremum distances between a
//! 2D point and an elementary 2D curve (line, circle, ellipse, hyperbola,
//! parabola).
//!
//! Result encoding: `mySqDist[4]` / `myIsMin[4]` / `myPoint[4]` fixed arrays
//! as in the OCCT header (L158-160).  The Ellipse perform runs
//! `math_TrigonometricFunctionRoots` (kernel `math::root::trig_function_roots`)
//! and the Hyperbola / Parabola performs run `math_DirectPolynomialRoots`
//! (kernel `math::direct_polynomial_roots`).
//!
//! gp payload mapping follows [`super::extrema_ext_elc2d`] (the gp_Parab2d
//! `Focal` reads `Parabola2d::focal_param / 2`, the Y directions are the
//! kernel 90-degree turns of the stored X directions).

use glam::DVec2;

use crate::base::extrema::POnCurve2d;
use crate::base::extrema_ext_elc::{elclib_adjust_periodic, elclib_circle2d_parameter, REAL_LAST};
use crate::base::extrema_ext_elc2d::{
    default_p_on_curve2d, elclib2d_circle_value, elclib2d_ellipse_value, elclib2d_hyperbola_value,
    elclib2d_parabola_value, hyperbola2d_y_dir, parabola2d_y_dir,
};
use crate::core::precision::{PCONFUSION, CONFUSION};
use crate::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};
use crate::math::direct_polynomial_roots::DirectPolynomialRoots;
use crate::math::root::trig_function_roots;

/// OCCT Extrema_ExtPElC2d (hxx L33-161).
#[derive(Debug, Clone)]
pub struct ExtremaExtPElC2d {
    /// hxx L156: bool myDone.
    my_done: bool,
    /// hxx L157: int myNbExt.
    my_nb_ext: usize,
    /// hxx L158: double mySqDist[4].
    my_sq_dist: [f64; 4],
    /// hxx L159: bool myIsMin[4].
    my_is_min: [bool; 4],
    /// hxx L160: Extrema_POnCurv2d myPoint[4].
    my_point: [POnCurve2d; 4],
}

impl ExtremaExtPElC2d {
    /// OCCT Extrema_ExtPElC2d() (cxx L34-44).
    pub fn new() -> Self {
        let mut this = ExtremaExtPElC2d {
            my_done: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 4],
            my_is_min: [false; 4],
            my_point: [
                default_p_on_curve2d(),
                default_p_on_curve2d(),
                default_p_on_curve2d(),
                default_p_on_curve2d(),
            ],
        };
        // cxx L39-43.
        for an_i in 0..4 {
            this.my_sq_dist[an_i] = REAL_LAST;
            this.my_is_min[an_i] = false;
        }
        this
    }

    /// OCCT Extrema_ExtPElC2d(const gp_Pnt2d& P, const gp_Lin2d& L, Tol, Uinf,
    /// Usup) (cxx L48-55).
    pub fn point_line(p: DVec2, l: &Line2d, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC2d::new();
        this.perform_point_line(p, l, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC2d(const gp_Pnt2d& P, const gp_Circ2d& C, Tol,
    /// Uinf, Usup) (cxx L84-91).
    pub fn point_circle(p: DVec2, c: &Circle2d, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC2d::new();
        this.perform_point_circle(p, c, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC2d(const gp_Pnt2d& P, const gp_Elips2d& E, Tol,
    /// Uinf, Usup) (cxx L158-165).
    pub fn point_ellipse(p: DVec2, e: &Ellipse2d, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC2d::new();
        this.perform_point_ellipse(p, e, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC2d(const gp_Pnt2d& P, const gp_Hypr2d& C, Tol,
    /// Uinf, Usup) (cxx L218-225).
    pub fn point_hyperbola(p: DVec2, c: &Hyperbola2d, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC2d::new();
        this.perform_point_hyperbola(p, c, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC2d(const gp_Pnt2d& P, const gp_Parab2d& C, Tol,
    /// Uinf, Usup) (cxx L289-296).
    pub fn point_parabola(p: DVec2, c: &Parabola2d, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC2d::new();
        this.perform_point_parabola(p, c, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC2d::Perform(P, gp_Lin2d, Tol, Uinf, Usup)
    /// (cxx L57-80).
    pub fn perform_point_line(&mut self, p: DVec2, l: &Line2d, tol: f64, u_inf: f64, u_sup: f64) {
        // cxx L63-65.
        self.my_done = true;
        self.my_nb_ext = 0;

        // cxx L67-70: gp_Vec2d V1 = L.Direction(); OR = L.Location();
        // gp_Vec2d V(OR, P); Mydist = V1.Dot(V).
        let a_v1 = l.direction;
        let a_or = l.origin;
        let a_v = p - a_or;
        let a_mydist = a_v1.dot(a_v);

        // cxx L71-79.
        if (a_mydist >= u_inf - tol) && (a_mydist <= u_sup + tol) {
            self.my_nb_ext = 1;
            let a_myp = a_or + a_mydist * a_v1;
            self.my_sq_dist[0] = p.distance_squared(a_myp);
            self.my_point[0] = POnCurve2d {
                param: a_mydist,
                point: a_myp,
            };
            self.my_is_min[0] = true;
        }
    }

    /// OCCT Extrema_ExtPElC2d::Perform(P, gp_Circ2d, Tol, Uinf, Usup)
    /// (cxx L93-154).
    pub fn perform_point_circle(&mut self, p: DVec2, c: &Circle2d, tol: f64, u_inf: f64, u_sup: f64) {
        // cxx L100-101.
        let a_oc = c.center;
        self.my_nb_ext = 0;

        // cxx L103-106: OC.IsEqual(P, Precision::Confusion()).
        if (a_oc - p).length() <= CONFUSION {
            self.my_done = false;
        } else {
            // cxx L112-118.
            self.my_done = true;
            // cxx L113: gp_Dir2d V(gp_Vec2d(P, OC)) — the direction from P
            // toward the center.
            let a_v = (a_oc - p).normalize();
            let a_radius = c.radius;
            let mut a_p1 = a_oc + a_radius * a_v;
            // cxx L116: U1 = ElCLib::Parameter(C, P1).
            let mut a_u1 = elclib_circle2d_parameter(c, a_p1);
            let mut a_u2 = a_u1 + std::f64::consts::PI;
            let mut a_p2 = a_oc - a_radius * a_v;

            // cxx L119-121: AdjustPeriodic(Uinf, Uinf + 2*PI, PConfusion,
            // myuinf, U1) and the same for U2 (myuinf is the throwaway
            // companion value of the OCCT 5-argument form).
            let mut a_myuinf = u_inf;
            elclib_adjust_periodic(
                u_inf,
                u_inf + std::f64::consts::PI + std::f64::consts::PI,
                PCONFUSION,
                &mut a_u1,
                &mut a_myuinf,
            );
            elclib_adjust_periodic(
                u_inf,
                u_inf + std::f64::consts::PI + std::f64::consts::PI,
                PCONFUSION,
                &mut a_u2,
                &mut a_myuinf,
            );

            // cxx L122-127: the snap of U1 back to Uinf when it lands one
            // period beyond the window; P1 is recomputed from the angle.
            let a_two_pi = std::f64::consts::PI + std::f64::consts::PI;
            if ((a_u1 - a_two_pi - u_inf) < tol) && ((a_u1 - a_two_pi - u_inf) > -tol) {
                a_u1 = u_inf;
                a_p1 = elclib2d_circle_value(a_u1, c);
            }

            // cxx L129-134.
            if ((a_u2 - a_two_pi - u_inf) < tol) && ((a_u2 - a_two_pi - u_inf) > -tol) {
                a_u2 = u_inf;
                a_p2 = elclib2d_circle_value(a_u2, c);
            }

            // cxx L136-143.
            if ((u_inf - a_u1) < tol) && ((a_u1 - u_sup) < tol) {
                self.my_sq_dist[0] = p.distance_squared(a_p1);
                self.my_point[0] = POnCurve2d {
                    param: a_u1,
                    point: a_p1,
                };
                self.my_is_min[0] = true;
                self.my_nb_ext += 1;
            }

            // cxx L145-152.
            if ((u_inf - a_u2) < tol) && ((a_u2 - u_sup) < tol) {
                self.my_sq_dist[self.my_nb_ext] = p.distance_squared(a_p2);
                self.my_point[self.my_nb_ext] = POnCurve2d {
                    param: a_u2,
                    point: a_p2,
                };
                self.my_is_min[self.my_nb_ext] = true;
                self.my_nb_ext += 1;
            }
        }
    }

    /// OCCT Extrema_ExtPElC2d::Perform(P, gp_Elips2d, Tol, Uinf, Usup)
    /// (cxx L167-214).
    pub fn perform_point_ellipse(
        &mut self,
        p: DVec2,
        e: &Ellipse2d,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        // cxx L173-177.
        self.my_done = false;
        self.my_nb_ext = 0;
        let a_or = e.center;

        // cxx L179-184: A = MajorRadius, B = MinorRadius; the degenerate
        // center-point / circle-like arm stays not-done.
        let a = e.major_radius;
        let b = e.minor_radius;
        let a_v = p - a_or;

        if (a_or - p).length() <= CONFUSION && (a - b).abs() <= tol {
            return;
        } else {
            // cxx L190-191: X = V.Dot(XAxis dir), Y = V.Dot(YAxis dir).
            let a_x = a_v.dot(e.major_dir);
            let a_y = a_v.dot(e.minor_dir);

            // cxx L193: math_TrigonometricFunctionRoots
            // Sol(0., (B*B - A*A)/2., -B*Y, A*X, 0., Uinf, Usup).
            let a_sol = trig_function_roots(0.0, (b * b - a * a) / 2.0, -b * a_y, a * a_x, 0.0, u_inf, u_sup);

            // cxx L195-198.
            if !a_sol.done {
                return;
            }

            // cxx L199-212.
            let a_nb_sol = a_sol.roots.len();
            self.my_nb_ext = 0;
            for a_no_sol in 1..=a_nb_sol {
                let a_us = a_sol.roots[a_no_sol - 1];
                let a_cu = elclib2d_ellipse_value(a_us, e);
                self.my_sq_dist[self.my_nb_ext] = a_cu.distance_squared(p);
                // cxx L208: myIsMin[myNbExt] = (NoSol == 0) — the OCCT index
                // quirk (NoSol runs from 1), preserved verbatim.
                self.my_is_min[self.my_nb_ext] = a_no_sol == 0;
                self.my_point[self.my_nb_ext] = POnCurve2d {
                    param: a_us,
                    point: a_cu,
                };
                self.my_nb_ext += 1;
            }
            self.my_done = true;
        }
    }

    /// OCCT Extrema_ExtPElC2d::Perform(P, gp_Hypr2d, Tol, Uinf, Usup)
    /// (cxx L227-285).
    pub fn perform_point_hyperbola(
        &mut self,
        p: DVec2,
        h: &Hyperbola2d,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        // cxx L233-235.
        let a_o = h.center;
        self.my_done = false;
        self.my_nb_ext = 0;

        // cxx L237-244.
        let a_r = h.semi_major;
        let a_r_minor = h.semi_minor;
        let a_opp = p - a_o;
        let a_tol2 = tol * tol;
        let a_x = a_opp.dot(h.major_dir);
        let a_y = a_opp.dot(hyperbola2d_y_dir(h));
        let a_c1 = (a_r * a_r + a_r_minor * a_r_minor) / 4.0;

        // cxx L244: math_DirectPolynomialRoots
        // Sol(C1, -(X*R + Y*r)/2., 0., (X*R - Y*r)/2., -C1).
        let a_sol = DirectPolynomialRoots::new_quartic(
            a_c1,
            -(a_x * a_r + a_y * a_r_minor) / 2.0,
            0.0,
            (a_x * a_r - a_y * a_r_minor) / 2.0,
            -a_c1,
        );

        // cxx L245-248.
        if !a_sol.is_done() {
            return;
        }

        // cxx L249-284.
        let mut a_tb_ext = [DVec2::ZERO; 4];
        let a_nb_sol = a_sol.nb_solutions();
        for a_no_sol in 1..=a_nb_sol {
            let a_vs = a_sol.value(a_no_sol);
            if a_vs > 0.0 {
                let a_us = a_vs.ln();
                if (a_us >= u_inf) && (a_us <= u_sup) {
                    let a_cu = elclib2d_hyperbola_value(a_us, h);
                    let mut a_deja_enr = false;
                    for a_no_ext in 0..self.my_nb_ext {
                        if a_tb_ext[a_no_ext].distance_squared(a_cu) < a_tol2 {
                            a_deja_enr = true;
                            break;
                        }
                    }
                    if !a_deja_enr {
                        // The quartic answers at most 4 roots, so myNbExt
                        // stays below the TbExt[4] bound of cxx L254.
                        a_tb_ext[self.my_nb_ext] = a_cu;
                        self.my_sq_dist[self.my_nb_ext] = a_cu.distance_squared(p);
                        // cxx L277: myIsMin[myNbExt] = (NoSol == 0) — the
                        // OCCT index quirk (NoSol runs from 1), preserved.
                        self.my_is_min[self.my_nb_ext] = a_no_sol == 0;
                        self.my_point[self.my_nb_ext] = POnCurve2d {
                            param: a_us,
                            point: a_cu,
                        };
                        self.my_nb_ext += 1;
                    }
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT Extrema_ExtPElC2d::Perform(P, gp_Parab2d, Tol, Uinf, Usup)
    /// (cxx L298-351).
    pub fn perform_point_parabola(
        &mut self,
        p: DVec2,
        c: &Parabola2d,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        // cxx L304-306.
        self.my_done = false;
        self.my_nb_ext = 0;
        let a_o = c.origin;

        // cxx L308-312: Tol2; F = C.Focal() (the rcad focal_param / 2);
        // X, Y the projections on the position X / Y directions.
        let a_tol2 = tol * tol;
        let a_f = c.focal_param / 2.0;
        let a_opp = p - a_o;
        let a_x = a_opp.dot(c.axis_dir);
        let a_y = a_opp.dot(parabola2d_y_dir(c));

        // cxx L314: math_DirectPolynomialRoots
        // Sol(1./(4.*F), 0., 2.*F - X, -2.*F*Y).
        let a_sol = DirectPolynomialRoots::new_cubic(1.0 / (4.0 * a_f), 0.0, 2.0 * a_f - a_x, -2.0 * a_f * a_y);

        // cxx L315-318.
        if !a_sol.is_done() {
            return;
        }

        // cxx L319-350.
        let mut a_tb_ext = [DVec2::ZERO; 3];
        let a_nb_sol = a_sol.nb_solutions();
        for a_no_sol in 1..=a_nb_sol {
            let a_us = a_sol.value(a_no_sol);
            if (a_us >= u_inf) && (a_us <= u_sup) {
                let a_cu = elclib2d_parabola_value(a_us, c);
                let mut a_deja_enr = false;
                for a_no_ext in 0..self.my_nb_ext {
                    if a_tb_ext[a_no_ext].distance_squared(a_cu) < a_tol2 {
                        a_deja_enr = true;
                        break;
                    }
                }
                if !a_deja_enr {
                    // The cubic answers at most 3 roots, so myNbExt stays
                    // below the TbExt[3] bound of cxx L324.
                    a_tb_ext[self.my_nb_ext] = a_cu;
                    self.my_sq_dist[self.my_nb_ext] = a_cu.distance_squared(p);
                    // cxx L344: myIsMin[myNbExt] = (NoSol == 0) — the OCCT
                    // index quirk (NoSol runs from 1), preserved.
                    self.my_is_min[self.my_nb_ext] = a_no_sol == 0;
                    self.my_point[self.my_nb_ext] = POnCurve2d {
                        param: a_us,
                        point: a_cu,
                    };
                    self.my_nb_ext += 1;
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT IsDone() (cxx L355-358).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT NbExt() (cxx L362-369).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_nb_ext
    }

    /// OCCT SquareDistance(N) (cxx L373-380) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT IsMin(N) (cxx L384-391) — 1-based.
    pub fn is_min(&self, n: usize) -> bool {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_is_min[n - 1]
    }

    /// OCCT Point(N) (cxx L395-402) — 1-based; the OCCT const& return is a
    /// copy in the rcad value encoding.
    pub fn point(&self, n: usize) -> POnCurve2d {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_point[n - 1].clone()
    }
}

impl Default for ExtremaExtPElC2d {
    fn default() -> Self {
        ExtremaExtPElC2d::new()
    }
}
