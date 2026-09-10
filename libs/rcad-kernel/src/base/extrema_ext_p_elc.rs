//! OCCT Extrema_ExtPElC (TKGeomBase/Extrema/Extrema_ExtPElC.hxx L35-161 and
//! Extrema_ExtPElC.cxx L31-533) — all the extremum distances between a point
//! and an elementary curve (line, circle, ellipse, hyperbola, parabola).
//!
//! Result encoding: `mySqDist[4]` / `myIsMin[4]` / `myPoint[4]` fixed arrays as
//! in the OCCT header (L158-160).

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema_ext_elc::{
    elclib_adjust_periodic, elclib_circle_value, elclib_ellipse_value, elclib_hyperbola_value,
    elclib_parabola_value, GP_RESOLUTION, PRECISION_INFINITE, REAL_LAST,
};
use crate::base::extrema_gp::vec_angle_with_ref;
use crate::core::precision::{ANGULAR, SQUARE_CONFUSION};
use crate::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};
use crate::math::direct_polynomial_roots::DirectPolynomialRoots;
use crate::math::root::trig_function_roots;

/// OCCT Extrema_ExtPElC (hxx L35-161).
#[derive(Debug, Clone)]
pub struct ExtremaExtPElC {
    /// hxx L156: bool myDone.
    my_done: bool,
    /// hxx L157: int myNbExt.
    my_nb_ext: usize,
    /// hxx L158: double mySqDist[4].
    my_sq_dist: [f64; 4],
    /// hxx L159: bool myIsMin[4].
    my_is_min: [bool; 4],
    /// hxx L160: Extrema_POnCurv myPoint[4].
    my_point: [POnCurve; 4],
}

/// The OCCT default-constructed Extrema_POnCurv.
fn default_p_on_curve() -> POnCurve {
    POnCurve {
        param: 0.0,
        point: DVec3::ZERO,
    }
}

impl ExtremaExtPElC {
    /// OCCT Extrema_ExtPElC() (cxx L33-43).
    pub fn new() -> Self {
        let mut this = ExtremaExtPElC {
            my_done: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 4],
            my_is_min: [false; 4],
            my_point: [default_p_on_curve(), default_p_on_curve(), default_p_on_curve(), default_p_on_curve()],
        };
        for i in 0..4 {
            this.my_sq_dist[i] = REAL_LAST;
            this.my_is_min[i] = false;
        }
        this
    }

    /// OCCT Extrema_ExtPElC(const gp_Pnt& P, const gp_Lin& L, Tol, Uinf, Usup)
    /// (cxx L47-54) with `Perform` (cxx L58-81).
    pub fn point_line(
        p: DVec3,
        l: &Line3,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) -> Self {
        let mut this = ExtremaExtPElC::new();
        this.perform_point_line(p, l, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC::Perform(P, gp_Lin, Tol, Uinf, Usup) (cxx L58-81).
    pub fn perform_point_line(
        &mut self,
        p: DVec3,
        l: &Line3,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        self.my_done = false;
        self.my_nb_ext = 0;
        let v1 = l.direction;
        let or = l.origin;
        let v = p - or;
        let my_dist = v1.dot(v);
        if (my_dist >= u_inf - tol) && (my_dist <= u_sup + tol) {
            let my_p = or + my_dist * v1;
            self.my_sq_dist[0] = (p - my_p).length_squared();
            self.my_point[0] = POnCurve {
                param: my_dist,
                point: my_p,
            };
            self.my_is_min[0] = true;
            self.my_nb_ext = 1;
            self.my_done = true;
        }
    }

    /// OCCT Extrema_ExtPElC(const gp_Pnt& P, const gp_Circ& C, Tol, Uinf,
    /// Usup) (cxx L83-90) with `Perform` (cxx L92-190).
    pub fn point_circle(p: DVec3, c: &Circle3, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC::new();
        this.perform_point_circle(p, c, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC::Perform(P, gp_Circ, Tol, Uinf, Usup) (cxx L92-190).
    pub fn perform_point_circle(
        &mut self,
        p: DVec3,
        c: &Circle3,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        self.my_done = false;
        self.my_nb_ext = 0;

        // 1- Projection of the point P in the plane of circle -> Pp (cxx L125-128).
        let o = c.center;
        let axe = c.normal;
        let trsl = axe * (-(p - o).dot(axe));
        let pp = p + trsl;

        // 2- u solutions in [0, 2*PI] (cxx L132-150).
        let o_pp = pp - o;
        if o_pp.length() < tol {
            return;
        }
        let mut u_sol = [0.0f64; 2];
        u_sol[0] = vec_angle_with_ref(c.x_dir, o_pp, axe); // -M_PI < U1 < M_PI

        // cxx L140-148.
        let a_ang_tol = ANGULAR;
        if u_sol[0] + std::f64::consts::PI < a_ang_tol {
            u_sol[0] = -std::f64::consts::PI;
        } else if u_sol[0] - std::f64::consts::PI > -a_ang_tol {
            u_sol[0] = std::f64::consts::PI;
        }
        u_sol[1] = u_sol[0] + std::f64::consts::PI;

        // cxx L152-171.
        let mut my_u_inf = u_inf;
        let a_r = c.radius;
        let tol_u = if a_r > GP_RESOLUTION {
            tol / a_r
        } else {
            PRECISION_INFINITE
        };
        elclib_adjust_periodic(
            u_inf,
            u_inf + 2.0 * std::f64::consts::PI,
            tol_u,
            &mut my_u_inf,
            &mut u_sol[0],
        );
        elclib_adjust_periodic(
            u_inf,
            u_inf + 2.0 * std::f64::consts::PI,
            tol_u,
            &mut my_u_inf,
            &mut u_sol[1],
        );
        if ((u_sol[0] - 2.0 * std::f64::consts::PI - u_inf) < tol_u)
            && ((u_sol[0] - 2.0 * std::f64::consts::PI - u_inf) > -tol_u)
        {
            u_sol[0] = u_inf;
        }
        if ((u_sol[1] - 2.0 * std::f64::consts::PI - u_inf) < tol_u)
            && ((u_sol[1] - 2.0 * std::f64::consts::PI - u_inf) > -tol_u)
        {
            u_sol[1] = u_inf;
        }

        // 3- extrema in [Umin, Umax] (cxx L175-188).
        for no_sol in 0..=1usize {
            let us = u_sol[no_sol];
            if ((u_inf - us) < tol_u) && ((us - u_sup) < tol_u) {
                let cu = elclib_circle_value(us, c);
                self.my_sq_dist[self.my_nb_ext] = (cu - p).length_squared();
                self.my_is_min[self.my_nb_ext] = no_sol == 0;
                self.my_point[self.my_nb_ext] = POnCurve { param: us, point: cu };
                self.my_nb_ext += 1;
            }
        }
        self.my_done = true;
    }

    /// OCCT Extrema_ExtPElC(const gp_Pnt& P, const gp_Elips& C, Tol, Uinf,
    /// Usup) (cxx L194-201) with `Perform` (cxx L203-281).
    pub fn point_ellipse(p: DVec3, c: &Ellipse3, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC::new();
        this.perform_point_ellipse(p, c, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC::Perform(P, gp_Elips, Tol, Uinf, Usup)
    /// (cxx L203-281).
    pub fn perform_point_ellipse(
        &mut self,
        p: DVec3,
        c: &Ellipse3,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        self.my_done = false;
        self.my_nb_ext = 0;

        // 1- Projection of P in the plane of the ellipse (cxx L232-235).
        let o = c.center;
        let axe = c.normal;
        let trsl = axe * (-(p - o).dot(axe));
        let pp = p + trsl;

        // 2- solutions (cxx L239-261).
        let a = c.major_radius;
        let b = c.minor_radius;
        let o_pp = pp - o;
        let o_pp_magn = o_pp.length();
        if o_pp_magn < tol && (a - b).abs() < tol {
            return;
        }
        let y_dir = c.normal.cross(c.major_dir);
        let x = o_pp.dot(c.major_dir);
        let y = o_pp.dot(y_dir);

        let ko2 = (b * b - a * a) / 2.0;
        let mut ko3 = -b * y;
        let ko4 = a * x;
        if ko3.abs() < 1.0e-16 * ko2.abs().max(ko3.abs()) {
            ko3 = 0.0;
        }

        // OCCT: math_TrigonometricFunctionRoots Sol(0., ko2, ko3, ko4, 0.,
        // Uinf, Usup).
        let sol = trig_function_roots(0.0, ko2, ko3, ko4, 0.0, u_inf, u_sup);
        if !sol.done {
            return;
        }
        let nb_sol = sol.roots.len();
        for no_sol in 1..=nb_sol {
            let us = sol.roots[no_sol - 1];
            let cu = elclib_ellipse_value(us, c);
            self.my_sq_dist[self.my_nb_ext] = (cu - p).length_squared();
            self.my_point[self.my_nb_ext] = POnCurve { param: us, point: cu };
            let cu1 = elclib_ellipse_value(us + 0.1, c);
            self.my_is_min[self.my_nb_ext] = self.my_sq_dist[self.my_nb_ext] < (cu1 - p).length_squared();
            self.my_nb_ext += 1;
        }
        self.my_done = true;
    }

    /// OCCT Extrema_ExtPElC(const gp_Pnt& P, const gp_Hypr& C, Tol, Uinf,
    /// Usup) (cxx L285-292) with `Perform` (cxx L294-389).
    pub fn point_hyperbola(p: DVec3, c: &Hyperbola3, tol: f64, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC::new();
        this.perform_point_hyperbola(p, c, tol, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC::Perform(P, gp_Hypr, Tol, Uinf, Usup)
    /// (cxx L294-389).
    pub fn perform_point_hyperbola(
        &mut self,
        p: DVec3,
        c: &Hyperbola3,
        tol: f64,
        u_inf: f64,
        u_sup: f64,
    ) {
        self.my_done = false;
        self.my_nb_ext = 0;

        // 1- Projection of P in the plane of the hyperbola (cxx L333-336).
        let o = c.center;
        let axe = c.normal;
        let trsl = axe * (-(p - o).dot(axe));
        let pp = p + trsl;

        // 2- solutions (cxx L340-348).
        let tol2 = tol * tol;
        let r = c.semi_major;
        let r_minor = c.semi_minor;
        let o_pp = pp - o;
        let y_dir = c.normal.cross(c.major_dir);
        let x = o_pp.dot(c.major_dir);
        let y = o_pp.dot(y_dir);

        let c1 = (r * r + r_minor * r_minor) / 4.0;
        let sol = DirectPolynomialRoots::new_quartic(
            c1,
            -(x * r + y * r_minor) / 2.0,
            0.0,
            (x * r - y * r_minor) / 2.0,
            -c1,
        );
        if !sol.is_done() {
            return;
        }
        let mut tb_ext = [DVec3::ZERO; 4];
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let vs = sol.value(no_sol);
            if vs > 0.0 {
                // cxx L364-365.
                let us = vs.ln();
                if (us >= u_inf) && (us <= u_sup) {
                    let cu = elclib_hyperbola_value(us, c);
                    // cxx L368-376: duplicate check.
                    let mut deja_enr = false;
                    for no_ext in 0..self.my_nb_ext {
                        if (tb_ext[no_ext] - cu).length_squared() < tol2 {
                            deja_enr = true;
                            break;
                        }
                    }
                    if !deja_enr {
                        tb_ext[self.my_nb_ext] = cu;
                        self.my_sq_dist[self.my_nb_ext] = (cu - p).length_squared();
                        let cu1 = elclib_hyperbola_value(us + 1.0, c);
                        self.my_is_min[self.my_nb_ext] =
                            self.my_sq_dist[self.my_nb_ext] < (cu1 - p).length_squared();
                        self.my_point[self.my_nb_ext] = POnCurv_or_default(us, cu);
                        self.my_nb_ext += 1;
                    }
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT Extrema_ExtPElC(const gp_Pnt& P, const gp_Parab& C, Tol, Uinf,
    /// Usup) (cxx L393-400) with `Perform` (cxx L402-480).
    pub fn point_parabola(p: DVec3, c: &Parabola3, u_inf: f64, u_sup: f64) -> Self {
        let mut this = ExtremaExtPElC::new();
        this.perform_point_parabola(p, c, u_inf, u_sup);
        this
    }

    /// OCCT Extrema_ExtPElC::Perform(P, gp_Parab, Tol, Uinf, Usup)
    /// (cxx L402-480).
    pub fn perform_point_parabola(
        &mut self,
        p: DVec3,
        c: &Parabola3,
        u_inf: f64,
        u_sup: f64,
    ) {
        self.my_done = false;
        self.my_nb_ext = 0;

        // 1- Projection of P in the plane of the parabola (cxx L430-435).
        let o = c.vertex;
        let axe = c.normal;
        let trsl = axe * (-(p - o).dot(axe));
        let pp = p + trsl;

        // 2- solutions (cxx L437-443).
        let f = c.focal_param;
        let o_pp = pp - o;
        let y_dir = c.normal.cross(c.axis_dir);
        let x = o_pp.dot(c.axis_dir);
        let y = o_pp.dot(y_dir);
        let sol = DirectPolynomialRoots::new_cubic(1.0 / (4.0 * f), 0.0, 2.0 * f - x, -2.0 * f * y);
        if !sol.is_done() {
            return;
        }
        let mut tb_ext = [DVec3::ZERO; 3];
        let nb_sol = sol.nb_solutions();
        for no_sol in 1..=nb_sol {
            let us = sol.value(no_sol);
            if (us >= u_inf) && (us <= u_sup) {
                let cu = elclib_parabola_value(us, c);
                // cxx L460-468: duplicate check.
                let mut deja_enr = false;
                for no_ext in 0..self.my_nb_ext {
                    if (tb_ext[no_ext] - cu).length_squared() < SQUARE_CONFUSION {
                        deja_enr = true;
                        break;
                    }
                }
                if !deja_enr {
                    tb_ext[self.my_nb_ext] = cu;
                    self.my_sq_dist[self.my_nb_ext] = (cu - p).length_squared();
                    let cu1 = elclib_parabola_value(us + 1.0, c);
                    self.my_is_min[self.my_nb_ext] =
                        self.my_sq_dist[self.my_nb_ext] < (cu1 - p).length_squared();
                    self.my_point[self.my_nb_ext] = POnCurv_or_default(us, cu);
                    self.my_nb_ext += 1;
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT IsDone() (cxx L484-487).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT NbExt() (cxx L491-498).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_nb_ext
    }

    /// OCCT SquareDistance(N) (cxx L502-509) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT IsMin(N) (cxx L513-520) — 1-based.
    pub fn is_min(&self, n: usize) -> bool {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_is_min[n - 1]
    }

    /// OCCT Point(N) (cxx L524-531) — 1-based.
    pub fn point(&self, n: usize) -> &POnCurve {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        &self.my_point[n - 1]
    }
}

/// The Extrema_POnCurv construction used by the storage loops.
fn POnCurv_or_default(param: f64, point: DVec3) -> POnCurve {
    POnCurve { param, point }
}

impl Default for ExtremaExtPElC {
    fn default() -> Self {
        ExtremaExtPElC::new()
    }
}
