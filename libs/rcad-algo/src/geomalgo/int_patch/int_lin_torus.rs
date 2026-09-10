//! IntAna_IntLinTorus — analytic line × torus intersection.
//!
//! 1:1 translation of OCCT `IntAna_IntLinTorus.cxx` (L29-132): the line is
//! reparametrised so its location is nearest to the torus location, expressed
//! in the torus reference frame, and substituted into the torus implicit
//! equation (a quartic).  Every root is validated by re-evaluating the torus
//! point (ElSLib::Parameters/Value round trip) and rejecting solutions that
//! miss the surface by more than 1e-5 (L103).

use glam::DVec3;
use rcad_kernel::geom::{Line3, ToroidalSurface};
use rcad_kernel::math::direct_polynomial_roots::DirectPolynomialRoots;
use rcad_kernel::math::el::{elslib_torus_parameters, elslib_torus_value};

/// OCCT IntAna_IntLinTorus (IntAna_IntLinTorus.hxx) — the intersection
/// results: the parameter on the line, the (fi, theta) parameters on the
/// torus and the 3D point, 1-based accessors.
#[derive(Debug, Clone)]
pub struct IntLinTorus {
    done: bool,
    nbpt: usize,
    the_fi: Vec<f64>,
    the_param: Vec<f64>,
    the_theta: Vec<f64>,
    the_point: Vec<DVec3>,
}

impl IntLinTorus {
    /// OCCT IntAna_IntLinTorus() (cxx L29-36) — empty.
    pub fn new() -> Self {
        IntLinTorus {
            done: false,
            nbpt: 0,
            the_fi: Vec::new(),
            the_param: Vec::new(),
            the_theta: Vec::new(),
            the_point: Vec::new(),
        }
    }

    /// OCCT IntAna_IntLinTorus(L, T) (cxx L38-41).
    pub fn new_line_torus(l: &Line3, t: &ToroidalSurface) -> Self {
        let mut r = IntLinTorus::new();
        r.perform(l, t);
        r
    }

    /// OCCT Perform(L, T) (cxx L43-132).
    pub fn perform(&mut self, l: &Line3, t: &ToroidalSurface) {
        let pl = l.origin;
        let dl = l.direction;

        // Reparametrize the line: set its location as nearest to the
        // location of torus.
        let tor_loc = t.center;
        let param_of_new_pl = (tor_loc - pl).dot(dl);
        let new_pl = pl + param_of_new_pl * dl;

        //--------------------------------------------------------------
        //-- Coefficients of the line in the torus reference frame
        //--
        // OCCT trsf.SetTransformation(T.Position()) — the local frame of
        // the torus (X = ref_dir orthogonal to the axis, Y = axis ^ X).
        let x_ax = (t.ref_dir - t.axis * t.ref_dir.dot(t.axis)).normalize_or_zero();
        let y_ax = t.axis.cross(x_ax).normalize_or_zero();
        let rel = new_pl - t.center;
        let x0 = rel.dot(x_ax);
        let y0 = rel.dot(y_ax);
        let z0 = rel.dot(t.axis);
        let x1 = dl.dot(x_ax);
        let y1 = dl.dot(y_ax);
        let z1 = dl.dot(t.axis);

        let r_major = t.major_radius;
        let r_major2 = r_major * r_major;
        let r_minor = t.minor_radius;
        let r_minor2 = r_minor * r_minor;

        let a = x1 * x1 + y1 * y1 + z1 * z1;
        let b = 2.0 * (x1 * x0 + y1 * y0 + z1 * z0);
        let c = x0 * x0 + y0 * y0 + z0 * z0 - (r_major2 + r_minor2);

        let a4 = a * a;
        let a3 = 2.0 * a * b;
        let a2 = 2.0 * a * c + 4.0 * r_major2 * z1 * z1 + b * b;
        let a1 = 2.0 * b * c + 8.0 * r_major2 * z1 * z0;
        let a0 = c * c + 4.0 * r_major2 * (z0 * z0 - r_minor2);

        let mdpr = DirectPolynomialRoots::new_quartic(a4, a3, a2, a1, a0);
        if mdpr.is_done() {
            let mut the_param: Vec<f64> = Vec::new();
            let mut the_fi: Vec<f64> = Vec::new();
            let mut the_theta: Vec<f64> = Vec::new();
            let mut the_point: Vec<DVec3> = Vec::new();
            let n = mdpr.nb_solutions();
            let mut a_nb_bad_sol = 0usize;
            for i in 1..=n {
                let mut tt = mdpr.value(i);
                tt += param_of_new_pl;
                // OCCT ElCLib::Value(t, L).
                let p_sol_l = l.origin + tt * l.direction;
                let (u, v) = elslib_torus_parameters(
                    p_sol_l,
                    t.center,
                    x_ax,
                    y_ax,
                    t.axis,
                    t.major_radius,
                    t.minor_radius,
                );
                let p_sol_t = elslib_torus_value(
                    u,
                    v,
                    t.center,
                    t.axis,
                    t.major_radius,
                    t.minor_radius,
                );
                let dist = p_sol_t.distance_squared(p_sol_l);

                if dist > 0.0000000001 {
                    a_nb_bad_sol += 1;
                } else {
                    the_param.push(tt);
                    the_fi.push(u);
                    the_theta.push(v);
                    the_point.push(p_sol_l);
                }
            }
            let nbsolvalid = the_param.len();
            if n > 0 && nbsolvalid == 0 && a_nb_bad_sol == n {
                self.nbpt = 0;
                self.done = false;
            } else {
                self.nbpt = nbsolvalid;
                self.done = true;
                self.the_param = the_param;
                self.the_fi = the_fi;
                self.the_theta = the_theta;
                self.the_point = the_point;
            }
        } else {
            self.nbpt = 0;
            self.done = false;
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.nbpt
    }

    /// OCCT ParamOnLine(Index) — 1-based.
    pub fn param_on_line(&self, index: usize) -> f64 {
        self.the_param[index - 1]
    }

    /// OCCT ParamOnTorus(Index, U, V) — 1-based; returns (fi, theta).
    pub fn param_on_torus(&self, index: usize) -> (f64, f64) {
        (self.the_fi[index - 1], self.the_theta[index - 1])
    }

    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> DVec3 {
        self.the_point[index - 1]
    }
}

impl Default for IntLinTorus {
    fn default() -> Self {
        IntLinTorus::new()
    }
}
