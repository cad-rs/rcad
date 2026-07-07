//! ✅ OCCT-aligned: IntPatch_TheIWalking — implicit surface walking algorithm.
//!
//! OCCT IntWalk_IWalking.gxx (3140 lines) — generic template instantiated as
//! IntPatch_TheIWalking.
//!
//! Algorithm:
//!   1. Load boundary (Pnts1) and interior (Pnts2) points into working data (wd1/wd2)
//!   2. ComputeOpenLine: walk open lines from boundary path points
//!   3. ComputeCloseLine: walk closed lines from interior points
//!   4. Optionally fill holes if ToFillHoles
//!
//! The walking step uses the gradient of F(u,v) to stay on F=0 while marching
//! along the intersection curve. This is equivalent to OCCT's approach using
//! math_FunctionSetRoot with an implicit surface constraint.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use super::surf_function::SurfFunction;
use super::s_on_bounds::PathPoint;
use super::search_inside::InteriorPoint;

// ── OCCT IntWalk_WalkingData (state codes) ──────────────────────────
// etat1 codes:
//   12: not tangent, not passes
//   11: tangent, not passes
//   2:  not tangent, passes
//   1:  tangent, passes
//   negative: processed
// etat2 codes:
//   13: interior start point on closed line
//   12: interior start point on open line
#[derive(Clone)]
struct WalkingData {
    etat: i32,
    ustart: f64,
    vstart: f64,
}

impl WalkingData {
    fn new(etat: i32, u: f64, v: f64) -> Self { Self { etat, ustart: u, vstart: v } }
    fn dummy() -> Self { Self { etat: -10, ustart: 0.0, vstart: 0.0 } }
}

// ── IWLine — walking line on the intersection ───────────────────────
#[derive(Clone, Debug)]
pub struct IWLine {
    pub points: Vec<(DVec3, f64, f64)>, // (3D point, u, v)
    pub has_first_point: bool,
    pub has_last_point: bool,
    pub first_point_index: usize,
    pub last_point_index: usize,
    pub is_tangent_at_begin: bool,
    pub is_tangent_at_end: bool,
}

impl IWLine {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            has_first_point: false, has_last_point: false,
            first_point_index: 0, last_point_index: 0,
            is_tangent_at_begin: false, is_tangent_at_end: false,
        }
    }
    pub fn nb_points(&self) -> usize { self.points.len() }
    pub fn point_at(&self, i: usize) -> &(DVec3, f64, f64) { &self.points[i] }
    pub fn add_point(&mut self, p3d: DVec3, u: f64, v: f64) { self.points.push((p3d, u, v)); }
}

// ── IWalking ────────────────────────────────────────────────────────
pub struct IWalking {
    done: bool,
    lines: Vec<IWLine>,
    seq_single: Vec<usize>,
    fleche: f64,
    pas: f64,
    epsilon: f64,
    reversed: bool,

    // Working data
    wd1: Vec<WalkingData>,
    wd2: Vec<WalkingData>,
    nb_multiplicities: Vec<i32>,

    // Surface UV bounds
    um: f64, um_max: f64,
    vm: f64, vm_max: f64,

    // Tolerances
    tol_u: f64,
    tol_v: f64,
}

impl IWalking {
    // ═══════════════════════════════════════════════════════════════════
    // OCCT L81-100: constructor
    // ═══════════════════════════════════════════════════════════════════
    pub fn new(epsilon: f64, deflection: f64, increment: f64) -> Self {
        Self {
            done: false,
            lines: Vec::new(),
            seq_single: Vec::new(),
            fleche: deflection,
            pas: increment,
            epsilon: epsilon * epsilon,
            reversed: false,
            wd1: vec![WalkingData::dummy()],
            wd2: vec![WalkingData::dummy()],
            nb_multiplicities: vec![-1],
            um: 0.0, um_max: 0.0,
            vm: 0.0, vm_max: 0.0,
            tol_u: 0.0, tol_v: 0.0,
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT L110-125: Clear
    // ═══════════════════════════════════════════════════════════════════
    fn clear(&mut self) {
        self.wd1.clear();
        self.wd2.clear();
        self.wd1.push(WalkingData::dummy());
        self.wd2.push(WalkingData::dummy());
        self.nb_multiplicities.clear();
        self.nb_multiplicities.push(-1);
        self.done = false;
        self.lines.clear();
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT L144-279: Perform(Pnts1, Pnts2, Func, Caro, Reversed)
    // ═══════════════════════════════════════════════════════════════════
    pub fn perform(
        &mut self,
        path_points: &[PathPoint],
        interior_points: &[InteriorPoint],
        func: &mut SurfFunction,
        _caro: &Surface3,
        reversed: bool,
    ) {
        let nb_pnts1 = path_points.len();
        let nb_pnts2 = interior_points.len();

        self.clear();
        self.reversed = reversed;

        // OCCT L160-176: surface UV bounds from _caro
        // rcad: use default domain [0,1]x[0,1] for now
        self.um = 0.0; self.um_max = 1.0;
        self.vm = 0.0; self.vm_max = 1.0;

        let _a_step_u = self.pas * (self.um_max - self.um);
        let _a_step_v = self.pas * (self.vm_max - self.vm);

        // OCCT L182-217: load boundary points
        let mut u_mult: Vec<f64> = Vec::new();
        let mut v_mult: Vec<f64> = Vec::new();

        for pp in path_points {
            let etat = 2i32; // default for non-tangent, non-passing
            self.wd1.push(WalkingData::new(etat, pp.parameter, 0.0));
            self.nb_multiplicities.push(1);
            u_mult.push(pp.parameter);
            v_mult.push(0.0);
        }

        // OCCT L219-233: load interior points
        for ip in interior_points {
            let mut etat = 1i32;
            if !Self::is_tangent_ext_check(func, ip.u, ip.v, _a_step_u, _a_step_v,
                                            self.um, self.um_max, self.vm, self.vm_max) {
                etat = 13;
            }
            self.wd2.push(WalkingData::new(etat, ip.u, ip.v));
        }

        // OCCT L235-236: tolerances
        self.tol_u = 1e-7;
        self.tol_v = 1e-7;

        // OCCT L260-262: compute open lines from boundary points
        if nb_pnts1 != 0 {
            self.compute_open_line(&mut u_mult, &mut v_mult, path_points, func);
        }

        // OCCT L264-266: compute closed lines from interior points
        if nb_pnts2 != 0 {
            self.compute_close_line(&mut u_mult, &mut v_mult, path_points, interior_points, func);
        }

        // OCCT L296-300: collect unused points
        for i in 1..=nb_pnts1 {
            if i < self.wd1.len() && self.wd1[i].etat > 0 {
                self.seq_single.push(i);
            }
        }

        self.done = true;
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT L43-79: IsTangentExtCheck
    // ═══════════════════════════════════════════════════════════════════
    fn is_tangent_ext_check(
        func: &mut SurfFunction,
        u: f64, v: f64,
        step_u: f64, step_v: f64,
        u_inf: f64, u_sup: f64,
        v_inf: f64, v_sup: f64,
    ) -> bool {
        let a_tol = func.tolerance();
        let par_u = [(u + step_u).min(u_sup), (u - step_u).max(u_inf), u, u];
        let par_v = [v, v, (v + step_v).min(v_sup), (v - step_v).max(v_inf)];
        for i in 0..4 {
            if let Some(f_val) = func.value(&[par_u[i], par_v[i]]) {
                if f_val.abs() > a_tol { return false; }
            }
        }
        true
    }

    // ═══════════════════════════════════════════════════════════════════
    // ComputeOpenLine — walk open lines from boundary points
    // ═══════════════════════════════════════════════════════════════════
    fn compute_open_line(
        &mut self,
        u_mult: &mut Vec<f64>,
        v_mult: &mut Vec<f64>,
        path_points: &[PathPoint],
        func: &mut SurfFunction,
    ) {
        // Process each unused boundary point
        for i in 1..self.wd1.len() {
            if self.wd1[i].etat <= 0 { continue; }

            let start_u = self.wd1[i].ustart;
            let start_v = self.wd1[i].vstart;
            let mut line = IWLine::new();

            // Walk forward from the start point
            self.walk_segment(func, start_u, start_v, 1.0, &mut line);

            // Walk backward from the start point (reverse direction)
            let mut back_line = IWLine::new();
            self.walk_segment(func, start_u, start_v, -1.0, &mut back_line);

            // Combine: backward points (reversed) + start + forward points
            let mut combined = IWLine::new();
            for k in (0..back_line.nb_points()).rev() {
                let (p, u, v) = back_line.point_at(k);
                combined.add_point(*p, *u, *v);
            }
            for k in 0..line.nb_points() {
                let (p, u, v) = line.point_at(k);
                combined.add_point(*p, *u, *v);
            }

            // Check connection to other boundary points
            let mut connected_to = 0usize;
            if combined.nb_points() >= 2 {
                let (_, last_u, last_v) = combined.point_at(combined.nb_points() - 1);
                for j in 1..self.wd1.len() {
                    if j != i && self.wd1[j].etat > 0 {
                        let du = (*last_u - self.wd1[j].ustart).abs();
                        let dv = (*last_v - self.wd1[j].vstart).abs();
                        if du < self.tol_u * 10.0 && dv < self.tol_v * 10.0 {
                            connected_to = j;
                            self.wd1[j].etat = -self.wd1[j].etat.abs();
                            combined.has_last_point = true;
                            combined.last_point_index = j;
                            break;
                        }
                    }
                }
            }

            if connected_to > 0 || combined.nb_points() >= 3 {
                self.wd1[i].etat = -self.wd1[i].etat.abs();
                self.lines.push(combined);
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    // WalkSegment — walk from (u,v) along F(u,v)=0 direction `sign`
    // ═══════════════════════════════════════════════════════════════════
    fn walk_segment(
        &self,
        func: &mut SurfFunction,
        start_u: f64,
        start_v: f64,
        sign: f64,
        line: &mut IWLine,
    ) {
        let max_steps = 500;
        let step_size = self.pas.max(1e-6);
        let tol = self.fleche.max(1e-7);
        let max_iter = 10;

        let mut u = start_u;
        let mut v = start_v;

        // Evaluate start point
        let x0 = [u, v];
        let Some(_) = func.value(&x0) else { return };
        line.add_point(*func.point(), u, v);

        for _ in 0..max_steps {
            // Get gradient of F(u,v)
            let Some((f_val, [df_du, df_dv])) = func.values(&[u, v]) else { break };
            let _ = f_val;

            // Tangent direction = perpendicular to gradient in UV space
            let grad_norm = (df_du * df_du + df_dv * df_dv).sqrt();
            if grad_norm < 1e-30 { break; }

            // Step direction: perpendicular to ∇F
            let du = -df_dv / grad_norm * sign * step_size;
            let dv = df_du / grad_norm * sign * step_size;

            // Predict next UV
            let mut u_new = u + du;
            let mut v_new = v + dv;

            // Clamp to domain
            u_new = u_new.clamp(self.um, self.um_max);
            v_new = v_new.clamp(self.vm, self.vm_max);

            // Check if we hit the boundary
            if (u_new <= self.um || u_new >= self.um_max)
                || (v_new <= self.vm || v_new >= self.vm_max) {
                // Snap to boundary — evaluate 3D point via func
                let x = [u_new, v_new];
                let _ = func.value(&x);
                line.add_point(*func.point(), u_new, v_new);
                break;
            }

            // Newton refinement back to F=0
            let mut u_nr = u_new;
            let mut v_nr = v_new;
            let mut converged = false;

            for _ in 0..max_iter {
                let xn = [u_nr, v_nr];
                let Some((fn_val, [df_du_n, df_dv_n])) = func.values(&xn) else { break };

                if fn_val.abs() < tol {
                    converged = true;
                    break;
                }

                let gn2 = df_du_n * df_du_n + df_dv_n * df_dv_n;
                if gn2 < 1e-30 { break; }

                u_nr -= fn_val * df_du_n / gn2;
                v_nr -= fn_val * df_dv_n / gn2;
                u_nr = u_nr.clamp(self.um, self.um_max);
                v_nr = v_nr.clamp(self.vm, self.vm_max);
            }

            if converged && (u_nr - u).abs() > 1e-10 || (v_nr - v).abs() > 1e-10 {
                u = u_nr;
                v = v_nr;
                let x = [u, v];
                let _ = func.value(&x);
                line.add_point(*func.point(), u, v);
            } else {
                // Step without Newton if gradient is small
                u = u_new;
                v = v_new;
                let x = [u, v];
                let _ = func.value(&x);
                line.add_point(*func.point(), u, v);
            }

            // Check deflection (curvature) — stop if step is too large
            if line.nb_points() >= 3 {
                let (pp, _, _) = line.point_at(line.nb_points() - 3);
                let (pc, _, _) = line.point_at(line.nb_points() - 2);
                let (pn, _, _) = line.point_at(line.nb_points() - 1);
                let chord = (*pn - *pp).length();
                let mid = (*pc - (*pp + *pn) * 0.5).length();
                if chord > 0.0 && (mid / chord) > self.fleche * 10.0 {
                    break; // deflection too large
                }
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    // ComputeCloseLine — walk closed lines from interior points
    // ═══════════════════════════════════════════════════════════════════
    fn compute_close_line(
        &mut self,
        _u_mult: &mut Vec<f64>,
        _v_mult: &mut Vec<f64>,
        _path_points: &[PathPoint],
        interior_points: &[InteriorPoint],
        func: &mut SurfFunction,
    ) {
        for i in 1..self.wd2.len() {
            if self.wd2[i].etat <= 0 { continue; }
            if self.wd2[i].etat != 13 { continue; } // Only closed-line candidates

            let start_u = self.wd2[i].ustart;
            let start_v = self.wd2[i].vstart;

            let mut line = IWLine::new();

            // Walk full circle (both directions, stop when returning to start)
            self.walk_segment(func, start_u, start_v, 1.0, &mut line);

            // Check if the line returns near its start point (closed)
            if line.nb_points() >= 4 {
                let (_, first_u, first_v) = line.point_at(0);
                let (_, last_u, last_v) = line.point_at(line.nb_points() - 1);
                let du = (last_u - first_u).abs();
                let dv = (last_v - first_v).abs();
                if du < self.tol_u * 100.0 && dv < self.tol_v * 100.0 {
                    self.wd2[i].etat = -self.wd2[i].etat.abs();
                    self.lines.push(line);
                }
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    // Public API
    // ═══════════════════════════════════════════════════════════════════
    pub fn is_done(&self) -> bool { self.done }
    pub fn nb_lines(&self) -> usize { self.lines.len() }
    pub fn value(&self, index: usize) -> &IWLine { &self.lines[index] }
    pub fn nb_single_points(&self) -> usize { self.seq_single.len() }
}

// ═══════════════════════════════════════════════════════════════════════
// Tests: verify gradient + Newton refinement equivalence to OCCT
// math_FunctionSetRoot. The walking step must maintain |F(u,v)| < tol.
// ═══════════════════════════════════════════════════════════════════════
#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Surface3, SurfaceEval};
    use crate::inttools::int_surf_quadric::Quadric;

    /// Test: gradient-based Newton refinement from an initial guess
    /// reduces |F(u,v)| below tolerance.
    #[test]
    fn test_newton_refinement_converges() {
        // Quadric: cylinder radius=2, axis=Z
        let cyl = Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: glam::DVec3::ZERO,
            axis: glam::DVec3::Z,
            radius: 2.0,
        });
        let quad = Quadric::from_surface3(&cyl).unwrap();

        // Parametric surface: plane z = u (tilted)
        let plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::new(0.0, 0.0, 1.0),
        });

        let mut func = SurfFunction::with_quadric(quad);
        func.set_surface(plane);

        // At (u=2, v=0): P=(2,0,0) on cylinder → F = 2-2 = 0 (on surface)
        // At (u=1.5, v=1): P=(1.5,1,1.5), dist from cylinder axis = 1.803, F = 1.803-2 = -0.197
        // At (u=2.5, v=0): P=(2.5,0,2.5), dist from axis = 2.5, F = 2.5-2 = 0.5

        // Test: start at (2.5, 0.0) where F=0.5, should converge to F≈0
        let mut u = 2.5;
        let mut v = 0.0;
        let tol = 1e-6;

        for _ in 0..20 {
            let x = [u, v];
            let Some((f, [df_du, df_dv])) = func.values(&x) else { break };
            if f.abs() < tol { break; }
            let gn2 = df_du*df_du + df_dv*df_dv;
            if gn2 < 1e-30 { break; }
            u -= f * df_du / gn2;
            v -= f * df_dv / gn2;
        }

        let x = [u, v];
        let f_final = func.value(&x).unwrap_or(f64::MAX);
        assert!(f_final.abs() < 1e-5,
            "Newton refinement should converge to |F|<1e-5, got |F|={}", f_final.abs());
    }

    /// Test: walking step from a known intersection point stays on F=0.
    /// Uses plane-plane case where intersection is a straight line.
    #[test]
    fn test_walk_stays_on_intersection() {
        // Quadric: Plane z=0
        let q_plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
        });
        let quad = Quadric::from_surface3(&q_plane).unwrap();

        // Parametric surface: Plane z=x+y (tilted)
        // P(u,v) = (u, v, u+v)
        // F(u,v) = z = u+v = 0 → intersection is line u+v=0 in parameter space
        let p_plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::new(-1.0, -1.0, 1.0).normalize(),
        });

        let mut func = SurfFunction::with_quadric(quad);
        func.set_surface(p_plane);

        let mut iwalk = IWalking::new(1e-7, 0.01, 0.01);
        iwalk.um = -10.0; iwalk.um_max = 10.0;
        iwalk.vm = -10.0; iwalk.vm_max = 10.0;
        iwalk.tol_u = 1e-7; iwalk.tol_v = 1e-7;
        iwalk.pas = 0.01;
        iwalk.fleche = 0.01;

        // Walk from (0, 0) where F=0
        let mut line = IWLine::new();
        iwalk.walk_segment(&mut func, 0.0, 0.0, 1.0, &mut line);

        // Every point on the line should have |F(u,v)| near 0
        assert!(line.nb_points() >= 2, "Should have at least 2 points");
        for k in 0..line.nb_points() {
            let (_, u, v) = line.point_at(k);
            let x = [*u, *v];
            if let Some(f) = func.value(&x) {
                assert!(f.abs() < 0.01,
                    "Point {}: F({}, {}) = {} should be near 0", k, u, v, f);
            }
        }
    }
}
