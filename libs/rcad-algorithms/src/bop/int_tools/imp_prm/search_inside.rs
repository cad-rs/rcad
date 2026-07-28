//! IntPatch_TheSearchInside — interior starting point search.
//!
//! OCCT IntStart_SearchInside.gxx (~300 lines) + .lxx.
//!
//! Finds points inside the parametric surface domain where F(u,v)=0.
//! Uses random/structured sampling + Gauss-Newton refinement via SurfFunction.
//! These interior points serve as starting points for closed intersection curves.

use super::surf_function::SurfFunction;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::Surface3;

// ── OCCT IntSurf_InteriorPoint ───────────────────────────────────────
#[derive(Clone, Debug)]
pub struct InteriorPoint {
    pub value: DVec3,
    pub u: f64,
    pub v: f64,
    pub direction: DVec3,
    pub direction_2d: DVec2,
}

impl InteriorPoint {
    pub fn new(value: DVec3, u: f64, v: f64, dir: DVec3, dir2d: DVec2) -> Self {
        Self {
            value,
            u,
            v,
            direction: dir,
            direction_2d: dir2d,
        }
    }
}

// ── OCCT IntPatch_TheSearchInside ────────────────────────────────────
pub struct SearchInside {
    done: bool,
    list: Vec<InteriorPoint>,
}

impl SearchInside {
    // OCCT L34: default constructor
    pub fn new() -> Self {
        Self {
            done: false,
            list: Vec::new(),
        }
    }

    // OCCT L41-44: Perform(F, Surf, T, Epsilon)
    // Search for interior points by sampling the parametric surface domain.
    // rcad: grid sampling + Newton-Raphson refinement (no TopolTool).
    pub fn perform(
        &mut self,
        func: &mut SurfFunction,
        u_min: f64,
        u_max: f64,
        v_min: f64,
        v_max: f64,
        epsilon: f64,
    ) {
        self.list.clear();
        self.done = false;

        // OCCT: samples based on TopolTool (NbSamples). rcad: 8x8 grid.
        let n_u = 8;
        let n_v = 8;

        for i in 0..=n_u {
            for j in 0..=n_v {
                let u = u_min + (i as f64 / n_u as f64) * (u_max - u_min);
                let v = v_min + (j as f64 / n_v as f64) * (v_max - v_min);
                let x = [u, v];

                let Some(f) = func.value(&x) else { continue };

                // Near zero → refine directly
                if f.abs() < epsilon * 10.0 {
                    self.refine_and_add(func, u, v, epsilon);
                    continue;
                }

                // Sign change in U direction → zero crossing
                if i < n_u {
                    let u2 = u_min + ((i + 1) as f64 / n_u as f64) * (u_max - u_min);
                    let x2 = [u2, v];
                    if let Some(f2) = func.value(&x2) {
                        if f * f2 < 0.0 {
                            let t = f.abs() / (f.abs() + f2.abs());
                            let u_zero = u + t * (u2 - u);
                            self.refine_and_add(func, u_zero, v, epsilon);
                        }
                    }
                }

                // Sign change in V direction → zero crossing
                if j < n_v {
                    let v2 = v_min + ((j + 1) as f64 / n_v as f64) * (v_max - v_min);
                    let x2 = [u, v2];
                    if let Some(f2) = func.value(&x2) {
                        if f * f2 < 0.0 {
                            let t = f.abs() / (f.abs() + f2.abs());
                            let v_zero = v + t * (v2 - v);
                            self.refine_and_add(func, u, v_zero, epsilon);
                        }
                    }
                }
            }
        }

        self.done = true;
    }

    // OCCT L46-49: Perform(F, Surf, UStart, VStart) — from a given start point
    pub fn perform_from_point(&mut self, func: &mut SurfFunction, u_start: f64, v_start: f64) {
        self.list.clear();
        self.done = false;

        self.refine_and_add(func, u_start, v_start, 1e-7);

        self.done = true;
    }

    // ── Gauss-Newton refinement (identical to OCCT math_FunctionSetRoot) ──
    fn refine_and_add(&mut self, func: &mut SurfFunction, u: f64, v: f64, eps: f64) {
        let tol = eps.max(1e-8);
        if let Some((un, vn)) = super::i_walking::IWalking::gauss_newton_root(func, u, v, tol) {
            if !self.is_duplicate(un, vn, tol) {
                let dir_3d = func.direction_3d();
                let dir_2d = func.direction_2d();
                self.list
                    .push(InteriorPoint::new(*func.point(), un, vn, dir_3d, dir_2d));
            }
        }
    }

    fn is_duplicate(&self, u: f64, v: f64, tol: f64) -> bool {
        self.list
            .iter()
            .any(|p| (p.u - u).abs() < tol && (p.v - v).abs() < tol)
    }

    // ── Public API ───────────────────────────────────────────────────
    pub fn is_done(&self) -> bool {
        self.done
    }
    pub fn nb_points(&self) -> usize {
        self.list.len()
    }
    pub fn value(&self, index: usize) -> &InteriorPoint {
        &self.list[index]
    }
}
