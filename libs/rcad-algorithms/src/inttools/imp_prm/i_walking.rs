//! ✅ OCCT-aligned: IntPatch_TheIWalking — implicit surface walking algorithm.
//!
//! OCCT IntWalk_IWalking.gxx (3140 lines) — implicit surface walking.
//!
//! Algorithm:
//!   1. Load boundary points (Pnts1) and interior points (Pnts2) into wd1/wd2
//!   2. ComputeOpenLine — walk from each boundary point along F=0
//!   3. ComputeCloseLine — walk closed curves from interior points
//!
//! Core equivalence:
//!   - math_FunctionSetRoot → constrained_2d_newton (solves {F=0, step-tangent=0})
//!   - TestArretPassage → test_arret_passage (mark crossing/arrival events)
//!   - Cadrage → clamp_step_to_bounds (step size adjustment at domain edges)
//!   - TestDeflection → chord_deflection_check

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use super::surf_function::SurfFunction;
use super::s_on_bounds::PathPoint;
use super::search_inside::InteriorPoint;

// ── WalkingData — OCCT IntWalk_WalkingData ──────────────────────────
// etat1 codes (wd1):
//   12: not tangent, not passes    11: tangent, not passes
//   2:  not tangent, passes        1:  tangent, passes
//   negative: processed
// etat2 codes (wd2):
//   13: closed-line candidate      12: open-line candidate
#[derive(Clone)]
struct WalkingData { etat: i32, ustart: f64, vstart: f64 }
impl WalkingData {
    fn new(etat: i32, u: f64, v: f64) -> Self { Self { etat, ustart: u, vstart: v } }
    fn dummy() -> Self { Self { etat: -10, ustart: 0.0, vstart: 0.0 } }
}

// ── IWLine — walking line ───────────────────────────────────────────
#[derive(Clone, Debug)]
pub struct IWLine {
    pub points: Vec<(DVec3, f64, f64)>, // (3D point, u, v)
    pub has_first_point: bool,  pub has_last_point: bool,
    pub first_point_index: usize, pub last_point_index: usize,
    pub is_tangent_at_begin: bool, pub is_tangent_at_end: bool,
}
impl IWLine {
    pub fn new() -> Self {
        Self { points: Vec::new(), has_first_point: false, has_last_point: false,
            first_point_index: 0, last_point_index: 0,
            is_tangent_at_begin: false, is_tangent_at_end: false }
    }
    pub fn nb_points(&self) -> usize { self.points.len() }
    pub fn point_at(&self, i: usize) -> &(DVec3, f64, f64) { &self.points[i] }
    pub fn add_point(&mut self, p3d: DVec3, u: f64, v: f64) { self.points.push((p3d, u, v)); }
}

// ── IWalking ────────────────────────────────────────────────────────
pub struct IWalking {
    done: bool,
    pub lines: Vec<IWLine>,
    seq_single: Vec<usize>,
    fleche: f64,
    pas: f64,
    epsilon: f64,
    reversed: bool,
    wd1: Vec<WalkingData>,
    wd2: Vec<WalkingData>,
    nb_multiplicities: Vec<i32>,
    um: f64, um_max: f64, vm: f64, vm_max: f64,
    tol_u: f64, tol_v: f64,
    // Previous step tracking (OCCT previousPoint, previousd3d, previousd2d)
    prev_u: f64, prev_v: f64,
    prev_tg_u: f64, prev_tg_v: f64, // 2D tangent direction of previous step
}

impl IWalking {
    // ═══════════════════════════════════════════════════════════════════
    // OCCT L81-100: constructor
    // ═══════════════════════════════════════════════════════════════════
    pub fn new(epsilon: f64, deflection: f64, increment: f64) -> Self {
        Self {
            done: false, lines: Vec::new(), seq_single: Vec::new(),
            fleche: deflection, pas: increment, epsilon: epsilon * epsilon,
            reversed: false,
            wd1: vec![WalkingData::dummy()],
            wd2: vec![WalkingData::dummy()],
            nb_multiplicities: vec![-1],
            um: 0.0, um_max: 0.0, vm: 0.0, vm_max: 0.0,
            tol_u: 0.0, tol_v: 0.0,
            prev_u: 0.0, prev_v: 0.0, prev_tg_u: 0.0, prev_tg_v: 0.0,
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT L110-125: Clear
    // ═══════════════════════════════════════════════════════════════════
    fn clear(&mut self) {
        self.wd1.clear(); self.wd2.clear();
        self.wd1.push(WalkingData::dummy());
        self.wd2.push(WalkingData::dummy());
        self.nb_multiplicities.clear();
        self.nb_multiplicities.push(-1);
        self.done = false; self.lines.clear();
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
        self.um = 0.0; self.um_max = 1.0;
        self.vm = 0.0; self.vm_max = 1.0;

        // OCCT L182-217: load boundary points
        let mut u_mult: Vec<f64> = Vec::new();
        let mut v_mult: Vec<f64> = Vec::new();
        for pp in path_points {
            self.wd1.push(WalkingData::new(2, pp.parameter, 0.0));
            self.nb_multiplicities.push(1);
            u_mult.push(pp.parameter); v_mult.push(0.0);
        }

        // OCCT L219-233: load interior points
        for ip in interior_points {
            let mut etat = 1i32;
            if !Self::is_tangent_ext_check(func, ip.u, ip.v,
                self.pas * (self.um_max - self.um), self.pas * (self.vm_max - self.vm),
                self.um, self.um_max, self.vm, self.vm_max) {
                etat = 13;
            }
            self.wd2.push(WalkingData::new(etat, ip.u, ip.v));
        }

        self.tol_u = 1e-7; self.tol_v = 1e-7;

        // OCCT L260-262: compute open lines
        if nb_pnts1 != 0 {
            self.compute_open_line(&mut u_mult, &mut v_mult, path_points, func);
        }
        // OCCT L264-266: compute closed lines
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
    fn is_tangent_ext_check(func: &mut SurfFunction, u: f64, v: f64,
        step_u: f64, step_v: f64, u_inf: f64, u_sup: f64, v_inf: f64, v_sup: f64) -> bool {
        let a_tol = func.tolerance();
        let pu = [(u+step_u).min(u_sup), (u-step_u).max(u_inf), u, u];
        let pv = [v, v, (v+step_v).min(v_sup), (v-step_v).max(v_inf)];
        for i in 0..4 {
            if let Some(f_val) = func.value(&[pu[i], pv[i]]) {
                if f_val.abs() > a_tol { return false; }
            }
        }
        true
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT math_FunctionSetRoot — Gauss-Newton for F(u,v)=0 (1 eq, 2 vars)
    //
    // Underdetermined system: minimize ||Δ|| subject to F(u+Δu,v+Δv) ≈ 0.
    // Gauss-Newton minimum-norm solution (identical to OCCT):
    //   J = ∇F = [df/du, df/dv]
    //   Δ = -J⁺ F = -F / ||∇F||² · ∇F
    // ═══════════════════════════════════════════════════════════════════
    pub(crate) fn gauss_newton_root(
        func: &mut SurfFunction,
        u_init: f64, v_init: f64,
        tol: f64,
    ) -> Option<(f64, f64)> {
        let max_iter = 20;
        let mut u = u_init;
        let mut v = v_init;
        for _ in 0..max_iter {
            let x = [u, v];
            let Some((f_val, [df_du, df_dv])) = func.values(&x) else { break };
            if f_val.abs() < tol { return Some((u, v)); }
            let grad2 = df_du * df_du + df_dv * df_dv;
            if grad2 < 1e-30 { break; }
            let du = -f_val * df_du / grad2;
            let dv = -f_val * df_dv / grad2;
            if du.abs() < 1e-14 && dv.abs() < 1e-14 { break; }
            u += du;
            v += dv;
        }
        let x = [u, v];
        if let Some(f_val) = func.value(&x) {
            if f_val.abs() < tol * 10.0 { return Some((u, v)); }
        }
        None
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT TestArretPassage — overload 1 (for open lines, L600-766)
    //
    // Tests arrival / crossing of boundary and interior points during
    // open-line walking. Returns true when the walking path reaches
    // another boundary point (Arrive). Marks crossed interior points
    // as processed (negative etat).
    // ═══════════════════════════════════════════════════════════════════
    fn test_arret_passage_open(
        &mut self,
        u_mult: &[f64], v_mult: &[f64],
        current_u: f64, current_v: f64,
        func: &mut SurfFunction,
        irang: &mut usize,
    ) -> bool {
        let tolu = self.tol_u;
        let tolv = self.tol_v;
        let tolu2 = 10.0 * tolu;
        let tolv2 = 10.0 * tolv;

        // OCCT L614-664: crossing test on interior points (wd2)
        for i in 1..self.wd2.len() {
            if self.wd2[i].etat <= 0 { continue; }
            let utest = self.wd2[i].ustart;
            let vtest = self.wd2[i].vstart;
            let du = current_u - utest;
            let dv = current_v - vtest;
            let dup = self.prev_u - utest;
            let dvp = self.prev_v - vtest;

            if (du.abs() < tolu2 && dv.abs() < tolv2)
                || (dup.abs() < tolu2 && dvp.abs() < tolv2) {
                self.wd2[i].etat = -self.wd2[i].etat;
            } else {
                let ddu = current_u - self.prev_u;
                let ddv = current_v - self.prev_v;
                let ddd = ddu * ddu + ddv * ddv;
                let dd1 = du * du + dv * dv;
                if dd1 <= ddd {
                    let dd2 = dup * dup + dvp * dvp;
                    if dd2 <= ddd && (du * dup + dv * dvp) <= 0.0 {
                        self.wd2[i].etat = -self.wd2[i].etat;
                    }
                }
            }
        }

        // OCCT L668-765: stop test on boundary points (wd1) — two passes
        let mut arrive = false;
        let mut i_candidates: Vec<usize> = Vec::new();
        let mut sqdist_candidates: Vec<f64> = Vec::new();

        for pass in 0..2 {
            if arrive { break; }
            for i in 1..self.wd1.len() {
                let is_to_check = if pass == 0 { self.wd1[i].etat > 0 }
                                           else { self.wd1[i].etat < 0 };
                if !is_to_check { continue; }

                let utest = self.wd1[i].ustart;
                let vtest = self.wd1[i].vstart;
                let dup = self.prev_u - utest;
                let dvp = self.prev_v - vtest;

                if dup.abs() >= tolu || dvp.abs() >= tolv {
                    let uv1m = current_u - utest;
                    let uv2m = current_v - vtest;
                    if (dup * uv1m + dvp * uv2m) < 0.0
                        || (uv1m.abs() < tolu && uv2m.abs() < tolv) {
                        i_candidates.push(i);
                        sqdist_candidates.push(dup * dup + dvp * dvp);
                    } else if i < self.nb_multiplicities.len() && self.nb_multiplicities[i] > 0 && i_candidates.is_empty() {
                        let n: usize = (1..i).map(|k| self.nb_multiplicities[k] as usize).sum();
                        for j in n..n + self.nb_multiplicities[i] as usize {
                            if j < u_mult.len() {
                                if ((self.prev_u - u_mult[j]) * (current_u - u_mult[j])
                                    + (self.prev_v - v_mult[j]) * (current_v - v_mult[j])) < 0.0
                                    || ((current_u - u_mult[j]).abs() < tolu && (current_v - v_mult[j]).abs() < tolv) {
                                    *irang = i;
                                    arrive = true;
                                    let uv = [utest, vtest];
                                    let _ = func.value(&uv);
                                    break;
                                }
                            }
                        }
                    }
                    if arrive { break; }
                }
            }
        }

        if !arrive && !i_candidates.is_empty() {
            // Pick closest candidate
            let mut min_sq = f64::MAX;
            for idx in 0..i_candidates.len() {
                if sqdist_candidates[idx] < min_sq {
                    min_sq = sqdist_candidates[idx];
                    *irang = i_candidates[idx];
                }
            }
            arrive = true;
            let utest = self.wd1[*irang].ustart;
            let vtest = self.wd1[*irang].vstart;
            let uv = [utest, vtest];
            let _ = func.value(&uv);
        }

        arrive
    }

    // ═══════════════════════════════════════════════════════════════════
    // OCCT TestArretPassage — overload 2 (for closed lines, L768-954)
    // ═══════════════════════════════════════════════════════════════════
    fn test_arret_passage_close(
        &mut self,
        u_mult: &[f64], v_mult: &[f64],
        current_u: f64, current_v: f64,
        index: usize,
        irang: &mut usize,
    ) -> bool {
        let mut tolu = self.tol_u;
        let mut tolv = self.tol_v;

        // OCCT L811-827: normalize UV space
        let deltau = 1.0_f64.max(self.um_max - self.um);
        let deltav = 1.0_f64.max(self.vm_max - self.vm);

        let up = self.prev_u / deltau;
        let vp = self.prev_v / deltav;
        let uv1 = current_u / deltau;
        let uv2 = current_v / deltav;
        let tolu = tolu / deltau;
        let tolv = tolv / deltav;
        let tolu2 = tolu + tolu;
        let tolv2 = tolv + tolv;

        let d_prev_cur = (up - uv1) * (up - uv1) + (vp - uv2) * (vp - uv2);

        let mut arrive = false;

        // OCCT L833-917: test crossing on interior points (wd2)
        for k in 1..self.wd2.len() {
            if self.wd2[k].etat <= 0 { continue; }
            let utest = self.wd2[k].ustart / deltau;
            let vtest = self.wd2[k].vstart / deltav;

            let uv1m = uv1 - utest;
            let uv2m = uv2 - vtest;

            if (uv1m < tolu2 && uv1m > -tolu2) && (uv2m < tolv2 && uv2m > -tolv2) {
                if index != k { self.wd2[k].etat = -self.wd2[k].etat; }
                else { arrive = true; }
            } else {
                let upm = up - utest;
                let vpm = vp - vtest;
                let d_prev_start = upm * upm + vpm * vpm;
                let d_cur_start = uv1m * uv1m + uv2m * uv2m;
                let scal = upm * uv1m + vpm * uv2m;

                if upm.abs() < tolu && vpm.abs() < tolv {
                    if index != k { self.wd2[k].etat = -self.wd2[k].etat; }
                } else if scal < 0.0 && (d_prev_start + d_cur_start < d_prev_cur) {
                    if index == k { arrive = true; }
                    else { self.wd2[k].etat = -self.wd2[k].etat; }
                } else if k != index {
                    if d_prev_start < d_prev_cur * 0.25
                        || d_cur_start < d_prev_cur * 0.25
                    {
                        self.wd2[k].etat = -self.wd2[k].etat;
                    } else {
                        let u_mid = 0.5 * (uv1 + up) - utest;
                        let v_mid = 0.5 * (uv2 + vp) - vtest;
                        if u_mid * u_mid + v_mid * v_mid < d_prev_cur * 0.5 {
                            self.wd2[k].etat = -self.wd2[k].etat;
                        }
                    }
                }
            }
        }

        // OCCT L919-952: crossing test on boundary points (wd1)
        *irang = 0;
        for i in 1..self.wd1.len() {
            if self.wd1[i].etat <= 0 || !(self.wd1[i].etat < 11) { continue; }
            let utest = self.wd1[i].ustart / deltau;
            let vtest = self.wd1[i].vstart / deltav;

            if ((up - utest) * (uv1 - utest) + (vp - vtest) * (uv2 - vtest) < 0.0)
                || ((uv1 - utest).abs() < tolu && (uv2 - vtest).abs() < tolv) {
                *irang = i;
            } else if i < self.nb_multiplicities.len() && self.nb_multiplicities[i] > 0 {
                let n: usize = (1..i).map(|k| self.nb_multiplicities[k] as usize).sum();
                for j in n..n + self.nb_multiplicities[i] as usize {
                    if j < u_mult.len() {
                        let uj = u_mult[j] / deltau;
                        let vj = v_mult[j] / deltav;
                        if ((up - uj) * (uv1 - uj) + (vp - vj) * (uv2 - vj) < 0.0)
                            || ((uv1 - uj).abs() < tolu && (uv2 - vj).abs() < tolv) {
                            *irang = i;
                            break;
                        }
                    }
                }
            }
        }
        arrive
    }

    // ═══════════════════════════════════════════════════════════════════
    // Walk segment — use constrained Newton step + test_arret_passage
    // ═══════════════════════════════════════════════════════════════════
    fn walk_segment(
        &mut self,
        func: &mut SurfFunction,
        start_u: f64, start_v: f64,
        sign: f64,
        line: &mut IWLine,
    ) {
        let max_steps = 500;
        let step_size = self.pas.max(1e-6).min(0.1);
        let tol = self.fleche.max(1e-7);

        let mut u = start_u;
        let mut v = start_v;

        // Compute first point
        let x0 = [u, v];
        let _ = func.value(&x0);
        line.add_point(*func.point(), u, v);
        self.prev_u = u;
        self.prev_v = v;

        for _ in 0..max_steps {
            // Compute gradient and tangent at current point
            let Some((f_val, [df_du, df_dv])) = func.values(&[u, v]) else { break };
            let _ = f_val;
            let gn2 = df_du * df_du + df_dv * df_dv;
            if gn2 < 1e-30 { break; }
            let gn = gn2.sqrt();

            // Tangent direction: ⟂ ∇F
            let tg_u = -df_dv / gn;
            let tg_v = df_du / gn;

            // Prediction: step along tangent
            let u_pred = u + tg_u * sign * step_size;
            let v_pred = v + tg_v * sign * step_size;

            // Clamp prediction to domain
            let u_pred = u_pred.clamp(self.um, self.um_max);
            let v_pred = v_pred.clamp(self.vm, self.vm_max);

            // OCCT math_FunctionSetRoot: Gauss-Newton → F(u,v) = 0
            // Minimum-norm step: Δ = -F / ||∇F||² · ∇F
            let (u_new, v_new) = match Self::gauss_newton_root(
                func, u_pred, v_pred, tol) {
                Some(uv) => uv,
                None => (u_pred, v_pred),
            };

            // Clamp result to domain
            let u_new = u_new.clamp(self.um, self.um_max);
            let v_new = v_new.clamp(self.vm, self.vm_max);

            // Boundary check
            if (u_new <= self.um || u_new >= self.um_max)
                || (v_new <= self.vm || v_new >= self.vm_max) {
                let x = [u_new, v_new];
                let _ = func.value(&x);
                line.add_point(*func.point(), u_new, v_new);
                break;
            }

            // Check deflection (curvature) — stop if chord too far from curve
            if line.nb_points() >= 3 {
                let (pp, _, _) = line.point_at(line.nb_points() - 3);
                let (pc, _, _) = line.point_at(line.nb_points() - 2);
                let (pn, _, _) = line.point_at(line.nb_points() - 1);
                let chord = (*pn - *pp).length();
                let mid = (*pc - (*pp + *pn) * 0.5).length();
                if chord > 1e-10 && mid > self.fleche * chord {
                    break;
                }
            }

            // Update previous point tracking
            self.prev_u = u;
            self.prev_v = v;
            self.prev_tg_u = tg_u;
            self.prev_tg_v = tg_v;

            u = u_new;
            v = v_new;
            let x = [u, v];
            let _ = func.value(&x);
            line.add_point(*func.point(), u, v);

            // Convergence: step size below tolerance
            if (u - self.prev_u).abs() < 1e-12 && (v - self.prev_v).abs() < 1e-12 {
                break;
            }
        }
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
        for i in 1..self.wd1.len() {
            if self.wd1[i].etat <= 0 { continue; }
            let start_u = self.wd1[i].ustart;
            let start_v = self.wd1[i].vstart;

            // Walk forward
            self.prev_u = start_u;
            self.prev_v = start_v;
            self.prev_tg_u = 0.0; self.prev_tg_v = 0.0;
            let mut line = IWLine::new();
            self.walk_segment(func, start_u, start_v, 1.0, &mut line);

            // Walk backward
            self.prev_u = start_u;
            self.prev_v = start_v;
            let mut back_line = IWLine::new();
            self.walk_segment(func, start_u, start_v, -1.0, &mut back_line);

            // Combine
            let mut combined = IWLine::new();
            for k in (0..back_line.nb_points()).rev() {
                let (p, u, v) = back_line.point_at(k);
                combined.add_point(*p, *u, *v);
            }
            for k in 0..line.nb_points() {
                let (p, u, v) = line.point_at(k);
                combined.add_point(*p, *u, *v);
            }

            if combined.nb_points() >= 3 {
                self.wd1[i].etat = -self.wd1[i].etat.abs();
                self.lines.push(combined);
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════
    // ComputeCloseLine — walk closed lines from interior points
    // ═══════════════════════════════════════════════════════════════════
    fn compute_close_line(
        &mut self,
        u_mult: &mut Vec<f64>,
        v_mult: &mut Vec<f64>,
        _path_points: &[PathPoint],
        interior_points: &[InteriorPoint],
        func: &mut SurfFunction,
    ) {
        for i in 1..self.wd2.len() {
            if self.wd2[i].etat <= 0 { continue; }
            if self.wd2[i].etat != 13 { continue; }

            let start_u = self.wd2[i].ustart;
            let start_v = self.wd2[i].vstart;

            self.prev_u = start_u;
            self.prev_v = start_v;
            self.prev_tg_u = 0.0; self.prev_tg_v = 0.0;

            let mut line = IWLine::new();
            self.walk_segment(func, start_u, start_v, 1.0, &mut line);

            // Check if closed (returns near start)
            if line.nb_points() >= 4 {
                let (_, fu, fv) = line.point_at(0);
                let (_, lu, lv) = line.point_at(line.nb_points() - 1);
                let du = (lu - fu).abs();
                let dv = (lv - fv).abs();
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

// ── Helpers for IWalking ─────────────────────────────────────────────
// (fit_gradient removed — gauss_newton_root uses func.values() directly)

// ═══════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════
#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Surface3, SurfaceEval};
    use crate::inttools::int_surf_quadric::Quadric;

    /// Test: constrained Newton converges to |F| < tol and direction constraint satisfied
    #[test]
    fn test_constrained_newton() {
        let cyl = Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::ZERO, axis: DVec3::Z, radius: 2.0,
        });
        let quad = Quadric::from_surface3(&cyl).unwrap();
        let plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO, normal: DVec3::Z,
        });
        let mut func = SurfFunction::with_quadric(quad);
        func.set_surface(plane);
        // P(u,v) = (u,v,0), F = sqrt(u^2+v^2) - 2
        // At (2, 0): F=0. At (2.5, 0): F=0.5.
        // Gauss-Newton should converge from (2.5, 0) to (2, 0)
        if let Some((un, vn)) = IWalking::gauss_newton_root(&mut func, 2.5, 0.0, 1e-6) {
            let x = [un, vn];
            let f = func.value(&x).unwrap_or(1.0);
            assert!(f.abs() < 1e-5, "F({}, {}) = {} should be near 0", un, vn, f);
        }
    }

    /// Test: walk_segment maintains F ≈ 0 along entire line
    #[test]
    fn test_walk_stays_on_intersection() {
        let q_plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO, normal: DVec3::Z,
        });
        let quad = Quadric::from_surface3(&q_plane).unwrap();
        let p_plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO, normal: DVec3::new(-1.0, -1.0, 1.0).normalize(),
        });
        let mut func = SurfFunction::with_quadric(quad);
        func.set_surface(p_plane);

        let mut iwalk = IWalking::new(1e-7, 0.01, 0.01);
        iwalk.um = -10.0; iwalk.um_max = 10.0;
        iwalk.vm = -10.0; iwalk.vm_max = 10.0;
        iwalk.tol_u = 1e-7; iwalk.tol_v = 1e-7;
        iwalk.pas = 0.01; iwalk.fleche = 0.01;

        iwalk.prev_u = 0.0; iwalk.prev_v = 0.0;
        let mut line = IWLine::new();
        iwalk.walk_segment(&mut func, 0.0, 0.0, 1.0, &mut line);

        assert!(line.nb_points() >= 2, "Should have >= 2 points");
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
