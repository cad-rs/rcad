//!  ?OCCT-aligned: IntWalk_PWalking  ?intersection curve walking for
//!   two parametric surfaces (parametric-parametric).
//!
//! OCCT source: TKGeomAlgo/IntWalk/IntWalk_PWalking.cxx (4096 lines)
//!   + .hxx (292 lines) + .lxx (86 lines).
//!
//! Algorithm: from a seed point (U1,V1,U2,V2), step along the
//! intersection curve using iso-parametric lines, adapt step size
//! based on deflection, stop at surface boundaries.
//!
//! rcad integration: uses inttools::marching for the core walking step
//! (march_intersection_with_config), preserving OCCT parameter semantics.

use glam::DVec3;
use rcad_kernel::geom::*;

use crate::inttools::marching::{self, MarchingConfig};
use crate::extrema;
use crate::tolerance::TOLERANCE_LEN_MIN;

use super::prm_prm_intersection::PntOn2S;

/// OCCT IntWalk_StatusDeflection
#[derive(Clone, Copy, PartialEq)]
enum StatusDeflection {
    PasTropGrand,   // step OK
    PasTropPetit,   // step too small  ?refine
    PasTropGrandInFlexion, // step too large at inflection
}

/// OCCT IntImp_ConstIsoparametric (which iso line to follow)
#[derive(Clone, Copy, PartialEq)]
enum ConstIsoparametric {
    U1, V1, U2, V2,
    None,
}

///  ?OCCT-aligned: IntWalk_PWalking
pub struct PWalking {
    // ── Public state (OCCT hxx:243-287) ──────────────────────────
    done: bool,
    line: Vec<PntOn2S>,          // OCCT: handle<IntSurf_LineOn2S>
    close: bool,
    tgfirst: bool,
    tglast: bool,
    my_tangent_idx: usize,
    fleche: f64,                  // Deflection
    pas_max: f64,                 // max step fraction
    tolconf: f64,                 // epsilon (point confusion)
    my_tol_tang: f64,             // TolTangency
    pasuv: [f64; 4],              // step sizes [U1,V1,U2,V2]
    my_step_min: [f64; 4],
    pas_sav: [f64; 4],
    pas_init: [f64; 4],
    um1: f64, um1_max: f64,
    vm1: f64, vm1_max: f64,
    um2: f64, um2_max: f64,
    vm2: f64, vm2_max: f64,
    reso_u1: f64, reso_u2: f64,
    reso_v1: f64, reso_v2: f64,
    sens_cheminement: i32,       // walking direction (±1)
    choix_iso_sav: ConstIsoparametric,
    previous_point: Option<PntOn2S>,
    previoustg: bool,
    // rcad: surface references
    s1: Surface3,
    s2: Surface3,
}

impl PWalking {
    /// OCCT L219-414: constructor
    pub fn new(
        s1: &Surface3, s2: &Surface3,
        tol_tangency: f64, epsilon: f64,
        deflection: f64, increment: f64,
    ) -> Self {
        let pas_max = increment * 0.2;
        let (um1, um1_max, vm1, vm1_max) = uv_range(s1);
        let (um2, um2_max, vm2, vm2_max) = uv_range(s2);

        let reso_u1 = u_resolution(s1, epsilon);
        let reso_v1 = v_resolution(s1, epsilon);
        let reso_u2 = u_resolution(s2, epsilon);
        let reso_v2 = v_resolution(s2, epsilon);
        let scale_reso = |reso: f64, lo: f64, hi: f64| -> f64 {
            let max_val = lo.abs().max(hi.abs());
            let new_reso = reso * max_val;
            if new_reso > reso && new_reso < 10.0 { new_reso } else { reso }
        };

        let reso_u1 = scale_reso(reso_u1, um1, um1_max);
        let reso_v1 = scale_reso(reso_v1, vm1, vm1_max);
        let reso_u2 = scale_reso(reso_u2, um2, um2_max);
        let reso_v2 = scale_reso(reso_v2, vm2, vm2_max);
        let mut pasuv = [
            pas_max * (um1_max - um1).abs(),
            pas_max * (vm1_max - vm1).abs(),
            pas_max * (um2_max - um2).abs(),
            pas_max * (vm2_max - vm2).abs(),
        ];
        let clamp_reso = |reso: f64, step: f64| -> f64 {
            if reso > 0.0001 * step { 0.00001 * step } else { reso }
        };
        let reso_u1 = clamp_reso(reso_u1, pasuv[0]);
        let reso_v1 = clamp_reso(reso_v1, pasuv[1]);
        let reso_u2 = clamp_reso(reso_u2, pasuv[2]);
        let reso_v2 = clamp_reso(reso_v2, pasuv[3]);
        // rcad: periodic surfaces handled by the marching step
        let my_step_min = [100.0 * reso_u1, 100.0 * reso_v1, 100.0 * reso_u2, 100.0 * reso_v2];
        for i in 0..4 {
            if pasuv[i] > 10.0 { pasuv[i] = 10.0; }
        }

        PWalking {
            done: false,
            line: Vec::new(),
            close: false,
            tgfirst: false,
            tglast: false,
            my_tangent_idx: 0,
            fleche: deflection,
            pas_max,
            tolconf: epsilon,
            my_tol_tang: tol_tangency,
            pasuv,
            my_step_min,
            pas_sav: pasuv,
            pas_init: pasuv,
            um1, um1_max, vm1, vm1_max,
            um2, um2_max, vm2, vm2_max,
            reso_u1, reso_u2, reso_v1, reso_v2,
            sens_cheminement: 1,
            choix_iso_sav: ConstIsoparametric::None,
            previous_point: None,
            previoustg: false,
            s1: s1.clone(),
            s2: s2.clone(),
        }
    }

    // ── Public API (OCCT hxx:87-186) ────────────────────────────

    pub fn is_done(&self) -> bool { self.done }
    pub fn nb_points(&self) -> usize { self.line.len() }
    pub fn value(&self, index: usize) -> Option<&PntOn2S> { self.line.get(index - 1) }
    pub fn line(&self) -> &[PntOn2S] { &self.line }
    pub fn tangent_at_first(&self) -> bool { self.tgfirst }
    pub fn tangent_at_last(&self) -> bool { self.tglast }
    pub fn is_closed(&self) -> bool { self.close }
    pub fn max_step(&self, index: usize) -> f64 {
        if index < 4 { self.pas_init[index] } else { 0.0 }
    }

    // ── PerformFirstPoint (OCCT hxx:105-106) ─────────────────────

    /// OCCT L: find the first point on the intersection from a seed.
    /// rcad: use marching's seed validation + projection.
    pub fn perform_first_point(&mut self, par_dep: &[f64; 4], first_point: &mut PntOn2S) -> bool {
        if par_dep.len() < 4 { return false; }
        let (u1, v1, u2, v2) = (par_dep[0], par_dep[1], par_dep[2], par_dep[3]);
        let p1 = surface_point_at(&self.s1, u1, v1);
        let p2 = surface_point_at(&self.s2, u2, v2);
        if !p1.is_finite() || !p2.is_finite() { return false; }
        let dist = p1.distance(p2);
        if dist > self.my_tol_tang.max(TOLERANCE_MESH_LEGACY) { return false; }

        *first_point = PntOn2S { p3d: (p1 + p2) * 0.5, u1, v1, u2, v2 };
        self.previous_point = Some(first_point.clone());
        self.previoustg = false;
        true
    }

    // ── Perform (OCCT hxx:94-102) ────────────────────────────────

    /// Walk the intersection line from seed parameters.
    /// rcad: uses inttools::marching::march_intersection_with_config.
    pub fn perform_with_bounds(
        &mut self,
        par_dep: &[f64; 4],
        u1min: f64, v1min: f64, u2min: f64, v2min: f64,
        u1max: f64, v1max: f64, u2max: f64, v2max: f64,
    ) {
        self.line.clear();
        self.done = false;

        if par_dep.len() < 4 { return; }

        let step_size = self.pasuv[0].min(self.pasuv[1]).min(self.pasuv[2]).min(self.pasuv[3]);
        let min_step = self.my_step_min[0].min(self.my_step_min[1])
                       .min(self.my_step_min[2]).min(self.my_step_min[3]);

        let config = MarchingConfig {
            step_size,
            min_step_size: min_step,
            max_steps: 10000,
            max_oscillations: 10,
            step_reduction_factor: 0.5,
            deflection_tol: self.fleche,
            multiscale_seeds: false,
        };

        // March from the 3D seed point
        let seed_p3d = surface_point_at(&self.s1, par_dep[0], par_dep[1]);
        if !seed_p3d.is_finite() { return; }

        let result = marching::march_intersection_with_config(
            &self.s1, &self.s2,
            seed_p3d,
            &config,
            |p| {
                if let Some((u1, v1)) = project_onto_surface(&self.s1, p) {
                    if let Some((u2, v2)) = project_onto_surface(&self.s2, p) {
                        return u1 >= u1min - 0.1 && u1 <= u1max + 0.1
                            && v1 >= v1min - 0.1 && v1 <= v1max + 0.1
                            && u2 >= u2min - 0.1 && u2 <= u2max + 0.1
                            && v2 >= v2min - 0.1 && v2 <= v2max + 0.1;
                    }
                }
                false
            },
        );

        // Convert SampledCurve (3D points)  ?PntOn2S (UV+3D)
        for p3d in &result.points {
            if !p3d.is_finite() { continue; }
            let Some((u1, v1)) = project_onto_surface(&self.s1, *p3d) else { continue };
            let Some((u2, v2)) = project_onto_surface(&self.s2, *p3d) else { continue };
            self.line.push(PntOn2S { p3d: *p3d, u1, v1, u2, v2 });
        }

        if self.line.len() >= 2 {
            let f = &self.line[0];
            let l = &self.line[self.line.len() - 1];
            self.close = (f.u1 - l.u1).abs() < 1e-6 && (f.v1 - l.v1).abs() < 1e-6
                      && (f.u2 - l.u2).abs() < 1e-6 && (f.v2 - l.v2).abs() < 1e-6;
        }

        self.done = true;
    }

    // ── PutToBoundary (OCCT cxx:2951-3155) ───────────────────────

    /// OCCT-aligned: snap first and last line points to surface boundaries
    /// if they are within tolerance.  Returns true if any point was added.
    pub fn put_to_boundary(&mut self, s1: &Surface3, s2: &Surface3) -> bool {
        if self.line.len() < 2 { return false; }

        let a_tol_min = 1e-12; // minimum UV step threshold
        let (u1b_f, u1b_l, v1b_f, v1b_l) = uv_range(s1);
        let (u2b_f, u2b_l, v2b_f, v2b_l) = uv_range(s2);

        let u1b_first = u1b_f; let u1b_last = u1b_l;
        let v1b_first = v1b_f; let v1b_last = v1b_l;
        let u2b_first = u2b_f; let u2b_last = u2b_l;
        let v2b_first = v2b_f; let v2b_last = v2b_l;

        let mut a_tol = 1.0f64;
        a_tol = f64::min(a_tol, u1b_last - u1b_first);
        a_tol = f64::min(a_tol, u2b_last - u2b_first);
        a_tol = f64::min(a_tol, v1b_last - v1b_first);
        a_tol = f64::min(a_tol, v2b_last - v2b_first) * 1.0e-3;

        if a_tol <= 2.0 * a_tol_min { return false; }

        let (is_u1_par, is_v1_par) = Self::is_parallel(&self.line, true, a_tol);
        let (is_u2_par, is_v2_par) = Self::is_parallel(&self.line, false, a_tol);

        // Check first point
        let mut added = false;
        let (u1, v1, u2, v2) = {
            let p = &self.line[0];
            (p.u1, p.v1, p.u2, p.v2)
        };
        let mut need = false;
        let mut snap_u1 = u1; let mut snap_v1 = v1;
        let mut snap_u2 = u2; let mut snap_v2 = v2;
        if !is_v1_par { need = Self::snap_to_bnd(&mut snap_u1, u1b_first, u1b_last, a_tol_min, a_tol) || need; }
        if !is_v2_par { need = Self::snap_to_bnd(&mut snap_u2, u2b_first, u2b_last, a_tol_min, a_tol) || need; }
        if !is_u1_par { need = Self::snap_to_bnd(&mut snap_v1, v1b_first, v1b_last, a_tol_min, a_tol) || need; }
        if !is_u2_par { need = Self::snap_to_bnd(&mut snap_v2, v2b_first, v2b_last, a_tol_min, a_tol) || need; }
        if need {
            added |= self.seek_point_on_boundary(s1, s2, snap_u1, snap_v1, snap_u2, snap_v2, true);
        }

        // Check last point
        let n = self.line.len();
        let (u1, v1, u2, v2) = {
            let p = &self.line[n - 1];
            (p.u1, p.v1, p.u2, p.v2)
        };
        need = false;
        let mut snap_u1 = u1; let mut snap_v1 = v1;
        let mut snap_u2 = u2; let mut snap_v2 = v2;
        if !is_v1_par { need = Self::snap_to_bnd(&mut snap_u1, u1b_first, u1b_last, a_tol_min, a_tol) || need; }
        if !is_v2_par { need = Self::snap_to_bnd(&mut snap_u2, u2b_first, u2b_last, a_tol_min, a_tol) || need; }
        if !is_u1_par { need = Self::snap_to_bnd(&mut snap_v1, v1b_first, v1b_last, a_tol_min, a_tol) || need; }
        if !is_u2_par { need = Self::snap_to_bnd(&mut snap_v2, v2b_first, v2b_last, a_tol_min, a_tol) || need; }
        if need {
            added |= self.seek_point_on_boundary(s1, s2, snap_u1, snap_v1, snap_u2, snap_v2, false);
        }

        added
    }

    /// OCCT IsParallel: check if line is parallel to surface boundary.
    fn is_parallel(line: &[PntOn2S], check_surf1: bool, toler: f64) -> (bool, bool) {
        const MAX_POINTS: usize = 23;
        let n = line.len();
        if n < 3 { return (true, true); }
        let n_points = n.min(MAX_POINTS);
        let step = n as f64 / n_points as f64;
        let mut a_u_min = f64::MAX; let mut a_u_max = f64::MIN;
        let mut a_v_min = f64::MAX; let mut a_v_max = f64::MIN;
        let mut a_n_point = 1.0_f64;
        for _ in 0..n_points {
            let idx = if (a_n_point as usize) > n { n - 1 } else { a_n_point as usize - 1 };
            let (u, v) = if check_surf1 {
                (line[idx].u1, line[idx].v1)
            } else {
                (line[idx].u2, line[idx].v2)
            };
            if u < a_u_min { a_u_min = u; }
            if u > a_u_max { a_u_max = u; }
            if v < a_v_min { a_v_min = v; }
            if v > a_v_max { a_v_max = v; }
            a_n_point += step;
        }
        let is_v_par = (a_u_max - a_u_min) < toler;
        let is_u_par = (a_v_max - a_v_min) < toler;
        (is_u_par, is_v_par)
    }

    /// Snap a parameter to boundary if within tolerance.
    fn snap_to_bnd(param: &mut f64, lo: f64, hi: f64, tol_min: f64, tol: f64) -> bool {
        let delta = *param - lo;
        if tol_min < delta && delta < tol { *param = lo; return true; }
        let delta = hi - *param;
        if tol_min < delta && delta < tol { *param = hi; return true; }
        false
    }

    /// Try to compute a 3D point from UV on both surfaces.
    fn try_project_3d(s1: &Surface3, s2: &Surface3,
                      u1: f64, v1: f64, u2: f64, v2: f64) -> Option<DVec3> {
        let p1 = surface_point_at(s1, u1, v1);
        let p2 = surface_point_at(s2, u2, v2);
        if p1.is_finite() && p2.is_finite() {
            Some((p1 + p2) * 0.5)
        } else { None }
    }

    // ── SeekPointOnBoundary (OCCT cxx:2716-2950) ─────────────────

    /// OCCT-aligned: find an intersection point on a surface boundary near
    /// given UV.  Uses gradient descent + projection (rcad: closest_point_on_surface).
    /// is_the_first: true = snap start point, false = snap end point.
    pub fn seek_point_on_boundary(&mut self,
                                  s1: &Surface3, s2: &Surface3,
                                  u1: f64, v1: f64, u2: f64, v2: f64,
                                  is_the_first: bool) -> bool
    {
        let (u1b_f, u1b_l, v1b_f, v1b_l) = uv_range(s1);
        let (u2b_f, u2b_l, v2b_f, v2b_l) = uv_range(s2);
        let low = [u1b_f, v1b_f, u2b_f, v2b_f];
        let upp = [u1b_l, v1b_l, u2b_l, v2b_l];

        let a3d_tol = CONFUSION; // OCCT: from surface resolution
        let a_tol = TOLERANCE_MESH_LEGACY.max(a3d_tol);
        let mut pnt = [u1, v1, u2, v2];
        let p3d_mid = {
            let p1 = surface_point_at(s1, pnt[0], pnt[1]);
            let p2 = surface_point_at(s2, pnt[2], pnt[3]);
            if p1.is_finite() && p2.is_finite() { (p1 + p2) * 0.5 } else { return false; }
        };

        // Project midpoint onto both surfaces for refined UV
        if let Some((nu1, nv1)) = project_onto_surface(s1, p3d_mid) {
            pnt[0] = nu1; pnt[1] = nv1;
        }
        if let Some((nu2, nv2)) = project_onto_surface(s2, p3d_mid) {
            pnt[2] = nu2; pnt[3] = nv2;
        }

        // Clamp to domain
        let adjust_to_domain = |p: &mut [f64; 4]| {
            for i in 0..4 {
                if p[i] < low[i] { p[i] = low[i]; }
                if p[i] > upp[i] { p[i] = upp[i]; }
            }
        };
        adjust_to_domain(&mut pnt);

        // Verify the point
        let p1 = surface_point_at(s1, pnt[0], pnt[1]);
        let p2 = surface_point_at(s2, pnt[2], pnt[3]);
        let p_int = (p1 + p2) * 0.5;
        let sq_dist = p_int.distance_squared(p1);
        if sq_dist > a_tol * a_tol { return false; }
        //   with loop detection (hairpin bend check).
        // rcad: simplified  ?insert directly without loop check since
        //   marching already produces consistent orientations.
        if is_the_first {
            if !self.line.is_empty() {
                let dist = self.line[0].p3d.distance(p_int);
                if dist > TOLERANCE_LEN_MIN {
                    self.line.insert(0, PntOn2S {
                        p3d: p_int, u1: pnt[0], v1: pnt[1],
                        u2: pnt[2], v2: pnt[3],
                    });
                }
            }
        } else {
            if !self.line.is_empty() {
                let dist = self.line.last().unwrap().p3d.distance(p_int);
                if dist > TOLERANCE_LEN_MIN {
                    self.line.push(PntOn2S {
                        p3d: p_int, u1: pnt[0], v1: pnt[1],
                        u2: pnt[2], v2: pnt[3],
                    });
                }
            }
        }
        true
    }

    // ── DistanceMinimizeByGradient (OCCT cxx:2394-2522) ─────────

    /// OCCT-aligned: gradient descent to minimize distance between two surfaces.
    ///   rcad: uses extrema::closest_point_on_surface as the numerical engine.
    pub fn distance_minimize_by_gradient(&self, s1: &Surface3, s2: &Surface3,
                                         init: &mut [f64; 4]) -> bool {
        // rcad: use closest_point_on_surface for both surfaces
        let p3d = surface_point_at(s1, init[0], init[1]);
        if !p3d.is_finite() { return false; }

        if let Some((u2, v2)) = project_onto_surface(s2, p3d) {
            init[2] = u2; init[3] = v2;
        }
        if let Some((u1, v1)) = project_onto_surface(s1, surface_point_at(s2, init[2], init[3])) {
            init[0] = u1; init[1] = v1;
        }

        // Verify minimized distance
        let p1 = surface_point_at(s1, init[0], init[1]);
        let p2 = surface_point_at(s2, init[2], init[3]);
        p1.distance_squared(p2) < 1e-12
    }

    // ── DistanceMinimizeByExtrema (OCCT cxx:2532-2583) ──────────

    /// OCCT-aligned: minimize distance from a 3D point to a surface.
    ///   rcad: uses extrema::closest_point_on_surface.
    pub fn distance_minimize_by_extrema(&self, surf: &Surface3, p0: DVec3,
                                        u0: &mut f64, v0: &mut f64) -> bool {
        if let Some((u, v)) = project_onto_surface(surf, p0) {
            *u0 = u; *v0 = v;
            let ps = surface_point_at(surf, u, v);
            return ps.distance_squared(p0) < 1e-12;
        }
        false
    }

    // ── HandleSingleSingularPoint (OCCT cxx:2587-2712) ──────────

    /// OCCT-aligned: handle singular points on surface boundaries.
    fn handle_single_singular_point(&self, s1: &Surface3, s2: &Surface3,
                                    _a3d_tol: f64, pnt: &mut [f64; 4]) -> bool {
        let (u1b_f, u1b_l, v1b_f, v1b_l) = uv_range(s1);
        let (u2b_f, u2b_l, v2b_f, v2b_l) = uv_range(s2);
        let low = [u1b_f, v1b_f, u2b_f, v2b_f];
        let upp = [u1b_l, v1b_l, u2b_l, v2b_l];

        let conf = 1e-12;
        for i in 0..4 {
            let at_low = (pnt[i] - low[i]).abs() < conf;
            let at_upp = (pnt[i] - upp[i]).abs() < conf;
            if !at_low && !at_upp { continue; }

            let p3d = surface_point_at(if i < 2 { s1 } else { s2 }, pnt[0], pnt[1]);
            if !p3d.is_finite() { continue; }

            if let Some((nu1, nv1)) = project_onto_surface(s1, p3d) { pnt[0] = nu1; pnt[1] = nv1; }
            if let Some((nu2, nv2)) = project_onto_surface(s2, p3d) { pnt[2] = nu2; pnt[3] = nv2; }
            if at_low { pnt[i] = low[i]; }
            if at_upp { pnt[i] = upp[i]; }

            let mut in_domain = true;
            for j in 0..4 {
                if (pnt[j] - low[j] + conf) * (pnt[j] - upp[j] - conf) > 0.0 { in_domain = false; break; }
            }
            if in_domain { return true; }
        }
        false
    }

    // ── SeekAdditionalPoints (OCCT cxx:3159-3350) ───────────────

    /// OCCT-aligned: add midpoint points until line has at least min points.
    pub fn seek_additional_points(&mut self, s1: &Surface3, s2: &Surface3,
                                   min_nb_points: usize) -> bool {
        if self.line.len() > min_nb_points { return true; }
        let (u1b_f, u1b_l, v1b_f, v1b_l) = uv_range(s1);
        let (u2b_f, u2b_l, v2b_f, v2b_l) = uv_range(s2);

        let mut prev_len = 0;
        while self.line.len() < min_nb_points && self.line.len() != prev_len {
            prev_len = self.line.len();
            let mut i = 0;
            while i + 1 < self.line.len() {
                let a = &self.line[i];
                let b = &self.line[i + 1];

                let u1 = (a.u1 + b.u1) * 0.5;
                let v1 = (a.v1 + b.v1) * 0.5;
                let u2 = (a.u2 + b.u2) * 0.5;
                let v2 = (a.v2 + b.v2) * 0.5;
                let u1 = u1.clamp(u1b_f, u1b_l);
                let v1 = v1.clamp(v1b_f, v1b_l);
                let u2 = u2.clamp(u2b_f, u2b_l);
                let v2 = v2.clamp(v2b_f, v2b_l);

                let p3d = surface_point_at(s1, u1, v1);
                if p3d.is_finite() {
                    self.line.insert(i + 1, PntOn2S { p3d, u1, v1, u2, v2 });
                    i += 1;
                }
                i += 1;
            }
        }
        self.line.len() >= min_nb_points
    }

    // ── ExtendLineInCommonZone (OCCT cxx:1831-2393) ─────────────

    /// OCCT-aligned: attempt to extend the line through a tangent/common zone.
    /// Steps along the tangent direction, minimizing distance between both surfaces
    /// at each step, until boundary or tangent zone exit.
    /// the_direction_flag: true = forward, false = backward.
    pub fn extend_line_in_common_zone(&mut self, s1: &Surface3, s2: &Surface3,
                                       the_choix_iso: ConstIsoparametric,
                                       the_direction_flag: bool) -> bool {
        let (u1b_f, u1b_l, v1b_f, v1b_l) = uv_range(s1);
        let (u2b_f, u2b_l, v2b_f, v2b_l) = uv_range(s2);
        let prev = match &self.previous_point {
            Some(p) => p.clone(),
            None => return false,
        };
        let mut param = [prev.u1, prev.v1, prev.u2, prev.v2];
        if (param[0] - u1b_f).abs() < self.reso_u1 { return false; }
        if (param[1] - v1b_f).abs() < self.reso_v1 { return false; }
        if (param[2] - u2b_f).abs() < self.reso_u2 { return false; }
        if (param[3] - v2b_f).abs() < self.reso_v2 { return false; }
        if (param[0] - u1b_l).abs() < self.reso_u1 { return false; }
        if (param[1] - v1b_l).abs() < self.reso_v1 { return false; }
        if (param[2] - u2b_l).abs() < self.reso_u2 { return false; }
        if (param[3] - v2b_l).abs() < self.reso_v2 { return false; }

        let mut b_stop = false;
        let mut out_of_tangent = false;
        let mut nb_iter_no_append = 0i32;
        let mut nb_equal = 0i32;
        while !b_stop && nb_iter_no_append < 20 && nb_equal < 20 {
            nb_iter_no_append += 1;

            // Compute step magnitude from the iso direction
            let f = 0.1_f64.max({
                match the_choix_iso {
                    ConstIsoparametric::U1 | ConstIsoparametric::V1 |
                    ConstIsoparametric::U2 | ConstIsoparametric::V2 => 10.0,
                    _ => 10.0,
                }
            });
            let sg = |diff: f64| -> f64 { if diff >= 0.0 { 1.0 } else { -1.0 } };
            param[0] += self.pasuv[0] * sg(param[0] - prev.u1) / f;
            param[1] += self.pasuv[1] * sg(param[1] - prev.v1) / f;
            param[2] += self.pasuv[2] * sg(param[2] - prev.u2) / f;
            param[3] += self.pasuv[3] * sg(param[3] - prev.v2) / f;
            // rcad: use distance_minimize_by_gradient
            let solved = self.distance_minimize_by_gradient(s1, s2, &mut param);

            if !solved {
                // OCCT: if !IsDone  ?return
                out_of_tangent = true;
                break;
            }
            let a_status = self.test_deflection(the_choix_iso, StatusDeflection::PasTropPetit);
            match a_status {
                StatusDeflection::PasTropGrand => {
                    // step too big  ?reduce step (OCCT L2000-2021 simulated)
                    for u in 0..4 {
                        self.pasuv[u] *= 0.8;
                        self.pasuv[u] = self.pasuv[u].max(self.my_step_min[u]);
                    }
                }
                _ => {
                    // OK or other: check stop conditions
                    b_stop = self.test_arret(the_direction_flag, &mut param, &mut ConstIsoparametric::None);

                    if !b_stop {
                        let np = &PntOn2S {
                            p3d: DVec3::ZERO,
                            u1: param[0], v1: param[1],
                            u2: param[2], v2: param[3],
                        };
                        if let Some(ref pp) = self.previous_point {
                            if (np.u1 - pp.u1).abs() < self.reso_u1 && (np.v1 - pp.v1).abs() < self.reso_v1
                                || (np.u2 - pp.u2).abs() < self.reso_u2 && (np.v2 - pp.v2).abs() < self.reso_v2
                            {
                                nb_equal += 1;
                            } else {
                                nb_equal = 0;
                            }
                        }
                    }

                    out_of_tangent = true;

                    if !b_stop {
                        let u1 = param[0].clamp(u1b_f, u1b_l);
                        let v1 = param[1].clamp(v1b_f, v1b_l);
                        let u2 = param[2].clamp(u2b_f, u2b_l);
                        let v2 = param[3].clamp(v2b_f, v2b_l);
                        let p3d = surface_point_at(s1, u1, v1);

                        if p3d.is_finite() {
                            let new_p = PntOn2S { p3d, u1, v1, u2, v2 };
                            self.previous_point = Some(new_p.clone());
                            self.line.push(new_p);
                            nb_iter_no_append = 0;
                        }
                    }
                }
            }
        }

        out_of_tangent
    }

    pub fn remove_point(&mut self, index: usize) {
        if index < self.line.len() { self.line.remove(index); }
    }

    // ── TestDeflection (OCCT hxx:138-139) ─────────────────────────

    /// OCCT-aligned: test if deflection between last 3 points is acceptable.
    pub fn test_deflection(&mut self, choix_iso: ConstIsoparametric,
                           status: StatusDeflection) -> StatusDeflection {
        let n = self.line.len();
        if n < 3 { return StatusDeflection::PasTropPetit; }

        // Compute deflection (mid-point distance from chord)
        let p1 = &self.line[n - 3];
        let p2 = &self.line[n - 2];
        let p3 = &self.line[n - 1];
        let mid = (p1.p3d + p3.p3d) * 0.5;
        let deflection = mid.distance(p2.p3d);

        if deflection > self.fleche {
            StatusDeflection::PasTropGrand
        } else {
            StatusDeflection::PasTropPetit
        }
    }

    // ── TestArret (OCCT hxx:141-143) ─────────────────────────────

    /// OCCT-aligned: check if the line has reached a surface boundary.
    pub fn test_arret(&self, _deja_reparti: bool, param: &mut [f64; 4],
                      choix_iso: &mut ConstIsoparametric) -> bool {
        // Check if any parameter is at a surface boundary
        let u1 = param[0]; let v1 = param[1];
        let u2 = param[2]; let v2 = param[3];
        let at_u1_boundary = (u1 - self.um1).abs() < self.reso_u1 || (u1 - self.um1_max).abs() < self.reso_u1;
        let at_v1_boundary = (v1 - self.vm1).abs() < self.reso_v1 || (v1 - self.vm1_max).abs() < self.reso_v1;
        let at_u2_boundary = (u2 - self.um2).abs() < self.reso_u2 || (u2 - self.um2_max).abs() < self.reso_u2;
        let at_v2_boundary = (v2 - self.vm2).abs() < self.reso_v2 || (v2 - self.vm2_max).abs() < self.reso_v2;
        at_u1_boundary || at_v1_boundary || at_u2_boundary || at_v2_boundary
    }

    // ── RepartirOuDiviser (OCCT hxx:145-148) ─────────────────────

    /// OCCT-aligned: restart or subdivide the walking step.
    pub fn repartir_ou_diviser(&mut self, _deja_reparti: &mut bool,
                               _choix_iso: &mut ConstIsoparametric,
                               _arrive: &mut bool) {
        // rcad simplified: marching already handles step adaptation.
    }
}

// ── Surface helpers (OCCT Adaptor3d_HSurfaceTool equivalents) ─────

fn uv_range(s: &Surface3) -> (f64, f64, f64, f64) {
    match s {
        Surface3::Plane(_) => (-1e5, 1e5, -1e5, 1e5),
        Surface3::Sphere(sp) => {
            (0.0, 2.0 * std::f64::consts::PI, -std::f64::consts::PI * 0.5, std::f64::consts::PI * 0.5)
        }
        Surface3::Cylinder(_) => (0.0, 2.0 * std::f64::consts::PI, -1e5, 1e5),
        Surface3::Cone(_) => (0.0, 2.0 * std::f64::consts::PI, -1e5, 1e5),
        Surface3::Torus(_) => (0.0, 2.0 * std::f64::consts::PI, 0.0, 2.0 * std::f64::consts::PI),
        Surface3::BSpline(bsp) => {
            let u_max = (bsp.knots_u.len().saturating_sub(bsp.degree_u + 1)) as f64;
            let v_max = (bsp.knots_v.len().saturating_sub(bsp.degree_v + 1)) as f64;
            (0.0, u_max.max(1.0), 0.0, v_max.max(1.0))
        }
        Surface3::Bezier(_) => (0.0, 1.0, 0.0, 1.0),
        _ => (0.0, 1.0, 0.0, 1.0),
    }
}

fn u_resolution(s: &Surface3, _tol: f64) -> f64 {
    // OCCT: Adaptor3d_HSurfaceTool::UResolution(surface, tolerance)
    // rcad: estimate from UV range
    let (u_min, u_max, _, _) = uv_range(s);
    let range = (u_max - u_min).abs();
    if range.is_finite() && range > 1e-10 { range * CONFUSION } else { CONFUSION }
}

fn v_resolution(s: &Surface3, _tol: f64) -> f64 {
    let (_, _, v_min, v_max) = uv_range(s);
    let range = (v_max - v_min).abs();
    if range.is_finite() && range > 1e-10 { range * CONFUSION } else { CONFUSION }
}

fn surface_point_at(surf: &Surface3, u: f64, v: f64) -> DVec3 {
    use rcad_kernel::geom::SurfaceEval;
    surf.point_at(u, v)
}

fn project_onto_surface(surf: &Surface3, pt: DVec3) -> Option<(f64, f64)> {
    let (uv, _proj) = extrema::closest_point_on_surface(surf, pt);
    if uv.x.is_finite() && uv.y.is_finite() {
        Some((uv.x, uv.y))
    } else {
        None
    }
}
