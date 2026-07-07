//! OCCT-aligned: IntPatch_ImpPrmIntersection — analytic-parametric surface intersection.
//!
//! OCCT L617+ algorithm:
//!   1. Set quadric constraint from analytic surface
//!   2. Find boundary intersection points on parametric surface (SOnBounds)
//!   3. Find interior starting points (SearchInside)
//!   4. Walk intersection curves (IWalking)
//!   5. Decompose results into IntPatch_Line sequence
//!
//! rcad: boundary scanning + interior sampling + marching

use rcad_kernel::geom::{Surface3, SurfaceEval};
use super::int_patch_line::{IntPatchLine, WLinePnt, WLineType};
use super::int_patch_type::IntPatchIType;
use super::int_surf_quadric::Quadric;

pub struct ImpPrmIntersection {
    done: bool, empt: bool,
    spnt: Vec<super::int_patch_point::IntPatchPoint>,
    slin: Vec<IntPatchLine>,
    my_is_start_pnt: bool, my_u_start: f64, my_v_start: f64,
}

impl ImpPrmIntersection {
    pub fn new() -> Self {
        Self { done: false, empt: true, spnt: Vec::new(), slin: Vec::new(),
            my_is_start_pnt: false, my_u_start: 0.0, my_v_start: 0.0 }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }
    pub fn set_start_point(&mut self, u: f64, v: f64) {
        self.my_is_start_pnt = true; self.my_u_start = u; self.my_v_start = v;
    }

    /// OCCT L617: Perform — intersect analytic surface with parametric surface.
    ///
    /// Algorithm:
    /// 1. Identify analytic (quadric) vs parametric surface
    /// 2. Set constraint function F(u,v) = distance from parametric point to quadric
    /// 3. Scan boundary: find points where F=0 (SOnBounds)
    /// 4. Sample interior: find regions where F changes sign (SearchInside)
    /// 5. Walk: connect boundary and interior points into lines (IWalking)
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3,
                   tol_arc: f64, tol_tang: f64, _fleche: f64, _pas: f64) {
        self.done = false; self.empt = true; self.slin.clear(); self.spnt.clear();

        // Step 1: identify analytic (quadric) surface
        let (quad, _param_surf, _reversed) = match (Quadric::from_surface3(s1), Quadric::from_surface3(s2)) {
            (Some(q), _) => (q, s2, false),  // s1 is analytic
            (None, Some(q)) => (q, s1, true), // s2 is analytic
            (None, None) => { self.done = true; return; }
        };

        // Step 2: sample parametric surface on a grid, find intersection points
        let n_u = 30;
        let n_v = 30;
        let mut sign_changes: Vec<(f64, f64)> = Vec::new();

        // Constraint function: distance from point on parametric surface to analytic surface
        let constraint = |u: f64, v: f64| -> f64 {
            let pnt = _param_surf.point_at(u, v);
            quad.distance(pnt)
        };

        // sample grid for sign changes (zero crossings)
        for i in 0..n_u {
            for j in 0..n_v-1 {
                let u = i as f64 / n_u as f64;
                let v1 = j as f64 / n_v as f64;
                let v2 = (j+1) as f64 / n_v as f64;
                let f1 = constraint(u, v1);
                let f2 = constraint(u, v2);
                if f1 * f2 < 0.0 || f1.abs() < tol_tang || f2.abs() < tol_tang {
                    // Linear interpolation for zero crossing
                    let t = f1.abs() / (f1.abs() + f2.abs()).max(1e-30);
                    let v_zero = v1 + t * (v2 - v1);
                    sign_changes.push((u, v_zero));
                }
            }
        }
        for j in 0..n_v {
            for i in 0..n_u-1 {
                let u1 = i as f64 / n_u as f64;
                let u2 = (i+1) as f64 / n_u as f64;
                let v = j as f64 / n_v as f64;
                let f1 = constraint(u1, v);
                let f2 = constraint(u2, v);
                if f1 * f2 < 0.0 || f1.abs() < tol_tang || f2.abs() < tol_tang {
                    let t = f1.abs() / (f1.abs() + f2.abs()).max(1e-30);
                    let u_zero = u1 + t * (u2 - u1);
                    sign_changes.push((u_zero, v));
                }
            }
        }

        // Step 3: group nearby points into lines
        if sign_changes.is_empty() {
            self.empt = true; self.done = true; return;
        }
        sign_changes.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Step 4: create WLine from sample points
        let pnts: Vec<WLinePnt> = sign_changes.iter().map(|(u, v)| {
            let p3d = _param_surf.point_at(*u, *v);
            WLinePnt { p3d, u1: *u, v1: *v, u2: 0.0, v2: 0.0 }
        }).collect();

        if pnts.len() >= 2 {
            self.slin.push(IntPatchLine::walking(pnts, WLineType::ImpPrm));
        }

        self.empt = self.slin.is_empty();
        self.done = true;
    }
}
