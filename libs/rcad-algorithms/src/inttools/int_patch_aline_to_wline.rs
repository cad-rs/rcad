//! IntPatch_ALineToWLine — convert analytic line to walking point sequence.
//!
//! OCCT IntPatch_ALineToWLine.hxx / .cxx (33K)
//!
//! Converts IntPatch_ALine (analytic intersection curves like circles/ellipses)
//! into IntPatch_WLine (point sequences) for marching refinement.
//! Used when the analytic intersection needs numerical refinement.

use rcad_kernel::geom::{Curve3, Surface3};
use rcad_kernel::CurveEval;
use super::int_patch_line::{IntPatchLine, WLinePnt, WLineType};

/// ALineToWLine — converts analytic lines to walking point sequences.
pub struct ALineToWLine {
    nb_points: usize,
    tol_open_domain: f64,
    tol_transition: f64,
    tol_3d: f64,
}

impl ALineToWLine {
    pub fn new(_s1: &Surface3, _s2: &Surface3, nb_points: usize) -> Self {
        Self { nb_points, tol_open_domain: 1e-5, tol_transition: 1e-5, tol_3d: 1e-7 }
    }

    pub fn set_tol_open_domain(&mut self, t: f64) { self.tol_open_domain = t; }
    pub fn tol_open_domain(&self) -> f64 { self.tol_open_domain }
    pub fn set_tol_transition(&mut self, t: f64) { self.tol_transition = t; }
    pub fn tol_transition(&self) -> f64 { self.tol_transition }
    pub fn set_tol_3d(&mut self, t: f64) { self.tol_3d = t; }
    pub fn tol_3d(&self) -> f64 { self.tol_3d }

    /// OCCT L53: MakeWLine — convert an analytic IntPatchLine to a walking point sequence.
    /// Returns a new IntPatchLine with wline_pnts filled, or None if conversion fails.
    pub fn make_wline(&self, line: &IntPatchLine) -> Option<IntPatchLine> {
        if line.is_wline() { return Some(line.clone()); }
        let n = self.nb_points.max(10);
        let t0 = line.t_range[0];
        let t1 = line.t_range[1];
        let dt = (t1 - t0) / (n as f64);
        if !dt.is_finite() || dt.abs() > 1e10 { return None; }

        let mut pnts = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = t0 + i as f64 * dt;
            let p3d = line.curve.point_at(t);
            pnts.push(WLinePnt { p3d, u1: 0.0, v1: 0.0, u2: 0.0, v2: 0.0 });
        }
        Some(IntPatchLine::walking(pnts, WLineType::ImpImp))
    }
}
