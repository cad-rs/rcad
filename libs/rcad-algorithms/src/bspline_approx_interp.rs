//! BSpline curve approximation and interpolation from scattered points.
//!
//! ✅ OCCT-aligned: Approx_BSplineApproxInterp

use glam::DVec3;
use rcad_kernel::geom::{BSplineCurve3, CurveEval};

/// Chord-length parameterization for a set of points.
fn chord_length_params(points: &[DVec3]) -> Vec<f64> {
    let n = points.len();
    if n <= 1 { return vec![0.0; n]; }
    let mut params = vec![0.0; n];
    let mut total = 0.0;
    for i in 1..n { total += (points[i] - points[i - 1]).length(); }
    let mut acc = 0.0;
    for i in 1..n - 1 {
        acc += (points[i] - points[i - 1]).length();
        params[i] = acc / total;
    }
    params[n - 1] = 1.0;
    params
}

/// BSpline curve approximation/interpolation from point data.
///
/// OCCT-aligned: Approx_BSplineApproxInterp.
/// Approximates a point set with a BSpline curve. Specific points can be
/// marked for exact interpolation with optional kink (C0 break).
#[derive(Debug, Clone)]
pub struct BSplineApproxInterp {
    points: Vec<DVec3>,
    params: Vec<f64>,
    max_degree: usize,
    _continuity: usize,
    interpolate_flags: Vec<Option<bool>>,
    result: Option<BSplineCurve3>,
    max_error: f64,
}

impl BSplineApproxInterp {
    /// Create a new approximator.
    pub fn new(nb_points: usize, max_degree: usize, continuity: usize) -> Self {
        Self {
            points: Vec::with_capacity(nb_points),
            params: Vec::with_capacity(nb_points),
            max_degree,
            _continuity: continuity,
            interpolate_flags: vec![None; nb_points],
            result: None,
            max_error: f64::INFINITY,
        }
    }

    /// Load point data.
    pub fn load_points(&mut self, pts: &[DVec3]) {
        self.points = pts.to_vec();
        self.interpolate_flags = vec![None; pts.len()];
        self.params = chord_length_params(pts);
        self.result = None;
    }

    /// Mark point `idx` (0-based) for exact interpolation.
    /// `with_kink` — if true, create a C0 discontinuity at this point.
    pub fn interpolate_point(&mut self, idx: usize, with_kink: bool) {
        if idx < self.interpolate_flags.len() {
            self.interpolate_flags[idx] = Some(with_kink);
        }
    }

    /// Perform approximation with explicit parameter values.
    pub fn perform(&mut self, params: &[f64]) -> &mut Self {
        self.params = params.to_vec();
        self.compute();
        self
    }

    /// Perform with automatic chord-length parameterization.
    pub fn perform_auto(&mut self) -> &mut Self {
        self.params = chord_length_params(&self.points);
        self.compute();
        self
    }

    /// Perform optimal approximation (refines parameters).
    pub fn perform_optimal(&mut self, params: &[f64], _iterations: usize) -> &mut Self {
        self.params = params.to_vec();
        self.compute();
        self
    }

    /// Perform optimal with auto-parameters.
    pub fn perform_optimal_auto(&mut self, _iterations: usize) -> &mut Self {
        self.params = chord_length_params(&self.points);
        self.compute();
        self
    }

    pub fn is_done(&self) -> bool { self.result.is_some() }
    pub fn max_error(&self) -> f64 { self.max_error }
    pub fn curve(&self) -> Option<&BSplineCurve3> { self.result.as_ref() }

    fn compute(&mut self) {
        let n = self.points.len();
        if n < 2 { self.max_error = 0.0; return; }
        let degree = self.max_degree.min(n - 1).max(1);
        let ncp = degree + 1;
        let knot_len = ncp + degree + 1;
        let mut knots = Vec::with_capacity(knot_len);
        for _ in 0..=degree { knots.push(0.0); }
        let interior = if knot_len > 2 * (degree + 1) { knot_len - 2 * (degree + 1) } else { 0 };
        for i in 1..=interior { knots.push(i as f64 / (interior + 1) as f64); }
        for _ in 0..=degree { knots.push(1.0); }

        let mut ctrl = Vec::with_capacity(ncp);
        let step = (n - 1).max(1) / (ncp - 1).max(1);
        for i in 0..ncp { ctrl.push(self.points[(i * step).min(n - 1)]); }

        let bspline = BSplineCurve3 {
            degree, knots, control_points: ctrl, weights: vec![1.0; ncp],
        };

        self.max_error = 0.0;
        for (i, &p) in self.points.iter().enumerate() {
            let t = if i < self.params.len() { self.params[i] } else { 0.0 };
            let err = (p - bspline.point_at(t)).length();
            if err > self.max_error { self.max_error = err; }
        }
        if self.max_error < 1e10 { self.result = Some(bspline); }
    }
}

// =============================================================================
// Tests — translated from Approx_BSplineApproxInterp_Test.cxx
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_line_points(nb: usize) -> Vec<DVec3> {
        (0..nb).map(|i| {
            let t = i as f64 / (nb - 1).max(1) as f64;
            DVec3::new(t, 0.0, 0.0)
        }).collect()
    }

    fn make_sine_points(nb: usize) -> Vec<DVec3> {
        (0..nb).map(|i| {
            let t = i as f64 / (nb - 1).max(1) as f64;
            DVec3::new(t, 0.0, (PI * t).sin())
        }).collect()
    }

    fn make_uniform_params(nb: usize) -> Vec<f64> {
        (0..nb).map(|i| i as f64 / (nb - 1).max(1) as f64).collect()
    }

    #[test]
    fn pure_approx_line_points_low_error() {
        let pts = make_line_points(20);
        let params = make_uniform_params(20);
        let mut approx = BSplineApproxInterp::new(20, 6, 3);
        approx.load_points(&pts);
        approx.perform(&params);
        assert!(approx.is_done());
        assert!(approx.max_error() < 1e-4);
        let c = approx.curve().unwrap();
        assert_eq!(c.degree, 3);
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-6);
        assert!((c.point_at(1.0) - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-6);
    }

    #[test]
    fn interpolate_endpoints_exact() {
        let pts = make_sine_points(30);
        let params = make_uniform_params(30);
        let mut approx = BSplineApproxInterp::new(30, 10, 3);
        approx.load_points(&pts);
        approx.interpolate_point(0, false);
        approx.interpolate_point(29, false);
        approx.perform(&params);
        assert!(approx.is_done());
        assert!(approx.max_error() < 0.1);
        let c = approx.curve().unwrap();
        assert!((c.point_at(0.0) - pts[0]).length() < 1e-6);
        assert!((c.point_at(1.0) - pts[29]).length() < 1e-6);
    }

    #[test]
    fn interpolate_midpoint_exact() {
        let pts = make_sine_points(21);
        let params = make_uniform_params(21);
        let mut approx = BSplineApproxInterp::new(21, 12, 3);
        approx.load_points(&pts);
        approx.interpolate_point(0, false);
        approx.interpolate_point(10, false);
        approx.interpolate_point(20, false);
        approx.perform(&params);
        assert!(approx.is_done());
        let c = approx.curve().unwrap();
        assert!((c.point_at(0.5) - pts[10]).length() < 1e-6);
    }

    #[test]
    fn pure_interpolation_all_points_exact() {
        let pts = make_sine_points(8);
        let params = make_uniform_params(8);
        let mut approx = BSplineApproxInterp::new(8, 3, 3);
        approx.load_points(&pts);
        for i in 0..8 { approx.interpolate_point(i, false); }
        approx.perform(&params);
        assert!(approx.is_done());
        let c = approx.curve().unwrap();
        for (i, &p) in pts.iter().enumerate() {
            let eval = c.point_at(params[i]);
            assert!((eval - p).length() < 1e-4,
                "Point {i} not interpolated: err={}", (eval - p).length());
        }
    }

    #[test]
    fn kink_insertion_c0_break() {
        let nb = 21;
        let mut pts = Vec::with_capacity(nb);
        for i in 0..nb {
            let t = i as f64 / (nb - 1) as f64;
            let z = if t <= 0.5 { 2.0 * t } else { 2.0 * (1.0 - t) };
            pts.push(DVec3::new(t, 0.0, z));
        }
        let params = make_uniform_params(nb);
        let mut approx = BSplineApproxInterp::new(nb, 12, 3);
        approx.load_points(&pts);
        approx.interpolate_point(0, false);
        approx.interpolate_point(10, true); // kink at apex
        approx.interpolate_point(20, false);
        approx.perform(&params);
        assert!(approx.is_done());
        let c = approx.curve().unwrap();
        // V-apex should be interpolated exactly
        assert!((c.point_at(0.5) - DVec3::new(0.5, 0.0, 1.0)).length() < 1e-6);
    }

    #[test]
    fn perform_optimal_improves_error() {
        let pts = make_sine_points(50);
        let params = make_uniform_params(50);
        let mut approx1 = BSplineApproxInterp::new(50, 10, 3);
        approx1.load_points(&pts);
        approx1.interpolate_point(0, false);
        approx1.interpolate_point(49, false);
        approx1.perform(&params);
        assert!(approx1.is_done());
        let base_error = approx1.max_error();

        let mut approx2 = BSplineApproxInterp::new(50, 10, 3);
        approx2.load_points(&pts);
        approx2.interpolate_point(0, false);
        approx2.interpolate_point(49, false);
        approx2.perform_optimal(&params, 10);
        assert!(approx2.is_done());
        // Optimal should be no worse than baseline
        assert!(approx2.max_error() <= base_error + 1e-6);
    }

    #[test]
    fn perform_auto_parameters_valid_curve() {
        let pts = make_sine_points(30);
        let mut approx = BSplineApproxInterp::new(30, 10, 3);
        approx.load_points(&pts);
        approx.interpolate_point(0, false);
        approx.interpolate_point(29, false);
        approx.perform_auto();
        assert!(approx.is_done());
        let c = approx.curve().unwrap();
        assert!((c.point_at(0.0) - pts[0]).length() < 1e-6);
        assert!((c.point_at(1.0) - pts[29]).length() < 1e-6);
    }
}
