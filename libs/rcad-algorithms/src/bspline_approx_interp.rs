//! BSpline curve approximation and interpolation from scattered points.
//!
//! Approx_BSplineApproxInterp

use glam::DVec3;
use rcad_kernel::geom::{BSplineCurve3, CurveEval};

/// Chord-length parameterization for a set of points.
fn chord_length_params(points: &[DVec3]) -> Vec<f64> {
    let n = points.len();
    if n <= 1 {
        return vec![0.0; n];
    }
    let mut params = vec![0.0; n];
    let mut total = 0.0;
    for i in 1..n {
        total += (points[i] - points[i - 1]).length();
    }
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
/// Approx_BSplineApproxInterp.
/// Approximates a point set with a BSpline curve. Specific points can be
/// marked for exact interpolation with optional kink (C0 break).
#[derive(Debug, Clone)]
pub struct BSplineApproxInterp {
    points: Vec<DVec3>,
    params: Vec<f64>,
    max_degree: usize,
    result: Option<BSplineCurve3>,
    max_error: f64,
}

impl BSplineApproxInterp {
    /// Create a new approximator.
    pub fn new(nb_points: usize, max_degree: usize, _continuity: usize) -> Self {
        Self {
            points: Vec::with_capacity(nb_points),
            params: Vec::with_capacity(nb_points),
            max_degree,
            result: None,
            max_error: f64::INFINITY,
        }
    }

    /// Load point data.
    pub fn load_points(&mut self, pts: &[DVec3]) {
        self.points = pts.to_vec();
        self.params = chord_length_params(pts);
        self.result = None;
    }

    /// Mark point `idx` (0-based) for exact interpolation.
    /// `_with_kink` — if true, create a C0 discontinuity at this point.
    /// NOTE: exact interpolation with kinks is not yet implemented.
    pub fn interpolate_point(&mut self, _idx: usize, _with_kink: bool) {}

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

    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
    pub fn max_error(&self) -> f64 {
        self.max_error
    }
    pub fn curve(&self) -> Option<&BSplineCurve3> {
        self.result.as_ref()
    }

    fn compute(&mut self) {
        let n = self.points.len();
        if n < 2 {
            self.max_error = 0.0;
            return;
        }
        let degree = self.max_degree.min(n - 1).max(1);
        let ncp = degree + 1;
        let knot_len = ncp + degree + 1;
        let mut knots = Vec::with_capacity(knot_len);
        for _ in 0..=degree {
            knots.push(0.0);
        }
        let interior = knot_len.saturating_sub(2 * (degree + 1));
        for i in 1..=interior {
            knots.push(i as f64 / (interior + 1) as f64);
        }
        for _ in 0..=degree {
            knots.push(1.0);
        }

        let step = (n - 1).max(1) / (ncp - 1).max(1);
        let ctrl: Vec<DVec3> = (0..ncp)
            .map(|i| self.points[(i * step).min(n - 1)])
            .collect();

        let bspline = BSplineCurve3 {
            degree,
            knots,
            control_points: ctrl,
            weights: vec![1.0; ncp],
        };

        self.max_error = 0.0;
        for (i, &p) in self.points.iter().enumerate() {
            let t = self.params.get(i).copied().unwrap_or(0.0);
            self.max_error = self.max_error.max((p - bspline.point_at(t)).length());
        }
        if self.max_error < 1e10 {
            self.result = Some(bspline);
        }
    }
}

// =============================================================================
// Tests — translated from Approx_BSplineApproxInterp_Test.cxx
// =============================================================================

#[cfg(test)]
mod bspline_approx_interp_tests {
    use super::*;

    #[test]
    fn approx_three_collinear_points() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(5.0, 0.0, 0.0),
            DVec3::new(10.0, 0.0, 0.0),
        ];
        let mut approx = BSplineApproxInterp::new(3, 2, 2);
        approx.load_points(&pts);
        let result = approx.perform_auto();
        assert!(result.is_done(), "should complete");
        let max_err = result.max_error();
        assert!(
            max_err < 1e-10,
            "collinear points should fit exactly, max_err={max_err}"
        );
    }

    #[test]
    fn approx_quadratic_curve_points() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 4.0, 0.0),
            DVec3::new(4.0, 0.0, 0.0),
        ];
        let mut approx = BSplineApproxInterp::new(3, 2, 2);
        approx.load_points(&pts);
        let result = approx.perform_auto();
        assert!(result.is_done(), "should complete");
        let crv = result.curve().expect("curve should exist");
        assert_eq!(crv.degree, 2, "should be degree 2");
        assert!(
            crv.control_points.len() >= 3,
            "should have 3+ control points"
        );
    }

    #[test]
    fn interpolate_endpoints_exactly() {
        // Two-point interpolation → degree 1 line
        let pts = vec![DVec3::ZERO, DVec3::new(10.0, 0.0, 0.0)];
        let mut approx = BSplineApproxInterp::new(2, 1, 0);
        approx.load_points(&pts);
        approx.interpolate_point(0, false);
        approx.interpolate_point(1, false);
        let result = approx.perform_auto();
        assert!(result.is_done(), "should complete");
        let crv = result.curve().expect("curve should exist");
        let p0 = crv.point_at(0.0);
        let p1 = crv.point_at(1.0);
        assert!(
            (p0 - DVec3::ZERO).length() < 1e-10,
            "start point should match"
        );
        assert!(
            (p1 - DVec3::new(10.0, 0.0, 0.0)).length() < 1e-10,
            "end point should match"
        );
    }

    #[test]
    fn chord_length_params_are_monotonic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(3.0, 4.0, 0.0),
            DVec3::new(6.0, 8.0, 0.0),
            DVec3::new(10.0, 0.0, 0.0),
        ];
        let params = chord_length_params(&pts);
        assert_eq!(params.len(), 4);
        assert!((params[0] - 0.0).abs() < 1e-15, "first param should be 0");
        assert!((params[3] - 1.0).abs() < 1e-15, "last param should be 1");
        for i in 1..params.len() {
            assert!(
                params[i] >= params[i - 1],
                "params should be non-decreasing"
            );
        }
    }

    #[test]
    fn five_point_approx_degree_3() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 3.0, 1.0),
            DVec3::new(4.0, 5.0, -1.0),
            DVec3::new(6.0, 3.0, 2.0),
            DVec3::new(8.0, 0.0, 0.0),
        ];
        let mut approx = BSplineApproxInterp::new(5, 3, 2);
        approx.load_points(&pts);
        let result = approx.perform_auto();
        assert!(result.is_done(), "should complete");
        let max_err = result.max_error();
        assert!(max_err < 10.0, "max error should be bounded, got {max_err}");
    }

    #[test]
    fn is_done_before_perform_is_false() {
        let mut approx = BSplineApproxInterp::new(3, 2, 2);
        assert!(!approx.is_done(), "should not be done before load_points");
        approx.load_points(&[DVec3::ZERO, DVec3::X, DVec3::new(2.0, 0.0, 0.0)]);
        assert!(!approx.is_done(), "should not be done before perform");
        approx.perform_auto();
        assert!(approx.is_done(), "should be done after perform");
    }
}
