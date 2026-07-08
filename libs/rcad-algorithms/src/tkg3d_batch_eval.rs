//! OCCT GeomGridEval — batch curve/surface evaluation.
//!
//! OCCT source: src/ModelingData/TKG3d/GeomGridEval/
//!
//! Takes an array of parameters and evaluates all points (and optionally D1/D2)
//! at once using per-type optimized tight loops. This matches the OCCT
//! GeomGridEval_Line/Circle/BSplineCurve/BezierCurve/OtherCurve pattern.

use glam::DVec3;
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;

/// Result of a batch point evaluation.
pub struct BatchPoints {
    pub points: Vec<DVec3>,
}

/// Result of a batch D1 evaluation.
pub struct BatchD1 {
    pub points: Vec<DVec3>,
    pub d1: Vec<DVec3>,
}

/// Result of a batch D2 evaluation.
pub struct BatchD2 {
    pub points: Vec<DVec3>,
    pub d1: Vec<DVec3>,
    pub d2: Vec<DVec3>,
}

/// Batch evaluate (D0) a line at multiple parameters.
pub fn batch_eval_line(line: &Line3, params: &[f64]) -> Vec<DVec3> {
    let o = line.origin;
    let d = line.direction;
    params.iter().map(|&t| o + t * d).collect()
}

/// Batch evaluate (D0) a circle at multiple parameters.
pub fn batch_eval_circle(circle: &Circle3, params: &[f64]) -> Vec<DVec3> {
    let c = circle.center;
    let xx = circle.x_dir;
    let yy = circle.y_dir;
    let r = circle.radius;
    params.iter().map(|&u| {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        c + r * (cos_u * xx + sin_u * yy)
    }).collect()
}

/// Batch evaluate (D0) a circle with position and derivatives.
pub fn batch_eval_circle_d1(circle: &Circle3, params: &[f64]) -> BatchD1 {
    let c = circle.center;
    let xx = circle.x_dir;
    let yy = circle.y_dir;
    let r = circle.radius;
    let mut points = Vec::with_capacity(params.len());
    let mut d1 = Vec::with_capacity(params.len());
    for &u in params {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        points.push(c + r * (cos_u * xx + sin_u * yy));
        d1.push(r * (-sin_u * xx + cos_u * yy));
    }
    BatchD1 { points, d1 }
}

/// Batch evaluate (D0) a circle with position, D1, D2.
pub fn batch_eval_circle_d2(circle: &Circle3, params: &[f64]) -> BatchD2 {
    let c = circle.center;
    let xx = circle.x_dir;
    let yy = circle.y_dir;
    let r = circle.radius;
    let mut points = Vec::with_capacity(params.len());
    let mut d1 = Vec::with_capacity(params.len());
    let mut d2 = Vec::with_capacity(params.len());
    for &u in params {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        points.push(c + r * (cos_u * xx + sin_u * yy));
        d1.push(r * (-sin_u * xx + cos_u * yy));
        d2.push(r * (-cos_u * xx - sin_u * yy));
    }
    BatchD2 { points, d1, d2 }
}

/// Batch evaluate (D0) a BSpline curve at multiple parameters.
pub fn batch_eval_bspline(bs: &BSplineCurve3, params: &[f64]) -> Vec<DVec3> {
    params.iter().map(|&t| bs.point_at(t)).collect()
}

/// Batch evaluate (D0) any curve via generic dispatch.
pub fn batch_eval_curve(curve: &Curve3, params: &[f64]) -> Vec<DVec3> {
    params.iter().map(|&t| curve.point_at(t)).collect()
}

/// GeomGridEval_Curve equivalent: per-type batch evaluation with optimal paths.
pub struct GridEvalCurve;

impl GridEvalCurve {
    /// Evaluate an array of parameters using the best per-type batch evaluator.
    /// Falls back to generic dispatch for unsupported types.
    pub fn evaluate_grid(curve: &Curve3, params: &[f64]) -> Vec<DVec3> {
        match curve {
            Curve3::Line(l) => Self::eval_line(l, params),
            Curve3::Circle(c) => Self::eval_circle(c, params),
            _ => params.iter().map(|&t| curve.point_at(t)).collect(),
        }
    }

    fn eval_line(line: &Line3, params: &[f64]) -> Vec<DVec3> {
        batch_eval_line(line, params)
    }

    fn eval_circle(circle: &Circle3, params: &[f64]) -> Vec<DVec3> {
        batch_eval_circle(circle, params)
    }
}

// =============================================================================
// Tests — GeomGridEval GTests
// =============================================================================

#[cfg(test)]
mod grideval_tests {
    use super::*;

    fn uniform_params(first: f64, last: f64, n: usize) -> Vec<f64> {
        let step = if n > 1 { (last - first) / (n - 1) as f64 } else { 0.0 };
        (0..n).map(|i| first + i as f64 * step).collect()
    }

    // ── Line ─────────────────────────────────────────────────────────

    #[test]
    fn grideval_line_basic() {
        let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
        let params = uniform_params(0.0, 10.0, 11);
        let pts = batch_eval_line(&line, &params);
        assert_eq!(pts.len(), 11);
        for (i, &t) in params.iter().enumerate() {
            assert!((pts[i] - DVec3::new(t, 0.0, 0.0)).length() < TOL);
        }
    }

    #[test]
    fn grideval_line_non_origin() {
        let line = Line3 {
            origin: DVec3::new(1.0, 2.0, 3.0),
            direction: DVec3::new(1.0, 1.0, 1.0).normalize(),
        };
        let params = uniform_params(0.0, 5.0, 6);
        let pts = batch_eval_line(&line, &params);
        assert!((pts[0] - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
        // At t=5, point moves 5 units along normalized (1,1,1)
        let d = 5.0 / 3.0_f64.sqrt();
        assert!((pts[5] - DVec3::new(1.0 + d, 2.0 + d, 3.0 + d)).length() < TOL);
    }

    // ── Circle ───────────────────────────────────────────────────────

    #[test]
    fn grideval_circle_basic() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let params = vec![0.0, std::f64::consts::PI / 2.0, std::f64::consts::PI,
                          3.0 * std::f64::consts::PI / 2.0, 2.0 * std::f64::consts::PI];
        let pts = batch_eval_circle(&circle, &params);
        // Circle3::new: normal=Z → x_dir=Y, y_dir=-X
        // P(0) = R*Y = (0,2,0)
        assert!((pts[0] - DVec3::new(0.0, 2.0, 0.0)).length() < TOL);
        assert!((pts[1] - DVec3::new(-2.0, 0.0, 0.0)).length() < TOL);
        assert!((pts[2] - DVec3::new(0.0, -2.0, 0.0)).length() < TOL);
        assert!((pts[3] - DVec3::new(2.0, 0.0, 0.0)).length() < TOL);
        assert!((pts[4] - pts[0]).length() < TOL);
    }

    #[test]
    fn grideval_circle_non_standard() {
        // Circle in YZ plane (normal = X), center at (1,0,0)
        let circle = Circle3::new(DVec3::new(1.0, 0.0, 0.0), DVec3::X, 3.0);
        // Circle3::new: normal=X → x_dir=any_perpendicular(X)=Z, y_dir=X×Z=-Y
        // P(0) = center + 3*Z = (1,0,3)
        let params = uniform_params(0.0, 2.0 * std::f64::consts::PI, 9);
        let pts = batch_eval_circle(&circle, &params);
        for pt in &pts {
            assert!((pt.x - 1.0).abs() < TOL);
            let r = (*pt - DVec3::new(1.0, 0.0, 0.0)).length();
            assert!((r - 3.0).abs() < TOL);
        }
    }

    #[test]
    fn grideval_circle_d1_accuracy() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let params = uniform_params(0.0, 2.0 * std::f64::consts::PI, 9);
        let result = batch_eval_circle_d1(&circle, &params);
        assert_eq!(result.points.len(), 9);
        assert_eq!(result.d1.len(), 9);
        // D1 magnitude should be R = 5
        for d in &result.d1 {
            assert!((d.length() - 5.0).abs() < TOL);
        }
    }

    #[test]
    fn grideval_circle_d2_accuracy() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let params = uniform_params(0.0, 2.0 * std::f64::consts::PI, 9);
        let result = batch_eval_circle_d2(&circle, &params);
        assert_eq!(result.points.len(), 9);
        assert_eq!(result.d2.len(), 9);
        // D2 magnitude should also be R = 5
        for d in &result.d2 {
            assert!((d.length() - 5.0).abs() < TOL);
        }
    }

    // ── BSpline ──────────────────────────────────────────────────────

    fn make_bspline() -> BSplineCurve3 {
        BSplineCurve3 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO, DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(3.0, 2.0, 0.0), DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        }
    }

    #[test]
    fn grideval_bspline_basic() {
        let bs = make_bspline();
        let bs_curve = Curve3::BSpline(bs.clone());
        let params = uniform_params(0.0, 1.0, 11);
        // Compare batch eval vs direct point_at
        for (&t, pt) in params.iter().zip(batch_eval_bspline(&bs, &params)) {
            let expected = bs_curve.point_at(t);
            assert!((pt - expected).length() < 1e-9);
        }
    }

    #[test]
    fn grideval_bspline_endpoints() {
        let bs = make_bspline();
        let pts = batch_eval_bspline(&bs, &[0.0, 1.0]);
        assert!((pts[0] - DVec3::ZERO).length() < TOL);
        assert!((pts[1] - DVec3::new(4.0, 0.0, 0.0)).length() < TOL);
    }

    // ── GridEvalCurve (dispatcher) ────────────────────────────────────

    #[test]
    fn grideval_curve_dispatcher_line() {
        let curve = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let params = uniform_params(0.0, 10.0, 5);
        let pts = GridEvalCurve::evaluate_grid(&curve, &params);
        assert_eq!(pts.len(), 5);
        assert!((pts[0] - DVec3::ZERO).length() < TOL);
        assert!((pts[4] - DVec3::new(10.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn grideval_curve_dispatcher_circle() {
        let curve = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 2.0));
        let params = uniform_params(0.0, 6.28318, 5);
        let pts = GridEvalCurve::evaluate_grid(&curve, &params);
        assert_eq!(pts.len(), 5);
        for pt in &pts {
            let r = pt.length();
            assert!((r - 2.0).abs() < 1e-4);
        }
    }

    #[test]
    fn grideval_curve_dispatcher_ellipse_fallback() {
        // Ellipse not directly optimized → falls through to generic dispatch
        let curve = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 3.0, minor_radius: 2.0,
        });
        let params = uniform_params(0.0, 6.28318, 9);
        let pts = GridEvalCurve::evaluate_grid(&curve, &params);
        assert_eq!(pts.len(), 9);
        for pt in &pts {
            // Should satisfy ellipse equation x²/a² + y²/b² = 1
            let val = pt.x*pt.x/9.0 + pt.y*pt.y/4.0;
            assert!((val - 1.0).abs() < 1e-6);
        }
    }
}
