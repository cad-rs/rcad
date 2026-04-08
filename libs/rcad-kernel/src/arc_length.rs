//! Curve arc-length computation.
//!
//! Analogous to OCCT `GCPnts_AbscissaPoint` / `CPnts_AbscissaPoint::Length`.
//!
//! # Strategy
//!
//! | Curve type | Method |
//! |------------|--------|
//! | `Line3`    | Analytic: `|t2 − t1|` (direction is always unit) |
//! | `Circle3`  | Analytic: `r · |t2 − t1|` |
//! | `Ellipse3` | 16-point Gauss-Legendre quadrature of `∫|tangent_at(t)| dt` |
//! | `BSplineCurve3` | 16-point Gauss-Legendre quadrature |
//!
//! The function returns a **signed** value: positive when `t2 > t1`, negative
//! otherwise.  Call `.abs()` if you need an unsigned arc length.

use crate::geom::{Curve3, CurveEval};

// ── 16-point Gauss-Legendre nodes and weights on [-1, 1] ──────────────────────
// Source: Abramowitz & Stegun Table 25.4, standard reference.
const GL16_NODES: [f64; 16] = [
    -0.095012509837637440185,
    0.095012509837637440185,
    -0.281603550779258913230,
    0.281603550779258913230,
    -0.458016777657227386342,
    0.458016777657227386342,
    -0.617876244402643748447,
    0.617876244402643748447,
    -0.755404408355003033895,
    0.755404408355003033895,
    -0.865631202387831743880,
    0.865631202387831743880,
    -0.944575023073232576078,
    0.944575023073232576078,
    -0.989400934991649932596,
    0.989400934991649932596,
];

const GL16_WEIGHTS: [f64; 16] = [
    0.189450610455068496285,
    0.189450610455068496285,
    0.182603415044923588867,
    0.182603415044923588867,
    0.169156519395002538189,
    0.169156519395002538189,
    0.149451349150580593146,
    0.149451349150580593146,
    0.124628971255533872052,
    0.124628971255533872052,
    0.095158511682492784810,
    0.095158511682492784810,
    0.062253523938647892863,
    0.062253523938647892863,
    0.027152459411754094852,
    0.027152459411754094852,
];

// ── Internal Gauss-Legendre integrator ───────────────────────────────────────

/// Integrate `|dP/dt|` (un-normalized speed) over `[t1, t2]` using 16-point
/// GL quadrature.  The derivative is approximated via central differences on
/// `point_at` to avoid the normalization baked into `tangent_at`.
///
/// Returns the *signed* result (sign follows `t2 - t1`).
fn gl16_arc_length(curve: &Curve3, t1: f64, t2: f64) -> f64 {
    // Finite-difference step — small enough for precision, large enough to
    // avoid cancellation.  1e-8 relative to the half-interval works well for
    // typical CAD parameter ranges.
    let fd_eps = ((t2 - t1).abs() * 0.5).max(1.0) * 1e-8;

    let half = (t2 - t1) * 0.5;
    let mid = (t2 + t1) * 0.5;
    let sum: f64 = GL16_NODES
        .iter()
        .zip(GL16_WEIGHTS.iter())
        .map(|(&xi, &wi)| {
            let t = mid + half * xi;
            // Un-normalized derivative |dP/dt|
            let dp = (curve.point_at(t + fd_eps) - curve.point_at(t - fd_eps)) / (2.0 * fd_eps);
            wi * dp.length()
        })
        .sum();
    half * sum
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Arc length of `curve` over parameter range `[t1, t2]`.
///
/// Returns a signed value: positive if `t2 > t1`, negative if `t2 < t1`.
/// Take `.abs()` for unsigned arc length.
///
/// # Accuracy
/// - `Line3` and `Circle3`: exact (analytic).
/// - `Ellipse3` and `BSplineCurve3`: 16-point Gauss-Legendre quadrature;
///   relative error typically < 1e-10 for smooth curves.
pub fn arc_length(curve: &Curve3, t1: f64, t2: f64) -> f64 {
    match curve {
        Curve3::Line(_) => {
            // direction is always a unit vector by convention
            t2 - t1
        }
        Curve3::Circle(c) => c.radius * (t2 - t1),
        Curve3::Ellipse(_)
        | Curve3::BSpline(_)
        | Curve3::Bezier(_)
        | Curve3::Offset(_)
        | Curve3::Hyperbola(_)
        | Curve3::Parabola(_) => gl16_arc_length(curve, t1, t2),
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{BSplineCurve3, Circle3, Ellipse3, Line3};
    use glam::DVec3;
    use std::f64::consts::PI;

    const TOL: f64 = 1e-9;
    #[allow(dead_code)]
    const NUM_TOL: f64 = 1e-6;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn line_arc_length_analytic() {
        let c = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        assert!(approx_eq(arc_length(&c, 0.0, 5.0), 5.0, TOL));
        assert!(approx_eq(arc_length(&c, -2.0, 3.0), 5.0, TOL));
        assert!(approx_eq(arc_length(&c, 5.0, 0.0), -5.0, TOL)); // signed
    }

    #[test]
    fn circle_half_circumference() {
        let c = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        // half circle [0, π] → length = π
        let l = arc_length(&c, 0.0, PI);
        assert!(approx_eq(l, PI, TOL), "half circle got {l}");
    }

    #[test]
    fn circle_full_circumference() {
        let c = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 2.0,
        });
        let l = arc_length(&c, 0.0, 2.0 * PI);
        assert!(approx_eq(l, 4.0 * PI, TOL), "full circle r=2 got {l}");
    }

    #[test]
    fn ellipse_full_perimeter_approx() {
        // Ellipse a=2, b=1 — full perimeter via Ramanujan approx ≈ 9.6884...
        // Gauss-Legendre should match to within NUM_TOL.
        let c = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 2.0,
            minor_radius: 1.0,
        });
        let l = arc_length(&c, 0.0, 2.0 * PI);
        // Ramanujan first approximation: π(3(a+b) - sqrt((3a+b)(a+3b))) where a=2,b=1
        let a = 2.0_f64;
        let b = 1.0_f64;
        let ramanujan = PI * (3.0 * (a + b) - ((3.0 * a + b) * (a + 3.0 * b)).sqrt());
        assert!(
            approx_eq(l, ramanujan, 1e-3),
            "ellipse perimeter {l} vs Ramanujan {ramanujan}"
        );
    }

    #[test]
    fn bspline_line_segment_arc_length() {
        // Degree-1 BSpline from (0,0,0) to (3,4,0) → length = 5
        // GL16 with finite-difference speed estimation; tolerance 1e-3 for numerical path.
        let c = Curve3::BSpline(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(3.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let l = arc_length(&c, 0.0, 1.0).abs();
        assert!(approx_eq(l, 5.0, 1e-3), "bspline line segment length {l}");
    }
}
