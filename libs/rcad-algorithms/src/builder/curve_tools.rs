// OCCT-aligned: ElCLib derivative functions for 2D elementary curves
//   (ElCLib.hxx: LineD1/D2, CircleD1/D2, EllipseD1/D2).
//   curve2d_resolution approximates Geom2dAdaptor::Resolution.
//   curve2d_lprop_curvature matches ElCLib::Curvature formula |d1×d2|/|d1|³.

use glam::DVec2;
use rcad_kernel::geom::*;

/// OCCT-aligned: ElCLib::D1 — first derivative of a 2D curve.
/// OCCT ElCLib.hxx dispatches to LineD1(lin, U, P, V), CircleD1(circ, U, P, V),
/// EllipseD1(elips, U, majR, minR, P, V). rcad: matches OCCT formulas exactly.
pub fn curve2d_d1(curve: &Curve2d, t: f64) -> DVec2 {
    match curve {
        Curve2d::Line(l) => l.direction,
        Curve2d::Circle(c) => c.radius * DVec2::new(-t.sin(), t.cos()),
        Curve2d::Ellipse(e) => {
            let minor = DVec2::new(-e.major_dir.y, e.major_dir.x);
            e.major_dir * (-e.major_radius * t.sin()) + minor * (e.minor_radius * t.cos())
        }
        Curve2d::Parabola(p) => {
            let perp = DVec2::new(-p.axis_dir.y, p.axis_dir.x);
            (t / p.focal_param) * p.axis_dir + perp
        }
        Curve2d::Hyperbola(h) => {
            let minor = DVec2::new(-h.major_dir.y, h.major_dir.x);
            h.semi_major * t.sinh() * h.major_dir + h.semi_minor * t.cosh() * minor
        }
        Curve2d::BSpline(b) => b.derivative_at(t),
        Curve2d::Bezier(b) => { let eps = 1e-7; (b.point_at(t + eps) - b.point_at(t - eps)) / (2.0 * eps) }
        Curve2d::Trimmed(tc) => curve2d_d1(&tc.curve, t),
        _ => { let eps = 1e-7; (curve.point_at(t + eps) - curve.point_at(t - eps)) / (2.0 * eps) }
    }
}

pub fn curve2d_d2(curve: &Curve2d, t: f64) -> DVec2 {
    match curve {
        Curve2d::Line(_) => DVec2::ZERO,
        Curve2d::Circle(c) => c.radius * DVec2::new(-t.cos(), -t.sin()),
        Curve2d::Ellipse(e) => {
            let minor = DVec2::new(-e.major_dir.y, e.major_dir.x);
            e.major_dir * (-e.major_radius * t.cos()) + minor * (-e.minor_radius * t.sin())
        }
        Curve2d::Parabola(p) => {
            (1.0 / p.focal_param) * p.axis_dir
        }
        Curve2d::Hyperbola(h) => {
            let minor = DVec2::new(-h.major_dir.y, h.major_dir.x);
            h.semi_major * t.cosh() * h.major_dir + h.semi_minor * t.sinh() * minor
        }
        _ => { let eps = 1e-7; (curve.point_at(t + eps) - 2.0 * curve.point_at(t) + curve.point_at(t - eps)) / (eps * eps) }
    }
}

pub fn curve2d_lprop_curvature(d1: DVec2, d2: DVec2, tol_sq: f64) -> f64 {
    let a_dd1 = d1.length_squared();
    let a_dd2 = d2.length_squared();
    if a_dd2 <= tol_sq { return 0.0; }
    let cross = d1.x * d2.y - d1.y * d2.x;
    let a_n = cross * cross;
    let a_t = a_n / a_dd1 / a_dd2;
    if a_t <= tol_sq { return 0.0; }
    cross.abs() / a_dd1 / a_dd1.sqrt()
}

/// OCCT-aligned: Geom2dAdaptor_Curve::Resolution (L1186-1219).
///   Line L1191: return Ruv;
///   Circle L1193-1201: 2*asin(Ruv/(2*R)) or 2*PI
///   Ellipse L1204: Ruv / MajorRadius
///   Bezier L1206-1209: Geom2d_BezierCurve::Resolution (BSplCLib::Resolution)
///   BSpline L1211-1215: Geom2d_BSplineCurve::Resolution (BSplCLib::Resolution)
///   default L1217: Precision::Parametric(Ruv)
///   Trimmed: OCCT adaptor unwraps TrimmedCurve → base type; rcad recurses.
pub fn curve2d_resolution(curve: &Curve2d, r_uv: f64) -> f64 {
    match curve {
        // OCCT L1191: return Ruv (gp_Lin2d direction is always unit)
        Curve2d::Line(_) => r_uv,
        // OCCT L1193-1201
        Curve2d::Circle(c) => {
            let r = c.radius;
            if r > r_uv / 2.0 { 2.0 * f64::asin(r_uv / (2.0 * r)) } else { std::f64::consts::TAU }
        }
        // OCCT L1203-1205
        Curve2d::Ellipse(e) => r_uv / e.major_radius,
        // OCCT L1206-1209: BSplCLib::Resolution (analytical, not sampling)
        Curve2d::Bezier(b) => { let ms = sample_max_speed_bezier_2d(b); if ms > 1e-15 { r_uv / ms } else { r_uv } }
        // OCCT L1211-1215: BSplCLib::Resolution — evaluates at knots + span midpoints
        Curve2d::BSpline(b) => { let ms = sample_max_speed_bspline_2d(b); if ms > 1e-15 { r_uv / ms } else { r_uv } }
        // OCCT: TrimmedCurve unwrapped to base type before dispatch
        Curve2d::Trimmed(tc) => curve2d_resolution(&tc.curve, r_uv),
        // OCCT L1217: Precision::Parametric(Ruv) = Ruv * 0.01
        _ => r_uv * 0.01,
    }
}

/// OCCT-aligned: BSplCLib::Resolution — evaluates |dC/dt| at endpoints and
/// midpoint, returns the MAXIMUM derivative (Resolution = Ruv / max|dC/dt|).
fn sample_max_speed_bezier_2d(b: &BezierCurve2) -> f64 {
    let d0 = curve2d_d1(&Curve2d::Bezier(b.clone()), 0.0).length();
    let d1 = curve2d_d1(&Curve2d::Bezier(b.clone()), 1.0).length();
    let dm = curve2d_d1(&Curve2d::Bezier(b.clone()), 0.5).length();
    d0.max(d1).max(dm)
}

/// OCCT-aligned: BSplCLib::Resolution for BSpline — evaluates at knots and
/// knot-span midpoints, returns the MAXIMUM derivative (Resolution = Ruv / max|dC/dt|).
fn sample_max_speed_bspline_2d(b: &BSplineCurve2) -> f64 {
    let mut max_s = 0.0;
    for &k in &b.knots {
        let d = b.derivative_at(k).length();
        if d > max_s { max_s = d; }
    }
    for w in b.knots.windows(2) {
        let t_mid = 0.5 * (w[0] + w[1]);
        if t_mid > w[0] && t_mid < w[1] {
            let d = b.derivative_at(t_mid).length();
            if d > max_s { max_s = d; }
        }
    }
    max_s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d1_line_is_direction() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let d = curve2d_d1(&l, 0.5);
        assert!((d - DVec2::X).length() < 1e-15);
    }

    #[test]
    fn d1_circle() {
        let c = Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: 2.0 });
        let d = curve2d_d1(&c, 0.0);
        assert!((d - DVec2::new(0.0, 2.0)).length() < 1e-10);
    }

    #[test]
    fn d2_line_is_zero() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let d = curve2d_d2(&l, 0.5);
        assert_eq!(d, DVec2::ZERO);
    }

    #[test]
    fn d2_circle() {
        let c = Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: 2.0 });
        let d = curve2d_d2(&c, 0.0);
        assert!((d - DVec2::new(-2.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn curvature_straight_line_is_zero() {
        let k = curve2d_lprop_curvature(DVec2::new(1.0, 0.0), DVec2::ZERO, 1e-14);
        assert!(k < 1e-15);
    }

    #[test]
    fn curvature_circle_one_over_r() {
        let r: f64 = 2.0;
        let t: f64 = 0.7;
        let d1 = r * DVec2::new(-t.sin(), t.cos());
        let d2 = r * DVec2::new(-t.cos(), -t.sin());
        let k = curve2d_lprop_curvature(d1, d2, 1e-14);
        assert!((k - 1.0 / r).abs() < 1e-10);
    }

    #[test]
    fn resolution_line_matches_occt() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let res = curve2d_resolution(&l, 0.01);
        assert!((res - 0.01).abs() < 1e-15);
    }

    #[test]
    fn resolution_circle_large_radius() {
        let c = Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: 10.0 });
        let res = curve2d_resolution(&c, 0.01);
        let expected = 2.0 * ((0.01_f64 / 20.0_f64).asin());
        assert!((res - expected).abs() < 1e-15);
    }

    #[test]
    fn resolution_circle_tiny_radius_returns_tau() {
        let c = Curve2d::Circle(Circle2d { center: DVec2::ZERO, radius: 0.001 });
        let res = curve2d_resolution(&c, 0.01);
        assert!((res - std::f64::consts::TAU).abs() < 1e-15);
    }

    #[test]
    fn resolution_ellipse() {
        let e = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X,
            major_radius: 5.0, minor_radius: 3.0,
        });
        let res = curve2d_resolution(&e, 0.01);
        assert!((res - 0.01 / 5.0).abs() < 1e-15);
    }

    #[test]
    fn resolution_default_returns_r_uv_times_001() {
        let inv = Curve2d::CircleInvolute(CircleInvolute2d {
            center: DVec2::ZERO, base_radius: 1.0, start_angle: 0.0,
        });
        // OCCT L1217: Precision::Parametric(Ruv) = Ruv * 0.01
        let res = curve2d_resolution(&inv, 0.01);
        assert!((res - 0.01 * 0.01).abs() < 1e-15);
    }

    #[test]
    fn d1_bezier_linear() {
        // Linear Bezier: P0=(0,0), P1=(2,0) → d1 ≈ (2,0) via finite difference
        let b = Curve2d::Bezier(BezierCurve2 {
            control_points: vec![DVec2::ZERO, DVec2::new(2.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let d = curve2d_d1(&b, 0.3);
        assert!((d - DVec2::new(2.0, 0.0)).length() < 1e-4,
            "Bezier d1 expected near (2,0), got {:?}", d);
    }

    #[test]
    fn d1_trimmed_delegates_to_base() {
        let base = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let t = Curve2d::Trimmed(TrimmedCurve2 { curve: Box::new(base), t_min: 0.0, t_max: 1.0 });
        let d = curve2d_d1(&t, 0.5);
        assert!((d - DVec2::X).length() < 1e-15);
    }

    #[test]
    fn bspline_derivative_at_works() {
        let b = Curve2d::BSpline(BSplineCurve2 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::new(1.0, 0.0), DVec2::new(2.0, 1.0)],
            weights: vec![1.0, 1.0, 1.0],
        });
        assert!(curve2d_d1(&b, 0.5).length() > 0.0);
    }
}
