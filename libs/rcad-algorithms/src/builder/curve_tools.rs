// �?OCCT-aligned: ElCLib derivative functions for 2D elementary curves
//   (ElCLib.hxx: LineD1/D2, CircleD1/D2, EllipseD1/D2).
//   curve2d_resolution approximates Geom2dAdaptor::Resolution.
//   curve2d_lprop_curvature matches ElCLib::Curvature formula |d1×d2|/|d1|³.

use glam::DVec2;
use rcad_kernel::geom::*;
use crate::tolerance::TOLERANCE_CLAMP_MIN;

/// OCCT-aligned: ElCLib::D1 �?first derivative of a 2D curve.
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

/// �?OCCT-aligned: Geom2dAdaptor_Curve::Resolution (L1186-1219).
///   Line: return Ruv (gp_Lin2d direction is always unit)
///   Circle: 2*asin(Ruv/(2*R)) or 2*PI
///   Ellipse: Ruv / MajorRadius
///   Bezier: BSplCLib::Resolution (via max derivative sampling)
///   BSpline: BSplCLib::Resolution (at knots + midpoints)
///   default: Precision::Parametric(Ruv) = Ruv * 0.01
///   Trimmed: unwrap to base type then dispatch.
pub fn curve2d_resolution(curve: &Curve2d, r_uv: f64) -> f64 {
    match curve {
        Curve2d::Line(_) => r_uv,
        Curve2d::Circle(c) => {
            let r = c.radius;
            if r > r_uv / 2.0 { 2.0 * f64::asin(r_uv / (2.0 * r)) } else { std::f64::consts::TAU }
        }
        Curve2d::Ellipse(e) => r_uv / e.major_radius,
        Curve2d::Bezier(b) => { let ms = sample_max_speed_bezier_2d(b); if ms > TOLERANCE_CLAMP_MIN { r_uv / ms } else { r_uv } }
        Curve2d::BSpline(b) => { let ms = sample_max_speed_bspline_2d(b); if ms > TOLERANCE_CLAMP_MIN { r_uv / ms } else { r_uv } }
        Curve2d::Trimmed(tc) => curve2d_resolution(&tc.curve, r_uv),
        _ => r_uv * 0.01,
    }
}

/// OCCT-aligned: BSplCLib::Resolution �?evaluates |dC/dt| at endpoints and
/// midpoint, returns the MAXIMUM derivative (Resolution = Ruv / max|dC/dt|).
fn sample_max_speed_bezier_2d(b: &BezierCurve2) -> f64 {
    let d0 = curve2d_d1(&Curve2d::Bezier(b.clone()), 0.0).length();
    let d1 = curve2d_d1(&Curve2d::Bezier(b.clone()), 1.0).length();
    let dm = curve2d_d1(&Curve2d::Bezier(b.clone()), 0.5).length();
    d0.max(d1).max(dm)
}

/// OCCT-aligned: BSplCLib::Resolution for BSpline �?evaluates at knots and
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


