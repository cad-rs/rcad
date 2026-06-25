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

pub fn curve2d_resolution(curve: &Curve2d, r_uv: f64) -> f64 {
    match curve {
        Curve2d::Line(_) => r_uv,
        Curve2d::Circle(c) => {
            let r = c.radius;
            if r > r_uv / 2.0 { 2.0 * (r_uv / (2.0 * r)).asin() } else { std::f64::consts::TAU }
        }
        Curve2d::Ellipse(e) => r_uv / e.major_radius,
        Curve2d::Bezier(b) => { let ms = sample_max_speed_2d(b, 16); if ms > 1e-15 { r_uv / ms } else { r_uv } }
        Curve2d::BSpline(b) => { let ms = sample_max_speed_bspline_2d(b, 32); if ms > 1e-15 { r_uv / ms } else { r_uv } }
        Curve2d::Trimmed(tc) => curve2d_resolution(&tc.curve, r_uv),
        _ => r_uv * 0.01,
    }
}

fn sample_max_speed_2d(b: &BezierCurve2, n: usize) -> f64 {
    let mut ms: f64 = 0.0;
    for i in 0..=n { let t = i as f64 / n as f64; let d = curve2d_d1(&Curve2d::Bezier(b.clone()), t); ms = ms.max(d.length()); }
    ms
}

fn sample_max_speed_bspline_2d(b: &BSplineCurve2, n: usize) -> f64 {
    let mut ms: f64 = 0.0;
    let t0 = b.knots.first().copied().unwrap_or(0.0);
    let t1 = b.knots.last().copied().unwrap_or(1.0);
    let span = (t1 - t0).max(1e-15);
    for i in 0..=n { let t = t0 + (i as f64 / n as f64) * span; ms = ms.max(b.derivative_at(t).length()); }
    ms
}
