use glam::DVec2;
use rcad_kernel::geom::*;
use super::curve_tools::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveGeomType { Line, Circle, Ellipse, Other }

pub fn curve_geom_type(curve: &Curve2d) -> CurveGeomType {
    match curve { Curve2d::Line(_) => CurveGeomType::Line, Curve2d::Circle(_) => CurveGeomType::Circle, Curve2d::Ellipse(_) => CurveGeomType::Ellipse, _ => CurveGeomType::Other }
}

pub fn dir_to_angle(dir: DVec2) -> f64 {
    let a = dir.y.atan2(dir.x);
    if a < 0.0 { a + std::f64::consts::TAU } else { a }
}

pub fn angle_2d(curve: &Curve2d, t: f64, domain: [f64; 2], b_is_in: bool) -> Option<f64> {
    let first = domain[0]; let last = domain[1];
    let range = (last - first).abs();
    if range < 1e-15 { return None; }
    let tol2d = 1e-5;
    let mut dt = curve2d_resolution(curve, tol2d).max(1e-7);
    let typ = curve_geom_type(curve);
    if typ != CurveGeomType::Line {
        let d1 = curve2d_d1(curve, t); let d2 = curve2d_d2(curve, t);
        if d1.length_squared() > 1e-14 {
            let r = curve2d_lprop_curvature(d1, d2, 1e-14);
            if r > 1e-7 {
                let r_curv = 1.0 / r; let cosphi = r_curv / (r_curv + tol2d);
                if cosphi < 1.0 { dt = dt.max(cosphi.acos()); }
            }
        }
    }
    let max_dt = 0.05 * range;
    let a_tx = if max_dt < 5e-5 { (5e-5_f64).min(range / 2.0) } else { max_dt };
    if dt > a_tx { dt = a_tx; }
    let t1 = if (t - first).abs() < (t - last).abs() { (t + dt).min(last) } else { (t - dt).max(first) };
    let dir = if b_is_in { curve.point_at(t) - curve.point_at(t1) } else { curve.point_at(t1) - curve.point_at(t) };
    if dir.length_squared() < 1e-40 { return None; }
    Some(dir_to_angle(dir))
}

pub fn clock_wise_angle(angle_in: f64, angle_out: f64) -> f64 {
    const TAU: f64 = std::f64::consts::TAU;
    let ai = if angle_in >= TAU { angle_in - TAU } else { angle_in };
    let ao = if angle_out >= TAU { angle_out - TAU } else { angle_out };
    let a1 = ai + std::f64::consts::PI;
    let a1n = if a1 >= TAU { a1 - TAU } else { a1 };
    let mut d = a1n - ao;
    if d <= 0.0 { d += TAU; }
    if d > 0.0 && d <= 1e-14 { d = TAU; }
    d
}
