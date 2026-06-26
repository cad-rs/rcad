// OCCT-aligned: Angle2D, dir_to_angle (-> Angle(gp_Dir2d)), ClockWiseAngle
//   WireSplitter_1.cxx L768-840 (Angle2D), L621-650 (ClockWiseAngle)
//
//   clock_wise_angle: ✓ matches OCCT formula a1n - ao, epsilon guard
//
//   tolerance_2d: ✓ 1.1x multiplier only for BSpline (OCCT L889-892)
//   Architectural gap: vertex parameter t is provided by caller via domain array
//     OCCT uses BRep_Tool::Parameter(aV, anEdge, myFace) — TopoDS stores
//     per-vertex-per-edge parameter in its shape/vertex/edge data model.
//     rcad DS lacks this per-vertex-on-edge parameter table.  Fix requires
//     adding a `param_on_edge: HashMap<(usize,usize), f64>` to DS or
//     storing the parameter in WireSegment/PaveBlock directly.

use glam::DVec2;
use rcad_kernel::geom::*;
use super::curve_tools::*;

/// OCCT-aligned: Tolerance2D — BOPAlgo_WireSplitter_1.cxx L859-881
///   aTol2D = max(UResolution(aTolV3D), VResolution(aTolV3D), aTolV3D)
///   For BSpline surface: multiplied by 1.1
pub(crate) fn tolerance_2d(vt: f64, surface: &Surface3, v_opt: Option<f64>) -> f64 {
    let u_res = u_resolution(vt, surface, v_opt);
    let v_res = v_resolution(vt, surface);
    let mut t2d = u_res.max(v_res).max(vt);
    if matches!(surface, Surface3::BSpline(_)) {
        t2d *= 1.1;
    }
    t2d
}

/// OCCT-aligned: BRepAdaptor_Surface::UResolution
///   For Cone: tol / radius_at(V) where radius_at(V) = radius + V * tan(half_angle).
///   When v_opt is None, falls back to apex radius (radius at V=0).
fn u_resolution(vt: f64, surface: &Surface3, v_opt: Option<f64>) -> f64 {
    match surface {
        Surface3::Sphere(s) => vt / s.radius.max(1e-15),
        Surface3::Cylinder(c) => vt / c.radius.max(1e-15),
        Surface3::Cone(c) => {
            let r = match v_opt {
                Some(v) => (c.radius + v * c.half_angle_rad.tan()).abs(),
                None => c.radius,
            };
            vt / r.max(1e-15)
        }
        Surface3::Torus(t) => vt / t.major_radius.max(1e-15),
        _ => vt,
    }
}

/// OCCT-aligned: BRepAdaptor_Surface::VResolution
fn v_resolution(vt: f64, surface: &Surface3) -> f64 {
    match surface {
        Surface3::Sphere(s) => vt / s.radius.max(1e-15),
        Surface3::Cylinder(_) => vt,
        Surface3::Cone(_) => vt,
        Surface3::Torus(t) => vt / t.minor_radius.max(1e-15),
        _ => vt,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveGeomType { Line, Circle, Ellipse, Parabola, Hyperbola, Other }

pub fn curve_geom_type(curve: &Curve2d) -> CurveGeomType {
    match curve { Curve2d::Line(_) => CurveGeomType::Line, Curve2d::Circle(_) => CurveGeomType::Circle, Curve2d::Ellipse(_) => CurveGeomType::Ellipse,
    Curve2d::Parabola(_) => CurveGeomType::Parabola,
    Curve2d::Hyperbola(_) => CurveGeomType::Hyperbola, _ => CurveGeomType::Other }
}

pub fn dir_to_angle(dir: DVec2) -> f64 {
    let a = dir.y.atan2(dir.x);
    if a < 0.0 { a + std::f64::consts::TAU } else { a }
}

pub fn angle_2d(curve: &Curve2d, t: f64, domain: [f64; 2], b_is_in: bool, surface: &Surface3, geom_tol: f64, v_opt: Option<f64>) -> Option<f64> {
    let first = domain[0]; let last = domain[1];
    let range = (last - first).abs();
    if range < 1e-15 { return None; }
    let a_tol_2d = 2.0 * tolerance_2d(geom_tol, surface, v_opt);
    let mut dt = curve2d_resolution(curve, a_tol_2d).max(1e-9);
    let typ = curve_geom_type(curve);
    if typ != CurveGeomType::Line {
        let d1 = curve2d_d1(curve, t); let d2 = curve2d_d2(curve, t);
        if d1.length_squared() > 1e-14 {
            let r = curve2d_lprop_curvature(d1, d2, 1e-14);
            if r > 1e-7 {
                let r_curv = 1.0 / r; let cosphi = r_curv / (r_curv + a_tol_2d);
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
