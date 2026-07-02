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
    let t1 = if (t - first).abs() < (t - last).abs() { t + dt } else { t - dt };
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

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    #[test]
    fn dir_to_angle_positive_x() {
        let a = dir_to_angle(DVec2::new(1.0, 0.0));
        assert!((a - 0.0).abs() < 1e-15);
    }

    #[test]
    fn dir_to_angle_positive_y() {
        let a = dir_to_angle(DVec2::new(0.0, 1.0));
        assert!((a - std::f64::consts::FRAC_PI_2).abs() < 1e-15);
    }

    #[test]
    fn dir_to_angle_negative_x() {
        let a = dir_to_angle(DVec2::new(-1.0, 0.0));
        assert!((a - std::f64::consts::PI).abs() < 1e-15);
    }

    #[test]
    fn dir_to_angle_wraps_negative() {
        let a = dir_to_angle(DVec2::new(-1.0, -1e-15));
        assert!(a >= 0.0 && a < std::f64::consts::TAU);
    }

    #[test]
    fn clock_wise_angle_basic() {
        // AIn=0, AOut=PI/2: A1=PI, dA=PI-PI/2=PI/2
        let d = clock_wise_angle(0.0, std::f64::consts::FRAC_PI_2);
        assert!((d - std::f64::consts::FRAC_PI_2).abs() < 1e-14);
    }

    #[test]
    fn clock_wise_angle_opposite() {
        // AIn=PI, AOut=0: A1=0 (PI+PI=TAU→0), dA=0-0=0 → dA=TAU
        let d = clock_wise_angle(std::f64::consts::PI, 0.0);
        assert!((d - std::f64::consts::TAU).abs() < 1e-14);
    }

    #[test]
    fn clock_wise_angle_wraps_large() {
        let d = clock_wise_angle(0.0, std::f64::consts::TAU);
        assert!((d - std::f64::consts::PI).abs() < 1e-14);
    }

    #[test]
    fn tolerance_2d_plane() {
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let t = tolerance_2d(1e-5, &p, None);
        assert!((t - 1e-5).abs() < 1e-10); // URes=VRes=vt → max=vt
    }

    #[test]
    fn tolerance_2d_sphere_max_is_vt() {
        let s = Surface3::Sphere(SphericalSurface { center: DVec3::ZERO, axis: DVec3::Z, radius: 10.0, ref_dir: DVec3::X });
        let t = tolerance_2d(1e-3, &s, None);
        // URes=VRes=1e-3/10=1e-4, vt=1e-3 → max=1e-3
        assert!((t - 1e-3).abs() < 1e-10);
    }

    #[test]
    fn tolerance_2d_cylinder_u_res_dominates() {
        let c = Surface3::Cylinder(CylindricalSurface { origin: DVec3::ZERO, axis: DVec3::Z, radius: 0.01, ref_dir: DVec3::X });
        let t = tolerance_2d(1e-3, &c, None);
        // URes=1e-3/0.01=0.1, VRes=1e-3, vt=1e-3 → max=0.1
        assert!((t - 0.1).abs() < 1e-10);
    }

    #[test]
    fn tolerance_2d_bspline_scaled() {
        let b = Surface3::BSpline(BSplineSurface {
            degree_u: 3, degree_v: 3,
            knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![vec![DVec3::ZERO; 4]; 4],
            weights: vec![vec![1.0; 4]; 4],
        });
        let t = tolerance_2d(1e-5, &b, None);
        // Base = max(vt, vt, vt) * 1.1 = 1.1e-5
        assert!((t - 1.1e-5).abs() < 1e-10);
    }

    #[test]
    fn u_resolution_cone_with_v() {
        let c = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO, axis: DVec3::Z,
            radius: 1.0, half_angle_rad: 0.5,
        });
        // V=2: r = 1 + 2*tan(0.5) ≈ 1 + 2*0.546 = 2.092, URes=1e-3/2.092≈4.78e-4
        let u = u_resolution(1e-3, &c, Some(2.0));
        assert!(u > 0.0 && u < 1e-3);
    }

    #[test]
    fn u_resolution_cone_without_v() {
        let c = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO, axis: DVec3::Z,
            radius: 2.0, half_angle_rad: 0.3,
        });
        let u = u_resolution(1e-3, &c, None);
        assert!((u - 1e-3 / 2.0).abs() < 1e-10);
    }

    #[test]
    fn curve_geom_type_classification() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let c = Curve2d::Circle(Circle2d { center: DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: 1.0  });
        let e = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X,
            major_radius: 2.0, minor_radius: 1.0,
        });
        assert_eq!(curve_geom_type(&l), CurveGeomType::Line);
        assert_eq!(curve_geom_type(&c), CurveGeomType::Circle);
        assert_eq!(curve_geom_type(&e), CurveGeomType::Ellipse);
    }

    #[test]
    fn angle_2d_returns_none_for_zero_range() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let a = angle_2d(&l, 0.5, [0.0, 0.0], false, &p, 1e-5, None);
        assert!(a.is_none());
    }

    #[test]
    fn angle_2d_line_returns_some() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::new(1.0, 0.0) });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let a = angle_2d(&l, 0.5, [0.0, 1.0], false, &p, 1e-5, None);
        assert!(a.is_some());
    }

    #[test]
    fn angle_2d_is_in_flips_direction() {
        let l = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let a_out = angle_2d(&l, 0.5, [0.0, 1.0], false, &p, 1e-5, None).unwrap();
        let a_in = angle_2d(&l, 0.5, [0.0, 1.0], true, &p, 1e-5, None).unwrap();
        let diff = (a_out - a_in).abs();
        assert!((diff - std::f64::consts::PI).abs() < 0.01);
    }

    #[test]
    fn angle_2d_circle_non_line_triggers_curvature() {
        let c = Curve2d::Circle(Circle2d { center: DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: 10.0  });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        // Circle at t=0: d1=(0,10), d2=(-10,0), curvature=1/10
        let a = angle_2d(&c, 0.0, [0.0, std::f64::consts::TAU], false, &p, 1e-5, None);
        assert!(a.is_some());
        // At t=0, direction should be upward (0, positive) → angle ~ PI/2
        assert!((a.unwrap() - std::f64::consts::FRAC_PI_2).abs() < 0.1);
    }
}
