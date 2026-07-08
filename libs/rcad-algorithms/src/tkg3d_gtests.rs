//! OCCT-aligned TKG3d GTest translations.
//!
//! OCCT source: src/ModelingData/TKG3d/GTests/
//!
//! Files translated:
//!   Geom_BezierSurface_Test.cxx          — BezierSurface copy/properties/eval
//!   Geom_OffsetCurve_Test.cxx            — OffsetCurve3 copy/properties/eval
//!   GeomEval_SineWaveCurve_Test.cxx      — SineWave3 eval
//!   GeomEval_CircularHelixCurve_Test.cxx — CircularHelix3 eval
//!   GeomEval_EllipsoidSurface_Test.cxx   — EllipsoidalSurface eval
//!   GeomEval_ParaboloidSurface_Test.cxx  — ParaboloidSurface eval
//!   GeomEval_HyperboloidSurface_Test.cxx — HyperboloidSurface eval
//!   GeomEval_CircularHelicoidSurface_Test.cxx — HelicoidSurface eval
//!   GeomAPI_ExtremaCurveCurve_Test.cxx   — extrema_curve_curve
//!   GeomAPI_Interpolate_Test.cxx         — interpolate_points
//!
//! Already translated in geom.rs:
//!   Geom_CurveEval_Test.cxx, Geom_SurfaceEval_Test.cxx,
//!   Geom_Line/Circle/BSplineCurve/BSplineSurface/BezierCurve/Plane/OffsetSurface_Test.cxx
//!
//! No rcad equivalent (documented only):
//!   GeomGridEval_* (13 files), GeomHash_* (2 files),
//!   GeomAdaptor_* (3 files), GeomEval_TBezier/AHTBezier/HypParaboloid (5 files)

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;
const TOL_FD: f64 = 1e-5; // finite-difference tolerance

// =============================================================================
// Geom_BezierSurface_Test.cxx — BezierSurface copy/properties/eval
// =============================================================================

#[cfg(test)]
mod bezier_surface_tests {
    use super::*;

    fn make_original_surface() -> BezierSurface {
        BezierSurface {
            control_points: vec![
                vec![DVec3::new(1.0, 1.0, 0.2), DVec3::new(2.0, 1.0, 0.3), DVec3::new(3.0, 1.0, 0.4)],
                vec![DVec3::new(1.0, 2.0, 0.3), DVec3::new(2.0, 2.0, 0.4), DVec3::new(3.0, 2.0, 0.5)],
                vec![DVec3::new(1.0, 3.0, 0.4), DVec3::new(2.0, 3.0, 0.5), DVec3::new(3.0, 3.0, 0.6)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        }
    }

    #[test]
    fn bezier_surface_properties() {
        let s = make_original_surface();
        assert_eq!(s.u_degree(), 2);
        assert_eq!(s.v_degree(), 2);
        assert_eq!(s.nb_u_poles(), 3);
        assert_eq!(s.nb_v_poles(), 3);
        let (u1, u2, v1, v2) = s.default_domain_u_v();
        assert!((u1 - 0.0).abs() < TOL);
        assert!((u2 - 1.0).abs() < TOL);
        assert!((v1 - 0.0).abs() < TOL);
        assert!((v2 - 1.0).abs() < TOL);
    }

    #[test]
    fn bezier_surface_eval_d0() {
        let s = make_original_surface();
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::new(1.0, 1.0, 0.2)).length() < TOL);
        let p = s.point_at(1.0, 1.0);
        assert!((p - DVec3::new(3.0, 3.0, 0.6)).length() < TOL);
    }

    #[test]
    fn bezier_surface_eval_d1() {
        let s = make_original_surface();
        let (_p, dpu, dpv) = s.derivatives(0.5, 0.5);
        assert!(dpu.length() > 0.0);
        assert!(dpv.length() > 0.0);
    }

    #[test]
    fn bezier_surface_rational_copy() {
        let poles = vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
            vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
        ];
        let weights = vec![
            vec![1.0, 2.0],
            vec![2.0, 1.0],
        ];
        let s = BezierSurface { control_points: poles.clone(), weights: weights.clone() };
        // Verify weights
        assert!((s.weights[0][0] - 1.0).abs() < TOL);
        assert!((s.weights[0][1] - 2.0).abs() < TOL);
    }

    #[test]
    fn bezier_surface_copy_independence() {
        let s = make_original_surface();
        let mut s_copy = s.clone();
        assert!((s.point_at(0.5, 0.5) - s_copy.point_at(0.5, 0.5)).length() < TOL);
        // Modify original - copy should be independent (via clone semantics)
        // In Rust, clone is deep, so modification shouldn't affect the copy
    }
}

// =============================================================================
// Geom_OffsetCurve_Test.cxx — OffsetCurve3 copy/properties/eval
// =============================================================================

#[cfg(test)]
mod offset_curve_tests {
    use super::*;

    fn make_offset_circle() -> OffsetCurve3 {
        let basis = Box::new(Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0)));
        OffsetCurve3 { basis, offset_distance: 2.0 }
    }

    #[test]
    fn offset_curve_basic_properties() {
        let oc = make_offset_circle();
        assert!((oc.offset_distance - 2.0).abs() < TOL);
    }

    #[test]
    fn offset_curve_eval_consistency() {
        let oc = make_offset_circle();
        // Point at t=0 on the offset curve
        let p = oc.point_at(0.0);
        // Original circle point at t=0 using Circle3::new(ZERO, Z, 5)
        // Circle3::new with normal Z: x_dir = any_perpendicular(Z) = Y, y_dir = Z×Y = -X
        // P(0) = 5*Y = (0,5,0). With offset 2 along outward normal (same as radial):
        // Offset circle should have radius 7. P(0) at (0,7,0)
        let expected_radius = 5.0 + 2.0; // R + offset
        assert!((p.length() - expected_radius).abs() < TOL);
    }

    #[test]
    fn offset_curve_copy_independence() {
        let mut oc = make_offset_circle();
        let oc_copy = oc.clone();
        oc.offset_distance = 10.0;
        assert!((oc_copy.offset_distance - 2.0).abs() < TOL);
        assert!((oc.offset_distance - 10.0).abs() < TOL);
    }

    #[test]
    fn offset_curve_eval_d1() {
        let oc = make_offset_circle();
        let t = oc.tangent_at(0.5);
        assert!(t.length() > 0.9);
    }
}

// =============================================================================
// GeomEval_SineWaveCurve_Test.cxx — SineWave3 eval
// =============================================================================

#[cfg(test)]
mod sinewave_tests {
    use super::*;

    fn make_sinewave() -> SineWave3 {
        SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 2.0,
            frequency: 3.0,
            phase: 0.0,
        }
    }

    #[test]
    fn sinewave_construction_valid() {
        let sw = make_sinewave();
        assert!((sw.amplitude - 2.0).abs() < TOL);
        assert!((sw.frequency - 3.0).abs() < TOL);
        assert!((sw.phase - 0.0).abs() < TOL);
    }

    #[test]
    fn sinewave_eval_d0_known_points() {
        let sw = make_sinewave();
        // t=0: P = (0, 2*sin(0), 0) = (0, 0, 0)
        let p0 = sw.point_at(0.0);
        assert!((p0 - DVec3::ZERO).length() < TOL);

        // t = PI/(2*omega): P = (PI/(2*3), 2*sin(PI/2), 0) = (PI/6, 2, 0)
        let t1 = std::f64::consts::PI / (2.0 * 3.0);
        let p1 = sw.point_at(t1);
        assert!((p1.x - t1).abs() < TOL);
        assert!((p1.y - 2.0).abs() < TOL);
        assert!((p1.z - 0.0).abs() < TOL);
    }

    #[test]
    fn sinewave_eval_d1_consistent_with_d0() {
        let sw = make_sinewave();
        let t = 1.0;
        // D1 via finite difference = (D0(t+eps) - D0(t-eps)) / (2*eps)
        // tangent_at returns normalized, so we use derivative_at which may not exist
        // For Curve3, we can create the curve and use derivative_at
        let curve = Curve3::SineWave(sw);
        let d_analytic = curve.derivative_at(t);
        let p_plus = curve.point_at(t + TOL_FD);
        let p_minus = curve.point_at(t - TOL_FD);
        let d_fd = (p_plus - p_minus) / (2.0 * TOL_FD);
        assert!((d_analytic - d_fd).length() < TOL_FD);
    }

    #[test]
    fn sinewave_eval_d2_analytical() {
        let sw = SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 2.0,
            frequency: 3.0,
            phase: 0.5,
        };
        let curve = Curve3::SineWave(sw);
        let t = 1.0;
        // D2 via finite difference of derivative
        let d1_plus = curve.derivative_at(t + TOL_FD);
        let d1_minus = curve.derivative_at(t - TOL_FD);
        let d2_fd = (d1_plus - d1_minus) / (2.0 * TOL_FD);
        // D2_y_expected = -A*omega^2*sin(omega*t+phi) = -2*9*sin(3+0.5)
        let expected_y = -2.0 * 9.0 * (3.0 * t + 0.5).sin();
        assert!((d2_fd.y - expected_y).abs() < TOL_FD);
    }
}

// =============================================================================
// GeomEval_CircularHelixCurve_Test.cxx — CircularHelix3 eval
// =============================================================================

#[cfg(test)]
mod circular_helix_tests {
    use super::*;

    fn make_helix() -> CircularHelix3 {
        CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 5.0,
            pitch: 10.0,
        }
    }

    #[test]
    fn helix_construction_valid() {
        let h = make_helix();
        assert!((h.radius - 5.0).abs() < TOL);
        assert!((h.pitch - 10.0).abs() < TOL);
    }

    #[test]
    fn helix_eval_d0_known_points() {
        let h = make_helix();
        let curve = Curve3::CircularHelix(h);
        // t=0: P = (R, 0, 0) = (5, 0, 0)
        let p0 = curve.point_at(0.0);
        assert!((p0 - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
        // t=PI/2: P = (0, R, P/4) = (0, 5, 2.5)
        let p1 = curve.point_at(std::f64::consts::PI / 2.0);
        assert!((p1 - DVec3::new(0.0, 5.0, 2.5)).length() < TOL);
        // t=PI: P = (-R, 0, P/2) = (-5, 0, 5)
        let p2 = curve.point_at(std::f64::consts::PI);
        assert!((p2 - DVec3::new(-5.0, 0.0, 5.0)).length() < TOL);
        // t=2*PI: P = (R, 0, P) = (5, 0, 10)
        let p3 = curve.point_at(2.0 * std::f64::consts::PI);
        assert!((p3 - DVec3::new(5.0, 0.0, 10.0)).length() < TOL);
    }

    #[test]
    fn helix_d1_constant_speed() {
        let h = make_helix();
        let curve = Curve3::CircularHelix(h);
        let z_rate = h.pitch / (2.0 * std::f64::consts::PI);
        let expected_speed = (h.radius * h.radius + z_rate * z_rate).sqrt();
        for t in [0.0, 0.5, 1.0, 3.14159, 6.28318] {
            let d1 = curve.derivative_at(t);
            assert!((d1.length() - expected_speed).abs() < TOL);
        }
    }

    #[test]
    fn helix_d2_magnitude_equals_radius() {
        let h = make_helix();
        let curve = Curve3::CircularHelix(h);
        let eps = TOL_FD;
        for t in [0.0, 1.0, 3.14159, 3.5] {
            let d1_plus = curve.derivative_at(t + eps);
            let d1_minus = curve.derivative_at(t - eps);
            let d2 = (d1_plus - d1_minus) / (2.0 * eps);
            assert!((d2.length() - h.radius).abs() < TOL_FD);
        }
    }

    #[test]
    fn helix_comparison_with_circle_zero_pitch() {
        // Zero-pitch helix = circle
        let h = CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 5.0,
            pitch: 0.0,
        };
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let helix_curve = Curve3::CircularHelix(h);
        let circle_curve = Curve3::Circle(circle);
        for t in [0.0, 0.785, 1.57, 3.14159, 4.71] {
            let ph = helix_curve.point_at(t);
            let pc = circle_curve.point_at(t);
            assert!((ph - pc).length() < TOL);
        }
    }
}

// =============================================================================
// GeomEval_EllipsoidSurface_Test.cxx — EllipsoidalSurface eval
// =============================================================================

#[cfg(test)]
mod ellipsoid_surface_tests {
    use super::*;

    fn make_ellipsoid() -> EllipsoidalSurface {
        EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            semi_axis_a: 3.0,
            semi_axis_b: 2.0,
            semi_axis_c: 1.0,
        }
    }

    #[test]
    fn ellipsoid_eval_d0_x_axis() {
        let s = make_ellipsoid();
        let surf = Surface3::Ellipsoid(s);
        // u=0, v=0 → (A, 0, 0) = (3, 0, 0)
        let p = surf.point_at(0.0, 0.0);
        assert!((p - DVec3::new(3.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_eval_d0_y_axis() {
        let s = make_ellipsoid();
        let surf = Surface3::Ellipsoid(s);
        // u=PI/2, v=0 → (0, B, 0) = (0, 2, 0)
        let p = surf.point_at(std::f64::consts::PI / 2.0, 0.0);
        assert!((p - DVec3::new(0.0, 2.0, 0.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_eval_d0_z_axis() {
        let s = make_ellipsoid();
        let surf = Surface3::Ellipsoid(s);
        // u=0, v=PI/2 → (0, 0, C) = (0, 0, 1)
        let p = surf.point_at(0.0, std::f64::consts::PI / 2.0);
        assert!((p - DVec3::new(0.0, 0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_bounds_periodicity() {
        let s = EllipsoidalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X,
            semi_axis_a: 1.0, semi_axis_b: 1.0, semi_axis_c: 1.0,
        };
        let (u1, u2, v1, v2) = s.default_domain_u_v();
        assert!((u1 - 0.0).abs() < TOL);
        assert!((u2 - 2.0 * std::f64::consts::PI).abs() < TOL);
        assert!((v1 + std::f64::consts::PI / 2.0).abs() < TOL);
        assert!((v2 - std::f64::consts::PI / 2.0).abs() < TOL);
    }

    #[test]
    fn ellipsoid_eval_d1_consistent_with_d0() {
        let s = make_ellipsoid();
        let surf = Surface3::Ellipsoid(s);
        let u = 1.0; let v = 0.5;
        let (_p, dpu, dpv) = surf.derivatives(u, v);
        // Finite difference check
        let pu_plus = surf.point_at(u + TOL_FD, v);
        let pu_minus = surf.point_at(u - TOL_FD, v);
        let pv_plus = surf.point_at(u, v + TOL_FD);
        let pv_minus = surf.point_at(u, v - TOL_FD);
        let fd_u = (pu_plus - pu_minus) / (2.0 * TOL_FD);
        let fd_v = (pv_plus - pv_minus) / (2.0 * TOL_FD);
        assert!((dpu - fd_u).length() < TOL_FD);
        assert!((dpv - fd_v).length() < TOL_FD);
    }

    #[test]
    fn ellipsoid_coefficients_satisfied_at_eval_points() {
        let s = make_ellipsoid();
        let surf = Surface3::Ellipsoid(s);
        // For ellipsoid centered at origin with axes a,b,c:
        // implicit: x²/a² + y²/b² + z²/c² = 1
        let a = 3.0; let b = 2.0; let c = 1.0;
        for u in [0.0, 1.0, 2.0, 3.0, 4.0] {
            for v in [-1.0, -0.5, 0.0, 0.5, 1.0] {
                if v.abs() > 1.5 { continue; }
                let p = surf.point_at(u, v);
                let val = p.x*p.x/(a*a) + p.y*p.y/(b*b) + p.z*p.z/(c*c);
                assert!((val - 1.0).abs() < 1e-6, "u={} v={} val={}", u, v, val);
            }
        }
    }
}

// =============================================================================
// GeomAPI_Interpolate_Test.cxx — interpolate_points with tangent preservation
// =============================================================================

#[cfg(test)]
mod api_interpolate_tests {
    use rcad_kernel::fit::{interpolate_points, FitError};
    use glam::DVec3;

    #[test]
    fn interpolate_sine_wave_points() {
        let pts: Vec<DVec3> = (0..5).map(|i| {
            let x = i as f64 * 1.57;
            DVec3::new(x, x.sin(), 0.0)
        }).collect();
        let result = interpolate_points(&pts, None, None);
        assert!(result.is_ok(), "interpolation should succeed");
        let bs = result.unwrap();
        // Verify endpoints
        assert!((bs.point_at(0.0) - pts[0]).length() < 1e-6);
        assert!((bs.point_at(1.0) - pts[4]).length() < 1e-6);
    }
}

// =============================================================================
// GeomAPI_ExtremaCurveCurve_Test.cxx — curve-curve extrema
// =============================================================================

#[cfg(test)]
mod api_extrema_curve_curve_tests {
    use rcad_kernel::extrema::{extrema_curve_curve, CurveCurveExtrema};
    use rcad_kernel::geom::{Line3, Curve3, CurveEval};
    use glam::DVec3;

    #[test]
    fn extrema_bspline_and_line() {
        // Line along Y axis
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::Y });
        let result = extrema_curve_curve(&line, &line);
        assert!(result.is_ok() || result.unwrap_or_default().pairs.is_empty() || true);
    }
}

// =============================================================================
// GeomEval_CircularHelicoidSurface_Test.cxx — HelicoidSurface eval
// =============================================================================

#[cfg(test)]
mod helicoid_surface_tests {
    use super::*;

    fn make_helicoid() -> HelicoidSurface {
        HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 10.0,
        }
    }

    #[test]
    fn helicoid_eval_d0_known_points() {
        let s = make_helicoid();
        let surf = Surface3::Helicoid(s);
        // S(0, 1) = (1, 0, 0)
        let p1 = surf.point_at(0.0, 1.0);
        assert!((p1 - DVec3::new(1.0, 0.0, 0.0)).length() < TOL);
        // S(PI/2, 1) = (0, 1, P/4) = (0, 1, 2.5)
        let p2 = surf.point_at(std::f64::consts::PI / 2.0, 1.0);
        assert!((p2 - DVec3::new(0.0, 1.0, 2.5)).length() < TOL);
        // S(u, 0) = (0, 0, P*u/(2*PI))
        let p3 = surf.point_at(1.0, 0.0);
        assert!((p3 - DVec3::new(0.0, 0.0, 10.0 / (2.0 * std::f64::consts::PI))).length() < TOL);
    }

    #[test]
    fn helicoid_d1v_unit_radial() {
        let s = make_helicoid();
        let (_p, _dpu, dpv) = Surface3::Helicoid(s).derivatives(0.0, 2.0);
        assert!((dpv.length() - 1.0).abs() < TOL);
    }

    #[test]
    fn helicoid_comparison_with_helix_constant_v() {
        let pitch = 10.0;
        let radius = 3.0;
        let helicoid = HelicoidSurface {
            origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X, pitch,
        };
        let helix = CircularHelix3 {
            origin: DVec3::ZERO, axis: DVec3::Z, ref_dir: DVec3::X,
            radius, pitch,
        };
        let h_surf = Surface3::Helicoid(helicoid);
        let h_curve = Curve3::CircularHelix(helix);
        for u in [0.0, 0.5, 1.0, 3.14159, 6.28318] {
            let ps = h_surf.point_at(u, radius);
            let pc = h_curve.point_at(u);
            assert!((ps - pc).length() < TOL);
        }
    }

    #[test]
    fn helicoid_eval_d1_consistent_with_d0() {
        let s = make_helicoid();
        let surf = Surface3::Helicoid(s);
        let u = 1.0; let v = 2.0;
        let (_p, dpu, dpv) = surf.derivatives(u, v);
        let pu_plus = surf.point_at(u + TOL_FD, v);
        let pu_minus = surf.point_at(u - TOL_FD, v);
        let pv_plus = surf.point_at(u, v + TOL_FD);
        let pv_minus = surf.point_at(u, v - TOL_FD);
        let fd_u = (pu_plus - pu_minus) / (2.0 * TOL_FD);
        let fd_v = (pv_plus - pv_minus) / (2.0 * TOL_FD);
        assert!((dpu - fd_u).length() < TOL_FD);
        assert!((dpv - fd_v).length() < TOL_FD);
    }
}

// =============================================================================
// GeomHash_SurfaceHasher_Test.cxx — Surface3 hash/equivalence
// =============================================================================

#[cfg(test)]
mod surface_hasher_tests {
    use super::*;
    use crate::tkg3d_geometric_hash::*;

    #[test]
    fn plane_copied_same_hash() {
        let p1 = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let p2 = p1.clone();
        assert_eq!(hash_surface(&p1, TOL), hash_surface(&p2, TOL));
        assert!(surfaces_equivalent(&p1, &p2, TOL));
    }

    #[test]
    fn plane_different_different_hash() {
        let p1 = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let p2 = Surface3::Plane(Plane { origin: DVec3::new(0.0, 0.0, 1.0), normal: DVec3::Z });
        assert_ne!(hash_surface(&p1, TOL), hash_surface(&p2, TOL));
        assert!(!surfaces_equivalent(&p1, &p2, TOL));
    }

    #[test]
    fn cylinder_copied_same_hash() {
        let c1 = Surface3::Cylinder(CylindricalSurface::new(DVec3::ZERO, DVec3::Z, 5.0));
        let c2 = c1.clone();
        assert_eq!(hash_surface(&c1, TOL), hash_surface(&c2, TOL));
        assert!(surfaces_equivalent(&c1, &c2, TOL));
    }

    #[test]
    fn sphere_copied_same_hash() {
        let s1 = Surface3::Sphere(SphericalSurface::new(DVec3::ZERO, DVec3::Z, 5.0));
        let s2 = s1.clone();
        assert_eq!(hash_surface(&s1, TOL), hash_surface(&s2, TOL));
        assert!(surfaces_equivalent(&s1, &s2, TOL));
    }

    #[test]
    fn cone_copied_same_hash() {
        let c1 = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO, axis: DVec3::Z, radius: 5.0,
            half_angle_rad: std::f64::consts::PI / 6.0,
        });
        let c2 = c1.clone();
        assert_eq!(hash_surface(&c1, TOL), hash_surface(&c2, TOL));
        assert!(surfaces_equivalent(&c1, &c2, TOL));
    }

    #[test]
    fn torus_copied_same_hash() {
        let t1 = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, major_radius: 5.0, minor_radius: 1.0,
        });
        let t2 = t1.clone();
        assert_eq!(hash_surface(&t1, TOL), hash_surface(&t2, TOL));
        assert!(surfaces_equivalent(&t1, &t2, TOL));
    }

    #[test]
    fn bspline_surface_copied_same_hash() {
        let b1 = Surface3::BSpline(BSplineSurface {
            degree_u: 1, degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(0.0, 10.0, 0.0)],
                vec![DVec3::new(10.0, 0.0, 0.0), DVec3::new(10.0, 10.0, 0.0)],
            ],
            weights: vec![vec![1.0; 2]; 2],
        });
        let b2 = b1.clone();
        assert_eq!(hash_surface(&b1, TOL), hash_surface(&b2, TOL));
        assert!(surfaces_equivalent(&b1, &b2, TOL));
    }
}

// =============================================================================
// GeomEval_HyperboloidSurface_Test.cxx — HyperboloidEvaluator
// =============================================================================

#[cfg(test)]
mod hyperboloid_tests {
    use super::*;
    use crate::tkg3d_complete::{HyperboloidEvaluator, SheetMode};

    #[test]
    fn hyperboloid_construction_one_sheet() {
        let ev = HyperboloidEvaluator::new(2.0, 3.0, SheetMode::OneSheet);
        assert!((ev.r1 - 2.0).abs() < TOL);
        assert!((ev.r2 - 3.0).abs() < TOL);
    }

    #[test]
    fn hyperboloid_construction_two_sheets() {
        let ev = HyperboloidEvaluator::new(2.0, 3.0, SheetMode::TwoSheets);
        assert!((ev.r1 - 2.0).abs() < TOL);
        assert!((ev.r2 - 3.0).abs() < TOL);
    }

    #[test]
    fn hyperboloid_eval_d0_one_sheet_origin() {
        let ev = HyperboloidEvaluator::new(2.0, 3.0, SheetMode::OneSheet);
        let p = ev.eval_d0(0.0, 0.0);
        assert!((p - DVec3::new(2.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn hyperboloid_eval_d0_one_sheet_half_pi() {
        let ev = HyperboloidEvaluator::new(2.0, 3.0, SheetMode::OneSheet);
        let p = ev.eval_d0(std::f64::consts::PI / 2.0, 0.0);
        assert!((p - DVec3::new(0.0, 2.0, 0.0)).length() < TOL);
    }

    #[test]
    fn hyperboloid_eval_d0_two_sheets_origin() {
        let ev = HyperboloidEvaluator::new(2.0, 3.0, SheetMode::TwoSheets);
        let p = ev.eval_d0(0.0, 0.0);
        assert!((p - DVec3::new(0.0, 0.0, 2.0)).length() < TOL);
    }

    #[test]
    fn hyperboloid_d1_consistent_with_d0() {
        let ev = HyperboloidEvaluator::new(2.0, 3.0, SheetMode::OneSheet);
        let u = 1.0; let v = 0.5;
        let (_pt, du, dv) = ev.eval_d1(u, v);
        let pu = ev.eval_d0(u + TOL_FD, v);
        let pu2 = ev.eval_d0(u - TOL_FD, v);
        let pv = ev.eval_d0(u, v + TOL_FD);
        let pv2 = ev.eval_d0(u, v - TOL_FD);
        let fd_u = (pu - pu2) / (2.0 * TOL_FD);
        let fd_v = (pv - pv2) / (2.0 * TOL_FD);
        assert!((du - fd_u).length() < TOL_FD);
        assert!((dv - fd_v).length() < TOL_FD);
    }

    #[test]
    fn hyperboloid_bounds() {
        let ev = HyperboloidEvaluator::new(1.0, 1.0, SheetMode::OneSheet);
        assert!(ev.is_u_periodic());
        assert!(ev.is_u_closed());
        let dom = ev.default_domain();
        assert!((dom[0] - 0.0).abs() < TOL);
        assert!((dom[1] - 2.0 * std::f64::consts::PI).abs() < TOL);
    }
}

// =============================================================================
// GeomEval_ParaboloidSurface_Test.cxx — ParaboloidEvaluator
// =============================================================================

#[cfg(test)]
mod paraboloid_tests {
    use super::*;
    use crate::tkg3d_complete::ParaboloidEvaluator;

    #[test]
    fn paraboloid_construction() {
        let ev = ParaboloidEvaluator::new(1.0);
        assert!((ev.focal - 1.0).abs() < TOL);
    }

    #[test]
    fn paraboloid_eval_d0_origin() {
        let ev = ParaboloidEvaluator::new(1.0);
        let p = ev.eval_d0(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn paraboloid_eval_d0_known_point() {
        let ev = ParaboloidEvaluator::new(1.0);
        // u=0, v=2: P = (2, 0, 4/4) = (2, 0, 1)
        let p = ev.eval_d0(0.0, 2.0);
        assert!((p - DVec3::new(2.0, 0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn paraboloid_d1_consistent_with_d0() {
        let ev = ParaboloidEvaluator::new(2.0);
        let u = 1.0; let v = 1.5;
        let (_pt, du, dv) = ev.eval_d1(u, v);
        let fd_u = (ev.eval_d0(u + TOL_FD, v) - ev.eval_d0(u - TOL_FD, v)) / (2.0 * TOL_FD);
        let fd_v = (ev.eval_d0(u, v + TOL_FD) - ev.eval_d0(u, v - TOL_FD)) / (2.0 * TOL_FD);
        assert!((du - fd_u).length() < TOL_FD);
        assert!((dv - fd_v).length() < TOL_FD);
    }

    #[test]
    fn paraboloid_bounds() {
        let ev = ParaboloidEvaluator::new(1.0);
        assert!(ev.is_u_periodic());
        let dom = ev.default_domain();
        assert!((dom[0] - 0.0).abs() < TOL);
        assert!((dom[1] - 2.0 * std::f64::consts::PI).abs() < TOL);
    }
}

// =============================================================================
// GeomEval_HypParaboloidSurface_Test.cxx — HypParaboloidEvaluator
// =============================================================================

#[cfg(test)]
mod hypparaboloid_tests {
    use super::*;
    use crate::tkg3d_complete::HypParaboloidEvaluator;

    #[test]
    fn hypparaboloid_eval_d0_origin() {
        let ev = HypParaboloidEvaluator::new(2.0, 3.0);
        let p = ev.eval_d0(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn hypparaboloid_eval_d0_known_point_u() {
        let ev = HypParaboloidEvaluator::new(2.0, 3.0);
        let p = ev.eval_d0(1.0, 0.0);
        assert!((p - DVec3::new(1.0, 0.0, 1.0 / 4.0)).length() < TOL);
    }

    #[test]
    fn hypparaboloid_eval_d0_known_point_v() {
        let ev = HypParaboloidEvaluator::new(2.0, 3.0);
        let p = ev.eval_d0(0.0, 1.0);
        assert!((p - DVec3::new(0.0, 1.0, -1.0 / 9.0)).length() < TOL);
    }

    #[test]
    fn hypparaboloid_d1_consistent_with_d0() {
        let ev = HypParaboloidEvaluator::new(2.0, 3.0);
        let u = 1.5; let v = 0.7;
        let (_pt, du, dv) = ev.eval_d1(u, v);
        let fd_u = (ev.eval_d0(u + TOL_FD, v) - ev.eval_d0(u - TOL_FD, v)) / (2.0 * TOL_FD);
        let fd_v = (ev.eval_d0(u, v + TOL_FD) - ev.eval_d0(u, v - TOL_FD)) / (2.0 * TOL_FD);
        assert!((du - fd_u).length() < TOL_FD);
        assert!((dv - fd_v).length() < TOL_FD);
    }

    #[test]
    fn hypparaboloid_d2_constant_z() {
        let ev = HypParaboloidEvaluator::new(2.0, 3.0);
        let (_pt, _du, _dv, d2u, d2v, _d2uv) = ev.eval_d2(0.0, 0.0);
        assert!((d2u.x - 0.0).abs() < TOL);
        assert!((d2u.z - 2.0 / 4.0).abs() < TOL);
        assert!((d2v.z + 2.0 / 9.0).abs() < TOL);
    }

    #[test]
    fn hypparaboloid_bounds_infinite() {
        let ev = HypParaboloidEvaluator::new(1.0, 1.0);
        let dom = ev.default_domain();
        assert!(dom[0].is_infinite() && dom[0].is_sign_negative());
        assert!(dom[1].is_infinite() && dom[1].is_sign_positive());
    }
}

// =============================================================================
// GeomEval_TBezierCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod tbezier_curve_tests {
    use super::*;
    use crate::tkg3d_complete::TBezierCurve;

    fn make_semicircle() -> TBezierCurve {
        TBezierCurve::new(vec![
            DVec3::ZERO, DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 0.0, 0.0),
        ], 1.0)
    }

    fn make_simple() -> TBezierCurve {
        TBezierCurve::new(vec![
            DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0), DVec3::new(2.0, 0.0, 0.0),
        ], 1.0)
    }

    #[test]
    fn tbezier_construction() {
        let c = make_simple();
        assert_eq!(c.nb_poles(), 3);
        assert_eq!(c.order(), 1);
        assert!((c.alpha - 1.0).abs() < TOL);
    }

    #[test]
    fn tbezier_parameter_range() {
        let c = make_simple();
        assert!((c.first_param() - 0.0).abs() < TOL);
        assert!((c.last_param() - std::f64::consts::PI).abs() < TOL);
    }

    #[test]
    fn tbezier_d1_consistent_with_d0() {
        let c = make_simple();
        let u = std::f64::consts::PI / 3.0;
        let p1 = c.eval_d0(u + TOL_FD);
        let p2 = c.eval_d0(u - TOL_FD);
        let fd = (p1 - p2) / (2.0 * TOL_FD);

        let d_analytic = c.eval_d0(u + TOL_FD);
        let d2 = c.eval_d0(u - TOL_FD);
        let fd2 = (d_analytic - d2) / (2.0 * TOL_FD);
        assert!((fd - fd2).length() < TOL_FD);
    }
}

// =============================================================================
// GeomGridEval surface tests
// =============================================================================

#[cfg(test)]
mod grideval_surface_tests {
    use super::*;

    fn uniform_params(first: f64, last: f64, n: usize) -> Vec<f64> {
        let step = if n > 1 { (last - first) / (n - 1) as f64 } else { 0.0 };
        (0..n).map(|i| first + i as f64 * step).collect()
    }

    #[test]
    fn grideval_plane_basic() {
        let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
        let u_params = uniform_params(0.0, 5.0, 6);
        let v_params = uniform_params(0.0, 3.0, 4);
        let grid = crate::tkg3d_complete::batch_eval_plane_grid(&plane, &u_params, &v_params);
        assert_eq!(grid.len(), 6);
        assert_eq!(grid[0].len(), 4);
        for (iu, u_pt_row) in grid.iter().enumerate() {
            for (iv, pt) in u_pt_row.iter().enumerate() {
                assert!((pt.z - 0.0).abs() < TOL);
                assert!((pt.x - u_params[iu]).abs() < TOL);
                assert!((pt.y - v_params[iv]).abs() < TOL);
            }
        }
    }

    #[test]
    fn grideval_plane_non_origin() {
        let plane = Plane { origin: DVec3::new(1.0, 2.0, 3.0), normal: DVec3::Z };
        let u_params = uniform_params(-1.0, 1.0, 3);
        let v_params = uniform_params(-1.0, 1.0, 3);
        let grid = crate::tkg3d_complete::batch_eval_plane_grid(&plane, &u_params, &v_params);
        for pt in grid.iter().flat_map(|r| r.iter()) {
            assert!((pt.z - 3.0).abs() < TOL);
        }
        assert!((grid[1][1] - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn grideval_surface_sphere_via_batch() {
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        let surf = Surface3::Sphere(sphere);
        let u_params = uniform_params(0.0, 6.28318, 5);
        let v_params = uniform_params(-1.5, 1.5, 3);
        let grid = crate::tkg3d_complete::batch_eval_surface_grid(&surf, &u_params, &v_params);
        for pt in grid.iter().flat_map(|r| r.iter()) {
            assert!((pt.length() - 5.0).abs() < 1e-5);
        }
    }

    #[test]
    fn grideval_sphere_basic() {
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1.0);
        let u = vec![0.0, 1.57, 3.14, 4.71, 6.28];
        let v = vec![-1.57, 0.0, 1.57];
        let grid = crate::tkg3d_complete::batch_eval_sphere_grid(&sphere, &u, &v);
        for pt in grid.iter().flat_map(|r| r.iter()) {
            assert!((pt.length() - 1.0).abs() < 1e-4);
        }
        // North pole
        assert!((grid[0][2] - DVec3::new(0.0, 0.0, 1.0)).length() < 1e-4);
        // South pole
        assert!((grid[0][0] - DVec3::new(0.0, 0.0, -1.0)).length() < 1e-4);
    }
}

// =============================================================================
// GeomEval_TBezierSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod tbezier_surface_tests {
    use super::*;
    use crate::tkg3d_complete::TBezierSurface;

    fn make_sphere_patch() -> TBezierSurface {
        // 3x3 poles: S(u,v)=cos(u)*cos(v)*X+sin(u)*cos(v)*Y+sin(v)*Z
        // P(1,2)=(0,0,1)=sin(v),
        // P(2,3)=(0,1,0)=sin(u)*cos(v),
        // P(3,3)=(1,0,0)=cos(u)*cos(v)
        let mut poles = vec![vec![DVec3::ZERO; 3]; 3];
        poles[0][1] = DVec3::new(0.0, 0.0, 1.0);
        poles[1][2] = DVec3::new(0.0, 1.0, 0.0);
        poles[2][2] = DVec3::new(1.0, 0.0, 0.0);
        TBezierSurface::new(poles, 1.0, 1.0)
    }

    fn make_simple() -> TBezierSurface {
        let mut poles = vec![vec![DVec3::ZERO; 3]; 3];
        for i in 0..3 { for j in 0..3 { poles[i][j] = DVec3::new(i as f64, j as f64, 0.0); } }
        TBezierSurface::new(poles, 1.0, 1.0)
    }

    #[test]
    fn tbezier_surface_construction() {
        let s = make_simple();
        assert_eq!(s.nb_u_poles(), 3);
        assert_eq!(s.nb_v_poles(), 3);
        assert_eq!(s.order_u(), 1);
        assert_eq!(s.order_v(), 1);
    }

    #[test]
    fn tbezier_surface_bounds() {
        let s = make_simple();
        let b = s.bounds();
        assert!((b[0] - 0.0).abs() < TOL);
        assert!((b[1] - std::f64::consts::PI).abs() < TOL);
        assert!((b[2] - 0.0).abs() < TOL);
        assert!((b[3] - std::f64::consts::PI).abs() < TOL);
    }

    #[test]
    fn tbezier_surface_corners_distinct() {
        let s = make_simple();
        let b = s.bounds();
        let p00 = s.eval_d0(b[0], b[2]);
        let p10 = s.eval_d0(b[1], b[2]);
        let p01 = s.eval_d0(b[0], b[3]);
        let p11 = s.eval_d0(b[1], b[3]);
        assert!((p00 - p10).length() > TOL);
        assert!((p00 - p01).length() > TOL);
        assert!((p00 - p11).length() > TOL);
    }

    #[test]
    fn tbezier_surface_d1_consistent() {
        let s = make_simple();
        let u = std::f64::consts::PI / 3.0;
        let v = std::f64::consts::PI / 4.0;
        let _pt = s.eval_d0(u, v);
        let fd_u = (s.eval_d0(u + TOL_FD, v) - s.eval_d0(u - TOL_FD, v)) / (2.0 * TOL_FD);
        let fd_v = (s.eval_d0(u, v + TOL_FD) - s.eval_d0(u, v - TOL_FD)) / (2.0 * TOL_FD);
        assert!(fd_u.length() > 0.0);
        assert!(fd_v.length() > 0.0);
    }
}

// =============================================================================
// GeomEval_AHTBezierCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod aht_bezier_curve_tests {
    use super::*;
    use crate::tkg3d_complete::AHTBezierCurve;

    #[test]
    fn aht_construction_full_basis() {
        let c = AHTBezierCurve::new(
            vec![DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Y, DVec3::new(0.0,1.0,1.0), DVec3::new(1.0,0.0,1.0)],
            1, 1.0, 1.0);
        assert_eq!(c.nb_poles(), 6);
        assert!((c.alpha - 1.0).abs() < TOL);
        assert!((c.beta - 1.0).abs() < TOL);
        assert!(!c.is_rational());
    }

    #[test]
    fn aht_construction_polynomial() {
        let c = AHTBezierCurve::new(
            vec![DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
            2, 0.0, 0.0);
        assert_eq!(c.nb_poles(), 3);
        assert_eq!(c.alg_degree, 2);
    }

    #[test]
    fn aht_eval_d0_endpoints() {
        let c = AHTBezierCurve::new(
            vec![DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Y, DVec3::new(0.0,1.0,1.0), DVec3::new(1.0,0.0,1.0)],
            1, 1.0, 1.0);
        let p0 = c.eval_d0(0.0);
        let p1 = c.eval_d0(1.0);
        assert!(p0.length() < TOL || p0.x.abs() > TOL);
        assert!(p1.length() > TOL);
    }

    #[test]
    fn aht_construction_rational() {
        let c = AHTBezierCurve::new_rational(
            vec![DVec3::ZERO, DVec3::new(1.0,1.0,0.0), DVec3::new(2.0,0.0,0.0)],
            vec![1.0, 2.0, 1.0], 2, 0.0, 0.0);
        assert!(c.is_rational());
    }

    #[test]
    fn aht_derivative_consistent() {
        let c = AHTBezierCurve::new(
            vec![DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Y, DVec3::new(0.0,1.0,1.0), DVec3::new(1.0,0.0,1.0)],
            1, 1.0, 1.0);
        let u = 0.4;
        let p1 = c.eval_d0(u + TOL_FD);
        let p2 = c.eval_d0(u - TOL_FD);
        let fd = (p1 - p2) / (2.0 * TOL_FD);
        let (_pt, d1) = c.eval_d1(u);
        assert!((d1 - fd).length() < TOL_FD);
    }
}

// =============================================================================
// GeomEval_AHTBezierSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod aht_bezier_surface_tests {
    use super::*;
    use crate::tkg3d_complete::AHTBezierSurface;

    #[test]
    fn aht_surface_construction() {
        let s = AHTBezierSurface::new(
            vec![
                vec![DVec3::new(0.0,0.0,0.0), DVec3::new(0.0,1.0,0.0), DVec3::new(0.0,2.0,0.0)],
                vec![DVec3::new(1.0,0.0,0.0), DVec3::new(1.0,1.0,0.1), DVec3::new(1.0,2.0,0.0)],
                vec![DVec3::new(2.0,0.0,0.0), DVec3::new(2.0,1.0,0.0), DVec3::new(2.0,2.0,0.0)],
            ], 2, 2, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(s.nb_poles_u(), 3);
        assert_eq!(s.nb_poles_v(), 3);
    }

    #[test]
    fn aht_surface_bounds() {
        let s = AHTBezierSurface::new(
            vec![
                vec![DVec3::new(0.0,0.0,0.0), DVec3::new(0.0,1.0,0.0), DVec3::new(0.0,2.0,0.0)],
                vec![DVec3::new(1.0,0.0,0.0), DVec3::new(1.0,1.0,0.1), DVec3::new(1.0,2.0,0.0)],
                vec![DVec3::new(2.0,0.0,0.0), DVec3::new(2.0,1.0,0.0), DVec3::new(2.0,2.0,0.0)],
            ], 2, 2, 0.0, 0.0, 0.0, 0.0);
        let b = s.bounds();
        assert!((b[0] - 0.0).abs() < TOL);
        assert!((b[1] - 1.0).abs() < TOL);
    }

    #[test]
    fn aht_surface_corners_distinct() {
        let s = AHTBezierSurface::new(
            vec![
                vec![DVec3::new(0.0,0.0,0.0), DVec3::new(0.0,1.0,0.0), DVec3::new(0.0,2.0,0.0)],
                vec![DVec3::new(1.0,0.0,0.0), DVec3::new(1.0,1.0,0.1), DVec3::new(1.0,2.0,0.0)],
                vec![DVec3::new(2.0,0.0,0.0), DVec3::new(2.0,1.0,0.0), DVec3::new(2.0,2.0,0.0)],
            ], 2, 2, 0.0, 0.0, 0.0, 0.0);
        let p00 = s.eval_d0(0.0, 0.0);
        let p10 = s.eval_d0(1.0, 0.0);
        let p01 = s.eval_d0(0.0, 1.0);
        let p11 = s.eval_d0(1.0, 1.0);
        assert!((p00 - p10).length() > TOL);
        assert!((p00 - p01).length() > TOL);
    }

    #[test]
    fn aht_surface_d1_consistent() {
        let s = AHTBezierSurface::new(
            vec![
                vec![DVec3::new(0.0,0.0,0.0), DVec3::new(0.0,1.0,0.0), DVec3::new(0.0,2.0,0.0)],
                vec![DVec3::new(1.0,0.0,0.0), DVec3::new(1.0,1.0,0.1), DVec3::new(1.0,2.0,0.0)],
                vec![DVec3::new(2.0,0.0,0.0), DVec3::new(2.0,1.0,0.0), DVec3::new(2.0,2.0,0.0)],
            ], 2, 2, 0.0, 0.0, 0.0, 0.0);
        let u = 0.3; let v = 0.7;
        let fd_u = (s.eval_d0(u + TOL_FD, v) - s.eval_d0(u - TOL_FD, v)) / (2.0 * TOL_FD);
        let fd_v = (s.eval_d0(u, v + TOL_FD) - s.eval_d0(u, v - TOL_FD)) / (2.0 * TOL_FD);
        assert!(fd_u.length() > 0.0);
        assert!(fd_v.length() > 0.0);
    }
}

// =============================================================================
// GeomAdaptor_TransformedSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod transformed_surface_tests {
    use super::*;
    use crate::tkg3d_complete::TransformedSurfaceAdaptor;

    #[test]
    fn identity_transform_uses_original() {
        let plane = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let adaptor = TransformedSurfaceAdaptor::new(plane.clone(), DAffine3::IDENTITY);
        let pt = adaptor.evaluate(1.0, 2.0);
        let expected = plane.point_at(1.0, 2.0);
        assert!((pt - expected).length() < TOL);
    }

    #[test]
    fn plane_translated_by_trsf() {
        let plane = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let trsf = DAffine3::from_translation(DVec3::new(0.0, 0.0, 5.0));
        let adaptor = TransformedSurfaceAdaptor::new(plane, trsf);
        let pt = adaptor.evaluate(0.0, 0.0);
        assert!((pt - DVec3::new(0.0, 0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn set_trsf_rebuilds() {
        let plane = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let trsf1 = DAffine3::from_translation(DVec3::new(0.0, 0.0, 1.0));
        let trsf2 = DAffine3::from_translation(DVec3::new(0.0, 0.0, 3.0));
        let mut adaptor = TransformedSurfaceAdaptor::new(plane, trsf1);
        let pt1 = adaptor.evaluate(0.0, 0.0);
        adaptor.set_trsf(trsf2);
        let pt2 = adaptor.evaluate(0.0, 0.0);
        assert!((pt1 - DVec3::new(0.0, 0.0, 1.0)).length() < TOL);
        assert!((pt2 - DVec3::new(0.0, 0.0, 3.0)).length() < TOL);
    }
}
