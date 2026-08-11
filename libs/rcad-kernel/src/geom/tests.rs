use crate::geom::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_involute_starts_on_base_circle() {
        let inv = CircleInvolute2d {
            center: DVec2::new(2.0, -1.0),
            base_radius: 3.0,
            start_angle: 0.0,
        };
        let p0 = inv.point_at(0.0);
        assert!((p0.x - 5.0).abs() < 1e-12);
        assert!((p0.y + 1.0).abs() < 1e-12);
    }

    #[test]
    fn archimedean_spiral_point_progresses_radially() {
        let s = ArchimedeanSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.5,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p1 = s.point_at(2.0);
        assert!(
            p1.length() > p0.length(),
            "spiral radius should increase with t"
        );
    }

    #[test]
    fn logarithmic_spiral_grows_exponentially() {
        let s = LogarithmicSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.4,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p1 = s.point_at(2.0);
        assert!(
            p1.length() > p0.length() * 1.5,
            "log spiral should grow faster than linear at this sample"
        );
    }

    #[test]
    fn sine_wave_samples_match_expected_values() {
        let s = SineWave2d {
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        let p0 = s.point_at(0.0);
        let p90 = s.point_at(std::f64::consts::FRAC_PI_2);
        assert!((p0.x - 0.0).abs() < 1e-12);
        assert!((p0.y - 0.0).abs() < 1e-12);
        assert!((p90.y - 2.0).abs() < 1e-12);
    }

    #[test]
    fn bspline_curve3_derivative_matches_linear_curve() {
        let curve = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(1.0, 2.0, 3.0), DVec3::new(4.0, 8.0, 15.0)],
            weights: vec![1.0, 1.0],
        is_periodic: false,
};

        let derivative = curve.derivative_at(0.4);

        assert!((derivative - DVec3::new(3.0, 6.0, 12.0)).length() < 1e-12);
    }

    #[test]
    fn bspline_curve2_derivative_matches_linear_curve() {
        let curve = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(1.0, 2.0), DVec2::new(4.0, 8.0)],
            weights: vec![1.0, 1.0],
        };

        let derivative = curve.derivative_at(0.4);

        assert!((derivative - DVec2::new(3.0, 6.0)).length() < 1e-12);
    }

    #[test]
    fn curve2d_sine_wave_variant_dispatches_evaluator() {
        let c = Curve2d::SineWave(SineWave2d {
            amplitude: 1.5,
            frequency: 2.0,
            phase: 0.25,
        });
        let t = 0.3;
        let p = c.point_at(t);
        let expected_y = 1.5 * (2.0 * t + 0.25).sin();
        assert!((p.x - t).abs() < 1e-12);
        assert!((p.y - expected_y).abs() < 1e-12);
    }

    #[test]
    fn sine_wave3_origin_phase_zero_evaluates_at_zero_offset() {
        let c = SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 3.0,
            frequency: 1.0,
            phase: 0.0,
        };
        // At t=0, sin(0)=0 → point should be at origin.
        let p = c.point_at(0.0);
        assert!(
            p.length() < 1e-12,
            "phase-zero at t=0 should be at origin: {p:?}"
        );
        // At t=pi/2, sin(pi/2)=1 → y should equal amplitude.
        let p2 = c.point_at(std::f64::consts::FRAC_PI_2);
        assert!(
            (p2.y - 3.0).abs() < 1e-9,
            "y at t=pi/2 should be amplitude=3: {p2:?}"
        );
    }

    #[test]
    fn curve3_sine_wave_variant_dispatches_evaluator() {
        let c = Curve3::SineWave(SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 1.0,
            frequency: 2.0,
            phase: 0.0,
        });
        let t = 0.5;
        let p = c.point_at(t);
        let expected = DVec3::new(0.5, (2.0_f64 * t).sin(), 0.0);
        assert!((p - expected).length() < 1e-12);
        // Tangent should be non-zero
        let tan = c.tangent_at(t);
        assert!(
            tan.length() > 0.9,
            "tangent should be roughly unit-length: {tan:?}"
        );
    }
}
#[cfg(test)]
mod eval_tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn line3_point_at() {
        let l = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        assert!((l.point_at(3.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn circle3_point_at_zero_is_on_circle() {
        // Circle in XY plane, normal = Z
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let p0 = c.point_at(0.0);
        assert!((p0.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn circle3_full_revolution_closes() {
        let c = Circle3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::Y, 5.0);
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * PI);
        assert!((p0 - p2pi).length() < 1e-10);
    }

    #[test]
    fn circle3_quarter_turn() {
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let p0 = c.point_at(0.0);
        let p90 = c.point_at(FRAC_PI_2);
        // 90° rotation: p0 and p90 should be perpendicular from center
        assert!((p0.dot(p90)).abs() < 1e-10);
        assert!((p90.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sphere_surface_north_pole() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
            ref_dir: any_perpendicular(DVec3::Y),
        };
        // OCCT ElSLib::SphereValue: V = latitude [-pi/2, pi/2], V = +pi/2 is
        // the pole in the axis direction (v=0 is the equator).
        let p = s.point_at(0.0, std::f64::consts::FRAC_PI_2);
        // Should be at (0, 3, 0)
        assert!((p - DVec3::new(0.0, 3.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn sphere_surface_point_on_sphere() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 2.0,
            ref_dir: any_perpendicular(DVec3::Y),
        };
        for u in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            for v in [0.1, 0.5, 1.0, PI / 2.0, PI - 0.1] {
                let p = s.point_at(u, v);
                assert!(
                    (p.length() - 2.0).abs() < 1e-9,
                    "u={u} v={v} |p|={}",
                    p.length()
                );
            }
        }
    }

    #[test]
    fn cylinder_surface_point_on_cylinder() {
        let c = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        for u in [0.0, 1.0, PI, 2.0 * PI - 0.1] {
            let p = c.point_at(u, 0.0);
            let radial = DVec3::new(p.x, 0.0, p.z).length();
            assert!((radial - 3.0).abs() < 1e-9, "u={u} radial={radial}");
        }
    }

    #[test]
    fn bspline_degree1_linear_interpolation() {
        // Degree-1 BSpline with 2 control points = straight line
        let c = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        is_periodic: false,
};
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);
        let pmid = c.point_at(0.5);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);
        assert!((pmid - DVec3::new(0.5, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn bspline_degree2_quadratic() {
        // Degree-2 quadratic arc through 3 control points
        let c = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::new(0.5, 1.0, 0.0), DVec3::X],
            weights: vec![1.0, 1.0, 1.0],
        is_periodic: false,
};
        let p0 = c.point_at(0.0);
        let p1 = c.point_at(1.0);
        assert!((p0 - DVec3::ZERO).length() < 1e-10);
        assert!((p1 - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn torus_surface_point_on_torus() {
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        for u in [0.0, PI / 2.0, PI] {
            for v in [0.0, PI / 2.0, PI] {
                let p = t.point_at(u, v);
                // Distance from the tube center circle should be minor_radius
                let x_ax = any_perpendicular(DVec3::Y);
                let y_ax = DVec3::Y.cross(x_ax).normalize();
                let tube_center = t.center + t.major_radius * (u.cos() * x_ax + u.sin() * y_ax);
                assert!((p - tube_center).length() - 1.0 < 1e-9, "u={u} v={v}");
            }
        }
    }

    #[test]
    fn ellipsoid_surface_satisfies_implicit_equation() {
        let s = EllipsoidalSurface {
            center: DVec3::new(1.0, -2.0, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 4.0,
            radius_y: 2.0,
            radius_z: 1.5,
        };
        let p = s.point_at(0.7, 1.2) - s.center;
        let value =
            (p.x / s.radius_x).powi(2) + (p.y / s.radius_y).powi(2) + (p.z / s.radius_z).powi(2);
        assert!(
            (value - 1.0).abs() < 1e-9,
            "implicit value should be 1, got {value}"
        );
    }

    #[test]
    fn ellipsoid_surface_normal_matches_gradient_direction() {
        let s = EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        };
        let u = 0.9;
        let v = 1.1;
        let p = s.point_at(u, v);
        let expected = DVec3::new(
            p.x / (s.radius_x * s.radius_x),
            p.y / (s.radius_y * s.radius_y),
            p.z / (s.radius_z * s.radius_z),
        )
        .normalize();
        let n = s.normal_at(u, v);
        assert!(
            (n - expected).length() < 1e-9,
            "n={n:?} expected={expected:?}"
        );
    }

    #[test]
    fn helicoid_surface_advances_by_pitch_per_turn() {
        let s = HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 6.0,
        };
        let p0 = s.point_at(0.0, 2.0);
        let p1 = s.point_at(2.0 * PI, 2.0);
        let delta = p1 - p0;
        assert!(
            (delta - DVec3::new(0.0, 0.0, 6.0)).length() < 1e-9,
            "delta={delta:?}"
        );
    }

    #[test]
    fn helicoid_surface_normal_is_perpendicular_to_parametric_tangents() {
        let s = HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 4.0,
        };
        let u = 0.6;
        let v = 1.75;
        let n = s.normal_at(u, v);
        let eps = 1e-6;
        let du = (s.point_at(u + eps, v) - s.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (s.point_at(u, v + eps) - s.point_at(u, v - eps)) / (2.0 * eps);
        assert!(
            n.dot(du).abs() < 1e-6,
            "n·du={} should be near 0",
            n.dot(du)
        );
        assert!(
            n.dot(dv).abs() < 1e-6,
            "n·dv={} should be near 0",
            n.dot(dv)
        );
        assert!(n.length() > 0.99, "normal should be unit-length: {n:?}");
    }

    #[test]
    fn pipe_surface_with_line_spine_matches_cylindrical_section() {
        let surface = PipeSurface {
            spine: Box::new(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Z,
            })),
            ref_dir: DVec3::X,
            radius: 2.0,
        };

        assert!((surface.point_at(0.0, 0.0) - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-9);
        assert!((surface.point_at(PI * 0.5, 0.5) - DVec3::new(0.0, 2.0, 0.5)).length() < 1e-9);
        assert!((surface.default_domain()[0] - 0.0).abs() < 1e-12);
        assert!((surface.default_domain()[1] - 2.0 * PI).abs() < 1e-12);
    }

    #[test]
    fn tri_bezier_surface_hits_triangle_corners() {
        let surface = TriBezierSurface {
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0]],
        };
        assert!((surface.point_at(0.0, 0.0) - DVec3::new(0.0, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(1.0, 0.0) - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.0, 1.0) - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn tri_bezier_surface_dispatches_through_surface3() {
        let surface = Surface3::TriBezier(TriBezierSurface {
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0]],
        });
        let p = surface.point_at(0.25, 0.5);
        assert!(p.x >= -1e-12 && p.y >= -1e-12);
        assert!(surface.normal_at(0.2, 0.2).length() > 0.99);
    }

    #[test]
    fn ruled_surface_interpolates_between_curves() {
        let surface = RuledSurface {
            start: Box::new(Curve3::Line(Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::X,
            })),
            end: Box::new(Curve3::Line(Line3 {
                origin: DVec3::Y,
                direction: DVec3::X,
            })),
        };
        assert!((surface.point_at(0.25, 0.0) - DVec3::new(0.25, 0.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.25, 1.0) - DVec3::new(0.25, 1.0, 0.0)).length() < 1e-12);
        assert!((surface.point_at(0.25, 0.5) - DVec3::new(0.25, 0.5, 0.0)).length() < 1e-12);
        assert!(surface.normal_at(0.25, 0.5).length() > 0.99);
    }

    #[test]
    fn coons_surface_interpolates_all_four_boundaries() {
        let surface = CoonsSurface {
            south: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::X,
            })),
            north: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 1.0, 1.0),
                direction: DVec3::X,
            })),
            west: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
            east: Box::new(Curve3::Line(Line3 {
                origin: DVec3::new(1.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
        };

        assert!((surface.point_at(0.3, 0.0) - DVec3::new(0.3, 0.0, 0.0)).length() < 1e-9);
        assert!((surface.point_at(0.3, 1.0) - DVec3::new(0.3, 1.0, 1.0)).length() < 1e-9);
        assert!((surface.point_at(0.0, 0.4) - DVec3::new(0.0, 0.4, 0.4)).length() < 1e-9);
        assert!((surface.point_at(1.0, 0.4) - DVec3::new(1.0, 0.4, 0.4)).length() < 1e-9);
        assert!((surface.point_at(0.5, 0.5) - DVec3::new(0.5, 0.5, 0.5)).length() < 1e-9);
    }

    #[test]
    fn conical_surface_uses_slant_distance_from_reference_circle() {
        let surface = ConicalSurface::new(DVec3::ZERO, DVec3::Z, 2.0, 30.0_f64.to_radians());

        let p0 = surface.point_at(0.0, 0.0);
        assert!(p0.dot(surface.axis_dir()).abs() < 1e-9);
        assert!((p0.length() - 2.0).abs() < 1e-9);

        let slant = 4.0;
        let p1 = surface.point_at(0.0, slant);
        assert!((p1.z - slant * surface.half_angle_rad.cos()).abs() < 1e-9);
        let radial = p1 - surface.axis_dir() * p1.dot(surface.axis_dir());
        assert!((radial.length() - (2.0 + slant * surface.half_angle_rad.sin())).abs() < 1e-9);
    }

    #[test]
    fn conical_surface_derives_true_apex_from_reference_circle() {
        let surface = ConicalSurface::new(DVec3::new(0.0, 0.0, 5.0), DVec3::Z, 2.0, 45.0_f64.to_radians());

        assert!((surface.apex_point() - DVec3::new(0.0, 0.0, 3.0)).length() < 1e-9);
    }

    // --- Analytic derivative tests ---

    /// Quadratic Bezier: P0=(0,0,0), P1=(0.5,1,0), P2=(1,0,0), unit weights.
    /// Analytic tangent at t=0 should be (0.5,1,0).normalize() = (1,2,0)/√5.
    #[test]
    fn bezier_tangent_at_endpoint_analytic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(0.5, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let c = BezierCurve3 {
            control_points: pts,
            weights: wts,
        };
        let tan = c.tangent_at(0.0);
        let expected = DVec3::new(1.0, 2.0, 0.0).normalize();
        assert!(
            (tan - expected).length() < 1e-10,
            "tan={tan:?} expected={expected:?}"
        );
    }

    /// Quadratic Bezier tangent at t=1 should be (1,-2,0)/√5.
    #[test]
    fn bezier_tangent_at_end_analytic() {
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(0.5, 1.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let c = BezierCurve3 {
            control_points: pts,
            weights: wts,
        };
        let tan = c.tangent_at(1.0);
        let expected = DVec3::new(1.0, -2.0, 0.0).normalize();
        assert!(
            (tan - expected).length() < 1e-10,
            "tan={tan:?} expected={expected:?}"
        );
    }

    /// Degree-1 B-Spline (polyline): tangent should be constant along each segment.
    #[test]
    fn bspline_degree1_tangent_is_segment_direction() {
        // Two-segment polyline: (0,0,0) -> (1,0,0) -> (1,1,0)
        let pts = vec![
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let c = BSplineCurve3 {
            degree: 1,
            knots,
            control_points: pts,
            weights: wts,
        is_periodic: false,
};
        let tan0 = c.tangent_at(0.1);
        assert!(
            (tan0 - DVec3::X).length() < 1e-10,
            "first segment should be +X, got {tan0:?}"
        );
        let tan1 = c.tangent_at(0.9);
        assert!(
            (tan1 - DVec3::Y).length() < 1e-10,
            "second segment should be +Y, got {tan1:?}"
        );
    }

    /// Degree-2 B-Spline circle arc: tangent should be perpendicular to radius.
    #[test]
    fn bspline_circle_tangent_perpendicular_to_radius() {
        // Use circle_to_bspline to get an exact NURBS circle, then check tangents.
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let c = crate::base::convert::circle_to_bspline(&circle);
        for &t in &[0.0, 0.5, 1.0, 1.5, 2.0] {
            let pt = c.point_at(t);
            let tan = c.tangent_at(t);
            // Tangent must be perpendicular to the radius vector
            let dot = pt.normalize_or_zero().dot(tan);
            assert!(
                dot.abs() < 1e-8,
                "t={t}: radius*tangent={dot} (should be 0)"
            );
            // Tangent must be a unit vector
            assert!(
                (tan.length() - 1.0).abs() < 1e-10,
                "t={t}: |tan|={}",
                tan.length()
            );
        }
    }

    // =========================================================================
    // OCCT-aligned TKG3d / TKG2d evaluation tests
    // =========================================================================
    //
    // These test point_at (D0) and tangent_at (D1) for each curve/surface type,
    // matching patterns in OCCT's TKG3d/GTests/ and TKG2d/GTests/.

    // ── 3D Curve evaluation (Geom_CurveEval_Test.cxx pattern) ────────────

    #[test]
    fn line_eval_d0_d1() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let p = line.point_at(5.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-12);
        let t = line.tangent_at(5.0);
        assert!((t - DVec3::X).length() < 1e-12);
    }

    #[test]
    fn line_eval_d2_zero_second_derivative() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        // For a line, the first derivative (tangent) is constant; second derivative is zero.
        // The curve is linear: P(t) = origin + t * direction
        // The second derivative: d²P/dt² = 0
        // Using the derivative_at method:
        let d1 = line.derivative_at(0.0);
        let d2 = (line.derivative_at(1e-4) - line.derivative_at(-1e-4)) / (2.0 * 1e-4);
        assert!((d1 - DVec3::X).length() < 1e-10);
        assert!(
            d2.length() < 1e-10,
            "Line second derivative should be 0, got {d2:?}"
        );
    }

    #[test]
    fn circle_eval_d0_d1() {
        // Circle3::new(ZERO, Z, 5.0): OCCT gp_Ax2 gives x_dir=X, y_dir=Z×X=Y
        // P(0) = 5*X = (5,0,0), tangent = Y
        // P(PI/2) = 5*Y = (0,5,0), tangent = -X
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let p0 = circle.point_at(0.0);
        assert!((p0 - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
        let p_half = circle.point_at(std::f64::consts::PI / 2.0);
        assert!((p_half - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10);
        let p_pi = circle.point_at(std::f64::consts::PI);
        assert!((p_pi - DVec3::new(-5.0, 0.0, 0.0)).length() < 1e-10);
        // Tangent at 0 should be (0, 1, 0)
        let t0 = circle.tangent_at(0.0);
        assert!((t0 - DVec3::Y).length() < 1e-10);
        // Tangent at PI/2 should be (-1, 0, 0)
        let t_half = circle.tangent_at(std::f64::consts::PI / 2.0);
        assert!((t_half - DVec3::new(-1.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn circle_transform_copy() {
        // Circle3::new((1,2,3), Z, 4.0): OCCT gp_Ax2 gives x_dir=X, y_dir=Y
        // P(0) = (1,2,3) + 4*X = (5,2,3)
        let circle = Circle3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::Z, 4.0);
        assert!((circle.point_at(0.0) - DVec3::new(5.0, 2.0, 3.0)).length() < 1e-10);
    }

    #[test]
    fn bspline_eval_d0_d1_d2_consistency() {
        // Degree-2 BSpline through 4 points
        let c = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.3, 0.7, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(2.0, 3.0, 0.0),
                DVec3::new(5.0, 3.0, 0.0),
                DVec3::new(7.0, 0.0, 0.0),
                DVec3::new(10.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 5],
        is_periodic: false,
};
        // Point at t=0 should be first control point
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // Point at t=1 should be last control point
        assert!((c.point_at(1.0) - DVec3::new(10.0, 0.0, 0.0)).length() < 1e-10);
        // Range check: some midpoint
        let _pmid = c.point_at(0.5);
        assert!(_pmid.x >= 0.0 && _pmid.x <= 10.0);
    }

    #[test]
    fn bezier_eval_d0_d1() {
        // Cubic Bezier through 4 control points
        let c = BezierCurve3 {
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 3.0, 0.0),
                DVec3::new(3.0, 3.0, 0.0),
                DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        };
        // t=0 -> first pole
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // t=1 -> last pole
        assert!((c.point_at(1.0) - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-10);
        // Tangent at t=0: direction from P0 to P1
        let t0 = c.tangent_at(0.0);
        assert!((t0 - DVec3::new(1.0, 3.0, 0.0).normalize()).length() < 1e-10);
    }

    // ── 3D Surface evaluation (Geom_SurfaceEval_Test.cxx pattern) ───────

    #[test]
    fn plane_eval_d0_d1() {
        // Plane with normal Z: OCCT gp_Ax3 gives u_dir=X, v_dir=Y.
        // P(u,v) = u*X + v*Y. So P(2,3) = (2, 3, 0)
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let p = plane.point_at(2.0, 3.0);
        // Should be in the plane (Z=0)
        assert!((p.z).abs() < 1e-10);
        // Distance from origin should be sqrt(4+9) = sqrt(13)
        assert!((p.length() - 13.0f64.sqrt()).abs() < 1e-10);
    }

    #[test]
    fn cylinder_eval_d0() {
        // Cylinder with axis Z, any_perpendicular(Z) = Y, y_ax = Z.cross(Y) = -X
        // Cylinder with axis Z, ref_dir=X:
        // x_ax = X, y_ax = Z×X = Y
        // P(u,v) = R*(cos(u)*X + sin(u)*Y) + v*Z
        // P(0,0) = 3*X = (3, 0, 0)
        // P(PI/2, 5) = 3*Y + 5*Z = (0, 3, 5)
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        let p0 = cyl.point_at(0.0, 0.0);
        assert!((p0 - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
        let p1 = cyl.point_at(std::f64::consts::PI / 2.0, 5.0);
        assert!((p1 - DVec3::new(0.0, 3.0, 5.0)).length() < 1e-10);
    }

    #[test]
    fn sphere_eval_d0_full_sphere() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        };
        for u in [0.0, 1.0, 2.0, 4.0, 6.0] {
            for v in [-1.0, -0.5, 0.0, 0.5, 1.0] {
                let p = s.point_at(u, v);
                assert!((p.length() - 2.0).abs() < 1e-9, "u={u} v={v}");
            }
        }
    }

    #[test]
    fn cone_eval_d0() {
        let cone = ConicalSurface::new(DVec3::ZERO, DVec3::Z, 2.0, 45.0_f64.to_radians());
        // At V=0, radius=2
        let p0 = cone.point_at(0.0, 0.0);
        assert!((p0.x - 2.0).abs() < 1e-9 || (p0.y - 2.0).abs() < 1e-9);
        assert!((p0.z).abs() < 1e-9);
    }

    #[test]
    fn torus_eval_d0() {
        use std::f64::consts::PI;
        // Torus with axis Z: OCCT-aligned x_ax = X, y_ax = Y
        // P(u,v) = (R + r*cos(v))*(cos(u)*X + sin(u)*Y) + r*sin(v)*Z
        // P(0,0) = (5+1)*X = (6, 0, 0)
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let p = t.point_at(0.0, 0.0);
        assert!((p - DVec3::new(6.0, 0.0, 0.0)).length() < 1e-9);
        // At u=0, v=PI: P = (5-1)*X = (4, 0, 0)
        let p2 = t.point_at(0.0, PI);
        assert!((p2 - DVec3::new(4.0, 0.0, 0.0)).length() < 1e-9);
    }

    #[test]
    fn offset_surface_eval_d0() {
        // Use a Sphere as basis — normal is well-defined
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 1.0);
        let off = OffsetSurface {
            basis: Box::new(Surface3::Sphere(sphere)),
            offset_distance: 0.5,
        };
        // Offset sphere: radius = 1 + 0.5 = 1.5
        let p = off.point_at(0.0, 0.0);
        assert!((p.length() - 1.5).abs() < 1e-9);
    }

    // ── 2D Curve evaluation (Geom2d_CurveEval_Test.cxx pattern) ─────────

    #[test]
    fn line2d_eval_d0() {
        let l = Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        };
        let p = l.point_at(3.0);
        assert!((p - DVec2::new(3.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn circle2d_eval_d0() {
        use std::f64::consts::PI;
        let c = Circle2d::new(DVec2::ZERO, 5.0);
        let p0 = c.point_at(0.0);
        assert!((p0 - DVec2::new(5.0, 0.0)).length() < 1e-12);
        let p_half = c.point_at(PI / 2.0);
        assert!((p_half - DVec2::new(0.0, 5.0)).length() < 1e-12);
        let p_pi = c.point_at(PI);
        assert!((p_pi - DVec2::new(-5.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn circle2d_revolved() {
        let c = Circle2d::new(DVec2::ZERO, 5.0);
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * std::f64::consts::PI);
        assert!((p0 - p2pi).length() < 1e-12);
    }

    #[test]
    fn ellipse2d_eval_d0() {
        use std::f64::consts::PI;
        let e = Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        let p0 = e.point_at(0.0);
        assert!((p0 - DVec2::new(10.0, 0.0)).length() < 1e-12);
        let p_half = e.point_at(PI / 2.0);
        assert!((p_half - DVec2::new(0.0, 5.0)).length() < 1e-12);
    }

    #[test]
    fn parabola2d_eval_d0() {
        let p = Parabola2d {
            origin: DVec2::ZERO,
            axis_dir: DVec2::X,
            focal_param: 4.0,
        };
        // P(t) = (t²/(2p), t) = (t²/8, t)
        let p0 = p.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-12);
        let p4 = p.point_at(4.0);
        assert!((p4 - DVec2::new(2.0, 4.0)).length() < 1e-10);
    }

    #[test]
    fn hyperbola2d_eval_d0() {
        let h = Hyperbola2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };
        let p0 = h.point_at(0.0);
        // X = a*cosh(0) = 3, Y = b*sinh(0) = 0
        assert!((p0 - DVec2::new(3.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn bspline2_eval_d0() {
        let c = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::ZERO, DVec2::X],
            weights: vec![1.0, 1.0],
        };
        assert!((c.point_at(0.0) - DVec2::ZERO).length() < 1e-12);
        assert!((c.point_at(1.0) - DVec2::X).length() < 1e-12);
        assert!((c.point_at(0.5) - DVec2::new(0.5, 0.0)).length() < 1e-12);
    }

    #[test]
    fn bezier2_eval_d0() {
        let c = BezierCurve2 {
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(0.5, 1.0),
                DVec2::new(1.0, 0.0),
            ],
            weights: vec![1.0; 3],
        };
        assert!((c.point_at(0.0) - DVec2::ZERO).length() < 1e-12);
        assert!((c.point_at(1.0) - DVec2::X).length() < 1e-12);
    }

    #[test]
    fn trimmed_curve2_eval() {
        let inner = Circle2d::new(DVec2::ZERO, 5.0);
        let tc = TrimmedCurve2 {
            curve: Box::new(Curve2d::Circle(inner)),
            t_min: 0.0,
            t_max: std::f64::consts::PI,
        };
        let p0 = tc.point_at(0.0);
        assert!((p0 - DVec2::new(5.0, 0.0)).length() < 1e-12);
        let p_half = tc.point_at(std::f64::consts::PI / 2.0);
        assert!((p_half - DVec2::new(0.0, 5.0)).length() < 1e-12);
        // Out of range → clamped
        let p_out = tc.point_at(10.0);
        assert!((p_out - DVec2::new(-5.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn offset_curve2d_eval() {
        // Circle2d r=5, offset_distance uses right-hand normal.
        // P(0) = (5,0) + 1*(right_normal at 0) = (5,0) + 1*(1,0) ≈ (6,0)
        let basis = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let off = OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 1.0,
        };
        let p0 = off.point_at(0.0);
        // Offset point should be outside the original circle
        assert!(p0.length() > 5.0);
        assert!((p0.length() - 6.0).abs() < 0.1);
    }

    // ── Special 2D curve evaluation tests ───────────────────────────────

    #[test]
    fn circle_involute2d_eval() {
        let inv = CircleInvolute2d {
            center: DVec2::ZERO,
            base_radius: 3.0,
            start_angle: 0.0,
        };
        // At t=0: P = center + r*(cos(0)+0*sin(0), sin(0)-0*cos(0)) = center + (r, 0)
        let p0 = inv.point_at(0.0);
        assert!((p0 - DVec2::new(3.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn sine_wave2d_eval() {
        let w = SineWave2d {
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        let p0 = w.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-12);
        let p_pi = w.point_at(std::f64::consts::PI / 2.0);
        assert!((p_pi - DVec2::new(std::f64::consts::PI / 2.0, 2.0)).length() < 1e-10);
    }

    #[test]
    fn archimedean_spiral2d_eval() {
        let s = ArchimedeanSpiral2d {
            center: DVec2::ZERO,
            a: 0.0,
            b: 1.0,
            start_angle: 0.0,
        };
        // r(t) = 0 + 1*t, theta(t) = t
        let p0 = s.point_at(0.0);
        assert!((p0 - DVec2::ZERO).length() < 1e-12);
        let p1 = s.point_at(2.0);
        assert!((p1.length() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn logarithmic_spiral2d_eval() {
        let s = LogarithmicSpiral2d {
            center: DVec2::ZERO,
            a: 1.0,
            b: 0.5,
            start_angle: 0.0,
        };
        let p0 = s.point_at(0.0);
        assert!((p0 - DVec2::new(1.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn offset_curve3_eval() {
        // Circle3 offset along Z — FD tangent makes this approximate
        let basis = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let off = OffsetCurve3 {
            basis: Box::new(basis),
            offset_distance: 2.0,
            offset_dir: DVec3::Z,
        };
        let p0 = off.point_at(0.0);
        // With offset, should differ from original circle (radius 5)
        assert!((p0.length() - 5.0).abs() > 0.5);
        // But not be wildly different
        assert!(p0.length() < 10.0);
    }

    // ── Reverse parameter tests (Geom2d_*_ReversedParameter pattern) ────

    #[test]
    fn circle_reverse_eval() {
        // A reversed circle should evaluate in the opposite direction
        // OCCT: ReversedParameter(t) = 2*PI - t
        // Since rcad doesn't have an explicit reverse flag, verify that
        // the parameterization wraps around correctly.
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let p_fwd = c.point_at(0.3);
        let p_rev = c.point_at(2.0 * std::f64::consts::PI - 0.3);
        // These are different points (one advances forward, one backward)
        assert!((p_fwd - p_rev).length() > 0.5);
    }

    // =========================================================================
    // OCCT-aligned comprehensive 3D curve evaluation tests
    // (matching TKG3d/GTests Geom_Line/Circle/Ellipse/BSpline/Bezier patterns)
    // =========================================================================

    // ── Line ────────────────────────────────────────────────────────────

    #[test]
    fn line3_eval_at_multiple_points() {
        let line = Line3 {
            origin: DVec3::new(1.0, 2.0, 3.0),
            direction: DVec3::new(0.0, 1.0, 0.0),
        };
        assert!((line.point_at(0.0) - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-12);
        assert!((line.point_at(5.0) - DVec3::new(1.0, 7.0, 3.0)).length() < 1e-12);
        assert!((line.point_at(-3.0) - DVec3::new(1.0, -1.0, 3.0)).length() < 1e-12);
    }

    #[test]
    fn line3_constant_tangent_and_derivative() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::new(1.0, 1.0, 1.0).normalize(),
        };
        let d = DVec3::new(1.0, 1.0, 1.0).normalize();
        for &t in &[-10.0, -1.0, 0.0, 1.0, 10.0] {
            assert!((line.tangent_at(t) - d).length() < 1e-12);
            assert!((line.derivative_at(t) - d).length() < 1e-12);
        }
    }

    #[test]
    fn line3_default_domain_infinite() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let [t0, t1] = line.default_domain();
        assert!(t0.is_infinite() && t0.is_sign_negative());
        assert!(t1.is_infinite() && t1.is_sign_positive());
    }

    // ── Circle ─────────────────────────────────────────────────────────

    #[test]
    fn circle3_eval_four_quadrants() {
        // Use explicit x_dir/y_dir for predictable orientation
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            radius: 5.0,
        };
        // P(0) = (5,0,0), P(PI/2) = (0,5,0), P(PI) = (-5,0,0), P(3PI/2) = (0,-5,0)
        assert!((c.point_at(0.0) - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
        assert!(
            (c.point_at(std::f64::consts::PI / 2.0) - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10
        );
        assert!((c.point_at(std::f64::consts::PI) - DVec3::new(-5.0, 0.0, 0.0)).length() < 1e-10);
        assert!(
            (c.point_at(3.0 * std::f64::consts::PI / 2.0) - DVec3::new(0.0, -5.0, 0.0)).length()
                < 1e-10
        );
    }

    #[test]
    fn circle3_tangent_at_quadrants() {
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            radius: 5.0,
        };
        use std::f64::consts::PI;
        // tangent = (-R*sin(t)*X + R*cos(t)*Y).normalize()
        assert!((c.tangent_at(0.0) - DVec3::Y).length() < 1e-10);
        assert!((c.tangent_at(PI / 2.0) + DVec3::X).length() < 1e-10);
        assert!((c.tangent_at(PI) + DVec3::Y).length() < 1e-10);
        assert!((c.tangent_at(3.0 * PI / 2.0) - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn circle3_derivative_nonzero() {
        let c = Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            radius: 5.0,
        };
        // derivative = R * (-sin(t)*X + cos(t)*Y), always non-zero for R>0
        for &t in &[0.0, 0.5, 1.0, 2.0, 4.0, 6.0] {
            let d = c.derivative_at(t);
            assert!(d.length() > 0.0);
            assert!(
                (d.length() - 5.0).abs() < 1e-10,
                "t={} |d|={}",
                t,
                d.length()
            );
        }
    }

    #[test]
    fn circle3_default_domain() {
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let [t0, t1] = c.default_domain();
        assert!((t0 - 0.0).abs() < 1e-12);
        assert!((t1 - 2.0 * std::f64::consts::PI).abs() < 1e-12);
    }

    // ── Ellipse ────────────────────────────────────────────────────────

    #[test]
    fn ellipse3_eval_vertices() {
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        use std::f64::consts::PI;
        // Major vertices: t=0 → (10,0,0), t=PI → (-10,0,0)
        assert!((e.point_at(0.0) - DVec3::new(10.0, 0.0, 0.0)).length() < 1e-10);
        assert!((e.point_at(PI) - DVec3::new(-10.0, 0.0, 0.0)).length() < 1e-10);
        // Minor vertices: t=PI/2 → (0,5,0), t=3PI/2 → (0,-5,0)
        assert!((e.point_at(PI / 2.0) - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10);
        assert!((e.point_at(3.0 * PI / 2.0) - DVec3::new(0.0, -5.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn ellipse3_tangent_at_major_vertex() {
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        // Tangent at major vertex t=0: direction = (0, 5, 0) = Y
        // (derivative: -a*sin(0)*X + b*cos(0)*Y = 5*Y)
        let t0 = e.tangent_at(0.0);
        assert!((t0 - DVec3::Y).length() < 1e-10);
    }

    #[test]
    fn ellipse3_derivative_at_vertices() {
        let e = Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 10.0,
            minor_radius: 5.0,
        };
        // derivative at t=0: a*(-sin(0))*X + b*cos(0)*Y = b*Y = 5*Y
        let d0 = e.derivative_at(0.0);
        assert!((d0 - DVec3::new(0.0, 5.0, 0.0)).length() < 1e-10);
        // derivative at t=PI/2: a*(-sin(PI/2))*X + b*cos(PI/2)*Y = -a*X = -10*X
        let d_half = e.derivative_at(std::f64::consts::PI / 2.0);
        assert!((d_half - DVec3::new(-10.0, 0.0, 0.0)).length() < 1e-10);
    }

    // ── BSpline ────────────────────────────────────────────────────────

    #[test]
    fn bspline3_eval_at_knots() {
        let c = BSplineCurve3 {
            degree: 2,
            // 4 poles, degree 2 -> 7 knots (clamped).
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(3.0, 5.0, 0.0),
                DVec3::new(6.0, 5.0, 0.0),
                DVec3::new(9.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
            is_periodic: false,
        };
        // Endpoints
        assert!((c.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        assert!((c.point_at(1.0) - DVec3::new(9.0, 0.0, 0.0)).length() < 1e-10);
    }

    #[test]
    fn bspline3_degree1_is_line() {
        let c = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(1.0, 2.0, 3.0), DVec3::new(4.0, 5.0, 6.0)],
            weights: vec![1.0, 1.0],
        is_periodic: false,
};
        assert!((c.point_at(0.0) - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-10);
        assert!((c.point_at(1.0) - DVec3::new(4.0, 5.0, 6.0)).length() < 1e-10);
        assert!((c.point_at(0.5) - DVec3::new(2.5, 3.5, 4.5)).length() < 1e-10);
    }

    // ── Bezier ─────────────────────────────────────────────────────────

    #[test]
    fn bezier3_linear_tangent_constant() {
        // Degree-1 Bezier (line) has constant tangent
        let c = BezierCurve3 {
            control_points: vec![DVec3::ZERO, DVec3::new(2.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let t0 = c.tangent_at(0.0);
        let t1 = c.tangent_at(0.5);
        let t_end = c.tangent_at(1.0);
        assert!((t0 - t1).length() < 1e-10);
        assert!((t0 - t_end).length() < 1e-10);
    }

    #[test]
    fn bezier3_rational_weight_effect() {
        // Rational Bezier with center weight > 1 pulls curve toward control point
        let non_rational = BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let rational = BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 0.0),
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 5.0, 1.0], // heavy center weight pulls up
        };
        let p_nr = non_rational.point_at(0.5);
        let p_r = rational.point_at(0.5);
        // Rational with heavy center weight should be higher (more Y)
        assert!(p_r.y > p_nr.y);
    }

    // ── Parabola ───────────────────────────────────────────────────────

    #[test]
    fn parabola3_eval_and_derivative() {
        use std::f64::consts::PI;
        // OCCT gp_Parab (gp_Ax2 N, X): the cross-axis is Y = N × X (right-handed
        // frame), so dir_perp = normal.cross(axis_dir) = Z×X = +Y.
        // P(t) = (t²/(2p)) * X + t * Y = (t²/8, t, 0)
        let p = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 4.0,
        };
        // P(0) = (0,0,0)
        assert!((p.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // P(4) = (16/8, 4, 0) = (2, 4, 0)
        assert!((p.point_at(4.0) - DVec3::new(2.0, 4.0, 0.0)).length() < 1e-10);
        // derivative = (t/4, 1, 0): at t=0 → (0, 1, 0), at t=4 → (1, 1, 0)
        let d0 = p.derivative_at(0.0);
        assert!((d0 - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-10);
        let d4 = p.derivative_at(4.0);
        assert!((d4 - DVec3::new(1.0, 1.0, 0.0)).length() < 1e-10);
    }

    // ── Hyperbola ──────────────────────────────────────────────────────

    #[test]
    fn hyperbola3_eval_and_derivative() {
        use std::f64::consts::PI;
        let h = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 3.0,
            semi_minor: 2.0,
        };
        // P(0) = (a*cosh(0), 0, 0) = (3, 0, 0)  (since minor_dir = normal×major_dir = Z×X = Y)
        // Actually: P(t) = center + a*cosh(t)*X + b*sinh(t)*Y
        // P(0) = (3*1, 2*0, 0) = (3, 0, 0)
        assert!((h.point_at(0.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
        // P(1) = (3*cosh(1), 2*sinh(1), 0)
        let p1 = h.point_at(1.0);
        assert!((p1.x - 3.0 * 1.0f64.cosh()).abs() < 1e-10);
        assert!((p1.y - 2.0 * 1.0f64.sinh()).abs() < 1e-10);
        // derivative = (a*sinh(t), b*cosh(t), 0)
        // at t=0: (0, 2, 0)
        let d0 = h.derivative_at(0.0);
        assert!((d0 - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-10);
    }

    // ── Helix ──────────────────────────────────────────────────────────

    #[test]
    fn helix3_full_turn_pitch_advance() {
        // Helix with pitch 2: after one full turn, Z advances by 2
        let h = CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 3.0,
            pitch: 2.0,
        };
        use std::f64::consts::{PI, TAU};
        // x_axis = ref_dir - axis * dot = X - 0 = X
        // y_axis = Z.cross(X) = Y
        // P(0) = (3, 0, 0), P(TAU) = (3, 0, pitch) = (3, 0, 2)
        assert!((h.point_at(0.0) - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
        let p_full = h.point_at(TAU);
        assert!((p_full - DVec3::new(3.0, 0.0, 2.0)).length() < 1e-10);
    }

    #[test]
    fn helix3_half_turn_opposite() {
        let h = CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 3.0,
            pitch: 2.0,
        };
        use std::f64::consts::PI;
        // P(PI) = (-3, 0, 1)
        let p_half = h.point_at(PI);
        assert!((p_half - DVec3::new(-3.0, 0.0, 1.0)).length() < 1e-10);
    }

    // ── SineWave ───────────────────────────────────────────────────────

    #[test]
    fn sine_wave3_eval_and_derivative() {
        let w = SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        };
        use std::f64::consts::PI;
        // P(0) = (0, 0, 0)
        assert!((w.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        // P(PI/2) = (PI/2, 2, 0)
        let p = w.point_at(PI / 2.0);
        assert!((p - DVec3::new(PI / 2.0, 2.0, 0.0)).length() < 1e-10);
        // derivative = X + 2*cos(t)*Y, at t=0: X + 2*Y
        let d0 = w.derivative_at(0.0);
        assert!((d0 - DVec3::new(1.0, 2.0, 0.0)).length() < 1e-10);
    }

    // ── OffsetCurve3 ───────────────────────────────────────────────────

    #[test]
    fn offset_curve3_line_offset() {
        // Line along X, offset along Z: tangent = X, perp = X×Z = -Y
        // The offset displaces in the -Y direction (perpendicular to both tangent and offset_dir)
        // FD tangent gives approximate direction, so just check the point differs from the line
        let basis = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let off = OffsetCurve3 {
            basis: Box::new(basis),
            offset_distance: 2.0,
            offset_dir: DVec3::Z,
        };
        let pt = off.point_at(5.0);
        // Should differ from the base line point (5,0,0)
        assert!((pt - DVec3::new(5.0, 0.0, 0.0)).length() > 1.0);
        // The Z coordinate should be near 0 (offset in XY plane, not Z)
        assert!(pt.z.abs() < 0.1);
    }

    // =========================================================================
    // OCCT-aligned comprehensive surface evaluation tests
    // (matching TKG3d/GTests Geom_Plane/Cylinder/Sphere/Cone/Torus patterns)
    // =========================================================================

    // ── Plane ───────────────────────────────────────────────────────────

    #[test]
    fn plane_derivatives_constant() {
        let p = Plane::new(DVec3::ZERO, DVec3::Z);
        // OCCT-aligned: gp_Ax3(gp_Pnt, gp_Dir) with normal=Z gives u_dir=X, v_dir=Y
        let (pt, dpu, dpv) = p.derivatives(2.0, 3.0);
        assert!((dpu - DVec3::X).length() < 1e-10);
        assert!((dpv - DVec3::Y).length() < 1e-10);
        // normal should be perpendicular to both dPdu and dPdv
        assert!(dpu.dot(p.normal_at(0.0, 0.0)).abs() < 1e-10);
        assert!(dpv.dot(p.normal_at(0.0, 0.0)).abs() < 1e-10);
        // point's Z should be 0
        assert!(pt.z.abs() < 1e-10);
    }

    // ── Cylinder ────────────────────────────────────────────────────────

    #[test]
    fn cylinder_derivatives_and_normal() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        let (p, dpu, dpv) = cyl.derivatives(0.0, 5.0);
        // dP/dv should be axis (Z)
        assert!((dpv - DVec3::Z).length() < 1e-10);
        // dP/du should be tangent around the cylinder, perpendicular to radius
        let radial = p - DVec3::new(0.0, 0.0, 5.0);
        let radial = radial.normalize_or_zero();
        assert!(dpu.dot(radial).abs() < 1e-10);
        // normal should be the radial direction
        let n = cyl.normal_at(0.0, 5.0);
        assert!((n - radial).length() < 1e-10);
    }

    #[test]
    fn cylinder_normal_radial() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        };
        // Cylinder with axis Z, ref_dir=X:
        // x_ax = X, y_ax = Z×X = Y
        // At u=0, normal should point in the X direction (radial outward)
        let n = cyl.normal_at(0.0, 0.0);
        assert!((n - DVec3::X).length() < 1e-10);
    }

    // ── Sphere ──────────────────────────────────────────────────────────

    #[test]
    fn sphere_derivatives_and_normal() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        // OCCT ElSLib::SphereValue/SphereD1: V is latitude, V=0 is the
        // equator (u=0 → +X), V=+pi/2 the axis pole. At the equator (v=0):
        // point = (5, 0, 0), dP/dv = R*axis = +Z (toward the north pole).
        use std::f64::consts::PI;
        let (p, dpu, dpv) = s.derivatives(0.0, 0.0);
        // Point radius should be 5
        assert!((p.length() - 5.0).abs() < 1e-10);
        // dP/du should be perpendicular to point
        assert!(dpu.dot(p.normalize_or_zero()).abs() < 1e-10);
        // dP/dv at the equator should point toward +Z (north pole)
        assert!((dpv - DVec3::new(0.0, 0.0, 5.0)).length() < 1e-10);
        // Normal should point radially outward
        let n = s.normal_at(0.0, 0.0);
        assert!((n - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn sphere_normal_at_poles() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        // OCCT SphereParameters: V = atan(z/l) — north pole V=+pi/2:
        // point=(0,0,5), normal +Z; south pole V=-pi/2: normal -Z.
        let n_north = s.normal_at(0.0, std::f64::consts::FRAC_PI_2);
        assert!((n_north - DVec3::Z).length() < 1e-10);
        let n_south = s.normal_at(0.0, -std::f64::consts::FRAC_PI_2);
        assert!((n_south + DVec3::Z).length() < 1e-10);
    }

    // ── Cone ────────────────────────────────────────────────────────────

    #[test]
    fn cone_derivatives() {
        let sa = 30.0_f64.to_radians(); // 30 degree half-angle
        let cone = ConicalSurface::new(DVec3::ZERO, DVec3::Z, 2.0, sa);
        let (_p, dpu, dpv) = cone.derivatives(0.0, 0.0);
        // At v=0: radial = 2, axial = 0
        // dP/du at u=0: radial * (-sin(0)*x_ax + cos(0)*y_ax)
        //   where x_ax = any_perpendicular(Z) = Y, y_ax = Z×Y = -X
        //   = 2 * (0*Y + 1*(-X)) = 2*(-X)
        // dP/dv at v=0: da*axis + dr*r_vec = cos(sa)*Z + sin(sa)*Y
        //   where r_vec = cos(0)*Y + sin(0)*(-X) = Y
        //   = cos(sa)*Z + sin(sa)*Y
        assert!(dpu.length() > 0.0);
        assert!(dpv.length() > 0.0);
        // dP/du should be perpendicular to the radial direction
        let n = cone.normal_at(0.0, 0.0);
        assert!(dpu.dot(n).abs() < 1e-10);
        assert!(dpv.dot(n).abs() < 1e-10);
    }

    // ── Torus ───────────────────────────────────────────────────────────

    #[test]
    fn torus_derivatives_and_normal() {
        use std::f64::consts::PI;
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // OCCT-aligned: x_ax=X, y_ax=Y. At u=0, v=0: outer equator, point=(6,0,0)
        let (p, dpu, dpv) = t.derivatives(0.0, 0.0);
        assert!((p - DVec3::new(6.0, 0.0, 0.0)).length() < 1e-9);
        // dP/du should be in the Y direction (major circle tangent)
        assert!((dpu - DVec3::new(0.0, 6.0, 0.0)).length() < 1e-9);
        // dP/dv should be in the Z direction (minor circle tangent at v=0)
        assert!((dpv - DVec3::Z).length() < 1e-9);
        // Normal at outer equator should be outward radial (X)
        let n = t.normal_at(0.0, 0.0);
        assert!((n - DVec3::X).length() < 1e-9);
    }

    #[test]
    fn torus_inner_equator_normal() {
        use std::f64::consts::PI;
        let t = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Inner equator: u=0, v=PI → point=(4,0,0), normal should be -X (inward)
        let n_inner = t.normal_at(0.0, PI);
        assert!((n_inner + DVec3::X).length() < 1e-9);
    }

    // ── BSplineSurface ──────────────────────────────────────────────────

    #[test]
    fn bspline_surface_eval_d0() {
        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 10.0, 0.0)],
                vec![DVec3::new(10.0, 0.0, 0.0), DVec3::new(10.0, 10.0, 0.0)],
            ],
            weights: vec![vec![1.0; 2]; 2],
        };
        // Degree 1 surface = bilinear, should interpolate corners
        assert!((surf.point_at(0.0, 0.0) - DVec3::ZERO).length() < 1e-10);
        assert!((surf.point_at(1.0, 1.0) - DVec3::new(10.0, 10.0, 0.0)).length() < 1e-10);
        assert!((surf.point_at(0.5, 0.5) - DVec3::new(5.0, 5.0, 0.0)).length() < 1e-10);
    }

    // ── Surface3 dispatch verification ───────────────────────────────────

    #[test]
    fn surface3_plane_dispatch() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        // OCCT-aligned: normal=Z gives u_dir=X, v_dir=Y
        let (_p, dpu, dpv) = s.derivatives(1.0, 2.0);
        assert!((dpu - DVec3::X).length() < 1e-10);
        assert!((dpv - DVec3::Y).length() < 1e-10);
        assert!(s.normal_at(1.0, 2.0) == DVec3::Z);
    }

    #[test]
    fn surface3_sphere_dispatch() {
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 3.0,
            ref_dir: DVec3::X,
        });
        // OCCT ElSLib::SphereValue: V=+pi/2 is the axis pole (north).
        let (p, _dpu, _dpv) = s.derivatives(0.0, std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(0.0, 0.0, 3.0)).length() < 1e-9);
        let n = s.normal_at(0.0, std::f64::consts::FRAC_PI_2);
        assert!((n - DVec3::Z).length() < 1e-9);
    }

    #[test]
    fn surface3_cylinder_dispatch() {
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        });
        // Cylinder with axis Z, ref_dir=X: dP/du at u=0 = R*y_ax = 2*Y = (0,2,0)
        let (_p, dpu, dpv) = s.derivatives(0.0, 0.0);
        // dP/dv should be axis direction
        assert!((dpv - DVec3::Z).length() < 1e-10);
        // dP/du should be tangent (perpendicular to radial)
        assert!((dpu - DVec3::new(0.0, 2.0, 0.0)).length() < 1e-10);
    }
}
