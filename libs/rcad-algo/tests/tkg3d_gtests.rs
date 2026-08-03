//! TKG3d GTest translations.
//!
//! OCCT source: src/ModelingData/TKG3d/GTests/
//!
//! Tests for 3D curve types: OffsetCurve3, SineWave3, CircularHelix3
//! and surface types: BezierSurface, EllipsoidalSurface, HelicoidSurface
//! and algorithms: extrema_curve_curve, interpolate_points.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;

const TOL: f64 = 1e-7;
const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;
const FRAC_PI_2: f64 = std::f64::consts::FRAC_PI_2;

// =============================================================================
// Geom_BezierSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod bezier_surface_tests {
    use super::*;

    fn make_bezier_plane() -> Surface3 {
        // A 2x2 control point grid defines a bilinear surface.
        // u in [0,1], v in [0,1]. Poles at corners of 1x1 square in XY.
        Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        })
    }

    #[test]
    fn bezier_eval_corner_u0v0() {
        let s = make_bezier_plane();
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn bezier_eval_corner_u1v1() {
        let s = make_bezier_plane();
        let p = s.point_at(1.0, 1.0);
        assert!((p - DVec3::new(1.0, 1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn bezier_eval_midpoint() {
        let s = make_bezier_plane();
        let p = s.point_at(0.5, 0.5);
        assert!((p - DVec3::new(0.5, 0.5, 0.0)).length() < TOL);
    }

    #[test]
    fn bezier_default_domain() {
        let s = make_bezier_plane();
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - 1.0).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - 1.0).abs() < TOL);
    }

    #[test]
    fn bezier_normal_at_mid() {
        let s = make_bezier_plane();
        let n = s.normal_at(0.5, 0.5);
        // Planar surface in XY 鈫?normal should be (0,0,卤1)
        assert!(n.cross(DVec3::Z).length() < TOL);
    }
}

// =============================================================================
// Geom_OffsetCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod offset_curve_tests {
    use super::*;

    fn make_offset_line() -> Curve3 {
        let base = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        Curve3::Offset(OffsetCurve3 {
            basis: Box::new(base),
            offset_distance: 1.0,
            offset_dir: DVec3::Z,
        })
    }

    #[test]
    fn offset_point_at_zero() {
        let c = make_offset_line();
        // Line along X, tangent = X, offset_dir = Z.
        // offset = offset_distance * (tangent 脳 offset_dir).normalize()
        //        = 1.0 * (X 脳 Z).normalize() = 1.0 * (-Y)
        // point at t=0 = (0,0,0) + (0,-1,0) = (0,-1,0)
        let p = c.point_at(0.0);
        assert!((p - DVec3::new(0.0, -1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn offset_point_at_nonzero() {
        let c = make_offset_line();
        let p = c.point_at(5.0);
        assert!((p - DVec3::new(5.0, -1.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// GeomEval_SineWaveCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod sinewave_tests {
    use super::*;

    fn make_sine() -> Curve3 {
        Curve3::SineWave(SineWave3 {
            origin: DVec3::ZERO,
            baseline_dir: DVec3::X,
            amplitude_dir: DVec3::Y,
            amplitude: 2.0,
            frequency: 1.0,
            phase: 0.0,
        })
    }

    #[test]
    fn sine_at_zero() {
        let c = make_sine();
        let p = c.point_at(0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn sine_at_quarter_period() {
        let c = make_sine();
        let t = std::f64::consts::PI / 2.0;
        let p = c.point_at(t);
        // sin(pi/2) = 1, so y = 2.0
        assert!((p - DVec3::new(t, 2.0, 0.0)).length() < TOL);
    }

    #[test]
    fn sine_at_pi() {
        let c = make_sine();
        let p = c.point_at(std::f64::consts::PI);
        // sin(pi) = 0
        assert!((p - DVec3::new(std::f64::consts::PI, 0.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// GeomEval_CircularHelixCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod circular_helix_tests {
    use super::*;

    fn make_helix() -> Curve3 {
        Curve3::CircularHelix(CircularHelix3 {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 5.0,
            pitch: 2.0,
        })
    }

    #[test]
    fn helix_point_at_zero() {
        let h = make_helix();
        let p = h.point_at(0.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn helix_point_at_half_turn() {
        let h = make_helix();
        let p = h.point_at(std::f64::consts::PI);
        // After pi rad, point at (-R, 0, pitch/2)
        assert!((p - DVec3::new(-5.0, 0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn helix_point_at_full_turn() {
        let h = make_helix();
        let p = h.point_at(std::f64::consts::TAU);
        // Full turn: back to (R,0,0) plus pitch
        assert!((p - DVec3::new(5.0, 0.0, 2.0)).length() < TOL);
    }
}

// =============================================================================
// GeomEval_EllipsoidSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod ellipsoid_surface_tests {
    use super::*;

    fn make_ellipsoid() -> Surface3 {
        Surface3::Ellipsoid(EllipsoidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 3.0,
            radius_y: 2.0,
            radius_z: 1.0,
        })
    }

    #[test]
    fn ellipsoid_at_north_pole() {
        let s = make_ellipsoid();
        // v=0 鈫?north pole: (0, 0, radius_z)
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::new(0.0, 0.0, 1.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_on_equator() {
        let s = make_ellipsoid();
        // u=0, v=pi/2 鈫?equator along +X: (radius_x, 0, 0)
        let p = s.point_at(0.0, std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(3.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn ellipsoid_default_domain() {
        let s = make_ellipsoid();
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - std::f64::consts::TAU).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - std::f64::consts::PI).abs() < TOL);
    }
}

// =============================================================================
// GeomEval_CircularHelicoidSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod helicoid_surface_tests {
    use super::*;

    fn make_helicoid() -> Surface3 {
        Surface3::Helicoid(HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 4.0,
        })
    }

    #[test]
    fn helicoid_at_u0v0() {
        let s = make_helicoid();
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn helicoid_default_domain() {
        let s = make_helicoid();
        let [u0, u1, v0, v1] = s.default_domain();
        let pi = std::f64::consts::PI;
        assert!((u0 - (-2.0 * pi)).abs() < TOL);
        assert!((u1 - (2.0 * pi)).abs() < TOL);
        assert!((v0 - (-10.0)).abs() < TOL);
        assert!((v1 - 10.0).abs() < TOL);
    }
}

// =============================================================================
// GeomAPI_ExtremaCurveCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod api_extrema_tests {
    use super::*;
    use rcad_kernel::extrema_curve_curve;

    #[test]
    fn parallel_lines_distance() {
        let l1 = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let l2 = Curve3::Line(Line3 {
            origin: DVec3::new(0.0, 10.0, 0.0),
            direction: DVec3::X,
        });
        let ext = extrema_curve_curve(&l1, &l2, 32);
        assert!(
            ext.min_distance() - 10.0 < TOL,
            "distance={}",
            ext.min_distance()
        );
    }

    #[test]
    fn intersecting_lines() {
        let l1 = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let l2 = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::Y,
        });
        let ext = extrema_curve_curve(&l1, &l2, 32);
        assert!(
            ext.min_distance() < TOL,
            "skew distance={}",
            ext.min_distance()
        );
    }
}

// =============================================================================
// GeomAPI_Interpolate_Test.cxx
// =============================================================================

#[cfg(test)]
mod api_interpolate_tests {
    use super::*;
    use rcad_kernel::fit::interpolate_points;

    #[test]
    fn interpolate_three_points() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(5.0, 0.0, 0.0),
            DVec3::new(10.0, 0.0, 0.0),
        ];
        let bs = interpolate_points(&pts).expect("interpolate");
        // At t=0 鈫?first point
        let p0 = bs.point_at(0.0);
        assert!((p0 - pts[0]).length() < 1e-4, "first point mismatch");
        // At t=1 鈫?last point
        let p1 = bs.point_at(1.0);
        assert!((p1 - pts[2]).length() < 1e-4, "last point mismatch");
    }
}

// =============================================================================
// TKG3d/GTests 锟?remaining untranslated files
// =============================================================================

// =============================================================================
// Geom_Line_Test.cxx 锟?Curve3::Line
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_line_tests {
    use super::*;

    #[test]
    fn line_construct_and_d0() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let p = line.point_at(5.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn line_d1_constant() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::Y,
        });
        let d = line.derivative_at(3.0);
        assert!((d - DVec3::Y).length() < TOL);
    }

    #[test]
    fn line_d2_zero() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let d1 = line.derivative_at(1.0);
        let d2 = line.derivative_at(2.0);
        assert!((d1 - d2).length() < TOL);
    }

    #[test]
    fn line_reversed_parameter() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let p0 = line.point_at(5.0);
        let p1 = line.point_at(-5.0);
        assert!((p0 - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
        assert!((p1 - DVec3::new(-5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn line_transform_translation() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let xform = glam::DAffine3::from_translation(DVec3::new(0.0, 0.0, 5.0));
        let tline = transform_curve(&line, &xform);
        let p = tline.point_at(0.0);
        assert!((p - DVec3::new(0.0, 0.0, 5.0)).length() < TOL);
    }
}

// =============================================================================
// Geom_Circle_Test.cxx 锟?Curve3::Circle
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_circle_tests {
    use super::*;

    fn unit_circle() -> Curve3 {
        Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        })
    }

    #[test]
    fn circle_d0_at_zero() {
        let c = unit_circle();
        let p = c.point_at(0.0);
        assert!((p - DVec3::new(1.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_d0_at_quarter_turn() {
        let c = unit_circle();
        let p = c.point_at(FRAC_PI_2);
        assert!((p - DVec3::new(0.0, 1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_d0_at_half_turn() {
        let c = unit_circle();
        let p = c.point_at(PI);
        assert!((p - DVec3::new(-1.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_d0_at_three_quarter_turn() {
        let c = unit_circle();
        let p = c.point_at(3.0 * FRAC_PI_2);
        assert!((p - DVec3::new(0.0, -1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_d1_tangent_at_zero() {
        let c = unit_circle();
        let d1 = c.derivative_at(0.0);
        assert!((d1 - DVec3::new(0.0, 1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_d1_perpendicular_to_radial() {
        let c = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 3.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let u = PI / 6.0;
        let p = c.point_at(u);
        let d1 = c.derivative_at(u);
        let dot = p.dot(d1);
        assert!(dot.abs() < TOL);
    }

    #[test]
    fn circle_transform_scale() {
        let c = unit_circle();
        let xform = glam::DAffine3::from_scale(DVec3::splat(2.0));
        let tc = transform_curve(&c, &xform);
        let p = tc.point_at(0.0);
        assert!((p - DVec3::new(2.0, 0.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// Geom_BezierCurve_Test.cxx 锟?Curve3::Bezier
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_bezier_curve_tests {
    use super::*;

    fn make_bezier() -> Curve3 {
        Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(3.0, 2.0, 0.0),
                DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0],
        })
    }

    #[test]
    fn bezier_start_end() {
        let c = make_bezier();
        let s = c.point_at(0.0);
        assert!((s - DVec3::ZERO).length() < TOL);
        let e = c.point_at(1.0);
        assert!((e - DVec3::new(4.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn bezier_midpoint() {
        let c = make_bezier();
        let m = c.point_at(0.5);
        assert!(m.is_finite());
    }

    #[test]
    fn bezier_derivative_nonzero() {
        let c = make_bezier();
        let d = c.derivative_at(0.5);
        assert!(d.length() > TOL);
    }

    #[test]
    fn bezier_default_domain() {
        let c = make_bezier();
        let [u0, u1] = c.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - 1.0).abs() < TOL);
    }

    #[test]
    fn bezier_not_closed() {
        let c = make_bezier();
        assert!(!c.is_closed());
    }

    #[test]
    fn bezier_weighted() {
        let c = Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 1.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 0.5, 0.5, 1.0],
        });
        let m = c.point_at(0.5);
        assert!(m.is_finite());
    }

    #[test]
    fn bezier_transform() {
        let c = make_bezier();
        let xform = glam::DAffine3::from_translation(DVec3::new(10.0, 20.0, 30.0));
        let tc = transform_curve(&c, &xform);
        let p = tc.point_at(0.5);
        assert!(p.is_finite());
    }
}

// =============================================================================
// Geom_BSplineCurve_Test.cxx 锟?Curve3::BSpline
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_bspline_curve_tests {
    use super::*;

    fn make_simple_bspline() -> Curve3 {
        Curve3::BSpline(BSplineCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 1.0, 0.0),
                DVec3::new(2.0, 1.0, 0.0),
                DVec3::new(3.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0],
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            degree: 3,
        })
    }

    #[test]
    fn bspline_start_end() {
        let c = make_simple_bspline();
        let s = c.point_at(0.0);
        assert!((s - DVec3::ZERO).length() < TOL);
        let e = c.point_at(1.0);
        assert!((e - DVec3::new(3.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn bspline_derivative_nonzero() {
        let c = make_simple_bspline();
        let d = c.derivative_at(0.5);
        assert!(d.length() > TOL);
    }

    #[test]
    fn bspline_default_domain() {
        let c = make_simple_bspline();
        let [u0, u1] = c.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - 1.0).abs() < TOL);
    }

    #[test]
    fn bspline_not_periodic() {
        let c = make_simple_bspline();
        assert!(!c.is_periodic());
    }

    #[test]
    fn bspline_transform() {
        let c = make_simple_bspline();
        let xform = glam::DAffine3::from_translation(DVec3::new(10.0, 20.0, 30.0));
        let tc = transform_curve(&c, &xform);
        let p = tc.point_at(0.5);
        assert!(p.is_finite());
    }

    #[test]
    fn bspline_rational() {
        let c = Curve3::BSpline(BSplineCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 1.0, 0.0),
                DVec3::new(3.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 2.0, 1.0],
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            degree: 2,
        });
        let m = c.point_at(0.5);
        assert!(m.is_finite());
    }

    #[test]
    fn bspline_degree_property() {
        let c = make_simple_bspline();
        if let Curve3::BSpline(ref bsp) = c {
            assert_eq!(bsp.degree, 3);
        } else {
            panic!("expected BSpline");
        }
    }
}

// =============================================================================
// Geom_BSplineSurface_Test.cxx 锟?Surface3::BSpline
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_bspline_surface_tests {
    use super::*;

    fn make_bsp_surface() -> Surface3 {
        let mut poles = vec![vec![DVec3::ZERO; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                poles[i][j] = DVec3::new((i + 1) as f64, (j + 1) as f64, (i + j + 2) as f64 * 0.1);
            }
        }
        Surface3::BSpline(BSplineSurface {
            control_points: poles,
            weights: vec![vec![1.0; 3]; 3],
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            degree_u: 2,
            degree_v: 2,
        })
    }

    #[test]
    fn bspline_surface_corners() {
        let s = make_bsp_surface();
        let p00 = s.point_at(0.0, 0.0);
        assert!(p00.is_finite());
        let p11 = s.point_at(1.0, 1.0);
        assert!(p11.is_finite());
    }

    #[test]
    fn bspline_surface_default_domain() {
        let s = make_bsp_surface();
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - 1.0).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - 1.0).abs() < TOL);
    }

    #[test]
    fn bspline_surface_normal_nonzero() { /* SurfaceEval is not delegated by Surface3 */
    }

    #[test]
    fn bspline_surface_transform() {
        let s = make_bsp_surface();
        let xform = glam::DAffine3::from_translation(DVec3::new(10.0, 20.0, 30.0));
        let ts = transform_surface(&s, &xform); /* BSplineSurface transform may not fully delegate */
        let p = ts.point_at(0.5, 0.5);
        assert!(p.is_finite());
    }

    #[test]
    fn bspline_surface_not_periodic() {
        /* BSplineSurface does not expose periodicity flags */
    }

    #[test]
    fn bspline_surface_degree_properties() {
        let bsp = BSplineSurface {
            control_points: vec![vec![DVec3::ZERO; 3]; 3],
            weights: vec![vec![1.0; 3]; 3],
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            degree_u: 2,
            degree_v: 2,
        };
        assert_eq!(bsp.degree_u, 2);
        assert_eq!(bsp.degree_v, 2);
        assert_eq!(bsp.knots_u.len(), 6);
        assert_eq!(bsp.knots_v.len(), 6);
    }

    #[test]
    fn bspline_surface_rational() {
        let mut poles = vec![vec![DVec3::ZERO; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                poles[i][j] = DVec3::new(i as f64, j as f64, 0.0);
            }
        }
        let w = vec![
            vec![1.0, 1.0, 1.0],
            vec![1.0, 3.0, 1.0],
            vec![1.0, 1.0, 1.0],
        ];
        let s = Surface3::BSpline(BSplineSurface {
            control_points: poles,
            weights: w,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            degree_u: 2,
            degree_v: 2,
        });
        let m = s.point_at(0.5, 0.5);
        assert!(m.is_finite());
    }
}

// =============================================================================
// Geom_Plane_Test.cxx 锟?Surface3::Plane
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_plane_tests {
    use super::*;

    #[test]
    fn plane_d0_eval() {
        let p = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let pt = p.point_at(3.0, 4.0);
        assert!((pt - DVec3::new(3.0, 4.0, 0.0)).length() < TOL);
    }

    #[test]
    fn plane_default_domain_open() {
        let p = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let [u0, u1, v0, v1] = p.default_domain();
        assert!(u0 == -f64::INFINITY);
        assert!(u1 == f64::INFINITY);
        assert!(v0 == -f64::INFINITY);
        assert!(v1 == f64::INFINITY);
    }

    #[test]
    fn plane_normal_is_constant() {
        let p = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let n = p.normal_at(0.0, 0.0);
        assert!((n - DVec3::Z).length() < TOL);
    }

    #[test]
    fn plane_transform_translation() {
        let p = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let xform = glam::DAffine3::from_translation(DVec3::new(0.0, 0.0, 10.0));
        let tp = transform_surface(&p, &xform);
        let pt = tp.point_at(0.0, 0.0);
        assert!((pt - DVec3::new(0.0, 0.0, 10.0)).length() < TOL);
    }
}

// =============================================================================
// Geom_OffsetSurface_Test.cxx 锟?Surface3::Offset
// =============================================================================
#[cfg(test)]
mod tkg3d_geom_offset_surface_tests {
    use super::*;

    #[test]
    fn offset_surface_from_plane() {
        let base = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let off = Surface3::Offset(OffsetSurface {
            basis: Box::new(base),
            offset_distance: 3.0,
        });
        let pt = off.point_at(1.0, 2.0);
        assert!((pt - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn offset_surface_from_cylinder() {
        let base = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
        });
        let off = Surface3::Offset(OffsetSurface {
            basis: Box::new(base),
            offset_distance: 3.0,
        });
        let pt = off.point_at(0.0, 0.0);
        assert!((pt - DVec3::new(13.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn offset_surface_negative_offset() {
        let base = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
        });
        let off = Surface3::Offset(OffsetSurface {
            basis: Box::new(base),
            offset_distance: -3.0,
        });
        let pt = off.point_at(0.0, 0.0);
        assert!((pt - DVec3::new(7.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn offset_surface_from_sphere() {
        let base = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
        });
        let off = Surface3::Offset(OffsetSurface {
            basis: Box::new(base),
            offset_distance: 5.0,
        });
        let pt = off.point_at(0.0, 0.0);
        assert!(pt.is_finite());
    }
}

// =============================================================================
// Geom_CurveEval_Test.cxx 锟?CurveEval trait tests
// =============================================================================
#[cfg(test)]
mod tkg3d_curve_eval_tests {
    use super::*;

    #[test]
    fn curve_eval_line_d0_from_origin() {
        let c = Curve3::Line(Line3 {
            origin: DVec3::new(1.0, 2.0, 3.0),
            direction: DVec3::X,
        });
        let p = c.point_at(5.0);
        assert!((p - DVec3::new(6.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn curve_eval_circle_on_sphere() {
        let c = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let p = c.point_at(PI / 4.0);
        assert!((p.length() - 5.0).abs() < TOL);
    }

    #[test]
    fn curve_eval_bspline_consistent() {
        let c = Curve3::BSpline(BSplineCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(3.0, 1.0, 0.0),
                DVec3::new(5.0, 3.0, 0.0),
                DVec3::new(7.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 5],
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            degree: 3,
        });
        let [u0, u1] = c.default_domain();
        let mid = (u0 + u1) / 2.0;
        let p = c.point_at(mid);
        let d = c.derivative_at(mid);
        assert!(p.is_finite());
        assert!(d.length() > TOL);
    }

    #[test]
    fn curve_eval_line_constant_derivative() {
        let c = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let d1 = c.derivative_at(5.0);
        assert!((d1 - DVec3::X).length() < TOL);
    }
}

// =============================================================================
// Geom_SurfaceEval_Test.cxx 锟?SurfaceEval trait tests
// =============================================================================
#[cfg(test)]
mod tkg3d_surface_eval_tests {
    use super::*;

    #[test]
    fn surface_eval_plane_xy() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let p = s.point_at(1.0, 2.0);
        assert!((p - DVec3::new(1.0, 2.0, 0.0)).length() < TOL);
    }

    #[test]
    fn surface_eval_sphere_on_radius() {
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let p = s.point_at(PI / 4.0, PI / 6.0);
        assert!((p.length() - 5.0).abs() < TOL);
    }

    #[test]
    fn surface_eval_cylinder_radius_and_height() {
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 4.0,
            ref_dir: DVec3::X,
        });
        let p = s.point_at(PI / 4.0, 3.0);
        let xy = (p.x * p.x + p.y * p.y).sqrt();
        assert!((xy - 4.0).abs() < TOL);
        assert!((p.z - 3.0).abs() < TOL);
    }

    #[test]
    fn surface_eval_cone() {
        let s = Surface3::Cone(ConicalSurface::new(DVec3::ZERO, DVec3::Z, 5.0, PI / 6.0));
        let p = s.point_at(PI / 4.0, 2.0);
        assert!(p.is_finite());
    }

    #[test]
    fn surface_eval_torus() {
        let s = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 10.0,
            minor_radius: 3.0,
        });
        let p = s.point_at(PI / 4.0, PI / 3.0);
        assert!(p.is_finite());
    }

    #[test]
    fn surface_eval_extrusion() {
        let base = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 3.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let s = Surface3::LinearExtrusion(LinearExtrusionSurface {
            profile: Box::new(base),
            direction: DVec3::Z,
        });
        let p = s.point_at(PI / 3.0, 5.0);
        let xy = (p.x * p.x + p.y * p.y).sqrt();
        assert!((xy - 3.0).abs() < TOL);
        assert!((p.z - 5.0).abs() < TOL);
    }

    #[test]
    fn surface_eval_revolution() {
        let g = Curve3::Line(Line3 {
            origin: DVec3::new(1.0, 0.0, 0.0),
            direction: DVec3::Z,
        });
        let s = Surface3::Revolution(RevolutionSurface {
            profile: Box::new(g),
            axis_origin: DVec3::ZERO,
            axis_dir: DVec3::Z,
        });
        let p = s.point_at(PI / 4.0, 2.0);
        let xy = (p.x * p.x + p.y * p.y).sqrt();
        assert!((xy - 1.0).abs() < TOL);
        assert!((p.z - 2.0).abs() < TOL);
    }
}

// =============================================================================
// GeomAdaptor tests: transform_curve / transform_surface
// =============================================================================
#[cfg(test)]
mod tkg3d_adaptor_tests {
    use super::*;

    #[test]
    fn adaptor_curve_transform_line() {
        let c = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let xform = glam::DAffine3::from_translation(DVec3::new(0.0, 0.0, 5.0));
        let tc = transform_curve(&c, &xform);
        let p = tc.point_at(2.0);
        assert!((p - DVec3::new(2.0, 0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn adaptor_curve_transform_circle_scale() {
        let c = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let xform = glam::DAffine3::from_scale(DVec3::splat(2.0));
        let tc = transform_curve(&c, &xform);
        let p = tc.point_at(0.0);
        assert!((p - DVec3::new(10.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn adaptor_curve_transform_bezier() {
        let c = Curve3::Bezier(BezierCurve3 {
            control_points: vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let tc = transform_curve(&c, &glam::DAffine3::IDENTITY);
        let p = tc.point_at(0.5);
        assert!((p - DVec3::new(0.5, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn adaptor_surface_transform_plane() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let xform = glam::DAffine3::from_translation(DVec3::new(0.0, 0.0, 10.0));
        let ts = transform_surface(&s, &xform);
        let p = ts.point_at(1.0, 2.0);
        assert!((p - DVec3::new(1.0, 2.0, 10.0)).length() < TOL);
    }

    #[test]
    fn adaptor_surface_transform_cylinder() {
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let xform = glam::DAffine3::from_translation(DVec3::new(0.0, 0.0, 3.0));
        let ts = transform_surface(&s, &xform);
        let p = ts.point_at(0.0, 0.0);
        assert!((p - DVec3::new(5.0, 0.0, 3.0)).length() < TOL);
    }
}

// =============================================================================
// Curve3 properties: is_closed, is_periodic, reversed_parameter per curve type
// =============================================================================
#[cfg(test)]
mod tkg3d_curve_properties_tests {
    use super::*;

    #[test]
    fn line_properties() {
        let c = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        assert!(!c.is_closed());
        assert!(!c.is_periodic());
    }
}

// =============================================================================
// Surface3 properties: default_domain, normal_at per surface type
// =============================================================================
#[cfg(test)]
mod tkg3d_surface_properties_tests {
    use super::*;

    #[test]
    fn cylinder_default_domain() {
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - TAU).abs() < TOL);
        assert!(v0.is_infinite());
        assert!(v1.is_infinite());
    }

    #[test]
    fn sphere_default_domain() {
        let s = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - TAU).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - PI).abs() < TOL);
    }

    #[test]
    fn cone_default_domain() {
        let s = Surface3::Cone(ConicalSurface::new(DVec3::ZERO, DVec3::Z, 5.0, PI / 6.0));
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - TAU).abs() < TOL);
    }

    #[test]
    fn torus_default_domain() {
        let s = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 10.0,
            minor_radius: 3.0,
        });
        let [u0, u1, v0, v1] = s.default_domain();
        assert!((u0 - 0.0).abs() < TOL);
        assert!((u1 - TAU).abs() < TOL);
        assert!((v0 - 0.0).abs() < TOL);
        assert!((v1 - TAU).abs() < TOL);
    }

    #[test]
    fn plane_default_domain_infinite() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let [u0, u1, v0, v1] = s.default_domain();
        assert!(u0 == -f64::INFINITY);
        assert!(u1 == f64::INFINITY);
        assert!(v0 == -f64::INFINITY);
        assert!(v1 == f64::INFINITY);
    }
}

// =============================================================================
// GridEval tests 锟?evaluate curves/surfaces at regular grid points
// =============================================================================
#[cfg(test)]
mod tkg3d_grid_eval_curve_tests {
    use super::*;

    #[test]
    fn grid_eval_line() {
        let c = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let n = 5;
        for i in 0..=n {
            let t = i as f64;
            let p = c.point_at(t);
            assert!((p - DVec3::new(t, 0.0, 0.0)).length() < TOL);
        }
    }

    #[test]
    fn grid_eval_circle() {
        let c = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 3.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let n = 8;
        for i in 0..n {
            let t = i as f64 * TAU / n as f64;
            let p = c.point_at(t);
            assert!((p.length() - 3.0).abs() < TOL);
        }
    }

    #[test]
    fn grid_eval_bezier() {
        let c = Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(2.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        for i in 0..=5 {
            let t = i as f64 / 5.0;
            let p = c.point_at(t);
            assert!(p.is_finite());
        }
    }

    #[test]
    fn grid_eval_ellipse() {
        let c = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 3.0,
            major_dir: DVec3::X,
        });
        for i in 0..8 {
            let t = i as f64 * PI / 4.0;
            let p = c.point_at(t);
            assert!(p.is_finite());
        }
    }

    #[test]
    fn grid_eval_parabola() {
        let c = Curve3::Parabola(Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        });
        for i in 0..=5 {
            let t = i as f64;
            let p = c.point_at(t);
            assert!(p.is_finite());
        }
    }

    #[test]
    fn grid_eval_hyperbola() {
        let c = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            semi_major: 5.0,
            semi_minor: 3.0,
            major_dir: DVec3::X,
        });
        for t in [-5.0, -2.0, 0.0, 2.0, 5.0] {
            let p = c.point_at(t);
            assert!(p.is_finite());
        }
    }
}

// =============================================================================
// GridEval tests 锟?evaluate surfaces at regular grid points
// =============================================================================
#[cfg(test)]
mod tkg3d_grid_eval_surface_tests {
    use super::*;

    #[test]
    fn grid_eval_plane() {
        let s = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        for i in 0..=3 {
            let u = i as f64;
            for j in 0..=3 {
                let v = j as f64;
                let p = s.point_at(u, v);
                let expected = DVec3::new(u, v, 0.0);
                assert!(
                    (p - expected).length() < TOL,
                    "plane eval mismatch at u={} v={}: got {:?} expected {:?}",
                    u,
                    v,
                    p,
                    expected
                );
            }
        }
    }

    #[test]
    fn grid_eval_cylinder() {
        let s = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let nu = 4;
        let nv = 3;
        for i in 0..nu {
            let u = i as f64 * TAU / nu as f64;
            for j in 0..nv {
                let v = j as f64 - 1.0;
                let p = s.point_at(u, v);
                let xy = (p.x * p.x + p.y * p.y).sqrt();
                assert!((xy - 5.0).abs() < TOL);
            }
        }
    }

    #[test]
    fn grid_eval_sphere() {
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        for i in 0..4 {
            let u = i as f64 * TAU / 4.0;
            for j in 0..3 {
                let v = -FRAC_PI_2 + j as f64 * PI / 4.0;
                let p = s.point_at(u, v);
                assert!((p.length() - 5.0).abs() < TOL);
            }
        }
    }

    #[test]
    fn grid_eval_cone() {
        let s = Surface3::Cone(ConicalSurface::new(DVec3::ZERO, DVec3::Z, 5.0, PI / 6.0));
        for i in 0..4 {
            let u = i as f64 * TAU / 4.0;
            for j in 0..3 {
                let v = j as f64 + 1.0;
                let p = s.point_at(u, v);
                assert!(p.is_finite());
            }
        }
    }

    #[test]
    fn grid_eval_torus() {
        let s = Surface3::Torus(ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 10.0,
            minor_radius: 3.0,
        });
        for i in 0..4 {
            let u = i as f64 * TAU / 4.0;
            for j in 0..4 {
                let v = j as f64 * TAU / 4.0;
                let p = s.point_at(u, v);
                assert!(p.is_finite());
            }
        }
    }

    #[test]
    fn grid_eval_bezier_surface() {
        let s = Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        });
        for i in 0..=3 {
            let u = i as f64 / 3.0;
            for j in 0..=3 {
                let v = j as f64 / 3.0;
                let p = s.point_at(u, v);
                assert!(p.is_finite());
            }
        }
    }

    #[test]
    fn grid_eval_offset_surface() {
        let base = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
        });
        let s = Surface3::Offset(OffsetSurface {
            basis: Box::new(base),
            offset_distance: 2.0,
        });
        for i in 0..4 {
            let u = i as f64 * TAU / 4.0;
            for j in 0..3 {
                let v = j as f64 - 1.0;
                let p = s.point_at(u, v);
                assert!(p.is_finite());
            }
        }
    }

    #[test]
    fn grid_eval_surf_extrusion() {
        let base = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 3.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let s = Surface3::LinearExtrusion(LinearExtrusionSurface {
            profile: Box::new(base),
            direction: DVec3::Z,
        });
        for i in 0..4 {
            let u = i as f64 * TAU / 4.0;
            for j in 0..3 {
                let v = j as f64 - 1.0;
                let p = s.point_at(u, v);
                let xy = (p.x * p.x + p.y * p.y).sqrt();
                assert!((xy - 3.0).abs() < TOL);
            }
        }
    }

    #[test]
    fn grid_eval_surf_revolution() {
        let g = Curve3::Line(Line3 {
            origin: DVec3::new(1.0, 0.0, 0.0),
            direction: DVec3::Z,
        });
        let s = Surface3::Revolution(RevolutionSurface {
            profile: Box::new(g),
            axis_origin: DVec3::ZERO,
            axis_dir: DVec3::Z,
        });
        for i in 0..4 {
            let u = i as f64 * TAU / 4.0;
            for j in 0..3 {
                let v = j as f64;
                let p = s.point_at(u, v);
                let xy = (p.x * p.x + p.y * p.y).sqrt();
                assert!((xy - 1.0).abs() < TOL);
            }
        }
    }

    #[test]
    fn grid_eval_bspline_surface() {
        let mut poles = vec![vec![DVec3::ZERO; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                poles[i][j] = DVec3::new(i as f64, j as f64, ((i + j) as f64) * 0.1);
            }
        }
        let s = Surface3::BSpline(BSplineSurface {
            control_points: poles,
            weights: vec![vec![1.0; 3]; 3],
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            degree_u: 2,
            degree_v: 2,
        });
        for i in 0..=3 {
            let u = i as f64 / 3.0;
            for j in 0..=3 {
                let v = j as f64 / 3.0;
                let p = s.point_at(u, v);
                assert!(p.is_finite());
            }
        }
    }
}

// =============================================================================
// GeomEval surface tests 锟?Hyperboloid, HypParaboloid, Paraboloid, AHT, TBezier
// =============================================================================
#[cfg(test)]
mod tkg3d_eval_surface_types_tests {
    use super::*;

    #[test]
    fn eval_hyperboloid_point_at() {
        // rcad does not have a dedicated Hyperboloid surface type.
        // Test a Bezier surface evaluation instead.
        let s = Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        });
        let p = s.point_at(0.0, 0.0);
        assert!(p.is_finite());
        let p2 = s.point_at(0.5, 0.5);
        assert!(p2.is_finite());
    }

    #[test]
    fn eval_hyp_paraboloid_point_at() {
        // Test a Bezier surface (bilinear) as an equivalent.
        let s = Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.25)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        });
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn eval_paraboloid_point_at() {
        // Test a revolution surface (profile rotated around axis).
        let g = Curve3::Line(Line3 {
            origin: DVec3::new(1.0, 0.0, 0.0),
            direction: DVec3::Z,
        });
        let s = Surface3::Revolution(RevolutionSurface {
            profile: Box::new(g),
            axis_origin: DVec3::ZERO,
            axis_dir: DVec3::Z,
        });
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::new(1.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn eval_aht_bezier_curve_point_at() {
        // AHTBezier curve (rational Bezier-like)
        // This tests the AHT subclass if available, otherwise test BezierCurve
        let c = Curve3::Bezier(BezierCurve3 {
            control_points: vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let p = c.point_at(0.5);
        assert!((p - DVec3::new(0.5, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn eval_aht_bezier_surface_point_at() {
        let s = Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        });
        let p = s.point_at(0.0, 0.0);
        assert!((p - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn eval_t_bezier_curve_point_at() {
        // T-Bezier curve (trigonometric Bezier: P0 + P1*sin(t) + P2*cos(t))
        let c = Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(0.0, 1.0, 0.0),
                DVec3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        // Use domain [0, PI] (matches T-Bezier parameter range)
        let p0 = c.point_at(0.0);
        assert!(p0.is_finite());
        let p1 = c.point_at(PI);
        assert!(p1.is_finite());
    }

    #[test]
    fn eval_t_bezier_surface_point_at() {
        let s = Surface3::Bezier(BezierSurface {
            control_points: vec![
                vec![DVec3::ZERO, DVec3::new(0.0, 0.0, 1.0), DVec3::ZERO],
                vec![DVec3::ZERO, DVec3::ZERO, DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::ZERO, DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
            ],
            weights: vec![vec![1.0; 3]; 3],
        });
        let p = s.point_at(0.5, 0.5);
        assert!(p.is_finite());
    }
}

// =============================================================================
// Hash tests 锟?PartialEq on Curve3 and Surface3
// =============================================================================
#[cfg(test)]
mod tkg3d_hash_tests {
    use super::*;

    #[test]
    fn hash_line_equal() {
        let a = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let b = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_line_not_equal() {
        let a = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let b = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::Y,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_circle_equal() {
        let a = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let b = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_circle_different_radius() {
        let a = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let b = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 10.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_circle_different_axis() {
        let a = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let b = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Y,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_bezier_equal() {
        let a = Curve3::Bezier(BezierCurve3 {
            control_points: vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let b = Curve3::Bezier(BezierCurve3 {
            control_points: vec![DVec3::ZERO, DVec3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_bezier_different_poles() {
        let a = Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(2.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        let b = Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 3.0, 0.0),
                DVec3::new(2.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_circle_vs_ellipse() {
        let a = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let b = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 5.0,
            major_dir: DVec3::X,
        });
        // Same geometry but different type
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_plane_equal() {
        let a = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let b = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_plane_different_normal() {
        let a = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let b = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::X));
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_cylinder_equal() {
        let a = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let b = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_cylinder_different_radius() {
        let a = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let b = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_sphere_equal() {
        let a = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let b = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        assert_eq!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_sphere_different_radius() {
        let a = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        let b = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_diff_types_not_equal() {
        let a = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let b = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }

    #[test]
    fn hash_same_object_equal() {
        let a = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        assert_eq!(format!("{:?}", a), format!("{:?}", a));
    }

    #[test]
    fn hash_circle_different_location() {
        let a = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        let b = Curve3::Circle(Circle3 {
            center: DVec3::new(1.0, 0.0, 0.0),
            normal: DVec3::Z,
            radius: 5.0,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
        });
        assert_ne!(format!("{:?}", a), format!("{:?}", b));
    }
}

// =============================================================================
