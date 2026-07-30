//! TKG2d GTest translations.
//!
//! OCCT source: src/ModelingData/TKG2d/GTests/
//!
//! Tests for 2D curve types: Line, Circle, Ellipse, Parabola, Hyperbola,
//! Bezier, BSpline, Offset, Trimmed, SineWave, ArchimedeanSpiral,
//! CircleInvolute, LogarithmicSpiral, and API functions.

use glam::DVec2;
use rcad_kernel::geom::*;

const TOL: f64 = 1e-10;
const PI: f64 = std::f64::consts::PI;

// =============================================================================
// Geom2d_Line_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_line_tests {
    use super::*;

    fn make_line() -> Curve2d {
        Curve2d::Line(Line2d::new(Point2::ZERO, Vec2::X))
    }

    #[test]
    fn line_construct_from_point_dir() {
        let l = make_line();
        let (d, o) = match &l {
            Curve2d::Line(l) => (l.direction, l.origin),
            _ => unreachable!(),
        };
        assert!((d - Vec2::X).length() < TOL);
        assert!((o - Point2::ZERO).length() < TOL);
    }

    #[test]
    fn line_eval_d0() {
        let l = make_line();
        let p = l.point_at(5.0);
        assert!((p - Point2::new(5.0, 0.0)).length() < TOL);
    }

    #[test]
    fn line_d1_constant() {
        let l = make_line();
        let d1 = l.derivative_at(5.0);
        assert!((d1 - Vec2::X).length() < TOL);
    }

    #[test]
    fn line_diagonal_value() {
        let d = Vec2::new(1.0, 1.0).normalize();
        let l = Curve2d::Line(Line2d::new(Point2::ZERO, d));
        let p = l.point_at(2.0f64.sqrt());
        assert!((p - Point2::new(1.0, 1.0)).length() < TOL);
    }

    #[test]
    fn line_distance() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        // Point (5, 3) has perpendicular distance 3 from X-axis line
        let d = l.distance(Point2::new(5.0, 3.0));
        assert!((d - 3.0).abs() < TOL, "Distance should be 3, got {d}");
        // Point on the line has distance 0
        let d = l.distance(Point2::new(0.0, 0.0));
        assert!(d < TOL, "Point on line should have distance 0, got {d}");
    }

    #[test]
    fn line_not_closed_not_periodic() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        assert!(!l.is_closed(), "Line should not be closed");
        assert!(!l.is_periodic(), "Line should not be periodic");
    }

    #[test]
    fn line_reversed_parameter() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        assert!((l.reversed_parameter(5.0) - (-5.0)).abs() < TOL);
    }

    #[test]
    fn line_set_direction() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        let l = l.with_direction(Vec2::Y);
        assert!((l.direction - Vec2::Y).length() < TOL);
    }

    #[test]
    fn line_set_location() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        let l = l.with_origin(Point2::new(3.0, 4.0));
        assert!((l.origin - Point2::new(3.0, 4.0)).length() < TOL);
    }

    #[test]
    fn line_copy() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        let c = l; // Copy (Line2d is Copy)
        assert!((c.direction - Vec2::X).length() < TOL);
        // Verify independence: modifying copy doesn't affect original
        let c = c.with_direction(Vec2::Y);
        assert!(
            (l.direction - Vec2::X).length() < TOL,
            "Original should be unchanged"
        );
    }

    #[test]
    fn line_transform_translation() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        let l = l.translate(Point2::new(5.0, 10.0));
        assert!((l.origin - Point2::new(5.0, 10.0)).length() < TOL);
    }

    #[test]
    fn line_transform_rotation() {
        let l = Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        };
        // Rotate 90 degrees around origin: horizontal becomes vertical
        let l = l.rotate(Point2::ZERO, PI / 2.0);
        assert!((l.direction.x).abs() < TOL, "Direction X should be ~0");
        assert!(
            (l.direction.y - 1.0).abs() < TOL,
            "Direction Y should be ~1"
        );
    }
}

// =============================================================================
// Geom2d_Circle_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_circle_tests {
    use super::*;

    fn make_circle() -> Curve2d {
        Curve2d::Circle(Circle2d::new(Point2::ZERO, 5.0))
    }

    #[test]
    fn circle_closed_periodic() {
        let c = make_circle();
        let p0 = c.point_at(0.0);
        let p2pi = c.point_at(2.0 * PI);
        assert!((p0 - p2pi).length() < TOL);
    }

    #[test]
    fn circle_point_at_quarter() {
        let c = make_circle();
        let p = c.point_at(PI / 2.0);
        assert!((p - Point2::new(0.0, 5.0)).length() < TOL);
    }

    #[test]
    fn circle_point_at_pi() {
        let c = make_circle();
        let p = c.point_at(PI);
        assert!((p - Point2::new(-5.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_origin_offset() {
        let c = Curve2d::Circle(Circle2d::new(Point2::new(10.0, 20.0), 3.0));
        let p = c.point_at(0.0);
        assert!((p - Point2::new(13.0, 20.0)).length() < TOL);
    }
}

// =============================================================================
// Geom2d_Ellipse_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_ellipse_tests {
    use super::*;

    fn make_ellipse() -> Curve2d {
        Curve2d::Ellipse(Ellipse2d {
            center: Point2::ZERO,
            major_dir: Vec2::X,
            major_radius: 5.0,
            minor_radius: 3.0,
        })
    }

    #[test]
    fn ellipse_point_at_major() {
        let e = make_ellipse();
        let p0 = e.point_at(0.0);
        assert!((p0 - Point2::new(5.0, 0.0)).length() < TOL);
    }

    #[test]
    fn ellipse_point_at_minor() {
        let e = make_ellipse();
        let p = e.point_at(PI / 2.0);
        assert!((p - Point2::new(0.0, 3.0)).length() < TOL);
    }
}

// =============================================================================
// Geom2d_Parabola_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_parabola_tests {
    use super::*;

    fn make_parabola() -> Curve2d {
        Curve2d::Parabola(Parabola2d {
            origin: Point2::ZERO,
            axis_dir: Vec2::X,
            focal_param: 2.0,
        })
    }

    #[test]
    fn parabola_point_at_zero() {
        let p = make_parabola();
        let pt = p.point_at(0.0);
        assert!((pt - Point2::ZERO).length() < TOL);
    }
}

// =============================================================================
// Geom2d_BSpline_Test.cxx (simplified, covers basic eval)
// =============================================================================

#[cfg(test)]
mod geom2d_bspline_tests {
    use super::*;

    #[test]
    fn bspline_eval_linear() {
        // BSpline with 2 poles = line segment
        let bs = Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point2::ZERO, Point2::new(10.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let mid = bs.point_at(0.5);
        assert!((mid - Point2::new(5.0, 0.0)).length() < 1e-7);
    }
}

// =============================================================================
// Geom2d_TrimmedCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_trimmed_tests {
    use super::*;

    #[test]
    fn trimmed_truncates_domain() {
        let base = Curve2d::Line(Line2d {
            origin: Point2::ZERO,
            direction: Vec2::X,
        });
        let t = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(base),
            t_min: 2.0,
            t_max: 5.0,
        });
        let pt = t.point_at(2.0);
        assert!((pt - Point2::new(2.0, 0.0)).length() < TOL);
        let pt = t.point_at(5.0);
        assert!((pt - Point2::new(5.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// Geom2d_SineWave_Test.cxx
// =============================================================================

#[cfg(test)]
mod sinewave_tests {
    use super::*;

    fn make_sine() -> Curve2d {
        Curve2d::SineWave(SineWave2d {
            amplitude: 1.0,
            frequency: 1.0,
            phase: 0.0,
        })
    }

    #[test]
    fn sine_at_zero() {
        let s = make_sine();
        let pt = s.point_at(0.0);
        assert!((pt - Point2::new(0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn sine_at_quarter_period() {
        let s = make_sine();
        let pt = s.point_at(PI / 2.0);
        assert!((pt - Point2::new(PI / 2.0, 1.0)).length() < TOL);
    }
}

// =============================================================================
// Geom2d_ArchimedeanSpiral_Test.cxx
// =============================================================================

#[cfg(test)]
mod archimedean_spiral_tests {
    use super::*;

    #[test]
    fn spiral_at_zero() {
        let s = Curve2d::ArchimedeanSpiral(ArchimedeanSpiral2d {
            center: Point2::ZERO,
            a: 0.0,
            b: 1.0,
            start_angle: 0.0,
        });
        let pt = s.point_at(0.0);
        assert!((pt - Point2::ZERO).length() < TOL);
    }
}

// =============================================================================
// Geom2d_CircleInvolute_Test.cxx
// =============================================================================

#[cfg(test)]
mod circle_involute_tests {
    use super::*;

    #[test]
    fn involute_at_zero() {
        let c = Curve2d::CircleInvolute(CircleInvolute2d {
            center: Point2::ZERO,
            base_radius: 1.0,
            start_angle: 0.0,
        });
        let pt = c.point_at(0.0);
        assert!((pt - Point2::new(1.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// Geom2d_LogarithmicSpiral_Test.cxx
// =============================================================================

#[cfg(test)]
mod logspiral_tests {
    use super::*;

    #[test]
    fn logspiral_at_zero() {
        let s = Curve2d::LogarithmicSpiral(LogarithmicSpiral2d {
            center: Point2::ZERO,
            a: 1.0,
            b: 0.5,
            start_angle: 0.0,
        });
        let pt = s.point_at(0.0);
        assert!((pt - Point2::new(1.0, 0.0)).length() < TOL);
    }
}

// =============================================================================
// Geom2dAPI_InterCurveCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod api_intercurve_tests {
    use super::*;
    use rcad_kernel::base::int_ana2d::AnaIntersection2d;

    #[test]
    fn line_line_intersect_origin() {
        let l1 = Line2d { origin: Point2::new(-10.0, 0.0), direction: Vec2::X };
        let l2 = Line2d { origin: Point2::new(0.0, -10.0), direction: Vec2::Y };
        let mut inter = AnaIntersection2d::new();
        inter.perform_lin_lin(&l1, &l2);
        assert!(inter.is_done());
        assert!(!inter.is_empty());
        let d = inter.point(0).value().distance(Point2::new(0.0, 0.0));
        assert!(d < 1e-6, "intersection at origin, got dist {d}");
    }
}
