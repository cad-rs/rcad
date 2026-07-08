//! OCCT ExtremaPC GTest translations — point-to-curve distance extrema.
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/
//!   ExtremaPC_Line_Test.cxx
//!   ExtremaPC_Circle_Test.cxx
//!   ExtremaPC_Ellipse_Test.cxx
//!   ExtremaPC_Parabola_Test.cxx
//!   ExtremaPC_Hyperbola_Test.cxx
//!   ExtremaPC_BezierCurve_Test.cxx
//!   ExtremaPC_BSplineCurve_Test.cxx
//!   ExtremaPC_OffsetCurve_Test.cxx
//!   ExtremaPC_Curve_Test.cxx (aggregator)
//!   ExtremaPC_Comparison_Test.cxx
//!   ExtremaPC_SearchMode_Test.cxx
//!   ExtremaPC_ExtendedGeometry_Test.cxx
//!   Extrema_ExtPC_Test.cxx
//!
//! Rust equivalent: rcad_algorithms::extrema::{distance_point_curve, closest_point_on_curve}

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Line3, Ellipse3, Parabola3, Hyperbola3,
                        BezierCurve3, BSplineCurve3, CurveEval, OffsetCurve3};

const TOL: f64 = 1e-6;

/// Minimal GTest-style helper: assert near with tolerance
fn assert_near(actual: f64, expected: f64, tol: f64, msg: &str) {
    assert!((actual - expected).abs() < tol, "{}: expected {}, got {}", msg, expected, actual);
}

// =============================================================================
// ExtremaPC_Line_Test.cxx — point-to-line distance
// =============================================================================

#[cfg(test)]
mod extremapc_line_tests {
    use super::*;
    use rcad_algorithms::extrema::{distance_point_curve, closest_point_on_curve};

    fn make_line() -> Curve3 {
        Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X })
    }

    #[test]
    fn line_point_on_line_at_origin() {
        let line = make_line();
        let (dist, _param) = distance_point_curve(DVec3::ZERO, &line);
        assert_near(dist, 0.0, TOL, "point on origin");
    }

    #[test]
    fn line_point_on_line_positive() {
        let line = make_line();
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &line);
        assert_near(dist, 0.0, TOL, "point on line");
    }

    #[test]
    fn line_point_off_line_y_offset() {
        let line = make_line();
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 3.0, 0.0), &line);
        assert_near(dist, 3.0, TOL, "y offset distance");
    }

    #[test]
    fn line_point_off_line_z_offset() {
        let line = make_line();
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 0.0, 4.0), &line);
        assert_near(dist, 4.0, TOL, "z offset distance");
    }

    #[test]
    fn line_point_off_line_yz_offset() {
        let line = make_line();
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 3.0, 4.0), &line);
        assert_near(dist, 5.0, TOL, "yz offset distance");
    }
}

// =============================================================================
// ExtremaPC_Circle_Test.cxx — point-to-circle distance
// =============================================================================

#[cfg(test)]
mod extremapc_circle_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    fn make_circle() -> Curve3 {
        Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 10.0))
    }

    #[test]
    fn circle_point_outside_on_x_axis() {
        let circle = make_circle();
        let (dist, _param) = distance_point_curve(DVec3::new(20.0, 0.0, 0.0), &circle);
        // Distance from (20,0,0) to circle center = 20, minus radius 10 = 10
        assert_near(dist, 10.0, TOL, "point outside on X axis");
    }

    #[test]
    fn circle_point_outside_on_y_axis() {
        let circle = make_circle();
        let (dist, _param) = distance_point_curve(DVec3::new(0.0, 25.0, 0.0), &circle);
        assert_near(dist, 15.0, TOL, "point outside on Y axis");
    }

    #[test]
    fn circle_point_inside_on_x_axis() {
        let circle = make_circle();
        let (dist, _param) = distance_point_curve(DVec3::new(3.0, 0.0, 0.0), &circle);
        assert_near(dist, 7.0, TOL, "point inside on X axis");
    }

    #[test]
    fn circle_point_inside_on_y_axis() {
        let circle = make_circle();
        let (dist, _param) = distance_point_curve(DVec3::new(0.0, 5.0, 0.0), &circle);
        assert_near(dist, 5.0, TOL, "point inside on Y axis");
    }

    #[test]
    fn circle_point_at_center() {
        let circle = make_circle();
        let (dist, _param) = distance_point_curve(DVec3::ZERO, &circle);
        assert_near(dist, 10.0, TOL, "point at center distance = radius");
    }
}

// =============================================================================
// ExtremaPC_Ellipse_Test.cxx — point-to-ellipse distance
// =============================================================================

#[cfg(test)]
mod extremapc_ellipse_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    fn make_ellipse() -> Curve3 {
        Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 20.0, minor_radius: 10.0,
        })
    }

    #[test]
    fn ellipse_point_on_major_axis_outside() {
        let e = make_ellipse();
        let (dist, _param) = distance_point_curve(DVec3::new(30.0, 0.0, 0.0), &e);
        assert_near(dist, 10.0, TOL, "point on major axis outside");
    }

    #[test]
    fn ellipse_point_on_major_axis_negative() {
        let e = make_ellipse();
        let (dist, _param) = distance_point_curve(DVec3::new(-30.0, 0.0, 0.0), &e);
        assert_near(dist, 10.0, TOL, "point on major axis negative");
    }

    #[test]
    fn ellipse_point_on_minor_axis_outside() {
        let e = make_ellipse();
        let (dist, _param) = distance_point_curve(DVec3::new(0.0, 20.0, 0.0), &e);
        assert_near(dist, 10.0, TOL, "point on minor axis outside");
    }

    #[test]
    fn ellipse_point_on_minor_axis_inside() {
        let e = make_ellipse();
        let (dist, _param) = distance_point_curve(DVec3::new(0.0, 5.0, 0.0), &e);
        assert_near(dist, 5.0, TOL, "point on minor axis inside");
    }
}

// =============================================================================
// ExtremaPC_Parabola_Test.cxx — point-to-parabola distance
// =============================================================================

#[cfg(test)]
mod extremapc_parabola_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    #[test]
    fn parabola_point_at_vertex() {
        // y² = 20x parabola (focal param = 10)
        let p = Curve3::Parabola(Parabola3 {
            center: DVec3::ZERO, normal: DVec3::Z, x_dir: DVec3::X,
            focal_param: 10.0,
        });
        let (dist, _param) = distance_point_curve(DVec3::ZERO, &p);
        assert_near(dist, 0.0, 0.1, "point at vertex");
    }

    #[test]
    fn parabola_point_on_axis_positive() {
        let p = Curve3::Parabola(Parabola3 {
            center: DVec3::ZERO, normal: DVec3::Z, x_dir: DVec3::X,
            focal_param: 10.0,
        });
        let (dist, _param) = distance_point_curve(DVec3::new(10.0, 0.0, 0.0), &p);
        // Distance to vertex should be close to 10 (depends on parabola shape)
        assert!(dist > 0.0);
    }

    #[test]
    fn parabola_point_above_plane() {
        let p = Curve3::Parabola(Parabola3 {
            center: DVec3::ZERO, normal: DVec3::Z, x_dir: DVec3::X,
            focal_param: 10.0,
        });
        let (dist, _param) = distance_point_curve(DVec3::new(0.0, 0.0, -5.0), &p);
        // Z component dominates
        assert!(dist > 4.0);
    }
}

// =============================================================================
// ExtremaPC_Hyperbola_Test.cxx — point-to-hyperbola distance
// =============================================================================

#[cfg(test)]
mod extremapc_hyperbola_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    #[test]
    fn hyperbola_point_at_vertex() {
        // x²/16 - y²/9 = 1
        let h = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            semi_major: 4.0, semi_minor: 3.0,
        });
        let (dist, _param) = distance_point_curve(DVec3::new(4.0, 0.0, 0.0), &h);
        assert_near(dist, 0.0, TOL, "point at vertex");
    }

    #[test]
    fn hyperbola_point_on_axis() {
        let h = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            semi_major: 4.0, semi_minor: 3.0,
        });
        // Point between two branches — distance should be positive
        let (dist, _param) = distance_point_curve(DVec3::new(0.0, 0.0, 0.0), &h);
        assert!(dist > 3.0, "distance from origin to hyperbola branches");
    }
}

// =============================================================================
// ExtremaPC_BezierCurve_Test.cxx — point-to-Bezier distance
// =============================================================================

#[cfg(test)]
mod extremapc_bezier_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    fn make_cubic_bezier() -> Curve3 {
        Curve3::Bezier(BezierCurve3 {
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(3.0, 2.0, 0.0),
                DVec3::new(4.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        })
    }

    fn make_linear_bezier() -> Curve3 {
        Curve3::Bezier(BezierCurve3 {
            control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(10.0, 0.0, 0.0)],
            weights: vec![1.0; 2],
        })
    }

    #[test]
    fn bezier_point_on_curve_start() {
        let b = make_cubic_bezier();
        let (dist, _param) = distance_point_curve(DVec3::ZERO, &b);
        assert_near(dist, 0.0, TOL, "point at start of Bezier");
    }

    #[test]
    fn bezier_point_on_curve_end() {
        let b = make_cubic_bezier();
        let (dist, _param) = distance_point_curve(DVec3::new(4.0, 0.0, 0.0), &b);
        assert_near(dist, 0.0, TOL, "point at end of Bezier");
    }

    #[test]
    fn bezier_point_on_curve_middle() {
        let b = make_cubic_bezier();
        let pt_mid = b.point_at(0.5);
        let (dist, _param) = distance_point_curve(pt_mid, &b);
        assert_near(dist, 0.0, 1e-3, "point at middle of Bezier");
    }

    #[test]
    fn bezier_linear_point_projection() {
        let b = make_linear_bezier();
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 3.0, 0.0), &b);
        assert_near(dist, 3.0, TOL, "point projected onto linear Bezier");
    }

    #[test]
    fn bezier_linear_point_before_start() {
        let b = make_linear_bezier();
        let (dist, _param) = distance_point_curve(DVec3::new(-5.0, 0.0, 0.0), &b);
        assert_near(dist, 5.0, TOL, "point before start of linear Bezier");
    }
}

// =============================================================================
// ExtremaPC_BSplineCurve_Test.cxx — point-to-BSpline distance
// =============================================================================

#[cfg(test)]
mod extremapc_bspline_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    fn make_cubic_bspline() -> Curve3 {
        Curve3::BSpline(BSplineCurve3 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::new(0.0, 0.0, 0.0),
                DVec3::new(1.0, 2.0, 0.0),
                DVec3::new(2.0, 2.0, 0.0),
                DVec3::new(3.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 4],
        })
    }

    fn make_linear_bspline() -> Curve3 {
        Curve3::BSpline(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(10.0, 0.0, 0.0)],
            weights: vec![1.0; 2],
        })
    }

    #[test]
    fn bspline_point_on_curve_start() {
        let b = make_cubic_bspline();
        let (dist, _param) = distance_point_curve(DVec3::ZERO, &b);
        assert_near(dist, 0.0, TOL, "point at start");
    }

    #[test]
    fn bspline_point_on_curve_end() {
        let b = make_cubic_bspline();
        let (dist, _param) = distance_point_curve(DVec3::new(3.0, 0.0, 0.0), &b);
        assert_near(dist, 0.0, TOL, "point at end");
    }

    #[test]
    fn bspline_point_near_curve_above() {
        let b = make_cubic_bspline();
        let (dist, _param) = distance_point_curve(DVec3::new(1.5, 3.0, 0.0), &b);
        assert!(dist > 0.0);
        assert!(dist < 2.0, "distance above curve");
    }

    #[test]
    fn bspline_linear_point_projection() {
        let b = make_linear_bspline();
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 3.0, 0.0), &b);
        assert_near(dist, 3.0, TOL, "point projected onto linear BSpline");
    }
}

// =============================================================================
// ExtremaPC_OffsetCurve_Test.cxx — point-to-offset-curve distance
// =============================================================================

#[cfg(test)]
mod extremapc_offset_curve_tests {
    use super::*;
    use rcad_algorithms::extrema::distance_point_curve;

    #[test]
    fn offset_curve_basic_projection() {
        // Offset a line by 2 units
        let base = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let offset = Curve3::Offset(OffsetCurve3 {
            basis: Box::new(base),
            offset_distance: 2.0,
        });
        let (dist, _param) = distance_point_curve(DVec3::new(5.0, 2.0, 0.0), &offset);
        assert_near(dist, 0.0, TOL, "point on offset curve");
    }
}

// =============================================================================
// IntAna_IntQuadQuad_Test.cxx — quadric-quadric intersection
// =============================================================================

#[cfg(test)]
mod intana_intquadquad_tests {
    use rcad_algorithms::int_ana::intersect_plane_plane_intana;

    #[test]
    fn plane_plane_intersection_exists() {
        // Two non-parallel planes intersect in a line
        let p1 = rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::Z };
        let p2 = rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::X };
        let result = intersect_plane_plane_intana(p1, p2);
        assert!(result.is_some(), "two non-parallel planes should intersect");
    }
}

// =============================================================================
// GeomConvert_Test.cxx — geometry conversion tests
// =============================================================================

#[cfg(test)]
mod geomconvert_tests {
    use rcad_kernel::geom::{Line3, Curve3, Circle3, CurveEval};
    use rcad_algorithms::geom_convert::{line_to_bspline, circle_to_bspline};

    #[test]
    fn line_conversion_exact() {
        let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
        let bs = line_to_bspline(&line, 1);
        assert_eq!(bs.degree, 1);
        assert_eq!(bs.control_points.len(), 2);
        assert!((bs.point_at(0.0) - DVec3::ZERO).length() < 1e-10);
        assert!((bs.point_at(1.0) - DVec3::X).length() < 1e-10);
    }

    #[test]
    fn circle_conversion_radius_constant() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let bs = circle_to_bspline(&circle, 2);
        for i in 0..8 {
            let t = i as f64 / 8.0;
            let p = bs.point_at(t);
            let r = (p - DVec3::ZERO).length();
            assert!((r - 5.0).abs() < 1e-4, "radius should be constant");
        }
    }
}
