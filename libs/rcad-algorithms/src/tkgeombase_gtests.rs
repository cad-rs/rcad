//! TKGeomBase GTest translations.
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/
//!
//! Tests for: ExtremaPC, GeomConvert, IntAna,
//! Hermit, CompCurveToBSplineCurve, ProjLib_Cone.

use crate::extrema::*;
use crate::geom_convert::*;
use crate::int_ana::*;
use crate::tkgeombase_algo::*;
use glam::DVec3;
use rcad_kernel::geom::*;

const TOL: f64 = 1e-7;

// =============================================================================
// ExtremaPC_Line_Test.cxx — point-line distance
// =============================================================================

#[cfg(test)]
mod extremapc_line_tests {
    use super::*;

    #[test]
    fn point_on_line_distance_zero() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &line);
        assert!(d < TOL, "point on line should have distance 0, got {d}");
    }

    #[test]
    fn point_off_line_by_3() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 3.0, 0.0), &line);
        assert!((d - 3.0).abs() < TOL, "distance should be 3, got {d}");
    }

    #[test]
    fn point_line_returns_projection_param() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let (_d, t) = distance_point_curve(DVec3::new(7.0, 0.0, 0.0), &line);
        assert!(
            (t - 7.0).abs() < TOL,
            "projection param should be 7, got {t}"
        );
    }
}

// =============================================================================
// ExtremaPC_Circle_Test.cxx — point-circle distance
// =============================================================================

#[cfg(test)]
mod extremapc_circle_tests {
    use super::*;

    #[test]
    fn point_at_circle_center() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let (d, _t) = distance_point_curve(DVec3::ZERO, &circle);
        assert!(
            (d - 5.0).abs() < TOL,
            "center distance should be radius 5, got {d}"
        );
    }

    #[test]
    fn point_on_circle_edge() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &circle);
        assert!(
            d < TOL,
            "point on circle edge should have distance 0, got {d}"
        );
    }
}

// =============================================================================
// ExtremaPC_Ellipse_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_ellipse_tests {
    use super::*;

    #[test]
    fn point_at_ellipse_major() {
        let ell = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 5.0,
            minor_radius: 3.0,
        });
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &ell);
        assert!(
            d < TOL,
            "point on ellipse edge should have distance 0, got {d}"
        );
    }
}

// =============================================================================
// ExtremaPC_BSplineCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_bspline_tests {
    use super::*;

    #[test]
    fn point_at_bspline_start() {
        let bs = Curve3::BSpline(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::new(10.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        });
        let (d, _t) = distance_point_curve(DVec3::ZERO, &bs);
        assert!(d < TOL, "point at start should have distance 0, got {d}");
    }
}

// =============================================================================
// GeomConvert_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom_convert_tests {
    use super::*;

    #[test]
    fn line_to_bspline_roundtrip() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let bs = line_to_bspline(&line, 1);
        let p0 = bs.point_at(0.0);
        let p5 = bs.point_at(5.0);
        assert!((p0 - DVec3::ZERO).length() < TOL);
        assert!((p5 - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_to_bspline_works() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let bs = circle_to_bspline(&circle, 3);
        // Just check the curve has valid control points
        assert!(!bs.control_points.is_empty(), "should have control points");
        assert!(bs.degree > 0, "should have positive degree");
    }
}

// =============================================================================
// IntAna_IntQuadQuad_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_ana_plane_plane_tests {
    use super::*;

    #[test]
    fn perpendicular_planes_intersect() {
        let p1 = Plane::new(DVec3::ZERO, DVec3::Z);
        let p2 = Plane::new(DVec3::ZERO, DVec3::X);
        let result = intersect_plane_plane_intana(&p1, &p2);
        match result {
            PlnPlnResult::Line(_) => {}
            _ => panic!(
                "perpendicular planes should intersect in a line, got {:?}",
                result
            ),
        }
    }

    #[test]
    fn parallel_offset_planes_no_intersection() {
        let p1 = Plane::new(DVec3::ZERO, DVec3::Z);
        let p2 = Plane::new(DVec3::new(0.0, 0.0, 10.0), DVec3::Z);
        let result = intersect_plane_plane_intana(&p1, &p2);
        assert!(matches!(result, PlnPlnResult::Parallel));
    }

    #[test]
    fn identical_planes_coincident() {
        let p1 = Plane::new(DVec3::ZERO, DVec3::Z);
        let p2 = Plane::new(DVec3::ZERO, DVec3::Z);
        let result = intersect_plane_plane_intana(&p1, &p2);
        assert!(matches!(result, PlnPlnResult::Coincident));
    }
}

// =============================================================================
// Hermit_Test.cxx — BSpline Hermite evaluation
// =============================================================================

#[cfg(test)]
mod hermit_tests {
    use super::*;

    fn make_bspline() -> BSplineCurve3 {
        BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::new(10.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn hermit_returns_2d_curve() {
        let bs = make_bspline();
        let result = hermit_solution(&bs);
        let p = result.point_at(0.0);
        assert!(p.y.is_finite());
    }

    #[test]
    fn hermit_solutionbis_returns_finite() {
        let bs = make_bspline();
        let (a, b) = hermit_solutionbis(&bs);
        assert!(a.is_finite());
        assert!(b.is_finite());
    }
}

// =============================================================================
// GeomConvert_CompCurveToBSplineCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod concat_bspline_tests {
    use super::*;

    fn make_bspline_seg(start: DVec3, end: DVec3) -> BSplineCurve3 {
        BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![start, end],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn concat_two_linear_bsplines() {
        let c1 = make_bspline_seg(DVec3::ZERO, DVec3::new(5.0, 0.0, 0.0));
        let c2 = make_bspline_seg(DVec3::new(5.0, 0.0, 0.0), DVec3::new(10.0, 0.0, 0.0));
        let result = concat_bsplines(&c1, &c2, 1e-7);
        assert!(result.is_some(), "should concatenate continuous curves");
    }
}

// =============================================================================
// ProjLib_Cone_Test.cxx — project circle onto cone
// =============================================================================

#[cfg(test)]
mod projlib_cone_tests {
    use super::*;

    #[test]
    fn project_circle_onto_cone_yields_2d_curve() {
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0);
        let cone = ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4,
        };
        let result = project_circle_onto_cone(&circle, &cone);
        assert!(result.is_some(), "circle on cone should project");
    }
}

// =============================================================================
// ExtremaPC_SearchMode_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_searchmode_tests {
    use super::*;

    fn make_line() -> Curve3 {
        Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        })
    }

    fn make_circle() -> Curve3 {
        Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 10.0))
    }

    #[test]
    fn line_min_mode() {
        let line = make_line();
        let pt = DVec3::new(25.0, 10.0, 0.0);
        let r = find_extrema_curve(&line, pt, 1e-6, SearchMode::Min);
        assert!(r.is_done());
        assert_eq!(r.nb_ext(), 1);
        assert!((r.min_square_distance().sqrt() - 10.0).abs() < TOL);
    }

    #[test]
    fn line_max_mode() {
        let pt = DVec3::new(0.0, 10.0, 0.0);
        let bounded = Curve3::Line(Line3 {
            origin: DVec3::new(-50.0, 0.0, 0.0),
            direction: DVec3::X,
        });
        let r = find_extrema_curve(&bounded, pt, 1e-6, SearchMode::Max);
        assert!(r.is_done());
        assert!(r.nb_ext() >= 1);
    }

    #[test]
    fn circle_min_mode_point_outside() {
        let circle = make_circle();
        let pt = DVec3::new(20.0, 0.0, 0.0);
        let r = find_extrema_curve(&circle, pt, 1e-6, SearchMode::Min);
        assert!(r.is_done());
        assert_eq!(r.nb_ext(), 1);
        assert!((r.min_square_distance().sqrt() - 10.0).abs() < TOL);
    }

    #[test]
    fn circle_max_mode_point_outside() {
        let circle = make_circle();
        let pt = DVec3::new(20.0, 0.0, 0.0);
        let r = find_extrema_curve(&circle, pt, 1e-6, SearchMode::Max);
        assert!(r.is_done());
        assert_eq!(r.nb_ext(), 1);
        assert!((r.max_square_distance().sqrt() - 30.0).abs() < TOL);
    }

    #[test]
    fn circle_min_mode_point_inside() {
        let circle = make_circle();
        let pt = DVec3::new(3.0, 0.0, 0.0);
        let r = find_extrema_curve(&circle, pt, 1e-6, SearchMode::Min);
        assert!(r.is_done());
        assert_eq!(r.nb_ext(), 1);
        assert!((r.min_square_distance().sqrt() - 7.0).abs() < TOL);
    }

    #[test]
    fn circle_minmax_mode_different_min_max() {
        let circle = make_circle();
        let pt = DVec3::new(15.0, 0.0, 0.0);
        let r = find_extrema_curve(&circle, pt, 1e-6, SearchMode::MinMax);
        assert!(r.is_done());
        assert!(r.nb_ext() >= 2);
        assert!((r.min_square_distance().sqrt() - 5.0).abs() < TOL);
        assert!((r.max_square_distance().sqrt() - 25.0).abs() < 2.0);
    }
}

// =============================================================================
// ExtremaPC_ExtendedGeometry_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_extended_geo_tests {
    use super::*;

    #[test]
    fn circle_translated() {
        let c = Curve3::Circle(Circle3::new(DVec3::new(100.0, 200.0, 50.0), DVec3::Z, 25.0));
        let pt = DVec3::new(150.0, 200.0, 50.0);
        let (dist, _) = distance_point_curve(pt, &c);
        assert!((dist - 25.0).abs() < TOL);
    }

    #[test]
    fn circle_rotated_xy_plane() {
        let normal = DVec3::new(1.0, 1.0, 0.0).normalize();
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, normal, 15.0));
        let pt = DVec3::new(20.0, 20.0, 5.0);
        let (dist, _) = distance_point_curve(pt, &c);
        assert!(dist > 0.0);
    }

    #[test]
    fn circle_very_small() {
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 0.001));
        let pt = DVec3::new(0.002, 0.0, 0.0);
        let (dist, _) = distance_point_curve(pt, &c);
        assert!((dist - 0.001).abs() < 1e-9);
    }

    #[test]
    fn circle_very_large() {
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1000.0));
        let pt = DVec3::new(1500.0, 0.0, 0.0);
        let (dist, _) = distance_point_curve(pt, &c);
        assert!((dist - 500.0).abs() < TOL);
    }

    #[test]
    fn ellipse_translated() {
        let e = Curve3::Ellipse(Ellipse3 {
            center: DVec3::new(50.0, 100.0, 25.0),
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 30.0,
            minor_radius: 15.0,
        });
        let pt = DVec3::new(90.0, 100.0, 25.0);
        let (dist, _) = distance_point_curve(pt, &e);
        assert!((dist - 10.0).abs() < TOL);
    }

    #[test]
    fn ellipse_high_eccentricity() {
        let e = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 50.0,
            minor_radius: 2.0,
        });
        let pt = DVec3::new(40.0, 5.0, 0.0);
        let (dist, _) = distance_point_curve(pt, &e);
        assert!(dist > 0.0);
    }
}

// =============================================================================
// ExtremaPC_Comparison_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_comparison_tests {
    use super::*;

    #[test]
    fn line_point_on_line() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let (d1, _) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &line);
        let (d2, _) = distance_point_curve(DVec3::new(5.0, 3.0, 4.0), &line);
        assert!(d1 < TOL);
        assert!((d2 - 5.0).abs() < TOL);
    }

    #[test]
    fn circle_point_outside() {
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 10.0));
        let (d1, _) = distance_point_curve(DVec3::new(20.0, 0.0, 0.0), &c);
        let (d2, _) = distance_point_curve(DVec3::new(3.0, 0.0, 0.0), &c);
        assert!((d1 - 10.0).abs() < TOL);
        assert!((d2 - 7.0).abs() < TOL);
    }

    #[test]
    fn circle_point_off_plane() {
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 10.0));
        let (dist, _) = distance_point_curve(DVec3::new(15.0, 0.0, 5.0), &c);
        assert!(dist > 0.0);
    }

    #[test]
    fn ellipse_point_on_major_axis() {
        let e = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 20.0,
            minor_radius: 10.0,
        });
        let (dist, _) = distance_point_curve(DVec3::new(30.0, 0.0, 0.0), &e);
        assert!((dist - 10.0).abs() < TOL);
    }

    #[test]
    fn parabola_point_near_vertex() {
        let p = Curve3::Parabola(Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 4.0,
        });
        let (dist, _) = distance_point_curve(DVec3::new(1.0, 2.0, 0.0), &p);
        assert!(dist > 0.0);
    }
}

// =============================================================================
// ExtremaPC_Parabola_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_parabola_tests {
    use super::*;

    #[test]
    fn parabola_point_at_vertex() {
        let p = Curve3::Parabola(Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 10.0,
        });
        let (dist, _) = distance_point_curve(DVec3::ZERO, &p);
        assert!(dist < 0.1, "point at vertex, dist={}", dist);
    }

    #[test]
    fn parabola_point_above_plane() {
        let p = Curve3::Parabola(Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 10.0,
        });
        let (dist, _) = distance_point_curve(DVec3::new(0.0, 0.0, -5.0), &p);
        assert!(dist > 4.0);
    }
}

// =============================================================================
// ExtremaPC_Hyperbola_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_hyperbola_tests {
    use super::*;

    #[test]
    fn hyperbola_point_at_vertex() {
        let h = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 4.0,
            semi_minor: 3.0,
        });
        let (dist, _) = distance_point_curve(DVec3::new(4.0, 0.0, 0.0), &h);
        assert!(dist < TOL, "point at vertex, dist={}", dist);
    }

    #[test]
    fn hyperbola_point_on_axis() {
        let h = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 4.0,
            semi_minor: 3.0,
        });
        let (dist, _) = distance_point_curve(DVec3::ZERO, &h);
        assert!(dist > 3.0, "distance from origin to hyperbola branches");
    }
}

// =============================================================================
// ExtremaPC_BezierCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_bezier_tests {
    use super::*;

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
    fn bezier_point_on_start() {
        let b = make_cubic_bezier();
        let (dist, _) = distance_point_curve(DVec3::ZERO, &b);
        assert!(dist < TOL, "point at start, dist={}", dist);
    }

    #[test]
    fn bezier_point_on_end() {
        let b = make_cubic_bezier();
        let (dist, _) = distance_point_curve(DVec3::new(4.0, 0.0, 0.0), &b);
        assert!(dist < TOL, "point at end, dist={}", dist);
    }

    #[test]
    fn bezier_linear_projection() {
        let b = make_linear_bezier();
        let (dist, _) = distance_point_curve(DVec3::new(5.0, 3.0, 0.0), &b);
        assert!(
            (dist - 3.0).abs() < TOL,
            "point projected onto linear Bezier, dist={}",
            dist
        );
    }

    #[test]
    fn bezier_linear_before_start() {
        let b = make_linear_bezier();
        let (dist, _) = distance_point_curve(DVec3::new(-5.0, 0.0, 0.0), &b);
        assert!(
            (dist - 5.0).abs() < TOL,
            "point before start, dist={}",
            dist
        );
    }
}

// =============================================================================
// Extrema_ExtPC_Test.cxx
// =============================================================================

#[cfg(test)]
mod extrema_extpc_tests {
    use super::*;

    #[test]
    fn bug24945_cylinder_parameter_normalization() {
        let pt = DVec3::new(-1725.97, 843.257, -4.22741e-13);
        let c = Curve3::Circle(Circle3::new(
            DVec3::new(0.0, 843.257, 0.0),
            -DVec3::Y,
            1725.9708621929999,
        ));
        let (dist, _param) = distance_point_curve(pt, &c);
        assert!(
            dist < 1.0,
            "bug24945: projected distance should be small, got {}",
            dist
        );
    }
}
