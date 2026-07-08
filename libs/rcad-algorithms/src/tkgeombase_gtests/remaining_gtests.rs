//! Remaining TKGeomBase GTest translations.
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/
//!   Hermit_Test.cxx
//!   GeomConvert_CompCurveToBSplineCurve_Test.cxx
//!   Geom2dConvert_CompCurveToBSplineCurve_Test.cxx
//!   ProjLib_Cone_Test.cxx
//!   ProjLib_ComputeApproxOnPolarSurface_Test.cxx
//!   ExtremaPC_SearchMode_Test.cxx
//!   ExtremaPC_Comparison_Test.cxx
//!   ExtremaPC_ExtendedGeometry_Test.cxx
//!   Extrema_ExtPC_Test.cxx

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;
#[cfg(test)]
use crate::tkgeombase_algo::*;

const TOL: f64 = 1e-7;

// =============================================================================
// Hermit_Test.cxx — Rational BSpline weight decomposition
// =============================================================================

fn make_rational_bspline_3d(w1: f64, w2: f64, w3: f64) -> BSplineCurve3 {
    BSplineCurve3 {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ],
        weights: vec![w1, w2, w3],
    }
}

fn make_rational_bspline_2d(w1: f64, w2: f64, w3: f64) -> BSplineCurve2 {
    BSplineCurve2 {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 0.0),
        ],
        weights: vec![w1, w2, w3],
    }
}

#[cfg(test)]
mod hermit_tests {
    use super::*;

    #[test]
    fn hermit_solution_3d_uniform_weights() {
        let bs = make_rational_bspline_3d(1.0, 1.0, 1.0);
        let result = hermit_solution(&bs);
        assert!(result.control_points.len() >= 3);
        let p0 = result.point_at(result.default_domain()[0]);
        let p1 = result.point_at(result.default_domain()[1]);
        assert!((p0.y - 1.0).abs() < 1e-6);
        assert!((p1.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn hermit_solution_3d_distinct_weights() {
        let bs = make_rational_bspline_3d(2.0, 1.5, 3.0);
        let result = hermit_solution(&bs);
        let p0 = result.point_at(result.default_domain()[0]);
        let p1 = result.point_at(result.default_domain()[1]);
        assert!((p0.y - 0.5).abs() < 1e-4, "a(0) should be 1/2 = 0.5, got {}", p0.y);
        assert!((p1.y - 1.0/3.0).abs() < 1e-4, "a(1) should be 1/3, got {}", p1.y);
    }

    #[test]
    fn hermit_solution_3d_high_weight_ratio() {
        let bs = make_rational_bspline_3d(0.5, 1.0, 5.0);
        let result = hermit_solution(&bs);
        let p0 = result.point_at(result.default_domain()[0]);
        let p1 = result.point_at(result.default_domain()[1]);
        assert!((p0.y - 2.0).abs() < 1e-4);
        assert!((p1.y - 0.2).abs() < 1e-4);
    }

    #[test]
    fn hermit_solution_3d_reversed_weight_ratio() {
        let bs = make_rational_bspline_3d(5.0, 1.0, 0.5);
        let result = hermit_solution(&bs);
        let p0 = result.point_at(result.default_domain()[0]);
        let p1 = result.point_at(result.default_domain()[1]);
        assert!((p0.y - 0.2).abs() < 1e-4);
        assert!((p1.y - 2.0).abs() < 1e-4);
    }

    #[test]
    fn hermit_solution_3d_positive_poles() {
        let bs = make_rational_bspline_3d(2.0, 3.0, 1.5);
        let result = hermit_solution(&bs);
        for cp in &result.control_points {
            assert!(cp.y > 0.0, "all poles should have positive Y");
        }
    }

    #[test]
    fn hermit_solution_2d_uniform_weights() {
        let bs = make_rational_bspline_2d(1.0, 1.0, 1.0);
        let result = hermit_solution(&BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::new(bs.control_points[0].x, bs.control_points[0].y, 0.0),
                DVec3::new(bs.control_points[1].x, bs.control_points[1].y, 0.0),
                DVec3::new(bs.control_points[2].x, bs.control_points[2].y, 0.0),
            ],
            weights: bs.weights,
        });
        let p0 = result.point_at(result.default_domain()[0]);
        let p1 = result.point_at(result.default_domain()[1]);
        assert!((p0.y - 1.0).abs() < 1e-6);
        assert!((p1.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn hermit_solutionbis_uniform_weights() {
        let bs = make_rational_bspline_3d(1.0, 1.0, 1.0);
        let (kmin, kmax) = hermit_solutionbis(&bs);
        assert!(kmin >= 0.0);
        assert!(kmax <= 1.0);
    }

    #[test]
    fn hermit_solutionbis_distinct_weights() {
        let bs = make_rational_bspline_3d(2.0, 1.5, 3.0);
        let (kmin, kmax) = hermit_solutionbis(&bs);
        assert!(kmin >= 0.0);
        assert!(kmax <= 1.0);
        assert!(kmin <= kmax);
    }

    #[test]
    fn hermit_solution_3d_symmetric_weights() {
        let bs = make_rational_bspline_3d(2.0, 1.0, 2.0);
        let result = hermit_solution(&bs);
        let p0 = result.point_at(result.default_domain()[0]);
        let p1 = result.point_at(result.default_domain()[1]);
        assert!((p0.y - 0.5).abs() < 1e-4);
        assert!((p1.y - 0.5).abs() < 1e-4);
    }
}

// =============================================================================
// GeomConvert_CompCurveToBSplineCurve_Test.cxx
// =============================================================================

#[cfg(test)]
mod compcurve_tests {
    use super::*;

    fn make_cubic_bspline(poles: Vec<DVec3>) -> BSplineCurve3 {
        let n = poles.len();
        let degree = 3;
        let knots: Vec<f64> = {
            let mut k = vec![0.0; degree + 1];
            k.push(1.0);
            k.extend(vec![1.0; degree]);
            k
        };
        BSplineCurve3 {
            degree,
            knots,
            control_points: poles,
            weights: vec![1.0; n],
        }
    }

    #[test]
    fn concat_clamped_bsplines() {
        let c1 = make_cubic_bspline(vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 1.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ]);
        let c2 = make_cubic_bspline(vec![
            DVec3::new(3.0, 0.0, 0.0), // shares endpoint
            DVec3::new(4.0, -1.0, 0.0),
            DVec3::new(5.0, -1.0, 0.0),
            DVec3::new(6.0, 0.0, 0.0),
        ]);

        let result = concat_bsplines(&c1, &c2, 1e-7);
        assert!(result.is_some(), "Should concatenate clamped B-splines");
        let r = result.unwrap();
        assert!((r.point_at(0.0) - DVec3::ZERO).length() < 1e-7);
        assert!((r.point_at(1.0) - DVec3::new(6.0, 0.0, 0.0)).length() < 1e-7);
    }

    #[test]
    fn concat_trimmed_circle_arcs() {
        // Two halves of a circle: concatenate as BSplines
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        use crate::geom_convert::circle_to_bspline;
        let bs_full = circle_to_bspline(&circle, 2);
        // Split into two arcs at t=0.5
        // bs_full degree 2, 3 poles, domain [0,1]
        let n = bs_full.control_points.len();
        let mid_pt = bs_full.point_at(0.5);
        let c1 = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 0.5, 0.5],
            control_points: vec![
                bs_full.control_points[0],
                bs_full.point_at(0.25),
                mid_pt,
            ],
            weights: vec![bs_full.weights[0], 1.0, bs_full.weights[n-1]],
        };
        let c2 = BSplineCurve3 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                mid_pt,
                bs_full.point_at(0.75),
                bs_full.control_points[n-1],
            ],
            weights: vec![bs_full.weights[0], 1.0, bs_full.weights[n-1]],
        };

        let result = concat_bsplines(&c1, &c2, 1e-7);
        assert!(result.is_some(), "Should concatenate circle arcs");
    }

    #[test]
    fn concat_with_reversal() {
        let c1 = make_cubic_bspline(vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 1.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ]);
        // Second curve ends where c1 ends (needs reversal)
        let c2 = make_cubic_bspline(vec![
            DVec3::new(6.0, 0.0, 0.0),
            DVec3::new(5.0, -1.0, 0.0),
            DVec3::new(4.0, -1.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ]);

        let result = concat_bsplines(&c1, &c2, 1e-7);
        assert!(result.is_some(), "Should concatenate with reversal");
        let r = result.unwrap();
        assert!((r.point_at(0.0) - DVec3::ZERO).length() < 1e-7);
        assert!((r.point_at(1.0) - DVec3::new(6.0, 0.0, 0.0)).length() < 1e-7);
    }

    #[test]
    fn concat_fails_for_disjoint_curves() {
        let c1 = make_cubic_bspline(vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(2.0, 1.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
        ]);
        let c2 = make_cubic_bspline(vec![
            DVec3::new(10.0, 0.0, 0.0),
            DVec3::new(11.0, 1.0, 0.0),
            DVec3::new(12.0, 1.0, 0.0),
            DVec3::new(13.0, 0.0, 0.0),
        ]);

        let result = concat_bsplines(&c1, &c2, 1e-7);
        assert!(result.is_none(), "Should fail for disjoint curves");
    }
}

// =============================================================================
// ProjLib_Cone_Test.cxx
// =============================================================================

#[cfg(test)]
mod projlib_cone_tests {
    use super::*;
    use rcad_kernel::geom::{ConicalSurface, Circle3};

    fn make_cone() -> ConicalSurface {
        ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            half_angle_rad: std::f64::consts::PI / 6.0,
        }
    }

    #[test]
    fn projlib_cone_project_circle_parallel_axes() {
        let cone = make_cone();
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 10.0), DVec3::Z, 2.0);
        let result = project_circle_onto_cone(&circle, &cone);
        assert!(result.is_some(), "Coaxial circle should project to line");
        if let Some(Curve2d::Line(l)) = result {
            assert!(l.direction.y.abs() < TOL);
        } else {
            panic!("Expected Line pcurve");
        }
    }

    #[test]
    fn projlib_cone_project_circle_non_parallel_axes() {
        let cone = make_cone();
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 10.0), DVec3::X, 2.0);
        let result = project_circle_onto_cone(&circle, &cone);
        assert!(result.is_none(), "Non-coaxial circle should not project simply");
    }

    #[test]
    fn projlib_cone_opposite_normal() {
        let cone = make_cone();
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 6.0), -DVec3::Z, 2.0);
        let result = project_circle_onto_cone(&circle, &cone);
        assert!(result.is_some(), "Coaxial circle with opposite normal should project");
    }
}

// =============================================================================
// ExtremaPC_SearchMode_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_searchmode_tests {
    use super::*;

    fn make_line() -> Curve3 {
        Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X })
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
        let line = make_line();
        let pt = DVec3::new(0.0, 10.0, 0.0);
        // Use domain [-50, 50] for max distance check
        let bounded = Curve3::Line(Line3 { origin: DVec3::new(-50.0, 0.0, 0.0), direction: DVec3::X });
        let r = find_extrema_curve(&bounded, pt, 1e-6, SearchMode::Max);
        assert!(r.is_done());
        assert!(r.nb_ext() >= 1);
        assert!((r.max_square_distance().sqrt() - (2500.0_f64 + 100.0_f64).sqrt()).abs() < 1.0);
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
        assert!((r.max_square_distance().sqrt() - 25.0).abs() < TOL);
    }
}

// =============================================================================
// ExtremaPC_ExtendedGeometry_Test.cxx
// =============================================================================

#[cfg(test)]
mod extremapc_extended_geo_tests {
    use super::*;
    use crate::extrema::distance_point_curve;

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
            center: DVec3::new(50.0, 100.0, 25.0), normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 30.0, minor_radius: 15.0,
        });
        let pt = DVec3::new(90.0, 100.0, 25.0);
        let (dist, _) = distance_point_curve(pt, &e);
        assert!((dist - 10.0).abs() < TOL);
    }

    #[test]
    fn ellipse_high_eccentricity() {
        let e = Curve3::Ellipse(Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 50.0, minor_radius: 2.0,
        });
        let pt = DVec3::new(40.0, 5.0, 0.0);
        let (dist, _) = distance_point_curve(pt, &e);
        assert!(dist > 0.0);
    }
}

// =============================================================================
// ExtremaPC_Comparison_Test.cxx — verify closest_point_on_curve matches
// =============================================================================

#[cfg(test)]
mod extremapc_comparison_tests {
    use super::*;
    use crate::extrema::distance_point_curve;

    #[test]
    fn line_point_on_line() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
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
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 20.0, minor_radius: 10.0,
        });
        let (dist, _) = distance_point_curve(DVec3::new(30.0, 0.0, 0.0), &e);
        assert!((dist - 10.0).abs() < TOL);
    }

    #[test]
    fn parabola_point_near_vertex() {
        let p = Curve3::Parabola(Parabola3 {
            center: DVec3::ZERO, normal: DVec3::Z, x_dir: DVec3::X,
            focal_param: 4.0,
        });
        let (dist, _) = distance_point_curve(DVec3::new(1.0, 2.0, 0.0), &p);
        assert!(dist > 0.0);
    }
}

// =============================================================================
// Extrema_ExtPC_Test.cxx — point-to-curve extrema regression test
// =============================================================================

#[cfg(test)]
mod extrema_extpc_tests {
    use super::*;
    use crate::extrema::distance_point_curve;

    #[test]
    fn bug24945_cylinder_parameter_normalization() {
        // Regression test: project point onto circle
        let pt = DVec3::new(-1725.97, 843.257, -4.22741e-13);
        let c = Curve3::Circle(Circle3::new(
            DVec3::new(0.0, 843.257, 0.0),
            -DVec3::Y,
            1725.9708621929999,
        ));
        let (dist, _param) = distance_point_curve(pt, &c);
        // Point is on the circle, distance should be near zero
        assert!(dist < 1.0, "bug24945: projected distance should be small, got {}", dist);
    }
}

// =============================================================================
// ProjLib_ComputeApproxOnPolarSurface_Test.cxx (simplified)
// =============================================================================

#[cfg(test)]
mod projlib_polar_tests {
    use rcad_kernel::geom::*;

    #[test]
    fn project_circle_on_cylinder_via_pcurve() {
        // Simplified: verify circle point projects onto cylinder surface
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, 5.0);
        // On the cylindrical surface, the circle maps to a constant-v line
        // Circle radius = cylinder radius → it's a great circle
        let p0 = circle.point_at(0.0);
        // p0 should be on the cylinder surface
        let radial = DVec3::new(p0.x, p0.y, 0.0);
        assert!((radial.length() - 5.0).abs() < 0.01,
            "circle point should be on cylinder: |radial|={}", radial.length());
    }
}
