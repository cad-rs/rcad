//! OCCT-aligned TKGeomBase GTest translations.
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/
//!
//! Tests for: ExtremaPC, GeomConvert, IntAna,
//! Hermit, CompCurveToBSplineCurve, ProjLib_Cone.

use glam::DVec3;
use rcad_kernel::geom::*;
use crate::extrema::*;
use crate::int_ana::*;
use crate::geom_convert::*;
use crate::tkgeombase_algo::*;

const TOL: f64 = 1e-7;

// =============================================================================
// ExtremaPC_Line_Test.cxx — point-line distance
// =============================================================================

#[cfg(test)]
mod extremapc_line_tests {
    use super::*;

    #[test]
    fn point_on_line_distance_zero() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &line);
        assert!(d < TOL, "point on line should have distance 0, got {d}");
    }

    #[test]
    fn point_off_line_by_3() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 3.0, 0.0), &line);
        assert!((d - 3.0).abs() < TOL, "distance should be 3, got {d}");
    }

    #[test]
    fn point_line_returns_projection_param() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let (_d, t) = distance_point_curve(DVec3::new(7.0, 0.0, 0.0), &line);
        assert!((t - 7.0).abs() < TOL, "projection param should be 7, got {t}");
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
        assert!((d - 5.0).abs() < TOL, "center distance should be radius 5, got {d}");
    }

    #[test]
    fn point_on_circle_edge() {
        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let (d, _t) = distance_point_curve(DVec3::new(5.0, 0.0, 0.0), &circle);
        assert!(d < TOL, "point on circle edge should have distance 0, got {d}");
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
        assert!(d < TOL, "point on ellipse edge should have distance 0, got {d}");
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
        let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
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
        let p1 = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
        let p2 = Plane { origin: DVec3::ZERO, normal: DVec3::X };
        let result = intersect_plane_plane_intana(&p1, &p2);
        match result {
            PlnPlnResult::Line(_) => {},
            _ => panic!("perpendicular planes should intersect in a line, got {:?}", result),
        }
    }

    #[test]
    fn parallel_offset_planes_no_intersection() {
        let p1 = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
        let p2 = Plane { origin: DVec3::new(0.0, 0.0, 10.0), normal: DVec3::Z };
        let result = intersect_plane_plane_intana(&p1, &p2);
        assert!(matches!(result, PlnPlnResult::Parallel));
    }

    #[test]
    fn identical_planes_coincident() {
        let p1 = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
        let p2 = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
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
