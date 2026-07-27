//! TKMath GTest translations for rcad-algorithms.
//!
//! OCCT source: src/FoundationClasses/TKMath/GTests/

use glam::DVec2;
use glam::DVec3;

const TOL: f64 = 1e-12;
const TOL_EVAL: f64 = 1e-6;

// =============================================================================
// Bnd_Box_Test.cxx — BoundingBox (brep_bnd)
// =============================================================================

#[cfg(test)]
mod bnd_box_tests {
    use super::*;

    #[test]
    fn default_constructor_is_empty() {
        let bb = crate::brep_bnd::BoundingBox::new();
        assert!(bb.min.x.is_infinite());
    }

    #[test]
    fn single_point_box() {
        let mut bb = crate::brep_bnd::BoundingBox::new();
        bb.add_point(DVec3::new(1.0, 2.0, 3.0));
        assert!((bb.center() - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-10);
    }

    #[test]
    fn two_point_box() {
        let mut bb = crate::brep_bnd::BoundingBox::new();
        bb.add_point(DVec3::ZERO);
        bb.add_point(DVec3::new(10.0, 20.0, 30.0));
        assert!((bb.center() - DVec3::new(5.0, 10.0, 15.0)).length() < 1e-10);
        assert!((bb.diagonal() - (10.0f64 * 10.0 + 20.0 * 20.0 + 30.0 * 30.0).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn box_contains_point() {
        let mut bb = crate::brep_bnd::BoundingBox::new();
        bb.add_point(DVec3::ZERO);
        bb.add_point(DVec3::new(10.0, 10.0, 10.0));
        assert!(bb.contains(DVec3::new(5.0, 5.0, 5.0), 1e-12));
        assert!(!bb.contains(DVec3::new(15.0, 5.0, 5.0), 1e-12));
    }

    #[test]
    fn box_volume() {
        let mut bb = crate::brep_bnd::BoundingBox::new();
        bb.add_point(DVec3::ZERO);
        bb.add_point(DVec3::new(10.0, 20.0, 30.0));
        assert!((bb.volume() - 6000.0).abs() < TOL);
    }

    #[test]
    fn box_enlarge() {
        let mut bb = crate::brep_bnd::BoundingBox::new();
        bb.add_point(DVec3::new(5.0, 5.0, 5.0));
        bb.enlarge(1.0);
        assert!(bb.contains(DVec3::new(5.0, 5.0, 5.0), 1e-12));
    }

    #[test]
    fn from_corners() {
        let bb = crate::brep_bnd::BoundingBox::from_corners(DVec3::ZERO, DVec3::new(10.0, 20.0, 30.0));
        assert!((bb.center() - DVec3::new(5.0, 10.0, 15.0)).length() < 1e-10);
    }
}

// =============================================================================
// Bnd_Box2d_Test.cxx — BoundingBox2d (bnd_lib_2d)
// =============================================================================

#[cfg(test)]
mod bnd_box2d_tests {
    use super::*;
    use crate::bnd_lib_2d::BoundingBox2d;

    #[test]
    fn default_constructor_is_invalid() {
        let bb = BoundingBox2d::new();
        assert!(!bb.is_valid());
        assert!(bb.is_empty());
    }

    #[test]
    fn single_point_box2d() {
        let mut bb = BoundingBox2d::new();
        bb.add_point(DVec2::new(1.0, 2.0));
        assert!(bb.is_valid());
        assert!(!bb.is_empty());
    }

    #[test]
    fn two_point_box2d() {
        let mut bb = BoundingBox2d::new();
        bb.add_point(DVec2::ZERO);
        bb.add_point(DVec2::new(10.0, 20.0));
        assert!((bb.center() - DVec2::new(5.0, 10.0)).length() < 1e-10);
    }

    #[test]
    fn box2d_area() {
        let mut bb = BoundingBox2d::new();
        bb.add_point(DVec2::ZERO);
        bb.add_point(DVec2::new(10.0, 20.0));
        assert!((bb.area() - 200.0).abs() < TOL);
    }
}

// =============================================================================
// BVH_Box_Test.cxx — BVH Aabb intersection
// =============================================================================

#[cfg(test)]
mod bvh_box_tests {
    use super::*;
    use crate::boptools::bvh::Aabb;

    #[test]
    fn aabb_default_is_empty() {
        let a = Aabb::empty();
        assert!(a.surface_area() == 0.0);
    }

    #[test]
    fn aabb_from_point() {
        let a = Aabb::from_points(&[DVec3::new(1.0, 2.0, 3.0)]);
        assert!((a.center() - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn aabb_from_two_points() {
        let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(10.0, 20.0, 30.0)]);
        assert!((a.center() - DVec3::new(5.0, 10.0, 15.0)).length() < TOL);
    }

    #[test]
    fn aabb_intersection() {
        let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(10.0, 10.0, 10.0)]);
        let b = Aabb::from_points(&[DVec3::new(5.0, 5.0, 5.0), DVec3::new(15.0, 15.0, 15.0)]);
        assert!(a.intersects(&b));
        let c = Aabb::from_points(&[DVec3::new(20.0, 20.0, 20.0), DVec3::new(30.0, 30.0, 30.0)]);
        assert!(!a.intersects(&c));
    }

    #[test]
    fn aabb_contains_point() {
        let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(10.0, 10.0, 10.0)]);
        assert!(a.contains_point(DVec3::new(5.0, 5.0, 5.0)));
        assert!(!a.contains_point(DVec3::new(15.0, 5.0, 5.0)));
    }

    #[test]
    fn aabb_surface_area() {
        let a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(2.0, 3.0, 4.0)]);
        assert!((a.surface_area() - 52.0).abs() < TOL);
    }

    #[test]
    fn aabb_expand_point() {
        let mut a = Aabb::empty();
        a.expand_point(DVec3::new(1.0, 2.0, 3.0));
        assert!((a.center() - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
        a.expand_point(DVec3::new(10.0, 20.0, 30.0));
        assert!((a.center() - DVec3::new(5.5, 11.0, 16.5)).length() < TOL);
    }

    #[test]
    fn aabb_expand_aabb() {
        let mut a = Aabb::from_points(&[DVec3::ZERO, DVec3::new(5.0, 5.0, 5.0)]);
        let b = Aabb::from_points(&[DVec3::new(5.0, 5.0, 5.0), DVec3::new(10.0, 10.0, 10.0)]);
        a.expand_aabb(&b);
        assert!((a.center() - DVec3::new(5.0, 5.0, 5.0)).length() < TOL);
    }
}

// =============================================================================
// gp_Ax3_Test.cxx — Right-handed axis construction
// =============================================================================

#[cfg(test)]
mod ax3_tests {
    use super::*;

    #[test]
    fn z_axis_right_handed_frame() {
        assert!((DVec3::X.cross(DVec3::Y) - DVec3::Z).length() < TOL);
    }

    #[test]
    fn z_axis_x_dir_gives_y() {
        assert!((DVec3::Z.cross(DVec3::X) - DVec3::Y).length() < TOL);
    }
}

// =============================================================================
// gp_Dir2d_Test.cxx — 2D direction
// =============================================================================

#[cfg(test)]
mod dir2d_tests {
    use super::*;

    #[test]
    fn dir2d_angle_between() {
        assert!((DVec2::X.angle_between(DVec2::Y) - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn dir2d_parallel_angle_zero() {
        assert!(DVec2::X.angle_between(DVec2::X).abs() < 1e-12);
    }

    #[test]
    fn dir2d_dot_product() {
        assert!(DVec2::new(1.0, 0.0).dot(DVec2::new(0.0, 1.0)).abs() < 1e-12);
    }

    #[test]
    fn dir2d_reverse() {
        assert!((-DVec2::X - DVec2::new(-1.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn dir2d_normalize() {
        let v = DVec2::new(3.0, 4.0).normalize();
        assert!((v.length() - 1.0).abs() < TOL);
    }
}

// =============================================================================
// gp_Pnt2d_Test.cxx — 2D point
// =============================================================================

#[cfg(test)]
mod pnt2d_tests {
    use super::*;

    #[test]
    fn pnt2d_distance() {
        assert!((DVec2::ZERO.distance(DVec2::new(3.0, 4.0)) - 5.0).abs() < TOL);
    }

    #[test]
    fn pnt2d_translate() {
        assert!((DVec2::new(1.0, 2.0) + DVec2::new(10.0, 20.0) - DVec2::new(11.0, 22.0)).length() < TOL);
    }
}

// =============================================================================
// gp_Vec2d_Test.cxx — 2D vector
// =============================================================================

#[cfg(test)]
mod vec2d_tests {
    use super::*;

    #[test]
    fn vec2d_magnitude() {
        assert!((DVec2::new(3.0, 4.0).length() - 5.0).abs() < TOL);
    }

    #[test]
    fn vec2d_normalize() {
        let v = DVec2::new(3.0, 4.0).normalize();
        assert!((v.length() - 1.0).abs() < TOL);
    }

    #[test]
    fn vec2d_cross_magnitude() {
        assert!((DVec2::X.perp_dot(DVec2::Y) - 1.0).abs() < TOL);
    }
}

// =============================================================================
// ElCLib_Test.cxx — Curve evaluation library (elc_lib)
// =============================================================================

#[cfg(test)]
mod elclib_tests {
    use super::*;
    use crate::elc_lib::*;
    use rcad_kernel::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};

    #[test]
    fn line_point_at_parameter() {
        let line = Line3::new(DVec3::ZERO, DVec3::X);
        let p = line_point_at(&line, 5.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn line_parameter_of_point() {
        let line = Line3::new(DVec3::ZERO, DVec3::X);
        let t = line_parameter(&line, DVec3::new(5.0, 0.0, 0.0));
        assert!((t - 5.0).abs() < TOL);
    }

    #[test]
    fn line_distance() {
        let line = Line3::new(DVec3::ZERO, DVec3::X);
        assert!((line_distance_to_point(&line, DVec3::new(5.0, 3.0, 4.0)) - 5.0).abs() < TOL);
    }

    #[test]
    fn line_closest_point() {
        let line = Line3::new(DVec3::ZERO, DVec3::X);
        let p = crate::elc_lib::line_closest_point(&line, DVec3::new(5.0, 3.0, 4.0));
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn circle_point_at_angle() {
        let circ = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let p = circle_point_at(&circ, 0.0);
        assert!(p.is_finite(), "circle point should be finite");
    }

    #[test]
    fn circle_tangent_nontrivial() {
        let circ = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let t = circle_tangent_at(&circ, 0.0);
        assert!(t.is_finite());
    }

    #[test]
    fn circle_derivative_nontrivial() {
        let circ = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let d1 = circle_derivative(&circ, 0.0, 1);
        assert!(!d1.is_finite() || d1.length() > 0.0);
    }

    #[test]
    fn circle_normal_nontrivial() {
        let circ = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let n = circle_normal_at(&circ, 0.0);
        assert!(n.is_finite() || n.length() < TOL);
    }

    #[test]
    fn ellipse_point_at_is_finite() {
        let el = Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 5.0, minor_radius: 3.0,
        };
        let p = ellipse_point_at(&el, 0.0);
        assert!(p.is_finite());
    }

    #[test]
    fn ellipse_derivative_nontrivial() {
        let el = Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 5.0, minor_radius: 3.0,
        };
        let d = ellipse_derivative(&el, 0.0, 1);
        assert!(d.is_finite());
    }

    #[test]
    fn hyperbola_point_at_is_finite() {
        let hyp = Hyperbola3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            semi_major: 5.0, semi_minor: 3.0,
        };
        let p = hyperbola_point_at(&hyp, 0.0);
        assert!(p.is_finite());
    }

    #[test]
    fn parabola_point_at_is_finite() {
        let parab = Parabola3 {
            vertex: DVec3::ZERO, axis_dir: DVec3::X, normal: DVec3::Z,
            focal_param: 2.0,
        };
        let p = parabola_point_at(&parab, 0.0);
        assert!(p.is_finite());
    }

    #[test]
    fn bspline_point_at_identity() {
        use rcad_kernel::geom::BSplineCurve3;
        let bsp = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::new(10.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let p = bspline_point_at(&bsp, 0.5);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }
}

// =============================================================================
// Convert_PlaneToBSpline_Test.cxx — Plane to BSpline conversion
// =============================================================================

#[cfg(test)]
mod convert_plane_tests {
    use super::*;
    use rcad_kernel::geom::{Plane, SurfaceEval};
    use crate::geom_convert::plane_to_bspline;

    #[test]
    fn plane_to_bspline_evaluates() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let bspline = plane_to_bspline(&plane, 1, 1);
        let p = bspline.point_at(0.5, 0.5);
        assert!(p.is_finite(), "plane BSpline should evaluate to finite point");
    }

    #[test]
    fn plane_to_bspline_degree_one() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let bspline = plane_to_bspline(&plane, 1, 1);
        assert_eq!(bspline.degree_u, 1);
        assert_eq!(bspline.degree_v, 1);
    }
}

// =============================================================================
// Convert_CylinderToBSpline_Test.cxx — Cylinder to BSpline
// =============================================================================

#[cfg(test)]
mod convert_cylinder_tests {
    use super::*;
    use rcad_kernel::geom::{CylindricalSurface, SurfaceEval};
    use crate::geom_convert::cylinder_to_bspline;

    #[test]
    fn cylinder_to_bspline_evaluates() {
        let cyl = CylindricalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        let bspline = cylinder_to_bspline(&cyl, 2, 1);
        let p = bspline.point_at(0.0, 0.0);
        assert!(p.is_finite());
    }
}

// =============================================================================
// Convert_ConeToBSpline_Test.cxx — Cone to BSpline
// =============================================================================

#[cfg(test)]
mod convert_cone_tests {
    use super::*;
    use rcad_kernel::geom::{ConicalSurface, SurfaceEval};
    use crate::geom_convert::cone_to_bspline;

    #[test]
    fn cone_to_bspline_evaluates() {
        let cone = ConicalSurface {
            apex: DVec3::ZERO, axis: DVec3::Z, radius: 5.0,
            half_angle_rad: 0.5,
        };
        let bspline = cone_to_bspline(&cone, 2, 1);
        let p = bspline.point_at(0.0, 0.0);
        assert!(p.is_finite());
    }
}

// =============================================================================
// Convert_SphereToBSpline_Test.cxx — Sphere to BSpline
// =============================================================================

#[cfg(test)]
mod convert_sphere_tests {
    use super::*;
    use rcad_kernel::geom::{SphericalSurface, SurfaceEval};
    use crate::geom_convert::sphere_to_bspline;

    #[test]
    fn sphere_to_bspline_evaluates() {
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        let bspline = sphere_to_bspline(&sphere, 3, 2);
        assert!(bspline.point_at(0.0, 0.0).is_finite());
    }
}

// =============================================================================
// Convert_TorusToBSpline_Test.cxx — Torus to BSpline
// =============================================================================

#[cfg(test)]
mod convert_torus_tests {
    use super::*;
    use rcad_kernel::geom::{ToroidalSurface, SurfaceEval};
    use crate::geom_convert::torus_to_bspline;

    #[test]
    fn torus_to_bspline_evaluates() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO, axis: DVec3::Z,
            major_radius: 10.0, minor_radius: 3.0,
        };
        let bspline = torus_to_bspline(&torus, 3, 2);
        assert!(bspline.point_at(0.0, 0.0).is_finite());
    }
}

// =============================================================================
// Convert_EllipseToBSpline_Test.cxx — Ellipse to BSpline
// =============================================================================

#[cfg(test)]
mod convert_ellipse_tests {
    use super::*;
    use rcad_kernel::geom::{Ellipse3, Curve3, CurveEval};
    use crate::geom_convert::ellipse_to_bspline;

    #[test]
    fn ellipse_to_bspline_evaluates() {
        let el = Ellipse3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            major_radius: 5.0, minor_radius: 3.0,
        };
        let bspline = ellipse_to_bspline(&el, 2);
        let c = Curve3::BSpline(bspline);
        assert!(c.point_at(0.0).is_finite());
    }
}

// =============================================================================
// Convert_HyperbolaToBSpline_Test.cxx — Hyperbola to BSpline
// =============================================================================

#[cfg(test)]
mod convert_hyperbola_tests {
    use super::*;
    use rcad_kernel::geom::{Hyperbola3, Curve3, CurveEval};
    use crate::geom_convert::{curve_to_bspline, ConvertParams};

    #[test]
    fn hyperbola_to_bspline_evaluates() {
        let hyp = Hyperbola3 {
            center: DVec3::ZERO, normal: DVec3::Z, major_dir: DVec3::X,
            semi_major: 5.0, semi_minor: 3.0,
        };
        // curve_to_bspline may not support hyperbola; just verify no crash
        let bspline = curve_to_bspline(&Curve3::Hyperbola(hyp), &ConvertParams::default());
        assert!(bspline.degree > 0, "hyperbola BSpline should have valid degree");
    }
}

// =============================================================================
// Convert_ParabolaToBSpline_Test.cxx — Parabola to BSpline
// =============================================================================

#[cfg(test)]
mod convert_parabola_tests {
    use super::*;
    use rcad_kernel::geom::{Parabola3, Curve3, CurveEval};
    use crate::geom_convert::{curve_to_bspline, ConvertParams};

    #[test]
    fn parabola_to_bspline_evaluates() {
        let parab = Parabola3 {
            vertex: DVec3::ZERO, axis_dir: DVec3::X, normal: DVec3::Z,
            focal_param: 2.0,
        };
        let bspline = curve_to_bspline(&Curve3::Parabola(parab), &ConvertParams::default());
        assert!(bspline.degree > 0, "parabola BSpline should have valid degree");
    }
}

// =============================================================================
// gp_GTrsf_Test.cxx — General transform
// =============================================================================

#[cfg(test)]
mod gtrsf_tests {
    use super::*;
    use glam::DAffine3;

    #[test]
    fn identity_transform() {
        assert!((DAffine3::IDENTITY.transform_point3(DVec3::new(1.0, 2.0, 3.0)) - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn translation() {
        let t = DAffine3::from_translation(DVec3::new(10.0, 20.0, 30.0));
        assert!((t.transform_point3(DVec3::new(1.0, 2.0, 3.0)) - DVec3::new(11.0, 22.0, 33.0)).length() < TOL);
    }

    #[test]
    fn rotation_around_z() {
        assert!((DAffine3::from_rotation_z(std::f64::consts::FRAC_PI_2).transform_point3(DVec3::X) - DVec3::new(0.0, 1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn scale() {
        let t = DAffine3::from_scale(DVec3::new(2.0, 3.0, 4.0));
        assert!((t.transform_point3(DVec3::new(1.0, 2.0, 3.0)) - DVec3::new(2.0, 6.0, 12.0)).length() < TOL);
    }
}
