//! OCCT-aligned TKMath GTest translations.
//!
//! OCCT source: src/FoundationClasses/TKMath/GTests/
//!
//! Files translated:
//!   gp_Lin_Test.cxx          — Line3
//!   gp_Circ_Test.cxx         — Circle3
//!   gp_Pln_Test.cxx          — Plane
//!   gp_Elips_Test.cxx        — Ellipse3
//!   gp_Hypr_Test.cxx         — Hyperbola3
//!   gp_Parab_Test.cxx        — Parabola3
//!   gp_Torus_Test.cxx        — ToroidalSurface
//!   gp_Cylinder_Test.cxx     — CylindricalSurface
//!   gp_Sphere_Test.cxx       — SphericalSurface
//!   gp_Cone_Test.cxx         — ConicalSurface
//!   gp_Dir_Test.cxx          — Direction (unit vector) operations
//!   gp_Mat_Test.cxx          — Matrix operations
//!   gp_Trsf_Test.cxx         — Transform
//!   gp_Quaternion_Test.cxx   — Quaternion
//!   gp_XYZ_Test.cxx          — DVec3 arithmetic
//!   Convert_CircleToBSplineCurve_Test.cxx — Circle to BSpline conversion
//!   Convert_SphereToBSplineSurface_Test.cxx — Sphere to BSpline conversion

use glam::DVec3;

const TOL: f64 = 1e-12;
const TOL_ANG: f64 = 1e-12;
const TOL_EVAL: f64 = 1e-6;

// =============================================================================
// gp_XYZ_Test.cxx — DVec3 (XYZ) arithmetic
// =============================================================================

#[cfg(test)]
mod xyz_tests {
    use super::*;

    #[test]
    fn default_constructor_is_zero() {
        let v = DVec3::ZERO;
        assert!((v.x - 0.0).abs() < TOL);
        assert!((v.y - 0.0).abs() < TOL);
        assert!((v.z - 0.0).abs() < TOL);
    }

    #[test]
    fn coordinate_constructor() {
        let v = DVec3::new(1.0, 2.0, 3.0);
        assert!((v.x - 1.0).abs() < TOL);
        assert!((v.y - 2.0).abs() < TOL);
        assert!((v.z - 3.0).abs() < TOL);
    }

    #[test]
    fn distance_and_square_distance() {
        let a = DVec3::ZERO;
        let b = DVec3::new(3.0, 4.0, 0.0);
        assert!((a.distance(b) - 5.0).abs() < TOL);
        assert!((a.distance_squared(b) - 25.0).abs() < TOL);
    }

    #[test]
    fn cross_product() {
        assert!((DVec3::X.cross(DVec3::Y) - DVec3::Z).length() < TOL);
    }

    #[test]
    fn dot_product() {
        assert!((DVec3::new(1.0, 2.0, 3.0).dot(DVec3::new(4.0, 5.0, 6.0)) - 32.0).abs() < TOL);
    }

    #[test]
    fn add_and_subtract() {
        let a = DVec3::new(1.0, 2.0, 3.0);
        let b = DVec3::new(4.0, 5.0, 6.0);
        assert!((a + b - DVec3::new(5.0, 7.0, 9.0)).length() < TOL);
        assert!((a - b - DVec3::new(-3.0, -3.0, -3.0)).length() < TOL);
    }

    #[test]
    fn scale_and_negate() {
        let a = DVec3::new(1.0, -2.0, 3.0);
        assert!((a * 2.0 - DVec3::new(2.0, -4.0, 6.0)).length() < TOL);
        assert!((-a - DVec3::new(-1.0, 2.0, -3.0)).length() < TOL);
    }

    #[test]
    fn normalize() {
        let n = DVec3::new(3.0, 4.0, 0.0).normalize();
        assert!((n.length() - 1.0).abs() < TOL);
        assert!((n - DVec3::new(0.6, 0.8, 0.0)).length() < TOL);
    }
}

// =============================================================================
// gp_Dir_Test.cxx — Direction (unit vector) operations
// =============================================================================

#[cfg(test)]
mod dir_tests {
    use super::*;

    #[test]
    fn direction_angle_between() {
        let angle = DVec3::X.angle_between(DVec3::Y);
        assert!((angle - std::f64::consts::FRAC_PI_2).abs() < TOL_ANG);
    }

    #[test]
    fn direction_parallel_angle_zero() {
        assert!(DVec3::X.angle_between(DVec3::X).abs() < TOL_ANG);
    }

    #[test]
    fn direction_opposite_angle_pi() {
        assert!((DVec3::X.angle_between(-DVec3::X) - std::f64::consts::PI).abs() < TOL_ANG);
    }

    #[test]
    fn direction_cross_product_orthogonal() {
        assert!((DVec3::X.cross(DVec3::Y) - DVec3::Z).length() < TOL);
    }

    #[test]
    fn direction_dot_perpendicular_is_zero() {
        assert!(DVec3::X.dot(DVec3::Y).abs() < TOL);
    }

    #[test]
    fn direction_is_normalized() {
        let v = DVec3::new(1.0, 2.0, 3.0).normalize();
        assert!((v.length() - 1.0).abs() < TOL);
    }

    #[test]
    fn direction_reverse() {
        assert!((-DVec3::X - DVec3::NEG_X).length() < TOL);
    }
}

// =============================================================================
// gp_Lin_Test.cxx — Line3
// =============================================================================

#[cfg(test)]
mod lin_tests {
    use super::*;
    use crate::geom::Line3;

    #[test]
    fn constructor_from_point_and_dir() {
        let line = Line3::new(DVec3::ZERO, DVec3::Z);
        assert!((line.origin - DVec3::ZERO).length() < TOL);
        assert!((line.direction - DVec3::Z).length() < TOL);
    }

    #[test]
    fn constructor_normalizes_direction() {
        let line = Line3::new(DVec3::ZERO, DVec3::new(0.0, 3.0, 4.0));
        assert!((line.direction.length() - 1.0).abs() < TOL);
        assert!((line.direction - DVec3::new(0.0, 0.6, 0.8)).length() < TOL);
    }

    #[test]
    fn distance_to_point() {
        let line = Line3::new(DVec3::ZERO, DVec3::X);
        assert!((line.distance(DVec3::new(5.0, 3.0, 4.0)) - 5.0).abs() < TOL);
    }

    #[test]
    fn distance_to_point_on_line_is_zero() {
        assert!(
            Line3::new(DVec3::ZERO, DVec3::X)
                .distance(DVec3::new(10.0, 0.0, 0.0))
                .abs()
                < TOL
        );
    }

    #[test]
    fn reversed_parameter_negates() {
        assert!((Line3::new(DVec3::ZERO, DVec3::X).reversed_parameter(5.0) - (-5.0)).abs() < TOL);
    }

    #[test]
    fn lin_is_not_closed() {
        assert!(!Line3::new(DVec3::ZERO, DVec3::X).is_closed());
    }

    #[test]
    fn lin_is_not_periodic() {
        assert!(!Line3::new(DVec3::ZERO, DVec3::X).is_periodic());
    }
}

// =============================================================================
// gp_Circ_Test.cxx — Circle3 (via CurveEval)
// =============================================================================

#[cfg(test)]
mod circ_tests {
    use super::*;
    use crate::CurveEval;
    use crate::geom::Circle3;

    fn make_circ() -> Circle3 {
        Circle3::new(DVec3::ZERO, DVec3::Z, 5.0)
    }

    #[test]
    fn constructor_with_radius() {
        let c = make_circ();
        assert!((c.radius - 5.0).abs() < TOL);
        assert!((c.center - DVec3::ZERO).length() < TOL);
    }

    #[test]
    fn constructor_normalizes_normal() {
        let c = Circle3::new(DVec3::ZERO, DVec3::new(0.0, 3.0, 4.0), 2.0);
        assert!((c.normal.length() - 1.0).abs() < TOL);
        assert!((c.normal - DVec3::new(0.0, 0.6, 0.8)).length() < TOL);
    }

    #[test]
    fn point_at_zero_angle() {
        // OCCT-aligned gp_Ax2: normal=Z gives x_dir=X, so P(0) = center + r*X = (5,0,0)
        let p = make_circ().point_at(0.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_quarter_turn() {
        // y_dir = Z×X = Y, P(π/2) = center + r*Y = (0,5,0)
        let p = make_circ().point_at(std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(0.0, 5.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_half_turn() {
        let p = make_circ().point_at(std::f64::consts::PI);
        assert!((p - DVec3::new(-5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_full_turn_returns_to_start() {
        let c = make_circ();
        assert!((c.point_at(std::f64::consts::TAU) - c.point_at(0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn center_offset() {
        let c = Circle3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::Z, 5.0);
        // OCCT-aligned: x_dir=X, P(0) = center + 5*X = (6, 2, 3)
        let p = c.point_at(0.0);
        assert!((p - DVec3::new(6.0, 2.0, 3.0)).length() < TOL_EVAL);
    }
}

// =============================================================================
// gp_Pln_Test.cxx — Plane
// =============================================================================

#[cfg(test)]
mod pln_tests {
    use super::*;
    use crate::geom::Plane;

    #[test]
    fn constructor_from_origin_and_normal() {
        let p = Plane::new(DVec3::ZERO, DVec3::Z);
        assert!((p.origin - DVec3::ZERO).length() < TOL);
        assert!((p.normal - DVec3::Z).length() < TOL);
    }

    #[test]
    fn constructor_normalizes_normal() {
        let p = Plane::new(DVec3::ZERO, DVec3::new(0.0, 3.0, 4.0));
        assert!((p.normal.length() - 1.0).abs() < TOL);
    }
}

// =============================================================================
// gp_Elips_Test.cxx — Ellipse3
// =============================================================================

#[cfg(test)]
mod elips_tests {
    use super::*;
    use crate::CurveEval;
    use crate::geom::Ellipse3;

    fn make_elips() -> Ellipse3 {
        Ellipse3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            major_radius: 5.0,
            minor_radius: 3.0,
        }
    }

    #[test]
    fn constructor_with_major_minor() {
        let e = make_elips();
        assert!((e.major_radius - 5.0).abs() < TOL);
        assert!((e.minor_radius - 3.0).abs() < TOL);
    }

    #[test]
    fn point_at_zero_angle_on_major_axis() {
        // P(0) = center + major_radius*major_dir = (5,0,0)
        let p = make_elips().point_at(0.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_quarter_on_minor_axis() {
        // P(π/2) = center + minor_radius*minor_dir
        // minor_dir = normal × major_dir = Z×X = Y
        let p = make_elips().point_at(std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(0.0, 3.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn foci_distance() {
        let e = make_elips();
        let c = (e.major_radius * e.major_radius - e.minor_radius * e.minor_radius).sqrt();
        assert!((c - 4.0).abs() < TOL);
    }
}

// =============================================================================
// gp_Hypr_Test.cxx — Hyperbola3
// =============================================================================

#[cfg(test)]
mod hypr_tests {
    use super::*;
    use crate::CurveEval;
    use crate::geom::Hyperbola3;

    #[test]
    fn constructor_with_semi_axes() {
        let h = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 5.0,
            semi_minor: 3.0,
        };
        assert!((h.semi_major - 5.0).abs() < TOL);
        assert!((h.semi_minor - 3.0).abs() < TOL);
    }

    #[test]
    fn point_at_zero_time() {
        let h = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 5.0,
            semi_minor: 3.0,
        };
        // P(0) = center + a*cosh(0)*major_dir + b*sinh(0)*minor_dir = (5,0,0)
        assert!((h.point_at(0.0) - DVec3::new(5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn eccentricity_formula() {
        let h = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 5.0,
            semi_minor: 3.0,
        };
        let e = (1.0 + (h.semi_minor * h.semi_minor) / (h.semi_major * h.semi_major)).sqrt();
        assert!((e - (34.0_f64).sqrt() / 5.0).abs() < TOL);
    }
}

// =============================================================================
// gp_Parab_Test.cxx — Parabola3
// =============================================================================

#[cfg(test)]
mod parab_tests {
    use super::*;
    use crate::CurveEval;
    use crate::geom::Parabola3;

    #[test]
    fn constructor_with_focal_length() {
        let p = Parabola3 {
            vertex: DVec3::ZERO,
            axis_dir: DVec3::X,
            normal: DVec3::Z,
            focal_param: 2.0,
        };
        assert!((p.focal_param - 2.0).abs() < TOL);
    }

    #[test]
    fn point_at_vertex_is_origin() {
        let p = Parabola3 {
            vertex: DVec3::ZERO,
            axis_dir: DVec3::X,
            normal: DVec3::Z,
            focal_param: 2.0,
        };
        assert!((p.point_at(0.0) - DVec3::ZERO).length() < TOL_EVAL);
    }
}

// =============================================================================
// gp_Torus_Test.cxx — ToroidalSurface (SurfaceEval)
// =============================================================================

#[cfg(test)]
mod torus_tests {
    use super::*;
    use crate::SurfaceEval;
    use crate::geom::ToroidalSurface;

    fn make_torus() -> ToroidalSurface {
        ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 10.0,
            minor_radius: 3.0,
        }
    }

    #[test]
    fn constructor_with_major_minor_radius() {
        let t = make_torus();
        assert!((t.major_radius - 10.0).abs() < TOL);
        assert!((t.minor_radius - 3.0).abs() < TOL);
    }

    #[test]
    fn point_at_uv_zero() {
        // OCCT-aligned: any_perpendicular(Z) = X as x_ax
        // P(0,0) = center + (R+r)*X = (13, 0, 0)
        let p = make_torus().point_at(0.0, 0.0);
        assert!((p - DVec3::new(13.0, 0.0, 0.0)).length() < TOL_EVAL);
    }
}

// =============================================================================
// gp_Cylinder_Test.cxx — CylindricalSurface (SurfaceEval)
// =============================================================================

#[cfg(test)]
mod cylinder_tests {
    use super::*;
    use crate::SurfaceEval;
    use crate::geom::CylindricalSurface;

    #[test]
    fn constructor_with_radius() {
        let c = CylindricalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        assert!((c.radius - 5.0).abs() < TOL);
    }

    #[test]
    fn point_at_u0_v0() {
        // OCCT-aligned gp_Ax2: axis=Z gives ref_dir=X, so P(0,0) = r*X = (5,0,0)
        let c = CylindricalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        let p = c.point_at(0.0, 0.0);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_u90_v0() {
        // P(π/2, 0) = r * axis×ref_dir = 5 * Z×X = 5*Y = (0,5,0)
        let c = CylindricalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        let p = c.point_at(std::f64::consts::FRAC_PI_2, 0.0);
        assert!((p - DVec3::new(0.0, 5.0, 0.0)).length() < TOL_EVAL);
    }
}

// =============================================================================
// gp_Sphere_Test.cxx — SphericalSurface (SurfaceEval)
// =============================================================================

#[cfg(test)]
mod sphere_tests {
    use super::*;
    use crate::SurfaceEval;
    use crate::geom::SphericalSurface;

    fn make_sphere() -> SphericalSurface {
        // axis=Z: OCCT-aligned ref_dir = X
        SphericalSurface::new(DVec3::ZERO, DVec3::Z, 5.0)
    }

    #[test]
    fn constructor_with_radius() {
        let s = make_sphere();
        assert!((s.radius - 5.0).abs() < TOL);
    }

    #[test]
    fn point_at_north_pole() {
        // v=0 (north pole): P = center + r*axis = (0,0,5)
        let p = make_sphere().point_at(0.0, 0.0);
        assert!((p - DVec3::new(0.0, 0.0, 5.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_equator() {
        // v=π/2 (equator), u=0: P = r*ref_dir = r*X = (5,0,0)
        let p = make_sphere().point_at(0.0, std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOL_EVAL);
    }

    #[test]
    fn point_at_90_degrees_along_equator() {
        // v=π/2 (equator), u=π/2: P = r*axis×ref_dir = 5*Z×X = (0,5,0)
        let p = make_sphere().point_at(std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        assert!((p - DVec3::new(0.0, 5.0, 0.0)).length() < TOL_EVAL);
    }
}

// =============================================================================
// gp_Cone_Test.cxx — ConicalSurface (SurfaceEval)
// =============================================================================

#[cfg(test)]
mod cone_tests {
    use super::*;
    use crate::SurfaceEval;
    use crate::geom::ConicalSurface;

    #[test]
    fn constructor_with_half_angle() {
        let c = ConicalSurface::new(DVec3::ZERO, DVec3::Z, 0.0, std::f64::consts::FRAC_PI_4);
        assert!((c.half_angle_rad - std::f64::consts::FRAC_PI_4).abs() < TOL);
    }

    #[test]
    fn point_at_v_zero_is_finite() {
        let c = ConicalSurface::new(DVec3::ZERO, DVec3::Z, 5.0, 0.5);
        assert!(c.point_at(0.0, 0.0).is_finite());
    }
}

// =============================================================================
// gp_Mat_Test.cxx — Matrix 3x3 operations (DMat3)
// =============================================================================

#[cfg(test)]
mod mat_tests {
    use super::*;
    use glam::DMat3;

    #[test]
    fn identity_multiply() {
        assert!(
            (DMat3::IDENTITY * DVec3::new(1.0, 2.0, 3.0) - DVec3::new(1.0, 2.0, 3.0)).length()
                < TOL
        );
    }

    #[test]
    fn scale_matrix() {
        let m = DMat3::from_diagonal(DVec3::new(2.0, 3.0, 4.0));
        assert!((m * DVec3::new(1.0, 2.0, 3.0) - DVec3::new(2.0, 6.0, 12.0)).length() < TOL);
    }

    #[test]
    fn matrix_multiply() {
        let s6 = DMat3::from_diagonal(DVec3::splat(2.0)) * DMat3::from_diagonal(DVec3::splat(3.0));
        assert!((s6 * DVec3::new(1.0, 1.0, 1.0) - DVec3::splat(6.0)).length() < TOL);
    }

    #[test]
    fn transpose() {
        let m = DMat3::from_cols(
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(4.0, 5.0, 6.0),
            DVec3::new(7.0, 8.0, 9.0),
        );
        assert!(
            (m.transpose() * DVec3::new(1.0, 0.0, 0.0) - DVec3::new(1.0, 4.0, 7.0)).length() < TOL
        );
    }

    #[test]
    fn determinant_identity() {
        assert!((DMat3::IDENTITY.determinant() - 1.0).abs() < TOL);
    }

    #[test]
    fn determinant_scale() {
        assert!((DMat3::from_diagonal(DVec3::new(2.0, 3.0, 4.0)).determinant() - 24.0).abs() < TOL);
    }
}

// =============================================================================
// gp_Trsf_Test.cxx — Transform (DAffine3)
// =============================================================================

#[cfg(test)]
mod trsf_tests {
    use super::*;
    use glam::DAffine3;

    #[test]
    fn identity_transform() {
        assert!(
            (DAffine3::IDENTITY.transform_point3(DVec3::new(1.0, 2.0, 3.0))
                - DVec3::new(1.0, 2.0, 3.0))
            .length()
                < TOL
        );
    }

    #[test]
    fn translation() {
        let t = DAffine3::from_translation(DVec3::new(10.0, 20.0, 30.0));
        assert!(
            (t.transform_point3(DVec3::new(1.0, 2.0, 3.0)) - DVec3::new(11.0, 22.0, 33.0)).length()
                < TOL
        );
    }

    #[test]
    fn rotation_around_z() {
        assert!(
            (DAffine3::from_rotation_z(std::f64::consts::FRAC_PI_2).transform_point3(DVec3::X)
                - DVec3::new(0.0, 1.0, 0.0))
            .length()
                < TOL
        );
    }

    #[test]
    fn scale() {
        let t = DAffine3::from_scale(DVec3::new(2.0, 3.0, 4.0));
        assert!(
            (t.transform_point3(DVec3::new(1.0, 2.0, 3.0)) - DVec3::new(2.0, 6.0, 12.0)).length()
                < TOL
        );
    }
}

// =============================================================================
// gp_Quaternion_Test.cxx — Quaternion (DQuat)
// =============================================================================

#[cfg(test)]
mod quat_tests {
    use super::*;
    use glam::DQuat;

    #[test]
    fn identity_quaternion_no_rotation() {
        assert!((DQuat::IDENTITY * DVec3::X - DVec3::X).length() < TOL);
    }

    #[test]
    fn rotation_90_around_z() {
        let q = DQuat::from_rotation_z(std::f64::consts::FRAC_PI_2);
        assert!((q * DVec3::X - DVec3::new(0.0, 1.0, 0.0)).length() < TOL);
    }

    #[test]
    fn rotation_180_around_z() {
        assert!(
            (DQuat::from_rotation_z(std::f64::consts::PI) * DVec3::X - DVec3::new(-1.0, 0.0, 0.0))
                .length()
                < TOL
        );
    }

    #[test]
    fn quaternion_slerp() {
        let mid = DQuat::IDENTITY.slerp(DQuat::from_rotation_z(std::f64::consts::FRAC_PI_2), 0.5);
        let r = mid * DVec3::X;
        assert!((r - DVec3::new(1.0, 1.0, 0.0).normalize()).length() < 1e-6);
    }
}

// =============================================================================
// Convert_CircleToBSplineCurve_Test.cxx — Circle to BSpline conversion
// =============================================================================

#[cfg(test)]
mod convert_circle_tests {
    use super::*;
    use crate::geom::{Circle3, Curve3, CurveEval};

    #[test]
    fn circle_to_bspline_degree_two() {
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let bspline = crate::base::convert::circle_to_bspline(&c);
        assert_eq!(bspline.degree, 2);
    }

    #[test]
    fn circle_to_bspline_eval_matches_circle() {
        let c = Circle3::new(DVec3::ZERO, DVec3::Z, 5.0);
        let bspline = crate::base::convert::circle_to_bspline(&c);
        let curve = Curve3::BSpline(bspline);
        // BSpline domain is [0, 1], circle domain is [0, TAU].
        // Map the circle parameter t ∈ [0, TAU] to BSpline parameter s ∈ [0, 1].
        let eval_at_circle_t = |t: f64| {
            let s = t / std::f64::consts::TAU; // normalized to [0, 1]
            curve.point_at(s)
        };
        // Only test a few sample points to check rough alignment
        for i in 0..4 {
            let t = i as f64 * std::f64::consts::FRAC_PI_2; // 0, π/2, π, 3π/2
            let p_circ = c.point_at(t);
            let p_bsp = eval_at_circle_t(t);
            assert!(
                (p_circ - p_bsp).length() < 1e-4,
                "mismatch at t={t}: {p_circ:?} vs {p_bsp:?}"
            );
        }
    }
}

// =============================================================================
// Convert_SphereToBSplineSurface_Test.cxx — Sphere to BSpline conversion
// =============================================================================

#[cfg(test)]
mod convert_sphere_tests {
    use super::*;
    use crate::geom::{SphericalSurface, Surface3, SurfaceEval};
    use crate::base::convert::surface_to_bspline;

    #[test]
    fn sphere_to_bspline_evaluates() {
        let sphere = SphericalSurface::new(DVec3::ZERO, DVec3::Z, 5.0);
        let bspline = surface_to_bspline(&Surface3::Sphere(sphere), 10, 10);
        for u in [0.0, 0.25, 0.5, 0.75, 1.0] {
            for v in [0.0, 0.25, 0.5, 0.75, 1.0] {
                assert!(
                    bspline.point_at(u, v).is_finite(),
                    "sphere BSpline point at ({u},{v}) should be finite"
                );
            }
        }
    }
}
