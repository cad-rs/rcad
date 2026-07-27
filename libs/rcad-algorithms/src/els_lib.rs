//! ElSLib-style elementary surface utilities.
//!
//! Provides utilities for elementary surfaces analogous to OCCT `ElSLib` package.
//! Includes evaluation, parameter computation, and differential properties
//! for Plane, Cylinder, Sphere, Cone, Torus, and BSplineSurface.
//!
//! # Overview
//!
//! For each elementary surface type, this module provides:
//! - `point_at(surf, u, v)`: Compute 3D point from surface parameters
//! - `parameters(surf, point)`: Compute (u, v) parameters from a 3D point
//! - `normal(surf, u, v)`: Surface normal at (u, v)
//! - `tangent_u(surf, u, v)`: Partial derivative dS/du (u-tangent)
//! - `tangent_v(surf, u, v)`: Partial derivative dS/dv (v-tangent)
//!
//! # Coordinate Conventions
//!
//! - Plane: u and v are Cartesian coordinates in the plane's local frame
//! - Cylinder: u = azimuth angle [0, 2*pi], v = height along axis
//! - Sphere: u = longitude [0, 2*pi], v = colatitude [0, pi] (0 = north pole)
//! - Cone: u = azimuth [0, 2*pi], v = slant distance from reference circle
//! - Torus: u = major angle [0, 2*pi], v = minor angle [0, 2*pi]

use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    BSplineSurface, ConicalSurface, CylindricalSurface, Plane, SphericalSurface, SurfaceEval,
    ToroidalSurface, any_perpendicular,
};
use std::f64::consts::PI;

// =============================================================================
// Plane Utilities
// =============================================================================

/// Compute the 3D point on a plane at parameters (u, v).
///
/// The plane is parameterized as: P(u, v) = origin + u * x_axis + v * y_axis
/// where x_axis and y_axis form an orthonormal basis with the normal.
pub fn plane_point_at(plane: &Plane, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(plane.normal);
    let y_ax = plane.normal.cross(x_ax);
    plane.origin + u * x_ax + v * y_ax
}

/// Compute the (u, v) parameters for a point on or near a plane.
///
/// Projects the point onto the plane and returns the local coordinates.
/// The result satisfies: `plane_point_at(plane, u, v) == point.projected_onto_plane`.
pub fn plane_parameters(plane: &Plane, point: DVec3) -> DVec2 {
    let x_ax = any_perpendicular(plane.normal);
    let y_ax = plane.normal.cross(x_ax);
    let d = point - plane.origin;
    DVec2::new(d.dot(x_ax), d.dot(y_ax))
}

/// Get the normal vector of a plane (constant across the surface).
pub fn plane_normal(plane: &Plane) -> DVec3 {
    plane.normal
}

/// Get the u-tangent vector of a plane (constant across the surface).
///
/// This is the partial derivative dS/du, pointing along the local x-axis.
pub fn plane_tangent_u(plane: &Plane) -> DVec3 {
    any_perpendicular(plane.normal)
}

/// Get the v-tangent vector of a plane (constant across the surface).
///
/// This is the partial derivative dS/dv, pointing along the local y-axis.
pub fn plane_tangent_v(plane: &Plane) -> DVec3 {
    let x_ax = any_perpendicular(plane.normal);
    plane.normal.cross(x_ax)
}

// =============================================================================
// Cylinder Utilities
// =============================================================================

/// Compute the 3D point on a cylinder at parameters (u, v).
///
/// - u: azimuth angle [0, 2*pi]
/// - v: height along the cylinder axis
///
/// P(u, v) = origin + radius * (cos(u) * x_axis + sin(u) * y_axis) + v * axis
pub fn cylinder_point_at(cyl: &CylindricalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = cyl.ref_dir.normalize();
    let y_ax = cyl.axis.cross(x_ax).normalize();
    cyl.origin + cyl.radius * (u.cos() * x_ax + u.sin() * y_ax) + v * cyl.axis
}

/// Compute the (u, v) parameters for a point on or near a cylinder.
///
/// Projects the point onto the cylinder surface and returns the corresponding
/// parameters. The u angle is normalized to [0, 2*pi).
pub fn cylinder_parameters(cyl: &CylindricalSurface, point: DVec3) -> DVec2 {
    let axis = cyl.axis.normalize_or_zero();
    let x_ax = cyl.ref_dir.normalize();
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from cylinder origin to point
    let d = point - cyl.origin;

    // Height along axis (v parameter)
    let v = d.dot(axis);

    // Radial component
    let radial = d - v * axis;
    // atan2(y, x) where y = radial.dot(y_ax), x = radial.dot(x_ax)
    let u = radial.dot(y_ax).atan2(radial.dot(x_ax));

    // Normalize u to [0, 2*pi)
    let u = if u < 0.0 { u + 2.0 * PI } else { u };

    DVec2::new(u, v)
}

/// Get the normal vector of a cylinder at parameters (u, v).
///
/// The normal points outward from the axis (radially).
pub fn cylinder_normal(cyl: &CylindricalSurface, u: f64, _v: f64) -> DVec3 {
    let x_ax = cyl.ref_dir.normalize();
    let y_ax = cyl.axis.cross(x_ax).normalize();
    (u.cos() * x_ax + u.sin() * y_ax).normalize()
}

/// Get the u-tangent vector of a cylinder at parameters (u, v).
///
/// This is the azimuthal tangent (along the circular cross-section).
pub fn cylinder_tangent_u(cyl: &CylindricalSurface, u: f64, _v: f64) -> DVec3 {
    let x_ax = cyl.ref_dir.normalize();
    let y_ax = cyl.axis.cross(x_ax).normalize();
    (-u.sin() * x_ax + u.cos() * y_ax).normalize()
}

/// Get the v-tangent vector of a cylinder at parameters (u, v).
///
/// This is the axial tangent (along the cylinder axis).
pub fn cylinder_tangent_v(cyl: &CylindricalSurface, _u: f64, _v: f64) -> DVec3 {
    cyl.axis.normalize_or_zero()
}

// =============================================================================
// Sphere Utilities
// =============================================================================

/// Compute the 3D point on a sphere at parameters (u, v).
///
/// - u: longitude angle [0, 2*pi]
/// - v: colatitude angle [0, pi] (0 = north pole, pi = south pole)
///
/// P(u, v) = center + radius * (sin(v) * (cos(u) * x + sin(u) * y) + cos(v) * axis)
pub fn sphere_point_at(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = sph.ref_dir.normalize();
    let y_ax = sph.axis.cross(x_ax).normalize();
    sph.center + sph.radius * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * sph.axis)
}

/// Compute the (u, v) parameters for a point on or near a sphere.
///
/// Projects the point onto the sphere surface and returns the corresponding
/// angles. At the poles, u is set to 0.0.
pub fn sphere_parameters(sph: &SphericalSurface, point: DVec3) -> DVec2 {
    let axis = sph.axis.normalize_or_zero();
    let x_ax = sph.ref_dir.normalize();
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from center to point (normalized)
    let d = (point - sph.center).normalize_or_zero();

    // Colatitude: angle from axis
    let cos_v = d.dot(axis);
    let v = cos_v.clamp(-1.0, 1.0).acos();

    // Longitude: angle in the equatorial plane
    let sin_v = v.sin();
    let u = if sin_v.abs() > TOLERANCE_LINEAR_ULTRA_STRICT {
        let radial = d - cos_v * axis;
        let u_raw = radial.dot(y_ax).atan2(radial.dot(x_ax));
        if u_raw < 0.0 { u_raw + 2.0 * PI } else { u_raw }
    } else {
        0.0 // At poles, u is undefined; use 0
    };

    DVec2::new(u, v)
}

/// Get the normal vector of a sphere at parameters (u, v).
///
/// The normal points outward from the center.
pub fn sphere_normal(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = sph.ref_dir.normalize();
    let y_ax = sph.axis.cross(x_ax).normalize();
    (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * sph.axis).normalize()
}

/// Get the u-tangent vector of a sphere at parameters (u, v).
///
/// This is the longitude tangent (along lines of latitude).
pub fn sphere_tangent_u(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = sph.ref_dir.normalize();
    let y_ax = sph.axis.cross(x_ax).normalize();
    v.sin() * (-u.sin() * x_ax + u.cos() * y_ax)
}

/// Get the v-tangent vector of a sphere at parameters (u, v).
///
/// This is the colatitude tangent (along meridians).
pub fn sphere_tangent_v(sph: &SphericalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = sph.ref_dir.normalize();
    let y_ax = sph.axis.cross(x_ax).normalize();
    sph.radius * (v.cos() * (u.cos() * x_ax + u.sin() * y_ax) - v.sin() * sph.axis)
}

// =============================================================================
// Cone Utilities
// =============================================================================

/// Compute the 3D point on a cone at parameters (u, v).
///
/// - u: azimuth angle [0, 2*pi]
/// - v: slant distance from the reference circle at apex
///
/// The reference circle has radius `cone.radius` at the apex point.
/// Positive v moves toward larger radius if half_angle > 0.
pub fn cone_point_at(cone: &ConicalSurface, u: f64, v: f64) -> DVec3 {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();
    let radial = cone.radius_at_slant(v);
    let axial = cone.axial_from_slant(v);
    cone.apex + axial * axis + radial * (u.cos() * x_ax + u.sin() * y_ax)
}

/// Compute the (u, v) parameters for a point on or near a cone.
///
/// Projects the point onto the cone surface and returns the corresponding
/// parameters.
pub fn cone_parameters(cone: &ConicalSurface, point: DVec3) -> DVec2 {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from apex to point
    let d = point - cone.apex;

    // Axial distance from apex
    let axial = d.dot(axis);

    // Radial component
    let radial_vec = d - axial * axis;
    let radial_dist = radial_vec.length();

    // Azimuth angle
    let u = if radial_dist > TOLERANCE_LINEAR_ULTRA_STRICT {
        let u_raw = radial_vec.dot(y_ax).atan2(radial_vec.dot(x_ax));
        if u_raw < 0.0 { u_raw + 2.0 * PI } else { u_raw }
    } else {
        0.0
    };

    // Slant distance from reference circle
    // At the reference circle (v=0), axial = 0 and radial = cone.radius
    // v is the distance along the cone surface from this reference
    let slant_from_apex = (axial * axial + radial_dist * radial_dist).sqrt();
    let ref_slant = if cone.half_angle_rad.tan().abs() > TOLERANCE_LINEAR_ULTRA_STRICT {
        cone.radius / cone.half_angle_rad.sin()
    } else {
        0.0
    };
    let v = slant_from_apex - ref_slant;

    DVec2::new(u, v)
}

/// Get the normal vector of a cone at parameters (u, v).
///
/// The normal is constant along lines of constant u (generators).
pub fn cone_normal(cone: &ConicalSurface, u: f64, _v: f64) -> DVec3 {
    let axis = cone.axis_dir();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();
    let radial = u.cos() * x_ax + u.sin() * y_ax;
    let half = cone.half_angle_rad;
    (radial * half.cos() - axis * half.sin()).normalize()
}

// =============================================================================
// Torus Utilities
// =============================================================================

/// Compute the 3D point on a torus at parameters (u, v).
///
/// - u: major angle [0, 2*pi] (angle around the main axis)
/// - v: minor angle [0, 2*pi] (angle around the tube)
///
/// The torus is centered at `center` with the main axis `axis`.
/// Major radius is the distance from center to tube center.
/// Minor radius is the tube radius.
pub fn torus_point_at(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();

    // Center of the tube cross-section at angle u
    let tube_center = torus.center + torus.major_radius * (u.cos() * x_ax + u.sin() * y_ax);

    // Radial direction from main axis to tube center
    let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();

    // Point on the tube surface
    tube_center + torus.minor_radius * (v.cos() * radial + v.sin() * torus.axis)
}

/// Compute the (u, v) parameters for a point on or near a torus.
///
/// Projects the point onto the torus surface and returns the corresponding
/// parameters.
pub fn torus_parameters(torus: &ToroidalSurface, point: DVec3) -> DVec2 {
    let axis = torus.axis.normalize_or_zero();
    let x_ax = any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();

    // Vector from torus center to point
    let d = point - torus.center;

    // Height above/below the equatorial plane
    let z = d.dot(axis);

    // Radial component in equatorial plane
    let radial_2d = d - z * axis;
    let r_2d = radial_2d.length();

    // Major angle u
    let u = if r_2d > TOLERANCE_LINEAR_ULTRA_STRICT {
        let u_raw = radial_2d.dot(y_ax).atan2(radial_2d.dot(x_ax));
        if u_raw < 0.0 { u_raw + 2.0 * PI } else { u_raw }
    } else {
        0.0
    };

    // Minor angle v
    // The tube center at angle u is at distance major_radius from axis
    // dv is the distance from the tube center in the radial direction
    let dv = r_2d - torus.major_radius;

    let v = if dv.abs() > TOLERANCE_LINEAR_ULTRA_STRICT || z.abs() > TOLERANCE_LINEAR_ULTRA_STRICT {
        let v_raw = z.atan2(dv);
        if v_raw < 0.0 { v_raw + 2.0 * PI } else { v_raw }
    } else {
        0.0
    };

    DVec2::new(u, v)
}

/// Get the normal vector of a torus at parameters (u, v).
///
/// The normal points outward from the tube surface.
pub fn torus_normal(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();
    let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
    (v.cos() * radial + v.sin() * torus.axis).normalize()
}

/// Get the u-tangent vector of a torus at parameters (u, v).
///
/// This is the tangent along the major circle (around the main axis).
pub fn torus_tangent_u(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();

    // Derivative of tube center w.r.t. u
    let dcenter_du = torus.major_radius * (-u.sin() * x_ax + u.cos() * y_ax);

    // Derivative of radial direction w.r.t. u
    let dradial_du = -u.sin() * x_ax + u.cos() * y_ax;

    dcenter_du + torus.minor_radius * v.cos() * dradial_du
}

/// Get the v-tangent vector of a torus at parameters (u, v).
///
/// This is the tangent along the minor circle (around the tube).
pub fn torus_tangent_v(torus: &ToroidalSurface, u: f64, v: f64) -> DVec3 {
    let x_ax = any_perpendicular(torus.axis);
    let y_ax = torus.axis.cross(x_ax).normalize();
    let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
    torus.minor_radius * (-v.sin() * radial + v.cos() * torus.axis)
}

// =============================================================================
// BSplineSurface Utilities
// =============================================================================

/// Compute the 3D point on a BSpline surface at parameters (u, v).
///
/// Uses tensor-product NURBS evaluation via the `SurfaceEval` trait.
pub fn bspline_surface_point_at(surf: &BSplineSurface, u: f64, v: f64) -> DVec3 {
    surf.point_at(u, v)
}

/// Compute the normal vector of a BSpline surface at parameters (u, v).
///
/// Uses finite differences to compute the cross product of partial derivatives.
pub fn bspline_surface_normal(surf: &BSplineSurface, u: f64, v: f64) -> DVec3 {
    surf.normal_at(u, v)
}

/// Compute the first partial derivatives of a BSpline surface at (u, v).
///
/// Returns `[du, dv, dudv]` where:
/// - `du` is the partial derivative with respect to u (dS/du)
/// - `dv` is the partial derivative with respect to v (dS/dv)
/// - `dudv` is the mixed second derivative (d2S/dudv)
///
/// Uses central finite differences for numerical stability.
pub fn bspline_surface_derivatives(surf: &BSplineSurface, u: f64, v: f64) -> [DVec3; 3] {
    let eps = TOLERANCE_RETRY_LADDER_MID;
    let [u0, u1, v0, v1] = surf.default_domain();

    // Clamp to domain bounds
    let u_minus = (u - eps).max(u0);
    let u_plus = (u + eps).min(u1);
    let v_minus = (v - eps).max(v0);
    let v_plus = (v + eps).min(v1);

    // First derivatives using central differences where possible
    let du = if u_plus > u_minus {
        (surf.point_at(u_plus, v) - surf.point_at(u_minus, v)) / (u_plus - u_minus)
    } else if u_plus > u0 {
        (surf.point_at(u_plus, v) - surf.point_at(u, v)) / (u_plus - u)
    } else if u_minus < u1 {
        (surf.point_at(u, v) - surf.point_at(u_minus, v)) / (u - u_minus)
    } else {
        DVec3::ZERO
    };

    let dv = if v_plus > v_minus {
        (surf.point_at(u, v_plus) - surf.point_at(u, v_minus)) / (v_plus - v_minus)
    } else if v_plus > v0 {
        (surf.point_at(u, v_plus) - surf.point_at(u, v)) / (v_plus - v)
    } else if v_minus < v1 {
        (surf.point_at(u, v) - surf.point_at(u, v_minus)) / (v - v_minus)
    } else {
        DVec3::ZERO
    };

    // Mixed second derivative using finite differences
    let dudv = if u_plus > u_minus && v_plus > v_minus {
        let p_pp = surf.point_at(u_plus, v_plus);
        let p_pm = surf.point_at(u_plus, v_minus);
        let p_mp = surf.point_at(u_minus, v_plus);
        let p_mm = surf.point_at(u_minus, v_minus);

        ((p_pp - p_pm) - (p_mp - p_mm)) / ((u_plus - u_minus) * (v_plus - v_minus))
    } else {
        DVec3::ZERO
    };

    [du, dv, dudv]
}

// =============================================================================
// Tests
// =============================================================================
