//! Surface construction algorithms.
//!
//! OCCT GC package: GC_MakePlane, GC_MakeCylindricalSurface,
//! GC_MakeConicalSurface, GC_MakeTrimmedCone, GC_MakeTrimmedCylinder.

#![allow(clippy::manual_clamp)]

use glam::DVec3;

use crate::geom::{
    ConicalSurface, CylindricalSurface, Plane, Point3, Vec3,
};

use super::{
    GceError, TOL_CONF, points_coincident,
};

// ============================================================================
// GC_MakePlane
// ============================================================================

/// Construct a plane from a point and normal direction.
///
/// OCCT: `GC_MakePlane(gp_Pnt, gp_Dir)`.
pub fn make_plane_pn(point: Point3, normal: Vec3) -> Result<Plane, GceError> {
    let n = normal.normalize_or_zero();
    if n.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Plane::new(point, n))
}

/// Construct a plane from three points.
///
/// OCCT: `GC_MakePlane(gp_Pnt, gp_Pnt, gp_Pnt)`.
/// Returns `GceError::ConfusedPoints` or `GceError::ColinearPoints` for degenerate input.
pub fn make_plane_3p(p1: Point3, p2: Point3, p3: Point3) -> Result<Plane, GceError> {
    if points_coincident(p1, p2) || points_coincident(p1, p3) || points_coincident(p2, p3) {
        return Err(GceError::ConfusedPoints);
    }
    let d1 = p2 - p1;
    let d2 = p3 - p1;
    let normal = d1.cross(d2);
    if normal.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ColinearPoints);
    }
    Ok(Plane::new(p1, normal))
}

/// Construct a plane from its cartesian equation: A*x + B*y + C*z + D = 0.
///
/// OCCT: `GC_MakePlane(double, double, double, double)`.
/// Returns `GceError::BadEquation` when sqrt(A² + B² + C²) is below resolution.
pub fn make_plane_abcd(a: f64, b: f64, c: f64, d: f64) -> Result<Plane, GceError> {
    let norm_sq = a * a + b * b + c * c;
    if norm_sq < TOL_CONF * TOL_CONF {
        return Err(GceError::BadEquation);
    }
    let inv_norm = 1.0 / norm_sq.sqrt();
    let normal = DVec3::new(a, b, c) * inv_norm;
    // D is signed distance from origin along normal (with sign convention)
    let origin = normal * (-d * inv_norm);
    Ok(Plane::new(origin, normal))
}

/// Construct a plane parallel to an existing plane passing through a point.
///
/// OCCT: `GC_MakePlane(gp_Pln, gp_Pnt)`.
pub fn make_plane_parallel_point(plane: &Plane, point: Point3) -> Result<Plane, GceError> {
    Ok(Plane::new(point, plane.normal))
}

/// Construct a plane parallel to an existing plane at a signed distance.
///
/// OCCT: `GC_MakePlane(gp_Pln, double)`.
/// Positive distance follows the normal direction.
pub fn make_plane_parallel_dist(plane: &Plane, dist: f64) -> Result<Plane, GceError> {
    if dist.abs() < TOL_CONF {
        return Err(GceError::ZeroDistance);
    }
    let origin = plane.origin + plane.normal * dist;
    Ok(Plane::new(origin, plane.normal))
}

/// Construct a plane through an axis (location point on plane + direction defines normal).
///
/// OCCT: `GC_MakePlane(gp_Ax1)`.
/// The axis location lies on the plane; the axis direction is the plane normal.
pub fn make_plane_axis(origin: Point3, axis_dir: Vec3) -> Result<Plane, GceError> {
    let normal = axis_dir.normalize_or_zero();
    if normal.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(Plane::new(origin, normal))
}

// ============================================================================
// GC_MakeCylindricalSurface
// ============================================================================

/// Construct a cylindrical surface from origin, axis, and radius.
///
/// OCCT: `GC_MakeCylindricalSurface(gp_Ax2, double)`.
/// radius must be >= 0.
pub fn make_cylindrical_surface(
    origin: Point3,
    axis: Vec3,
    radius: f64,
    ref_dir: Vec3,
) -> Result<CylindricalSurface, GceError> {
    if radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let axis = axis.normalize_or_zero();
    if axis.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    let ref_dir = ref_dir.normalize_or_zero();
    if ref_dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(CylindricalSurface {
        origin,
        axis,
        radius,
        ref_dir,
        y_dir: None,
    })
}

/// Construct a cylindrical surface from an existing one, parallel at a point.
///
/// OCCT: `GC_MakeCylindricalSurface(gp_CylindricalSurface, gp_Pnt)`.
/// Creates a cylinder with same axis and radius passing through the point.
pub fn make_cylindrical_surface_point(
    cyl: &CylindricalSurface,
    point: Point3,
) -> Result<CylindricalSurface, GceError> {
    let d = point - cyl.origin;
    let along_axis = d.dot(cyl.axis);
    let radial = d - cyl.axis * along_axis;
    let radius = radial.length();
    if radius < TOL_CONF {
        return Err(GceError::NullRadius);
    }
    Ok(CylindricalSurface {
        origin: cyl.origin,
        axis: cyl.axis,
        radius,
        ref_dir: cyl.ref_dir,
        y_dir: cyl.y_dir,
    })
}

// ============================================================================
// GC_MakeConicalSurface
// ============================================================================

/// Construct a conical surface from apex, axis, reference radius, and half-angle.
///
/// OCCT: `GC_MakeConicalSurface(gp_Ax2, double, double)`.
/// The `apex` is the point where radius = 0; `ref_radius` is the radius at `apex` + 1 unit
/// along the axis (i.e., semi-angle = atan(ref_radius / 1.0)).
pub fn make_conical_surface(
    apex: Point3,
    axis: Vec3,
    radius: f64,
    half_angle_rad: f64,
) -> Result<ConicalSurface, GceError> {
    if half_angle_rad < 0.0 {
        return Err(GceError::NullAngle);
    }
    if radius < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    let axis = axis.normalize_or_zero();
    if axis.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(ConicalSurface::new(apex, axis, radius, half_angle_rad))
}

// ============================================================================
// GC_MakeTrimmedCylinder
// ============================================================================

/// Result of trimming a cylinder: stores the cylinder with axial bounds.
#[derive(Debug, Clone, Copy)]
pub struct TrimmedCylinder {
    /// The underlying cylinder surface.
    pub cylinder: CylindricalSurface,
    /// Lower axial bound (v parameter).
    pub v_min: f64,
    /// Upper axial bound (v parameter).
    pub v_max: f64,
}

/// Construct a trimmed cylinder from an existing cylinder and axial range.
///
/// OCCT: `GC_MakeTrimmedCylinder(gp_CylindricalSurface, double, double)`.
/// Creates a cylinder restricted to the axial range [v_min, v_max].
pub fn make_trimmed_cylinder(
    cyl: &CylindricalSurface,
    v_min: f64,
    v_max: f64,
) -> Result<TrimmedCylinder, GceError> {
    if (v_max - v_min).abs() < TOL_CONF {
        return Err(GceError::NullLength);
    }
    let (v1, v2) = if v_min < v_max { (v_min, v_max) } else { (v_max, v_min) };
    Ok(TrimmedCylinder {
        cylinder: *cyl,
        v_min: v1,
        v_max: v2,
    })
}

// ============================================================================
// GC_MakeTrimmedCone
// ============================================================================

/// Result of trimming a cone: stores the cone with axial bounds.
#[derive(Debug, Clone, Copy)]
pub struct TrimmedCone {
    /// The underlying conical surface.
    pub cone: ConicalSurface,
    /// Lower axial bound (v parameter, along axis from apex).
    pub v_min: f64,
    /// Upper axial bound (v parameter, along axis from apex).
    pub v_max: f64,
}

/// Construct a trimmed cone from an existing cone and axial range.
///
/// OCCT: `GC_MakeTrimmedCone(gp_ConicalSurface, double, double)`.
/// Creates a cone restricted to the axial range [v_min, v_max].
/// v_min and v_max are distances from apex along the axis.
pub fn make_trimmed_cone(
    cone: &ConicalSurface,
    v_min: f64,
    v_max: f64,
) -> Result<TrimmedCone, GceError> {
    if (v_max - v_min).abs() < TOL_CONF {
        return Err(GceError::NullLength);
    }
    let (v1, v2) = if v_min < v_max { (v_min, v_max) } else { (v_max, v_min) };
    Ok(TrimmedCone {
        cone: *cone,
        v_min: v1,
        v_max: v2,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_plane_pn() {
        let p = make_plane_pn(DVec3::new(1.0, 2.0, 3.0), DVec3::Z).unwrap();
        assert!((p.normal - DVec3::Z).length() < 1e-12);
        assert!((p.origin - DVec3::new(1.0, 2.0, 3.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_plane_3p() {
        let p = make_plane_3p(
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        )
        .unwrap();
        assert!((p.normal - DVec3::Z).length() < 1e-12);
    }

    #[test]
    fn test_make_plane_3p_collinear() {
        assert_eq!(
            make_plane_3p(
                DVec3::ZERO,
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(2.0, 0.0, 0.0),
            )
            .unwrap_err(),
            GceError::ColinearPoints
        );
    }

    #[test]
    fn test_make_plane_abcd() {
        // z = 1 => 1*z - 1 = 0 => A=0, B=0, C=1, D=-1
        let p = make_plane_abcd(0.0, 0.0, 1.0, -1.0).unwrap();
        assert!((p.normal - DVec3::Z).length() < 1e-12);
        assert!((p.origin - DVec3::new(0.0, 0.0, 1.0)).length() < 1e-10);
    }

    #[test]
    fn test_make_plane_bad_equation() {
        assert_eq!(
            make_plane_abcd(0.0, 0.0, 0.0, 1.0).unwrap_err(),
            GceError::BadEquation
        );
    }

    #[test]
    fn test_make_cylindrical_surface() {
        let c = make_cylindrical_surface(
            DVec3::ZERO,
            DVec3::Z,
            5.0,
            DVec3::X,
        )
        .unwrap();
        assert!((c.radius - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_conical_surface() {
        let c = make_conical_surface(
            DVec3::ZERO,
            DVec3::Z,
            2.0,
            0.5,
        )
        .unwrap();
        assert!((c.half_angle_rad - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_make_trimmed_cylinder() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 5.0,
            ref_dir: DVec3::X,
        };
        let tc = make_trimmed_cylinder(&cyl, 0.0, 10.0).unwrap();
        assert!((tc.v_min - 0.0).abs() < 1e-12);
        assert!((tc.v_max - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_make_trimmed_cone() {
        let cone = ConicalSurface::new(DVec3::ZERO, DVec3::Z, 0.0, 0.25);
        let tc = make_trimmed_cone(&cone, 1.0, 5.0).unwrap();
        assert!((tc.v_min - 1.0).abs() < 1e-12);
        assert!((tc.v_max - 5.0).abs() < 1e-12);
    }
}
