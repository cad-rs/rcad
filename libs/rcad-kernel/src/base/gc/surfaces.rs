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
// gce_MakeCone
// ============================================================================

/// OCCT gce_MakeCone::gce_MakeCone(const gp_Pnt& P1, const gp_Pnt& P2,
/// const double R1, const double R2) (gce_MakeCone.cxx L228-280): the cone
/// defined by two points (axis) and two radii (section radii at each point).
///
/// The error ladder keeps the OCCT order (NullAxis on a sub-Resolution
/// distance — note the STRICT `<` here, the gce_MakeDir form uses `<=` —
/// then NegativeRadius, then NullAngle).  The D2 reference direction is the
/// OCCT L258-270 perpendicular selection over (x, y, z); for a normalized
/// D1 at least one component exceeds gp::Resolution, so the trailing
/// unreachable arm mirrors the OCCT uninitialized-D2 UB.
pub fn make_cone_2p_2r(p1: Point3, p2: Point3, r1: f64, r2: f64) -> Result<ConicalSurface, GceError> {
    const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;
    let dist = p1.distance(p2);
    if dist < GP_RESOLUTION {
        return Err(GceError::NullAxis);
    }
    if r1 < 0.0 || r2 < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    // OCCT L241: Angle = std::abs(atan((R1 - R2) / dist)).
    let mut angle = ((r1 - r2) / dist).atan().abs();
    if (std::f64::consts::FRAC_PI_2 - angle).abs() < GP_RESOLUTION || angle.abs() < GP_RESOLUTION
    {
        return Err(GceError::NullAngle);
    }
    // OCCT L243: D1(P2.XYZ() - P1.XYZ()) — normalized by the gp_Dir ctor.
    let d1 = (p2 - p1).normalize_or_zero();
    // OCCT L245-270: the D2 perpendicular selection.
    let x = d1.x;
    let y = d1.y;
    let z = d1.z;
    let d2 = if x.abs() > GP_RESOLUTION {
        Vec3::new(-y, x, 0.0).normalize_or_zero()
    } else if y.abs() > GP_RESOLUTION {
        Vec3::new(0.0, -z, y).normalize_or_zero()
    } else if z.abs() > GP_RESOLUTION {
        Vec3::new(z, 0.0, -x).normalize_or_zero()
    } else {
        unreachable!("a normalized D1 has a component above gp::Resolution");
    };
    // OCCT L272-275: R1 > R2 flips the (signed) semi-angle.
    if r1 > r2 {
        angle = -angle;
    }
    // OCCT L276: gp_Cone(gp_Ax2(P1, D1, D2), Angle, R1) — the Ax2 location
    // P1 is the reference point at radius R1 (the rcad `apex` field), D1 the
    // axis, D2 the X direction.
    Ok(ConicalSurface {
        apex: p1,
        axis: d1,
        radius: r1,
        half_angle_rad: angle,
        ref_dir: d2,
    })
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod gce_cone_tests {
    use super::*;

    // The OCCT (P1, P2, R1, R2) ctor: the cone of axis P1->P2, R1 at P1 and
    // R2 at P2; the semi-angle atan((R1-R2)/dist) is SIGNED POSITIVE when
    // R1 < R2 (the cone opens along +axis) and flips when R1 > R2 (cxx
    // L272-275).  The D2 reference direction is the (z,0,-x) arm for a
    // +Z axis.
    #[test]
    fn make_cone_2p_2r_opens_along_axis_and_selects_d2() {
        let cone = make_cone_2p_2r(DVec3::ZERO, DVec3::new(0.0, 0.0, 2.0), 1.0, 2.0).unwrap();
        let expected: f64 = (0.5f64).atan().abs();
        assert!((cone.half_angle_rad - expected).abs() < 1e-15);
        assert!((cone.axis - DVec3::Z).length() < 1e-15);
        assert!((cone.ref_dir - DVec3::X).length() < 1e-15);
        assert_eq!(cone.radius, 1.0);
        assert_eq!(cone.apex, DVec3::ZERO);

        // R1 > R2: the semi-angle flips (cxx L272-275).
        let cone2 = make_cone_2p_2r(DVec3::ZERO, DVec3::new(0.0, 0.0, 2.0), 2.0, 1.0).unwrap();
        assert!((cone2.half_angle_rad + expected).abs() < 1e-15);

        // Error ladder: sub-Resolution distance -> NullAxis; negative
        // radius -> NegativeRadius; degenerate angle -> NullAngle.
        assert!(matches!(
            make_cone_2p_2r(DVec3::ZERO, DVec3::ZERO, 1.0, 2.0),
            Err(GceError::NullAxis)
        ));
        assert!(matches!(
            make_cone_2p_2r(DVec3::ZERO, DVec3::new(0.0, 0.0, 2.0), -1.0, 2.0),
            Err(GceError::NegativeRadius)
        ));
        assert!(matches!(
            make_cone_2p_2r(DVec3::ZERO, DVec3::new(0.0, 0.0, 2.0), 1.0, 1.0),
            Err(GceError::NullAngle)
        ));
    }
}

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
            y_dir: None,
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
