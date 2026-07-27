//! GC_Make/gce_Make-style geometry construction utilities.
//!
//! Provides simple geometry construction functions analogous to OCCT's
//! GC (Geom_ Curve/Surface construction) and gce (gp_ primitive construction)
//! packages.
//!
//! GC_MakeSegment2d, GC_MakePlane, GC_MakeCircle2d,
//!                  GC_MakeArcOfCircle, GC_MakeParabola2d, GC_MakeConicalSurface,
//!                  gce_MakeHypr, gce_MakeElips, gce_MakeCylinder,
//!                  gce_MakeCone, gce_MakeCirc2d

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;

/// Error status for geometry construction (equivalent to OCCT gce_ErrorType).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MakeError {
    Ok,
    ConfusedPoints, // points too close
    BadAngle,       // invalid angle
    NullAxis,
    NullRadius,
    Other,
}

// =============================================================================
// 2D Segment (GC_MakeSegment2d)
// =============================================================================

/// Create a 2D line segment from two points. Returns error if points are coincident.
pub fn make_segment2d(p1: DVec2, p2: DVec2) -> Result<Line2d, MakeError> {
    let dir = p2 - p1;
    let len = dir.length();
    if len < 1e-15 {
        return Err(MakeError::ConfusedPoints);
    }
    Ok(Line2d {
        origin: p1,
        direction: dir / len,
    })
}

// =============================================================================
// Plane (GC_MakePlane)
// =============================================================================

/// Create an offset plane from a base plane by a distance along the normal.
pub fn make_plane_offset(plane: &Plane, distance: f64) -> Plane {
    Plane::new(plane.origin + distance * plane.normal, plane.normal)
}

// =============================================================================
// 2D Circle (gce_MakeCirc2d / GC_MakeCircle2d)
// =============================================================================

/// Create a 2D circle from three points. Returns error if points are collinear.
pub fn make_circle2d_3points(p1: DVec2, p2: DVec2, p3: DVec2) -> Result<Circle2d, MakeError> {
    let a = p2 - p1;
    let b = p3 - p1;
    let cross = a.x * b.y - a.y * b.x;
    if cross.abs() < 1e-15 {
        return Err(MakeError::ConfusedPoints);
    }
    let a2 = a.length_squared();
    let b2 = b.length_squared();
    let inv = 0.5 / cross;
    let cx = p1.x + (b2 * a.y - a2 * b.y) * inv;
    let cy = p1.y + (a2 * b.x - b2 * a.x) * inv;
    let center = DVec2::new(cx, cy);
    let radius = (p1 - center).length();
    if radius < 1e-15 {
        return Err(MakeError::NullRadius);
    }
    Ok(Circle2d::new(center, radius))
}

/// Create a 2D circle concentric with a base circle passing through a point.
pub fn make_circle2d_concentric(base: &Circle2d, point: DVec2) -> Result<Circle2d, MakeError> {
    let r = (point - base.center).length();
    if r < 1e-15 {
        return Err(MakeError::NullRadius);
    }
    Ok(Circle2d::new(base.center, r))
}

// =============================================================================
// 2D Parabola (GC_MakeParabola2d)
// =============================================================================

/// Create a 2D parabola from an axis placement and focus length.
/// focal = distance from vertex to focus.
pub fn make_parabola2d(origin: DVec2, axis_dir: DVec2, focal: f64) -> Parabola2d {
    Parabola2d {
        origin,
        axis_dir: axis_dir.normalize(),
        focal_param: 2.0 * focal,
    }
}

// =============================================================================
// Conical Surface (GC_MakeConicalSurface)
// =============================================================================

/// Create a conical surface from an axis and angle.
/// Returns BadAngle if angle is below resolution.
pub fn make_conical_surface(
    apex: DVec3,
    axis: DVec3,
    radius: f64,
    half_angle_rad: f64,
) -> Result<ConicalSurface, MakeError> {
    if half_angle_rad.abs() < 1e-15 {
        return Err(MakeError::BadAngle);
    }
    Ok(ConicalSurface {
        apex,
        axis,
        radius,
        half_angle_rad,
    })
}

// =============================================================================
// 3D Circle Arc (GC_MakeArcOfCircle)
// =============================================================================

/// Create a circular arc from first point `p1`, tangent direction at p1, and end point p2.
/// Returns the path center and radius.
pub fn make_arc_of_circle(
    p1: DVec3,
    tangent_at_p1: DVec3,
    p2: DVec3,
) -> Result<(DVec3, f64, DVec3), MakeError> {
    // Find center by intersecting perpendicular bisectors
    let mid = (p1 + p2) * 0.5;
    let chord = p2 - p1;
    let mid_perp = chord.cross(tangent_at_p1);
    if mid_perp.length_squared() < 1e-30 {
        return Err(MakeError::ConfusedPoints);
    }
    // Distance along perpendicular from midpoint to center
    let d = mid_perp.dot(p1 - mid) / mid_perp.dot(mid_perp);
    let center = mid + d * mid_perp;
    let radius = (p1 - center).length();
    let normal = tangent_at_p1.cross(p2 - p1).normalize_or_zero();
    Ok((center, radius, normal))
}

// =============================================================================
// Ellipse (gce_MakeElips)
// =============================================================================

/// Create an ellipse from three points: two foci and a point on the ellipse.
pub fn make_ellipse_from_foci(
    focus1: DVec3,
    focus2: DVec3,
    point_on: DVec3,
) -> Result<Ellipse3, MakeError> {
    let center = (focus1 + focus2) * 0.5;
    let major_dir = (focus2 - focus1).normalize_or_zero();
    if major_dir.length_squared() < 1e-30 {
        return Err(MakeError::ConfusedPoints);
    }
    let c = (focus2 - focus1).length() * 0.5; // focal half-distance
    let a = 0.5 * ((point_on - focus1).length() + (point_on - focus2).length()); // semi-major
    if a <= c {
        return Err(MakeError::Other);
    }
    let b = (a * a - c * c).sqrt(); // semi-minor
    let normal = (focus2 - focus1)
        .cross(point_on - focus1)
        .normalize_or_zero();
    Ok(Ellipse3 {
        center,
        normal,
        major_dir,
        major_radius: a,
        minor_radius: b,
    })
}

// =============================================================================
// Hyperbola (gce_MakeHypr)
// =============================================================================

/// Create a hyperbola from axis and radii.
pub fn make_hyperbola(
    center: DVec3,
    normal: DVec3,
    major_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Hyperbola3 {
    Hyperbola3 {
        center,
        normal,
        major_dir,
        semi_major: major_radius,
        semi_minor: minor_radius,
    }
}

/// Create a 2D hyperbola.
pub fn make_hyperbola2d(
    center: DVec2,
    major_dir: DVec2,
    major_radius: f64,
    minor_radius: f64,
) -> Hyperbola2d {
    Hyperbola2d {
        center,
        major_dir,
        semi_major: major_radius,
        semi_minor: minor_radius,
    }
}

// =============================================================================
// Cylinder (gce_MakeCylinder)
// =============================================================================

/// Create a cylindrical surface from an axis and radius.
pub fn make_cylinder(origin: DVec3, axis: DVec3, radius: f64) -> CylindricalSurface {
    let ax = axis.normalize_or_zero();
    CylindricalSurface {
        origin,
        axis: ax,
        radius,
        ref_dir: any_perpendicular(ax),
    }
}

// =============================================================================
// Cone (gce_MakeCone)
// =============================================================================

/// Create an offset cone surface from a base cone and offset distance.
pub fn make_cone_offset(cone: &ConicalSurface, distance: f64) -> ConicalSurface {
    // OCCT: offset along the axis direction shifts the apex
    let axis = cone.axis.normalize_or_zero();
    let shift = distance / cone.half_angle_rad.sin();
    let new_apex = cone.apex + shift * axis;
    let new_radius = cone.radius + distance * cone.half_angle_rad.cos();
    ConicalSurface {
        apex: new_apex,
        axis,
        radius: new_radius,
        half_angle_rad: cone.half_angle_rad,
    }
}

// =============================================================================
// Tests — translated from OCCT GTests
// =============================================================================
