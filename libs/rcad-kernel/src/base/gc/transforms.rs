//! Geometric transformation construction algorithms.
//!
//! OCCT GC package: GC_MakeMirror, GC_MakeRotation, GC_MakeScale,
//! GC_MakeTranslation (3D and 2D variants).

use glam::{DAffine2, DAffine3, DMat2, DMat3, DVec2, DVec3};

use super::{GceError, TOL_ANG, TOL_CONF};

// ============================================================================
// Transformation helpers
// ============================================================================

/// Constructs a 3D rotation from a rotation axis (origin + direction) and angle.
fn rotation_3d(origin: DVec3, axis_dir: DVec3, angle_rad: f64) -> DAffine3 {
    // OCCT gp_Trsf::SetRotation(gp_Ax1, angle)
    // Use Rodrigues' rotation formula
    let k = axis_dir.normalize();
    let cos_a = angle_rad.cos();
    let sin_a = angle_rad.sin();
    let one_minus_cos = 1.0 - cos_a;

    // Rotation matrix components
    let kx = k.x;
    let ky = k.y;
    let kz = k.z;

    let rot = DMat3::from_cols(
        DVec3::new(
            cos_a + kx * kx * one_minus_cos,
            kx * ky * one_minus_cos + kz * sin_a,
            kx * kz * one_minus_cos - ky * sin_a,
        ),
        DVec3::new(
            ky * kx * one_minus_cos - kz * sin_a,
            cos_a + ky * ky * one_minus_cos,
            ky * kz * one_minus_cos + kx * sin_a,
        ),
        DVec3::new(
            kz * kx * one_minus_cos + ky * sin_a,
            kz * ky * one_minus_cos - kx * sin_a,
            cos_a + kz * kz * one_minus_cos,
        ),
    );

    // Translate to origin, rotate, translate back
    let t = origin;
    let translation = t - rot * t;
    DAffine3 {
        matrix3: rot,
        translation,
    }
}

/// Constructs a 2D rotation around a center point.
fn rotation_2d(center: DVec2, angle_rad: f64) -> DAffine2 {
    let cos_a = angle_rad.cos();
    let sin_a = angle_rad.sin();
    let rot = DMat2::from_cols(
        DVec2::new(cos_a, sin_a),
        DVec2::new(-sin_a, cos_a),
    );
    // Translate to origin, rotate, translate back
    let t = center;
    let translation = t - rot * t;
    DAffine2 {
        matrix2: rot,
        translation,
    }
}

/// Constructs a 3D mirror (reflection) across a plane defined by point and normal.
fn mirror_3d(origin: DVec3, normal: DVec3) -> DAffine3 {
    let n = normal.normalize();
    // Householder reflection: M = I - 2*n*n^T
    let nx = n.x;
    let ny = n.y;
    let nz = n.z;
    let two_nn = 2.0;

    let reflection = DMat3::from_cols(
        DVec3::new(1.0 - two_nn * nx * nx, -two_nn * nx * ny, -two_nn * nx * nz),
        DVec3::new(-two_nn * ny * nx, 1.0 - two_nn * ny * ny, -two_nn * ny * nz),
        DVec3::new(-two_nn * nz * nx, -two_nn * nz * ny, 1.0 - two_nn * nz * nz),
    );

    // Translate origin to the mirror plane
    let t = origin;
    let translation = t - reflection * t;
    DAffine3 {
        matrix3: reflection,
        translation,
    }
}

/// Constructs a 2D mirror (reflection) across a line defined by origin and direction.
fn mirror_2d(origin: DVec2, direction: DVec2) -> DAffine2 {
    let d = direction.normalize();
    let dx = d.x;
    let dy = d.y;
    // Reflection matrix across the line: R = 2*(d⊗d) - I
    let reflection = DMat2::from_cols(
        DVec2::new(2.0 * dx * dx - 1.0, 2.0 * dx * dy),
        DVec2::new(2.0 * dy * dx, 2.0 * dy * dy - 1.0),
    );

    let t = origin;
    let translation = t - reflection * t;
    DAffine2 {
        matrix2: reflection,
        translation,
    }
}

// ============================================================================
// GC_MakeMirror (3D)
// ============================================================================

/// Construct a mirror transformation across a point (central symmetry).
///
/// OCCT: `GC_MakeMirror(gp_Pnt)`.
/// Returns a transformation that maps P → 2*point - P.
pub fn make_mirror_point(point: DVec3) -> DAffine3 {
    DAffine3::from_translation(point)
        * DAffine3::from_scale(DVec3::splat(-1.0))
        * DAffine3::from_translation(-point)
}

/// Construct a mirror transformation across a line (axis).
///
/// OCCT: `GC_MakeMirror(gp_Ax1)` — mirror across the line defined by origin and direction.
pub fn make_mirror_axis(origin: DVec3, direction: DVec3) -> Result<DAffine3, GceError> {
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    // Mirror across a line: reflect through 180° rotation around the line
    // Equivalent to R(v) = 2*proj_line(v) - v
    Ok(rotation_3d(origin, dir, std::f64::consts::PI))
}

/// Construct a mirror transformation across a plane.
///
/// OCCT: `GC_MakeMirror(gp_Pnt, gp_Dir)` — mirror across the plane through point with normal.
pub fn make_mirror_plane(origin: DVec3, normal: DVec3) -> Result<DAffine3, GceError> {
    let n = normal.normalize_or_zero();
    if n.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(mirror_3d(origin, n))
}

// ============================================================================
// GC_MakeRotation (3D)
// ============================================================================

/// Construct a rotation transformation around an axis through the origin.
///
/// OCCT: `GC_MakeRotation(gp_Ax1, double)`.
pub fn make_rotation(
    origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
) -> Result<DAffine3, GceError> {
    let dir = axis_dir.normalize_or_zero();
    if dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    if angle_rad.abs() < TOL_ANG {
        return Err(GceError::NullAngle);
    }
    Ok(rotation_3d(origin, dir, angle_rad))
}

// ============================================================================
// GC_MakeScale (3D)
// ============================================================================

/// Construct a scaling transformation relative to a center point.
///
/// OCCT: `GC_MakeScale(gp_Pnt, double)`.
pub fn make_scale(center: DVec3, scale: f64) -> Result<DAffine3, GceError> {
    if scale < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    if scale.abs() < TOL_CONF {
        return Err(GceError::NullRadius);
    }
    Ok(DAffine3::from_translation(center)
        * DAffine3::from_scale(DVec3::splat(scale))
        * DAffine3::from_translation(-center))
}

// ============================================================================
// GC_MakeTranslation (3D)
// ============================================================================

/// Construct a translation from a vector.
///
/// OCCT: `GC_MakeTranslation(gp_Vec)`.
pub fn make_translation_vec(offset: DVec3) -> DAffine3 {
    DAffine3::from_translation(offset)
}

/// Construct a translation from two points (from P1 to P2).
///
/// OCCT: `GC_MakeTranslation(gp_Pnt, gp_Pnt)`.
pub fn make_translation_2p(p1: DVec3, p2: DVec3) -> Result<DAffine3, GceError> {
    let offset = p2 - p1;
    if offset.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    Ok(DAffine3::from_translation(offset))
}

// ============================================================================
// 2D variants
// ============================================================================

// --- GC_MakeMirror2d ---

/// Construct a 2D mirror transformation across a point (central symmetry).
///
/// OCCT: `GC_MakeMirror2d(gp_Pnt2d)`.
pub fn make_mirror2d_point(point: DVec2) -> DAffine2 {
    DAffine2::from_translation(point) * DAffine2::from_scale(DVec2::splat(-1.0)) * DAffine2::from_translation(-point)
}

/// Construct a 2D mirror transformation across a line.
///
/// OCCT: `GC_MakeMirror2d(gp_Ax2d)`.
pub fn make_mirror2d_axis(origin: DVec2, direction: DVec2) -> Result<DAffine2, GceError> {
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < 0.5 {
        return Err(GceError::NullAxis);
    }
    Ok(mirror_2d(origin, dir))
}

// --- GC_MakeRotation2d ---

/// Construct a 2D rotation transformation around a center point.
///
/// OCCT: `GC_MakeRotation2d(gp_Pnt2d, double)`.
pub fn make_rotation2d(center: DVec2, angle_rad: f64) -> Result<DAffine2, GceError> {
    if angle_rad.abs() < TOL_ANG {
        return Err(GceError::NullAngle);
    }
    Ok(rotation_2d(center, angle_rad))
}

// --- GC_MakeScale2d ---

/// Construct a 2D scaling transformation relative to a center point.
///
/// OCCT: `GC_MakeScale2d(gp_Pnt2d, double)`.
pub fn make_scale2d(center: DVec2, scale: f64) -> Result<DAffine2, GceError> {
    if scale < 0.0 {
        return Err(GceError::NegativeRadius);
    }
    if scale.abs() < TOL_CONF {
        return Err(GceError::NullRadius);
    }
    Ok(DAffine2::from_translation(center)
        * DAffine2::from_scale(DVec2::splat(scale))
        * DAffine2::from_translation(-center))
}

// --- GC_MakeTranslation2d ---

/// Construct a 2D translation from a vector.
///
/// OCCT: `GC_MakeTranslation2d(gp_Vec2d)`.
pub fn make_translation2d_vec(offset: DVec2) -> DAffine2 {
    DAffine2::from_translation(offset)
}

/// Construct a 2D translation from two points.
///
/// OCCT: `GC_MakeTranslation2d(gp_Pnt2d, gp_Pnt2d)`.
pub fn make_translation2d_2p(p1: DVec2, p2: DVec2) -> Result<DAffine2, GceError> {
    let offset = p2 - p1;
    if offset.length_squared() < TOL_CONF * TOL_CONF {
        return Err(GceError::ConfusedPoints);
    }
    Ok(DAffine2::from_translation(offset))
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_mirror_plane() {
        let mirror = make_mirror_plane(DVec3::ZERO, DVec3::Z).unwrap();
        // Mirror (0,0,1) across z=0 plane -> (0,0,-1)
        let p = mirror.transform_point3(DVec3::new(0.0, 0.0, 1.0));
        assert!((p - DVec3::new(0.0, 0.0, -1.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_rotation() {
        let rot = make_rotation(DVec3::ZERO, DVec3::Z, std::f64::consts::FRAC_PI_2).unwrap();
        // Rotate (1,0,0) by 90° around Z -> (0,1,0)
        let p = rot.transform_point3(DVec3::new(1.0, 0.0, 0.0));
        assert!((p - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_scale() {
        let s = make_scale(DVec3::ZERO, 2.0).unwrap();
        let p = s.transform_point3(DVec3::new(1.0, 2.0, 3.0));
        assert!((p - DVec3::new(2.0, 4.0, 6.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_translation_vec() {
        let t = make_translation_vec(DVec3::new(1.0, 2.0, 3.0));
        let p = t.transform_point3(DVec3::new(5.0, 7.0, 9.0));
        assert!((p - DVec3::new(6.0, 9.0, 12.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_rotation2d() {
        let rot = make_rotation2d(DVec2::ZERO, std::f64::consts::FRAC_PI_2).unwrap();
        let p = rot.transform_point2(DVec2::new(1.0, 0.0));
        assert!((p - DVec2::new(0.0, 1.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_scale2d() {
        let s = make_scale2d(DVec2::ZERO, 2.0).unwrap();
        let p = s.transform_point2(DVec2::new(1.0, 2.0));
        assert!((p - DVec2::new(2.0, 4.0)).length() < 1e-12);
    }

    #[test]
    fn test_make_mirror2d_axis() {
        let mirror = make_mirror2d_axis(DVec2::ZERO, DVec2::X).unwrap();
        // Mirror (0,1) across X axis -> (0,-1)
        let p = mirror.transform_point2(DVec2::new(0.0, 1.0));
        assert!((p - DVec2::new(0.0, -1.0)).length() < 1e-12);
    }
}
