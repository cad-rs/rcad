use glam::DVec3;

/// Absolute tolerance for point coincidence.
///
/// Matches `rcad_kernel::tolerance::CONFUSION` = `Precision::Confusion()` in OCCT.
/// Two points are considered coincident when their distance is below this value.
pub const TOLERANCE_ABS: f64 = 1e-7;

/// Angular tolerance for parallel/perpendicular checks (radians, as cross-product magnitude).
///
/// This is intentionally **looser** than `rcad_kernel::tolerance::ANGULAR` (1e-12):
/// the algorithms layer needs to tolerate slightly imperfect parallelism that
/// arises from floating-point accumulation during intersection computation.
/// Used in [`vectors_parallel`] as `cross(a,b).length_squared() < TOLERANCE_ANG²`.
pub const TOLERANCE_ANG: f64 = 1e-9;

/// Tolerance squared — avoids `sqrt` in distance checks.
pub const TOLERANCE_ABS_SQ: f64 = TOLERANCE_ABS * TOLERANCE_ABS;

#[inline]
pub fn points_coincide(a: DVec3, b: DVec3) -> bool {
    (a - b).length_squared() < TOLERANCE_ABS_SQ
}

#[inline]
pub fn is_zero_vec(v: DVec3) -> bool {
    v.length_squared() < TOLERANCE_ABS_SQ
}

/// Returns true if two unit vectors are parallel (or anti-parallel).
#[inline]
pub fn vectors_parallel(a: DVec3, b: DVec3) -> bool {
    a.cross(b).length_squared() < TOLERANCE_ANG * TOLERANCE_ANG
}

#[inline]
pub fn params_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < TOLERANCE_ABS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coincident_points() {
        let a = DVec3::new(1.0, 2.0, 3.0);
        let b = DVec3::new(1.0, 2.0, 3.0 + 1e-8);
        assert!(points_coincide(a, b));
    }

    #[test]
    fn non_coincident_points() {
        let a = DVec3::ZERO;
        let b = DVec3::new(1e-6, 0.0, 0.0);
        assert!(!points_coincide(a, b));
    }

    #[test]
    fn parallel_vectors() {
        assert!(vectors_parallel(DVec3::X, DVec3::X));
        assert!(vectors_parallel(DVec3::X, -DVec3::X));
        assert!(!vectors_parallel(DVec3::X, DVec3::Y));
    }
}
