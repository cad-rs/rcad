//! gp_Vec / gp_Dir helpers consumed by the Extrema translations.
//!
//! These are the exact bodies of the OCCT inline/compiled members the Extrema
//! package calls (`gp_Dir::Angle`, `gp_Dir::AngleWithRef`, `gp_Vec::Angle`,
//! `gp_Vec::AngleWithRef`, `gp_Dir::IsParallel`, `gp_Dir::IsNormal`,
//! `gp_Vec::IsNormal`).  They live here rather than in `math/gp.rs` because the
//! E3-Q file domain is restricted to `base/extrema*.rs` (interface request:
//! fold them back into the kernel gp module once the domain opens).

use glam::DVec3;

/// OCCT gp_Dir::Angle (gp_Dir.cxx L27-50) — the angle in [0, PI] computed with
/// the arccos/arcsin switch that keeps precision near 0 and PI.
pub(crate) fn dir_angle(a: DVec3, other: DVec3) -> f64 {
    let cosinus = a.dot(other);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = a.cross(other).length();
        if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT gp_Dir::AngleWithRef (gp_Dir.cxx L55-81) — the signed angle, the sign
/// taken from the reference direction.
pub(crate) fn dir_angle_with_ref(a: DVec3, other: DVec3, v_ref: DVec3) -> f64 {
    let xyz = a.cross(other);
    let cosinus = a.dot(other);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(v_ref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488-494) — the vectors are normalized first
/// (`gp_Dir(coord)`), so the result is `gp_Dir::Angle` of the directions.
pub(crate) fn vec_angle(a: DVec3, other: DVec3) -> f64 {
    dir_angle(a.normalize_or_zero(), other.normalize_or_zero())
}

/// OCCT gp_Vec::AngleWithRef (gp_Vec.hxx L498-505).
pub(crate) fn vec_angle_with_ref(a: DVec3, other: DVec3, v_ref: DVec3) -> f64 {
    dir_angle_with_ref(a.normalize_or_zero(), other, v_ref)
}

/// OCCT gp_Dir::IsParallel (gp_Dir.hxx L191-195).
pub(crate) fn dir_is_parallel(a: DVec3, other: DVec3, ang_tol: f64) -> bool {
    let an_ang = dir_angle(a, other);
    an_ang <= ang_tol || std::f64::consts::PI - an_ang <= ang_tol
}

/// OCCT gp_Vec::IsNormal (gp_Vec.hxx L136-141) applied to unit inputs.
pub(crate) fn dir_is_normal(a: DVec3, other: DVec3, ang_tol: f64) -> bool {
    let an_ang = (std::f64::consts::FRAC_PI_2 - vec_angle(a, other)).abs();
    an_ang <= ang_tol
}
