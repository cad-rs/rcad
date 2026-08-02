//! OCCT IntPatch_PointLine — curvature radius of the surface/surface
//! intersection line (IntPatch_PointLine.cxx L50-141).
//!
//! 1:1 Rust translation of IntPatch_PointLine::CurvatureRadiusOfIntersLine.

use glam::DVec3;
use rcad_kernel::geom::{Surface3, SurfaceEval};

/// OCCT IntPatch_PointLine::CurvatureRadiusOfIntersLine (L50-141).
///
/// Returns the radius of curvature of the intersection line of theS1/theS2 at
/// the point given by its UV parameters on both surfaces.  Returns a negative
/// value when the computation is not possible.
pub fn curvature_radius_of_inters_line(
    s1: &Surface3,
    s2: &Surface3,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
) -> f64 {
    // constexpr double aSmallValue   = 1.0 / Precision::Infinite();
    let a_small_value = 1.0 / f64::INFINITY;
    let a_sq_small_value = a_small_value * a_small_value;

    let (a_pt, a_du1, a_dv1, a_duu1, a_duv1, a_dvv1) = s1.derivatives2(u1, v1);
    let (_, a_du2, a_dv2, a_duu2, a_duv2, a_dvv2) = s2.derivatives2(u2, v2);

    let a_n1 = a_du1.cross(a_dv1);
    let a_n2 = a_du2.cross(a_dv2);
    // Tangent vector to the intersection curve
    let a_c_tan = a_n1.cross(a_n2);
    let a_sq_magn_f_der = a_c_tan.length_squared();

    if a_sq_magn_f_der < 1.0e-8 {
        return -1.0;
    }

    let mut a_du_s1 = 0.0f64;
    let mut a_dv_s1 = 0.0f64;
    let mut a_du_s2 = 0.0f64;
    let mut a_dv_s2 = 1.0f64;

    // This algorithm is described in NonSingularProcessing() function
    // in ApproxInt_ImpPrmSvSurfaces.gxx file.
    let mut a_sq_n_magn = a_n1.length_squared();
    let a_tg_u = a_c_tan.cross(a_du1);
    let a_tg_v = a_c_tan.cross(a_dv1);
    let a_delta_u = a_tg_v.length_squared() / a_sq_n_magn;
    let a_delta_v = a_tg_u.length_squared() / a_sq_n_magn;

    a_du_s1 = a_delta_u.sqrt().copysign(a_tg_v.dot(a_n1));
    a_dv_s1 = -a_delta_v.sqrt().copysign(a_tg_u.dot(a_n1));

    a_sq_n_magn = a_n2.length_squared();
    let a_tg_u2 = a_c_tan.cross(a_du2);
    let a_tg_v2 = a_c_tan.cross(a_dv2);
    let a_delta_u2 = a_tg_v2.length_squared() / a_sq_n_magn;
    let a_delta_v2 = a_tg_u2.length_squared() / a_sq_n_magn;

    a_du_s2 = a_delta_u2.sqrt().copysign(a_tg_v2.dot(a_n2));
    a_dv_s2 = -a_delta_v2.sqrt().copysign(a_tg_u2.dot(a_n2));

    // According to "Marching along surface/surface intersection curves
    // with an adaptive step length" by Tz.E.Stoyagov, the system:
    //   { A*a + B*b = F1
    //   { B*a + C*b = F2
    // where a and b should be found.  After that, the 2nd derivative of the
    // intersection curve is r''(t) = a*aN1 + b*aN2.
    let a_a = a_n1.dot(a_n1);
    let a_b = a_n1.dot(a_n2);
    let a_c = a_n2.dot(a_n2);
    let a_det_syst = a_b * a_b - a_a * a_c;

    if a_det_syst.abs() < a_small_value {
        // Undetermined system solution
        return -1.0;
    }

    let a_f1 = a_du_s1 * a_du_s1 * a_duu1.dot(a_n1) + 2.0 * a_du_s1 * a_dv_s1 * a_duv1.dot(a_n1)
        + a_dv_s1 * a_dv_s1 * a_dvv1.dot(a_n1);
    let a_f2 = a_du_s2 * a_du_s2 * a_duu2.dot(a_n2) + 2.0 * a_du_s2 * a_dv_s2 * a_duv2.dot(a_n2)
        + a_dv_s2 * a_dv_s2 * a_dvv2.dot(a_n2);

    // Principal normal to the intersection curve
    let a_c_norm = DVec3::ZERO
        + (a_f1 * a_c - a_f2 * a_b) / a_det_syst * a_n1
        + (a_a * a_f2 - a_f1 * a_b) / a_det_syst * a_n2;
    // CrossSquareMagnitude(aCNorm, aCTan) = |aCNorm x aCTan|^2
    let a_cross = a_c_norm.cross(a_c_tan);
    let a_sq_magn_s_der = a_cross.length_squared();

    if a_sq_magn_s_der < a_sq_small_value {
        // Intersection curve has null curvature in observed point
        return f64::INFINITY;
    }

    // square of curvature radius
    let a_fact_sq_rad = a_sq_magn_f_der * a_sq_magn_f_der * a_sq_magn_f_der / a_sq_magn_s_der;

    a_fact_sq_rad.sqrt()
}
