//! ComputationMethods — cylinder-cylinder coefficient tables — 1:1 translation
//! of OCCT `IntPatch_ImpImpIntersection.cxx` class ComputationMethods
//! (L3949-4062), its stCoeffsValue constructor (L4063-4252), and the
//! CylCylMonotonicity / CylCylComputeParameters methods (L5251-5452).
//!
//! rcad data-model notes:
//!   - `gp_Cylinder` -> `rcad_kernel::geom::CylindricalSurface` (origin, axis,
//!     radius, ref_dir).
//!   - `math_Vector`(3) -> `glam::DVec3`; component (1)/(2)/(3) -> .x/.y/.z.
//!   - OCCT `throw Standard_Failure` -> `panic!` (rcad convention).

use glam::DVec3;
use rcad_kernel::geom::CylindricalSurface;
use rcad_kernel::precision::ANGULAR;

use super::cycy_common::{A_NUL_VALUE, inscribe_point, short_cos_form};

/// OCCT ComputationMethods::stCoeffsValue — stores the equations coefficients.
#[derive(Debug, Clone)]
pub struct StCoeffsValue {
    pub m_vec_a1: DVec3,
    pub m_vec_a2: DVec3,
    pub m_vec_b1: DVec3,
    pub m_vec_b2: DVec3,
    pub m_vec_c1: DVec3,
    pub m_vec_c2: DVec3,
    pub m_vec_d: DVec3,

    pub m_k21: f64, // sinU2
    pub m_k11: f64, // sinU1
    pub m_l21: f64, // cosU2
    pub m_l11: f64, // cosU1
    pub m_m1: f64,  // Free member

    pub m_k22: f64, // sinU2
    pub m_k12: f64, // sinU1
    pub m_l22: f64, // cosU2
    pub m_l12: f64, // cosU1
    pub m_m2: f64,  // Free member

    pub m_k1: f64,
    pub m_l1: f64,
    pub m_k2: f64,
    pub m_l2: f64,

    pub m_fiv1: f64,
    pub m_psiv1: f64,
    pub m_fiv2: f64,
    pub m_psiv2: f64,

    pub m_b: f64,
    pub m_c: f64,
    pub m_fi1: f64,
    pub m_fi2: f64,
}

/// OCCT ComputationMethods::stCoeffsValue::stCoeffsValue(theCyl1, theCyl2)
/// (L4063-4252).
#[allow(unused_assignments)]
pub fn new_coeffs(cyl1: &CylindricalSurface, cyl2: &CylindricalSurface) -> StCoeffsValue {
    let x1 = cyl1.ref_dir.normalize_or_zero();
    let y1 = cyl1.axis.normalize_or_zero().cross(x1).normalize_or_zero();
    let z1 = cyl1.axis.normalize_or_zero();
    let x2 = cyl2.ref_dir.normalize_or_zero();
    let y2 = cyl2.axis.normalize_or_zero().cross(x2).normalize_or_zero();
    let z2 = cyl2.axis.normalize_or_zero();

    let mut c = StCoeffsValue {
        m_vec_a1: -cyl1.radius * x1,
        m_vec_a2: cyl2.radius * x2,
        m_vec_b1: -cyl1.radius * y1,
        m_vec_b2: cyl2.radius * y2,
        m_vec_c1: z1,
        m_vec_c2: -z2,
        m_vec_d: cyl2.origin - cyl1.origin,
        m_k21: 0.0,
        m_k11: 0.0,
        m_l21: 0.0,
        m_l11: 0.0,
        m_m1: 0.0,
        m_k22: 0.0,
        m_k12: 0.0,
        m_l22: 0.0,
        m_l12: 0.0,
        m_m2: 0.0,
        m_k1: 0.0,
        m_l1: 0.0,
        m_k2: 0.0,
        m_l2: 0.0,
        m_fiv1: 0.0,
        m_psiv1: 0.0,
        m_fiv2: 0.0,
        m_psiv2: 0.0,
        m_b: 0.0,
        m_c: 0.0,
        m_fi1: 0.0,
        m_fi2: 0.0,
    };

    // enum CoupleOfEquation
    const COENONE: i32 = 0;
    const COE12: i32 = 1;
    const COE23: i32 = 2;
    const COE13: i32 = 3;
    let mut a_found_couple = COENONE;
    let mut a_det_v1v2 = 0.0;

    let a_delta1 = c.m_vec_c1.x * c.m_vec_c2.y - c.m_vec_c1.y * c.m_vec_c2.x; // 1-2
    let a_delta2 = c.m_vec_c1.y * c.m_vec_c2.z - c.m_vec_c1.z * c.m_vec_c2.y; // 2-3
    let a_delta3 = c.m_vec_c1.x * c.m_vec_c2.z - c.m_vec_c1.z * c.m_vec_c2.x; // 1-3
    let an_abs_d1 = a_delta1.abs(); // 1-2
    let an_abs_d2 = a_delta2.abs(); // 2-3
    let an_abs_d3 = a_delta3.abs(); // 1-3

    if an_abs_d1 >= an_abs_d2 {
        if an_abs_d3 > an_abs_d1 {
            a_found_couple = COE13;
            a_det_v1v2 = a_delta3;
        } else {
            a_found_couple = COE12;
            a_det_v1v2 = a_delta1;
        }
    } else {
        if an_abs_d3 > an_abs_d2 {
            a_found_couple = COE13;
            a_det_v1v2 = a_delta3;
        } else {
            a_found_couple = COE23;
            a_det_v1v2 = a_delta2;
        }
    }

    // If sine of the angle between the axes is too small then the axes are
    // considered parallel (compare with angular tolerance, see
    // AxeOperator::AxeOperator in IntAna_QuadQuadGeo.cxx).
    if a_det_v1v2.abs() < ANGULAR {
        panic!("Error. Exception in divide by zerro (IntCyCyTrim)!!!!");
    }

    match a_found_couple {
        COE12 => {}
        COE23 => {
            // rotate every vector: (1,2,3) -> (2,3,1)
            c.m_vec_a1 = DVec3::new(c.m_vec_a1.y, c.m_vec_a1.z, c.m_vec_a1.x);
            c.m_vec_a2 = DVec3::new(c.m_vec_a2.y, c.m_vec_a2.z, c.m_vec_a2.x);
            c.m_vec_b1 = DVec3::new(c.m_vec_b1.y, c.m_vec_b1.z, c.m_vec_b1.x);
            c.m_vec_b2 = DVec3::new(c.m_vec_b2.y, c.m_vec_b2.z, c.m_vec_b2.x);
            c.m_vec_c1 = DVec3::new(c.m_vec_c1.y, c.m_vec_c1.z, c.m_vec_c1.x);
            c.m_vec_c2 = DVec3::new(c.m_vec_c2.y, c.m_vec_c2.z, c.m_vec_c2.x);
            c.m_vec_d = DVec3::new(c.m_vec_d.y, c.m_vec_d.z, c.m_vec_d.x);
        }
        COE13 => {
            // swap components 2 and 3
            c.m_vec_a1 = DVec3::new(c.m_vec_a1.x, c.m_vec_a1.z, c.m_vec_a1.y);
            c.m_vec_a2 = DVec3::new(c.m_vec_a2.x, c.m_vec_a2.z, c.m_vec_a2.y);
            c.m_vec_b1 = DVec3::new(c.m_vec_b1.x, c.m_vec_b1.z, c.m_vec_b1.y);
            c.m_vec_b2 = DVec3::new(c.m_vec_b2.x, c.m_vec_b2.z, c.m_vec_b2.y);
            c.m_vec_c1 = DVec3::new(c.m_vec_c1.x, c.m_vec_c1.z, c.m_vec_c1.y);
            c.m_vec_c2 = DVec3::new(c.m_vec_c2.x, c.m_vec_c2.z, c.m_vec_c2.y);
            c.m_vec_d = DVec3::new(c.m_vec_d.x, c.m_vec_d.z, c.m_vec_d.y);
        }
        _ => {}
    }

    //------- For V1 (begin)
    // sinU2
    c.m_k21 = (c.m_vec_c2.y * c.m_vec_b2.x - c.m_vec_c2.x * c.m_vec_b2.y) / a_det_v1v2;
    // sinU1
    c.m_k11 = (c.m_vec_c2.y * c.m_vec_b1.x - c.m_vec_c2.x * c.m_vec_b1.y) / a_det_v1v2;
    // cosU2
    c.m_l21 = (c.m_vec_c2.y * c.m_vec_a2.x - c.m_vec_c2.x * c.m_vec_a2.y) / a_det_v1v2;
    // cosU1
    c.m_l11 = (c.m_vec_c2.y * c.m_vec_a1.x - c.m_vec_c2.x * c.m_vec_a1.y) / a_det_v1v2;
    // Free member
    c.m_m1 = (c.m_vec_c2.y * c.m_vec_d.x - c.m_vec_c2.x * c.m_vec_d.y) / a_det_v1v2;
    //------- For V1 (end)

    //------- For V2 (begin)
    // sinU2
    c.m_k22 = (c.m_vec_c1.x * c.m_vec_b2.y - c.m_vec_c1.y * c.m_vec_b2.x) / a_det_v1v2;
    // sinU1
    c.m_k12 = (c.m_vec_c1.x * c.m_vec_b1.y - c.m_vec_c1.y * c.m_vec_b1.x) / a_det_v1v2;
    // cosU2
    c.m_l22 = (c.m_vec_c1.x * c.m_vec_a2.y - c.m_vec_c1.y * c.m_vec_a2.x) / a_det_v1v2;
    // cosU1
    c.m_l12 = (c.m_vec_c1.x * c.m_vec_a1.y - c.m_vec_c1.y * c.m_vec_a1.x) / a_det_v1v2;
    // Free member
    c.m_m2 = (c.m_vec_c1.x * c.m_vec_d.y - c.m_vec_c1.y * c.m_vec_d.x) / a_det_v1v2;
    //------- For V2 (end)

    let (m_k1, m_fiv1) = short_cos_form(c.m_l11, c.m_k11);
    c.m_k1 = m_k1;
    c.m_fiv1 = m_fiv1;
    let (m_l1, m_psiv1) = short_cos_form(c.m_l21, c.m_k21);
    c.m_l1 = m_l1;
    c.m_psiv1 = m_psiv1;
    let (m_k2, m_fiv2) = short_cos_form(c.m_l12, c.m_k12);
    c.m_k2 = m_k2;
    c.m_fiv2 = m_fiv2;
    let (m_l2, m_psiv2) = short_cos_form(c.m_l22, c.m_k22);
    c.m_l2 = m_l2;
    c.m_psiv2 = m_psiv2;

    let a_a1 = c.m_vec_c1.z * c.m_k21 + c.m_vec_c2.z * c.m_k22 - c.m_vec_b2.z; // sinU2
    let a_a2 = c.m_vec_c1.z * c.m_l21 + c.m_vec_c2.z * c.m_l22 - c.m_vec_a2.z; // cosU2
    let a_b1 = c.m_vec_b1.z - c.m_vec_c1.z * c.m_k11 - c.m_vec_c2.z * c.m_k12; // sinU1
    let a_b2 = c.m_vec_a1.z - c.m_vec_c1.z * c.m_l11 - c.m_vec_c2.z * c.m_l12; // cosU1

    c.m_c = c.m_vec_d.z - c.m_vec_c1.z * c.m_m1 - c.m_vec_c2.z * c.m_m2; // Free

    let mut a_a = 0.0;
    let (m_b, m_fi1) = short_cos_form(a_b2, a_b1);
    c.m_b = m_b;
    c.m_fi1 = m_fi1;
    let (a_a_out, m_fi2) = short_cos_form(a_a2, a_a1);
    a_a = a_a_out;
    c.m_fi2 = m_fi2;

    c.m_b /= a_a;
    c.m_c /= a_a;

    c
}

/// OCCT ComputationMethods::CylCylMonotonicity (L5251-5322) — determines if
/// U2(U1) function is increasing.
pub fn cyl_cyl_monotonicity(
    u1par: f64,
    wl_index: usize,
    coeffs: &StCoeffsValue,
    period: f64,
    is_increasing: &mut bool,
) -> bool {
    // For "+/-" sign. If isPlus == TRUE, "+" is chosen, otherwise "-".
    let is_plus = match wl_index {
        0 => true,
        1 => false,
        _ => return false,
    };

    let mut u1_temp = u1par - coeffs.m_fi1;
    inscribe_point(0.0, period, &mut u1_temp, 0.0, period, false);

    *is_increasing = true;

    if ((std::f64::consts::PI - u1_temp) < f64::MIN_POSITIVE) && (u1_temp < period) {
        *is_increasing = false;
    }

    if coeffs.m_b < 0.0 {
        *is_increasing = !*is_increasing;
    }

    if !is_plus {
        *is_increasing = !*is_increasing;
    }

    true
}

/// OCCT ComputationMethods::CylCylComputeParameters (L5323-5405) — computes U2
/// (U-parameter of the 2nd cylinder) and, if theDelta != 0, estimates the
/// tolerance of U2-computing.
pub fn cyl_cyl_compute_parameters_u2(
    u1par: f64,
    wl_index: i32,
    coeffs: &StCoeffsValue,
    u2: &mut f64,
    mut delta: Option<&mut f64>,
) -> bool {
    // This formula is got from some experience and can be changed.
    let a_tol0 = (10.0 * f64::EPSILON * coeffs.m_b).min(A_NUL_VALUE);
    let a_tol = 1.0 - a_tol0;

    if wl_index < 0 || wl_index > 1 {
        return false;
    }

    let a_sign = if wl_index != 0 { -1.0 } else { 1.0 };

    let mut an_arg = (u1par - coeffs.m_fi1).cos();
    an_arg = coeffs.m_b * an_arg + coeffs.m_c;

    if an_arg >= a_tol {
        if let Some(d) = delta.as_deref_mut() {
            *d = 0.0;
        }
        an_arg = 1.0;
    } else if an_arg <= -a_tol {
        if let Some(d) = delta.as_deref_mut() {
            *d = 0.0;
        }
        an_arg = -1.0;
    } else if let Some(d) = delta.as_deref_mut() {
        let a_delta = (1.0 - an_arg).min(1.0 + an_arg);
        if (a_delta * a_delta < f64::MIN_POSITIVE) || (a_delta >= 2.0) {
            panic!("IntPatch_ImpImpIntersection_4.gxx, CylCylComputeParameters()");
        }
        *d = a_tol0 / (a_delta * (2.0 - a_delta)).sqrt();
    }

    *u2 = an_arg.acos();
    *u2 = coeffs.m_fi2 + a_sign * *u2;

    true
}

/// OCCT ComputationMethods::CylCylComputeParameters (L5406-5425) — computes V1
/// and V2 (V-parameters of the 1st and 2nd cylinder respectively).
pub fn cyl_cyl_compute_parameters_v(
    u1: f64,
    u2: f64,
    coeffs: &StCoeffsValue,
    v1: &mut f64,
    v2: &mut f64,
) -> bool {
    *v1 = coeffs.m_k21 * u2.sin() + coeffs.m_k11 * u1.sin() + coeffs.m_l21 * u2.cos()
        + coeffs.m_l11 * u1.cos() + coeffs.m_m1;

    *v2 = coeffs.m_k22 * u2.sin() + coeffs.m_k12 * u1.sin() + coeffs.m_l22 * u2.cos()
        + coeffs.m_l12 * u1.cos() + coeffs.m_m2;

    true
}

/// OCCT ComputationMethods::CylCylComputeParameters (L5426-5452) — computes U2,
/// V1 and V2.
pub fn cyl_cyl_compute_parameters(
    u1par: f64,
    wl_index: i32,
    coeffs: &StCoeffsValue,
    u2: &mut f64,
    v1: &mut f64,
    v2: &mut f64,
) -> bool {
    if !cyl_cyl_compute_parameters_u2(u1par, wl_index, coeffs, u2, None) {
        return false;
    }
    if !cyl_cyl_compute_parameters_v(u1par, *u2, coeffs, v1, v2) {
        return false;
    }
    true
}
