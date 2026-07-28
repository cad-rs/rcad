//! Analytic intersection of a cylinder and a cone.
//!
//! # Case classification
//!
//! ## Coaxial (axes coincide)
//!
//! When the cylinder axis and cone axis are the same line, the intersection is
//! one or two **circles** at the height(s) where the cone radius equals the
//! cylinder radius (one per nappe of the infinite double cone):
//!
//! ```text
//! r_cone(h) = (h 鈭?h_apex) 路 tan(尾)  where h_apex is the height of the apex
//! r_cone(h) = r_cyl  鈫? h = h_apex 卤 r_cyl / tan(尾)
//! ```
//!
//! Returns [`CoaxialCircle`](CylinderConeResult::CoaxialCircle) if only one
//! nappe produces a real circle, or
//! [`CoaxialTwoCircles`](CylinderConeResult::CoaxialTwoCircles) when both
//! nappes intersect the cylinder.  Returns
//! [`NoIntersection`](CylinderConeResult::NoIntersection) when the apex is
//! on the positive-axis side (the circle is real).
//!
//! ## Parallel axes (non-coaxial)
//!
//! When the axes are parallel but distinct (cross-product 鈮?0), the
//! cylinder-cone intersection is a quartic curve in general.  We perform a
//! radial distance test to detect obvious non-intersections and fall back to
//! marching otherwise.
//!
//! ## General / skew axes
//!
//! For all other configurations we return [`General`](CylinderConeResult::General)
//! so the caller can fall back to numeric marching.

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Circle3, ConicalSurface, CylindricalSurface, any_perpendicular};
use std::f64::consts::TAU;

use super::pcurve_derive::refine_polyline;


// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Result type
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Analytic result of cylinder 脳 cone intersection.
#[derive(Debug, Clone)]
pub enum CylinderConeResult {
    /// The surfaces do not intersect.
    NoIntersection,
    /// Coaxial configuration: single intersection circle (one nappe).
    CoaxialCircle(Circle3),
    /// Coaxial configuration: two intersection circles (one per nappe
    /// of the double cone).  OCCT IntAna_QuadQuadGeo returns both.
    CoaxialTwoCircles(Circle3, Circle3),
    /// Parallel-offset axes (parallel but not coaxial): intersection is one or
    /// two polylines (branches) on the cylinder surface.  Each branch runs from
    /// the near-side tangent point (h_min) to the far-side tangent point (h_max).
    ///
    /// Two branches occur when the intersection crosses both the near and far
    /// sides of the cylinder; one branch (or a single point) occurs at tangent
    /// transitions.
    ParallelOffsetPolyline(Vec<Vec<DVec3>>),
    /// General case (skew axes or oblique angle not handled analytically).
    /// The caller should fall back to numeric marching.
    General,
    /// Skew axes (non-parallel, non-coaxial): analytic quartic solution.
    ///
    /// The cylinder and cone intersect in a quartic space curve.  For each
    /// cylinder azimuth u 鈭?[0, 2蟺) the cone equation reduces to a quadratic
    /// in the cylinder height v, solved analytically.  Two branches (卤 sqrt)
    /// are returned as polylines.
    ///
    /// Used when the axes are skew (neither parallel nor coincident), avoiding
    /// expensive numeric marching on a dense grid.
    SkewQuartic(Vec<Vec<DVec3>>),
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Main function
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the analytic intersection of `cyl` and `cone`.
pub fn intersect_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
) -> CylinderConeResult {
    let a_cyl = cyl.axis.normalize();
    let a_cone = cone.axis.normalize();

    let cross = a_cyl.cross(a_cone);
    let sin_angle = cross.length(); // |sin 胃| between the two axes

    // 鈹€鈹€ Parallel axes (including coaxial) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if sin_angle < rcad_kernel::ANGULAR {
        return intersect_parallel_cylinder_cone(cyl, cone, a_cyl, a_cone);
    }

    // 鈹€鈹€ Skew axes 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Try analytic quartic solver first; fall back to numeric marching if it
    // produces no result (e.g. surfaces barely graze each other).

    let skew_result = intersect_skew_cylinder_cone(cyl, cone);
    if !skew_result.is_empty() {
        return CylinderConeResult::SkewQuartic(skew_result);
    }

    CylinderConeResult::General
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Parallel (and coaxial) axes
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn intersect_parallel_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
    a_cyl: DVec3,
    a_cone: DVec3,
) -> CylinderConeResult {
    let apex = cone.apex_point();
    // Make sure a_cone points the same direction as a_cyl for height arithmetic.
    // (The cross product is ~0 so they are parallel; they may anti-parallel.)
    let _a_cone = if a_cyl.dot(a_cone) >= 0.0 {
        a_cone
    } else {
        -a_cone
    };

    let r_cyl = cyl.radius;
    let tan_beta = cone.half_angle_rad.tan();

    // Perpendicular distance between the two axes.
    let delta = apex - cyl.origin;
    let delta_perp = delta - a_cyl * delta.dot(a_cyl);
    let d_perp = delta_perp.length();

    // 鈹€鈹€ Coaxial 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if d_perp < rcad_kernel::rcad_kernel::CONFUSION {
        let h_apex = (apex - cyl.origin).dot(a_cyl);

        if tan_beta.abs() < 1e-14 {
            // Degenerate cone (half_angle = 0), no lateral surface.
            return CylinderConeResult::NoIntersection;
        }

        // OCCT IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cylinder), coaxial case:
        // the intersection is one or two circles at the height where the cone
        // radius equals the cylinder radius.  For a double-napped infinite cone
        // there are two solutions: one on each nappe.
        let h_offset = r_cyl / tan_beta;
        let h_upper = h_apex + h_offset; // circle on the upper nappe (h > h_apex)
        let h_lower = h_apex - h_offset; // circle on the lower nappe (h < h_apex)

        let mut first: Option<Circle3> = None;
        let mut second: Option<Circle3> = None;

        if h_upper > h_apex + rcad_kernel::rcad_kernel::CONFUSION {
            first = Some(Circle3::new(cyl.origin + a_cyl * h_upper, a_cyl, r_cyl));
        }
        if h_lower < h_apex - rcad_kernel::rcad_kernel::CONFUSION {
            second = Some(Circle3::new(cyl.origin + a_cyl * h_lower, a_cyl, r_cyl));
        }

        return match (first, second) {
            (Some(c1), Some(c2)) => CylinderConeResult::CoaxialTwoCircles(c1, c2),
            (Some(c1), None) => CylinderConeResult::CoaxialCircle(c1),
            (None, Some(c2)) => CylinderConeResult::CoaxialCircle(c2),
            (None, None) => CylinderConeResult::NoIntersection,
        };
    }

    // 鈹€鈹€ Parallel but offset axes 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // At each axial height h the cone has radius r_cone(h) = |h - h_apex|路tan(尾).
    // The cylinder cross-section is a circle of radius r_cyl centred at a
    // perpendicular offset d_perp from the cone axis.  Intersection reduces to
    // a circle-circle test in the plane perpendicular to the axes, which is
    // solved analytically via axial sweep.
    if d_perp < rcad_kernel::rcad_kernel::CONFUSION {
        // Should not reach here (coaxial case above), but guard against
        // numerical near-zero.
        return CylinderConeResult::NoIntersection;
    }
    return intersect_parallel_offset_cylinder_cone(cyl, cone, a_cyl, d_perp, delta_perp);
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Parallel-offset axial sweep
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the intersection of a cylinder and a cone with parallel but
/// non-coaxial axes via axial sweep.
///
/// At each height `h` along the common axis direction `a`, the cone
/// cross-section is a circle of radius `r_cone(h) = |h - h_apex|路tan(尾)`
/// centred on the cone axis, and the cylinder cross-section is a circle of
/// radius `r_cyl` centred on the cylinder axis.  The two axes are separated
/// by perpendicular distance `d_perp` (vector `delta_perp` from cylinder axis
/// to cone apex in the perpendicular plane).
///
/// The intersection exists when the two circles overlap:
///   `|r_cyl - r_cone(h)| 鈮?d_perp 鈮?r_cyl + r_cone(h)`
///
/// Both the upper nappe (`h > h_apex`) and lower nappe (`h < h_apex`) are
/// swept; each nappe may produce 0, 1, or 2 branches.
fn intersect_parallel_offset_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
    a: DVec3,
    d_perp: f64,
    delta_perp: DVec3,
) -> CylinderConeResult {
    let r_cyl = cyl.radius;
    let tan_beta = cone.half_angle_rad.tan();
    if tan_beta.abs() < 1e-14 {
        return CylinderConeResult::NoIntersection;
    }

    let apex = cone.apex_point();
    let h_apex = (apex - cyl.origin).dot(a);

    // Geometry constants for the sweep.
    let abs_diff = (d_perp - r_cyl).abs();
    let d_sum = d_perp + r_cyl;

    // Perpendicular basis: u points from cylinder axis toward cone axis,
    // v = a 脳 u completes the right-handed frame.
    let u = delta_perp / d_perp;
    let v = a.cross(u).normalize_or_zero();

    const N_SAMPLES: usize = 128;
    let mut branches: Vec<Vec<DVec3>> = Vec::new();

    // Sweep both nappes (upper: sign = +1, lower: sign = -1).
    for &sign in &[1.0, -1.0] {
        let h_min = h_apex + sign * abs_diff / tan_beta;
        let h_max = h_apex + sign * d_sum / tan_beta;

        let (h_lo, h_hi) = if sign > 0.0 {
            (h_min, h_max)
        } else {
            (h_max, h_min)
        };

        if h_hi - h_lo <= rcad_kernel::rcad_kernel::CONFUSION * 10.0 {
            continue; // degenerately narrow intersection band
        }

        let mut branch_plus: Vec<DVec3> = Vec::with_capacity(N_SAMPLES + 1);
        let mut branch_minus: Vec<DVec3> = Vec::with_capacity(N_SAMPLES + 1);

        for i in 0..=N_SAMPLES {
            let t = i as f64 / N_SAMPLES as f64;
            let h = h_lo + (h_hi - h_lo) * t;
            let r_cone = (h - h_apex).abs() * tan_beta;

            let center = cyl.origin + a * h;

            // cos(蠁) where 蠁 is the angle on the cylinder cross-section circle
            // from the u-direction.  Derived from |P(h,蠁) - P_cone(h)|虏 = r_cone虏.
            let numer = r_cyl * r_cyl + d_perp * d_perp - r_cone * r_cone;
            let denom = 2.0 * r_cyl * d_perp;
            let cos_phi = (numer / denom).clamp(-1.0, 1.0);

            if cos_phi.abs() >= 1.0 - rcad_kernel::ANGULAR {
                // Tangent: both branches meet at the same point.
                let pt = center + r_cyl * u * cos_phi;
                branch_plus.push(pt);
                branch_minus.push(pt);
            } else {
                let phi = cos_phi.acos();
                let sin_phi = phi.sin();
                let pt_plus = center + r_cyl * (u * cos_phi + v * sin_phi);
                let pt_minus = center + r_cyl * (u * cos_phi - v * sin_phi);
                branch_plus.push(pt_plus);
                branch_minus.push(pt_minus);
            }
        }

        // Register non-degenerate branches.
        // Deduplicate: skip trailing points that are nearly identical to the first
        // point to avoid zero-length edge chains.
        let dedup = |pts: &mut Vec<DVec3>| {
            while pts.len() >= 3 {
                let n = pts.len();
                if (pts[n - 1] - pts[0]).length_squared() < 1e-24 {
                    pts.pop();
                } else {
                    break;
                }
            }
        };

        if branch_plus.len() >= 2 {
            dedup(&mut branch_plus);
            if branch_plus.len() >= 2 {
                branches.push(branch_plus);
            }
        }
        if branch_minus.len() >= 2 {
            dedup(&mut branch_minus);
            if branch_minus.len() >= 2 {
                branches.push(branch_minus);
            }
        }
    }

    if branches.is_empty() {
        CylinderConeResult::NoIntersection
    } else {
        CylinderConeResult::ParallelOffsetPolyline(branches)
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Skew-axis analytic solver
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the intersection of a cylinder and cone with skew (non-parallel)
/// axes using an analytic quartic solver.
///
/// # Theory
///
/// Cylinder parametrization (u = azimuth [0, 2蟺), v = height along axis):
///
/// ```text
/// P(u,v) = O_cyl + v路a_cyl + r_cyl路(cos(u)路x_cyl + sin(u)路y_cyl)
/// ```
///
/// The cone surface satisfies:
///
/// ```text
/// |P - O_cone|虏路cos虏(伪) = ((P - O_cone)路a_cone)虏
/// ```
///
/// Substituting P(u,v) gives F(v) = 0 where
///
/// ```text
/// a_v路v虏 + b_v(u)路v + c_v(u) = 0
///
/// a_v   = (a_cyl路a_cone)虏 - cos虏(伪)                          (constant)
/// b_v(u) = 2路(D0路a_cone)路(a_cyl路a_cone) - 2路cos虏(伪)路(D0路a_cyl)
/// c_v(u) = (D0路a_cone)虏 - cos虏(伪)路|D0|虏
/// D0(u)  = O_cyl - O_cone + r_cyl路(cos(u)路x_cyl + sin(u)路y_cyl)
/// ```
///
/// For each u, we solve the quadratic for v, giving two branches (卤 sqrt).
fn intersect_skew_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
) -> Vec<Vec<DVec3>> {
    let a_cyl = cyl.axis.normalize();
    let a_cone = cone.axis_dir();
    let o_cyl = cyl.origin;
    let o_cone = cone.apex_point();
    let r_cyl = cyl.radius;
    let cos_alpha = cone.half_angle_rad.cos();
    let cos2 = cos_alpha * cos_alpha; // cos虏(伪)

    // Perpendicular basis for cylinder (must match CylindricalSurface::point_at).
    let x_cyl = any_perpendicular(a_cyl);
    let y_cyl = a_cyl.cross(x_cyl).normalize();

    // Constant coefficient a_v = (a_cyl路a_cone)虏 - cos虏(伪).
    let a_dot = a_cyl.dot(a_cone);
    let a_v = a_dot * a_dot - cos2;

    // Pre-computed constants for b_v and c_v.
    let a_dot2 = 2.0 * a_dot; // 2路(a_cyl路a_cone)
    let two_cos2 = 2.0 * cos2; // 2路cos虏(伪)

    let delta_o = o_cyl - o_cone; // O_cyl - O_cone

    const N_SAMPLES: usize = 128;
    const CHORD_TOL: f64 = crate::bop::int_tools::CHORD_TOLERANCE;
    const REFINE_DEPTH: usize = crate::bop::int_tools::CHORD_REFINE_DEPTH;
    let mut branch_plus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);
    let mut branch_minus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cos_u, sin_u) = (u.cos(), u.sin());

        // D0(u) = (O_cyl - O_cone) + r_cyl路(cos(u)路x_cyl + sin(u)路y_cyl)
        let d0 = delta_o + r_cyl * (cos_u * x_cyl + sin_u * y_cyl);

        let d0_a = d0.dot(a_cone); // D0路a_cone
        let d0_cyl = d0.dot(a_cyl); // D0路a_cyl
        let d0_sq = d0.length_squared(); // |D0|虏

        // b_v(u) = 2路(D0路a_cone)路(a_cyl路a_cone) - 2路cos虏(伪)路(D0路a_cyl)
        let b_v = a_dot2 * d0_a - two_cos2 * d0_cyl;

        // c_v(u) = (D0路a_cone)虏 - cos虏(伪)路|D0|虏
        let c_v = d0_a * d0_a - cos2 * d0_sq;

        if a_v.abs() > 1e-12 {
            // Quadratic: a_v路v虏 + b_v路v + c_v = 0
            let disc = b_v * b_v - 4.0 * a_v * c_v;

            if disc < 0.0 {
                // No intersection at this u azimuth.
                continue;
            }

            let sqrt_disc = disc.sqrt();
            let two_a_v = 2.0 * a_v;

            let v_plus = (-b_v + sqrt_disc) / two_a_v;
            let v_minus = (-b_v - sqrt_disc) / two_a_v;

            // Only push valid (finite) points.
            if v_plus.is_finite() {
                let p = cyl.point_at(u, v_plus);
                if p.is_finite() {
                    branch_plus.push((u, p));
                }
            }
            if v_minus.is_finite() {
                let p = cyl.point_at(u, v_minus);
                if p.is_finite() {
                    branch_minus.push((u, p));
                }
            }
        } else if b_v.abs() > 1e-12 {
            // Near-degenerate a_v 鈮?0: solve linear b_v路v + c_v = 0
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cyl.point_at(u, v);
                if p.is_finite() {
                    branch_plus.push((u, p));
                    branch_minus.push((u, p));
                }
            }
        }
    }

    // Adaptive refinement: subdivide segments where chord error exceeds tolerance
    let eval_plus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d0 = delta_o + r_cyl * (cos_u * x_cyl + sin_u * y_cyl);
        let d0_a = d0.dot(a_cone);
        let d0_cyl = d0.dot(a_cyl);
        let d0_sq = d0.length_squared();
        let b_v = a_dot2 * d0_a - two_cos2 * d0_cyl;
        let c_v = d0_a * d0_a - cos2 * d0_sq;
        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                return None;
            }
            let v = (-b_v + disc.sqrt()) / (2.0 * a_v);
            if v.is_finite() {
                let p = cyl.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else if b_v.abs() > 1e-12 {
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cyl.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else {
            None
        }
    };
    let eval_minus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d0 = delta_o + r_cyl * (cos_u * x_cyl + sin_u * y_cyl);
        let d0_a = d0.dot(a_cone);
        let d0_cyl = d0.dot(a_cyl);
        let d0_sq = d0.length_squared();
        let b_v = a_dot2 * d0_a - two_cos2 * d0_cyl;
        let c_v = d0_a * d0_a - cos2 * d0_sq;
        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                return None;
            }
            let v = (-b_v - disc.sqrt()) / (2.0 * a_v);
            if v.is_finite() {
                let p = cyl.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else if b_v.abs() > 1e-12 {
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cyl.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else {
            None
        }
    };

    let (mut branch_plus, mut branch_minus): (Vec<DVec3>, Vec<DVec3>) = (
        refine_polyline(&branch_plus, eval_plus, CHORD_TOL, REFINE_DEPTH)
            .into_iter()
            .map(|(_, p)| p)
            .collect(),
        refine_polyline(&branch_minus, eval_minus, CHORD_TOL, REFINE_DEPTH)
            .into_iter()
            .map(|(_, p)| p)
            .collect(),
    );

    // Dedup: remove trailing points that nearly duplicate the first point
    // (closed curve degeneracy).
    let dedup = |pts: &mut Vec<DVec3>| {
        while pts.len() >= 3 {
            let n = pts.len();
            if (pts[n - 1] - pts[0]).length_squared() < 1e-24 {
                pts.pop();
            } else {
                break;
            }
        }
    };

    let mut branches = Vec::new();
    if branch_plus.len() >= 2 {
        dedup(&mut branch_plus);
        if branch_plus.len() >= 2 {
            branches.push(branch_plus);
        }
    }
    if branch_minus.len() >= 2 {
        // Check the minus branch is distinct from the plus branch by comparing
        // their first points. If they're nearly the same, the branches collapsed
        // (tangent intersection) 鈥?keep only one.
        let is_distinct = branches.is_empty()
            || (branch_minus[0] - branches[0][0]).length_squared() > 1e-24;
        if is_distinct {
            dedup(&mut branch_minus);
            if branch_minus.len() >= 2 {
                branches.push(branch_minus);
            }
        }
    }

    branches
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Tests
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
