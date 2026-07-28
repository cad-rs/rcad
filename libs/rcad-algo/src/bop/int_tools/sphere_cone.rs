//! Analytic intersection of a sphere and a cone.
//!
//! # Case classification
//!
//! ## Sphere center on cone axis (axis-aligned case)
//!
//! When the sphere center lies on the cone axis, the intersection consists of
//! circles at axial heights where the sphere's cross-section radius equals
//! the cone's radius at that height. Solve:
//!
//!   r_sphere(z) = sqrt(R虏 - (z - z_c)虏)
//!   r_cone(z) = r_ref + (z - z_ref) * tan(half_angle)
//!
//! The intersection points satisfy r_sphere(z) = r_cone(z), which leads to
//! a quartic equation. We use numerical root-finding to locate solutions.
//!
//! ## General case
//!
//! For all other configurations the intersection is a space curve of degree <= 4.
//! We return `General` so the caller falls back to numeric marching.

use glam::DVec3;
use rcad_kernel::any_perpendicular;
use rcad_kernel::geom::{Circle3, ConicalSurface, SphericalSurface};



// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Result type
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Analytic result of sphere x cone intersection.
#[derive(Debug, Clone)]
pub enum SphereConeResult {
    /// The sphere and cone do not intersect.
    NoIntersection,
    /// Sphere center is on cone axis; intersection is a single circle.
    SingleCircle(Circle3),
    /// Sphere center is on cone axis; intersection consists of two circles.
    TwoCircles(Circle3, Circle3),
    /// The intersection is a tangent point.
    TangentPoint(DVec3),
    /// Off-axis intersection: one or more polyline branches on the cone surface.
    Polyline(Vec<Vec<DVec3>>),
    /// General case. Caller should fall back to marching.
    General,
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Main function
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the analytic intersection of `sphere` and `cone`.
pub fn intersect_sphere_cone(sphere: &SphericalSurface, cone: &ConicalSurface) -> SphereConeResult {
    intersect_sphere_cone_with_tolerance(sphere, cone, 0.0)
}

/// Compute sphere-cone intersection with additional fuzzy tolerance.
///
/// This relaxes axis-aligned and distance early-out checks by `fuzzy_tol` so
/// near-coincident cases can still classify into analytic branches.
pub fn intersect_sphere_cone_with_tolerance(
    sphere: &SphericalSurface,
    cone: &ConicalSurface,
    fuzzy_tol: f64,
) -> SphereConeResult {
    let tol = rcad_kernel::rcad_kernel::CONFUSION + fuzzy_tol.max(0.0);
    let axis = cone.axis_dir();
    let apex = cone.apex_point();

    // Project sphere center onto cone axis
    let d = sphere.center - apex;
    let z_c = d.dot(axis); // axial distance from cone apex to sphere center
    let foot = apex + axis * z_c;
    let d_perp = (sphere.center - foot).length();

    // 鈹€鈹€ Axis-aligned case: sphere center on cone axis 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if d_perp < tol * 10.0 {
        return intersect_sphere_cone_on_axis(sphere, cone, z_c, tol);
    }

    // 鈹€鈹€ Off-axis: try 胃-parameterized solver 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let result = intersect_sphere_cone_off_axis(sphere, cone);
    if !matches!(result, SphereConeResult::General) {
        return result;
    }

    // 鈹€鈹€ General case: numerical fallback 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    SphereConeResult::General
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Axis-aligned case
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn intersect_sphere_cone_on_axis(
    sphere: &SphericalSurface,
    cone: &ConicalSurface,
    z_c: f64,
    tol: f64,
) -> SphereConeResult {
    let axis = cone.axis_dir();
    let apex = cone.apex_point();
    let tan_half = cone.half_angle_rad.tan();
    let big_r = sphere.radius;
    let r_ref = cone.radius;

    // Cone radius at axial distance z from apex: r_cone(z) = r_ref + z * tan_half
    // Sphere cross-section radius at axial distance z from sphere center:
    //   r_sphere(z) = sqrt(big_r虏 - (z - z_c)虏)  when |z - z_c| <= big_r
    //
    // Intersection: r_sphere(z) = r_cone(z)
    // Let u = z - z_c (offset from sphere center along axis)
    //   sqrt(big_r虏 - u虏) = r_ref + (z_c + u) * tan_half
    // Square both sides:
    //   big_r虏 - u虏 = (r_ref + z_c * tan_half + u * tan_half)虏
    // This is a quartic in u. We solve by sampling and bisection.

    // Sampling range: u in [-big_r, +big_r]
    let n = 128usize;
    let mut roots: Vec<f64> = Vec::new();

    let f = |u: f64| -> f64 {
        if u.abs() > big_r {
            return f64::NAN;
        }
        let r_sphere_sq = big_r * big_r - u * u;
        if r_sphere_sq < 0.0 {
            return f64::NAN;
        }
        let r_sphere = r_sphere_sq.sqrt();
        let z = z_c + u;
        let r_cone = r_ref + z * tan_half;
        r_sphere - r_cone
    };

    let mut prev_u = -big_r;
    let mut prev_f = f(prev_u);

    for i in 1..=n {
        let u = -big_r + 2.0 * big_r * i as f64 / n as f64;
        let curr_f = f(u);

        if prev_f.is_nan() || curr_f.is_nan() {
            prev_u = u;
            prev_f = curr_f;
            continue;
        }

        // Sign change indicates a root
        if prev_f * curr_f < 0.0 {
            // Bisection
            let mut lo = prev_u;
            let mut hi = u;
            for _ in 0..64 {
                let mid = (lo + hi) * 0.5;
                let f_mid = f(mid);
                if f_mid.is_nan() {
                    break;
                }
                if f(lo) * f_mid < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            roots.push((lo + hi) * 0.5);
        } else if curr_f.abs() < tol {
            // Near-zero: check if this is a tangent point
            roots.push(u);
        }

        prev_u = u;
        prev_f = curr_f;
    }

    // Remove duplicate roots
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|a, b| (*a - *b).abs() < tol);

    // Convert roots to circles
    let mut circles: Vec<Circle3> = Vec::new();
    for u in roots {
        let z = z_c + u;
        let r_cone = r_ref + z * tan_half;

        // Skip if radius is negative or too small
        if r_cone < tol {
            // Tangent point at apex or near apex
            let pt = apex + axis * z;
            // Check if this is actually on the sphere
            if (pt - sphere.center).length() < big_r + tol {
                continue; // Will be handled as a point later if needed
            }
            continue;
        }

        // Verify the solution
        let r_sphere_sq = big_r * big_r - u * u;
        if r_sphere_sq < -tol {
            continue;
        }

        let center = apex + axis * z;
        circles.push(Circle3::new(center, axis, r_cone.max(0.0)));
    }

    match circles.len() {
        0 => SphereConeResult::NoIntersection,
        1 => SphereConeResult::SingleCircle(circles[0]),
        2 => SphereConeResult::TwoCircles(circles[0], circles[1]),
        _ => SphereConeResult::General, // More than 2 circles is unusual
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Off-axis 胃-parameterized solver
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Off-axis sphere-cone intersection using a 胃-parameterized quadratic solver.
///
/// Parameterizes the cone by 胃 (angle around axis) and s (axial distance from
/// true apex).  For each 胃, solves `A路s虏 + B(胃)路s + D = 0` where the sphere
/// constraint |P(胃,s) 锟?C|虏 = R虏 gives:
///
/// | Term | Expression |
/// |------|------------|
/// | A    | `1 + tan虏(尾)` (constant, 锟?1) |
/// | B(胃) | `锟?路(Cz + tan(尾)路(Cx路cos(胃) + Cy路sin(胃)))` |
/// | D    | `|C_local|虏 锟?R虏` (constant) |
///
/// Returns [`SphereConeResult::Polyline`] with one or more branches.
fn intersect_sphere_cone_off_axis(
    sphere: &SphericalSurface,
    cone: &ConicalSurface,
) -> SphereConeResult {
    let axis = cone.axis_dir();
    let tan_beta = cone.half_angle_rad.tan();
    if tan_beta.abs() < 1e-12 {
        return SphereConeResult::General; // degenerate (cylinder)
    }
    let apex_true = cone.apex_point();
    let sphere_r = sphere.radius;
    if sphere_r < rcad_kernel::rcad_kernel::CONFUSION {
        return SphereConeResult::General;
    }

    // Local basis: u, v 锟?axis
    let u = any_perpendicular(axis);
    let v = axis.cross(u).normalize();

    // Decompose sphere-center offset in the (u, v, axis) basis
    let cl = sphere.center - apex_true;
    let cx = cl.dot(u);
    let cy = cl.dot(v);
    let cz = cl.dot(axis);
    let local_sq = cl.length_squared();

    // Constant coefficients
    let a_coef = 1.0 + tan_beta * tan_beta; // A 锟?1
    let d_coef = local_sq - sphere_r * sphere_r;

    const N_THETA: usize = 256;
    let mut lower_branch: Vec<Option<(f64, DVec3)>> = Vec::with_capacity(N_THETA + 1);
    let mut upper_branch: Vec<Option<(f64, DVec3)>> = Vec::with_capacity(N_THETA + 1);

    for i in 0..=N_THETA {
        let theta = std::f64::consts::TAU * i as f64 / N_THETA as f64;
        let (cos_t, sin_t) = (theta.cos(), theta.sin());

        let b_theta = -2.0 * (cz + tan_beta * (cx * cos_t + cy * sin_t));
        let delta = b_theta * b_theta - 4.0 * a_coef * d_coef;

        if delta < 0.0 {
            lower_branch.push(None);
            upper_branch.push(None);
            continue;
        }

        let sqrt_delta = delta.sqrt();
        // Stable quadratic: "far" root via standard formula, "near" root via Vieta
        let s_far = (-b_theta - b_theta.signum() * sqrt_delta) / (2.0 * a_coef);
        let s_near = if s_far.abs() > 1e-15 {
            d_coef / (a_coef * s_far)
        } else {
            // s_far 锟?0 (D 锟?0, near-tangent case), compute directly
            (-b_theta + b_theta.signum() * sqrt_delta) / (2.0 * a_coef)
        };

        // Order so that s_lower 锟?s_upper
        let (s_lower, s_upper) = if s_far <= s_near {
            (s_far, s_near)
        } else {
            (s_near, s_far)
        };

        // Eval 3D point at given s on the 胃 ray
        let pt_at_s = |s: f64| -> DVec3 {
            let radial = s * tan_beta;
            apex_true + axis * s + radial * (u * cos_t + v * sin_t)
        };

        // Lower branch 锟?closer to apex
        if s_lower >= 0.0 {
            lower_branch.push(Some((theta, pt_at_s(s_lower))));
        } else {
            lower_branch.push(None);
        }

        // Upper branch 锟?further from apex (skip if coincident with lower)
        if s_upper >= 0.0 && (s_upper - s_lower).abs() > rcad_kernel::rcad_kernel::CONFUSION * 0.1 {
            upper_branch.push(Some((theta, pt_at_s(s_upper))));
        } else {
            upper_branch.push(None);
        }
    }

    // Extract contiguous valid runs from a branch array, handling 胃=0/2蟺 wrap.
    let extract_runs = |branch: &[Option<(f64, DVec3)>]| -> Vec<Vec<(f64, DVec3)>> {
        let n = branch.len();
        // Find first gap so we can rotate to avoid wrap issues
        let gap_at = branch.iter().position(|x| x.is_none()).unwrap_or(n);
        let mut rotated: Vec<Option<(f64, DVec3)>> = Vec::with_capacity(n);
        rotated.extend_from_slice(&branch[gap_at..]);
        rotated.extend_from_slice(&branch[..gap_at]);

        let mut curves: Vec<Vec<(f64, DVec3)>> = Vec::new();
        let mut current: Vec<(f64, DVec3)> = Vec::new();
        for pt in &rotated {
            match pt {
                Some(p) => current.push(*p),
                None => {
                    if current.len() >= 2 {
                        curves.push(current.clone());
                    }
                    current.clear();
                }
            }
        }
        if current.len() >= 2 {
            curves.push(current);
        }

        // Fix 胃 monotonicity in runs that wrap the 胃=0/2蟺 boundary
        // (rotation may place TAU before 0 in a single run)
        for curve in &mut curves {
            for j in 1..curve.len() {
                if curve[j].0 < curve[j - 1].0 - 1e-12 {
                    for k in j..curve.len() {
                        curve[k].0 += std::f64::consts::TAU;
                    }
                    break;
                }
            }
        }

        curves
    };

    let lower_param_branches = extract_runs(&lower_branch);
    let upper_param_branches = extract_runs(&upper_branch);

    // 鈹€鈹€ Adaptive chord-error refinement 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    const CHORD_TOL: f64 = crate::bop::int_tools::CHORD_TOLERANCE;
    const REFINE_DEPTH: usize = crate::bop::int_tools::CHORD_REFINE_DEPTH;

    // Lower-branch eval: recompute 3D point at any 胃, picking the near-apex root
    let lower_eval = {
        move |theta: f64| -> Option<DVec3> {
            let (cos_t, sin_t) = (theta.cos(), theta.sin());
            let b_theta = -2.0 * (cz + tan_beta * (cx * cos_t + cy * sin_t));
            let delta = b_theta * b_theta - 4.0 * a_coef * d_coef;
            if delta < 0.0 {
                return None;
            }
            let sqrt_delta = delta.sqrt();
            let s_far = (-b_theta - b_theta.signum() * sqrt_delta) / (2.0 * a_coef);
            let s_near = if s_far.abs() > 1e-15 {
                d_coef / (a_coef * s_far)
            } else {
                (-b_theta + b_theta.signum() * sqrt_delta) / (2.0 * a_coef)
            };
            let (s_lower, _s_upper) = if s_far <= s_near {
                (s_far, s_near)
            } else {
                (s_near, s_far)
            };
            if s_lower < 0.0 {
                return None;
            }
            let radial = s_lower * tan_beta;
            Some(apex_true + axis * s_lower + radial * (u * cos_t + v * sin_t))
        }
    };

    // Upper-branch eval: recompute 3D point at any 胃, picking the far-from-apex root
    let upper_eval = {
        move |theta: f64| -> Option<DVec3> {
            let (cos_t, sin_t) = (theta.cos(), theta.sin());
            let b_theta = -2.0 * (cz + tan_beta * (cx * cos_t + cy * sin_t));
            let delta = b_theta * b_theta - 4.0 * a_coef * d_coef;
            if delta < 0.0 {
                return None;
            }
            let sqrt_delta = delta.sqrt();
            let s_far = (-b_theta - b_theta.signum() * sqrt_delta) / (2.0 * a_coef);
            let s_near = if s_far.abs() > 1e-15 {
                d_coef / (a_coef * s_far)
            } else {
                (-b_theta + b_theta.signum() * sqrt_delta) / (2.0 * a_coef)
            };
            let (s_lower, s_upper) = if s_far <= s_near {
                (s_far, s_near)
            } else {
                (s_near, s_far)
            };
            if s_upper < 0.0 || (s_upper - s_lower).abs() <= rcad_kernel::rcad_kernel::CONFUSION * 0.1 {
                return None;
            }
            let radial = s_upper * tan_beta;
            Some(apex_true + axis * s_upper + radial * (u * cos_t + v * sin_t))
        }
    };

    let mut result: Vec<Vec<DVec3>> = Vec::new();

    for branch in &lower_param_branches {
        if branch.len() >= 4 {
            let refined = crate::bop::int_tools::pcurve_derive::refine_polyline(
                branch,
                &lower_eval,
                CHORD_TOL,
                REFINE_DEPTH,
            );
            result.push(refined.into_iter().map(|(_, p)| p).collect());
        } else if branch.len() >= 2 {
            result.push(branch.iter().map(|&(_, p)| p).collect());
        }
    }

    for branch in &upper_param_branches {
        if branch.len() >= 4 {
            let refined = crate::bop::int_tools::pcurve_derive::refine_polyline(
                branch,
                &upper_eval,
                CHORD_TOL,
                REFINE_DEPTH,
            );
            result.push(refined.into_iter().map(|(_, p)| p).collect());
        } else if branch.len() >= 2 {
            result.push(branch.iter().map(|&(_, p)| p).collect());
        }
    }

    // Closed-curve check: drop duplicate endpoint for runs where first 锟?last
    for branch in &mut result {
        if branch.len() >= 3 {
            let d = (branch[0] - branch[branch.len() - 1]).length();
            if d < rcad_kernel::rcad_kernel::CONFUSION * 10.0 {
                branch.pop();
            }
        }
    }

    // Filter very short branches and empty
    let result: Vec<Vec<DVec3>> = result.into_iter().filter(|b| b.len() >= 3).collect();

    if result.is_empty() {
        SphereConeResult::NoIntersection
    } else {
        SphereConeResult::Polyline(result)
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Tests
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
