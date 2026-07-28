//! Analytic intersection of two tori.
//!
//! # Case classification
//!
//! ## Coaxial tori (same axis line)
//!
//! When both tori share the same axis line, the intersection consists of circles
//! at axial heights where the torus tube circles (in the rho-z half-plane)
//! intersect each other.
//!
//! In the (rho, z) half-plane:
//! - Torus1 tube: (rho - R1)² + (z - z1)² = r1²  (circle centered at (R1, z1))
//! - Torus2 tube: (rho - R2)² + (z - z2)² = r2²  (circle centered at (R2, z2))
//!
//! The intersection of two circles can be 0, 1, or 2 points in the (rho, z) plane,
//! which by rotational symmetry gives 0, 1, or 2 circles in 3D.
//!
//! ## Tangent case
//!
//! When the two tube circles are tangent (touching), the 3D intersection is a
//! single tangent circle.
//!
//! ## Skew case (non-coaxial, non-parallel axes)
//!
//! Solved analytically by parameterizing one torus P₁(u,v) and substituting
//! into the other torus's implicit equation:
//!
//! ```text
//! (|P|² + R₂² - r₂²)² = 4R₂²(|P|² - (P·a₂)²)
//! ```
//!
//! At each u ∈ [0, 2π) the substitution yields a quadratic trigonometric
//! equation in v (the tube angle).  Via t = tan(v/2) this becomes a monic
//! quartic in t, solved via Ferrari's method.
//!
//! ## General case
//!
//! For all other configurations (skew axes, offset axes) the intersection is a
//! complex space curve. We return `General` so the caller falls back to numeric
//! marching.

use std::f64::consts::TAU;

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Circle3, ToroidalSurface, any_perpendicular};

use crate::tolerance::*;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of torus x torus intersection.
#[derive(Debug, Clone)]
pub enum TorusTorusResult {
    /// The tori do not intersect.
    NoIntersection,
    /// Coaxial case: single intersection circle.
    SingleCircle(Circle3),
    /// Coaxial case: two intersection circles.
    TwoCircles(Circle3, Circle3),
    /// The tori are tangent, giving one tangent circle.
    TangentCircle(Circle3),
    /// Coaxial tori with identical geometry (same axis, same radii, same center).
    Coaxial,
    /// Skew (non-coaxial) intersection: analytic quartic solution via
    /// torus parameterization.  Returns one or more 3D polylines.
    SkewQuartic(Vec<Vec<DVec3>>),
    /// General case. Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `t1` and `t2`.
pub fn intersect_torus_torus(t1: &ToroidalSurface, t2: &ToroidalSurface) -> TorusTorusResult {
    intersect_torus_torus_with_tolerance(t1, t2, 0.0)
}

/// Compute torus-torus intersection with additional fuzzy tolerance.
///
/// This relaxes coaxial detection by `fuzzy_tol` so near-coaxial cases
/// can still classify into analytic branches.
pub fn intersect_torus_torus_with_tolerance(
    t1: &ToroidalSurface,
    t2: &ToroidalSurface,
    fuzzy_tol: f64,
) -> TorusTorusResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    let a1 = t1.axis.normalize();
    let a2 = t2.axis.normalize();
    let cross = a1.cross(a2);
    let sin_angle = cross.length();

    // Project t2 center onto t1 axis
    let delta = t2.center - t1.center;
    let d_perp = (delta - a1 * delta.dot(a1)).length();

    // ── Coaxial: same axis line ───────────────────────────────────────────────
    if sin_angle < TOLERANCE_ANG && d_perp < tol {
        return intersect_torus_torus_coaxial(t1, t2, a1);
    }

    // ── Skew (non-coaxial): try analytic quartic solver ──────────────────────
    // Parameterize t1 via (u,v) and substitute into t2's implicit.
    let skew_result = intersect_skew_torus_torus(t1, t2);
    if !skew_result.is_empty() {
        return TorusTorusResult::SkewQuartic(skew_result);
    }

    // ── General case: numerical fallback ─────────────────────────────────────
    TorusTorusResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Coaxial case
// ─────────────────────────────────────────────────────────────────────────────

#[allow(non_snake_case)]
fn intersect_torus_torus_coaxial(
    t1: &ToroidalSurface,
    t2: &ToroidalSurface,
    axis: DVec3,
) -> TorusTorusResult {
    let R1 = t1.major_radius;
    let r1 = t1.minor_radius;
    let R2 = t2.major_radius;
    let r2 = t2.minor_radius;

    // Check for identical tori
    let dz_centers = (t2.center - t1.center).dot(axis);
    if dz_centers.abs() < TOLERANCE_ABS
        && (R1 - R2).abs() < TOLERANCE_ABS
        && (r1 - r2).abs() < TOLERANCE_ABS
    {
        return TorusTorusResult::Coaxial;
    }

    // Two circles in (rho, z) plane:
    // Circle 1: (rho - R1)² + (z - 0)² = r1²   (center at (R1, 0))
    // Circle 2: (rho - R2)² + (z - dz)² = r2²  (center at (R2, dz))
    //
    // Solve for intersection of two circles.
    // Let u = rho, v = z. Then:
    //   (u - R1)² + v² = r1²
    //   (u - R2)² + (v - dz)² = r2²
    //
    // Expand both:
    //   u² - 2*R1*u + R1² + v² = r1²
    //   u² - 2*R2*u + R2² + v² - 2*dz*v + dz² = r2²
    //
    // Subtract: -2*(R1 - R2)*u + R1² - R2² + 2*dz*v - dz² = r1² - r2²
    //   v = (r1² - r2² + 2*(R1 - R2)*u + dz² - R1² + R2²) / (2*dz)  [if dz != 0]
    //
    // If dz = 0: circles are concentric in (rho, z) → intersection only if
    // tube circles touch.

    if dz_centers.abs() < TOLERANCE_ABS {
        // Concentric tori: intersection only if tubes touch
        // Distance between tube centers in (rho, z) plane is |R1 - R2|
        let d_tube = (R1 - R2).abs();

        // Check for tangent tubes
        if (d_tube - (r1 + r2)).abs() < TOLERANCE_ABS {
            // Tubes touch at one circle at z = 0, rho = midpoint
            let rho = (R1 + R2) / 2.0;
            if rho > TOLERANCE_ABS {
                return TorusTorusResult::TangentCircle(Circle3::new(t1.center, axis, rho));
            }
        }

        // Check for one tube inside the other (no intersection)
        if d_tube + r1.min(r2) < r1.max(r2) - TOLERANCE_ABS {
            return TorusTorusResult::NoIntersection;
        }

        // Overlapping tubes with same major radius
        if (R1 - R2).abs() < TOLERANCE_ABS && (r1 - r2).abs() < TOLERANCE_ABS {
            return TorusTorusResult::Coaxial;
        }

        // General concentric case: tubes may intersect at multiple heights
        // Fall back to numerical for now
        return TorusTorusResult::General;
    }

    // Linear relation: v = A*u + B
    let A = (R1 - R2) / dz_centers;
    let B = (r1 * r1 - r2 * r2 + dz_centers * dz_centers - R1 * R1 + R2 * R2) / (2.0 * dz_centers);

    // Substitute into circle 1: (u - R1)² + (A*u + B)² = r1²
    // (1 + A²)*u² + (-2*R1 + 2*A*B)*u + (R1² + B² - r1²) = 0
    let a_q = 1.0 + A * A;
    let b_q = -2.0 * R1 + 2.0 * A * B;
    let c_q = R1 * R1 + B * B - r1 * r1;

    let disc = b_q * b_q - 4.0 * a_q * c_q;

    if disc < -TOLERANCE_ABS {
        return TorusTorusResult::NoIntersection;
    }

    if disc.abs() < TOLERANCE_ABS {
        // Tangent: one solution
        let u = -b_q / (2.0 * a_q);
        if u < TOLERANCE_ABS {
            return TorusTorusResult::NoIntersection;
        }
        let v = A * u + B;
        let center = t1.center + axis * v;
        return TorusTorusResult::TangentCircle(Circle3::new(center, axis, u));
    }

    // Two solutions
    let sqrt_disc = disc.sqrt();
    let u1 = (-b_q - sqrt_disc) / (2.0 * a_q);
    let u2 = (-b_q + sqrt_disc) / (2.0 * a_q);

    let mut circles: Vec<Circle3> = Vec::new();

    for u in [u1, u2] {
        if u > TOLERANCE_ABS {
            let v = A * u + B;
            let center = t1.center + axis * v;
            circles.push(Circle3::new(center, axis, u));
        }
    }

    match circles.len() {
        0 => TorusTorusResult::NoIntersection,
        1 => TorusTorusResult::SingleCircle(circles[0]),
        2 => TorusTorusResult::TwoCircles(circles[0], circles[1]),
        _ => TorusTorusResult::General,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Skew (non-coaxial) case — analytic quartic solver
// ─────────────────────────────────────────────────────────────────────────────

/// Skew torus × torus intersection via torus-parameterized quartic solver.
///
/// Parameterize torus 1 P₁(u,v) and substitute into torus 2's implicit:
///
/// ```text
/// (|P|² + R₂² - r₂²)² = 4R₂²(|P|² - (P·a₂)²)
/// ```
///
/// At each u the substitution yields a quadratic trigonometric equation in v:
///
/// ```text
/// A_const + A_cos·cos(v) + A_sin·sin(v) + A_cos2·cos²(v)
///     + A_sin2·sin²(v) + A_cossin·cos(v)·sin(v) = 0
/// ```
///
/// Via t = tan(v/2):
/// ```text
/// a₄·t⁴ + a₃·t³ + a₂·t² + a₁·t + a₀ = 0
/// ```
/// Solved via Ferrari's method (crate::solve_quartic).
fn intersect_skew_torus_torus(t1: &ToroidalSurface, t2: &ToroidalSurface) -> Vec<Vec<DVec3>> {
    let a1 = t1.axis.normalize();
    let a2 = t2.axis.normalize();
    let o1 = t1.center;
    let o2 = t2.center;

    let d = o1 - o2; // t1 center in t2-centered frame

    let r_major1 = t1.major_radius;
    let r_minor1 = t1.minor_radius;
    let r_major2 = t2.major_radius;
    let r_minor2 = t2.minor_radius;

    let r_major1_sq = r_major1 * r_major1;
    let r_minor1_sq = r_minor1 * r_minor1;
    let r_major2_sq = r_major2 * r_major2;
    let r_minor2_sq = r_minor2 * r_minor2;
    let c2_sq = r_major2_sq - r_minor2_sq; // R₂² - r₂²
    let c2_sq_sq = c2_sq * c2_sq; // (R₂² - r₂²)²

    // Perpendicular basis for torus 1
    let x1 = any_perpendicular(a1);
    let y1 = a1.cross(x1).normalize();

    // ── Pre-computed constants (independent of u) ──────────────────────────
    let d_sq = d.length_squared();
    let d_dot_a1 = d.dot(a1);
    let d_dot_a2 = d.dot(a2);
    let a1_dot_a2 = a1.dot(a2);

    // Components of D and a2 in the plane of torus 1
    let d_dot_x1 = d.dot(x1);
    let d_dot_y1 = d.dot(y1);
    let a2_dot_x1 = a2.dot(x1);
    let a2_dot_y1 = a2.dot(y1);

    // Constant (u-independent) parts of the trigonometric coefficients
    // sin_S = coefficient of sin(v) in |P₁ - O₂|²
    let sin_s = 2.0 * r_minor1 * d_dot_a1;

    // sin_T = coefficient of sin(v) in (P₁ - O₂)·a₂
    let sin_t = r_minor1 * a1_dot_a2;

    // Constant sub-expressions for A_sin2 and parts of A_const, A_sin
    let term_r2_sq_minus_m2_sq = r_minor2_sq - r_major2_sq;

    let sin_s_sq = sin_s * sin_s;
    let sin_t_sq = sin_t * sin_t;
    let four_r2_sq = 4.0 * r_major2_sq;

    // A_sin2 and the sin_T contribution to A_cossin are constant:
    let a_sin2_const = sin_s_sq + four_r2_sq * sin_t_sq;
    let a_cossin_sin_t_const = four_r2_sq * sin_t;

    const N_SAMPLES: usize = 128;
    let mut samples: Vec<(f64, f64, DVec3)> = Vec::new();

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cu, su) = (u.cos(), u.sin());

        // r_xy(u) = cos(u)·x₁ + sin(u)·y₁
        let drxy = d_dot_x1 * cu + d_dot_y1 * su; // D · r_xy(u)
        let rxy_a2 = a2_dot_x1 * cu + a2_dot_y1 * su; // r_xy(u) · a₂

        // ── Compute trigonometric coefficients ────────────────────────────
        // S = |P₁ - O₂|² = M + N·cos(v) + sin_S·sin(v)
        let m = d_sq + r_major1_sq + r_minor1_sq + 2.0 * r_major1 * drxy;
        let n = 2.0 * r_minor1 * (drxy + r_major1);

        // T = (P₁ - O₂)·a₂ = Q + R_cos·cos(v) + sin_T·sin(v)
        let q = d_dot_a2 + r_major1 * rxy_a2;
        let r_cos = r_minor1 * rxy_a2;

        let m_sq = m * m;
        let n_sq = n * n;
        let q_sq = q * q;
        let r_cos_sq = r_cos * r_cos;

        // A_const = M² + 2(r₂² - R₂²)M + 4R₂²Q² + (R₂² - r₂²)²
        let a_const = m_sq + 2.0 * term_r2_sq_minus_m2_sq * m + four_r2_sq * q_sq + c2_sq_sq;

        // A_cos = 2MN + 2(r₂² - R₂²)N + 8R₂²·Q·R_cos
        let a_cos = 2.0 * m * n + 2.0 * term_r2_sq_minus_m2_sq * n + 8.0 * r_major2_sq * q * r_cos;

        // A_sin = 2M·sin_S + 2(r₂² - R₂²)·sin_S + 8R₂²·Q·sin_T
        let a_sin =
            2.0 * (m * sin_s + term_r2_sq_minus_m2_sq * sin_s) + 8.0 * r_major2_sq * q * sin_t;

        // A_cos2 = N² + 4R₂²·R_cos²
        let a_cos2 = n_sq + four_r2_sq * r_cos_sq;

        // A_sin2 = sin_S² + 4R₂²·sin_T²  (constant)
        let a_sin2 = a_sin2_const;

        // A_cossin = 2N·sin_S + 8R₂²·R_cos·sin_T
        let a_cossin = 2.0 * n * sin_s + a_cossin_sin_t_const * r_cos;

        // ── Quartic in t = tan(v/2) ───────────────────────────────────────
        let a4 = a_const - a_cos + a_cos2;
        let a3 = 2.0 * (a_sin - a_cossin);
        let a2 = 2.0 * a_const - 2.0 * a_cos2 + 4.0 * a_sin2;
        let a1 = 2.0 * (a_sin + a_cossin);
        let a0 = a_const + a_cos + a_cos2;

        let t_roots = crate::solve_quartic(a4, a3, a2, a1, a0);

        // Convert t = tan(v/2) → v ∈ [0, 2π)
        let mut v_vals: Vec<f64> = Vec::new();
        for &t in &t_roots {
            if t.is_finite() {
                let v = 2.0 * t.atan();
                let v_norm = if v < 0.0 { v + TAU } else { v };
                v_vals.push(v_norm);
            }
        }

        // Check for root at/near v = π (t → ∞).
        // The value of the trigonometric equation at v = π equals a₄.
        if a4.abs() < 1e-6 {
            v_vals.push(std::f64::consts::PI);
        }

        // Deduplicate v values
        v_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v_vals.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

        for &v in &v_vals {
            if v.is_finite() {
                let p = t1.point_at(u, v);
                if p.is_finite() {
                    samples.push((u, v, p));
                }
            }
        }
    }

    if samples.len() < 2 {
        return vec![];
    }

    // ── Branch extraction ──────────────────────────────────────────────────
    samples.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });

    let mut branches: Vec<Vec<(f64, DVec3)>> = Vec::new();

    let mut i = 0;
    while i < samples.len() {
        let u_cur = samples[i].0;
        let mut cur_pts: Vec<(f64, DVec3)> = Vec::new();
        while i < samples.len() && (samples[i].0 - u_cur).abs() < 1e-14 {
            cur_pts.push((samples[i].1, samples[i].2));
            i += 1;
        }
        cur_pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        if branches.is_empty() {
            for &(_v, p) in &cur_pts {
                branches.push(vec![(u_cur, p)]);
            }
        } else {
            let v_threshold = if cur_pts.len() >= 2 {
                let v_range = cur_pts.last().unwrap().0 - cur_pts[0].0;
                (v_range / cur_pts.len() as f64) * 2.5
            } else {
                (TAU / N_SAMPLES as f64) * 20.0
            };

            let mut assigned = vec![false; cur_pts.len()];
            for branch in &mut branches {
                let last_v = branch.last().unwrap().0;
                let mut best_idx = None;
                let mut best_dist = f64::MAX;
                for (j, &(v, _)) in cur_pts.iter().enumerate() {
                    if !assigned[j] {
                        let dist = (v - last_v).abs();
                        if dist < best_dist && dist < v_threshold {
                            best_dist = dist;
                            best_idx = Some(j);
                        }
                    }
                }
                if let Some(idx) = best_idx {
                    branch.push((u_cur, cur_pts[idx].1));
                    assigned[idx] = true;
                }
            }
            for (j, &(_, p)) in cur_pts.iter().enumerate() {
                if !assigned[j] {
                    branches.push(vec![(u_cur, p)]);
                }
            }
        }
    }

    // ── Adaptive chord-error refinement ────────────────────────────────────
    const CHORD_TOL: f64 = crate::bop::int_tools::CHORD_TOLERANCE;
    const REFINE_DEPTH: usize = crate::bop::int_tools::CHORD_REFINE_DEPTH;

    let refined_3d: Vec<Vec<DVec3>> = branches
        .into_iter()
        .filter(|b| b.len() >= 4)
        .map(|branch| {
            let pts_for_eval = branch.clone();
            let eval_fn = move |u_mid: f64| -> Option<DVec3> {
                let idx = pts_for_eval.partition_point(|&(u, _)| u < u_mid);
                if idx == 0 || idx >= pts_for_eval.len() {
                    return None;
                }
                let (u0, p0) = pts_for_eval[idx - 1];
                let (u1, p1) = pts_for_eval[idx];
                let t_frac = ((u_mid - u0) / (u1 - u0)).clamp(0.0, 1.0);
                let expected = p0.lerp(p1, t_frac);

                let (cu, su) = (u_mid.cos(), u_mid.sin());
                let drxy = d_dot_x1 * cu + d_dot_y1 * su;
                let rxy_a2 = a2_dot_x1 * cu + a2_dot_y1 * su;

                let m_val = d_sq + r_major1_sq + r_minor1_sq + 2.0 * r_major1 * drxy;
                let n_val = 2.0 * r_minor1 * (drxy + r_major1);
                let q_val = d_dot_a2 + r_major1 * rxy_a2;
                let r_cos_val = r_minor1 * rxy_a2;

                let m_sq = m_val * m_val;
                let n_sq = n_val * n_val;
                let q_sq = q_val * q_val;
                let r_cos_sq = r_cos_val * r_cos_val;

                let a_const =
                    m_sq + 2.0 * term_r2_sq_minus_m2_sq * m_val + four_r2_sq * q_sq + c2_sq_sq;
                let a_cos = 2.0 * m_val * n_val
                    + 2.0 * term_r2_sq_minus_m2_sq * n_val
                    + 8.0 * r_major2_sq * q_val * r_cos_val;
                let a_sin = 2.0 * (m_val * sin_s + term_r2_sq_minus_m2_sq * sin_s)
                    + 8.0 * r_major2_sq * q_val * sin_t;
                let a_cos2 = n_sq + four_r2_sq * r_cos_sq;
                let a_sin2 = a_sin2_const;
                let a_cossin = 2.0 * n_val * sin_s + a_cossin_sin_t_const * r_cos_val;

                let a4 = a_const - a_cos + a_cos2;
                let a3_v = 2.0 * (a_sin - a_cossin);
                let a2_v = 2.0 * a_const - 2.0 * a_cos2 + 4.0 * a_sin2;
                let a1_v = 2.0 * (a_sin + a_cossin);
                let a0_v = a_const + a_cos + a_cos2;

                let t_roots = crate::solve_quartic(a4, a3_v, a2_v, a1_v, a0_v);

                // Convert t = tan(v/2) → v ∈ [0, 2π)
                let mut v_candidates: Vec<f64> = t_roots
                    .iter()
                    .filter(|t| t.is_finite())
                    .map(|&t| {
                        let v = 2.0 * t.atan();
                        if v < 0.0 { v + TAU } else { v }
                    })
                    .collect();

                // Check for root at/near v = π (t → ∞)
                if a4.abs() < 1e-6 {
                    v_candidates.push(std::f64::consts::PI);
                }

                v_candidates
                    .iter()
                    .filter(|v| v.is_finite())
                    .map(|&v| (v, t1.point_at(u_mid, v)))
                    .filter(|(_, p)| p.is_finite())
                    .min_by(|(_, pa), (_, pb)| {
                        (pa - expected)
                            .length_squared()
                            .partial_cmp(&(pb - expected).length_squared())
                            .unwrap()
                    })
                    .map(|(_, p)| p)
            };

            let refined = crate::bop::int_tools::pcurve_derive::refine_polyline(
                &branch,
                eval_fn,
                CHORD_TOL,
                REFINE_DEPTH,
            );
            refined.into_iter().map(|(_, p)| p).collect()
        })
        .collect();

    // Dedup near-duplicate trailing points
    let mut result = Vec::new();
    for mut branch in refined_3d {
        while branch.len() >= 3 {
            let n = branch.len();
            if (branch[n - 1] - branch[0]).length_squared() < TOLERANCE_VEC_SQ_MIN {
                branch.pop();
            } else {
                break;
            }
        }
        if branch.len() >= 2 {
            result.push(branch);
        }
    }

    result
}
