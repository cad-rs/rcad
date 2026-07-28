//! Analytic intersection of a torus and a cone.
//!
//! # Case classification
//!
//! ## Coaxial case (cone apex on torus axis, cone axis = torus axis)
//!
//! When both surfaces share the same axis line and the cone apex lies on the
//! torus axis, the intersection consists of circles at axial heights where
//! the torus tube circle (in the rho-z half-plane) intersects the cone line.
//!
//! In the (rho, z) half-plane:
//! - Torus tube: (rho - R)² + z² = r²  (circle centered at (R, 0))
//! - Cone: rho = (z - z_apex) * tan(half_angle)  (line from apex)
//!
//! Substituting gives a quadratic equation in z.
//!
//! ## Skew case (non-coaxial, non-parallel axes)
//!
//! Solved analytically by parameterizing the cone surface P(u,v) and
//! substituting into the torus implicit equation:
//!
//! ```text
//! (|P|² + R² - r²)² = 4R²(|P|² - (P·a_tor)²)
//! ```
//!
//! The cone P(u,v) is linear in v (slant distance), so the substitution
//! yields a monic quartic in v at each azimuth u ∈ [0, 2π). Solved via
//! Ferrari's method (same approach as cylinder-torus).
//!
//! ## General case
//!
//! For all other configurations the intersection is a complex space curve.
//! We return `General` so the caller falls back to numeric marching.

use std::f64::consts::TAU;

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Circle3, ConicalSurface, ToroidalSurface, any_perpendicular};

use crate::tolerance::*;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of torus x cone intersection.
#[derive(Debug, Clone)]
pub enum TorusConeResult {
    /// The torus and cone do not intersect.
    NoIntersection,
    /// Coaxial case: single intersection circle.
    SingleCircle(Circle3),
    /// Coaxial case: two intersection circles.
    TwoCircles(Circle3, Circle3),
    /// The intersection is a tangent circle.
    TangentCircle(Circle3),
    /// Skew (non-coaxial) intersection: analytic quartic solution via
    /// cone parameterization.  Returns one or more 3D polylines.
    SkewQuartic(Vec<Vec<DVec3>>),
    /// General case. Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `torus` and `cone`.
pub fn intersect_torus_cone(torus: &ToroidalSurface, cone: &ConicalSurface) -> TorusConeResult {
    intersect_torus_cone_with_tolerance(torus, cone, 0.0)
}

/// Compute torus-cone intersection with additional fuzzy tolerance.
///
/// This relaxes coaxial detection by `fuzzy_tol` so near-coaxial cases
/// can still classify into analytic branches.
pub fn intersect_torus_cone_with_tolerance(
    torus: &ToroidalSurface,
    cone: &ConicalSurface,
    fuzzy_tol: f64,
) -> TorusConeResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    let t_axis = torus.axis.normalize();
    let c_axis = cone.axis_dir();
    let cross = t_axis.cross(c_axis);
    let sin_angle = cross.length();
    let apex = cone.apex_point();

    // Project cone apex onto torus axis
    let t = (apex - torus.center).dot(t_axis);
    let foot = torus.center + t_axis * t;
    let d_apex = (apex - foot).length();

    // ── Coaxial: same axis line ───────────────────────────────────────────────
    if sin_angle < TOLERANCE_ANG && d_apex < tol {
        return intersect_torus_cone_coaxial(torus, cone, t_axis);
    }

    // ── Skew (non-coaxial, non-parallel): try analytic quartic solver ────────
    // Parameterize the cone and substitute into the torus implicit.
    // Returns polyline branches when the quartic has real roots.
    let skew_result = intersect_skew_torus_cone(torus, cone);
    if !skew_result.is_empty() {
        return TorusConeResult::SkewQuartic(skew_result);
    }

    // ── General case: numerical fallback ─────────────────────────────────────
    TorusConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Coaxial case
// ─────────────────────────────────────────────────────────────────────────────

#[allow(non_snake_case)]
fn intersect_torus_cone_coaxial(
    torus: &ToroidalSurface,
    cone: &ConicalSurface,
    axis: DVec3,
) -> TorusConeResult {
    let R = torus.major_radius;
    let r = torus.minor_radius;
    let ta = cone.half_angle_rad.tan();

    // Determine cone orientation relative to torus axis
    let sigma = if cone.axis_dir().dot(axis) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let apex = cone.apex_point();
    let r_ref = cone.radius;

    // z_apex: axial coordinate of cone apex relative to torus center
    let z_apex = (apex - torus.center).dot(axis);

    // In the (rho, z) half-plane:
    // Torus tube: (rho - R)² + z² = r²   (tube center at (R, 0))
    // Cone: rho = r_ref + sigma * (z - z_apex) * ta
    //
    // Let's set up the equation:
    // Let a_cone = sigma * ta  (cone slope in rho-z plane, signed)
    // rho_cone(z) = r_ref + a_cone * (z - z_apex)
    //
    // Substitute into torus equation:
    // (r_ref + a_cone*(z - z_apex) - R)² + z² = r²
    //
    // Let A = a_cone, B = r_ref - R - A*z_apex
    // Then rho_cone = A*z + B
    //
    // (A*z + B - R)² + z² = r²
    // A²*z² + 2*A*(B-R)*z + (B-R)² + z² = r²
    // (A² + 1)*z² + 2*A*(B-R)*z + (B-R)² - r² = 0

    let A = sigma * ta;
    let B = r_ref - A * z_apex;
    let rho_offset = B - R;

    let a_q = A * A + 1.0;
    let b_q = 2.0 * A * rho_offset;
    let c_q = rho_offset * rho_offset - r * r;

    let disc = b_q * b_q - 4.0 * a_q * c_q;

    if disc < -TOLERANCE_ABS {
        return TorusConeResult::NoIntersection;
    }

    if disc.abs() < TOLERANCE_ABS {
        // Tangent: one solution
        let z = -b_q / (2.0 * a_q);
        let rho = A * z + B;

        if rho < TOLERANCE_ABS {
            return TorusConeResult::NoIntersection;
        }

        let center = torus.center + axis * z;
        return TorusConeResult::TangentCircle(Circle3::new(center, axis, rho));
    }

    // Two solutions
    let sqrt_disc = disc.sqrt();
    let z1 = (-b_q - sqrt_disc) / (2.0 * a_q);
    let z2 = (-b_q + sqrt_disc) / (2.0 * a_q);

    let mut circles: Vec<Circle3> = Vec::new();

    for z in [z1, z2] {
        let rho = A * z + B;
        if rho > TOLERANCE_ABS {
            let center = torus.center + axis * z;
            circles.push(Circle3::new(center, axis, rho));
        }
    }

    match circles.len() {
        0 => TorusConeResult::NoIntersection,
        1 => TorusConeResult::SingleCircle(circles[0]),
        2 => TorusConeResult::TwoCircles(circles[0], circles[1]),
        _ => TorusConeResult::General,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Skew (non-coaxial) case — analytic quartic solver
// ─────────────────────────────────────────────────────────────────────────────

/// Skew torus × cone intersection via cone-parameterized quartic solver.
///
/// Parameterize the cone surface P(u,v) and substitute into the torus implicit
/// equation.  At each cone azimuth u ∈ [0, 2π) the equation becomes a monic
/// quartic in v (slant distance).  Solve via Ferrari's method and extract
/// polyline branches by proximity clustering across consecutive u values.
///
/// # Cone parameterization
///
/// ```text
/// P(u,v) = apex + v·cos(α)·a_cone + (r_ref + v·sin(α))·r_dir(u)
///
/// where   r_dir(u) = cos(u)·x_cone + sin(u)·y_cone
/// ```
///
/// This is linear in v: P(u,v) = B0(u) + v·B1(u).
///
/// # Torus implicit equation
///
/// ```text
/// (|P|² + R² - r²)² = 4R²(|P|² - (P·a_tor)²)
/// ```
///
/// Substituting P(u,v) and grouping by powers of v gives a monic quartic:
///   v⁴ + a₃(u)·v³ + a₂(u)·v² + a₁(u)·v + a₀(u) = 0
///
/// The coefficients are derived from:
/// ```
/// |P|² = S0(u) + v·S1(u) + v²     (since |B1|² = 1 for the cone)
/// P·a_tor = T0(u) + v·T1(u)
/// ```
fn intersect_skew_torus_cone(torus: &ToroidalSurface, cone: &ConicalSurface) -> Vec<Vec<DVec3>> {
    let a_tor = torus.axis.normalize();
    let a_cone = cone.axis_dir();
    let apex = cone.apex_point();
    let r_ref = cone.radius;

    let o = apex - torus.center; // cone apex in torus-centered frame

    let r_major = torus.major_radius;
    let r_minor = torus.minor_radius;
    let r_major_sq = r_major * r_major;
    let r_minor_sq = r_minor * r_minor;
    let c_sq = r_major_sq - r_minor_sq; // R² - r²
    let c_sq_sq = c_sq * c_sq; // (R² - r²)²

    let k_sin = cone.half_angle_rad.sin();
    let k_cos = cone.half_angle_rad.cos();

    // Perpendicular basis for cone (matches ConicalSurface::point_at)
    let x_cone = any_perpendicular(a_cone);
    let y_cone = a_cone.cross(x_cone).normalize();

    // ── Pre-computed constants (independent of u) ──────────────────────────
    let apex_sq = o.length_squared();
    let a_dot_at = a_cone.dot(a_tor);
    let apex_o_at = o.dot(a_tor); // apex · a_tor
    let apex_o_ac = o.dot(a_cone); // apex · a_cone

    // Components of apex-in-torus-frame and a_tor in the cone's perpendicular plane
    let ox = o.dot(x_cone);
    let oy = o.dot(y_cone);
    let cx = a_tor.dot(x_cone);
    let cy = a_tor.dot(y_cone);

    // Constant (u-independent) part of S1:
    //   S1 = 2·(k_cos·apex_o_ac + k_sin·(apex_r + r_ref))
    //       = 2·k_cos·apex_o_ac + 2·k_sin·r_ref + 2·k_sin·apex_r
    //   The u-dependent part is 2·k_sin·apex_r.
    let s1_const = 2.0 * (k_cos * apex_o_ac + k_sin * r_ref);

    const N_SAMPLES: usize = 128;
    let mut samples: Vec<(f64, f64, DVec3)> = Vec::new();

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cu, su) = (u.cos(), u.sin());

        // r_dir(u) = cos(u)·x_cone + sin(u)·y_cone
        let apex_r = ox * cu + oy * su; // apex · r_dir(u)
        let at_r = cx * cu + cy * su; // a_tor · r_dir(u)

        // ── Compute quartic coefficients ───────────────────────────────────
        // S0 = |B0|² = |o|² + r_ref² + 2·r_ref·(o·r_dir)
        let s0 = apex_sq + r_ref * r_ref + 2.0 * r_ref * apex_r;

        // S1 = 2·B0·B1 = 2·(k_cos·apex_o_ac + k_sin·(apex_r + r_ref))
        let s1 = s1_const + 2.0 * k_sin * apex_r;

        // T0 = B0·a_tor = apex_o_at + r_ref·(r_dir·a_tor)
        let t0 = apex_o_at + r_ref * at_r;

        // T1 = B1·a_tor = k_cos·(a_cone·a_tor) + k_sin·(r_dir·a_tor)
        let t1 = k_cos * a_dot_at + k_sin * at_r;

        // Monic quartic: v⁴ + a₃·v³ + a₂·v² + a₁·v + a₀ = 0
        let a3 = 2.0 * s1;
        let a2 =
            2.0 * s0 - 2.0 * r_major_sq - 2.0 * r_minor_sq + s1 * s1 + 4.0 * r_major_sq * t1 * t1;
        let a1 = 2.0 * s1 * (s0 - r_major_sq - r_minor_sq) + 8.0 * r_major_sq * t0 * t1;
        let a0 =
            s0 * s0 - 2.0 * s0 * (r_major_sq + r_minor_sq) + c_sq_sq + 4.0 * r_major_sq * t0 * t0;

        let v_roots = crate::solve_quartic(1.0, a3, a2, a1, a0);

        for &v in &v_roots {
            if v.is_finite() {
                let p = cone.point_at(u, v);
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
    // Sort by u, then v.
    samples.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });

    // Cluster points into branches by proximity across consecutive u values.
    let mut branches: Vec<Vec<(f64, DVec3)>> = Vec::new();

    let mut i = 0;
    while i < samples.len() {
        let u_cur = samples[i].0;
        // Collect all points at this u
        let mut cur_pts: Vec<(f64, DVec3)> = Vec::new();
        while i < samples.len() && (samples[i].0 - u_cur).abs() < 1e-14 {
            cur_pts.push((samples[i].1, samples[i].2)); // (v, point)
            i += 1;
        }
        // Sort within u by v
        cur_pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        if branches.is_empty() {
            // Start a branch for each point at this u
            for &(_v, p) in &cur_pts {
                branches.push(vec![(u_cur, p)]);
            }
        } else {
            // Determine v-matching threshold
            let v_threshold = if cur_pts.len() >= 2 {
                let v_range = cur_pts.last().unwrap().0 - cur_pts[0].0;
                (v_range / cur_pts.len() as f64) * 2.5
            } else {
                (TAU / N_SAMPLES as f64) * 20.0
            };

            // For each existing branch, find the closest v at the current u
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
            // Unassigned points start new branches
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
            // branch: Vec<(f64, DVec3)> = (u, point)
            let pts_for_eval = branch.clone();
            let eval_fn = move |u_mid: f64| -> Option<DVec3> {
                // Find neighboring points by binary search on u
                let idx = pts_for_eval.partition_point(|&(u, _)| u < u_mid);
                if idx == 0 || idx >= pts_for_eval.len() {
                    return None;
                }
                let (u0, p0) = pts_for_eval[idx - 1];
                let (u1, p1) = pts_for_eval[idx];
                let t = ((u_mid - u0) / (u1 - u0)).clamp(0.0, 1.0);
                let expected = p0.lerp(p1, t);

                // Recompute quartic at u_mid
                let (cu, su) = (u_mid.cos(), u_mid.sin());
                let apex_r = ox * cu + oy * su;
                let at_r = cx * cu + cy * su;
                let s0_val = apex_sq + r_ref * r_ref + 2.0 * r_ref * apex_r;
                let s1_val = s1_const + 2.0 * k_sin * apex_r;
                let t0_val = apex_o_at + r_ref * at_r;
                let t1_val = k_cos * a_dot_at + k_sin * at_r;
                let a3 = 2.0 * s1_val;
                let a2 = 2.0 * s0_val - 2.0 * r_major_sq - 2.0 * r_minor_sq
                    + s1_val * s1_val
                    + 4.0 * r_major_sq * t1_val * t1_val;
                let a1 = 2.0 * s1_val * (s0_val - r_major_sq - r_minor_sq)
                    + 8.0 * r_major_sq * t0_val * t1_val;
                let a0 = s0_val * s0_val - 2.0 * s0_val * (r_major_sq + r_minor_sq)
                    + c_sq_sq
                    + 4.0 * r_major_sq * t0_val * t0_val;

                let v_roots = crate::solve_quartic(1.0, a3, a2, a1, a0);

                v_roots
                    .iter()
                    .filter(|v| v.is_finite())
                    .map(|&v| (v, cone.point_at(u_mid, v)))
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

    // Dedup near-duplicate trailing points (closed curve degeneracy)
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
