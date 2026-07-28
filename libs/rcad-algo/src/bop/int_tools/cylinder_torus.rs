//! Analytic intersection of a cylinder with a torus.
//!
//! # Cases
//!
//! - **Coaxial**: Two circles when cylinder radius intersects torus tube
//! - **Skew** (non-coaxial, non-parallel): Quartic curve solved analytically
//!   via Ferrari's method at each u azimuth
//! - **General**: Numerical marching fallback

use std::f64::consts::TAU;

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Circle3, CylindricalSurface, ToroidalSurface, any_perpendicular};



/// Result of cylinder 脳 torus intersection.
#[derive(Debug, Clone)]
pub enum CylinderTorusResult {
    /// No intersection.
    NoIntersection,
    /// Single tangent circle.
    TangentCircle(Circle3),
    /// Two circles (coaxial case).
    TwoCircles(Circle3, Circle3),
    /// Skew/quartic intersection: one or more polyline branches sampled on the
    /// cylinder parameterization.  For each cylinder azimuth u 锟?[0, 2蟺), solve
    /// the quartic equation in v derived from the torus implicit equation.
    SkewQuartic(Vec<Vec<DVec3>>),
    /// Complex intersection, fall back to numerical marching.
    General,
}

/// Compute the analytic intersection of `cylinder` and `torus`.
pub fn intersect_cylinder_torus(
    cylinder: &CylindricalSurface,
    torus: &ToroidalSurface,
) -> CylinderTorusResult {
    intersect_cylinder_torus_with_tolerance(cylinder, torus, 0.0)
}

/// Cylinder 脳 torus intersection with fuzzy tolerance.
pub fn intersect_cylinder_torus_with_tolerance(
    cylinder: &CylindricalSurface,
    torus: &ToroidalSurface,
    fuzzy_tol: f64,
) -> CylinderTorusResult {
    let tol = rcad_kernel::rcad_kernel::CONFUSION + fuzzy_tol.max(0.0);

    let a_cyl = cylinder.axis.normalize();
    let a_tor = torus.axis.normalize();

    // Check for coaxial case
    let cross = a_cyl.cross(a_tor);
    let sin_angle = cross.length();

    // Project cylinder origin onto torus axis
    let delta = cylinder.origin - torus.center;
    let d_perp = (delta - a_tor * delta.dot(a_tor)).length();

    // Coaxial: same axis line
    if sin_angle < rcad_kernel::ANGULAR && d_perp < tol {
        return intersect_cylinder_torus_coaxial(cylinder, torus, tol);
    }

    // Skew: try quartic solver
    let skew_result = intersect_skew_cylinder_torus(cylinder, torus);
    if !skew_result.is_empty() {
        return CylinderTorusResult::SkewQuartic(skew_result);
    }

    // General case: numerical fallback
    CylinderTorusResult::General
}

/// Solve a鈧劼穢锟?+ a鈧兟穢鲁 + a鈧偮穢虏 + a鈧伮穢 + a鈧€ = 0, return sorted real roots.
fn solve_quartic_real(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    if a4.abs() < 1e-15 {
        return solve_cubic_real(a3, a2, a1, a0);
    }

    // Normalize to monic: x锟?+ b路x鲁 + c路x虏 + d路x + e = 0
    let (b, c, d, e) = (a3 / a4, a2 / a4, a1 / a4, a0 / a4);

    // Depress: x = y - b/4 锟?y锟?+ p路y虏 + q路y + r = 0
    let bb = b * b;
    let p = c - 3.0 * bb / 8.0;
    let q = d - b * c / 2.0 + b * bb / 8.0;
    let r = e - b * d / 4.0 + bb * c / 16.0 - 3.0 * bb * bb / 256.0;

    // Handle special: depressed quartic is biquadratic (q 锟?0)
    if q.abs() < 1e-14 {
        let mut roots = Vec::new();
        let y2_roots = solve_quadratic_real(1.0, p, r);
        for &y2 in &y2_roots {
            if y2 >= 0.0 {
                roots.push(y2.sqrt());
                if y2 > 0.0 {
                    roots.push(-y2.sqrt());
                }
            }
        }
        for x in &mut roots {
            *x -= b / 4.0;
        }
        roots.sort_by(|a, bb| a.partial_cmp(bb).unwrap());
        return roots;
    }

    // Resolvent cubic: m鲁 + 2p路m虏 + (p虏 - 4r)路m - q虏 = 0
    let m_roots = solve_cubic_real(1.0, 2.0 * p, p * p - 4.0 * r, -q * q);
    let m = m_roots
        .into_iter()
        .find(|&m| m >= 0.0 && m.is_finite())
        .unwrap_or_else(|| {
            // Try Newton to find a non-negative root
            newton_root(1.0, 2.0 * p, p * p - 4.0 * r, -q * q)
        });

    // If all resolvent roots are negative, we can still proceed with
    // complex arithmetic, but for practical CAD purposes this means the
    // quartic has 0 or 4 complex roots 锟?no real solutions.
    if m < 0.0 {
        return Vec::new();
    }

    let sqrt_m = m.sqrt();

    // Factor: y锟?+ p路y虏 + q路y + r = (y虏 + A路y + B)(y虏 - A路y + C)
    // where A = sqrt_m, B + C = p + m, C - B = q / sqrt_m
    let (b_val, c_val) = if sqrt_m > 1e-15 {
        let half_q_over_sqrt_m = 0.5 * q / sqrt_m;
        (
            0.5 * (p + m) - half_q_over_sqrt_m,
            0.5 * (p + m) + half_q_over_sqrt_m,
        )
    } else {
        (0.5 * (p + m), 0.5 * (p + m))
    };

    let mut roots = Vec::new();

    // y虏 + sqrt_m路y + b_val = 0
    let d1 = sqrt_m * sqrt_m - 4.0 * b_val;
    if d1 >= 0.0 {
        let sd = d1.sqrt();
        roots.push((-sqrt_m + sd) / 2.0);
        roots.push((-sqrt_m - sd) / 2.0);
    }

    // y虏 - sqrt_m路y + c_val = 0
    let d2 = sqrt_m * sqrt_m - 4.0 * c_val;
    if d2 >= 0.0 {
        let sd = d2.sqrt();
        roots.push((sqrt_m + sd) / 2.0);
        roots.push((sqrt_m - sd) / 2.0);
    }

    // De-depress
    for x in &mut roots {
        *x -= b / 4.0;
    }

    roots.sort_by(|a, bb| a.partial_cmp(bb).unwrap());
    roots
}

/// Solve a鈧兟穢鲁 + a鈧偮穢虏 + a鈧伮穢 + a鈧€ = 0, return sorted real roots.
fn solve_cubic_real(a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    if a3.abs() < 1e-15 {
        return solve_quadratic_real(a2, a1, a0);
    }

    // Normalize: x鲁 + p路x虏 + q路x + r = 0
    let (p, q, r) = (a2 / a3, a1 / a3, a0 / a3);

    // Depress: x = t - p/3 锟?t鲁 + a路t + b = 0
    let pp = p * p;
    let a = q - pp / 3.0;
    let b = 2.0 * p * pp / 27.0 - p * q / 3.0 + r;

    if a.abs() < 1e-15 {
        // t鲁 = -b
        return vec![b.signum() * (-b).abs().powf(1.0 / 3.0) - p / 3.0];
    }

    let disc = b * b / 4.0 + a * a * a / 27.0;

    if disc > 0.0 {
        // One real root (Cardano)
        let sqrt_d = disc.sqrt();
        let u = (-b / 2.0 + sqrt_d).cbrt();
        let v = (-b / 2.0 - sqrt_d).cbrt();
        vec![u + v - p / 3.0]
    } else if disc.abs() < 1e-15 {
        // Multiple real roots (discriminant 锟?0)
        let u = (-b / 2.0).cbrt();
        let t1 = 2.0 * u;
        let t2 = -u;
        let mut roots = vec![t1 - p / 3.0, t2 - p / 3.0];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    } else {
        // Three real roots (trigonometric / Vieta)
        let r_val = (-a / 3.0).sqrt();
        let cos_phi = (-b / (2.0 * r_val * r_val * r_val)).clamp(-1.0, 1.0);
        let phi = cos_phi.acos();
        let two_r = 2.0 * r_val;
        let pi = std::f64::consts::PI;
        let t1 = two_r * (phi / 3.0).cos();
        let t2 = two_r * ((phi + 2.0 * pi) / 3.0).cos();
        let t3 = two_r * ((phi + 4.0 * pi) / 3.0).cos();
        let mut roots = vec![t1 - p / 3.0, t2 - p / 3.0, t3 - p / 3.0];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    }
}

/// Solve a鈧偮穢虏 + a鈧伮穢 + a鈧€ = 0, return sorted real roots.
fn solve_quadratic_real(a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    if a2.abs() < 1e-15 {
        if a1.abs() < 1e-15 {
            return Vec::new();
        }
        return vec![-a0 / a1];
    }

    let disc = a1 * a1 - 4.0 * a2 * a0;
    if disc < 0.0 {
        return Vec::new();
    }
    if disc.abs() < 1e-15 {
        return vec![-a1 / (2.0 * a2)];
    }
    let sd = disc.sqrt();
    let x1 = (-a1 - sd) / (2.0 * a2);
    let x2 = (-a1 + sd) / (2.0 * a2);
    vec![x1.min(x2), x1.max(x2)]
}

/// Newton-iteration fallback to find one real root of a鈧兟穢鲁 + a鈧偮穢虏 + a鈧伮穢 + a鈧€.
fn newton_root(a3: f64, a2: f64, a1: f64, a0: f64) -> f64 {
    let f = |x: f64| -> (f64, f64) {
        let fx = ((a3 * x + a2) * x + a1) * x + a0;
        let dfx = (3.0 * a3 * x + 2.0 * a2) * x + a1;
        (fx, dfx)
    };

    for &start in &[0.0, 1.0, -1.0, 2.0, -2.0] {
        let mut x = start;
        for _ in 0..32 {
            let (fx, dfx) = f(x);
            if dfx.abs() < 1e-15 {
                break;
            }
            let step = fx / dfx;
            x -= step;
            if step.abs() < 1e-14 {
                return x;
            }
        }
        if x.is_finite() {
            return x;
        }
    }
    0.0
}

/// Skew (non-coaxial, non-parallel) cylinder 脳 torus intersection via quartic
/// solver.  Parameterize the cylinder as P(u,v), substitute into the torus
/// implicit equation (x虏+y虏+z虏-R虏-r虏)虏 = 4R虏(r虏-z虏), obtaining a quartic in v
/// at each u.  Solve via Ferrari's method, extract branches.
fn intersect_skew_cylinder_torus(
    cyl: &CylindricalSurface,
    torus: &ToroidalSurface,
) -> Vec<Vec<DVec3>> {
    let a = cyl.axis.normalize();
    let a_tor = torus.axis.normalize();
    let o = cyl.origin - torus.center; // cylinder origin in torus frame
    let rho = cyl.radius;
    let r_major = torus.major_radius;
    let r_minor = torus.minor_radius;

    let r_major_sq = r_major * r_major;
    let r_minor_sq = r_minor * r_minor;
    let r_sum_sq = r_major_sq + r_minor_sq;
    let r_diff_sq = r_major_sq - r_minor_sq; // R虏 - r虏

    // Perpendicular basis for cylinder radial direction
    let x_dir = any_perpendicular(a);
    let y_dir = a.cross(x_dir).normalize();

    // Constants for the quartic coefficients
    let o_sq = o.length_squared();
    let d1 = a.dot(a_tor); // a路a_tor 锟?constant
    let d1_sq = d1 * d1;
    let o_dot_a = o.dot(a);
    let o_dot_a_tor = o.dot(a_tor);

    const N_SAMPLES: usize = 128;
    // Store (u, v, point) tuples for proximity-based branch extraction
    let mut samples: Vec<(f64, f64, DVec3)> = Vec::new();

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cu, su) = (u.cos(), u.sin());

        // Radial direction on cylinder at this u
        let r_dir = cu * x_dir + su * y_dir;

        // C锟?u) = a路O + 蟻路(a路r_dir)
        let c1 = o_dot_a + rho * a.dot(r_dir);

        // C鈧€(u) = |O|虏 + 蟻虏 + 2蟻路O路r_dir
        let c0 = o_sq + rho * rho + 2.0 * rho * o.dot(r_dir);

        // D锟?= a路a_tor (constant)
        // D鈧€(u) = O路a_tor + 蟻路(r_dir路a_tor)
        let d0 = o_dot_a_tor + rho * r_dir.dot(a_tor);

        // Quartic: v锟?+ a鈧兟穠鲁 + a鈧偮穠虏 + a鈧伮穠 + a鈧€ = 0
        let a3 = 4.0 * c1;
        let a2 = 4.0 * c1 * c1 + 2.0 * c0 - 2.0 * r_sum_sq + 4.0 * r_major_sq * d1_sq;
        let a1 = 4.0 * c1 * c0 - 4.0 * r_sum_sq * c1 + 8.0 * r_major_sq * d1 * d0;
        let a0 = c0 * c0 - 2.0 * r_sum_sq * c0 + r_diff_sq * r_diff_sq + 4.0 * r_major_sq * d0 * d0;

        let v_roots = solve_quartic_real(1.0, a3, a2, a1, a0);

        for &v in &v_roots {
            if v.is_finite() {
                let p = cyl.point_at(u, v);
                if p.is_finite() {
                    samples.push((u, v, p));
                }
            }
        }
    }

    if samples.len() < 2 {
        return vec![];
    }

    // Extract branches by proximity: group (u, v) pairs into contiguous curves.
    // Sort by u, then v.
    samples.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    });

    // Clustering: at each distinct u, we have a set of v values.
    // Group them into branches by proximity across consecutive u values.
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
            // Determine v-matching threshold based on typical spacing between
            // adjacent roots at this u.  v-values can span orders-of-magnitude
            // more than u-values, so a u-based threshold is meaningless.
            let v_threshold = if cur_pts.len() >= 2 {
                let v_range = cur_pts.last().unwrap().0 - cur_pts[0].0;
                (v_range / cur_pts.len() as f64) * 2.5 // 2.5脳 typical gap
            } else {
                (TAU / N_SAMPLES as f64) * 20.0 // fallback: generous bound
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

    // 鈹€鈹€ Adaptive chord-error refinement 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
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
                let r_dir = cu * x_dir + su * y_dir;
                let c1 = o_dot_a + rho * a.dot(r_dir);
                let c0 = o_sq + rho * rho + 2.0 * rho * o.dot(r_dir);
                let d0 = o_dot_a_tor + rho * r_dir.dot(a_tor);
                let a3 = 4.0 * c1;
                let a2 = 4.0 * c1 * c1 + 2.0 * c0 - 2.0 * r_sum_sq + 4.0 * r_major_sq * d1_sq;
                let a1 = 4.0 * c1 * c0 - 4.0 * r_sum_sq * c1 + 8.0 * r_major_sq * d1 * d0;
                let a0 = c0 * c0 - 2.0 * r_sum_sq * c0
                    + r_diff_sq * r_diff_sq
                    + 4.0 * r_major_sq * d0 * d0;

                let v_roots = solve_quartic_real(1.0, a3, a2, a1, a0);

                v_roots
                    .iter()
                    .filter(|v| v.is_finite())
                    .map(|&v| (v, cyl.point_at(u_mid, v)))
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
            if (branch[n - 1] - branch[0]).length_squared() < 1e-24 {
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

fn intersect_cylinder_torus_coaxial(
    cylinder: &CylindricalSurface,
    torus: &ToroidalSurface,
    tol: f64,
) -> CylinderTorusResult {
    // In the (rho, z) half-plane:
    // Torus tube: (rho - major_r)^2 + z^2 = minor_r^2
    // Cylinder: rho = r_cyl

    let major_r = torus.major_radius;
    let minor_r = torus.minor_radius;
    let r_cyl = cylinder.radius;

    // Cylinder radius must intersect the tube circle
    // Tube circle center at (major_r, 0) in (rho, z) plane
    // Distance from tube center to cylinder: |major_r - r_cyl|
    let d = (major_r - r_cyl).abs();

    if d > minor_r + tol {
        return CylinderTorusResult::NoIntersection;
    }

    if (d - minor_r).abs() < tol {
        // Tangent: one circle
        let z = 0.0;
        let center = torus.center + torus.axis * z;
        return CylinderTorusResult::TangentCircle(Circle3::new(center, torus.axis, r_cyl));
    }

    // Two intersection points in (rho, z): z = +/-sqrt(minor_r^2 - d^2)
    let z_offset = (minor_r * minor_r - d * d).sqrt();

    let center1 = torus.center + torus.axis * z_offset;
    let center2 = torus.center - torus.axis * z_offset;

    CylinderTorusResult::TwoCircles(
        Circle3::new(center1, torus.axis, r_cyl),
        Circle3::new(center2, torus.axis, r_cyl),
    )
}
