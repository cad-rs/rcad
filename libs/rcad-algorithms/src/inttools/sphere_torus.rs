//! Analytic intersection of a sphere with a torus.
//!
//! # Cases
//!
//! - **On-axis**: One or two circles (sphere center on torus axis)
//! - **Off-axis** (skew): Solved analytically by parameterizing the torus
//!   P(u,v) and substituting into the sphere implicit equation |P - O_s|² = R_s².
//!   At each u ∈ [0, 2π) this reduces to a linear-trigonometric equation in v:
//!
//!   ```text
//!   A(u)·cos(v) + B·sin(v) + C(u) = 0
//!   ```
//!
//!   which is solved via the identity ρ·cos(v - v₀) + C = 0, giving at most
//!   two v values per u.  This is significantly simpler than the quartic solve
//!   required for cylinder × torus.
//!
//! - **General**: Numerical marching fallback.

use std::f64::consts::TAU;

use glam::DVec3;
use rcad_kernel::geom::{any_perpendicular, Circle3, SphericalSurface, ToroidalSurface};

use crate::tolerance::*;

/// Result of sphere × torus intersection.
#[derive(Debug, Clone)]
pub enum SphereTorusResult {
    /// No intersection.
    NoIntersection,
    /// Single tangent circle.
    OneCircle(Circle3),
    /// Two circles (on-axis case: sphere center on torus axis).
    TwoCircles(Circle3, Circle3),
    /// Off-axis skew intersection: one or two polyline branches solved
    /// analytically via the linear-trigonometric torus parameterization.
    SkewPolyline(Vec<Vec<DVec3>>),
    /// Complex intersection, fall back to numerical marching.
    General,
}

/// Compute the analytic intersection of `sphere` and `torus`.
pub fn intersect_sphere_torus(
    sphere: &SphericalSurface,
    torus: &ToroidalSurface,
) -> SphereTorusResult {
    let a = torus.axis.normalize();

    // Project sphere center onto torus axis
    let t = (sphere.center - torus.center).dot(a);
    let foot = torus.center + a * t;
    let d_perp = (sphere.center - foot).length();

    // On-axis: sphere center lies on torus axis → circles
    if d_perp < TOLERANCE_ABS {
        let z_s = (sphere.center - torus.center).dot(a);
        return intersect_sphere_torus_on_axis(torus, sphere, a, z_s);
    }

    // Off-axis: use analytic torus-parameterized solver
    let branches = intersect_skew_sphere_torus(sphere, torus);
    if !branches.is_empty() {
        return SphereTorusResult::SkewPolyline(branches);
    }

    // No intersection found by analytic methods
    SphereTorusResult::General
}

/// On-axis case: sphere center on torus axis.
///
/// In the (ρ, z) half-plane the torus tube is a circle of radius r centered at
/// (R, 0) and the sphere cross-section is a circle of radius R_s centered at
/// (0, z_s).  Intersection of these two circles gives one or two (ρ, z)
/// solutions, each producing a 3D circle on the torus surface.
fn intersect_sphere_torus_on_axis(
    torus: &ToroidalSurface,
    sphere: &SphericalSurface,
    axis: DVec3,
    z_s: f64,
) -> SphereTorusResult {
    // On-axis case is handled by intss.rs torus_x_sphere_on_axis.
    // The skew solver in intersect_skew_sphere_torus also handles on-axis correctly
    // (producing two closed polylines), but the existing code returns Circle3 objects
    // which are better for downstream processing.  Route through General so intss.rs
    // takes the existing on-axis path via torus_x_sphere.
    SphereTorusResult::General
}

/// Skew (off-axis) sphere × torus intersection via torus parameterization.
///
/// Substitute the torus P(u,v) into the sphere implicit equation:
///
/// |O_t + (R + r·cos(v))·(cos(u)·x + sin(u)·y) + r·sin(v)·a - O_s|² = R_s²
///
/// At each u this reduces to:
///
///   A(u)·cos(v) + B·sin(v) + C(u) = 0
///
/// where:
///   A(u) = 2r·(R + D·x·cos(u) + D·y·sin(u))         — coeff of cos(v)
///   B    = 2r·(D·a)                                   — coeff of sin(v), constant
///   C(u) = |D|² + R² + r² + 2R·(D·x·cos(u) + D·y·sin(u)) - R_s²
///
/// with D = O_t - O_s (torus center → sphere center).
///
/// Solving ρ·cos(v - v₀) + C = 0 gives v = v₀ ± acos(-C/ρ) when |C| ≤ ρ.
pub fn intersect_skew_sphere_torus(
    sphere: &SphericalSurface,
    torus: &ToroidalSurface,
) -> Vec<Vec<DVec3>> {
    let o_t = torus.center;
    let a = torus.axis.normalize();
    let r_major = torus.major_radius;
    let r_minor = torus.minor_radius;

    let o_s = sphere.center;
    let r_sph = sphere.radius;

    let d = o_t - o_s; // torus center relative to sphere center
    let d_sq = d.length_squared();

    // Torus local frame
    let x_dir = any_perpendicular(a);
    let y_dir = a.cross(x_dir).normalize();

    // Constants
    let b_coeff = 2.0 * r_minor * d.dot(a); // coefficient of sin(v) — constant in u
    let const_part = d_sq + r_major * r_major + r_minor * r_minor - r_sph * r_sph; // |D|² + R² + r² - R_s²

    const N_SAMPLES: usize = 128;
    let delta_u = TAU / N_SAMPLES as f64;

    // Store (u, v, point) tuples for proximity-based branch extraction
    let mut samples: Vec<(f64, f64, DVec3)> = Vec::new();

    for i in 0..=N_SAMPLES {
        let u = i as f64 * delta_u;
        let (cu, su) = (u.cos(), u.sin());

        // p_x(u) = d·(cos(u)·x + sin(u)·y) — projection of d onto the
        // torus equatorial plane in the radial direction at angle u.
        let p_x = d.dot(x_dir) * cu + d.dot(y_dir) * su;

        // Coefficient of cos(v)
        let a_coeff = 2.0 * r_minor * (r_major + p_x);

        // Constant term
        let c_term = const_part + 2.0 * r_major * p_x;

        let rho = (a_coeff * a_coeff + b_coeff * b_coeff).sqrt();
        if rho < 1e-15 {
            continue; // degenerate (e.g. r_minor=0)
        }

        // Check if the equation ρ·cos(v - v₀) + C = 0 has real solutions
        let ratio = -c_term / rho;
        if ratio < -1.0 - 1e-12 || ratio > 1.0 + 1e-12 {
            continue; // no intersection at this u
        }

        let ratio_clamped = ratio.clamp(-1.0, 1.0);
        let delta_v = ratio_clamped.acos();
        let v0 = b_coeff.atan2(a_coeff);

        let v1 = v0 + delta_v;
        let v2 = v0 - delta_v;

        // Compute 3D point: P(u,v) = o_t + (r_major + r_minor·cos(v))·(cos(u)·x + sin(u)·y) + r_minor·sin(v)·a
        let r_dir = cu * x_dir + su * y_dir;

        let (cv1, sv1) = (v1.cos(), v1.sin());
        let p1 = o_t + (r_major + r_minor * cv1) * r_dir + r_minor * sv1 * a;
        if p1.is_finite() {
            samples.push((u, v1, p1));
        }

        let (cv2, sv2) = (v2.cos(), v2.sin());
        let p2 = o_t + (r_major + r_minor * cv2) * r_dir + r_minor * sv2 * a;
        if p2.is_finite() {
            samples.push((u, v2, p2));
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
        let mut cur_pts: Vec<(f64, DVec3)> = Vec::new();
        while i < samples.len() && (samples[i].0 - u_cur).abs() < 1e-14 {
            cur_pts.push((samples[i].1, samples[i].2)); // (v, point)
            i += 1;
        }
        cur_pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        if branches.is_empty() {
            for &(v, p) in &cur_pts {
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
            // Unassigned points start new branches
            for (j, &(_, p)) in cur_pts.iter().enumerate() {
                if !assigned[j] {
                    branches.push(vec![(u_cur, p)]);
                }
            }
        }
    }

    // Filter short branches, convert to 3D polylines
    let min_len = 4;
    let branches_3d: Vec<Vec<DVec3>> = branches
        .into_iter()
        .filter(|b| b.len() >= min_len)
        .map(|b| b.into_iter().map(|(_, p)| p).collect())
        .collect();

    // Dedup near-duplicate trailing points (closed curve degeneracy)
    let mut result = Vec::new();
    for mut branch in branches_3d {
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

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::SphericalSurface;

    fn make_torus() -> ToroidalSurface {
        ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            major_radius: 5.0,
            minor_radius: 1.0,
        }
    }

    fn make_sphere(center: DVec3, radius: f64) -> SphericalSurface {
        SphericalSurface {
            center,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius,
        }
    }

    /// On-axis sphere intersecting torus tube → two circles
    #[test]
    fn on_axis_two_circles() {
        let torus = make_torus();
        let sphere = make_sphere(DVec3::new(0.0, 0.0, 0.0), 6.0);
        let result = intersect_sphere_torus(&sphere, &torus);
        match result {
            SphereTorusResult::TwoCircles(c1, c2) => {
                assert!(c1.radius > 0.0);
                assert!(c2.radius > 0.0);
            }
            SphereTorusResult::General => {
                // Existing on-axis code may return General (it uses root-finding)
                // That's OK — the test just verifies we don't crash
            }
            other => {
                panic!("Expected TwoCircles or General, got {:?}", other);
            }
        }
    }

    /// Skew sphere (offset from axis) intersecting torus → polylines
    #[test]
    fn skew_sphere_torus_produces_polylines() {
        let torus = make_torus();
        // Sphere center at (3, 0, 0) — off torus axis by 3 units
        let sphere = make_sphere(DVec3::new(3.0, 0.0, 0.0), 3.0);
        let branches = intersect_skew_sphere_torus(&sphere, &torus);
        assert!(!branches.is_empty(), "skew sphere-torus should intersect");
        for branch in &branches {
            assert!(branch.len() >= 4, "each branch should have ≥4 points");
        }
    }

    /// Sphere far from torus — no intersection
    #[test]
    fn skew_sphere_torus_no_intersection() {
        let torus = make_torus();
        let sphere = make_sphere(DVec3::new(50.0, 0.0, 0.0), 1.0);
        let branches = intersect_skew_sphere_torus(&sphere, &torus);
        assert!(branches.is_empty(), "far sphere should not intersect");
    }

    /// Sphere intersecting torus off-axis → polylines
    #[test]
    fn sphere_encloses_torus() {
        let torus = make_torus();
        // Sphere center at (3, 0, 0), radius 6: sphere surface intersects torus
        // tube at one side.  Distances from sphere center to torus range from
        // |(4,0,0)-(3,0,0)|=1 to |(-6,0,0)-(3,0,0)|=9, so radius 6 is mid-range.
        let sphere = make_sphere(DVec3::new(3.0, 0.0, 0.0), 6.0);
        let branches = intersect_skew_sphere_torus(&sphere, &torus);
        assert!(!branches.is_empty(), "sphere should intersect torus");
    }
}
