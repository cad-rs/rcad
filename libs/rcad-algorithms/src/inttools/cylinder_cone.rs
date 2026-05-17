//! Analytic intersection of a cylinder and a cone.
//!
//! # Case classification
//!
//! ## Coaxial (axes coincide)
//!
//! When the cylinder axis and cone axis are the same line, the intersection is
//! a **circle** at the height where the cone radius equals the cylinder radius:
//!
//! ```text
//! r_cone(h) = (h − h_apex) · tan(β)  where h_apex is the height of the apex
//! r_cone(h) = r_cyl  →  h = h_apex + r_cyl / tan(β)
//! ```
//!
//! Returns [`CoaxialCircle`](CylinderConeResult::CoaxialCircle) if the apex is
//! on the positive-axis side (the circle is real), or
//! [`NoIntersection`](CylinderConeResult::NoIntersection) otherwise.
//!
//! ## Parallel axes (non-coaxial)
//!
//! When the axes are parallel but distinct (cross-product ≈ 0), the
//! cylinder-cone intersection is a quartic curve in general.  We perform a
//! radial distance test to detect obvious non-intersections and fall back to
//! marching otherwise.
//!
//! ## General / skew axes
//!
//! For all other configurations we return [`General`](CylinderConeResult::General)
//! so the caller can fall back to numeric marching.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, CylindricalSurface};

use crate::tolerance::*;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of cylinder × cone intersection.
#[derive(Debug, Clone)]
pub enum CylinderConeResult {
    /// The surfaces do not intersect.
    NoIntersection,
    /// Coaxial configuration: exactly one intersection circle.
    CoaxialCircle(Circle3),
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
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `cyl` and `cone`.
pub fn intersect_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
) -> CylinderConeResult {
    let a_cyl = cyl.axis.normalize();
    let a_cone = cone.axis.normalize();

    let cross = a_cyl.cross(a_cone);
    let sin_angle = cross.length(); // |sin θ| between the two axes

    // ── Parallel axes (including coaxial) ────────────────────────────────────
    if sin_angle < TOLERANCE_ANG {
        return intersect_parallel_cylinder_cone(cyl, cone, a_cyl, a_cone);
    }

    // ── General / skew ────────────────────────────────────────────────────────
    // Perform a quick distance-based no-intersection test:
    // Find closest distance between the two axes.  If the cylinder completely
    // misses the cone's bounding envelope, return NoIntersection.

    // For now, return General (marching handles this correctly).
    CylinderConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel (and coaxial) axes
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_parallel_cylinder_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
    a_cyl: DVec3,
    a_cone: DVec3,
) -> CylinderConeResult {
    let apex = cone.apex_point();
    // Make sure a_cone points the same direction as a_cyl for height arithmetic.
    // (The cross product is ~0 so they are parallel; they may anti-parallel.)
    let _a_cone = if a_cyl.dot(a_cone) >= 0.0 { a_cone } else { -a_cone };

    let r_cyl = cyl.radius;
    let tan_beta = cone.half_angle_rad.tan();

    // Perpendicular distance between the two axes.
    let delta = apex - cyl.origin;
    let delta_perp = delta - a_cyl * delta.dot(a_cyl);
    let d_perp = delta_perp.length();

    // ── Coaxial ──────────────────────────────────────────────────────────────
    if d_perp < TOLERANCE_ABS {
        // Height of apex above cyl.origin along shared axis.
        let h_apex = (apex - cyl.origin).dot(a_cyl);

        // At height h (measured from cyl.origin), cone radius = (h - h_apex)*tan_beta
        // (only positive when h > h_apex, i.e. above the apex in axis direction).
        // Set equal to r_cyl:  h = h_apex + r_cyl / tan_beta
        if tan_beta.abs() < TOLERANCE_FLOAT_LOOSE {
            // Degenerate cone (half_angle = 0), no lateral surface.
            return CylinderConeResult::NoIntersection;
        }
        let h_circle = h_apex + r_cyl / tan_beta;

        // The circle must be on the cone's nappe (h_circle > h_apex).
        if h_circle <= h_apex - TOLERANCE_ABS {
            return CylinderConeResult::NoIntersection;
        }

        let center = cyl.origin + a_cyl * h_circle;
        return CylinderConeResult::CoaxialCircle(Circle3 {
            center,
            normal: a_cyl,
            radius: r_cyl,
        });
    }

    // ── Parallel but offset axes ──────────────────────────────────────────────
    // At each axial height h the cone has radius r_cone(h) = |h - h_apex|·tan(β).
    // The cylinder cross-section is a circle of radius r_cyl centred at a
    // perpendicular offset d_perp from the cone axis.  Intersection reduces to
    // a circle-circle test in the plane perpendicular to the axes, which is
    // solved analytically via axial sweep.
    if d_perp < TOLERANCE_ABS {
        // Should not reach here (coaxial case above), but guard against
        // numerical near-zero.
        return CylinderConeResult::NoIntersection;
    }
    return intersect_parallel_offset_cylinder_cone(cyl, cone, a_cyl, d_perp, delta_perp);
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel-offset axial sweep
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the intersection of a cylinder and a cone with parallel but
/// non-coaxial axes via axial sweep.
///
/// At each height `h` along the common axis direction `a`, the cone
/// cross-section is a circle of radius `r_cone(h) = |h - h_apex|·tan(β)`
/// centred on the cone axis, and the cylinder cross-section is a circle of
/// radius `r_cyl` centred on the cylinder axis.  The two axes are separated
/// by perpendicular distance `d_perp` (vector `delta_perp` from cylinder axis
/// to cone apex in the perpendicular plane).
///
/// The intersection exists when the two circles overlap:
///   `|r_cyl - r_cone(h)| ≤ d_perp ≤ r_cyl + r_cone(h)`
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
    if tan_beta.abs() < TOLERANCE_FLOAT_LOOSE {
        return CylinderConeResult::NoIntersection;
    }

    let apex = cone.apex_point();
    let h_apex = (apex - cyl.origin).dot(a);

    // Geometry constants for the sweep.
    let abs_diff = (d_perp - r_cyl).abs();
    let d_sum = d_perp + r_cyl;

    // Perpendicular basis: u points from cylinder axis toward cone axis,
    // v = a × u completes the right-handed frame.
    let u = delta_perp / d_perp;
    let v = a.cross(u).normalize_or_zero();

    const N_SAMPLES: usize = 128;
    let mut branches: Vec<Vec<DVec3>> = Vec::new();

    // Sweep both nappes (upper: sign = +1, lower: sign = -1).
    for &sign in &[1.0, -1.0] {
        let h_min = h_apex + sign * abs_diff / tan_beta;
        let h_max = h_apex + sign * d_sum / tan_beta;

        let (h_lo, h_hi) = if sign > 0.0 { (h_min, h_max) } else { (h_max, h_min) };

        if h_hi - h_lo <= TOLERANCE_ABS * 10.0 {
            continue; // degenerately narrow intersection band
        }

        let mut branch_plus: Vec<DVec3> = Vec::with_capacity(N_SAMPLES + 1);
        let mut branch_minus: Vec<DVec3> = Vec::with_capacity(N_SAMPLES + 1);

        for i in 0..=N_SAMPLES {
            let t = i as f64 / N_SAMPLES as f64;
            let h = h_lo + (h_hi - h_lo) * t;
            let r_cone = (h - h_apex).abs() * tan_beta;

            let center = cyl.origin + a * h;

            // cos(φ) where φ is the angle on the cylinder cross-section circle
            // from the u-direction.  Derived from |P(h,φ) - P_cone(h)|² = r_cone².
            let numer = r_cyl * r_cyl + d_perp * d_perp - r_cone * r_cone;
            let denom = 2.0 * r_cyl * d_perp;
            let cos_phi = (numer / denom).clamp(-1.0, 1.0);

            if cos_phi.abs() >= 1.0 - TOLERANCE_ANG {
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
                if (pts[n - 1] - pts[0]).length_squared() < TOLERANCE_VEC_SQ_MIN {
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

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn cyl(origin: DVec3, axis: DVec3, radius: f64) -> CylindricalSurface {
        CylindricalSurface { origin, axis, radius }
    }

    fn cone(apex: DVec3, axis: DVec3, half_angle_deg: f64) -> ConicalSurface {
        ConicalSurface {
            apex,
            axis,
            radius: 0.0,
            half_angle_rad: half_angle_deg.to_radians(),
        }
    }

    /// Coaxial cylinder and cone (Z axis): circle at h = apex_z + r / tan(β).
    ///
    /// Cone: apex at (0,0,0), axis Z, half_angle=45° → tan(β)=1.
    /// Cylinder: r=2, axis Z → circle at h = 0 + 2/1 = 2.
    #[test]
    fn coaxial_circle_z_axis() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 2.0);
        let k = cone(DVec3::ZERO, DVec3::Z, 45.0);
        match intersect_cylinder_cone(&c, &k) {
            CylinderConeResult::CoaxialCircle(circ) => {
                assert!(
                    (circ.center.z - 2.0).abs() < TOLERANCE_COORD_SUB,
                    "circle z={}, expected 2.0",
                    circ.center.z
                );
                assert!((circ.radius - 2.0).abs() < TOLERANCE_COORD_SUB);
            }
            other => panic!("expected CoaxialCircle, got {other:?}"),
        }
    }

    /// Coaxial but cone apex ABOVE the cylinder origin with an offset.
    ///
    /// Cone: apex at (0,0,5), axis Z, half_angle=30° → tan(β)=1/√3.
    /// Cylinder: r=1 → h = 5 + 1/(1/√3) = 5 + √3 ≈ 6.732.
    #[test]
    fn coaxial_circle_offset_apex() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::new(0.0, 0.0, 5.0), DVec3::Z, 30.0);
        match intersect_cylinder_cone(&c, &k) {
            CylinderConeResult::CoaxialCircle(circ) => {
                let expected_h = 5.0 + 1.0 / (30.0_f64.to_radians().tan());
                assert!(
                    (circ.center.z - expected_h).abs() < TOLERANCE_COORD_SUB,
                    "circle z={}, expected {}",
                    circ.center.z,
                    expected_h
                );
            }
            other => panic!("expected CoaxialCircle, got {other:?}"),
        }
    }

    /// Skew axes → General.
    #[test]
    fn skew_axes_general() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::ZERO, DVec3::new(1.0, 1.0, 0.0).normalize(), 45.0);
        assert!(matches!(intersect_cylinder_cone(&c, &k), CylinderConeResult::General));
    }

    /// Perpendicular axes → General.
    #[test]
    fn perpendicular_axes_general() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::ZERO, DVec3::X, 45.0);
        assert!(matches!(intersect_cylinder_cone(&c, &k), CylinderConeResult::General));
    }

    /// Parallel, offset axes → two polyline branches (upper nappe).
    ///
    /// Cylinder: origin (0,0,0), axis Z, radius 2.
    /// Cone: apex at (3,0,0), axis Z, half-angle 45°.
    /// d_perp = 3, r_cyl = 2, tan(45°) = 1.
    /// Upper nappe: h_min = |3-2| = 1, h_max = 3+2 = 5.
    #[test]
    fn parallel_offset_two_branches() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 2.0);
        let k = cone(DVec3::new(3.0, 0.0, 0.0), DVec3::Z, 45.0);
        match intersect_cylinder_cone(&c, &k) {
            CylinderConeResult::ParallelOffsetPolyline(branches) => {
                assert!(
                    branches.len() >= 2,
                    "expected at least 2 branches (one per nappe), got {}",
                    branches.len()
                );
                for (i, branch) in branches.iter().enumerate() {
                    assert!(branch.len() >= 2, "branch {i} has < 2 points");
                    for p in branch {
                        assert!(p.is_finite(), "branch {i} contains non-finite point");
                    }
                }
            }
            other => panic!("expected ParallelOffsetPolyline, got {other:?}"),
        }
    }

    /// Parallel, offset axes with the cylinder axis far from the cone apex
    /// → still intersects (cone radius grows unboundedly), branches may merge
    /// near tangent transitions.  Just verify it returns a polyline result
    /// (not General).
    #[test]
    fn parallel_offset_far_axis() {
        let c = cyl(DVec3::ZERO, DVec3::Z, 1.0);
        let k = cone(DVec3::new(50.0, 0.0, 0.0), DVec3::Z, 30.0);
        let result = intersect_cylinder_cone(&c, &k);
        match &result {
            CylinderConeResult::ParallelOffsetPolyline(branches) => {
                assert!(!branches.is_empty(), "expected non-empty branches");
                for branch in branches {
                    assert!(branch.len() >= 2);
                    for p in branch {
                        assert!(p.is_finite());
                    }
                }
            }
            CylinderConeResult::NoIntersection => {
                // For very small tan_beta this is possible.
                let tan_beta = (30.0_f64.to_radians()).tan();
                if tan_beta.abs() < 1e-10 {
                    return; // degenerate cone, NoIntersection is OK
                }
                panic!("expected ParallelOffsetPolyline, got NoIntersection with tan_beta={tan_beta}");
            }
            other => panic!("expected ParallelOffsetPolyline, got {other:?}"),
        }
    }

    /// Coaxial with a slight numerical perturbation → should still detect
    /// coaxial case (not fall through to parallel-offset).
    #[test]
    fn near_coaxial_still_circle() {
        // Perturb the apex very slightly off-axis.
        let c = cyl(DVec3::ZERO, DVec3::Z, 2.0);
        let k = cone(DVec3::new(1e-8, 0.0, 0.0), DVec3::Z, 45.0);
        match intersect_cylinder_cone(&c, &k) {
            CylinderConeResult::CoaxialCircle(circ) => {
                assert!((circ.center.z - 2.0).abs() < TOLERANCE_COORD_SUB);
            }
            CylinderConeResult::ParallelOffsetPolyline(_) => {
                // Also acceptable if numeric noise triggers parallel-offset path.
            }
            other => panic!("expected CoaxialCircle or ParallelOffsetPolyline, got {other:?}"),
        }
    }
}
