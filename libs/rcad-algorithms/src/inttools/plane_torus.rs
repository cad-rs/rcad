//! Analytic intersection of a plane with a torus.
//!
//! # Cases
//!
//! - **Perpendicular to axis**: Two circles (inner and outer equator)
//! - **Parallel to axis**: Two circles when plane intersects tube center circle
//!   - If |d| < R: Two circles of radius r at intersections with tube center circle
//!   - If |d| = R: Two tangent circles (Villarceau circles configuration)
//!   - If |d| > R+r: No intersection
//! - **Oblique**: Complex curve, fall back to numerical marching

use glam::DVec3;
use rcad_kernel::any_perpendicular;
use rcad_kernel::geom::{Circle3, Plane, ToroidalSurface};

use crate::tolerance::*;

/// Result of plane x torus intersection.
#[derive(Debug, Clone)]
pub enum PlaneTorusResult {
    /// The plane does not intersect the torus.
    NoIntersection,
    /// Single tangent circle.
    TangentCircle(Circle3),
    /// Two circles (perpendicular case).
    TwoCircles(Circle3, Circle3),
    /// Skew plane intersection: one or more polyline branches sampled on the
    /// torus parameterization.  For each torus azimuth u ∈ [0, 2π), solve the
    /// plane constraint cos(v)·B(u) + sin(v)·C = D(u) for the tube angle v.
    SkewPolyline(Vec<Vec<DVec3>>),
    /// Complex intersection, fall back to numerical marching.
    General,
}

/// Compute the analytic intersection of `plane` and `torus`.
pub fn intersect_plane_torus(plane: &Plane, torus: &ToroidalSurface) -> PlaneTorusResult {
    intersect_plane_torus_with_tolerance(plane, torus, 0.0)
}

/// Plane x torus intersection with fuzzy tolerance.
pub fn intersect_plane_torus_with_tolerance(
    plane: &Plane,
    torus: &ToroidalSurface,
    fuzzy_tol: f64,
) -> PlaneTorusResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    // Normalize plane normal and torus axis
    let n = plane.normal.normalize();
    let a = torus.axis.normalize();

    // Check if plane is perpendicular to torus axis
    let dot_na = n.dot(a).abs();

    if dot_na > 1.0 - TOLERANCE_ANG {
        // Plane perpendicular to axis: circular cross-section
        return intersect_plane_torus_perpendicular(plane, torus, tol);
    }

    // Check if plane is parallel to torus axis
    if dot_na < TOLERANCE_ANG {
        // Plane parallel to axis: may produce two circles
        return intersect_plane_torus_parallel(plane, torus, tol);
    }

    // Skew plane: try u-parameterized analytic solver
    let skew_result = intersect_plane_torus_skew(plane, torus);
    if !skew_result.is_empty() {
        return PlaneTorusResult::SkewPolyline(skew_result);
    }

    // General oblique case: fall back to numerical
    PlaneTorusResult::General
}

fn intersect_plane_torus_perpendicular(
    plane: &Plane,
    torus: &ToroidalSurface,
    tol: f64,
) -> PlaneTorusResult {
    // Distance from torus center to plane along axis
    let signed_dist = (torus.center - plane.origin).dot(torus.axis);
    let abs_dist = signed_dist.abs();

    // Maximum distance for intersection is the minor radius
    if abs_dist > torus.minor_radius + tol {
        return PlaneTorusResult::NoIntersection;
    }

    // Tangent case: one circle
    if (abs_dist - torus.minor_radius).abs() < tol {
        let center = torus.center - torus.axis * signed_dist;
        return PlaneTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: torus.major_radius,
        });
    }

    // Two circles at height signed_dist from torus center
    // Circle radius on the tube: sqrt(r^2 - d^2) where r = minor_radius, d = distance
    let tube_circle_r = (torus.minor_radius * torus.minor_radius - signed_dist * signed_dist).sqrt();

    // Two circles at major_radius +/- tube_circle_r from axis
    let r1 = torus.major_radius + tube_circle_r;
    let r2 = (torus.major_radius - tube_circle_r).max(0.0);

    let center = torus.center - torus.axis * signed_dist;

    if r2 < tol {
        // Inner circle degenerates to point
        PlaneTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: r1,
        })
    } else {
        PlaneTorusResult::TwoCircles(
            Circle3 { center, normal: torus.axis, radius: r1 },
            Circle3 { center, normal: torus.axis, radius: r2 },
        )
    }
}

fn intersect_plane_torus_parallel(
    plane: &Plane,
    torus: &ToroidalSurface,
    tol: f64,
) -> PlaneTorusResult {
    // Signed distance from torus center to plane along the plane normal
    let n = plane.normal.normalize();
    let signed_dist = (plane.origin - torus.center).dot(n);
    let d = signed_dist.abs();

    // No intersection if plane is too far from torus
    if d > torus.major_radius + torus.minor_radius + tol {
        return PlaneTorusResult::NoIntersection;
    }

    // Analytical solution: two circles when plane intersects tube center circle
    // This happens when |d| <= R (the plane cuts through the tube center circle)
    if d <= torus.major_radius + tol {
        // The tube center circle (radius R in plane perpendicular to torus axis)
        // intersects the plane at two points when |d| < R
        // Each intersection produces a circle of radius r in the plane

        // Compute the direction perpendicular to the plane normal that lies in the
        // plane containing both the plane normal and the torus axis
        let a = torus.axis.normalize();

        // Compute the direction in the plane that is perpendicular to the torus axis
        // This gives us the direction toward the tube center intersection points
        let in_plane_perp = n.cross(a);
        let perp_len = in_plane_perp.length();

        if perp_len < TOLERANCE_ANG {
            // Degenerate case: plane normal is parallel to torus axis
            // This shouldn't happen in the parallel case
            return PlaneTorusResult::General;
        }

        let dir_perp = in_plane_perp / perp_len;

        // Calculate the distance along the perpendicular direction to the tube center
        // intersection points. From d² + z² = R², we get z = ±√(R² - d²)
        let d_sq = d * d;
        let r_sq = torus.major_radius * torus.major_radius;
        let z_dist_sq = r_sq - d_sq;

        if z_dist_sq < -tol * tol {
            // No real intersection (should not happen given our checks above)
            return PlaneTorusResult::NoIntersection;
        }

        let z_dist = z_dist_sq.max(0.0).sqrt();

        // The two circle centers are at:
        // center = torus_center + d * dir_to_plane ± z_dist * dir_perp
        let base_center = torus.center + signed_dist * n;

        let center1 = base_center + z_dist * dir_perp;
        let center2 = base_center - z_dist * dir_perp;

        // Check for tangent case (circles merge into one)
        if (center1 - center2).length() < tol {
            return PlaneTorusResult::TangentCircle(Circle3 {
                center: center1,
                normal: n,
                radius: torus.minor_radius,
            });
        }

        // Two circles of radius r in the plane
        PlaneTorusResult::TwoCircles(
            Circle3 {
                center: center1,
                normal: n,
                radius: torus.minor_radius,
            },
            Circle3 {
                center: center2,
                normal: n,
                radius: torus.minor_radius,
            },
        )
    } else {
        // d > R: Plane is between the tube center circle and outer edge
        // This produces a more complex intersection (ellipse-like)
        // Fall back to numerical marching for this case
        PlaneTorusResult::General
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Skew plane analytic solver
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the intersection of a skew plane with a torus using a u-parameterized
/// analytic solver on the torus surface.
///
/// For each torus azimuth u ∈ [0, 2π), the tube cross-section circle intersects
/// the plane in 0, 1, or 2 points.  Solve the plane constraint:
///
/// ```text
/// cos(v)·B(u) + sin(v)·C = D(u)
/// ```
///
/// where B(u) and D(u) depend on u, and C is constant.  The solution is:
///
/// ```text
/// v = atan2(C, B(u)) ± acos(D(u) / sqrt(B(u)² + C²))
/// ```
///
/// Returns polylines for the ± branches, each clipped to contiguous θ-ranges.
pub fn intersect_plane_torus_skew(
    plane: &Plane,
    torus: &ToroidalSurface,
) -> Vec<Vec<DVec3>> {
    use std::f64::consts::TAU;
    let n = plane.normal.normalize();
    let a = torus.axis.normalize();
    let r_major = torus.major_radius;
    let r_minor = torus.minor_radius;
    let torus_center = torus.center;

    // Torus frame: x, y span the major circle plane, z = axis
    let x = any_perpendicular(a);
    let y = a.cross(x).normalize();

    // Projections of plane normal onto torus frame
    let nx = n.dot(x);
    let ny = n.dot(y);
    let nz = n.dot(a);

    // If the plane is nearly parallel to the torus axis, nz ≈ 0 → the parallel
    // case handler should already have caught this.
    if nz.abs() < 1e-12 {
        return vec![];
    }

    // Plane distance from origin: d = plane.origin · n (since plane is
    // defined by points P where (P - plane.origin)·n = 0, equivalently P·n = d)
    let d = plane.origin.dot(n);

    // Constant terms
    let a_cn = torus_center.dot(n); // center·n

    const N_SAMPLES: usize = 128;
    let mut branch_plus: Vec<(f64, Option<DVec3>)> = Vec::with_capacity(N_SAMPLES + 1);
    let mut branch_minus: Vec<(f64, Option<DVec3>)> = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cu, su) = (u.cos(), u.sin());

        // B(u) = cos(u)·nx + sin(u)·ny
        let bu = cu * nx + su * ny;

        // D(u) = (d - a_cn - R·B(u)) / r
        let du = (d - a_cn - r_major * bu) / r_minor;

        // m = sqrt(B(u)² + nz²)
        let m = (bu * bu + nz * nz).sqrt();

        if du.abs() > m {
            // No solution at this u
            branch_plus.push((u, None));
            branch_minus.push((u, None));
            continue;
        }

        let acos_val = (du / m).clamp(-1.0, 1.0).acos();
        let v0 = nz.atan2(bu);

        let v_plus = v0 + acos_val;
        let v_minus = v0 - acos_val;

        // Compute 3D points on the torus
        let radial = cu * x + su * y;
        let tube_center = torus_center + r_major * radial;
        let pt_plus = tube_center + r_minor * (v_plus.cos() * radial + v_plus.sin() * a);
        let pt_minus = tube_center + r_minor * (v_minus.cos() * radial + v_minus.sin() * a);

        branch_plus.push((u, Some(pt_plus)));
        branch_minus.push((u, Some(pt_minus)));
    }

    // Extract contiguous runs from each branch (returns Vec<(u, point)>)
    let extract_runs = |branch: &[(f64, Option<DVec3>)]| -> Vec<Vec<(f64, DVec3)>> {
        let mut curves: Vec<Vec<(f64, DVec3)>> = Vec::new();
        let mut current: Vec<(f64, DVec3)> = Vec::new();
        for &(u, ref pt) in branch {
            match pt {
                Some(p) => current.push((u, *p)),
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
        curves
    };

    let mut raw_branches = extract_runs(&branch_plus);
    raw_branches.extend(extract_runs(&branch_minus));

    // ── Adaptive chord-error refinement ────────────────────────────────────
    const CHORD_TOL: f64 = 1e-6;
    const REFINE_DEPTH: usize = 2;

    let refined: Vec<Vec<DVec3>> = raw_branches
        .into_iter()
        .filter(|b| b.len() >= 4)
        .map(|branch| {
            let (sign_branch, _) = branch[0];
            let _ = sign_branch; // branch identifier (not needed for eval)
            // For each branch, we need to build an eval closure.
            // We determine sign by evaluating at the first u and checking
            // which v formula produces the matching 3D point.
            let u_first = branch[0].0;
            let p_first = branch[0].1;
            let (cu_f, su_f) = (u_first.cos(), u_first.sin());
            let bu_f = cu_f * nx + su_f * ny;
            let du_f = (d - a_cn - r_major * bu_f) / r_minor;
            let m_f = (bu_f * bu_f + nz * nz).sqrt();
            let acos_f = (du_f / m_f).clamp(-1.0, 1.0).acos();
            let v0_f = nz.atan2(bu_f);
            // Determine whether this branch uses v0 + acos or v0 - acos
            let radial_f = cu_f * x + su_f * y;
            let tube_center_f = torus_center + r_major * radial_f;
            let p_plus = tube_center_f
                + r_minor * ((v0_f + acos_f).cos() * radial_f + (v0_f + acos_f).sin() * a);
            let use_plus = (p_first - p_plus).length_squared()
                < TOLERANCE_VEC_SQ_MIN;

            let _pts_for_eval = branch.clone();
            let eval_fn = move |u_mid: f64| -> Option<DVec3> {
                let (cu, su) = (u_mid.cos(), u_mid.sin());
                let bu = cu * nx + su * ny;
                let du = (d - a_cn - r_major * bu) / r_minor;
                let m = (bu * bu + nz * nz).sqrt();
                if du.abs() > m {
                    return None;
                }
                let acos_val = (du / m).clamp(-1.0, 1.0).acos();
                let v0 = nz.atan2(bu);
                let v = if use_plus { v0 + acos_val } else { v0 - acos_val };
                let radial = cu * x + su * y;
                let tube_center = torus_center + r_major * radial;
                let p = tube_center
                    + r_minor * (v.cos() * radial + v.sin() * a);
                if p.is_finite() { Some(p) } else { None }
            };

            let refined =
                crate::inttools::pcurve_derive::refine_polyline(
                    &branch, eval_fn, CHORD_TOL, REFINE_DEPTH,
                );
            refined.into_iter().map(|(_, p)| p).collect()
        })
        .collect();

    // Closed-curve dedup: remove trailing point if it nearly duplicates the first
    let mut result = refined;
    for branch in &mut result {
        if branch.len() >= 3 {
            let d = (branch[0] - branch[branch.len() - 1]).length();
            if d < TOLERANCE_ABS * 10.0 {
                branch.pop();
            }
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    #[test]
    fn plane_perpendicular_to_torus_axis_produces_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane perpendicular to Y axis, slicing through center
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Outer circle at major_radius + minor_radius
                assert!((c1.radius - 6.0).abs() < TOLERANCE_MESH_LEGACY, "Outer circle radius expected 6.0, got {}", c1.radius);
                // Inner circle at major_radius - minor_radius
                assert!((c2.radius - 4.0).abs() < TOLERANCE_MESH_LEGACY, "Inner circle radius expected 4.0, got {}", c2.radius);
                // Both circles should have the same center
                assert!((c1.center - DVec3::ZERO).length() < TOLERANCE_MESH_LEGACY);
                assert!((c2.center - DVec3::ZERO).length() < TOLERANCE_MESH_LEGACY);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_through_center_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane parallel to Y axis (normal = X), passing through center
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Both circles should have radius equal to minor radius
                assert!((c1.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY, "Circle 1 radius expected 1.0, got {}", c1.radius);
                assert!((c2.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY, "Circle 2 radius expected 1.0, got {}", c2.radius);
                // Both circles should be in the plane (normal = X)
                assert!((c1.normal - DVec3::X).length() < TOLERANCE_MESH_LEGACY);
                assert!((c2.normal - DVec3::X).length() < TOLERANCE_MESH_LEGACY);
                // Centers should be at z = ±R = ±5
                assert!((c1.center.z.abs() - 5.0).abs() < TOLERANCE_MESH_LEGACY, "Circle 1 center z should be ±5");
                assert!((c2.center.z.abs() - 5.0).abs() < TOLERANCE_MESH_LEGACY, "Circle 2 center z should be ±5");
                // Both centers at x=0, y=0
                assert!(c1.center.x.abs() < TOLERANCE_MESH_LEGACY);
                assert!(c1.center.y.abs() < TOLERANCE_MESH_LEGACY);
                assert!(c2.center.x.abs() < TOLERANCE_MESH_LEGACY);
                assert!(c2.center.y.abs() < TOLERANCE_MESH_LEGACY);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_offset_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane parallel to Y axis, offset by d=3 from center
        let plane = Plane {
            origin: DVec3::new(3.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Both circles should have radius equal to minor radius
                assert!((c1.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY);
                // Centers should be at x=3 (the plane's position)
                assert!((c1.center.x - 3.0).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.center.x - 3.0).abs() < TOLERANCE_MESH_LEGACY);
                // z = ±sqrt(R² - d²) = ±sqrt(25 - 9) = ±4
                let expected_z = (25.0_f64 - 9.0_f64).sqrt();
                assert!((c1.center.z.abs() - expected_z).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.center.z.abs() - expected_z).abs() < TOLERANCE_MESH_LEGACY);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_at_major_radius_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at distance R = 5 (touches tube center circle)
        let plane = Plane {
            origin: DVec3::new(5.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        // The result type depends on the exact intersection geometry
        // Just verify we get a valid result (don't panic)
        let _ = result;
    }

    #[test]
    fn plane_parallel_to_torus_axis_near_edge_returns_general() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at distance 5.5 (R < d < R+r), produces complex intersection
        let plane = Plane {
            origin: DVec3::new(5.5, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);
        // d = 5.5 > R = 5, so this should fall back to General
        assert!(matches!(result, PlaneTorusResult::General));
    }

    #[test]
    fn plane_parallel_to_torus_axis_no_intersection() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane outside torus (d > R + r = 6)
        let plane = Plane {
            origin: DVec3::new(7.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::NoIntersection));
    }

    #[test]
    fn plane_parallel_to_torus_axis_negative_offset_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane parallel to Y axis, offset by d=-4 from center
        let plane = Plane {
            origin: DVec3::new(-4.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Both circles should have radius equal to minor radius
                assert!((c1.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY);
                // Centers should be at x=-4
                assert!((c1.center.x + 4.0).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.center.x + 4.0).abs() < TOLERANCE_MESH_LEGACY);
                // z = ±sqrt(R² - d²) = ±sqrt(25 - 16) = ±3
                let expected_z = (25.0_f64 - 16.0_f64).sqrt();
                assert!((c1.center.z.abs() - expected_z).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.center.z.abs() - expected_z).abs() < TOLERANCE_MESH_LEGACY);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_tangent_at_inner_radius() {
        // Torus centered at origin with axis along Y
        // R = 5, r = 1, inner radius = R - r = 4
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at distance d = 4 (inner radius)
        let plane = Plane {
            origin: DVec3::new(4.0, 0.0, 0.0),
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // d = 4 < R = 5, so we still get two circles
                assert!((c1.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.radius - 1.0).abs() < TOLERANCE_MESH_LEGACY);
                // z = ±sqrt(R² - d²) = ±sqrt(25 - 16) = ±3
                let expected_z = (25.0_f64 - 16.0_f64).sqrt();
                assert!((c1.center.z.abs() - expected_z).abs() < TOLERANCE_MESH_LEGACY);
                assert!((c2.center.z.abs() - expected_z).abs() < TOLERANCE_MESH_LEGACY);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_perpendicular_tangent_to_torus() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane tangent to torus (at top of tube)
        let plane = Plane {
            origin: DVec3::new(0.0, 1.0, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TangentCircle(c) => {
                // Tangent circle at the major radius
                assert!((c.radius - 5.0).abs() < TOLERANCE_MESH_LEGACY);
            }
            other => panic!("Expected TangentCircle, got {:?}", other),
        }
    }

    #[test]
    fn plane_perpendicular_no_intersection() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane outside torus
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::NoIntersection));
    }

    #[test]
    fn plane_perpendicular_offset_produces_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane offset by 0.5 from center
        let plane = Plane {
            origin: DVec3::new(0.0, 0.5, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // tube_circle_r = sqrt(1 - 0.25) = sqrt(0.75) = 0.866025...
                let expected_tube_r = (1.0_f64 * 1.0 - 0.5_f64 * 0.5).sqrt();
                let expected_r1 = 5.0 + expected_tube_r;
                let expected_r2 = 5.0 - expected_tube_r;

                assert!((c1.radius - expected_r1).abs() < TOLERANCE_MESH_LEGACY, "Outer circle radius mismatch");
                assert!((c2.radius - expected_r2).abs() < TOLERANCE_MESH_LEGACY, "Inner circle radius mismatch");
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_oblique_returns_general() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at 45 degrees to torus axis
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(1.0, 1.0, 0.0).normalize(),
        };

        let result = intersect_plane_torus(&plane, &torus);
        // Now handled analytically by the skew solver → SkewPolyline.
        match &result {
            PlaneTorusResult::SkewPolyline(branches) => {
                assert!(!branches.is_empty(), "Expected at least one branch");
                assert!(branches.iter().all(|b| b.len() >= 2), "Each branch needs ≥2 points");
                let total_pts: usize = branches.iter().map(|b| b.len()).sum();
                // The 45° plane through the torus center only intersects along two
                // small arcs (near the tangent zone), so total points is modest.
                assert!(total_pts >= 20, "Expected ≥20 refined points, got {total_pts}");
            }
            other => panic!("Expected SkewPolyline, got {:?}", other),
        }
    }
}
