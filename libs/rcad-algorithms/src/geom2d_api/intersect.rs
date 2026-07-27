use super::helpers::{curve2d_domain, refine_curve2d_intersection};
use super::*;

// =============================================================================
// InterCurveCurve - 2D curve-curve intersection
// =============================================================================

/// Find intersection points between two 2D curves.
///
/// Uses sampling to find initial candidates, then Newton refinement for accuracy.
/// Returns all intersection points within the given tolerance.
///
/// # Arguments
/// * `curve1` - First 2D curve
/// * `curve2` - Second 2D curve
/// * `tol` - Tolerance for considering points as coincident
///
/// # Returns
/// Vector of intersection points with parameters on each curve.
pub fn intersect_curves2d(
    curve1: &Curve2d,
    curve2: &Curve2d,
    tol: f64,
) -> Vec<Curve2dIntersection> {
    let domain1 = curve2d_domain(curve1);
    let domain2 = curve2d_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 64;
    let mut candidates: Vec<(f64, f64, f64)> = Vec::new(); // (dist, t1, t2)

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < tol * 10.0 {
                candidates.push((dist, t1, t2));
            }
        }
    }

    // Sort by distance and refine candidates
    candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut intersections: Vec<Curve2dIntersection> = Vec::new();

    for (_, t1, t2) in candidates {
        // Newton refinement
        let (refined_t1, refined_t2) =
            refine_curve2d_intersection(curve1, curve2, domain1, domain2, t1, t2);

        let p1 = curve1.point_at(refined_t1);
        let p2 = curve2.point_at(refined_t2);
        let dist = (p2 - p1).length();

        if dist < tol {
            // Check if this intersection is already found
            let is_duplicate = intersections.iter().any(|int| {
                (int.param1 - refined_t1).abs() < tol * 10.0
                    && (int.param2 - refined_t2).abs() < tol * 10.0
            });

            if !is_duplicate {
                intersections.push(Curve2dIntersection {
                    point: (p1 + p2) * 0.5,
                    param1: refined_t1,
                    param2: refined_t2,
                });
            }
        }
    }

    intersections
}
