use super::helpers::{
    curve2d_derivative, curve2d_domain, curve2d_second_derivative, curve2d_tangent,
    refine_curve2d_distance, refine_point_curve2d_distance,
};
use super::*;

// =============================================================================
// ProjectPointOnCurve - Project point on 2D curve
// =============================================================================

/// Project a point onto a 2D curve, finding the closest point.
///
/// Uses sampling to find initial candidates, then Newton refinement.
///
/// # Arguments
/// * `point` - The point to project
/// * `curve` - The 2D curve to project onto
///
/// # Returns
/// A tuple (closest_point, parameter) where closest_point is on the curve
/// and parameter is the curve parameter at that point.
pub fn project_point_on_curve2d(point: DVec2, curve: &Curve2d) -> (DVec2, f64) {
    let domain = curve2d_domain(curve);

    // Sample the curve to find initial candidates
    let n_samples = 100;
    let mut best_t = domain[0];
    let mut best_dist = f64::INFINITY;

    for i in 0..=n_samples {
        let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n_samples as f64;
        let p = curve.point_at(t);
        let dist = (p - point).length();
        if dist < best_dist {
            best_dist = dist;
            best_t = t;
        }
    }

    // Newton refinement
    let refined_t = refine_point_curve2d_distance(curve, domain, point, best_t);
    let closest = curve.point_at(refined_t);

    (closest, refined_t)
}

// =============================================================================
// ExtremaCurveCurve - Distance between 2D curves
// =============================================================================

/// Compute the minimum distance between two 2D curves.
///
/// Uses sampling to find initial candidates, then Newton refinement.
///
/// # Arguments
/// * `curve1` - First 2D curve
/// * `curve2` - Second 2D curve
///
/// # Returns
/// A tuple (distance, param1, param2) where distance is the minimum Euclidean
/// distance between the curves, and param1, param2 are the parameters at the
/// closest points.
pub fn distance_between_curves2d(curve1: &Curve2d, curve2: &Curve2d) -> (f64, f64, f64) {
    let domain1 = curve2d_domain(curve1);
    let domain2 = curve2d_domain(curve2);

    // Sample both curves to find initial candidates
    let n_samples = 48;
    let mut best_dist = f64::INFINITY;
    let mut best_t1 = domain1[0];
    let mut best_t2 = domain2[0];

    for i in 0..=n_samples {
        let t1 = domain1[0] + (domain1[1] - domain1[0]) * i as f64 / n_samples as f64;
        let p1 = curve1.point_at(t1);

        for j in 0..=n_samples {
            let t2 = domain2[0] + (domain2[1] - domain2[0]) * j as f64 / n_samples as f64;
            let p2 = curve2.point_at(t2);
            let dist = (p2 - p1).length();

            if dist < best_dist {
                best_dist = dist;
                best_t1 = t1;
                best_t2 = t2;
            }
        }
    }

    // Newton refinement
    let (refined_t1, refined_t2) =
        refine_curve2d_distance(curve1, curve2, domain1, domain2, best_t1, best_t2);
    let p1 = curve1.point_at(refined_t1);
    let p2 = curve2.point_at(refined_t2);
    let final_dist = (p2 - p1).length();

    (final_dist, refined_t1, refined_t2)
}

// =============================================================================
// ExtremaCurvePoint - Distance from point to 2D curve
// =============================================================================

/// Compute the distance from a point to a 2D curve.
///
/// # Arguments
/// * `point` - The query point
/// * `curve` - The 2D curve
///
/// # Returns
/// A tuple (distance, parameter) where distance is the minimum Euclidean
/// distance from the point to the curve, and parameter is the curve parameter
/// at the closest point.
pub fn distance_point_to_curve2d(point: DVec2, curve: &Curve2d) -> (f64, f64) {
    let (closest, param) = project_point_on_curve2d(point, curve);
    let distance = (closest - point).length();
    (distance, param)
}

// =============================================================================
// Angle and Curvature Analysis
// =============================================================================

/// Compute the angle of the tangent vector at a parameter on a 2D curve.
///
/// The angle is measured from the positive X-axis, in radians, in the
/// counter-clockwise direction.
///
/// # Arguments
/// * `curve` - The 2D curve
/// * `t` - Parameter value
///
/// # Returns
/// The angle in radians of the tangent vector at parameter t.
pub fn curve2d_angle_at(curve: &Curve2d, t: f64) -> f64 {
    let tangent = curve2d_tangent(curve, t);
    tangent.y.atan2(tangent.x)
}

/// Compute the curvature at a parameter on a 2D curve.
///
/// Curvature is defined as |dT/ds| where T is the unit tangent and s is
/// the arc length. For a parametric curve C(t), this is:
///   kappa = |C' x C''| / |C'|^3
///
/// # Arguments
/// * `curve` - The 2D curve
/// * `t` - Parameter value
///
/// # Returns
/// The curvature value (positive for counter-clockwise turning, negative for
/// clockwise turning).
pub fn curve2d_curvature_at(curve: &Curve2d, t: f64) -> f64 {
    let d1 = curve2d_derivative(curve, t);
    let d2 = curve2d_second_derivative(curve, t);

    // In 2D, the cross product magnitude is |x1*y2 - y1*x2|
    let cross = d1.x * d2.y - d1.y * d2.x;
    let speed = d1.length();

    if speed < TOLERANCE_FLOAT_DEDUP {
        return 0.0;
    }

    cross / speed.powi(3)
}
