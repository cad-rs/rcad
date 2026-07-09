use super::*;
use super::helpers::{chord_length_params_2d, clamped_knots_from_params, solve_interpolation_2d};

// =============================================================================
// PointsToBSpline - Fit BSpline to 2D points
// =============================================================================

/// Fit a B-spline curve through a set of 2D points with specified degree.
///
/// Uses chord-length parameterization and builds a clamped B-spline.
/// This is a convenience wrapper around the kernel's interpolate_points_2d.
///
/// # Arguments
/// * `points` - Slice of 2D points to fit
/// * `degree` - Desired degree (will be clamped to n-1 for n points)
///
/// # Returns
/// A BSplineCurve2 that approximates the input points.
pub fn points_to_bspline2d(points: &[DVec2], degree: usize) -> BSplineCurve2 {
    let n = points.len();
    if n < 2 {
        return BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: points.to_vec(),
            weights: vec![1.0; n.max(1)],
        };
    }

    let actual_degree = degree.min(n - 1);

    // Use chord-length parameterization
    let params = chord_length_params_2d(points);
    let knots = clamped_knots_from_params(&params, actual_degree);

    // Build collocation matrix and solve
    let control_points = solve_interpolation_2d(&params, &knots, actual_degree, points);

    BSplineCurve2 {
        degree: actual_degree,
        knots,
        control_points,
        weights: vec![1.0; n],
    }
}

/// Fit a B-spline curve through a set of 2D points with cubic interpolation.
///
/// Equivalent to calling `points_to_bspline2d(points, 3)`.
/// The curve passes exactly through all input points.
///
/// # Arguments
/// * `points` - Slice of 2D points to interpolate
///
/// # Returns
/// A cubic BSplineCurve2 passing through all points.
pub fn points_to_bspline2d_interpolate(points: &[DVec2]) -> BSplineCurve2 {
    points_to_bspline2d(points, 3)
}
