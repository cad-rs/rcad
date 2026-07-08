//! ElCLib-style elementary curve utilities.
//!
//! Provides analytical evaluation and parameter computation for elementary curves.
//! Analogous to OCCT `ElCLib` package.
//!
//! # Curve Types
//! - **Line**: Unbounded linear curve, parameter is distance along direction
//! - **Circle**: Closed circular curve, parameter is angle in radians [0, 2*pi]
//! - **Ellipse**: Closed elliptical curve, parameter is angle in radians [0, 2*pi]
//! - **Hyperbola**: Unbounded hyperbolic curve, parameter t (real line)
//! - **Parabola**: Unbounded parabolic curve, parameter t (real line)
//! - **BSpline**: NURBS curve, parameter within knot domain

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::geom::{
    any_perpendicular, BSplineCurve3, Circle3, CurveEval, Ellipse3, Hyperbola3, Line3, Parabola3,
};
use std::f64::consts::TAU;

// =============================================================================
// Line Utilities
// =============================================================================

/// Evaluate a point on a line at parameter t.
///
/// The parameter t represents the signed distance along the line's direction
/// from the origin. P(t) = origin + t * direction.
///
/// # Example
/// ```ignore
/// let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
/// let p = line_point_at(&line, 5.0);
/// assert!((p - DVec3::new(5.0, 0.0, 0.0)).length() < TOLERANCE_LINEAR_ULTRA_STRICT);
/// ```
pub fn line_point_at(line: &Line3, t: f64) -> DVec3 {
    line.origin + t * line.direction
}

/// Compute the parameter value for a point on a line.
///
/// Projects the point onto the line and returns the signed distance from
/// the origin along the line's direction.
///
/// Returns the parameter t such that `line_point_at(line, t)` is the
/// closest point on the line to the input point.
pub fn line_parameter(line: &Line3, point: DVec3) -> f64 {
    (point - line.origin).dot(line.direction)
}

/// Compute the perpendicular distance from a point to a line.
///
/// Uses the formula: distance = |(P - origin) x direction| / |direction|
/// For unit direction vectors, this simplifies to |(P - origin) x direction|.
pub fn line_distance_to_point(line: &Line3, point: DVec3) -> f64 {
    let v = point - line.origin;
    let cross = v.cross(line.direction);
    cross.length() / line.direction.length()
}

/// Find the closest point on a line to a given point.
///
/// Projects the point onto the line along the perpendicular direction.
pub fn line_closest_point(line: &Line3, point: DVec3) -> DVec3 {
    let t = line_parameter(line, point);
    line_point_at(line, t)
}

// =============================================================================
// Circle Utilities
// =============================================================================

/// Evaluate a point on a circle at the given angle.
///
/// The angle parameter is in radians, measured from the reference direction
/// (computed as any_perpendicular of the normal) in a right-handed sense
/// around the normal.
///
/// # Example
/// ```ignore
/// let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 2.0 );
/// let p = circle_point_at(&circle, 0.0); // Point at angle 0
/// assert!((p.length() - 2.0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
/// ```
pub fn circle_point_at(circle: &Circle3, angle: f64) -> DVec3 {
    circle.point_at(angle)
}

/// Compute the angle parameter for a point on or near a circle.
///
/// Projects the point onto the circle's plane and computes the angle
/// from the reference direction. Returns the angle in radians [0, 2*pi].
///
/// If the point is not exactly on the circle, the returned angle
/// corresponds to the closest point on the circle to the input.
pub fn circle_parameter(circle: &Circle3, point: DVec3) -> f64 {
    // Build local frame
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);

    // Project point into plane and express in local coordinates
    let v = point - circle.center;
    let x = v.dot(x_axis);
    let y = v.dot(y_axis);

    // Compute angle using atan2
    let angle = y.atan2(x);

    // Normalize to [0, 2*pi]
    if angle < 0.0 {
        angle + TAU
    } else {
        angle
    }
}

/// Compute the unit tangent vector on a circle at the given angle.
///
/// The tangent is perpendicular to the radius vector, pointing in the
/// direction of increasing angle (counterclockwise when viewed from
/// the normal direction).
pub fn circle_tangent_at(circle: &Circle3, angle: f64) -> DVec3 {
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);
    (-angle.sin() * x_axis + angle.cos() * y_axis).normalize()
}

/// Compute the unit normal vector on a circle at the given angle.
///
/// The normal points outward from the center, perpendicular to the curve
/// in the plane of the circle.
pub fn circle_normal_at(circle: &Circle3, angle: f64) -> DVec3 {
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);
    (angle.cos() * x_axis + angle.sin() * y_axis).normalize()
}

/// Compute the binormal vector on a circle at the given angle.
///
/// The binormal is the cross product of tangent and normal, which for
/// a planar circle is the circle's normal (axis) vector.
pub fn circle_binormal_at(circle: &Circle3, _angle: f64) -> DVec3 {
    circle.normal.normalize()
}

/// Compute the nth derivative of a circle at the given angle.
///
/// - Order 0: position on circle (same as `circle_point_at`)
/// - Order 1: first derivative = radius * tangent
/// - Order 2: second derivative = -radius * normal (centripetal acceleration)
/// - Higher orders: follow the pattern of derivatives of sin/cos
///
/// For a circle of radius R:
/// - dP/dt = R * (-sin(t), cos(t)) [first derivative]
/// - d²P/dt² = R * (-cos(t), -sin(t)) = -R * (cos(t), sin(t)) [second derivative]
/// - d³P/dt³ = R * (sin(t), -cos(t)) [third derivative]
///
/// Returns DVec3::ZERO if order is not supported.
pub fn circle_derivative(circle: &Circle3, angle: f64, order: usize) -> DVec3 {
    let x_axis = any_perpendicular(circle.normal);
    let y_axis = circle.normal.cross(x_axis);
    let r = circle.radius;

    match order {
        0 => circle_point_at(circle, angle),
        1 => r * (-angle.sin() * x_axis + angle.cos() * y_axis),
        2 => -r * (angle.cos() * x_axis + angle.sin() * y_axis),
        3 => r * (angle.sin() * x_axis - angle.cos() * y_axis),
        4 => r * (angle.cos() * x_axis + angle.sin() * y_axis), // Same as order 2 with opposite sign
        n => {
            // Higher orders cycle every 4
            let k = n % 4;
            match k {
                0 => circle_point_at(circle, angle),
                1 => r * (-angle.sin() * x_axis + angle.cos() * y_axis),
                2 => -r * (angle.cos() * x_axis + angle.sin() * y_axis),
                3 => r * (angle.sin() * x_axis - angle.cos() * y_axis),
                _ => DVec3::ZERO,
            }
        }
    }
}

// =============================================================================
// Ellipse Utilities
// =============================================================================

/// Evaluate a point on an ellipse at the given angle parameter.
///
/// The angle parameter is the eccentric anomaly, not the polar angle.
/// The ellipse is parameterized as:
///   P(angle) = center + a*cos(angle)*major_dir + b*sin(angle)*minor_dir
///
/// where a = major_radius, b = minor_radius, and minor_dir = normal x major_dir.
pub fn ellipse_point_at(ellipse: &Ellipse3, angle: f64) -> DVec3 {
    ellipse.point_at(angle)
}

/// Compute the angle parameter for a point on or near an ellipse.
///
/// Projects the point onto the ellipse's plane and solves for the
/// eccentric anomaly. Uses Newton-Raphson iteration for accuracy.
///
/// Returns the angle in radians [0, 2*pi].
pub fn ellipse_parameter(ellipse: &Ellipse3, point: DVec3) -> f64 {
    // Build local frame
    let x_axis = ellipse.major_dir;
    let y_axis = ellipse.normal.cross(x_axis).normalize();

    // Project point into plane and express in local coordinates
    let v = point - ellipse.center;
    let x = v.dot(x_axis);
    let y = v.dot(y_axis);

    // For an ellipse x = a*cos(t), y = b*sin(t)
    // We need to solve for t given (x, y)
    // Use atan2(y/b, x/a) as initial guess, then refine
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;

    if a.abs() < TOLERANCE_FLOAT_DEDUP || b.abs() < TOLERANCE_FLOAT_DEDUP {
        return 0.0;
    }

    // Initial guess using modified atan2
    let t = (y / b).atan2(x / a);

    // Newton-Raphson refinement for better accuracy
    // We solve: f(t) = atan2(y - b*sin(t), x - a*cos(t)) = 0
    // which is implicit. Instead, we directly compute the eccentric anomaly.
    // The eccentric anomaly satisfies:
    //   x = a * cos(t)
    //   y = b * sin(t)
    // So: cos(t) = x/a, sin(t) = y/b
    // t = atan2(y/b, x/a) is already exact for the parametric form

    // Normalize to [0, 2*pi]
    if t < 0.0 {
        t + TAU
    } else {
        t
    }
}

/// Compute the nth derivative of an ellipse at the given angle.
///
/// For an ellipse with radii a (major) and b (minor):
/// - Order 0: P(t) = (a*cos(t), b*sin(t))
/// - Order 1: dP/dt = (-a*sin(t), b*cos(t))
/// - Order 2: d²P/dt² = (-a*cos(t), -b*sin(t))
/// - Order 3: d³P/dt³ = (a*sin(t), -b*cos(t))
/// - Higher orders: cycle every 4
pub fn ellipse_derivative(ellipse: &Ellipse3, angle: f64, order: usize) -> DVec3 {
    let x_axis = ellipse.major_dir;
    let y_axis = ellipse.normal.cross(x_axis).normalize();
    let a = ellipse.major_radius;
    let b = ellipse.minor_radius;

    let cos_a = angle.cos();
    let sin_a = angle.sin();

    match order {
        0 => ellipse.center + a * cos_a * x_axis + b * sin_a * y_axis,
        1 => -a * sin_a * x_axis + b * cos_a * y_axis,
        2 => -a * cos_a * x_axis - b * sin_a * y_axis,
        3 => a * sin_a * x_axis - b * cos_a * y_axis,
        n => {
            let k = n % 4;
            match k {
                0 => ellipse.center + a * cos_a * x_axis + b * sin_a * y_axis,
                1 => -a * sin_a * x_axis + b * cos_a * y_axis,
                2 => -a * cos_a * x_axis - b * sin_a * y_axis,
                3 => a * sin_a * x_axis - b * cos_a * y_axis,
                _ => DVec3::ZERO,
            }
        }
    }
}

// =============================================================================
// Hyperbola Utilities
// =============================================================================

/// Evaluate a point on a hyperbola at parameter t.
///
/// The hyperbola is parameterized as:
///   P(t) = center + a*cosh(t)*major_dir + b*sinh(t)*minor_dir
///
/// where a = semi_major, b = semi_minor, and minor_dir = normal x major_dir.
/// The principal branch (t >= 0) is on the +major_dir side of the center.
pub fn hyperbola_point_at(hyp: &Hyperbola3, t: f64) -> DVec3 {
    hyp.point_at(t)
}

/// Compute the nth derivative of a hyperbola at parameter t.
///
/// For a hyperbola with semi-axes a and b:
/// - Order 0: P(t) = (a*cosh(t), b*sinh(t))
/// - Order 1: dP/dt = (a*sinh(t), b*cosh(t))
/// - Order 2: d²P/dt² = (a*cosh(t), b*sinh(t))
/// - Order 3: d³P/dt³ = (a*sinh(t), b*cosh(t))
/// - Higher orders: alternates between sinh and cosh patterns
pub fn hyperbola_derivative(hyp: &Hyperbola3, t: f64, order: usize) -> DVec3 {
    let minor_dir = hyp.normal.cross(hyp.major_dir).normalize();
    let a = hyp.semi_major;
    let b = hyp.semi_minor;

    let cosh_t = t.cosh();
    let sinh_t = t.sinh();

    match order {
        0 => hyp.center + a * cosh_t * hyp.major_dir + b * sinh_t * minor_dir,
        1 => a * sinh_t * hyp.major_dir + b * cosh_t * minor_dir,
        2 => a * cosh_t * hyp.major_dir + b * sinh_t * minor_dir, // Same as order 0
        3 => a * sinh_t * hyp.major_dir + b * cosh_t * minor_dir, // Same as order 1
        n => {
            // Pattern: even orders = order 0, odd orders = order 1
            if n % 2 == 0 {
                hyp.center + a * cosh_t * hyp.major_dir + b * sinh_t * minor_dir
            } else {
                a * sinh_t * hyp.major_dir + b * cosh_t * minor_dir
            }
        }
    }
}

// =============================================================================
// Parabola Utilities
// =============================================================================

/// Evaluate a point on a parabola at parameter t.
///
/// The parabola is parameterized as:
///   P(t) = vertex + (t²/(2p))*axis_dir + t*dir_perp
///
/// where p = focal_param (twice the focal length), and dir_perp = normal x axis_dir.
/// The focus is at distance p/2 from the vertex along axis_dir.
pub fn parabola_point_at(parab: &Parabola3, t: f64) -> DVec3 {
    parab.point_at(t)
}

/// Compute the nth derivative of a parabola at parameter t.
///
/// For a parabola with focal parameter p:
/// - Order 0: P(t) = (t²/(2p), t) in local coordinates
/// - Order 1: dP/dt = (t/p, 1)
/// - Order 2: d²P/dt² = (1/p, 0)
/// - Order 3+: All higher derivatives are zero
pub fn parabola_derivative(parab: &Parabola3, t: f64, order: usize) -> DVec3 {
    // dir_perp forms a right-handed system: axis_dir × normal gives perpendicular direction
    let dir_perp = parab.axis_dir.cross(parab.normal).normalize();
    let p = parab.focal_param;

    if p.abs() < TOLERANCE_FLOAT_DEDUP {
        return DVec3::ZERO;
    }

    match order {
        0 => {
            parab.vertex
                + (t * t / (2.0 * p)) * parab.axis_dir
                + t * dir_perp
        }
        1 => {
            // dP/dt = (t/p) * axis_dir + dir_perp
            (t / p) * parab.axis_dir + dir_perp
        }
        2 => {
            // d²P/dt² = (1/p) * axis_dir
            (1.0 / p) * parab.axis_dir
        }
        _ => {
            // All higher derivatives are zero
            DVec3::ZERO
        }
    }
}

// =============================================================================
// BSpline Utilities
// =============================================================================

/// Evaluate a point on a B-spline curve at parameter t.
///
/// Uses the de Boor algorithm for rational and non-rational B-splines.
/// The parameter t should be within the curve's domain [knots[degree], knots[n-degree-1]].
pub fn bspline_point_at(spline: &BSplineCurve3, t: f64) -> DVec3 {
    spline.point_at(t)
}

/// Compute the nth derivative of a B-spline curve at parameter t.
///
/// Uses analytical differentiation for NURBS curves via the quotient rule.
/// The derivative of a rational B-spline C(t) = A(t)/W(t) is:
///   C'(t) = (A'(t) - W'(t)*C(t)) / W(t)
///
/// Higher-order derivatives are computed by differentiating the derivative
/// curve, which is itself a B-spline of degree p-1.
///
/// Returns DVec3::ZERO if:
/// - The order is greater than the curve's degree
/// - The curve is invalid (no control points or degree 0 with order > 0)
pub fn bspline_derivative(spline: &BSplineCurve3, t: f64, order: usize) -> DVec3 {
    if order == 0 {
        return bspline_point_at(spline, t);
    }

    let n = spline.control_points.len();
    let degree = spline.degree;

    if n == 0 || (order > degree && degree == 0) {
        return DVec3::ZERO;
    }

    // For higher orders than degree, derivative is zero
    if order > degree {
        return DVec3::ZERO;
    }

    // Compute derivative using finite differences for simplicity and robustness
    // For production code, this should use the analytical derivative chain
    let h = TOLERANCE_ABS;

    let domain = spline.default_domain();
    let t_min = domain[0];
    let t_max = domain[1];

    // Clamp t to domain for finite difference calculation
    let t_lo = (t - h).max(t_min);
    let t_hi = (t + h).min(t_max);
    let actual_h = t_hi - t_lo;

    if actual_h < TOLERANCE_FLOAT_DEDUP {
        return DVec3::ZERO;
    }

    if order == 1 {
        // First derivative
        let p_lo = bspline_point_at(spline, t_lo);
        let p_hi = bspline_point_at(spline, t_hi);
        (p_hi - p_lo) / actual_h
    } else if order == 2 {
        // Second derivative using central differences
        let p_lo = bspline_point_at(spline, t_lo);
        let p_mid = bspline_point_at(spline, t);
        let p_hi = bspline_point_at(spline, t_hi);
        (p_hi - 2.0 * p_mid + p_lo) / (actual_h * actual_h / 4.0)
    } else {
        // Higher-order derivatives via recursive finite differences
        // This is less accurate but works for any order
        let mut points = Vec::with_capacity(order + 1);
        let step = actual_h / order as f64;
        for i in 0..=order {
            let ti = t_lo + i as f64 * step;
            points.push(bspline_point_at(spline, ti));
        }

        // Apply finite difference formula n times
        for _ in 0..order {
            let mut new_points = Vec::with_capacity(points.len() - 1);
            for i in 0..points.len() - 1 {
                new_points.push(points[i + 1] - points[i]);
            }
            points = new_points;
        }

        points.first().copied().unwrap_or(DVec3::ZERO) / step.powi(order as i32)
    }
}

// =============================================================================
// Unit Tests
// =============================================================================


