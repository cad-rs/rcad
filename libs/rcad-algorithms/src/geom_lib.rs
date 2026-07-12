//! GeomLib-style geometry utilities.
//!
//! This module provides additional geometry utilities analogous to OpenCASCADE's GeomLib package.
//!
//! # Closure Checking
//! - [`is_curve_closed`] - Check if a curve is closed within tolerance
//! - [`is_surface_u_closed`] - Check if a surface is closed in U direction
//! - [`is_surface_v_closed`] - Check if a surface is closed in V direction
//!
//! # Degeneracy Removal
//! - [`remove_degenerate_curve_sections`] - Remove degenerate portions from a curve
//!
//! # Normal Estimation
//! - [`estimate_normal`] - Estimate surface normal at a UV point
//! - [`estimate_normal_by_neighbors`] - Estimate normal using neighboring samples
//!
//! # Curve Tools
//! - [`reverse_curve`] - Reverse the parametric direction of a curve
//! - [`trim_curve`] - Trim a curve to a parameter range
//! - [`transform_curve`] - Apply an affine transformation to a curve
//!
//! # Surface Tools
//! - [`reverse_surface_u`] - Reverse the U parametric direction of a surface
//! - [`reverse_surface_v`] - Reverse the V parametric direction of a surface
//! - [`trim_surface`] - Trim a surface to UV bounds
//! - [`transform_surface`] - Apply an affine transformation to a surface
//!
//! # Continuity Checking
//! - [`check_curve_continuity`] - Check the continuity order of a curve
//! - [`check_surface_continuity`] - Check the continuity orders of a surface

use crate::tolerance::*;
use glam::{DAffine3, DVec3};

use rcad_kernel::geom::{
    any_perpendicular, BSplineCurve3, BSplineSurface, BezierCurve3, BezierSurface, Circle3,
    ConicalSurface, Curve3, CurveEval, CylindricalSurface, Ellipse3, Line3, Plane,
    SphericalSurface, Surface3, SurfaceEval, ToroidalSurface, TrimmedSurface,
};

// =============================================================================
// Closure Checking
// =============================================================================

/// Check if a curve is closed within the given tolerance.
///
/// A curve is considered closed if the distance between its start and end points
/// is less than the tolerance.
///
/// # Arguments
/// * `curve` - The curve to check
/// * `tol` - The tolerance for closure checking
///
/// # Returns
/// `true` if the curve is closed within tolerance
///
/// # Example
/// ```rust
/// # use rcad_algorithms::tolerance::*;
/// use glam::DVec3;
/// use rcad_kernel::geom::{Circle3, Curve3};
/// use rcad_algorithms::geom_lib::is_curve_closed;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
///
/// let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
/// assert!(is_curve_closed(&Curve3::Circle(circle), TOLERANCE_MESH_LEGACY));
/// ```
pub fn is_curve_closed(curve: &Curve3, tol: f64) -> bool {
    let [t0, t1] = curve.default_domain();

    // For infinite domains, the curve is not closed in the geometric sense
    if t0.is_infinite() || t1.is_infinite() {
        return false;
    }

    let p0 = curve.point_at(t0);
    let p1 = curve.point_at(t1);

    p0.distance(p1) < tol
}

/// Check if a surface is closed in the U direction within the given tolerance.
///
/// A surface is considered closed in U if the iso-V curves at U min and U max
/// are coincident within tolerance.
///
/// # Arguments
/// * `surface` - The surface to check
/// * `tol` - The tolerance for closure checking
///
/// # Returns
/// `true` if the surface is closed in U within tolerance
///
/// # Example
/// ```rust
/// # use rcad_algorithms::tolerance::*;
/// use glam::DVec3;
/// use rcad_kernel::geom::{CylindricalSurface, Surface3};
/// use rcad_algorithms::geom_lib::is_surface_u_closed;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
///
/// let cylinder = CylindricalSurface {
///     origin: DVec3::ZERO,
///     axis: DVec3::Z,
///     radius: 1.0,
/// };
/// assert!(is_surface_u_closed(&Surface3::Cylinder(cylinder), TOLERANCE_MESH_LEGACY));
/// ```
pub fn is_surface_u_closed(surface: &Surface3, tol: f64) -> bool {
    let [u0, u1, v0, v1] = surface.default_domain();

    // For infinite domains, the surface is not closed
    if u0.is_infinite() || u1.is_infinite() {
        return false;
    }

    // Sample the boundary curves at U min and U max
    let num_samples = 21;
    let v_step = (v1 - v0) / (num_samples - 1) as f64;

    for i in 0..num_samples {
        let v = v0 + i as f64 * v_step;
        let p0 = surface.point_at(u0, v);
        let p1 = surface.point_at(u1, v);

        if p0.distance(p1) >= tol {
            return false;
        }
    }

    true
}

/// Check if a surface is closed in the V direction within the given tolerance.
///
/// A surface is considered closed in V if the iso-U curves at V min and V max
/// are coincident within tolerance.
///
/// # Arguments
/// * `surface` - The surface to check
/// * `tol` - The tolerance for closure checking
///
/// # Returns
/// `true` if the surface is closed in V within tolerance
///
/// # Example
/// ```rust
/// # use rcad_algorithms::tolerance::*;
/// use glam::DVec3;
/// use rcad_kernel::geom::{SphericalSurface, Surface3};
/// use rcad_algorithms::geom_lib::is_surface_v_closed;
/// use rcad_algorithms::tolerance::TOLERANCE_MESH_LEGACY;
///
/// let sphere = SphericalSurface {
///     center: DVec3::ZERO,
///     axis: DVec3::Z,
///     radius: 1.0,
///     ref_dir: any_perpendicular(DVec3::Z),
/// };
/// // Sphere is closed in U but not in V (has poles)
/// assert!(!is_surface_v_closed(&Surface3::Sphere(sphere), TOLERANCE_MESH_LEGACY));
/// ```
pub fn is_surface_v_closed(surface: &Surface3, tol: f64) -> bool {
    let [u0, u1, v0, v1] = surface.default_domain();

    // For infinite domains, the surface is not closed
    if v0.is_infinite() || v1.is_infinite() {
        return false;
    }

    // Sample the boundary curves at V min and V max
    let num_samples = 21;
    let u_step = (u1 - u0) / (num_samples - 1) as f64;

    for i in 0..num_samples {
        let u = u0 + i as f64 * u_step;
        let p0 = surface.point_at(u, v0);
        let p1 = surface.point_at(u, v1);

        if p0.distance(p1) >= tol {
            return false;
        }
    }

    true
}

// =============================================================================
// Degeneracy Removal
// =============================================================================

/// Remove degenerate sections from a curve.
///
/// Degenerate sections are portions of a curve where consecutive points are
/// within the tolerance distance of each other, indicating a collapsed segment.
///
/// # Arguments
/// * `curve` - The curve to process
/// * `tol` - The tolerance for detecting degenerate sections
///
/// # Returns
/// * `Some(Curve3)` - A new curve with degenerate sections removed, if any were found
/// * `None` - If no degenerate sections were found or removal is not applicable
///
/// # Note
/// Currently only implemented for BSpline curves. Other curve types return `None`.
pub fn remove_degenerate_curve_sections(curve: &Curve3, tol: f64) -> Option<Curve3> {
    match curve {
        Curve3::BSpline(bspline) => remove_degenerate_bspline_sections(bspline, tol)
            .map(Curve3::BSpline),
        Curve3::Bezier(bezier) => remove_degenerate_bezier_sections(bezier, tol)
            .map(Curve3::Bezier),
        _ => None, // Other curve types don't typically have degenerate sections
    }
}

fn remove_degenerate_bspline_sections(
    curve: &BSplineCurve3,
    tol: f64,
) -> Option<BSplineCurve3> {
    let n = curve.control_points.len();
    if n < 2 {
        return None;
    }

    // Find non-degenerate control point ranges
    let mut new_ctrl = Vec::new();
    let mut new_weights = Vec::new();
    let mut any_removed = false;

    for i in 0..n {
        let is_degenerate = if i > 0 {
            curve.control_points[i]
                .distance(curve.control_points[i - 1]) < tol
        } else {
            false
        };

        if !is_degenerate {
            new_ctrl.push(curve.control_points[i]);
            new_weights.push(curve.weights[i]);
        } else {
            any_removed = true;
        }
    }

    if !any_removed || new_ctrl.len() < 2 {
        return None;
    }

    // Rebuild knot vector for the new control point count
    let degree = curve.degree;
    let new_n = new_ctrl.len();
    let mut new_knots = Vec::with_capacity(new_n + degree + 1);

    // Create a uniform clamped knot vector for the reduced curve
    for _ in 0..=degree {
        new_knots.push(0.0);
    }
    for i in 1..(new_n - degree) {
        new_knots.push(i as f64 / (new_n - degree) as f64);
    }
    for _ in 0..=degree {
        new_knots.push(1.0);
    }

    Some(BSplineCurve3 {
        degree,
        knots: new_knots,
        control_points: new_ctrl,
        weights: new_weights,
    })
}

fn remove_degenerate_bezier_sections(
    curve: &BezierCurve3,
    tol: f64,
) -> Option<BezierCurve3> {
    let n = curve.control_points.len();
    if n < 2 {
        return None;
    }

    let mut new_ctrl = Vec::new();
    let mut new_weights = Vec::new();
    let mut any_removed = false;

    for i in 0..n {
        let is_degenerate = if i > 0 {
            curve.control_points[i]
                .distance(curve.control_points[i - 1]) < tol
        } else {
            false
        };

        if !is_degenerate {
            new_ctrl.push(curve.control_points[i]);
            new_weights.push(curve.weights[i]);
        } else {
            any_removed = true;
        }
    }

    if !any_removed || new_ctrl.len() < 2 {
        return None;
    }

    Some(BezierCurve3 {
        control_points: new_ctrl,
        weights: new_weights,
    })
}

// =============================================================================
// Normal Estimation
// =============================================================================

/// Estimate the surface normal at a given UV parameter.
///
/// This computes the cross product of the partial derivatives with respect to
/// U and V, normalized to a unit vector.
///
/// # Arguments
/// * `surface` - The surface to evaluate
/// * `u` - The U parameter
/// * `v` - The V parameter
///
/// # Returns
/// The estimated unit normal vector at the given point
///
/// # Note
/// For surfaces with well-defined analytical normals (planes, spheres, etc.),
/// prefer using `surface.normal_at(u, v)` directly.
pub fn estimate_normal(surface: &Surface3, u: f64, v: f64) -> DVec3 {
    // First try to use the surface's built-in normal computation
    let normal = surface.normal_at(u, v);
    if normal.length_squared() > TOLERANCE_LEN_MIN {
        return normal;
    }

    // Fallback to numerical differentiation if the analytical normal is degenerate
    estimate_normal_numerical(surface, u, v)
}

fn estimate_normal_numerical(surface: &Surface3, u: f64, v: f64) -> DVec3 {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Use a small step for numerical differentiation
    let h_u = ((u1 - u0) * TOLERANCE_MESH_LEGACY).max(TOLERANCE_LINEAR_RELAX_8);
    let h_v = ((v1 - v0) * TOLERANCE_MESH_LEGACY).max(TOLERANCE_LINEAR_RELAX_8);

    // Compute partial derivatives numerically
    let p = surface.point_at(u, v);
    let pu = surface.point_at(u + h_u, v);
    let pv = surface.point_at(u, v + h_v);

    let du = (pu - p) / h_u;
    let dv = (pv - p) / h_v;

    let normal = du.cross(dv);
    if normal.length_squared() > TOLERANCE_LEN_MIN {
        normal.normalize()
    } else {
        // Attempt to find any valid normal using cross products with coordinate axes
        any_perpendicular(du)
    }
}

/// Estimate the surface normal using neighboring samples for robustness.
///
/// This method samples the surface at neighboring points and uses the average
/// normal direction, which can be more robust for surfaces with local irregularities.
///
/// # Arguments
/// * `surface` - The surface to evaluate
/// * `u` - The U parameter
/// * `v` - The V parameter
/// * `step` - The step size for sampling neighbors (in parameter space)
///
/// # Returns
/// The estimated unit normal vector, averaged from neighboring samples
pub fn estimate_normal_by_neighbors(surface: &Surface3, u: f64, v: f64, step: f64) -> DVec3 {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Clamp step to valid range
    let step_u = step.min((u1 - u0) * 0.1_f64).max(TOLERANCE_LINEAR_RELAX_8);
    let step_v = step.min((v1 - v0) * 0.1_f64).max(TOLERANCE_LINEAR_RELAX_8);

    // Sample a 3x3 grid around the point
    let mut normals: Vec<DVec3> = Vec::new();

    for du in [-step_u, 0.0, step_u] {
        for dv in [-step_v, 0.0, step_v] {
            let u_s = (u + du).clamp(u0, u1);
            let v_s = (v + dv).clamp(v0, v1);

            let n = estimate_normal(surface, u_s, v_s);
            if n.length_squared() > TOLERANCE_LEN_MIN {
                normals.push(n);
            }
        }
    }

    if normals.is_empty() {
        return DVec3::Z; // Fallback
    }

    // Average the normals
    let sum: DVec3 = normals.iter().sum();
    if sum.length_squared() > TOLERANCE_LEN_MIN {
        sum.normalize()
    } else {
        normals[0]
    }
}

// =============================================================================
// Curve Tools
// =============================================================================

/// Reverse the parametric direction of a curve.
///
/// This creates a new curve that traces the same geometry but in the opposite
/// parametric direction. The parameter domain remains the same, but evaluation
/// at parameter `t` now returns what was previously at `t_max - t`.
///
/// # Arguments
/// * `curve` - The curve to reverse
///
/// # Returns
/// A new curve with reversed parametric direction
pub fn reverse_curve(curve: Curve3) -> Curve3 {
    match curve {
        Curve3::Line(mut line) => {
            line.direction = -line.direction;
            Curve3::Line(line)
        }
        Curve3::Circle(mut circle) => {
            circle.normal = -circle.normal;
            Curve3::Circle(circle)
        }
        Curve3::Ellipse(mut ellipse) => {
            ellipse.normal = -ellipse.normal;
            Curve3::Ellipse(ellipse)
        }
        Curve3::BSpline(bspline) => {
            Curve3::BSpline(reverse_bspline_curve(&bspline))
        }
        Curve3::Bezier(bezier) => {
            Curve3::Bezier(reverse_bezier_curve(&bezier))
        }
        Curve3::Offset(mut offset) => {
            offset.offset_distance = -offset.offset_distance;
            offset.basis = Box::new(reverse_curve(*offset.basis));
            Curve3::Offset(offset)
        }
        Curve3::Hyperbola(h) => Curve3::Hyperbola(h), // Reversing hyperbola has no simple representation
        Curve3::Parabola(h) => Curve3::Parabola(h), // Reversing parabola has no simple representation
        Curve3::CircularHelix(mut h) => {
            h.pitch = -h.pitch;
            Curve3::CircularHelix(h)
        }
        Curve3::SineWave(mut s) => {
            s.phase = -s.phase;
            Curve3::SineWave(s)
        }
    }
}

fn reverse_bspline_curve(curve: &BSplineCurve3) -> BSplineCurve3 {
    let degree = curve.degree;

    // Reverse control points and weights
    let new_ctrl: Vec<DVec3> = curve.control_points.iter().rev().cloned().collect();
    let new_weights: Vec<f64> = curve.weights.iter().rev().cloned().collect();

    // Reverse and renormalize the knot vector
    let knots = &curve.knots;
    let k0 = knots[0];
    let k1 = knots[knots.len() - 1];

    let new_knots: Vec<f64> = knots
        .iter()
        .rev()
        .map(|&k| k0 + (k1 - k))
        .collect();

    BSplineCurve3 {
        degree,
        knots: new_knots,
        control_points: new_ctrl,
        weights: new_weights,
    }
}

fn reverse_bezier_curve(curve: &BezierCurve3) -> BezierCurve3 {
    BezierCurve3 {
        control_points: curve.control_points.iter().rev().cloned().collect(),
        weights: curve.weights.iter().rev().cloned().collect(),
    }
}

/// Trim a curve to the specified parameter range.
///
/// This creates a new curve that represents the portion of the original curve
/// between parameters `t1` and `t2`.
///
/// # Arguments
/// * `curve` - The curve to trim
/// * `t1` - The start parameter
/// * `t2` - The end parameter
///
/// # Returns
/// A new curve representing the trimmed portion
///
/// # Note
/// For BSpline curves, this uses exact knot insertion. For other curve types,
/// this creates a trimmed representation using appropriate techniques.
pub fn trim_curve(curve: &Curve3, t1: f64, t2: f64) -> Curve3 {
    match curve {
        Curve3::BSpline(bspline) => {
            let trimmed = trim_bspline_curve(bspline, t1, t2);
            Curve3::BSpline(trimmed)
        }
        Curve3::Bezier(bezier) => {
            // For Bezier curves, subdivision using de Casteljau
            let trimmed = trim_bezier_curve(bezier, t1, t2);
            Curve3::Bezier(trimmed)
        }
        Curve3::Line(line) => {
            // For lines, create a new line segment
            let p1 = line.point_at(t1);
            let p2 = line.point_at(t2);
            Curve3::Line(Line3 {
                origin: p1,
                direction: (p2 - p1).normalize(),
            })
        }
        Curve3::Circle(circle) => {
            // For circles, create an arc (still represented as a circle for now)
            // A proper implementation would create a BSpline arc
            Curve3::Circle(*circle)
        }
        Curve3::Ellipse(ellipse) => {
            Curve3::Ellipse(*ellipse)
        }
        _ => curve.clone(), // Fallback: return unchanged
    }
}

fn trim_bspline_curve(curve: &BSplineCurve3, t1: f64, t2: f64) -> BSplineCurve3 {
    // Use the kernel's trim_curve from extend module
    rcad_kernel::extend::trim_curve(curve, t1, t2)
}

fn trim_bezier_curve(curve: &BezierCurve3, t1: f64, t2: f64) -> BezierCurve3 {
    // Subdivide using de Casteljau algorithm
    // First subdivide at t2, then take the first part and subdivide at t1/t2
    let upper = subdivide_bezier_at(curve, t2).1;
    let t1_adj = t1 / t2;
    subdivide_bezier_at(&upper, t1_adj).0
}

fn subdivide_bezier_at(curve: &BezierCurve3, t: f64) -> (BezierCurve3, BezierCurve3) {
    let n = curve.control_points.len();
    if n < 2 {
        return (curve.clone(), curve.clone());
    }

    let one_minus_t = 1.0 - t;

    // de Casteljau subdivision
    let mut lower = Vec::with_capacity(n);
    let mut upper = Vec::with_capacity(n);

    // We need to track weights separately for rational Bezier
    let mut points = curve.control_points.clone();
    let mut weights = curve.weights.clone();

    lower.push(points[0]);
    let mut lower_weights = vec![weights[0]];
    upper.push(points[n - 1]);
    let mut upper_weights = vec![weights[n - 1]];

    for i in 1..n {
        let mut new_points = Vec::with_capacity(n - i);
        let mut new_weights = Vec::with_capacity(n - i);

        for j in 0..(n - i) {
            let w0 = weights[j];
            let w1 = weights[j + 1];
            let p0 = points[j];
            let p1 = points[j + 1];

            // Weighted interpolation for rational Bezier
            let w = one_minus_t * w0 + t * w1;
            let p = if w.abs() > TOLERANCE_FLOAT_LOOSE {
                (one_minus_t * w0 * p0 + t * w1 * p1) / w
            } else {
                p0
            };

            new_points.push(p);
            new_weights.push(w);
        }

        points = new_points;
        weights = new_weights;

        lower.push(points[0]);
        lower_weights.push(weights[0]);
        upper.push(points[points.len() - 1]);
        upper_weights.push(weights[weights.len() - 1]);
    }

    upper.reverse();
    upper_weights.reverse();

    (
        BezierCurve3 {
            control_points: lower,
            weights: lower_weights,
        },
        BezierCurve3 {
            control_points: upper,
            weights: upper_weights,
        },
    )
}

/// Apply an affine transformation to a curve.
///
/// # Arguments
/// * `curve` - The curve to transform
/// * `transform` - The affine transformation to apply
///
/// # Returns
/// A new curve with the transformation applied
///
/// # Note
/// For analytic curves (lines, circles, etc.), this transforms the defining
/// parameters. For BSpline curves, this transforms the control points.
pub fn transform_curve(curve: &Curve3, transform: DAffine3) -> Curve3 {
    match curve {
        Curve3::Line(line) => {
            Curve3::Line(Line3 {
                origin: transform.transform_point3(line.origin),
                direction: transform.transform_vector3(line.direction).normalize(),
            })
        }
        Curve3::Circle(circle) => {
            Curve3::Circle(Circle3::new(transform.transform_point3(circle.center), transform.transform_vector3(circle.normal).normalize(), circle.radius * transform.matrix3.x_axis.length(),
            ))
        }
        Curve3::Ellipse(ellipse) => {
            let scale = transform.matrix3.x_axis.length();
            Curve3::Ellipse(Ellipse3 {
                center: transform.transform_point3(ellipse.center),
                normal: transform.transform_vector3(ellipse.normal).normalize(),
                major_dir: transform.transform_vector3(ellipse.major_dir).normalize(),
                major_radius: ellipse.major_radius * scale,
                minor_radius: ellipse.minor_radius * scale,
            })
        }
        Curve3::BSpline(bspline) => {
            Curve3::BSpline(transform_bspline_curve(bspline, transform))
        }
        Curve3::Bezier(bezier) => {
            Curve3::Bezier(transform_bezier_curve(bezier, transform))
        }
        Curve3::Offset(offset) => {
            let mut new_offset = offset.clone();
            new_offset.basis = Box::new(transform_curve(&new_offset.basis, transform));
            new_offset.offset_dir = transform.transform_vector3(new_offset.offset_dir).normalize();
            Curve3::Offset(new_offset)
        }
        Curve3::Hyperbola(h) => {
            let scale = transform.matrix3.x_axis.length();
            Curve3::Hyperbola(rcad_kernel::geom::Hyperbola3 {
                center: transform.transform_point3(h.center),
                normal: transform.transform_vector3(h.normal).normalize(),
                major_dir: transform.transform_vector3(h.major_dir).normalize(),
                semi_major: h.semi_major * scale,
                semi_minor: h.semi_minor * scale,
            })
        }
        Curve3::Parabola(p) => {
            let scale = transform.matrix3.x_axis.length();
            Curve3::Parabola(rcad_kernel::geom::Parabola3 {
                vertex: transform.transform_point3(p.vertex),
                normal: transform.transform_vector3(p.normal).normalize(),
                axis_dir: transform.transform_vector3(p.axis_dir).normalize(),
                focal_param: p.focal_param * scale,
            })
        }
        Curve3::CircularHelix(h) => {
            let scale = transform.matrix3.x_axis.length();
            Curve3::CircularHelix(rcad_kernel::geom::CircularHelix3 {
                origin: transform.transform_point3(h.origin),
                axis: transform.transform_vector3(h.axis).normalize(),
                ref_dir: transform.transform_vector3(h.ref_dir).normalize(),
                radius: h.radius * scale,
                pitch: h.pitch * scale,
            })
        }
        Curve3::SineWave(s) => {
            let scale = transform.matrix3.x_axis.length();
            Curve3::SineWave(rcad_kernel::geom::SineWave3 {
                origin: transform.transform_point3(s.origin),
                baseline_dir: transform.transform_vector3(s.baseline_dir).normalize(),
                amplitude_dir: transform.transform_vector3(s.amplitude_dir).normalize(),
                amplitude: s.amplitude * scale,
                frequency: s.frequency,
                phase: s.phase,
            })
        }
    }
}

fn transform_bspline_curve(curve: &BSplineCurve3, transform: DAffine3) -> BSplineCurve3 {
    BSplineCurve3 {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: curve
            .control_points
            .iter()
            .map(|&p| transform.transform_point3(p))
            .collect(),
        weights: curve.weights.clone(),
    }
}

fn transform_bezier_curve(curve: &BezierCurve3, transform: DAffine3) -> BezierCurve3 {
    BezierCurve3 {
        control_points: curve
            .control_points
            .iter()
            .map(|&p| transform.transform_point3(p))
            .collect(),
        weights: curve.weights.clone(),
    }
}

// =============================================================================
// Surface Tools
// =============================================================================

/// Reverse the U parametric direction of a surface.
///
/// This creates a new surface where the U parameter runs in the opposite direction.
///
/// # Arguments
/// * `surface` - The surface to reverse
///
/// # Returns
/// A new surface with reversed U direction
pub fn reverse_surface_u(surface: Surface3) -> Surface3 {
    match surface {
        Surface3::BSpline(bspline) => {
            Surface3::BSpline(reverse_bspline_surface_u(&bspline))
        }
        Surface3::Bezier(bezier) => {
            Surface3::Bezier(reverse_bezier_surface_u(&bezier))
        }
        Surface3::Trimmed(mut trimmed) => {
            // Swap u bounds
            let [u1, u2, v1, v2] = trimmed.trim;
            trimmed.trim = [u2, u1, v1, v2];
            Surface3::Trimmed(trimmed)
        }
        _ => surface, // Analytic surfaces: reversing U typically has no geometric effect
    }
}

/// Reverse the V parametric direction of a surface.
///
/// This creates a new surface where the V parameter runs in the opposite direction.
///
/// # Arguments
/// * `surface` - The surface to reverse
///
/// # Returns
/// A new surface with reversed V direction
pub fn reverse_surface_v(surface: Surface3) -> Surface3 {
    match surface {
        Surface3::BSpline(bspline) => {
            Surface3::BSpline(reverse_bspline_surface_v(&bspline))
        }
        Surface3::Bezier(bezier) => {
            Surface3::Bezier(reverse_bezier_surface_v(&bezier))
        }
        Surface3::Trimmed(mut trimmed) => {
            // Swap v bounds
            let [u1, u2, v1, v2] = trimmed.trim;
            trimmed.trim = [u1, u2, v2, v1];
            Surface3::Trimmed(trimmed)
        }
        _ => surface, // Analytic surfaces: reversing V typically has no geometric effect
    }
}

fn reverse_bspline_surface_u(surface: &BSplineSurface) -> BSplineSurface {
    // Reverse control points and weights in U direction
    let new_ctrl: Vec<Vec<DVec3>> = surface.control_points.iter().rev().cloned().collect();
    let new_weights: Vec<Vec<f64>> = surface.weights.iter().rev().cloned().collect();

    // Reverse and renormalize U knot vector
    let knots_u = &surface.knots_u;
    let k0 = knots_u[0];
    let k1 = knots_u[knots_u.len() - 1];

    let new_knots_u: Vec<f64> = knots_u
        .iter()
        .rev()
        .map(|&k| k0 + (k1 - k))
        .collect();

    BSplineSurface {
        degree_u: surface.degree_u,
        degree_v: surface.degree_v,
        knots_u: new_knots_u,
        knots_v: surface.knots_v.clone(),
        control_points: new_ctrl,
        weights: new_weights,
    }
}

fn reverse_bspline_surface_v(surface: &BSplineSurface) -> BSplineSurface {
    // Reverse control points and weights in V direction
    let new_ctrl: Vec<Vec<DVec3>> = surface
        .control_points
        .iter()
        .map(|row| row.iter().rev().cloned().collect())
        .collect();
    let new_weights: Vec<Vec<f64>> = surface
        .weights
        .iter()
        .map(|row| row.iter().rev().cloned().collect())
        .collect();

    // Reverse and renormalize V knot vector
    let knots_v = &surface.knots_v;
    let k0 = knots_v[0];
    let k1 = knots_v[knots_v.len() - 1];

    let new_knots_v: Vec<f64> = knots_v
        .iter()
        .rev()
        .map(|&k| k0 + (k1 - k))
        .collect();

    BSplineSurface {
        degree_u: surface.degree_u,
        degree_v: surface.degree_v,
        knots_u: surface.knots_u.clone(),
        knots_v: new_knots_v,
        control_points: new_ctrl,
        weights: new_weights,
    }
}

fn reverse_bezier_surface_u(surface: &BezierSurface) -> BezierSurface {
    BezierSurface {
        control_points: surface.control_points.iter().rev().cloned().collect(),
        weights: surface.weights.iter().rev().cloned().collect(),
    }
}

fn reverse_bezier_surface_v(surface: &BezierSurface) -> BezierSurface {
    BezierSurface {
        control_points: surface
            .control_points
            .iter()
            .map(|row| row.iter().rev().cloned().collect())
            .collect(),
        weights: surface
            .weights
            .iter()
            .map(|row| row.iter().rev().cloned().collect())
            .collect(),
    }
}

/// Trim a surface to the specified UV bounds.
///
/// This creates a new surface restricted to the parameter range
/// `[u1, u2] x [v1, v2]`.
///
/// # Arguments
/// * `surface` - The surface to trim
/// * `u1` - The minimum U parameter
/// * `u2` - The maximum U parameter
/// * `v1` - The minimum V parameter
/// * `v2` - The maximum V parameter
///
/// # Returns
/// A new surface representing the trimmed region
pub fn trim_surface(surface: &Surface3, u1: f64, u2: f64, v1: f64, v2: f64) -> Surface3 {
    // For BSpline surfaces, we could do exact trimming via knot insertion
    // For now, use the TrimmedSurface wrapper
    Surface3::Trimmed(TrimmedSurface::new(surface.clone(), u1, u2, v1, v2))
}

/// Apply an affine transformation to a surface.
///
/// # Arguments
/// * `surface` - The surface to transform
/// * `transform` - The affine transformation to apply
///
/// # Returns
/// A new surface with the transformation applied
pub fn transform_surface(surface: &Surface3, transform: DAffine3) -> Surface3 {
    match surface {
        Surface3::Plane(plane) => {
            Surface3::Plane(Plane::new(transform.transform_point3(plane.origin), transform.transform_vector3(plane.normal).normalize()))
        }
        Surface3::Cylinder(cyl) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Cylinder(CylindricalSurface {
                origin: transform.transform_point3(cyl.origin),
                axis: transform.transform_vector3(cyl.axis).normalize(),
                radius: cyl.radius * scale,
                ref_dir: transform.transform_vector3(cyl.ref_dir).normalize(),
            })
        }
        Surface3::Sphere(sphere) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Sphere(SphericalSurface {
                center: transform.transform_point3(sphere.center),
                axis: transform.transform_vector3(sphere.axis).normalize(),
                radius: sphere.radius * scale,
                ref_dir: any_perpendicular(transform.transform_vector3(sphere.axis).normalize()),
            })
        }
        Surface3::Cone(cone) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Cone(ConicalSurface {
                apex: transform.transform_point3(cone.apex),
                axis: transform.transform_vector3(cone.axis).normalize(),
                radius: cone.radius * scale,
                half_angle_rad: cone.half_angle_rad,
            })
        }
        Surface3::Torus(torus) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Torus(ToroidalSurface {
                center: transform.transform_point3(torus.center),
                axis: transform.transform_vector3(torus.axis).normalize(),
                major_radius: torus.major_radius * scale,
                minor_radius: torus.minor_radius * scale,
            })
        }
        Surface3::BSpline(bspline) => {
            Surface3::BSpline(transform_bspline_surface(bspline, transform))
        }
        Surface3::Bezier(bezier) => {
            Surface3::Bezier(transform_bezier_surface(bezier, transform))
        }
        Surface3::Trimmed(trimmed) => {
            Surface3::Trimmed(TrimmedSurface {
                basis: Box::new(transform_surface(&trimmed.basis, transform)),
                trim: trimmed.trim,
            })
        }
        Surface3::Offset(offset) => {
            let mut new_offset = offset.clone();
            new_offset.basis = Box::new(transform_surface(&new_offset.basis, transform));
            Surface3::Offset(new_offset)
        }
        Surface3::LinearExtrusion(extrusion) => {
            Surface3::LinearExtrusion(rcad_kernel::geom::LinearExtrusionSurface {
                profile: Box::new(transform_curve(&extrusion.profile, transform)),
                direction: transform.transform_vector3(extrusion.direction).normalize(),
            })
        }
        Surface3::Revolution(revolution) => {
            Surface3::Revolution(rcad_kernel::geom::RevolutionSurface {
                profile: Box::new(transform_curve(&revolution.profile, transform)),
                axis_origin: transform.transform_point3(revolution.axis_origin),
                axis_dir: transform.transform_vector3(revolution.axis_dir).normalize(),
            })
        }
        Surface3::Ruled(ruled) => {
            Surface3::Ruled(rcad_kernel::geom::RuledSurface {
                start: Box::new(transform_curve(&ruled.start, transform)),
                end: Box::new(transform_curve(&ruled.end, transform)),
            })
        }
        Surface3::Coons(coons) => {
            Surface3::Coons(rcad_kernel::geom::CoonsSurface {
                south: Box::new(transform_curve(&coons.south, transform)),
                north: Box::new(transform_curve(&coons.north, transform)),
                west: Box::new(transform_curve(&coons.west, transform)),
                east: Box::new(transform_curve(&coons.east, transform)),
            })
        }
        Surface3::Ellipsoid(ell) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Ellipsoid(rcad_kernel::geom::EllipsoidalSurface {
                center: transform.transform_point3(ell.center),
                axis: transform.transform_vector3(ell.axis).normalize(),
                ref_dir: transform.transform_vector3(ell.ref_dir).normalize(),
                radius_x: ell.radius_x * scale,
                radius_y: ell.radius_y * scale,
                radius_z: ell.radius_z * scale,
            })
        }
        Surface3::Helicoid(h) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Helicoid(rcad_kernel::geom::HelicoidSurface {
                origin: transform.transform_point3(h.origin),
                axis: transform.transform_vector3(h.axis).normalize(),
                ref_dir: transform.transform_vector3(h.ref_dir).normalize(),
                pitch: h.pitch * scale,
            })
        }
        Surface3::Pipe(pipe) => {
            let scale = transform.matrix3.x_axis.length();
            Surface3::Pipe(rcad_kernel::geom::PipeSurface {
                spine: Box::new(transform_curve(&pipe.spine, transform)),
                ref_dir: transform.transform_vector3(pipe.ref_dir).normalize(),
                radius: pipe.radius * scale,
            })
        }
        Surface3::TriBezier(tri) => {
            Surface3::TriBezier(transform_tri_bezier_surface(tri, transform))
        }
    }
}

fn transform_bspline_surface(surface: &BSplineSurface, transform: DAffine3) -> BSplineSurface {
    BSplineSurface {
        degree_u: surface.degree_u,
        degree_v: surface.degree_v,
        knots_u: surface.knots_u.clone(),
        knots_v: surface.knots_v.clone(),
        control_points: surface
            .control_points
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&p| transform.transform_point3(p))
                    .collect()
            })
            .collect(),
        weights: surface.weights.clone(),
    }
}

fn transform_bezier_surface(surface: &BezierSurface, transform: DAffine3) -> BezierSurface {
    BezierSurface {
        control_points: surface
            .control_points
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&p| transform.transform_point3(p))
                    .collect()
            })
            .collect(),
        weights: surface.weights.clone(),
    }
}

fn transform_tri_bezier_surface(
    surface: &rcad_kernel::geom::TriBezierSurface,
    transform: DAffine3,
) -> rcad_kernel::geom::TriBezierSurface {
    rcad_kernel::geom::TriBezierSurface {
        control_points: surface
            .control_points
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&p| transform.transform_point3(p))
                    .collect()
            })
            .collect(),
        weights: surface.weights.clone(),
    }
}

// =============================================================================
// Continuity Checking
// =============================================================================

/// Check the continuity order of a curve.
///
/// Returns the highest continuity level achieved:
/// - 0: C0 (positional continuity)
/// - 1: C1 (tangent continuity)
/// - 2: C2 (curvature continuity)
/// - Higher values for smoother curves
///
/// # Arguments
/// * `curve` - The curve to check
/// * `tol` - The tolerance for checking continuity
///
/// # Returns
/// The continuity order (0, 1, 2, etc.)
pub fn check_curve_continuity(curve: &Curve3, tol: f64) -> usize {
    match curve {
        Curve3::Line(_) => usize::MAX, // Lines are infinitely smooth
        Curve3::Circle(_) => usize::MAX, // Circles are infinitely smooth
        Curve3::Ellipse(_) => usize::MAX, // Ellipses are infinitely smooth
        Curve3::Hyperbola(_) => usize::MAX,
        Curve3::Parabola(_) => usize::MAX,
        Curve3::BSpline(bspline) => check_bspline_curve_continuity(bspline, tol),
        Curve3::Bezier(_) => usize::MAX, // Bezier curves are C-infinity within their domain
        Curve3::Offset(offset) => check_curve_continuity(&offset.basis, tol),
        Curve3::CircularHelix(_) => usize::MAX,
        Curve3::SineWave(_) => usize::MAX,
    }
}

fn check_bspline_curve_continuity(curve: &BSplineCurve3, tol: f64) -> usize {
    // Check the maximum continuity based on knot multiplicities
    // A BSpline curve is C^(degree - multiplicity) at each knot
    let degree = curve.degree;
    let knots = &curve.knots;

    if knots.len() < 2 {
        return 0;
    }

    // Find internal knots and their multiplicities
    let mut min_continuity = degree; // Best case

    let mut i = 0;
    while i < knots.len() {
        let knot = knots[i];
        let mut multiplicity = 1;
        while i + multiplicity < knots.len() && (knots[i + multiplicity] - knot).abs() < tol {
            multiplicity += 1;
        }

        // Only consider internal knots (not the endpoints)
        if i > degree && i + multiplicity < knots.len() - degree {
            let continuity = degree.saturating_sub(multiplicity);
            min_continuity = min_continuity.min(continuity);
        }

        i += multiplicity;
    }

    min_continuity
}

/// Check the continuity orders of a surface.
///
/// Returns a tuple `(u_continuity, v_continuity)` where each value represents
/// the highest continuity level achieved in the respective parametric direction.
///
/// # Arguments
/// * `surface` - The surface to check
/// * `tol` - The tolerance for checking continuity
///
/// # Returns
/// A tuple of (U continuity, V continuity)
pub fn check_surface_continuity(surface: &Surface3, tol: f64) -> (usize, usize) {
    match surface {
        Surface3::Plane(_) => (usize::MAX, usize::MAX),
        Surface3::Cylinder(_) => (usize::MAX, usize::MAX),
        Surface3::Sphere(_) => (usize::MAX, usize::MAX),
        Surface3::Cone(_) => (usize::MAX, usize::MAX),
        Surface3::Torus(_) => (usize::MAX, usize::MAX),
        Surface3::Ellipsoid(_) => (usize::MAX, usize::MAX),
        Surface3::BSpline(bspline) => check_bspline_surface_continuity(bspline, tol),
        Surface3::Bezier(_) => (usize::MAX, usize::MAX),
        Surface3::Trimmed(trimmed) => check_surface_continuity(&trimmed.basis, tol),
        Surface3::Offset(offset) => check_surface_continuity(&offset.basis, tol),
        Surface3::LinearExtrusion(extrusion) => {
            let u_cont = check_curve_continuity(&extrusion.profile, tol);
            (u_cont, usize::MAX)
        }
        Surface3::Revolution(revolution) => {
            let v_cont = check_curve_continuity(&revolution.profile, tol);
            (usize::MAX, v_cont)
        }
        Surface3::Ruled(ruled) => {
            let cont1 = check_curve_continuity(&ruled.start, tol);
            let cont2 = check_curve_continuity(&ruled.end, tol);
            (cont1.min(cont2), usize::MAX)
        }
        Surface3::Coons(_) => (0, 0), // Coons patches are typically C0 at boundaries
        Surface3::Helicoid(_) => (usize::MAX, usize::MAX),
        Surface3::Pipe(pipe) => {
            let v_cont = check_curve_continuity(&pipe.spine, tol);
            (usize::MAX, v_cont)
        }
        Surface3::TriBezier(_) => (usize::MAX, usize::MAX),
    }
}

fn check_bspline_surface_continuity(surface: &BSplineSurface, tol: f64) -> (usize, usize) {
    let degree_u = surface.degree_u;
    let degree_v = surface.degree_v;

    let u_cont = check_knot_continuity(&surface.knots_u, degree_u, tol);
    let v_cont = check_knot_continuity(&surface.knots_v, degree_v, tol);

    (u_cont, v_cont)
}

fn check_knot_continuity(knots: &[f64], degree: usize, tol: f64) -> usize {
    if knots.len() < 2 {
        return 0;
    }

    let mut min_continuity = degree;

    let mut i = 0;
    while i < knots.len() {
        let knot = knots[i];
        let mut multiplicity = 1;
        while i + multiplicity < knots.len() && (knots[i + multiplicity] - knot).abs() < tol {
            multiplicity += 1;
        }

        // Only consider internal knots
        if i > degree && i + multiplicity < knots.len() - degree {
            let continuity = degree.saturating_sub(multiplicity);
            min_continuity = min_continuity.min(continuity);
        }

        i += multiplicity;
    }

    min_continuity
}

// =============================================================================
// Tests
// =============================================================================


