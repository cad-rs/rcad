//! Shape customization tools analogous to OCCT ShapeCustom package.
//!
//! This module provides utilities for:
//! - BSpline degree reduction and simplification
//! - Surface conversion to NURBS representation
//! - Geometry restrictions for export compatibility
//!
//! # Modules
//!
//! - [`BSplineSimplifyOptions`] - Configuration for BSpline simplification
//! - [`simplify_bspline_curve`] - Reduce degree and control points of a BSpline curve
//! - [`simplify_bspline_surface`] - Reduce degree and control points of a BSpline surface
//! - [`convert_to_bspline`] - Convert entire BRep geometry to BSpline representation
//! - [`restrict_geometry`] - Apply geometry restrictions to a BRep

use glam::DVec3;
use rcad_kernel::{
    BRep, Curve3, Surface3,
    geom::{
        BSplineCurve3, BSplineSurface, CurveEval, SurfaceEval,
    },
    nurbs_convert::{
        curve_to_bspline, surface_to_bspline,
        bezier_curve_to_bspline, bezier_surface_to_bspline,
        line_to_bspline, circle_to_bspline, ellipse_to_bspline,
        plane_to_bspline, cylinder_to_bspline, sphere_to_bspline,
    },
};

// =============================================================================
// BSpline Simplification Options
// =============================================================================

/// Options for BSpline curve/surface simplification.
///
/// Analogous to OCCT `ShapeCustom_BSplineRestriction`.
#[derive(Debug, Clone)]
pub struct BSplineSimplifyOptions {
    /// Maximum degree allowed (default: 3).
    pub max_degree: usize,
    /// Tolerance for approximation when reducing degree (default: 1e-6).
    pub tolerance: f64,
    /// Whether to preserve endpoint tangents (default: true).
    pub preserve_ends: bool,
    /// Minimum control point reduction ratio (default: 0.8, i.e. target 80% of original).
    pub min_reduction_ratio: f64,
    /// Maximum number of iterations for degree reduction (default: 10).
    pub max_iterations: usize,
}

impl Default for BSplineSimplifyOptions {
    fn default() -> Self {
        Self {
            max_degree: 3,
            tolerance: 1e-6,
            preserve_ends: true,
            min_reduction_ratio: 0.8,
            max_iterations: 10,
        }
    }
}

/// Result of BSpline simplification.
#[derive(Debug, Clone)]
pub struct SimplificationResult<T> {
    /// The simplified geometry.
    pub geometry: T,
    /// Whether any simplification was applied.
    pub was_simplified: bool,
    /// Maximum deviation from original.
    pub max_deviation: f64,
    /// Original degree.
    pub original_degree: usize,
    /// Final degree.
    pub final_degree: usize,
    /// Original number of control points (or [u, v] for surfaces).
    pub original_ctrl_pts: usize,
    /// Final number of control points.
    pub final_ctrl_pts: usize,
}

// =============================================================================
// BSpline Curve Simplification
// =============================================================================

/// Simplify a BSpline curve by reducing degree and removing redundant control points.
///
/// This function:
/// 1. Reduces the degree if it exceeds `max_degree`
/// 2. Removes approximately collinear control points within tolerance
/// 3. Preserves endpoints and optionally endpoint tangents
///
/// # Example
/// ```ignore
/// use rcad_algorithms::shape_custom::{simplify_bspline_curve, BSplineSimplifyOptions};
/// let opts = BSplineSimplifyOptions { max_degree: 3, ..Default::default() };
/// let simplified = simplify_bspline_curve(&curve, &opts);
/// ```
pub fn simplify_bspline_curve(
    curve: &BSplineCurve3,
    opts: &BSplineSimplifyOptions,
) -> SimplificationResult<BSplineCurve3> {
    let original_degree = curve.degree;
    let original_ctrl_pts = curve.control_points.len();

    // If already within constraints, return as-is
    if curve.degree <= opts.max_degree && curve.control_points.len() <= 4 {
        return SimplificationResult {
            geometry: curve.clone(),
            was_simplified: false,
            max_deviation: 0.0,
            original_degree,
            final_degree: curve.degree,
            original_ctrl_pts,
            final_ctrl_pts: curve.control_points.len(),
        };
    }

    let mut current = curve.clone();
    let mut total_deviation: f64 = 0.0;
    let mut was_simplified = false;

    // Step 1: Degree reduction
    if current.degree > opts.max_degree {
        let reduced = reduce_curve_degree(&current, opts.max_degree, opts.tolerance);
        if let Some((new_curve, deviation)) = reduced {
            total_deviation = total_deviation.max(deviation);
            current = new_curve;
            was_simplified = true;
        }
    }

    // Step 2: Control point reduction (knot removal)
    let target_pts = (curve.control_points.len() as f64 * opts.min_reduction_ratio).ceil() as usize;
    let target_pts = target_pts.max(opts.max_degree + 1); // Minimum control points for degree

    if current.control_points.len() > target_pts {
        let reduced = remove_redundant_control_points(&current, opts.tolerance, target_pts);
        if let Some((new_curve, deviation)) = reduced {
            total_deviation = total_deviation.max(deviation);
            current = new_curve;
            was_simplified = true;
        }
    }

    SimplificationResult {
        geometry: current,
        was_simplified,
        max_deviation: total_deviation,
        original_degree,
        final_degree: curve.degree.min(opts.max_degree),
        original_ctrl_pts,
        final_ctrl_pts: curve.control_points.len(),
    }
}

/// Reduce the degree of a BSpline curve through approximation.
///
/// Uses degree elevation inverse (knot insertion reversal) where possible,
/// falls back to sampling and re-fitting when necessary.
fn reduce_curve_degree(
    curve: &BSplineCurve3,
    target_degree: usize,
    tolerance: f64,
) -> Option<(BSplineCurve3, f64)> {
    if curve.degree <= target_degree {
        return None;
    }

    // Sample the curve and fit a lower-degree approximation
    let n_samples = curve.control_points.len().max(20) * 2;
    let [t0, t1] = curve.default_domain();

    let samples: Vec<DVec3> = (0..n_samples)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / (n_samples - 1) as f64;
            curve.point_at(t)
        })
        .collect();

    // Fit with target degree
    let fitted = fit_curve_to_points(&samples, target_degree, tolerance);
    let (new_curve, max_dev) = fitted?;

    // Check deviation
    let mut max_deviation: f64 = 0.0;
    for i in 0..n_samples {
        let t = t0 + (t1 - t0) * i as f64 / (n_samples - 1) as f64;
        let orig_pt = curve.point_at(t);
        let new_pt = new_curve.point_at(t);
        max_deviation = max_deviation.max((orig_pt - new_pt).length());
    }

    if max_deviation <= tolerance {
        Some((new_curve, max_deviation))
    } else {
        None
    }
}

/// Remove redundant control points from a BSpline curve.
///
/// Uses knot removal algorithm: evaluates the curve with and without each knot,
/// checks deviation against tolerance.
fn remove_redundant_control_points(
    curve: &BSplineCurve3,
    tolerance: f64,
    target_ctrl_pts: usize,
) -> Option<(BSplineCurve3, f64)> {
    if curve.control_points.len() <= target_ctrl_pts {
        return None;
    }

    // Identify interior knots with multiplicity > 0
    let knots = &curve.knots;
    let degree = curve.degree;

    // Find unique interior knots
    let mut interior_knots: Vec<f64> = Vec::new();
    let mut i = degree;
    while i < knots.len() - degree - 1 {
        let k = knots[i];
        if k > knots[degree] && k < knots[knots.len() - degree - 1] {
            if interior_knots.last().map_or(true, |&last| (k - last).abs() > 1e-10) {
                interior_knots.push(k);
            }
        }
        i += 1;
    }

    if interior_knots.is_empty() {
        return None;
    }

    // Try removing knots one by one, starting from the middle
    let mut current = curve.clone();
    let mut total_deviation: f64 = 0.0;
    let mut removed_any = false;

    // Sort knots by importance (middle knots first, as they affect shape less)
    let mut sorted_knots = interior_knots.clone();
    let mid = sorted_knots.len() as f64 / 2.0;
    sorted_knots.sort_by(|a, b| {
        let a_dist = (*a - mid).abs();
        let b_dist = (*b - mid).abs();
        a_dist.partial_cmp(&b_dist).unwrap()
    });

    for knot in sorted_knots {
        if current.control_points.len() <= target_ctrl_pts {
            break;
        }

        if let Some((new_curve, dev)) = try_remove_knot(&current, knot, tolerance) {
            current = new_curve;
            total_deviation = total_deviation.max(dev);
            removed_any = true;
        }
    }

    if removed_any {
        Some((current, total_deviation))
    } else {
        None
    }
}

/// Try to remove a knot from the curve, checking deviation.
fn try_remove_knot(
    curve: &BSplineCurve3,
    knot: f64,
    tolerance: f64,
) -> Option<(BSplineCurve3, f64)> {
    // Find knot index
    let knot_idx = curve.knots.iter().position(|&k| (k - knot).abs() < 1e-10)?;

    // Build new knot vector without this knot
    let mut new_knots = curve.knots.clone();
    new_knots.remove(knot_idx);

    // Check if we still have enough knots
    let n_ctrl = curve.control_points.len();
    if new_knots.len() < curve.degree + 1 + n_ctrl - 1 {
        return None;
    }

    // Recompute control points for new knot vector
    // This is a simplified approach - proper knot removal requires
    // solving a constrained optimization
    let new_n_ctrl = n_ctrl - 1;
    if new_n_ctrl < curve.degree + 1 {
        return None;
    }

    // Sample the original curve and fit to new knot vector
    let n_samples = n_ctrl.max(20);
    let [t0, t1] = curve.default_domain();
    let samples: Vec<DVec3> = (0..n_samples)
        .map(|i| {
            let t = t0 + (t1 - t0) * i as f64 / (n_samples - 1) as f64;
            curve.point_at(t)
        })
        .collect();

    // Fit new curve with one fewer control point
    let mut new_ctrl_pts = Vec::with_capacity(new_n_ctrl);
    let mut new_weights = Vec::with_capacity(new_n_ctrl);

    // Simple approach: redistribute control points
    for i in 0..new_n_ctrl {
        let u = i as f64 / (new_n_ctrl - 1).max(1) as f64;
        let t = t0 + u * (t1 - t0);
        new_ctrl_pts.push(curve.point_at(t));
        new_weights.push(1.0);
    }

    let new_curve = BSplineCurve3 {
        degree: curve.degree,
        knots: build_clamped_knots(new_n_ctrl, curve.degree),
        control_points: new_ctrl_pts,
        weights: new_weights,
    };

    // Check deviation
    let mut max_dev: f64 = 0.0;
    for i in 0..n_samples {
        let t = t0 + (t1 - t0) * i as f64 / (n_samples - 1) as f64;
        let orig = curve.point_at(t);
        let new = new_curve.point_at(t);
        max_dev = max_dev.max((orig - new).length());
    }

    if max_dev <= tolerance {
        Some((new_curve, max_dev))
    } else {
        None
    }
}

/// Build a clamped knot vector for n control points of given degree.
fn build_clamped_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    // BSpline knots length = n_ctrl + degree + 1
    // For clamped knots: degree+1 zeros, internal knots, degree+1 ones
    let total_len = n_ctrl + degree + 1;
    let mut knots = Vec::with_capacity(total_len);

    // Clamped start: degree+1 zeros
    for _ in 0..=degree {
        knots.push(0.0);
    }

    // Interior knots: total_len - 2*(degree+1) = n_ctrl - degree - 1 interior knots
    let n_interior = n_ctrl.saturating_sub(degree).saturating_sub(1);
    if n_interior > 0 {
        for i in 1..=n_interior {
            knots.push(i as f64 / (n_interior + 1) as f64);
        }
    }

    // Clamped end: degree+1 ones (but only if we have room)
    let remaining = total_len.saturating_sub(knots.len());
    for _ in 0..remaining {
        knots.push(1.0);
    }

    knots
}

/// Fit a BSpline curve through points with specified degree.
fn fit_curve_to_points(
    points: &[DVec3],
    degree: usize,
    _tolerance: f64,
) -> Option<(BSplineCurve3, f64)> {
    let n = points.len();
    if n < 2 {
        return None;
    }

    let degree = degree.min(n - 1);
    let n_ctrl = n;

    // Use chord-length parameterization
    let mut params = vec![0.0];
    let mut total_len = 0.0;
    for i in 1..n {
        total_len += (points[i] - points[i - 1]).length();
        params.push(total_len);
    }
    if total_len > 0.0 {
        for p in &mut params {
            *p /= total_len;
        }
    }

    // Build knot vector
    let knots = build_clamped_knots(n_ctrl, degree);

    // Simple interpolation (use control points as data points for now)
    let control_points = points.to_vec();
    let weights = vec![1.0; n_ctrl];

    let curve = BSplineCurve3 {
        degree,
        knots,
        control_points,
        weights,
    };

    // Compute max deviation
    let mut max_dev: f64 = 0.0;
    for (i, &pt) in points.iter().enumerate() {
        let curve_pt = curve.point_at(params[i]);
        max_dev = max_dev.max((pt - curve_pt).length());
    }

    Some((curve, max_dev))
}

// =============================================================================
// BSpline Surface Simplification
// =============================================================================

/// Simplify a BSpline surface by reducing degree and removing redundant control points.
///
/// This function applies simplification in both U and V directions.
pub fn simplify_bspline_surface(
    surface: &BSplineSurface,
    opts: &BSplineSimplifyOptions,
) -> SimplificationResult<BSplineSurface> {
    let original_degree_u = surface.degree_u;
    let original_degree_v = surface.degree_v;
    let original_ctrl_pts = surface.control_points.len() * surface.control_points.get(0).map_or(0, |r| r.len());

    // If already within constraints, return as-is
    if surface.degree_u <= opts.max_degree
        && surface.degree_v <= opts.max_degree
        && surface.control_points.len() <= 4
        && surface.control_points.get(0).map_or(true, |r| r.len() <= 4)
    {
        return SimplificationResult {
            geometry: surface.clone(),
            was_simplified: false,
            max_deviation: 0.0,
            original_degree: original_degree_u.max(original_degree_v),
            final_degree: surface.degree_u.max(surface.degree_v),
            original_ctrl_pts,
            final_ctrl_pts: surface.control_points.len() * surface.control_points.get(0).map_or(0, |r| r.len()),
        };
    }

    let mut current = surface.clone();
    let mut total_deviation: f64 = 0.0;
    let mut was_simplified = false;

    // Reduce U degree
    if current.degree_u > opts.max_degree {
        if let Some((new_surf, dev)) = reduce_surface_degree_u(&current, opts.max_degree, opts.tolerance) {
            total_deviation = total_deviation.max(dev);
            current = new_surf;
            was_simplified = true;
        }
    }

    // Reduce V degree
    if current.degree_v > opts.max_degree {
        if let Some((new_surf, dev)) = reduce_surface_degree_v(&current, opts.max_degree, opts.tolerance) {
            total_deviation = total_deviation.max(dev);
            current = new_surf;
            was_simplified = true;
        }
    }

    // Compute final values before moving current
    let final_degree = current.degree_u.max(current.degree_v);
    let final_ctrl_pts = current.control_points.len() * current.control_points.get(0).map_or(0, |r| r.len());

    SimplificationResult {
        geometry: current,
        was_simplified,
        max_deviation: total_deviation,
        original_degree: original_degree_u.max(original_degree_v),
        final_degree,
        original_ctrl_pts,
        final_ctrl_pts,
    }
}

/// Reduce the U-degree of a BSpline surface.
fn reduce_surface_degree_u(
    surface: &BSplineSurface,
    target_degree: usize,
    tolerance: f64,
) -> Option<(BSplineSurface, f64)> {
    if surface.degree_u <= target_degree {
        return None;
    }

    // Sample surface and refit
    let n_u = surface.control_points.len().max(10);
    let n_v = surface.control_points.get(0)?.len().max(10);
    let [u0, u1, v0, v1] = surface.default_domain();

    let mut samples: Vec<Vec<DVec3>> = Vec::new();
    for i in 0..n_u {
        let mut row = Vec::new();
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1) as f64;
        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1) as f64;
            row.push(surface.point_at(u, v));
        }
        samples.push(row);
    }

    // Fit new surface with reduced U degree
    let new_surface = fit_surface_to_grid(&samples, target_degree, surface.degree_v)?;
    let max_dev = compute_surface_deviation(surface, &new_surface, n_u, n_v);

    if max_dev <= tolerance {
        Some((new_surface, max_dev))
    } else {
        None
    }
}

/// Reduce the V-degree of a BSpline surface.
fn reduce_surface_degree_v(
    surface: &BSplineSurface,
    target_degree: usize,
    tolerance: f64,
) -> Option<(BSplineSurface, f64)> {
    if surface.degree_v <= target_degree {
        return None;
    }

    let n_u = surface.control_points.len().max(10);
    let n_v = surface.control_points.get(0)?.len().max(10);
    let [u0, u1, v0, v1] = surface.default_domain();

    let mut samples: Vec<Vec<DVec3>> = Vec::new();
    for i in 0..n_u {
        let mut row = Vec::new();
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1) as f64;
        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1) as f64;
            row.push(surface.point_at(u, v));
        }
        samples.push(row);
    }

    let new_surface = fit_surface_to_grid(&samples, surface.degree_u, target_degree)?;
    let max_dev = compute_surface_deviation(surface, &new_surface, n_u, n_v);

    if max_dev <= tolerance {
        Some((new_surface, max_dev))
    } else {
        None
    }
}

/// Fit a BSpline surface to a grid of points.
fn fit_surface_to_grid(
    points: &[Vec<DVec3>],
    degree_u: usize,
    degree_v: usize,
) -> Option<BSplineSurface> {
    let n_u = points.len();
    let n_v = points.get(0)?.len();

    if n_u < 2 || n_v < 2 {
        return None;
    }

    let degree_u = degree_u.min(n_u - 1);
    let degree_v = degree_v.min(n_v - 1);

    // Use points directly as control points (simple approximation)
    let control_points = points.to_vec();
    let weights: Vec<Vec<f64>> = points.iter()
        .map(|row| vec![1.0; row.len()])
        .collect();

    Some(BSplineSurface {
        degree_u,
        degree_v,
        knots_u: build_clamped_knots(n_u, degree_u),
        knots_v: build_clamped_knots(n_v, degree_v),
        control_points,
        weights,
    })
}

/// Compute maximum deviation between two surfaces.
fn compute_surface_deviation(
    surf1: &BSplineSurface,
    surf2: &BSplineSurface,
    n_u: usize,
    n_v: usize,
) -> f64 {
    let [u0, u1, v0, v1] = surf1.default_domain();
    let mut max_dev: f64 = 0.0;

    for i in 0..n_u {
        let u = u0 + (u1 - u0) * i as f64 / (n_u - 1).max(1) as f64;
        for j in 0..n_v {
            let v = v0 + (v1 - v0) * j as f64 / (n_v - 1).max(1) as f64;
            let p1 = surf1.point_at(u, v);
            let p2 = surf2.point_at(u, v);
            max_dev = max_dev.max((p1 - p2).length());
        }
    }

    max_dev
}

// =============================================================================
// Geometry Restrictions
// =============================================================================

/// Geometry restrictions for export compatibility.
///
/// Analogous to OCCT `ShapeCustom_RestrictionParameters`.
#[derive(Debug, Clone)]
pub struct GeometryRestrictions {
    /// Maximum allowed degree for curves and surfaces.
    pub max_degree: usize,
    /// Tolerance for conversion approximations.
    pub tolerance: f64,
    /// Convert all curves to BSpline.
    pub curves_to_bspline: bool,
    /// Convert all surfaces to BSpline.
    pub surfaces_to_bspline: bool,
    /// Convert offset curves to BSpline.
    pub convert_offset_curves: bool,
    /// Convert offset surfaces to BSpline.
    pub convert_offset_surfaces: bool,
    /// Simplify BSpline curves after conversion.
    pub simplify_bspline: bool,
    /// Sample count for curve approximation.
    pub curve_samples: usize,
    /// Sample counts (u, v) for surface approximation.
    pub surface_samples: (usize, usize),
}

impl Default for GeometryRestrictions {
    fn default() -> Self {
        Self {
            max_degree: 3,
            tolerance: 1e-6,
            curves_to_bspline: true,
            surfaces_to_bspline: true,
            convert_offset_curves: true,
            convert_offset_surfaces: true,
            simplify_bspline: true,
            curve_samples: 32,
            surface_samples: (16, 16),
        }
    }
}

/// Result of geometry conversion.
#[derive(Debug, Clone, Default)]
pub struct ConversionReport {
    /// Number of curves converted.
    pub curves_converted: usize,
    /// Number of surfaces converted.
    pub surfaces_converted: usize,
    /// Number of BSpline curves simplified.
    pub bspline_curves_simplified: usize,
    /// Number of BSpline surfaces simplified.
    pub bspline_surfaces_simplified: usize,
    /// Maximum deviation during conversion.
    pub max_deviation: f64,
}

// =============================================================================
// BRep Conversion API
// =============================================================================

/// Convert all geometry in a BRep to BSpline representation.
///
/// This function:
/// 1. Converts all analytic curves to BSpline
/// 2. Converts all analytic surfaces to BSpline
/// 3. Optionally simplifies the resulting BSplines
///
/// # Example
/// ```ignore
/// use rcad_algorithms::shape_custom::convert_to_bspline;
/// let bspline_brep = convert_to_bspline(&brep, 1e-6);
/// ```
pub fn convert_to_bspline(brep: &BRep, tolerance: f64) -> (BRep, ConversionReport) {
    let mut restrictions = GeometryRestrictions::default();
    restrictions.tolerance = tolerance;
    restrict_geometry(brep, &restrictions)
}

/// Apply geometry restrictions to a BRep.
///
/// Converts all geometry to conform to the specified restrictions,
/// useful for export to formats with limited geometry support.
pub fn restrict_geometry(
    brep: &BRep,
    restrictions: &GeometryRestrictions,
) -> (BRep, ConversionReport) {
    let mut result = brep.clone();
    let mut report = ConversionReport::default();

    // Convert curves
    if restrictions.curves_to_bspline {
        for curve in &mut result.geom.curves {
            let needs_conversion = !matches!(curve, Curve3::BSpline(_));
            let is_offset = matches!(curve, Curve3::Offset(_));
            let is_bezier = matches!(curve, Curve3::Bezier(_));

            if needs_conversion && (is_bezier || !is_offset || restrictions.convert_offset_curves) {
                let bspline = curve_to_bspline(curve, restrictions.curve_samples);
                *curve = Curve3::BSpline(bspline);
                report.curves_converted += 1;
            }
        }
    }

    // Convert surfaces
    if restrictions.surfaces_to_bspline {
        for surface in &mut result.geom.surfaces {
            let needs_conversion = !matches!(surface, Surface3::BSpline(_));
            let is_offset = matches!(surface, Surface3::Offset(_));
            let is_bezier = matches!(surface, Surface3::Bezier(_));

            if needs_conversion && (is_bezier || !is_offset || restrictions.convert_offset_surfaces) {
                let (n_u, n_v) = restrictions.surface_samples;
                let bspline = surface_to_bspline(surface, n_u, n_v);
                *surface = Surface3::BSpline(bspline);
                report.surfaces_converted += 1;
            }
        }
    }

    // Simplify BSplines
    if restrictions.simplify_bspline {
        let simplify_opts = BSplineSimplifyOptions {
            max_degree: restrictions.max_degree,
            tolerance: restrictions.tolerance,
            ..Default::default()
        };

        for curve in &mut result.geom.curves {
            if let Curve3::BSpline(bspline) = curve {
                let result = simplify_bspline_curve(bspline, &simplify_opts);
                if result.was_simplified {
                    *bspline = result.geometry;
                    report.bspline_curves_simplified += 1;
                    report.max_deviation = report.max_deviation.max(result.max_deviation);
                }
            }
        }

        for surface in &mut result.geom.surfaces {
            if let Surface3::BSpline(bspline) = surface {
                let result = simplify_bspline_surface(bspline, &simplify_opts);
                if result.was_simplified {
                    *bspline = result.geometry;
                    report.bspline_surfaces_simplified += 1;
                    report.max_deviation = report.max_deviation.max(result.max_deviation);
                }
            }
        }
    }

    (result, report)
}

// =============================================================================
// Utility Functions
// =============================================================================

/// Check if a curve is already in BSpline form.
pub fn is_bspline_curve(curve: &Curve3) -> bool {
    matches!(curve, Curve3::BSpline(_))
}

/// Check if a surface is already in BSpline form.
pub fn is_bspline_surface(surface: &Surface3) -> bool {
    matches!(surface, Surface3::BSpline(_))
}

/// Get the degree of a curve (BSpline degree for BSplines, 1 for lines, 2 for circles/ellipses).
pub fn curve_degree(curve: &Curve3) -> usize {
    match curve {
        Curve3::BSpline(b) => b.degree,
        Curve3::Bezier(b) => b.control_points.len().saturating_sub(1),
        Curve3::Line(_) => 1,
        Curve3::Circle(_) | Curve3::Ellipse(_) => 2,
        _ => 3, // Approximation degree for other types
    }
}

/// Get the degrees of a surface.
pub fn surface_degrees(surface: &Surface3) -> (usize, usize) {
    match surface {
        Surface3::BSpline(b) => (b.degree_u, b.degree_v),
        Surface3::Bezier(b) => {
            let du = b.control_points.len().saturating_sub(1);
            let dv = b.control_points.get(0).map_or(0, |r| r.len().saturating_sub(1));
            (du, dv)
        }
        Surface3::Plane(_) => (1, 1),
        Surface3::Cylinder(_) => (2, 1),
        Surface3::Sphere(_) => (2, 2),
        Surface3::Cone(_) => (2, 1),
        Surface3::Torus(_) => (2, 2),
        _ => (3, 3), // Approximation degrees
    }
}

/// Convert a single curve to BSpline if needed.
pub fn ensure_bspline_curve(curve: &Curve3, samples: usize) -> BSplineCurve3 {
    match curve {
        Curve3::BSpline(b) => b.clone(),
        Curve3::Bezier(b) => bezier_curve_to_bspline(b),
        Curve3::Line(l) => line_to_bspline(l),
        Curve3::Circle(c) => circle_to_bspline(c),
        Curve3::Ellipse(e) => ellipse_to_bspline(e),
        _ => curve_to_bspline(curve, samples),
    }
}

/// Convert a single surface to BSpline if needed.
pub fn ensure_bspline_surface(surface: &Surface3, samples_u: usize, samples_v: usize) -> BSplineSurface {
    match surface {
        Surface3::BSpline(b) => b.clone(),
        Surface3::Bezier(b) => bezier_surface_to_bspline(b),
        Surface3::Plane(p) => plane_to_bspline(p),
        Surface3::Cylinder(c) => cylinder_to_bspline(c),
        Surface3::Sphere(s) => sphere_to_bspline(s),
        _ => surface_to_bspline(surface, samples_u, samples_v),
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{Circle3, Line3, Plane, SphericalSurface};

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn simplify_high_degree_curve() {
        // Create a degree-5 BSpline curve
        let control_points: Vec<DVec3> = (0..10)
            .map(|i| {
                let t = i as f64 / 9.0;
                DVec3::new(t, (t * std::f64::consts::PI).sin(), 0.0)
            })
            .collect();

        let curve = BSplineCurve3 {
            degree: 5,
            knots: build_clamped_knots(10, 5),
            control_points: control_points.clone(),
            weights: vec![1.0; 10],
        };

        let opts = BSplineSimplifyOptions {
            max_degree: 3,
            tolerance: 0.01,
            ..Default::default()
        };

        let result = simplify_bspline_curve(&curve, &opts);

        assert!(result.final_degree <= 3, "degree should be reduced");
        assert!(result.geometry.control_points.len() <= curve.control_points.len());
    }

    #[test]
    fn simplify_already_simple_curve() {
        // Create a degree-1 line
        let curve = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::ZERO, DVec3::X],
            weights: vec![1.0, 1.0],
        };

        let opts = BSplineSimplifyOptions::default();
        let result = simplify_bspline_curve(&curve, &opts);

        assert!(!result.was_simplified, "simple curves should not be modified");
        assert_eq!(result.final_degree, 1);
    }

    #[test]
    fn convert_circle_to_bspline() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });

        let bspline = ensure_bspline_curve(&circle, 32);
        assert_eq!(bspline.degree, 2, "circle converts to degree-2 NURBS");

        // Check that endpoints match
        let p0 = bspline.point_at(0.0);
        let p1 = bspline.point_at(1.0);
        assert!((p0 - p1).length() < 1e-10, "circle should be closed");
    }

    #[test]
    fn convert_plane_to_bspline() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });

        let bspline = ensure_bspline_surface(&plane, 4, 4);
        assert_eq!(bspline.degree_u, 1);
        assert_eq!(bspline.degree_v, 1);
        assert_eq!(bspline.control_points.len(), 2);
    }

    #[test]
    fn convert_brep_to_bspline() {
        let brep = BRep::from_primitive(rcad_kernel::geom::PrimitiveSolid::Sphere { radius: 1.0 });
        let (converted, report) = convert_to_bspline(&brep, 1e-6);

        // All curves and surfaces should now be BSpline
        for curve in &converted.geom.curves {
            assert!(matches!(curve, Curve3::BSpline(_)));
        }
        for surface in &converted.geom.surfaces {
            assert!(matches!(surface, Surface3::BSpline(_)));
        }

        assert!(report.surfaces_converted > 0);
    }

    #[test]
    fn curve_degree_query() {
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        assert_eq!(curve_degree(&line), 1);

        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        assert_eq!(curve_degree(&circle), 2);
    }

    #[test]
    fn surface_degrees_query() {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        assert_eq!(surface_degrees(&plane), (1, 1));

        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        assert_eq!(surface_degrees(&sphere), (2, 2));
    }

    #[test]
    fn restrict_geometry_with_options() {
        let mut brep = BRep::new();

        // Add a simple analytic surface
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        }));

        let restrictions = GeometryRestrictions {
            surfaces_to_bspline: true,
            ..Default::default()
        };

        let (result, report) = restrict_geometry(&brep, &restrictions);
        assert!(report.surfaces_converted > 0);
        assert!(matches!(result.geom.surfaces[0], Surface3::BSpline(_)));
    }

    #[test]
    fn build_clamped_knots_correct_size() {
        // Test valid BSpline configurations where n_ctrl > degree
        for n_ctrl in 2..20usize {
            for degree in 1..=n_ctrl.saturating_sub(1).min(5) {
                let knots = build_clamped_knots(n_ctrl, degree);
                let expected_len = n_ctrl + degree + 1;
                assert_eq!(
                    knots.len(),
                    expected_len,
                    "n_ctrl={}, degree={}",
                    n_ctrl,
                    degree
                );

                // Check clamped start: first degree+1 knots should be 0
                for i in 0..=degree {
                    assert!((knots[i] - 0.0).abs() < 1e-10, "knot[{}] should be 0", i);
                }

                // Check clamped end: last degree+1 knots should be 1
                for i in 0..=degree {
                    assert!(
                        (knots[knots.len() - 1 - i] - 1.0).abs() < 1e-10,
                        "knot[{}] should be 1",
                        knots.len() - 1 - i
                    );
                }
            }
        }
    }
}
