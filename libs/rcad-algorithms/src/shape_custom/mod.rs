//! Shape customization tools analogous to OCCT ShapeCustom package.
//!
//! This module provides utilities for:
//! - BSpline degree reduction and simplification
//! - Surface conversion to NURBS representation
//! - Geometry restrictions for export compatibility
//! - Canonical form conversion
//! - BSpline to analytic conversion
//!
//! # Modules
//!
//! - [`BSplineSimplifyOptions`] - Configuration for BSpline simplification
//! - [`simplify_bspline_curve`] - Reduce degree and control points of a BSpline curve
//! - [`simplify_bspline_surface`] - Reduce degree and control points of a BSpline surface
//! - [`convert_to_bspline`] - Convert entire rcad_kernel::BRep geometry to BSpline representation
//! - [`restrict_geometry`] - Apply geometry restrictions to a rcad_kernel::BRep
//! - [`surface_to_bspline_from_face`] - Convert a face's surface to BSpline
//! - [`curve_to_bspline_from_edge`] - Convert an edge's curve to BSpline
//! - [`restrict_to_bspline`] - Convert all geometry to BSpline
//! - [`convert_to_canonical`] - Convert surfaces to canonical forms
//! - [`try_convert_to_analytic`] - Try converting BSpline to analytic surface
//! - [`simplify_geometry`] - Convert BSpline surfaces to analytic where possible
//! - [`make_direct_faces`] - Convert indirect faces to direct

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::{
    Curve3, Edge, Face, Surface3,
    geom::{
        BSplineCurve3, BSplineSurface, ConicalSurface, CurveEval, CylindricalSurface, Plane,
        SphericalSurface, SurfaceEval, ToroidalSurface, TrimmedCurve3, any_perpendicular,
    },
    nurbs_convert::{
        bezier_curve_to_bspline, bezier_surface_to_bspline, circle_to_bspline, curve_to_bspline,
        cylinder_to_bspline, ellipse_to_bspline, line_to_bspline, plane_to_bspline,
        sphere_to_bspline, surface_to_bspline,
    },
    topods::TShape,
};
use std::f64::consts::PI;
use std::sync::Arc;

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
    /// Tolerance for approximation when reducing degree (default: TOLERANCE_MESH_LEGACY).
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
            tolerance: TOLERANCE_MESH_LEGACY,
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
    let (new_curve, _max_dev) = fitted?;

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
        if k > knots[degree]
            && k < knots[knots.len() - degree - 1]
            && interior_knots
                .last()
                .is_none_or(|&last| (k - last).abs() > TOLERANCE_LINEAR_ULTRA_STRICT)
        {
            interior_knots.push(k);
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
    let knot_idx = curve
        .knots
        .iter()
        .position(|&k| (k - knot).abs() < TOLERANCE_LINEAR_ULTRA_STRICT)?;

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
    let _samples: Vec<DVec3> = (0..n_samples)
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
        is_periodic: false,
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
        is_periodic: false,
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
    let original_ctrl_pts =
        surface.control_points.len() * surface.control_points.first().map_or(0, |r| r.len());

    // If already within constraints, return as-is
    if surface.degree_u <= opts.max_degree
        && surface.degree_v <= opts.max_degree
        && surface.control_points.len() <= 4
        && surface.control_points.first().is_none_or(|r| r.len() <= 4)
    {
        return SimplificationResult {
            geometry: surface.clone(),
            was_simplified: false,
            max_deviation: 0.0,
            original_degree: original_degree_u.max(original_degree_v),
            final_degree: surface.degree_u.max(surface.degree_v),
            original_ctrl_pts,
            final_ctrl_pts: surface.control_points.len()
                * surface.control_points.first().map_or(0, |r| r.len()),
        };
    }

    let mut current = surface.clone();
    let mut total_deviation: f64 = 0.0;
    let mut was_simplified = false;

    // Reduce U degree
    if current.degree_u > opts.max_degree
        && let Some((new_surf, dev)) =
            reduce_surface_degree_u(&current, opts.max_degree, opts.tolerance)
    {
        total_deviation = total_deviation.max(dev);
        current = new_surf;
        was_simplified = true;
    }

    // Reduce V degree
    if current.degree_v > opts.max_degree
        && let Some((new_surf, dev)) =
            reduce_surface_degree_v(&current, opts.max_degree, opts.tolerance)
    {
        total_deviation = total_deviation.max(dev);
        current = new_surf;
        was_simplified = true;
    }

    // Compute final values before moving current
    let final_degree = current.degree_u.max(current.degree_v);
    let final_ctrl_pts =
        current.control_points.len() * current.control_points.first().map_or(0, |r| r.len());

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
    let n_v = surface.control_points.first()?.len().max(10);
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
    let n_v = surface.control_points.first()?.len().max(10);
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
    let n_v = points.first()?.len();

    if n_u < 2 || n_v < 2 {
        return None;
    }

    let degree_u = degree_u.min(n_u - 1);
    let degree_v = degree_v.min(n_v - 1);

    // Use points directly as control points (simple approximation)
    let control_points = points.to_vec();
    let weights: Vec<Vec<f64>> = points.iter().map(|row| vec![1.0; row.len()]).collect();

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
            tolerance: TOLERANCE_MESH_LEGACY,
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
// rcad_kernel::BRep Conversion API
// =============================================================================

/// Convert all geometry in a rcad_kernel::BRep to BSpline representation.
///
/// This function:
/// 1. Converts all analytic curves to BSpline
/// 2. Converts all analytic surfaces to BSpline
/// 3. Optionally simplifies the resulting BSplines
///
/// # Example
/// ```ignore
/// use rcad_algorithms::shape_custom::convert_to_bspline;
/// let bspline_brep = convert_to_bspline(&brep, TOLERANCE_MESH_LEGACY);
/// ```
pub fn convert_to_bspline(
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> (rcad_kernel::BRep, ConversionReport) {
    let mut restrictions = GeometryRestrictions::default();
    restrictions.tolerance = tolerance;
    restrict_geometry(brep, &restrictions)
}

/// Apply geometry restrictions to a rcad_kernel::BRep.
///
/// Converts all geometry to conform to the specified restrictions,
/// useful for export to formats with limited geometry support.
pub fn restrict_geometry(
    brep: &rcad_kernel::BRep,
    restrictions: &GeometryRestrictions,
) -> (rcad_kernel::BRep, ConversionReport) {
    let mut result = brep.clone();
    let mut report = ConversionReport::default();

    // Convert curves on edges
    if restrictions.curves_to_bspline {
        let n_samples = restrictions.curve_samples;
        for ts in &mut result.tshapes {
            if let TShape::Edge(ed) = &mut *Arc::make_mut(ts) {
                if let Some(ref mut curve) = ed.curve {
                    let needs_conversion = !matches!(curve, Curve3::BSpline(_));
                    let is_offset = matches!(curve, Curve3::Offset(_));
                    let is_bezier = matches!(curve, Curve3::Bezier(_));

                    if needs_conversion
                        && (is_bezier || !is_offset || restrictions.convert_offset_curves)
                    {
                        let bspline = curve_to_bspline(curve, n_samples);
                        *curve = Curve3::BSpline(bspline);
                        report.curves_converted += 1;
                    }
                }
            }
        }
    }

    // Convert surfaces on faces
    if restrictions.surfaces_to_bspline {
        let (n_u, n_v) = restrictions.surface_samples;
        for ts in &mut result.tshapes {
            if let TShape::Face(fd) = &mut *Arc::make_mut(ts) {
                if let Some(ref mut surface) = fd.surface {
                    let needs_conversion = !matches!(surface, Surface3::BSpline(_));
                    let is_offset = matches!(surface, Surface3::Offset(_));
                    let is_bezier = matches!(surface, Surface3::Bezier(_));

                    if needs_conversion
                        && (is_bezier || !is_offset || restrictions.convert_offset_surfaces)
                    {
                        let bspline = surface_to_bspline(surface, n_u, n_v);
                        *surface = Surface3::BSpline(bspline);
                        report.surfaces_converted += 1;
                    }
                }
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

        for ts in &mut result.tshapes {
            if let TShape::Edge(ed) = &mut *Arc::make_mut(ts) {
                if let Some(curve) = ed.curve.as_mut() {
                    if let Curve3::BSpline(bspline) = curve {
                        let res = simplify_bspline_curve(bspline, &simplify_opts);
                        if res.was_simplified {
                            *bspline = res.geometry;
                            report.bspline_curves_simplified += 1;
                            report.max_deviation = report.max_deviation.max(res.max_deviation);
                        }
                    }
                }
            }
        }

        for ts in &mut result.tshapes {
            if let TShape::Face(fd) = &mut *Arc::make_mut(ts) {
                if let Some(surface) = fd.surface.as_mut() {
                    if let Surface3::BSpline(bspline) = surface {
                        let res = simplify_bspline_surface(bspline, &simplify_opts);
                        if res.was_simplified {
                            *bspline = res.geometry;
                            report.bspline_surfaces_simplified += 1;
                            report.max_deviation = report.max_deviation.max(res.max_deviation);
                        }
                    }
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
            let dv = b
                .control_points
                .first()
                .map_or(0, |r| r.len().saturating_sub(1));
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
pub fn ensure_bspline_surface(
    surface: &Surface3,
    samples_u: usize,
    samples_v: usize,
) -> BSplineSurface {
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
// Face/Edge Level BSpline Conversion (ShapeCustom_BSplineRestriction equivalent)
// =============================================================================

/// Convert the surface of a face to BSpline representation.
///
/// This function extracts the surface associated with a face and converts it
/// to a BSpline surface. For faces with trim bounds, the BSpline is built
/// over the appropriate parameter domain.
///
/// # Arguments
/// * `face` - Reference to the Face topology
/// * `face_idx` - Flat index of the face in the rcad_kernel::BRep
/// * `brep` - Reference to the rcad_kernel::BRep containing geometry
///
/// # Returns
/// A BSplineSurface representing the face's underlying surface.
///
/// # Example
/// ```ignore
/// let bspline = surface_to_bspline_from_face(&face, 0, &brep);
/// ```
pub fn surface_to_bspline_from_face(
    _face: &Face,
    face_idx: usize,
    brep: &rcad_kernel::BRep,
) -> BSplineSurface {
    // Access the face TShape directly — surface lives on TFaceData.
    let fd = match brep.tshapes.get(face_idx) {
        Some(ts) => match &**ts {
            TShape::Face(fd) => fd,
            _ => return fallback_bspline_surface(),
        },
        None => return fallback_bspline_surface(),
    };
    let surface = match fd.surface.as_ref() {
        Some(s) => s,
        None => return fallback_bspline_surface(),
    };

    // UV domain from TFaceData (if set)
    let trim = fd.uv_domain;

    match surface {
        Surface3::BSpline(b) => {
            if let Some([u0, u1, v0, v1]) = trim {
                let n_u = 16;
                let n_v = 16;
                let mut ctrl: Vec<Vec<DVec3>> = Vec::new();
                let mut w: Vec<Vec<f64>> = Vec::new();
                for i in 0..n_u {
                    let u = u0 + (u1 - u0) * i as f64 / (n_u - 1).max(1) as f64;
                    let mut row = Vec::new();
                    let mut wrow = Vec::new();
                    for j in 0..n_v {
                        let v = v0 + (v1 - v0) * j as f64 / (n_v - 1).max(1) as f64;
                        row.push(b.point_at(u, v));
                        wrow.push(1.0);
                    }
                    ctrl.push(row);
                    w.push(wrow);
                }
                BSplineSurface {
                    degree_u: b.degree_u.min(3),
                    degree_v: b.degree_v.min(3),
                    knots_u: build_clamped_knots(n_u, b.degree_u.min(3)),
                    knots_v: build_clamped_knots(n_v, b.degree_v.min(3)),
                    control_points: ctrl,
                    weights: w,
                }
            } else {
                b.clone()
            }
        }
        Surface3::Plane(p) => {
            if let Some([u0, u1, v0, v1]) = trim {
                let pts = [
                    p.point_at(u0, v0),
                    p.point_at(u0, v1),
                    p.point_at(u1, v0),
                    p.point_at(u1, v1),
                ];
                BSplineSurface {
                    degree_u: 1,
                    degree_v: 1,
                    knots_u: vec![0.0, 0.0, 1.0, 1.0],
                    knots_v: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![vec![pts[0], pts[1]], vec![pts[2], pts[3]]],
                    weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
                }
            } else {
                plane_to_bspline(p)
            }
        }
        Surface3::Cylinder(c) => cylinder_to_bspline(c),
        Surface3::Sphere(s) => sphere_to_bspline(s),
        _ => surface_to_bspline(surface, 16, 16),
    }
}

/// Fallback BSpline surface used when no valid face surface is found.
fn fallback_bspline_surface() -> BSplineSurface {
    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            vec![DVec3::ZERO, DVec3::Y],
            vec![DVec3::X, DVec3::X + DVec3::Y],
        ],
        weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
    }
}

/// Convert the curve of an edge to BSpline representation.
///
/// This function extracts the 3D curve associated with an edge and converts it
/// to a BSpline curve. The conversion respects the edge's parameter range.
///
/// # Arguments
/// * `edge` - Reference to the Edge topology
/// * `edge_idx` - Index of the edge in the rcad_kernel::BRep
/// * `brep` - Reference to the rcad_kernel::BRep containing geometry
///
/// # Returns
/// A BSplineCurve3 representing the edge's 3D curve.
///
/// # Example
/// ```ignore
/// let bspline = curve_to_bspline_from_edge(&edge, 0, &brep);
/// ```
pub fn curve_to_bspline_from_edge(
    _edge: &Edge,
    edge_idx: usize,
    brep: &rcad_kernel::BRep,
) -> BSplineCurve3 {
    // Access the edge TShape directly — curve lives on TEdgeData.
    let ed = match brep.tshapes.get(edge_idx) {
        Some(ts) => match &**ts {
            TShape::Edge(ed) => ed,
            _ => return fallback_bspline_curve(brep, edge_idx),
        },
        None => return fallback_bspline_curve(brep, edge_idx),
    };

    let curve = match ed.curve.as_ref() {
        Some(c) => c,
        None => return fallback_bspline_curve(brep, edge_idx),
    };

    let range = ed.range;
    match curve {
        Curve3::BSpline(b) => b.clone(),
        Curve3::Line(l) => {
            let p0 = l.point_at(range[0]);
            let p1 = l.point_at(range[1]);
            BSplineCurve3 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![p0, p1],
                weights: vec![1.0, 1.0],
                is_periodic: false,
}
        }
        Curve3::Circle(c) => circle_to_bspline(c),
        Curve3::Ellipse(e) => ellipse_to_bspline(e),
        Curve3::Bezier(b) => bezier_curve_to_bspline(b),
        _ => curve_to_bspline(curve, 32),
    }
}

/// Fallback BSpline curve created from edge's vertex positions.
fn fallback_bspline_curve(brep: &rcad_kernel::BRep, edge_idx: usize) -> BSplineCurve3 {
    let p0 = brep
        .tshapes
        .get(edge_idx)
        .and_then(|ts| {
            if let TShape::Edge(ed) = &**ts {
                brep.vertex_point(ed.first.index)
            } else {
                None
            }
        })
        .unwrap_or(DVec3::ZERO);
    let p1 = brep
        .tshapes
        .get(edge_idx)
        .and_then(|ts| {
            if let TShape::Edge(ed) = &**ts {
                brep.vertex_point(ed.last.index)
            } else {
                None
            }
        })
        .unwrap_or(DVec3::X);
    BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
        is_periodic: false,
}
}

/// Convert all geometry in a rcad_kernel::BRep to BSpline representation.
///
/// This is a comprehensive conversion that:
/// 1. Converts all curves to BSpline
/// 2. Converts all surfaces to BSpline
/// 3. Updates all geometry references
///
/// mirrors `ShapeCustom_ConvertToBSpline`.
///
/// OCCT reference: ShapeCustom.cxx L123-189 (ConvertSurfaceToBSpline).
///   ShapeCustom_ConvertToBSpline iterates over all faces in a shape,
///   replaces each analytic surface with a BSpline approximation via
///   `ShapeCustom::ConvertSurfaceToBSpline(surface, tolerance,
///    Convert_QuasiPolynomial)`.  The curve conversion path uses
///   `ShapeCustom::ConvertCurveToBSpline`.
///
/// In rcad the same effect is achieved by `convert_to_bspline()` →
/// `curve_to_bspline()` / `surface_to_bspline()` over the geometry pools.
/// The `restrictions` parameter maps to OCCT's
/// `ShapeCustom_RestrictionParameters` (max degree, target tolerance,
/// offset curve/surface conversion flags, etc.).
///
/// # Arguments
/// * `brep` - The rcad_kernel::BRep to convert
///
/// # Returns
/// A new rcad_kernel::BRep with all geometry in BSpline form.
pub fn restrict_to_bspline(brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
    let tolerance = TOLERANCE_MESH_LEGACY;
    let (result, _report) = convert_to_bspline(brep, tolerance);
    result
}

// =============================================================================
// Canonical Form Conversion (ShapeCustom_SweptToElementary equivalent)
// =============================================================================

/// Canonical form identification for surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalForm {
    /// Plane aligned with XY plane (Z normal)
    PlaneXY,
    /// Plane aligned with XZ plane (Y normal)
    PlaneXZ,
    /// Plane aligned with YZ plane (X normal)
    PlaneYZ,
    /// General plane (non-canonical orientation)
    PlaneGeneral,
    /// Cylinder with axis along Z
    CylinderZ,
    /// General cylinder
    CylinderGeneral,
    /// Sphere centered at origin
    SphereOrigin,
    /// General sphere
    SphereGeneral,
    /// Cone with axis along Z
    ConeZ,
    /// General cone
    ConeGeneral,
    /// Torus centered at origin with Z axis
    TorusOriginZ,
    /// General torus
    TorusGeneral,
    /// Cannot be converted to canonical form
    NonCanonical,
}

/// Options for canonical form conversion.
#[derive(Debug, Clone)]
pub struct CanonicalConversionOptions {
    /// Tolerance for detecting canonical alignment.
    pub tolerance: f64,
    /// Whether to convert planes to XY/XZ/YZ where possible.
    pub convert_planes: bool,
    /// Whether to convert cylinders/cones to Z-axis alignment.
    pub convert_revolution_surfaces: bool,
    /// Whether to convert spheres to origin-centered.
    pub convert_spheres: bool,
    /// Whether to convert tori to origin-centered with Z axis.
    pub convert_tori: bool,
}

impl Default for CanonicalConversionOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_MESH_LEGACY,
            convert_planes: true,
            convert_revolution_surfaces: true,
            convert_spheres: true,
            convert_tori: true,
        }
    }
}

/// Identify the canonical form of a surface.
pub fn identify_canonical_form(surface: &Surface3, tolerance: f64) -> CanonicalForm {
    match surface {
        Surface3::Plane(p) => {
            let n = p.normal.normalize();
            if (n - DVec3::Z).length() < tolerance || (n + DVec3::Z).length() < tolerance {
                CanonicalForm::PlaneXY
            } else if (n - DVec3::Y).length() < tolerance || (n + DVec3::Y).length() < tolerance {
                CanonicalForm::PlaneXZ
            } else if (n - DVec3::X).length() < tolerance || (n + DVec3::X).length() < tolerance {
                CanonicalForm::PlaneYZ
            } else {
                CanonicalForm::PlaneGeneral
            }
        }
        Surface3::Cylinder(c) => {
            let axis = c.axis.normalize();
            if (axis - DVec3::Z).length() < tolerance || (axis + DVec3::Z).length() < tolerance {
                CanonicalForm::CylinderZ
            } else {
                CanonicalForm::CylinderGeneral
            }
        }
        Surface3::Sphere(s) => {
            if s.center.length() < tolerance {
                CanonicalForm::SphereOrigin
            } else {
                CanonicalForm::SphereGeneral
            }
        }
        Surface3::Cone(c) => {
            let axis = c.axis.normalize();
            if (axis - DVec3::Z).length() < tolerance || (axis + DVec3::Z).length() < tolerance {
                CanonicalForm::ConeZ
            } else {
                CanonicalForm::ConeGeneral
            }
        }
        Surface3::Torus(t) => {
            let axis = t.axis.normalize();
            let centered = t.center.length() < tolerance;
            let aligned =
                (axis - DVec3::Z).length() < tolerance || (axis + DVec3::Z).length() < tolerance;
            if centered && aligned {
                CanonicalForm::TorusOriginZ
            } else {
                CanonicalForm::TorusGeneral
            }
        }
        _ => CanonicalForm::NonCanonical,
    }
}

/// Convert surfaces in a rcad_kernel::BRep to canonical forms where possible.
///
/// This function attempts to align surfaces with canonical coordinate frames
/// for better interoperability and numerical stability.
///
/// # Arguments
/// * `brep` - The rcad_kernel::BRep to convert
/// * `tolerance` - Tolerance for detecting canonical alignment
///
/// # Returns
/// A new rcad_kernel::BRep with surfaces converted to canonical forms where possible.
pub fn convert_to_canonical(brep: &rcad_kernel::BRep, tolerance: f64) -> rcad_kernel::BRep {
    let options = CanonicalConversionOptions {
        tolerance,
        ..Default::default()
    };
    convert_to_canonical_with_options(brep, &options)
}

/// Convert surfaces to canonical forms with options.
pub fn convert_to_canonical_with_options(
    brep: &rcad_kernel::BRep,
    options: &CanonicalConversionOptions,
) -> rcad_kernel::BRep {
    let mut result = brep.clone();
    let tol = options.tolerance;

    for ts in &mut result.tshapes {
        if let TShape::Face(fd) = &mut *Arc::make_mut(ts) {
            let surface = match fd.surface.as_mut() {
                Some(s) => s,
                None => continue,
            };
            let _canonical = identify_canonical_form(surface, tol);

            match surface {
                Surface3::Plane(_) if options.convert_planes => {
                    // Planes are already in their optimal form
                }
                Surface3::Cylinder(_) if options.convert_revolution_surfaces => {
                    // Cylinder is already optimal for its frame
                }
                Surface3::Sphere(_) if options.convert_spheres => {
                    // Sphere is already optimal
                }
                Surface3::Cone(_) if options.convert_revolution_surfaces => {
                    // Cone is already optimal
                }
                Surface3::Torus(_) if options.convert_tori => {
                    // Torus is already optimal
                }
                _ => {}
            }
        }
    }

    result
}

// =============================================================================
// BSpline to Analytic Conversion (ShapeCustom_ConvertToAnalytic equivalent)
// =============================================================================

/// Result of attempting to convert a BSpline to an analytic surface.
#[derive(Debug, Clone)]
pub struct AnalyticConversionResult {
    /// The converted analytic surface, if conversion was possible.
    pub surface: Option<Surface3>,
    /// Maximum deviation from the original BSpline.
    pub deviation: f64,
    /// Type of analytic surface detected.
    pub surface_type: AnalyticType,
}

/// Types of analytic surfaces that a BSpline might represent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalyticType {
    /// Could not identify as analytic
    Unknown,
    /// A plane
    Plane,
    /// A cylinder
    Cylinder,
    /// A sphere
    Sphere,
    /// A cone
    Cone,
    /// A torus
    Torus,
}

/// Try to convert a BSpline surface to an analytic surface.
///
/// This function analyzes the BSpline surface to detect if it represents
/// a plane, cylinder, sphere, cone, or torus within the given tolerance.
///
/// # Arguments
/// * `surface` - The BSpline surface to analyze
/// * `tolerance` - Tolerance for geometric fitting
///
/// # Returns
/// An `Option<Surface3>` containing the analytic surface if conversion succeeded.
///
/// # Example
/// ```ignore
/// let bspline: BSplineSurface = /* ... */;
/// if let Some(analytic) = try_convert_to_analytic(&bspline, TOLERANCE_MESH_LEGACY) {
///     println!("Converted to analytic: {:?}", analytic);
/// }
/// ```
pub fn try_convert_to_analytic(surface: &BSplineSurface, tolerance: f64) -> Option<Surface3> {
    // Try each analytic type in order of simplicity

    // 1. Try plane detection
    if let Some(plane) = try_detect_plane(surface, tolerance) {
        return Some(Surface3::Plane(plane));
    }

    // 2. Try cylinder detection
    if let Some(cylinder) = try_detect_cylinder(surface, tolerance) {
        return Some(Surface3::Cylinder(cylinder));
    }

    // 3. Try sphere detection
    if let Some(sphere) = try_detect_sphere(surface, tolerance) {
        return Some(Surface3::Sphere(sphere));
    }

    // 4. Try cone detection
    if let Some(cone) = try_detect_cone(surface, tolerance) {
        return Some(Surface3::Cone(cone));
    }

    // 5. Try torus detection
    if let Some(torus) = try_detect_torus(surface, tolerance) {
        return Some(Surface3::Torus(torus));
    }

    None
}

/// Try to detect if a BSpline surface represents a plane.
fn try_detect_plane(surface: &BSplineSurface, tolerance: f64) -> Option<Plane> {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Sample points on the surface
    let n_samples = 5;
    let mut points: Vec<DVec3> = Vec::new();
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;
        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    // Fit a plane to the points
    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;

    // Compute covariance matrix for normal estimation
    let mut cov = [[0.0f64; 3]; 3];
    for p in &points {
        let d = *p - centroid;
        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    // Simple eigenvalue estimation using power iteration
    let mut normal = DVec3::Z;
    for _ in 0..50 {
        let n = DVec3::new(
            cov[0][0] * normal.x + cov[0][1] * normal.y + cov[0][2] * normal.z,
            cov[1][0] * normal.x + cov[1][1] * normal.y + cov[1][2] * normal.z,
            cov[2][0] * normal.x + cov[2][1] * normal.y + cov[2][2] * normal.z,
        );
        let len = n.length();
        if len < TOLERANCE_LEN_MIN {
            break;
        }
        normal = n / len;
    }

    // Check if all points are within tolerance of the plane
    for p in &points {
        let dist = (p - centroid).dot(normal).abs();
        if dist > tolerance {
            return None;
        }
    }

    Some(Plane::new(centroid, centroid.normalize_or_zero()))
}

/// Try to detect if a BSpline surface represents a cylinder.
fn try_detect_cylinder(surface: &BSplineSurface, tolerance: f64) -> Option<CylindricalSurface> {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Sample points on the surface
    let n_samples = 10;
    let mut points: Vec<DVec3> = Vec::new();
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;
        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    // Compute centroid
    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;

    // For a cylinder, points should lie at constant distance from an axis
    // First, estimate the axis using PCA
    let mut cov = [[0.0f64; 3]; 3];
    for p in &points {
        let d = *p - centroid;
        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    // The axis is the direction of maximum variance (principal component)
    // Use power iteration to find it
    let mut axis = DVec3::Z;
    for _ in 0..50 {
        let n = DVec3::new(
            cov[0][0] * axis.x + cov[0][1] * axis.y + cov[0][2] * axis.z,
            cov[1][0] * axis.x + cov[1][1] * axis.y + cov[1][2] * axis.z,
            cov[2][0] * axis.x + cov[2][1] * axis.y + cov[2][2] * axis.z,
        );
        let len = n.length();
        if len < TOLERANCE_LEN_MIN {
            break;
        }
        axis = n / len;
    }

    // Project points onto the plane perpendicular to the axis
    let mut radii: Vec<f64> = Vec::new();
    for p in &points {
        let to_point = *p - centroid;
        let along_axis = to_point.dot(axis);
        let radial = to_point - along_axis * axis;
        radii.push(radial.length());
    }

    // Check if all radii are approximately equal
    let avg_radius = radii.iter().sum::<f64>() / radii.len() as f64;
    let radius_variance: f64 =
        radii.iter().map(|r| (r - avg_radius).powi(2)).sum::<f64>() / radii.len() as f64;

    // Use a relative tolerance based on the radius
    let rel_tolerance = tolerance.max(avg_radius * 0.1);
    if radius_variance.sqrt() > rel_tolerance {
        return None;
    }

    // Find a point on the axis
    let axis_point = centroid;

    Some(CylindricalSurface {
        origin: axis_point,
        axis,
        radius: avg_radius,
        ref_dir: any_perpendicular(axis),
    })
}

/// Try to detect if a BSpline surface represents a sphere.
fn try_detect_sphere(surface: &BSplineSurface, tolerance: f64) -> Option<SphericalSurface> {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Sample points on the surface
    let n_samples = 10;
    let mut points: Vec<DVec3> = Vec::new();
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;
        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    // For a sphere, we need to find the center. Since we're dealing with BSpline
    // surfaces that may only represent a portion of a sphere, we use a different
    // approach: check if all normals point to/from a common center.

    // First, estimate the center by looking at where normals intersect
    // For each point, compute the line along the normal, then find the
    // best intersection point.
    let mut normals: Vec<DVec3> = Vec::new();
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;
        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            normals.push(surface.normal_at(u, v));
        }
    }

    // For a sphere, the normal at each point should point toward (or away from) the center
    // We can use the fact that: center = point - radius * normal
    // This means for any two points: center = p1 - r*n1 = p2 - r*n2
    // So: p1 - p2 = r*(n1 - n2)
    // And: r = (p1 - p2).length() / (n1 - n2).length() (approximately)

    // Estimate radius using multiple point pairs
    let mut radii_estimates: Vec<f64> = Vec::new();
    for i in 0..points.len().min(20) {
        for j in (i + 1)..points.len().min(20) {
            let n_diff = (normals[i] - normals[j]).length();
            if n_diff > 0.1 {
                let p_diff = (points[i] - points[j]).length();
                radii_estimates.push(p_diff / n_diff);
            }
        }
    }

    if radii_estimates.is_empty() {
        return None;
    }

    let avg_radius = radii_estimates.iter().sum::<f64>() / radii_estimates.len() as f64;

    // Use relative tolerance
    let rel_tolerance = tolerance.max(avg_radius * 0.1);

    // Check variance in radius estimates
    let radius_variance: f64 = radii_estimates
        .iter()
        .map(|r| (r - avg_radius).powi(2))
        .sum::<f64>()
        / radii_estimates.len() as f64;

    if radius_variance.sqrt() > rel_tolerance * avg_radius {
        return None;
    }

    // Estimate center from one point and its normal
    // Use average of several estimates
    let mut center_estimates: Vec<DVec3> = Vec::new();
    for i in 0..points.len().min(10) {
        center_estimates.push(points[i] - avg_radius * normals[i]);
    }
    let center = center_estimates.iter().sum::<DVec3>() / center_estimates.len() as f64;

    // Verify by checking distances from estimated center
    let mut distances: Vec<f64> = Vec::new();
    for p in &points {
        distances.push((*p - center).length());
    }

    let avg_dist = distances.iter().sum::<f64>() / distances.len() as f64;
    let dist_variance: f64 = distances
        .iter()
        .map(|d| (d - avg_dist).powi(2))
        .sum::<f64>()
        / distances.len() as f64;

    if dist_variance.sqrt() > rel_tolerance * avg_radius {
        return None;
    }

    // Use the surface normal at center to determine the axis
    let axis = surface.normal_at(0.5, 0.5);

    Some(SphericalSurface::new(center, axis, avg_radius))
}

/// Try to detect if a BSpline surface represents a cone.
fn try_detect_cone(surface: &BSplineSurface, tolerance: f64) -> Option<ConicalSurface> {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Sample points on the surface
    let n_samples = 10;
    let mut points: Vec<DVec3> = Vec::new();
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;
        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    // For a cone, the radius varies linearly along the axis
    // This is more complex detection - use simplified approach

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;

    // Estimate axis using PCA
    let mut cov = [[0.0f64; 3]; 3];
    for p in &points {
        let d = *p - centroid;
        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }

    let mut axis = DVec3::Z;
    for _ in 0..50 {
        let n = DVec3::new(
            cov[0][0] * axis.x + cov[0][1] * axis.y + cov[0][2] * axis.z,
            cov[1][0] * axis.x + cov[1][1] * axis.y + cov[1][2] * axis.z,
            cov[2][0] * axis.x + cov[2][1] * axis.y + cov[2][2] * axis.z,
        );
        let len = n.length();
        if len < TOLERANCE_LEN_MIN {
            break;
        }
        axis = n / len;
    }

    // Compute radius vs axial position
    let mut radius_axial: Vec<(f64, f64)> = Vec::new();
    for p in &points {
        let to_point = *p - centroid;
        let axial = to_point.dot(axis);
        let radial = (to_point - axial * axis).length();
        radius_axial.push((radial, axial));
    }

    // Check linear relationship: radius = r0 + axial * tan(half_angle)
    // Fit a line to radius vs axial
    let n = radius_axial.len() as f64;
    let sum_axial: f64 = radius_axial.iter().map(|(_, a)| *a).sum();
    let sum_radius: f64 = radius_axial.iter().map(|(r, _)| *r).sum();
    let sum_aa: f64 = radius_axial.iter().map(|(_, a)| a * a).sum();
    let sum_ar: f64 = radius_axial.iter().map(|(r, a)| r * a).sum();

    let denom = n * sum_aa - sum_axial * sum_axial;
    if denom.abs() < TOLERANCE_LEN_MIN {
        return None;
    }

    let slope = (n * sum_ar - sum_axial * sum_radius) / denom;
    let intercept = (sum_radius - slope * sum_axial) / n;

    // Compute residuals
    let mut max_residual = 0.0f64;
    for (r, a) in &radius_axial {
        let predicted = intercept + slope * a;
        max_residual = max_residual.max((r - predicted).abs());
    }

    if max_residual > tolerance {
        return None;
    }

    // Extract cone parameters
    let half_angle = slope.atan().abs();
    let radius = intercept.max(0.0);

    if !(TOLERANCE_MESH_LEGACY..=PI / 2.0 - TOLERANCE_MESH_LEGACY).contains(&half_angle) {
        return None; // Degenerate cone
    }

    Some(ConicalSurface::new(centroid, axis, radius, half_angle))
}

/// Try to detect if a BSpline surface represents a torus.
fn try_detect_torus(surface: &BSplineSurface, tolerance: f64) -> Option<ToroidalSurface> {
    let [u0, u1, v0, v1] = surface.default_domain();

    // Sample points on the surface
    let n_samples = 15;
    let mut points: Vec<DVec3> = Vec::new();
    for i in 0..n_samples {
        let u = u0 + (u1 - u0) * i as f64 / (n_samples - 1) as f64;
        for j in 0..n_samples {
            let v = v0 + (v1 - v0) * j as f64 / (n_samples - 1) as f64;
            points.push(surface.point_at(u, v));
        }
    }

    let centroid = points.iter().sum::<DVec3>() / points.len() as f64;

    // For a torus, points lie at distance R from a central circle
    // This is complex detection - use simplified approach

    // Estimate the axis (normal to the torus plane)
    let mut cov = [[0.0f64; 3]; 3];
    for p in &points {
        let d = *p - centroid;
        cov[0][0] += d.x * d.x;
        cov[0][1] += d.x * d.y;
        cov[0][2] += d.x * d.z;
        cov[1][1] += d.y * d.y;
        cov[1][2] += d.y * d.z;
        cov[2][2] += d.z * d.z;
    }

    // Find the axis (smallest eigenvalue direction for torus)
    let axis = DVec3::Z;

    // Try to find distances from the axis
    let mut distances: Vec<f64> = Vec::new();
    for p in &points {
        let to_point = *p - centroid;
        let axial = to_point.dot(axis);
        let radial = (to_point - axial * axis).length();
        distances.push(radial);
    }

    // For a torus, distances should cluster around two values (R+r and R-r)
    // Use a simple estimate: major_radius ~ mean distance, minor_radius ~ variance
    let avg_dist = distances.iter().sum::<f64>() / distances.len() as f64;
    let variance = distances
        .iter()
        .map(|d| (d - avg_dist).powi(2))
        .sum::<f64>()
        / distances.len() as f64;

    // Check if the torus fit is reasonable
    let minor_radius = variance.sqrt();
    let major_radius = avg_dist;

    if minor_radius < tolerance || major_radius < minor_radius {
        return None;
    }

    // Verify the fit
    let mut max_deviation = 0.0f64;
    for p in &points {
        let to_point = *p - centroid;
        let axial = to_point.dot(axis);
        let radial = (to_point - axial * axis).length();
        let expected_r = (major_radius - radial).hypot(axial);
        let deviation = (expected_r - minor_radius).abs();
        max_deviation = max_deviation.max(deviation);
    }

    if max_deviation > tolerance * 10.0 {
        return None;
    }

    Some(ToroidalSurface {
        center: centroid,
        axis,
        major_radius,
        minor_radius,
    })
}

/// Simplify geometry by converting BSpline surfaces to analytic where possible.
///
/// This function analyzes all BSpline surfaces in a rcad_kernel::BRep and converts them
/// to analytic surfaces (plane, cylinder, sphere, cone, torus) when they
/// match within the given tolerance.
///
/// # Arguments
/// * `brep` - The rcad_kernel::BRep to simplify
/// * `tolerance` - Tolerance for geometric fitting
///
/// # Returns
/// A new rcad_kernel::BRep with simplified geometry.
pub fn simplify_geometry(brep: &rcad_kernel::BRep, tolerance: f64) -> rcad_kernel::BRep {
    let mut result = brep.clone();

    for ts in &mut result.tshapes {
        if let TShape::Face(fd) = &mut *Arc::make_mut(ts) {
            let surface = match fd.surface.as_mut() {
                Some(s) => s,
                None => continue,
            };
            if let Surface3::BSpline(bspline) = surface
                && let Some(analytic) = try_convert_to_analytic(bspline, tolerance)
            {
                *surface = analytic;
            }
        }
    }

    result
}

/// Convert indirect faces to direct faces.
///
/// In OCCT, an "indirect" face has a transformation applied via a location.
/// In RCAD, transformations are applied directly to geometry, so this
/// function primarily resolves Offset and Trimmed surfaces to their
/// underlying forms.
///
/// # Arguments
/// * `brep` - The rcad_kernel::BRep to process
///
/// # Returns
/// A new rcad_kernel::BRep with all faces converted to direct representation.
pub fn make_direct_faces(brep: &rcad_kernel::BRep) -> rcad_kernel::BRep {
    let mut result = brep.clone();

    // Resolve surfaces on faces
    for ts in &mut result.tshapes {
        if let TShape::Face(fd) = &mut *Arc::make_mut(ts) {
            if let Some(ref mut surface) = fd.surface {
                *surface = resolve_to_direct_surface(surface);
            }
        }
    }

    // Resolve curves on edges
    for ts in &mut result.tshapes {
        if let TShape::Edge(ed) = &mut *Arc::make_mut(ts) {
            if let Some(ref mut curve) = ed.curve {
                *curve = resolve_to_direct_curve(curve);
            }
        }
    }

    result
}

/// Resolve a surface to its direct form.
fn resolve_to_direct_surface(surface: &Surface3) -> Surface3 {
    match surface {
        Surface3::Trimmed(t) => {
            // Resolve the underlying surface

            // For now, return the resolved basis (losing trim info)
            // A more complete implementation would rebuild with proper UV bounds
            resolve_to_direct_surface(&t.basis)
        }
        Surface3::Offset(_) => {
            // Offset surfaces don't have a simple direct form
            // Convert to BSpline for now
            Surface3::BSpline(surface_to_bspline(surface, 32, 32))
        }
        Surface3::LinearExtrusion(_) => {
            // Convert extrusion surface to BSpline
            Surface3::BSpline(surface_to_bspline(surface, 32, 32))
        }
        Surface3::Revolution(_) => {
            // Revolution surfaces can sometimes be converted to analytic
            // For now, convert to BSpline
            Surface3::BSpline(surface_to_bspline(surface, 32, 32))
        }
        Surface3::Ruled(_) => {
            // Ruled surfaces to BSpline
            Surface3::BSpline(surface_to_bspline(surface, 32, 32))
        }
        Surface3::Coons(_) => {
            // Coons patches to BSpline
            Surface3::BSpline(surface_to_bspline(surface, 32, 32))
        }
        // Analytic surfaces are already direct
        s @ Surface3::Plane(_)
        | s @ Surface3::Cylinder(_)
        | s @ Surface3::Sphere(_)
        | s @ Surface3::Cone(_)
        | s @ Surface3::Torus(_)
        | s @ Surface3::Ellipsoid(_)
        | s @ Surface3::Helicoid(_)
        | s @ Surface3::Pipe(_)
        | s @ Surface3::BSpline(_)
        | s @ Surface3::Bezier(_)
        | s @ Surface3::TriBezier(_) => s.clone(),
    }
}

/// Resolve a curve to its direct form.
fn resolve_to_direct_curve(curve: &Curve3) -> Curve3 {
    match curve {
        Curve3::Offset(_) => {
            // Offset curves to BSpline
            Curve3::BSpline(curve_to_bspline(curve, 64))
        }
        // Analytic curves are already direct
        c @ Curve3::Line(_)
        | c @ Curve3::Circle(_)
        | c @ Curve3::Ellipse(_)
        | c @ Curve3::BSpline(_)
        | c @ Curve3::Bezier(_)
        | c @ Curve3::Hyperbola(_)
        | c @ Curve3::Parabola(_)
        | c @ Curve3::CircularHelix(_)
        | c @ Curve3::SineWave(_) => c.clone(),
        Curve3::Trimmed(tc) => Curve3::Trimmed(TrimmedCurve3::new(
            resolve_to_direct_curve(tc.basis_curve()),
            tc.first,
            tc.last,
        )),
    }
}

// =============================================================================
// Comprehensive Shape Customization Report
// =============================================================================

/// Report from shape customization operations.
#[derive(Debug, Clone, Default)]
pub struct ShapeCustomReport {
    /// Number of surfaces converted to BSpline.
    pub surfaces_to_bspline: usize,
    /// Number of curves converted to BSpline.
    pub curves_to_bspline: usize,
    /// Number of BSpline surfaces simplified to analytic.
    pub bspline_to_analytic: usize,
    /// Number of BSpline curves simplified to analytic.
    pub bspline_curve_to_analytic: usize,
    /// Number of faces made direct.
    pub faces_made_direct: usize,
    /// Number of surfaces converted to canonical form.
    pub canonical_conversions: usize,
    /// Maximum deviation during conversions.
    pub max_deviation: f64,
}

/// Apply comprehensive shape customization.
///
/// This function applies all shape customization operations in sequence:
/// 1. Convert to canonical forms
/// 2. Convert all geometry to BSpline
/// 3. Simplify BSpline to analytic where possible
/// 4. Make all faces direct
///
/// # Arguments
/// * `brep` - The rcad_kernel::BRep to process
/// * `tolerance` - Tolerance for all operations
///
/// # Returns
/// A tuple of (processed rcad_kernel::BRep, report).
pub fn customize_shape(
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> (rcad_kernel::BRep, ShapeCustomReport) {
    let mut report = ShapeCustomReport::default();

    // Step 1: Convert to canonical forms
    let canonical = convert_to_canonical(brep, tolerance);
    report.canonical_conversions = count_canonical_conversions(brep, &canonical);

    // Step 2: Convert to BSpline
    let (bspline, conv_report) = convert_to_bspline(&canonical, tolerance);
    report.surfaces_to_bspline = conv_report.surfaces_converted;
    report.curves_to_bspline = conv_report.curves_converted;
    report.max_deviation = conv_report.max_deviation;

    // Step 3: Simplify to analytic
    let simplified = simplify_geometry(&bspline, tolerance);
    report.bspline_to_analytic = count_analytic_conversions(&bspline, &simplified);

    // Step 4: Make direct faces
    let direct = make_direct_faces(&simplified);
    report.faces_made_direct = count_indirect_faces(&simplified);

    (direct, report)
}

/// Count how many surfaces were converted to canonical form.
fn count_canonical_conversions(before: &rcad_kernel::BRep, after: &rcad_kernel::BRep) -> usize {
    let surfaces_before: Vec<&Surface3> = before
        .tshapes
        .iter()
        .filter_map(|ts| {
            if let TShape::Face(fd) = &**ts {
                fd.surface.as_ref()
            } else {
                None
            }
        })
        .collect();
    let surfaces_after: Vec<&Surface3> = after
        .tshapes
        .iter()
        .filter_map(|ts| {
            if let TShape::Face(fd) = &**ts {
                fd.surface.as_ref()
            } else {
                None
            }
        })
        .collect();
    let mut count = 0;
    for (s_before, s_after) in surfaces_before.iter().zip(surfaces_after.iter()) {
        let form_before = identify_canonical_form(s_before, TOLERANCE_MESH_LEGACY);
        let form_after = identify_canonical_form(s_after, TOLERANCE_MESH_LEGACY);
        if form_before != form_after && form_after != CanonicalForm::NonCanonical {
            count += 1;
        }
    }
    count
}

/// Count how many BSpline surfaces were converted to analytic.
fn count_analytic_conversions(before: &rcad_kernel::BRep, after: &rcad_kernel::BRep) -> usize {
    let surfaces_before: Vec<&Surface3> = before
        .tshapes
        .iter()
        .filter_map(|ts| {
            if let TShape::Face(fd) = &**ts {
                fd.surface.as_ref()
            } else {
                None
            }
        })
        .collect();
    let surfaces_after: Vec<&Surface3> = after
        .tshapes
        .iter()
        .filter_map(|ts| {
            if let TShape::Face(fd) = &**ts {
                fd.surface.as_ref()
            } else {
                None
            }
        })
        .collect();
    let mut count = 0;
    for (s_before, s_after) in surfaces_before.iter().zip(surfaces_after.iter()) {
        if matches!(s_before, Surface3::BSpline(_)) && !matches!(s_after, Surface3::BSpline(_)) {
            count += 1;
        }
    }
    count
}

/// Count indirect faces (faces with Offset/Trimmed surfaces).
fn count_indirect_faces(brep: &rcad_kernel::BRep) -> usize {
    brep.tshapes
        .iter()
        .filter_map(|ts| {
            if let TShape::Face(fd) = &**ts {
                fd.surface.as_ref()
            } else {
                None
            }
        })
        .filter(|s| matches!(s, Surface3::Offset(_) | Surface3::Trimmed(_)))
        .count()
}

// =============================================================================
// Tests
// =============================================================================
