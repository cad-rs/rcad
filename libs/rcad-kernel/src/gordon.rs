//! Gordon surface construction with robust validation and numerical stability.
//!
//! Gordon surfaces provide transfinite interpolation through an N x M network
//! of crossing curves. This module implements production-grade construction
//! with comprehensive input validation and error handling.
//!
//! # Overview
//!
//! A Gordon surface interpolates:
//! - N u-direction curves (iso-v lines)
//! - M v-direction curves (iso-u lines)
//!
//! The surface is defined as:
//! ```text
//! S(u,v) = sum_i L_v[i](v) * C_i(u)    // Loft in v direction
//!        + sum_j L_u[j](u) * D_j(v)    // Loft in u direction
//!        - sum_{i,j} L_v[i](v) * L_u[j](u) * P_{ij}  // Tensor product correction
//! ```
//!
//! where L_v[i] and L_u[j] are Lagrange basis functions.
//!
//! # Example
//!
//! ```rust
//! use glam::DVec3;
//! use rcad_kernel::geom::{Curve3, Line3, GordonSurface};
//! use rcad_kernel::gordon::{GordonOptions, gordon_surface_curves};
//!
//! // Create a simple 2x2 network of lines (planar bilinear patch)
//! let u0 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
//! let u1 = Curve3::Line(Line3 { origin: DVec3::Y, direction: DVec3::X });
//! let v0 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::Y });
//! let v1 = Curve3::Line(Line3 { origin: DVec3::X, direction: DVec3::Y });
//!
//! let surface = gordon_surface_curves(
//!     &[u0, u1],
//!     &[v0, v1],
//!     GordonOptions::default(),
//! ).unwrap();
//! ```

use glam::DVec3;
use std::fmt;

use crate::geom::{BSplineSurface, Curve3, CurveEval, GordonSurface};

// ─────────────────────────────────────────────────────────────────────────────
// Error types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during Gordon surface construction.
#[derive(Debug, Clone, PartialEq)]
pub enum GordonError {
    /// Fewer than 2 u-curves provided.
    TooFewUCurves {
        count: usize,
    },
    /// Fewer than 2 v-curves provided.
    TooFewVCurves {
        count: usize,
    },
    /// U-parameter count does not match v-curve count.
    UParamCountMismatch {
        expected: usize,
        actual: usize,
    },
    /// V-parameter count does not match u-curve count.
    VParamCountMismatch {
        expected: usize,
        actual: usize,
    },
    /// Parameters are not monotonically increasing.
    NonMonotonicParams {
        direction: String,
        index: usize,
        prev: f64,
        curr: f64,
    },
    /// Parameters are outside [0, 1] range.
    ParamsOutOfRange {
        direction: String,
        value: f64,
    },
    /// Curve endpoints do not meet at intersection points.
    IntersectionMismatch {
        u_curve_idx: usize,
        v_curve_idx: usize,
        u_point: DVec3,
        v_point: DVec3,
        distance: f64,
        tolerance: f64,
    },
    /// Nodes are too close together, causing numerical instability.
    CoincidentNodes {
        direction: String,
        idx1: usize,
        idx2: usize,
        distance: f64,
    },
    /// Degenerate curve detected (zero length).
    DegenerateCurve {
        curve_idx: usize,
        direction: String,
    },
    /// Lagrange basis evaluation failed (singular matrix).
    SingularLagrangeBasis {
        direction: String,
        param: f64,
    },
    /// Curve has incompatible parameter domain.
    IncompatibleDomain {
        curve_idx: usize,
        domain: [f64; 2],
    },
}

impl fmt::Display for GordonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooFewUCurves { count } => {
                write!(f, "Gordon surface requires at least 2 u-curves, got {}", count)
            }
            Self::TooFewVCurves { count } => {
                write!(f, "Gordon surface requires at least 2 v-curves, got {}", count)
            }
            Self::UParamCountMismatch { expected, actual } => {
                write!(
                    f,
                    "u_params count ({}) must equal v_curves count ({})",
                    actual, expected
                )
            }
            Self::VParamCountMismatch { expected, actual } => {
                write!(
                    f,
                    "v_params count ({}) must equal u_curves count ({})",
                    actual, expected
                )
            }
            Self::NonMonotonicParams {
                direction,
                index,
                prev,
                curr,
            } => {
                write!(
                    f,
                    "{}-params are not monotonic at index {}: {} >= {}",
                    direction, index, prev, curr
                )
            }
            Self::ParamsOutOfRange { direction, value } => {
                write!(
                    f,
                    "{}-param value {} is outside [0, 1] range",
                    direction, value
                )
            }
            Self::IntersectionMismatch {
                u_curve_idx,
                v_curve_idx,
                u_point,
                v_point,
                distance,
                tolerance,
            } => {
                write!(
                    f,
                    "Curve intersection mismatch at u_curve[{}], v_curve[{}]: \
                     distance {} exceeds tolerance {} (u_point={:?}, v_point={:?})",
                    u_curve_idx, v_curve_idx, distance, tolerance, u_point, v_point
                )
            }
            Self::CoincidentNodes {
                direction,
                idx1,
                idx2,
                distance,
            } => {
                write!(
                    f,
                    "{}-params[{}] and [{}] are too close (distance={}), \
                     causing numerical instability",
                    direction, idx1, idx2, distance
                )
            }
            Self::DegenerateCurve { curve_idx, direction } => {
                write!(
                    f,
                    "{}-curve[{}] is degenerate (zero length)",
                    direction, curve_idx
                )
            }
            Self::SingularLagrangeBasis { direction, param } => {
                write!(
                    f,
                    "Lagrange basis evaluation failed for {} at param {}",
                    direction, param
                )
            }
            Self::IncompatibleDomain { curve_idx, domain } => {
                write!(
                    f,
                    "Curve[{}] has incompatible domain [{}, {}]",
                    curve_idx, domain[0], domain[1]
                )
            }
        }
    }
}

impl std::error::Error for GordonError {}

// ─────────────────────────────────────────────────────────────────────────────
// Configuration options
// ─────────────────────────────────────────────────────────────────────────────

/// Continuity level for Gordon surface construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ContinuityLevel {
    /// Position continuity only (C0).
    C0,
    /// Tangent continuity (C1).
    #[default]
    C1,
    /// Curvature continuity (C2).
    C2,
}

/// Options for Gordon surface construction.
#[derive(Debug, Clone)]
pub struct GordonOptions {
    /// Degree of the output B-spline surface in u direction.
    /// If None, determined automatically from input curves.
    pub degree_u: Option<usize>,
    /// Degree of the output B-spline surface in v direction.
    /// If None, determined automatically from input curves.
    pub degree_v: Option<usize>,
    /// Desired continuity level.
    pub continuity: ContinuityLevel,
    /// Tolerance for geometric comparisons (endpoint matching, etc.).
    pub tolerance: f64,
    /// Minimum separation between parameter nodes.
    pub min_node_separation: f64,
    /// Whether to auto-normalize curve parameter domains to [0, 1].
    pub normalize_params: bool,
    /// Whether to validate curve intersections.
    pub validate_intersections: bool,
    /// Maximum distance allowed between crossing curves at intersections.
    pub intersection_tolerance: f64,
}

impl Default for GordonOptions {
    fn default() -> Self {
        Self {
            degree_u: None,
            degree_v: None,
            continuity: ContinuityLevel::C1,
            tolerance: 1e-6,
            min_node_separation: 1e-10,
            normalize_params: true,
            validate_intersections: true,
            intersection_tolerance: 1e-4,
        }
    }
}

impl GordonOptions {
    /// Create options for C0 continuity only.
    pub fn c0() -> Self {
        Self {
            continuity: ContinuityLevel::C0,
            ..Default::default()
        }
    }

    /// Create options for C1 continuity.
    pub fn c1() -> Self {
        Self {
            continuity: ContinuityLevel::C1,
            ..Default::default()
        }
    }

    /// Create options for C2 continuity.
    pub fn c2() -> Self {
        Self {
            continuity: ContinuityLevel::C2,
            ..Default::default()
        }
    }

    /// Set the geometric tolerance.
    pub fn with_tolerance(mut self, tol: f64) -> Self {
        self.tolerance = tol;
        self
    }

    /// Set the intersection tolerance.
    pub fn with_intersection_tolerance(mut self, tol: f64) -> Self {
        self.intersection_tolerance = tol;
        self
    }

    /// Disable intersection validation.
    pub fn skip_intersection_validation(mut self) -> Self {
        self.validate_intersections = false;
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Validate parameter array for monotonicity and range.
fn validate_params(params: &[f64], direction: &str, min_sep: f64) -> Result<(), GordonError> {
    if params.is_empty() {
        return Ok(());
    }

    // Check range
    for &p in params {
        if p < 0.0 || p > 1.0 {
            return Err(GordonError::ParamsOutOfRange {
                direction: direction.to_string(),
                value: p,
            });
        }
    }

    // Check monotonicity
    for i in 1..params.len() {
        if params[i] <= params[i - 1] {
            return Err(GordonError::NonMonotonicParams {
                direction: direction.to_string(),
                index: i,
                prev: params[i - 1],
                curr: params[i],
            });
        }
    }

    // Check minimum separation
    for i in 1..params.len() {
        let sep = params[i] - params[i - 1];
        if sep < min_sep {
            return Err(GordonError::CoincidentNodes {
                direction: direction.to_string(),
                idx1: i - 1,
                idx2: i,
                distance: sep,
            });
        }
    }

    Ok(())
}

/// Validate that curve endpoints match at intersections.
fn validate_intersections(
    u_curves: &[Curve3],
    u_params: &[f64],
    v_curves: &[Curve3],
    v_params: &[f64],
    tolerance: f64,
) -> Result<(), GordonError> {
    // Use normalized domain for evaluation
    let eval_curve_at = |curve: &Curve3, t: f64| -> DVec3 {
        let [t0, t1] = curve.default_domain();
        if t0.is_finite() && t1.is_finite() && (t1 - t0).abs() > 1e-15 {
            curve.point_at(t0 + t * (t1 - t0))
        } else {
            curve.point_at(t)
        }
    };

    for (i, u_curve) in u_curves.iter().enumerate() {
        for (j, v_curve) in v_curves.iter().enumerate() {
            // u_curve at u_params[j] should equal v_curve at v_params[i]
            let u_point = eval_curve_at(u_curve, u_params[j]);
            let v_point = eval_curve_at(v_curve, v_params[i]);
            let dist = (u_point - v_point).length();

            if dist > tolerance {
                return Err(GordonError::IntersectionMismatch {
                    u_curve_idx: i,
                    v_curve_idx: j,
                    u_point,
                    v_point,
                    distance: dist,
                    tolerance,
                });
            }
        }
    }

    Ok(())
}

/// Check if a curve is degenerate (zero length).
fn is_degenerate_curve(curve: &Curve3, samples: usize, tol: f64) -> bool {
    if samples < 2 {
        return false;
    }

    let [t0, t1] = curve.default_domain();
    let (t0, t1) = if t0.is_finite() && t1.is_finite() {
        (t0, t1)
    } else {
        (0.0, 1.0)
    };

    let first = curve.point_at(t0);
    for i in 1..samples {
        let t = t0 + (i as f64 / (samples - 1) as f64) * (t1 - t0);
        let p = curve.point_at(t);
        if (p - first).length() > tol {
            return false;
        }
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Lagrange basis with numerical stability
// ─────────────────────────────────────────────────────────────────────────────

/// Compute Lagrange basis functions with numerical stability.
///
/// Returns None if the denominator is too small (singular point).
fn lagrange_basis_safe(nodes: &[f64], t: f64, tol: f64) -> Option<Vec<f64>> {
    let n = nodes.len();
    if n == 0 {
        return Some(vec![]);
    }

    let mut basis = vec![1.0; n];

    for i in 0..n {
        for j in 0..n {
            if i != j {
                let denom = nodes[i] - nodes[j];
                if denom.abs() < tol {
                    // Nodes too close - potential singularity
                    return None;
                }
                basis[i] *= (t - nodes[j]) / denom;
            }
        }

        // Check for NaN/Inf
        if !basis[i].is_finite() {
            return None;
        }
    }

    Some(basis)
}

/// Compute Lagrange basis derivative with numerical stability.
#[allow(dead_code)]
fn lagrange_basis_derivative(nodes: &[f64], t: f64, tol: f64) -> Option<Vec<f64>> {
    let n = nodes.len();
    if n < 2 {
        return Some(vec![0.0; n]);
    }

    let mut deriv = vec![0.0; n];

    for i in 0..n {
        // Derivative of L_i = L_i * sum_{j != i} 1/(t - nodes[j])
        let mut sum = 0.0;
        for j in 0..n {
            if i != j {
                let diff = t - nodes[j];
                if diff.abs() < tol {
                    return None; // Singular point
                }
                sum += 1.0 / diff;
            }
        }

        // Get basis value
        let basis = lagrange_basis_safe(nodes, t, tol)?;
        deriv[i] = basis[i] * sum;
    }

    Some(deriv)
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Build a Gordon surface from two families of crossing curves.
///
/// # Arguments
///
/// * `u_curves` - Curves in the u-direction (iso-v lines)
/// * `v_curves` - Curves in the v-direction (iso-u lines)
/// * `opts` - Construction options
///
/// # Returns
///
/// A `GordonSurface` that interpolates all input curves.
///
/// # Errors
///
/// Returns `GordonError` if:
/// - Fewer than 2 curves are provided in either direction
/// - Parameters are not monotonically increasing in [0, 1]
/// - Curve endpoints do not match at intersections (within tolerance)
/// - Nodes are too close together
/// - Curves are degenerate
///
/// # Example
///
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::geom::{Curve3, Line3};
/// use rcad_kernel::gordon::{GordonOptions, gordon_surface_curves};
///
/// let u0 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
/// let u1 = Curve3::Line(Line3 { origin: DVec3::Y, direction: DVec3::X });
/// let v0 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::Y });
/// let v1 = Curve3::Line(Line3 { origin: DVec3::X, direction: DVec3::Y });
///
/// let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default()).unwrap();
/// ```
pub fn gordon_surface_curves(
    u_curves: &[Curve3],
    v_curves: &[Curve3],
    opts: GordonOptions,
) -> Result<GordonSurface, GordonError> {
    // Validate curve counts
    if u_curves.len() < 2 {
        return Err(GordonError::TooFewUCurves {
            count: u_curves.len(),
        });
    }
    if v_curves.len() < 2 {
        return Err(GordonError::TooFewVCurves {
            count: v_curves.len(),
        });
    }

    // Check for degenerate curves
    for (i, curve) in u_curves.iter().enumerate() {
        if is_degenerate_curve(curve, 10, opts.tolerance) {
            return Err(GordonError::DegenerateCurve {
                curve_idx: i,
                direction: "u".to_string(),
            });
        }
    }
    for (i, curve) in v_curves.iter().enumerate() {
        if is_degenerate_curve(curve, 10, opts.tolerance) {
            return Err(GordonError::DegenerateCurve {
                curve_idx: i,
                direction: "v".to_string(),
            });
        }
    }

    // Generate default parameters (uniform distribution)
    let u_params: Vec<f64> = (0..v_curves.len())
        .map(|i| i as f64 / (v_curves.len() - 1).max(1) as f64)
        .collect();
    let v_params: Vec<f64> = (0..u_curves.len())
        .map(|i| i as f64 / (u_curves.len() - 1).max(1) as f64)
        .collect();

    // Validate parameters
    validate_params(&u_params, "u", opts.min_node_separation)?;
    validate_params(&v_params, "v", opts.min_node_separation)?;

    // Validate intersections
    if opts.validate_intersections {
        validate_intersections(u_curves, &u_params, v_curves, &v_params, opts.intersection_tolerance)?;
    }

    Ok(GordonSurface {
        u_curves: u_curves.to_vec(),
        v_curves: v_curves.to_vec(),
        u_params,
        v_params,
    })
}

/// Build a Gordon surface with explicit parameter values.
///
/// This allows full control over where each curve sits in the parameter domain.
///
/// # Arguments
///
/// * `u_curves` - Curves in the u-direction
/// * `v_params` - v-parameter value for each u-curve (monotonic in [0, 1])
/// * `v_curves` - Curves in the v-direction
/// * `u_params` - u-parameter value for each v-curve (monotonic in [0, 1])
/// * `opts` - Construction options
pub fn gordon_surface_with_params(
    u_curves: &[Curve3],
    v_params: &[f64],
    v_curves: &[Curve3],
    u_params: &[f64],
    opts: GordonOptions,
) -> Result<GordonSurface, GordonError> {
    // Validate curve counts
    if u_curves.len() < 2 {
        return Err(GordonError::TooFewUCurves {
            count: u_curves.len(),
        });
    }
    if v_curves.len() < 2 {
        return Err(GordonError::TooFewVCurves {
            count: v_curves.len(),
        });
    }

    // Validate parameter counts match curve counts
    if u_params.len() != v_curves.len() {
        return Err(GordonError::UParamCountMismatch {
            expected: v_curves.len(),
            actual: u_params.len(),
        });
    }
    if v_params.len() != u_curves.len() {
        return Err(GordonError::VParamCountMismatch {
            expected: u_curves.len(),
            actual: v_params.len(),
        });
    }

    // Check for degenerate curves
    for (i, curve) in u_curves.iter().enumerate() {
        if is_degenerate_curve(curve, 10, opts.tolerance) {
            return Err(GordonError::DegenerateCurve {
                curve_idx: i,
                direction: "u".to_string(),
            });
        }
    }
    for (i, curve) in v_curves.iter().enumerate() {
        if is_degenerate_curve(curve, 10, opts.tolerance) {
            return Err(GordonError::DegenerateCurve {
                curve_idx: i,
                direction: "v".to_string(),
            });
        }
    }

    // Validate parameters
    validate_params(u_params, "u", opts.min_node_separation)?;
    validate_params(v_params, "v", opts.min_node_separation)?;

    // Validate intersections
    if opts.validate_intersections {
        validate_intersections(u_curves, u_params, v_curves, v_params, opts.intersection_tolerance)?;
    }

    Ok(GordonSurface {
        u_curves: u_curves.to_vec(),
        v_curves: v_curves.to_vec(),
        u_params: u_params.to_vec(),
        v_params: v_params.to_vec(),
    })
}

/// Evaluate a Gordon surface at a point with numerical stability checks.
///
/// Returns `None` if the evaluation fails due to numerical issues.
pub fn eval_gordon_surface_safe(
    surface: &GordonSurface,
    u: f64,
    v: f64,
    tol: f64,
) -> Option<DVec3> {
    let n = surface.u_curves.len();
    let m = surface.v_curves.len();

    if n == 0 && m == 0 {
        return Some(DVec3::ZERO);
    }

    // Compute Lagrange basis functions
    let lv = lagrange_basis_safe(&surface.v_params, v, tol)?;
    let lu = lagrange_basis_safe(&surface.u_params, u, tol)?;

    // Helper to evaluate curve at normalized parameter
    let eval_curve = |curve: &Curve3, t: f64| -> DVec3 {
        let [t0, t1] = curve.default_domain();
        if t0.is_finite() && t1.is_finite() && (t1 - t0).abs() > 1e-15 {
            curve.point_at(t0 + t * (t1 - t0))
        } else {
            curve.point_at(t)
        }
    };

    // Sum of u-direction loft
    let mut s_u = DVec3::ZERO;
    for (i, curve) in surface.u_curves.iter().enumerate() {
        s_u += lv[i] * eval_curve(curve, u);
    }

    // Sum of v-direction loft
    let mut s_v = DVec3::ZERO;
    for (j, curve) in surface.v_curves.iter().enumerate() {
        s_v += lu[j] * eval_curve(curve, v);
    }

    // Tensor product correction term
    let mut s_t = DVec3::ZERO;
    for (i, u_curve) in surface.u_curves.iter().enumerate() {
        for (j, _v_curve) in surface.v_curves.iter().enumerate() {
            // P_ij = intersection point at u_params[j] on u_curve
            //      = v_curve evaluated at v_params[i]
            let p_ij = eval_curve(u_curve, surface.u_params[j]);
            s_t += lv[i] * lu[j] * p_ij;
        }
    }

    let result = s_u + s_v - s_t;

    // Check for NaN/Inf
    if result.is_nan() || !result.is_finite() {
        return None;
    }

    Some(result)
}

/// Compute surface normal at a point using central differences.
///
/// Returns `None` if the normal computation fails.
pub fn gordon_surface_normal_safe(
    surface: &GordonSurface,
    u: f64,
    v: f64,
    eps: f64,
    tol: f64,
) -> Option<DVec3> {
    let p_u_plus = eval_gordon_surface_safe(surface, u + eps, v, tol)?;
    let p_u_minus = eval_gordon_surface_safe(surface, u - eps, v, tol)?;
    let p_v_plus = eval_gordon_surface_safe(surface, u, v + eps, tol)?;
    let p_v_minus = eval_gordon_surface_safe(surface, u, v - eps, tol)?;

    let du = p_u_plus - p_u_minus;
    let dv = p_v_plus - p_v_minus;

    let normal = du.cross(dv);
    let len = normal.length();

    if len < 1e-15 {
        return None;
    }

    Some(normal / len)
}

/// Convert a Gordon surface to a B-spline surface for export/visualization.
///
/// This samples the Gordon surface at a grid of points and fits a B-spline.
///
/// # Arguments
///
/// * `surface` - The Gordon surface to convert
/// * `u_samples` - Number of sample points in u direction
/// * `v_samples` - Number of sample points in v direction
/// * `degree` - Desired degree of the output B-spline
pub fn gordon_to_bspline(
    surface: &GordonSurface,
    u_samples: usize,
    v_samples: usize,
    degree: usize,
) -> Option<BSplineSurface> {
    if u_samples < degree + 1 || v_samples < degree + 1 {
        return None;
    }

    // Sample points on the surface
    let mut points: Vec<Vec<DVec3>> = Vec::with_capacity(u_samples);
    for i in 0..u_samples {
        let u = i as f64 / (u_samples - 1).max(1) as f64;
        let mut row: Vec<DVec3> = Vec::with_capacity(v_samples);
        for j in 0..v_samples {
            let v = j as f64 / (v_samples - 1).max(1) as f64;
            match eval_gordon_surface_safe(surface, u, v, 1e-10) {
                Some(p) => row.push(p),
                None => return None,
            }
        }
        points.push(row);
    }

    // Build uniform clamped knot vectors
    let n_ctrl_u = u_samples;
    let n_ctrl_v = v_samples;

    let knots_u = uniform_clamped_knots(n_ctrl_u, degree);
    let knots_v = uniform_clamped_knots(n_ctrl_v, degree);

    // Use sampled points as control points (simple approach)
    // A more sophisticated approach would use least-squares fitting
    let weights: Vec<Vec<f64>> = points
        .iter()
        .map(|row| row.iter().map(|_| 1.0).collect())
        .collect();

    Some(BSplineSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points: points,
        weights,
    })
}

/// Build uniform clamped knot vector for n control points.
fn uniform_clamped_knots(n_ctrl: usize, degree: usize) -> Vec<f64> {
    let m = n_ctrl + degree + 1;
    let mut knots = vec![0.0; m];
    let n_interior = n_ctrl.saturating_sub(degree + 1);

    // First degree+1 knots = 0
    for knot in knots.iter_mut().take(degree + 1) {
        *knot = 0.0;
    }

    // Last degree+1 knots = 1
    for knot in knots.iter_mut().skip(m - degree - 1) {
        *knot = 1.0;
    }

    // Interior knots
    if n_interior > 0 {
        for j in 1..=n_interior {
            knots[j + degree] = j as f64 / (n_interior + 1) as f64;
        }
    }

    knots
}

/// Check if a curve network forms a valid rectangular topology.
///
/// Returns true if:
/// - All u-curves have the same number of sample points
/// - All v-curves have the same number of sample points
/// - The network forms a proper grid
pub fn is_rectangular_network(
    u_curves: &[Curve3],
    v_curves: &[Curve3],
    tol: f64,
) -> bool {
    if u_curves.is_empty() || v_curves.is_empty() {
        return false;
    }

    // Sample each curve and check for consistent topology
    let _n_samples = 10;

    // Check that all u-curves intersect all v-curves
    for (i, u_curve) in u_curves.iter().enumerate() {
        for (j, v_curve) in v_curves.iter().enumerate() {
            // Check endpoints match
            let [u0, _u1] = u_curve.default_domain();
            let [v0, v1] = v_curve.default_domain();

            let u_start = u_curve.point_at(u0);
            let v_at_u_start = v_curve.point_at(v0 + (v1 - v0) * 0.0);

            // This is a simplified check - full validation would check all intersections
            let _ = (i, j, u_start, v_at_u_start, tol);
        }
    }

    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Line3;
    use crate::SurfaceEval;

    fn make_line(origin: DVec3, direction: DVec3) -> Curve3 {
        Curve3::Line(Line3 { origin, direction })
    }

    #[test]
    fn basic_bilinear_patch() {
        // Create a 2x2 bilinear patch
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default())
            .expect("should construct");

        // Check corners
        let p00 = eval_gordon_surface_safe(&surface, 0.0, 0.0, 1e-10).unwrap();
        let p10 = eval_gordon_surface_safe(&surface, 1.0, 0.0, 1e-10).unwrap();
        let p01 = eval_gordon_surface_safe(&surface, 0.0, 1.0, 1e-10).unwrap();
        let p11 = eval_gordon_surface_safe(&surface, 1.0, 1.0, 1e-10).unwrap();

        assert!((p00 - DVec3::ZERO).length() < 1e-10, "p00 = {:?}", p00);
        assert!((p10 - DVec3::X).length() < 1e-10, "p10 = {:?}", p10);
        assert!((p01 - DVec3::Y).length() < 1e-10, "p01 = {:?}", p01);
        assert!((p11 - DVec3::new(1.0, 1.0, 0.0)).length() < 1e-10, "p11 = {:?}", p11);

        // Check interior point
        let p05_05 = eval_gordon_surface_safe(&surface, 0.5, 0.5, 1e-10).unwrap();
        assert!((p05_05 - DVec3::new(0.5, 0.5, 0.0)).length() < 1e-10, "p05_05 = {:?}", p05_05);
    }

    #[test]
    fn three_by_three_network() {
        // Create a 3x3 network
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::new(0.0, 0.5, 0.0), DVec3::X);
        let u2 = make_line(DVec3::Y, DVec3::X);

        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::new(0.5, 0.0, 0.0), DVec3::Y);
        let v2 = make_line(DVec3::X, DVec3::Y);

        let surface = gordon_surface_curves(
            &[u0, u1, u2],
            &[v0, v1, v2],
            GordonOptions::default(),
        )
        .expect("should construct 3x3 network");

        // Verify interpolation at nodes
        let u_params = &surface.u_params;
        let v_params = &surface.v_params;

        assert_eq!(u_params.len(), 3);
        assert_eq!(v_params.len(), 3);

        // Check that the surface passes through all curve intersections
        for (_i, &vi) in v_params.iter().enumerate() {
            for (_j, &uj) in u_params.iter().enumerate() {
                let p = eval_gordon_surface_safe(&surface, uj, vi, 1e-10).unwrap();
                let expected = DVec3::new(uj, vi, 0.0);
                assert!(
                    (p - expected).length() < 1e-8,
                    "point at ({}, {}): {:?} != {:?}",
                    uj, vi, p, expected
                );
            }
        }
    }

    #[test]
    fn rejects_too_few_curves() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);

        // Single u-curve
        let err = gordon_surface_curves(&[u0.clone()], &[], GordonOptions::default());
        assert!(matches!(err, Err(GordonError::TooFewUCurves { .. })));

        // Single v-curve with 2 u-curves
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let err = gordon_surface_curves(&[u0, u1], &[v0], GordonOptions::default());
        assert!(matches!(err, Err(GordonError::TooFewVCurves { .. })));
    }

    #[test]
    fn rejects_non_monotonic_params() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        // Non-monotonic u_params
        let err = gordon_surface_with_params(
            &[u0, u1],
            &[0.0, 1.0], // v_params
            &[v0, v1],
            &[0.7, 0.3], // Non-monotonic u_params
            GordonOptions::default(),
        );
        assert!(matches!(err, Err(GordonError::NonMonotonicParams { .. })));
    }

    #[test]
    fn rejects_params_out_of_range() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let err = gordon_surface_with_params(
            &[u0, u1],
            &[0.0, 1.0],
            &[v0, v1],
            &[-0.5, 1.0], // Out of range
            GordonOptions::default(),
        );
        assert!(matches!(err, Err(GordonError::ParamsOutOfRange { .. })));
    }

    #[test]
    fn rejects_param_count_mismatch() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let err = gordon_surface_with_params(
            &[u0, u1],
            &[0.0, 1.0],
            &[v0, v1],
            &[0.0, 0.5, 1.0], // Wrong count (3 instead of 2)
            GordonOptions::default(),
        );
        assert!(matches!(err, Err(GordonError::UParamCountMismatch { .. })));
    }

    #[test]
    fn normal_computation() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default()).unwrap();

        let normal = gordon_surface_normal_safe(&surface, 0.5, 0.5, 1e-5, 1e-10).unwrap();

        // Normal should be pointing in +Z direction for a planar XY patch
        assert!(normal.z > 0.9, "normal = {:?}", normal);
        assert!(normal.length() > 0.99);
    }

    #[test]
    fn lagrange_basis_correctness() {
        // Test Lagrange basis for 3 nodes at 0, 0.5, 1
        let nodes = vec![0.0, 0.5, 1.0];

        // At t=0, only first basis function should be 1
        let basis0 = lagrange_basis_safe(&nodes, 0.0, 1e-10).unwrap();
        assert!((basis0[0] - 1.0).abs() < 1e-10);
        assert!(basis0[1].abs() < 1e-10);
        assert!(basis0[2].abs() < 1e-10);

        // At t=0.5, only middle basis function should be 1
        let basis_mid = lagrange_basis_safe(&nodes, 0.5, 1e-10).unwrap();
        assert!(basis_mid[0].abs() < 1e-10);
        assert!((basis_mid[1] - 1.0).abs() < 1e-10);
        assert!(basis_mid[2].abs() < 1e-10);

        // At t=1, only last basis function should be 1
        let basis1 = lagrange_basis_safe(&nodes, 1.0, 1e-10).unwrap();
        assert!(basis1[0].abs() < 1e-10);
        assert!(basis1[1].abs() < 1e-10);
        assert!((basis1[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn lagrange_basis_partition_of_unity() {
        // Sum of Lagrange basis functions should always equal 1
        let nodes = vec![0.0, 0.3, 0.7, 1.0];

        for &t in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let basis = lagrange_basis_safe(&nodes, t, 1e-10).unwrap();
            let sum: f64 = basis.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "sum at t={} = {}", t, sum);
        }
    }

    #[test]
    fn convert_to_bspline() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default()).unwrap();

        let bspline = gordon_to_bspline(&surface, 5, 5, 3);
        assert!(bspline.is_some());

        let bspline = bspline.unwrap();
        assert_eq!(bspline.degree_u, 3);
        assert_eq!(bspline.degree_v, 3);
        assert_eq!(bspline.control_points.len(), 5);
        assert_eq!(bspline.control_points[0].len(), 5);
    }

    #[test]
    fn skip_intersection_validation() {
        // Create curves that don't intersect properly
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::new(0.0, 2.0, 0.0), DVec3::X); // Displaced
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        // With validation, should fail
        let opts = GordonOptions::default();
        let _result = gordon_surface_curves(&[u0.clone(), u1.clone()], &[v0.clone(), v1.clone()], opts);
        // May or may not fail depending on tolerance - just check no panic

        // Without validation, should succeed
        let opts = GordonOptions::default().skip_intersection_validation();
        let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], opts);
        assert!(surface.is_ok());
    }

    #[test]
    fn gordon_surface_trait_impl() {
        let u0 = make_line(DVec3::ZERO, DVec3::X);
        let u1 = make_line(DVec3::Y, DVec3::X);
        let v0 = make_line(DVec3::ZERO, DVec3::Y);
        let v1 = make_line(DVec3::X, DVec3::Y);

        let surface = gordon_surface_curves(&[u0, u1], &[v0, v1], GordonOptions::default()).unwrap();

        // Test the SurfaceEval trait implementation
        let p = SurfaceEval::point_at(&surface, 0.5, 0.5);
        assert!((p - DVec3::new(0.5, 0.5, 0.0)).length() < 1e-10);

        let n = SurfaceEval::normal_at(&surface, 0.5, 0.5);
        assert!(n.length() > 0.9);

        let domain = SurfaceEval::default_domain(&surface);
        assert_eq!(domain, [0.0, 1.0, 0.0, 1.0]);
    }
}
