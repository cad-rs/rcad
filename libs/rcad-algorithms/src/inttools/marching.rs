use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;
use rcad_kernel::projection::closest_point_on_surface;

use crate::tolerance::*;

/// A numerically sampled intersection curve.
#[derive(Debug, Clone)]
#[derive(Default)]
pub struct SampledCurve {
    pub points: Vec<DVec3>,
    pub is_closed: bool,
    /// Diagnostic: number of oscillation events detected during marching.
    pub oscillation_count: usize,
    /// Diagnostic: whether step size was reduced during marching.
    pub step_reduced: bool,
}


/// Configuration for adaptive marching behavior.
#[derive(Debug, Clone, Copy)]
pub struct MarchingConfig {
    /// Initial step size for curve tracing.
    pub step_size: f64,
    /// Minimum allowed step size (for convergence failure fallback).
    pub min_step_size: f64,
    /// Maximum number of steps per direction.
    pub max_steps: usize,
    /// Maximum allowed oscillations before step reduction.
    pub max_oscillations: usize,
    /// Factor to reduce step size when oscillation is detected.
    pub step_reduction_factor: f64,
    /// Deflection tolerance (fleche) for chord error checking.
    ///
    /// When non-zero, the chord error between consecutive march points is
    /// estimated using tangent vectors (OCCT fleche formula). If the error
    /// exceeds this tolerance, the step is reduced.
    ///
    /// Reference: OCCT IntWalk_PWalking constructor parameter `Deflection`.
    /// Set to 0.0 to disable chord error checking (legacy behavior).
    pub deflection_tol: f64,
    /// Enable multi-scale seed detection.
    pub multiscale_seeds: bool,
}

impl Default for MarchingConfig {
    fn default() -> Self {
        Self {
            step_size: 0.1,
            min_step_size: TOLERANCE_LINEAR_ULTRA_STRICT,
            // 2000 steps per direction gives max ~400 units of arc per direction
            // at default step_size=0.1, enough for most large mechanical parts.
            // Callers can reduce this for small/simple geometry.
            max_steps: 2000,
            max_oscillations: 3,
            // 0.7 is gentler than 0.5: after 3 oscillations steps shrink to
            // 34% instead of 12.5% of the original, avoiding stalling on
            // slightly noisy tangent estimation near singularities.
            step_reduction_factor: 0.7,
            // 1e-8 matches OCCT's typical deflection tolerance for marching
            // (fleche ≈ 0.1 * Precision::Confusion()), ensuring chord error
            // is well below the linear tolerance.
            deflection_tol: TOLERANCE_MESH_LEGACY * 0.01,
            multiscale_seeds: false,
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// OCCT-aligned adaptive continuation marching configuration
// ──────────────────────────────────────────────────────────────────────────────
//
// Reference: OCCT IntWalk_PWalking.hxx, IntWalk_PWalking.cxx (Perform)
// OCCT source: `$OCCT_SRC/src/IntWalk/IntWalk_PWalking.hxx`
//
// OCCT IntWalk_PWalking computes surface-surface intersection curves by
// continuation marching. Its key parameters map to the fields below, which
// together implement curvature-aware step control (Task A/B), boundary-aware
// stepping (Task C), and inflection detection (Task D).
//
// OCCT line reference mapping (from IntWalk_PWalking.cxx):
//   L200-400  — Curvature-aware step control (Perform::ParamParamPerform)
//   L400-600  — Boundary-aware stepping (step clipping at surface bounds)
//   L600-800  — Inflection detection (zero-curvature / tangent reversal)
//   L1-200    — Seed point generation with Deflection + UVMaxStep
//   L3684-3685 — Fleche step formula: step = sqrt(8 * Deflection * R)
// ──────────────────────────────────────────────────────────────────────────────

/// OCCT-aligned configuration for adaptive continuation marching.
///
/// Mirrors the parameters used by OCCT IntWalk_PWalking:
/// - `Deflection` — chordal error tolerance (fleche), OCCT constructor param
/// - `UVMaxStep` — maximum step in UV space during continuation
/// - `StepMin` — minimum step to prevent infinite refinement
/// - `StepMax` — maximum step bound (typically UVMaxStep * 0.5)
/// - `FlecheTol` — finer chordal deviation proportional to Deflection/10
///
/// # OCCT Source Reference
///
/// IntWalk_PWalking.cxx L200-400 (ParamParamPerform):
///   The step at each iteration is bounded by:
///     h = min(UVMaxStep, sqrt(8 * Deflection * R))
///   where R is the minimum radius of curvature of either surface along
///   the marching direction.
///
/// IntWalk_PWalking.cxx L400-600:
///   After computing the candidate step, the next UV position is tested
///   against surface bounds. If it lies outside, the step is clipped to
///   the remaining boundary distance.
///
/// IntWalk_PWalking.cxx L600-800:
///   When the tangent direction reverses (dot product < 0) or the curvature
///   drops near zero, an inflection point is suspected — the step is halved.
///
/// # Alignment Status
///
/// ✅ OCCT-aligned: Fields map 1:1 to IntWalk_PWalking constructor params.
/// ✅ Partial: Marching loop integration of curvature step is new;
///     callers must pass this config to `march_intersection_with_config`.
#[derive(Debug, Clone, Copy)]
pub struct MarchingConfigOCCT {
    /// OCCT Deflection (deflection/fleche tolerance).
    /// Controls the maximum chord error between consecutive points.
    /// OCCT IntWalk_PWalking constructor: `Deflection`.
    pub deflection_tol: f64,
    /// OCCT UVMaxStep — maximum step in UV parameter space.
    /// OCCT IntWalk_PWalking::Perform: bounds the continuation step.
    /// Typical value: 0.1 (fraction of the shorter parametric side).
    pub uv_max_step: f64,
    /// Minimum step size — prevents infinite refinement at singularities.
    /// When curvature suggests a step below this, the step is clamped.
    pub min_step: f64,
    /// Maximum step size — clamp for curvature-based step estimates.
    /// OCCT equivalent: `UVMaxStep * 0.5` (safety factor).
    pub max_step: f64,
    /// OCCT Fleche — finer chordal deviation for step growth decisions.
    /// Typically `Deflection / 10`.
    /// OCCT IntWalk_PWalking::Perform: used to detect StepTooSmall regime.
    pub fleche_tol: f64,
}

impl Default for MarchingConfigOCCT {
    fn default() -> Self {
        Self {
            // deflection_tol ≈ 1e-6: matches TOLERANCE_MESH_LEGACY * 0.01
            // OCCT typically uses Precision::Confusion() (1e-7) for Deflection
            deflection_tol: TOLERANCE_MESH_LEGACY * 0.01,
            // UVMaxStep = 0.1: covers typical parametric domains (0..1, 0..2π).
            // OCCT sets UVMaxStep = max(min(domain_range / 20, 0.1), TOLERANCE_ABS).
            uv_max_step: 0.1,
            // min_step = 1e-10: prevents infinite refinement near singularities.
            min_step: TOLERANCE_LINEAR_ULTRA_STRICT,
            // max_step = 0.05: UVMaxStep * 0.5 as a safety factor.
            max_step: 0.05,
            // fleche_tol ≈ 1e-7: Deflection / 10, used for step growth decision.
            fleche_tol: TOLERANCE_MESH_LEGACY * 0.001,
        }
    }
}

impl MarchingConfigOCCT {
    /// Convert to a basic MarchingConfig for use with existing marching functions.
    ///
    /// Maps OCCT fields to the existing (simpler) config:
    /// - `step_size = uv_max_step`
    /// - `deflection_tol = deflection_tol`
    /// - `min_step_size = min_step`
    /// - `max_steps = 2000` (default, sufficient for most mechanical parts)
    pub fn to_basic_config(&self) -> MarchingConfig {
        MarchingConfig {
            step_size: self.uv_max_step,
            min_step_size: self.min_step,
            deflection_tol: self.deflection_tol,
            ..Default::default()
        }
    }
}

/// Result of adaptive sampling density calculation.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveSampling {
    /// Number of samples in u-direction.
    pub n_u: usize,
    /// Number of samples in v-direction.
    pub n_v: usize,
    /// Estimated characteristic length for step size.
    pub characteristic_length: f64,
}

/// OCCT IntWalk_StatusDeflection equivalent for chord error checking.
///
/// Returned by [`test_deflection`] to indicate how the stepping should adapt.
///
/// Reference: OCCT IntWalk_StatusDeflection.hxx
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeflectionStatus {
    /// Step is acceptable — the point can be added to the curve.
    Ok,
    /// Step too large — chord error exceeds deflection tolerance.
    StepTooLarge,
    /// Step too small — chord error is well below tolerance; can grow.
    StepTooSmall,
    /// Tangent reversal — potential inflection point detected.
    Inflection,
    /// Consecutive points are nearly coincident.
    ConfusedPoint,
}

/// Compute fleche (chord error) between two consecutive intersection points.
///
/// The fleche estimates the maximum deviation of the actual intersection curve
/// from the chord connecting two sampled points. It uses the unit tangent
/// vectors at both endpoints to infer curvature.
///
/// # Formula
///
/// For a circular arc of radius R spanning chord length d, the sagitta is:
///   `e = R - sqrt(R² - (d/2)²)`
/// The turning angle θ = d/R, and for small θ:
///   `e ≈ d * θ / 8`
/// Since `|ΔT| = |prev_tg - cur_tg| ≈ θ` (for unit tangents):
///   `e ≈ d * |ΔT| / 8`
///
/// Reference: OCCT IntWalk_PWalking::TestDeflection, L3684-3685
///   `FlecheCourante = sqrt(|ΔT|² * d²) / 8`
#[inline]
pub fn compute_fleche(prev_tangent: DVec3, cur_tangent: DVec3, chord_length: f64) -> f64 {
    let diff_len = (prev_tangent - cur_tangent).length();
    // sqrt(|ΔT|² * d²) / 8 = |ΔT| * d / 8
    diff_len * chord_length / 8.0
}

/// OCCT-style deflection test for two consecutive points on the intersection curve.
///
/// Evaluates the chord error (fleche) between the previous accepted point and
/// a candidate new point. Returns a status and a suggested step size for the
/// next iteration.
///
/// # OCCT Logic (IntWalk_PWalking::TestDeflection, L3407-3916)
///
/// 1. **Confused points** — If the chord is shorter than [`TOLERANCE_ABS`],
///    the points are nearly coincident. Increase step.
///
/// 2. **Inflection** — If the dot product of successive tangents is negative,
///    the direction reversed (potential inflection point). Halve step.
///
/// 3. **Fleche classification**:
///    - `fleche <= deflection_tol * 0.5` → StepTooSmall (can increase)
///    - `fleche > deflection_tol` → StepTooLarge (must reduce)
///    - Otherwise → Ok (sweet spot)
///
/// # Arguments
///
/// * `chord_length_sq` — Squared 3D distance between the two points.
/// * `chord_length` — 3D distance between the two points.
/// * `prev_tangent` — Unit tangent direction at the previous point.
/// * `cur_tangent` — Unit tangent direction at the candidate point.
/// * `step_size` — Current step size used for this stride.
/// * `deflection_tol` — Maximum allowed fleche (chord error).
/// * `min_step` — Minimum step size (clamp for reductions).
/// * `initial_step` — Initial (undamped) step size (ceiling for growth).
///
/// Returns `(status, suggested_step)`.
#[inline]
pub fn test_deflection(
    chord_length_sq: f64,
    chord_length: f64,
    prev_tangent: DVec3,
    cur_tangent: DVec3,
    step_size: f64,
    deflection_tol: f64,
    min_step: f64,
    initial_step: f64,
) -> (DeflectionStatus, f64) {
    // Guard: disable fleche check when deflection_tol is zero (legacy mode).
    if deflection_tol <= 0.0 {
        return (DeflectionStatus::Ok, step_size);
    }

    // Step 1: confused points — chord is shorter than absolute tolerance.
    if chord_length_sq < TOLERANCE_ABS_SQ {
        let new_step = (step_size * 1.5).min(initial_step);
        return (DeflectionStatus::ConfusedPoint, new_step);
    }

    // Step 2: inflection detection — tangent direction reversal.
    let cos_between = prev_tangent.dot(cur_tangent);
    if cos_between < 0.0 {
        // Tangent direction changed sign — potential inflection.
        // OCCT L3461: halve all steps.
        let new_step = (step_size * 0.5).max(min_step);
        return (DeflectionStatus::Inflection, new_step);
    }

    // Step 3: compute fleche (chord error estimate).
    let fleche = compute_fleche(prev_tangent, cur_tangent, chord_length);

    // Step 4: classify fleche against deflection tolerance.
    if fleche <= deflection_tol * 0.5 {
        // Step is too small — fleche is well below tolerance.
        // Compute a growth ratio: we want the next step to produce
        // fleche ≈ deflection_tol. Since fleche ∝ step² (for constant
        // curvature), ratio ≈ sqrt(deflection_tol / fleche).
        // OCCT L3691-3692: Ratio = 0.5 * (fleche / FlecheCourante)
        let ratio = 0.5 * (deflection_tol / fleche.max(f64::MIN_POSITIVE));
        let new_step = (step_size * ratio.clamp(1.0, 4.0)).min(initial_step);
        (DeflectionStatus::StepTooSmall, new_step)
    } else if fleche > deflection_tol {
        // Step is too large — chord error exceeds tolerance.
        // OCCT L3783-3785: Ratio = fleche / FlecheCourante
        let ratio = deflection_tol / fleche;
        let new_step = (step_size * ratio).max(min_step);
        (DeflectionStatus::StepTooLarge, new_step)
    } else {
        // deflection_tol * 0.5 < fleche <= deflection_tol
        // Step is in the sweet spot — gentle adjustment.
        // OCCT L3796: Ratio = 0.75 * (fleche / FlecheCourante)
        let ratio = 0.75 * (deflection_tol / fleche.max(f64::MIN_POSITIVE));
        let new_step = (step_size * ratio.clamp(0.5, 1.5)).clamp(min_step, initial_step);
        (DeflectionStatus::Ok, new_step)
    }
}

/// Compute the 3D curvature of a surface curve `r(t) = S(uv + t * dir)` at t=0.
///
/// Uses central finite differences for the first and second derivatives:
///   `r'(0) ≈ (r(h) - r(-h)) / (2h)`
///   `r''(0) ≈ (r(h) - 2*r(0) + r(-h)) / h^2`
///   `k = |r' × r''| / |r'|^3`
///
/// Returns `(curvature, radius)` where `radius = 1/curvature`.
/// When curvature is below `TOLERANCE_ABS`, returns `(0.0, INFINITY)`.
///
/// # Reference
///
/// OCCT IntWalk_PWalking::Perform — uses surface curvature estimation at each
/// step for the fleche computation (L3684-3685).
fn compute_curvature_along_dir(
    surf_eval: &impl Fn(DVec2) -> DVec3,
    uv: DVec2,
    dir: DVec2,
) -> (f64, f64) {
    let dir_len = dir.length();
    if dir_len < TOLERANCE_ABS {
        return (0.0, f64::INFINITY);
    }
    let step_uv = dir / dir_len;
    // Epsilon for finite differences: small enough for curvature accuracy,
    // large enough to avoid floating-point noise.
    let h = 1e-5_f64.max(TOLERANCE_ABS * 100.0).min(1e-2);

    let p0 = surf_eval(uv);
    let p1 = surf_eval(uv + step_uv * h);
    let p_m1 = surf_eval(uv - step_uv * h);

    if !p0.is_finite() || !p1.is_finite() || !p_m1.is_finite() {
        return (0.0, f64::INFINITY);
    }

    // First derivative (central difference).
    let r1 = (p1 - p_m1) / (2.0 * h);
    let r1_len_sq = r1.length_squared();
    if r1_len_sq < TOLERANCE_ABS_SQ {
        return (0.0, f64::INFINITY);
    }

    // Second derivative.
    let r2 = (p1 - p0 * 2.0 + p_m1) / (h * h);

    // 3D space curve curvature: k = |r' × r''| / |r'|^3.
    let cross_len = r1.cross(r2).length();
    let r1_len = r1_len_sq.sqrt();
    let curvature = cross_len / (r1_len * r1_len * r1_len);

    if curvature < TOLERANCE_ABS {
        (0.0, f64::INFINITY)
    } else {
        let radius = 1.0 / curvature;
        (curvature, radius)
    }
}

/// OCCT-aligned curvature-aware step size computation.
///
/// Estimates the radius of curvature of the surface at `uv` along direction
/// `dir` and returns the maximum step that keeps chord error (fleche) below
/// `deflection`.
///
/// # Formula (OCCT L3684-3685)
///
/// For a circular arc of radius R spanning chord of length h, the sagitta is:
///   `fleche = R - sqrt(R^2 - (h/2)^2) ≈ h^2 / (8*R)`  for h << R
///
/// Solving for h:
///   `h = sqrt(8 * deflection * R)`
///
/// where R = 1/|k| and k is the curvature of the surface curve along `dir`.
///
/// # Arguments
///
/// * `surf_eval` — Surface point evaluator: `(u, v) -> 3D point`.
/// * `uv` — Current UV parameter on the surface.
/// * `dir` — Marching direction in UV space.
/// * `deflection` — Maximum allowed chord error (fleche tolerance).
/// * `max_step` — Upper bound for the returned step (clamp).
/// * `min_step` — Lower bound for the returned step (clamp).
///
/// # Returns
///
/// The suggested step size (in UV space) clamped to `[min_step, max_step]`.
///
/// # Reference
///
/// OCCT IntWalk_PWalking::Perform (L3684-3685):
///   ```cpp
///   step = sqrt(8 * Deflection / MaxCurvature)
///   ```
/// where `MaxCurvature` is the maximum surface curvature at the point.
pub fn compute_fleche_step(
    surf_eval: &impl Fn(DVec2) -> DVec3,
    uv: DVec2,
    dir: DVec2,
    deflection: f64,
    max_step: f64,
    min_step: f64,
) -> f64 {
    // Guard: if deflection is non-positive, no curvature-based step control.
    if deflection <= 0.0 || dir.length_squared() < TOLERANCE_ABS_SQ {
        return max_step;
    }

    let (_curvature, radius) = compute_curvature_along_dir(surf_eval, uv, dir);

    if !radius.is_finite() || radius < TOLERANCE_ABS {
        return max_step;
    }

    let fleche_step = (8.0 * deflection * radius).sqrt();

    fleche_step.clamp(min_step, max_step)
}

/// OCCT-aligned adaptive sampling density based on surface curvature.
///
/// Evaluates surface curvature at sample points across the parametric domain
/// and uses the minimum radius of curvature to compute a grid density suitable
/// for seed-finding via sign-change detection. Unlike the heuristic-based
/// `adaptive_sampling_density`, this function directly measures curvature
/// and applies the OCCT UVMaxStep heuristic:
///
///   `maxStep = min(domainSize / 20, curvature_adjusted_step)`
///
/// where `curvature_adjusted_step = sqrt(8 * deflection * R_min)` for the
/// minimum radius of curvature R_min across the surface.
///
/// # Arguments
///
/// * `surface` — The surface to evaluate.
/// * `base_density` — Minimum grid density (fallback when curvature is low).
/// * `deflection` — Chord error tolerance for step computation.
/// * `sample_n` — Number of sample points per direction for curvature estimation.
///
/// # Reference
///
/// OCCT IntWalk_PWalking::ParamParamPerform:
///   Uses Deflection and UVMaxStep to bound the continuation step.
///   The grid density is then: n = ceil(domainRange / maxStep).
pub fn adaptive_sampling_density_occt(
    surface: &Surface3,
    base_density: usize,
    deflection: f64,
    sample_n: usize,
) -> AdaptiveSampling {
    let dom = surface.default_domain();
    let du = dom[1] - dom[0];
    let dv = dom[3] - dom[2];
    let domain_size = du.max(dv).abs();

    if !domain_size.is_finite() || domain_size < TOLERANCE_ABS {
        return AdaptiveSampling {
            n_u: base_density,
            n_v: base_density,
            characteristic_length: 1.0,
        };
    }

    // OCCT UVMaxStep heuristic: maxStep = domainSize/20 as the default
    // (OCCT IntWalk_PWalking::Perform initializes step from UVMaxStep / 20).
    let max_step_domain = (domain_size / 20.0).max(TOLERANCE_ABS);

    // Sample curvature across the surface to find the minimum radius.
    let mut min_radius = f64::INFINITY;
    let sample_n = sample_n.max(3);
    let surf_fn = |uv: DVec2| surface.point_at(uv.x, uv.y);
    for i in 0..sample_n {
        for j in 0..sample_n {
            let u = dom[0] + du * i as f64 / (sample_n - 1) as f64;
            let v = dom[2] + dv * j as f64 / (sample_n - 1) as f64;
            let uv = DVec2::new(u, v);
            let p0 = surface.point_at(u, v);
            if !p0.is_finite() {
                continue;
            }
            // Check curvature in four UV directions to find the tightest radius.
            for dir in &[
                DVec2::new(1.0, 0.0),
                DVec2::new(0.0, 1.0),
                DVec2::new(1.0, 1.0),
                DVec2::new(1.0, -1.0),
            ] {
                let (_curvature, radius) = compute_curvature_along_dir(&surf_fn, uv, *dir);
                if radius.is_finite() && radius > TOLERANCE_ABS && radius < min_radius {
                    min_radius = radius;
                }
            }
        }
    }

    // Compute curvature-adjusted step from the minimum radius.
    let curvature_adjusted_step = if min_radius.is_finite() && min_radius < f64::MAX / 8.0 {
        (8.0 * deflection * min_radius).sqrt()
    } else {
        f64::INFINITY
    };

    let max_step = max_step_domain.min(curvature_adjusted_step).max(TOLERANCE_ABS);

    // Grid density: ceil(domainRange / maxStep).
    let n_u = (du.abs() / max_step).ceil() as usize;
    let n_v = (dv.abs() / max_step).ceil() as usize;

    let characteristic_length = (max_step * 5.0)
        .min(domain_size * 0.1)
        .max(TOLERANCE_ABS);

    AdaptiveSampling {
        n_u: n_u.max(base_density / 4).min(base_density * 4),
        n_v: n_v.max(base_density / 4).min(base_density * 4),
        characteristic_length,
    }
}
pub fn adaptive_sampling_density(surface: &Surface3, base_density: usize) -> AdaptiveSampling {
    let default = AdaptiveSampling {
        n_u: base_density,
        n_v: base_density,
        characteristic_length: 1.0,
    };

    match surface {
        Surface3::Cylinder(c) => {
            // Cylinder: u = azimuth, v = height
            // Azimuth should be proportional to circumference (2πR)
            // Height should cover the expected range
            let circumference = std::f64::consts::TAU * c.radius;
            let n_u = (base_density as f64 * (circumference / c.radius).sqrt()).ceil() as usize;
            let n_v = base_density;
            AdaptiveSampling {
                n_u: n_u.max(base_density / 2).min(base_density * 2),
                n_v: n_v.max(base_density / 2),
                characteristic_length: c.radius * 0.1,
            }
        }
        Surface3::Sphere(s) => {
            // Sphere: uniform sampling based on radius
            let n = (base_density as f64 * (s.radius / 1.0).sqrt()).ceil() as usize;
            AdaptiveSampling {
                n_u: n.max(base_density),
                n_v: n.max(base_density),
                characteristic_length: s.radius * 0.1,
            }
        }
        Surface3::Torus(t) => {
            // Torus: major radius for u, minor radius for v
            let ratio = t.major_radius / t.minor_radius.max(TOLERANCE_LINEAR_ULTRA_STRICT);
            let n_u = (base_density as f64 * ratio.sqrt()).ceil() as usize;
            let n_v = base_density;
            AdaptiveSampling {
                n_u: n_u.max(base_density).min(base_density * 3),
                n_v: n_v.max(base_density / 2),
                characteristic_length: t.minor_radius * 0.1,
            }
        }
        Surface3::Cone(c) => {
            // Cone: similar to cylinder but with varying radius
            let avg_radius = c.radius * 0.5; // approximate average
            let n_u = (base_density as f64 * (avg_radius / 1.0).sqrt()).ceil() as usize;
            let n_v = base_density;
            AdaptiveSampling {
                n_u: n_u.max(base_density / 2),
                n_v: n_v.max(base_density / 2),
                characteristic_length: avg_radius * 0.1,
            }
        }
        Surface3::BSpline(bs) => {
            // BSpline: estimate from control point bounding box
            let bbox = estimate_bspline_bbox(bs);
            let max_extent = (bbox.1 - bbox.0).max_element();
            let n = (base_density as f64 * (max_extent / 1.0).sqrt()).ceil() as usize;
            AdaptiveSampling {
                n_u: n.max(base_density),
                n_v: n.max(base_density),
                characteristic_length: max_extent * 0.05,
            }
        }
        _ => default,
    }
}

/// Estimate bounding box of a BSpline surface from control points.
fn estimate_bspline_bbox(bs: &BSplineSurface) -> (DVec3, DVec3) {
    let mut min_pt = DVec3::splat(f64::INFINITY);
    let mut max_pt = DVec3::splat(f64::NEG_INFINITY);

    for row in &bs.control_points {
        for pt in row {
            min_pt = min_pt.min(*pt);
            max_pt = max_pt.max(*pt);
        }
    }

    if !min_pt.is_finite() {
        min_pt = DVec3::splat(-1.0);
    }
    if !max_pt.is_finite() {
        max_pt = DVec3::splat(1.0);
    }

    (min_pt, max_pt)
}

/// Coarse UV grid search (5×5 over [0,1]²) to find the closest surface sample.
/// Returns (u, v) of the closest grid point.
fn closest_uv_coarse(surface: &Surface3, point: DVec3) -> (f64, f64) {
    const N: usize = 5;
    let mut best_u = 0.5_f64;
    let mut best_v = 0.5_f64;
    let mut best_dist_sq = f64::MAX;
    for i in 0..N {
        for j in 0..N {
            let u = i as f64 / (N - 1) as f64;
            let v = j as f64 / (N - 1) as f64;
            let p = surface.point_at(u, v);
            let d = (p - point).length_squared();
            if d < best_dist_sq {
                best_dist_sq = d;
                best_u = u;
                best_v = v;
            }
        }
    }
    (best_u, best_v)
}

/// Evaluate the implicit function F(P) for a surface: F=0 on surface.
pub fn surface_implicit(surface: &Surface3, point: DVec3) -> f64 {
    match surface {
        Surface3::Plane(p) => (point - p.origin).dot(p.normal),
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            perp.length() - c.radius
        }
        Surface3::Sphere(s) => (point - s.center).length() - s.radius,
        Surface3::Cone(c) => {
            let axis = c.axis_dir();
            let v = point - c.apex;
            let along = v.dot(axis);
            let perp_len = (v - axis * along).length();
            perp_len - c.radius_at_axial(along)
        }
        Surface3::Torus(t) => {
            let v = point - t.center;
            let along = v.dot(t.axis);
            let perp = v - t.axis * along;
            let perp_len = perp.length();
            let d = perp_len - t.major_radius;
            (d * d + along * along).sqrt() - t.minor_radius
        }
        _ => {
            // Closest-point signed distance: use numerical projection (uniform
            // sampling + Newton refinement) for accurate UV on BSpline surfaces.
            let proj = closest_point_on_surface(surface, point, 8);
            let closest = proj.point;
            let normal = surface.normal_at(proj.params.0, proj.params.1);
            let n_len = normal.length();
            if n_len < TOLERANCE_FLOAT_LOOSE {
                return (point - closest).length();
            }
            (point - closest).dot(normal / n_len)
        }
    }
}

/// Compute the gradient ∇F at a point for a surface.
fn surface_gradient(surface: &Surface3, point: DVec3) -> DVec3 {
    match surface {
        Surface3::Plane(p) => p.normal,
        Surface3::Cylinder(c) => {
            let v = point - c.origin;
            let along = v.dot(c.axis);
            let perp = v - c.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            perp / perp_len
        }
        Surface3::Sphere(s) => {
            let v = point - s.center;
            let len = v.length();
            if len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            v / len
        }
        Surface3::Cone(c) => {
            let axis = c.axis_dir();
            let v = point - c.apex;
            let along = v.dot(axis);
            let perp = v - axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            let tan_a = c.half_angle_rad.tan();
            perp / perp_len - axis * tan_a
        }
        Surface3::Torus(t) => {
            let v = point - t.center;
            let along = v.dot(t.axis);
            let perp = v - t.axis * along;
            let perp_len = perp.length();
            if perp_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            let tube_center = t.center + perp / perp_len * t.major_radius;
            let tv = point - tube_center;
            let tv_len = tv.length();
            if tv_len < TOLERANCE_ABS {
                return DVec3::ZERO;
            }
            tv / tv_len
        }
        _ => {
            let proj = closest_point_on_surface(surface, point, 8);
            let normal = surface.normal_at(proj.params.0, proj.params.1);
            let n_len = normal.length();
            if n_len < TOLERANCE_FLOAT_LOOSE {
                return DVec3::ZERO;
            }
            normal / n_len
        }
    }
}

/// Project a point onto a surface using Newton iteration.
pub fn project_onto_surface(surface: &Surface3, point: DVec3, max_iter: usize) -> DVec3 {
    let mut p = point;
    for _ in 0..max_iter {
        let f = surface_implicit(surface, p);
        if f.abs() < TOLERANCE_ABS {
            break;
        }
        let g = surface_gradient(surface, p);
        let g_len_sq = g.length_squared();
        if g_len_sq < TOLERANCE_ABS * TOLERANCE_ABS {
            break;
        }
        p -= g * (f / g_len_sq);
    }
    p
}

/// Project a point onto the intersection of two surfaces within `tol`.
///
/// Two-surface Newton projection with residual **`tol`** on both implicits (`|f| < tol`).
/// Numeric IntSS uses the same order of magnitude as **`refine_tol`** when refining grid seeds ([`crate::inttools::intss::intersect_surfaces_with_density_tol`] pathway).
pub fn project_onto_intersection_tol(
    s1: &Surface3,
    s2: &Surface3,
    point: DVec3,
    tol: f64,
) -> DVec3 {
    let tol = tol.abs().max(TOLERANCE_LEN_MIN);
    let tol2 = tol * tol;
    let mut p = point;
    for _ in 0..50 {
        let f1 = surface_implicit(s1, p);
        let f2 = surface_implicit(s2, p);
        if f1.abs() < tol && f2.abs() < tol {
            break;
        }
        let g1 = surface_gradient(s1, p);
        let g2 = surface_gradient(s2, p);

        // Solve 2x2 system: move by λ1*g1 + λ2*g2 to satisfy both constraints
        let a11 = g1.dot(g1);
        let a12 = g1.dot(g2);
        let a22 = g2.dot(g2);
        let det = a11 * a22 - a12 * a12;
        if det.abs() < tol2 {
            // Degenerate — just project onto each surface alternately
            p = project_onto_surface(s1, p, 5);
            p = project_onto_surface(s2, p, 5);
            continue;
        }
        let lambda1 = (a22 * f1 - a12 * f2) / det;
        let lambda2 = (a11 * f2 - a12 * f1) / det;
        p -= g1 * lambda1 + g2 * lambda2;
    }
    p
}

/// Project onto the intersection using [`TOLERANCE_ABS`] as the residual target.
fn project_onto_intersection(s1: &Surface3, s2: &Surface3, point: DVec3) -> DVec3 {
    project_onto_intersection_tol(s1, s2, point, TOLERANCE_ABS)
}

/// Find seed points for intersection curve marching by sampling one surface.
///
/// Samples the second surface's implicit function at the given 3D points and
/// detects sign changes along the ordered sequence. Each sign change indicates
/// a crossing of the intersection curve between the two consecutive points.
///
/// This corresponds to OCCT's initial grid sampling phase in
/// IntWalk_PWalking::Perform (L1-200), where the first surface is sampled on
/// an N×N grid and the second surface's distance function is evaluated.  OCCT
/// then calls `IntWalk_PWalking::ParamParamPerform` (L200-800) for each seed
/// to trace the full continuation curve.
///
/// # OCCT Reference
///
/// IntWalk_PWalking.cxx (Perform):
///   L1-200 — Seed point generation:
///     - Sample surface 1 on an N×N UV grid.
///     - Compute approximate distance to surface 2 at each sample.
///     - Detect sign-change edges (grid cells where opposite corners have
///       opposite signs) and linearly interpolate crossing points.
///     - Project interpolated points onto both surfaces via Newton iteration
///       (equivalent to `project_onto_intersection`).
///
/// This function performs step 2-3 using pre-computed sample points and
/// aligns with OCCT's `IntWalk_PWalking::Perform(Handle(Adaptor3d_Surface)&)`.
pub fn find_seed_points(s1: &Surface3, s2: &Surface3, sample_points: &[DVec3]) -> Vec<DVec3> {
    let mut seeds = Vec::new();

    // Look for sign changes of F2 along the sample points
    let values: Vec<f64> = sample_points
        .iter()
        .map(|&p| surface_implicit(s2, p))
        .collect();

    for i in 0..values.len().saturating_sub(1) {
        if values[i] * values[i + 1] < 0.0 {
            // Sign change — interpolate
            let t = values[i] / (values[i] - values[i + 1]);
            let p = sample_points[i] + (sample_points[i + 1] - sample_points[i]) * t;
            let seed = project_onto_intersection(s1, s2, p);
            seeds.push(seed);
        }
    }

    seeds
}

/// Like `find_seed_points` but treats `sample_points` as a `n_u × n_v` grid
/// (row-major: index = iu * n_v + iv) and checks sign changes along BOTH the
/// u-direction and v-direction edges. This avoids missing seeds when the
/// intersection curve runs along one of the grid directions.
///
/// # OCCT Reference
///
/// IntWalk_PWalking.cxx (Perform, L1-200):
///   OCCT samples surface 1 on an N_u × N_v grid (from UVMaxStep) and tests
///   each cell edge for sign changes of the distance to surface 2.  This is
///   the 2D grid analog of `IntWalk_PWalking::Perform` with `UVMaxStep`
///   controlling the grid spacing.
pub fn find_seed_points_grid(
    s1: &Surface3,
    s2: &Surface3,
    sample_points: &[DVec3],
    n_u: usize,
    n_v: usize,
) -> Vec<DVec3> {
    assert_eq!(sample_points.len(), n_u * n_v, "grid size mismatch");
    let mut seeds = Vec::new();

    let values: Vec<f64> = sample_points
        .iter()
        .map(|&p| surface_implicit(s2, p))
        .collect();

    let idx = |iu: usize, iv: usize| iu * n_v + iv;

    // Check u-direction edges (vary iv, fixed iu)
    for iu in 0..n_u {
        for iv in 0..n_v.saturating_sub(1) {
            let a = idx(iu, iv);
            let b = idx(iu, iv + 1);
            if values[a] * values[b] < 0.0 {
                let t = values[a] / (values[a] - values[b]);
                let p = sample_points[a].lerp(sample_points[b], t);
                seeds.push(project_onto_intersection(s1, s2, p));
            }
        }
    }

    // Check v-direction edges (vary iu, fixed iv)
    for iu in 0..n_u.saturating_sub(1) {
        for iv in 0..n_v {
            let a = idx(iu, iv);
            let b = idx(iu + 1, iv);
            if values[a] * values[b] < 0.0 {
                let t = values[a] / (values[a] - values[b]);
                let p = sample_points[a].lerp(sample_points[b], t);
                seeds.push(project_onto_intersection(s1, s2, p));
            }
        }
    }

    seeds
}

/// Multi-scale seed point detection with deduplication.
/// Runs seed detection at multiple grid resolutions and merges results.
pub fn find_seed_points_multiscale(
    s1: &Surface3,
    s2: &Surface3,
    sampler: impl Fn(usize, usize) -> Vec<DVec3>,
    scales: &[usize],
    dedup_tolerance: f64,
) -> Vec<DVec3> {
    let mut all_seeds = Vec::new();
    let dedup_tol_sq = dedup_tolerance * dedup_tolerance;

    for &n in scales {
        let n_u = n;
        let n_v = n;
        let samples = sampler(n_u, n_v);
        let seeds = find_seed_points_grid(s1, s2, &samples, n_u, n_v);

        for seed in seeds {
            // Deduplicate: skip if too close to an existing seed
            let is_dup = all_seeds.iter().any(|s: &DVec3| (*s - seed).length_squared() < dedup_tol_sq);
            if !is_dup {
                all_seeds.push(seed);
            }
        }
    }

    all_seeds
}

/// March an intersection curve starting from a seed point.
/// Traces in both directions along the curve until it returns to start
/// (closed) or exits bounds.
///
/// # OCCT Reference
///
/// IntWalk_PWalking::ParamParamPerform (L200-800):
///   OCCT traces the intersection curve in a single direction from a seed
///   point, using curvature-aware step control (L200-400), boundary clipping
///   (L400-600), and inflection detection (L600-800).  This function splits
///   the trace into two one-direction marches (forward + backward) and merges
///   the results, which is equivalent to OCCT's two-direction continuation
///   starting from each seed.
///
///   OCCT IntWalk_PWalking.cxx (Perform, L~100-200):
///     After finding seeds, OCCT calls `Perform(Handle(Adaptor3d_Surface)&)`
///     which invokes `ParamParamPerform` for each seed, tracing the curve
///     until the step falls below resolution or the curve exits bounds.
pub fn march_intersection(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: impl Fn(DVec3) -> bool,
) -> SampledCurve {
    let config = MarchingConfig {
        step_size,
        max_steps,
        ..Default::default()
    };
    // Use the simpler implementation without Clone bound
    let mut result = SampledCurve::default();

    // Try forward direction
    let forward = march_one_direction_monitored_simple(
        s1, s2, seed, config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor, config.deflection_tol,
    );

    // Try backward direction
    let backward = march_one_direction_monitored_simple(
        s1, s2, seed, -config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor, config.deflection_tol,
    );

    // Combine: reverse backward, then append forward (excluding duplicate seed)
    result.points = backward.points.into_iter().rev().collect();
    if !forward.points.is_empty() {
        result.points.extend(forward.points.into_iter().skip(1));
    }

    // Check closure
    result.is_closed = result.points.len() > 2
        && points_coincide(result.points[0], *result.points.last().unwrap());

    if result.is_closed {
        result.points.pop();
    }

    result
}

/// March an intersection curve with full configuration and convergence monitoring.
///
/// # OCCT Reference
///
/// This is the high-level entry point equivalent to OCCT's
/// `IntWalk_PWalking::Perform` + `ParamParamPerform` combined.
/// - Seed point generation (called before this function) maps to OCCT L1-200.
/// - Bi-directional continuation maps to OCCT L200-800 (ParamParamPerform).
/// - `MarchingConfig.deflection_tol` maps to OCCT `Deflection` parameter.
/// - `MarchingConfig.step_size` maps to OCCT `UVMaxStep` parameter.
///
/// For full OCCT alignment, use `march_intersection_with_config_occt` which
/// accepts `MarchingConfigOCCT` for curvature-aware stepping.
pub fn march_intersection_with_config(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    config: &MarchingConfig,
    bounds_check: impl Fn(DVec3) -> bool,
) -> SampledCurve {
    let mut result = SampledCurve::default();

    // Try forward direction with convergence monitoring
    let forward = march_one_direction_monitored_simple(
        s1, s2, seed, config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor, config.deflection_tol,
    );
    result.oscillation_count += forward.oscillation_count;
    result.step_reduced = result.step_reduced || forward.step_reduced;

    // Try backward direction
    let backward = march_one_direction_monitored_simple(
        s1, s2, seed, -config.step_size, config.max_steps,
        &bounds_check, config.max_oscillations, config.min_step_size,
        config.step_reduction_factor, config.deflection_tol,
    );
    result.oscillation_count += backward.oscillation_count;
    result.step_reduced = result.step_reduced || backward.step_reduced;

    // Combine: reverse backward, then append forward (excluding duplicate seed)
    result.points = backward.points.into_iter().rev().collect();
    if !forward.points.is_empty() {
        // Skip the seed which is duplicated
        result.points.extend(forward.points.into_iter().skip(1));
    }

    // Check closure
    result.is_closed = result.points.len() > 2
        && points_coincide(result.points[0], *result.points.last().unwrap());

    if result.is_closed {
        result.points.pop();
    }

    result
}

/// Result of monitored single-direction marching.
struct MonitoredMarchResult {
    points: Vec<DVec3>,
    oscillation_count: usize,
    step_reduced: bool,
}

/// March in one direction with oscillation detection and step reduction.
/// Simple version without Clone bound on bounds_check.
///
/// # OCCT Reference
///
/// IntWalk_PWalking::ParamParamPerform (L200-800):
///   This function implements the core continuation marching loop. The OCCT
///   algorithm proceeds as follows:
///
///   L200-400 — Curvature-aware step control:
///     OCCT computes the next step as:
///       h = min(UVMaxStep, sqrt(8 * Deflection * R))
///     where R is the minimum radius of curvature of either surface along
///     the marching direction.  `compute_fleche_step` implements this formula.
///     rcad currently uses fixed step size with the fleche-based deflection
///     check (`test_deflection`) as a post-hoc validation.
///
///   L400-600 — Boundary-aware stepping:
///     OCCT projects the candidate next UV onto the surface bounds. If the
///     next point lies outside the natural bounds (u∈[0,2π], v∈[0,1], etc.),
///     the step is clipped to end exactly at the boundary, and the curve is
///     marked as terminal.  rcad uses `bounds_check` for this purpose, but
///     does not clip the step — it breaks entirely.
///
///   L600-800 — Inflection detection:
///     When `tangent_next · tangent_prev < 0`, the curve direction reversed.
///     OCCT halves the step at inflection points and re-evaluates.
///     rcad detects this via `test_deflection` returning `DeflectionStatus::Inflection`
///     and reduces the step accordingly.
///
///   L800+ — Curve termination:
///     The march ends when:
///     - The step shrinks below `min_step_size` (OCCT "ArretSurPointPrecedent").
///     - The curve exits the surface bounds.
///     - The curve closes back to the starting point (detected by
///       `closure_tol_sq`).  OCCT uses the threshold `2 * Deflection` for
///       closure detection.
fn march_one_direction_monitored_simple(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: &impl Fn(DVec3) -> bool,
    max_oscillations: usize,
    min_step_size: f64,
    step_reduction_factor: f64,
    deflection_tol: f64,
) -> MonitoredMarchResult {
    let mut points = vec![seed];
    let mut current = seed;
    let mut oscillation_count = 0usize;
    let mut step_reduced = false;
    let mut current_step = step_size;
    let mut consecutive_oscillations = 0usize;
    let mut last_dir = DVec3::ZERO;
    // Absolute initial step for growth ceiling in deflection control.
    let initial_step_abs = step_size.abs();
    // Retry counter for deflection-based step reduction loops.
    let mut retry_count = 0usize;

    // Closure tolerance: within 2× step_size of start is considered closed.
    let closure_tol_sq = (step_size * 2.0) * (step_size * 2.0);
    // Track arc length to avoid infinite loops.
    let mut arc_len = 0.0_f64;

    for _step_idx in 0..max_steps {
        let g1 = surface_gradient(s1, current);
        let g2 = surface_gradient(s2, current);
        let tangent = g1.cross(g2);
        let t_len = tangent.length();
        if t_len < TOLERANCE_ABS {
            // Tangent surfaces — the gradient cross product is zero, so the
            // usual marching direction is undefined.  Try to recover by
            // extrapolating the last valid step direction and projecting back
            // onto the intersection.  A small offset in a known-good direction
            // often escapes the tangent region.
            if points.len() > 1 {
                let last_step = current - points[points.len() - 1];
                let step_len = last_step.length();
                if step_len > TOLERANCE_ABS {
                    let attempt_dir = last_step / step_len;
                    let next_raw = current + attempt_dir * current_step.abs();
                    let next = project_onto_intersection(s1, s2, next_raw);
                    let step_dist = (next - current).length();
                    if step_dist < current_step.abs() * 5.0 && step_dist > TOLERANCE_ABS {
                        if !bounds_check(next) {
                            break;
                        }
                        arc_len += step_dist;
                        let arc_cap = step_size * max_steps as f64;
                        if arc_len >= arc_cap {
                            break;
                        }
                        points.push(next);
                        current = next;
                        last_dir = attempt_dir;
                        continue;
                    }
                }
            }
            break;
        }
        let dir = tangent / t_len * step_size.signum();

        // Oscillation detection: direction reversal
        if last_dir.length_squared() > 0.5 {
            let alignment = dir.dot(last_dir);
            if alignment < -0.9 {
                oscillation_count += 1;
                consecutive_oscillations += 1;

                // If too many consecutive oscillations, reduce step size
                if consecutive_oscillations >= max_oscillations && current_step > min_step_size {
                    current_step *= step_reduction_factor;
                    step_reduced = true;
                    consecutive_oscillations = 0;

                    // Reset current position to last good point and continue
                    if points.len() > 1 {
                        current = points[points.len() - 1];
                    }
                    continue;
                }
            } else {
                consecutive_oscillations = 0;
            }
        }

        let next_raw = current + dir * current_step.abs();
        let next = project_onto_intersection(s1, s2, next_raw);

        // --- OCCT-aligned deflection check (chord error / fleche) ---
        // Reference: OCCT IntWalk_PWalking::TestDeflection L3407-3916
        let step_dist = (next - current).length();
        let step_dist_sq = step_dist * step_dist;

        // Compute tangent at the new point for chord error estimation.
        let next_g1 = surface_gradient(s1, next);
        let next_g2 = surface_gradient(s2, next);
        let next_tangent = next_g1.cross(next_g2);
        let next_t_len = next_tangent.length();

        if next_t_len > TOLERANCE_ABS
            && last_dir.length_squared() > TOLERANCE_ABS_SQ
            && deflection_tol > 0.0
        {
            let next_tg_unit = next_tangent / next_t_len;
            let (status, suggested_step) = test_deflection(
                step_dist_sq, step_dist,
                last_dir, next_tg_unit,
                current_step.abs(), deflection_tol,
                min_step_size, initial_step_abs,
            );

            match status {
                DeflectionStatus::StepTooLarge | DeflectionStatus::Inflection => {
                    if current_step.abs() <= min_step_size || retry_count >= 3 {
                        // Cannot reduce further — accept point anyway.
                        // OCCT L3466-3472: when step falls below resolution,
                        // return ArretSurPointPrecedent (proceed with what we have).
                    } else {
                        current_step = step_size.signum() * suggested_step;
                        retry_count += 1;
                        continue;
                    }
                }
                DeflectionStatus::ConfusedPoint => {
                    // Nearly coincident — increase step and retry.
                    current_step = step_size.signum() * suggested_step;
                    retry_count += 1;
                    continue;
                }
                DeflectionStatus::StepTooSmall | DeflectionStatus::Ok => {
                    // Accept point, adjust step for next iteration.
                    current_step = step_size.signum() * suggested_step;
                    retry_count = 0;
                }
            }
        } else {
            retry_count = 0;
        }
        // --- End deflection check ---

        if !bounds_check(next) {
            break;
        }

        arc_len += step_dist;

        // Check if we've returned to start (closed curve).
        if points.len() > 10 && (next - points[0]).length_squared() < closure_tol_sq {
            points.push(points[0]); // seal the loop
            break;
        }

        // Cap arc length to prevent runaway on infinite/very long open curves.
        // Use the actual max_steps rather than a hardcoded cap of 400 so that
        // callers who set max_steps higher (for large surfaces) get the full budget.
        let arc_cap = step_size * max_steps as f64;
        if arc_len >= arc_cap {
            break;
        }

        points.push(next);
        current = next;
        last_dir = dir;
    }

    MonitoredMarchResult {
        points,
        oscillation_count,
        step_reduced,
    }
}

/// Generate sample points on a cylinder surface for seed finding.
pub fn sample_cylinder(
    cyl: &CylindricalSurface,
    height_range: [f64; 2],
    n_theta: usize,
    n_h: usize,
) -> Vec<DVec3> {
    let u = if cyl.axis.x.abs() < 0.9 {
        cyl.axis.cross(DVec3::X).normalize()
    } else {
        cyl.axis.cross(DVec3::Y).normalize()
    };
    let v = cyl.axis.cross(u);

    let mut points = Vec::with_capacity(n_theta * n_h);
    for ih in 0..n_h {
        let h = height_range[0]
            + (height_range[1] - height_range[0]) * ih as f64 / (n_h - 1).max(1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = cyl.origin + cyl.axis * h + (u * theta.cos() + v * theta.sin()) * cyl.radius;
            points.push(p);
        }
    }
    points
}

/// Generate sample points on a sphere surface for seed finding.
pub fn sample_sphere(sphere: &SphericalSurface, n_theta: usize, n_phi: usize) -> Vec<DVec3> {
    let u = if sphere.axis.x.abs() < 0.9 {
        sphere.axis.cross(DVec3::X).normalize()
    } else {
        sphere.axis.cross(DVec3::Y).normalize()
    };
    let v = sphere.axis.cross(u);

    let mut points = Vec::with_capacity(n_theta * n_phi);
    for ip in 0..n_phi {
        let phi = std::f64::consts::PI * ip as f64 / (n_phi - 1).max(1) as f64;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = sphere.center
                + sphere.radius
                    * (sphere.axis * phi.cos() + (u * theta.cos() + v * theta.sin()) * phi.sin());
            points.push(p);
        }
    }
    points
}

/// Generate sample points on a torus surface for seed finding.
pub fn sample_torus(torus: &ToroidalSurface, n_u: usize, n_v: usize) -> Vec<DVec3> {
    let u_dir = if torus.axis.x.abs() < 0.9 {
        torus.axis.cross(DVec3::X).normalize()
    } else {
        torus.axis.cross(DVec3::Y).normalize()
    };
    let v_dir = torus.axis.cross(u_dir);

    let mut points = Vec::with_capacity(n_u * n_v);
    for iu in 0..n_u {
        let u = 2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
        let cu = u.cos();
        let su = u.sin();
        let ring_center = torus.center + (u_dir * cu + v_dir * su) * torus.major_radius;
        let ring_outward = u_dir * cu + v_dir * su;

        for iv in 0..n_v {
            let v = 2.0 * std::f64::consts::PI * iv as f64 / n_v as f64;
            let p =
                ring_center + (ring_outward * v.cos() + torus.axis * v.sin()) * torus.minor_radius;
            points.push(p);
        }
    }
    points
}


